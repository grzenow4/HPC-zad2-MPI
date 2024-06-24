#include <cmath>
#include <fstream>
#include <iostream>

#include <mpi.h>

#include "grid.hh"

void Grid::createMatrixElementType() {
    const int count = 3;
    int lengths[count] = {1, 1, 1};
    MPI_Aint displacements[count];
    MPI_Datatype types[count] = {MPI_UINT32_T, MPI_UINT32_T, MPI_DOUBLE};

    displacements[0] = offsetof(MatrixElement, row);
    displacements[1] = offsetof(MatrixElement, col);
    displacements[2] = offsetof(MatrixElement, val);

    MPI_Type_create_struct(count, lengths, displacements, types, &MPI_MATRIX_ELEMENT);
    MPI_Type_commit(&MPI_MATRIX_ELEMENT);
}

int Grid::getElementOwnerRank(uint32_t i, uint32_t j, int idx) {
    if (idx == 1) {
        std::swap(i, j);
    }

    uint32_t colsPerStrip = matrixSize / strips;
    uint32_t extraCols = matrixSize % strips;

    int strip = j < (colsPerStrip + 1) * extraCols
              ? j / (colsPerStrip + 1)
              : extraCols + (j - (colsPerStrip + 1) * extraCols) / colsPerStrip;
    int layer = strip % layers;

    int layerCol = strip / layers;

    uint32_t rowsPerProc = matrixSize / sqrtLayerNumProcesses;
    uint32_t extraRows = matrixSize % sqrtLayerNumProcesses;

    int layerRow = i < (rowsPerProc + 1) * extraRows
                 ? i / (rowsPerProc + 1)
                 : extraRows + (i - (rowsPerProc + 1) * extraRows) / rowsPerProc;

    if (idx == 1) {
        std::swap(layerRow, layerCol);
    }

    return layer * layerNumProcesses + layerRow * sqrtLayerNumProcesses + layerCol;
}

uint32_t Grid::getFirstRowOfProcessA() {
    uint32_t rowsPerProc = matrixSize / sqrtLayerNumProcesses;
    int extraRows = matrixSize % sqrtLayerNumProcesses;
    return procRow * rowsPerProc + std::min(procRow, extraRows);
}

uint32_t Grid::getFirstColOfProcessA() {
    int procStrip = layers * procCol + procLayer;
    uint32_t colsPerStrip = matrixSize / strips;
    int extraCols = matrixSize % strips;
    return procStrip * colsPerStrip + std::min(procStrip, extraCols);
}

uint32_t Grid::getFirstRowOfProcessB() {
    int procStrip = layers * procRow + procLayer;
    uint32_t rowsPerStrip = matrixSize / strips;
    int extraRows = matrixSize % strips;
    return procStrip * rowsPerStrip + std::min(procStrip, extraRows);
}

uint32_t Grid::getFirstColOfProcessB() {
    uint32_t colsPerProc = matrixSize / sqrtLayerNumProcesses;
    int extraCols = matrixSize % sqrtLayerNumProcesses;
    return procCol * colsPerProc + std::min(procCol, extraCols);
}

Grid::Grid(int numProcesses, int myRank, int layers)
    : numProcesses(numProcesses)
    , myRank(myRank)
    , layers(layers) {
    layerNumProcesses = numProcesses / layers;
    sqrtLayerNumProcesses = std::sqrt(layerNumProcesses);
    strips = sqrtLayerNumProcesses * layers;

    procLayer = myRank / layerNumProcesses;
    layerRank = myRank % layerNumProcesses;
    procRow = layerRank / sqrtLayerNumProcesses;
    procCol = layerRank % sqrtLayerNumProcesses;
    createMatrixElementType();
}

Grid::~Grid() {
    MPI_Type_free(&MPI_MATRIX_ELEMENT);
}

void Grid::readMatrix(const std::string& file, int idx) {
    std::ifstream matrixFile(file);

    if (!matrixFile.is_open()) {
        std::cerr << "Error opening a file\n";
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    uint32_t numRows, numCols, maxRowNonZero;
    uint64_t numNonZero;
    matrixFile >> numRows >> numCols >> numNonZero >> maxRowNonZero;
    matrixSize = numRows;

    std::vector<MatrixElement> elements;

    if (myRank == 0) {
        std::vector<uint32_t> counts(numProcesses);
        std::vector<double> values(numNonZero);
        std::vector<uint32_t> columns(numNonZero);
        std::vector<uint32_t> rowIndex(numRows + 1);

        for (uint64_t i = 0; i < numNonZero; i++) {
            matrixFile >> values[i];
        }
        for (uint64_t i = 0; i < numNonZero; i++) {
            matrixFile >> columns[i];
        }
        for (uint64_t i = 0; i <= numRows; i++) {
            matrixFile >> rowIndex[i];
        }

        for (uint32_t row = 0; row < numRows; row++) {
            for (uint32_t j = rowIndex[row]; j < rowIndex[row + 1]; j++) {
                uint32_t col = columns[j];
                int rank = getElementOwnerRank(row, col, idx);
                counts[rank]++;
            }
        }

        for (int p = 1; p < numProcesses; p++) {
            MPI_Send(&counts[p], 1, MPI_UINT32_T, p, 0, MPI_COMM_WORLD);
        }

        for (uint32_t row = 0; row < numRows; row++) {
            for (uint32_t j = rowIndex[row]; j < rowIndex[row + 1]; j++) {
                uint32_t col = columns[j];
                double val = values[j];
                int rank = getElementOwnerRank(row, col, idx);

                MatrixElement elem = {row, col, val};

                if (rank == 0) {
                    elements.push_back(elem);
                } else {
                    MPI_Send(&elem, 1, MPI_MATRIX_ELEMENT, rank, 0, MPI_COMM_WORLD);
                }
            }
        }
    } else {
        uint32_t count;
        MPI_Recv(&count, 1, MPI_UINT32_T, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (uint32_t i = 0; i < count; i++) {
            MatrixElement elem;
            MPI_Recv(&elem, 1, MPI_MATRIX_ELEMENT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            elements.push_back(elem);
        }
    }

    if (idx == 0) {
        matrixA = Matrix(elements);
        matrixA.scaleElements(getFirstRowOfProcessA(), getFirstColOfProcessA(), false);
    } else {
        matrixB = Matrix(elements);
        matrixB.scaleElements(getFirstRowOfProcessB(), getFirstColOfProcessB(), false);
    }

    matrixFile.close();
}

void Grid::readMatrices(const std::string& fileA, const std::string& fileB) {
    readMatrix(fileA, 0);
    readMatrix(fileB, 1);
}

Matrix Grid::broadcastRow(Matrix matrix, int stage) {
    MPI_Comm row_comm;
    MPI_Comm_split(
        MPI_COMM_WORLD,
        procLayer * sqrtLayerNumProcesses + procRow,
        myRank,
        &row_comm
    );

    uint32_t numElements = matrixA.getElements().size();
    MPI_Request req;
    MPI_Ibcast(
        &numElements,
        1,
        MPI_UINT32_T,
        stage,
        row_comm,
        &req
    );
    MPI_Wait(&req, MPI_STATUS_IGNORE);

    std::vector<MatrixElement> buffer(numElements);

    if (procCol == stage) {
        buffer = matrix.getElements();
    }

    MPI_Ibcast(
        buffer.data(),
        numElements,
        MPI_MATRIX_ELEMENT,
        stage,
        row_comm,
        &req
    );
    MPI_Wait(&req, MPI_STATUS_IGNORE);

    MPI_Comm_free(&row_comm);

    if (procCol != stage) {
        matrix = Matrix(buffer);
    }

    return matrix;
}

Matrix Grid::broadcastCol(Matrix matrix, int stage) {
    MPI_Comm col_comm;
    MPI_Comm_split(
        MPI_COMM_WORLD,
        procLayer * sqrtLayerNumProcesses + procCol,
        myRank,
        &col_comm
    );

    uint32_t numElements = matrixB.getElements().size();
    MPI_Request req;
    MPI_Ibcast(
        &numElements,
        1,
        MPI_UINT32_T,
        stage,
        col_comm,
        &req
    );
    MPI_Wait(&req, MPI_STATUS_IGNORE);

    std::vector<MatrixElement> buffer(numElements);

    if (procRow == stage) {
        buffer = matrix.getElements();
    }

    MPI_Ibcast(
        buffer.data(),
        numElements,
        MPI_MATRIX_ELEMENT,
        stage,
        col_comm,
        &req
    );
    MPI_Wait(&req, MPI_STATUS_IGNORE);

    MPI_Comm_free(&col_comm);

    if (procRow != stage) {
        matrix = Matrix(buffer);
    }

    return matrix;
}

std::vector<std::vector<MatrixElement>> Grid::colSplit() {
    uint32_t q = matrixSize / sqrtLayerNumProcesses;
    int r = matrixSize % sqrtLayerNumProcesses;
    uint32_t cols = q + (procCol < r);
    uint32_t colsPerLayer = cols / layers;
    int extraCols = cols % layers;

    std::vector<std::vector<MatrixElement>> res(layers);

    for (auto elem: matrixC.getElements()) {
        int which = elem.col < (colsPerLayer + 1) * extraCols
                  ? elem.col / (colsPerLayer + 1)
                  : extraCols + (elem.col - (colsPerLayer + 1) * extraCols) / colsPerLayer;
        res[which].push_back(elem);
    }

    return res;
}

void Grid::allToAll(std::vector<std::vector<MatrixElement>> matricesD) {
    MPI_Comm fiber_comm;
    MPI_Comm_split(MPI_COMM_WORLD, layerRank, myRank, &fiber_comm);

    std::vector<int> sendCounts(layers);
    for (int i = 0; i < layers; i++) {
        sendCounts[i] = matricesD[i].size();
    }

    std::vector<int> recvCounts(layers);
    MPI_Alltoall(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, fiber_comm);

    std::vector<int> sendDisplacements(layers, 0), recvDisplacements(layers, 0);
    for (int i = 1; i < layers; i++) {
        sendDisplacements[i] = sendDisplacements[i - 1] + sendCounts[i - 1];
        recvDisplacements[i] = recvDisplacements[i - 1] + recvCounts[i - 1];
    }

    uint32_t totalSendElements = 0;
    for (auto count: sendCounts) {
        totalSendElements += count;
    }
    std::vector<MatrixElement> sendBuffer;
    sendBuffer.reserve(totalSendElements);
    for (auto elements: matricesD) {
        sendBuffer.insert(sendBuffer.end(), elements.begin(), elements.end());
    }

    uint32_t totalRecvElements = 0;
    for (auto count: recvCounts) {
        totalRecvElements += count;
    }
    std::vector<MatrixElement> recvBuffer(totalRecvElements);

    MPI_Alltoallv(
        sendBuffer.data(), sendCounts.data(), sendDisplacements.data(), MPI_MATRIX_ELEMENT,
        recvBuffer.data(), recvCounts.data(), recvDisplacements.data(), MPI_MATRIX_ELEMENT,
        fiber_comm
    );

    std::vector<std::vector<MatrixElement>> allElements(layers);
    for (int i = 0; i < layers; ++i) {
        allElements[i].insert(
            allElements[i].end(),
            recvBuffer.begin() + recvDisplacements[i],
            recvBuffer.begin() + recvDisplacements[i] + recvCounts[i]
        );
    }

    matrixC = Matrix();
    for (auto elements: allElements) {
        matrixC = matrixC.add(elements);
    }

    MPI_Comm_free(&fiber_comm);
}

void Grid::scaleMatrices() {
    uint32_t q = matrixSize / sqrtLayerNumProcesses;
    int r = matrixSize % sqrtLayerNumProcesses;
    uint32_t rowOffset = procRow * q + std::min(procRow, r);
    uint32_t colOffset = procCol * q + std::min(procCol, r);
    matrixC.scaleElements(rowOffset, colOffset, true);
}

void Grid::SUMMA_2D(bool only2D) {
    int stages = sqrtLayerNumProcesses;
    std::vector<Matrix> partialResults(stages);
    for (int i = 0; i < stages; i++) {
        Matrix A_recv = broadcastRow(matrixA, i);
        Matrix B_recv = broadcastCol(matrixB, i);
        partialResults[i] = A_recv.multiply(B_recv);
    }
    for (auto partial: partialResults) {
        matrixC = matrixC.add(partial);
    }
    if (only2D) {
        scaleMatrices();
    }
}

void Grid::SUMMA_3D() {
    SUMMA_2D(false);
    std::vector<std::vector<MatrixElement>> matricesD = colSplit();
    allToAll(matricesD);
    scaleMatrices();
}

void Grid::printGValue(double gValue) {
    double gValGlob;
    double gValLoc = matrixC.countGreater(gValue);
    MPI_Reduce(&gValLoc, &gValGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myRank == 0) {
        std::cout << gValGlob << "\n";
    }
}

void Grid::printMatrix() {
    std::vector<MatrixElement> allElements;
    std::vector<int> recvCounts(numProcesses);
    std::vector<int> displacements(numProcesses);

    int numElements = matrixC.getElements().size();
    MPI_Gather(&numElements, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (myRank == 0) {
        displacements[0] = 0;
        for (int i = 1; i < numProcesses; i++) {
            displacements[i] = displacements[i - 1] + recvCounts[i - 1];
        }

        allElements.resize((uint32_t) displacements[numProcesses - 1] + (uint32_t) recvCounts[numProcesses - 1]);
    }

    MPI_Gatherv(
        matrixC.getElements().data(), numElements, MPI_MATRIX_ELEMENT,
        allElements.data(), recvCounts.data(), displacements.data(), MPI_MATRIX_ELEMENT,
        0, MPI_COMM_WORLD
    );

    if (myRank == 0) {
        Matrix m(allElements);
        m.sort();
        size_t idx = 0;

        std::cout << matrixSize << " " << matrixSize << "\n";

        for (uint32_t i = 0; i < matrixSize; i++) {
            for (uint32_t j = 0; j < matrixSize; j++) {
                if (idx < m.getElements().size() &&
                    m.getElements()[idx].row == i &&
                    m.getElements()[idx].col == j) {
                    std::cout << m.getElements()[idx].val << " ";
                    idx++;
                } else {
                    std::cout << "0 ";
                }
            }
            std::cout << "\n";
        }
    }
}
