#include <cmath>
#include <fstream>
#include <iostream>

#include <mpi.h>

#include "grid.hh"

void Grid::createMatrixElementType() {
    const int count = 3;
    int lengths[count] = {1, 1, 1};
    MPI_Aint displacements[count];
    MPI_Datatype types[count] = {MPI_INT, MPI_INT, MPI_DOUBLE};

    displacements[0] = offsetof(MatrixElement, row);
    displacements[1] = offsetof(MatrixElement, col);
    displacements[2] = offsetof(MatrixElement, val);

    MPI_Type_create_struct(count, lengths, displacements, types, &MPI_MATRIX_ELEMENT);
    MPI_Type_commit(&MPI_MATRIX_ELEMENT);
}

int Grid::getElementOwnerRank(int i, int j, int idx) {
    if (idx == 1) {
        std::swap(i, j);
    }

    int colsPerStrip = matrixSize / strips;
    int extraCols = matrixSize % strips;

    int strip = j < (colsPerStrip + 1) * extraCols
                ? j / (colsPerStrip + 1)
                : extraCols + (j - (colsPerStrip + 1) * extraCols) / colsPerStrip;
    int layer = strip % layers;

    int layerCol = strip / layers;

    int rowsPerProc = matrixSize / sqrtLayerNumProcesses;
    int extraRows = matrixSize % sqrtLayerNumProcesses;

    int layerRow = i < (rowsPerProc + 1) * extraRows
                   ? i / (rowsPerProc + 1)
                   : extraRows + (i - (rowsPerProc + 1) * extraRows) / rowsPerProc;

    if (idx == 1) {
        std::swap(layerRow, layerCol);
    }

    return layer * layerNumProcesses + layerRow * sqrtLayerNumProcesses + layerCol;
}

int Grid::getFirstRowOfProcessA() {
    int rowsPerProc = matrixSize / sqrtLayerNumProcesses;
    int extraRows = matrixSize % sqrtLayerNumProcesses;
    return procRow * rowsPerProc + std::min(procRow, extraRows);
}

int Grid::getFirstColOfProcessA() {
    int procStrip = layers * procCol + procLayer;
    int colsPerStrip = matrixSize / strips;
    int extraCols = matrixSize % strips;
    return procStrip * colsPerStrip + std::min(procStrip, extraCols);
}

int Grid::getFirstRowOfProcessB() {
    int procStrip = layers * procRow + procLayer;
    int rowsPerStrip = matrixSize / strips;
    int extraRows = matrixSize % strips;
    return procStrip * rowsPerStrip + std::min(procStrip, extraRows);
}

int Grid::getFirstColOfProcessB() {
    int colsPerProc = matrixSize / sqrtLayerNumProcesses;
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

void Grid::readMatrix(const std::string& file, Matrix &matrix, int idx) {
    std::ifstream matrixFile(file);

    if (!matrixFile.is_open()) {
        std::cerr << "Error opening a file\n";
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    int numRows, numCols, numNonZero, maxRowNonZero;
    matrixFile >> numRows >> numCols >> numNonZero >> maxRowNonZero;
    matrixSize = numRows;

    std::vector<MatrixElement> elements;

    if (myRank == 0) {
        std::vector<int> counts(numProcesses);
        std::vector<double> values(numNonZero);
        std::vector<int> columns(numNonZero);
        std::vector<int> rowIndex(numRows + 1);

        for (int i = 0; i < numNonZero; i++) {
            matrixFile >> values[i];
        }
        for (int i = 0; i < numNonZero; i++) {
            matrixFile >> columns[i];
        }
        for (int i = 0; i <= numRows; i++) {
            matrixFile >> rowIndex[i];
        }

        for (int row = 0; row < numRows; row++) {
            for (int j = rowIndex[row]; j < rowIndex[row + 1]; j++) {
                int col = columns[j];
                int rank = getElementOwnerRank(row, col, idx);
                counts[rank]++;
            }
        }

        for (int p = 1; p < numProcesses; p++) {
            MPI_Send(&counts[p], 1, MPI_INT, p, 0, MPI_COMM_WORLD);
        }

        for (int row = 0; row < numRows; row++) {
            for (int j = rowIndex[row]; j < rowIndex[row + 1]; j++) {
                int col = columns[j];
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
        int count;
        MPI_Recv(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < count; i++) {
            MatrixElement elem;
            MPI_Recv(&elem, 1, MPI_MATRIX_ELEMENT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            elements.push_back(elem);
        }
    }

    matrix = Matrix(elements);
    if (idx == 0) {
        matrix.scaleElements(-getFirstRowOfProcessA(), -getFirstColOfProcessA());
    } else {
        matrix.scaleElements(-getFirstRowOfProcessB(), -getFirstColOfProcessB());
    }

    matrixFile.close();
}

void Grid::readMatrices(const std::string& fileA, const std::string& fileB) {
    readMatrix(fileA, matrixA, 0);
    readMatrix(fileB, matrixB, 1);
}

Matrix Grid::broadcastRow(Matrix matrix, int stage) {
    MPI_Comm row_comm;
    MPI_Comm_split(
        MPI_COMM_WORLD,
        procLayer * sqrtLayerNumProcesses + procRow,
        myRank,
        &row_comm
    );

    int root = procLayer * layerNumProcesses + procRow * sqrtLayerNumProcesses + stage;

    int numElements = matrixA.getElements().size();
    MPI_Request req;
    MPI_Ibcast(
        &numElements,
        1,
        MPI_INT,
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

    int root = procLayer * layerNumProcesses + stage * sqrtLayerNumProcesses + procCol;

    int numElements = matrixB.getElements().size();
    MPI_Request req;
    MPI_Ibcast(
        &numElements,
        1,
        MPI_INT,
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
    int q = matrixSize / sqrtLayerNumProcesses;
    int r = matrixSize % sqrtLayerNumProcesses;
    int cols = q + (procCol < r);
    int colsPerLayer = cols / layers;
    int extraCols = cols % layers;

    std::vector<std::vector<MatrixElement>> res(layers);

    for (auto elem: matrixC.getElements()) {
        int which;
        if (elem.col < (colsPerLayer + 1) * extraCols) {
            which = elem.col / (colsPerLayer + 1);
        } else {
            which = extraCols + (elem.col - (colsPerLayer + 1) * extraCols) / colsPerLayer;
        }
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

    int totalSendElements = 0;
    for (auto count: sendCounts) {
        totalSendElements += count;
    }
    std::vector<MatrixElement> sendBuffer;
    sendBuffer.reserve(totalSendElements);
    for (auto elements: matricesD) {
        sendBuffer.insert(sendBuffer.end(), elements.begin(), elements.end());
    }

    int totalRecvElements = 0;
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

    int q = matrixSize / sqrtLayerNumProcesses;
    int r = matrixSize % sqrtLayerNumProcesses;
    int rowOffset = procRow * q + std::min(procRow, r);
    int colOffset = procCol * q + std::min(procCol, r);
    matrixC.scaleElements(rowOffset, colOffset);

    MPI_Comm_free(&fiber_comm);
}

void Grid::SUMMA_2D() {
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
}

void Grid::SUMMA_3D() {
    SUMMA_2D();
    std::vector<std::vector<MatrixElement>> matricesD = colSplit();
    allToAll(matricesD);
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

        allElements.resize(displacements[numProcesses - 1] + recvCounts[numProcesses - 1]);
    }

    MPI_Gatherv(
        matrixC.getElements().data(), numElements, MPI_MATRIX_ELEMENT,
        allElements.data(), recvCounts.data(), displacements.data(), MPI_MATRIX_ELEMENT,
        0, MPI_COMM_WORLD
    );

    if (myRank == 0) {
        Matrix m(allElements);
        m.sort();
        int idx = 0;

        std::cout << matrixSize << " " << matrixSize << "\n";

        for (int i = 0; i < matrixSize; i++) {
            for (int j = 0; j < matrixSize; j++) {
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
