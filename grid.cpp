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

int Grid::getElementOwnerRank(int i, int j) {
    int q = matrixSize / sqrtNumProcesses;
    int r = matrixSize % sqrtNumProcesses;

    int pRow = (i < (q + 1) * r) ? i / (q + 1) : r + (i - (q + 1) * r) / q;
    int pCol = (j < (q + 1) * r) ? j / (q + 1) : r + (j - (q + 1) * r) / q;

    return pRow * sqrtNumProcesses + pCol;
}

int Grid::getFirstRowOfProcess() {
    int q = matrixSize / sqrtNumProcesses;
    int r = matrixSize % sqrtNumProcesses;
    if (procRow < r) {
        return procRow * (q + 1);
    }
    return procRow * q + r;
}

int Grid::getFirstColOfProcess() {
    int q = matrixSize / sqrtNumProcesses;
    int r = matrixSize % sqrtNumProcesses;
    if (procCol < r) {
        return procCol * (q + 1);
    }
    return procCol * q + r;
}

Grid::Grid(int numProcesses, int myRank)
    : numProcesses(numProcesses)
    , myRank(myRank)
    , noElements(numProcesses, std::vector<int>(3, 0)) {
    sqrtNumProcesses = std::sqrt(numProcesses);
    procRow = myRank / sqrtNumProcesses;
    procCol = myRank % sqrtNumProcesses;
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

    matrixFile.close();

    std::vector<MatrixElement> elements;

    for (int row = 0; row < numRows; row++) {
        for (int j = rowIndex[row]; j < rowIndex[row + 1]; j++) {
            int col = columns[j];
            double val = values[j];
            int rank = getElementOwnerRank(row, col);

            MatrixElement elem = {row, col, val};
            noElements[rank][idx]++;

            if (myRank == 0) {
                if (rank == 0) {
                    elements.push_back(elem);
                } else {
                    MPI_Send(&elem, 1, MPI_MATRIX_ELEMENT, rank, 0, MPI_COMM_WORLD);
                }
            } else if (myRank == rank) {
                MPI_Recv(&elem, 1, MPI_MATRIX_ELEMENT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                elem.row -= getFirstRowOfProcess();
                elem.col -= getFirstColOfProcess();
                elements.push_back(elem);
            }
        }
    }

    matrix = Matrix(elements);
    matrix.sort();
}

void Grid::readMatrices(const std::string& fileA, const std::string& fileB) {
    readMatrix(fileA, matrixA, 0);
    readMatrix(fileB, matrixB, 1);
}

Matrix Grid::broadcastRow(Matrix matrix, int stage) {
    MPI_Comm row_comm;
    MPI_Comm_split(MPI_COMM_WORLD, procRow, myRank, &row_comm);

    int root = procRow * sqrtNumProcesses + stage;
    int numElements = noElements[root][0];

    std::vector<MatrixElement> buffer(numElements);

    if (procCol == stage) {
        buffer = matrix.elements;
    }

    MPI_Request req;
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
    MPI_Comm_split(MPI_COMM_WORLD, procCol, myRank, &col_comm);

    int root = stage * sqrtNumProcesses + procCol;
    int numElements = noElements[root][1];

    std::vector<MatrixElement> buffer(numElements);

    if (procRow == stage) {
        buffer = matrix.elements;
    }

    MPI_Request req;
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

void Grid::SUMMA_2D() {
    int stages = sqrtNumProcesses;
    std::vector<Matrix> partialResults(stages);
    for (int i = 0; i < stages; i++) {
        Matrix A_recv = broadcastRow(matrixA, i);
        Matrix B_recv = broadcastCol(matrixB, i);
        partialResults[i] = A_recv.multiply(B_recv);
    }
    for (auto partial: partialResults) {
        matrixC = matrixC.add(partial);
    }
    for (auto &elem: matrixC.elements) {
        elem.row += getFirstRowOfProcess();
        elem.col += getFirstColOfProcess();
    }
}

// GPT-generated, only for debug now.
void Grid::printMatrix() {
    std::vector<MatrixElement> allElements;
    std::vector<int> recvCounts(numProcesses);
    std::vector<int> displacements(numProcesses);

    // Get the number of elements to be sent by each process
    int numElements = matrixC.elements.size();
    MPI_Gather(&numElements, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (myRank == 0) {
        // Calculate displacements for gather operation
        displacements[0] = 0;
        for (int i = 1; i < numProcesses; i++) {
            displacements[i] = displacements[i - 1] + recvCounts[i - 1];
        }

        // Resize the vector to hold all elements
        allElements.resize(displacements[numProcesses - 1] + recvCounts[numProcesses - 1]);
    }

    // Gather elements from all processes to process 0
    MPI_Gatherv(matrixC.elements.data(), numElements, MPI_MATRIX_ELEMENT, allElements.data(), recvCounts.data(), displacements.data(), MPI_MATRIX_ELEMENT, 0, MPI_COMM_WORLD);

    // Print the gathered matrix in process 0
    if (myRank == 0) {
        Matrix m(allElements);
        m.sort();
        int idx = 0;

        for (int i = 0; i < matrixSize; i++) {
            for (int j = 0; j < matrixSize; j++) {
                if (idx < m.elements.size() && m.elements[idx].row == i && m.elements[idx].col == j) {
                    std::cout << m.elements[idx].val << " ";
                    idx++;
                } else {
                    std::cout << "0 ";
                }
            }
            std::cout << "\n";
        }
    }
}
