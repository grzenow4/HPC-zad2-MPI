#include <cmath>
#include <fstream>
#include <iostream>

#include <mpi.h>

#include "grid.hh"

static int getElementOwnerRank(int numProcesses, int n, int i, int j) {
    int sqrtNumProcesses = std::sqrt(numProcesses);
    int q = n / sqrtNumProcesses;
    int r = n % sqrtNumProcesses;

    int procRow = (i < (q + 1) * r) ? i / (q + 1) : r + (i - (q + 1) * r) / q;
    int procCol = (j < (q + 1) * r) ? j / (q + 1) : r + (j - (q + 1) * r) / q;

    return procRow * sqrtNumProcesses + procCol;
}

Grid::Grid(int numProcesses, int myRank) : numProcesses(numProcesses), myRank(myRank) {
    sqrtNumProcesses = std::sqrt(numProcesses);
    procRow = myRank / sqrtNumProcesses;
    procCol = myRank % sqrtNumProcesses;
}

void Grid::readMatrix(const std::string& file, Matrix &matrix) {
    std::ifstream matrixFile(file);

    if (!matrixFile.is_open()) {
        std::cerr << "Error opening a file\n";
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    int numRows, numCols, numNonZero, maxRowNonZero;
    matrixFile >> numRows >> numCols >> numNonZero >> maxRowNonZero;

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
            int rank = getElementOwnerRank(numProcesses, numRows, row, col);

            if (myRank == 0) {
                if (rank == 0) {
                    elements.push_back({row, col, val});
                } else {
                    MPI_Send(&row, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
                    MPI_Send(&col, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
                    MPI_Send(&val, 1, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD);
                }
            } else if (myRank == rank) {
                int r, c;
                double v;
                MPI_Recv(&r, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&c, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&v, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                elements.push_back({r, c, v});
            }
        }
    }

    matrix = Matrix(elements);
    matrix.sort();
}

void Grid::readMatrices(const std::string& fileA, const std::string& fileB) {
    readMatrix(fileA, matrixA);
    readMatrix(fileB, matrixB);
}

void Grid::SUMMA_2D() {
    int stages = sqrtNumProcesses;
    std::vector<Matrix> partialResults(stages);
    for (int i = 0; i < stages; i++) {
        // A_recv <- BCAST(A, P_2D(i, s), P_2D(s, j));
        // B_recv <- BCAST(B, P_2D(i, s), P_2D(s, j));
        // C[s]   <- LocalMultiply(A_recv, B_recv);
    }
    // Todo: merge
    // C <- Merge(C[0..(stages-1)]);
}

// For now it prints matrices A and B, but later it should print the matrix C.
void Grid::printMatrix() {
    std::string msg = "P2D(" + std::to_string(procRow) + ", " + std::to_string(procCol) + ")\nA = [";
    for (auto el: matrixA.elements) {
        msg += "{" + std::to_string(el.row) + ", " + std::to_string(el.col) + ", " + std::to_string(el.val) + "} ";
    }
    msg += "]\nB = [";
    for (auto el: matrixB.elements) {
        msg += "{" + std::to_string(el.row) + ", " + std::to_string(el.col) + ", " + std::to_string(el.val) + "} ";
    }
    msg += "]\n";
    std::cerr << msg;
}
