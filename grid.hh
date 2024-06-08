#pragma once

#include "matrix.hh"

class Grid {
private:
    int numProcesses;
    int myRank;
    int sqrtNumProcesses;
    int procRow;
    int procCol;
    int matrixSize;
    Matrix matrixA;
    Matrix matrixB;
    Matrix matrixC;
    std::vector<std::vector<int>> noElements;
    MPI_Datatype MPI_MATRIX_ELEMENT;

    void createMatrixElementType();

    int getElementOwnerRank(int i, int j);
    int getFirstRowOfProcess();
    int getFirstColOfProcess();

    void readMatrix(const std::string& file, Matrix &matrix, int idx);

    Matrix broadcastRow(Matrix matrix, int stage);
    Matrix broadcastCol(Matrix matrix, int stage);

public:
    Grid(int numProcesses, int myRank);

    ~Grid();

    void readMatrices(const std::string& fileA, const std::string& fileB);

    void SUMMA_2D();

    void printMatrix();
};
