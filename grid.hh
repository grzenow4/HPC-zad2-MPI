#pragma once

#include "matrix.hh"

class Grid {
private:
    int numProcesses;
    int myRank;
    int sqrtNumProcesses;
    int procRow;
    int procCol;
    Matrix matrixA;
    Matrix matrixB;
    Matrix matrixC;

public:
    Grid(int numProcesses, int myRank);

    void readMatrix(const std::string& file, Matrix &matrix);

    void readMatrices(const std::string& fileA, const std::string& fileB);

    void SUMMA_2D();

    void printMatrix();
};
