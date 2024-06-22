#pragma once

#include "matrix.hh"

class Grid {
private:
    int numProcesses;
    int myRank;
    int layers;
    int strips;
    int layerNumProcesses;
    int sqrtLayerNumProcesses;
    int procLayer;
    int layerRank;
    int procRow;
    int procCol;
    int matrixSize;
    Matrix matrixA;
    Matrix matrixB;
    Matrix matrixC;
    MPI_Datatype MPI_MATRIX_ELEMENT;

    void createMatrixElementType();

    int getElementOwnerRank(int i, int j, int which);
    int getFirstRowOfProcessA();
    int getFirstColOfProcessA();
    int getFirstRowOfProcessB();
    int getFirstColOfProcessB();

    void readMatrix(const std::string& file, Matrix &matrix, int idx);

    Matrix broadcastRow(Matrix matrix, int stage);
    Matrix broadcastCol(Matrix matrix, int stage);

    std::vector<std::vector<MatrixElement>> colSplit();
    void allToAll(std::vector<std::vector<MatrixElement>> matricesD);

public:
    Grid(int numProcesses, int myRank, int layers);

    ~Grid();

    void readMatrices(const std::string& fileA, const std::string& fileB);

    void SUMMA_2D();
    void SUMMA_3D();

    void printGValue(double gValue);
    void printMatrix();
};
