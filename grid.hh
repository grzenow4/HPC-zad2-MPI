#pragma once

#include "matrix.hh"

class Grid {
private:
    int numProcesses;          // how many processes
    int myRank;                // rank of a process
    int layers;                // how many layers
    int strips;                // layers * sqrt(numProcesses / layers)
    int layerNumProcesses;     // numProcesses / layers
    int sqrtLayerNumProcesses; // sqrt(numProcesses / layers)
    int procLayer;             // layer of a process
    int layerRank;             // rank of a process in a layer
    int procRow;               // row of a process in a layer grid
    int procCol;               // col of a process in a layer grid
    uint32_t matrixSize;
    Matrix matrixA;
    Matrix matrixB;
    Matrix matrixC;
    MPI_Datatype MPI_MATRIX_ELEMENT;

    /*
     * Creates a helper MPI_MATRIX_ELEMENT type for sending matrix elements.
     */
    void createMatrixElementType();

    /*
     * Returns the rank of a process that receives an (i,j) matrix element.
     * If idx == 0 look at the matrixA, otherwise at the matrixB.
     */
    int getElementOwnerRank(uint32_t i, uint32_t j, int idx);

    /*
     * Functions that return top left corner indices
     * of a submatrices A and B of a process.
     */
    uint32_t getFirstRowOfProcessA();
    uint32_t getFirstColOfProcessA();
    uint32_t getFirstRowOfProcessB();
    uint32_t getFirstColOfProcessB();

    /*
     * Reads a matrix from the input file.
     * If idx == 0 save the result in the matrixA, otherwise in the matrixB.
     */
    void readMatrix(const std::string& file, int idx);

    /*
     * Broadcasts row/col in the 2D-SUMMA algorithm.
     * Detailed description in paper.
     */
    Matrix broadcastRow(Matrix matrix, int stage);
    Matrix broadcastCol(Matrix matrix, int stage);

    /*
     * Splits matrixC (D in paper) into l parts to perform AllToall
     * operation later. Detailed description in paper.
     */
    std::vector<std::vector<MatrixElement>> colSplit();

    /*
     * Sends and receives parts of a matrixC (D in paper) to/from other
     * processes in the same fiber. Detailed description in paper.
     */
    void allToAll(std::vector<std::vector<MatrixElement>> matricesD);

    /*
     * After completion of the SUMMA algorithm, each process has to scale his
     * elements by his top left corner, adding it to all indices.
     */
    void scaleMatrices();

public:
    Grid(int numProcesses, int myRank, int layers);

    /*
     * Frees allocated MPI_MATRIX_ELEMENT at the beginning.
     */
    ~Grid();

    /*
     * Reads input matrices.
     */
    void readMatrices(const std::string& fileA, const std::string& fileB);

    /*
     * Implementation of SUMMA algorithms as described in paper.
     * Detailed description in paper.
     * 
     * Note: if only2D is set, then `scaleMatrices` is invoked a bit earlier.
     * We need that flag in case of `-t 2D` option is set.
     */
    void SUMMA_2D(bool only2D);
    void SUMMA_3D();

    /*
     * Prints gValue and the whole result matrix.
     */
    void printGValue(double gValue);
    void printMatrix();
};
