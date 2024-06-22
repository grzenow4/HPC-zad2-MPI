#pragma once

#include <vector>

struct MatrixElement {
    int row;
    int col;
    double val;

    bool operator<(const MatrixElement& other) {
        return row < other.row || (row == other.row && col < other.col);
    }
};

class Matrix {
public:
    std::vector<MatrixElement> elements;

    Matrix transpose();

public:
    Matrix() = default;

    Matrix(std::vector<MatrixElement> elems);

    std::vector<MatrixElement> getElements();

    void sort();

    void scaleElements(int r, int c);

    int countGreater(double g);

    Matrix add(Matrix b);

    Matrix multiply(Matrix b);
};
