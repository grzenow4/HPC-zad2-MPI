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
public: // private:
    std::vector<MatrixElement> elements;

    void sort();

    Matrix transpose();

public:
    Matrix() = default;

    Matrix(std::vector<MatrixElement> elems);

    Matrix add(Matrix b);

    Matrix multiply(Matrix b);
};
