#pragma once

#include <vector>

struct MatrixElement {
    int row;
    int col;
    double val;
};

class Matrix {
public: // private:
    std::vector<MatrixElement> elements;

public:
    Matrix() = default;

    Matrix(std::vector<MatrixElement> elems);

    void sort();

    Matrix transpose();

    Matrix add(Matrix other);

    Matrix multiply(Matrix other);

    void print();
};
