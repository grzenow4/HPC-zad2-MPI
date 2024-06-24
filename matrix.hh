#pragma once

#include <cstdint>
#include <vector>

struct MatrixElement {
    uint32_t row;
    uint32_t col;
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

    void scaleElements(uint32_t r, uint32_t c, bool sign);

    uint32_t countGreater(double g);

    Matrix add(Matrix b);

    Matrix multiply(Matrix b);
};
