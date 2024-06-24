#include <algorithm>

#include "matrix.hh"

Matrix::Matrix(std::vector<MatrixElement> elems) : elements(elems) {}

std::vector<MatrixElement> Matrix::getElements() {
    return elements;
}

void Matrix::sort() {
    std::sort(elements.begin(), elements.end());
}

void Matrix::scaleElements(uint32_t r, uint32_t c, bool sign) {
    for (auto &elem: elements) {
        if (sign) {
            elem.row += r;
            elem.col += c;
        } else {
            elem.row -= r;
            elem.col -= c;
        }
    }
}

u_int32_t Matrix::countGreater(double g) {
    u_int32_t res = 0;
    for (auto elem: elements) {
        if (elem.val > g) {
            res++;
        }
    }
    return res;
}

Matrix Matrix::transpose() {
    std::vector<MatrixElement> res;
    for (auto elem: elements) {
        res.push_back({elem.col, elem.row, elem.val});
    }
    return Matrix(res);
}

Matrix Matrix::add(Matrix b) {
    uint32_t i = 0, j = 0;
    uint32_t lenA = elements.size(), lenB = b.elements.size();
    std::vector<MatrixElement> res;

    while (i < lenA && j < lenB) {
        if (elements[i] < b.elements[j]) {
            res.push_back(elements[i]);
            i++;
        } else if (b.elements[j] < elements[i]) {
            res.push_back(b.elements[j]);
            j++;
        } else {
            res.push_back({elements[i].row,
                           elements[i].col,
                           elements[i].val + b.elements[j].val});
            i++;
            j++;
        }
    }

    while (i < lenA) {
        res.push_back(elements[i]);
        i++;
    }

    while (j < lenB) {
        res.push_back(b.elements[j]);
        j++;
    }

    return Matrix(res);
}

Matrix Matrix::multiply(Matrix b) {
    b = b.transpose();
    b.sort();

    uint32_t i, j;
    uint32_t lenA = elements.size(), lenB = b.elements.size();
    std::vector<MatrixElement> res;

    for (i = 0; i < lenA;) {
        uint32_t r = elements[i].row; // current row

        for (j = 0; j < lenB;) {
            uint32_t c = b.elements[j].row; // current column

            uint32_t tempA = i, tempB = j;
            double sum = 0;

            while (tempA < lenA && elements[tempA].row == r &&
                   tempB < lenB && b.elements[tempB].row == c) {
                if (elements[tempA].col < b.elements[tempB].col) {
                    tempA++;
                } else if (elements[tempA].col > b.elements[tempB].col) {
                    tempB++;
                } else {
                    sum += elements[tempA++].val * b.elements[tempB++].val;
                }
            }

            if (sum != 0) {
                res.push_back({r, c, sum});
            }

            while (j < lenB && b.elements[j].row == c) {
                j++;
            }
        }

        while (i < lenA && elements[i].row == r) {
            i++;
        }
    }

    return Matrix(res);
}
