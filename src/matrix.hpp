#pragma once

#include <iostream>
#include <vector>
#include <stdexcept>
#include <cstdint>
#include <algorithm>

typedef uint64_t u64;

// Forward declarations
struct Matrix;
struct Vector;
Matrix multiply(const Matrix& A, const Matrix& B);
Vector multiply(const Matrix& A, const Vector& x);

struct Vector {
    std::vector<u64> data;

    Vector(size_t n) : data(n, 0) {}
    Vector(const std::vector<u64>& x): data(x) {}

    u64& operator[](size_t i) {
        return data[i];
    }

    const u64& operator[](size_t i) const {
        return data[i];
    }

    // Addition operator for vectors
    Vector operator+(const Vector& other) const {
        if (data.size() != other.data.size()) {
            throw std::invalid_argument("Vector dimensions must match for addition.");
        }
        Vector result(data.size());
        for (size_t i = 0; i < data.size(); ++i) {
            result[i] = data[i] + other[i];
        }
        return result;
    }

    // Subtraction operator for vectors
    Vector operator-(const Vector& other) const {
        if (data.size() != other.data.size()) {
            throw std::invalid_argument("Vector dimensions must match for subtraction.");
        }
        Vector result(data.size());
        for (size_t i = 0; i < data.size(); ++i) {
            result[i] = data[i] - other[i];
        }
        return result;
    }

    // In-place addition for vectors (+=)
    Vector& operator+=(const Vector& other) {
        if (data.size() != other.data.size()) {
            throw std::invalid_argument("Vector dimensions must match for addition.");
        }
        for (size_t i = 0; i < data.size(); ++i) {
            data[i] += other[i];
        }
        return *this;
    }

    // In-place subtraction for vectors (-=)
    Vector& operator-=(const Vector& other) {
        if (data.size() != other.data.size()) {
            throw std::invalid_argument("Vector dimensions must match for subtraction.");
        }
        for (size_t i = 0; i < data.size(); ++i) {
            data[i] -= other[i];
        }
        return *this;
    }

    // Scalar multiplication for vectors
    Vector operator*(u64 scalar) const {
        Vector result(data.size());
        for (size_t i = 0; i < data.size(); ++i) {
            result[i] = data[i] * scalar;
        }
        return result;
    }

    // Scalar multiplication from left side (u64 * Vector)
    friend Vector operator*(u64 scalar, const Vector& vec) {
        return vec * scalar;
    }

    // In-place scalar multiplication for vectors (*=)
    Vector& operator*=(u64 scalar) {
        for (size_t i = 0; i < data.size(); ++i) {
            data[i] *= scalar;
        }
        return *this;
    }

    const size_t size() const {
        return data.size();
    }

    void print() const {
        for (auto val : data) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
};

struct Matrix {
    std::vector<std::vector<u64>> data;
    size_t rows, cols;

    Matrix(size_t r, size_t c) : rows(r), cols(c), data(r, std::vector<u64>(c, 0)) {}
    Matrix(const std::vector<std::vector<u64>>& A) : data(A), rows(A.size()), cols(A[0].size()) {} // A must be non-empty 

    std::vector<u64>& operator[](size_t i) {
        return data[i];
    }

    const std::vector<u64>& operator[](size_t i) const {
        return data[i];
    }

    // Addition operator for matrices
    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Matrix dimensions must match for addition.");
        }
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result[i][j] = data[i][j] + other[i][j];
            }
        }
        return result;
    }

    // Subtraction operator for matrices
    Matrix operator-(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Matrix dimensions must match for subtraction.");
        }
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result[i][j] = data[i][j] - other[i][j];
            }
        }
        return result;
    }

    // In-place addition for matrices (+=)
    Matrix& operator+=(const Matrix& other) {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Matrix dimensions must match for addition.");
        }
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                data[i][j] += other[i][j];
            }
        }
        return *this;
    }

    // In-place subtraction for matrices (-=)
    Matrix& operator-=(const Matrix& other) {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Matrix dimensions must match for subtraction.");
        }
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                data[i][j] -= other[i][j];
            }
        }
        return *this;
    }

    // Scalar multiplication for matrices (Matrix * u64)
    Matrix operator*(u64 scalar) const {
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result[i][j] = data[i][j] * scalar;
            }
        }
        return result;
    }

    // Scalar multiplication from left side (u64 * Matrix)
    friend Matrix operator*(u64 scalar, const Matrix& mat) {
        return mat * scalar;
    }

    // In-place scalar multiplication for matrices (*=)
    Matrix& operator*=(u64 scalar) {
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                data[i][j] *= scalar;
            }
        }
        return *this;
    }

    // Matrix multiplication (*)
    Matrix operator*(const Matrix& other) const {
        return multiply(*this, other);
    }

    // In-place multiplication for matrices (*=)
    Matrix& operator*=(const Matrix& other) {
        *this = (*this) * other;
        return *this;
    }

    // Matrix-vector multiplication (Matrix * Vector)
    Vector operator*(const Vector& vec) const {
        return multiply(*this, vec);
    }

    const size_t size() const {
        return data.size();
    }

    void print() const {
        for (const auto& row : data) {
            for (auto val : row) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    }
};

// Function to multiply two matrices
Matrix multiply(const Matrix& A, const Matrix& B) {
    if (A.cols != B.rows) {
        throw std::invalid_argument("Matrix dimensions do not match for multiplication.");
    }

    Matrix result(A.rows, B.cols);
    for (size_t i = 0; i < A.rows; ++i) {
        for (size_t j = 0; j < B.cols; ++j) {
            u64 sum = 0;
            for (size_t k = 0; k < A.cols; ++k) {
                sum += A[i][k] * B[k][j];
            }
            result[i][j] = sum;
        }
    }
    return result;
}

// Function to multiply a matrix and a vector
Vector multiply(const Matrix& A, const Vector& x) {
    if (A.cols != x.data.size()) {
        throw std::invalid_argument("Matrix and vector dimensions do not match for multiplication.");
    }

    Vector result(A.rows);
    for (size_t i = 0; i < A.rows; ++i) {
        u64 sum = 0;
        for (size_t j = 0; j < A.cols; ++j) {
            sum += A[i][j] * x[j];
        }
        result[i] = sum;
    }
    return result;
}
