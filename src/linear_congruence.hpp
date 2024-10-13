#pragma once

#include <iostream>
#include <vector>
#include <stdexcept>
#include <cstdint>
#include <algorithm>
#include <tuple>

#include "matrix.hpp"

typedef uint64_t u64;
typedef __int128_t i128;

std::tuple<i128, i128, i128> extended_gcd(i128 a, i128 b) {
    if (a == 0) {
        return std::make_tuple(b, (i128)0, (i128)1); // We must do this in order to avoid GCC throwing errors
    }

    auto [gcd, x1, y1] = extended_gcd(b % a, a);
    i128 x = y1 - (b / a) * x1;  
    i128 y = x1;
    return {gcd, x, y};
}

/*
Solves a single linear congruence: cx = b (mod r)
*/
i128 solve_congruence(i128 c, i128 b, i128 r) {
    auto [gcd, x, _ignore_1] = extended_gcd(c, r);

    if (b % gcd != 0) {
        throw std::invalid_argument("No solution exists for the congruence.");
    }

    c /= gcd;
    b /= gcd;
    r /= gcd;

    auto [gcd_inv, inv_c, _ignore_2] = extended_gcd(c, r);

    if (gcd_inv != 1) {
        throw std::invalid_argument("Modular inverse doesn't exist.");
    }

    auto x0 = (inv_c * b) % r;
    return x0;
}

/*
"Diagonalizes" matrix A.
Generates matrices S and T such that D = SAT is diagonal.
Returns S and T. 
*/
std::tuple<Matrix, Matrix> diagonalize(Matrix A) {
    size_t n_rows = A.rows, n_cols = A.cols;
    Matrix S(n_rows, n_rows), T(n_cols, n_cols);

    // Initialize S and T to identity matrices
    for (size_t i = 0; i < n_rows; ++i) S[i][i] = 1;
    for (size_t i = 0; i < n_cols; ++i) T[i][i] = 1;

    for (size_t i = 0; i < std::min(n_rows, n_cols); ++i) {
        // std::cerr << "working on row/col " << i << std::endl;
        while (true) {
            bool col_nonzero = false;
            for (size_t k = i + 1; k < n_rows; ++k) {
                if (A[k][i] != 0) {
                    col_nonzero = true;
                    break;
                }
            }

            if (col_nonzero) {
                size_t pivot_row = i;
                for (size_t k = i + 1; k < n_rows; ++k) {
                    if ((A[k][i] < A[pivot_row][i] || A[pivot_row][i] == 0) && A[k][i] != 0) {
                        pivot_row = k;
                    }
                }
                std::swap(A[i], A[pivot_row]);
                std::swap(S[i], S[pivot_row]);

                for (size_t k = i + 1; k < n_rows; ++k) {
                    if (A[k][i] != 0) {
                        u64 multiplier = -(A[k][i] / A[i][i]);
                        for (size_t j = 0; j < n_cols; ++j) {
                            A[k][j] = (A[k][j] + multiplier * A[i][j]);
                        }
                        for (size_t j = 0; j < n_rows; ++j) {
                            S[k][j] = (S[k][j] + multiplier * S[i][j]);
                        }
                    }
                }
            } else {
                bool row_nonzero = false;
                for (size_t k = i + 1; k < n_cols; ++k) {
                    if (A[i][k] != 0) {
                        row_nonzero = true;
                        break;
                    }
                }

                if (!row_nonzero) break;

                size_t pivot_col = i;
                for (size_t k = i + 1; k < n_cols; ++k) {
                    if ((A[i][k] < A[i][pivot_col] || A[i][pivot_col] == 0) && A[i][k] != 0) {
                        pivot_col = k;
                    }
                }

                for (size_t j = 0; j < std::max(n_rows, n_cols); ++j) {
                    if (j < n_rows) {
                        std::swap(A[j][i], A[j][pivot_col]);
                    }
                    std::swap(T[j][i], T[j][pivot_col]);
                }

                for (size_t k = i + 1; k < n_cols; ++k) {
                    if (A[i][k] != 0) {
                        u64 multiplier = -(A[i][k] / A[i][i]);
                        for (size_t j = 0; j < n_rows; ++j) {
                            A[j][k] = (A[j][k] + multiplier * A[j][i]);
                        }
                        for (size_t j = 0; j < n_cols; ++j) {
                            T[j][k] = (T[j][k] + multiplier * T[j][i]);
                        }
                    }
                }
            }
        }
    }

    // std::cerr << "Final matrix: " << std::endl;
    // A.print();

    return {S, T};
}

/*
Solves Dx' = Sb, where D is a diagonal matrix and Sb is a vector.
This is simply a series of single linear congruences.
*/
Vector solve_diagonal(const Matrix& D, const Vector& Sb) {
    size_t n = Sb.size();
    Vector x_prime(D.cols);

    for (size_t i = 0; i < n; ++i) {
        if (D[i][i] != 0) {
            // Solve D[i][i] * x'_i = S b[i] (mod 2^64)
            x_prime[i] = solve_congruence(D[i][i], Sb[i], __int128_t(1) << 64);
        } else {
            // If D[i][i] is 0, check if Sb[i] is 0
            if (Sb[i] == 0) {
                x_prime[i] = 0;  // Arbitrary solution, choose x'_i = 0
            } else {
                throw std::invalid_argument("No solution exists for the equation 0 * x' = S b[i] mod 2^64");
            }
        }
    }

    return x_prime;
}

// Solves a system of linear congruences.
Vector solve_linear_system(Matrix A, Vector b) {
    // Step 1: "Diagonalize" A into D = S A T
    auto [S, T] = diagonalize(A);
    Matrix D = S * A * T;  // D = S * A * T

    /*
    std::cerr << "Diagonalization done. \n";
    std::cerr << "S: " << S.rows << ' ' << S.cols << "\n";
    S.print();
    std::cerr << "T: " << T.rows << ' ' << T.cols << "\n";
    T.print();
    std::cerr << "D: " << D.rows << ' ' << D.cols << "\n";
    D.print();
    */

    // Step 2: Solve D x' = S b
    Vector Sb = S * b;
    Vector x_prime = solve_diagonal(D, Sb);

    // std::cerr << "Done solving diagonal" << std::endl;

    // Step 3: Compute x = T x'
    Vector x = T * x_prime;

    return x;
}