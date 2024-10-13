#include <iostream>
#include <vector>
#include <stdexcept>
#include <cstdint>
#include <algorithm>
#include <tuple>

#include "matrix.hpp"
#include "linear_congruence.hpp"

typedef uint64_t u64;

int main() {
    Matrix A(2, 2);
    A[0][0] = 3; A[0][1] = 4;
    A[1][0] = 8; A[1][1] = 2;

    Vector b(2);
    b[0] = 69;
    b[1] = 278;

    Vector solution = solve_linear_system(A, b);

    std::cout << "Solution:" << std::endl;
    solution.print();
    std::cout << "Checking A * x: " << std::endl;
    multiply(A, solution).print();

    return 0;
}
