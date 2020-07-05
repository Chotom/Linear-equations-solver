/* *********************************************************************************************************************
 * author: Tomasz Czochański
 *
 * Main program
 **********************************************************************************************************************/

#include <iostream>
#include "main_methods.h"


int main() {

    auto matrix = create_matrix(5 + E, -1, -1, SIZE_N);
    auto vector_b = create_vector_b(SIZE_N);
    auto N = std::vector<size_t>({100, 500, 1000, 2000, 3000, 5000});

    // Time for Jacobi method ~50s for N = 5000
    for (auto &size_n: N) {
        matrix = create_matrix(5 + E, -1, -1, size_n);
        vector_b = create_vector_b(size_n);
        time_method_measure(jacobi_method, matrix, vector_b);
    }
    std::cout << std::endl;

    // Time for Gauss-Seidel method ~30s for N = 5000
    for (auto &size_n: N) {
        matrix = create_matrix(5 + E, -1, -1, size_n);
        vector_b = create_vector_b(size_n);
        time_method_measure(gauss_seidel_method, matrix, vector_b);
    }
    std::cout << std::endl;

    // Time for LU factorization method ~10min for N = 5000
    for (auto &size_n: N) {
        matrix = create_matrix(5 + E, -1, -1, size_n);
        vector_b = create_vector_b(size_n);
        time_method_measure(lu_factory, matrix, vector_b);
    }

    return 0;
}
