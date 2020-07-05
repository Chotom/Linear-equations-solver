/* *********************************************************************************************************************
 * Author: Tomasz Czochański
 *
 * Methods for calculate system of linear equations and print them.
 **********************************************************************************************************************/

#include <vector>
#include <iostream>
#include <chrono>
#include "vector_functions.h"

// PRINTS --------------------------------------------------------------------------------------------------------------
template<class F, class Arg1, class Arg2>
void time_method_measure(F func, Arg1 arg1, Arg2 arg2) {

    // measure time
    auto start = std::chrono::steady_clock::now();
    func(arg1, arg2);
    auto end = std::chrono::steady_clock::now();
    auto delta_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    // print time
    std::cout << "Time measure in milliseconds: " << delta_time << " ms" << std::endl;
}

void print_methods_values(char *name, std::vector<double> &result, std::vector<double> &residuum, int iter,
                          bool print_vectors) {
    // print method
    std::cout << std::endl;
    std::cout << "method: " << name << std::endl;;

    if (print_vectors) {
        //print result
        std::cout << "result: ";
        print_vector(result);
        std::cout << std::endl;

        //print residuum
        std::cout << "residuum: ";
        print_vector(residuum);
        std::cout << std::endl;

        //print residuum norm
        std::cout << "residuum norm: ";
        std::cout << calc_norm_residuum_vector(residuum);
        std::cout << std::endl;
    }
    //print iterations
    if (iter == 0)
        std::cout << "iterations: " << "this method has no iterations" << std::endl;
    else
        std::cout << "iterations: " << iter << std::endl;
}


// CALC METHODS --------------------------------------------------------------------------------------------------------
auto jacobi_method(std::vector<std::vector<double>> &matrix, std::vector<double> &vector_b) {

    //Init
    const double NORM = 10e-10;
    const size_t N = vector_b.size();
    auto vector_x = create_vector_with_value(1, N);
    auto vector_next_x = create_vector_with_value(1, N);
    auto vector_residuum = create_vector_with_value(1, N);
    int iter = 0;

    // Calc jacobi methods
    while (calc_norm_residuum_vector(vector_residuum) > NORM) {
        for (int i = 0; i < N; ++i) {
            vector_next_x[i] = (vector_b[i] - sum_multiply_vector(0, i, matrix[i], vector_x)
                                - sum_multiply_vector(i + 1, N, matrix[i], vector_x)) / matrix[i][i];
        }

        vector_x = vector_next_x;
        vector_residuum = matrix_multiply_vector(matrix, vector_x);
        vector_residuum = vector_subtract_vector(vector_residuum, vector_b);
        iter++;
    }

    //print
    print_methods_values("Jacobi", vector_x, vector_residuum, iter, PRINT_VECTORS);

    return vector_x;
}

auto gauss_seidel_method(std::vector<std::vector<double>> &matrix, std::vector<double> &vector_b) {

    // Init
    const double NORM = 1e-9;
    const size_t N = vector_b.size();
    auto vector_x = create_vector_with_value(1, N);
    auto vector_next_x = create_vector_with_value(1, N);
    auto vector_residuum = create_vector_with_value(1, N);
    int iter = 0;

    // Calc Gauss-Seidel method
    while (calc_norm_residuum_vector(vector_residuum) > NORM) {
        for (int i = 0; i < N; ++i) {
            vector_next_x[i] = (vector_b[i] - sum_multiply_vector(0, i, matrix[i], vector_next_x)
                                - sum_multiply_vector(i + 1, N, matrix[i], vector_x)) / matrix[i][i];
        }

        vector_x = vector_next_x;
        vector_residuum = matrix_multiply_vector(matrix, vector_x);
        vector_residuum = vector_subtract_vector(vector_residuum, vector_b);
        iter++;
    }

    // print values
    print_methods_values("Gauss-Seidel", vector_x, vector_residuum, iter, PRINT_VECTORS);

    return vector_x;
}


auto lu_factory(std::vector<std::vector<double>> &matrix, std::vector<double> &vector_b) {

    // Init
    const size_t N = vector_b.size();
    auto matrix_U = matrix;
    auto matrix_L = create_matrix(1, 0, 0, N);
    auto vector_y = create_vector_with_value(1, N);
    auto vector_x = create_vector_with_value(1, N);
    auto vector_residuum = create_vector_with_value(1, N);
    int iter = 0;

    // Calc L and U
    for (int i = 0; i < N - 1; ++i) {
        for (int j = i + 1; j < N; ++j) {
            matrix_L[j][i] = matrix_U[j][i] / matrix_U[i][i];
            // U[j][i:N] = U[j][i:N] - L[j][i]*U[i][i:N]
            for (int k = i; k < N; k++)
                matrix_U[j][k] = matrix_U[j][k] - matrix_L[j][i] * matrix_U[i][k];
        }
    }

    // resolve Ly = b
    vector_y[0] = vector_b[0] / matrix_L[0][0];
    for (int i = 1; i < N; ++i)
        vector_y[i] = (1 / matrix_L[i][i]) * (vector_b[i] - sum_multiply_vector(0, i, matrix_L[i], vector_y));

    // resolve Ux = y
    vector_x[N - 1] = vector_y[N - 1] / matrix_U[N - 1][N - 1];
    for (int i = N - 1; i >= 0; --i)
        vector_x[i] = (1 / matrix_U[i][i]) * (vector_y[i] - sum_multiply_vector(i + 1, N, matrix_U[i], vector_x));

    // print values
    vector_residuum = matrix_multiply_vector(matrix, vector_x);
    vector_residuum = vector_subtract_vector(vector_residuum, vector_b);
    print_methods_values("LU factorization", vector_x, vector_residuum, iter, PRINT_VECTORS);

    return vector_x;
}
