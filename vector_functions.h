/* *********************************************************************************************************************
 * author: Tomasz Czochański
 *
 * Functions for vectors and matrix for general uses:
 *  create_vector_b(size_t)
 *  create_vector_with_value(double, size_t)
 *  matrix_multiply_vector(vector<vector<double>>, vector<double>)
 *  vector_subtract_vector(vector<double>, vector<double>)
 *
 *  calc_norm_residuum_vector(vector<double>)
 *  sum_multiply_vector(int, int, vector<double>, vector<double>)
 *
 *  create_matrix(double, double, double, size_t)
 *
 *  print_vector(vector<double>)
 *  print_matrix(vector<vector<double>>)
 **********************************************************************************************************************/

#include <vector>
#include <cmath>
#include <iostream>
#include "project_consts.h"

// --- RETURNS VECTOR --------------------------------------------------------------------------------------------------

/**
 * Create vector b
 *
 * @param N size of vector
 * @return vector filled with values n*(F + 1) for n-element
 */
auto create_vector_b(size_t N) {

    std::vector<double> vector_b;
    vector_b.reserve(N);

    for (int i = 0; i < N; ++i)
        vector_b.push_back(sin(i * (F + 1)));

    return vector_b;
}

/**
 * Create vector filled with given value
 *
 * @param x value for vector
 * @param N size of vector
 * @return vector filled with x for n-element
 */
auto create_vector_with_value(double x, size_t N) {

    std::vector<double> vector_x;
    vector_x.reserve(N);

    // Fill vector with same given values
    for (int i = 0; i < N; ++i)
        vector_x.push_back(x);

    return vector_x;
}

/**
 *  Multiply Matrix and vector
 *
 * @param matrix
 * @param vec
 * @return vector
 */
auto matrix_multiply_vector(std::vector<std::vector<double>> &matrix, std::vector<double> &vec) {

    const size_t N = vec.size();
    std::vector<double> vector_ret = create_vector_with_value(0, N);

    // Calc
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            vector_ret[i] += (matrix[i][j] * vec[j]);

    return vector_ret;
}

/**
 *  Subtract two vectors
 *
 * @param vec1
 * @param vec2
 * @return vector
 */
auto vector_subtract_vector(std::vector<double> &vec1, std::vector<double> &vec2) {

    const size_t N = vec1.size();
    std::vector<double> vector_ret = create_vector_with_value(0, N);

    // Calc
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            vector_ret[i] = vec1[i] - vec2[i];

    return vector_ret;
}

// --- RETURNS DOUBLE --------------------------------------------------------------------------------------------------

/**
 * Return value of the norm for residuum vector for iterative methods
 *
 * @param vec residuum vector
 * @return double value of sqrt of sum of the squares
 */
auto calc_norm_residuum_vector(std::vector<double> &vec) {

    double norm = 0;

    // Calc sum of the squares
    for (auto &val: vec)
        norm += val * val;

    return sqrt(norm);
}

/**
 * Return sum value of multiply vectors indexes from given start to end of interval
 *
 * @param start index of vector to begin interval
 * @param end   index of vector to end interval
 * @param first_vector  first vector
 * @param second_vector  second vector
 * @return double value of sum of first_vector[i]*second_vector[i]
 */
auto sum_multiply_vector(int start, int end, std::vector<double> &first_vector, std::vector<double> &second_vector) {

    double sum = 0;

    // Calc sum of multiplying
    for (int j = start; j < end; ++j)
        sum += first_vector[j] * second_vector[j];

    return sum;
}

// --- RETURNS MATRIX (VECTOR OF VECTORS) ------------------------------------------------------------------------------

/**
 * Create vector of vectors and fill with given values
 *
 * @param a1 value for main diagonal in matrix
 * @param a2 value for two neighbour diagonals
 * @param a3 value for two next neighbour diagonals
 * @param N size of matrix
 * @return vector<vector<double>> with given values and NxN size
 */
auto create_matrix(double a1, double a2, double a3, size_t N) {

    std::vector<std::vector<double>> matrix;

    // Init matrix
    for (int i = 0; i < N; ++i) {
        std::vector<double> row;
        row.reserve(N);

        for (int j = 0; j < N; ++j)
            row.push_back(0);
        matrix.push_back(row);
    }

    // Set values for matrix
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j)
                // Set diagonal values in matrix
                matrix[i][j] = a1;
            else if ((i + 1 < N || j + 1 < N) && (i + 1 == j || j + 1 == i))
                // Set values for two neighbour diagonals
                matrix[i][j] = a2;
            else if ((i + 2 < N || j + 2 < N) && (i + 2 == j || j + 2 == i))
                // Set values for next two neighbour diagonals
                matrix[i][j] = a3;
        }
    }

    return matrix;
}

// --- NO RETURNS ------------------------------------------------------------------------------------------------------

/**
 * Print values of vector
 *
 * @param vector vector<double>to print
 */
void print_vector(std::vector<double> &v) {
    std::cout << "Vector: ";
    for (auto &val: v)
        std::cout << "\t" << val;
}

/**
 * Print values of matrix
 *
 * @param matrix vector<vector<double>> to print
 */
void print_matrix(std::vector<std::vector<double>> &matrix) {
    std::cout << "Matrix: " << std::endl;
    for (auto &vec: matrix) {
        for (auto &val: vec)
            std::cout << val << "\t";
        std::cout << std::endl;
    }
}
