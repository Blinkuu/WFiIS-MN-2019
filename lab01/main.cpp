#include <iostream>
#include <iomanip>
#include <fstream>

#include "gsl/gsl_math.h"
#include "gsl/gsl_linalg.h"
//#include "usr/include/gsl/gsl_math.h"
//#include "usr/include/gsl/gsl_linalg.h"

gsl_matrix* construct_matrix(size_t size);
gsl_permutation* construct_permutation_vector(size_t size);
gsl_vector* construct_vector(size_t size);

void copy_matrix(gsl_matrix* dest, const gsl_matrix* src);
void fill_matrix_A(gsl_matrix* mat);

void lu_decomposition(gsl_matrix* mat, gsl_permutation* perm, int* signum);
double calculate_det_lu(gsl_matrix* lu, int signum);
void calculate_inverse_matrix(gsl_matrix* dest, const gsl_matrix* lu, gsl_permutation* p);
void calculate_matrix_product(gsl_matrix* dest, const gsl_matrix* A, const gsl_matrix* B);

void print_matrix(const gsl_matrix* mat);
void print_vector(const gsl_vector* vec);

double get_norm(const gsl_matrix* A);

void write_to_file(const gsl_matrix* lu, const double det_A, const gsl_matrix* matrix_multiplication_result);

int main() {

    // Definitions
    const int n = 4;
    int signum = 0;
    gsl_matrix* A = construct_matrix(n);
    gsl_matrix* A_copy = construct_matrix(n);
    gsl_matrix* B = construct_matrix(n);
    gsl_matrix* C = construct_matrix(n);
    gsl_permutation* p = construct_permutation_vector(n);

    // Matrix setup
    fill_matrix_A(A);
    copy_matrix(A_copy, A);

    // LU decomposition
    lu_decomposition(A, p, &signum);

    // LU determinant
    double det_A = calculate_det_lu(A, signum);
    std::cout << "Wyznacznik macierzy A: " << det_A << "\n\n";

    // Calculating inverse matrix from LU and storing it in B
    calculate_inverse_matrix(B, A, p);

    // Calculating matrix product and storing it in C
    calculate_matrix_product(C, A_copy, B);

    // Printing C matrix
    print_matrix(C);

    // Writing important stuff to files
    write_to_file(A, det_A, C);

    // Calculating condition number
    double condition_number = get_norm(A_copy) * get_norm(B);
    std::cout << "\nWskaznik uwarunkowania macierzy: " << condition_number << "\n";

    // Freeing memory
    gsl_matrix_free(A);
    gsl_matrix_free(A_copy);
    gsl_matrix_free(B);
    gsl_matrix_free(C);
    gsl_permutation_free(p);

    return 0;
}

gsl_matrix* construct_matrix(const size_t size) {
    return gsl_matrix_calloc(size, size);
}

gsl_permutation* construct_permutation_vector(const size_t size) {
    return gsl_permutation_alloc(size);
}

gsl_vector* construct_vector(const size_t size) {
    return gsl_vector_calloc(size);
}

void copy_matrix(gsl_matrix* dest, const gsl_matrix* src) {
    gsl_matrix_memcpy(dest, src);
}

void fill_matrix_A(gsl_matrix* mat) {
    for (size_t  i = 0; i < mat->size1; ++i)
        for (size_t  j = 0; j < mat->size2; ++j) {
            double val = 1.0 / ((double) i + (double) j + 2.0);
            gsl_matrix_set(mat, i, j, val);
        }
}

void lu_decomposition(gsl_matrix* mat, gsl_permutation* perm, int* signum){
    gsl_linalg_LU_decomp(mat, perm, signum);
}

double calculate_det_lu(gsl_matrix* lu, int signum) {
    return gsl_linalg_LU_det(lu, signum);
}

void calculate_inverse_matrix(gsl_matrix* dest, const gsl_matrix* lu, gsl_permutation* p) {
    size_t n = lu->size1;
    gsl_vector* b_vec = construct_vector(n);
    gsl_vector* x_vec = construct_vector(n);

    for (size_t k = 0; k < n; ++k) {
        for (size_t i = 0; i < n; ++i)
            gsl_vector_set(b_vec, i, 0.0);
        gsl_vector_set(b_vec, k, 1.0);

        gsl_linalg_LU_solve(lu, p, b_vec, x_vec);

        for (size_t i = 0; i < n; ++i)
            gsl_matrix_set(dest, i, k, gsl_vector_get(x_vec, i));
    }

    gsl_vector_free(b_vec);
    gsl_vector_free(x_vec);
}

void calculate_matrix_product(gsl_matrix* dest, const gsl_matrix* A, const gsl_matrix* B) {
    size_t n = A->size1;

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < n; ++k) {
                sum += (gsl_matrix_get(A, i, k) * gsl_matrix_get(B, k, j));
            }
            gsl_matrix_set(dest, i, j, sum);
        }
    }
}

void print_matrix(const gsl_matrix* mat) {
    std::cout << "Macierz\n";
    for (size_t i = 0; i < mat->size1; ++i) {
        for (size_t j = 0; j < mat->size2; ++j) {
            double val = gsl_matrix_get(mat, i, j);
            std::cout << std::setw(20) << val;
        }
        std::cout << "\n";
    }
}

void print_vector(const gsl_vector* vec) {
    for (size_t i = 0; i < vec->size; ++i)
        std::cout << std::setw(5) << vec->data[i];
    std::cout << "\n";
}

double get_norm(const gsl_matrix* A) {
    double max = gsl_matrix_get(A, 0, 0);
    for(size_t i = 0; i < A->size1; ++i)
        for(size_t j = 0; j < A->size2; ++j)
            if(gsl_matrix_get(A, i, j) > max)
                max = gsl_matrix_get(A, i, j);
    return max;
}

void write_to_file(const gsl_matrix* lu, const double det_A, const gsl_matrix* matrix_multiplication_result) {
    std::ofstream file;
    file.open("macierz_lu.txt");
    for(size_t i = 0; i < lu->size1; ++i) {
        for (size_t j = 0; j < lu->size2; ++j)
            if (j == i)
                file << std::setw(20)<<  gsl_matrix_get(lu, i, j);
    }
    file.close();

    file.open("macierz_A_wyznacznik.txt");
    file << det_A;
    file.close();

    file.open("macierz_A_invA.txt");
    for(size_t i = 0; i < matrix_multiplication_result->size1; ++i) {
        for (size_t j = 0; j < matrix_multiplication_result->size2; ++j)
            file << std::setw(20) <<  gsl_matrix_get(matrix_multiplication_result, i, j);
        file << "\n";
    }
    file.close();
}