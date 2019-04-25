#include <iostream>
#include <array>
#include <cmath>
#include <fstream>
#include <chrono>

#define N 1000
#define m 5

#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define ABS(X)   ((X) > 0   ? (X) :-(X))

#define ARRAY std::array<double, N>
#define MATRIX std::array<std::array<double, N>, N>

void print_matrix(const MATRIX& A);
void print_vec(const ARRAY& vec);
double dot_product(const ARRAY& vec1, const ARRAY& vec2);
double euclides_norm(const ARRAY& vec1);
ARRAY scalar_vec_mul(double scalar, ARRAY& vec);
ARRAY vec_vec_sum(const ARRAY& vec1, const ARRAY& vec2);
ARRAY vec_vec_diff(const ARRAY& vec1, const ARRAY& vec2);
ARRAY mat_vec_mul(const MATRIX& mat, const ARRAY& vec);

int main() {

    MATRIX A;
    ARRAY b_vec;
    ARRAY x_vec;
    ARRAY r_vec;

    for(int i = 0; i < N; ++i) {
        b_vec[i] = i;
        x_vec[i] = 0.0;
        
        for(int j = 0; j < N; ++j)
            ABS(i - j) <= m ? A[i][j] = 1.0 / (1.0 + ABS(i - j)) : A[i][j] = 0.0;
    }

    // przyklad b
    //x_vec[0] = 1.0;

    std::ofstream file;
    file.open("data.txt");
    std::ofstream file2;
    file2.open("data2.txt");

    auto t1 = std::chrono::high_resolution_clock::now();

    unsigned k = 0;
    double e_norm;
    double x_norm;
    do {
        r_vec = vec_vec_diff(b_vec, mat_vec_mul(A, x_vec));
        e_norm = sqrt(dot_product(r_vec, r_vec));

        double alfa = dot_product(r_vec, r_vec)/dot_product(r_vec, mat_vec_mul(A, r_vec));

        x_vec = vec_vec_sum(x_vec, scalar_vec_mul(alfa, r_vec));
        x_norm = sqrt(dot_product(x_vec, x_vec));

        file << k++ << " " << e_norm << "\n";
        file2 << k++ << " " << x_norm << "\n";
    } while (e_norm > 1.E-6);

    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout <<"It took me: " << std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count() << " s\n\n";

    std::cout << k;

    file.close();
    file2.close();

    return 0;
}

void print_matrix(const MATRIX& A) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            std::cout << A[i][j];
        }
        std::cout << "\n";
    }
}

void print_vec(const ARRAY& vec) {
    for(int i = 0; i < N; ++i)
        std::cout << vec[i] << " ";
    std::cout << "\n";
}

double dot_product(const ARRAY& vec1, const ARRAY& vec2) {
    double dot = 0.0;
    for(int i = 0; i < N; ++i) {
        dot += vec1[i] * vec2[i];
    }
    return dot;
}

ARRAY scalar_vec_mul(double scalar, ARRAY& vec) {
    for (int i = 0; i < N; ++i)
        vec[i] *= scalar;
    return vec;
}

ARRAY vec_vec_sum(const ARRAY& vec1, const ARRAY& vec2) {
    ARRAY y;

    for(int i = 0; i < N; ++i)
        y[i] = vec1[i] + vec2[i];

    return y;
}

ARRAY vec_vec_diff(const ARRAY& vec1, const ARRAY& vec2) {
    ARRAY y;

    for(int i = 0; i < N; ++i)
        y[i] = vec1[i] - vec2[i];

    return y;
}

ARRAY mat_vec_mul(const MATRIX& mat, const ARRAY& vec) {
    ARRAY y;
    for(int i = 0; i < N; ++i) {
        int jmin = MAX(0, i - m);
        int jmax = MIN(i+m, N-1);
        y[i] = 0.0;
        for(int j = jmin; j <= jmax; ++j)
            y[i] += mat[i][j] * vec[j];
    }
    return y;
}