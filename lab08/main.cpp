#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <utility>
#include <iomanip>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>

#define NODE_NR 10
#define N 1000

#define NODE_VECTOR std::array<double, NODE_NR>
#define VECTOR std::array<double, N>

double fun(double x) {
    return 1.0 / (1.0 + x * x);
    //return std::cos(2 * x);
}

double fun_second_der(double x) {
    double dx = 0.01;
    return (fun(x - dx) - 2 * fun(x) + fun(x + dx)) / (dx * dx);
}

gsl_matrix* construct_matrix(size_t size);

gsl_permutation* construct_permutation_vector(size_t size);

gsl_vector* construct_vector(size_t size);

void fill_matrix(gsl_matrix* mat, double lambda, double mi);

void print_matrix(const gsl_matrix* mat);

void print_vector(const gsl_vector* vec);

gsl_vector* wyznacz_M(const NODE_VECTOR& x_nodes, const NODE_VECTOR& y_nodes, double h, double alfa, double beta);

double wyznacz_Sx(const NODE_VECTOR& x_nodes, const NODE_VECTOR& y_nodes, gsl_vector* m, int n, double x, double h);

int main() {

    NODE_VECTOR nodes{0};
    NODE_VECTOR nodes_values{0};
    VECTOR net{0};
    VECTOR fun_values{0};
    VECTOR inter_fun_values{0};

    std::pair<double, double> range{-5, 5};
    double delta_x = (range.second - range.first) / static_cast<double>(N - 1);
    double h = (range.second - range.first) / static_cast<double>(NODE_NR - 1);

    for(auto i = 0; i < NODE_NR; ++i) {
        nodes[i] = range.first + i * h;
        nodes_values[i] = fun(nodes[i]);
    }

    for(auto i = 0; i < N; ++i) {
        net[i] = range.first + i * delta_x;
        fun_values[i] = fun(net[i]);
    }

    gsl_vector* m = wyznacz_M(nodes, nodes_values, h, fun_second_der(nodes[0]), fun_second_der(nodes[NODE_NR - 1]));

    std::ofstream wezly;
    wezly.open("wezly.txt");
    for(auto i = 0; i < NODE_NR; ++i) {
        wezly << nodes[i] << "\t" << fun_second_der(nodes[i]) << "\t" << gsl_vector_get(m, i) << std::endl;
        std::cout << fun_second_der(nodes[i]) << "\t" << gsl_vector_get(m, i) << "\t" << fabs(fun_second_der(nodes[i]) - gsl_vector_get(m, i)) << std::endl;
    }
    wezly.close();

    for(auto i = 0; i < N; ++i) {
        inter_fun_values[i] = wyznacz_Sx(nodes, nodes_values, m, NODE_NR, net[i], h);
    }

    std::ofstream plik;
    plik.open("wyniki.txt");
    for(int i = 0; i < N; ++i) {
        plik << net[i] << "\t" << inter_fun_values[i] << "\t" << fun_values[i] << std::endl;
    }
    plik.close();

    gsl_vector_free(m);

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

void lu_decomposition(gsl_matrix* mat, gsl_permutation* perm, int* signum) {
    gsl_linalg_LU_decomp(mat, perm, signum);
}

void fill_matrix(gsl_matrix* mat, double lambda, double mi) {
    for(size_t i = 0; i < mat->size1; ++i)
        for(size_t j = 0; j < mat->size2; ++j) {
            double val{0.0};
            if((i == 0 && j == 0) || ((i == mat->size1 - 1 && j == mat->size2 - 1))) {
                val = 1.0;
            } else if(i == j) {
                val = 2.0;
            } else if(i - 1 == j) {
                i == mat->size1 - 1 ? val = 0.0 : val = mi;
            } else if(i == j - 1) {
                i == 0 ? val = 0.0 : val = lambda;
            } else {
                val = 0.0;
            }

            gsl_matrix_set(mat, i, j, val);
        }
}

void print_matrix(const gsl_matrix* mat) {
    std::cout << "Macierz\n";
    for(size_t i = 0; i < mat->size1; ++i) {
        for(size_t j = 0; j < mat->size2; ++j) {
            double val = gsl_matrix_get(mat, i, j);
            std::cout << std::setw(20) << val;
        }
        std::cout << "\n";
    }
}

void print_vector(const gsl_vector* vec) {
    for(size_t i = 0; i < vec->size; ++i)
        std::cout << std::setw(20) << vec->data[i];
    std::cout << "\n";
}

gsl_vector* wyznacz_M(const NODE_VECTOR& x_nodes, const NODE_VECTOR& y_nodes, double h, double alfa, double beta) {
    double lambda = h / (h + h);
    double mi = 1.0 - lambda;

    gsl_matrix* mat = construct_matrix(NODE_NR);
    fill_matrix(mat, lambda, mi);

    gsl_vector* x = construct_vector(NODE_NR);
    gsl_vector* b = construct_vector(NODE_NR);
    for(auto i = 1; i < NODE_NR - 1; ++i) {
        double val = (6.0 / (h + h)) * ((y_nodes[i + 1] - y_nodes[i]) / (h) - (y_nodes[i] - y_nodes[i - 1]) / (h));
        gsl_vector_set(b, i, val);
    }
    gsl_vector_set(b, 0, alfa);
    gsl_vector_set(b, NODE_NR - 1, beta);

    gsl_permutation* p = construct_permutation_vector(NODE_NR);
    int s;
    gsl_linalg_LU_decomp(mat, p, &s);
    gsl_linalg_LU_solve(mat, p, b, x);

    gsl_matrix_free(mat);
    gsl_vector_free(b);
    gsl_permutation_free(p);

    return x;
}

double wyznacz_Sx(const NODE_VECTOR& x_nodes, const NODE_VECTOR& y_nodes, gsl_vector* m, int n, double x, double h) {
    auto i = 0;
    for(auto j = 1; j < n; ++j) {
        ++i;
        if(x >= x_nodes[j - 1] && x <= x_nodes[j]) {
            break;
        }
    }

    double A = (y_nodes[i] - y_nodes[i - 1]) / h - (gsl_vector_get(m, i) - gsl_vector_get(m, i - 1)) * h / 6.0;
    double B = y_nodes[i - 1] - gsl_vector_get(m, i - 1) * h * h / 6.0;

    return   gsl_vector_get(m, i - 1) * std::pow(x_nodes[i] - x, 3) / (6 * h)
           + gsl_vector_get(m, i) * std::pow(x - x_nodes[i - 1], 3) / (6 * h)
           + A * (x - x_nodes[i - 1])
           + B;
}