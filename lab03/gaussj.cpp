//
// Created by blinku on 15.03.19.
//
#include <iostream>
#include <array>
#include <cmath>
#include <fstream>
#include <chrono>

#include "/opt/NR/numerical_recipes.c/nrutil.h"
#include "/opt/NR/numerical_recipes.c/nrutil.c"
#include "/opt/NR/numerical_recipes.c/gaussj.c"

#define N 10000
#define m 5

int main() {

    float** MatA = matrix(1, N, 1, N);
    float** VecB = matrix(1, N, 1, 1);

    for(int i = 1; i <= N; i++){
        VecB[i][1] = i + 1.0f;
        for(int j = 1; j <= N; j++){
            MatA[i][j] = 0.0;
            if( abs(i - j) <= m){
                MatA[i][j] = 1.0 / (1 + abs(i - j));
            }
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    gaussj(MatA, N, VecB, 1);
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout <<"It took me: " << std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count() << " s\n\n";

    free_matrix(MatA, 1, N, 1, N);
    free_matrix(VecB, 1, N, 1, 1);

    return 0;
}
