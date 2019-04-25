#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <array>

#define N 500

int main() {

    double b = 2.0;
    double a = 1.0/2.0;
    double h = 2.0 * b / (N - 1);

    std::array<double, N+1> x_i;

    std::array<double, N+1> a_i;
    std::array<double, N+1> d_i;
    std::array<double, N+1> c_i;

    std::array<double, N+1> ro_i;

    std::array<double, N+1> analitic_i;

    for (size_t i = 1; i <= N; ++i) {
        x_i[i] = -b + h * (i-1);

        if ( x_i[i] >= - b && x_i[i] < -a) {
            ro_i[i] = 0.0;
            analitic_i[i] = x_i[i] / 16.0 + 1.0/8.0;
        } else if ( x_i[i] >= -a && x_i[i] < 0.0) {
            ro_i[i] = 1.0;
            analitic_i[i] = -pow(x_i[i],2) / 2.0 - (7.0 * x_i[i])/16.0;
        } else if ( x_i[i] == 0.0) {
            ro_i[i] = 0.0;
        } else if ( x_i[i] > 0.0 && x_i[i] <= a ) {
            ro_i[i] = -1.0;
            analitic_i[i] = pow(x_i[i],2) / 2.0 - (7.0 * x_i[i])/16.0;
        } else if ( x_i[i] > a && x_i[i] <= b) {
            ro_i[i] = 0.0;
            analitic_i[i] = x_i[i] / 16.0 - 1.0/8.0;
        }

        d_i[i] = -2.0 / pow(h, 2);
        a_i[i] = 1.0 / pow(h, 2);
        c_i[i] = a_i[i];
    }

    // warunki brzegowe
    d_i[1] = 1.0;
    c_i[1] = 0.0;
    ro_i[1] = 0.0;


    d_i[N] = 1.0;
    c_i[N] = 0.0;
    ro_i[N] = 0.0;

    std::array<double, N+1> u_i;
    std::array<double, N+1> l_i;

    u_i[1] = d_i[1];

    for (size_t i = 2; i <= N; ++i) {
        l_i[i] = a_i[i] / u_i[i-1];
        u_i[i] = d_i[i] - l_i[i] * c_i[i-1];
    }

    std::array<double, N+1> y_i;
    y_i[1] = -ro_i[1];

    for(size_t i = 2; i <= N; ++i) {
        y_i[i] = -ro_i[i] - l_i[i] * y_i[i - 1];
    }

    std::array<double, N+1> v_i;

    v_i[N] = y_i[N] / u_i[N];

    for(size_t i = N - 1; i >= 1; --i) {
        v_i[i] = (y_i[i] - c_i[i] * v_i[i + 1])/u_i[i];
    }

    std::ofstream file;
    file.open("output.txt");
    for(size_t i = 1; i <= N; ++i) {
        file << x_i[i] << " " << v_i[i] << std::endl;
    }
    file.close();

    file.open("analitic.txt");
    for(size_t i = 1; i <= N; ++i) {
        file << x_i[i] << " " << analitic_i[i] << std::endl;
    }
    file.close();

    return 0;
}
