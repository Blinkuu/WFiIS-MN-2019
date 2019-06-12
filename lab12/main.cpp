#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <utility>
#include <iomanip>
#include <functional>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

constexpr int n = 7;
constexpr double delta = 0.0000001;

double f1(double x) {
	return std::sin(x + delta)/(x + delta);
}

double f2(double x) {
	return (std::cos(x + delta) - std::exp(x + delta))/std::sin(x + delta);
}

double f3(double x) {
	return -1.0/((x + delta) * std::exp(1.0/(x + delta)));
}


double h(const std::pair<double, double>& r, int iter) {
	return (r.second - r.first)/std::pow(2, iter);
}

int main() 
{

	std::function<double(double)> fun = f3;
	std::pair<double, double> range(1.0, 0.0);


	std::array<std::array<double, n + 1>, n + 1> d_arr{0.0};

	double h_val = h(range, 0);

	d_arr[0][0] = h_val*(0.5*fun(range.first) + 0.5*fun(range.second));

	for(int i = 1; i <= n; ++i) {

		auto sum = [&]() {
			double s = 0.0;
			for(int j = 1; j <= std::pow(2,i - 1); ++j) {
				s += fun(range.first + (2*j - 1)*h(range, i));
			}
			return s;
		};

		d_arr[i][0] = 0.5*d_arr[i-1][0] + h(range, i)*sum();
	}

	for(int k = 1; k <= n; ++k) {
		for(int i = k; i <= n; ++i) {
			d_arr[i][k] = (std::pow(4, k)*d_arr[i][k-1] - d_arr[i-1][k-1])/(std::pow(4,k) - 1.0);
		}
	}

	
	for(int k = 0; k <= n; ++k) {
		for(int i = 0; i <= n; ++i) {
			std::cout << d_arr[k][i] << " ";
		}
		std::cout << std::endl;
	}

	std::ofstream file;
	file.open("wyniki3.txt");

	for(int k = 0; k <= n; ++k) {
		for(int i = 0; i <= n; ++i) {
			file << d_arr[k][i] << " ";
		}
		file << std::endl;
	}

	file.close();

    return 0;
}
