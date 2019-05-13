#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <utility>
#include <iomanip>
#include <functional>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>

#define M_S 30
#define M_C 30

#define N 100

double f1(double x)
{
	return 2*std::sin(x) + std::sin(2*x) + 2*std::sin(3*x);
}

double f1_noise(double x)
{
	double alfa = []() { return (rand()/(RAND_MAX + 1.0)) - 0.5; }();
	return 2*std::sin(x) + std::sin(2*x) + 2*std::sin(3*x) + alfa;
}

double f2(double x)
{
	return 2*std::sin(x) + std::sin(2*x) + 2*std::cos(x) + std::cos(2*x);
}

double f3(double x)
{
	return 2*std::sin(1.1*x) + std::sin(2.1*x) + 2*std::cos(3.1*x);
}

double F(double x, const std::array<double, M_S + 1>& a, const std::array<double, M_C>& b)
{
	auto sum_sin = [&](const std::array<double, M_S + 1>& arr, int m) {
		double s = 0.0;
		for(int k = 1; k <= m; ++k)
		{
			s += arr[k] * sin(k*x);
		}
		return s;
	};

	auto sum_cos = [&](const std::array<double, M_C>& arr, int m) {
		double s = 0.0;
		for(int k = 0; k < m; ++k)
		{
			s += arr[k] * cos(k*x);
		}
		return s;
	};

	double sum_a = sum_sin(a, M_S);
	double sum_b = sum_cos(b, M_C);

	return sum_a + sum_b;
}

int main() 
{
	std::function<double(double)> fun = f1_noise;

	std::array<double, N> x{0};
	std::array<double, N> y{0};

	std::array<double, M_S + 1> a{0};
	std::array<double, M_C> b;

	std::pair<double, double> range{0.0, 2.0*M_PI};
	double step = (range.second - range.first)/(N - 1);

	for(int i = 0; i < N; ++i)
	{
		x[i] = i * step;
		y[i] = fun(x[i]);
	}

	for(int k = 1; k <= M_S; ++k)
	{
		a[k] = 0.0;
		for(int i = 0; i < N; ++i)
		{
			a[k] +=  y[i] * sin(k * x[i])/(N / 2.0);
		}
	}

	if(M_C)
	{
		b[0] = (1.0/N) * [&]() {
		double s = 0.0;
		for(auto x : y)
		{
			s += x;
		}
		return s;
		}();
	}


	for(int k = 1; k < M_C; ++k)
	{
		b[k] = 0.0;
		for(int i = 0; i < N; ++i)
		{
			b[k] +=  y[i] * cos(k * x[i])/(N / 2.0);
		}
	}
  
	std::array<double, N> f_vals{0};

	for(int i = 0; i < N; ++i)
	{
		f_vals[i] = F(x[i], a, b);
	}

	std::ofstream file;
	file.open("wyniki5.txt");

	for(int i = 0; i < N; ++i)
	{
		file << x[i] << '\t' << y[i] << '\t' << f_vals[i] << std::endl;
	}

	file.close();

    return 0;
}
