#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <utility>
#include <iomanip>
#include <functional>

double fun1(double x)
{
	return std::log(std::pow(x, 5) + 3*std::pow(x, 2) + x + 9);
}

double fun2(double x)
{
	return std::pow(x, 6);
}

int main() 
{
	std::function<double(double)> fun = fun2;
	std::pair<double, double> range{-4.0, 1.0};

	double x_a = range.first;
	double x_b = range.second;
	double r = (std::sqrt(5) - 1.0)/2.0;
	double l_1 = 1.0/3.0;
	double l_2 = 2.0/3.0;
	unsigned counter = 0;
	while(std::abs(x_a - x_b) > std::pow(10,-6))
	{
		double x1 = x_a + l_1*(x_b - x_a);
		double x2 = x_a + l_2*(x_b - x_a);

		fun(x1) < fun(x2) ? x_b = x2 : x_a = x1;
	
		++counter;
	}

	double xmin = (x_a + x_b)/2.0;

	std::ofstream file;
	file.open("zaleznos2.txt");

	file << counter << '\t' << xmin << '\t' << std::abs(0.0 - xmin) << std::endl;

	file.close();

    return 0;
}
