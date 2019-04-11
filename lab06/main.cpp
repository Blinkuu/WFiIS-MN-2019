#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <iomanip>
#include <functional>

#define MAX_ITER 20

void newton(std::function<double(double)> func, std::function<double(double)> func_der, double real_zero, double start);
void metoda_siecznych(std::function<double(double)> func, std::function<double(double)> func_der, double real_zero, double start1, double start2);

double f_function(double x);
double f_derivative_function(double x);

double g_function(double x);
double g_derivative_function(double x);

double epsilon(double x1, double x2);
double p(double e0, double e1, double e2);

int main()
{

	newton(f_function, f_derivative_function, 1.0, 3.0);
	std::cout << "\n\n";
	metoda_siecznych(f_function, f_derivative_function, 1.0, 3.0, 3.01);
	std::cout << "\n\n";

	//newton(g_function, g_derivative_function, -3.28427753730695, -20.0);
	//std::cout << "\n\n";
	//metoda_siecznych(g_function, g_derivative_function, -3.28427753730695, -20.0, -20.1);
	//std::cout << "\n\n";

	return 0;
}

void newton(std::function<double(double)> func, std::function<double(double)> func_der, double real_zero, double start)
{
	std::ofstream newton;
	newton.open("newton.txt");

	int i = 0;
	double x_zero = start;
	while(i++ < MAX_ITER)
	{
		double e_before = epsilon(x_zero, real_zero);
		x_zero = x_zero - func(x_zero)/func_der(x_zero);
		double e_now = epsilon(x_zero, real_zero);

		double for_next = x_zero - func(x_zero)/func_der(x_zero);

		double e_next = epsilon(for_next, real_zero);

		newton << i << "\t" << std::setprecision(15) << func(x_zero) << "\t" << x_zero << "\t" << e_now << "\t" << p(e_before, e_now, e_next) << std::endl;
		std::cout << x_zero << '\n'; 
	}

	newton.close();
}

void metoda_siecznych(std::function<double(double)> func, std::function<double(double)> func_der, double real_zero, double start1, double start2)
{
	std::ofstream sieczne;
	sieczne.open("sieczne.txt");

	int i = 0;
	double x_after = start1;
	double x_before = start2;
	while(i++ < MAX_ITER)
	{
		double e_before = epsilon(x_after, real_zero);
		double buf = x_after;
		x_after = x_after - (func(x_after)*(x_after - x_before))/(func(x_after) - func(x_before));
		x_before = buf;
		double e_now = epsilon(x_after, real_zero);

		double for_next = x_after - (func(x_after)*(x_after - x_before))/(func(x_after) - func(x_before));

		double e_next = epsilon(for_next, real_zero);

		sieczne << i << "\t" << std::setprecision(15) << func(x_after) << "\t" << x_after << "\t" << e_now << "\t" << p(e_before, e_now, e_next) << std::endl;
		std::cout << x_after << '\n'; 
	}

	sieczne.close();
}

double f_function(double x)
{
	return pow(log(x) - x, 6) - 1.0;
}

double f_derivative_function(double x)
{
	return 6.0*(1.0/x - 1.0)*pow(log(x) - x, 5);
}

double g_function(double x)
{
	return x*x*x + 2*x*x - 3*x + 4;
}

double g_derivative_function(double x)
{
	return 3*x*x + 4*x - 3;
}

double epsilon(double x1, double x2)
{
	return fabs(x1 - x2);
}

double p(double e0, double e1, double e2) 
{
	return log(e1/e2)/log(e0/e1);
}