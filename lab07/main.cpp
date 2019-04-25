#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <utility>

#define NODE_NR 20
#define STEPS 1000

double f_fun(double x)
{
	return exp(-(x*x));
}

double lagrange_interpolation(double x, const std::array<double, NODE_NR + 1>& nodes, const std::array<double, NODE_NR + 1>& f_values)
{
	double w = 0.0;

	for(int j = 0; j <= NODE_NR; ++j)
	{
		double u = 1.0;
		for(int i = 0; i <= NODE_NR; ++i)
		{
			if(i!=j)
			{
				u *= (x - nodes[i])/(nodes[j] - nodes[i]);
			}
		}
		w += f_values[j] * u;
	}

	return w;
}

int main()
{

	const std::pair<double, double> range{-5.0, 5.0};
	std::array<double, STEPS> domain{0};
	std::array<double, STEPS> f_values{0};

	std::array<double, STEPS> lagrange_values{0};

	std::array<double, NODE_NR + 1> nodes{};
	std::array<double, NODE_NR + 1> f_values_in_nodes{0};

	double step = (range.second- range.first)/STEPS;
	double node_step = (range.second - range.first)/NODE_NR;
	
	for(int i = 0; i < STEPS; ++i)
	{
		domain[i] = range.first + i * step;
		f_values[i] = f_fun(domain[i]);
	}

	for(int i = 0; i <= NODE_NR; ++i)
	{
		// Czebyszew:
		//nodes[i] = 0.5*((range.second - range.first)*cosf(M_PI*(2*i + 1)/(2*NODE_NR + 2)) + range.first + range.second);
		// Ordinary:
		nodes[i] = range.first + i * node_step;
		f_values_in_nodes[i] = f_fun(nodes[i]);
	}

	for(int i = 0; i < STEPS; ++i)
	{
		lagrange_values[i] = lagrange_interpolation(domain[i], nodes, f_values_in_nodes);
		std::cout << lagrange_values[i] << std::endl;
	}

	std::ofstream file;
	file.open("data.txt");
	for(int i = 0; i < STEPS; ++i)
	{
		file << i << '\t' << f_values[i] << '\t' << lagrange_values[i] << std::endl;
	}
	file.close();

	return 0;
}