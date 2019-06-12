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

constexpr int N = 2;
constexpr double T = 1.0;
constexpr double t_max = 3.0 * T;
constexpr double s = T/20.0;
constexpr double w = 2.0*M_PI/T;

double f1_noise(double t) {
	double noise = std::rand()/(RAND_MAX + 1.0) - 0.5;
	return std::sin(w*t) + std::sin(2*w*t) + std::sin(3*w*t) + noise;
}

double f1(double t) {
	return std::sin(w*t) + std::sin(2*w*t) + std::sin(3*w*t);
}

double g1(double t) {
	return 1.0/(s*std::sqrt(2*M_PI)) * std::exp(-(t*t)/(2*s*s));
}

int main() 
{

	for(int k = 8; k <= 12; k += 2)
	{
		int arr_size = 2 * std::pow(N,k);
		double* f1_vals = new double[arr_size];
		double* f1_noise_vals = new double[arr_size];
		double* f1_noise_vals_file = new double[arr_size];
		double* g1_vals_1 = new double[arr_size];
		double* g1_vals_2= new double[arr_size];
		double dt = t_max/(arr_size/2.0);
		double t = 0.0;
		for(int i = 0; i < arr_size; ++i) {
			if(i % 2 == 0) {
				f1_vals[i] = f1(t);
				f1_noise_vals[i] = f1_noise(t);
				f1_noise_vals_file[i] = f1_noise(t);
				g1_vals_1[i] = g1(t);
				g1_vals_2[i] = g1(t);

				t += dt;
			} else {
				f1_vals[i] = 0.0;
				f1_noise_vals[i] = 0.0;
				g1_vals_1[i] = 0.0;
				g1_vals_2[i] = 0.0;
			}
		}

		int success_f = gsl_fft_complex_radix2_forward(f1_noise_vals, 1, arr_size/2);
		int success_g1 = gsl_fft_complex_radix2_forward(g1_vals_1, 1, arr_size/2);
		int success_g2 = gsl_fft_complex_radix2_backward(g1_vals_2, 1, arr_size/2);

		double* f_final = new double[arr_size];

		for(int i = 0; i < arr_size/2; ++i) {
			double a1 = f1_noise_vals[2*i];
			double b1 = f1_noise_vals[2*i + 1];

			double a2 = g1_vals_1[2*i] + g1_vals_2[2*i];
			double b2 = g1_vals_1[2*i + 1] + g1_vals_2[2*i + 1];

			f_final[2*i] = a1 * a1 - b1 * b2;
			f_final[2*i + 1] = a1 * b2 + a2 * b1;
		}

		int success_f_backward = gsl_fft_complex_radix2_backward(f_final, 1, arr_size/2);

		double f_max = std::abs(f_final[0]);

		for(int i = 2; i < arr_size; i+=2) {
			if( std::abs(f_final[i]) > f_max)
				f_max =  std::abs(f_final[i]);
		}

		std::ofstream file;
		file.open(std::string("wyniki_k_").append(std::to_string(k)).append(".txt"));
		
		for(int i = 0; i < arr_size; i += 2) {
			file << i << '\t' << f1_vals[i] << '\t' << f1_noise_vals_file[i] << '\t' << f_final[i]*2.5/f_max << std::endl;
		}
		
		file.close();

		delete[] f1_vals;
		delete[] f1_noise_vals;
		delete[] f1_noise_vals_file;
		delete[] g1_vals_1;
		delete[] g1_vals_2;
		delete[] f_final;
	}

    return 0;
}
