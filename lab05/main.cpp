#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include </usr/include/gsl/gsl_eigen.h>

#define n 7
#define K_val 7
#define IT_MAX 12

#define MATRIX std::array<std::array<double, n>,n>
#define VECTOR std::array<double, n>

void print_matrix(const MATRIX& mat);
VECTOR mul_mat_vec(const MATRIX& A, const VECTOR& x_vec);
VECTOR mul_vec_scalar(const VECTOR& vec, double scalar);
double dot_product(const VECTOR& vec1, const VECTOR& vec2);
MATRIX tensor_product(const VECTOR& vec1, const VECTOR& vec2, double lambda);
double euclides_norm(const VECTOR& vec);
MATRIX mul_mat_mat(const MATRIX& mat1, const MATRIX& mat2);
MATRIX transpose(const MATRIX& mat);

int main()
{

	MATRIX A;
	MATRIX W;
	MATRIX X;
	VECTOR x_vec;
	VECTOR x_vec_next;

	double lambda;

	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			A[i][j] = 1.0/(sqrt(2.0 + abs(i-j)));
			W[i][j] = A[i][j];
		}
	}
		
	print_matrix(A);
	std::cout << "\n\n";

	std::ofstream plik1;
	plik1.open("lambdas.txt");


	for(int k = 0; k < K_val; ++k) 
	{
		std::cout << "k = " << k << ":\t";
		plik1 << "k = " << k << ":\t";

		for(int ini = 0; ini < n; ++ini) x_vec[ini] = 1.0;
		for(int i = 0; i <= IT_MAX; ++i)
		{
			x_vec_next = mul_mat_vec(A, x_vec);
			double l = dot_product(x_vec_next, x_vec);
			double m = dot_product(x_vec, x_vec);
			lambda = l/m;

			std::cout << lambda << " ";
			plik1 << lambda << " ";

			x_vec = mul_vec_scalar(x_vec_next, 1.0/euclides_norm(x_vec_next));
		}

		for(int x = 0; x < n; ++x)
		{
			X[x][k] = x_vec[x];
		}


		std::cout << "\n";
		plik1 << "\n";
		
		for(int x = 0; x < n; ++x)
		{
			for(int y = 0; y < n; ++y)
			{
				A[x][y] = A[x][y] - tensor_product(x_vec, x_vec, lambda)[x][y]; 	
			}
		}
		
	}
	plik1.close();

	MATRIX D = mul_mat_mat(mul_mat_mat(transpose(X), W), X);

	std::cout << "\n";
	std::cout << "\n";
	print_matrix(D);


	std::ofstream plik2;
	plik2.open("d_matrix.txt");
	for(int i = 0; i < n; ++i)
		{
		for(int j = 0; j < n; ++j)
		{
			plik2 << D[i][j] << "\t";
		}
		plik2 << "\n";
	}
	plik2.close();

    return 0;
}

void print_matrix(const MATRIX& mat)
{
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			std::cout << mat[i][j] << "\t";
		}
		std::cout << "\n";
	}
}

VECTOR mul_mat_vec(const MATRIX& A, const VECTOR& vec)
{
	VECTOR result;
	for(int i = 0; i < n; ++i)
	{
		result[i] = 0.0;
		for(int j = 0; j < n; ++j)
		{
			result[i] += A[i][j] * vec[j];
		}
	}

	return result;
}

VECTOR mul_vec_scalar(const VECTOR& vec, double scalar)
{
	VECTOR result;
	for(int i = 0; i < n; ++i)
	{
		result[i] = vec[i] * scalar;
	}

	return result;
}

double dot_product(const VECTOR& vec1, const VECTOR& vec2)
{
	double sum = 0.0;
	for(int i = 0; i < n; ++i)
		sum += vec1[i] * vec2[i];

	return sum;
}

MATRIX tensor_product(const VECTOR& vec1, const VECTOR& vec2, double lambda)
{
	MATRIX result;
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			result[i][j] = lambda * vec1[i] * vec2[j];
		}
	}

	return result;
}

double euclides_norm(const VECTOR& vec)
{
	return sqrt(dot_product(vec, vec));
}

MATRIX mul_mat_mat(const MATRIX& mat1, const MATRIX& mat2)
{
	MATRIX result;
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			result[i][j] = 0.0;
			for(int k = 0; k < n; ++k) result[i][j] += mat1[i][k] * mat2[k][j];
		}
	}
	return result;
}

MATRIX transpose(const MATRIX& mat)
{
	MATRIX result;
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			result[i][j] = mat[j][i];
		}
	}
	return result;
}