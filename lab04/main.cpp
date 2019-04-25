#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_eigen.h>

#define L 10
#define n 200
#define N 1

double del_kro(int i, int j) {
	double tmp = 0.0;
	if(i==j)
		++tmp;
	return tmp;
}

int main() {
	double del_x = static_cast<double>(L)/static_cast<double>(n+1);
	gsl_vector* x_vec = gsl_vector_alloc(n);
	for(int i = 0; i < n; ++i) {
		gsl_vector_set(x_vec, i, - L / 2.0 + del_x * (i + 1));
	}

 	gsl_matrix* A =  gsl_matrix_alloc(n, n);
 	gsl_matrix* B =  gsl_matrix_alloc(n, n);
 	gsl_vector* ro = gsl_vector_alloc(n);
 	gsl_vector* eval = gsl_vector_alloc(n);
 	gsl_matrix* evec = gsl_matrix_alloc(n,n);
 	gsl_eigen_gensymmv_workspace* w = gsl_eigen_gensymmv_alloc(n);

    std::ofstream file1;
    file1.open("data1.txt", std::ios_base::app);

	for(int alfa = 0; alfa <= 100; alfa+= 2) {
		for(int i = 0; i < n; ++i) {
			gsl_vector_set(ro, i, 1 + 4 * alfa * gsl_vector_get(x_vec, i) * gsl_vector_get(x_vec, i));
			for(int j = 0; j < n; ++j) {
				gsl_matrix_set(A, i, j, (-del_kro(i, j+1) + 2*del_kro(i, j) - del_kro(i, j-1))/del_x);
				gsl_matrix_set(B, i, j, (gsl_vector_get(ro, i)/static_cast<double>(N) * del_kro(i, j)));
			}
		}
		gsl_eigen_gensymmv(A, B, eval, evec, w);
		gsl_eigen_gensymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

		file1 << alfa;
		for(int i = 0; i < 6; ++i) {
		   file1 << " " << sqrt(gsl_vector_get(eval, i)) << " ";	
		}
		file1 << std::endl;

		if(alfa == 0) {
			std::ofstream file2;
			file2.open("data2.txt");
			for(int j = 0; j < n; ++j) {
				file2 << j << " ";
				for(int i = 0; i < 6; ++i) {
					file2 <<  gsl_matrix_get(evec, j, i) << " ";	
				}
				file2 << std::endl;
			}
			file2.close();
		} else if (alfa == 100) {
			std::ofstream file3;
			file3.open("data3.txt");
			for(int j = 0; j < n; ++j) {
				file3 << j << " ";
				for(int i = 0; i < 6; ++i) {
					file3 <<  gsl_matrix_get(evec, j, i) << " ";	
				}
				file3 << std::endl;
			}
			file3.close();
		}
	}

    file1.close();

    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_matrix_free(evec);

    gsl_vector_free(ro);
    gsl_vector_free(eval);

    return 0;
}
