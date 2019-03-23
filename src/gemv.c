#include "complexe.h"
#include "mnblas.h"
#include <stdio.h>
#include <stdlib.h>

float* copy_float (int size, const float *a) {
	float* copy = (float*) malloc(size * sizeof(float)) ;

	for (int i = 0; i < size; i++) {
		copy[i] = a[i] ;
	}

	return copy ;
}

double* copy_double (int size, const double *a) {
	double* copy = (double*) malloc(size * sizeof(double)) ;

	for (int i = 0; i < size; i++) {
		copy[i] = a[i] ;
	}

	return copy ;
}

complexe_float_t* copy_float_complexe (int size, const complexe_float_t *a) {
	complexe_float_t* copy = (complexe_float_t*) malloc(size * sizeof(complexe_float_t)) ;

	for (int i = 0; i < size; i++) {
		copy[i] = a[i] ;
	}

	return copy ;
}

complexe_double_t* copy_double_complexe (int size, const complexe_double_t *a) {
	complexe_double_t* copy = (complexe_double_t*) malloc(size * sizeof(complexe_double_t)) ;

	for (int i = 0; i < size; i++) {
		copy[i] = a[i] ;
	}

	return copy ;
}

float *transposee_s(const MNCBLAS_LAYOUT layout, const int m, const int n,
		    const float *a)
{
	float *abis = (float *)malloc(n * m * sizeof(float));

	if (layout == MNCblasRowMajor) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				abis[m * i + j] = a[i + n * j];
			}
		}
	} else {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				abis[n * i + j] = a[i + m * j];
			}
		}
	}

	return abis;
}

double *transposee_d(const MNCBLAS_LAYOUT layout, const int m, const int n,
		     const double *a)
{
	double *abis = (double *)malloc(n * m * sizeof(double));

	if (layout == MNCblasRowMajor) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				abis[m * i + j] = a[i + n * j];
			}
		}
	} else {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				abis[n * i + j] = a[i + m * j];
			}
		}
	}

	return abis;
}

void *transposee_c(const MNCBLAS_LAYOUT layout, const int m, const int n,
		   const void *a)
{
	complexe_float_t *abis =
	    (complexe_float_t *) malloc(n * m * sizeof(complexe_float_t *));

	if (layout == MNCblasRowMajor) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				abis[m * i + j] =
				    ((complexe_float_t *) a)[i + n * j];
			}
		}
	} else {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				abis[n * i + j] =
				    ((complexe_float_t *) a)[i + m * j];
			}
		}
	}

	return abis;
}

void *transposee_z(const MNCBLAS_LAYOUT layout, const int m, const int n,
		   const double *a)
{
	complexe_double_t *abis =
	    (complexe_double_t *) malloc(n * m * sizeof(complexe_double_t *));

	if (layout == MNCblasRowMajor) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				abis[m * i + j] =
				    ((complexe_double_t *) a)[i + n * j];
			}
		}
	} else {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				abis[n * i + j] =
				    ((complexe_double_t *) a)[i + m * j];
			}
		}
	}

	return abis;
}

void
mncblas_sgemv(const MNCBLAS_LAYOUT layout,
	      const MNCBLAS_TRANSPOSE TransA, const int M, const int N,
	      const float alpha, const float *A, const int lda,
	      const float *X, const int incX, const float beta,
	      float *Y, const int incY)
{
	int Nbis = N;
	int Mbis = M;
	float *Abis = (float *)malloc(N * M * sizeof(float));
	if (TransA != MNCblasNoTrans) {
		Nbis = M;
		Mbis = N;
	}
	int lenY = (1 + (Mbis - 1) * incY);
	if ((M == 0) || (N == 0)) {	//si la mtrice est nulle on a juste y=beta*y
		for (int i = 0; i < lenY; i += incY) {
			Y[i] = beta * Y[i];
		}
	} else {
		float *copyY = copy_float (Nbis, Y) ;
		if (TransA != MNCblasNoTrans) {	//On fait la transposée si demandée
			Abis = transposee_s(layout, M, N, A);
		} else {
			for (int k = 0; k < N * M; k++) {	//sinon abis est une copie de A
				Abis[k] = A[k];
			}
		}
		for (int j = 0; j < M * N; j++) {	//calcul de alpha*matrice
			Abis[j] *= alpha;
		}
		if (layout == MNCblasRowMajor) {
			for (int p = 0; p < Mbis; p++) {
				for (int q = 0; q < Nbis; q++) {
					Y[p * incY] +=
					    X[q * incX] * Abis[N * p + q];
				}
			}
		} else {
			for (int p = 0; p < Mbis; p++) {
				for (int q = 0; q < Nbis; q++) {
					Y[p * incY] +=
					    X[q * incX] * Abis[p + M * q];
				}
			}
		}
		mncblas_saxpy(Mbis, beta, copyY, incY, Y, incY);
	}
	return;
}

void mncblas_dgemv(MNCBLAS_LAYOUT layout,
		   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
		   const double alpha, const double *A, const int lda,
		   const double *X, const int incX, const double beta,
		   double *Y, const int incY)
{
	int Nbis = N;
	int Mbis = M;
	double *Abis = (double *)malloc(N * M * sizeof(double));
	if (TransA != MNCblasNoTrans) {
		Nbis = M;
		Mbis = N;
	}
	int lenY = (1 + (Mbis - 1) * incY);
	if ((M == 0) || (N == 0)) {	//si la mtrice est nulle on a juste y=beta*y
		for (int i = 0; i < lenY; i += incY) {
			Y[i] = beta * Y[i];
		}
	} else {
		double *copyY = copy_double (Nbis, Y) ;
		if (TransA != MNCblasNoTrans) {	//On fait la transposée si demandée
			Abis = transposee_d(layout, M, N, A);
		} else {
			for (int k = 0; k < N * M; k++) {	//sinon abis est une copie de A
				Abis[k] = A[k];
			}
		}
		for (int j = 0; j < M * N; j++) {	//calcul de alpha*matrice
			Abis[j] *= alpha;
		}
		if (layout == MNCblasRowMajor) {
			for (int p = 0; p < Mbis; p++) {
				for (int q = 0; q < Nbis; q++) {
					Y[p * incY] +=
					    X[q * incX] * Abis[N * p + q];
				}
			}
		} else {
			for (int p = 0; p < Mbis; p++) {
				for (int q = 0; q < Nbis; q++) {
					Y[p * incY] +=
					    X[q * incX] * Abis[p + M * q];
				}
			}
		}
		mncblas_daxpy(Mbis, beta, copyY, incY, Y, incY);
	}
	return;
}

void mncblas_cgemv(MNCBLAS_LAYOUT layout,
		   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
		   const void *alpha, const void *A, const int lda,
		   const void *X, const int incX, const void *beta,
		   void *Y, const int incY)
{
	int Nbis = N;
	int Mbis = M;
	complexe_float_t *Abis = (complexe_float_t *) malloc(N * M *
						sizeof(complexe_float_t));
	if (TransA != MNCblasNoTrans) {
		Nbis = M;
		Mbis = N;
	}
	int lenY = (1 + (Mbis - 1) * incY);
	if ((M == 0) || (N == 0)) {	//si la mtrice est nulle on a juste y=beta*y
		for (int i = 0; i < lenY; i += incY) {
			((complexe_float_t *) Y)[i] =
			    mult_complexe_float(*((complexe_float_t *) beta),
						((complexe_float_t *) Y)[i]);
		}
	} else {
		complexe_float_t *copyY = copy_float_complexe (Nbis, Y) ;
		if (TransA != MNCblasNoTrans) {	//On fait la transposée si demandée
			Abis = transposee_c(layout, M, N, A);
		} else {
			for (int k = 0; k < N * M; k++) {	//sinon abis est une copie de A
				Abis[k] = ((complexe_float_t *) A)[k];
			}
		}
		if (TransA == MNCblasConjTrans) {
			for (int incr = 0; incr < N * M; incr++) {
				Abis[incr] = conjg_float(Abis[incr]);
			}
		}
		for (int j = 0; j < M * N; j++) {	//calcul de alpha*matrice
			Abis[j] =
			    mult_complexe_float(*((complexe_float_t *) alpha),
						((complexe_float_t *) Abis)[j]);
		}
		if (layout == MNCblasRowMajor) {
			for (int p = 0; p < Mbis; p++) {
				for (int q = 0; q < Nbis; q++) {
					((complexe_float_t *)Y)[p * incY] =
					    add_complexe_float
					    (mult_complexe_float
					     (((complexe_float_t *) X)[q *
								       incX],
					      Abis[N * p + q]), ((complexe_float_t *)Y)[p * incY]);
				}
			}
		} else {
			for (int p = 0; p < Mbis; p++) {
				for (int q = 0; q < Nbis; q++) {
					((complexe_float_t *)Y)[p * incY] =
					    add_complexe_float
					    (mult_complexe_float
					     (((complexe_float_t *) X)[q *
								       incX],
					      Abis[p + M * q]), ((complexe_float_t *)Y)[p * incY]);
				}
			}
		}
		mncblas_caxpy(Mbis, beta, copyY, incY, Y, incY);
	}
	return;
}

void mncblas_zgemv(MNCBLAS_LAYOUT layout,
		   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
		   const void *alpha, const void *A, const int lda,
		   const void *X, const int incX, const void *beta,
		   void *Y, const int incY) {
	int Nbis = N;
	int Mbis = M;
		complexe_double_t *Abis = (complexe_double_t *) malloc(N * M *
						sizeof(complexe_double_t));
	if (TransA != MNCblasNoTrans) {
		Nbis = M;
		Mbis = N;
	}
	int lenY = (1 + (Mbis - 1) * incY);
	if ((M == 0) || (N == 0)) {	//si la mtrice est nulle on a juste y=beta*y
		for (int i = 0; i < lenY; i += incY) {
			((complexe_double_t *) Y)[i] =
			    mult_complexe_double(*((complexe_double_t *) beta),
						((complexe_double_t *) Y)[i]);
		}
	} else {
		complexe_double_t *copyY = copy_double_complexe (Nbis, Y) ;
		if (TransA != MNCblasNoTrans) {	//On fait la transposée si demandée
			Abis = transposee_z(layout, M, N, A);
		} else {
			for (int k = 0; k < N * M; k++) {	//sinon abis est une copie de A
				Abis[k] = ((complexe_double_t *) A)[k];
			}
		}
		if (TransA == MNCblasConjTrans) {
			for (int incr = 0; incr < N * M; incr++) {
				Abis[incr] = conjg_double(Abis[incr]);
			}
		}
		for (int j = 0; j < M * N; j++) {	//calcul de alpha*matrice
			Abis[j] =
			    mult_complexe_double(*((complexe_double_t *) alpha),
						((complexe_double_t *) Abis)[j]);
		}
		if (layout == MNCblasRowMajor) {
			for (int p = 0; p < Mbis; p++) {
				for (int q = 0; q < Nbis; q++) {
					((complexe_double_t *)Y)[p * incY] =
					    add_complexe_double
					    (mult_complexe_double
					     (((complexe_double_t *) X)[q *
								       incX],
					      Abis[N * p + q]), ((complexe_double_t *)Y)[p * incY]);
				}
			}
		} else {
			for (int p = 0; p < Mbis; p++) {
				for (int q = 0; q < Nbis; q++) {
					((complexe_double_t *)Y)[p * incY] =
					    add_complexe_double
					    (mult_complexe_double
					     (((complexe_double_t *) X)[q *
								       incX],
					      Abis[p + M * q]), ((complexe_double_t *)Y)[p * incY]);
				}
			}
		}
		mncblas_zaxpy(Mbis, beta, copyY, incY, Y, incY);
	}
	return;
}



