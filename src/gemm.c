#include "aux_blas3.h"

#include <stdio.h>
#include <stdlib.h>


void mncblas_sgemm (const MNCBLAS_LAYOUT Layout, const MNCBLAS_TRANSPOSE transa, 
        const MNCBLAS_TRANSPOSE transb, const int m, const int n, const int k, 
        const float alpha, const float *a, const float *b, const float beta, 
		float *c)
{
	float* opA, *opB ;
	//Création de op(A)
	if (transa == MNCblasNoTrans) {
		opA = copie_float (Layout, m, k, a) ;
	} else if (transa == MNCblasTrans) {
		opA = transposee_float (Layout, m, k, a) ;
	} else {
		float* tempA = transposee_float (Layout, m, k, a) ;
		opA = conjuguee_float (Layout, m, k, tempA) ;
		free (tempA) ;
	}
	//Création de op(B)
	if (transb == MNCblasNoTrans) {
		opB = copie_float (Layout, k, n, b) ;
	} else if (transb == MNCblasTrans) {
		opB = transposee_float (Layout, k, n, b) ;
	} else {
		float* tempB = transposee_float (Layout, k, n, b) ;
		opB = conjuguee_float (Layout, k, n, tempB) ;
		free (tempB) ;
	}

	//Calcul de C
	float* copyC = copie_float (Layout, k, n, c) ;
	if (Layout == MNCblasRowMajor) {
		for (int i = 0 ; i < m ; i++) {
			for (int j = 0 ; j < n ; j++) {
				float sum = 0 ;
				for (int l = 0 ; l < k ; l++) {
					sum += opA[k * i + l] * opB[n * l + j] ;
				}
				c[n * i + j] = alpha * sum + beta * copyC[n * i + j] ;
			}
		}
	} else {	//Layout == MNCblasColMajor
		for (int i = 0 ; i < n ; i++) {
			for (int j = 0 ; j < m ; j++) {
				float sum = 0 ;
				for (int l = 0 ; l < k ; l++) {
					sum += opA[m * l + j] * opB[k * i + l] ;
				}
				c[m * i + j] = alpha * sum + beta * copyC[m * i + j] ;
			}
		}
	}
	free (opA) ;
	free (opB) ;
	free (copyC) ;
}

void mncblas_dgemm (const MNCBLAS_LAYOUT Layout, const MNCBLAS_TRANSPOSE transa, 
        const MNCBLAS_TRANSPOSE transb, const int m, const int n, const int k, 
        const double alpha, const double *a, const double *b, const double beta, 
		double *c)
{
	double* opA, *opB ;
	//Création de op(A)
	if (transa == MNCblasNoTrans) {
		opA = copie_double (Layout, m, k, a) ;
	} else if (transa == MNCblasTrans) {
		opA = transposee_double (Layout, m, k, a) ;
	} else {
		double* tempA = transposee_double (Layout, m, k, a) ;
		opA = conjuguee_double (Layout, m, k, tempA) ;
		free (tempA) ;
	}
	//Création de op(B)
	if (transb == MNCblasNoTrans) {
		opB = copie_double (Layout, k, n, b) ;
	} else if (transb == MNCblasTrans) {
		opB = transposee_double (Layout, k, n, b) ;
	} else {
		double* tempB = transposee_double (Layout, k, n, b) ;
		opB = conjuguee_double (Layout, k, n, tempB) ;
		free (tempB) ;
	}

	//Calcul de C
	double* copyC = copie_double (Layout, k, n, c) ;
	if (Layout == MNCblasRowMajor) {
		for (int i = 0 ; i < m ; i++) {
			for (int j = 0 ; j < n ; j++) {
				double sum = 0 ;
				for (int l = 0 ; l < k ; l++) {
					sum += opA[k * i + l] * opB[n * l + j] ;
				}
				c[n * i + j] = alpha * sum + beta * copyC[n * i + j] ;
			}
		}
	} else {	//Layout == MNCblasColMajor
		for (int i = 0 ; i < n ; i++) {
			for (int j = 0 ; j < m ; j++) {
				float sum = 0 ;
				for (int l = 0 ; l < k ; l++) {
					sum += opA[m * l + j] * opB[k * i + l] ;
				}
				c[m * i + j] = alpha * sum + beta * copyC[m * i + j] ;
			}
		}
	}
	free (opA) ;
	free (opB) ;
	free (copyC) ;
}

void mncblas_cgemm (const MNCBLAS_LAYOUT Layout, const MNCBLAS_TRANSPOSE transa, 
        const MNCBLAS_TRANSPOSE transb, const int m, const int n, const int k, 
        const void* alpha, const void *a, const void *b, const void* beta, 
		void *c)
{
	complexe_float_t* opA, *opB ;
	//Création de op(A)
	if (transa == MNCblasNoTrans) {
		opA = copie_float_complexe (Layout, m, k, a) ;
	} else if (transa == MNCblasTrans) {
		opA = transposee_float_complexe (Layout, m, k, a) ;
	} else {
		complexe_float_t* tempA = transposee_float_complexe (Layout, m, k, a) ;
		opA = conjuguee_float (Layout, m, k, tempA) ;
		free (tempA) ;
	}
	//Création de op(B)
	if (transb == MNCblasNoTrans) {
		opB = copie_float_complexe (Layout, k, n, b) ;
	} else if (transb == MNCblasTrans) {
		opB = transposee_float_complexe (Layout, k, n, b) ;
	} else {
		complexe_float_t* tempB = transposee_float_complexe (Layout, k, n, b) ;
		opB = conjuguee_float (Layout, k, n, tempB) ;
		free (tempB) ;
	}

	//Calcul de C
	complexe_float_t* copyC = copie_float_complexe (Layout, k, n, c) ;
	if (Layout == MNCblasRowMajor) {
		for (int i = 0 ; i < m ; i++) {
			for (int j = 0 ; j < n ; j++) {
				complexe_float_t sum ;
				sum.real = 0 ;
				sum.imaginary = 0 ;
				for (int l = 0 ; l < k ; l++) {
					sum = add_complexe_float (mult_complexe_float (opA[k * i + l], opB[n * l + j]), sum) ;
				}
				((complexe_float_t*)c)[n * i + j] = add_complexe_float (mult_complexe_float (*((complexe_float_t*)alpha), sum), mult_complexe_float (*((complexe_float_t*)beta), copyC[n * i + j])) ;
			}
		}
	} else {	//Layout == MNCblasColMajor
		for (int i = 0 ; i < n ; i++) {
			for (int j = 0 ; j < m ; j++) {
				complexe_float_t sum ;
				sum.real = 0 ;
				sum.imaginary = 0 ;
				for (int l = 0 ; l < k ; l++) {
					sum = add_complexe_float (mult_complexe_float (opA[m * l + j], opB[k * i + l]), sum) ;
				}
				((complexe_float_t*)c)[m * i + j] = add_complexe_float (mult_complexe_float (*((complexe_float_t*)alpha), sum), mult_complexe_float (*((complexe_float_t*)beta), copyC[m * i + j])) ;
			}
		}
	}
	free (opA) ;
	free (opB) ;
	free (copyC) ;
}

void mncblas_zgemm (const MNCBLAS_LAYOUT Layout, const MNCBLAS_TRANSPOSE transa, 
        const MNCBLAS_TRANSPOSE transb, const int m, const int n, const int k, 
        const void* alpha, const void *a, const void *b, const void* beta, 
		void *c)
{
	complexe_double_t* opA, *opB ;
	//Création de op(A)
	if (transa == MNCblasNoTrans) {
		opA = copie_double_complexe (Layout, m, k, a) ;
	} else if (transa == MNCblasTrans) {
		opA = transposee_double_complexe (Layout, m, k, a) ;
	} else {
		complexe_double_t* tempA = transposee_double_complexe (Layout, m, k, a) ;
		opA = conjuguee_double (Layout, m, k, tempA) ;
		free (tempA) ;
	}
	//Création de op(B)
	if (transb == MNCblasNoTrans) {
		opB = copie_double_complexe (Layout, k, n, b) ;
	} else if (transb == MNCblasTrans) {
		opB = transposee_double_complexe (Layout, k, n, b) ;
	} else {
		complexe_double_t* tempB = transposee_double_complexe (Layout, k, n, b) ;
		opB = conjuguee_double (Layout, k, n, tempB) ;
		free (tempB) ;
	}

	//Calcul de C
	complexe_double_t* copyC = copie_double_complexe (Layout, k, n, c) ;
	if (Layout == MNCblasRowMajor) {
		for (int i = 0 ; i < m ; i++) {
			for (int j = 0 ; j < n ; j++) {
				complexe_double_t sum ;
				sum.real = 0 ;
				sum.imaginary = 0 ;
				for (int l = 0 ; l < k ; l++) {
					sum = add_complexe_double (mult_complexe_double (opA[k * i + l], opB[n * l + j]), sum) ;
				}
				((complexe_double_t*)c)[n * i + j] = add_complexe_double (mult_complexe_double (*((complexe_double_t*)alpha), sum), mult_complexe_double (*((complexe_double_t*)beta), copyC[n * i + j])) ;
			}
		}
	} else {	//Layout == MNCblasColMajor
		for (int i = 0 ; i < n ; i++) {
			for (int j = 0 ; j < m ; j++) {
				complexe_double_t sum ;
				sum.real = 0 ;
				sum.imaginary = 0 ;
				for (int l = 0 ; l < k ; l++) {
					sum = add_complexe_double (mult_complexe_double (opA[m * l + j], opB[k * i + l]), sum) ;
				}
				((complexe_double_t*)c)[m * i + j] = add_complexe_double (mult_complexe_double (*((complexe_double_t*)alpha), sum), mult_complexe_double (*((complexe_double_t*)beta), copyC[m * i + j])) ;
			}
		}
	}
	free (opA) ;
	free (opB) ;
	free (copyC) ;
}