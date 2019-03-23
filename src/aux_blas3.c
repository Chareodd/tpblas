#include "aux_blas3.h"

#include <stdio.h>
#include <stdlib.h>

float* copie_float (const MNCBLAS_LAYOUT layout, const int m, const int n, const float *a) {
	float* abis = (float*) malloc (n * m * sizeof(float)) ;

	if (layout == MNCblasRowMajor) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				abis[n * i + j] = a[n * i + j] ;
			}
		}
	} else {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				abis[m * i + j] = a[m * i + j] ;
			}
		}
	}

	return abis ;
}

double* copie_double (const MNCBLAS_LAYOUT layout, const int m, const int n, const double *a) {
	double* abis = (double*) malloc (n * m * sizeof(double)) ;

	if (layout == MNCblasRowMajor) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				abis[n * i + j] = a[n * i + j] ;
			}
		}
	} else {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				abis[m * i + j] = a[m * i + j] ;
			}
		}
	}

	return abis ;
}

complexe_float_t* copie_float_complexe (const MNCBLAS_LAYOUT layout, const int m, const int n, const complexe_float_t *a) {
	complexe_float_t* abis = (complexe_float_t*) malloc (n * m * sizeof(complexe_float_t)) ;

	if (layout == MNCblasRowMajor) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				abis[n * i + j] = a[n * i + j] ;
			}
		}
	} else {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				abis[m * i + j] = a[m * i + j] ;
			}
		}
	}

	return abis ;
}

complexe_double_t* copie_double_complexe (const MNCBLAS_LAYOUT layout, const int m, const int n, const complexe_double_t *a) {
	complexe_double_t* abis = (complexe_double_t*) malloc (n * m * sizeof(complexe_double_t)) ;

	if (layout == MNCblasRowMajor) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				abis[n * i + j] = a[n * i + j] ;
			}
		}
	} else {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				abis[m * i + j] = a[m * i + j] ;
			}
		}
	}

	return abis ;
}

void* conjuguee_float (const MNCBLAS_LAYOUT layout, const int m, const int n, const void *a) {
	complexe_float_t* abis = (complexe_float_t*) malloc (n * m * sizeof(complexe_float_t)) ;

	if (layout == MNCblasRowMajor) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				abis[n * i + j] = conjg_float (((complexe_float_t*)a)[n * i + j]) ;
			}
		}
	} else {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				abis[m * i + j] = conjg_float (((complexe_float_t*)a)[m * i + j]) ;
			}
		}
	}

    return (void*) abis ;
}

void* conjuguee_double (const MNCBLAS_LAYOUT layout, const int m, const int n, const void *a) {
	complexe_double_t* abis = (complexe_double_t*) malloc (n * m * sizeof(complexe_double_t)) ;

	if (layout == MNCblasRowMajor) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				abis[n * i + j] = conjg_double (((complexe_double_t*)a)[n * i + j]) ;
			}
		}
	} else {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				abis[m * i + j] = conjg_double (((complexe_double_t*)a)[m * i + j]) ;
			}
		}
	}

    return (void*) abis ;
}

float* transposee_float(const MNCBLAS_LAYOUT layout, const int m, const int n,
	   const float *a)
{
	float* abis = (float*) malloc (n * m * sizeof(float)) ;

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

    return abis ;
}

double* transposee_double(const MNCBLAS_LAYOUT layout, const int m, const int n,
	   const double *a)
{
	double* abis = (double*) malloc (n * m * sizeof(double)) ;

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

    return abis ;
}

complexe_float_t* transposee_float_complexe(const MNCBLAS_LAYOUT layout, const int m, const int n,
	   const complexe_float_t *a)
{
	complexe_float_t* abis = (complexe_float_t*) malloc (n * m * sizeof(complexe_float_t)) ;

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

    return abis ;
}

complexe_double_t* transposee_double_complexe(const MNCBLAS_LAYOUT layout, const int m, const int n,
	   const complexe_double_t *a)
{
	complexe_double_t* abis = (complexe_double_t*) malloc (n * m * sizeof(complexe_double_t)) ;

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

    return abis ;
}
