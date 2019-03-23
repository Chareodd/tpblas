#include "complexe.h"
#include "mnblas.h"

float* copie_float (const MNCBLAS_LAYOUT layout, const int m, const int n, const float *a) ;
double* copie_double (const MNCBLAS_LAYOUT layout, const int m, const int n, const double *a) ;
complexe_float_t* copie_float_complexe (const MNCBLAS_LAYOUT layout, const int m, const int n, const complexe_float_t *a) ;
complexe_double_t* copie_double_complexe (const MNCBLAS_LAYOUT layout, const int m, const int n, const complexe_double_t *a) ;

void* conjuguee_float (const MNCBLAS_LAYOUT layout, const int m, const int n, const void *a) ;
void* conjuguee_double (const MNCBLAS_LAYOUT layout, const int m, const int n, const void *a) ;

float* transposee_float(const MNCBLAS_LAYOUT layout, const int m, const int n, const float *a) ;
double* transposee_double(const MNCBLAS_LAYOUT layout, const int m, const int n, const double *a) ;
complexe_float_t* transposee_float_complexe(const MNCBLAS_LAYOUT layout, const int m, const int n, const complexe_float_t *a) ;
complexe_double_t* transposee_double_complexe(const MNCBLAS_LAYOUT layout, const int m, const int n, const complexe_double_t *a) ;