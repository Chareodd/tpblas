#include "mnblas.h"
#include "complexe.h"

#include <stdlib.h>

float f_abs (float f) {
    if (f < 0) {
        return -f ;
    } else {
        return f ;
    }
}

double d_abs (double d) {
    if (d < 0) {
        return -d ;
    } else {
        return d ;
    }
}

float mnblas_sasum (const int n, const float *x, const int incx) {
    float sum = 0 ;
    for (int i = 0 ; i < n ; i++) {
        sum += f_abs(x[i]) ;
    }
    return sum ;
}

float mnblas_scasum (const int n, const void *x, const int incx) {
    float sum = 0 ;
    for (int i = 0 ; i < n ; i++) {
        sum += f_abs(((complexe_float_t*)x)[i].real) + f_abs(((complexe_float_t*)x)[i].imaginary) ;
    }
    return sum ;
}

double mnblas_dasum (const int n, const double *x, const int incx) {
    double sum = 0 ;
    for (int i = 0 ; i < n ; i++) {
        sum += d_abs(x[i]) ;
    }
    return sum ;
}

double mnblas_dzasum (const int n, const void *x, const int incx) {
    double sum = 0 ;
    for (int i = 0 ; i < n ; i++) {
        sum += d_abs(((complexe_double_t*)x)[i].real) + d_abs(((complexe_double_t*)x)[i].imaginary) ;
    }
    return sum ;
}