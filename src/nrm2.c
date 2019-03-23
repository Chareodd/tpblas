#include "mnblas.h"
#include "complexe.h"
#include <math.h>
#include <stdlib.h>

float mnblas_snrm2 (const int n, const float *x, const int incx) {
    float res = 0;
    for (int i = 0; i<n; i+=incx) {
        res += x[i]*x[i];
    }
    res=sqrt(res);
    return res;
}

double mnblas_dnrm2 (const int n, const double *x, const int  incx) {
    double res = 0;
    for (int i = 0; i<n; i+=incx) {
        res += x[i]*x[i];
    }
    res=sqrt(res);
    return res;
}

float mnblas_scnrm2 (const int  n, const void *x, const int incx) {
    float* res = (float*) malloc (sizeof(float)) ;
    mncblas_cdotu_sub(n, x, incx, x, incx, res) ;
    return sqrt(*res);
}

double mnblas_dznrm2 (const int n, const void *x, const int incx){
    float* res = (float*) malloc (sizeof(float)) ;
    mncblas_zdotu_sub(n, x, incx, x, incx, res) ;
    return sqrt(*res);
}
