#include "mnblas.h"
#include "complexe.h"

#include <stdlib.h>

float mncblas_sdot(const int N, const float *X, const int incX,
                 const float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float dot = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      dot = dot + X [i] * Y [j] ;
    }

  return dot ;
}

double mncblas_ddot(const int N, const double *X, const int incX,
                 const double *Y, const int incY)
{
  register unsigned int i = 0;
  register unsigned int j = 0 ;
  register double dot = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY) {
      dot += X [i] * Y [j] ;
  }

  return dot;
}

void   mncblas_cdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  ((complexe_float_t*)dotu)->real = 0 ;
  ((complexe_float_t*)dotu)->imaginary = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY) {
      *((complexe_float_t*)dotu) = add_complexe_float( *((complexe_float_t*)dotu), mult_complexe_float( ((complexe_float_t*)X)[i], ((complexe_float_t*)Y)[j])) ;
  }
}

void   mncblas_cdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0;
  register unsigned int j = 0 ;
  ((complexe_float_t*)dotc)->real = 0 ;
  ((complexe_float_t*)dotc)->imaginary = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY) {
      *((complexe_float_t*)dotc) = add_complexe_float( *((complexe_float_t*)dotc), mult_complexe_float( conjg_float(((complexe_float_t*)X)[i]), ((complexe_float_t*)Y)[j])) ;
  }
}

void   mncblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  ((complexe_double_t*)dotu)->real = 0 ;
  ((complexe_double_t*)dotu)->imaginary = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY) {
      *((complexe_double_t*)dotu) = add_complexe_double( *((complexe_double_t*)dotu), mult_complexe_double( ((complexe_double_t*)X)[i], ((complexe_double_t*)Y)[j])) ;
  }

}

void   mncblas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0;
  register unsigned int j = 0 ;
  ((complexe_double_t*)dotc)->real = 0 ;
  ((complexe_double_t*)dotc)->imaginary = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY) {
      *((complexe_double_t*)dotc) = add_complexe_double( *((complexe_double_t*)dotc), mult_complexe_double( conjg_double(((complexe_double_t*)X)[i]), ((complexe_double_t*)Y)[j])) ;
  }
}
