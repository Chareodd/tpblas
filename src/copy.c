#include "mnblas.h"
#include "complexe.h"

void mncblas_scopy(const int N, const float *X, const int incX,
                 float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      Y [j] = X [i] ;
    }

  return ;
}

void mncblas_dcopy(const int N, const double *X, const int incX,
                 double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY) {
      Y [j] = X [i] ;
  }

  return ;
}

void mncblas_ccopy(const int N, const void *X, const int incX,
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY) {
      ((complexe_float_t *)Y)[j].real = ((complexe_float_t *)X)[i].real ;
      ((complexe_float_t *)Y)[j].imaginary = ((complexe_float_t *)X)[i].imaginary ;
  }

  return ;
}

void mncblas_zcopy(const int N, const void *X, const int incX,
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY) {
      ((complexe_double_t *)Y)[j].real = ((complexe_double_t *)X)[i].real ;
      ((complexe_double_t *)Y)[j].imaginary = ((complexe_double_t *)X)[i].imaginary ;
  }

  return ;
}
