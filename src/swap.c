#include "mnblas.h"
#include "complexe.h"

void mncblas_sswap(const int N, float *X, const int incX,
                 float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float save ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_dswap(const int N, double *X, const int incX,
                 double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register double save ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_cswap(const int N, void *X, const int incX,
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register complexe_float_t * Xc = X ;
  register complexe_float_t * Yc = Y ;
  register float saver;
  register float savei;
  for (; ((i < N) && (j < N)) ; i += incX, j += incY) {
      saver=Yc[j].real;
      savei=Yc[j].imaginary;
      Yc[j].real = Xc[i].real ;
      Yc[j].imaginary = Xc[i].imaginary ;
      Xc[i].real=saver;
      Xc[i].imaginary=savei;
  }

  return ;
}

void mncblas_zswap(const int N, void *X, const int incX,
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register complexe_double_t * Xc = X ;
  register complexe_double_t * Yc = Y ;
  register double saver;
  register double savei;
  for (; ((i < N) && (j < N)) ; i += incX, j += incY) {
      saver=Yc[j].real;
      savei=Yc[j].imaginary;
      Yc[j].real = Xc[i].real ;
      Yc[j].imaginary = Xc[i].imaginary ;
      Xc[i].real=saver;
      Xc[i].imaginary=savei;
  }

  return ;
}
