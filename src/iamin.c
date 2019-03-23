#include "mnblas.h"
#include "complexe.h"
#include "absolute_value.h"

int mncblas_isamin(const int N, const float *X, const int incX){
  register unsigned int i = 0 ;
  register unsigned int index_min = 0 ;
  register float min = 0 ;

  for (; (i < N) ; i += incX){
      if(absolute_value_float(X[i])<min){
        min = absolute_value_float(X[i]);
        index_min = i;
      }
  }
  return index_min;
}

int mncblas_idamin(const int N, const double *X, const int incX){
  register unsigned int i = 0 ;
  register unsigned int index_min = 0 ;
  register double min = 0 ;

  for (; (i < N) ; i += incX){
      if(absolute_value_double(X[i])<min){
        min = absolute_value_double(X[i]);
        index_min = i;
      }
  }
  return index_min;
}

int mncblas_icamin(const int N, const void *X, const int incX){
  register unsigned int i = 0 ;
  register unsigned int index_min = 0 ;
  register float min = 0 ;

  for (; (i < N) ; i += incX){
      if((absolute_value_float(((complexe_float_t*)X)[i].real)+absolute_value_float(((complexe_float_t*)X)[i].imaginary))<min){
        min = (absolute_value_float(((complexe_float_t*)X)[i].real)+absolute_value_float(((complexe_float_t*)X)[i].imaginary));
        index_min = i;
      }
  }
  return index_min;
}

int mncblas_izamin(const int N, const void *X, const int incX){
  register unsigned int i = 0 ;
  register unsigned int index_min = 0 ;
  register double min = 0 ;

  for (; (i < N) ; i += incX){
      if((absolute_value_double(((complexe_double_t*)X)[i].real)+absolute_value_double(((complexe_double_t*)X)[i].imaginary))<min){
        min = (absolute_value_double(((complexe_double_t*)X)[i].real)+absolute_value_double(((complexe_double_t*)X)[i].imaginary));
        index_min = i;
      }
  }
  return index_min;
}

