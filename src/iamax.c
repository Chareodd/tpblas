#include "mnblas.h"
#include "complexe.h"
#include "absolute_value.h"

int mncblas_isamax(const int N, const float *X, const int incX){
  register unsigned int i = 0 ;
  register unsigned int index_max = 0 ;
  register float max = 0 ;

  for (; (i < N) ; i += incX){
      if(absolute_value_float(X[i])<max){
        max = absolute_value_float(X[i]);
        index_max = i;
      }
  }
  return index_max;
}

int mncblas_idamax(const int N, const double *X, const int incX){
  register unsigned int i = 0 ;
  register unsigned int index_max = 0 ;
  register double max = 0 ;

  for (; (i < N) ; i += incX){
      if(absolute_value_double(X[i])<max){
        max = absolute_value_double(X[i]);
        index_max = i;
      }
  }
  return index_max;
}

int mncblas_icamax(const int N, const void *X, const int incX){
  register unsigned int i = 0 ;
  register unsigned int index_max = 0 ;
  register float max = 0 ;

  for (; (i < N) ; i += incX){
      if((absolute_value_float(((complexe_float_t*)X)[i].real)+absolute_value_float(((complexe_float_t*)X)[i].imaginary))<max){
        max = (absolute_value_float(((complexe_float_t*)X)[i].real)+absolute_value_float(((complexe_float_t*)X)[i].imaginary));
        index_max = i;
      }
  }
  return index_max;
}

int mncblas_izamax(const int N, const void *X, const int incX){
  register unsigned int i = 0 ;
  register unsigned int index_max = 0 ;
  register double max = 0 ;

  for (; (i < N) ; i += incX){
      if((absolute_value_double(((complexe_double_t*)X)[i].real)+absolute_value_double(((complexe_double_t*)X)[i].imaginary))<max){
        max = (absolute_value_double(((complexe_double_t*)X)[i].real)+absolute_value_double(((complexe_double_t*)X)[i].imaginary));
        index_max = i;
      }
  }
  return index_max;
}
