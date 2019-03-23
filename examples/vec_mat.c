#include <stdio.h>
#include <stdlib.h>

#include "vec_mat.h"

//VECTEUR

vfloat vector_init_float (int size, float x)
{
  register unsigned int i ;
  vfloat res = (vfloat) malloc (size * sizeof(float)) ;
  for (i = 0; i < size; i++)
    res [i] = x ;

  return res ;
}

vdouble vector_init_double (int size, double x)
{
  register unsigned int i ;
  vdouble res = (vdouble) malloc (size * sizeof(double)) ;
  for (i = 0; i < size; i++)
    res [i] = x ;

  return res ;
}

vcfloat vector_init_float_complexe (int size, complexe_float_t x)
{
  register unsigned int i ;
  vcfloat res = (vcfloat) malloc (size * sizeof(complexe_float_t)) ;
  for (i = 0; i < size; i++)
    res [i] = x ;

  return res ;
}

vcdouble vector_init_double_complexe (int size, complexe_double_t x)
{
  register unsigned int i ;
  vcdouble res = (vcdouble) malloc (size * sizeof(complexe_double_t)) ;
  for (i = 0; i < size; i++)
    res [i] = x ;

  return res ;
}

void vector_print_float (int size, vfloat V)
{
  register unsigned int i ;

  for (i = 0; i < size; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;

  return ;
}

void vector_print_double (int size, vdouble V)
{
  register unsigned int i ;

  for (i = 0; i < size; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;

  return ;
}

void vector_print_float_complexe (int size, vcfloat V)
{
  register unsigned int i ;

  for (i = 0; i < size; i++)
    printf ("%f + i%f  ", V[i].real, V[i].imaginary) ;
  printf ("\n") ;

  return ;
}

void vector_print_double_complexe (int size, vcdouble V)
{
  register unsigned int i ;

  for (i = 0; i < size; i++)
    printf ("%f + i%f  ", V[i].real, V[i].imaginary) ;
  printf ("\n") ;

  return ;
}

//MATRICE

mfloat matrix_init_float (int row, int col, float x) 
{
  register unsigned int i, j ;
  mfloat res = (mfloat) malloc (row * col * sizeof(float)) ;
  for (i = 0; i < row; i++) {
      for (j = 0; j < col; j++) {
          res [col * i + j] = x ;
      }
  }

  return res ;
}

mdouble matrix_init_double (int row, int col, double x) 
{
  register unsigned int i, j ;
  mdouble res = (mdouble) malloc (row * col * sizeof(double)) ;
  for (i = 0; i < row; i++) {
      for (j = 0; j < col; j++) {
          res [col * i + j] = x ;
      }
  }

  return res ;
}

mcfloat matrix_init_float_complexe (int row, int col, complexe_float_t x) 
{
  register unsigned int i, j ;
  mcfloat res = (mcfloat) malloc (row * col * sizeof(complexe_float_t)) ;
  for (i = 0; i < row; i++) {
      for (j = 0; j < col; j++) {
          res [col * i + j] = x ;
      }
  }

  return res ;
}

mcdouble matrix_init_double_complexe (int row, int col, complexe_double_t x)
{
  register unsigned int i, j ;
  mcdouble res = (mcdouble) malloc (row * col * sizeof(complexe_float_t)) ;
  for (i = 0; i < row; i++) {
      for (j = 0; j < col; j++) {
          res [col * i + j] = x ;
      }
  }

  return res ;
}

void matrix_print_float (MNCBLAS_LAYOUT Layout, int row, int col, mfloat M) 
{
  register unsigned int i, j ;

  printf ("\n") ;
  if (Layout == MNCblasRowMajor) {
    for (i = 0 ; i < row ; i++) {
      for (j = 0 ; j < col ; j++) {
        printf ("%f  ", M[i * col + j]) ;
      }
      printf ("\n") ;
    }
  } else {
    for (i = 0 ; i < col ; i++) {
      for (j = 0 ; j < row ; j++) {
        printf ("%f  ", M[i * row + j]) ;
      }
      printf ("\n") ;
    }
  }

  return ;
}

void matrix_print_double (MNCBLAS_LAYOUT Layout, int row, int col, mdouble M)
{
  register unsigned int i, j ;

  printf ("\n") ;
  if (Layout == MNCblasRowMajor) {
    for (i = 0 ; i < row ; i++) {
      for (j = 0 ; j < col ; j++) {
        printf ("%f  ", M[i * col + j]) ;
      }
      printf ("\n") ;
    }
  } else {
    for (i = 0 ; i < col ; i++) {
      for (j = 0 ; j < row ; j++) {
        printf ("%f  ", M[i * row + j]) ;
      }
      printf ("\n") ;
    }
  }


  return ;
}

void matrix_print_float_complexe (MNCBLAS_LAYOUT Layout, int row, int col, mcfloat M)
{
  register unsigned int i, j ;

  printf ("\n") ;
  if (Layout == MNCblasRowMajor) {
    for (i = 0 ; i < row ; i++) {
      for (j = 0 ; j < col ; j++) {
        printf ("%f + i%f  ", M[i * col + j].real, M[i * col + j].imaginary) ;
      }
      printf ("\n") ;
    }
  } else {
    for (i = 0 ; i < col ; i++) {
      for (j = 0 ; j < row ; j++) {
        printf ("%f + i%f  ", M[i * col + j].real, M[i * col + j].imaginary) ;
      }
      printf ("\n") ;
    }
  }


  return ;
}

void matrix_print_double_complexe (MNCBLAS_LAYOUT Layout, int row, int col, mcdouble M)
{
  register unsigned int i, j ;

  printf ("\n") ;
  if (Layout == MNCblasRowMajor) {
    for (i = 0 ; i < row ; i++) {
      for (j = 0 ; j < col ; j++) {
        printf ("%f + i%f  ", M[i * col + j].real, M[i * col + j].imaginary) ;
      }
      printf ("\n") ;
    }
  } else {
    for (i = 0 ; i < col ; i++) {
      for (j = 0 ; j < row ; j++) {
        printf ("%f + i%f  ", M[i * col + j].real, M[i * col + j].imaginary) ;
      }
      printf ("\n") ;
    }
  }


  return ;
}
