#include <stdio.h>
#include <x86intrin.h>

#include "vec_mat.h"

#include "flop.h"

#define NB_FOIS     1
#define ROW         1024
#define COL         1024
#define DISP_MAT    0

mfloat a1, b1, c1 ;
mdouble a2, b2, c2 ;
mcfloat a3, b3, c3 ;
mcdouble a4, b4, c4 ;

int main (int argc, char **argv)
{
 unsigned long long start, end ;
 int i ;
 int k = ROW ;
 MNCBLAS_LAYOUT Layout = MNCblasRowMajor ;
 MNCBLAS_TRANSPOSE transa = MNCblasNoTrans ;
 MNCBLAS_TRANSPOSE transb = MNCblasNoTrans ;  
 complexe_float_t init_float_complexe0 ;
 init_float_complexe0.real = 0 ;
 init_float_complexe0.imaginary = 0 ;
 complexe_float_t init_float_complexe1 ;
 init_float_complexe1.real = 1 ;
 init_float_complexe1.imaginary = 1 ;
 complexe_float_t init_float_complexe2 ;
 init_float_complexe2.real = 2 ;
 init_float_complexe2.imaginary = 2 ;
 complexe_double_t init_double_complexe0 ;
 init_double_complexe0.real = 0 ;
 init_double_complexe0.imaginary = 0 ;
 complexe_double_t init_double_complexe1 ;
 init_double_complexe1.real = 1 ;
 init_double_complexe1.imaginary = 1 ;
 complexe_double_t init_double_complexe2 ;
 init_double_complexe2.real = 2 ;
 init_double_complexe2.imaginary = 2 ;
 complexe_float_t *alpha_float = (complexe_float_t*) malloc(sizeof(complexe_float_t)) ;
 alpha_float->real = 1 ;
 alpha_float->imaginary = 0 ;
 complexe_double_t *alpha_double = (complexe_double_t*) malloc(sizeof(complexe_double_t)) ;
 alpha_double->real = 1 ;
 alpha_double->imaginary = 0 ;
 complexe_float_t *beta_float = (complexe_float_t*) malloc(sizeof(complexe_float_t)) ;
 beta_float->real = 1 ;
 beta_float->imaginary = 0 ;
 complexe_double_t *beta_double = (complexe_double_t*) malloc(sizeof(complexe_double_t)) ;
 beta_double->real = 1 ;
 beta_double->imaginary = 0 ;

 for (i = 0 ; i < NB_FOIS; i++)
   {
     a1 = matrix_init_float (ROW, COL, 1.0) ;
     b1 = matrix_init_float (ROW, COL, 2.0) ;
     c1 = matrix_init_float (ROW, COL, 0) ;
     a2 = matrix_init_double (ROW, COL, 1.0) ;
     b2 = matrix_init_double (ROW, COL, 2.0) ;
     c2 = matrix_init_double (ROW, COL, 0) ;
     a3 = matrix_init_float_complexe (ROW, COL, init_float_complexe1) ;
     b3 = matrix_init_float_complexe (ROW, COL, init_float_complexe2) ;
     c3 = matrix_init_float_complexe (ROW, COL, init_float_complexe0) ;
     a4 = matrix_init_double_complexe (ROW, COL, init_double_complexe1) ;
     b4 = matrix_init_double_complexe (ROW, COL, init_double_complexe2) ;
     c4 = matrix_init_double_complexe (ROW, COL, init_double_complexe0) ;

     //Produit/somme matrice-matrice float
     if (DISP_MAT) {
         printf("a1 = ") ;
         matrix_print_float(Layout, ROW, COL, a1) ;
         printf("b1 = ") ;
         matrix_print_float(Layout, ROW, COL, b1) ;
     }
     start = _rdtsc () ;
         mncblas_sgemm (Layout, transa, transb, ROW, COL, k, 1, a1, b1, 1, c1) ;
     end = _rdtsc () ;
     if (DISP_MAT) {
         printf("c1 = ") ;
         matrix_print_float(Layout, ROW, COL, c1) ;
     }
     printf ("Nombre de cycles: %Ld \n", end-start) ;
     calcul_flop ("sgemm ", ROW * COL * k, end-start) ;

     //Produit/somme matrice-matrice double
     if (DISP_MAT) {
         printf("a2 = ") ;
         matrix_print_double(Layout, ROW, COL, a2) ;
         printf("b2 = ") ;
         matrix_print_double(Layout, ROW, COL, b2) ;
     }
     start = _rdtsc () ;
         mncblas_dgemm (Layout, transa, transb, ROW, COL, k, 1, a2, b2, 1, c2) ;
     end = _rdtsc () ;
     if (DISP_MAT) {
         printf("c2 = ") ;
         matrix_print_double(Layout, ROW, COL, c2) ;
     }
     printf ("Nombre de cycles: %Ld \n", end-start) ;
     calcul_flop ("dgemm ", ROW * COL * k, end-start) ;

     ////Produit/somme matrice-matrice float complexe
     if (DISP_MAT) {
         printf("a3 = ") ;
         matrix_print_float_complexe(Layout, ROW, COL, a3) ;
         printf("b3 = ") ;
         matrix_print_float_complexe(Layout, ROW, COL, b3) ;
     }
     start = _rdtsc () ;
        mncblas_cgemm (Layout, transa, transb, ROW, COL, k, alpha_float, a3, b3, beta_float, c3) ;
     end = _rdtsc () ;
     if (DISP_MAT) {
         printf("c3 = ") ;
         matrix_print_float_complexe(Layout, ROW, COL, c3) ;
     }
     printf ("Nombre de cycles: %Ld \n", end-start) ;
     calcul_flop ("cgemm ", ROW * COL * k, end-start) ;

     //Produit/somme matrice-matrice double
     if (DISP_MAT) {
         printf("a4 = ") ;
         matrix_print_double_complexe(Layout, ROW, COL, a4) ;
         printf("b4= ") ;
         matrix_print_double_complexe(Layout, ROW, COL, b4) ;
     }
     start = _rdtsc () ;
        mncblas_zgemm (Layout, transa, transb, ROW, COL, k, alpha_double, a4, b4, beta_double, c4) ;
     end = _rdtsc () ;
     if (DISP_MAT) {
         printf("c4 = ") ;
         matrix_print_double_complexe(Layout, ROW, COL, c4) ;
     }
     printf ("Nombre de cycles: %Ld \n", end-start) ;
     calcul_flop ("zgemm ", ROW * COL * k, end-start) ;
   }
   return 0 ;
}