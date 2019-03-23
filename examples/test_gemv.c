#include <stdio.h>
#include <x86intrin.h>

#include "vec_mat.h"

#include "flop.h"

#define NB_FOIS         1
#define ROW             1024
#define COL             1024   //qui est aussi la taille des vecteurs
#define DISP_MAT_VEC    0

mfloat a1 ;
mdouble a2 ;
mcfloat a3 ;
mcdouble a4 ;

vfloat x1, r1 ;
vdouble x2, r2 ;
vcfloat x3, r3 ;
vcdouble x4, r4 ;

int main (int argc, char **argv)
{
 unsigned long long start, end ;
 int i ;
 MNCBLAS_LAYOUT Layout = MNCblasRowMajor ;
 MNCBLAS_TRANSPOSE transa = MNCblasNoTrans ;
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
     x1 = vector_init_float (COL, 2.0) ;
     r1 = vector_init_float (COL, 0) ;
     a2 = matrix_init_double (ROW, COL, 1.0) ;
     x2 = vector_init_double (COL, 2.0) ;
     r2 = vector_init_double (COL, 0) ;
     a3 = matrix_init_float_complexe (ROW, COL, init_float_complexe1) ;
     x3 = vector_init_float_complexe (COL, init_float_complexe2) ;
     r3 = vector_init_float_complexe (COL, init_float_complexe0) ;
     a4 = matrix_init_double_complexe (ROW, COL, init_double_complexe1) ;
     x4 = vector_init_double_complexe (COL, init_double_complexe2) ;
     r4 = vector_init_double_complexe (COL, init_double_complexe0) ;

     //Produit/somme matrice-matrice float
     if (DISP_MAT_VEC) {
         printf("a1 = ") ;
         matrix_print_float(Layout, ROW, COL, a1) ;
         printf("x1 = ") ;
         vector_print_float(COL, x1) ;
     }
     start = _rdtsc () ;
         mncblas_sgemv (Layout, transa, ROW, COL, 1, a1, 0, x1, 1, 1, r1, 1) ;
     end = _rdtsc () ;
     if (DISP_MAT_VEC) {
         printf("r1 = ") ;
         vector_print_float(COL, r1) ;
     }
     printf ("Nombre de cycles: %Ld \n", end-start) ;
     calcul_flop ("sgemv ", 2 * ROW * COL, end-start) ;

     //Produit/somme matrice-matrice double
     if (DISP_MAT_VEC) {
         printf("a2 = ") ;
         matrix_print_double(Layout, ROW, COL, a2) ;
         printf("x2 = ") ;
         vector_print_double(COL, x2) ;
     }
     start = _rdtsc () ;
         mncblas_dgemv (Layout, transa, ROW, COL, 1, a2, 0, x2, 1, 1, r2, 1) ;
     end = _rdtsc () ;
     if (DISP_MAT_VEC) {
         printf("r2 = ") ;
         vector_print_double(COL, r2) ;
     }
     printf ("Nombre de cycles: %Ld \n", end-start) ;
     calcul_flop ("dgemv ", 2 * ROW * COL, end-start) ;

     ////Produit/somme matrice-matrice float complexe
     if (DISP_MAT_VEC) {
         printf("a3 = ") ;
         matrix_print_float_complexe(Layout, ROW, COL, a3) ;
         printf("x3 = ") ;
         vector_print_float_complexe(COL, x3) ;
     }
     start = _rdtsc () ;
        mncblas_cgemv (Layout, transa, ROW, COL, alpha_float, a3, 0, x3, 1, beta_float, r3, 1) ;
     end = _rdtsc () ;
     if (DISP_MAT_VEC) {
         printf("r3 = ") ;
         vector_print_float_complexe(COL, r3) ;
     }
     printf ("Nombre de cycles: %Ld \n", end-start) ;
     calcul_flop ("cgemv ", 2 * ROW * COL, end-start) ;

     //Produit/somme matrice-matrice double
     if (DISP_MAT_VEC) {
         printf("a4 = ") ;
         matrix_print_double_complexe(Layout, ROW, COL, a4) ;
         printf("x4= ") ;
         vector_print_double_complexe(COL, x4) ;
     }
     start = _rdtsc () ;
        mncblas_zgemv (Layout, transa, ROW, COL, alpha_double, a4, 0, x4, 1, beta_double, r4, 1) ;
     end = _rdtsc () ;
     if (DISP_MAT_VEC) {
         printf("r4 = ") ;
         vector_print_double_complexe(COL, r4) ;
     }
     printf ("Nombre de cycles: %Ld \n", end-start) ;
     calcul_flop ("zgemv ", 2 * ROW * COL, end-start) ;
   }
   return 0 ;
}