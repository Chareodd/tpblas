#include <stdio.h>
#include <x86intrin.h>

#include "vec_mat.h"

#include "flop.h"

#define NB_FOIS     1
#define VECSIZE    1024
#define DISP_VEC    0

vfloat vec1, vec2 ;
vdouble vec3, vec4 ;
vcfloat vec5, vec6 ;
vcdouble vec7, vec8 ;

int main (int argc, char **argv)
{
 unsigned long long start, end ;
 float res_float ;
 double res_double ;
 int i ;
 complexe_float_t init_float_complexe1 ;
 init_float_complexe1.real = 1 ;
 init_float_complexe1.imaginary = 1 ;
 complexe_float_t init_float_complexe2 ;
 init_float_complexe2.real = 2 ;
 init_float_complexe2.imaginary = 2 ;
 complexe_double_t init_double_complexe1 ;
 init_double_complexe1.real = 1 ;
 init_double_complexe1.imaginary = 1 ;
 complexe_double_t init_double_complexe2 ;
 init_double_complexe2.real = 2 ;
 init_double_complexe2.imaginary = 2 ;
 complexe_float_t *dotu = (complexe_float_t*) malloc(sizeof(complexe_float_t)) ;
 complexe_double_t *dotc = (complexe_double_t*) malloc(sizeof(complexe_double_t)) ;

 for (i = 0 ; i < NB_FOIS; i++)
   {
     vec1 = vector_init_float (VECSIZE, 1.0) ;
     vec2 = vector_init_float (VECSIZE, 2.0) ;
     vec3 = vector_init_double (VECSIZE, 1.0) ;
     vec4 = vector_init_double (VECSIZE, 2.0) ;
     vec5 = vector_init_float_complexe (VECSIZE, init_float_complexe1) ;
     vec6 = vector_init_float_complexe (VECSIZE, init_float_complexe2) ;
     vec7 = vector_init_double_complexe (VECSIZE, init_double_complexe1) ;
     vec8 = vector_init_double_complexe (VECSIZE, init_double_complexe2) ;

     //Produit scalaire float
     if (DISP_VEC) {
         printf("vec1 = ") ;
         vector_print_float(VECSIZE, vec1) ;
         printf("vec2 = ") ;
         vector_print_float(VECSIZE, vec2) ;
     }
     start = _rdtsc () ;
        res_float = mncblas_sdot (VECSIZE, vec1, 1, vec2, 1) ;
     end = _rdtsc () ;
     printf ("mncblas_sdot %d : res = %3.2f nombre de cycles: %Ld \n", i, res_float, end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;

     //Produit scalaire double
     if (DISP_VEC) {
         printf("vec3 = ") ;
         vector_print_double(VECSIZE, vec3) ;
         printf("vec4 = ") ;
         vector_print_double(VECSIZE, vec4) ;
     }
     start = _rdtsc () ;
        res_double = mncblas_ddot (VECSIZE, vec3, 1, vec4, 1) ;
     end = _rdtsc () ;
     printf ("mncblas_ddot %d : res = %3.2f nombre de cycles: %Ld \n", i, res_double, end-start) ;
     calcul_flop ("ddot ", 2 * VECSIZE, end-start) ;

     //Produit scalaire float complexe
     if (DISP_VEC) {
         printf("vec5 = ") ;
         vector_print_float_complexe(VECSIZE, vec5) ;
         printf("vec6 = ") ;
         vector_print_float_complexe(VECSIZE, vec6) ;
     }
     start = _rdtsc () ;
        mncblas_cdotu_sub (VECSIZE, vec5, 1, vec6, 1, dotu) ;
     end = _rdtsc () ;
     printf ("mncblas_cdotu_sub %d : res = %3.2f + i %3.2f  nombre de cycles: %Ld \n", i, dotu->real, dotu->imaginary, end-start) ;
     calcul_flop ("cdotu ", 2 * VECSIZE, end-start) ;

     //Produit scalaire float complexe conjugué
     if (DISP_VEC) {
         printf("vec5 = ") ;
         vector_print_float_complexe(VECSIZE, vec5) ;
         printf("vec6= ") ;
         vector_print_float_complexe(VECSIZE, vec6) ;
     }
     start = _rdtsc () ;
        mncblas_cdotc_sub (VECSIZE, vec5, 1, vec6, 1, dotc) ;
     end = _rdtsc () ;
     printf ("mncblas_cdotc_sub %d : res = %3.2f + i %3.2f  nombre de cycles: %Ld \n", i, dotc->real, dotc->imaginary, end-start) ;
     calcul_flop ("cdotc ", 2 * VECSIZE, end-start) ;

     //Produit scalaire double complexe
     if (DISP_VEC) {
         printf("vec7 = ") ;
         vector_print_double_complexe(VECSIZE, vec7) ;
         printf("vec8 = ") ;
         vector_print_double_complexe(VECSIZE, vec8) ;
     }
     start = _rdtsc () ;
        mncblas_zdotu_sub (VECSIZE, vec7, 1, vec8, 1, dotu) ;
     end = _rdtsc () ;
     printf ("mncblas_zdotu_sub %d : res = %3.2f + i %3.2f  nombre de cycles: %Ld \n", i, dotu->real, dotu->imaginary, end-start) ;
     calcul_flop ("zdotu ", 2 * VECSIZE, end-start) ;

     //Produit scalaire double complexe conjugué
     if (DISP_VEC) {
         printf("vec7 = ") ;
         vector_print_double_complexe(VECSIZE, vec7) ;
         printf("vec8 = ") ;
         vector_print_double_complexe(VECSIZE, vec8) ;
     }
     start = _rdtsc () ;
        mncblas_zdotc_sub (VECSIZE, vec7, 1, vec8, 1, dotc) ;
     end = _rdtsc () ;
     printf ("mncblas_zdotc_sub %d : res = %3.2f + i %3.2f  nombre de cycles: %Ld \n", i, dotc->real, dotc->imaginary, end-start) ;
     calcul_flop ("zdotc ", 2 * VECSIZE, end-start) ;
   }
}
