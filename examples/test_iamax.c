#include <stdio.h>
#include <x86intrin.h>

#include "vec_mat.h"

#include "flop.h"

#define NB_FOIS     1
#define VECSIZE    10
#define DISP_VEC    0

vfloat vec1;
vdouble vec3;
vcfloat vec5;
vcdouble vec7;

int main (int argc, char **argv)
{
 unsigned long long start, end ;
 int res_index;
 int i ;
 complexe_float_t init_float_complexe1 ;
 init_float_complexe1.real = 1 ;
 init_float_complexe1.imaginary = 1 ;
 complexe_double_t init_double_complexe1 ;
 init_double_complexe1.real = 1 ;
 init_double_complexe1.imaginary = 1 ;
 complexe_double_t init_double_complexe2 ;

 for (i = 0 ; i < NB_FOIS; i++)
   {
     vec1 = vector_init_float (VECSIZE, 1.0) ;
     vec3 = vector_init_double (VECSIZE, 1.0) ;
     vec5 = vector_init_float_complexe (VECSIZE, init_float_complexe1) ;
     vec7 = vector_init_double_complexe (VECSIZE, init_double_complexe1) ;

     //Index maximal float
     if (DISP_VEC) {
         printf("vec1 = ") ;
         vector_print_float(VECSIZE, vec1) ;
     }
     start = _rdtsc () ;
        res_index = mncblas_isamax(VECSIZE, vec1, 1) ;
     end = _rdtsc () ;
     printf ("mncblas_isamax %d : res = %3.2f nombre de cycles: %Ld \n", i, res_index, end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;

     //Index maximal double
     if (DISP_VEC) {
         printf("vec3 = ") ;
         vector_print_double(VECSIZE, vec3) ;
     }
     start = _rdtsc () ;
        res_index = mncblas_idamax(VECSIZE, vec3, 1) ;
     end = _rdtsc () ;
     printf ("mncblas_idamax %d : res = %3.2f nombre de cycles: %Ld \n", i, res_index, end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;

     //Index maximal simple complexe
     if (DISP_VEC) {
         printf("vec5 = ") ;
         vector_print_float_complexe(VECSIZE, vec5) ;
     }
     start = _rdtsc () ;
        res_index = mncblas_icamax(VECSIZE, vec5, 1) ;
     end = _rdtsc () ;
     printf ("mncblas_icamax %d : res = %3.2f nombre de cycles: %Ld \n", i, res_index, end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;

     //Index maximal double complexe
     if (DISP_VEC) {
         printf("vec7 = ") ;
         vector_print_float_complexe(VECSIZE, vec7) ;
     }
     start = _rdtsc () ;
        res_index = mncblas_izamax(VECSIZE, vec7, 1) ;
     end = _rdtsc () ;
     printf ("mncblas_izamax %d : res = %3.2f nombre de cycles: %Ld \n", i, res_index, end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;
   }
}
