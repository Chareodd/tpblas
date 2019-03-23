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

        //Copie float
        if (DISP_VEC) {
                printf("vec1 = ") ;
                vector_print_float(VECSIZE, vec1) ;
                printf("vec2 (initial) = ") ;
                vector_print_float(VECSIZE, vec2) ;
        }
        start = _rdtsc () ;
                mncblas_scopy(VECSIZE, vec1, 1, vec2, 1) ;
        end = _rdtsc () ;
        printf("Test de la copie (float) %d : vec2 = ", i) ;
        if (DISP_VEC) vector_print_float(VECSIZE, vec2) ;
        printf("Nombre de cycles: %Ld \n", end-start) ;
        calcul_flop ("Copie float ", 1, end-start) ;

        //Copie double
        if (DISP_VEC) {
                printf("\nvec3 = ") ;
                vector_print_double(VECSIZE, vec3) ;
                printf("vec4 (initial) = ") ;
                vector_print_double(VECSIZE, vec4) ;
        }
        start = _rdtsc () ;
                mncblas_dcopy(VECSIZE, vec3, 1, vec4, 1) ;
        end = _rdtsc () ;
        printf("Test de la copie (double) %d : vec4 = ", i) ;
        if (DISP_VEC) vector_print_double(VECSIZE, vec4) ;
        printf("Nombre de cycles: %Ld \n", end-start) ;
        calcul_flop ("Copie double ", 1, end-start) ;

        //Copie float complexe
        if (DISP_VEC) {
                printf("\nvec5 = ") ;
                vector_print_float_complexe(VECSIZE, vec5) ;
                printf("vec6 (initial) = ") ;
                vector_print_float_complexe(VECSIZE, vec6) ;
        }
        start = _rdtsc () ;
                mncblas_ccopy(VECSIZE, vec5, 1, vec6, 1) ;
        end = _rdtsc () ;
        printf("Test de la copie (float complexe) %d : vec6 = ", i) ;
        if (DISP_VEC) vector_print_float_complexe(VECSIZE, vec6) ;
        printf("Nombre de cycles: %Ld \n", end-start) ;
        calcul_flop ("Copie float complexe ", 1, end-start) ;

        //Copie double complexe
        if (DISP_VEC) {
                printf("\nvec7 = ") ;
                vector_print_double_complexe(VECSIZE, vec7) ;
                printf("vec8 (initial) = ") ;
                vector_print_double_complexe(VECSIZE, vec8) ;
        }
        start = _rdtsc () ;
                mncblas_zcopy(VECSIZE, vec7, 1, vec8, 1) ;
        end = _rdtsc () ;
        printf("Test de la copie (double complexe) %d : vec8 = ", i) ;
        if (DISP_VEC) vector_print_double_complexe(VECSIZE, vec8) ;
        printf("Nombre de cycles: %Ld \n", end-start) ;
        calcul_flop ("Copie double complexe ", 1, end-start) ;
   }
}