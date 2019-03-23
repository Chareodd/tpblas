#include <stdio.h>
#include <x86intrin.h>

#include "vec_mat.h"

#include "flop.h"

#define NB_FOIS     1
#define VECSIZE    1024
#define DISP_VEC    0

vfloat vec1 ;
vdouble vec2 ;
vcfloat vec3 ;
vcdouble vec4 ;

int main (int argc, char **argv)
{
    unsigned long long start, end ;
    int i ;
    complexe_float_t init_float_complexe1 ;
    init_float_complexe1.real = 2 ;
    init_float_complexe1.imaginary = 2 ;
    complexe_double_t init_double_complexe1 ;
    init_double_complexe1.real = 2 ;
    init_double_complexe1.imaginary = 2 ;
    float res1, res3 ;
    double res2, res4 ;

    for (i = 0 ; i < NB_FOIS; i++)
   {
        vec1 = vector_init_float (VECSIZE, 2.0) ;
        vec2 = vector_init_double (VECSIZE, 2.0) ;
        vec3 = vector_init_float_complexe (VECSIZE, init_float_complexe1) ;
        vec4 = vector_init_double_complexe (VECSIZE, init_double_complexe1) ;

        //Somme magnitude float
        if (DISP_VEC) {
                printf("vec1 = ") ;
                vector_print_float(VECSIZE, vec1) ;
        }
        start = _rdtsc () ;
                res1 = mnblas_sasum(VECSIZE, vec1, 1) ;
        end = _rdtsc () ;
        printf("Test de la somme des magnitudes (float) %d : res = %f\n", i, res1) ;
        printf("Nombre de cycles: %Ld \n", end-start) ;
        calcul_flop ("Somme magnitude float ", VECSIZE, end-start) ;

        //Somme magnitude double
        if (DISP_VEC) {
                printf("\nvec2 = ") ;
                vector_print_double(VECSIZE, vec2) ;
        }
        start = _rdtsc () ;
                res2 = mnblas_dasum(VECSIZE, vec2, 1) ;
        end = _rdtsc () ;
        printf("Test de la somme des magnitudes (double) %d : res = %f\n", i, res2) ;
        printf("Nombre de cycles: %Ld \n", end-start) ;
        calcul_flop ("Somme magnitude double ", VECSIZE, end-start) ;

        //Somme magnitude float complexe
        if (DISP_VEC) {
                printf("\nvec3 = ") ;
                vector_print_float_complexe(VECSIZE, vec3) ;
        }
        start = _rdtsc () ;
                res3 = mnblas_scasum(VECSIZE, vec3, 1) ;
        end = _rdtsc () ;
        printf("Test de la somme des magnitudes (float complexe) %d : res = %f\n", i, res3) ;
        printf("Nombre de cycles: %Ld \n", end-start) ;
        calcul_flop ("Somme magnitude float complexe ", 2 * VECSIZE, end-start) ;

        //Somme magnitude double complexe
        if (DISP_VEC) {
                printf("\nvec4 = ") ;
                vector_print_double_complexe(VECSIZE, vec4) ;
        }
        start = _rdtsc () ;
                res4 = mnblas_dzasum(VECSIZE, vec4, 1) ;
        end = _rdtsc () ;
        printf("Test de la somme des magnitudes (double complexe) %d : res = %f\n", i, res4) ;
        printf("Nombre de cycles: %Ld \n", end-start) ;
        calcul_flop ("Somme magnitude double complexe ", 2 * VECSIZE, end-start) ;
   }
}