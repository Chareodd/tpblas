#include <stdio.h>
#include <x86intrin.h>

#include "vec_mat.h"

#include "flop.h"

#define NB_FOIS    1

int main (int argc, char **argv)
{
    unsigned long long start, end ;

    complexe_float_t res_float_complexe ;
    complexe_double_t res_double_complexe ;

    complexe_float_t init_float_complexe1 ;
    complexe_float_t init_float_complexe2 ;
    complexe_double_t init_double_complexe1 ;
    complexe_double_t init_double_complexe2 ;

    for (int i = 0; i < NB_FOIS; i++)
    {
        init_float_complexe1.real = 1 ;
        init_float_complexe1.imaginary = 1 ;
        
        init_float_complexe2.real = 2 ;
        init_float_complexe2.imaginary = 2 ;
        
        init_double_complexe1.real = 1 ;
        init_double_complexe1.imaginary = 1 ;
        
        init_double_complexe2.real = 2 ;
        init_double_complexe2.imaginary = 2 ;

        //Conjugué float
        start = _rdtsc () ;
            res_float_complexe = conjg_float(init_float_complexe2) ;
        end = _rdtsc () ;
        printf("Test du conjugué (float) %d : res = %3.2f + i %3.2f \nNombre de cycles: %Ld \n\n", i, res_float_complexe.real, res_float_complexe.imaginary, end-start) ;
        calcul_flop ("Conjugué float ", 1, end-start) ;

        //Conjugué double
        start = _rdtsc () ;
            res_double_complexe = conjg_double(init_double_complexe2) ;
        end = _rdtsc () ;
        printf("Test du conjugué (double) %d : res = %3.2f + i %3.2f \nNombre de cycles: %Ld\n\n", i, res_float_complexe.real, res_float_complexe.imaginary, end-start) ;
        calcul_flop ("Conjugué double ", 1, end-start) ;

        //Addition float
        start = _rdtsc () ;
            res_float_complexe = add_complexe_float(init_float_complexe1, init_float_complexe2) ;
        end = _rdtsc () ;
        printf("Test de la somme (float) %d : res = %3.2f + i %3.2f \nNombre de cycles: %Ld\n\n", i, res_float_complexe.real, res_float_complexe.imaginary, end-start) ;
        calcul_flop ("Addition float ", 2, end-start) ;

        //Addition double
        start = _rdtsc () ;
            res_double_complexe = add_complexe_double(init_double_complexe1, init_double_complexe2) ;
        end = _rdtsc () ;
        printf("Test de la somme (double) %d : res = %3.2f + i %3.2f \nNombre de cycles: %Ld\n\n", i, res_double_complexe.real, res_double_complexe.imaginary, end-start) ;
        calcul_flop ("Addition double ", 2, end-start) ;

        //Multiplication float
        start = _rdtsc () ;
            res_float_complexe = mult_complexe_float(init_float_complexe1, init_float_complexe2) ;
        end = _rdtsc () ;
        printf("Test de la multiplication (float) %d : res = %3.2f + i %3.2f \nNombre de cycles: %Ld\n\n", i, res_float_complexe.real, res_float_complexe.imaginary, end-start) ;
        calcul_flop ("Multiplication float ", 6, end-start) ;

        //Multiplication double
        start = _rdtsc () ;
            res_double_complexe = mult_complexe_double(init_double_complexe1, init_double_complexe2) ;
        end = _rdtsc () ;
        printf("Test de la multiplication (double) %d : res = %3.2f + i %3.2f \nNombre de cycles: %Ld \n\n", i, res_double_complexe.real, res_double_complexe.imaginary, end-start) ;
        calcul_flop ("Multiplication double ", 6, end-start) ;
    }
}       