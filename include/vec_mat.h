#include "complexe.h"
#include "mnblas.h"

//VECTEUR

typedef float* vfloat ;
typedef double* vdouble ;
typedef complexe_float_t* vcfloat ;
typedef complexe_double_t* vcdouble ;

//Initialise toutes les composantes d'un vecteur avec une unique valeur
vfloat vector_init_float (int size, float x) ;

vdouble vector_init_double (int size, double x) ;

vcfloat vector_init_float_complexe (int size, complexe_float_t x) ;

vcdouble vector_init_double_complexe (int size, complexe_double_t x) ;

//Affiche toutes les composantes d'un vecteur
void vector_print_float (int size, vfloat V) ;

void vector_print_double (int size, vdouble V) ;

void vector_print_float_complexe (int size, vcfloat V) ;

void vector_print_double_complexe (int size, vcdouble V) ;

//MATRICE

typedef float* mfloat ;
typedef double* mdouble ;
typedef complexe_float_t* mcfloat ;
typedef complexe_double_t* mcdouble ;

//Initialise toutes les composantes d'un vecteur avec une unique valeur
mfloat matrix_init_float (int row, int col, float x) ;

mdouble matrix_init_double (int row, int col, double x) ;

mcfloat matrix_init_float_complexe (int row, int col, complexe_float_t x) ;

mcdouble matrix_init_double_complexe (int row, int col, complexe_double_t x) ;

//Affiche toutes les composantes d'un vecteur
void matrix_print_float (MNCBLAS_LAYOUT Layout, int row, int col, mfloat M) ;

void matrix_print_double (MNCBLAS_LAYOUT Layout, int row, int col, mdouble M) ;

void matrix_print_float_complexe (MNCBLAS_LAYOUT Layout, int row, int col, mcfloat M) ;

void matrix_print_double_complexe (MNCBLAS_LAYOUT Layout, int row, int col, mcdouble M) ;