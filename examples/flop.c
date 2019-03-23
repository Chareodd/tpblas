#include <stdio.h>
#include <x86intrin.h>

static const float duree_cycle = (float) 1 / (float) 2.6 ;
// duree du cycle en nano seconde 10^-9

void calcul_flop (char *message, int nb_operations_flottantes, unsigned long long int cycles)
{
  printf ("%s %d operations %5.3f GFLOP/s\n", message, nb_operations_flottantes, ((float) nb_operations_flottantes) / (((float) cycles) * duree_cycle)) ;
  return ;
}
