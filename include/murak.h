#ifndef MURAK
#define MURAK
#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include "fp_def.h"

const int max_elem = 36;

const double alpha_k[max_elem] = 
 {
   5., 5.,
   7., 7., 5., 5., 5., 5., 5., 5., //through Ne
   7., 7., 5., 5., 5., 5., 5., 5., //through Ar
   7., 7., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., //through Zn
   5., 5., 5., 5., 5., 5.
 };

void get_murak_grid_f(int size, FP1* r, FP1* w, int Z, const int m);
void get_murak_grid_f(int size, FP1* r, FP1* w, FP1* er, int Z, FP1 zeta, const int m);
void get_murak_grid(int size, double* r, double* w, double* er, int Z, double zeta, const int m);

#endif
