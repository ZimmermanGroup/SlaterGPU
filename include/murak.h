#ifndef MURAK
#define MURAK
#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>

const int max_elem = 36;

const double alpha_k[max_elem] = 
 {
   5., 5.,
   7., 7., 5., 5., 5., 5., 5., 5., //through Ne
   7., 7., 5., 5., 5., 5., 5., 5., //through Ar
   7., 7., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., //through Zn
   5., 5., 5., 5., 5., 5.
 };


void get_eumac_grid(int size, double* r, double* w, const double rmax, const int m);

void get_murak_grid_f(int size, float* r, float* w, int Z, const int m);
void get_murak_grid(int size, double* r, double* w, int Z, const int m);
void get_murak_grid_f(int size, float* r, float* w, float* er, int Z, float zeta, const int m);
void get_murak_grid(int size, double* r, double* w, double* er, int Z, double zeta, const int m);

//grid specific to a particular exponential
void get_murak_grid_zeta(int size, double* r, double* w, const double zeta, const int m);

#endif
