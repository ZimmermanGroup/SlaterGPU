#ifndef SPHARM
#define SPHARM

#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <vector>

using namespace std;

void eval_sh_3r(int gs, float* grid, float* val, int n1, int l1, int m1);
void eval_sh_3r(int tid, int gs, float* grid, float* val, int n1, int l1, int m1);

void eval_sh_3rd(int gs, double* grid, double* val, int n1, int l1, int m1);
void eval_sh_3rd(int tid, int gs, double* grid, double* val, int n1, int l1, int m1);


void eval_sh_s(int gs, float* grid, float* val, int l1, int m1, float A1, float B1, float C1);
void eval_sh_4r(int gs, float* grid, float* val, int l1, int m1, float A1, float B1, float C1);

void eval_sh(int tid, int gs, float* grid, float* val, int n1, int l1, int m1, float zeta);
void eval_shd(int tid, int gs, double* grid, double* val, int n1, int l1, int m1, double zeta1);

void eval_sh_s(int tid, int gs, float* grid, float* val, int l1, int m1, float A1, float B1, float C1);
void eval_sh_4r(int tid, int gs, float* grid, float* val, int l1, int m1, float A1, float B1, float C1);

#endif
