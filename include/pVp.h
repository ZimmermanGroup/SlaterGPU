#ifndef PVPH
#define PVPH

#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void acc_assign(int size, float* vec, float v1);
void acc_assign(int size, double* vec, double v1);

void eval_p(int gs, float* grid1, float* val, int n1, int l1, int m1, float zeta1);
void eval_dp_3r(int gs, float* grid, float* val, int n1, int l1, int m1);

//double precision versions
void eval_pd(int tid, int gs, double* grid1, double* val, int n1, int l1, int m1, double zeta1);
void eval_dp_3rd(int tid, int gs, double* grid, double* val, int n1, int l1, int m1);

//Hessian elements (hess.cpp)
void eval_h(int gs, double* grid, double* val, int n1, int l1, int m1, double zeta1);
void eval_h(int gs, float* grid, float* val, int n1, int l1, int m1, float zeta1);

void eval_ss_pd(int gs, double* grid, double* val, double* tmp, int n1, int l1, int m1, const double zeta1, double Rc);

#endif
