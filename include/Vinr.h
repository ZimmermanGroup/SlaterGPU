#ifndef VINRH
#define VINRH

#include "integrals.h"

void eval_ne(int gs, float* grid1, float** val, int s1, int s2, int natoms, int* atno, float* coords, float A0, float B0, float C0);
void eval_ke(int gs, float* grid, float* val, int n, int l, float zeta);
void eval_ke(int gs, double* grid, double* val, int n, int l, double zeta);
void eval_ke_erf(int gs, float* grid, float* val, int n, int l, float zeta, float ef1);
void eval_ke3(int gs, float* grid, float* val, int n, int l, float zeta);
void eval_dke(int gs, float* grid, float* val, int n, int l, float zeta);
void eval_ne_3(int gs, float* grid, float** val, int s1, int s2, int natoms, int* atno, float* coords, float A0, float B0, float C0);

void eval_inr_r12(int gs, float* grid, float* val, int n1, int l1, float zeta1, int index);
void eval_inr_r12(int gs, double* grid, double* val, int n1, int l1, double zeta1, int index);

void eval_inr_d(int gs, float* grid, float* val, int n1, int l1, float zeta1);
void eval_inr_d(int gs, double* grid, double* val, int n1, int l1, double zeta1);

void eval_inr_yukawa(int gs, double* grid, double* val, int n1, int l1, double zeta1, double gam1);
void eval_inr_yukawa_mt(int gs, double* grid, double* val, int n1, int l1, double zt1, double gam1);

void get_inr(int n1, int l1, double zeta1, int nrad, float* r, float* inr);

double norm_sh(int l, int m); //cartesian functions
double norm_sh_theta_phi(int l, int m); //theta/phi functions
double norm_sv(int n, int l, int m, double zeta);
double norm(int n, int l, int m, double zeta);
size_t fact(size_t N);

#endif
