#ifndef VINRH
#define VINRH

#include "fp_def.h"
#include "integrals.h"

void eval_ne(int gs, FP1* grid1, FP1** val, int s1, int s2, int natoms, int* atno, FP1* coords, FP1 A0, FP1 B0, FP1 C0);
void eval_ke(int gs, FP1* grid, FP1* val, int n, int l, FP1 zeta);
void eval_ke3(int gs, FP1* grid, FP1* val, int n, int l, FP1 zeta);
void eval_dke(int gs, FP1* grid, FP1* val, int n, int l, FP1 zeta);
void eval_ne_3(int gs, FP1* grid, FP1** val, int s1, int s2, int natoms, int* atno, FP1* coords, FP1 A0, FP1 B0, FP1 C0);

void eval_inr_r12(int gs, FP1* grid, FP1* val, int n1, int l1, FP1 zeta1);
void eval_inr_d(int gs, FP1* grid, FP1* val, int n1, int l1, FP1 zeta1);

void get_inr(int n1, int l1, FP2 zeta1, int nrad, FP1* r, FP1* inr);
void eval_exp_r1(int tid, int gs, FP1* grid, FP1* val, FP1 zeta);
void eval_exp_r2(int tid, int gs, FP1* grid, FP1* val, FP1 zeta);
void eval_exp_r1(int gs, FP1* grid, FP1* val, FP1 zeta);
void eval_exp_r2(int gs, FP1* grid, FP1* val, FP1 zeta);
FP2 norm_sh(int l, int m);
FP2 norm_sv(int n, int l, int m, FP2 zeta);
FP2 norm(int n, int l, int m, FP2 zeta);
int fact(int N);

#endif
