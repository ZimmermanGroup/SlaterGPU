#ifndef PVPH
#define PVPH

#include <stdio.h>
#include <math.h>
#include "fp_def.h"

void acc_assign(int size, FP1* vec, FP1 v1);

void eval_p(int	gs, FP1* grid1, FP1* val, int n1, int l1, int m1, FP1 zeta1);
void eval_dp_3r(int gs, FP1* grid, FP1* val, int n1, int l1, int m1);

#endif
