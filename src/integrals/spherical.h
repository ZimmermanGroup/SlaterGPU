#ifndef SPHARM
#define SPHARM

#include "fp_def.h"
#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "fp_def.h"

using namespace std;

void eval_sh_3r(int gs, FP1* grid, FP1* val, int n1, int l1, int m1);
void eval_sh_3r(int tid, int gs, FP1* grid, FP1* val, int n1, int l1, int m1);

void eval_sh_s(int gs, FP1* grid, FP1* val, int l1, int m1, FP1 A1, FP1 B1, FP1 C1);
void eval_sh_4r(int gs, FP1* grid, FP1* val, int l1, int m1, FP1 A1, FP1 B1, FP1 C1);

void eval_sh(int tid, int gs, FP1* grid, FP1* val, int n1, int l1, int m1, FP1 zeta);

void eval_sh_s(int tid, int gs, FP1* grid, FP1* val, int l1, int m1, FP1 A1, FP1 B1, FP1 C1);
void eval_sh_4r(int tid, int gs, FP1* grid, FP1* val, int l1, int m1, FP1 A1, FP1 B1, FP1 C1);

#endif
