#ifndef OPTH
#define OPTH

#include <math.h>
#include "fp_def.h"

FP2 sd_step(FP1 scalar, FP1 maxstep, int natoms, FP2* grad, FP1* xyz);
FP2 cg_step(FP1 scalar, FP1 maxstep, FP2& grmsp, int natoms, FP2* step, FP2* grad, FP1* xyz);

#endif
