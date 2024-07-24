#ifndef OPTH
#define OPTH

#include <math.h>
//#include "fp_def.h"

double sd_step(float scalar, float maxstep, int natoms, double* grad, float* xyz);
double cg_step(float scalar, float maxstep, double& grmsp, int natoms, double* step, double* grad, float* xyz);

#endif
