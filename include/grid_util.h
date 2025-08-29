#ifndef _GRID_UTIL_H_
#define _GRID_UTIL_H_

#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <vector>

using std::vector;

void compare_pao_12(bool gbasis, int natoms, int* atno, double* coords, int nrad, int nang, double* ang_g, double* ang_w, vector<vector<double> >& basis);

void save_grid_ao_basis(bool gbasis, int natoms, int* atno, double* coords, int nrad, int nang, double* ang_g, double* ang_w, vector<vector<double> >& basis);

void save_grid_rho(bool gbasis, int natoms, int* atno, double* coords, int nrad, int nang, double* ang_g, double* ang_w, vector<vector<double> >& basis);

#endif // _GRID_UTIL_H_
