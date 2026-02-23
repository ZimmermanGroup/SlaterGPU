#ifndef _GRID_UTIL_H_
#define _GRID_UTIL_H_

#include "read.h"
#include "write.h"

#include "becke.h"
#include <accel.h>

void compare_pao_12(int ngpu, bool gbasis, int natoms, int* atno, double* coords, int nrad, int nang, double* ang_g, double* ang_w, vector<vector<double> >& basis, cusolverDnHandle_t cu_hdl);
void save_grid_rho(bool gbasis, int natoms, int* atno, double* coords, int nrad, int nang, double* ang_g, double* ang_w, vector<vector<double> >& basis);
void save_grid_ao_basis(bool gbasis, int natoms, int* atno, double* coords, int nrad, int nang, double* ang_g, double* ang_w, vector<vector<double> >& basis);


//void get_angular_grid(int size_ang, double* ang_g, double* ang_w);
//int electron_count(int charge, int natoms, int* atno, int& Nc, int& Na, int& Nb);

#endif // _GRID_UTIL_H_
