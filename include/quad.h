#ifndef QUADG
#define QUADG

#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>

//1D Gaussian quadrature. Generated in Mathematica
void get_quad(int npts, double* Qi);
void test_2d_quad(int npts);

void quad_grid_munuphi(int tid, const int nptsx, const int nptsy, const int nptsz, int qosp, double a, double* Qx, double* Qy, double* Qz, int gs, int offset, double* gridm, double* grid, double* wt);
//batched version
void quad_grid_munuphi(int tid, int wb, int nb, const int nptsx, const int nptsy, const int nptsz, int qosp, double a, double* Qx, double* Qy, double* Qz, int gs, int offset, double* gridm, double* grid, double* wt);

void quad_grid_munu(const int nptsx, const int nptsy, double a, double* Qx, double* Qy, int gs, double* gridm, double* grid, double* wt);
void quad_grid_munu_one_pt(const int nptsx, const int nptsy, double a, const double mu1, const double dmu, const double nu1, const double dnu, const double phi,
                    double* Qx, double* Qy, double* grid, double* wt);

#endif
