#ifndef PRINT_H
#define PRINT_H

#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>

#include <omp.h>

using std::vector;
using std::string;

#define A2B 1.8897261

void print_mos_col(int M, int N, vector<vector<double> > basis, double* jCA);
void print_square_diff(int N, double* S1, double* S2);
void print_square(int N, double* S);
void print_square(int N, float* S);
void print_square(int M, int N, float* S);
void print_square(int M, int N, double* S);
void print_square_col(int M, int N, double* S);
void print_square_col_sm(int M, int N, double* S);
void print_square_fine(int N, float* S);
void print_square_fine(int M, int N, double* S);
void print_square_fine(int N, double* S);
void print_square_ss(int N, double* S);
void print_square_sm(int N, float* S);
void print_square_sm(int N, double* S);
void print_square_ss_sm(int N, double* S);
void print_square_nxn(int No, int N, float* S);
void print_square_nxn(int No, int N, double* S);
void print_coords(int natoms, double* coords);
void print_gradient(int natoms, double* grad);

void print_rectangle(int N1, int N2, double* S);
void print_rectangle_e(int N1, int N2, double* S);
void print_rectangle_sm(int N1, int N2, double* S);
void print_rdm(int M, double* rdm);
void print_rdm(int M, float* rdm);
void print_vec(int gsa, float* grid, float* vxc);
void print_vec(int gsa, float* grid, double* vxc);
void print_vec(int gsa, double* grid, double* vxc);
void print_vec_fine(int gsa, double* grid, double* vxc);
void print_vec(int gsa, float* grid, float* A, float* B);
void print_vec(int gsa, float* grid, double* A, double* B);

void print_dft_vals(int natoms, int gs, double* grid, double* rho, double* drho, double* Td, double* ei, double* vc, int zpos, bool scirep);
void print_vc(int natoms, int gs, float* grid, float* rho, float* vc, int mode);
void print_vc(int natoms, int gs, float* grid, double* rho, double* vc, int mode);
void print_vc_shift(int natoms, int gs, double* grid, float* rho, double* vc, int mode);
void print_vc_shift(int natoms, int gs, float* grid, double* rho, double* vc, int mode);
void print_vc_shift(int natoms, int gs, double* grid, double* rho, double* vc, int mode);
void print_wt(int natoms, int gs, float* grid, float* rho, float* wt, int mode);
void print_vc(int natoms, int gs, float* grid, float* rho, double* vc, int mode);
void print_vc_shift(int natoms, int gs, float* grid, float* rho, float* vc, int mode);
void print_vc_shift(int natoms, int gs, float* grid, float* rho, double* vc, int mode);
void print_vc_shift(int natoms, int gs, float* grid, double* rho, float* vc, int mode);
void print_vxc(int nrad, int nang, int natoms, float* grid, double* vxc, string name);
void print_vxc(int nrad, int nang, int natoms, double* grid, double* vxc, string name);
void print_vxc(int nrad, int nang, int natoms, float* grid, float* vxcf, string name);
void print_vxc(int nrad, int nang, int natoms, float* grid, double* vxc);
void print_vxc(int nrad, int nang, int natoms, double* grid, double* vxc);
void print_vxc(int nrad, int nang, int natoms, float* grid, float* vxcf);
void print_sphere(int nrad, int nang, int natoms, float* grid, double* vxc);
void print_axes(int nrad, int nang, int natoms, float* grid, float* r1, float* rho, float* T, float* gf, float* vxch);
void print_axes(int nrad, int nang, int natoms, float* grid, float* r1, float* rho, float* T, float* gf, double* vxc);

#endif

