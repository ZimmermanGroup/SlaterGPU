#ifndef BECKE_H
#define BECKE_H

#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <vector>

using std::vector;

void print_gradient(int	natoms,	double*	grad);

void atomic_domain_cell_wt(const int ta, const float alpha, const float beta, const int natoms, const int gc, const int gs, float* grid1, float* wt1, int* atno, float* coords);

void becke_weight_nc(int natoms, float* grid1, float* wt1, float* coords);

void becke_weight_3c(int gs, float* grid1, float* wt1, float* grid2, float* wt2, float* grid3, float* wt3, 
                     int Z1, int Z2, int Z3, float A2, float B2, float C2, float A3, float B3, float C3);
void becke_weight_2c(int gs, float* grid1, float* wt1, float* grid2, float* wt2,
                     int Z1, int Z2, float A2, float B2, float C2);
void becke_weight_2d(int gs, double* grid1, double* wt1, double* grid2, double* wt2,
                     double zeta1, double zeta2, double A2, double B2, double C2);

void get_becke_grid_full(int natoms, int* atno, float* coords, int nrad, int nang, float* ang_g, float* ang_w, const int gc, float* grid, float* wt);
void get_becke_grid_full(int natoms, int* atno, float* coords, int nrad, int nang, double* ang_g, double* ang_w, const int gc, float* grid, float* wt);
void get_becke_grid_full(int natoms, int* atno, float* coords, int nrad, int nang, double* ang_g, double* ang_w, const int gc, double* grid, double* wt);
void get_becke_grid_full(int natoms, int* atno, double* coords, int nrad, int nang, double* ang_g, double* ang_w, const int gc, double* grid, double* wt);
void compute_rho(int natoms, int* atno, float* coords, vector<vector<double> > &basis, double* Pao, int gsa, float* grid, double* rho, double* drho, int prl);
void compute_rhod(int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* Pao, int gsa, double* grid, double* rho, double* drho, int prl);
void compute_rhod(int natoms, int* atno, float* coords, vector<vector<double> > &basis, double* Pao, int gsa, float* grid, double* rho, double* drho, int prl);
void compute_rho(bool gbasis, int natoms, int* atno, float* coords, vector<vector<double> > &basis, double* Pao, int gsa, float* grid, double* rho, double* drho, int prl);

void compute_fxcd(int natoms, int* atno, double* coords, vector<vector<double> > &basis, bool gga, bool tau, bool need_wt, double* Pao, double* vxc, double* vxcs, int gsa, double* grid, double* wt, double* fxc, int prl);
void compute_fxc(int natoms, int* atno, float* coords, vector<vector<double> > &basis, bool gga, bool tau, bool need_wt, double* Pao, double* vxc, double* vxcs, int gsa, float* grid, float* wt, double* fxc, int prl);
void compute_fxc(int natoms, int* atno, float* coords, vector<vector<double> > &basis, bool gga, bool tau, bool need_wt, double* Pao, float* vxc, float* vxcs, int gs, float* grid, float* wt, double* fxc, int prl);
void compute_fxc(bool gbasis, int natoms, int* atno, float* coords, vector<vector<double> > &basis, bool need_wt, int gsa, float* grid, float* wt, double* vxc, double* fxc, int prl);
void compute_fxc(bool gbasis, int natoms, int* atno, float* coords, vector<vector<double> > &basis, bool need_wt, int gsa, float* grid, float* wt, float* vxc, double* fxc, int prl);

void compute_delta(int natoms, int* atno, float* coords, vector<vector<double> > &basis1, vector<vector<double> > &basis2, int No, double* jCA, bool gga, bool tau, float* rho, int gsa, float* grid, float* wt, double* diff, int prl);

void compute_dft_grad(int natoms, int* atno, float* coords, vector<vector<double> > &basis, double* Pao, bool is_gga, double* vxc, double* vxcs, int gs, float* grid, double* grad, int prl);

//void atomic_charges(int natoms, int gs, float* chg, double* rho, float** wta, int prl);
void atomic_charges(int natoms, int gs, float* chg, double* rho, float* zta, float** gridall, float** wta, int prl);
void batomic_charges(float alpha, int natoms, int* atno, float* coords, int nrad, int nang, float* ang_g, float* ang_w, float* chg, double* rho, int prl);

void becke_charges(int natoms, int gs, float* chg, double* rho, float* wt, int prl);
void becke_charges(int natoms, int gs, double* chg, double* rho, double* wt, int prl);
void density_in_basis2(int natoms, int* atno, float* coords, vector<vector<double> > &basis1, vector<vector<double> > &basis2, int No, double* jCA, int gsa, float* grid, float* wt, double* Paom, int prl);

float becke_a(int Z1, int Z2);
float bf3(float f1);
double bf3d(double f1);

//stretch atom 1's radius
float becke_a(float alpha, int Z1, int Z2);

//ratio based on exponents
double becke_a_zeta(double zeta1, double zeta2);

#endif
