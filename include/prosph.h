#ifndef PROSPH
#define PROSPH

#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <vector>

using namespace std;

#include "integrals.h"
#include "read.h"
#include "quad.h"

void print_square_fine(int N, double* S);

//
void test_prosph();

void generate_ps_quad_grid(double Z1, int natoms, double* coordn, int quad_order, int quad_r_order, int nmu, int nnu, int nphi, double* grid, double* wt);
void generate_ps_quad_grid(int wb, int nb, double Z1, int natoms, double* coordn, int quad_order, int quad_r_order, int nmu, int nnu, int nphi, double* grid, double* wt);
void generate_ps_quad_grid(double cfn, double Z1, int natoms, double* coordn, int quad_order, int quad_r_order, int nmu, int nnu, int nphi, double* grid, double* wt);
void generate_ps_quad_grid(double cfn, int wb, int nb, double Z1, int natoms, double* coordn, int quad_order, int quad_r_order, int nmu, int nnu, int nphi, double* grid, double* wt);
void generate_ps_quad_grid_3c_refine(double ztm1, double ztm2, int nsplit, int natoms, double* coordn, int quad_order, int quad_r_order, int nmu, int nnu, int nphi, double* grid, double* wt);
void generate_ps_quad_grid_3c_refine(int wb, int nb, double ztm1, double ztm2, int nsplit, int natoms, double* coordn, int quad_order, int quad_r_order, int nmu, int nnu, int nphi, double* grid, double* wt);
void generate_ps_quad_grid_3c_refine(double cfn, int wb, int nb, double ztm1, double ztm2, int nsplit, int natoms, double* coordn, int quad_order, int quad_r_order, int nmu, int nnu, int nphi, double* grid, double* wt);


void initialize_ps_coords_2c(double a, double cf, int nmu, int nnu, int nphi, double phi0, double* grid, double* gridm, double* wt, int prl);
void initialize_ps_coords_3c(double cf, int nmu, int nnu, int nphi, double phi0, double* coordn, double* grid, double* gridm, double* wt, double* rot, int prl);
void initialize_ps_coords_batch(int wb, int nbatch, double a, double cf, int nmu, int nnu, int nphi, double phi0, double* grid, double* gridm, double* wt, int prl);
void get_2c_position(double* coordn, double* rot);
void get_3c_position(double* coordn, double* rot);
void get_4c_position(double* coordn, double* rot);
void gen_total_rot_n2(int natoms, double* coords, double* trot);
void gen_total_rot_n3(int natoms, double* coords, double* trot);
void rotate_3x3(double* rot_mat, double* A, double* B);
void reorient_grid(const double z0, int gs, double* grid, double* grid2, double* rot);

void integrate_STEnAC_2c(int natoms, int* atno, float* coordsf, vector<vector<double> > basis, vector<vector<double> > basis_aux, double epsilon, int nmu, int nnu, int nphi, bool do_coulomb, double* S, double* T, double* En, double* A, double* C, int prl);
void integrate_STEnAC_2c(int natoms, int* atno, double* coords, vector<vector<double> > basis, vector<vector<double> > basis_aux, double epsilon, int nmu, int nnu, int nphi, bool do_coulomb, double* S, double* T, double* En, double* A, double* C, int prl);
bool integrate_ol_4c(int natoms, int* atno, double* coords, vector<vector<double> > basis, double epsilon, int nmu, int nnu, int nphi, double* g, int prl);

double get_idp(int l1, int m1, int l2, int m2);
double get_idp_one_atom(bool do_coulomb, int l1, int m1, int l2, int m2);
//double get_idp13(int l1, int m1, int l2, int m2, int l3, int m3);
double get_idp_m4(int m1, int m2, int m3);
double get_idp_m5(int m1, int m2, int m3);
double get_idp_4b_m3(int m1, int m2, int m3, int m4);


//integrals expanded to second order in exp(r1)exp(r2)
double second_order_fordV(double Z1, double Z2, double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double a, double a0, double b1, double b2, double c1, double c2, double bc);
double second_order_fx2ordV(double Z1, double Z2, double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double a, double a0, double b1, double b2, double c1, double c2, double bc);

double second_order_fzzdV(double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double a, double a0, double b1, double b2, double c1, double c2, double bc);
double second_order_fzdV(double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double a, double zc, double a0, double b1, double b2, double c1, double c2, double bc);
double second_order_fx2dV(double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double a, double a0, double b1, double b2, double c1, double c2, double bc);
double second_order_fdV(double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double a, double a0, double b1, double b2, double c1, double c2, double bc);

double first_order_fdV(double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double phi0, double phi1, double phi2, double a, double b, double c, double d);
double first_order_fdV(double mu0, double mu1, double mu2, double nu1, double nu2, double a, double b);

#pragma acc routine seq
double ps_dV(double mu1, double mu2, double nu1, double nu2);

void get_ab_r1r2(int gs, const double z0, const double Z1, const double Z2, double* val1, double* val2, double* val1p, double* val2p, double* grid, double* gridm, double* wt, double* val);
void get_ab_2d(int gs, const double z0, double* val1, double* val2, double* val1p, double* val2p, double* grid, double* gridm, double* wt, double* val);
void get_ab_3d(int gs, const double z0, double* val1, double* val2, double* val1p, double* val2p, double* grid, double* gridm, double* wt, double* val);
void get_ab_mnp(int gs, const double a, bool do_r12, double Z1, double Z2, int n1, int l1, int m1, int n2, int l2, int m2, double norm1, double norm2, double zt1, double zt2, double* grid, double* gridm, double* wt, double* val);


#endif

int sign(double val);
