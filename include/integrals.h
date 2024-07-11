#ifndef INTEGRALS
#define INTEGRALS

#define DOUBLE 0
#define DEBUG 0
//KE is slightly sensitive to this cutoff

//was 1.e-10
#define WT_THRESH 1.e-10
#define WT_THRESH_D 1.e-10
#define CL_THRESH 1.e-7
#define REDUCEV 0


#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>

#include <chrono>

#include "accel.h"
#include <omp.h>

using namespace std;

#include "spherical.h"
#include "Vinr.h"
#include "pVp.h"
#include "reduce.h"
#include "murak.h"
#include "becke.h"
#include "sortUtil.h"

extern void sortByKeys( fp_t * A, xyz<fp_t> * Av, size_t N);

#define PI 3.141592653589793238


void print_square(int N, double* S);

double simple_test(int size1, int size2);

float compute_2c_trial(float zeta1, float zeta2, int nrad, int nang, double* ang_g, double* ang_w);


void compute_d_ST(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* GFao, double* Pao, double* xyz_grad, int prl);
void compute_d_ST(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* GFao, double* Pao, float* xyz_grad, int prl);

void compute_d_En(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* Pao, double* xyz_grad, int prl);
void compute_d_En(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* Pao, float* xyz_grad, int prl);

void compute_d_2c(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* dpq, double* xyz_grad, int prl);
void compute_d_2c(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* dpq, float* xyz_grad, int prl);

void compute_d_3c(int natoms, int* atno, float* coords, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, int nrad, int nang, double* ang_g0, double* ang_w0, double* dC, double* xyz_grad, int prl);
void compute_d_3c(int natoms, int* atno, float* coords, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, int nrad, int nang, double* ang_g0, double* ang_w0, double* dC, float* xyz_grad, int prl);

void compute_d_3c_para(int ngpu, int natoms, int* atno, float* coords, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, int nrad, int nang, double* ang_g0, double* ang_w0, double* dC, double* xyz_grad, int prl);
void compute_d_3c_para(int ngpu, int natoms, int* atno, float* coords, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, int nrad, int nang, double* ang_g0, double* ang_w0, double* dC, float* xyz_grad, int prl);


void compute_Enp(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* En, double* pVp, int prl);
void compute_Enp(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, float* En, float* pVp, int prl);

void compute_Enp_para(int ngpu, int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* En, double* pVp, int prl);
void compute_Enp_para(int ngpu, int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, float* En, float* pVp, int prl);

//electric fields in x,y,z (centered at origin)
void compute_Exyz(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g, double* ang_w, double* E, int prl);

void compute_ST(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* S, double* T, int prl);
void compute_ST(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, float* S, float* T, int prl);

void compute_cusp(int natoms, int* atno, float* coords, vector<vector<double> > &basis, double* pb, int prl);

void compute_all_2c(int natoms, int* atno, float* coordsf, vector<vector<double> > &basis_aux, int nrad, int nang, double* ang_g0, double* ang_w0, float* A, int prl);
//void compute_all_2c(int natoms, int* atno, float* coordsf, vector<vector<double> > &basis_aux, int nrad, int nang, double* ang_g0, double* ang_w0, double* A, int prl);

void compute_all_2c_v2(bool do_overlap, int natoms, int* atno, float* coordsf, vector<vector<double> > &basis_aux, int nrad, int nang, double* ang_g0, double* ang_w0, double* A, int prl);
void compute_all_2c_v2(bool do_overlap, int natoms, int* atno, float* coordsf, vector<vector<double> > &basis_aux, int nrad, int nang, double* ang_g0, double* ang_w0, float* A, int prl);

//fully double precision
void compute_Sd(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* S, int prl);
void compute_all_2c_v2d(bool do_overlap, int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* An, int prl);

void compute_all_3c(int natoms, int* atno, float* coords, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, int nrad, int nang, double* ang_g0, double* ang_w0, float* C, int prl);
//void compute_all_3c(int natoms, int* atno, float* coords, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, int nrad, int nang, double* ang_g0, double* ang_w0, double* C, int prl);

void compute_all_3c_v2(bool do_overlap, int natoms, int* atno, float* coords, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, int nrad, int nang, double* ang_g0, double* ang_w0, double* C, int prl);
void compute_all_3c_v2(bool do_overlap, int natoms, int* atno, float* coords, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, int nrad, int nang, double* ang_g0, double* ang_w0, float* C, int prl);

void compute_all_3c_para(int ngpu, bool do_overlap, int natoms, int* atno, float* coords, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, int nrad, int nang, double* ang_g0, double* ang_w0, double* C, int prl);
void compute_all_3c_para(int ngpu, bool do_overlap, int natoms, int* atno, float* coords, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, int nrad, int nang, double* ang_g0, double* ang_w0, float* C, int prl);

void compute_VdV(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, int nc, float* coordsc, double* Pao, double* V, double* dV, int prl);
void compute_VdV(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, int nc, float* coordsc, double* Pao, float* V, float* dV, int prl);

void compute_all_4c(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* g, int prl);
void compute_all_4c(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, float* g, int prl);

void compute_all_4c_gen(int ngpu, int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* g, int prl);
void compute_all_4c_gen(int ngpu, int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, float* g, int prl);

void compute_all_4c_ol_gend(int ngpu, int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g, double* ang_w, double* ol, int prl);

void compute_all_4c_ol_gen(int ngpu, int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* ol, int prl);
void compute_all_4c_ol_gen(int ngpu, int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, float* ol, int prl);

//Nate's Edits
void compute_all_4c_v2(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* g, int prl);
void compute_all_4c_v2(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, float* g, int prl);

double compute_2c(int Z1, int Z2, float zeta10, float zeta20, float A20, float B20, float C20, int nrad, int nang, double* ang_g0, double* ang_w0, int prl);
double compute_1s_1s(int Z1, double zeta1, double zeta2, double A2, double B2, double C2, int nrad, int nang, double* ang_g, double* ang_w);

//float compute_1s_1s(int Z1, double zeta1, double zeta2, double A2, double B2, double C2, int nrad, int nang, double* ang_g, double* ang_w);


///// Prolate spheroidal coordinates /////
void compute_STEn_ps(int natoms, int* atno, double* coords, vector<vector<double> > &basis, int quad_order, int nmu, int nnu, int nphi, double* S, double* T, double* En, int prl);
void compute_pVp_ps(int natoms, int* atno, double* coords, vector<vector<double> > &basis, int quad_order, int nmu, int nnu, int nphi, double* pVp, int prl);
void compute_pVp_3c_ps(int natoms, int* atno, double* coords, vector<vector<double> > &basis, int quad_order, int quad_r_order, int nsplit, int nmu, int nnu, int nphi, double* pVp, int prl);
void compute_2c_ps(bool do_overlap, bool do_yukawa, double gamma, int natoms, int* atno, double* coords, vector<vector<double> > &basis, int quad_order, int nmu, int nnu, int nphi, double* A, int prl);
void compute_3c_ps(bool do_yukawa, double gamma, int natoms, int* atno, double* coords, vector<vector<double> > &basis, vector<vector<double> >& basis_aux, int quad_order, int quad_r_order, int nsplit, int nmu, int nnu, int nphi, double* En, double* C, int prl);
void compute_4c_ol_ps(int natoms, int* atno, double* coords, vector<vector<double> > &basis, int quad_order, int quad_r_order, int nmu, int nnu, int nphi, double* ol, int prl);
/////////////////////////////////////////

void get_inr_1s(int nrad, double zeta, double* r, double* inr);
void get_inr_1s_f(int nrad, float zeta, float* r, float* inr);
void get_inr_2s_f(int nrad, float zeta, float* r, float* inr);
void get_inr_2p_f(int nrad, float zeta, float* r, float* inr);
double norm_sv(int n, int l, int m, double zeta);
double norm(int n, int l, int m, double zeta);
size_t fact(size_t N);


int get_imax_n2i(int natoms, int N, vector<vector<double> >& basis, int* n2i);
int get_imax_n2ip(int Nmax, int natoms, int N, vector<vector<double> >& basis, vector<vector<int> >& n2ip);


int find_center_of_grid(float Z1, int nrad);
void generate_central_grid_2(float* grid1, float* wt1, float Z1, int nrad, int nang, float* ang_g, float* ang_w);
void acc_assign(int size, float* vec, float v1);
void acc_assign(int size, double* vec, double v1);
void acc_copy(int size, double* v1, double* v2);
void acc_copyf(int size, float* v1, float* v2);
void acc_copyf(int size, float* v1, float* v2, float* v3, float* v4);

void copy_grid(int gs, float* grid1, float* grid2);
void copy_grid(int gs, double* grid1, double* grid2);
void copy_grid(int gs, float* grid1, float* wt1, float* grid2, float* wt2);
void eliminate_small_wt(int size, float* wt1);
void eliminate_small_wt_3(int size, float* wt1, float* wt2, float* wt3);
void eliminate_small_wt(int s1, int size, float* wt1);
void eliminate_small_wt_3(int s1, int size, float* wt1, float* wt2, float* wt3);
void recenter_grid_zero(int gs, float* grid, float x2, float y2, float z2);
void recenter_grid_zero(int gs, double* grid, double x2, double y2, double z2);
void recenter_grid(int gs, float* grid, float x2, float y2, float z2);
void recenter_grid(int gs, double* grid, double x2, double y2, double z2);

void add_r1_to_grid(int gs, float* grid1, float A2, float B2, float C2);
void add_r1_to_grid(int gs, double* grid1, double A2, double B2, double C2);
void add_r2_to_grid(int gs, float* grid1, float A2, float B2, float C2);
void add_r2_to_grid(int gs, double* grid1, double A2, double B2, double C2);
void add_r3_to_grid(int gs, float* grid1, float A3, float B3, float C3);
void add_r123_to_grid(int gs, float* grid1, float A1, float B1, float C1, float A2, float B2, float C2, float A3, float B3, float C3);
void add_r1_to_grid_6z(int gs, float* grid1, float* grid2, float* grid3, float* grid4, float* grid5, float* grid6);

void becke_weight_2c(int gs, float* grid1, float* wt1, float* grid2, float* wt2, int Z1, int Z2, float A2, float B2, float C2);
void becke_weight_3c(int gs, float* grid1, float* wt1, float* grid2, float* wt2, float* grid3, float* wt3, int Z1, int Z2, int Z3, float A2, float B2, float C2, float A3, float B3, float C3);


//integrals_aux.cpp:
void print_array(int size, float* vec);
void clean_small_values(int N, float* S);
void clean_small_values(int N, double* S);
void acc_assign(int size, float* vec, float v1);
void acc_assign(int size, double* vec, double v1);
void acc_assign(int size, float* vec1, float* vec2, float v1);
void acc_assign(int size, float* vec1, float* vec2, float* vec3, float v1);
void acc_copyf(int tid, int size, float* v1, float* v2);
void acc_copyf(int size, float* v1, float* v2);
void acc_copyf(int size, float* v1, float* v2, float* v3, float* v4);
void acc_copyf(int size, float* v1, float* v2, float* v3, float* v4, float* v5, float* v6);
void eliminate_small_wt_3(int size, float* wt1, float* wt2, float* wt3);
void eliminate_small_wt_3(int s1, int size, float* wt1, float* wt2, float* wt3);
void eliminate_small_wt(int size, float* wt1);
void eliminate_small_wt(int s1, int size, float* wt1);
void copy_grid(int gs, float* grid1, float* grid2);
void copy_grid(int gs, float* grid1, float* wt1, float* grid2, float* wt2);
void recenter_grid_zero(int gs, float* grid, float x2, float y2, float z2);
void recenter_grid(int gs, float* grid, float x2, float y2, float z2);
void recenter_grid_exp(int gs, float* grid, float* wt, float* val, float x2, float y2, float z2, float zeta2);
void add_r123_to_grid(int gs, float* grid1, float A1, float B1, float C1, float A2, float B2, float C2, float A3, float B3, float C3);
void add_r1_to_grid_6z(int gs, float* grid1, float* grid2, float* grid3, float* grid4, float* grid5, float* grid6);
void add_r1_to_grid(int gs, float* grid1, float A2, float B2, float C2);
void add_r2_to_grid(int gs, float* grid1, float A2, float B2, float C2);
void add_r3_to_grid(int gs, float* grid1, float A3, float B3, float C3);
void generate_central_grid_2d(bool use_murak, double* grid1, double* wt1, float Z1, int nrad, int nang, double* ang_g, double* ang_w);
void generate_central_grid_2(float* grid1, float* wt1, float Z1, int nrad, int nang, float* ang_g, float* ang_w);
void generate_central_grid(float* grid1, float* wt1, float* val1, int need_inr, float Z1, int n1, int l1, float zeta1, int nrad, int nang, float* ang_g, float* ang_w);
void transpose_C(int Naux, int N, float* C);
void transpose_C(int Naux, int N, double* C);
void copy_symm(int natoms, int N, int Naux, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, float* C, int type);
void copy_symm(int natoms, int N, int Naux, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, double* C, int type);
void copy_symm_3c_ps(int natoms, int N, int Naux, int* n2i, int* na2i, double* C);
void copy_symm_4c_ps(int natoms, int* n2i, int N, double* olp);
int get_natoms_with_basis(int natoms, int* atno, vector<vector<double> >& basis);

#endif
