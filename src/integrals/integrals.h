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

#include "fp_def.h"
#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include <chrono>

#if USE_ACC
#include "accel.h"
#endif

#include <omp.h>

using namespace std;

#include "spherical.h"
#include "Vinr.h"
#include "pVp.h"
#include "reduce.h"
#include "murak.h"

#if TEST_SORT
#include "lib/sortUtil.h"
extern void sortByKeys( fp_t * A, xyz<fp_t> * Av, size_t N);
#endif

#define PI 3.1415926535

FP2 simple_test(int size1, int size2);

FP1 compute_2c_trial(FP1 zeta1, FP1 zeta2, int nrad, int nang, FP2* ang_g, FP2* ang_w);

void compute_all_3c_para2(int ngpu, int natoms, int* atno, FP1* coords, 
                         vector<vector<FP2> > &basis, 
                         vector<vector<FP2> > &basis_aux, 
                         int nrad, int nang, FP2* ang_g0, FP2* ang_w0, 
                         FP2* C, int prl);
void compute_all_3c_para3(int ngpu, int natoms, int* atno, FP1* coords, 
                         vector<vector<FP2> > &basis, 
                         vector<vector<FP2> > &basis_aux, 
                         int nrad, int nang, FP2* ang_g0, FP2* ang_w0, 
                         FP2* C, int prl);
void compute_d_3c_para2(
  int ngpu, int natoms, int* atno, FP1* coords,
  vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, 
  int nrad, int nang, FP2* ang_g0, FP2* ang_w0, 
  FP2* dC, FP2* xyz_grad, int prl);

void compute_d_ST(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* GFao, FP2* Pao, FP2* xyz_grad, int prl);
// void compute_d_ST(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* GFao, FP2* Pao, FP1* xyz_grad, int prl);

void compute_d_En(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* Pao, FP2* xyz_grad, int prl);
// void compute_d_En(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* Pao, FP1* xyz_grad, int prl);

void compute_d_2c(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* dpq, FP2* xyz_grad, int prl);
// void compute_d_2c(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* dpq, FP1* xyz_grad, int prl);

void compute_d_3c(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* dC, FP2* xyz_grad, int prl);
// void compute_d_3c(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* dC, FP1* xyz_grad, int prl);

void compute_d_3c_para(int ngpu, int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* dC, FP2* xyz_grad, int prl);
// void compute_d_3c_para(int ngpu, int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* dC, FP1* xyz_grad, int prl);


void compute_Enp(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* En, FP2* pVp, int prl);
// void compute_Enp(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP1* En, FP1* pVp, int prl);

void compute_Enp_para(int ngpu, int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* En, FP2* pVp, int prl);
// void compute_Enp_para(int ngpu, int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP1* En, FP1* pVp, int prl);

void compute_Sd(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* S, int prl);

void compute_ST(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* S, FP2* T, int prl);
// void compute_ST(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP1* S, FP1* T, int prl);

void compute_all_2c(int natoms, int* atno, FP1* coordsf, vector<vector<FP2> > &basis_aux, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP1* A, int prl);

void compute_all_2c_v2(int natoms, int* atno, FP1* coordsf, vector<vector<FP2> > &basis_aux, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* A, int prl);
// void compute_all_2c_v2(int natoms, int* atno, FP1* coordsf, vector<vector<FP2> > &basis_aux, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP1* A, int prl);

void compute_all_3c(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP1* C, int prl);

void compute_all_3c_v2(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* C, int prl);
// void compute_all_3c_v2(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP1* C, int prl);

void compute_all_3c_para(int ngpu, int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* C, int prl);
// void compute_all_3c_para(int ngpu, int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP1* C, int prl);

void compute_VdV(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, int nc, FP1* coordsc, FP2* Pao, FP2* V, FP2* dV, int prl);
// void compute_VdV(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, int nc, FP1* coordsc, FP2* Pao, FP1* V, FP1* dV, int prl);



FP2 compute_2c(int Z1, int Z2, FP1 zeta10, FP1 zeta20, FP1 A20, FP1 B20, FP1 C20, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, int prl);
FP2 compute_1s_1s(int Z1, FP2 zeta1, FP2 zeta2, FP2 A2, FP2 B2, FP2 C2, int nrad, int nang, FP2* ang_g, FP2* ang_w);
//FP1 compute_1s_1s(int Z1, FP2 zeta1, FP2 zeta2, FP2 A2, FP2 B2, FP2 C2, int nrad, int nang, FP2* ang_g, FP2* ang_w);


void get_inr_1s(int nrad, FP2 zeta, FP2* r, FP2* inr);
void get_inr_1s_f(int nrad, FP1 zeta, FP1* r, FP1* inr);
void get_inr_2s_f(int nrad, FP1 zeta, FP1* r, FP1* inr);
void get_inr_2p_f(int nrad, FP1 zeta, FP1* r, FP1* inr);
FP2 norm_sv(int n, int l, int m, FP2 zeta);
FP2 norm(int n, int l, int m, FP2 zeta);
int fact(int N);


int get_imax_n2i(int natoms, int N, vector<vector<FP2> >& basis, int* n2i);


int find_center_of_grid(FP1 Z1, int nrad);
void generate_central_grid_2(FP1* grid1, FP1* wt1, FP1 Z1, int nrad, int nang, FP1* ang_g, FP1* ang_w);
void acc_assign(int size, FP1* vec, FP1 v1);
void acc_assign(int size, FP2* vec, FP2 v1);
void acc_copyf(int size, FP1* v1, FP1* v2);
void acc_copyf(int size, FP1* v1, FP1* v2, FP1* v3, FP1* v4);

void copy_grid(int gs, FP1* grid1, FP1* grid2);
void copy_grid(int gs, FP1* grid1, FP1* wt1, FP1* grid2, FP1* wt2);
void eliminate_small_wt(int size, FP1* wt1);
void eliminate_small_wt_3(int size, FP1* wt1, FP1* wt2, FP1* wt3);
void eliminate_small_wt(int s1, int size, FP1* wt1);
void eliminate_small_wt_3(int s1, int size, FP1* wt1, FP1* wt2, FP1* wt3);
void recenter_grid_zero(int gs, FP1* grid, FP1 x2, FP1 y2, FP1 z2);

void add_r1_to_grid(int gs, FP1* grid1, FP1 A2, FP1 B2, FP1 C2);
void add_r2_to_grid(int gs, double* grid1, double A2, double B2, double C2);
void add_r3_to_grid(int gs, FP1* grid1, FP1 A3, FP1 B3, FP1 C3);
void add_r123_to_grid(int gs, FP1* grid1, FP1 A1, FP1 B1, FP1 C1, FP1 A2, FP1 B2, FP1 C2, FP1 A3, FP1 B3, FP1 C3);
void add_r1_to_grid_6z(int gs, FP1* grid1, FP1* grid2, FP1* grid3, FP1* grid4, FP1* grid5, FP1* grid6);

void becke_weight_2c(int gs, FP1* grid1, FP1* wt1, FP1* grid2, FP1* wt2, int Z1, int Z2, FP1 A2, FP1 B2, FP1 C2);
void becke_weight_3c(int gs, FP1* grid1, FP1* wt1, FP1* grid2, FP1* wt2, FP1* grid3, FP1* wt3, int Z1, int Z2, int Z3, FP1 A2, FP1 B2, FP1 C2, FP1 A3, FP1 B3, FP1 C3);
FP1 becke_a(int Z1, int Z2);

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
void copy_grid(int gs, double* grid1, double* grid2);
void copy_grid(int gs, float* grid1, float* wt1, float* grid2, float* wt2);
void recenter_grid_zero(int gs, float* grid, float x2, float y2, float z2);
void recenter_grid_zero(int gs, double* grid, double x2, double y2, double z2);
void recenter_grid(int gs, float* grid, float x2, float y2, float z2);
void recenter_grid(int gs, double* grid, double x2, double y2, double z2);
void recenter_grid_exp(int gs, float* grid, float* wt, float* val, float x2, float y2, float z2, float zeta2);
void add_r123_to_grid(int gs, float* grid1, float A1, float B1, float C1, float A2, float B2, float C2, float A3, float B3, float C3);
void add_r1_to_grid_6z(int gs, float* grid1, float* grid2, float* grid3, float* grid4, float* grid5, float* grid6);
void add_r1_to_grid(int gs, float* grid1, float A2, float B2, float C2);
void add_r1_to_grid(int gs, double* grid1, double A2, double B2, double C2);
void add_r2_to_grid(int gs, float* grid1, float A2, float B2, float C2);
//void add_r2_to_grid(int gs, double* grid1, float A2, float B2, float C2);
void add_r3_to_grid(int gs, float* grid1, float A3, float B3, float C3);
void generate_central_grid_2d(bool use_murak, double* grid1, double* wt1, float Z1, int nrad, int nang, double* ang_g, double* ang_w);
void generate_central_grid_2(float* grid1, float* wt1, float Z1, int nrad, int nang, float* ang_g, float* ang_w);
void generate_central_grid(float* grid1, float* wt1, float* val1, int need_inr, float Z1, int n1, int l1, float zeta1, int nrad, int nang, float* ang_g, float* ang_w);
void transpose_C(int Naux, int N, float* C);
void transpose_C(int Naux, int N, double* C);
void copy_symm(int natoms, int N, int Naux, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, float* C, int type);
void copy_symm(int natoms, int N, int Naux, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, double* C, int type);
void becke_weight_2d(int gs, double* grid1, double* wt1, double* grid2, double* wt2,
                     double zeta1, double zeta2, double A2, double B2, double C2);
double bf3d(double f1);

#endif
