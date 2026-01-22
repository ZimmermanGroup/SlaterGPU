#ifndef CPU_UTILH
#define CPU_UTILH

#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>

#include <algorithm>
#include <math.h>
#include <complex>

#include <ctime>
#include <sys/time.h>

using namespace std;

//#define lapack_complex_float complex< float >
//#define lapack_complex_double complex< double >

//#include <lapacke.h>
//#include <cblas.h>

double randomf(double a, double b);

void solve_axeb(int dim, double* A, double* b);

void expmat_complex_cpu(int N, double* theta, double* thetai, double* etheta);
void expmat_cpu(int N, double *theta, double *etheta);

int la_diagR(int neig, double* A, double* eigen, double* eigeni);
void la_diag(int neig, int s1, double* A, double* Ae);
int invert_stable_cpu(double* A, int size, double delta);
int invert_stable_cpu(double* A, int size, double delta, int prl);
int invert_stable_cpu_2(double* A, int size, double delta);
int invert_stable_cpu_3(double* A, int size, double delta, bool root);
int invert_cpu(double* A, int size, int mode);
int mat_root_cpu(double* A, int size);
int mat_root_inv_cpu(double* A, int size);
int mat_root_inv_stable_cpu(double* A, int size, double inv_cutoff, int prl);
int LU_inv_stable_cpu(double* A, int size);

void trans_cpu(float* Bt, float* B, int m, int n);
void trans_cpu(double* Bt, double* B, int m, int n);
void cross(double* m, double* r1, double* r2);
int sign(double val);

void mat_times_mat_cpu(float* C, float* A, float* B, int M, int N, int K);
void mat_times_mat_cpu(float* C, float* A, float* B, int N);
void mat_times_mat_cpu(double* C, double* A, double* B, int M, int N, int K);
void mat_times_mat_ct_cpu(double* C, double* A, double* B, int M, int N, int K);
void mat_times_mat_cpu(double* C, double* A, double* B, int N);
void mat_times_mat_at_cpu(double* C, double* A, double* B, int N);
void mat_times_mat_bt_cpu(float* C, float* A, float* B, int M, int N, int K);
void mat_times_mat_bt_cpu(double* C, double* A, double* B, int M, int N, int K, int LDAB);
void mat_times_mat_bt_cpu(float* C, float* A, float* B, int N);
void mat_times_mat_bt_cpu(double* C, double* A, double* B, int M, int N, int K);
void mat_times_mat_bt_cpu(double* C, double* A, double* B, int N);
void mat_times_mat_at_cpu(double* C, double* A, double* B, int M, int N, int K);

#endif
