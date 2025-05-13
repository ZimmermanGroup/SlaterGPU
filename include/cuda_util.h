#ifndef CUDA_UTILH
#define CUDA_UTILH

//#include <cublas.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

int invert_eigen_cusolver(int size, double* A, double eig_max, cusolverDnHandle_t& cu_hdl);
int invert_stable_cusolver(int size, double* A, double delta, cusolverDnHandle_t& cu_hdl);
int invert_stable_cusolver_2(int size, double* A, double delta, cusolverDnHandle_t& cu_hdl);
int invert_cusolver(int N, float* A, cusolverDnHandle_t& cu_hdl);
int invert_cusolver(int N, double* A, cusolverDnHandle_t& cu_hdl);
void diagonalize_cusolver(int N, float* A, float* Ae, cusolverDnHandle_t& cu_hdl);
void diagonalize_cusolver(int Ne, int N, double* A, double* Ae, cusolverDnHandle_t& cu_hdl);

int mat_root_inv_cusolver(double* A, int size, cusolverDnHandle_t& cu_hdl);
int mat_root_inv_stable_cusolver(double* A, int size, double delta, cusolverDnHandle_t& cu_hdl);

void expmat(int N, double* theta1, double* U, cusolverDnHandle_t cu_hdl, cublasHandle_t cublasH);
//double expmat_complex(int N, double* theta1, double* theta1i, double* U, cusolverDnHandle_t cu_hdl, cublasHandle_t cublasH);
double expmat_complex(int N, double* theta1, double* theta1i, double* jCA, double* U, cusolverDnHandle_t cu_hdl, cublasHandle_t cublasH);

void mat_times_mat(float* C, float* A, float* B, int N);
void mat_times_mat(float* C, float* A, float* B, int M, int N, int K);
void mat_times_mat_at(float* C, float* A, float* B, int M, int N, int K);
void mat_times_mat_bt(float* C, float* A, float* B, int N);

void mat_times_mat(double* C, double* A, double* B, int N);
void mat_times_mat(double* C, double* A, double* B, int M, int N, int K);
void mat_times_mat_at(double* C, double* A, double* B, int N);
void mat_times_mat_at(double* C, double* A, double* B, int M, int N, int K);
void mat_times_mat_bt(double* C, double* A, double* B, int N);
void mat_times_mat_bt(double* C, double* A, double* B, int M, int N, int K);

#pragma acc routine seq
void trans(float* Bt, float* B, int m, int n);
#pragma acc routine seq
void trans(double* Bt, double* B, int m, int n);

void copy_to_all_gpu(int ngpu, int s1, double* A, int include_first);
void copy_to_all_gpu(int ngpu, int s1, float* A, int include_first);
void copy_to_all_gpu(int ngpu, int s1, int s2, double** A, int include_first);

#endif
