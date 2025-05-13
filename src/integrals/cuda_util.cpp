#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <algorithm>
#include <math.h>
#include <chrono>
#include <ctime>
#include "accel.h"
#include "cuda_util.h"
#include "cpu_util.h"

#include "cuda_runtime.h"

//these functions assume the matrices are already on the device, via ACC
//hoping to reduce overhead by getting cusolver handle just once


//CPMZ: incorporate collapse into mat_times_mat calls


void print_square_sm(int N, float* A);
void print_square_sm(int N, double* A);
double read_float(string filename);

#include <cuComplex.h>

void print_square_complex(int N, cuDoubleComplex* A)
{
  for (int i=0;i<N;i++)
  {
    for (int j=0;j<N;j++)
      printf("  %8.5f + %8.5fi",A[i*N+j].x,A[i*N+j].y);
    printf("\n");
  }
}

//void expmat_complex(int N, double* theta1, double* theta1i, double* U, cusolverDnHandle_t cu_hdl, cublasHandle_t cublasH)
double expmat_complex(int N, double* theta1, double* theta1i, double* jCA, double* U, cusolverDnHandle_t cu_hdl, cublasHandle_t cublasH)
{
 //this is a bit of a mess
 //takes exp(-i theta)

  int N2 = N*N;

  cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR;
  cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
  //cusolverEigRange_t range = CUSOLVER_EIG_RANGE_ALL;

  for (int i=0;i<N;i++)
  for (int j=0;j<i;j++)
  if (fabs(theta1[i*N+j]+theta1[j*N+i])>1.e-8)
    printf(" ERROR: asymmetry: %12.10f %12.10f \n",theta1[i*N+j],theta1[j*N+i]);

  int Lwork = 1;

 //MO coeffs onto device
  double* jCAi = &jCA[N2];
  cuDoubleComplex* CA = new cuDoubleComplex[N2];
  for (int i=0;i<N2;i++)
  {
    CA[i].x = jCA[i];
    CA[i].y = jCAi[i];
  }
  cuDoubleComplex* d_CA;
  cudaMalloc((void**)&d_CA,N2*sizeof(cuDoubleComplex));
  cudaMemcpy(d_CA,CA,sizeof(cuDoubleComplex)*N2,cudaMemcpyHostToDevice);

 //get theta1 onto device as complex matrix
  cuDoubleComplex* A1 = new cuDoubleComplex[N2];
  cuDoubleComplex* A;
  cudaMalloc((void**)&A,N2*sizeof(cuDoubleComplex));
  double* W;
  cudaMalloc ((void**)&W,sizeof(double)*N2);
  cuDoubleComplex* work;
  int* info;
  cudaMalloc((void**) &info,sizeof(int));

  for (int i=0;i<N2;i++)
  {
    A1[i].x = theta1i[i];
    A1[i].y = -theta1[i];
  }
  //printf("\n theta(complex) \n");
  //print_square_complex(N,A1);

  cudaMemcpy(A,A1,sizeof(cuDoubleComplex)*N2,cudaMemcpyHostToDevice);

  cusolverDnZheevd_bufferSize(cu_hdl,jobz,uplo,N,A,N,W,&Lwork);
  //cusolverDnZheevdx_bufferSize(cu_hdl,jobz,range,uplo,N,A,N,0,0,0,0,neig,W,&Lwork);

  //printf(" Lwork: %i \n",Lwork);
  cudaMalloc((void**)&work,Lwork*sizeof(cuDoubleComplex));

  cusolverDnZheevd(cu_hdl,jobz,uplo,N,A,N,W,work,Lwork,info);
  //cusolverDnZheevdx(cu_hdl,jobz,range,uplo,N,A,N,0,0,0,0,neig,W,work,Lwork,info);
  cudaDeviceSynchronize();

  double* Ue = new double[N];
  cudaMemcpy(Ue,W,sizeof(double)*N,cudaMemcpyDeviceToHost);

  cuDoubleComplex* d_eix;
  cudaMalloc((void**)&d_eix,N2*sizeof(cuDoubleComplex));

  cuDoubleComplex* eix = new cuDoubleComplex[N2];
  for (int i=0;i<N2;i++)
    eix[i].x = eix[i].y = 0.;
  for (int i=0;i<N;i++)
  {
    double c1 = cos(Ue[i]);
    double c1i = sin(Ue[i]);
    eix[i*N+i].x = c1;
    eix[i*N+i].y = c1i;
  }
  cudaMemcpy(d_eix,eix,sizeof(cuDoubleComplex)*N2,cudaMemcpyHostToDevice);

  cuDoubleComplex alpha; alpha.x = 1.; alpha.y = 0.;
  cuDoubleComplex beta; beta.x = 0.; beta.y = 0.;

  cublasOperation_t cu_notrans = CUBLAS_OP_N;
  cublasOperation_t cu_transC = CUBLAS_OP_C; //conjugate transpose
  cublasOperation_t cu_trans = CUBLAS_OP_T; //transpose

  cuDoubleComplex* Aeix;
  cudaMalloc((void**)&Aeix,N2*sizeof(cuDoubleComplex));

  //need A*eix*A(T)
  cublasZgemm(cublasH,cu_notrans,cu_notrans,N,N,N,&alpha,A,N,d_eix,N,&beta,Aeix,N);
  cublasZgemm(cublasH,cu_notrans,cu_transC,N,N,N,&alpha,Aeix,N,A,N,&beta,d_eix,N);

 #if 0
  printf("\n UU*(complex) \n");
  cublasZgemm(cublasH,cu_notrans,cu_transC,N,N,N,&alpha,d_eix,N,d_eix,N,&beta,A,N);
  cudaMemcpy(CA,A,sizeof(cuDoubleComplex)*N2,cudaMemcpyDeviceToHost);
  print_square_complex(N,CA);
 #endif

  cuDoubleComplex* U1 = new cuDoubleComplex[N2];
  cudaMemcpy(U1,d_eix,sizeof(cuDoubleComplex)*N2,cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();

 //CPMZ cannot use this in current form
  double ss = 0.;
  for (int i=0;i<N;i++)
  for (int j=0;j<i;j++)
  {
    double x1 = U1[i*N+j].x;
    double y1 = U1[i*N+j].y;
    ss += x1*x1+y1*y1;
  }
  ss = 0.;

 #if 1
  printf(" Zheevd eigenvalues:");
  for (int i=0;i<N;i++)
    printf(" %8.5f",Ue[i]);
  printf("\n");
 #endif

  for (int i=0;i<N2;i++)
    U[i] = U1[i].x;

 #if 0
  printf("\n U(complex) \n");
  print_square_complex(N,U1);
 #endif


 //reuse A to store rotated CA
  cudaDeviceSynchronize();
  cublasZgemm(cublasH,cu_trans,cu_notrans,N,N,N,&alpha,d_CA,N,d_eix,N,&beta,A,N);
  cudaDeviceSynchronize();

 //copy back MO coeffs
  cudaMemcpy(CA,A,sizeof(cuDoubleComplex)*N2,cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
  {
    jCA[j*N+i] = CA[i*N+j].x;
    jCAi[j*N+i] = CA[i*N+j].y;
  }

  //cudaFree(exi);
  cudaFree(d_CA);
  cudaFree(A);
  cudaFree(work);
  cudaFree(info);
  cudaFree(W);
  cudaFree(Aeix);
  cudaFree(d_eix);
  delete [] CA;
  delete [] A1;
  delete [] Ue;
  delete [] eix;
  delete [] U1;

  return ss;
}

void expmat(int N, double* theta1, double* U, cusolverDnHandle_t cu_hdl, cublasHandle_t cublasH)
{
 //this is a bit of a mess
 //takes exp(-i theta)

  int N2 = N*N;

  cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR;
  cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
  //cusolverEigRange_t range = CUSOLVER_EIG_RANGE_ALL;

  for (int i=0;i<N;i++)
  for (int j=0;j<i;j++)
  if (fabs(theta1[i*N+j]+theta1[j*N+i])>1.e-8)
    printf(" ERROR: asymmetry: %12.10f %12.10f \n",theta1[i*N+j],theta1[j*N+i]);

  int Lwork = 1;

 //get theta1 onto device as complex matrix
  cuDoubleComplex* A1 = new cuDoubleComplex[N2];
  cuDoubleComplex* A;
  cudaMalloc((void**)&A,N2*sizeof(cuDoubleComplex));
  double* W;
  cudaMalloc ((void**)&W,sizeof(double)*N2);
  cuDoubleComplex* work;
  int* info;
  cudaMalloc((void**) &info,sizeof(int));

  for (int i=0;i<N2;i++)
  {
    A1[i].x = 0.;
    A1[i].y = -theta1[i];
  }
  cudaMemcpy(A,A1,sizeof(cuDoubleComplex)*N2,cudaMemcpyHostToDevice);

  cusolverDnZheevd_bufferSize(cu_hdl,jobz,uplo,N,A,N,W,&Lwork);
  //cusolverDnZheevdx_bufferSize(cu_hdl,jobz,range,uplo,N,A,N,0,0,0,0,neig,W,&Lwork);

  //printf(" Lwork: %i \n",Lwork);
  cudaMalloc((void**)&work,Lwork*sizeof(cuDoubleComplex));

  cusolverDnZheevd(cu_hdl,jobz,uplo,N,A,N,W,work,Lwork,info);
  //cusolverDnZheevdx(cu_hdl,jobz,range,uplo,N,A,N,0,0,0,0,neig,W,work,Lwork,info);
  cudaDeviceSynchronize();

  double* Ue = new double[N];
  cudaMemcpy(Ue,W,sizeof(double)*N,cudaMemcpyDeviceToHost);

  cuDoubleComplex* d_eix;
  cudaMalloc((void**)&d_eix,N2*sizeof(cuDoubleComplex));

  cuDoubleComplex* eix = new cuDoubleComplex[N2];
  for (int i=0;i<N2;i++)
    eix[i].x = eix[i].y = 0.;
  for (int i=0;i<N;i++)
  {
    double c1 = cos(Ue[i]);
    double c1i = sin(Ue[i]);
    eix[i*N+i].x = c1;
    eix[i*N+i].y = c1i;
  }
  cudaMemcpy(d_eix,eix,sizeof(cuDoubleComplex)*N2,cudaMemcpyHostToDevice);

  cuDoubleComplex alpha; alpha.x = 1.; alpha.y = 0.;
  cuDoubleComplex beta; beta.x = 0.; beta.y = 0.;

  cublasOperation_t cu_notrans = CUBLAS_OP_N;
  cublasOperation_t cu_trans = CUBLAS_OP_C;
  //cublasOperation_t cu_transC = CUBLAS_OP_T; //conjugate trans

  cuDoubleComplex* Aeix;
  cudaMalloc((void**)&Aeix,N2*sizeof(cuDoubleComplex));

  //need A*eix*A(T)
  cublasZgemm(cublasH,cu_notrans,cu_notrans,N,N,N,&alpha,A,N,d_eix,N,&beta,Aeix,N);
  cublasZgemm(cublasH,cu_notrans,cu_trans,N,N,N,&alpha,Aeix,N,A,N,&beta,d_eix,N);

  cuDoubleComplex* U1 = new cuDoubleComplex[N2];
  cudaMemcpy(U1,d_eix,sizeof(cuDoubleComplex)*N2,cudaMemcpyDeviceToHost);

 #if 0
  printf(" Zheevd eigenvalues:");
  for (int i=0;i<N;i++)
    printf(" %8.5f",Ue[i]);
  printf("\n");
 #endif

  for (int i=0;i<N2;i++)
    U[i] = U1[i].x;

#if 0
  printf("\n U1(vec) \n");
  print_square_complex(N,U1);
#endif

  cudaFree(A);
  cudaFree(work);
  cudaFree(info);
  cudaFree(W);
  cudaFree(Aeix);
  cudaFree(d_eix);
  delete [] A1;
  delete [] Ue;
  delete [] eix;
  delete [] U1;

  return;
}

void mat_times_mat(float* C, float* A, float* B, int M, int N, int K)
{
  int MN = M*N;
  int MK = M*K;
  int NK = N*K;

 //#pragma acc parallel loop present(C[0:MN])
 // for (int i=0;i<MN;i++)
 //   C[i] = 0.;

 #pragma acc parallel loop collapse(2) present(A[0:MK],B[0:NK],C[0:MN])
  for (int i=0;i<M;i++)
  {
   //#pragma acc loop
    for (int j=0;j<N;j++)
    {
      float val = 0.;
     #pragma acc loop reduction(+:val)
      for (int k=0;k<K;k++)
        val += A[i*K+k]*B[k*N+j];
      C[i*N+j] = val;
    }
  }

  return;
}

void mat_times_mat(float* C, float* A, float* B, int N)
{
  return mat_times_mat(C,A,B,N,N,N);
}

void mat_times_mat(double* C, double* A, double* B, int M, int N, int K)
{
  int MN = M*N;
  int MK = M*K;
  int NK = N*K;

 //#pragma acc parallel loop present(C[0:MN])
 // for (int i=0;i<MN;i++)
 //   C[i] = 0.;

 #pragma acc parallel loop present(A[0:MK],B[0:NK],C[0:MN])
  for (int i=0;i<M;i++)
  {
   #pragma acc loop
    for (int j=0;j<N;j++)
    {
      double val = 0.;
     #pragma acc loop reduction(+:val)
      for (int k=0;k<K;k++)
        val += A[i*K+k]*B[k*N+j];
      C[i*N+j] = val;
    }
  }

  return;
}

void mat_times_mat(double* C, double* A, double* B, int N)
{
 #if !USE_ACC
  return mat_times_mat_cpu(C,A,B,N);
 #endif
  return mat_times_mat(C,A,B,N,N,N);
}

void mat_times_mat_at(double* C, double* A, double* B, int M, int N, int K)
{
  int MN = M*N;
  int MK = M*K;
  int NK = N*K;

 //#pragma acc parallel loop present(C[0:MN])
 // for (int i=0;i<MN;i++)
 //   C[i] = 0.;

 #pragma acc parallel loop collapse(2) present(A[0:MK],B[0:NK],C[0:MN])
  for (int i=0;i<M;i++)
  {
   //#pragma acc loop
    for (int j=0;j<N;j++)
    {
      double val = 0.;
     #pragma acc loop reduction(+:val)
      for (int k=0;k<K;k++)
        val += A[k*M+i]*B[k*N+j];
      C[i*N+j] = val;
    }
  }

  return;
}

void mat_times_mat_at(double* C, double* A, double* B, int N)
{
 #if !USE_ACC
  return mat_times_mat_at_cpu(C,A,B,N);
 #endif
  return mat_times_mat_at(C,A,B,N,N,N);
}

void mat_times_mat_bt(double* C, double* A, double* B, int M, int N, int K)
{
  int MN = M*N;
  int MK = M*K;
  int NK = N*K;

 //#pragma acc parallel loop present(C[0:MN])
 // for (int i=0;i<MN;i++)
 //   C[i] = 0.;

 #pragma acc parallel loop collapse(2) present(A[0:MK],B[0:NK],C[0:MN])
  for (int i=0;i<M;i++)
  {
   //#pragma acc loop
    for (int j=0;j<N;j++)
    {
      double val = 0.;
     #pragma acc loop reduction(+:val)
      for (int k=0;k<K;k++)
        val += A[i*K+k]*B[j*K+k];
      C[i*N+j] = val;
    }
  }

  return;
}

void mat_times_mat_bt(double* C, double* A, double* B, int N)
{
 #if !USE_ACC
  return mat_times_mat_bt_cpu(C,A,B,N);
 #endif
  return mat_times_mat_bt(C,A,B,N,N,N);
}

#pragma acc routine seq
void trans(float* Bt, float* B, int m, int n)
{
  for (int i=0;i<m;i++)
  for (int j=0;j<n;j++)
    Bt[i*n+j] = B[j*m+i];

  return;
}

#pragma acc routine seq
void trans(double* Bt, double* B, int m, int n)
{
  for (int i=0;i<m;i++)
  for (int j=0;j<n;j++)
    Bt[i*n+j] = B[j*m+i];

  return;
}

int mat_root_inv_stable_cusolver(double* A, int size, double delta, cusolverDnHandle_t& cu_hdl)
{
  int s2 = size*size;

  double* B = new double[s2];
  double* Bi = new double[s2];
  double* Beigen = new double[size];
  double* tmp = new double[s2]();

  #pragma acc enter data create(B[0:s2],Bi[0:s2],Beigen[0:size],tmp[0:s2])

  #pragma acc serial present(B[0:s2],A[0:s2],Beigen[0:size])
  {
    for (int i=0;i<s2;i++) B[i] = A[i];
    for (int i=0;i<size;i++) Beigen[i] = 0.;
  }

  diagonalize_cusolver(size,size,B,Beigen,cu_hdl);

 #if 1
  #pragma acc update self(Beigen[0:size])
  printf(" mat_root_inv_stable eigenvalues:");
  for (int i=0;i<std::min(10,size);i++)
  if (fabs(Beigen[i])<1e10)
    printf(" %5.2e",Beigen[i]);
  else
    printf(" large");
  printf("\n");
 #endif

  int nsmall = 0;
 #pragma acc parallel loop present(Beigen[0:size]) reduction(+:nsmall)
  for (int i=0;i<size;i++)
  if (Beigen[i]<delta)
  {
    Beigen[i] = 1.e20;
    nsmall++;
  }
  printf("  found %i small eigen \n",nsmall);

  if (nsmall>0)
  {
    printf(" lowest vector:");
    #pragma acc update self(B[0:s2])

    if (size<40)
    {
      for (int i=0;i<size;i++)
        printf(" %6.3f",B[i]);
    }
    else
    {
      for (int i=0;i<20;i++)
        printf(" %6.3f",B[i]);
      printf(" ... ");
      for (int i=size-20;i<size;i++)
        printf(" %6.3f",B[i]);
    }

    printf("\n");
  }


 #pragma acc serial present(Bi[0:s2],B[0:s2])
  trans(Bi,B,size,size);

 #pragma acc parallel loop present(tmp[0:s2])
  for (int i=0;i<s2;i++) tmp[i] = 0.;

 #pragma acc parallel loop present(tmp[0:s2],Bi[0:s2],Beigen[0:size])
  for (int i=0;i<size;i++)
 #pragma acc loop
  for (int j=0;j<size;j++)
    tmp[i*size+j] += Bi[i*size+j] / sqrt(Beigen[j]);
  mat_times_mat(A,tmp,B,size);

  #pragma acc exit data delete(B[0:s2],Bi[0:s2],Beigen[0:size],tmp[0:s2])

  delete [] B;
  delete [] Beigen;
  delete [] Bi;
  delete [] tmp;

  return nsmall;
}

int mat_root_inv_cusolver(double* A, int size, cusolverDnHandle_t& cu_hdl)
{
  int s2 = size*size;

  double* B = new double[s2];
  double* Bi = new double[s2];
  double* Beigen = new double[size];
  double* tmp = new double[s2]();

  #pragma acc enter data create(B[0:s2],Bi[0:s2],Beigen[0:size],tmp[0:s2])

  #pragma acc serial present(B[0:s2],A[0:s2],Beigen[0:size])
  {
    for (int i=0;i<s2;i++) B[i] = A[i];
    for (int i=0;i<size;i++) Beigen[i] = 0.;
  }

  diagonalize_cusolver(size,size,B,Beigen,cu_hdl);

 #if 1
  #pragma acc update self(Beigen[0:size])
  printf(" mat_root_inv eigenvalues:");
  for (int i=0;i<std::min(10,size);i++)
  if (fabs(Beigen[i])<1e8)
    printf(" %5.2e",Beigen[i]);
  else
    printf(" large");
  printf("\n");
 #endif

 #pragma acc serial present(Bi[0:s2],B[0:s2])
  trans(Bi,B,size,size);

 #pragma acc parallel loop present(tmp[0:s2])
  for (int i=0;i<s2;i++) tmp[i] = 0.;

  double thresh = 1.e-8;
  int nsmall = 0;
  for (int i=0;i<size;i++)
  if (Beigen[i]<thresh)
    nsmall++;

 #pragma acc parallel loop present(tmp[0:s2],Bi[0:s2],Beigen[0:size])
  for (int i=0;i<size;i++)
 #pragma acc loop
  for (int j=0;j<size;j++)
    tmp[i*size+j] += Bi[i*size+j] / sqrt(Beigen[j]);
  mat_times_mat(A,tmp,B,size);

  #pragma acc exit data delete(B[0:s2],Bi[0:s2],Beigen[0:size],tmp[0:s2])

  delete [] B;
  delete [] Beigen;
  delete [] Bi;
  delete [] tmp;

  return nsmall;
}


void diagonalize_cusolver(int N, float* A, float* Ae, cusolverDnHandle_t& cu_hdl)
{
  cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR;
  cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;

  int Lwork = 1;

 //is the &Lwork okay?
 #pragma acc host_data use_device(A,Ae)
  cusolverDnSsyevd_bufferSize(cu_hdl,jobz,uplo,N,A,N,Ae,&Lwork);

  int* info = new int[1];
  float* work = new float[Lwork];
  #pragma acc enter data create(work[0:Lwork],info[0:1])

 #pragma acc host_data use_device(A,Ae,work,info)
  cusolverDnSsyevd(cu_hdl,jobz,uplo,N,A,N,Ae,work,Lwork,info);

  cudaDeviceSynchronize();

 #pragma acc exit data delete(work[0:Lwork],info[0:1])

#if 0
  #pragma acc update self(A[0:N*N],Ae[0:N])
  printf("  eigenvalues: ");
  for (int i=0;i<N;i++)
    printf(" %6.3f",Ae[i]);
  printf("\n");
  printf("  eigenvectors: \n");
  print_square_sm(N,A);
#endif

  delete [] work;
  delete [] info;

  return;
}

void diagonalize_cusolver(int Ne, int N, double* A, double* Ae, cusolverDnHandle_t& cu_hdl)
{
  cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR;
  cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;

 //for syevdx
  cusolverEigRange_t rtype = CUSOLVER_EIG_RANGE_I;
  double vl = 0; double vu = 0;
  int il = 1;
  int iu = std::min(N,Ne);
  int nf;

  int Lwork = 1;

 //#pragma acc host_data use_device(A,Ae)
  //cusolverDnDsyevd_bufferSize(cu_hdl,jobz,uplo,N,A,N,Ae,&Lwork);

 #pragma acc host_data use_device(A,Ae)
  cusolverDnDsyevdx_bufferSize(cu_hdl,jobz,rtype,uplo,N,A,N,vl,vu,il,iu,&nf,Ae,&Lwork);

  //printf(" Lwork: %2i \n",Lwork);

  int* info = new int[1];
  double* work = new double[Lwork];
  #pragma acc enter data create(work[0:Lwork],info[0:1])

 #pragma acc host_data use_device(A,Ae,work,info)
  cusolverDnDsyevdx(cu_hdl,jobz,rtype,uplo,N,A,N,vl,vu,il,iu,&nf,Ae,work,Lwork,info);
  //cusolverDnDsyevd(cu_hdl,jobz,uplo,N,A,N,Ae,work,Lwork,info);

  //printf("  after syevdx \n");

  cudaDeviceSynchronize();

 #pragma acc exit data delete(work[0:Lwork],info[0:1])

#if 0
  #pragma acc update self(A[0:N*N])
  printf("  here update \n");
  #pragma acc update self(Ae[0:N])
  printf("  eigenvalues: ");
  for (int i=0;i<N;i++)
    printf(" %6.3f",Ae[i]);
  printf("\n");
  printf("  eigenvectors: \n");
  print_square_sm(N,A);
#endif

  delete [] work;
  delete [] info;

  return;
}

int invert_eigen_cusolver(int size, double* A, double eig_max, cusolverDnHandle_t& cu_hdl)
{
  int s2 = size*size;

  double* B = new double[s2];
  double* Bi = new double[s2];
  double* Ae = new double[size];
  double* tmp = new double[s2]();

  #pragma acc enter data create(B[0:s2],Bi[0:s2],Ae[0:size],tmp[0:s2])

  #pragma acc serial present(B[0:s2],A[0:s2])
   for (int i=0;i<s2;i++)
     B[i] = A[i];

  diagonalize_cusolver(size,size,B,Ae,cu_hdl);

 #pragma acc serial present(Bi[0:s2],B[0:s2])
  trans(Bi,B,size,size);

 #if 0
  #pragma acc update self(Ae[0:size])
  printf(" invert_eigen eigenvalues:");
  if (size>8)
  {
    for (int i=0;i<std::min(4,size);i++)
      printf(" %6.2f",Ae[i]);
    printf(" ...");
    for (int i=std::max(0,size-4);i<size;i++)
      printf(" %6.2f",Ae[i]);
  }
  else
  {
    for (int i=0;i<size;i++)
      printf(" %6.2f",Ae[i]);
  }
  printf("\n");
 #endif

  int reset = 0;
 #pragma acc parallel loop independent present(Ae[0:size]) reduction(+:reset)
  for (int i=0;i<size;i++)
  if (Ae[i]>eig_max)
  {
    reset++;
  }

  if (reset>0)
    printf("  found %i large eigen \n",reset);

  if (!reset)
  {
   #pragma acc parallel loop present(tmp[0:s2])
    for (int i=0;i<s2;i++)
      tmp[i] = 0.;

   #pragma acc parallel loop present(tmp[0:s2],Bi[0:s2],Ae[0:size])
    for (int i=0;i<size;i++)
   #pragma acc loop
    for (int j=0;j<size;j++)
      tmp[i*size+j] += Bi[i*size+j] / Ae[j];
    mat_times_mat(A,tmp,B,size);
  }

  #pragma acc exit data delete(B[0:s2],Bi[0:s2],Ae[0:size],tmp[0:s2])

  delete [] B;
  delete [] Ae;
  delete [] Bi;
  delete [] tmp;

  return reset;
}

int invert_stable_cusolver_2(int size, double* A, double delta, cusolverDnHandle_t& cu_hdl)
{
  int s2 = size*size;

  double* B = new double[s2];
  double* Bi = new double[s2];
  double* Ae = new double[size];
  double* tmp = new double[s2]();

  #pragma acc enter data create(B[0:s2],Bi[0:s2],Ae[0:size],tmp[0:s2])

  #pragma acc serial present(B[0:s2],A[0:s2])
   for (int i=0;i<s2;i++)
     B[i] = A[i];

  diagonalize_cusolver(size,size,B,Ae,cu_hdl);

 #pragma acc serial present(Bi[0:s2],B[0:s2])
  trans(Bi,B,size,size);

 #if 1
  #pragma acc update self(Ae[0:size])
  printf(" invert_stable eigenvalues:");
  for (int i=0;i<std::min(10,size);i++)
    printf(" %5.2e",Ae[i]);
  printf("\n");
 #endif

  int nneg = 0;
  int nsmall = 0;
 #pragma acc parallel loop independent present(Ae[0:size]) reduction(+:nsmall,nneg)
  for (int i=0;i<size;i++)
  {
    if (Ae[i]<0.)
    {
      Ae[i] = 1.e10;
      nneg++;
    }
    else if (Ae[i]<delta)
    {
      Ae[i] = delta;
      //Ae[i] = 1.e10;
      nsmall++;
    }
  }

  printf("  found %i small eigen \n",nsmall);

 #pragma acc parallel loop present(tmp[0:s2])
  for (int i=0;i<s2;i++)
    tmp[i] = 0.;

 #pragma acc parallel loop present(tmp[0:s2],Bi[0:s2],Ae[0:size])
  for (int i=0;i<size;i++)
 #pragma acc loop
  for (int j=0;j<size;j++)
    tmp[i*size+j] += Bi[i*size+j] / Ae[j];
  mat_times_mat(A,tmp,B,size);

  #pragma acc exit data delete(B[0:s2],Bi[0:s2],Ae[0:size],tmp[0:s2])

  delete [] B;
  delete [] Ae;
  delete [] Bi;
  delete [] tmp;

  return 0;
}

int invert_stable_cusolver(int size, double* A, double delta, cusolverDnHandle_t& cu_hdl)
{
  int s2 = size*size;

  double* B = new double[s2];
  double* Bi = new double[s2];
  double* Ae = new double[size];
  double* tmp = new double[s2]();

  #pragma acc enter data create(B[0:s2],Bi[0:s2],Ae[0:size],tmp[0:s2])

  #pragma acc serial present(B[0:s2],A[0:s2])
   for (int i=0;i<s2;i++)
     B[i] = A[i];

  diagonalize_cusolver(size,size,B,Ae,cu_hdl);

 #pragma acc serial present(Bi[0:s2],B[0:s2])
  trans(Bi,B,size,size);

 #if 0
  #pragma acc update self(Ae[0:size])
  printf(" invert_stable eigenvalues:");
  for (int i=0;i<size;i++)
    printf(" %5.2e",Ae[i]);
  printf("\n");

  #pragma acc update self(B[0:s2])

  double eig_thresh = read_float("ETHRESH");
  if (eig_thresh<=0.) eig_thresh = 0.01;

  double* vals = new double[size]();
  for (int i=0;i<size;i++)
  if (Ae[i]<-eig_thresh)
  {
    for (int j=0;j<size;j++)
      vals[j] += B[i*size+j]*B[i*size+j];
  }
  printf("  lowev components: ");
  for (int j=0;j<size;j++)
  if (vals[j]>0.03)
    printf(" %2i/%4.2f",j,sqrt(vals[j]));
  printf("\n");

  delete [] vals;
 #endif
 #if 1
  #pragma acc update self(Ae[0:size])
  printf(" invert_stable eigenvalues:");
  for (int i=0;i<std::min(10,size);i++)
    printf(" %5.2e",Ae[i]);
  printf("\n");
 #endif

  int nsmall = 0;
 #pragma acc parallel loop independent present(Ae[0:size]) reduction(+:nsmall)
  for (int i=0;i<size;i++)
  if (Ae[i]<delta)
  {
    Ae[i] = 1.e10;
    nsmall++;
  }

  printf("  found %i small eigen \n",nsmall);

 #pragma acc parallel loop present(tmp[0:s2])
  for (int i=0;i<s2;i++)
    tmp[i] = 0.;

 #pragma acc parallel loop present(tmp[0:s2],Bi[0:s2],Ae[0:size])
  for (int i=0;i<size;i++)
 #pragma acc loop
  for (int j=0;j<size;j++)
    tmp[i*size+j] += Bi[i*size+j] / Ae[j];
  mat_times_mat(A,tmp,B,size);

  #pragma acc exit data delete(B[0:s2],Bi[0:s2],Ae[0:size],tmp[0:s2])

  delete [] B;
  delete [] Ae;
  delete [] Bi;
  delete [] tmp;

  return 0;
}

int invert_cusolver(int N, float* A, cusolverDnHandle_t& cu_hdl)
{
  int reset = 0;

  int N2 = N*N;
  float* B = new float[N2];
  #pragma acc enter data create(B[0:N2])

  #pragma acc serial present(B[0:N2])
  {
    for (int i=0;i<N2;i++) B[i] = 0.;
    for (int i=0;i<N;i++) B[i*N+i] = 1.;
  }

  //int* pivots = new int[N];
  int* info = new int[1];
  //#pragma acc enter data create(pivots[0:N],info[1])
  #pragma acc enter data create(info[1])

  int Lwork = 0;
 #pragma acc host_data use_device(A)
  cusolverDnSgetrf_bufferSize(cu_hdl,N,N,A,N,&Lwork);

  float* work = new float[Lwork];
  #pragma acc enter data create(work[0:Lwork])

  #pragma acc wait
  //printf(" worksize: %2i \n",Lwork);

#if 1
 //check invert eigenstructure
  float* A2 = new float[N2]();
  float* Ae = new float[N]();
  #pragma acc enter data create(A2[0:N2],Ae[0:N])

  #pragma acc parallel loop present(A2[0:N2],A[0:N2])
  for (int i=0;i<N2;i++)
    A2[i] = A[i];

  diagonalize_cusolver(N,A2,Ae,cu_hdl);
  #pragma acc update self(Ae[0:N])

  printf("  debug invert. eigenval:");
  for (int i=0;i<std::min(10,N);i++)
    printf(" %10.8f",Ae[i]);
  printf("\n");

  for (int i=0;i<N;i++)
  if (Ae[i]>100.)
    reset = 1;

  #pragma acc exit data delete(A2[0:N2],Ae[0:N])

  delete [] A2;
  delete [] Ae;
#endif

 #pragma acc host_data use_device(A,B,work,info) //pivots absent
  {
   //solve AX=B with B=I. B overwritten by X.
    cusolverDnSgetrf(cu_hdl,N,N,A,N,work,NULL,info);
    cudaDeviceSynchronize();
    cusolverDnSgetrs(cu_hdl,CUBLAS_OP_N,N,N,A,N,NULL,B,N,info);
    cudaDeviceSynchronize();
  }

  #pragma acc parallel loop present(A[0:N2],B[0:N2])
  for (int i=0;i<N2;i++)
    A[i] = B[i];

  //#pragma acc update self(dLUInfo[1])
  //printf(" after getrf: info_gpu = %d\n", dLUInfo[0]);

  #pragma acc exit data delete(B[0:N2])
  #pragma acc exit data delete(info[0:1],work[0:Lwork])
  //#pragma acc exit data delete(pivots[0:N],info[0:1],work[0:Lwork])


  delete [] B;
  //delete [] pivots;
  delete [] info;
  delete [] work;

  return reset;
}

int invert_cusolver(int N, double* A, cusolverDnHandle_t& cu_hdl)
{
  int reset = 0;

  int N2 = N*N;
  double* B = new double[N2];
  #pragma acc enter data create(B[0:N2])

  #pragma acc serial present(B[0:N2])
  {
    for (int i=0;i<N2;i++) B[i] = 0.;
    for (int i=0;i<N;i++) B[i*N+i] = 1.;
  }

  int* info;
  cudaMalloc((void**) &info,sizeof(int));
  //int* pivots = new int[N];
  //int* info = new int[1];
  //#pragma acc enter data create(pivots[0:N],info[1])
  //#pragma acc enter data create(info[1])

#if 1
 //check invert eigenstructure
  double* A2 = new double[N2]();
  double* Ae = new double[N]();
  #pragma acc enter data create(A2[0:N2],Ae[0:N])

  #pragma acc parallel loop present(A2[0:N2],A[0:N2])
  for (int i=0;i<N2;i++)
    A2[i] = A[i];

  diagonalize_cusolver(N,N,A2,Ae,cu_hdl);
  #pragma acc update self(Ae[0:N])
  //#pragma acc update self(A2[0:N2],Ae[0:N])

  for (int i=0;i<N;i++)
  if (Ae[i]>100.)
    reset = 1;

  printf("\n debug invert. eigenval:");
  for (int i=0;i<std::min(10,N);i++)
    printf(" %10.8f",Ae[i]);
  printf("\n");

  //print_square_sm(N,A2);

  //printf(" lowest eigenvec: ");
  //for (int i=0;i<N;i++)
  //  printf(" %5.2f",A2[i]);
  //printf("\n");

  #pragma acc exit data delete(A2[0:N2],Ae[0:N])

  delete [] A2;
  delete [] Ae;
#endif

  int Lwork = 0;
 #pragma acc host_data use_device(A)
  cusolverDnDgetrf_bufferSize(cu_hdl,N,N,A,N,&Lwork);

  double* work;
  cudaMalloc((void**) &work,Lwork*sizeof(double));
  //#pragma acc enter data create(work[0:Lwork])

  #pragma acc wait
  //printf(" worksize: %2i \n",Lwork);

 #pragma acc host_data use_device(A,B)
  {
   //solve AX=B with B=I. B overwritten by X.
    cusolverDnDgetrf(cu_hdl,N,N,A,N,work,NULL,info);
    cudaDeviceSynchronize();
    cusolverDnDgetrs(cu_hdl,CUBLAS_OP_N,N,N,A,N,NULL,B,N,info);
    cudaDeviceSynchronize();
  }

  #pragma acc parallel loop present(A[0:N2],B[0:N2])
  for (int i=0;i<N2;i++)
    A[i] = B[i];

  //#pragma acc update self(dLUInfo[1])
  //printf(" after getrf: info_gpu = %d\n", dLUInfo[0]);

  #pragma acc exit data delete(B[0:N2])
  //#pragma acc exit data delete(work[0:Lwork])
  //#pragma acc exit data delete(pivots[0:N],info[0:1],work[0:Lwork])


  delete [] B;
  //delete [] pivots;
  //delete [] info;
  //delete [] work;
  cudaFree(info);
  cudaFree(work);

  return reset;
}

