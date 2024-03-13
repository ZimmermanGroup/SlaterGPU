#include "utils.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <chrono>

extern "C" void dgemm_(char* TA, char* TB, int* M, int* N, int* K, 
 double* ALPHA, double* A, int* LDA, 
 double* B, int* LDB, double* BETA, double* C, int* LDC);

extern "C" void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
                      double *w, double * work, int *lwork, int* info);

int get_NOs(double thresh, int N, double* rho1, double* rho1e, double* jCA, double* jCA1, int prl)
{
  int N2 = N*N;
  double rot[N2];
  double rho2[N2];
  for (int i=0;i<N2;i++)
    rho2[i] = rho1[i];

  DiagonalizeP(rho2,rho1e,N);
  if (prl>1)
  {
    printf("\n NO eigenvalues:");
    for (int i=N-1;i>=0;i--)
      printf(" %12.10f",rho1e[i]);
    printf("\n");
  }

  int nsig = 0;
  for (int i=N-1;i>=0;i--)
  if (fabs(rho1e[i])>thresh && rho1e[i]!=-1.)
    nsig++;

  for (int i=0;i<N;i++)
  {
    for (int j=0;j<N;j++)
      rot[(N-1-i)*N+j] = rho2[i*N+j];
  }
  double tmp[N];
  for (int i=0;i<N;i++)
    tmp[i] = rho1e[N-1-i];
  for (int i=0;i<N;i++)
    rho1e[i] = tmp[i];

  mat_times_mat_bt_cpu(jCA1,jCA,rot,N);

  //printf("\n jCA: \n"); print_square(N,jCA);
  //printf("\n jCA1: \n"); print_square(N,jCA1);

  return nsig;
}

void ao_to_mo_cpu(int N, double* jAO, double* jMO, double* jCA, double* tmp)
{
  mat_times_mat_cpu(tmp,jAO,jCA,N);
  mat_times_mat_at_cpu(jMO,jCA,tmp,N);
   
  return;
}

void mo_to_ao_cpu(int N, double* Pmo, double* Pao, double* jCA, double* tmp)
{
  mat_times_mat_bt_cpu(tmp,Pmo,jCA,N);
  mat_times_mat_cpu(Pao,jCA,tmp,N);

  return;
}

void mat_times_mat_cpu(double* C, double* A, double* B, int N)
{
  mat_times_mat_cpu(C,A,B,N,N,N);

  return;
}

void mat_times_mat_cpu(double* C, double* A, double* B, int M, int N, int K)
{
 #if 0
  for (int i=0;i<M*N;i++) C[i] = 0.;

  for (int i=0;i<M;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<K;k++)
    C[i*N+j] += A[i*K+k]*B[k*N+j];
 #endif

 #if 1
 //probably correct
  int LDA = K;
  int LDB = N;
  int LDC = N;
  double ALPHA = 1.0;
  double BETA = 0.0;
  char TA = 'N';
  char TB = 'N';
  dgemm_(&TB,&TA,&N,&M,&K,&ALPHA,B,&LDB,A,&LDA,&BETA,C,&LDC);
 #endif

  return;
}

void mat_times_mat_at_cpu(double* C, double* A, double* B, int N)
{
  int M = N; int K = N;
  int LDA = N;
  int LDB = N;
  int LDC = N;
   
  double ALPHA = 1.0;
  double BETA = 0.0;
  
 //C := alpha*op( A )*op( B ) + beta*C (op means A or B, possibly transposed)
 //CBlas version
  //cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC);
  
 //LIBLAS version
 // note: funny order of A/B to get C ordered
  char TB = 'N';
  char TA = 'T';
  dgemm_(&TB,&TA,&M,&N,&K,&ALPHA,B,&LDB,A,&LDA,&BETA,C,&LDC);

 #if 0
  printf("\n C (gemm): \n");
  for (int m=0;m<N;m++)
  {
    for (int n=0;n<N;n++)
      printf(" %8.5f",C[m*N+n]);
    printf("\n");
  }

  for (int i=0;i<N*N;i++) C[i] = 0.;

  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    C[i*N+j] += A[k*N+i]*B[k*N+j];

  printf("\n C (loops): \n");
  for (int m=0;m<N;m++)
  {
    for (int n=0;n<N;n++)
      printf(" %8.5f",C[m*N+n]);
    printf("\n");
  }
 #endif

  return;
}

void mat_times_mat_bt_cpu(double* C, double* A, double* B, int M, int N, int K, int LDAB)
{
 #if 0
  for (int i=0;i<M*N;i++) C[i] = 0.;

  for (int i=0;i<M;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<K;k++)
    C[i*N+j] += A[i*K+k]*B[j*K+k];

  printf("\n C (loops): \n");
  for (int m=0;m<M;m++)
  {
    for (int n=0;n<N;n++)
      printf(" %8.5f",C[m*N+n]);
    printf("\n");
  }
 #endif

  int LDA = LDAB; //K
  int LDB = LDAB; //K
  int LDC = N;
   
  double ALPHA = 1.0;
  double BETA = 0.0;
  
 //C := alpha*op( A )*op( B ) + beta*C (op means A or B, possibly transposed)
 //CBlas version
  //cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC);
  
 //LIBLAS version
 // note: funny order of A/B to get transpose right
  char TB = 'T';
  char TA = 'N';
  dgemm_(&TB,&TA,&N,&M,&K,&ALPHA,B,&LDB,A,&LDA,&BETA,C,&LDC);

 #if 0
  printf("\n C (gemm): \n");
  for (int m=0;m<M;m++)
  {
    for (int n=0;n<N;n++)
      printf(" %8.5f",C[m*N+n]);
    printf("\n");
  }
 #endif

  return;
}

void mat_times_mat_bt_cpu(double* C, double* A, double* B, int M, int N, int K)
{
  for (int i=0;i<M*N;i++) C[i] = 0.;

  for (int i=0;i<M;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<K;k++)
    C[i*N+j] += A[i*N+k]*B[j*N+k];

  return;
}

void mat_times_mat_bt_cpu(double* C, double* A, double* B, int N)
{
  //return mat_times_mat_bt_cpu(C,A,B,N,N,N);

  for (int i=0;i<N*N;i++) C[i] = 0.;

  int M = N; int K = N;
  int LDA = N;
  int LDB = N;
  int LDC = N;
   
  double ALPHA = 1.0;
  double BETA = 0.0;
  
 //C := alpha*op( A )*op( B ) + beta*C (op means A or B, possibly transposed)
 //CBlas version
  //cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC);
  
 //LIBLAS version
 // note: funny order of A/B to get transpose right
  char TB = 'T';
  char TA = 'N';
  dgemm_(&TB,&TA,&M,&N,&K,&ALPHA,B,&LDB,A,&LDA,&BETA,C,&LDC);

 #if 0
  printf("\n C (gemm): \n");
  for (int m=0;m<N;m++)
  {
    for (int n=0;n<N;n++)
      printf(" %8.5f",C[m*N+n]);
    printf("\n");
  }

  mat_times_mat_bt_cpu(C,A,B,N,N,N);
  printf("\n C (loops): \n");
  for (int m=0;m<N;m++)
  {
    for (int n=0;n<N;n++)
      printf(" %8.5f",C[m*N+n]);
    printf("\n");
  }
 #endif

  return;
}

void DiagonalizeP(double* A, double* Ae, int N)
{
  int N2 = N*N;
  char v = 'V';
  char u = 'L';
  int info;
  double * work = new double[N2];
  int worksize = max(128,N2);
  dsyev_(&v,&u,&N,A,&N,Ae,work,&worksize,&info);
  delete [] work;
}

void print_duration(chrono::high_resolution_clock::time_point t1, chrono::high_resolution_clock::time_point t2, string name)
{
  //return;

  chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
  if (time_span.count()<1.)
    printf(" %s took %5.1f ms \n",name.c_str(),time_span.count()*1000.);
  else
    printf(" %s took %6.3f s \n",name.c_str(),time_span.count());
}

void ao_to_mo_4c(double ** eri, double ** jMOI, double * jCA, int N)
{
  int N2 = N * N;
  double ** tmp_eri1 = new double*[N2];
  double ** tmp_eri2 = new double*[N2];
  double ** tmp_eri3 = new double*[N2];
  for (int i = 0 ; i < N2; i++) {
    tmp_eri1[i] = new double[N2]();
    tmp_eri2[i] = new double[N2]();
    tmp_eri3[i] = new double[N2]();
  }
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N2 * N; j++) {
      jMOI[i][j] = 0;
    }
  }
  #pragma omp parallel for
  for (int p = 0; p < N; p++) {
    for (int mu = 0; mu < N; mu++) {
      for (int q = 0; q < N; q++) {
        for (int r = 0; r < N; r++) {
          for (int s = 0; s < N; s++) {
            tmp_eri1[p*N+q][r*N+s] += jCA[mu*N+p] * eri[mu*N+q][r*N+s];
          }
        }
      }
    } // for mu
    for (int q = 0; q < N; q++) {
      for (int nu = 0; nu < N; nu++) {
        for (int r = 0; r < N; r++) {
          for (int s = 0; s < N; s++) {
            tmp_eri2[p*N+q][r*N+s] += jCA[nu*N+q] * tmp_eri1[p*N+nu][r*N+s];
          }
        }
      } // for nu
      for (int r = 0; r < N; r++) {
        for (int lam = 0; lam < N; lam++) {
          for (int s = 0; s < N; s++) {
            tmp_eri3[p*N+q][r*N+s] += jCA[lam*N+r] * tmp_eri2[p*N+q][lam*N+s];
          }
        } // for lam
        for (int s = 0; s < N; s++) {
          for (int sig = 0; sig < N; sig++) {
            jMOI[p][q*N2+r*N+s] += jCA[sig*N+s] * tmp_eri3[p*N+q][r*N+sig];
          } // for sig
        } // for s
      } // for r
    } // for q
  } // for i
  for (int i = 0 ; i < N2; i++) {
    delete [] tmp_eri1[i];
    delete [] tmp_eri2[i];
    delete [] tmp_eri3[i];
  }
  delete [] tmp_eri1;
  delete [] tmp_eri2;
  delete [] tmp_eri3;
}
