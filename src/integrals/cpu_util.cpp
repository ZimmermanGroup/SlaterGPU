#include "cpu_util.h"


///need to convert mat_times functions to gemm
// (for general M,N,K)


extern "C" {

extern void zgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int* K, complex<double>* ALPHA, complex<double>* A, int* LDA,
 complex<double>* B, int* LDB, complex<double>* BETA, complex<double>* C, int* LDC);

extern void zheev_(char* JOBZ, char* UPLO, int* N, complex<double>* A, int* LDA,
 double* W, complex<double>* WORK, int* LWORK, double* RWORK, int* INFO);

extern void dgeev_(char* JOBVL, char* JOBVR, int* N, double* A, int* LDA,
 double* WR, double* WI, double* VL, int* LDVL, double* VR, int* LDVR,
double* WORK, int* LWORK, int* INFO);

extern void dgemm_(char* TA, char* TB, int* M, int* N, int* K,
 double* ALPHA, double* A, int* LDA,
 double* B, int* LDB, double* BETA, double* C, int* LDC);

extern void dsyevx_(
    char * jobz, char * range, char * uplo,
    int * n,
    double* A, int * lda,
    double * vl,
    double * vu, int * il, int * iu,
    double * abstol, int* m,
    double* W,
    double* Z, int * ldz,
    double* work, int * lwork,
    int* iwork, int* IFAIL,
    int* info );

extern void dgetrf_(int* M, int* N, double* A, int* LDA, int* IPIV, int* INFO);

extern void dgetri_(int* N, double* A, int* LDA, int* IPIV, double* WORK, int* LWORK, int* INFO);
}

extern "C" {
    void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv,
                double* b, int* ldb, int* info);
}

void print_square(int N, double* A);
void print_square_sm(int N, double* A);

double randomf(double a, double b)
{
  timeval t1;
  gettimeofday(&t1, NULL);
  srand(t1.tv_usec*t1.tv_sec);

  double range = b-a;
  double randn = a + double(range*rand()/(RAND_MAX));
  return randn;
}

void solve_axeb(int dim, double* A, double* b)
{
  int info;
  int one = 1;
  int z[dim];
  dgesv_(&dim,&one,A,&dim,z,b,&dim,&info);

  return;
}

void expmat_complex_cpu(int N, double* theta, double* thetai, double* etheta)
{
  printf("\n ERROR expmat_complex_cpu broken by new compiler \n"); exit(-1);

#if 0
 //exp(-i*theta)
  int N2 = N*N;
  char JOBZ = 'V';
  char UPLO = 'U';
  int LWORK = N2;
  int INFO = 0;

  double evals[N];
  complex<double> itheta[N2];
  complex<double> WORK[LWORK];
  double RWORK[3*N-2];

  complex<double> tmp1[N2];
  complex<double> tmp2[N2];
  complex<double> tmp3[N2];

  for (int j = 0; j < N2; j++)
    itheta[j] = thetai[j] - 1.i * theta[j];

  zheev_(&JOBZ,&UPLO,&N,itheta,&N,evals,WORK,&LWORK,RWORK,&INFO);

  for (int i = 0; i < N2; i++)
    tmp3[i] = 0.;

  for (int j = 0; j < N; j++)
    tmp3[j*N+j] = cos(evals[j]) + 1.i*sin(evals[j]);

  for (int i = 0; i < N; i++)
  for (int j = 0; j < N; j++)
    tmp1[i*N+j] = conj(itheta[j*N+i]);

  char TRANSA = 'N';
  char TRANSB = 'N';

  complex<double> one = 1.0 + 0.i;
  complex<double> zero = 0. + 0.i;
  zgemm_(&TRANSA,&TRANSB,&N,&N,&N,&one,itheta,&N,tmp3,&N,&zero,tmp2,&N);
  zgemm_(&TRANSA,&TRANSB,&N,&N,&N,&one,tmp2,  &N,tmp1,&N,&zero,tmp3,&N);

 //CPMZ need to update this
  for (int i=0;i<N2;i++)
    etheta[i] = tmp3[i].real();

  printf("\n etheta(real): \n");
  print_square(N,etheta);

 //CPMZ need to update this
  for (int i=0;i<N2;i++)
    etheta[i] = tmp3[i].imag();

  printf("\n etheta(imag): \n");
  print_square(N,etheta);

  printf("\n\n ending early \n");
  fflush(stdout);
  exit(1);

 #if 0
  printf("Testing if matrix is unitary\n");
  double tmp4[N2];

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
              N, N, N, 1., etheta, N, etheta, N, 0., tmp4, N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      printf("%7.4f ", tmp4[i*N+j]);
    }
    printf("\n");
  }
 #endif
#endif

  return;
}

void expmat_cpu(int N, double *theta, double *etheta)
{
 //computes exp(-i theta)

  printf("\n WARNING: expmat_cpu broken by new compiler \n");

 #if 1
  int N2 = N*N;
  char JOBZ = 'V';
  char UPLO = 'U';
  int LWORK = N2;
  int INFO = 0;

  double evals[N];
  complex<double> itheta[N2];
  complex<double> WORK[LWORK];
  double RWORK[3*N-2];

  complex<double> tmp1[N2];
  complex<double> tmp2[N2];
  complex<double> tmp3[N2];

  for (int i = 0; i < N2; i++)
    itheta[i] = 0. - 1.i * theta[i];

  zheev_(&JOBZ,&UPLO,&N,itheta,&N,evals,WORK,&LWORK,RWORK,&INFO);

  for (int i = 0; i < N2; i++)
    tmp3[i] = 0.;

  for (int i = 0; i < N; i++)
    tmp3[i*N+i] = cos(evals[i]) + 1.i*sin(evals[i]);

  for (int i = 0; i < N; i++)
  for (int j = 0; j < N; j++)
    tmp1[i*N+j] = conj(itheta[j*N+i]);

  char TRANSA = 'N';
  char TRANSB = 'N';

  complex<double> one = 1.0 + 0.i;
  complex<double> zero = 0. + 0.i;
  zgemm_(&TRANSA,&TRANSB,&N,&N,&N,&one,itheta,&N,tmp3,&N,&zero,tmp2,&N);
  zgemm_(&TRANSA,&TRANSB,&N,&N,&N,&one,tmp2,  &N,tmp1,&N,&zero,tmp3,&N);

  for (int i=0;i<N2;i++)
    etheta[i] = tmp3[i].real();

  //printf("\n etheta: \n");
  //print_square(N,etheta);

  #if 0
  printf("Testing if matrix is unitary\n");
  double tmp4[N2];

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
              N, N, N, 1., etheta, N, etheta, N, 0., tmp4, N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      printf("%7.4f ", tmp4[i*N+j]);
    }
    printf("\n");
  }
  #endif

 #endif

  return;
}

int la_diagR(int size, double* A, double* eigen, double* eigeni)
{
  printf(" in la_diagR call, size: %i \n",size);

  int N = size;
  int LDA = size;

  //double AbsTol = 0.0;

  int N2 = N*N;
  double EVecR[N2];
  //double eigeni[N];
  for (int i=0;i<size;i++) eigen[i] = eigeni[i] = 0.;

  // Give dgeev extra work space, 4*N is minimum
  int LenWork = 32*N; //8*N min for dsyevx
  double Work[LenWork];

  int Info = 0;

  char JobVL = 'N';
  char JobVR = 'V';

  dgeev_(&JobVL,&JobVR,&N,A,&LDA,eigen,eigeni,NULL,&N,EVecR,&N,Work,&LenWork,&Info);

  if (Info != 0) {
    printf(" Info = %d\n",Info);
    printf("Call to dgeev failed in Diagonalize");
  }

  for (int i=0;i<size;i++)
  for (int j=0;j<size;j++)
    A[i*size+j] = EVecR[i*size+j];

 #if 1
  printf(" real part of eigenvalues \n");
  for (int i=0;i<size;i++)
    printf(" %6.3f",eigen[i]);
  printf("\n");
  printf(" imag part of eigenvalues \n");
  for (int i=0;i<size;i++)
    printf(" %6.3f",eigeni[i]);
  printf("\n");
 #endif

 //copy back just imaginary part
  //for (int i=0;i<size;i++)
  //  eigen[i] = eigeni[i];

  return 0;
}

void la_diag(int neig, int s1, double* A, double* Ae)
{
  #define USE_DSYEVX 1

  int s1a = s1;
  int LDA = s1;

  char jobz = 'V';
  char uplo = 'u';
 #if USE_DSYEVX
  char range = 'I';
  if (neig==s1)
    range = 'A';

  int il = 1;
  int iu = neig;

  double abstol = 0.;
  double vl = 1.; double vu = -1.;
  int neval = s1;
 #endif

  double* EVec = new double[LDA*s1];

 #if USE_DSYEVX
  int LenWork = 32*s1; //8*N min for dsyevx
 #else
  int LenWork = 1+6*s1+2*s1*s1; //1+6*N+2*N*N min for dsyevd
 #endif
  double* Work = new double[LenWork];

  int LenIWork = 10*s1;
  int* IWork = new int[LenIWork];
  int* IFail = new int[s1];

  int Info = 0;

 #if USE_DSYEVX
  dsyevx_(&jobz,&range,&uplo,&s1a,A,&LDA,&vl,&vu,&il,&iu,&abstol,
       &neval,Ae,EVec,&LDA,Work,&LenWork,IWork,IFail,&Info);
 #else
  dsyevd_(&jobz,&uplo,&s1a,A,&LDA,Ae, Work,&LenWork,IWork,&LenIWork, &Info);
 #endif

  #if 1
  if (Info != 0)
    printf(" Info = %d\n",Info);
  #endif

 #if USE_DSYEVX
  for (int i=0;i<s1;i++)
  for (int j=0;j<s1;j++)
    A[i*s1+j] = EVec[i*s1+j];
 #endif

 #if 0
  printf("\n printing eigenvectors \n");
  for (int i=0;i<s1;i++)
  {
    for (int j=0;j<s1;j++)
      printf(" %8.5f",A[i*s1+j]);
    printf("\n");
  }
  printf("\n");
 #endif
 #if 0
  printf(" eigenvalues: ");
  for (int i=0;i<min(4,s1);i++)
    printf(" %8.5f",Ae[i]);
  printf("\n");
 #endif

  delete [] EVec;
  delete [] Work;
  delete [] IWork;
  delete [] IFail;

}

void test_diag()
{
  int N = 4;
  int N2 = N*N;

  double* A = new double[N2];
  double* Ae = new double[N];
  for (int i=0;i<N2;i++)
    A[i] = 0.1;
  for (int i=0;i<N;i++)
    A[i*N+i] = 1.;

  la_diag(N,N,A,Ae);


  return;
}


int invert_stable_cpu_3(double* A, int size, double delta, bool root)
{
  if (A==NULL) return -1;
  //if (delta>0.) printf(" using Invert with %12.10f shift factor \n",delta);

  double* B = new double[size*size];
  for (int i=0;i<size*size;i++) B[i] = A[i];
  double* Beigen = new double[size];
  for (int i=0;i<size;i++) Beigen[i] = 0.;

  for (int i=0;i<size;i++)
    B[i*size+i] += delta;

  la_diag(size,size,B,Beigen);

  double* Bi = new double[size*size];

  trans_cpu(Bi,B,size,size);

  double* tmp = new double[size*size];
  for (int i=0;i<size*size;i++) tmp[i] = 0.;
  for (int i=0;i<size*size;i++) A[i] = 0.;

  int nf = 0;
  double mineig = 1000.;
  for (int i=0;i<size;i++)
  {
    if (Beigen[i]<mineig) mineig = Beigen[i];
    if (Beigen[i]<0.)
    {
      Beigen[i] = 1.e10;
      nf++;
    }
    else if (root)
      Beigen[i] = sqrt(Beigen[i]);
  }
  //printf(" found %2i small eigenvalues. lowest: %4.2g \n",nf,mineig);


#if 0
  printf(" B(eigen):");
  for (int i=0;i<size;i++)
    printf(" %6.4f",Beigen[i]);
  printf("\n");
#endif

  for (int i=0;i<size;i++)
  for (int j=0;j<size;j++)
    tmp[i*size+j] += Bi[i*size+j] / Beigen[j];
 //need MM

  mat_times_mat_cpu(A,tmp,B,size);
  //for (int i=0;i<size;i++)
  //for (int j=0;j<size;j++)
  //for (int k=0;k<size;k++)
  //  A[i*size+j] += tmp[i*size+k] * B[k*size+j];


  delete [] B;
  delete [] Beigen;
  delete [] Bi;
  delete [] tmp;

  return nf;
}


int invert_stable_cpu_2(double* A, int size, double delta)
{
  if (A==NULL) return -1;
  //if (delta>0.) printf(" using Invert with %12.10f minimum eigenvalue \n",delta);

  double* B = new double[size*size];
  for (int i=0;i<size*size;i++) B[i] = A[i];
  double* Beigen = new double[size];
  for (int i=0;i<size;i++) Beigen[i] = 0.;

  la_diag(size,size,B,Beigen);

  double* Bi = new double[size*size];
  bool* zeroes = new bool[size]();

  trans_cpu(Bi,B,size,size);

  double* tmp = new double[size*size];
  for (int i=0;i<size*size;i++) tmp[i] = 0.;
  for (int i=0;i<size*size;i++) A[i] = 0.;

  int nf = 0;
  double mineig = 1000.;
  for (int i=0;i<size;i++)
  {
    if (Beigen[i]<=0.)
    {
      zeroes[i] = 1;
      Beigen[i] = 1.e10;
      nf++;
    }
    else if (Beigen[i]<delta)
    {
      //printf("  screening: %14.12f \n",Beigen[i]);
      //Beigen[i] = 1.e10;
      Beigen[i] = delta;
      nf++;
    }
    if (Beigen[i]<mineig) mineig = Beigen[i];
  }
  //printf(" found %2i small eigenvalues. lowest: %4.2g \n",nf,mineig);

#if 0
  printf(" B(eigen):");
  for (int i=0;i<size;i++)
    printf(" %6.4f",Beigen[i]);
  printf("\n");
#endif

  for (int j=0;j<size;j++)
  if (!zeroes[j])
  for (int i=0;i<size;i++)
    tmp[i*size+j] += Bi[i*size+j] / Beigen[j];
 //need MM
  for (int i=0;i<size;i++)
  for (int j=0;j<size;j++)
  for (int k=0;k<size;k++)
    A[i*size+j] += tmp[i*size+k] * B[k*size+j];


  delete [] B;
  delete [] Beigen;
  delete [] Bi;
  delete [] tmp;
  delete [] zeroes;

  return nf;
}

int invert_stable_cpu(double* A, int size, double delta)
{
  return invert_stable_cpu(A,size,delta,0);
}

int invert_stable_cpu(double* A, int size, double delta, int prl)
{
  if (A==NULL) return -1;
  //double delta = INV_CUTOFF;
  if (delta>0. && prl) printf("   invert w/ %4.1e minimum eigenvalue \n",delta);

  double* B = new double[size*size];
  for (int i=0;i<size*size;i++) B[i] = A[i];
  double* Beigen = new double[size];
  for (int i=0;i<size;i++) Beigen[i] = 0.;

  la_diag(size,size,B,Beigen);

  double* Bi = new double[size*size];

  trans_cpu(Bi,B,size,size);

  double* tmp = new double[size*size];
  for (int i=0;i<size*size;i++) tmp[i] = 0.;
  for (int i=0;i<size*size;i++) A[i] = 0.;

 #if 0
  printf(" printing Bi \n");
  print_square_ss_sm(size,Bi);
 #endif

  int nf = 0;
  double mineig = 1000.;
  for (int i=0;i<size;i++)
  {
    if (Beigen[i]<delta)
    {
      if (prl>1) printf("    screening: %14.12f \n",Beigen[i]);
      Beigen[i] = 1.e10;
      //Beigen[i] = INV_MIN;
     // for (int j=0;j<size;j++)
     //   Bi[j*size+i] = B[i*size+j] = 0.;
      nf++;
    }
    if (Beigen[i]<mineig) mineig = Beigen[i];
  }
  if (prl>0) printf("   found %2i small eigenvalues. lowest: %4.2g \n",nf,mineig);

  if (prl>1)
  {
    printf(" B(eigen):");
    for (int i=0;i<size;i++)
      printf(" %5.2e",Beigen[i]);
    printf("\n");
  }

  for (int i=0;i<size;i++)
  for (int j=0;j<size;j++)
    tmp[i*size+j] += Bi[i*size+j] / Beigen[j];
 //need MM
  for (int i=0;i<size;i++)
  for (int j=0;j<size;j++)
  for (int k=0;k<size;k++)
    A[i*size+j] += tmp[i*size+k] * B[k*size+j];


  delete [] B;
  delete [] Beigen;
  delete [] Bi;
  delete [] tmp;

  return nf;
}

int invert_cpu(double* A, int size, int mode)
{
  printf("\n Invert not yet available \n");
 #if 0
  if (size<1)
  {
    printf("  WARNING: cannot invert, size: %i \n",size);
    return 1;
  }

  int print_inverse = 0;
  if (mode>0)
  {
    double icut = INV_CUTOFF;
    printf(" size(A): %3i \n",size);
    double* A2 = new double[size*size];
    for (int i=0;i<size*size;i++) A2[i] = A[i];
    double* A2e = new double[size]();
    la_diag(size,size,A2,A2e);
    //printf("  first/last eigenvalue of A:  %8g  %8g \n",A2e[0],A2e[size-1]);
    printf("  eigenvalues of A:");
    for (int i=0;i<size;i++) printf(" %2g",A2e[i]);
    printf("\n");
    int nf = 0;
    int nzf = 0;
    for (int i=0;i<size;i++)
    {
      if (std::isnan(A2e[i]))
        nf++;
      if (A2e[i]<icut)
        nzf++;
    }
    if (nf>0)
    {
      printf("  found %2i NaN's, exiting \n",nf);
      exit(1);
    }
    delete [] A2;
    delete [] A2e;

    if (nzf>2) print_inverse = 1;

    if (nzf>0)
      return Invert_stable(A,size,icut);
  }

  int LenWork = 4*size;
  double* Work = new double[LenWork];

  int Info = 0;

  //printf(" LenWork: %i \n",LenWork);

  int* IPiv = new int[size];

  dgetrf_(&size,&size,A,&size,IPiv,&Info);
  if (Info!=0)
  {
    printf("\n  after dgetrf, Info error is: %i \n",Info);
    delete [] IPiv;
    delete [] Work;

    for (int i=0;i<size*size;i++) A[i] = 0.;
    for (int i=0;i<size;i++)
       A[i*size+i] = 1.;

    exit(1);
    return 1;
  }

  dgetri_(&size,A,&size,IPiv,Work,&LenWork,&Info);
  if (Info!=0 || print_inverse)
  {
    printf(" after invert, Info error is: %i \n",Info);
    printf(" A-1: \n");
    for (int i=0;i<size;i++)
    {
      for (int j=0;j<size;j++)
        printf(" %6.3f",A[i*size+j]);
      printf("\n");
    }
  }

  delete [] IPiv;
  delete [] Work;
 #endif

  return 0;
}

int mat_root_cpu(double* A, int size)
{
  double* B = new double[size*size];
  for (int i=0;i<size*size;i++) B[i] = A[i];
  double* Beigen = new double[size];
  for (int i=0;i<size;i++) Beigen[i] = 0.;

  la_diag(size,size,B,Beigen);

  double* Bi = new double[size*size];
  //for (int i=0;i<size*size;i++) Bi[i] = B[i];

  trans_cpu(Bi,B,size,size);

  double* tmp = new double[size*size];
  for (int i=0;i<size*size;i++) tmp[i] = 0.;
  for (int i=0;i<size*size;i++) A[i] = 0.;

  for (int i=0;i<size;i++)
  for (int j=0;j<size;j++)
    tmp[i*size+j] += Bi[i*size+j] * sqrt(Beigen[j]);
  for (int i=0;i<size;i++)
  for (int j=0;j<size;j++)
  for (int k=0;k<size;k++)
    A[i*size+j] += tmp[i*size+k] * B[k*size+j];


  delete [] B;
  delete [] Beigen;
  delete [] Bi;
  delete [] tmp;

  return 0;
}

int mat_root_inv_cpu(double* A, int size)
{
  int size2 = size*size;

  double* B = new double[size2];
  for (int i=0;i<size2;i++) B[i] = A[i];
  double* Beigen = new double[size];
  for (int i=0;i<size;i++) Beigen[i] = 0.;

  la_diag(size,size,B,Beigen);

  double* Bi = new double[size2];
  //for (int i=0;i<size2;i++) Bi[i] = B[i];

  trans_cpu(Bi,B,size,size);

  double* tmp = new double[size2]();
  for (int i=0;i<size2;i++) A[i] = 0.;

 #if 0
  printf("  eigenvalues for mat_root_inv:");
  for (int i=0;i<min(size,8);i++) printf(" %5g",Beigen[i]);
  printf("\n");
 #endif

  int nsmall = 0;
  for (int i=0;i<size;i++)
  if (Beigen[i]<1.e-10)
    nsmall++;

  for (int i=0;i<size;i++)
  for (int j=0;j<size;j++)
    tmp[i*size+j] += Bi[i*size+j] / sqrt(Beigen[j]);
  for (int i=0;i<size;i++)
  for (int j=0;j<size;j++)
  for (int k=0;k<size;k++)
    A[i*size+j] += tmp[i*size+k] * B[k*size+j];


  delete [] B;
  delete [] Beigen;
  delete [] Bi;
  delete [] tmp;

  return nsmall;
}

int mat_root_inv_stable_cpu(double* A, int size, double inv_cutoff, int prl)
{
  int size2 = size*size;

  double* B = new double[size2];
  for (int i=0;i<size2;i++) B[i] = A[i];
  double* Beigen = new double[size];
  for (int i=0;i<size;i++) Beigen[i] = 0.;

  la_diag(size,size,B,Beigen);

  double* Bi = new double[size2];
  //for (int j=0;j<size2;j++) Bi[j] = B[j];
  trans_cpu(Bi,B,size,size);

  if (prl>1)
  {
    printf("  eigenvalues for mat_root_inv (cutoff: %10.8f):",inv_cutoff);
    for (int i=0;i<min(10,size);i++) printf(" %5g",Beigen[i]);
    printf("\n");
  }

  int nlow = 0;
  for (int i=0;i<size;i++)
  if (Beigen[i]<inv_cutoff)
  {
    Beigen[i] = 1.e20;
    nlow++;
  }
  if (prl>0)
    printf("  found %i low eigenvalues \n",nlow);

  double* tmp = new double[size2]();
  for (int i=0;i<size;i++)
  for (int j=0;j<size;j++)
    tmp[i*size+j] = Bi[i*size+j] / sqrt(Beigen[j]);

  for (int i=0;i<size2;i++) A[i] = 0.;
  for (int i=0;i<size;i++)
  for (int j=0;j<size;j++)
  for (int k=0;k<size;k++)
    A[i*size+j] += tmp[i*size+k] * B[k*size+j];


  delete [] B;
  delete [] Beigen;
  delete [] Bi;
  delete [] tmp;

  return nlow;
}

int LU_inv_stable_cpu(double* A, int size)
{
 //Need to check if debug is necessary, but need ZEST incorporation first
    int size2 = size*size;
    int infolu;
    int infoinv;
    int* piv = new int[size];
    double* work = new double[size2];

    dgetrf_(&size,&size,A,&size,piv,&infolu);
    if(infolu != 0)
    {    
      printf("LU factorization failed %i\n",infolu);
      delete [] piv; 
      delete [] work;
      exit(1);
    }    

    dgetri_(&size,A,&size,piv,work,&size2,&infoinv);
    if(infoinv != 0)
    {    
      printf("LU inversion failed %i\n",infoinv);
      exit(1);
    }    

    delete [] piv; 
    delete [] work;

    return infoinv;
}

void trans_cpu(float* Bt, float* B, int m, int n)
{
  for (int i=0;i<m;i++)
  for (int j=0;j<n;j++)
    Bt[i*n+j] = B[j*m+i];

  return;
}

void trans_cpu(double* Bt, double* B, int m, int n)
{
  for (int i=0;i<m;i++)
  for (int j=0;j<n;j++)
    Bt[i*n+j] = B[j*m+i];

  return;
}

void mat_times_mat_cpu(float* C, float* A, float* B, int M, int N, int K)
{
  for (int i=0;i<M*N;i++) C[i] = 0.;

  for (int i=0;i<M;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<K;k++)
    C[i*N+j] += A[i*K+k]*B[k*N+j];

  return;
}

void mat_times_mat_cpu(float* C, float* A, float* B, int N)
{
  mat_times_mat_cpu(C,A,B,N,N,N);

  return;
}

void mat_times_mat_ct_cpu(double* C, double* A, double* B, int M, int N, int K)
{
  int LDA = K; //M by K
  int LDB = N; //K by N
  int LDC = M; //M by N

  double ALPHA = 1.0;
  double BETA = 0.0;

 //LIBLAS version
 //need to transpose to get major right
  char TA = 'T';
  char TB = 'T';
  dgemm_(&TA,&TB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);

 // note: funny order of A/B to get C ordered (no _ct)
  //LDA = K; LDB = N; LDC = N;
  //char TA = 'N';
  //char TB = 'N';
  //dgemm_(&TB,&TA,&N,&M,&K,&ALPHA,B,&LDB,A,&LDA,&BETA,C,&LDC);

 #if 0
  printf("\n C (gemm): \n");
 //transposed
  for (int m=0;m<M;m++)
  {
    for (int n=0;n<N;n++)
      printf(" %8.5f",C[n*M+m]);
    printf("\n");
  }
 #endif

 #if 0
  for (int i=0;i<M*N;i++) C[i] = 0.;

  for (int i=0;i<M;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<K;k++)
    C[j*M+i] += A[i*K+k]*B[k*N+j];
 #endif

 #if 0
  printf("\n C (loops): \n");
  for (int m=0;m<M;m++)
  {
    for (int n=0;n<N;n++)
      printf(" %8.5f",C[m*N+n]);
    printf("\n");
  }
 #endif

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

void mat_times_mat_cpu(double* C, double* A, double* B, int N)
{
  int M = N; int K = N;
  int LDA = N;
  int LDB = N;
  int LDC = N;

  double ALPHA = 1.0;
  double BETA = 0.0;

 //LIBLAS version
 // note: funny order of A/B to get C ordered
  char TB = 'N';
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

  mat_times_mat_cpu(C,A,B,N,N,N);

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

void mat_times_mat_bt_cpu(double* C, double* A, double* B, int M, int N, int K)
{
  return mat_times_mat_bt_cpu(C,A,B,M,N,K,K);
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

void mat_times_mat_bt_cpu(float* C, float* A, float* B, int M, int N, int K)
{
  for (int i=0;i<M*N;i++) C[i] = 0.;
  printf("\n this mat_times_mat not ready \n"); exit(1);

 #if 0
  for (int i=0;i<M;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<K;k++)
    C[i*N+j] += A[i*N+k]*B[j*N+k];
 #endif

  return;
}

void mat_times_mat_bt_cpu(float* C, float* A, float* B, int N)
{
  return mat_times_mat_bt_cpu(C,A,B,N,N,N);
}

int sign(double x)
{
  if (x>0) return 1;
  else if (x<=0) return -1;

  return 0;
}

void cross(double* m, double* r1, double* r2)
{
  m[0] =  r1[1]*r2[2] - r2[1]*r1[2];
  m[1] = -r1[0]*r2[2] + r2[0]*r1[2];
  m[2] =  r1[0]*r2[1] - r2[0]*r1[1];

  return;
}

// Extended version of mat_times_mat_at_cpu with M, N, K parameters
void mat_times_mat_at_cpu(double* C, double* A, double* B, int M, int N, int K)
{
  // C = A^T * B
  // A is K x M, B is K x N, C is M x N
  
  int LDA = M;
  int LDB = N;
  int LDC = N;

  double ALPHA = 1.0;
  double BETA = 0.0;

  char TA = 'T';
  char TB = 'N';
  dgemm_(&TB,&TA,&N,&M,&K,&ALPHA,B,&LDB,A,&LDA,&BETA,C,&LDC);

  return;
}
