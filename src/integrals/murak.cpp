#include "murak.h"

using namespace std;
#include "fp_def.h"

 //working this in double precision, otherwise limit on grid size
void get_murak_grid_f(int size, FP1* r, FP1* w, int Z, const int m)
{
 //assuming alpha scalar == 1
  double alpha = alpha_k[Z-1];
  const int mm1 = m-1;
  const double w0 = 1./size;
  const double mal = m*alpha;

  if (size>165)
  {
    printf(" WARNING: single precision cannot handle this many radial grid points \n");
    //exit(1.);
  }

#if USE_ACC
 #pragma acc parallel loop independent present(r[0:size],w[0:size])
#endif
  for (int n=0;n<size;n++)
  {
   #if 1
    double i1 = (n+0.5)/size;
    double a1 = 1.0 - pow(i1,m);
    double l1 = log(a1);

    double r1 = -alpha*l1;

    double w1 = mal*w0*r1*r1;
    w1 *= pow(i1,mm1) / a1;
   #else
    FP1 i1 = (n+0.5f)/size;
    FP1 a1 = 1.0f - powf(i1,m);
    FP1 l1 = logf(a1);

    FP1 r1 = -alpha*l1;

    FP1 w1 = mal*w0*r1*r1;
    w1 *= powf(i1,mm1) / a1;
   #endif

    r[n] = r1;
    w[n] = w1;
  }

  if (size>165)
  {
   #pragma acc update self(r[0:size])
    printf("\n r:");
    for (int m=0;m<size;m++)
      printf(" %14.12f",r[m]);
  }

  return;
}

void get_murak_grid_f(int size, FP1* r, FP1* w, FP1* er, int Z, FP1 zeta, const int m)
{
 //assuming alpha scalar == 1
  FP1 alpha = alpha_k[Z-1];
  const int mm1 = m-1;
  const FP1 w0 = 1.f/size;
  const FP1 mal = m*alpha;

  if (size>165)
  {
    printf(" ERROR: single precision cannot handle this many radial grid points \n");
    exit(1.);
  }

#if USE_ACC
 #pragma acc parallel loop independent present(r[0:size],w[0:size],er[0:size])
#endif
  for (int n=0;n<size;n++)
  {
    FP1 i1 = (n+0.5f)/size;
    FP1 a1 = 1.0f - powf(i1,m);
    FP1 l1 = logf(a1);

    FP1 r1 = -alpha*l1;

    FP1 w1 = mal*w0*r1*r1;
    w1 *= powf(i1,mm1) / a1;

    r[n] = r1;
    w[n] = w1;
    er[n] = expf(-zeta*r1);
  }

 #if 0
 #pragma acc update self(r[0:size])
  printf("\n r:");
  for (int m=0;m<size;m++)
    printf(" %14.12f",r[m]);
 #endif 

  return;
}

void get_murak_grid(int size, double* r, double* w, double* er, int Z, double zeta, const int m)
{
  //printf("\n debug: Mura-Knowles grid \n");

 //assuming alpha scalar == 1
  double alpha = alpha_k[Z-1];
  const int mm1 = m-1;
  const double w0 = 1./size;
  const double mal = m*alpha;

#if USE_ACC
 //#pragma acc kernels
#endif
  for (int n=0;n<size;n++)
  {
    double i1 = (n+0.5)/size;
    double a1 = 1.0 - pow(i1,m);
    double l1 = log(a1);

    double r1 = -alpha*l1;
    r[n] = r1;

    double w1 = mal*w0*r1*r1;
    w1 *= pow(i1,mm1) / a1;
    w[n] = w1;

    er[n] = exp(-zeta*r1);
  }

 #if 0
  printf("\n r:");
  for (int m=0;m<size;m++)
    printf(" %8.5f",r[m]);
  printf("\n w:");
  for (int m=0;m<size;m++)
    printf(" %8.5f",w[m]);
  printf("\n exp(-zeta*r):");
 #endif
 #if 0
  for (int m=0;m<size;m++)
    printf(" %8.5f",er[m]);
  printf("\n");
 #endif

  return;
}
