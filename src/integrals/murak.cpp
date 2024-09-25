#include "murak.h"

using namespace std;

void get_eumac_grid(int size, double* r, double* w, const double rmax, const int m)
{
  printf("\n Euler Mac grid not ready \n");
  exit(-1);
  double rm = rmax/100.;
  double rm3 = rm*rm*rm;
  double w0 = rm3*(size+1);

  if (m==1)
 #pragma acc parallel loop present(r[0:size],w[0:size])
  for (int i=0;i<size;i++)
  {
    double u = i+1;
    double v = size-i;

    double r1 = rm*u/v;
    double w1 = w0*u*u*pow(v,-4.);

    r[i] = r1;
    w[i] = w1;
  }

  if (m==2)
 #pragma acc parallel loop present(r[0:size],w[0:size])
  for (int i=0;i<size;i++)
  {
    double u = i+1;
    double v = size-i;

    double r1 = rm * u*u / v / v;
    double w1 = 2.*w0*pow(u,5.)*pow(v,-7);

    r[i] = r1;
    w[i] = w1;
  }

  if (1)
  {
   #pragma acc update self(r[0:size])
    printf("\n reu:");
    for (int m=0;m<size;m++)
      printf(" %6.3e",r[m]);
    printf("\n");
  }
}

void get_murak_grid_zeta(int size, double* r, double* w, const double zeta, const int m)
{
 //assuming alpha scalar == 1
  double alpha = 1./zeta;
  const int mm1 = m-1;
  //const double w0 = 1./(size-1);
  const double w0 = 1./size;
  const double mal = m*alpha;

  //printf("  murak_grid zeta/alpha: %8.5f %8.5f \n",zeta,alpha);

#if USE_ACC
 #pragma acc parallel loop independent present(r[0:size],w[0:size])
#endif
  for (int n=0;n<size;n++)
  {
    //double i0 = (n-0.5)/size; if (i0<0.) i0 = 0.;
    //double a0 = 1.0 - pow(i0,m);
    //double l0 = log(a0);
    //double r0 = -alpha*l0;

    double i1 = (n+0.5)/size;
    double a1 = 1.0 - pow(i1,m);
    double l1 = log(a1);

    double r1 = -alpha*l1;

    double w1 = mal*w0*r1*r1;
    w1 *= pow(i1,mm1) / a1;
    //double w1 = r1*r1*(r1-r0);

    //printf("  r/w: %5.3e %5.3e \n",r1,w1);
    //printf("  rjw: %5.3e %5.3e \n",r1,mal*w0*r1*r1*pow(i1,mm1)/a1);

    r[n] = r1;
    w[n] = w1;
  }

  if (0)
  {
   #pragma acc update self(r[0:size])
    printf("\n r:");
    for (int m=0;m<size;m++)
      printf(" %6.3e",r[m]);
    printf("\n");
  }

  return;
}

 //working this in double precision, otherwise limit on grid size
void get_murak_grid_f(int size, float* r, float* w, int Z, const int m)
{
  if (Z==0) Z = 1;

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
    float i1 = (n+0.5f)/size;
    float a1 = 1.0f - powf(i1,m);
    float l1 = logf(a1);

    float r1 = -alpha*l1;

    float w1 = mal*w0*r1*r1;
    w1 *= powf(i1,mm1) / a1;
   #endif

    r[n] = r1;
    w[n] = w1;
  }

  if (size>165 && 0)
  {
   #pragma acc update self(r[0:size])
    printf("\n r:");
    for (int m=0;m<size;m++)
      printf(" %14.12f",r[m]);
  }

  return;
}

void get_murak_grid(int size, double* r, double* w, int Z, const int m)
{
  if (Z==0) Z = 1;

 //assuming alpha scalar == 1
  double alpha = alpha_k[Z-1];
  const int mm1 = m-1;
  const double w0 = 1./size;
  const double mal = m*alpha;

#if USE_ACC
 #pragma acc parallel loop independent present(r[0:size],w[0:size])
#endif
  for (int n=0;n<size;n++)
  {
    double i1 = (n+0.5)/size;
    double a1 = 1.0 - pow(i1,m);
    double l1 = log(a1);

    double r1 = -alpha*l1;

    double w1 = mal*w0*r1*r1;
    w1 *= pow(i1,mm1) / a1;
    r[n] = r1;
    w[n] = w1;
  }

  if (size>165 && 0)
  {
   #pragma acc update self(r[0:size])
    printf("\n r:");
    for (int m=0;m<size;m++)
      printf(" %14.12f",r[m]);
  }

  return;
}

void get_murak_grid_f(int size, float* r, float* w, float* er, int Z, float zeta, const int m)
{
  if (Z==0) Z = 1;

 //assuming alpha scalar == 1
  float alpha = alpha_k[Z-1];
  const int mm1 = m-1;
  const float w0 = 1.f/size;
  const float mal = m*alpha;

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
    float i1 = (n+0.5f)/size;
    float a1 = 1.0f - powf(i1,m);
    float l1 = logf(a1);

    float r1 = -alpha*l1;

    float w1 = mal*w0*r1*r1;
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
  if (Z==0) Z = 1;

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
