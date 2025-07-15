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

void get_murak_grid_zeta(int tid, int size, double* r, double* w, const double zeta, const int m)
{
 //assuming alpha scalar == 1
  double alpha = 1./zeta;
  const int mm1 = m-1;
  //const double w0 = 1./(size-1);
  const double w0 = 1./size;
  const double mal = m*alpha;

  //printf("    murak_grid zeta/alpha: %8.5f %8.5f \n",zeta,alpha);

 #pragma acc parallel loop independent present(r[0:size],w[0:size]) //async(tid+1)
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

  #pragma acc wait
  //acc_wait_all();

  return;
}

void get_murak_grid_zeta(int size, double* r, double* w, const double zeta, const int m)
{
  return get_murak_grid_zeta(-1,size,r,w,zeta,m);
}

//working this in double precision, otherwise limit on grid size
void get_murak_grid_f(int size, float* r, float* w, int Z, const int m)
{
  {
    //printf(" WARNING: single precision murak called for large radial grid \n");
    double rd[size];
    double wd[size];
    #pragma acc enter data create(rd[0:size],wd[0:size])

    get_murak_grid(size,rd,wd,Z,m);
   #pragma acc parallel loop present(r[0:size],rd[0:size])
    for (int j=0;j<size;j++)
      r[j] = rd[j];
   #pragma acc parallel loop present(w[0:size],wd[0:size])
    for (int j=0;j<size;j++)
      w[j] = wd[j];

    #pragma acc exit data delete(rd[0:size],wd[0:size])
    return;
  }

  if (size>165 && 0)
  {
   #pragma acc update self(r[0:size],w[0:size])
    printf("\n r:");
    for (int m=0;m<size;m++)
      printf(" %14.12f",r[m]);

    printf("\n w:");
    for (int m=0;m<size;m++)
      printf(" %14.12f",w[m]);
    printf("\n");
  }

  return;
}

void get_murak_grid(int tid, int size, double* r, double* w, int Z, const int m)
{
  if (Z==0) Z = 1;

 //assuming alpha scalar == 1
  double alpha = alpha_k[Z-1];
  const int mm1 = m-1;
  const double w0 = 1./size;
  const double mal = m*alpha;

 #pragma acc parallel loop independent present(r[0:size],w[0:size]) //async(tid+1)
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
    int npr = 10;
   #pragma acc update self(r[0:size])
    printf("\n r:");
    for (int m=0;m<npr;m++)
      printf(" %14.12f",r[m]);
    printf(" ... ");
    for (int m=size-npr;m<size;m++)
      printf(" %14.12f",r[m]);
  }

  #pragma acc wait
  //acc_wait_all();

  return;
}

void get_murak_grid(int size, double* r, double* w, int Z, const int m)
{
  return get_murak_grid(-1,size,r,w,Z,m);
}
