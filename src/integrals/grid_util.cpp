#include "grid_util.h"
#include "lebedev2.h"
#include <cstdio>
#include <algorithm> 

using namespace std;

void get_angular_grid(int size_ang, double*& ang_g, double*& ang_w)
{
 #if 1
  double* xa = new double[size_ang];
  double* ya = new double[size_ang];
  double* za = new double[size_ang];
  ld_by_order(size_ang,xa,ya,za,ang_w);
  for (int i=0;i<size_ang;i++)
  {
    ang_g[3*i+0] = xa[i];
    ang_g[3*i+1] = ya[i];
    ang_g[3*i+2] = za[i];
  }
  delete [] xa;
  delete [] ya;
  delete [] za;
 #elif 1
  printf("\n\n\n WARNING: using libgrid \n\n\n");
 //this version is not correct
  double* ang = new double[3*size_ang];
  libgrid::lebedev(ang,ang_w,size_ang);
  for (int i=0;i<3;i++)
  for (int j=0;j<size_ang;j++)
    ang_g[j*3+i] = ang[i*size_ang+j];
  //printf(" ang_w: ");
  //for (int j=0;j<size_ang;j++)
  //  printf(" %8.5f",ang_w[j]);
  //printf("\n");
 #else
  arma::mat ang(3,size_ang);
  arma::vec wang(size_ang);
  libgrid::lebedev(ang, wang);
  ang_g = ang.memptr();
  ang_w = wang.memptr();
 #endif

#if 0
  float sumawt = 0.;
  for (int j=0;j<size_ang;j++)
    sumawt += ang_w[j];
  printf(" ang_w: ");
  for (int j=0;j<size_ang;j++)
    printf(" %8.5f",ang_w[j]);
  printf("\n");
  printf(" sum(ang_w): %8.5f \n",sumawt);
#endif

  return;
}

int electron_count(int charge, int natoms, int* atno, int& Nc, int& Na, int& Nb)
{
  int No = 0;
  for (int m=0;m<natoms;m++)
    No += atno[m];
  No -= charge;
  if (No%2==0)
    Na = Nb = No/2;
  else
  {
    Nb = (No-1)/2;
    Na = Nb+1;
  }
  No = min(Na,Nb);

  Nc = 0;
  for (int m=0;m<natoms;m++)
  {
    if (atno[m]>36) 
    {
      printf(" WARNING: is the grid updated for Z>36? \n");
      Nc += 19;
    }
    else if (atno[m]>18)
      Nc += 9;
    else if (atno[m]>10)
      Nc += 5;
    else if (atno[m]>2)
      Nc++;
  }
  if (Nc>No) Nc = 0;

  return No;
}

