#include "prosph.h"

#define ZERO_TOL 1.e-10
#define QUAD_ORDER 8

//using Cartesian product quadrature (2D) in quad.cpp

//GPU:
// 1. accelerate half_split_ps?
// 2. asynchronous?

// get_cfn turned off

void copy_to_all_gpu(int ngpu, int s1, double* A, int include_first);
void swap_grids_gpu(int gs, double*& grid, double*& gridm, double*& wt, double* gridn, double* gridmn, double* wtn);
void auto_crash();

//splits every grid point into octrees
void full_split_ps(int gs, double* grid, double* gridm, double* wt, double a, int gs2, double* grid2, double* gridm2, double* wt2, int prl)
{
  const int div = 9*6; //1+8 grid points
  double a3 = a*a*a;

  int ng = 0;
  for (int j=0;j<gs;j++)
  {
    double mu0  = gridm[6*j+0]; double nu0  = gridm[6*j+1]; double phi0  = gridm[6*j+2];
    double dmu0 = gridm[6*j+3]; double dnu0 = gridm[6*j+4]; double dphi0 = gridm[6*j+5];

    double dmu = dmu0*0.5; double dnu = dnu0*0.5; double dphi = dphi0*0.5;
    //double mu1 = mu0-dmu; double mu2 = mu0;
    double mu1  = mu0-dmu*0.5;   double mu2  = mu0+dmu*0.5;
    double nu1  = nu0-dnu*0.5;   double nu2  = nu0+dnu*0.5;
    double phi1 = phi0-dphi*0.5; double phi2 = phi0+dphi*0.5;

    int i1 = j*div;
    grid2[i1+0] = grid[6*j+0]; grid2[i1+1] = grid[6*j+1]; grid2[i1+2] = grid[6*j+2]; //copy xyz from original grid
    gridm2[i1+0] = gridm[6*j+0]; gridm2[i1+1] = gridm[6*j+1]; gridm2[i1+2] = gridm[6*j+2]; //copy mnp/dmnp
    gridm2[i1+3] = gridm[6*j+3]; gridm2[i1+4] = gridm[6*j+4]; gridm2[i1+5] = gridm[6*j+5];
    wt2[9*j] = wt[j];

   //new mnp values
    gridm2[i1+6]  = mu1; gridm2[i1+7]  = nu1; gridm2[i1+8]  = phi1;
    gridm2[i1+12] = mu1; gridm2[i1+13] = nu1; gridm2[i1+14] = phi2;
    gridm2[i1+18] = mu1; gridm2[i1+19] = nu2; gridm2[i1+20] = phi1;
    gridm2[i1+24] = mu1; gridm2[i1+25] = nu2; gridm2[i1+26] = phi2;
    gridm2[i1+30] = mu2; gridm2[i1+31] = nu1; gridm2[i1+32] = phi1;
    gridm2[i1+36] = mu2; gridm2[i1+37] = nu1; gridm2[i1+38] = phi2;
    gridm2[i1+42] = mu2; gridm2[i1+43] = nu2; gridm2[i1+44] = phi1;
    gridm2[i1+48] = mu2; gridm2[i1+49] = nu2; gridm2[i1+50] = phi2;

    //printf("   mnp0: %8.5f %8.5f %8.5f \n",mu0,nu0,phi0);

    //double dmnp = dmu*dnu*dphi;
    for (int k=1;k<9;k++)
    {
      int i2 = i1+6*k;
      double mu = gridm2[i2+0]; double nu = gridm2[i2+1]; double phi = gridm2[i2+2];

      double sinhm = sinh(mu); double coshm = cosh(mu);
      double sinn  = sin(nu);  double cosn  = cos(nu);

      double x = a*sinhm*sinn*cos(phi);
      double y = a*sinhm*sinn*sin(phi);
      double z = a*coshm*cosn-a;
      //double wt1 = sinhm*sinn*(sinhm*sinhm+sinn*sinn)*dmnp;
      double wt1 = a3*ps_dV(mu-0.5*dmu,mu+0.5*dmu,nu-0.5*dnu,nu+0.5*dnu) * dphi;

      gridm2[i2+3] = dmu; gridm2[i2+4] = dnu; gridm2[i2+5] = dphi;
      grid2[i2+0] = x; grid2[i2+1] = y; grid2[i2+2] = z; //xyz of new pts
      wt2[9*j+k] = wt1;

      //printf("     mnp: %8.5f %8.5f %8.5f \n",mu,nu,phi);
    }

    ng += 9;
  }
  //printf("   full_split new size: %4i \n",ng);

  if (ng<1000 && prl>1)
  {
    printf("\n full_split grid: \n");
    for (int j=0;j<ng;j++)
      printf("  %8.5f %8.5f %8.5f  wt: %8.5f \n",grid2[6*j+0],grid2[6*j+1],grid2[6*j+2],wt2[j]);
      //printf("  %8.5f %8.5f %8.5f  wt: %8.5f \n",gridm2[6*j+0],gridm2[6*j+1],gridm2[6*j+2],wt2[j]);
  }
  //printf("  grid2 size: %3i \n",ng);


  return;
}

void half_split_ps(int ws, int gs, double* grid, double* gridm, double* wt, double a, double* grid2, double* gridm2, double* wt2, int prl)
{
  const int div = 3*6; //1+2 grid points
  double a3 = a*a*a;
  int gs6 = 6*gs;
  int gsb = 3*gs;
  int gsb6 = 6*gsb;

  if (ws==1)
 #pragma acc parallel loop present(grid[0:gs6],gridm[0:gs6],wt[0:gs],grid2[0:gsb6],gridm2[0:gsb6],wt2[0:gsb])
  for (int j=0;j<gs;j++)
  {
    int i1 = j*div;

    double mu0  = gridm[6*j+0]; double nu0  = gridm[6*j+1]; double phi0  = gridm[6*j+2];
    double dmu0 = gridm[6*j+3]; double dnu0 = gridm[6*j+4]; double dphi0 = gridm[6*j+5];

    double dmu = dmu0;  double dnu = dnu0; double dphi = dphi0;
    double mu1  = mu0;  double mu2  = mu0;
    double nu1  = nu0;  double nu2  = nu0;
    double phi1 = phi0; double phi2 = phi0;
    {
      dmu = dmu0*0.5;
      mu1 = mu0-dmu*0.5; mu2 = mu0+dmu*0.5;
    }

    grid2[i1+0] = grid[6*j+0]; grid2[i1+1] = grid[6*j+1]; grid2[i1+2] = grid[6*j+2]; //copy xyz from original grid
    gridm2[i1+0] = gridm[6*j+0]; gridm2[i1+1] = gridm[6*j+1]; gridm2[i1+2] = gridm[6*j+2]; //copy mnp/dmnp
    gridm2[i1+3] = gridm[6*j+3]; gridm2[i1+4] = gridm[6*j+4]; gridm2[i1+5] = gridm[6*j+5];
    wt2[3*j] = wt[j];

   //new mnp values
    gridm2[i1+6]  = mu1; gridm2[i1+7]  = nu1; gridm2[i1+8]  = phi1;
    gridm2[i1+12] = mu2; gridm2[i1+13] = nu2; gridm2[i1+14] = phi2;

   #pragma acc loop independent
    for (int k=1;k<3;k++)
    {
      int i2 = i1+6*k;
      double mu = gridm2[i2+0]; double nu = gridm2[i2+1]; double phi = gridm2[i2+2];

      double sinhm = sinh(mu); double coshm = cosh(mu);
      double sinn  = sin(nu);  double cosn  = cos(nu);

      double x = a*sinhm*sinn*cos(phi);
      double y = a*sinhm*sinn*sin(phi);
      double z = a*coshm*cosn-a;
      double wt1 = a3*ps_dV(mu-0.5*dmu,mu+0.5*dmu,nu-0.5*dnu,nu+0.5*dnu) * dphi;

      gridm2[i2+3] = dmu; gridm2[i2+4] = dnu; gridm2[i2+5] = dphi;
      grid2[i2+0] = x; grid2[i2+1] = y; grid2[i2+2] = z; //xyz of new pts
      wt2[3*j+k] = wt1;
    }
  }

  if (ws==2)
 #pragma acc parallel loop present(grid[0:gs6],gridm[0:gs6],wt[0:gs],grid2[0:gsb6],gridm2[0:gsb6],wt2[0:gsb])
  for (int j=0;j<gs;j++)
  {
    int i1 = j*div;

    double mu0  = gridm[6*j+0]; double nu0  = gridm[6*j+1]; double phi0  = gridm[6*j+2];
    double dmu0 = gridm[6*j+3]; double dnu0 = gridm[6*j+4]; double dphi0 = gridm[6*j+5];

    double dmu = dmu0;  double dnu = dnu0; double dphi = dphi0;
    double mu1  = mu0;  double mu2  = mu0;
    double nu1  = nu0;  double nu2  = nu0;
    double phi1 = phi0; double phi2 = phi0;
    {
      dnu = dnu0*0.5;
      nu1 = nu0-dnu*0.5; nu2 = nu0+dnu*0.5;
    }

    grid2[i1+0] = grid[6*j+0]; grid2[i1+1] = grid[6*j+1]; grid2[i1+2] = grid[6*j+2]; //copy xyz from original grid
    gridm2[i1+0] = gridm[6*j+0]; gridm2[i1+1] = gridm[6*j+1]; gridm2[i1+2] = gridm[6*j+2]; //copy mnp/dmnp
    gridm2[i1+3] = gridm[6*j+3]; gridm2[i1+4] = gridm[6*j+4]; gridm2[i1+5] = gridm[6*j+5];
    wt2[3*j] = wt[j];

   //new mnp values
    gridm2[i1+6]  = mu1; gridm2[i1+7]  = nu1; gridm2[i1+8]  = phi1;
    gridm2[i1+12] = mu2; gridm2[i1+13] = nu2; gridm2[i1+14] = phi2;

   #pragma acc loop independent
    for (int k=1;k<3;k++)
    {
      int i2 = i1+6*k;
      double mu = gridm2[i2+0]; double nu = gridm2[i2+1]; double phi = gridm2[i2+2];

      double sinhm = sinh(mu); double coshm = cosh(mu);
      double sinn  = sin(nu);  double cosn  = cos(nu);

      double x = a*sinhm*sinn*cos(phi);
      double y = a*sinhm*sinn*sin(phi);
      double z = a*coshm*cosn-a;
      double wt1 = a3*ps_dV(mu-0.5*dmu,mu+0.5*dmu,nu-0.5*dnu,nu+0.5*dnu) * dphi;

      gridm2[i2+3] = dmu; gridm2[i2+4] = dnu; gridm2[i2+5] = dphi;
      grid2[i2+0] = x; grid2[i2+1] = y; grid2[i2+2] = z; //xyz of new pts
      wt2[3*j+k] = wt1;
    }
  }

  if (ws==3)
 #pragma acc parallel loop present(grid[0:gs6],gridm[0:gs6],wt[0:gs],grid2[0:gsb6],gridm2[0:gsb6],wt2[0:gsb])
  for (int j=0;j<gs;j++)
  {
    int i1 = j*div;

    double mu0  = gridm[6*j+0]; double nu0  = gridm[6*j+1]; double phi0  = gridm[6*j+2];
    double dmu0 = gridm[6*j+3]; double dnu0 = gridm[6*j+4]; double dphi0 = gridm[6*j+5];

    double dmu = dmu0;  double dnu = dnu0; double dphi = dphi0;
    double mu1  = mu0;  double mu2  = mu0;
    double nu1  = nu0;  double nu2  = nu0;
    double phi1 = phi0; double phi2 = phi0;
    {
      dphi = dphi0*0.5;
      phi1 = phi0-dphi*0.5; phi2 = phi0+dphi*0.5;
    }

    grid2[i1+0] = grid[6*j+0]; grid2[i1+1] = grid[6*j+1]; grid2[i1+2] = grid[6*j+2]; //copy xyz from original grid
    gridm2[i1+0] = gridm[6*j+0]; gridm2[i1+1] = gridm[6*j+1]; gridm2[i1+2] = gridm[6*j+2]; //copy mnp/dmnp
    gridm2[i1+3] = gridm[6*j+3]; gridm2[i1+4] = gridm[6*j+4]; gridm2[i1+5] = gridm[6*j+5];
    wt2[3*j] = wt[j];

   //new mnp values
    gridm2[i1+6]  = mu1; gridm2[i1+7]  = nu1; gridm2[i1+8]  = phi1;
    gridm2[i1+12] = mu2; gridm2[i1+13] = nu2; gridm2[i1+14] = phi2;

   #pragma acc loop independent
    for (int k=1;k<3;k++)
    {
      int i2 = i1+6*k;
      double mu = gridm2[i2+0]; double nu = gridm2[i2+1]; double phi = gridm2[i2+2];

      double sinhm = sinh(mu); double coshm = cosh(mu);
      double sinn  = sin(nu);  double cosn  = cos(nu);

      double x = a*sinhm*sinn*cos(phi);
      double y = a*sinhm*sinn*sin(phi);
      double z = a*coshm*cosn-a;
      double wt1 = a3*ps_dV(mu-0.5*dmu,mu+0.5*dmu,nu-0.5*dnu,nu+0.5*dnu) * dphi;

      gridm2[i2+3] = dmu; gridm2[i2+4] = dnu; gridm2[i2+5] = dphi;
      grid2[i2+0] = x; grid2[i2+1] = y; grid2[i2+2] = z; //xyz of new pts
      wt2[3*j+k] = wt1;
    }
  }

  int ng = gsb;

  #pragma acc update self(grid2[0:gsb6],gridm2[0:gsb6],wt2[0:gsb])

  if (ng<300 && prl>1)
  {
    printf("\n half_split grid: \n");
    for (int j=0;j<ng;j++)
      printf("  %8.5f %8.5f %8.5f  wt: %8.5f \n",grid2[6*j+0],grid2[6*j+1],grid2[6*j+2],wt2[j]);
      //printf("  %8.5f %8.5f %8.5f  wt: %8.5f \n",gridm2[6*j+0],gridm2[6*j+1],gridm2[6*j+2],wt2[j]);
  }

  return;
}

double evaluate_integral_and_split_2ss_gpu(const double eps, int gso, double* val, double* wt, int* split)
{
  //int ngn = 0;
  double vt = 0.;

 #pragma acc parallel loop present(val[0:3*gso],split[0:gso]) reduction(+:vt)
  for (int j=0;j<gso;j++)
  {
    int i1 = 3*j;

    double v0 = val[i1+0];
    double v1 = val[i1+1]+val[i1+2];

    double diff = fabs(v1-v0);

    if (diff>eps)
      split[j] = 1;
    vt += v1;
  } //calculate error

  return vt;
}

int evaluate_integral_and_split(const double eps, int ss, int gso, double* val, double* wt, int* split, double& inactsum, double& actsum, double& maxd)
{
  double vtn = 0.;
  double vt0 = 0.;
  double vt1 = 0.;

  double mind = 1000.;
  maxd = 0.;

 //also need to account for nonsplit terms?
  int ngn = 0;
  for (int j=0;j<gso;j++)
  {
    int i1 = ss*j;

    double v0 = val[i1+0];
    double v1 = 0.;
    for (int k=1;k<ss;k++)
      v1 += val[i1+k];

    double w0 = wt[i1];
    double w1 = 0.;
    for (int k=1;k<ss;k++)
      w1 += wt[i1+k];

    double diff = fabs(v1-v0);
    if (diff>maxd) maxd = diff;
    if (diff<mind) mind = diff;

    if (0)
    //if (prl>1 || gso<500) 
    //if (j%100000==0 || (gso<10000 && j%100==0))
    //if ((j%1000)==0)
    {
      printf("        point %8i  v0/1: %11.8f %11.8f  diff: %5.2e  wt0/1: %6.3e %6.3e  diff: %5.2e ",j,v0,v1,diff,w0,w1,fabs(w1-w0));

      if (diff<eps) printf("*");
      printf("\n");
    }

    if (diff>eps)
    {
      split[j] = 1;
      ngn++;
    }
    else
      vtn += v1; //accumulate "inactive" terms

    vt0 += v0;
    vt1 += v1;
  } //calculate octree error

  inactsum = vtn;
  actsum = vt1-vtn;

  printf("  ints: %12.8f  %12.8f  max/mind: %5.2e %5.2e \n",vt0,vt1,maxd,mind);

  return ngn;
}

//input old grid+val, output new grid
//change ordering of grid to save time on copy step?
int adaptive_split_ps(const double eps, int ws, int ss, int gs, double a, double* grid, double* gridm, double* wt, double* val, bool& grid_updated, double& inactsum, double& actsum, double& maxd, double*& grid2, double*& gridm2, double*& wt2, int prl)
{
  //printf("  as_ps. eps: %5.1e \n",eps);

  int ss1 = ss-1;
  int gso = gs/ss; //assuming octree structure
  const int div = ss*6; //1+8 grid points
  const int divsp = ss1*6;

 //determine what dV elements to split up
  int* split = new int[gs]();
  int ngn = evaluate_integral_and_split(eps,ss,gso,val,wt,split,inactsum,actsum,maxd);
 //cannot inactivate, since pt may split over a different DOF
  if (ss<9) { actsum += inactsum; inactsum = 0.; }

  if (ngn==0) grid_updated = 0; else grid_updated = 1;

  if (eps<=0.) return 0;
  if (ss==9 && ngn<1) return 0;
  //if (ss<9 && ngn<1) return 0;

 //split up the grid
  int gs2 = ngn*ss1;
  if (ss<9) //keep "inactive" grid pts
    gs2 = (gso-ngn)+ngn*ss1;
  int gssp = gs2;
  double* gridsp = new double[6*gs2];
  double* gridmsp = new double[6*gs2];
  double* wtsp = new double[gs2];

  if (prl>1) printf("     gs2: %4i  gssp: %4i \n",gs2,gssp);

 //copy over pts from prior grid
  if (ss<9)
  {
    int i0 = 0;
    for (int i=0;i<gso;i++)
    if (!split[i])
    {
      for (int k=0;k<6;k++)
      {
        gridsp[i0*6+k] = grid[div*i+k];
        gridmsp[i0*6+k] = gridm[div*i+k];
      }
      wtsp[i0] = wt[ss*i];
      i0++;
    }
    //printf("            i0(inact): %3i \n",i0);
    for (int i=0;i<gso;i++)
    if (split[i])
    {
      for (int k=6;k<6*ss;k++)
      {
        gridsp[i0*6+k-6] = grid[div*i+k];
        gridmsp[i0*6+k-6] = gridm[div*i+k];
      }
      for (int k=1;k<ss;k++)
        wtsp[i0+k-1] = wt[ss*i+k];
      i0 += ss1;
    }
    //printf("          i0(  act): %3i \n",i0);
  }

  if (ss==9)
  for (int i=0;i<ngn;i++)
  {
    int j = split[i];
    int i0 = i*divsp;
    int i1 = j*div;
    for (int k=6;k<div;k++)
    {
      gridsp[i0+k-6] = grid[i1+k];
      gridmsp[i0+k-6] = gridm[i1+k];
    }
    for (int k=1;k<ss;k++)
      wtsp[ss1*i+k-1] = wt[ss*j+k];
  }

 //precise amount of space allocated here
  grid2 = new double[div*gssp];
  gridm2 = new double[div*gssp];
  wt2 = new double[ss*gssp];

  if (ss==9)
    full_split_ps(gssp,gridsp,gridmsp,wtsp,a,gs2,grid2,gridm2,wt2,prl);
  else
    half_split_ps(ws,gssp,gridsp,gridmsp,wtsp,a,grid2,gridm2,wt2,prl);

  delete [] gridsp;
  delete [] gridmsp;
  delete [] wtsp;
  delete [] split;

  return gssp*ss;
}

//not obvious whether serial loops should be on cpu or gpu
int adaptive_split_ps_2ss_gpu(const double eps, int ws, int gs, double a, double* grid, double* gridm, double* wt, double* val, bool& grid_updated, double& actsum, double*& grid2, double*& gridm2, double*& wt2, int prl)
{
  const int ss = 3;
  const int ss1 = ss-1;
  int gso = gs/ss;
  const int div = ss*6; //1+split grid points
  const int divsp = ss1*6; //for octree (ss==9)

 //determine what dV elements to split up
  int* split = new int[gso];
  #pragma acc enter data create(split[0:gso])
 #pragma acc parallel loop present(split[0:gso])
  for (int j=0;j<gso;j++)
    split[j] = 0;

  actsum = evaluate_integral_and_split_2ss_gpu(eps,gso,val,wt,split);
  //#pragma acc exit data copyout(split[0:gso])

  int ngn = 0;
 #pragma acc parallel loop present(split[0:gso]) reduction(+:ngn)
  for (int j=0;j<gso;j++)
    ngn += split[j];

  if (ngn==0) grid_updated = 0; else grid_updated = 1;

  if (eps<=0.) return 0;
  if (ss==9 && ngn<1) return 0;

 //split up the grid
  int gs2 = ngn*ss1;
  if (ss<9) //keep "inactive" grid pts
    gs2 = (gso-ngn)+ngn*ss1;
  int gssp = gs2;
  double* gridsp = new double[6*gs2];
  double* gridmsp = new double[6*gs2];
  double* wtsp = new double[gs2];
  #pragma acc enter data create(gridsp[0:6*gs2],gridmsp[0:6*gs2],wtsp[0:gs2])

  if (prl>1) printf("     gso: %4i  gs2: %4i  gssp: %4i  \n",gso,gs2,gssp);

 //copy over pts from prior grid
  {
    int i0 = 0;
   #pragma acc serial loop present(split[0:gso],gridsp[0:6*gs2],wtsp[0:gs2],gridmsp[0:6*gs2],grid[0:6*gs],gridm[0:6*gs],wt[0:gs])
    for (int i=0;i<gso;i++)
    if (!split[i])
    {
     #pragma acc loop
      for (int k=0;k<6;k++)
      {
        gridsp[i0*6+k] = grid[div*i+k];
        gridmsp[i0*6+k] = gridm[div*i+k];
      }
      wtsp[i0] = wt[ss*i];
      i0++;
    }
    //printf("            i0(inact): %3i \n",i0);
   #pragma acc serial loop present(split[0:gso],gridsp[0:6*gs2],gridmsp[0:6*gs2],wtsp[gs2],grid[0:6*gs],gridm[0:6*gs],wt[0:gs])
    for (int i=0;i<gso;i++)
    if (split[i])
    {
      //printf(" - %i: %i ",i,i0);
     #pragma acc loop
      for (int k=6;k<6*ss;k++)
      {
        gridsp[i0*6+k-6] = grid[div*i+k];
        gridmsp[i0*6+k-6] = gridm[div*i+k];
      }
     #pragma acc loop
      for (int k=1;k<ss;k++)
        wtsp[i0+k-1] = wt[ss*i+k];
      i0 += ss1;
    }
    //printf("          i0(  act): %3i \n",i0);
  }

  #pragma acc exit data delete(split[0:gso])
  delete [] split;

 //precise amount of space allocated here
  double* grid2p = new double[div*gssp];
  double* gridm2p = new double[div*gssp];
  double* wt2p = new double[ss*gssp];

  #pragma acc enter data create(grid2p[0:div*gssp],gridm2p[0:div*gssp],wt2p[0:ss*gssp])

  grid2 = grid2p; gridm2 = gridm2p; wt2 = wt2p;

  if (ss==9)
    full_split_ps(gssp,gridsp,gridmsp,wtsp,a,gs2,grid2,gridm2,wt2,prl);
  else
    half_split_ps(ws,gssp,gridsp,gridmsp,wtsp,a,grid2,gridm2,wt2,prl);

  #pragma acc exit data delete(gridsp[0:6*gs2],gridmsp[0:6*gs2],wtsp[0:gs2])

  delete [] gridsp;
  delete [] gridmsp;
  delete [] wtsp;

  return gssp*ss;
}

void get_two_grids(int gs, double* grid1, double* grid2, double* grid, double A1, double B1, double C1, double A2, double B2, double C2)
{
  copy_grid(gs,grid1,grid);
  recenter_grid_zero(gs,grid1,-A1,-B1,-C1);
  copy_grid(gs,grid2,grid);
  recenter_grid_zero(gs,grid2,-A2,-B2,-C2);
}

void get_three_grids(int gs, double* grid1, double* grid2, double* grid3, double* grid, double A1, double B1, double C1, double A2, double B2, double C2, double A3, double B3, double C3)
{
  copy_grid(gs,grid1,grid);
  recenter_grid_zero(gs,grid1,-A1,-B1,-C1);
  copy_grid(gs,grid2,grid);
  recenter_grid_zero(gs,grid2,-A2,-B2,-C2);
  copy_grid(gs,grid3,grid);
  recenter_grid_zero(gs,grid3,-A3,-B3,-C3);
}

void swap_grids_gpu(int gs, double*& grid, double*& gridm, double*& wt, double* gridn, double* gridmn, double* wtn)
{
  int gs6 = 6*gs;

  double* gridd = grid; double* gridmd = gridm; double* wtd = wt;
  #pragma acc exit data delete(gridd[0:gs6],gridmd[0:gs6],wtd[0:gs])

  delete [] grid; delete [] gridm; delete [] wt;

  grid = gridn; gridm = gridmn; wt = wtn;
}

void swap_grids(double*& grid, double*& gridm, double*& wt, double* gridn, double* gridmn, double* wtn)
{
  delete [] grid; delete [] gridm; delete [] wt;
  grid = gridn; gridm = gridmn; wt = wtn;
}

double get_cfn(int n1, int l1, double zt1, int n2, int l2, double zt2, double Z1, double Z2)
{
 //CPMZ disabled this
  return 1.;

  const double cm = 0.1;
  double cr = 1. + (n1-l1-1)*cm + (n2-l2-1)*cm;
  if (zt1*zt2<=1.) cr *= 1.1;
  return cr;
}

void assemble_sh_vals(int n1, int l1, int m1, double zt1, double norm1, int n2, int l2, int m2, double zt2, double norm2, int gs, double* grid1, double* grid2, double* val1, double* val2, double* val1p, double* val2p)
{
  int gs3 = gs*3;
 #pragma acc parallel loop present(val1[0:gs])
  for (int j=0;j<gs;j++)
    val1[j] = norm1;
 #pragma acc parallel loop present(val2[0:gs])
  for (int j=0;j<gs;j++)
    val2[j] = norm2;
 #pragma acc parallel loop present(val1p[0:gs3])
  for (int j=0;j<gs3;j++)
    val1p[j] = norm1;
 #pragma acc parallel loop present(val2p[0:gs3])
  for (int j=0;j<gs3;j++)
    val2p[j] = norm2;

  eval_shd(-1,gs,grid1,val1,n1,l1,m1,zt1); //nu=0  ->+z
  eval_shd(-1,gs,grid2,val2,n2,l2,m2,zt2); //nu=180->-z
  eval_pd(-1,gs,grid1,val1p,n1,l1,m1,zt1);
  eval_pd(-1,gs,grid2,val2p,n2,l2,m2,zt2);

  return;
}

//for 2c Coulomb integrals
void assemble_vc_vals(int n1, int l1, int m1, double zt1, double norm1v, int gs, double* grid1, double* val1, double* val1p, double* valt)
{
  int gs3 = gs*3;
 #pragma acc parallel loop present(valt[0:gs3])
  for (int j=0;j<gs;j++)
    valt[j] = norm1v;
 #pragma acc parallel loop present(val1[0:gs])
  for (int j=0;j<gs;j++)
    val1[j] = norm1v;
 #pragma acc parallel loop present(val1p[0:gs3])
  for (int j=0;j<gs3;j++)
    val1p[j] = norm1v;

 //S in valt, V in val1
  eval_sh_3rd(gs,grid1,valt,n1,l1,m1);
  eval_inr_r12(-1,gs,grid1,val1,n1,l1,zt1);

   //first term: SdV
 #pragma acc parallel loop present(val1p[0:gs3],valt[0:gs3])
  for (int j=0;j<gs;j++)
    val1p[3*j+0] = val1p[3*j+1] = val1p[3*j+2] = valt[j];
  eval_inr_d(gs,grid1,val1p,n1,l1,zt1);

   //second term: VdS
 #pragma acc parallel loop present(valt[0:gs3],val1[0:gs])
  for (int j=0;j<gs;j++)
    valt[3*j+0] = valt[3*j+1] = valt[3*j+2] = val1[j];
  eval_dp_3rd(-1,gs,grid1,valt,n1,l1,m1);

 #pragma acc parallel loop present(val1p[0:gs3],valt[0:gs3])
  for (int j=0;j<gs3;j++)
    val1p[j] += valt[j];

   //complete the SV term
  eval_sh_3rd(gs,grid1,val1,n1,l1,m1);

  return;
}

void assemble_vc_vals(int a12, int b12, int c12, int n1, int l1, int m1, double zt1, double norm1, int n2, int l2, int m2, double zt2, double norm2, int n3, int l3, int m3, double zt3, double norm3v, int gs, double* grid1, double* grid2, double* val1, double* val1p, double* val2, double* val2p, double* valt, double* valtp)
{
 //a12/b12/c12 --> which center the basis/potential is on
  double* grida = grid1; if (a12==2) grida = grid2;
  double* gridb = grid1; if (b12==2) gridb = grid2;
  double* gridc = grid1; if (c12==2) gridc = grid2;

  //v1 --> 1+2 ftns
  //v2 --> 3 ftns
  int gs3 = gs*3;
  double n12 = norm1*norm2;
 #pragma acc parallel loop present(val1[0:gs])
  for (int j=0;j<gs;j++)
    val1[j] = n12;
 #pragma acc parallel loop present(val2[0:gs])
  for (int j=0;j<gs;j++)
    val2[j] = n12;

  eval_shd(-1,gs,grida,val1,n1,l1,m1,zt1);

 //mu dnu
 #pragma acc parallel loop present(valtp[0:gs3],val1[0:gs])
  for (int j=0;j<gs;j++)
    valtp[3*j+0] = valtp[3*j+1] = valtp[3*j+2] = val1[j];
  eval_pd(-1,gs,gridb,valtp,n2,l2,m2,zt2);

 //nu dmu
  eval_shd(-1,gs,gridb,val2,n2,l2,m2,zt2);
 #pragma acc parallel loop present(valt[0:gs3],val2[0:gs])
  for (int j=0;j<gs;j++)
    valt[3*j+0] = valt[3*j+1] = valt[3*j+2] = val2[j];
  eval_pd(-1,gs,grida,valt,n1,l1,m1,zt1);

 //mu dnu + nu dmu
 #pragma acc parallel loop present(val1p[0:gs3],valt[0:gs3],valtp[0:gs3])
  for (int j=0;j<gs3;j++)
    val1p[j] = valt[j] + valtp[j];

 //mu nu
  eval_shd(-1,gs,gridb,val1,n2,l2,m2,zt2);

 #pragma acc parallel loop present(val2[0:gs])
  for (int j=0;j<gs;j++)
    val2[j] = norm3v;
  eval_sh_3rd(gs,gridc,val2,n3,l3,m3);

 //SdV
 #pragma acc parallel loop present(val2p[0:gs3],val2[0:gs])
  for (int j=0;j<gs;j++)
    val2p[3*j+0] = val2p[3*j+1] = val2p[3*j+2] = val2[j];
  eval_inr_d(gs,gridc,val2p,n3,l3,zt3);

 //VdS
 #pragma acc parallel loop present(valt[0:gs3])
  for (int j=0;j<gs;j++)
    valt[j] = norm3v;
  eval_inr_r12(-1,gs,gridc,valt,n3,l3,zt3);
 #pragma acc parallel loop present(valtp[0:gs3],valt[0:gs3])
  for (int j=0;j<gs;j++)
    valtp[3*j+0] = valtp[3*j+1] = valtp[3*j+2] = valt[j];
  eval_dp_3rd(-1,gs,gridc,valtp,n3,l3,m3);

 #pragma acc parallel loop present(val2p[0:gs3],valtp[0:gs3])
  for (int j=0;j<gs3;j++)
    val2p[j] += valtp[j];

 //finishing VS
  eval_inr_r12(-1,gs,gridc,val2,n3,l3,zt3);

  return;
}

//3c integrals, deal with special cases for phase
bool get_triple_phi(int m1, int m2, int m3, double& v1, double& v2)
{
  vector<int> ms;
  ms.push_back(m1);
  ms.push_back(m2);
  ms.push_back(m3);
  sort(ms.begin(),ms.end());

  bool found = 0;
  if (ms[0]==-4 && ms[1]==-1 && ms[2]==3)
  {
    v1 = 1.147923073128459347;
    v2 = 0.8642620985957808795;
    found = 1;
  }
  if (ms[0]==-4 && ms[1]==-3 && ms[2]==1)
  {
    v1 = 0.422873253666437272;
    v2 = 0.8642620985957808795;
    found = 1;
  }
  if (ms[0]==-5 && ms[1]==-2 && ms[2]==3)
  {
    v1 = 0.9507384372688643;
    v2 = 0.9057228709674523;
    found = 1;
  }
  if (ms[0]==-5 && ms[1]==-3 && ms[2]==2)
  {
    v1 = PI*0.5;
    v2 = 1.;
    found = 1;
  }
  return found;
}

void evaluate_over_grid_quad_3c(int qo, int n1, int l1, int m1, double zt1, double norm1, int n2, int l2, int m2, double zt2, double norm2, int n3, int l3, int m3, double zt3, double norm3v,
                        double z0, double A1, double B1, double C1, double A2, double B2, double C2, double A3, double B3, double C3,
                        int gs, double* grid, double* gridm, double* wt, double* rot, double* val)
{
 //grid comes in in 1+2 order, but doesn't matter here

  int qo2 = 2*qo;
  int qos = qo*qo*qo;
  double* Qx = new double[qo2];
  double* Qy = new double[qo2];
  double* Qz = new double[qo2];

  get_quad(qo,Qx);
  get_quad(qo,Qy);
  get_quad(qo,Qz);

  double norm123 = norm1*norm2*norm3v;

  int gsq = qos*gs;
  int gsq6 = 6*gsq;

  double* gridq = new double[6*gsq];
  double* gridqr = new double[6*gsq];
  double* wtq = new double[gsq];
  double* valq = new double[gsq];

  double* grid1 = new double[gsq6];
  double* grid2 = new double[gsq6];
  double* grid3 = new double[gsq6];

  #pragma acc enter data copyin(Qx[0:qo2],Qy[0:qo2],Qz[0:qo2])
  #pragma acc enter data create(gridq[0:gsq6],gridqr[0:gsq6],wtq[0:gsq],valq[0:gsq],grid1[0:gsq6],grid2[0:gsq6],grid3[0:gsq6])

  //#pragma acc update self(gridm[0:6*gs],wt[0:gs])
  //#pragma acc update self(val[0:gs])

  double* grida = grid1;
  double* gridb = grid2;
  double* gridc = grid3;

  quad_grid_munuphi(-1,qo,qo,qo,0,z0,Qx,Qy,Qz,gs,0,gridm,gridq,wtq);

  reorient_grid(-1,z0,gsq,gridq,gridqr,rot);

  get_three_grids(gsq,grid1,grid2,grid3,gridqr,A1,B1,C1,A2,B2,C2,A3,B3,C3);

 #pragma acc parallel loop present(valq[0:gsq])
  for (int k=0;k<gsq;k++)
    valq[k] = norm123;

  eval_shd(-1,gsq,grida,valq,n1,l1,m1,zt1);
  eval_shd(-1,gsq,gridb,valq,n2,l2,m2,zt2);

  eval_sh_3rd(gsq,gridc,valq,n3,l3,m3);
  eval_inr_r12(-1,gsq,gridc,valq,n3,l3,zt3);

 #pragma acc parallel loop present(val[0:gs],valq[0:gsq],wtq[0:gsq])
  for (int j=0;j<gs;j++)
  {
    int i1 = j*qos;

    double v1 = 0.;
   #pragma acc loop reduction(+:v1)
    for (int k=0;k<qos;k++)
      v1 += valq[i1+k]*wtq[i1+k];

    val[j] = v1;
  }

  #pragma acc exit data delete(Qx[0:qo2],Qy[0:qo2],Qz[0:qo2])
  #pragma acc exit data delete(gridq[0:gsq6],gridqr[0:gsq6],wtq[0:gsq],valq[0:gsq],grid1[0:gsq6],grid2[0:gsq6],grid3[0:gsq6])

  delete [] grid1;
  delete [] grid2;
  delete [] grid3;

  delete [] Qx;
  delete [] Qy;
  delete [] Qz;
  delete [] gridq;
  delete [] gridqr;
  delete [] wtq;
  delete [] valq;

  return;
}

bool set_4c_pass_idp(int m, int mmin, vector<int> ms, double& idp)
{
  bool pass = 0;

  if (idp!=0. && (mmin==-m || mmin>=0)) pass = 1;

  if (ms[0]==-1 && ms[1]==-1 && ms[2]==1 && ms[3]==1) { pass = 1; idp *= 4.; }
  if (ms[0]==-2 && ms[1]==-2 && ms[2]==-1 && ms[3]==-1) { idp *= 2.; }
  if (ms[0]==-2 && ms[1]==-2 && ms[2]==1 && ms[3]==1) { idp *= 2.; }

  if (m==4 && ms[0]==-2 && ms[1]==-2 && ms[2]==2 && ms[3]==2) { pass = 1; idp *= 4.; }

  if (m==2 && ms[0]==-3 && ms[1]==-3 && ms[2]==3 && ms[3]==3) { pass = 1; idp *= 4.; }
  if (m==4 && ms[0]==-1 && ms[1]==-1 && ms[2]==1 && ms[3]==3) { pass = 1; idp *= 19.31370849898476; }
  if (m==4 && ms[0]==-2 && ms[1]==-1 && ms[2]==2 && ms[3]==3) { pass = 1; idp *= 13.65685424949238; }
  if (m==2 && ms[0]==-1 && ms[1]==-1 && ms[2]==3 && ms[3]==3) { pass = 1; idp *= 4.; }

  return pass;
}

void do_4c_integrals_ps(const double epsilon, double cf, int nomp, int ss, int maxsteps, double z0, int natoms, double* coords, int* atno, vector<vector<double> > basis, int nmu, int nnu, int nphi, double* g, int prl)
{
  int quad_order = read_int("QUAD");
  if (quad_order<=0) quad_order = QUAD_ORDER;

  //if (nomp>1) { printf("\n ERROR: do_4c_integrals is one gpu only \n"); exit(-1); }

  const double zero_tol = ZERO_TOL;
  int N = basis.size();
  int N2 = N*N;
  int N3 = N*N2;

  printf("  do_4c basis.size: %2i \n",N);
  if (basis.size()<1) { exit(-1); }

  int nphi1 = nphi;

  int lmax = 0;
  for (int i=0;i<N;i++)
  {
    int l1 = (int)basis[i][1];
    if (lmax<l1)
      lmax = l1;
  }
  if (lmax>3)
  {
    printf("\n ERROR: 4c integrals require l<4 (l: %i) \n",lmax); exit(-1);
  }

  int gs = nmu*nnu*nphi1;
  int qo = quad_order;
  int qo2 = qo*2;
  int qos = qo*qo;
  int gsq = qos*gs;
  int gsq6 = 6*gsq;

  double* grid = new double[6*gs];
  double* gridm = new double[6*gs];
  double* wt = new double[gs];

 //AO evaluated on gridq
  double* valab = new double[gsq];
  double** val1 = new double*[N];
  for (int i=0;i<N;i++)
    val1[i] = new double[gsq];

  for (int i1=0;i1<N2*N2;i1++) g[i1] = 0.;

 //4c overlap integrals
  int at1 = 0; int at2 = 1;

  double A1 = coords[3*at1+0]; double B1 = coords[3*at1+1]; double C1 = coords[3*at1+2];
  double A2 = coords[3*at2+0]; double B2 = coords[3*at2+1]; double C2 = coords[3*at2+2];

 //create a quadrature grid
  double* Qx = new double[qo2];
  double* Qy = new double[qo2];
  get_quad(qo,Qx);
  get_quad(qo,Qy);

 //fixed size grid (q) centered at 0. grid1/2 on atom 1/2
  double* gridq = new double[gsq6];
  double* gridqr = new double[gsq6];
  double* wtq = new double[gsq];

  double* grid1 = new double[gsq6];
  double* grid2 = new double[gsq6];

  printf("\n nomp in do_4c_ol: %i \n",nomp);

 //alloc on all gpus
  for (int n=0;n<nomp;n++)
  {
    int tid = n; //omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc enter data copyin(Qx[0:qo2],Qy[0:qo2])
    #pragma acc enter data create(gridq[0:gsq6],gridqr[0:gsq6],wtq[0:gsq],grid1[0:gsq6],grid2[0:gsq6])

    #pragma acc enter data create(grid[0:6*gs],gridm[0:6*gs],wt[0:gs])
    #pragma acc enter data create(val1[0:N][0:gsq],valab[0:gsq])
  }
  acc_set_device_num(0,acc_device_nvidia);

  //auto_crash();

  double phi_phase = 0.;

 //create grid and make a copy to rotate xy later
 #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);

   //CPMZ hard-coded
    double cfn = 1.1*cf;
    initialize_ps_coords_batch(-1,0,1,z0,cfn,nmu,nnu,nphi1,phi_phase,grid,gridm,wt,prl);
    quad_grid_munu(qo,qo,z0,Qx,Qy,gs,gridm,gridq,wtq);

    #pragma acc parallel loop present(gridq[0:gsq6],gridqr[0:gsq6])
    for (int j=0;j<gsq6;j++)
      gridqr[j] = gridq[j];
  }
  acc_set_device_num(0,acc_device_nvidia);

 //this will be uneven effort distribution
  int mmax = 4;
  if (prl>0) printf("  outer m loop");
  for (int m=mmax;m>=0;m--)
  {
    double phi_phase = 0.;
    if (m>0) phi_phase = PI/m/2.;

    if (prl>0) { printf("."); fflush(stdout); }
    double cospp = cos(phi_phase); double sinpp = sin(phi_phase);
   #pragma omp parallel for schedule(static,1) num_threads(nomp)
    for (int n=0;n<nomp;n++)
    {
      int tid = omp_get_thread_num();
      acc_set_device_num(tid,acc_device_nvidia);

     #pragma acc parallel loop present(gridqr[0:gsq6],gridq[0:gsq6])
      for (int j=0;j<gsq;j++)
      {
        double x1 = gridq[6*j+0];
        double xn = cospp*x1;
        double yn = sinpp*x1;
        gridqr[6*j+0] = xn;
        gridqr[6*j+1] = yn;
      }
      get_two_grids(gsq,grid1,grid2,gridqr,A1,B1,C1,A2,B2,C2);
    }
    acc_set_device_num(0,acc_device_nvidia);

   #pragma omp parallel for schedule(static,1) num_threads(nomp)
    for (int n=0;n<nomp;n++)
    for (int i1=0;i1<N;i1++)
    {
      int tid = omp_get_thread_num();
      acc_set_device_num(tid,acc_device_nvidia);

     //evaluate all basis ftns on grid
      double* grida = grid1; if (basis[i1][9]==1.) grida = grid2; //basis is a double. comparing to int

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zt1 = basis1[3];
      double norm1 = basis1[4];

     #pragma acc parallel loop present(val1[0:N][0:gsq])
      for (int j=0;j<gsq;j++)
        val1[i1][j] = norm1;
      eval_shd(-1,gsq,grida,val1[i1],n1,l1,m1,zt1);
    }

   #pragma omp parallel for schedule(static,1) num_threads(nomp)
    for (int i1=0;i1<N;i1++)
    for (int i2=0;i2<=i1;i2++)
    {
      int tid = omp_get_thread_num();
      acc_set_device_num(tid,acc_device_nvidia);

      double* vala = val1[i1];
      double* valb = val1[i2];

     #pragma acc parallel loop present(valab[0:gsq],vala[0:gsq],valb[0:gsq],wtq[0:gsq])
      for (int j=0;j<gsq;j++)
        valab[j] = vala[j]*valb[j]*wtq[j];

      for (int i3=0;i3<N;i3++)
      for (int i4=0;i4<=i3;i4++)
      {
        int m1 = (int)basis[i1][2]; int m2 = (int)basis[i2][2]; int m3 = (int)basis[i3][2]; int m4 = (int)basis[i4][2];
        double idp = get_idp_4b_m3(m1,m2,m3,m4);
        vector<int> ms; ms.push_back(m1); ms.push_back(m2); ms.push_back(m3); ms.push_back(m4);
        sort(ms.begin(),ms.end());
        int mmin = ms[0];

        bool pass = set_4c_pass_idp(m,mmin,ms,idp);

        if (pass)
        {
          double* valc = val1[i3];
          double* vald = val1[i4];

          double v1 = 0.;
         #pragma acc parallel loop present(valab[0:gsq],valc[0:gsq],vald[0:gsq]) reduction(+:v1)
          for (int j=0;j<gsq;j++)
            v1 += valab[j]*valc[j]*vald[j];

          if (prl>0) printf("  %2i %2i %2i %2i  mmin/m: %2i %2i pp: %5.3f idp: %5.3f  v1: %8.5f \n",i1,i2,i3,i4,mmin,m,phi_phase,idp,v1);

          if (fabs(v1)>1.e-14)
          {
            v1 *= idp;

            g[i1*N3+i2*N2+i3*N+i4] = v1;
            g[i1*N3+i2*N2+i4*N+i3] = v1;
            g[i2*N3+i1*N2+i3*N+i4] = v1;
            g[i2*N3+i1*N2+i4*N+i3] = v1;

            g[i3*N3+i4*N2+i1*N+i2] = v1;
            g[i3*N3+i4*N2+i2*N+i1] = v1;
            g[i4*N3+i3*N2+i1*N+i2] = v1;
            g[i4*N3+i3*N2+i1*N+i2] = v1;
          }
        } //idp!=.0

      } //loop i4<=i3
    } //loop i2<=i1 (main integral work loop)
    acc_set_device_num(0,acc_device_nvidia);

  } //outer rotation phi
  if (prl>0) printf("\n");

  for (int n=0;n<nomp;n++)
  {
    int tid = n; //omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc exit data delete(val1[0:N][0:gsq],valab[0:gsq])
    #pragma acc exit data delete(gridq[0:gsq6],gridqr[0:gsq6],wtq[0:gsq],grid1[0:gsq6],grid2[0:gsq6])
    #pragma acc exit data delete(Qx[0:qo2],Qy[0:qo2])
    #pragma acc exit data delete(grid[0:6*gs],gridm[0:6*gs],wt[0:gs])
  }
  acc_set_device_num(0,acc_device_nvidia);

  delete [] valab;
  for (int i=0;i<N;i++)
    delete [] val1[i];
  delete [] val1;

  delete [] grid; delete [] gridm; delete [] wt;
  delete [] gridq; delete [] gridqr; delete [] wtq;
  delete [] grid1; delete [] grid2;

  delete [] Qx;
  delete [] Qy;


  //auto_crash();

  return;
}

void evaluate_over_grid_3c(int n, int m, int p, int n1, int l1, int m1, double zt1, double norm1, int n2, int l2, int m2, double zt2, double norm2, int n3, int l3, int m3, double zt3, double norm3v,
                        double z0, double A1, double B1, double C1, double A2, double B2, double C2, int gs, double* grid, double* gridm, double* wt, double* val)
{
  return;
}

void evaluate_over_grid_quad_En_3c(int qo, int n1, int l1, int m1, double zt1, double norm1, int n2, int l2, int m2, double zt2, double norm2, 
                        double z0, double A1, double B1, double C1, double A2, double B2, double C2, double A3, double B3, double C3, double Z3,
                        int gs, double* grid, double* gridm, double* wt, double* rot, double* val)
{
 //grid comes in in 1+2 order, but doesn't matter here

  int qo2 = 2*qo;
  int qos = qo*qo*qo;
  double* Qx = new double[qo2];
  double* Qy = new double[qo2];
  double* Qz = new double[qo2];

  get_quad(qo,Qx);
  get_quad(qo,Qy);
  get_quad(qo,Qz);

  double norm12 = norm1*norm2;

  int gsq = qos*gs;
  int gsq6 = 6*gsq;

  double* gridq = new double[6*gsq];
  double* gridqr = new double[6*gsq];
  double* wtq = new double[gsq];
  double* valq = new double[gsq];

  double* grid1 = new double[gsq6];
  double* grid2 = new double[gsq6];

  #pragma acc enter data copyin(Qx[0:qo2],Qy[0:qo2],Qz[0:qo2])
  #pragma acc enter data create(gridq[0:gsq6],gridqr[0:gsq6],wtq[0:gsq],valq[0:gsq],grid1[0:gsq6],grid2[0:gsq6])

  #pragma acc update self(gridm[0:6*gs],wt[0:gs])

  quad_grid_munuphi(-1,qo,qo,qo,0,z0,Qx,Qy,Qz,gs,0,gridm,gridq,wtq);

  reorient_grid(-1,z0,gsq,gridq,gridqr,rot);

  get_two_grids(gsq,grid1,grid2,gridqr,A1,B1,C1,A2,B2,C2);


 #pragma acc parallel loop present(valq[0:gsq])
  for (int k=0;k<gsq;k++)
    valq[k] = norm12;

  //nuclear attraction integrals
  {
    eval_shd(-1,gsq,grid1,valq,n1,l1,m1,zt1);
    eval_shd(-1,gsq,grid2,valq,n2,l2,m2,zt2);

   #pragma acc parallel loop present(valq[0:gsq],grid1[0:gsq],grid2[0:gsq]) //note grid1/2 here
    for (int j=0;j<gsq;j++)
    {
      double x1 = grid1[6*j+0] - A3;
      double y1 = grid1[6*j+1] - B3;
      double z1 = grid1[6*j+2] - C3;
      double Rn3 = sqrt(x1*x1+y1*y1+z1*z1);
      double ne = -Z3/Rn3;
      valq[j] *= ne;
    }
  }

 #pragma acc parallel loop present(val[0:gs],valq[0:gsq],wtq[0:gsq])
  for (int j=0;j<gs;j++)
  {
    int i1 = j*qos;

    double v1 = 0.;
   #pragma acc loop reduction(+:v1)
    for (int k=0;k<qos;k++)
      v1 += valq[i1+k]*wtq[i1+k];

    val[j] = v1;
  }

  #pragma acc exit data delete(Qx[0:qo2],Qy[0:qo2],Qz[0:qo2])
  #pragma acc exit data delete(gridq[0:gsq6],gridqr[0:gsq6],wtq[0:gsq],valq[0:gsq],grid1[0:gsq6],grid2[0:gsq6])

  delete [] grid1;
  delete [] grid2;

  delete [] Qx;
  delete [] Qy;
  delete [] Qz;
  delete [] gridq;
  delete [] gridqr;
  delete [] wtq;
  delete [] valq;

  return;
}

void get_2c_position(double* coordn, double* rot)
{
 //this ftn assumes two distinct atoms in coordn (on cpu)

  //#pragma acc enter data copyin(coordn[0:9])
  gen_total_rot_n2(2,coordn,rot);

  int prl = 0;
  if (prl>1)
  {
    //#pragma acc update self(rot[0:9])
    printf("\n rot(0): \n");
    print_square(3,rot);
  }

 //inverse of rot via transpose
  //#pragma acc parallel loop present(rot[0:9])
  for (int j=0;j<3;j++)
  for (int k=0;k<j;k++)
  {
    double r1 = rot[j*3+k];
    double r2 = rot[k*3+j];

    rot[j*3+k] = r2;
    rot[k*3+j] = r1;
  }

  //#pragma acc update self(rot[0:9])

  if (prl>1)
  {

    printf("\n rot(i): \n");
    print_square(3,rot);
  }

  //#pragma acc exit data delete(coordn[0:9])

  return;
}

void get_2c_position(double& A1, double& B1, double& C1, double& A2, double& B2, double& C2, double* rot)
{
 //center at atom 1
  double coordn[9];
  coordn[0] = 0.; coordn[1] = 0.; coordn[2] = 0.;
  coordn[3] = A2-A1; coordn[4] = B2-B1; coordn[5] = C2-C1;

 //save centers at atom 1
  A1 = coordn[0]; B1 = coordn[1]; C1 = coordn[2];
  A2 = coordn[3]; B2 = coordn[4]; C2 = coordn[5];

  int prl = 0;
  if (prl>1)
  {
    printf(" at1: %8.5f %8.5f %8.5f \n",A1,B1,C1);
    printf(" at2: %8.5f %8.5f %8.5f \n",A2,B2,C2);
  }

  if (prl>1)
  {
    printf(" XYZ: \n");
    for (int k=0;k<2;k++)
      printf("  %8.5f %8.5f %8.5f \n",coordn[3*k+0],coordn[3*k+1],coordn[3*k+2]);
  }

  bool ab_match = 0;
  if (coordn[3]==0. && coordn[4]==0. && coordn[5]==0.)
    ab_match = 1;

  if (ab_match)
  {
   //#pragma acc parallel loop present(rot[0:9])
    for (int j=0;j<9;j++)
      rot[j] = 0.;
   //#pragma acc parallel loop present(rot[0:9])
    for (int j=0;j<3;j++)
      rot[3*j+j] = 1.;

    //#pragma acc update self(rot[0:9])

    return;
  }

  return get_2c_position(coordn,rot);
}

void get_3c_position(double* coordn, double* rot)
{
  bool ab_match = 0; if (coordn[3]==0. && coordn[4]==0. && coordn[5]==0.) ab_match = 1;
  bool abc_match = 0; if (ab_match && coordn[6]==0. && coordn[7]==0. && coordn[8]==0.) abc_match = 1;

  if (abc_match)
  {
   //#pragma acc parallel loop present(rot[0:9])
    for (int j=0;j<9;j++)
      rot[j] = 0.;
   //#pragma acc parallel loop present(rot[0:9])
    for (int j=0;j<3;j++)
      rot[3*j+j] = 1.;

    //#pragma acc update self(rot[0:9])

    return;
  }
  else if (ab_match)
  {
    coordn[3] = coordn[6];
    coordn[4] = coordn[7];
    coordn[5] = coordn[8];

    coordn[6] = 0.;
    coordn[7] = 0.;
    coordn[8] = 0.;
  }

  int prl = 0;
  if (prl>1)
  {
    printf(" XYZ(rearranged. match? %i %i): \n",(int)ab_match,(int)abc_match);
    for (int k=0;k<3;k++)
      printf("  %8.5f %8.5f %8.5f \n",coordn[3*k+0],coordn[3*k+1],coordn[3*k+2]);
  }

  //#pragma acc enter data copyin(coordn[0:9])
  gen_total_rot_n3(3-ab_match-abc_match,coordn,rot);

  if (prl>0)
  {
    //#pragma acc update self(rot[0:9])
    printf("\n rot(0): \n");
    print_square(3,rot);
  }

 //inverse of rot via transpose
  //#pragma acc parallel loop present(rot[0:9])
  for (int j=0;j<3;j++)
  //#pragma acc loop independent
  for (int k=0;k<j;k++)
  {
    double r1 = rot[j*3+k];
    double r2 = rot[k*3+j];

    rot[j*3+k] = r2;
    rot[k*3+j] = r1;
  }

 //CPMZ check?
  //#pragma acc update self(rot[0:9])

  if (prl>0)
  {
    printf("\n rot(i): \n");
    print_square(3,rot);
  }

  //#pragma acc exit data delete(coordn[0:9])

  return;
}

void get_4c_position(double* coordn, double* rot)
{
 //assumes all 4 centers are unique
  int prl = 0;

  //#pragma acc enter data copyin(coordn[0:9])
  gen_total_rot_n3(3,coordn,rot);

 //inverse of rot via transpose
  //#pragma acc parallel loop present(rot[0:9])
  for (int j=0;j<3;j++)
  //#pragma acc loop independent
  for (int k=0;k<j;k++)
  {
    double r1 = rot[j*3+k];
    double r2 = rot[k*3+j];

    rot[j*3+k] = r2;
    rot[k*3+j] = r1;
  }

  //#pragma acc update self(rot[0:9])

  if (prl>0)
  {
    printf("\n rot(i): \n");
    print_square(3,rot);
  }

  //#pragma acc exit data delete(coordn[0:9])

  return;
}

void get_3c_position(double& A1, double& B1, double& C1, double& A2, double& B2, double& C2, double& A3, double& B3, double& C3, double* rot)
{
  int prl = 0;
  if (prl>1)
  {
    printf(" at1: %8.5f %8.5f %8.5f \n",A1,B1,C1);
    printf(" at2: %8.5f %8.5f %8.5f \n",A2,B2,C2);
    printf(" at3: %8.5f %8.5f %8.5f \n",A3,B3,C3);
  }

 //center at atom 1
  double coordn[9];
  coordn[0] = 0.;    coordn[1] = 0.;    coordn[2] = 0.;
  coordn[3] = A2-A1; coordn[4] = B2-B1; coordn[5] = C2-C1;
  coordn[6] = A3-A1; coordn[7] = B3-B1; coordn[8] = C3-C1;

 //save centers at atom 1
  A1 = coordn[0]; B1 = coordn[1]; C1 = coordn[2];
  A2 = coordn[3]; B2 = coordn[4]; C2 = coordn[5];
  A3 = coordn[6]; B3 = coordn[7]; C3 = coordn[8];

  return get_3c_position(coordn,rot);
}


int count_unique(int a, int b, int c)
{
  vector<int> all;
  all.push_back(a);
  all.push_back(b);
  all.push_back(c);
  sort(all.begin(),all.end());

  all.erase(unique(all.begin(),all.end()), all.end());

  return all.size();
}

void do_3c_integrals_ps(const double epsilon, double cf, int nomp, int ss, int maxsteps, int natoms, double* coords, int* atno, int* n2i, int* a2i, vector<vector<double> > basis, vector<vector<double> > basis_aux, int nmu, int nnu, int nphi, double* En, double* C, int prl)
{
  if (natoms>3) { printf("\n TESTING: natoms>3 in do_3c_integrals_ps \n"); }
  //prl = 1;

  int quad_order = read_int("QUAD");
  if (quad_order<=0) quad_order = QUAD_ORDER;

  const double zero_tol = ZERO_TOL;
  int N = basis.size();
  int N2 = N*N;
  int Naux = basis_aux.size();
  int nna = N*Naux;

  printf("\n  do_3c basis.size: %2i  aux: %2i  nsteps: %2i  eps: %4.1e \n",N,Naux,maxsteps,epsilon);
  if (basis.size()<1) { exit(-1); }

 //number of extra cells for third atom
  int nx3 = 0;
  int gs0 = nmu*nnu*nphi + nx3;
  double* grid0 = new double[6*gs0];
  double* gridm0 = new double[6*gs0];
  double* wt0 = new double[gs0];
  double rot[9];

  for (int n=0;n<nomp;n++)
  {
    int tid = n; //omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);
    #pragma acc enter data create(rot[0:9])
    #pragma acc enter data create(grid0[0:6*gs0],gridm0[0:6*gs0],wt0[0:gs0])
  }
  acc_set_device_num(0,acc_device_nvidia);

  double gpumem = (double)acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  double togb = 1./1024./1024./1024.;
  printf("  initial gpu memory available: %6.3f GB \n",gpumem*togb);

  if (C!=NULL) for (int i1=0;i1<N2*Naux;i1++) C[i1] = 0.;

 //3c Coulomb integrals
  if (C!=NULL)
  for (int n=0;n<natoms;n++)
  for (int m=n;m<natoms;m++) //m>=n
  for (int p=0;p<natoms;p++)
  {
    int s1 = 0; if (n>0) s1 = n2i[n-1]; int s2 = n2i[n];
    int s3 = 0; if (m>0) s3 = n2i[m-1]; int s4 = n2i[m];
    int s5 = 0; if (p>0) s5 = a2i[p-1]; int s6 = a2i[p];

   //avoiding redundant integrals (i1>=i2 if same atom)
    int s1r = s1; if (n==m) s1r = s2;

    double Z1 = (double)atno[n]; double Z2 = (double)atno[m]; double Z3 = (double)atno[p];
    double A1 = coords[3*n+0]; double B1 = coords[3*n+1]; double C1 = coords[3*n+2];
    double A2 = coords[3*m+0]; double B2 = coords[3*m+1]; double C2 = coords[3*m+2];
    double A3 = coords[3*p+0]; double B3 = coords[3*p+1]; double C3 = coords[3*p+2];

    int nphi1 = nphi;

    get_3c_position(A1,B1,C1,A2,B2,C2,A3,B3,C3,rot);
    copy_to_all_gpu(nomp,9,rot,1);
    double coordn[9];
    coordn[0] = A1; coordn[1] = B1; coordn[2] = C1;
    coordn[3] = A2; coordn[4] = B2; coordn[5] = C2;
    coordn[6] = A3; coordn[7] = B3; coordn[8] = C3;

    int ncenters = count_unique(n,m,p);

    double z0 = sqrt(A2*A2+B2*B2+C2*C2)*0.5;
    if (m==n) z0 = sqrt(A3*A3+B3*B3+C3*C3)*0.5;
    if (m==n && n==p) z0 = 0.5;

   #if 0
    int at1 = min(n,min(m,p)); int at2 = max(n,max(m,p));
   //CPMZ fix this
    int nphi1 = nphi;

    double Z1 = (double)atno[at1]; double Z2 = (double)atno[at2];
    double A1 = coords[3*at1+0]; double B1 = coords[3*at1+1]; double C1 = coords[3*at1+2];
    double A2 = coords[3*at2+0]; double B2 = coords[3*at2+1]; double C2 = coords[3*at2+2];
   #endif

    if (prl>1) printf("  3c Coulomb. s1-6: %2i %2i %2i %2i %3i %3i  (at123: %i %i %i) \n",s1,s2,s3,s4,s5,s6,n,m,p);

   #pragma omp parallel for schedule(dynamic,1) num_threads(nomp)
    for (int i2=s3;i2<s4;i2++)
    for (int i1=min(s1r,i2);i1<s2;i1++)
    for (int i3=s5;i3<s6;i3++)
    if (C[i1*nna+i2*Naux+i3]==0.)
    {
      int tid = omp_get_thread_num();
      acc_set_device_num(tid,acc_device_nvidia);

      double* grid; double* gridm; double* wt;
      double* val = NULL;
      double sum = 0.;

      vector<double> basis1 = basis[i1];
      vector<double> basis2 = basis[i2];
      vector<double> basis3 = basis_aux[i3];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zt1 = basis1[3];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zt2 = basis2[3];
      int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zt3 = basis3[3];
      double norm1 = basis1[4]; double norm2 = basis2[4]; double norm3 = basis3[4];
      double norm3v = 4.*PI/(2.*l3+1)*norm3;

      //double idp = get_idp_m4(m1,m2,m3);
      double idp = 1.;
      //if (natoms==2) idp = get_idp_m5(m1,m2,m3);

      double phi_phase = 0.;
     #if 0
      if (natoms==2 && idp!=0.)
      if (m1<0 || m2<0 || m3<0)
      {
       //y is zero, so align with maximum value of spherical harmonic
        int me = min(m1,min(m2,m3));
        phi_phase = PI/fabs(me)/2.;

        double v1 = phi_phase; double v2 = 1.;
        bool update_phi = get_triple_phi(m1,m2,m3,v1,v2);
        if (update_phi)
        {
          phi_phase = v1;
          idp /= v2;
        }
      }
     #endif

      int nbatch = 1;
     #if 0
      if (maxsteps>12) nbatch = 4;
      if (maxsteps>16) nbatch = 12;
      if (maxsteps>20) nbatch = 16;
     #endif

     //deal with long-range functions
      double cfn = cf*get_cfn(n1,l1,zt1,n2,l2,zt2,Z1,Z2);

      if (prl>0 && tid==0 && idp!=0.)
      {
        printf("\n GPU %i working on %2i-%2i-%2i \n",tid,i1,i2,i3);
        printf("  nlmz: %i %i %2i  %8.5f   idp/pi: %5.3f  phase: %8.5f \n",n1,l1,m1,zt1,idp/PI,phi_phase/PI);
        printf("  nlmz: %i %i %2i  %8.5f   idp/pi: %5.3f  phase: %8.5f \n",n2,l2,m2,zt2,idp/PI,phi_phase/PI);
        printf("  nlmz: %i %i %2i  %8.5f   idp/pi: %5.3f  phase: %8.5f \n",n3,l3,m3,zt3,idp/PI,phi_phase/PI);
      }

      if (idp!=0.)
      for (int wb=0;wb<nbatch;wb++)
      {
        int gs = ss*gs0; if (ncenters<3) gs = ss*(gs0-nx3); int gs6 = 6*gs;
        if (maxsteps>1)
        {
          grid = new double[gs6]; gridm = new double[gs6]; wt = new double[gs];
          #pragma acc enter data create(grid[0:gs6],gridm[0:gs6],wt[0:gs])
        }

        if (ncenters==3)
          initialize_ps_coords_3c(-1,cfn,nmu,nnu,nphi1,0.,coordn,grid0,gridm0,wt0,rot,prl-tid);
        else
          initialize_ps_coords_batch(-1,wb,nbatch,z0,cfn,nmu,nnu,nphi1,0.,grid0,gridm0,wt0,prl-tid);

        if (maxsteps>1)
          half_split_ps(1,gs0,grid0,gridm0,wt0,z0,grid,gridm,wt,prl);
        else
        {
         //handle non-adaptive case
          val = new double[gs0];
          #pragma acc enter data create(val[0:gs0])
          evaluate_over_grid_quad_3c(quad_order,n1,l1,m1,zt1,norm1,n2,l2,m2,zt2,norm2,n3,l3,m3,zt3,norm3v,z0,A1,B1,C1,A2,B2,C2,A3,B3,C3,gs0,grid0,gridm0,wt0,rot,val);

          sum = 0.;
          #pragma acc parallel loop present(val[0:gs0]) reduction(+:sum)
          for (int j=0;j<gs0;j++)
            sum += val[j];

          #pragma acc exit data delete(val[0:gs0])
          delete [] val; val = NULL;
        }

       //adapt the grid specifically for 3c Coulomb
        int wsc = 0; int nzero = 0; int nconv = 0; double sump = 1000.; int gsp = 1000000;
        if (maxsteps>1)
        for (int ns=0;ns<maxsteps;ns++)
        {
          int ws = (wsc%3)+1;
          if (ws==4) { ws = 1; wsc = 0; }
          wsc++;

          if (val!=NULL) delete [] val;
          val = new double[gs];

          #pragma acc enter data create(val[0:gs])

          //evaluate_over_grid_3c(n,m,p,n1,l1,m1,zt1,norm1,n2,l2,m2,zt2,norm2,n3,l3,m3,zt3,norm3v,z0,A1,B1,C1,A2,B2,C2,gs,grid,gridm,wt,val);
          evaluate_over_grid_quad_3c(quad_order,n1,l1,m1,zt1,norm1,n2,l2,m2,zt2,norm2,n3,l3,m3,zt3,norm3v,z0,A1,B1,C1,A2,B2,C2,A3,B3,C3,gs,grid,gridm,wt,rot,val);

          double* gridn = NULL; double* gridmn = NULL; double* wtn = NULL;
          double asum = 0.; bool gup = 0;
          int gs1 = adaptive_split_ps_2ss_gpu(epsilon,ws,gs,z0,grid,gridm,wt,val,gup,asum,gridn,gridmn,wtn,prl);
          sum = asum*idp;

          #pragma acc exit data delete(val[0:gs])

          swap_grids_gpu(gs,grid,gridm,wt,gridn,gridmn,wtn);
          gs = gs1; //gs6 = 6*gs;

          if (tid==0 && prl>0) printf("  step %2i  ws: %i  sum: %14.10f  gs: %7i \n",ns,ws,sum,gs);

          //if (fabs(sum-sump)<zero_tol && ws==1) nconv++;
          if (gs==gsp) nconv++; else nconv = 0;
          if (nconv>2) { if (tid==0 && prl>0) printf(" ** converged ** \n"); break; }
          if (fabs(sum)<zero_tol) nzero++; else nzero = 0;
          if (nzero>3) { if (tid==0 && prl>0) printf(" ** converged to zero ** \n"); break; }
          sump = sum;
          gsp = gs;

        } //loop ns, adapting grid
        if (tid==0 && prl>0) printf("\n");

        if (prl>0) printf("  final 3c(co). sum: %14.10f  tid: %i \n",sum,tid);

       //i1/i2->regular basis, i3->aux basis
       #pragma omp atomic
        C[i1*nna+i2*Naux+i3] += sum;

        if (maxsteps>1)
        {
          #pragma acc exit data delete(grid[0:gs6],gridm[0:gs6],wt[0:gs])

          delete [] grid;
          delete [] gridm;
          delete [] wt;
        }

      } //nonzero integral

     //cleanup
      if (val!=NULL) delete [] val;

    } //loop i1,i2
  } //outer loop over i1,i2
  acc_set_device_num(0,acc_device_nvidia);


 //3c e-n attraction integrals
  if (En!=NULL)
  for (int n=0;n<natoms;n++)
  for (int m=n;m<natoms;m++) //m>=n
  for (int p=0;p<natoms;p++)
  if (p!=n && p!=m)
  {
    int s1 = 0; if (n>0) s1 = n2i[n-1]; int s2 = n2i[n];
    int s3 = 0; if (m>0) s3 = n2i[m-1]; int s4 = n2i[m];

   //avoiding redundant integrals (i1>=i2 if same atom)
    int s1r = s1; if (n==m) s1r = s2;

    double Z1 = (double)atno[n]; double Z2 = (double)atno[m]; double Z3 = (double)atno[p];
    double A1 = coords[3*n+0]; double B1 = coords[3*n+1]; double C1 = coords[3*n+2];
    double A2 = coords[3*m+0]; double B2 = coords[3*m+1]; double C2 = coords[3*m+2];
    //if (n==m) { A2 = A1; B2 = B1; C2 = C1 + z0; }
    double A3 = coords[3*p+0]; double B3 = coords[3*p+1]; double C3 = coords[3*p+2];

    int nphi1 = nphi;

    get_3c_position(A1,B1,C1,A2,B2,C2,A3,B3,C3,rot);
    copy_to_all_gpu(nomp,9,rot,1);
    double coordn[9];
    coordn[0] = A1; coordn[1] = B1; coordn[2] = C1;
    coordn[3] = A2; coordn[4] = B2; coordn[5] = C2;
    coordn[6] = A3; coordn[7] = B3; coordn[8] = C3;

    int ncenters = count_unique(n,m,p);

    double z0 = sqrt(A2*A2+B2*B2+C2*C2)*0.5;
    if (m==n) z0 = sqrt(A3*A3+B3*B3+C3*C3)*0.5;

    if (prl>1) printf("  3c En. s1-4: %2i %2i %2i %2i (at: %i %i %i) z0: %8.5f \n",s1,s2,s3,s4,n,m,p,z0);

   #pragma omp parallel for schedule(dynamic,1) num_threads(nomp)
    for (int i2=s3;i2<s4;i2++)
    for (int i1=min(s1r,i2);i1<s2;i1++)
    {
      int tid = omp_get_thread_num();
      acc_set_device_num(tid,acc_device_nvidia);

      double* grid; double* gridm; double* wt;
      double* val = NULL;
      double sum = 0.;

      vector<double> basis1 = basis[i1];
      vector<double> basis2 = basis[i2];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zt1 = basis1[3];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zt2 = basis2[3];
      double norm1 = basis1[4]; double norm2 = basis2[4];

      int nbatch = 1;
     #if 0
      if (maxsteps>20) nbatch = 4;
      if (maxsteps>24) nbatch = 12;
      if (maxsteps>27) nbatch = 16;
     #endif

     //deal with long-range functions
      double cfn = cf*get_cfn(n1,l1,zt1,n2,l2,zt2,Z1,Z2);

      if (prl>0 && tid==0)
      {
        printf("\n GPU %i working on %2i-%2i (%2i) \n",tid,i1,i2,p);
        printf("  nlmz: %i %i %2i  %8.5f \n",n1,l1,m1,zt1);
        printf("  nlmz: %i %i %2i  %8.5f \n",n2,l2,m2,zt2);
      }

      for (int wb=0;wb<nbatch;wb++)
      {
        int gs = ss*gs0; if (ncenters<3) gs = ss*(gs0-nx3); int gs6 = 6*gs;
        if (maxsteps>1)
        {
          grid = new double[gs6]; gridm = new double[gs6]; wt = new double[gs];
          #pragma acc enter data create(grid[0:gs6],gridm[0:gs6],wt[0:gs])
        }

        if (ncenters==3)
          initialize_ps_coords_3c(-1,cfn,nmu,nnu,nphi1,0.,coordn,grid0,gridm0,wt0,rot,prl-tid);
        else
          initialize_ps_coords_batch(-1,wb,nbatch,z0,cfn,nmu,nnu,nphi1,0.,grid0,gridm0,wt0,prl-tid);

        if (maxsteps>1)
          half_split_ps(1,gs0,grid0,gridm0,wt0,z0,grid,gridm,wt,prl);
        else
        {
         //handle non-adaptive case
          val = new double[gs0];
          #pragma acc enter data create(val[0:gs0])
          evaluate_over_grid_quad_En_3c(quad_order,n1,l1,m1,zt1,norm1,n2,l2,m2,zt2,norm2,z0,A1,B1,C1,A2,B2,C2,A3,B3,C3,Z3,gs0,grid0,gridm0,wt0,rot,val);

          sum = 0.;
          #pragma acc parallel loop present(val[0:gs0]) reduction(+:sum)
          for (int j=0;j<gs0;j++)
            sum += val[j];

          #pragma acc exit data delete(val[0:gs0])
          delete [] val; val = NULL;
        }

       #if 0
        double* grid1 = new double[gs6]; //just for debugging
        #pragma acc enter data create(grid1[0:gs6])

        #pragma acc update self(grid[0:gs6])
        reorient_grid(-1,z0,gs,grid,grid1,rot);

        printf(" grid: \n");
        for (int k=0;k<gs;k++)
          printf("   %8.5f %8.5f %8.5f \n",grid[6*k+0],grid[6*k+1],grid[6*k+2]);

        #pragma acc update self(grid1[0:gs6])
        printf(" grid1: \n");
        for (int k=0;k<gs;k++)
          printf("   %8.5f %8.5f %8.5f \n",grid1[6*k+0],grid1[6*k+1],grid1[6*k+2]);

        #pragma acc exit data delete(grid1[0:gs6])
        delete [] grid1;
       #endif

       //adapt the grid specifically for e-n attraction
        int wsc = 0; int nzero = 0; int nconv = 0; double sump = 1000.; int gsp = 1000000;
        if (maxsteps>1)
        for (int ns=0;ns<maxsteps;ns++)
        {
          int ws = (wsc%3)+1;
          if (ws==4) { ws = 1; wsc = 0; }
          wsc++;

          if (val!=NULL) delete [] val;
          val = new double[gs];

          #pragma acc enter data create(val[0:gs])

          evaluate_over_grid_quad_En_3c(quad_order,n1,l1,m1,zt1,norm1,n2,l2,m2,zt2,norm2,z0,A1,B1,C1,A2,B2,C2,A3,B3,C3,Z3,gs,grid,gridm,wt,rot,val);

          double* gridn = NULL; double* gridmn = NULL; double* wtn = NULL;
          double asum = 0.; bool gup = 0;
          int gs1 = adaptive_split_ps_2ss_gpu(epsilon,ws,gs,z0,grid,gridm,wt,val,gup,asum,gridn,gridmn,wtn,prl);
          sum = asum;

          #pragma acc exit data delete(val[0:gs])

          swap_grids_gpu(gs,grid,gridm,wt,gridn,gridmn,wtn);
          gs = gs1; //gs6 = 6*gs;

          if (tid==0 && prl>0) printf("  step %2i  ws: %i  sum: %14.10f  gs: %7i \n",ns,ws,sum,gs);

          //if (fabs(sum-sump)<zero_tol && ws==1) nconv++;
          if (gs==gsp) nconv++; else nconv = 0;
          if (nconv>2) { if (tid==0 && prl>0) printf(" ** converged ** \n"); break; }
          if (fabs(sum)<zero_tol) nzero++; else nzero = 0;
          if (nzero>3) { if (tid==0 && prl>0) printf(" ** converged to zero ** \n"); break; }
          sump = sum;
          gsp = gs;

        } //loop ns, adapting grid
        if (tid==0 && prl>0) printf("\n");

        if (prl>0) printf("  final 3c(en). sum: %14.10f  tid: %i \n",sum,tid);

       //i1/i2->regular basis
       #pragma omp atomic
        En[i1*N+i2] += sum; 
        if (i1!=i2)
         #pragma omp atomic
          En[i2*N+i1] += sum;

        if (maxsteps>1)
        {
          #pragma acc exit data delete(grid[0:gs6],gridm[0:gs6],wt[0:gs])

          delete [] grid;
          delete [] gridm;
          delete [] wt;
        }
      } //nonzero integral

     //cleanup
      if (val!=NULL) delete [] val;

    } //loop i1,i2
  } //outer loop over i1,i2
  acc_set_device_num(0,acc_device_nvidia);


 //copy for symmetric indices
  if (C!=NULL)
  for (int n=0;n<natoms;n++)
  for (int m=n;m<natoms;m++)
  for (int p=0;p<natoms;p++)
  {
    int s1 = 0; if (n>0) s1 = n2i[n-1]; int s2 = n2i[n];
    int s3 = 0; if (m>0) s3 = n2i[m-1]; int s4 = n2i[m];
    int s5 = 0; if (p>0) s5 = a2i[p-1]; int s6 = a2i[p];

   //i1/i2 swap
    if (m==n)
    for (int i2=s3;i2<s4;i2++)
    for (int i1=i2+1;i1<s2;i1++)
    for (int i3=s5;i3<s6;i3++)
      C[i2*nna+i1*Naux+i3] = C[i1*nna+i2*Naux+i3];

    if (m!=n)
    for (int i2=s3;i2<s4;i2++)
    for (int i1=s1;i1<s2;i1++)
    for (int i3=s5;i3<s6;i3++)
      C[i2*nna+i1*Naux+i3] = C[i1*nna+i2*Naux+i3];
  }

  if (C!=NULL)
  if (prl>2)
  {
    printf("\n C(ps): \n");
    for (int i=0;i<Naux;i++)
    {
      //printf(" i: %i \n",i);
      for (int j=0;j<N;j++)
      for (int k=0;k<N;k++)
        printf("  %12.6f",C[j*nna+k*Naux+i]);
      printf("\n");
    }
  }

  for (int n=0;n<nomp;n++)
  {
    int tid = n; //omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);
    #pragma acc exit data delete(rot[0:9])
    #pragma acc exit data delete(grid0[0:6*gs0],gridm0[0:6*gs0],wt0[0:gs0])
  }
  acc_set_device_num(0,acc_device_nvidia);
  delete [] grid0; delete [] gridm0; delete [] wt0;

  //auto_crash();

  return;
}

void evaluate_over_grid_quad(int qo, bool do_coulomb, int type, int n1, int l1, int m1, double zt1, double norm1, double norm1v, int n2, int l2, int m2, double zt2, double norm2, 
                        double z0, int a12, int b12, double A1, double B1, double C1, double A2, double B2, double C2, double Z1, double Z2, bool do_rot, int gs, double* grid, double* gridm, double* wt, double* rot, double* val)
{
 //grid comes in in 1+2 order, but doesn't matter here

  int qo2 = 2*qo;
  int qos = qo*qo;
  double* Qx = new double[qo2];
  double* Qy = new double[qo2];
  double* Qz = new double[qo2];

  get_quad(qo,Qx);
  get_quad(qo,Qy);
  get_quad(qo,Qz);
  if (do_rot) //if z-aligned, can do 2d integral, otherwise 3d
    qos *= qo;

  double norm12 = norm1*norm2; if (do_coulomb) norm12 = norm1v*norm2;

  int gsq = qos*gs;
  int gsq6 = 6*gsq;

  double* gridq = new double[6*gsq];
  double* gridqr = new double[6*gsq];
  double* wtq = new double[gsq];
  double* valq = new double[gsq];

  double* grid1 = new double[gsq6];
  double* grid2 = new double[gsq6];

  #pragma acc enter data copyin(Qx[0:qo2],Qy[0:qo2],Qz[0:qo2])
  #pragma acc enter data create(gridq[0:gsq6],gridqr[0:gsq6],wtq[0:gsq],valq[0:gsq],grid1[0:gsq6],grid2[0:gsq6])

  #pragma acc update self(gridm[0:6*gs],wt[0:gs])
  //#pragma acc update self(val[0:gs])

  if (do_rot)
  {
    quad_grid_munuphi(-1,qo,qo,qo,0,z0,Qx,Qy,Qz,gs,0,gridm,gridq,wtq);
    reorient_grid(-1,z0,gsq,gridq,gridqr,rot);
    get_two_grids(gsq,grid1,grid2,gridqr,A1,B1,C1,A2,B2,C2);
  }
  else
  {
    quad_grid_munu(qo,qo,z0,Qx,Qy,gs,gridm,gridq,wtq);
    get_two_grids(gsq,grid1,grid2,gridq,A1,B1,C1,A2,B2,C2);
  }

  double* grida = grid1;
  double* gridb = grid2;

 #pragma acc parallel loop present(valq[0:gsq])
  for (int k=0;k<gsq;k++)
    valq[k] = norm12;

  if (type==1)
  {
   //overlap
    if (!do_coulomb)
      eval_shd(-1,gsq,grida,valq,n1,l1,m1,zt1);
    eval_shd(-1,gsq,gridb,valq,n2,l2,m2,zt2);

    if (do_coulomb)
    {
      eval_sh_3rd(gsq,grida,valq,n1,l1,m1);
      eval_inr_r12(-1,gsq,grida,valq,n1,l1,zt1);
    }
  } //S or 2c Coulomb

  if (type==2)
  {
    if (do_coulomb) { printf("\n ERROR: type==2 integrals cannot be Coulomb \n"); exit(-1); }
    if (Z1==0. || Z2==0.) { printf("\n ERROR: type==2 integrals cannot have Z=0. \n"); exit(-1); }

    //nuclear attraction integrals
    eval_shd(-1,gsq,grida,valq,n1,l1,m1,zt1);
    eval_shd(-1,gsq,gridb,valq,n2,l2,m2,zt2);

    double fn = 1.; if (a12==b12) fn = 0.5;
   #pragma acc parallel loop present(valq[0:gsq],grid1[0:gsq],grid2[0:gsq]) //note grid1/2 here
    for (int j=0;j<gsq;j++)
    {
      double Rn1 = grid1[6*j+3]; double Rn2 = grid2[6*j+3];
      double ne = -Z1/Rn1 - Z2/Rn2;
      valq[j] *= fn*ne;
    }

  } //En

  if (type==3)
  {
    //double* valqt = new double[gsq];
    //#pragma acc enter data create(valqt[0:gsq])

    //kinetic energy integrals
    eval_shd(-1,gsq,grida,valq,n1,l1,m1,zt1);
    eval_shd(-1,gsq,gridb,valq,n2,l2,m2,zt2);

  #if 0
   #pragma acc parallel loop present(valqt[0:gsq],valq[0:gsq])
    for (int j=0;j<gsq;j++)
      valqt[j] = valq[j];
  #endif

    eval_ked(-1,gsq,grida,valq,n1,l1,zt1);

  #if 0
    eval_ke(gsq,gridb,valqt,n2,l2,zt2);

   #pragma acc parallel loop present(valqt[0:gsq],valq[0:gsq])
    for (int j=0;j<gsq;j++)
      valqt[j] += valqt[j];

   #pragma acc parallel loop present(valq[0:gsq])
    for (int j=0;j<gsq;j++)
      valq[j] *= 0.5;

    #pragma acc exit data delete(valqt[0:gsq])
    delete [] valqt;
   #endif
  } //T


 #pragma acc parallel loop present(val[0:gs],valq[0:gsq],wtq[0:gsq])
  for (int j=0;j<gs;j++)
  {
    int i1 = j*qos;

    double v1 = 0.;
   #pragma acc loop reduction(+:v1)
    for (int k=0;k<qos;k++)
      v1 += valq[i1+k]*wtq[i1+k];

    val[j] = v1;
  }

  #pragma acc exit data delete(Qx[0:qo2],Qy[0:qo2],Qz[0:qo2])
  #pragma acc exit data delete(gridq[0:gsq6],gridqr[0:gsq6],wtq[0:gsq],valq[0:gsq],grid1[0:gsq6],grid2[0:gsq6])

  delete [] grid1;
  delete [] grid2;

  delete [] Qx;
  delete [] Qy;
  delete [] Qz;
  delete [] gridq;
  delete [] gridqr;
  delete [] wtq;
  delete [] valq;

  return;
}

void evaluate_over_grid(bool do_coulomb, int n1, int l1, int m1, double zt1, double norm1, double norm1v, int n2, int l2, int m2, double zt2, double norm2, 
                        double z0, double A1, double B1, double C1, double A2, double B2, double C2, int gs, double* grid, double* gridm, double* wt, double* val)
{
  return;
}

bool check_grid_for_batching(int nmu, int nnu, int nphi)
{
  int mult = nmu*nnu*nphi;
  if (mult%4>0)
  {
    printf("\n ERROR: cannot use batching \n");
    return 0;
  }

  return 1;
}

void integrate_STEnAC_2c(int natoms, int* atno, double* coords, vector<vector<double> > basis, vector<vector<double> > basis_aux, double epsilon, int nmu, int nnu, int nphi, bool do_coulomb, double* S, double* T, double* En, double* A, double* C, int prl)
{
  if (prl>0) printf("\n\n ---------- integrate all 2c using PS coordinates ---------- \n\n");

  //prl = 1;

 //this code only partially set up for more than 2 atoms.
  if (natoms>3) { printf("\n TESTING: natoms>3 in integrate_STEnAC_2c \n"); }

  int quad_order = read_int("QUAD");
  if (quad_order<=0) quad_order = QUAD_ORDER;

  int nomp_max = 1;
 #pragma omp parallel
  nomp_max = omp_get_num_threads();

  int ngpu = 0;
 #if USE_ACC
  ngpu = acc_get_num_devices(acc_device_nvidia);
  int nomp = ngpu;
 #else
  int nomp = nomp_max;
 #endif

  int gs_quad = nmu*nnu*nphi*quad_order*quad_order*quad_order;
  printf("  basis sizes: %2i  %2i  initial grid size: %7i \n",basis.size(),basis_aux.size(),gs_quad);
  if (basis.size()<1 || (basis_aux.size()<1 && do_coulomb))
  {
    return;
  }

  size_t gpumem = acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  printf("  gpu memory: %3zu GB \n",gpumem/1024/1024/1024);

  bool can_batch = check_grid_for_batching(nmu,nnu,nphi);

 //initialize the grid and other related quantities
  const int ss = 3;

  const double zero_tol = ZERO_TOL;
  double cf = 1.;

  vector<vector<double> > basis_phase1;
  if (do_coulomb) basis_phase1 = basis_aux;
  else basis_phase1 = basis;

  int N = basis_phase1.size();
  int Naux = basis_aux.size();

  int n2i[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis_phase1,n2i);

  int a2i[natoms];
  int imaxNa = get_imax_n2i(natoms,Naux,basis_aux,a2i);

  int gsmn = nmu*nnu;
  int gs0m = gsmn*nphi;
  double* grid0 = new double[6*gs0m];
  double* gridm0 = new double[6*gs0m];
  double* wt0 = new double[gs0m];
  double rot[9];

  for (int n=0;n<nomp;n++)
  {
    int tid = n; //omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);
    #pragma acc enter data create(rot[0:9])
    #pragma acc enter data create(grid0[0:6*gs0m],gridm0[0:6*gs0m],wt0[0:gs0m])
  }

  if (!do_coulomb)
  for (int j=0;j<N*N;j++)
    S[j] = En[j] = T[j] = 0.;
  if (do_coulomb)
  for (int j=0;j<Naux*Naux;j++)
    A[j] = 0.;


  //may want to evaluate T+N simultaneously

  int maxsteps = read_int("NSTEPS_PS"); //grid refinement steps
  if (maxsteps<1) maxsteps = 9;
  printf("  PS grid adapt steps: %2i \n",maxsteps);

 //2c integrals
  for (int n=0;n<natoms;n++)
  for (int m=n;m<natoms;m++)
  {
    int s1 = 0; if (n>0) s1 = n2i[n-1]; int s2 = n2i[n];
    int s3 = 0; if (m>0) s3 = n2i[m-1]; int s4 = n2i[m];

    double Z1 = atno[n]; double Z2 = atno[m];
    double A1 = coords[3*n+0]; double B1 = coords[3*n+1]; double C1 = coords[3*n+2];
    double A2 = coords[3*m+0]; double B2 = coords[3*m+1]; double C2 = coords[3*m+2];
   #if 0
   //always two atom
    int at1 = 0; int at2 = 1;
    double Z1 = atno[at1]; double Z2 = atno[at2];
    double A1 = coords[3*at1+0]; double B1 = coords[3*at1+1]; double C1 = coords[3*at1+2];
    double A2 = coords[3*at2+0]; double B2 = coords[3*at2+1]; double C2 = coords[3*at2+2];
   #endif

    bool do_rot = 0;
    if (m!=n)
    {
      do_rot = 1;
      get_2c_position(A1,B1,C1,A2,B2,C2,rot);
      copy_to_all_gpu(nomp,9,rot,1);
    }
    else
      A1 = B1 = C1 = A2 = B2 = C2 = 0.;

    double z0 = sqrt(A2*A2+B2*B2+C2*C2)*0.5; if (m==n) z0 = 1.;

   //CPMZ check this
    int nphi1 = nphi;
    if (m==n) nphi1 = -1;

   #pragma omp parallel for schedule(dynamic,1) num_threads(nomp)
    for (int i1=s1;i1<s2;i1++)
    for (int i2=s3;i2<s4;i2++)
    if (m!=n || i2>=i1)
    {
      int tid = omp_get_thread_num();
      acc_set_device_num(tid,acc_device_nvidia);

      double sum = 0.;

      vector<double> basis1 = basis_phase1[i1];
      vector<double> basis2 = basis_phase1[i2];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zt1 = basis1[3];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zt2 = basis2[3];
      double norm1 = basis1[4]; double norm2 = basis2[4];
      double norm1v = 4.*PI/(2.*l1+1)*norm1;

      //double idp = get_idp(l1,m1,l2,m2);
      double idp = 1.;
      if (m==n) idp = get_idp_one_atom(do_coulomb,l1,m1,l2,m2);

     //deal with long-range functions
      double cfn = cf*get_cfn(n1,l1,zt1,n2,l2,zt2,Z1,Z2);

      int nbatch = 1;
     #if 0
      if (maxsteps>20) nbatch = 4;
      if (maxsteps>24) nbatch = 24;
      if (maxsteps>27) nbatch = 16;
     #endif

      if (tid==0 && prl>0 && idp!=0.)
      {
        printf("\n\n GPU-%i  working on %2i-%2i \n",tid,i1,i2);
        printf("  nlmz: %i %i %2i  %8.5f   idp/pi: %5.3f \n",n1,l1,m1,zt1,idp/PI);
        printf("  nlmz: %i %i %2i  %8.5f   idp/pi: %5.3f \n",n2,l2,m2,zt2,idp/PI);
      }

      double phi_phase = 0.;
      if ((m1<0 || m2<0) && m==n)
      {
       //y is zero, so align with maximum value of spherical harmonic
       //m1==m2 for 2c integral to be nonzero
        phi_phase = PI/fabs(m1)/2.;
      }

     //adaptive grid, sampling over mu,nu or mu,nu,phi
      int wsm = 4; if (m==n) wsm = 3;

     //real work in here
      if (idp!=0.)
      for (int wb=0;wb<nbatch;wb++)
      {
        int gs0 = gsmn*abs(nphi1);
        int gs = ss*gs0; int gs6 = 6*gs;

        double* grid = new double[gs6];
        double* gridm = new double[gs6];
        double* wt = new double[gs];

        #pragma acc enter data create(grid[0:gs6],gridm[0:gs6],wt[0:gs])

        initialize_ps_coords_batch(-1,wb,nbatch,z0,cfn,nmu,nnu,nphi1,phi_phase,grid0,gridm0,wt0,prl-tid);
        half_split_ps(1,gs0,grid0,gridm0,wt0,z0,grid,gridm,wt,prl);

       //adapt the grid for overlap (or Coulomb)
        int wsc = 0; int nzero = 0; int nconv = 0; double sump = 1000.; int gsp = 1000000;
        for (int ns=0;ns<maxsteps;ns++)
        {
          int ws = (wsc%3)+1;
          if (ws==wsm) { ws = 1; wsc = 0; }
          wsc++;

          double* val = new double[gs];
          #pragma acc enter data create(val[0:gs])

         //integration work in here, except final sum
          //evaluate_over_grid(do_coulomb,n1,l1,m1,zt1,norm1,norm1v,n2,l2,m2,zt2,norm2,z0,A1,B1,C1,A2,B2,C2,gs,grid,gridm,wt,val);
         //new quadrature
          evaluate_over_grid_quad(quad_order,do_coulomb,1,n1,l1,m1,zt1,norm1,norm1v,n2,l2,m2,zt2,norm2,z0,n+1,m+1,A1,B1,C1,A2,B2,C2,Z1,Z2,do_rot,gs,grid,gridm,wt,rot,val);

          double* gridn = NULL; double* gridmn = NULL; double* wtn = NULL;
          double asum = 0.; bool gup = 0;
          int gs1 = adaptive_split_ps_2ss_gpu(epsilon,ws,gs,z0,grid,gridm,wt,val,gup,asum,gridn,gridmn,wtn,prl);
          sum = asum*idp;

          #pragma acc exit data delete(val[0:gs])
          delete [] val;

          swap_grids_gpu(gs,grid,gridm,wt,gridn,gridmn,wtn);
          gs = gs1; gs6 = 6*gs;

          if (tid==0 && prl>0) printf("  step %2i  ws: %i  sum: %14.10f  gs: %7i \n",ns,ws,sum,gs);

          //if (fabs(sum-sump)<zero_tol && ws==1) nconv++;
          if (gs==gsp) nconv++; else nconv = 0;
          if (fabs(sum)<zero_tol) nzero++; else nzero = 0;
          if (nzero>3) { if (tid==0 && prl>0) printf(" ** converged to zero ** \n"); break; }
          if (nconv>2) { if (tid==0 && prl>0) printf(" ** converged ** \n"); break; }
          sump = sum;
          gsp = gs;

        } //loop ns, adapting grid
        if (tid==0 && prl>0) printf("\n");

        if (!do_coulomb)
        {
         #pragma omp atomic
          S[i1*N+i2] += sum;
          if (i1!=i2)
           #pragma omp atomic
            S[i2*N+i1] += sum;
        }
        else
        {
         #pragma omp atomic
          A[i1*Naux+i2] += sum; 
          if (i1!=i2)
           #pragma omp atomic
            A[i2*Naux+i1] += sum;
        }

       ///////////////////// grid refinement and first 2c integrals done ////////////////////////////////////////

        //#pragma acc update self(grid[0:6*gs],gridm[0:6*gs],wt[0:gs])

       //reduce to "regular" part of the grid
        int gso = gs/ss;
        for (int j=0;j<gso;j++)
        {
          int j1 = 18*j;
          for (int k=0;k<6;k++)
          {
            grid[6*j+k] = grid[j1+k];
            gridm[6*j+k] = gridm[j1+k];
          }
          wt[j] = wt[ss*j];
        }
        gs = gso; gs6 = 6*gs;
        #pragma acc update device(grid[0:gs6],gridm[0:gs6],wt[0:gs])

        double* val = new double[gs];
        #pragma acc enter data create(val[0:gs])

        double vtn = 0.; double vtk = 0.;
        if (!do_coulomb)
        {
         //En
          evaluate_over_grid_quad(quad_order,do_coulomb,2,n1,l1,m1,zt1,norm1,norm1v,n2,l2,m2,zt2,norm2,z0,n+1,m+1,A1,B1,C1,A2,B2,C2,Z1,Z2,do_rot,gs,grid,gridm,wt,rot,val);
         #pragma acc parallel loop present(val[0:gs]) reduction(+:vtn)
          for (int j=0;j<gs;j++)
            vtn += val[j];
          vtn *= idp;

         //KE
          evaluate_over_grid_quad(quad_order,do_coulomb,3,n1,l1,m1,zt1,norm1,norm1v,n2,l2,m2,zt2,norm2,z0,n+1,m+1,A1,B1,C1,A2,B2,C2,Z1,Z2,do_rot,gs,grid,gridm,wt,rot,val);
          double vtk = 0.;
         #pragma acc parallel loop present(val[0:gs]) reduction(+:vtk)
          for (int j=0;j<gs;j++)
            vtk += val[j];
          vtk *= -0.5*idp;

         #pragma omp atomic
          En[i1*N+i2] += vtn; 
          if (i1!=i2)
           #pragma omp atomic
            En[i2*N+i1] += vtn;
         #pragma omp atomic
          T[i1*N+i2] += vtk;
          if (i1!=i2)
           #pragma omp atomic
            T[i2*N+i1] += vtk;
        } //En+T integrals

        if (prl>0 && nbatch==1 && tid==0) printf(" final(%2i-%2i) S: %11.8f En: %11.8f T: %11.8f (%4.1e) \n",i1,i2,sum,vtn,vtk,epsilon);
        if (prl>0 && tid==0 && nbatch>1) printf(" batch(%2i-%2i) S: %11.8f En: %11.8f T: %11.8f (%4.1e) \n",i1,i2,sum,vtn,vtk,epsilon);

        #pragma acc exit data delete(val[0:gs])
        #pragma acc exit data delete(grid[0:gs6],gridm[0:gs6],wt[0:gs])

        delete [] grid; delete [] gridm; delete [] wt;
        delete [] val;
      } //idp!=0.

    } //loop i1,i2

  } //outer loop over atom pairs

  for (int n=0;n<nomp;n++)
  {
    int tid = n; //omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);
    #pragma acc exit data delete(rot[0:9])
    #pragma acc exit data delete(grid0[0:6*gs0m],gridm0[0:6*gs0m],wt0[0:gs0m])
  }
  acc_set_device_num(0,acc_device_nvidia);
  delete [] grid0; delete [] gridm0; delete [] wt0;



  printf("\n");

  N = basis.size();
  imaxN = get_imax_n2i(natoms,N,basis,n2i);

  if (prl>1)
  {
    printf("  n2i:");
    for (int n=0;n<natoms;n++)
      printf(" %i",n2i[n]);
    printf("\n");
  }

  if (prl>1)
  {
    printf("\n En(2c): \n");
    print_square_fine(N,En);
  }


 /////////////////////////////////////////////////////////////////////////
  if (!do_coulomb)
  {
    do_3c_integrals_ps(epsilon,cf,nomp,ss,maxsteps,natoms,coords,atno,
      n2i,a2i,basis,basis_aux,nmu,nnu,nphi,En,NULL,prl);
  }
 /////////////////////////////////////////////////////////////////////////

 /////////////////////////////////////////////////////////////////////////
  if (do_coulomb)
  {
    do_3c_integrals_ps(epsilon,cf,nomp,ss,maxsteps,natoms,coords,atno,
      n2i,a2i,basis,basis_aux,nmu,nnu,nphi,NULL,C,prl);
  }
 /////////////////////////////////////////////////////////////////////////


 #if 0
  printf("\n S (1/4): \n");
  for (int n=N/2;n<N;n++)
  {
    for (int m=0;m<N/2;m++)
      printf(" %13.10f",S[n*N+m]);
    printf("\n");
  }

  printf("\n En (1/4): \n");
  for (int n=N/2;n<N;n++)
  {
    for (int m=0;m<N/2;m++)
      printf(" %13.10f",En[n*N+m]);
    printf("\n");
  }

  printf("\n T (1/4): \n");
  for (int n=N/2;n<N;n++)
  {
    for (int m=0;m<N/2;m++)
      printf(" %13.10f",T[n*N+m]);
    printf("\n");
  }
 #endif

  if (!do_coulomb && prl>0)
  {
    printf("\n S: \n");
    print_square_fine(N,S);
    printf("\n En: \n");
    print_square_fine(N,En);
    printf("\n T: \n");
    print_square_fine(N,T);
  }

  if (!do_coulomb && prl>-1)
  {
    printf("\n S diagonals accuracy: \n");
    for (int j=0;j<N;j++)
    {
      double v1 = log10(fabs(1.-S[j*N+j])+1.e-16);
      printf(" %5.2f",v1);
    }
    printf("\n");
  }

 //not really necessary, but might as well
  for (int j=0;j<N;j++)
    S[j*N+j] = 1.;

  if (do_coulomb && prl>0)
  {
    printf("\n A: \n");
    print_square_fine(Naux,A);

    printf("\n C: \n ");
    for (int k=0;k<Naux;k++)
      printf("%i%i%i ",(int)basis_aux[k][0],(int)basis_aux[k][1],(int)basis_aux[k][2]);
    printf("\n");
    for (int i=0;i<N;i++)
    for (int j=0;j<N;j++)
    {
      printf("%i%i/%i%i ",(int)basis[i][1],(int)basis[i][2],(int)basis[j][1],(int)basis[j][2]);
      int i1 = i*N+j;
      for (int k=0;k<Naux;k++)
        printf(" %12.8f",C[i1*Naux+k]);
      printf("\n");
    }
  }

  return;
}

void integrate_STEnAC_2c(int natoms, int* atno, float* coordsf, vector<vector<double> > basis, vector<vector<double> > basis_aux, double epsilon, int nmu, int nnu, int nphi, bool do_coulomb, double* S, double* T, double* En, double* A, double* C, int prl)
{
  double coords[3*natoms];
  for (int j=0;j<3*natoms;j++)
    coords[j] = coordsf[j];
  return integrate_STEnAC_2c(natoms,atno,coords,basis,basis_aux,epsilon,nmu,nnu,nphi,do_coulomb,S,T,En,A,C,prl);
}

bool integrate_ol_4c(int natoms, int* atno, double* coords, vector<vector<double> > basis, double epsilon, int nmu, int nnu, int nphi, double* g, int prl)
{
  if (prl>0) printf("\n\n ---------- integrate 4c OL using PS coordinates ---------- \n\n");

 //this code only partially set up for more than 2 atoms. For now, exit if natoms>2
  if (natoms>2) { printf("\n ERROR: natoms must be two in integrate_ol_4c \n"); return 0; }
  if (natoms>3) { printf("\n ERROR: natoms>3 in integrate_ol_4c \n"); return 0; }

  int nomp_max = 1;
 #pragma omp parallel
  nomp_max = omp_get_num_threads();

  int ngpu = 0;
 #if USE_ACC
  ngpu = acc_get_num_devices(acc_device_nvidia);
  int nomp = ngpu;
 #else
  int nomp = nomp_max;
 #endif

  printf("  basis size: %2i  ngpu/nomp: %2i \n",basis.size(),nomp);
  if (basis.size()<1)
    return 0;

  size_t gpumem = acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  printf("  gpu memory: %3zu GB \n",gpumem/1024/1024/1024);

 //initialize the grid and other related quantities
  const int ss = 3;

  if (coords[0]!=0. || coords[1]!=0. || coords[3]!=0. || coords[4]!=0.)
  { printf("\n\n ERROR: integrate_ol_4c atoms must be aligned on z axis \n"); return 0; }

 //CPMZ needs updating
  double zdist = coords[2]-coords[5];
  double z0 = zdist*0.5;
  if (zdist<0)
  { printf("\n\n ERROR: second atom must be along -z \n"); exit(-1); }

  double cf = 1.;
  int maxsteps = read_int("NSTEPS_PS"); //grid refinement steps
  if (maxsteps<1) maxsteps = 9;

  do_4c_integrals_ps(epsilon,cf,nomp,ss,maxsteps,z0,natoms,coords,atno,
      basis,nmu,nnu,nphi,g,prl);

  return 1;
}

vector<int> para_division(int ss, int gs, int nbatch)
{
 //batches must be in multiples of ss
  vector<int> batches;
  batches.push_back(0);

  int bs = gs/nbatch;
  if (bs%ss!=0)
    bs += ss-bs%ss;

  for (int j=0;j<nbatch;j++)
  {
    int i1 = (j+1)*bs;
    if (j+1==nbatch) i1 = gs;
    if (i1>=gs) { batches.push_back(gs); break; }

    batches.push_back(i1);
  }

  return batches;
}


//test ftns assume a=1
double test_1s_integrate_v2(int ss, int maxsteps, double epsilon, int gs1, double* grid1, double* gridm1, double* wt1, int prl)
{
  double z0 = 1.; //parameter

  int nbatch = 4;
  if (epsilon<1.e-10) nbatch = 64;
  if (epsilon<1.e-12) nbatch = 256;
  vector<int> batches = para_division(ss,gs1,nbatch);
  nbatch = batches.size()-1;

  printf("\n batches: ");
  for (int j=0;j<nbatch;j++)
    printf(" %3i",batches[j]);
  printf("\n");


  printf("\n\n 1s_integrate_v2 with nbatch: %2i  gs1: %6i \n",nbatch,gs1);
 

  double gsum = 0.;
  for (int ns=0;ns<nbatch;ns++)
  {
    int s0 = batches[ns];
    int s06 = 6*s0;
    int bs = batches[ns+1]-s0;
    int bs6 = 6*bs;

    printf("  working on batch %2i of initial size %3i. ws: %i \n",ns+1,bs);
    printf("    s0/6: %5i %5i \n",s0,s06);

    double* grid = new double[bs6];
    double* gridm = new double[bs6];
    double* wt = new double[bs];
    for (int j=0;j<bs6;j++)
    {
      grid[j] = grid1[s06+j];
      gridm[j] = gridm1[s06+j];
    }
    for (int j=0;j<bs;j++)
      wt[j] = wt1[s0+j];


    double tisum = 0.; //inactive terms
    double sum = 0.; //total sum
    for (int n=0;n<maxsteps;n++)
    {
      int ws = ((n+1)%3)+1;
      //ws = 1;

      double eps1 = epsilon; if (n+1==maxsteps) eps1 = 0.; //stops grid split at last step
      printf("    iter %2i with eps: %4.1e ws: %i \n",n+1,eps1,ws);

      double* val = new double[bs];
      for (int j=0;j<bs;j++)
      {
       //two 1s AOs
        double x1 = grid[6*j+0]; double y1 = grid[6*j+1]; double z1 = grid[6*j+2]-z0;
        double r1 = sqrt(x1*x1+y1*y1+z1*z1);
        double x2 = grid[6*j+0]; double y2 = grid[6*j+1]; double z2 = grid[6*j+2]+z0;
        double r2 = sqrt(x2*x2+y2*y2+z2*z2);

        double wt1 = wt[j];
        val[j] = 32.*exp(-2.*r1)*exp(-2.*r2)*wt1;
      }

      double* gridn = NULL;
      double* gridmn = NULL;
      double* wtn = NULL;

      double isum = 0.; double asum = 0.; bool gup = 0; double maxd;
      bs = adaptive_split_ps(eps1,ws,ss,bs,z0,grid,gridm,wt,val,gup,isum,asum,maxd,gridn,gridmn,wtn,prl);
      tisum += isum;
      sum = asum+tisum;

      printf("    asum: %11.8f tisum: %11.8f sum: %11.8f \n",asum,tisum,sum);

     //replace prior grid
      delete [] grid; delete [] gridm; delete [] wt;
      grid = gridn; gridm = gridmn; wt = wtn;

      delete [] val;

      if (bs<1 || n+1==maxsteps) { grid = NULL; gridm = NULL; wt = NULL; break; }
    }

    if (bs>0) printf("  this batch did not converge (bs: %4i) \n",bs);

    if (grid!=NULL) delete [] grid;
    if (gridm!=NULL) delete [] gridm;
    if (wt!=NULL) delete [] wt;

    gsum += sum;
    printf("\n");
  }
  printf("\n final sum: %11.8f  (%4.1e) \n",gsum,epsilon);

  return gsum;
}

double test_1s_integrate(int ss, int maxsteps, double epsilon, int gs, double*& grid, double*& gridm, double*& wt, int prl)
{
  double sum = 0.;
  return sum;
}


void test_prosph()
{
  double epsilon = read_float("EPSILON");
  if (epsilon<=0) epsilon = 1.e-12;

  int maxsteps = read_int("NSTEPS_PS");
  if (maxsteps<1) maxsteps = 4;
  int prl = 1+2;

  double z0 = 1.;

  int nmu = 10;
  int nnu = 10;
  int nphi = read_int("PHI");
  if (nphi<1) nphi = 4;

  if (epsilon<1.e-12) { nmu = 24; nnu = 20; }

  int gs0 = nmu*nnu*nphi;
  double* grid0 = new double[6*gs0];
  double* gridm0 = new double[6*gs0];
  double* wt0 = new double[gs0];

  initialize_ps_coords_2c(-1,z0,1.,nmu,nnu,nphi,0.,grid0,gridm0,wt0,prl);

  if (prl>0)
  {
    printf("\n grid0: \n");
    for (int j=0;j<gs0;j++)
      printf("  %8.5f %8.5f %8.5f  wt: %8.5f \n",grid0[6*j+0],grid0[6*j+1],grid0[6*j+2],wt0[j]);
      //printf("  %8.5f %8.5f %8.5f  wt: %8.5f \n",gridm0[6*j+0],gridm0[6*j+1],gridm0[6*j+2],wt0[j]);
    printf("  grid0 size: %3i \n",gs0);
  }

  int ss = 3;

 //first grid
  int gs = gs0*ss;
  double* grid = new double[6*gs];
  double* gridm = new double[6*gs];
  double* wt = new double[gs];

  if (ss==9) //octree grid
    full_split_ps(gs0,grid0,gridm0,wt0,z0,gs,grid,gridm,wt,prl+1);
  else
    half_split_ps(1,gs0,grid0,gridm0,wt0,z0,grid,gridm,wt,prl+1);

  printf("  initial grid size: %5i \n",gs);

  delete [] grid0;
  delete [] gridm0;
  delete [] wt0;

 //adaptive integration
  //double sum = test_1s_integrate_v2(ss,maxsteps,epsilon,gs,grid,gridm,wt,prl);
  double sum = test_1s_integrate(ss,maxsteps,epsilon,gs,grid,gridm,wt,prl);


 //cleanup
  delete [] grid;
  delete [] gridm;
  delete [] wt;

  return;
}

