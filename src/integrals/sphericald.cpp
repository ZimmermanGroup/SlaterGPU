#include "spherical.h"


#define SS_V1 0
#define SS_V3 1

// if integrating in PS coordinates, make sure nphi > 2l
 // this is based on tests on dihydrogen with higher angular momentum


// start spherical * exp //
//  double precision     //


#pragma acc routine seq
int facto(int N)
{
  int v = 1;
  //#pragma acc reduction(*:v)
  for (int i=2;i<=N;i++)
    v *= i;

  return v;
}

void get_1s_expd(int tid, int gs, double* grid, double* val, double zeta1)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double ezr = exp(-zeta1*r);
    val[i] *= ezr;
  }
  return;
}

void get_2s_expd(int tid, int gs, double* grid, double* val, double zeta1)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double rezr = r*exp(-zeta1*r);
    val[i] *= rezr;
  }
  return;
}

void get_3s_expd(int tid, int gs, double* grid, double* val, double zeta1)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double r2ezr = r*r*exp(-zeta1*r);
    val[i] *= r2ezr;
  }
  return;
}

void get_4s_expd(int tid, int gs, double* grid, double* val, double zeta1)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double r2 = r*r;
    double r3ezr = r2*r*exp(-zeta1*r);
    val[i] *= r3ezr;
  }
  return;
}

void get_5s_expd(int tid, int gs, double* grid, double* val, double zeta1)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double r4 = pow(r,4.);
    double r4ezr = r4*exp(-zeta1*r);
    val[i] *= r4ezr;
  }
  return;
}

void get_6s_expd(int tid, int gs, double* grid, double* val, double zeta1)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double r5 = pow(r,5.);
    double r5ezr = r5*exp(-zeta1*r);
    val[i] *= r5ezr;
  }
  return;
}

void get_7s_expd(int tid, int gs, double* grid, double* val, double zeta1)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double r6 = pow(r,6.);
    double r6ezr = r6*exp(-zeta1*r);
    val[i] *= r6ezr;
  }
  return;
}

void get_8s_expd(int tid, int gs, double* grid, double* val, double zeta1)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double r7 = pow(r,7.);
    double r7ezr = r7*exp(-zeta1*r);
    val[i] *= r7ezr;
  }
  return;
}

void get_9s_expd(int tid, int gs, double* grid, double* val, double zeta1)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double r8 = pow(r,8.);
    double r8ezr = r8*exp(-zeta1*r);
    val[i] *= r8ezr;
  }
  return;
}

void get_10s_expd(int tid, int gs, double* grid, double* val, double zeta1)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double r9 = pow(r,9.);
    double r9ezr = r9*exp(-zeta1*r);
    val[i] *= r9ezr;
  }
  return;
}

void get_11s_expd(int tid, int gs, double* grid, double* val, double zeta1)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double r10 = pow(r,10.);
    double r10ezr = r10*exp(-zeta1*r);
    val[i] *= r10ezr;
  }
  return;
}

void get_12s_expd(int tid, int gs, double* grid, double* val, double zeta1)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double r11 = pow(r,11.);
    double r11ezr = r11*exp(-zeta1*r);
    val[i] *= r11ezr;
  }
  return;
}

void get_2pxd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= x*ezr;
  }
}

void get_2pyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= y*ezr;
  }
}

void get_2pzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= z*ezr;
  }
}

void get_3pxd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double r = grid[6*i+3];
    double rezr = r*exp(-zeta*r);
    val[i] *= x*rezr;
  }
}

void get_3pyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double rezr = r*exp(-zeta*r);
    val[i] *= y*rezr;
  }
}

void get_3pzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double rezr = r*exp(-zeta*r);
    val[i] *= z*rezr;
  }
}

void get_4pxd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double r = grid[6*i+3];
    double rezr = r*r*exp(-zeta*r);
    val[i] *= x*rezr;
  }
}

void get_4pyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double rezr = r*r*exp(-zeta*r);
    val[i] *= y*rezr;
  }
}

void get_4pzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double rezr = r*r*exp(-zeta*r);
    val[i] *= z*rezr;
  }
}

void get_5pxd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double r = grid[6*i+3];
    double rezr = r*r*r*exp(-zeta*r);
    val[i] *= x*rezr;
  }
}

void get_5pyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double rezr = r*r*r*exp(-zeta*r);
    val[i] *= y*rezr;
  }
}

void get_5pzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double rezr = r*r*r*exp(-zeta*r);
    val[i] *= z*rezr;
  }
}

void get_6pxd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double r = grid[6*i+3];
    double r4 = pow(r,4.);
    double rezr = r4*exp(-zeta*r);
    val[i] *= x*rezr;
  }
}

void get_6pyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double r4 = pow(r,4.);
    double rezr = r4*exp(-zeta*r);
    val[i] *= y*rezr;
  }
}

void get_6pzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r4 = pow(r,4.);
    double rezr = r4*exp(-zeta*r);
    val[i] *= z*rezr;
  }
}

void get_7pxd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double r = grid[6*i+3];
    double r5 = pow(r,5.);
    double rezr = r5*exp(-zeta*r);
    val[i] *= x*rezr;
  }
}

void get_7pyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double r5 = pow(r,5.);
    double rezr = r5*exp(-zeta*r);
    val[i] *= y*rezr;
  }
}

void get_7pzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r5 = pow(r,5.);
    double rezr = r5*exp(-zeta*r);
    val[i] *= z*rezr;
  }
}

void get_8pxd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double r = grid[6*i+3];
    double r6 = pow(r,6.);
    double rezr = r6*exp(-zeta*r);
    val[i] *= x*rezr;
  }
}

void get_8pyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double r6 = pow(r,6.);
    double rezr = r6*exp(-zeta*r);
    val[i] *= y*rezr;
  }
}

void get_8pzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r6 = pow(r,6.);
    double rezr = r6*exp(-zeta*r);
    val[i] *= z*rezr;
  }
}

void get_9pxd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double r = grid[6*i+3];
    double r7 = pow(r,7.);
    double rezr = r7*exp(-zeta*r);
    val[i] *= x*rezr;
  }
}

void get_9pyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double r7 = pow(r,7.);
    double rezr = r7*exp(-zeta*r);
    val[i] *= y*rezr;
  }
}

void get_9pzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r7 = pow(r,7.);
    double rezr = r7*exp(-zeta*r);
    val[i] *= z*rezr;
  }
}

void get_3dxyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= x*y*ezr;
  }
}

void get_3dyzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= y*z*ezr;
  }
}

void get_3dz2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= (2.*z*z-x*x-y*y)*ezr;

   //ambiguous for which r to use
    //val[i] *= 3.*z*z-r*r; 
  }
}

void get_3dxzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= x*z*ezr;
  }
}

void get_3dx2y2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= (x*x-y*y)*ezr;
  }
}

//for Cartesian expansion of d functions
void get_3dxxd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= x*x*ezr;
  }
}

void get_3dyyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= y*y*ezr;
  }
}

void get_3dzzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= z*z*ezr;
  }
}

void get_4dxyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = r*exp(-zeta*r);
    val[i] *= x*y*ezr;
  }
}

void get_4dyzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*exp(-zeta*r);
    val[i] *= y*z*ezr;
  }
}

void get_4dz2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*exp(-zeta*r);
    val[i] *= (2.*z*z-x*x-y*y)*ezr;

   //ambiguous for which r to use
    //val[i] *= 3.*z*z-r*r; 
  }
}

void get_4dxzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*exp(-zeta*r);
    val[i] *= x*z*ezr;
  }
}

void get_4dx2y2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = r*exp(-zeta*r);
    val[i] *= (x*x-y*y)*ezr;
  }
}

void get_5dxyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = r*r*exp(-zeta*r);
    val[i] *= x*y*ezr;
  }
}

void get_5dyzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*r*exp(-zeta*r);
    val[i] *= y*z*ezr;
  }
}

void get_5dz2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*r*exp(-zeta*r);
    val[i] *= (2.*z*z-x*x-y*y)*ezr;

   //ambiguous for which r to use
    //val[i] *= 3.*z*z-r*r; 
  }
}

void get_5dxzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*r*exp(-zeta*r);
    val[i] *= x*z*ezr;
  }
}

void get_5dx2y2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = r*r*exp(-zeta*r);
    val[i] *= (x*x-y*y)*ezr;
  }
}

void get_6dxyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double r3 = pow(r,3.);
    double ezr = r3*exp(-zeta*r);
    val[i] *= x*y*ezr;
  }
}

void get_6dyzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r3 = pow(r,3.);
    double ezr = r3*exp(-zeta*r);
    val[i] *= y*z*ezr;
  }
}

void get_6dz2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r3 = pow(r,3.);
    double ezr = r3*exp(-zeta*r);
    val[i] *= (2.*z*z-x*x-y*y)*ezr;

   //ambiguous for which r to use
    //val[i] *= 3.*z*z-r*r; 
  }
}

void get_6dxzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r3 = pow(r,3.);
    double ezr = r3*exp(-zeta*r);
    val[i] *= x*z*ezr;
  }
}

void get_6dx2y2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double r3 = pow(r,3.);
    double ezr = r3*exp(-zeta*r);
    val[i] *= (x*x-y*y)*ezr;
  }
}

void get_7dxyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double r4 = pow(r,4.);
    double ezr = r4*exp(-zeta*r);
    val[i] *= x*y*ezr;
  }
}

void get_7dyzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r4 = pow(r,4.);
    double ezr = r4*exp(-zeta*r);
    val[i] *= y*z*ezr;
  }
}

void get_7dz2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r4 = pow(r,4.);
    double ezr = r4*exp(-zeta*r);
    val[i] *= (2.*z*z-x*x-y*y)*ezr;

   //ambiguous for which r to use
    //val[i] *= 3.*z*z-r*r; 
  }
}

void get_7dxzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r4 = pow(r,4.);
    double ezr = r4*exp(-zeta*r);
    val[i] *= x*z*ezr;
  }
}

void get_7dx2y2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double r4 = pow(r,4.);
    double ezr = r4*exp(-zeta*r);
    val[i] *= (x*x-y*y)*ezr;
  }
}

void get_4fm3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= (3.*x*x-y*y)*y*ezr;
  }
}

void get_4fm2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= x*y*z*ezr;
  }
}

void get_4fm1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= (4.*z*z-x*x-y*y)*y*ezr;
  }
}

void get_4f0d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= (2.*z*z-3.*x*x-3.*y*y)*z*ezr;
  }
}

void get_4fp1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= (4.*z*z-x*x-y*y)*x*ezr;
  }
}

void get_4fp2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= (x*x-y*y)*z*ezr;
  }
}

void get_4fp3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= (x*x-3.*y*y)*x*ezr;
  }
}

void get_5fm3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = r*exp(-zeta*r);
    val[i] *= (3.*x*x-y*y)*y*ezr;
  }
}

void get_5fm2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*exp(-zeta*r);
    val[i] *= x*y*z*ezr;
  }
}

void get_5fm1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*exp(-zeta*r);
    val[i] *= (4.*z*z-x*x-y*y)*y*ezr;
  }
}

void get_5f0d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*exp(-zeta*r);
    val[i] *= (2.*z*z-3.*x*x-3.*y*y)*z*ezr;
  }
}

void get_5fp1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*exp(-zeta*r);
    val[i] *= (4.*z*z-x*x-y*y)*x*ezr;
  }
}

void get_5fp2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*exp(-zeta*r);
    val[i] *= (x*x-y*y)*z*ezr;
  }
}

void get_5fp3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = r*exp(-zeta*r);
    val[i] *= (x*x-3.*y*y)*x*ezr;
  }
}

void get_6fm3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = r*r*exp(-zeta*r);
    val[i] *= (3.*x*x-y*y)*y*ezr;
  }
}

void get_6fm2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*r*exp(-zeta*r);
    val[i] *= x*y*z*ezr;
  }
}

void get_6fm1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*r*exp(-zeta*r);
    val[i] *= (4.*z*z-x*x-y*y)*y*ezr;
  }
}

void get_6f0d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*r*exp(-zeta*r);
    val[i] *= (2.*z*z-3.*x*x-3.*y*y)*z*ezr;
  }
}

void get_6fp1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*r*exp(-zeta*r);
    val[i] *= (4.*z*z-x*x-y*y)*x*ezr;
  }
}

void get_6fp2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*r*exp(-zeta*r);
    val[i] *= (x*x-y*y)*z*ezr;
  }
}

void get_6fp3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = r*r*r*exp(-zeta*r);
    val[i] *= (x*x-3.*y*y)*x*ezr;
  }
}

void get_7fm3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = r*r*r*exp(-zeta*r);
    val[i] *= (3.*x*x-y*y)*y*ezr;
  }
}

void get_7fm2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*r*r*exp(-zeta*r);
    val[i] *= x*y*z*ezr;
  }
}

void get_7fm1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*r*r*exp(-zeta*r);
    val[i] *= (4.*z*z-x*x-y*y)*y*ezr;
  }
}

void get_7f0d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*r*r*exp(-zeta*r);
    val[i] *= (2.*z*z-3.*x*x-3.*y*y)*z*ezr;
  }
}

void get_7fp1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*r*r*exp(-zeta*r);
    val[i] *= (4.*z*z-x*x-y*y)*x*ezr;
  }
}

void get_7fp2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = r*r*r*exp(-zeta*r);
    val[i] *= (x*x-y*y)*z*ezr;
  }
}

void get_7fp3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = r*r*r*exp(-zeta*r);
    val[i] *= (x*x-3.*y*y)*x*ezr;
  }
}

void get_5gm4d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    //double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= x*y * (x*x - y*y) * ezr;
  }
}

void get_5gm3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= y*z * (3.*x*x - y*y) * ezr;
  }
}

void get_5gm2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= x*y * (6.*z*z - x*x - y*y) * ezr;
  }
}

void get_5gm1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= y*z * (4.*z*z - 3.*x*x - 3.*y*y) * ezr;
  }
}

void get_5g0d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double z2 = z*z;
    double x2py2 = x*x+y*y;
    double r2 = r*r;
    //val[i] *= (35.*z2*z2 - 30.*z2*r2 + 3.*r2*r2) * ezr; //equal via algebra
    val[i] *= (35.*z2*z2 - 30.*z2*r2 + 3.*r2*r2) * ezr;
    //val[i] *= (3.*x2py2*x2py2 - 24.*x2py2*z2 + 8.*z2*z2) * ezr;
  }
}

void get_5gp1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= x*z * (4.*z*z - 3.*x*x - 3.*y*y) * ezr;
  }
}

void get_5gp2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    val[i] *= (x2 - y2) * (6.*z*z - x2 - y2) * ezr;
  }
}

void get_5gp3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    val[i] *= x*z * (x*x - 3.*y*y) * ezr;
  }
}

void get_5gp4d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    //double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    val[i] *= (x2 * (x2 - 3.*y2) - y2 * (3.*x2 - y2)) * ezr;
  }
}


void get_6hm5d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    //double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;

    val[i] *= -(5.*x2*x2 - 10.*x2*y2 + y2*y2)*y * ezr;
    //val[i] *= (5.*x2*x2 - 10.*x2*y2 + y2*y2)*y * ezr;
  }
}

void get_6hm4d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;

    val[i] *= (y2 - x2)*x*y*z * ezr;
    //val[i] *= (y2 - x2)*x*y*z * ezr;
  }
}

void get_6hm3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= y * (r2*(3.*x2-y2) + 9.*z2*(y2-3.*x2)) * ezr;
    //val[i] *= (x2*x2 - 10.*x2*y2 + 5.*y2*y2)*x * ezr;
  }
}

void get_6hm2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= x*y*z*(3.*z2-r2) * ezr;
    //val[i] *= y * (y2-3.*x2)*(x2 + y2 - 8.*z2) * ezr;
  }
}

void get_6hm1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= y * (14.*r2*z2 - 21.*z2*z2 - r2*r2) * ezr;
    //val[i] *= z * (x2*x2 - 6.*x2*y2 + y2*y2) * ezr;
  }
}

void get_6h0d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= z*(15.*r2*r2 - 70.*r2*z2 + 63.*z2*z2) * ezr;
    //val[i] *= (x2 + y2 - 2.*z2)*x*y*z * ezr;
  }
}

void get_6hp1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= x*(-r2*r2 + 14.*r2*z2 - 21.*z2*z2) * ezr;
    //val[i] *= x * (x2-3.*y2)*(x2 + y2 - 8.*z2) * ezr;
  }
}

void get_6hp2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= -z * (r2 - 3.*z2)*(x-y)*(x+y) * ezr;
    //val[i] *= y * (x2*x2 + y2*y2 - 12.*y2*z2 + 8.*z2*z2 + 2.*x2 * (y2-6.*z2)) * ezr;
  }
}

void get_6hp3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= x * (r2*(x2-3.*y2) + 9.*z2*(3.*y2-x2)) * ezr;
    //val[i] *= (x2 - y2) * (x2 + y2 - 2.*z2) * z * ezr;
  }
}

void get_6hp4d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    //double z2 = z*z;

    val[i] *= z * (x2*x2 - 6.*x2*y2 + y2*y2) * ezr;
    //val[i] *= z * (15.*x2*x2 + 15.*y2*y2 - 40.*y2*z2 + 8.*z2*z2 + 10.*x2 * (3.*y2 - 4.*z2)) * ezr;
  }
}

void get_6hp5d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    //double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    //double z2 = z*z;

    val[i] *= x * (10.*x2*y2 - x2*x2 - 5.*y2*y2) * ezr;
    //val[i] *= x * (x2*x2 + y2*y2 - 12.*y2*z2 + 8.*z2*z2 + 2.*x2* (y2 - 6.*z2)) * ezr;
  }
}

void get_7im6d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  { 
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    double xy = x*y;
 
    val[i] *= xy*(3.*x2*x2 - 10.*x2*y2 + 3.*y2*y2) * ezr;
  }
}

void get_7im5d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    double xz = x*z;

    val[i] *= xz * (x2*x2 - 10.*x2*y2 + 5.*y2*y2) * ezr;
  }
}

void get_7im4d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= (x*(x - y)*y*(x + y)*(r2 - 11.*z2)) * ezr;
  }
}

void get_7im3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= x*((x2 - 3.*y2)*z*(3.*r2 - 11.*z2)) * ezr;
  }
}

void get_7im2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= x*y * (r2*r2 - 18.*r2*z2 + 33.*z2*z2) * ezr;
  }
}

void get_7im1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= x*z * (5.*r2*r2 - 30.*r2*z2 + 33.*z2*z2) * ezr;
  }
}

void get_7i0d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= (-5.*r2*r2*r2 + 105.*r2*r2*z2 - 315.*r2*z2*z2 + 231.*z2*z2*z2) * ezr;
  }
}

void get_7ip1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= y*z * (5.*r2*r2 - 30.*r2*z2 + 33.*z2*z2) * ezr;
  }
}

void get_7ip2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= (x - y)*(x + y) * (r2*r2 - 18.*r2*z2 + 33.*z2*z2) * ezr;
  }
}

void get_7ip3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= y * (-3.*x2 + y2)*z*(-3.*r2 + 11.*z2) * ezr;
  }
}

void get_7ip4d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= (x2*x2 - 6.*x2*y2 + y2*y2)*(r2 - 11.*z2) * ezr;
  }
}

void get_7ip5d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    double yz = y*z;

    val[i] *= yz * (5.*x2*x2 - 10.*x2*y2 + y2*y2) * ezr;
  }
}

void get_7ip6d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;

    val[i] *= (x2*x2*x2 - 15.*x2*x2*y2 + 15.*x2*y2*y2 - y2*y2*y2) * ezr;
  }
}

void get_8jm7d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;

    val[i] *= x*(x2*x2*x2 - 21.*x2*x2*y2 + 35.*x2*y2*y2 - 7.*y2*y2*y2) * ezr;
  }
}

void get_8jm6d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x; 
    double y2 = y*y; 
    double xyz = x*y*z;

    val[i] *= xyz * (3.*x2*x2 - 10.*x2*y2 + 3.*y2*y2) * ezr; 
  }
}

void get_8jm5d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= x*(x2*x2 - 10.*x2*y2 + 5.*y2*y2) * (r2 - 13.*z2) * ezr;
  }
}

void get_8jm4d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= x * (x - y) * y * (x + y) * z * (-3.*r2 + 13.*z2) * ezr;
  }
}

void get_8jm3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= x*(x2 - 3.*y2) * (3.*r2*r2 - 66.*r2*z2 + 143.*z2*z2) * ezr;
  }
}

void get_8jm2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= x * y * z * (15.*r2*r2 - 110.*r2*z2 + 143.*z2*z2) * ezr;
  }
}

void get_8jm1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= x * (5.*r2*r2*r2 - 135.*r2*r2*z2 + 495.*r2*z2*z2 - 429.*z2*z2*z2) * ezr;
  }
}

void get_8j0d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= z*(-35.*r2*r2*r2 + 315.*r2*r2*z2 - 693.*r2*z2*z2 + 429.*z2*z2*z2) * ezr;
  }
}

void get_8jp1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= y * (5.*r2*r2*r2 - 135.*r2*r2*z2 + 495.*r2*z2*z2 - 429.*z2*z2*z2) * ezr;
  }
}

void get_8jp2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= (x - y) * (x + y) * z * (15.*r2*r2 - 110.*r2*z2 + 143.*z2*z2) * ezr;
  }
}

void get_8jp3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= y * (-3.*x2 + y2) * (3.*r2*r2 - 66.*r2*z2 + 143.*z2*z2) * ezr;
  }
}

void get_8jp4d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x; 
    double y2 = y*y; 
    double z2 = z*z; 
    double r2 = r*r; 

    val[i] *= (x2*x2 - 6.*x2*y2 + y2*y2) * z * (-3.*r2 + 13.*z2) * ezr; 
  }
}

void get_8jp5d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= y*(5.*x2*x2 - 10.*x2*y2 + y2*y2) * (r2 - 13.*z2) * ezr;     
  }
}

void get_8jp6d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;

    val[i] *= z * (x2*x2*x2 - 15.*x2*x2*y2 + 15.*x2*y2*y2 - y2*y2*y2) * ezr;
  }
}

void get_8jp7d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double x2 = x*x;
    double y2 = y*y;

    val[i] *= y * (-7.*x2*x2*x2 + 35.*x2*x2*y2 - 21.*x2*y2*y2 + y2*y2*y2) * ezr;
  }
}

// end spherical * exp //


// spherical w/r factors //

void get_px_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double r = grid[6*i+3];
    val[i] *= x/r;
  }
}

void get_py_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    val[i] *= y/r;
  }
}

void get_pz_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    val[i] *= z/r;
  }
}

void get_dxy_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    val[i] *= x*y/r/r;
  }
}

void get_dyz_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    val[i] *= y*z/r/r;
  }
}

void get_dz2_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    val[i] *= 3.*z*z/r/r-1.;
  }
}

void get_dxz_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    val[i] *= x*z/r/r;
  }
}

void get_dx2y2_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    val[i] *= (x*x-y*y)/r/r;
  }
}

//Cartesian d functions
void get_dxx_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double r = grid[6*i+3];
    val[i] *= x*x/r/r;
  }
}

void get_dyy_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    val[i] *= y*y/r/r;
  }
}

void get_dzz_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    val[i] *= z*z/r/r;
  }
}

void get_fm3_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double or3 = pow(r,-3.);
    val[i] *= (3.*x*x-y*y)*y;
    val[i] *= or3;
  }
}

void get_fm2_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or3 = pow(r,-3.);
    val[i] *= x*y*z*or3;
  }
}

void get_fm1_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or3 = pow(r,-3.);
    val[i] *= (4.*z*z-x*x-y*y)*y*or3;
  }
}

void get_f0_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or3 = pow(r,-3.);
    val[i] *= (2.*z*z-3.*x*x-3.*y*y)*z*or3;
  }
}

void get_fp1_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or3 = pow(r,-3.);
    val[i] *= (4.*z*z-x*x-y*y)*x*or3;
  }
}

void get_fp2_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or3 = pow(r,-3.);
    val[i] *= (x*x-y*y)*z*or3;
  }
}

void get_fp3_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double or3 = pow(r,-3.);
    val[i] *= (x*x-3.*y*y)*x*or3;
  }
}

void get_gm4_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double or4 = pow(r,-4.);
    val[i] *= x*y * (x*x - y*y) * or4;
  }
}

void get_gm3_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or4 = pow(r,-4.);
    val[i] *= y*z * (3.*x*x - y*y) * or4;
  }
}

void get_gm2_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or4 = pow(r,-4.);
    val[i] *= x*y * (6.*z*z - x*x - y*y) * or4;
  }
}

void get_gm1_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or4 = pow(r,-4.);
    val[i] *= y*z * (4.*z*z - 3.*x*x - 3.*y*y) * or4;
  }
}

void get_g0_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or4 = pow(r,-4.);
    double z2 = z*z;
    double x2py2 = x*x+y*y;
    double r2 = r*r;
    val[i] *= (35.*z2*z2 - 30.*z2*r2) * or4 + 3.;
  }
}

void get_gp1_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or4 = pow(r,-4.);
    val[i] *= x*z * (4.*z*z - 3.*x*x - 3.*y*y) * or4;
  }
}

void get_gp2_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or4 = pow(r,-4.);
    double x2 = x*x;
    double y2 = y*y;
    val[i] *= (x2 - y2) * (6.*z*z - x2 - y2) * or4;
  }
}

void get_gp3_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or4 = pow(r,-4.);
    val[i] *= x*z * (x*x - 3.*y*y) * or4;
  }
}

void get_gp4_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    //double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or4 = pow(r,-4.);
    double x2 = x*x;
    double y2 = y*y;
    val[i] *= (x2 * (x2 - 3.*y2) - y2 * (3.*x2 - y2)) * or4;
  }
}

void get_hm5_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    //double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or5 = pow(r,-5.);
    double x2 = x*x;
    double y2 = y*y;

    val[i] *= -(5.*x2*x2 - 10.*x2*y2 + y2*y2)*y*or5;
  }
}

void get_hm4_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or5 = pow(r,-5.);
    double x2 = x*x;
    double y2 = y*y;

    val[i] *= (y2 - x2)*x*y*z*or5;
  }
}

void get_hm3_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or5 = pow(r,-5.);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= y * (r2*(3.*x2-y2) + 9.*z2*(y2-3.*x2))*or5;
  }
}

void get_hm2_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or5 = pow(r,-5.);
    //double x2 = x*x;
    //double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= x*y*z*(3.*z2-r2)*or5;
  }
}

void get_hm1_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    //double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or5 = pow(r,-5.);
    //double x2 = x*x;
    //double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= y * (14.*r2*z2 - 21.*z2*z2 - r2*r2)*or5;
  }
}

void get_h0_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or5 = pow(r,-5.);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= z*(15.*r2*r2 - 70.*r2*z2 + 63.*z2*z2)*or5;
  }
}

void get_hp1_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    //double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or5 = pow(r,-5.);
    //double x2 = x*x;
    //double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= x*(-r2*r2 + 14.*r2*z2 - 21.*z2*z2)*or5;
  }
}

void get_hp2_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or5 = pow(r,-5.);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= -z * (r2 - 3.*z2)*(x-y)*(x+y)*or5;
  }
}

void get_hp3_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or5 = pow(r,-5.);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= x * (r2*(x2-3.*y2) + 9.*z2*(3.*y2-x2))*or5;
  }
}

void get_hp4_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or5 = pow(r,-5.);
    double x2 = x*x;
    double y2 = y*y;
    //double z2 = z*z;

    val[i] *= z * (x2*x2 - 6.*x2*y2 + y2*y2)*or5;
  }
}

void get_hp5_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    //double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or5 = pow(r,-5.);
    double x2 = x*x;
    double y2 = y*y;
    //double z2 = z*z;

    val[i] *= x * (10.*x2*y2 - x2*x2 - 5.*y2*y2)*or5;
  }
}

void get_im6_3rd(int tid, int gs, double* grid, double* val) 
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double or6 = pow(r,-6.);
    double x2 = x*x;
    double y2 = y*y; 
    double xy = x*y;

    val[i] *= xy*(3.*x2*x2 - 10.*x2*y2 + 3.*y2*y2)*or6;
  }
}   

void get_im5_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or6 = pow(r,-6.);
    double x2 = x*x;
    double y2 = y*y;
    double xz = x*z;

    val[i] *= xz * (x2*x2 - 10.*x2*y2 + 5.*y2*y2)*or6;
  }
}        

void get_im4_3rd(int tid, int gs, double* grid, double* val) 
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or6 = pow(r,-6.);
    double z2 = z*z;
    double r2 = r*r;
  
    val[i] *= (x*(x - y)*y*(x + y)*(r2 - 11.*z2))*or6;
  }
}     

void get_im3_3rd(int tid, int gs, double* grid, double* val) 
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or6 = pow(r,-6.);
    double x2 = x*x; 
    double y2 = y*y; 
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= (x*(x2 - 3.*y2)*z*(3.*r2 - 11.*z2))*or6;
  }
}     

void get_im2_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or6 = pow(r,-6.);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= x*y * (r2*r2 - 18.*r2*z2 + 33.*z2*z2)*or6;
  }
}

void get_im1_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or6 = pow(r,-6.);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= x*z * (5.*r2*r2 - 30.*r2*z2 + 33.*z2*z2)*or6;
  }
}

void get_i0_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or6 = pow(r,-6.);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= -5. + ((105.*r2*r2*z2 - 315.*r2*z2*z2 + 231.*z2*z2*z2)*or6);
  }
}

void get_ip1_3rd(int tid, int gs, double* grid, double* val) 
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or6 = pow(r,-6.);
    double z2 = z*z; 
    double r2 = r*r;

    val[i] *= y*z * (5.*r2*r2 - 30.*r2*z2 + 33.*z2*z2)*or6;
  }
}

void get_ip2_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or6 = pow(r,-6.);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= (x - y)*(x + y) * (r2*r2 - 18.*r2*z2 + 33.*z2*z2)*or6;
  }
}

void get_ip3_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or6 = pow(r,-6.);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;
    double yz = y*z;

    val[i] *= yz * (-3.*x2 + y2)*(-3.*r2 + 11.*z2)*or6;
  }
}

void get_ip4_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or6 = pow(r,-6.);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= (x2*x2 - 6.*x2*y2 + y2*y2)*(r2 - 11.*z2)*or6;
  }
}

void get_ip5_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or6 = pow(r,-6.);
    double x2 = x*x;
    double y2 = y*y;
    double yz = y*z;

    val[i] *= yz * (5.*x2*x2 - 10.*x2*y2 + y2*y2)*or6;
  }
}

void get_ip6_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double or6 = pow(r,-6.);
    double x2 = x*x;
    double y2 = y*y;

    val[i] *= (x2*x2*x2 - 15.*x2*x2*y2 + 15.*x2*y2*y2 - y2*y2*y2)*or6;
  }
}

void get_jm7_3rd(int tid, int gs, double* grid, double* val) 
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double or7 = pow(r,-7.);
    double x2 = x*x; 
    double y2 = y*y; 

    val[i] *= x*(x2*x2*x2 - 21.*x2*x2*y2 + 35.*x2*y2*y2 - 7.*y2*y2*y2)*or7;
  }
}

void get_jm6_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or7 = pow(r,-7.);
    double x2 = x*x;
    double y2 = y*y;
    double xyz = x*y*z;

    val[i] *= xyz * (3.*x2*x2 - 10.*x2*y2 + 3.*y2*y2)*or7;
  }
}

void get_jm5_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or7 = pow(r,-7.);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;
      
    val[i] *= x*(x2*x2 - 10.*x2*y2 + 5.*y2*y2) * (r2 - 13.*z2)*or7;
  }
}


void get_jm4_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or7 = pow(r,-7.);
    double z2 = z*z;
    double r2 = r*r;
      
    val[i] *= x * (x - y) * y * (x + y) * z * (-3.*r2 + 13.*z2)*or7;
  }
}

void get_jm3_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or7 = pow(r,-7.);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= x*(x2 - 3.*y2) * (3.*r2*r2 - 66.*r2*z2 + 143.*z2*z2)*or7;
  }
}

void get_jm2_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or7 = pow(r,-7.);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= x * y * z * (15.*r2*r2 - 110.*r2*z2 + 143.*z2*z2)*or7;
  }
}

void get_jm1_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or7 = pow(r,-7.);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= x * (5.*r2*r2*r2 - 135.*r2*r2*z2 + 495.*r2*z2*z2 - 429.*z2*z2*z2)*or7;
  }
}

void get_j0_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or7 = pow(r,-7.);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= z*(-35.*r2*r2*r2 + 315.*r2*r2*z2 - 693.*r2*z2*z2 + 429.*z2*z2*z2)*or7;
  }
}

void get_jp1_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or7 = pow(r,-7.);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= y * (5.*r2*r2*r2 - 135.*r2*r2*z2 + 495.*r2*z2*z2 - 429.*z2*z2*z2)*or7;
  }
}

void get_jp2_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or7 = pow(r,-7.);
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= (x - y) * (x + y) * z * (15.*r2*r2 - 110.*r2*z2 + 143.*z2*z2)*or7;
  }
}

void get_jp3_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or7 = pow(r,-7.);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= y * (-3.*x2 + y2) * (3.*r2*r2 - 66.*r2*z2 + 143.*z2*z2)*or7;
  }
}

void get_jp4_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or7 = pow(r,-7.);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= (x2*x2 - 6.*x2*y2 + y2*y2) * z * (-3.*r2 + 13.*z2)*or7;
  }
}

void get_jp5_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or7 = pow(r,-7.);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double r2 = r*r;

    val[i] *= y*(5.*x2*x2 - 10.*x2*y2 + y2*y2) * (r2 - 13.*z2)*or7;
  }
}

void get_jp6_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double or7 = pow(r,-7.);
    double x2 = x*x;
    double y2 = y*y;

    val[i] *= z * (x2*x2*x2 - 15.*x2*x2*y2 + 15.*x2*y2*y2 - y2*y2*y2)*or7;
  }
}

void get_jp7_3rd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double or7 = pow(r,-7.);
    double x2 = x*x;
    double y2 = y*y;

    val[i] *= y * (-7.*x2*x2*x2 + 35.*x2*x2*y2 - 21.*x2*y2*y2 + y2*y2*y2)*or7;
  }
}

// spherical w/r factors //


void eval_sh_3rd(int tid, int gs, double* grid, double* val, int n1, int l1, int m1)
{
  if (l1==0) //1s
  {
    return;
  }
  else if (l1==1) //p
  {
    if (m1==1)
      return get_px_3rd(tid,gs,grid,val);
    else if (m1==-1)
      return get_py_3rd(tid,gs,grid,val);
    else
      return get_pz_3rd(tid,gs,grid,val);
  }
  else if (l1==2) //d
  {
   #if CART_D
    if (m1==0)
      return get_dxx_3rd(tid,gs,grid,val);
    else if (m1==1)
      return get_dyy_3rd(tid,gs,grid,val);
    else if (m1==2)
      return get_dzz_3rd(tid,gs,grid,val);
    else if (m1==3)
      return get_dxy_3rd(tid,gs,grid,val);
    else if (m1==4)
      return get_dxz_3rd(tid,gs,grid,val);
    else if (m1==5)
      return get_dyz_3rd(tid,gs,grid,val);
   #else
    if (m1==-2)
      return get_dxy_3rd(tid,gs,grid,val);
    else if (m1==-1)
      return get_dyz_3rd(tid,gs,grid,val);
    else if (m1==0)
      return get_dz2_3rd(tid,gs,grid,val);
    else if (m1==1)
      return get_dxz_3rd(tid,gs,grid,val);
    else
      return get_dx2y2_3rd(tid,gs,grid,val);
   #endif
  }
  else if (l1==3) //f
  {
    if (m1==-3)
      return get_fm3_3rd(tid,gs,grid,val);
    else if (m1==-2)
      return get_fm2_3rd(tid,gs,grid,val);
    else if (m1==-1)
      return get_fm1_3rd(tid,gs,grid,val);
    else if (m1== 0)
      return get_f0_3rd(tid,gs,grid,val);
    else if (m1== 1)
      return get_fp1_3rd(tid,gs,grid,val);
    else if (m1== 2)
      return get_fp2_3rd(tid,gs,grid,val);
    else if (m1== 3)
      return get_fp3_3rd(tid,gs,grid,val);
  }
  else if (l1==4) //g
  {
    if (m1==-4)
      return get_gm4_3rd(tid,gs,grid,val);
    else if (m1==-3)
      return get_gm3_3rd(tid,gs,grid,val);
    else if (m1==-2)
      return get_gm2_3rd(tid,gs,grid,val);
    else if (m1==-1)
      return get_gm1_3rd(tid,gs,grid,val);
    else if (m1== 0)
      return get_g0_3rd(tid,gs,grid,val);
    else if (m1== 1)
      return get_gp1_3rd(tid,gs,grid,val);
    else if (m1== 2)
      return get_gp2_3rd(tid,gs,grid,val);
    else if (m1== 3)
      return get_gp3_3rd(tid,gs,grid,val);
    else if (m1== 4)
      return get_gp4_3rd(tid,gs,grid,val);
  }
  else if (l1==5) //h
  {
    if (m1==-5)
      return get_hm5_3rd(tid,gs,grid,val);
    else if (m1==-4)
      return get_hm4_3rd(tid,gs,grid,val);
    else if (m1==-3)
      return get_hm3_3rd(tid,gs,grid,val);
    else if (m1==-2)
      return get_hm2_3rd(tid,gs,grid,val);
    else if (m1==-1)
      return get_hm1_3rd(tid,gs,grid,val);
    else if (m1== 0)
      return get_h0_3rd(tid,gs,grid,val);
    else if (m1== 1)
      return get_hp1_3rd(tid,gs,grid,val);
    else if (m1== 2)
      return get_hp2_3rd(tid,gs,grid,val);
    else if (m1== 3)
      return get_hp3_3rd(tid,gs,grid,val);
    else if (m1== 4)
      return get_hp4_3rd(tid,gs,grid,val);
    else if (m1== 5)
      return get_hp5_3rd(tid,gs,grid,val);
  }
  else if (l1==6) //i
  {
    if (m1==-6)
      return get_im6_3rd(tid,gs,grid,val);
    else if (m1==-5)
      return get_im5_3rd(tid,gs,grid,val);
    else if (m1==-4)
      return get_im4_3rd(tid,gs,grid,val);
    else if (m1==-3)
      return get_im3_3rd(tid,gs,grid,val);
    else if (m1==-2)
      return get_im2_3rd(tid,gs,grid,val);
    else if (m1==-1)
      return get_im1_3rd(tid,gs,grid,val);
    else if (m1== 0)
      return get_i0_3rd(tid,gs,grid,val);
    else if (m1== 1)
      return get_ip1_3rd(tid,gs,grid,val);
    else if (m1== 2)
      return get_ip2_3rd(tid,gs,grid,val);
    else if (m1== 3)
      return get_ip3_3rd(tid,gs,grid,val);
    else if (m1== 4)
      return get_ip4_3rd(tid,gs,grid,val);
    else if (m1== 5)
      return get_ip5_3rd(tid,gs,grid,val);
    else if (m1== 6)
      return get_ip6_3rd(tid,gs,grid,val);
  }
  else if (l1==7) //j
  {
    if (m1==-7)
      return get_jm7_3rd(tid,gs,grid,val);
    else if (m1==-6)
      return get_jm6_3rd(tid,gs,grid,val);
    else if (m1==-5)
      return get_jm5_3rd(tid,gs,grid,val);
    else if (m1==-4)
      return get_jm4_3rd(tid,gs,grid,val);
    else if (m1==-3)
      return get_jm3_3rd(tid,gs,grid,val);
    else if (m1==-2)
      return get_jm2_3rd(tid,gs,grid,val);
    else if (m1==-1)
      return get_jm1_3rd(tid,gs,grid,val);
    else if (m1== 0)
      return get_j0_3rd(tid,gs,grid,val);
    else if (m1== 1)
      return get_jp1_3rd(tid,gs,grid,val);
    else if (m1== 2)
      return get_jp2_3rd(tid,gs,grid,val);
    else if (m1== 3)
      return get_jp3_3rd(tid,gs,grid,val);
    else if (m1== 4)
      return get_jp4_3rd(tid,gs,grid,val);
    else if (m1== 5)
      return get_jp5_3rd(tid,gs,grid,val);
    else if (m1== 6)
      return get_jp6_3rd(tid,gs,grid,val);
    else if (m1== 7)
      return get_jp7_3rd(tid,gs,grid,val);
  }

  return;
}

void eval_sh_3rd(int gs, double* grid, double* val, int n1, int l1, int m1)
{
  return eval_sh_3rd(0,gs,grid,val,n1,l1,m1);
}

void eval_shd(int tid, int gs, double* grid, double* val, int n1, int l1, int m1, double zeta1)
{
 //1s - 12s
 //2p - 9p
 //3-7d,4-7f,5g,6h

  if (l1==0)
  {
    if (n1==1)
      get_1s_expd(tid,gs,grid,val,zeta1);
    else if (n1==2)
      get_2s_expd(tid,gs,grid,val,zeta1);
    else if (n1==3)
      get_3s_expd(tid,gs,grid,val,zeta1);
    else if (n1==4)
      get_4s_expd(tid,gs,grid,val,zeta1);
    else if (n1==5)
      get_5s_expd(tid,gs,grid,val,zeta1);
    else if (n1==6)
      get_6s_expd(tid,gs,grid,val,zeta1);
    else if (n1==7)
      get_7s_expd(tid,gs,grid,val,zeta1);
    else if (n1==8)
      get_8s_expd(tid,gs,grid,val,zeta1);
    else if (n1==9)
      get_9s_expd(tid,gs,grid,val,zeta1);
    else if (n1==10)
      get_10s_expd(tid,gs,grid,val,zeta1);
    else if (n1==11)
      get_11s_expd(tid,gs,grid,val,zeta1);
    else if (n1==12)
      get_12s_expd(tid,gs,grid,val,zeta1);
  }
  else if (l1>0)
  {
    if (l1==1) //p
    {
      if (n1==2)
      {
        if (m1==1) //2p
          get_2pxd(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_2pyd(tid,gs,grid,val,zeta1);
        else
          get_2pzd(tid,gs,grid,val,zeta1);
      }
      else if (n1==3)
      {
        if (m1==1) //3p
          get_3pxd(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_3pyd(tid,gs,grid,val,zeta1);
        else
          get_3pzd(tid,gs,grid,val,zeta1);
      }
      else if (n1==4)
      {
        if (m1==1) //4p
          get_4pxd(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_4pyd(tid,gs,grid,val,zeta1);
        else
          get_4pzd(tid,gs,grid,val,zeta1);
      }
      else if (n1==5)
      {
        if (m1==1) //5p
          get_5pxd(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_5pyd(tid,gs,grid,val,zeta1);
        else
          get_5pzd(tid,gs,grid,val,zeta1);
      }
      else if (n1==6)
      {
        if (m1==1) //6p
          get_6pxd(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_6pyd(tid,gs,grid,val,zeta1);
        else
          get_6pzd(tid,gs,grid,val,zeta1);
      }
      else if (n1==7)
      {
        if (m1==1) //7p
          get_7pxd(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_7pyd(tid,gs,grid,val,zeta1);
        else
          get_7pzd(tid,gs,grid,val,zeta1);
      }
      else if (n1==8)
      {
        if (m1==1) //8p
          get_8pxd(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_8pyd(tid,gs,grid,val,zeta1);
        else
          get_8pzd(tid,gs,grid,val,zeta1);
      }
      else if (n1==9)
      {
        if (m1==1) //9p
          get_9pxd(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_9pyd(tid,gs,grid,val,zeta1);
        else
          get_9pzd(tid,gs,grid,val,zeta1);
      }
    }
    else if (l1==2) //d
    {
      if (n1==3)
      {
       #if CART_D
        if (m1==0)
          get_3dxxd(tid,gs,grid,val,zeta1);
        else if (m1==1)
          get_3dyyd(tid,gs,grid,val,zeta1);
        else if (m1==2)
          get_3dzzd(tid,gs,grid,val,zeta1);
        else if (m1==3)
          get_3dxyd(tid,gs,grid,val,zeta1);
        else if (m1==4)
          get_3dxzd(tid,gs,grid,val,zeta1);
        else
          get_3dyzd(tid,gs,grid,val,zeta1);
       #else
        if (m1==-2)
          get_3dxyd(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_3dyzd(tid,gs,grid,val,zeta1);
        else if (m1== 0)
          get_3dz2d(tid,gs,grid,val,zeta1);
        else if (m1== 1)
          get_3dxzd(tid,gs,grid,val,zeta1);
        else if (m1== 2)
          get_3dx2y2d(tid,gs,grid,val,zeta1);
       #endif
      }
      else if (n1==4)
      {
        if (m1==-2)
          get_4dxyd(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_4dyzd(tid,gs,grid,val,zeta1);
        else if (m1== 0)
          get_4dz2d(tid,gs,grid,val,zeta1);
        else if (m1== 1)
          get_4dxzd(tid,gs,grid,val,zeta1);
        else if (m1== 2)
          get_4dx2y2d(tid,gs,grid,val,zeta1);
      }
      else if (n1==5)
      {
        if (m1==-2)
          get_5dxyd(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_5dyzd(tid,gs,grid,val,zeta1);
        else if (m1== 0)
          get_5dz2d(tid,gs,grid,val,zeta1);
        else if (m1== 1)
          get_5dxzd(tid,gs,grid,val,zeta1);
        else if (m1== 2)
          get_5dx2y2d(tid,gs,grid,val,zeta1);
      }
      else if (n1==6)
      {
        if (m1==-2)
          get_6dxyd(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_6dyzd(tid,gs,grid,val,zeta1);
        else if (m1== 0)
          get_6dz2d(tid,gs,grid,val,zeta1);
        else if (m1== 1)
          get_6dxzd(tid,gs,grid,val,zeta1);
        else if (m1== 2)
          get_6dx2y2d(tid,gs,grid,val,zeta1);
      }
      else if (n1==7)
      {
        if (m1==-2)
          get_7dxyd(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_7dyzd(tid,gs,grid,val,zeta1);
        else if (m1== 0)
          get_7dz2d(tid,gs,grid,val,zeta1);
        else if (m1== 1)
          get_7dxzd(tid,gs,grid,val,zeta1);
        else if (m1== 2)
          get_7dx2y2d(tid,gs,grid,val,zeta1);
      }
    }
    else if (l1==3) //f
    {
      if (n1==4)
      {
        if (m1==-3)
          get_4fm3d(tid,gs,grid,val,zeta1);
        else if (m1==-2)
          get_4fm2d(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_4fm1d(tid,gs,grid,val,zeta1);
        else if (m1== 0)
          get_4f0d(tid,gs,grid,val,zeta1);
        else if (m1== 1)
          get_4fp1d(tid,gs,grid,val,zeta1);
        else if (m1== 2)
          get_4fp2d(tid,gs,grid,val,zeta1);
        else if (m1== 3)
          get_4fp3d(tid,gs,grid,val,zeta1);
      }
      else if (n1==5)
      {
        if (m1==-3)
          get_5fm3d(tid,gs,grid,val,zeta1);
        else if (m1==-2)
          get_5fm2d(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_5fm1d(tid,gs,grid,val,zeta1);
        else if (m1== 0)
          get_5f0d(tid,gs,grid,val,zeta1);
        else if (m1== 1)
          get_5fp1d(tid,gs,grid,val,zeta1);
        else if (m1== 2)
          get_5fp2d(tid,gs,grid,val,zeta1);
        else if (m1== 3)
          get_5fp3d(tid,gs,grid,val,zeta1);
      }
      else if (n1==6)
      {
        if (m1==-3)
          get_6fm3d(tid,gs,grid,val,zeta1);
        else if (m1==-2)
          get_6fm2d(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_6fm1d(tid,gs,grid,val,zeta1);
        else if (m1== 0)
          get_6f0d(tid,gs,grid,val,zeta1);
        else if (m1== 1)
          get_6fp1d(tid,gs,grid,val,zeta1);
        else if (m1== 2)
          get_6fp2d(tid,gs,grid,val,zeta1);
        else if (m1== 3)
          get_6fp3d(tid,gs,grid,val,zeta1);
      }
      else if (n1==7)
      {
        if (m1==-3)
          get_7fm3d(tid,gs,grid,val,zeta1);
        else if (m1==-2)
          get_7fm2d(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_7fm1d(tid,gs,grid,val,zeta1);
        else if (m1== 0)
          get_7f0d(tid,gs,grid,val,zeta1);
        else if (m1== 1)
          get_7fp1d(tid,gs,grid,val,zeta1);
        else if (m1== 2)
          get_7fp2d(tid,gs,grid,val,zeta1);
        else if (m1== 3)
          get_7fp3d(tid,gs,grid,val,zeta1);
      }
      else if (n1>7)
      {
        printf(" ERROR: 8f not available \n");
      }
    }
    else if (l1==4) //g
    {
      if (n1==6)
      {
        printf(" ERROR: 6g not available \n");
        exit(1);
      }

      if (m1==-4)
        get_5gm4d(tid,gs,grid,val,zeta1);
      else if (m1==-3)
        get_5gm3d(tid,gs,grid,val,zeta1);
      else if (m1==-2)
        get_5gm2d(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        get_5gm1d(tid,gs,grid,val,zeta1);
      else if (m1== 0)
        get_5g0d(tid,gs,grid,val,zeta1);
      else if (m1== 1)
        get_5gp1d(tid,gs,grid,val,zeta1);
      else if (m1== 2)
        get_5gp2d(tid,gs,grid,val,zeta1);
      else if (m1== 3)
        get_5gp3d(tid,gs,grid,val,zeta1);
      else if (m1== 4)
        get_5gp4d(tid,gs,grid,val,zeta1);
    }
    else if (l1==5) //h
    {
      if (m1==-5)
        get_6hm5d(tid,gs,grid,val,zeta1);
      else if (m1==-4)
        get_6hm4d(tid,gs,grid,val,zeta1);
      else if (m1==-3)
        get_6hm3d(tid,gs,grid,val,zeta1);
      else if (m1==-2)
        get_6hm2d(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        get_6hm1d(tid,gs,grid,val,zeta1);
      else if (m1== 0)
        get_6h0d(tid,gs,grid,val,zeta1);
      else if (m1== 1)
        get_6hp1d(tid,gs,grid,val,zeta1);
      else if (m1== 2)
        get_6hp2d(tid,gs,grid,val,zeta1);
      else if (m1== 3)
        get_6hp3d(tid,gs,grid,val,zeta1);
      else if (m1== 4)
        get_6hp4d(tid,gs,grid,val,zeta1);
      else if (m1== 5)
        get_6hp5d(tid,gs,grid,val,zeta1);
    }
    else if (l1==6) //i
    {
      if (m1==-6)
        get_7im6d(tid,gs,grid,val,zeta1);
      else if (m1==-5)
        get_7im5d(tid,gs,grid,val,zeta1);
      else if (m1==-4)
        get_7im4d(tid,gs,grid,val,zeta1);
      else if (m1==-3)
        get_7im3d(tid,gs,grid,val,zeta1);
      else if (m1==-2)
        get_7im2d(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        get_7im1d(tid,gs,grid,val,zeta1);
      else if (m1== 0)
        get_7i0d(tid,gs,grid,val,zeta1);
      else if (m1== 1)
        get_7ip1d(tid,gs,grid,val,zeta1);
      else if (m1== 2)
 	get_7ip2d(tid,gs,grid,val,zeta1);
      else if (m1== 3)
        get_7ip3d(tid,gs,grid,val,zeta1);
      else if (m1== 4)
        get_7ip4d(tid,gs,grid,val,zeta1);
      else if (m1== 5)
        get_7ip5d(tid,gs,grid,val,zeta1);
      else if (m1== 6)
        get_7ip6d(tid,gs,grid,val,zeta1);
    }
    else if (l1==7) //j
    {
      if (m1==-7)
        get_8jm7d(tid,gs,grid,val,zeta1);
      else if (m1==-6)
        get_8jm6d(tid,gs,grid,val,zeta1);
      else if (m1==-5)
        get_8jm5d(tid,gs,grid,val,zeta1);
      else if (m1==-4)
        get_8jm4d(tid,gs,grid,val,zeta1);
      else if (m1==-3)
        get_8jm3d(tid,gs,grid,val,zeta1);
      else if (m1==-2)
        get_8jm2d(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        get_8jm1d(tid,gs,grid,val,zeta1);
      else if (m1== 0)
        get_8j0d(tid,gs,grid,val,zeta1);
      else if (m1== 1)
        get_8jp1d(tid,gs,grid,val,zeta1);
      else if (m1== 2)
 	get_8jp2d(tid,gs,grid,val,zeta1);
      else if (m1== 3)
        get_8jp3d(tid,gs,grid,val,zeta1);
      else if (m1== 4)
        get_8jp4d(tid,gs,grid,val,zeta1);
      else if (m1== 5)
        get_8jp5d(tid,gs,grid,val,zeta1);
      else if (m1== 6)
        get_8jp6d(tid,gs,grid,val,zeta1);
      else if (m1== 7)
        get_8jp7d(tid,gs,grid,val,zeta1);
    }

  } //if l1>0

  return;
}

//Split Gaussian-Slater basis set
void eval_sgsd(int tid, int gs1, int gs2, double* grid, double* val, int n1, int l1, int m1, double zeta1, double Rc)
{
  //printf("  eval_sgsd gs12: %4i %4i \n",gs1,gs2);

  int gs = gs2;
  int gs6 = 6*gs;
  //int nml = n1-l1-1; //eval_sh_3rd has factor of r^-l
  int nm1 = n1-1;
  double zeta2 = 2.*zeta1*Rc;
  double norm2 = exp(zeta1*Rc*Rc);

  #pragma acc parallel loop present(grid[0:gs6],val[0:gs])
  for (int j=0;j<gs1;j++)
  {
    double r = grid[6*j+3];
    double rp = pow(r,nm1);

    val[j] *= rp*exp(-zeta1*r*r);
  }

  #pragma acc parallel loop present(grid[0:gs6],val[0:gs])
  for (int j=gs1;j<gs2;j++)
  {
    double r = grid[6*j+3];
    double rp = pow(r,nm1);

    val[j] *= norm2*rp*exp(-zeta2*r);
  }

  eval_sh_3rd(tid,gs,grid,val,n1,l1,m1);

  return;
}

void eval_sgsd(int tid, int gs1, int gs2, float* gridf, double* val, int n1, int l1, int m1, double zeta1, double Rc)
{
  int gs6 = 6*gs2;
  double* grid = new double[gs6];
  #pragma acc enter data create(grid[0:gs6])

 #pragma acc parallel loop present(grid[0:gs6],gridf[0:gs6])
  for (int j=0;j<gs6;j++)
    grid[j] = gridf[j];

  eval_sgsd(tid,gs1,gs2,grid,val,n1,l1,m1,zeta1,Rc);

  #pragma acc exit data delete(grid[0:gs6])
  delete [] grid;

  return;
}

void eval_sgs(int tid, int gs1, int gs2, float* gridf, float* valf, int n1, int l1, int m1, double zeta1, double Rc)
{
  int gs = gs2;
  int gs6 = 6*gs;
  double* val = new double[gs];
  double* grid = new double[gs6];
  #pragma acc enter data create(val[0:gs],grid[0:gs6])

 #pragma acc parallel loop present(grid[0:gs6],gridf[0:gs6])
  for (int j=0;j<gs6;j++)
    grid[j] = gridf[j];

 #pragma acc parallel loop present(val[0:gs],valf[0:gs])
  for (int j=0;j<gs;j++)
    val[j] = valf[j];

  eval_sgsd(tid,gs1,gs2,grid,val,n1,l1,m1,zeta1,Rc);

 #pragma acc parallel loop present(valf[0:gs],val[0:gs])
  for (int j=0;j<gs;j++)
    valf[j] = val[j];

  #pragma acc exit data delete(val[0:gs2],grid[0:gs6])
  delete [] val;
  delete [] grid;

  return;
}

void eval_sgs_ked(int tid, int gs1, int gs2, double* grid, double* val, int n1, int l1, int m1, double zeta1, double Rc)
{
  //printf("  using eval_sgs_ked (gs12: %4i %4i) \n",gs1,gs2);

  int gs = gs2;
  int gs6 = 6*gs;

  //int nm1 = n1-1;
  double zeta2 = 2.*zeta1*Rc;
  //double norm2 = exp(zeta1*Rc*Rc);

  int nnm = n1*(n1-1);
  int llp = l1*(l1+1);
  double fz2 = 4.*zeta1*zeta1;
  double ta = 2.*zeta1*(1.+2.*n1);

  #pragma acc parallel loop present(grid[0:gs6],val[0:gs])
  for (int j=0;j<gs1;j++)
  {
    double r = grid[6*j+3];
    double r2 = r*r;

    val[j] *= -ta + (nnm-llp)/r2 + fz2*r2;
  }

  double f1 = nnm - llp;
  double f2 = 2.*n1*zeta2;
  double f3 = zeta2*zeta2;

  #pragma acc parallel loop present(grid[0:gs6],val[0:gs])
  for (int j=gs1;j<gs2;j++)
  {
    double r = grid[6*j+3];
    double or1 = 1./r;
    double or2 = or1*or1;

    val[j] *= f1*or2 - f2*or1 + f3;
  }

  return;
}

//eval_ss calls eval_ssd
void eval_ss(int tid, int gs, float* gridf, float* valf, int n1, int l1, int m1, double zeta1)
{
  int gs6 = 6*gs;
  double* val = new double[gs];
  double* grid = new double[gs6];
  #pragma acc enter data create(val[0:gs],grid[0:gs6])

 #pragma acc parallel loop present(grid[0:gs6],gridf[0:gs6])
  for (int j=0;j<gs6;j++)
    grid[j] = gridf[j];

 #pragma acc parallel loop present(val[0:gs],valf[0:gs])
  for (int j=0;j<gs;j++)
    val[j] = valf[j];

  eval_ssd(tid,gs,grid,val,n1,l1,m1,zeta1);

 #pragma acc parallel loop present(valf[0:gs],val[0:gs])
  for (int j=0;j<gs;j++)
    valf[j] = val[j];

  #pragma acc exit data delete(val[0:gs],grid[0:gs6])
  delete [] val;
  delete [] grid;

  return;
}

//eval_ssd (grid float) calls eval_ssd
void eval_ssd(int tid, int gs, float* gridf, double* val, int n1, int l1, int m1, double zeta1)
{
  int gs6 = 6*gs;
  double* grid = new double[gs6];
  #pragma acc enter data create(grid[0:gs6])

 #pragma acc parallel loop present(grid[0:gs6],gridf[0:gs6])
  for (int j=0;j<gs6;j++)
    grid[j] = gridf[j];

  eval_ssd(tid,gs,grid,val,n1,l1,m1,zeta1);

  #pragma acc exit data delete(grid[0:gs6])
  delete [] grid;

  return;
}

#if SS_V3
//new version
void eval_ssd(int tid, int gs, double* grid, double* val, int n1, int l1, int m1, const double zeta1)
{
  int gs6 = 6*gs;

  //n1 = l1+1;
  int nl1 = l1;

 //"1s"
  double zt1 = zeta1;
  double zt2 = zt1*zt1;
  double ztos = 0.4*zt2;
  double zt3 = zt2*zt1/15.;
  double zt4 = 0.;
  double zt5 = 0.;
  double zt6 = 0.;
  if (n1-l1>1)
  {
    ztos = 0.5*zt2;
    zt3 = zt2*zt1/6.;
    zt4 = zt2*zt2/24.;
    zt5 = zt2*zt2*zt1/120.;
    zt6 = zt2*zt2*zt2/840.;
    //ztos = 0.45*zt2;
    //zt3 = 7./60. * zt2*zt1;
    //zt4 = zt2*zt2/60.;
    //zt6 = -zt2*zt2*zt2/2100.;
  }

 #pragma acc parallel loop present(grid[0:gs6],val[0:gs])
  for (int j=0;j<gs;j++)
  {
    double r = grid[6*j+3];
    double ezr = exp(-zeta1*r);
    double r2 = r*r;
    double r4 = r2*r2;
    double rnl = pow(r,nl1); //to compensate r^-l in eval_sh_3rd

    double v1 = (1. + zt1*r + ztos*r2 + zt3*r*r2 + zt4*r4 + zt5*r4*r + zt6*r4*r2);

    val[j] *= v1*rnl*ezr;
  }

  eval_sh_3rd(tid,gs,grid,val,n1,l1,m1);

  return;
}
#endif

#if SS_V1
void eval_ssd(int tid, int gs, double* grid, double* val, int n1, int l1, int m1, const double zeta1)
{
  int gs6 = 6*gs;

  //if (n1-l1>1) { printf(" ERROR: cannot use radial nodes in eval_ssd \n"); exit(-1); }

  int nl1 = n1-1;

  double zt1 = zeta1;
  double zt2 = zeta1*zeta1/2.;
  //double zt3 = 0.;
  //if (n1-l1-1==1) zt3 = zt2*zeta1/3.;
  //if (n1-l1-1==1) zt2 = 0.;
  //if (n1-l1-1==2) zt1 = zt2 = 0.;

  if (zt2==0.) printf(" eval_ss. nl: %i %i  zt12: %8.5f %8.5f  zeta: %8.5f \n",n1,l1,zt1,zt2,zeta1);

  //n1 = l1+1;
  nl1 = l1; //r power

  if (n1-l1-1==0)
 #pragma acc parallel loop present(grid[0:gs6],val[0:gs])
  for (int j=0;j<gs;j++)
  {
    double r = grid[6*j+3];

    double rp1 = r; double rp2 = r*r; //double rp3 = rp2*r;

    double ezr = exp(-zeta1*r);
    double rnl = pow(r,nl1); //to compensate r^-l in eval_sh_3rd
    double v1 = 1. + zt1*rp1 + zt2*rp2;

    val[j] *= v1*rnl*ezr;
  }

  int nrp = n1-l1+2;
  if (n1-l1-1>0)
 #pragma acc parallel loop present(grid[0:gs6],val[0:gs])
  for (int j=0;j<gs;j++)
  {
    double r = grid[6*j+3];

    double ezr = exp(-zeta1*r);
    double rnl = pow(r,nl1); //to compensate r^-l in eval_sh_3rd

    double vt = 0.;
    #pragma acc loop reduction(+:vt)
    for (int k=0;k<nrp;k++)
    {
      double dn = 1./facto(k);
      double rzp = pow(r*zeta1,k);

      vt += rzp*dn;
    }

    val[j] *= vt*rnl*ezr;
  }

  eval_sh_3rd(tid,gs,grid,val,n1,l1,m1);

  return;
}
#endif

void eval_ss_ked(int tid, int gs, double* grid, double* val, int n1, int l1, int m1, const double zeta1)
{
  int gs6 = 6*gs;

  bool first_loop = (n1-l1==1);
 #if SS_V3
  first_loop = 1;
 #endif
  //printf("  eval_ss_ked  first_loop: %i \n",(int)first_loop);

  double zt1 = zeta1;
  double zt2 = zeta1*zeta1;
  double ztos = zt2/2.;

 #if SS_V3
  ztos = 0.4*zt2;
  double zt3 = zt2*zt1/15.;
  double zt4 = 0.;
  double zt5 = 0.;
  double zt6 = 0.;
  if (n1-l1>1)
  {
    ztos = 0.5*zt2;
    zt3 = zt2*zt1/6.;
    zt4 = zt2*zt2/24.;
    zt5 = zt2*zt2*zt1/120.;
    zt6 = zt2*zt2*zt2/840.;
    //ztos = 0.45*zt2;
    //zt3 = 7./60. * zt2*zt1;
    //zt4 = zt2*zt2/60.;
    //zt6 = -zt2*zt2*zt2/2100.;
  }
 #endif

  int nl1 = l1; //r power

      n1 = l1+1;
  int n2 = n1+1;
  int n3 = n1+2;
  int n4 = n1+3;
  int n5 = n1+4;
  int n6 = n1+5;
  int n7 = n1+6;
  int llp = l1*(l1+1);

  double f11 = n1*(n1-1) - llp;
  double f12 = n2*(n2-1) - llp;
  double f13 = n3*(n3-1) - llp;
  double f21 = 2.*n1*zeta1;
  double f22 = 2.*n2*zeta1;
  double f23 = 2.*n3*zeta1;
  double f3 = zeta1*zeta1;

 #if SS_V3
  double f14 = n4*(n4-1) - llp;
  double f15 = n5*(n5-1) - llp;
  double f16 = n6*(n6-1) - llp;
  double f17 = n7*(n7-1) - llp;
  double f24 = 2.*n4*zeta1;
  double f25 = 2.*n5*zeta1;
  double f26 = 2.*n6*zeta1;
  double f27 = 2.*n7*zeta1;
 #endif

  if (first_loop)
 #pragma acc parallel loop present(grid[0:gs6],val[0:gs])
  for (int j=0;j<gs;j++)
  {
    double r = grid[6*j+3];
    double ezr = exp(-zeta1*r);
    double rnl = pow(r,nl1); //to compensate r^-l in eval_sh_3rd
    double or1 = 1./r;
    double or2 = or1*or1;

    double rp1 = r; //pow(r,nr2);
    double rp2 = r*r; //pow(r,nr3);

    double v1 = 1. *      (f11*or2 - f21*or1 + f3);
    double v2 = zt1 * (f12*or1 - f22 + f3*rp1);
    double v3 = ztos* (f13 - f23*rp1 + f3*rp2);

   #if SS_V3

    double rp3 = rp2*rp1;
    double rp4 = rp2*rp2;
    double rp5 = rp3*rp2;
    double rp6 = rp3*rp3;
    double v4 = zt3 * (f14*rp1 - f24*rp2 + f3*rp3);
    double v5 = zt4 * (f15*rp2 - f25*rp3 + f3*rp4);
    double v6 = zt5 * (f16*rp3 - f26*rp4 + f3*rp5);
    double v7 = zt6 * (f17*rp4 - f27*rp5 + f3*rp6);
    val[j] *= (v1+v2+v3+v4+v5+v6+v7)*rnl*ezr;

   #else
    val[j] *= (v1+v2+v3)*rnl*ezr;
   #endif
  }

 #if SS_V1
  int nrp = n1-l1+2;
  //printf(" n1/l1: %i %i  nrp: %i \n",n1,l1,nrp);
  if (!first_loop)
 #pragma acc parallel loop present(grid[0:gs6],val[0:gs])
  for (int j=0;j<gs;j++)
  {
    double r = grid[6*j+3];
    double ezr = exp(-zeta1*r);
    double rnl = pow(r,nl1); //to compensate r^-l in eval_sh_3rd
    double or1 = 1./r;
    double or2 = or1*or1;

    double vt = 0.;
    #pragma acc loop reduction(+:vt)
    for (int k=0;k<nrp;k++)
    {
      double dn = 1./facto(k);
      double rzp = pow(r*zeta1,k);
      int nn = 1+l1+k;
      double f1 = nn*(nn-1) - llp;
      double f2 = 2.*nn*zeta1;

      vt += rzp*dn*(f1*or2 - f2*or1 + f3);
    }

    val[j] *= vt*rnl*ezr;
  }
 #endif

  eval_sh_3rd(tid,gs,grid,val,-1,l1,m1);

  return;
}

void eval_pd(int gs, double* grid1, double* val, int n1, int l1, int m1, double zeta1);

void eval_ss_pd(int gs, double* grid, double* val, double* tmp, int n1, int l1, int m1, const double zeta1)
{
  int gs3 = 3*gs;

  double zt1 = zeta1;
  double zt2 = zt1*zt1;
  double zt3 = zt2*zt1;
  double zt4 = zt2*zt2;
  double zt6 = zt4*zt2;

  int ns = 3;
  if (n1-l1>1) ns = 5;
 #if SS_V3
  ns = 4;
  if (n1-l1>1) ns = 7;
 #endif

  #pragma acc parallel loop present(val[0:gs3])
  for (int j=0;j<gs3;j++)
    val[j] = 0.;

  //printf("  nl: %i %i  nml1: %i ns:  %i \n",n1,l1,nml1,ns);

  for (int n=0;n<ns;n++)
  {
    int n2 = l1+1+n;
   #if SS_V1
    int nml = n;
    double ztn = pow(zeta1,nml);
    double f1 = 1.; for (int i=1;i<=n;i++) f1 *= i;
    ztn /= f1;
   #endif

   #if SS_V3
    double ztn = 1.;
    if (n==1) ztn = zt1;
    if (n==2) ztn = 0.4*zt2;
    if (n==3) ztn = zt3/15.;
    if (ns==7)
    {
      if (n==2) ztn = 0.5*zt2;
      if (n==3) ztn = zt2*zt1/6.;
      if (n==4) ztn = zt2*zt2/24.;
      if (n==5) ztn = zt2*zt2*zt1/120.;
      if (n==6) ztn = zt2*zt2*zt2/840.;
      //if (n==2) ztn = 0.45*zt2;
      //if (n==3) ztn = 7./60.*zt3;
      //if (n==4) ztn = zt4/60.;
      //if (n==5) { ztn = -zt6/2100.; n2++; }
    }
   #endif

    if (ztn!=0.)
    {
      #pragma acc parallel loop present(tmp[0:gs3])
      for (int j=0;j<gs3;j++)
        tmp[j] = ztn;

      eval_pd(gs,grid,tmp,n2,l1,m1,zeta1);

      #pragma acc parallel loop present(val[0:gs3],tmp[0:gs3])
      for (int j=0;j<gs;j++)
      {
        val[3*j+0] += tmp[3*j+0];
        val[3*j+1] += tmp[3*j+1];
        val[3*j+2] += tmp[3*j+2];
      }
    } //if ztn

  }

  return;
}
