#include "spherical.h"


// if integrating in PS coordinates, make sure nphi > 2l
 // this is based on tests on dihydrogen with higher angular momentum


// start spherical * exp //
//  double precision     //


void get_1s_expd(int tid, int gs, double* grid, double* val, double zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double ezr = r*r*exp(-zeta*r);
    val[i] *= (x*x-3.*y*y)*x*ezr;
  }
}

void get_5gm4d(int tid, int gs, double* grid, double* val, double zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double r = grid[6*i+3];
    val[i] *= x/r;
  }
}

void get_py_3rd(int tid, int gs, double* grid, double* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    val[i] *= y/r;
  }
}

void get_pz_3rd(int tid, int gs, double* grid, double* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    val[i] *= z/r;
  }
}

void get_dxy_3rd(int tid, int gs, double* grid, double* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    val[i] *= 3.*z*z/r/r-1.;
  }
}

void get_dxz_3rd(int tid, int gs, double* grid, double* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double r = grid[6*i+3];
    val[i] *= x*x/r/r;
  }
}

void get_dyy_3rd(int tid, int gs, double* grid, double* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    val[i] *= y*y/r/r;
  }
}

void get_dzz_3rd(int tid, int gs, double* grid, double* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    val[i] *= z*z/r/r;
  }
}

void get_fm3_3rd(int tid, int gs, double* grid, double* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
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
 //3-7d,4-6f,5g,6h
 
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
      else if (n1>6)
      {
        printf(" ERROR: 7f not available \n");
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

