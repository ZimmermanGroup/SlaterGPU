#include "spherical.h"


// start spherical * exp //
//  double precision     //


void get_1s_expd(int tid, int gs, double* grid, double* val, float zeta1)
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

void get_2s_expd(int tid, int gs, double* grid, double* val, float zeta1)
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

void get_3s_expd(int tid, int gs, double* grid, double* val, float zeta1)
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

void get_4s_expd(int tid, int gs, double* grid, double* val, float zeta1)
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

void get_5s_expd(int tid, int gs, double* grid, double* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double r4 = pow(r,4.f);
    double r4ezr = r4*exp(-zeta1*r);
    val[i] *= r4ezr;
  }
  return;
}

void get_6s_expd(int tid, int gs, double* grid, double* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double r5 = pow(r,5.f);
    double r5ezr = r5*exp(-zeta1*r);
    val[i] *= r5ezr;
  }
  return;
}

void get_7s_expd(int tid, int gs, double* grid, double* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double r6 = pow(r,6.f);
    double r6ezr = r6*exp(-zeta1*r);
    val[i] *= r6ezr;
  }
  return;
}

void get_8s_expd(int tid, int gs, double* grid, double* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double r7 = pow(r,7.f);
    double r7ezr = r7*exp(-zeta1*r);
    val[i] *= r7ezr;
  }
  return;
}

void get_9s_expd(int tid, int gs, double* grid, double* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double r8 = pow(r,8.f);
    double r8ezr = r8*exp(-zeta1*r);
    val[i] *= r8ezr;
  }
  return;
}

void get_10s_expd(int tid, int gs, double* grid, double* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double r9 = pow(r,9.f);
    double r9ezr = r9*exp(-zeta1*r);
    val[i] *= r9ezr;
  }
  return;
}

void get_11s_expd(int tid, int gs, double* grid, double* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double r10 = pow(r,10.f);
    double r10ezr = r10*exp(-zeta1*r);
    val[i] *= r10ezr;
  }
  return;
}

void get_12s_expd(int tid, int gs, double* grid, double* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double r11 = pow(r,11.f);
    double r11ezr = r11*exp(-zeta1*r);
    val[i] *= r11ezr;
  }
  return;
}

void get_2pxd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_2pyd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_2pzd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_3pxd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_3pyd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_3pzd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_4pxd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_4pyd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_4pzd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_5pxd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_5pyd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_5pzd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_6pxd(int tid, int gs, double* grid, double* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double r = grid[6*i+3];
    double r4 = pow(r,4.f);
    double rezr = r4*exp(-zeta*r);
    val[i] *= x*rezr;
  }
}

void get_6pyd(int tid, int gs, double* grid, double* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double r4 = pow(r,4.f);
    double rezr = r4*exp(-zeta*r);
    val[i] *= y*rezr;
  }
}

void get_6pzd(int tid, int gs, double* grid, double* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r4 = pow(r,4.f);
    double rezr = r4*exp(-zeta*r);
    val[i] *= z*rezr;
  }
}

void get_7pxd(int tid, int gs, double* grid, double* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double r = grid[6*i+3];
    double r5 = pow(r,5.f);
    double rezr = r5*exp(-zeta*r);
    val[i] *= x*rezr;
  }
}

void get_7pyd(int tid, int gs, double* grid, double* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double r5 = pow(r,5.f);
    double rezr = r5*exp(-zeta*r);
    val[i] *= y*rezr;
  }
}

void get_7pzd(int tid, int gs, double* grid, double* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r5 = pow(r,5.f);
    double rezr = r5*exp(-zeta*r);
    val[i] *= z*rezr;
  }
}

void get_8pxd(int tid, int gs, double* grid, double* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double r = grid[6*i+3];
    double r6 = pow(r,6.f);
    double rezr = r6*exp(-zeta*r);
    val[i] *= x*rezr;
  }
}

void get_8pyd(int tid, int gs, double* grid, double* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double r6 = pow(r,6.f);
    double rezr = r6*exp(-zeta*r);
    val[i] *= y*rezr;
  }
}

void get_8pzd(int tid, int gs, double* grid, double* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r6 = pow(r,6.f);
    double rezr = r6*exp(-zeta*r);
    val[i] *= z*rezr;
  }
}

void get_9pxd(int tid, int gs, double* grid, double* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double r = grid[6*i+3];
    double r7 = pow(r,7.f);
    double rezr = r7*exp(-zeta*r);
    val[i] *= x*rezr;
  }
}

void get_9pyd(int tid, int gs, double* grid, double* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double y = grid[6*i+1];
    double r = grid[6*i+3];
    double r7 = pow(r,7.f);
    double rezr = r7*exp(-zeta*r);
    val[i] *= y*rezr;
  }
}

void get_9pzd(int tid, int gs, double* grid, double* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r7 = pow(r,7.f);
    double rezr = r7*exp(-zeta*r);
    val[i] *= z*rezr;
  }
}

void get_3dxyd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_3dyzd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_3dz2d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (2.f*z*z-x*x-y*y)*ezr;

   //ambiguous for which r to use
    //val[i] *= 3.f*z*z-r*r; 
  }
}

void get_3dxzd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_3dx2y2d(int tid, int gs, double* grid, double* val, float zeta)
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
void get_3dxxd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_3dyyd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_3dzzd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_4dxyd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_4dyzd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_4dz2d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (2.f*z*z-x*x-y*y)*ezr;

   //ambiguous for which r to use
    //val[i] *= 3.f*z*z-r*r; 
  }
}

void get_4dxzd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_4dx2y2d(int tid, int gs, double* grid, double* val, float zeta)
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

void get_5dxyd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_5dyzd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_5dz2d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (2.f*z*z-x*x-y*y)*ezr;

   //ambiguous for which r to use
    //val[i] *= 3.f*z*z-r*r; 
  }
}

void get_5dxzd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_5dx2y2d(int tid, int gs, double* grid, double* val, float zeta)
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

void get_6dxyd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_6dyzd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_6dz2d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (2.f*z*z-x*x-y*y)*ezr;

   //ambiguous for which r to use
    //val[i] *= 3.f*z*z-r*r; 
  }
}

void get_6dxzd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_6dx2y2d(int tid, int gs, double* grid, double* val, float zeta)
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

void get_7dxyd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_7dyzd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_7dz2d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (2.f*z*z-x*x-y*y)*ezr;

   //ambiguous for which r to use
    //val[i] *= 3.f*z*z-r*r; 
  }
}

void get_7dxzd(int tid, int gs, double* grid, double* val, float zeta)
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

void get_7dx2y2d(int tid, int gs, double* grid, double* val, float zeta)
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

void get_4fm3d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (3.f*x*x-y*y)*y*ezr;
  }
}

void get_4fm2d(int tid, int gs, double* grid, double* val, float zeta)
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

void get_4fm1d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (4.f*z*z-x*x-y*y)*y*ezr;
  }
}

void get_4f0d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (2.f*z*z-3.f*x*x-3.f*y*y)*z*ezr;
  }
}

void get_4fp1d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (4.f*z*z-x*x-y*y)*x*ezr;
  }
}

void get_4fp2d(int tid, int gs, double* grid, double* val, float zeta)
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

void get_4fp3d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (x*x-3.f*y*y)*x*ezr;
  }
}

void get_5fm3d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (3.f*x*x-y*y)*y*ezr;
  }
}

void get_5fm2d(int tid, int gs, double* grid, double* val, float zeta)
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

void get_5fm1d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (4.f*z*z-x*x-y*y)*y*ezr;
  }
}

void get_5f0d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (2.f*z*z-3.f*x*x-3.f*y*y)*z*ezr;
  }
}

void get_5fp1d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (4.f*z*z-x*x-y*y)*x*ezr;
  }
}

void get_5fp2d(int tid, int gs, double* grid, double* val, float zeta)
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

void get_5fp3d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (x*x-3.f*y*y)*x*ezr;
  }
}

void get_6fm3d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (3.f*x*x-y*y)*y*ezr;
  }
}

void get_6fm2d(int tid, int gs, double* grid, double* val, float zeta)
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

void get_6fm1d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (4.f*z*z-x*x-y*y)*y*ezr;
  }
}

void get_6f0d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (2.f*z*z-3.f*x*x-3.f*y*y)*z*ezr;
  }
}

void get_6fp1d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (4.f*z*z-x*x-y*y)*x*ezr;
  }
}

void get_6fp2d(int tid, int gs, double* grid, double* val, float zeta)
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

void get_6fp3d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (x*x-3.f*y*y)*x*ezr;
  }
}

void get_5gm4d(int tid, int gs, double* grid, double* val, float zeta)
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

void get_5gm3d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= y*z * (3.f*x*x - y*y) * ezr;
  }
}

void get_5gm2d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= x*y * (6.f*z*z - x*x - y*y) * ezr;
  }
}

void get_5gm1d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= y*z * (4.f*z*z - 3.f*x*x - 3.f*y*y) * ezr;
  }
}

void get_5g0d(int tid, int gs, double* grid, double* val, float zeta)
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
    //val[i] *= (35.f*z2*z2 - 30.f*z2*r2 + 3.f*r2*r2) * ezr; //equal via algebra
    val[i] *= (35.*z2*z2 - 30.*z2*r2 + 3.*r2*r2) * ezr;
    //val[i] *= (3.f*x2py2*x2py2 - 24.f*x2py2*z2 + 8.f*z2*z2) * ezr;
  }
}

void get_5gp1d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= x*z * (4.*z*z - 3.f*x*x - 3.f*y*y) * ezr;
  }
}

void get_5gp2d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (x2 - y2) * (6.f*z*z - x2 - y2) * ezr;
  }
}

void get_5gp3d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= x*z * (x*x - 3.f*y*y) * ezr;
  }
}

void get_5gp4d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (x2 * (x2 - 3.f*y2) - y2 * (3.f*x2 - y2)) * ezr;
  }
}

void get_6hm5d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (5.f*x2*x2 - 10.f*x2*y2 + y2*y2)*y * ezr;
  }
}

void get_6hm4d(int tid, int gs, double* grid, double* val, float zeta)
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
  }
}

void get_6hm3d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (x2*x2 - 10.f*x2*y2 + 5.f*y2*y2)*x * ezr;
  }
}

void get_6hm2d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= y * (y2-3.f*x2)*(x2 + y2 - 8.f*z2) * ezr;
  }
}

void get_6hm1d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= z * (x2*x2 - 6.f*x2*y2 + y2*y2) * ezr;
  }
}

void get_6h0d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (x2 + y2 - 2.f*z2)*x*y*z * ezr;
  }
}

void get_6hp1d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= x * (x2-3.f*y2)*(x2 + y2 - 8.f*z2) * ezr;
  }
}

void get_6hp2d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= y * (x2*x2 + y2*y2 - 12.f*y2*z2 + 8.f*z2*z2 + 2.f*x2 * (y2-6.f*z2)) * ezr;
  }
}

void get_6hp3d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= (x2 - y2) * (x2 + y2 - 2.f*z2) * z * ezr;
  }
}

void get_6hp4d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= z * (15.f*x2*x2 + 15.f*y2*y2 - 40.f*y2*z2 + 8.f*z2*z2 + 10.f*x2 * (3.f*y2 - 4.f*z2)) * ezr;
  }
}

void get_6hp5d(int tid, int gs, double* grid, double* val, float zeta)
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
    val[i] *= x * (x2*x2 + y2*y2 - 12.f*y2*z2 + 8.f*z2*z2 + 2.f*x2* (y2 - 6.f*z2)) * ezr;
  }
}

// end spherical * exp //

void eval_shd(int tid, int gs, double* grid, double* val, int n1, int l1, int m1, float zeta1)
{
 //1s - 12s
 //2p - 9p
 //3d,4f,5g,6h
 
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
        if (m1==0) //2p
          get_2pxd(tid,gs,grid,val,zeta1);
        else if (m1==1)
          get_2pyd(tid,gs,grid,val,zeta1);
        else
          get_2pzd(tid,gs,grid,val,zeta1);
      }
      else if (n1==3)
      {
        if (m1==0) //3p
          get_3pxd(tid,gs,grid,val,zeta1);
        else if (m1==1)
          get_3pyd(tid,gs,grid,val,zeta1);
        else
          get_3pzd(tid,gs,grid,val,zeta1);
      }
      else if (n1==4)
      {
        if (m1==0) //4p
          get_4pxd(tid,gs,grid,val,zeta1);
        else if (m1==1)
          get_4pyd(tid,gs,grid,val,zeta1);
        else
          get_4pzd(tid,gs,grid,val,zeta1);
      }
      else if (n1==5)
      {
        if (m1==0) //5p
          get_5pxd(tid,gs,grid,val,zeta1);
        else if (m1==1)
          get_5pyd(tid,gs,grid,val,zeta1);
        else
          get_5pzd(tid,gs,grid,val,zeta1);
      }
      else if (n1==6)
      {
        if (m1==0) //6p
          get_6pxd(tid,gs,grid,val,zeta1);
        else if (m1==1)
          get_6pyd(tid,gs,grid,val,zeta1);
        else
          get_6pzd(tid,gs,grid,val,zeta1);
      }
      else if (n1==7)
      {
        if (m1==0) //7p
          get_7pxd(tid,gs,grid,val,zeta1);
        else if (m1==1)
          get_7pyd(tid,gs,grid,val,zeta1);
        else
          get_7pzd(tid,gs,grid,val,zeta1);
      }
      else if (n1==8)
      {
        if (m1==0) //8p
          get_8pxd(tid,gs,grid,val,zeta1);
        else if (m1==1)
          get_8pyd(tid,gs,grid,val,zeta1);
        else
          get_8pzd(tid,gs,grid,val,zeta1);
      }
      else if (n1==9)
      {
        if (m1==0) //9p
          get_9pxd(tid,gs,grid,val,zeta1);
        else if (m1==1)
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
  } //if l1>0

  return;
}

