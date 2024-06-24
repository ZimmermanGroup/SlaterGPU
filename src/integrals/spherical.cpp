#include "spherical.h"

// start spherical * exp //

/* 
1s-7s
2p-7p
3d-7d
4f-6f
5g
6h
*/

void get_1s_exp(int tid, int gs, float* grid, float* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float ezr = expf(-zeta1*r);
    val[i] *= ezr;
  }
  return;
}

void get_2s_exp(int tid, int gs, float* grid, float* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float rezr = r*expf(-zeta1*r);
    val[i] *= rezr;
  }
  return;
}

void get_3s_exp(int tid, int gs, float* grid, float* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float r2ezr = r*r*expf(-zeta1*r);
    val[i] *= r2ezr;
  }
  return;
}

void get_4s_exp(int tid, int gs, float* grid, float* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float r2 = r*r;
    float r3ezr = r2*r*expf(-zeta1*r);
    val[i] *= r3ezr;
  }
  return;
}

void get_5s_exp(int tid, int gs, float* grid, float* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float r4 = powf(r,4.f);
    float r4ezr = r4*expf(-zeta1*r);
    val[i] *= r4ezr;
  }
  return;
}

void get_6s_exp(int tid, int gs, float* grid, float* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float r5 = powf(r,5.f);
    float r5ezr = r5*expf(-zeta1*r);
    val[i] *= r5ezr;
  }
  return;
}

void get_7s_exp(int tid, int gs, float* grid, float* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float r6 = powf(r,6.f);
    float r6ezr = r6*expf(-zeta1*r);
    val[i] *= r6ezr;
  }
  return;
}

void get_8s_exp(int tid, int gs, float* grid, float* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float r7 = powf(r,7.f);
    float r7ezr = r7*expf(-zeta1*r);
    val[i] *= r7ezr;
  }
  return;
}

void get_9s_exp(int tid, int gs, float* grid, float* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float r8 = powf(r,8.f);
    float r8ezr = r8*expf(-zeta1*r);
    val[i] *= r8ezr;
  }
  return;
}

void get_10s_exp(int tid, int gs, float* grid, float* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float r9 = powf(r,9.f);
    float r9ezr = r9*expf(-zeta1*r);
    val[i] *= r9ezr;
  }
  return;
}

void get_11s_exp(int tid, int gs, float* grid, float* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float r10 = powf(r,10.f);
    float r10ezr = r10*expf(-zeta1*r);
    val[i] *= r10ezr;
  }
  return;
}

void get_12s_exp(int tid, int gs, float* grid, float* val, float zeta1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float r11 = powf(r,11.f);
    float r11ezr = r11*expf(-zeta1*r);
    val[i] *= r11ezr;
  }
  return;
}

void get_2px(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= x*ezr;
  }
}

void get_2py(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= y*ezr;
  }
}

void get_2pz(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= z*ezr;
  }
}

void get_3px(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float r = grid[6*i+3];
    float rezr = r*expf(-zeta*r);
    val[i] *= x*rezr;
  }
}

void get_3py(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float rezr = r*expf(-zeta*r);
    val[i] *= y*rezr;
  }
}

void get_3pz(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float rezr = r*expf(-zeta*r);
    val[i] *= z*rezr;
  }
}

void get_4px(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float r = grid[6*i+3];
    float rezr = r*r*expf(-zeta*r);
    val[i] *= x*rezr;
  }
}

void get_4py(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float rezr = r*r*expf(-zeta*r);
    val[i] *= y*rezr;
  }
}

void get_4pz(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float rezr = r*r*expf(-zeta*r);
    val[i] *= z*rezr;
  }
}

void get_5px(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float r = grid[6*i+3];
    float rezr = r*r*r*expf(-zeta*r);
    val[i] *= x*rezr;
  }
}

void get_5py(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float rezr = r*r*r*expf(-zeta*r);
    val[i] *= y*rezr;
  }
}

void get_5pz(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float rezr = r*r*r*expf(-zeta*r);
    val[i] *= z*rezr;
  }
}

void get_6px(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float r = grid[6*i+3];
    float r4 = powf(r,4.f);
    float rezr = r4*expf(-zeta*r);
    val[i] *= x*rezr;
  }
}

void get_6py(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float r4 = powf(r,4.f);
    float rezr = r4*expf(-zeta*r);
    val[i] *= y*rezr;
  }
}

void get_6pz(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float r4 = powf(r,4.f);
    float rezr = r4*expf(-zeta*r);
    val[i] *= z*rezr;
  }
}

void get_7px(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float r = grid[6*i+3];
    float r5 = powf(r,5.f);
    float rezr = r5*expf(-zeta*r);
    val[i] *= x*rezr;
  }
}

void get_7py(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float r5 = powf(r,5.f);
    float rezr = r5*expf(-zeta*r);
    val[i] *= y*rezr;
  }
}

void get_7pz(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float r5 = powf(r,5.f);
    float rezr = r5*expf(-zeta*r);
    val[i] *= z*rezr;
  }
}

void get_8px(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float r = grid[6*i+3];
    float r6 = powf(r,6.f);
    float rezr = r6*expf(-zeta*r);
    val[i] *= x*rezr;
  }
}

void get_8py(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float r6 = powf(r,6.f);
    float rezr = r6*expf(-zeta*r);
    val[i] *= y*rezr;
  }
}

void get_8pz(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float r6 = powf(r,6.f);
    float rezr = r6*expf(-zeta*r);
    val[i] *= z*rezr;
  }
}

void get_9px(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float r = grid[6*i+3];
    float r7 = powf(r,7.f);
    float rezr = r7*expf(-zeta*r);
    val[i] *= x*rezr;
  }
}

void get_9py(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float r7 = powf(r,7.f);
    float rezr = r7*expf(-zeta*r);
    val[i] *= y*rezr;
  }
}

void get_9pz(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float r7 = powf(r,7.f);
    float rezr = r7*expf(-zeta*r);
    val[i] *= z*rezr;
  }
}

void get_3dxy(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= x*y*ezr;
  }
}

void get_3dyz(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= y*z*ezr;
  }
}

void get_3dz2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= (2.f*z*z-x*x-y*y)*ezr;

   //ambiguous for which r to use
    //val[i] *= 3.f*z*z-r*r; 
  }
}

void get_3dxz(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= x*z*ezr;
  }
}

void get_3dx2y2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= (x*x-y*y)*ezr;
  }
}

//for Cartesian expansion of d functions
void get_3dxx(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= x*x*ezr;
  }
}

void get_3dyy(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= y*y*ezr;
  }
}

void get_3dzz(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= z*z*ezr;
  }
}

void get_4dxy(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float ezr = r*expf(-zeta*r);
    val[i] *= x*y*ezr;
  }
}

void get_4dyz(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = r*expf(-zeta*r);
    val[i] *= y*z*ezr;
  }
}

void get_4dz2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = r*expf(-zeta*r);
    val[i] *= (2.f*z*z-x*x-y*y)*ezr;

   //ambiguous for which r to use
    //val[i] *= 3.f*z*z-r*r; 
  }
}

void get_4dxz(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = r*expf(-zeta*r);
    val[i] *= x*z*ezr;
  }
}

void get_4dx2y2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float ezr = r*expf(-zeta*r);
    val[i] *= (x*x-y*y)*ezr;
  }
}

void get_5dxy(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float ezr = r*r*expf(-zeta*r);
    val[i] *= x*y*ezr;
  }
}

void get_5dyz(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = r*r*expf(-zeta*r);
    val[i] *= y*z*ezr;
  }
}

void get_5dz2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = r*r*expf(-zeta*r);
    val[i] *= (2.f*z*z-x*x-y*y)*ezr;

   //ambiguous for which r to use
    //val[i] *= 3.f*z*z-r*r; 
  }
}

void get_5dxz(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = r*r*expf(-zeta*r);
    val[i] *= x*z*ezr;
  }
}

void get_5dx2y2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float ezr = r*r*expf(-zeta*r);
    val[i] *= (x*x-y*y)*ezr;
  }
}

void get_6dxy(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float r3 = powf(r,3.);
    float ezr = r3*expf(-zeta*r);
    val[i] *= x*y*ezr;
  }
}

void get_6dyz(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float r3 = powf(r,3.);
    float ezr = r3*expf(-zeta*r);
    val[i] *= y*z*ezr;
  }
}

void get_6dz2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float r3 = powf(r,3.);
    float ezr = r3*expf(-zeta*r);
    val[i] *= (2.f*z*z-x*x-y*y)*ezr;

   //ambiguous for which r to use
    //val[i] *= 3.f*z*z-r*r; 
  }
}

void get_6dxz(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float r3 = powf(r,3.);
    float ezr = r3*expf(-zeta*r);
    val[i] *= x*z*ezr;
  }
}

void get_6dx2y2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float r3 = powf(r,3.);
    float ezr = r3*expf(-zeta*r);
    val[i] *= (x*x-y*y)*ezr;
  }
}

void get_7dxy(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float r4 = powf(r,4.);
    float ezr = r4*expf(-zeta*r);
    val[i] *= x*y*ezr;
  }
}

void get_7dyz(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float r4 = powf(r,4.);
    float ezr = r4*expf(-zeta*r);
    val[i] *= y*z*ezr;
  }
}

void get_7dz2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float r4 = powf(r,4.);
    float ezr = r4*expf(-zeta*r);
    val[i] *= (2.f*z*z-x*x-y*y)*ezr;

   //ambiguous for which r to use
    //val[i] *= 3.f*z*z-r*r; 
  }
}

void get_7dxz(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float r4 = powf(r,4.);
    float ezr = r4*expf(-zeta*r);
    val[i] *= x*z*ezr;
  }
}

void get_7dx2y2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float r4 = powf(r,4.);
    float ezr = r4*expf(-zeta*r);
    val[i] *= (x*x-y*y)*ezr;
  }
}

void get_4fm3(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= (3.f*x*x-y*y)*y*ezr;
  }
}

void get_4fm2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= x*y*z*ezr;
  }
}

void get_4fm1(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= (4.f*z*z-x*x-y*y)*y*ezr;
  }
}

void get_4f0(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= (2.f*z*z-3.f*x*x-3.f*y*y)*z*ezr;
  }
}

void get_4fp1(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= (4.f*z*z-x*x-y*y)*x*ezr;
  }
}

void get_4fp2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= (x*x-y*y)*z*ezr;
  }
}

void get_4fp3(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= (x*x-3.f*y*y)*x*ezr;
  }
}

void get_5fm3(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float ezr = r*expf(-zeta*r);
    val[i] *= (3.f*x*x-y*y)*y*ezr;
  }
}

void get_5fm2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = r*expf(-zeta*r);
    val[i] *= x*y*z*ezr;
  }
}

void get_5fm1(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = r*expf(-zeta*r);
    val[i] *= (4.f*z*z-x*x-y*y)*y*ezr;
  }
}

void get_5f0(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = r*expf(-zeta*r);
    val[i] *= (2.f*z*z-3.f*x*x-3.f*y*y)*z*ezr;
  }
}

void get_5fp1(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = r*expf(-zeta*r);
    val[i] *= (4.f*z*z-x*x-y*y)*x*ezr;
  }
}

void get_5fp2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = r*expf(-zeta*r);
    val[i] *= (x*x-y*y)*z*ezr;
  }
}

void get_5fp3(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float ezr = r*expf(-zeta*r);
    val[i] *= (x*x-3.f*y*y)*x*ezr;
  }
}

void get_6fm3(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float ezr = r*r*expf(-zeta*r);
    val[i] *= (3.f*x*x-y*y)*y*ezr;
  }
}

void get_6fm2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = r*r*expf(-zeta*r);
    val[i] *= x*y*z*ezr;
  }
}

void get_6fm1(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = r*r*expf(-zeta*r);
    val[i] *= (4.f*z*z-x*x-y*y)*y*ezr;
  }
}

void get_6f0(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = r*r*expf(-zeta*r);
    val[i] *= (2.f*z*z-3.f*x*x-3.f*y*y)*z*ezr;
  }
}

void get_6fp1(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = r*r*expf(-zeta*r);
    val[i] *= (4.f*z*z-x*x-y*y)*x*ezr;
  }
}

void get_6fp2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = r*r*expf(-zeta*r);
    val[i] *= (x*x-y*y)*z*ezr;
  }
}

void get_6fp3(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float ezr = r*r*expf(-zeta*r);
    val[i] *= (x*x-3.f*y*y)*x*ezr;
  }
}

void get_5gm4(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    //float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= x*y * (x*x - y*y) * ezr;
  }
}

void get_5gm3(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= y*z * (3.f*x*x - y*y) * ezr;
  }
}

void get_5gm2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= x*y * (6.f*z*z - x*x - y*y) * ezr;
  }
}

void get_5gm1(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= y*z * (4.f*z*z - 3.f*x*x - 3.f*y*y) * ezr;
  }
}

void get_5g0(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    float z2 = z*z;
    float x2py2 = x*x+y*y;
    float r2 = r*r;
    //val[i] *= (35.f*z2*z2 - 30.f*z2*r2 + 3.f*r2*r2) * ezr; //equal via algebra
    val[i] *= (35.*z2*z2 - 30.*z2*r2 + 3.*r2*r2) * ezr;
    //val[i] *= (3.f*x2py2*x2py2 - 24.f*x2py2*z2 + 8.f*z2*z2) * ezr;
  }
}

void get_5gp1(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= x*z * (4.*z*z - 3.f*x*x - 3.f*y*y) * ezr;
  }
}

void get_5gp2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    float x2 = x*x;
    float y2 = y*y;
    val[i] *= (x2 - y2) * (6.f*z*z - x2 - y2) * ezr;
  }
}

void get_5gp3(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    val[i] *= x*z * (x*x - 3.f*y*y) * ezr;
  }
}

void get_5gp4(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    //float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    float x2 = x*x;
    float y2 = y*y;
    val[i] *= (x2 * (x2 - 3.f*y2) - y2 * (3.f*x2 - y2)) * ezr;
  }
}

void get_6hm5(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    //float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    float x2 = x*x;
    float y2 = y*y;

    val[i] *= -(5.f*x2*x2 - 10.f*x2*y2 + y2*y2)*y * ezr;
    //val[i] *= (5.f*x2*x2 - 10.f*x2*y2 + y2*y2)*y * ezr;
  }
}

void get_6hm4(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    float x2 = x*x;
    float y2 = y*y;

    val[i] *= (y2 - x2)*x*y*z * ezr;
    //val[i] *= (y2 - x2)*x*y*z * ezr;
  }
}

void get_6hm3(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float r2 = r*r;

    val[i] *= y * (r2*(3.f*x2-y2) + 9.f*z2*(y2-3.f*x2)) * ezr;
    //val[i] *= (x2*x2 - 10.f*x2*y2 + 5.f*y2*y2)*x * ezr;
  }
}

void get_6hm2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float r2 = r*r;

    val[i] *= x*y*z*(3.f*z2-r2) * ezr;
    //val[i] *= y * (y2-3.f*x2)*(x2 + y2 - 8.f*z2) * ezr;
  }
}

void get_6hm1(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float r2 = r*r;

    val[i] *= y * (14.f*r2*z2 - 21.f*z2*z2 - r2*r2) * ezr;
    //val[i] *= z * (x2*x2 - 6.f*x2*y2 + y2*y2) * ezr;
  }
}

void get_6h0(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float r2 = r*r;

    val[i] *= z*(15.f*r2*r2 - 70.f*r2*z2 + 63.f*z2*z2) * ezr;
    //val[i] *= (x2 + y2 - 2.f*z2)*x*y*z * ezr;
  }
}

void get_6hp1(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float r2 = r*r;

    val[i] *= x*(-r2*r2 + 14.f*r2*z2 - 21.f*z2*z2) * ezr;
    //val[i] *= x * (x2-3.f*y2)*(x2 + y2 - 8.f*z2) * ezr;
  }
}

void get_6hp2(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float r2 = r*r;

    val[i] *= -z * (r2 - 3.f*z2)*(x-y)*(x+y) * ezr;
    //val[i] *= y * (x2*x2 + y2*y2 - 12.f*y2*z2 + 8.f*z2*z2 + 2.f*x2 * (y2-6.f*z2)) * ezr;
  }
}

void get_6hp3(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float r2 = r*r;

    val[i] *= x * (r2*(x2-3.f*y2) + 9.f*z2*(3.f*y2-x2)) * ezr;
    //val[i] *= (x2 - y2) * (x2 + y2 - 2.f*z2) * z * ezr;
  }
}

void get_6hp4(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    float x2 = x*x;
    float y2 = y*y;
    //float z2 = z*z;

    val[i] *= z * (x2*x2 - 6.f*x2*y2 + y2*y2) * ezr;
    //val[i] *= z * (15.f*x2*x2 + 15.f*y2*y2 - 40.f*y2*z2 + 8.f*z2*z2 + 10.f*x2 * (3.f*y2 - 4.f*z2)) * ezr;
  }
}

void get_6hp5(int tid, int gs, float* grid, float* val, float zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    //float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    float x2 = x*x;
    float y2 = y*y;
    //float z2 = z*z;

    val[i] *= x * (10.f*x2*y2 - x2*x2 - 5.f*y2*y2) * ezr;
    //val[i] *= x * (x2*x2 + y2*y2 - 12.f*y2*z2 + 8.f*z2*z2 + 2.f*x2* (y2 - 6.f*z2)) * ezr;
  }
}

// end spherical * exp //


// spherical w/r factors //

void get_px_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float r = grid[6*i+3];
    val[i] *= x/r;
  }
}

void get_py_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    val[i] *= y/r;
  }
}

void get_pz_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    val[i] *= z/r;
  }
}

void get_dxy_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    val[i] *= x*y/r/r;
  }
}

void get_dyz_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    val[i] *= y*z/r/r;
  }
}

void get_dz2_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    val[i] *= 3.f*z*z/r/r-1.f;
  }
}

void get_dxz_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    val[i] *= x*z/r/r;
  }
}

void get_dx2y2_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    val[i] *= (x*x-y*y)/r/r;
  }
}

//Cartesian d functions
void get_dxx_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float r = grid[6*i+3];
    val[i] *= x*x/r/r;
  }
}

void get_dyy_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    val[i] *= y*y/r/r;
  }
}

void get_dzz_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    val[i] *= z*z/r/r;
  }
}

void get_fm3_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float or3 = powf(r,-3.f);
    //val[i] *= (3.f*x*x-y*y)*y*or3;
    //val[i] *= (3.f*x*x-y*y)*y/r/r/r;
    val[i] *= (3.f*x*x-y*y)*y;
    val[i] *= or3;
  }
}

void get_fm2_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or3 = powf(r,-3.f);
    //val[i] *= x*y*z*or3;
    //val[i] *= x*y*z/r/r/r;
    val[i] *= x*y*z;
    val[i] *= or3;
  }
}

void get_fm1_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or3 = powf(r,-3.f);
    //val[i] *= (4.f*z*z-x*x-y*y)*y*or3;
    //val[i] *= (4.f*z*z-x*x-y*y)*y/r/r/r;
    val[i] *= (4.f*z*z-x*x-y*y)*y;
    val[i] *= or3;
  }
}

void get_f0_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or3 = powf(r,-3.f);
    //double valx = (2.*z*z-3.*x*x-3.*y*y)*z*or3;
    //val[i] *= valx;
    val[i] *= (2.f*z*z-3.f*x*x-3.f*y*y)*z*or3;
    //val[i] *= (2.f*z*z-3.f*x*x-3.f*y*y)*z/r/r/r;
    //val[i] *= (2.f*z*z-3.f*x*x-3.f*y*y)*z;
    //val[i] *= or3;
  }
}

void get_fp1_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or3 = powf(r,-3.f);
    //val[i] *= (4.f*z*z-x*x-y*y)*x*or3;
    //val[i] *= (4.f*z*z-x*x-y*y)*x/r/r/r;
    val[i] *= (4.f*z*z-x*x-y*y)*x;
    val[i] *= or3;
  }
}

void get_fp2_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or3 = powf(r,-3.f);
    //val[i] *= (x*x-y*y)*z*or3;
    //val[i] *= (x*x-y*y)*z/r/r/r;
    val[i] *= (x*x-y*y)*z;
    val[i] *= or3;
  }
}

void get_fp3_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float r = grid[6*i+3];
    float or3 = powf(r,-3.f);
    //val[i] *= (x*x-3.f*y*y)*x*or3;
    //val[i] *= (x*x-3.f*y*y)*x/r/r/r;
    val[i] *= (x*x-3.f*y*y)*x;
    val[i] *= or3;
  }
}

void get_gm4_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    //float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or4 = powf(r,-4.f);
    val[i] *= x*y * (x*x - y*y) * or4;
  }
}

void get_gm3_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or4 = powf(r,-4.f);
    val[i] *= y*z * (3.f*x*x - y*y) * or4;
  }
}

void get_gm2_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or4 = powf(r,-4.f);
    val[i] *= x*y * (6.f*z*z - x*x - y*y) * or4;
  }
}

void get_gm1_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or4 = powf(r,-4.f);
    val[i] *= y*z * (4.f*z*z - 3.f*x*x - 3.f*y*y) * or4;
  }
}

void get_g0_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or4 = powf(r,-4.f);
    float z2 = z*z;
    float x2py2 = x*x+y*y;
    float r2 = r*r;
    //float r2 = x*x+y*y+z2;
    //val[i] *= (35.f*z2*z2 - 30.f*z2*r2 + 3.f*r2*r2) * or4;
    val[i] *= (35.f*z2*z2 - 30.f*z2*r2) * or4 + 3.f;
    //val[i] *= (3.f*x2py2*x2py2 - 24.f*x2py2*z2 + 8.f*z2*z2) * or4;
  }
}

void get_gp1_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or4 = powf(r,-4.f);
    val[i] *= x*z * (4.*z*z - 3.f*x*x - 3.f*y*y) * or4;
  }
}

void get_gp2_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or4 = powf(r,-4.f);
    float x2 = x*x;
    float y2 = y*y;
    val[i] *= (x2 - y2) * (6.f*z*z - x2 - y2) * or4;
  }
}

void get_gp3_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or4 = powf(r,-4.f);
    val[i] *= x*z * (x*x - 3.f*y*y) * or4;
  }
}

void get_gp4_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    //float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or4 = powf(r,-4.f);
    float x2 = x*x;
    float y2 = y*y;
    val[i] *= (x2 * (x2 - 3.f*y2) - y2 * (3.f*x2 - y2)) * or4;
  }
}

void get_hm5_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    //float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or5 = powf(r,-5.f);
    float x2 = x*x;
    float y2 = y*y;

    val[i] *= -(5.f*x2*x2 - 10.f*x2*y2 + y2*y2)*y;
    //val[i] *= (5.f*x2*x2 - 10.f*x2*y2 + y2*y2)*y;
    val[i] *= or5;
  }
}

void get_hm4_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or5 = powf(r,-5.f);
    float x2 = x*x;
    float y2 = y*y;

    val[i] *= (y2 - x2)*x*y*z;
    //val[i] *= (y2 - x2)*x*y*z;
    val[i] *= or5;
  }
}

void get_hm3_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or5 = powf(r,-5.f);
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float r2 = r*r;

    val[i] *= y * (r2*(3.f*x2-y2) + 9.f*z2*(y2-3.f*x2));
    //val[i] *= (x2*x2 - 10.f*x2*y2 + 5.f*y2*y2)*x;
    val[i] *= or5;
  }
}

void get_hm2_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or5 = powf(r,-5.f);
    //float x2 = x*x;
    //float y2 = y*y;
    float z2 = z*z;
    float r2 = r*r;

    val[i] *= x*y*z*(3.f*z2-r2);
    //val[i] *= y * (y2-3.f*x2)*(x2 + y2 - 8.f*z2);
    val[i] *= or5;
  }
}

void get_hm1_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    //float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or5 = powf(r,-5.f);
    //float x2 = x*x;
    //float y2 = y*y;
    float z2 = z*z;
    float r2 = r*r;

    val[i] *= y * (14.f*r2*z2 - 21.f*z2*z2 - r2*r2);
    //val[i] *= z * (x2*x2 - 6.f*x2*y2 + y2*y2);
    val[i] *= or5;
  }
}

void get_h0_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or5 = powf(r,-5.f);
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float r2 = r*r;

    val[i] *= z*(15.f*r2*r2 - 70.f*r2*z2 + 63.f*z2*z2);
    //val[i] *= (x2 + y2 - 2.f*z2)*x*y*z;
    val[i] *= or5;
  }
}

void get_hp1_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    //float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or5 = powf(r,-5.f);
    //float x2 = x*x;
    //float y2 = y*y;
    float z2 = z*z;
    float r2 = r*r;

    val[i] *= x*(-r2*r2 + 14.f*r2*z2 - 21.f*z2*z2);
    //val[i] *= x * (x2-3.f*y2)*(x2 + y2 - 8.f*z2);
    val[i] *= or5;
  }
}

void get_hp2_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or5 = powf(r,-5.f);
    //float x2 = x*x;
    //float y2 = y*y;
    float z2 = z*z;
    float r2 = r*r;

    val[i] *= -z * (r2 - 3.f*z2)*(x-y)*(x+y);
    //val[i] *= y * (x2*x2 + y2*y2 - 12.f*y2*z2 + 8.f*z2*z2 + 2.f*x2 * (y2-6.f*z2));
    val[i] *= or5;
  }
}

void get_hp3_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or5 = powf(r,-5.f);
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float r2 = r*r;

    val[i] *= x * (r2*(x2-3.f*y2) + 9.f*z2*(3.f*y2-x2));
    //val[i] *= (x2 - y2) * (x2 + y2 - 2.f*z2) * z;
    val[i] *= or5;
  }
}

void get_hp4_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or5 = powf(r,-5.f);
    float x2 = x*x;
    float y2 = y*y;
    //float z2 = z*z;

    val[i] *= z * (x2*x2 - 6.f*x2*y2 + y2*y2);
    //val[i] *= z * (15.f*x2*x2 + 15.f*y2*y2 - 40.f*y2*z2 + 8.f*z2*z2 + 10.f*x2 * (3.f*y2 - 4.f*z2));
    val[i] *= or5;
  }
}

void get_hp5_3r(int tid, int gs, float* grid, float* val)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    //float z = grid[6*i+2];
    float r = grid[6*i+3];
    float or5 = powf(r,-5.f);
    float x2 = x*x;
    float y2 = y*y;
    //float z2 = z*z;

    val[i] *= x * (10.f*x2*y2 - x2*x2 - 5.f*y2*y2);
    //val[i] *= x * (x2*x2 + y2*y2 - 12.f*y2*z2 + 8.f*z2*z2 + 2.f*x2* (y2 - 6.f*z2));
    val[i] *= or5;
  }
}

// spherical w/r factors //


void eval_sh_3r(int tid, int gs, float* grid, float* val, int n1, int l1, int m1)
{
  if (l1==0) //1s
  {
    return;
  }
  else if (l1==1) //p
  { 
    if (m1==1)
      return get_px_3r(tid,gs,grid,val);
    else if (m1==-1)
      return get_py_3r(tid,gs,grid,val);
    else
      return get_pz_3r(tid,gs,grid,val);
  }
  else if (l1==2) //d
  {
   #if CART_D
    if (m1==0)
      return get_dxx_3r(tid,gs,grid,val);
    else if (m1==1)
      return get_dyy_3r(tid,gs,grid,val);
    else if (m1==2)
      return get_dzz_3r(tid,gs,grid,val);
    else if (m1==3)
      return get_dxy_3r(tid,gs,grid,val);
    else if (m1==4)
      return get_dxz_3r(tid,gs,grid,val);
    else if (m1==5)
      return get_dyz_3r(tid,gs,grid,val);
   #else
    if (m1==-2)
      return get_dxy_3r(tid,gs,grid,val);
    else if (m1==-1)
      return get_dyz_3r(tid,gs,grid,val);
    else if (m1==0)
      return get_dz2_3r(tid,gs,grid,val);
    else if (m1==1)
      return get_dxz_3r(tid,gs,grid,val);
    else
      return get_dx2y2_3r(tid,gs,grid,val);
   #endif
  }
  else if (l1==3) //f
  {
    if (m1==-3)
      return get_fm3_3r(tid,gs,grid,val);
    else if (m1==-2)
      return get_fm2_3r(tid,gs,grid,val);
    else if (m1==-1)
      return get_fm1_3r(tid,gs,grid,val);
    else if (m1== 0)
      return get_f0_3r(tid,gs,grid,val);
    else if (m1== 1)
      return get_fp1_3r(tid,gs,grid,val);
    else if (m1== 2)
      return get_fp2_3r(tid,gs,grid,val);
    else if (m1== 3)
      return get_fp3_3r(tid,gs,grid,val);
  }
  else if (l1==4) //g
  {
    if (m1==-4)
      return get_gm4_3r(tid,gs,grid,val);
    else if (m1==-3)
      return get_gm3_3r(tid,gs,grid,val);
    else if (m1==-2)
      return get_gm2_3r(tid,gs,grid,val);
    else if (m1==-1)
      return get_gm1_3r(tid,gs,grid,val);
    else if (m1== 0)
      return get_g0_3r(tid,gs,grid,val);
    else if (m1== 1)
      return get_gp1_3r(tid,gs,grid,val);
    else if (m1== 2)
      return get_gp2_3r(tid,gs,grid,val);
    else if (m1== 3)
      return get_gp3_3r(tid,gs,grid,val);
    else if (m1== 4)
      return get_gp4_3r(tid,gs,grid,val);
  }
  else if (l1==5) //h
  {
    if (m1==-5)
      return get_hm5_3r(tid,gs,grid,val);
    else if (m1==-4)
      return get_hm4_3r(tid,gs,grid,val);
    else if (m1==-3)
      return get_hm3_3r(tid,gs,grid,val);
    else if (m1==-2)
      return get_hm2_3r(tid,gs,grid,val);
    else if (m1==-1)
      return get_hm1_3r(tid,gs,grid,val);
    else if (m1== 0)
      return get_h0_3r(tid,gs,grid,val);
    else if (m1== 1)
      return get_hp1_3r(tid,gs,grid,val);
    else if (m1== 2)
      return get_hp2_3r(tid,gs,grid,val);
    else if (m1== 3)
      return get_hp3_3r(tid,gs,grid,val);
    else if (m1== 4)
      return get_hp4_3r(tid,gs,grid,val);
    else if (m1== 5)
      return get_hp5_3r(tid,gs,grid,val);
  }

  return;
}


void eval_sh(int tid, int gs, float* grid, float* val, int n1, int l1, int m1, float zeta1)
{
 //1s - 12s
 //2p - 9p
 //3d,4f,5g,6h
 
  if (l1==0)
  {
    if (n1==1)
      get_1s_exp(tid,gs,grid,val,zeta1);
    else if (n1==2)
      get_2s_exp(tid,gs,grid,val,zeta1);
    else if (n1==3)
      get_3s_exp(tid,gs,grid,val,zeta1);
    else if (n1==4)
      get_4s_exp(tid,gs,grid,val,zeta1);
    else if (n1==5)
      get_5s_exp(tid,gs,grid,val,zeta1);
    else if (n1==6)
      get_6s_exp(tid,gs,grid,val,zeta1);
    else if (n1==7)
      get_7s_exp(tid,gs,grid,val,zeta1);
    else if (n1==8)
      get_8s_exp(tid,gs,grid,val,zeta1);
    else if (n1==9)
      get_9s_exp(tid,gs,grid,val,zeta1);
    else if (n1==10)
      get_10s_exp(tid,gs,grid,val,zeta1);
    else if (n1==11)
      get_11s_exp(tid,gs,grid,val,zeta1);
    else if (n1==12)
      get_12s_exp(tid,gs,grid,val,zeta1);
  }
  else if (l1>0)
  {
    if (l1==1) //p
    {
      if (n1==2)
      {
        if (m1==1) //2p
          get_2px(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_2py(tid,gs,grid,val,zeta1);
        else
          get_2pz(tid,gs,grid,val,zeta1);
      }
      else if (n1==3)
      {
        if (m1==1) //3p
          get_3px(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_3py(tid,gs,grid,val,zeta1);
        else
          get_3pz(tid,gs,grid,val,zeta1);
      }
      else if (n1==4)
      {
        if (m1==1) //4p
          get_4px(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_4py(tid,gs,grid,val,zeta1);
        else
          get_4pz(tid,gs,grid,val,zeta1);
      }
      else if (n1==5)
      {
        if (m1==1) //5p
          get_5px(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_5py(tid,gs,grid,val,zeta1);
        else
          get_5pz(tid,gs,grid,val,zeta1);
      }
      else if (n1==6)
      {
        if (m1==1) //6p
          get_6px(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_6py(tid,gs,grid,val,zeta1);
        else
          get_6pz(tid,gs,grid,val,zeta1);
      }
      else if (n1==7)
      {
        if (m1==1) //7p
          get_7px(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_7py(tid,gs,grid,val,zeta1);
        else
          get_7pz(tid,gs,grid,val,zeta1);
      }
      else if (n1==8)
      {
        if (m1==1) //8p
          get_8px(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_8py(tid,gs,grid,val,zeta1);
        else
          get_8pz(tid,gs,grid,val,zeta1);
      }
      else if (n1==9)
      {
        if (m1==1) //9p
          get_9px(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_9py(tid,gs,grid,val,zeta1);
        else
          get_9pz(tid,gs,grid,val,zeta1);
      }
    }
    else if (l1==2) //d
    {
      if (n1==3)
      {
       #if CART_D
        if (m1==0)
          get_3dxx(tid,gs,grid,val,zeta1);
        else if (m1==1)
          get_3dyy(tid,gs,grid,val,zeta1);
        else if (m1==2)
          get_3dzz(tid,gs,grid,val,zeta1);
        else if (m1==3)
          get_3dxy(tid,gs,grid,val,zeta1);
        else if (m1==4)
          get_3dxz(tid,gs,grid,val,zeta1);
        else
          get_3dyz(tid,gs,grid,val,zeta1);
       #else
        if (m1==-2)
          get_3dxy(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_3dyz(tid,gs,grid,val,zeta1);
        else if (m1== 0)
          get_3dz2(tid,gs,grid,val,zeta1);
        else if (m1== 1)
          get_3dxz(tid,gs,grid,val,zeta1);
        else if (m1== 2)
          get_3dx2y2(tid,gs,grid,val,zeta1);
       #endif
      }
      else if (n1==4)
      {
        if (m1==-2)
          get_4dxy(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_4dyz(tid,gs,grid,val,zeta1);
        else if (m1== 0)
          get_4dz2(tid,gs,grid,val,zeta1);
        else if (m1== 1)
          get_4dxz(tid,gs,grid,val,zeta1);
        else if (m1== 2)
          get_4dx2y2(tid,gs,grid,val,zeta1);
      }
      else if (n1==5)
      {
        if (m1==-2)
          get_5dxy(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_5dyz(tid,gs,grid,val,zeta1);
        else if (m1== 0)
          get_5dz2(tid,gs,grid,val,zeta1);
        else if (m1== 1)
          get_5dxz(tid,gs,grid,val,zeta1);
        else if (m1== 2)
          get_5dx2y2(tid,gs,grid,val,zeta1);
      }
      else if (n1==6)
      {
        if (m1==-2)
          get_6dxy(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_6dyz(tid,gs,grid,val,zeta1);
        else if (m1== 0)
          get_6dz2(tid,gs,grid,val,zeta1);
        else if (m1== 1)
          get_6dxz(tid,gs,grid,val,zeta1);
        else if (m1== 2)
          get_6dx2y2(tid,gs,grid,val,zeta1);
      }
      else if (n1==7)
      {
        if (m1==-2)
          get_7dxy(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_7dyz(tid,gs,grid,val,zeta1);
        else if (m1== 0)
          get_7dz2(tid,gs,grid,val,zeta1);
        else if (m1== 1)
          get_7dxz(tid,gs,grid,val,zeta1);
        else if (m1== 2)
          get_7dx2y2(tid,gs,grid,val,zeta1);
      }
    }
    else if (l1==3) //f
    {
      if (n1==4)
      {
        if (m1==-3)
          get_4fm3(tid,gs,grid,val,zeta1);
        else if (m1==-2)
          get_4fm2(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_4fm1(tid,gs,grid,val,zeta1);
        else if (m1== 0)
          get_4f0(tid,gs,grid,val,zeta1);
        else if (m1== 1)
          get_4fp1(tid,gs,grid,val,zeta1);
        else if (m1== 2)
          get_4fp2(tid,gs,grid,val,zeta1);
        else if (m1== 3)
          get_4fp3(tid,gs,grid,val,zeta1);
      }
      else if (n1==5)
      {
        if (m1==-3)
          get_5fm3(tid,gs,grid,val,zeta1);
        else if (m1==-2)
          get_5fm2(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_5fm1(tid,gs,grid,val,zeta1);
        else if (m1== 0)
          get_5f0(tid,gs,grid,val,zeta1);
        else if (m1== 1)
          get_5fp1(tid,gs,grid,val,zeta1);
        else if (m1== 2)
          get_5fp2(tid,gs,grid,val,zeta1);
        else if (m1== 3)
          get_5fp3(tid,gs,grid,val,zeta1);
      }
      else if (n1==6)
      {
        if (m1==-3)
          get_6fm3(tid,gs,grid,val,zeta1);
        else if (m1==-2)
          get_6fm2(tid,gs,grid,val,zeta1);
        else if (m1==-1)
          get_6fm1(tid,gs,grid,val,zeta1);
        else if (m1== 0)
          get_6f0(tid,gs,grid,val,zeta1);
        else if (m1== 1)
          get_6fp1(tid,gs,grid,val,zeta1);
        else if (m1== 2)
          get_6fp2(tid,gs,grid,val,zeta1);
        else if (m1== 3)
          get_6fp3(tid,gs,grid,val,zeta1);
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
        get_5gm4(tid,gs,grid,val,zeta1);
      else if (m1==-3)
        get_5gm3(tid,gs,grid,val,zeta1);
      else if (m1==-2)
        get_5gm2(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        get_5gm1(tid,gs,grid,val,zeta1);
      else if (m1== 0)
        get_5g0(tid,gs,grid,val,zeta1);
      else if (m1== 1)
        get_5gp1(tid,gs,grid,val,zeta1);
      else if (m1== 2)
        get_5gp2(tid,gs,grid,val,zeta1);
      else if (m1== 3)
        get_5gp3(tid,gs,grid,val,zeta1);
      else if (m1== 4)
        get_5gp4(tid,gs,grid,val,zeta1);
    }
    else if (l1==5) //h
    {
     //the ordering here might not be the canonical one
      if (m1==-5)
        get_6hm5(tid,gs,grid,val,zeta1);
      else if (m1==-4)
        get_6hm4(tid,gs,grid,val,zeta1);
      else if (m1==-3)
        get_6hm3(tid,gs,grid,val,zeta1);
      else if (m1==-2)
        get_6hm2(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        get_6hm1(tid,gs,grid,val,zeta1);
      else if (m1== 0)
        get_6h0(tid,gs,grid,val,zeta1);
      else if (m1== 1)
        get_6hp1(tid,gs,grid,val,zeta1);
      else if (m1== 2)
        get_6hp2(tid,gs,grid,val,zeta1);
      else if (m1== 3)
        get_6hp3(tid,gs,grid,val,zeta1);
      else if (m1== 4)
        get_6hp4(tid,gs,grid,val,zeta1);
      else if (m1== 5)
        get_6hp5(tid,gs,grid,val,zeta1);
    }
  } //if l1>0

  return;
}

void eval_sh_3r(int gs, float* grid, float* val, int n1, int l1, int m1)
{
  eval_sh_3r(0,gs,grid,val,n1,l1,m1);
 #if USE_ACC
  //#pragma acc wait
 #endif
}

