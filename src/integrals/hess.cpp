#include "pVp.h"

void get_h_1s(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3]; float r3 = r*r*r;
    float ezr = expf(-zeta*r);
    float ezor = ezr*zeta/r3;

    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float xy = x*y; float xz = x*z; float yz = y*z;
    float zr = zeta*r; float ozr = 1.f+zr;

   //xx, xy, xz, yy, yz, zz 
    val[6*i]   *= (-y2-z2+x2*zr)*ezor; 
    val[6*i+1] *= xy*ozr*ezor;
    val[6*i+2] *= xz*ozr*ezor;
    val[6*i+3] *= (-x2-z2+y2*zr)*ezor;
    val[6*i+4] *= yz*ozr*ezor;
    val[6*i+5] *= (-x2-y2+z2*zr)*ezor;
  }
  return;
}

void get_h_2s(int gs, float* grid, float* val, float zeta)
{
  float zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3]; float r2 = r*r; float r3 = r2*r;
    float ezr = expf(-zeta*r);
    float ezor = ezr/r3;

    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float xy = x*y; float xz = x*z; float yz = y*z;
    float zr = zeta*r; float ozr = 1.f+zr;

   //xx, xy, xz, yy, yz, zz
    val[6*i]   *= (y2+z2+zt2*x2*r2 - zr*(x2+r2))*ezor;
    val[6*i+1] *= xy*(-1.f-zr+zt2*r2)*ezor;
    val[6*i+2] *= xz*(-1.f-zr+zt2*r2)*ezor;
    val[6*i+3] *= (z2+zt2*y2*(y2+z2)-zr*(2.f*y2+z2)+x2*(1.f+zt2*y2-zr))*ezor;
    val[6*i+4] *= yz*(-1.f-zr+zt2*r2)*ezor;
    val[6*i+5] *= (y2+x2*(1.f+zt2*z2-zr) + zt2*z2*(y2+z2) - zr*(y2+2.f*z2))*ezor;
  }
  return;
}

void get_h_2px(int gs, float* grid, float* val, float zeta)
{
  float zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3]; float r3 = r*r*r;
    float ezr = expf(-zeta*r);
    float ezor = ezr*zeta/r3;

    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float xy = x*y; float xz = x*z; float yz = y*z;
    float zr = zeta*r; float ozr = 1.f+zr;

   //xx, xy, xz, yy, yz, zz
    val[6*i]   *= x*(-3.f*(y2+z2)+x2*(-2.f+zr))*ezor;
    val[6*i+1] *= y*(-y2-z2+x2*zr)*ezor;
    val[6*i+2] *= z*(-y2-z2+x2*zr)*ezor;
    val[6*i+3] *= x*(-x2-z2 + y2*zr)*ezor;
    val[6*i+4] *= x*y*z*ozr*ezor;
    val[6*i+5] *= x*(-x2-y2+z2*zr)*ezor;
  }
  return;
}

void get_h_2py(int gs, float* grid, float* val, float zeta)
{
  float zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3]; float r3 = r*r*r;
    float ezr = expf(-zeta*r);
    float ezor = ezr*zeta/r3;

    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float xy = x*y; float xz = x*z; float yz = y*z;
    float zr = zeta*r; float ozr = 1.f+zr;

   //xx, xy, xz, yy, yz, zz
    val[6*i]   *= y*(-y2-z2+x2*zr)*ezor;
    val[6*i+1] *= x*(-x2-z2+y2*zr)*ezor;
    val[6*i+2] *= x*y*z*ozr*ezor;
    val[6*i+3] *= y*(-3.f*x2-3.f*z2+y2*(-2.f+zr))*ezor;
    val[6*i+4] *= z*(-x2-z2+y2*zr)*ezor;
    val[6*i+5] *= y*(-x2-y2+z2*zr)*ezor;
  }
  return;
}

void get_h_2pz(int gs, float* grid, float* val, float zeta)
{
  float zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) //async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3]; float r3 = r*r*r;
    float ezr = expf(-zeta*r);
    float ezor = ezr*zeta/r3;

    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float xy = x*y; float xz = x*z; float yz = y*z;
    float zr = zeta*r; float ozr = 1.f+zr;

   //xx, xy, xz, yy, yz, zz
    val[6*i]   *= z*(-y2-z2+x2*zr)*ezor;
    val[6*i+1] *= x*y*z*ozr*ezor;
    val[6*i+2] *= x*(-x2-y2+z2*zr)*ezor;
    val[6*i+3] *= z*(-x2-z2+y2*zr)*ezor;
    val[6*i+4] *= y*(-x2-y2+z2*zr)*ezor;
    val[6*i+5] *= z*(-3.f*(x2+y2)+z2*(-2.f+zr))*ezor;
  }
  return;
}

void eval_h(int gs, float* grid, float* val, int n1, int l1, int m1, float zeta1)
{
  if (n1==1)
    return get_h_1s(gs,grid,val,zeta1);
  else if (n1==2)
  {
    if (l1==0)
      return get_h_2s(gs,grid,val,zeta1);
    else
    {
      if (m1==0)
        return get_h_2px(gs,grid,val,zeta1);
      else if (m1==1)
        return get_h_2py(gs,grid,val,zeta1);
      else
        return get_h_2pz(gs,grid,val,zeta1);
    }
  }

  return;
}
