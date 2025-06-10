#include "pVp.h"

//sign convention:
// plain derivatives (d/dx,d/dy,d/dz)

//to do//
//1. finish n=4,5
//2. test n=6,7 functions 
//2. 5f (needs testing)

//2. 6h --> need dp (needs testing)

//3. check double vs float

//4. 6f isn't quite right

//not available: 6g,(6h),7f,7g,7h

void get_p_1s(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    float ezr = expf(-zeta*r);
    float ezor = ezr*zeta/r;

    val[3*i]   *= -x*ezor;
    val[3*i+1] *= -y*ezor;
    val[3*i+2] *= -z*ezor;
  }
  return;
}

void get_p_2s(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ezor = ezr*(1.f-zr)/r;

    val[3*i]   *= x*ezor;
    val[3*i+1] *= y*ezor;
    val[3*i+2] *= z*ezor;
  }
  return;
}

void get_p_3s(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ezor = ezr*(2.f-zr);

    val[3*i]   *= x*ezor;
    val[3*i+1] *= y*ezor;
    val[3*i+2] *= z*ezor;
  }
  return;
}

void get_p_4s(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ezor = ezr*(3.f-zr)*r;

    val[3*i]   *= x*ezor;
    val[3*i+1] *= y*ezor;
    val[3*i+2] *= z*ezor;
  }
  return;
}

void get_p_5s(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ezor = ezr*(4.f-zr)*r*r;

    val[3*i]   *= x*ezor;
    val[3*i+1] *= y*ezor;
    val[3*i+2] *= z*ezor;
  }
  return;
}

void get_p_6s(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ezor = ezr*(5.f-zr)*r*r*r;

    val[3*i]   *= x*ezor;
    val[3*i+1] *= y*ezor;
    val[3*i+2] *= z*ezor;
  }
  return;
}

void get_p_7s(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double r2 = r*r;
    double ezor = ezr*(6.f-zr)*r2*r2;

    val[3*i]   *= x*ezor;
    val[3*i+1] *= y*ezor;
    val[3*i+2] *= z*ezor;
  }
  return;
}

void get_p_8s(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double r2 = r*r; double r3 = r2*r;
    double ezor = ezr*(7.f-zr)*r3*r2;

    val[3*i]   *= x*ezor;
    val[3*i+1] *= y*ezor;
    val[3*i+2] *= z*ezor;
  }
  return;
}

void get_p_2px(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    float ezor = ezr/r;
    double rmzx2 = r-zeta*x*x;

    val[3*i]   *= rmzx2*ezor;
    val[3*i+1] *= -zeta*x*y*ezor;
    val[3*i+2] *= -zeta*x*z*ezor;
  }
  return;
}

void get_p_2py(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    float ezor = ezr/r;
    double rmzy2 = r-zeta*y*y;

    val[3*i]   *= -zeta*x*y*ezor;
    val[3*i+1] *= rmzy2*ezor;
    val[3*i+2] *= -zeta*y*z*ezor;
  }
  return;
}

void get_p_2pz(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    float ezor = ezr/r;
    double rmzz2 = r-zeta*z*z;

    val[3*i]   *= -zeta*x*z*ezor;
    val[3*i+1] *= -zeta*y*z*ezor;
    val[3*i+2] *= rmzz2*ezor;
  }
  return;
}

void get_p_3px(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    float ezor = ezr/r;
    double omzr = 1.-zr;

    val[3*i]   *= (y*y+z*z+x*x*(2.-zr))*ezor;
    val[3*i+1] *= x*y*omzr*ezor;
    val[3*i+2] *= x*z*omzr*ezor;
  }
  return;
}

void get_p_3py(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    float ezor = ezr/r;
    double omzr = 1.-zr;

    val[3*i]   *= x*y*omzr*ezor;
    val[3*i+1] *= (x*x+z*z+y*y*(2.-zr))*ezor;
    val[3*i+2] *= y*z*omzr*ezor;
  }
  return;
}

void get_p_3pz(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    float ezor = ezr/r;
    double omzr = 1.-zr;

    val[3*i]   *= x*z*omzr*ezor;
    val[3*i+1] *= y*z*omzr*ezor;
    val[3*i+2] *= (x*x+y*y+z*z*(2.-zr))*ezor;
  }
  return;
}

void get_p_4px(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double tmzr = 2.-zr;

    val[3*i]   *= (y*y+z*z+x*x*(3.-zr))*ezr;
    val[3*i+1] *= x*y*tmzr*ezr;
    val[3*i+2] *= x*z*tmzr*ezr;
  }
  return;
}

void get_p_4py(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double tmzr = 2.-zr;

    val[3*i]   *= x*y*tmzr*ezr;
    val[3*i+1] *= (x*x+z*z+y*y*(3.-zr))*ezr;
    val[3*i+2] *= y*z*tmzr*ezr;
  }
  return;
}

void get_p_4pz(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double tmzr = 2.-zr;

    val[3*i]   *= x*z*tmzr*ezr;
    val[3*i+1] *= y*z*tmzr*ezr;
    val[3*i+2] *= (x*x+y*y+z*z*(3.-zr))*ezr;
  }
  return;
}

void get_p_5px(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double x2 = x*x;
    double r2 = r*r;
    double zr2 = zeta*r2;

    val[3*i]   *= (-x2*zr2+r*(3.*x2+r2))*ezr;
    val[3*i+1] *= -x*y*(zr2-3.*r)*ezr;
    val[3*i+2] *= -x*z*(zr2-3.*r)*ezr;
  }
  return;
}

void get_p_5py(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    float ezor = ezr/r;
    double y2 = y*y;
    double r2 = r*r;
    double zr2 = zeta*r2;

    val[3*i]   *= -x*y*(zr2-3.*r)*ezr;
    val[3*i+1] *= (-y2*zr2+r*(3.*y2+r2))*ezr;
    val[3*i+2] *= -y*z*(zr2-3.*r)*ezr;
  }
  return;
}

void get_p_5pz(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double z2 = z*z;
    double r2 = r*r;
    double zr2 = zeta*r2;

    val[3*i]   *= -x*z*(zr2-3.*r)*ezr;
    val[3*i+1] *= -y*z*(zr2-3.*r)*ezr;
    val[3*i+2] *= (-z2*zr2+r*(3.*z2+r2))*ezr;
  }
  return;
}

void get_p_6px(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double r2 = r*r;
    double zr2 = zeta*r2;

    val[3*i]   *= r2*(y2+z2+x2*(5.-zr))*ezr;
    val[3*i+1] *= x*y*r2*(4.-zr)*ezr;
    val[3*i+2] *= x*z*r2*(4.-zr)*ezr;
  }
  return;
}

void get_p_6py(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    float ezor = ezr/r;
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double r2 = r*r;
    double zr2 = zeta*r2;

    val[3*i]   *= x*y*r2*(4.-zr)*ezr;
    val[3*i+1] *= r2*(x2+z2+y2*(5.-zr))*ezr;
    val[3*i+2] *= y*z*r2*(4.-zr)*ezr;
  }
  return;
}

void get_p_6pz(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double r2 = r*r;
    double zr2 = zeta*r2;

    val[3*i]   *= x*z*r2*(4.-zr)*ezr;
    val[3*i+1] *= y*z*r2*(4.-zr)*ezr;
    val[3*i+2] *= r2*(x2+y2+z2*(5.-zr))*ezr;
  }
  return;
}

void get_p_7px(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double x2 = x*x;
    double r2 = r*r;
    double zr2 = zeta*r2;

    val[3*i]   *= r2*(-zeta*x2*r2 + r*(5.*x2+r2))*ezr;
    val[3*i+1] *= x*y*r2*(5.*r-zr2)*ezr;
    val[3*i+2] *= x*z*r2*(5.*r-zr2)*ezr;
  }
  return;
}

void get_p_7py(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    float ezor = ezr/r;
    double y2 = y*y;
    double r2 = r*r;
    double zr2 = zeta*r2;

    val[3*i]   *= x*y*r2*(5.*r-zr2)*ezr;
    val[3*i+1] *= r2*(-y2*zr2+r*(5.*y2+r2))*ezr;
    val[3*i+2] *= y*z*r2*(5.*r-zr2)*ezr;
  }
  return;
}

void get_p_7pz(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double z2 = z*z;
    double r2 = r*r;
    double zr2 = zeta*r2;

    val[3*i]   *= x*z*r2*(5.*r-zr2)*ezr;
    val[3*i+1] *= y*z*r2*(5.*r-zr2)*ezr;
    val[3*i+2] *= r2*(-z2*zr2+r*(5.*z2+r2))*ezr;
  }
  return;
}

void get_p_3dxy(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)                 
  for (int i=0;i<gs;i++)                                                                  
  {
    float x = grid[6*i];  
    float y = grid[6*i+1]; 
    float z = grid[6*i+2]; 
    float r = grid[6*i+3]; 
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ezor = ezr/r; 

    val[3*i]   *= y*ezor*(r-zeta*x*x);
    val[3*i+1] *= x*ezor*(r-zeta*y*y);
    val[3*i+2] *= -zeta*x*y*z*ezor;
  }
  return;
}

void get_p_3dxz(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ezor = ezr/r;

    val[3*i]   *= z*ezor*(r-zeta*x*x);
    val[3*i+1] *= -zeta*x*y*z*ezor;
    val[3*i+2] *= x*ezor*(r-zeta*z*z);
  }
  return;
}

void get_p_3dyz(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ezor = ezr/r;

    val[3*i]   *= -zeta*x*y*z*ezor;
    val[3*i+1] *= z*ezor*(r-zeta*y*y);
    val[3*i+2] *= y*ezor*(r-zeta*z*z);
  }
  return;
}

void get_p_3dx2y2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ezor = ezr/r;
    double zx2my2 = zeta*(x*x-y*y);

    val[3*i]   *= x*ezor*(2.f*r-zx2my2);
    val[3*i+1] *= -y*ezor*(2.f*r+zx2my2);
    val[3*i+2] *= -z*zx2my2*ezor;
  }
  return;
}

void get_p_3dz2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ezor = ezr/r;
    double ztz2r2 = zeta*(3.f*z*z-r*r);

    val[3*i]   *= -x*ezor*(ztz2r2+2.f*r);
    val[3*i+1] *= -y*ezor*(ztz2r2+2.f*r);
    val[3*i+2] *= z*ezor*(-ztz2r2+4.f*r);
  }
  return;
}

void get_p_4dxy(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)                 
  for (int i=0;i<gs;i++)                                                                  
  {
    float x = grid[6*i];  
    float y = grid[6*i+1]; 
    float z = grid[6*i+2]; 
    float r = grid[6*i+3]; 
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ezor = ezr/r; 
    float x2 = x*x; float y2 = y*y; float z2 = z*z;

    val[3*i]   *= y*ezor*(y2 + z2 + x2*(2.f-zr));
    val[3*i+1] *= x*ezor*(x2 + z2 + y2*(2.f-zr));
    val[3*i+2] *= x*y*z*(1.f-zeta*r)*ezor;
  }
  return;
}

void get_p_4dxz(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ezor = ezr/r;
    float x2 = x*x; float y2 = y*y; float z2 = z*z;

    val[3*i]   *= z*(y2 + z2 + x2*(2.f-zr))*ezor;
    val[3*i+1] *= x*y*z*(1.f-zr)*ezor;
    val[3*i+2] *= x*(x2 + y2 + z2*(2.f-zeta*r))*ezor;
  }
  return;
}

void get_p_4dyz(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ezor = ezr/r;
    float x2 = x*x; float y2 = y*y; float z2 = z*z;

    val[3*i]   *= x*y*z*(1.f-zr)*ezor;
    val[3*i+1] *= z*(x2 + z2 + y2*(2.f-zr))*ezor;
    val[3*i+2] *= y*(x2 + y2 + z2*(2.f-zr))*ezor;
  }
  return;
}

void get_p_4dx2y2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ezor = ezr/r;
    float x2 = x*x; float y2 = y*y; float z2 = z*z;

    val[3*i]   *=  x*(2.f*z2 + x2*(3.f-zr) + y2*(1.f+zr))*ezor;
    val[3*i+1] *= -y*(2.f*z2 + y2*(3.f-zr) + x2*(1.f+zr))*ezor;
    val[3*i+2] *= -z*(x-y)*(x+y)*(-1.f+zr)*ezor;
  }
  return;
}

void get_p_4dz2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ezor = ezr/r;
    double ztz2r2 = zeta*(3.f*z*z-r*r);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float f0 = y2-2.f*z2;
    float f1 = zr*f0;

    val[3*i]   *= x*(-3.f*y2 + f1 + x2*(-3.f+zr))*ezor;
    val[3*i+1] *= y*(-3.f*y2 + f1 + x2*(-3.f+zr))*ezor;
    val[3*i+2] *= z*(f1 + 3.f*(y2+2.f*z2) + x2*(3.f+zr))*ezor;
  }
  return;
}

void get_p_5dxy(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)                 
  for (int i=0;i<gs;i++)                                                                  
  {
    float x = grid[6*i];  
    float y = grid[6*i+1]; 
    float z = grid[6*i+2]; 
    float r = grid[6*i+3]; 
    float zr = zeta*r;
    float ezr = expf(-zr);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;

    val[3*i]   *= y*ezr*(y2 + z2 + x2*(3.f-zr));
    val[3*i+1] *= x*ezr*(x2 + z2 + y2*(3.f-zr));
    val[3*i+2] *= x*y*z*(2.f-zeta*r)*ezr;
  }
  return;
}

void get_p_5dxz(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;

    val[3*i]   *= z*(y2 + z2 + x2*(3.f-zr))*ezr;
    val[3*i+1] *= x*y*z*(2.f-zr)*ezr;
    val[3*i+2] *= x*(x2 + y2 + z2*(3.f-zeta*r))*ezr;
  }
  return;
}

void get_p_5dyz(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;

    val[3*i]   *= x*y*z*(2.f-zr)*ezr;
    val[3*i+1] *= z*(x2 + z2 + y2*(3.f-zr))*ezr;
    val[3*i+2] *= y*(x2 + y2 + z2*(3.f-zr))*ezr;
  }
  return;
}

void get_p_5dx2y2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    //double ezor = ezr/r;
    float x2 = x*x; float y2 = y*y; float z2 = z*z;

    val[3*i]   *=  x*(2.f*z2 + y2*zr + x2*(4.-zr))*ezr;
    val[3*i+1] *= -y*(2.f*z2 + x2*zr + y2*(4.-zr))*ezr;
    val[3*i+2] *=  z*(x-y)*(x+y)*(2.f-zr)*ezr;
  }
  return;
}

void get_p_5dz2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ztz2r2 = zeta*(3.f*z*z-r*r);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float f0 = -4.f*y2 + 2.f*z2;
    float f1 = zr*(y2-2.f*z2);

    val[3*i]   *= x*(f0 + f1 + x2*(-4.f+zr))*ezr;
    val[3*i+1] *= y*(f0 + f1 + x2*(-4.f+zr))*ezr;
    val[3*i+2] *= z*(f1 + 2.f*(y2+4.f*z2) + x2*(2.f+zr))*ezr;
  }
  return;
}

void get_p_6dxy(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)                 
  for (int i=0;i<gs;i++)                                                                  
  {
    float x = grid[6*i];  
    float y = grid[6*i+1]; 
    float z = grid[6*i+2]; 
    float r = grid[6*i+3]; 
    float zr = zeta*r;
    float ezr = expf(-zr);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    double r2 = r*r;

    val[3*i]   *= y*ezr*(-zeta*x2*r2 + r*(r2+3.f*x2));
    val[3*i+1] *= x*ezr*(-zeta*y2*r2 + r*(r2+3.f*y2));
    val[3*i+2] *= x*y*z*(3.f*r-zeta*r2)*ezr;
  }
  return;
}

void get_p_6dxz(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    double r2 = r*r;

    val[3*i]   *= z*(-zeta*x2*r2 + r*(r2+3.f*x2))*ezr;
    val[3*i+1] *= x*y*z*(3.f*r - zeta*r2)*ezr;
    val[3*i+2] *= x*(-zeta*z2*r2 + r*(r2+3.f*z2))*ezr;
  }
  return;
}

void get_p_6dyz(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ezor = ezr/r;
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    double r2 = r*r;

    val[3*i]   *= x*y*z*(3.f*r-zeta*r2)*ezr;
    val[3*i+1] *= z*(-zeta*y2*r2 + r*(r2+3.f*y2))*ezr;
    val[3*i+2] *= y*(-zeta*z2*r2 + r*(r2+3.f*z2))*ezr;
  }
  return;
}

void get_p_6dx2y2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    double r2 = r*r;

    val[3*i]   *= x*((5.f*x2-y2+2.f*z2)*r - zeta*(x-y)*(x+y)*r2)*ezr;
    val[3*i+1] *= y*((x2-5.f*y2-2.f*z2)*r - zeta*(x-y)*(x+y)*r2)*ezr;
    val[3*i+2] *= z*(x-y)*(x+y)*(3.f*r - zeta*r2)*ezr;
  }
  return;
}

void get_p_6dz2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ztz2r2 = zeta*(3.f*z*z-r*r);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    double r2 = r*r;
    float f0 = -5.f*(x2+y2)+4.f*z2;
    float f1 = zeta*(x2+y2-2.f*z2)*r2;

    val[3*i]   *= x*(f1 + r*f0)*ezr;
    val[3*i+1] *= y*(f1 + r*f0)*ezr;
    val[3*i+2] *= z*(f1 + r*(x2+y2+10.f*z2))*ezr;
  }
  return;
}

void get_p_7dxy(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)                 
  for (int i=0;i<gs;i++)                                                                  
  {
    float x = grid[6*i];  
    float y = grid[6*i+1]; 
    float z = grid[6*i+2]; 
    float r = grid[6*i+3]; 
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ezor = ezr/r; 
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    double r2 = r*r;
    float f1 = 5.f-zr;

    val[3*i]   *= y*r2*ezr*(y2 + z2 + x2*f1);
    val[3*i+1] *= x*r2*ezr*(x2 + z2 + y2*f1);
    val[3*i+2] *= x*y*z*r2*(4.f*r-zeta*r2)*ezr;
  }
  return;
}

void get_p_7dxz(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ezor = ezr/r;
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    double r2 = r*r;

    val[3*i]   *= z*r2*(y2 + z2 + x2*(5.f-zr))*ezr;
    val[3*i+1] *= x*y*z*r2*(4.f*r - zeta*r2)*ezr;
    val[3*i+2] *= x*r2*(x2 + y2 + z2*(5.-zr))*ezr;
  }
  return;
}

void get_p_7dyz(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ezor = ezr/r;
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    double r2 = r*r;

    val[3*i]   *= x*y*z*r2*(4.f*r-zeta*r2)*ezr;
    val[3*i+1] *= z*r2*(x2 + z2 + y2*(5.f-zr))*ezr;
    val[3*i+2] *= y*r2*(x2 + y2 + z2*(5.f-zr))*ezr;
  }
  return;
}

void get_p_7dx2y2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    double r2 = r*r;
    float f0 = 2.f-zr;
    float f1 = -6.f+zr;

    val[3*i]   *= x*r2*(2.f*z2 - y2*f0 - x2*f1)*ezr;
    val[3*i+1] *= y*r2*(-2.f*z2 + x2*f0 + y2*f1)*ezr;
    val[3*i+2] *= z*r2*(x-y)*(x+y)*(4.f-zr)*ezr;
  }
  return;
}

void get_p_7dz2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float zr = zeta*r;
    float ezr = expf(-zr);
    double ztz2r2 = zeta*(3.f*z*z-r*r);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    double r2 = r*r;
    float f0 = -6.f+zr;

    val[3*i]   *= x*r2*(6.f*(z-y)*(y+z) + (y2-2.f*z2)*zr + x2*f0)*ezr;
    val[3*i+1] *= y*r2*(6.f*(z-y)*(y+z) + (y2-2.f*z2)*zr + x2*f0)*ezr;
    val[3*i+2] *= z*r2*(zeta*(x2+y2)*r - 2.f*z2*f0)*ezr;
  }
  return;
}

void get_p_4fm3(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x;
    float y2 = y*y;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float zty2m3x2 = (y2-3.f*x2)*zeta;

    val[3*i]   *= x*y*ezor*(zty2m3x2 + 6.f*r);
    val[3*i+1] *= ezor*(zty2m3x2*y2 + 3.f*(x-y)*(x+y)*r);
    val[3*i+2] *= y*z*ezor*zty2m3x2;
  }
  return;
}

void get_p_4fm2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float y2m3x2 = y2-3.f*x2;

    val[3*i]   *= y*z*ezor*(r-zeta*x2);
    val[3*i+1] *= x*z*ezor*(r-zeta*y2);
    val[3*i+2] *= x*y*ezor*(r-zeta*z2);
  }
  return;
}

void get_p_4fm1(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float f1 = x2+y2-4.f*z2;

    val[3*i]   *= x*y*ezor*(zeta*f1-2.f*r);
    val[3*i+1] *= ezor*(zeta*y2*f1 - (x2+3.f*y2-4.f*z2)*r);
    val[3*i+2] *= y*z*ezor*(zeta*f1+8.f*r);
  }
  return;
}

void get_p_4f0(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float f1 = x2+y2;
    float f2 = 3.f*zeta*f1-2.f*zeta*z2;

    val[3*i]   *= x*z*ezor*(f2 - 6.f*r);
    val[3*i+1] *= y*z*ezor*(3.f*zeta*f1 - 2.f*zeta*z2 - 6.f*r);
    val[3*i+2] *= ezor*(f2*z2 - 3.f*(f1-2.f*z2)*r);
  }
  return;
}

void get_p_4fp1(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float f1 = x2+y2-4.f*z2;

    val[3*i]   *= ezor*(zeta*x2*f1 - (3.f*x2+y2-4.f*z2)*r);
    val[3*i+1] *= x*y*ezor*(zeta*f1 - 2.f*r);
    val[3*i+2] *= x*z*ezor*(zeta*f1 + 8.f*r);
  }
  return;
}

void get_p_4fp2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float f1 = zeta*(y2-x2);

    val[3*i]   *= x*z*ezor*(f1+2.f*r);
    val[3*i+1] *= y*z*ezor*(f1-2.f*r);
    val[3*i+2] *= (x-y)*(x+y)*ezor*(-zeta*z2+r);
  }
  return;
}

void get_p_4fp3(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float f1 = x2-3.f*y2;

    val[3*i]   *= ezor*(-zeta*(x2*x2-3.f*x2*y2) + 3.f*(x-y)*(x+y)*r);
    val[3*i+1] *= -x*y*ezor*(zeta*f1 + 6.f*r);
    val[3*i+2] *= -x*z*ezor*(zeta*f1);
  }
  return;
}

void get_p_5fm3(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float zty2m3x2 = (y2-3.f*x2)*zeta;

    val[3*i]   *= x*y*ezor*(6.*z2 + x2*(9.-3.*zr) + y2*(5.+zr));
    val[3*i+1] *= ezor*(3.*x2*x2 - 3.*y2*z2 + y2*y2*(zr-4.) + 3.*x2*(z2+y2*(1.-zr)));
    val[3*i+2] *= y*z*ezor*zty2m3x2*(zr-1.);
  }
  return;
}

void get_p_5fm2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float y2m3x2 = y2-3.f*x2;

    val[3*i]   *= y*z*ezor*(y2+z2 + x2*(2.-zr));
    val[3*i+1] *= x*z*ezor*(x2+z2 + y2*(2.-zr));
    val[3*i+2] *= x*y*ezor*(x2+y2 + z2*(2.-zr));
  }
  return;
}

void get_p_5fm1(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float f1 = x2+y2-4.f*z2;

    val[3*i]   *= x*y*ezor*(-3.*y2 + 2.*z2 + zeta*(y2-4.*z2)*r + x2*(-3.+zr));
    val[3*i+1] *= ezor*(-x2*x2+4.*z2*z2+y2*z2*(5.-4.*zr) + y2*y2*(zr-4.) + x2*(3.*z2+y2*(zr-5.)));
    val[3*i+2] *= y*z*ezor*(7.*y2 + 12.*z2 + zeta*(y2-4.*z2)*r + x2*(7.+zr));
  }
  return;
}

void get_p_5f0(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float f1 = x2+y2;
    float f2 = 3.f*zeta*f1-2.f*zeta*z2;

    val[3*i]   *= x*z*ezor*(-9.*y2 - 4.*z2 + zeta*(3.*y2-2.*z2)*r + 3.*x2*(zr-3.));
    val[3*i+1] *= y*z*ezor*(-9.*y2 - 4.*z2 + zeta*(3.*y2-2.*z2)*r + 3.*x2*(zr-3.));
    val[3*i+2] *= ezor*(-3.*x2*x2-3.*y2*y2 + 3.*zeta*y2*z2*r - 2.*z2*z2*(-4.+zr) + x2*(-6.*y2+3.*zeta*z2*r));
  }
  return;
}

void get_p_5fp1(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float f1 = x2+y2-4.f*z2;

    val[3*i]   *= ezor*(-y2*y2 + 3.*y2*z2 + 4.*z2*z2 + x2*x2*(-4.+zr) + x2*(z2*(5-4.*zr) + y2*(-5.+zr)));
    val[3*i+1] *= x*y*ezor*(-3.*y2 + 2.*z2 + zeta*(y2-4.*z2)*r + x2*(zr-3.));
    val[3*i+2] *= x*z*ezor*(7.*y2 + 12.*z2 + zeta*(y2-4.*z2)*r + x2*(7.+zr));
  }
  return;
}

void get_p_5fp2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float f1 = zeta*(y2-x2);

    val[3*i]   *=  x*z*ezor*(2.*z2 + x2*(3.-zr) + y2*(1.+zr));
    val[3*i+1] *= -y*z*ezor*(2.*z2 + y2*(3.-zr) + x2*(1.+zr));
    val[3*i+2] *=  (x-y)*(x+y)*ezor*(x2+y2 + z2*(2.-zr));
  }
  return;
}

void get_p_5fp3(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float f1 = x2-3.f*y2;

    val[3*i]   *= ezor*(-3.*y2*(y2+z2) + x2*x2*(4.-zr) + 3.*x2*(z2+y2*(zr-1.)));
    val[3*i+1] *= -x*y*ezor*(6.*z2+y2*(9.-3.*zr) + x2*(5.+zr));
    val[3*i+2] *= -x*z*ezor*(x2-3.*y2)*(zr-1.);
  }
  return;
}

void get_p_6fm3(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float zr = zeta*r;
    float ezr = expf(-zr);

    val[3*i]   *= x*y*(6.f*z2 - 3.f*x2*(zr-4.f) + y2*(4+zr))*ezr;
    val[3*i+1] *= (3.f*x2*x2 - 3.f*y2*z2 + y2*y2*(zr-5.f) + 3.f*x2*(z2+y*(2.f-zr)))*ezr;
    val[3*i+2] *= y*z*(y2-3.f*x2)*(zr-2.f)*ezr;
  }
  return;
}

void get_p_6fm2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float zr = zeta*r;
    float ezr = expf(-zr);
    float f1 = 3.f-zr;

    val[3*i]   *= y*z*(y2+z2+x2*f1)*ezr;
    val[3*i+1] *= x*z*(x2+z2+y2*f1)*ezr;
    val[3*i+2] *= x*y*(x2+y2+z2*f1)*ezr;
  }
  return;
}

void get_p_6fm1(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float zr = zeta*r;
    float ezr = expf(-zr);
    float f1 = zeta*(y2-4.f*z2)*r;

    val[3*i]   *= x*y*(-4.f*y2+6.f*z2+f1 + x2*(zr-4.f))*ezr;
    val[3*i+1] *= (-x2*x2 + 4.f*z2*z2 + y2*z2*(9.f-4.*zr) + y2*y2*(zr-5.f) + x2*(3.f*z2+y2*(zr-6.f)))*ezr;
    val[3*i+2] *= y*z*(6.f*y2+16.f*z2 + f1 + x2*(6.f+zr))*ezr;
  }
  return;
}

void get_p_6f0(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float zr = zeta*r;
    float ezr = expf(-zr);
    float f1 = zeta*(3.f*y2-2.f*z2)*r;
    float f2 = -2.f*(6.f*y2+z2) + 3.f*x2*(zr-4.f);

    val[3*i]   *= x*y*(f1 + f2)*ezr;
    val[3*i+1] *= y*z*(f1 + f2)*ezr;
    val[3*i+2] *= (-3.f*x2*x2 - 3.f*y2*y2 - 2.f*z2*z2*(zr-5.f) + 3.f*y2*z2*(zr-1.f) - 3.f*x2*(2.f*y2+z2*(1.-zr)))*ezr;
  }
  return;
}

void get_p_6fp1(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float zr = zeta*r;
    float ezr = expf(-zr);
    float f1 = zeta*(y2-4.f*z2)*r;

    val[3*i]   *= (-y2*y2+3.f*y2*z2+4.f*z2*z2+ x2*x2*(zr-5.f) + x2*(z2*(9.f-4.f*zr)+y2*(zr-6.f)))*ezr;
    val[3*i+1] *= x*y*(-4.f*y2 + 6.f*z2 + f1 + x2*(zr-4.f))*ezr;
    val[3*i+2] *= x*z*(6.f*y2 + 16.f*z2 + f1 + x2*(6.f+zr))*ezr;
  }
  return;
}

void get_p_6fp2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float zr = zeta*r;
    float ezr = expf(-zr);

    val[3*i]   *= x*z*(2.f*z2 + zr*y2 + x2*(4.f-zr))*ezr;
    val[3*i+1] *= y*z*(-2.f*z2 - zr*x2 + y2*(zr-4.f))*ezr;
    val[3*i+2] *= (x-y)*(x+y)*(x2 + y2 + z2*(3.f-zr))*ezr;
  }
  return;
}

void get_p_6fp3(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float zr = zeta*r;
    float ezr = expf(-zr);

    val[3*i]   *= (-3.f*y2*(y2+z2) + x2*x2*(5.f-zr) + 3.f*x2*(z2+y2*(zr-2.f)))*ezr;
    val[3*i+1] *= x*y*(-6.f*z2+3.f*y2*(zr-4.f) + x2*(4.f+zr))*ezr;
    val[3*i+2] *= x*z*(x2-3.f*y2)*(2.f-zr)*ezr;
  }
  return;
}

void get_p_5gm4(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float f0 = zeta*(y2-x2);
    float xmy = 3.f*x2-y2;
    float ymx = 3.f*y2-x2;

    val[3*i]   *= y*(x2*f0 + xmy*r)*ezor;
    val[3*i+1] *= x*(y2*f0 - ymx*r)*ezor;
    val[3*i+2] *= x*y*z*f0*ezor;
  }
  return;
}

void get_p_5gm3(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float f0 = (x-y)*(x+y);
    float xmy = 3.f*x2-y2;

    val[3*i]   *= x*y*z*(-zeta*xmy + 6.f*r)*ezor;
    val[3*i+1] *= z*(-zeta*y2*xmy + 3.*f0*r)*ezor;
    val[3*i+2] *= y*xmy*(-zeta*z2+r)*ezor;
  }
  return;
}

void get_p_5gm2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float f0 = 6.f*z2-x2-y2;
    float f1 = f0-2.f*x2;
    float f2 = f0-2.f*y2;

    val[3*i]   *= y*(-zeta*x2*f0 + f1*r)*ezor;
    val[3*i+1] *= x*(-zeta*y2*f0 + f2*r)*ezor;
    val[3*i+2] *= x*y*z*(-zeta*f0 + 12.f*r)*ezor;
  }
  return;
}

void get_p_5gm1(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float f0 = zeta*(3.f*(x2+y2)-4.f*z2);

    val[3*i]   *= x*y*z*(f0 - 6.f*r)*ezor;
    val[3*i+1] *= z*(y2*f0 + (4.f*z2-3.f*x2-9.f*y2)*r)*ezor;
    val[3*i+2] *= y*(z2*f0 - 3.f*(x2+y2-4.f*z2)*r)*ezor;
  }
  return;
}

void get_p_5g0(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    double x2py2 = x2+y2;
    double xy2 = x2py2*x2py2;
    double f0 = zeta*(-3.f*xy2+24.f*x2py2*z2-8.f*z2*z2);

    //val[3*i]   *= x*(12.f*x2py2*r - 48.f*z2*r + f0)*ezor;
    //val[3*i+1] *= y*(12.f*x2py2*r - 48.f*z2*r + f0)*ezor;
    //val[3*i+2] *= z*(-48.f*x2py2*r + 32.f*z2*r + f0)*ezor;

    val[3*i]   *= x*(12.f*(x2py2-4.f*z2)*r + f0)*ezor;
    val[3*i+1] *= y*(12.f*(x2py2-4.f*z2)*r + f0)*ezor;
    val[3*i+2] *= z*(-16.f*(3.f*x2py2-2.f*z2)*r + f0)*ezor;
  }
  return;
}

void get_p_5gp1(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float x2py2 = x2+y2;
    float f0 = 4.f*z2-9.f*x2-3.*y2;
    float f1 = 3.f*x2py2 - 4.f*z2;

    val[3*i]   *= z*(zeta*x2*f1+r*f0)*ezor;
    val[3*i+1] *= x*y*z*(zeta*f1 - 6.f*r)*ezor;
    val[3*i+2] *= x*(zeta*f1*z2 - 3.f*r*(x2py2-4.f*z2))*ezor;
  }
  return;
}

void get_p_5gp2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float x2py2 = x2+y2;
    float xy2 = x2py2*x2py2;
    float f0 = x2py2-6.f*z2;
    float f1 = (x-y)*(x+y);

    val[3*i]   *= x*(zeta*f1*f0 - 4.f*(x2-3.f*z2)*r)*ezor;
    val[3*i+1] *= y*(zeta*f1*f0 + 4.f*(y2-3.f*z2)*r)*ezor;
    val[3*i+2] *= z*f1*(zeta*f0 + 12.f*r)*ezor;
  }
  return;
}

void get_p_5gp3(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float f0 = x2-3.f*y2;

    val[3*i]   *= z*(-zeta*x2*f0 + 3.f*(x-y)*(x+y)*r)*ezor;
    val[3*i+1] *= -x*y*z*(zeta*f0 + 6.f*r)*ezor;
    val[3*i+2] *= x*f0*(-zeta*z2 + r)*ezor;
  }
  return;
}

void get_p_5gp4(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;
    float f0 = x2*x2-6.f*x2*y2+y2*y2;

    val[3*i]   *= x*(-zeta*f0 + 4.f*(x2-3.f*y2)*r)*ezor;
    val[3*i+1] *= -y*(zeta*f0 + 4.f*(3.f*x2-y2)*r)*ezor;
    val[3*i+2] *= z*zeta*f0*ezor;
  }
  return;
}

void get_p_6hm5(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;
    float zr = zeta*r;
    double ezr = exp(-zr);
    //double ezor = exp(-zr)/r;
    float f1 = 5.f*x4-10.f*x2*y+y4;

    val[3*i]   *= x*y*(20.f*(y2-x2) - (10.f*x2*y2-5.f*x4-y4)*zeta/r) * ezr;
    val[3*i+1] *= (y2*(20.f*x2-4.f*y2) + (-5.f*x4+10.f*x2*y2-y4) - y2*(-5.f*x4+10.f*x2*y2-y4)*zeta/r) * ezr;
    val[3*i+2] *= y*z*(-5.f*x4+10.f*x2*y2-y4)*zeta/r * ezr;
    //val[3*i]   *= -x*y*(zeta*f1 + 20.f*(y2-x2)*r)*ezor;
    //val[3*i+1] *= (-zeta*y2*f1 + 5.f*(x4-6.f*x2*y2+y4)*r)*ezor;
    //val[3*i+2] *= -y*z*f1*ezor;
  }
  return;
}

void get_p_6hm4(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;
    float zr = zeta*r;
    double ezr = exp(-zr);
    //float f0 = (x-y)*(x+y);

    val[3*i]   *= (-2.f*x2*y*z + y*z*(y2-x2) - x2*y*z*(y2-x2)*zeta/r) * ezr;
    val[3*i+1] *= ( 2.f*x2*y*z + x*z*(y2-x2) - x*y2*z*(y2-x2)*zeta/r) * ezr;
    val[3*i+2] *= (x*y*(y2-x2) - x*y*(y2-x2)*z2*zeta/r) * ezr;
    //val[3*i]   *= y*z*(zeta*x2*f0 + (y2-3.f*x2)*r)*ezor;
    //val[3*i+1] *= x*z*(zeta*y2*f0 + (3.f*y2-x2)*r)*ezor;
    //val[3*i+2] *= x*y*f0*(r-zeta*z2)*ezor;
  }
  return;
}

void get_p_6hm3(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z; float r2 = r*r;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;
    float zr = zeta*r;
    double ezr = exp(-zr);
    float f0 = (9.f*(y2-3.f*x2)*z2 + (3.f*x2-y2)*r2)*zeta/r;

    val[3*i]   *= y*(2.f*x*(3.f*x2-y2) - 54.f*x*z2 + 6.f*x*r2 - x*f0) * ezr;
    val[3*i+1] *= (y*(2.f*y*(3.f*x2-y2) + 18.f*y*z2 - 2.f*y*r2) + 9.f*(y2-3.f*x2)*z2 + (3.f*x2-y2)*r2 - y2*f0) * ezr;
    val[3*i+2] *= z*(y*(2.f*(3.f*x2-y2) + 18.f*(y2-3.f*x2)) - y*f0) * ezr;
    //val[3*i]   *= (-zeta*x2*f0 + 5.f*(x4-6.f*x2*y2+y4)*r)*ezor;
    //val[3*i+1] *= x*y*(zeta*f0 + 20.f*(x-y)*(x+y)*r)*ezor;
    //val[3*i+2] *= x*z*f0*ezor;
  }
  return;
}

void get_p_6hm2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;
    float zr = zeta*r;
    double ezr = exp(-zr);
    float f1 = 2.f*z2-x2-y2;
    float f0 = x*y*z*f1*zeta/r;

    val[3*i]   *= (-2.f*x2*y*z + y*z*f1 - x*f0) * ezr;
    val[3*i+1] *= (-2.f*x*y2*z + x*z*f1 - y*f0) * ezr;
    val[3*i+2] *= ( 4.f*x*y*z2 + x*y*f1 - z*f0) * ezr;
    //val[3*i]   *= x*y*(zeta*f2*f1 - 3.f*f0*r)*ezor;
    //val[3*i+1] *= (zeta*y2*f2*f1 + r*(-3.f*x4-6.f*x2*y2+5.f*y4+24.f*(x-y)*(x+y)*z2))*ezor;
    //val[3*i+2] *= y*z*f2*(zeta*f1 + 16.f*r)*ezor;
  }
  return;
}

void get_p_6hm1(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z; float r2 = r*r;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2; float r4 = r2*r2;
    float zr = zeta*r;
    double ezr = exp(-zr);
    float f0 = 28.f*z2-4.*r2;
    float f1 = (-21.f*z4 + 14.f*z2*r2 - r4)*zeta/r;

    val[3*i]   *= (x*f0 - x*f1) * ezr;
    val[3*i+1] *= (y*f0 - y*f1) * ezr;
    val[3*i+2] *= (-56.f*z*z2 + 24.f*z*r2 - z*f1) * ezr;
    //val[3*i]   *= x*x*(zeta*f0-4.f*(x2-3.f*y2)*r)*ezor;
    //val[3*i+1] *= y*z*(zeta*f0+4.f*(3.f*x2-y2)*r)*ezor;
    //val[3*i+2] *= f0*(-zeta*z2+r)*ezor;
  }
  return;
}

void get_p_6h0(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z; float r2 = r*r;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2; float r4 = r2*r2;
    float zr = zeta*r;
    double ezr = exp(-zr);
    float f0 = -140.f*z2 + 60.f*r2;
    float f1 = 112.f*z2 - 80.f*r2;
    float f2 = (63.f*z4 - 70.f*z2*r2 + 15.f*r4)*zeta/r;

    val[3*i]   *= (x*f0 - x*f2)*ezr;
    val[3*i+1] *= (y*f0 - y*f2)*ezr;
    val[3*i+2] *= (z*f1 - z*f2)*ezr;
    //val[3*i]   *= y*z*(-zeta*x2*f0 + (2.f*x2+f0)*r)*ezor;
    //val[3*i+1] *= x*z*(-zeta*y2*f0 + (2.f*y2+f0)*r)*ezor;
    //val[3*i+2] *= x*y*(-zeta*z2*f0 + (-4.f*z2+f0)*r)*ezor;
  }
  return;
}

void get_p_6hp1(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z; float r2 = r*r;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2; float r4 = r2*r2;
    float zr = zeta*r;
    double ezr = exp(-zr);
    float f0 = 28.f*z2 - 4.f*r2;
    float f1 = -56.f*z2 + 24.f*r2;
    float f2 = (-21.f*z4 + 14.f*z2*r2 - r4)*zeta/r;

    val[3*i]   *= (x*f0 - x*f2)*ezr;
    val[3*i+1] *= (y*f0 - y*f2)*ezr;
    val[3*i+2] *= (z*f1 - z*f2)*ezr;
    //val[3*i]   *= (-zeta*x2*f1*f0 + r*(5.f*x4-3.f*y2*(y2-8.f*z2)-6.f*x2*(y2+4.f*z2)))*ezor;
    //val[3*i+1] *= -x*y*(zeta*f1*f0 + 4.f*r*(x2+3.f*(y2-4.f*z2)))*ezor;
    //val[3*i+2] *= -x*z*f1*(zeta*f0+16.f*r)*ezor;
  }
  return;
}

void get_p_6hp2(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z; float r2 = r*r;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;
    float zr = zeta*r;
    double ezr = exp(-zr);
    float f0 = (x+y)*(x-y);
    float f1 = (3.f*z2 - r2);
    float f2 = z*f0*f1*zeta/r;

    val[3*i]   *= (-2.f*x*f0 - (x-y)*z*f1 - (x+y)*z*f1 + x*f2)*ezr;
    val[3*i+1] *= (-2.f*y*f0 - (x-y)*z*f1 + (x+y)*z*f1 + y*f2)*ezr;
    val[3*i+2] *= (f0*z2 - f0*f1 + z*f2)*ezr;
    //val[3*i]   *= -x*y*(-4.f*(f0-6.f*z2)*r + zeta*(f1-12.f*f0*z2+8.f*z4))*ezor;
    //val[3*i+1] *= (zeta*y2*(-f1+12.f*f0*z2-8.f*z4) + r*(x4+6.f*x2*y2+5.f*y4-12.f*(x2+3.f*y2)*z2+8.f*z4))*ezor;
    //val[3*i+2] *= -y*z*((24.f*f0-32.f*z2)*r + zeta*(f1-12.f*f0*z2+8.f*z4))*ezor;
  }
  return;
}

void get_p_6hp3(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z; float r2 = r*r;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;
    float zr = zeta*r;
    double ezr = exp(-zr);
    float f0 = x2-3.f*y2;
    float f1 = (-9.f*f0*z2+f0*r2)*zeta/r;

    val[3*i]   *= (x*(2.f*x*f0 - 18.f*x*z2 + 2.f*x*r2) - 9.f*f0*z2 + f0*r2 - x2*f1) * ezr;
    val[3*i+1] *= (x*(2.f*y*f0 + 54.f*y*z2 - 6.f*y*r2) - x*y*f1) * ezr;
    val[3*i+2] *= (x*(2.f*f0*z - 18.f*f0*z) - x*z*f1) * ezr;
    //val[3*i]   *= x*z*(-zeta*f0*f1 + 4.f*(x-z)*(x+z)*r)*ezor;
    //val[3*i+1] *= y*z*(-zeta*f0*f1 + 4.f*(z2-y2)*r)*ezor;
    //val[3*i+2] *= f0*(-zeta*z2*f1 + (f1-4.f*z2)*r)*ezor;
  }
  return;
}

void get_p_6hp4(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;
    float zr = zeta*r;
    double ezr = exp(-zr);
    float f0 = (x4 - 6.f*x2*y2 + y4)*z*zeta/r;

    val[3*i]   *= ((4.f*x2 - 12.f*y2)*x*z - x*f0) * ezr;
    val[3*i+1] *= ((-12.f*x2 + 4.f*y2)*y*z - y*f0) * ezr;
    val[3*i+2] *= (x4 - 6.f*x2*y2 + y4 - z*f0) * ezr;
    //val[3*i]   *= x*z*(20.f*f1*r + zeta*f2)*ezor;
    //val[3*i+1] *= y*z*(20.f*f1*r + zeta*f2)*ezor;
    //val[3*i+2] *= (zeta*z2*f2 + 5.f*r*(3.f*f0*f0-24.f*f0*z2+8.f*z4))*ezor;
  }
  return;
}

void get_p_6hp5(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;
    float zr = zeta*r;
    double ezr = exp(-zr);
    float f0 = (-x4 + 10.f*x2*y2 - 5.f*y4);
    float f1 = x*f0*zeta/r;

    val[3*i]   *= (x2*(-4.f*x2 + 20.f*y2) + f0 - x*f1) * ezr;
    val[3*i+1] *= (x*y*(20.f*x2 - 20.f*y2) - y*f1) * ezr;
    val[3*i+2] *= -(z*f1) * ezr;
    //val[3*i]   *= (-x2*f1 + r*(5.f*x4+6.f*x2*y2+y4-12.f*(3.f*x2+y2)*z2 + 8.f*z4))*ezor;
    //val[3*i+1] *= x*y*(4.f*(f0-6.f*z2)*r - f1)*ezor;
    //val[3*i+2] *= -x*z*((24.f*f0-32.f*z2)*r + f1)*ezor;
  }
  return;
}

#if 0
//incomplete
void get_p_7hm5(int gs, float* grid, float* val, float zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;
    float zr = zeta*r;
    double ezor = exp(-zr)/r;

    val[3*i]   *= -x*y*(20.f*y2*z2+5.f*x2*x2*(zr-5.f) + y4*(19.f+zr) - 10.f*x2*(2.f*z2+y2*(zr-1.f)))*ezor;
    val[3*i+1] *= (5.f*x2*x4+ 5.f*y4*z2 + y4*y2*(6.f-zr) + 5.f*x4*(z2-y2*(4.+zr)) + 5.f*x2*(-6.y2*z2 + y4*(2.f*zr-6.f*y2*z2)))*ezor;
    val[3*i+2] *= -y*z*(5.f*x4-10.f*x2*y2+y4)*(zr-1.f))*ezor;
  }
  return;
}
#endif

void get_dp_px(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    double xor3 = x/r/r/r;

    val[3*i+0] *= 1./r - x*xor3;
    val[3*i+1] *= -y*xor3;
    val[3*i+2] *= -z*xor3;
  }
  return;
}

void get_dp_py(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    double yor3 = y/r/r/r;

    val[3*i+0] *= -x*yor3;
    val[3*i+1] *= 1./r - y*yor3;
    val[3*i+2] *= -z*yor3;
  }
  return;
}

void get_dp_pz(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    double zor3 = z/r/r/r;

    val[3*i+0] *= -x*zor3;
    val[3*i+1] *= -y*zor3;
    val[3*i+2] *= 1./r - z*zor3;
  }
  return;
}

void get_dp_dxy(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r; double txyrm4 = 2.f*x*y/r2/r2;

    val[3*i+0] *= -x*txyrm4 + y/r2;
    val[3*i+1] *= -y*txyrm4 + x/r2;
    val[3*i+2] *= -z*txyrm4;
  }
  return;
}

void get_dp_dyz(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r; double tyzrm4 = 2.*y*z/r2/r2;

    val[3*i+0] *= -x*tyzrm4;
    val[3*i+1] *= -y*tyzrm4 + z/r2;
    val[3*i+2] *= -z*tyzrm4 + y/r2;
  }
  return;
}

void get_dp_dz2(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r; double szrm4 = 6.*z/r2/r2;

    val[3*i+0] *= -x*z*szrm4;
    val[3*i+1] *= -y*z*szrm4;
    val[3*i+2] *= szrm4*(x*x+y*y);
  }
  return;
}

void get_dp_dxz(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r; double txzrm4 = 2.*x*z/r2/r2;

    val[3*i+0] *= -x*txzrm4 + z/r2;
    val[3*i+1] *= -y*txzrm4;
    val[3*i+2] *= -z*txzrm4 + x/r2;
  }
  return;
}

void get_dp_dx2y2(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r; double tx2y2rm4 = 2.*(x*x-y*y)/r2/r2;

    val[3*i+0] *= -x*tx2y2rm4 + 2.f*x/r2;
    val[3*i+1] *= -y*tx2y2rm4 - 2.f*y/r2;
    val[3*i+2] *= -z*tx2y2rm4;
  }
  return;
}

void get_dp_fm3(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y;
    float rm3 = 1.f/r/r/r;
    float rm5 = rm3/r/r;
    double f0 = 3.f*x2-y2;
    double f1 = 3.f*y*f0*rm5;

    val[3*i+0] *= -x*f1 + 6.f*x*y*rm3;
    val[3*i+1] *= -y*f1 - 2.f*y2*rm3 + f0*rm3;
    val[3*i+2] *= -z*f1;
  }
  return;
}

void get_dp_fm2(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y;
    float rm3 = 1.f/r/r/r;
    float rm5 = rm3/r/r;
    double f1 = 3.f*x*y*z*rm5;

    val[3*i+0] *= -x*f1 + y*z*rm3;
    val[3*i+1] *= -y*f1 + x*z*rm3;
    val[3*i+2] *= -z*f1 + x*y*rm3;
  }
  return;
}

void get_dp_fm1(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y;
    float rm3 = 1.f/r/r/r;
    float rm5 = rm3/r/r;
    double f0 = 4.f*z*z-x*x-y*y;
    double f1 = 3.f*y*f0*rm5;

    val[3*i+0] *= -x*f1 - 2.f*x*y*rm3;
    val[3*i+1] *= -y*f1 - 2.f*y*y*rm3 + f0*rm3;
    val[3*i+2] *= -z*f1 + 8.f*y*z*rm3;
  }
  return;
}

void get_dp_f0(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y;
    float rm3 = 1.f/r/r/r;
    float rm5 = rm3/r/r;
    double f0 = 2.f*z*z-3.f*x*x-3.f*y*y;
    double f1 = 3.f*z*f0*rm5;

    val[3*i+0] *= -x*f1 - 6.f*x*z*rm3;
    val[3*i+1] *= -y*f1 - 6.f*y*z*rm3;
    val[3*i+2] *= -z*f1 + 4.*z*z*rm3 + f0*rm3;
  }
  return;
}

void get_dp_fp1(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y;
    float rm3 = 1.f/r/r/r;
    float rm5 = rm3/r/r;
    double f0 = 4.f*z*z-x*x-y*y;
    double f1 = 3.f*x*f0*rm5;

    val[3*i+0] *= -x*f1 - 2.f*x*x*rm3 + f0*rm3;
    val[3*i+1] *= -y*f1 - 2.f*x*y*rm3;
    val[3*i+2] *= -z*f1 + 8.*x*z*rm3;
  }
  return;
}

void get_dp_fp2(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y;
    float rm3 = 1.f/r/r/r;
    float rm5 = rm3/r/r;
    double f0 = x*x-y*y;
    double f1 = 3.f*z*f0*rm5;

    val[3*i+0] *= -x*f1 + 2.f*x*z*rm3;
    val[3*i+1] *= -y*f1 - 2.f*y*z*rm3;
    val[3*i+2] *= -z*f1 + f0*rm3;
  }
  return;
}

void get_dp_fp3(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    float x2 = x*x; float y2 = y*y;
    float rm3 = 1.f/r/r/r;
    float rm5 = rm3/r/r;
    double f0 = x*x-3.f*y*y;
    double f1 = 3.f*x*f0*rm5;

    val[3*i+0] *= -x*f1 + f0*rm3 + 2.f*x*x*rm3;
    val[3*i+1] *= -y*f1 - 6.f*x*y*rm3;
    val[3*i+2] *= -z*f1;
  }
  return;
}

void get_dp_gm4(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    float rm2 = 1.f/r/r;
    double f0 = (x2-y2)*rm2;
    double xyr = x*y*rm2;
    double f1 = 4.*xyr*f0*rm2;

    val[3*i+0] *= -x*f1 + 2.f*x*xyr*rm2 + y*f0*rm2;
    val[3*i+1] *= -y*f1 - 2.f*y*xyr*rm2 + x*f0*rm2;
    val[3*i+2] *= -z*f1;
  }
  return;
}

void get_dp_gm3(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    float rm2 = 1.f/r/r;
    double f0 = (3.*x2-y2)*rm2;
    double yzr = y*z*rm2;
    double f1 = yzr*f0*rm2;

    val[3*i+0] *= -4.*x*f1 + 6.f*x*yzr*rm2;
    val[3*i+1] *= -4.*y*f1 - 2.f*y*yzr*rm2 + z*f0*rm2;
    val[3*i+2] *= -4.*z*f1 + y*f0*rm2;
  }
  return;
}

void get_dp_gm2(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    float rm2 = 1.f/r/r;
    double f0 = (6.*z2-x2-y2)*rm2;
    double xyr = x*y*rm2;
    double f1 = xyr*f0*rm2;

    val[3*i+0] *= -4.f*x*f1 + 2.f*x*xyr*rm2 + y*f0*rm2;
    val[3*i+1] *= -4.f*y*f1 - 2.f*y*xyr*rm2 + x*f0*rm2;
    val[3*i+2] *= -4.f*z*f1 + 12.f*xyr*z*rm2;
  }
  return;
}

void get_dp_gm1(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double r2 = r*r;
    float rm2 = 1.f/r/r;
    double x2r = x2*rm2; double y2r = y2*rm2; double z2r = z2*rm2;
    double x2y2r = (x2+y2)*rm2;
    //double f0 = 4.*z2-3.*x2-3.*y2;
    double f0 = 7.*z2r-3.;
    double yzr = y*z*rm2;
    double f1 = 4.*yzr*f0*rm2;

    //val[3*i+0] *= -x*f1 - 6.f*x*yz*rm4;
    //val[3*i+1] *= -y*f1 - 6.f*y*yz*rm4 + z*f0*rm4;
    //val[3*i+2] *= -z*f1 + 8.f*z*yz*rm4 + y*f0*rm4;

    val[3*i+0] *= 2.*x*yzr*(3.*x2r+3.*y2r-11.*z2r)*rm2;
    val[3*i+1] *= z*(-3.*(-x2r*x2r+y2r*y2r) + z2r*(x2r-21.*y2r+4.*z2r))*rm2;
    val[3*i+2] *= y*(-3.*x2y2r*x2y2r + 21.*x2y2r*z2r - 4.*z2r*z2r)*rm2;
  }
  return;
}

void get_dp_g0(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    float rm2 = 1.f/r/r;
    double z2r = z2*rm2;
    double x2py2 = x2+y2;
    double x2py2r = x2py2*rm2;
   // double f0 = 3.*x2py2*x2py2 - 24.*x2py2*z2 + 8.*z2*z2;
   // double f1 = 4.*f0*rm6;

   // val[3*i+0] *= -x*f1 + x*(12.f*x2py2 - 48.f*z2)*rm4;
   // val[3*i+1] *= -y*f1 + y*(12.f*x2py2 - 48.f*z2)*rm4;
   // val[3*i+2] *= -z*f1 + (-48.f*z*x2py2 + 32.f*z*z2)*rm4;

    double f2 = (3.*x2py2-4.*z2)*rm2;
    double f3 = 20.*z2r*f2*rm2;

    val[3*i+0] *= x*f3;
    val[3*i+1] *= y*f3;
    val[3*i+2] *= -20.*z*x2py2r*f2*rm2;
  }
  return;
}

void get_dp_gp1(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    float rm2 = 1.f/r/r;
    double xzr = x*z*rm2;
    double f0 = (4.*z2-3.*x2-3.*y2)*rm2;
    double f1 = 4.*xzr*f0*rm2;

    val[3*i+0] *= -x*f1 - 6.f*x*xzr*rm2 + z*f0*rm2;
    val[3*i+1] *= -y*f1 - 6.f*y*xzr*rm2;
    val[3*i+2] *= -z*f1 + 8.f*z*xzr*rm2 + x*f0*rm2;
  }
  return;
}

void get_dp_gp2(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    float rm2 = 1.f/r/r;
    double x2my2r = (x2-y2)*rm2;
    double sz2x2y2r = (6.*z2-x2-y2)*rm2;
    double f0 = x2my2r*sz2x2y2r;
    double f1 = 4.*f0*rm2;

    val[3*i+0] *= -x*f1 - 2.f*x*x2my2r*rm2 + 2.f*x*sz2x2y2r*rm2;
    val[3*i+1] *= -y*f1 - 2.f*y*x2my2r*rm2 - 2.f*y*sz2x2y2r*rm2;
    val[3*i+2] *= -z*f1 + 12.f*z*x2my2r*rm2;
  }
  return;
}

void get_dp_gp3(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    float rm2 = 1.f/r/r;
    double x2r = x2*rm2;
    double xzr = x*z*rm2;
    double x2my2r = (x2-3.*y2)*rm2;
    double f1 = 4.*xzr*x2my2r*rm2;

    val[3*i+0] *= -x*f1 + 2.f*x2r*z*rm2 + z*x2my2r*rm2;
    val[3*i+1] *= -y*f1 - 6.f*y*xzr*rm2;
    val[3*i+2] *= -z*f1 + x*x2my2r*rm2;
  }
  return;
}

void get_dp_gp4(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    float r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    float rm2 = 1.f/r/r;
    double x2m3y2r = (x2-3.*y2)*rm2;
    double y2m3x2r = (y2-3.*x2)*rm2;
    double f0 = x2*x2m3y2r + y2*y2m3x2r;
    double f1 = 4.*f0*rm2*rm2;

  //double check this
    val[3*i+0] *= -x*f1 + 4.f*x*x2m3y2r*rm2;
    val[3*i+1] *= -y*f1 + 4.f*y*y2m3x2r*rm2;
    val[3*i+2] *= -z*f1;
  }
  return;
}


void get_dp_hm5(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;

    val[3*i+0] *= 5.f*x*y*(x4+5.f*y4+4.f*y2*z2-2.f*x2*(5.f*y2+2.f*z2))*rm7;
    val[3*i+1] *= -5.f*(x4*x2-10.f*x4*y2+5.f*x2*y4+(x4-6.f*x2*y2+y4)*z2)*rm7;
    val[3*i+2] *= 5.f*y*(5.f*x4-10.f*x2*y2+y4)*z*rm7;
  }
  return;
}

void get_dp_hm4(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;

    val[3*i+0] *= y*z*(2.f*x4+y2*(y2+z2)-x2*(7.f*y2+3.f*z2))*rm7;
    val[3*i+1] *= -x*z*(x4+2.f*y4-3.f*y2*z2+x2*(-7.f*y2+z2))*rm7;
    val[3*i+2] *= -x*(x-y)*y*(x+y)*(x2+y2-4.f*z2)*rm7;
  }
  return;
}

void get_dp_hm3(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;

    val[3*i+0] *= -3.f*x*y*(x4-3.f*y4+28.f*y2*z2+16.f*z4-2*x2*(y2+14.f*z2))*rm7;
    val[3*i+1] *= 3.f*(x2*(x2-3.f*y2)*(x2+y2)-7.f*(x4-6.f*x2*y2+y4)*z2+8.f*(-x2+y2)*z4)*rm7;
    val[3*i+2] *= 3.f*y*(-3.f*x2+y2)*z*(7.f*(x2+y2)-8.f*z2)*rm7;
  }
  return;
}

void get_dp_hm2(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;

    val[3*i+0] *= y*z*(2*x4-y4+y2*z2+2.f*z4+x2*(y2-11.f*z2))*rm7;
    val[3*i+1] *= x*z*(-x4+2.f*y4-11.f*y2*z2+2.f*z4+x2*(y2+z2))*rm7;
    val[3*i+2] *= x*y*(-(x2+y2)*(x2+y2)+10.f*(x2+y2)*z2-4.f*z4)*rm7;
  }
  return;
}

void get_dp_hm1(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;

    val[3*i+0] *= x*y*((x2+y2)*(x2+y2)-40*(x2+y2)*z2+64.f*z4)*rm7;
    val[3*i+1] *= -x2*(x2+y2)*(x2+y2)+(11.f*x2-29.f*y2)*(x2+y2)*z2+4.f*(x2+17.f*y2)*z4-8.f*z4*z2*rm7;
    val[3*i+2] *= y*z*(29.f*(x2+y2)*(x2+y2)-68.f*(x2+y2)*z2+8.f*z4)*rm7;
  }
  return;
}

void get_dp_h0(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;

    val[3*i+0] *= 15.f*x*z*(-(x2+y2)*(x2+y2)+12.f*(x2+y2)*z2-8.f*z4)*rm7;
    val[3*i+1] *= 15.f*y*z*(-(x2+y2)*(x2+y2)+12.f*(x2+y2)*z2-8.f*z4)*rm7;
    val[3*i+2] *= 15.f*(x2+y2)*((x2+y2)*(x2+y2)-12.f*(x2+y2)*z2+8.f*z4)*rm7;
  }
  return;
}

void get_dp_hp1(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;

    val[3*i+0] *= -y2*(x2+y2)*(x2+y2)+(x2+y2)*(-29.f*x2+11.f*y2)*z2+4.f*(17.f*x2+y2)*z4-8.f*z4*z2*rm7;
    val[3*i+1] *= x*y*((x2+y2)*(x2+y2)-40*(x2+y2)*z2+64.f*z4)*rm7;
    val[3*i+2] *= x*z*(29.f*(x2+y2)*(x2+y2)-68.f*(x2+y2)*z2+8.f*z4)*rm7;
  }
  return;
}

void get_dp_hp2(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;

    val[3*i+0] *= x*z*(x4-5.f*y4+14.f*y2*z2+4.f*z4-2*x2*(2*y2+5.f*z2))*rm7;
    val[3*i+1] *= -y*z*(-5.f*x4-4.f*x2*y2+y4+2.f*(7.f*x2-5.f*y2)*z2+4.f*z4)*rm7;
    val[3*i+2] *= -((x-y)*(x+y)*((x2+y2)*(x2+y2)-10.f*(x2+y2)*z2+4.f*z4))*rm7;
  }
  return;
}

void get_dp_hp3(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;

    val[3*i+0] *= 9.f*x4*y2+6.f*x2*y4-3.f*y4*y2+21*(x4-6.f*x2*y2+y4)*z2+24.f*(-x2+y2)*z4*rm7;
    val[3*i+1] *= 3.f*x*y*(-3.f*x4-2*x2*y2+y4+28.f*(x-y)*(x+y)*z2+16.f*z4)*rm7;
    val[3*i+2] *= -3.f*x*(x2-3.f*y2)*z*(7.f*(x2+y2)-8.f*z2)*rm7;
  }
  return;
}

void get_dp_hp4(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;

    val[3*i+0] *= -x*z*(x4-22.f*x2*y2+17.f*y4-4.f*(x2-3.f*y2)*z2)*rm7;
    val[3*i+1] *= -y*z*(17.f*x4+y4-4.f*y2*z2+x2*(-22.f*y2+12.f*z2))*rm7;
    val[3*i+2] *= (x4-6.f*x2*y2+y4)*(x2+y2-4.f*z2)*rm7;
  }
  return;
}

void get_dp_hp5(int gs, float* grid, float* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i];
    float y = grid[6*i+1];
    float z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    float x2 = x*x; float y2 = y*y; float z2 = z*z;
    float x4 = x2*x2; float y4 = y2*y2; float z4 = z2*z2;

    val[3*i+0] *= -5.f*(5.f*x4*y2-10.f*x2*y4+y4*y2+(x4-6.f*x2*y2+y4)*z2)*rm7;
    val[3*i+1] *= 5.f*x*y*(5.f*x4-10.f*x2*y2+y4+4.f*(x-y)*(x+y)*z2)*rm7;
    val[3*i+2] *= 5.f*x*(x4-10.f*x2*y2+5.f*y4)*z*rm7;
  }
  return;
}

void eval_dp_3r(int gs, float* grid, float* val, int n1, int l1, int m1)
{
  if (l1==0)
  {
    acc_assign(3*gs,val,0.);
  }
  else if (l1==1)
  {
    if (m1==1)
      return get_dp_px(gs,grid,val);
    else if (m1==-1)
      return get_dp_py(gs,grid,val);
    else
      return get_dp_pz(gs,grid,val);
  }
  else if (l1==2)
  {
    if (m1==-2)
      return get_dp_dxy(gs,grid,val);
    else if (m1==-1)
      return get_dp_dyz(gs,grid,val);
    else if (m1== 0)
      return get_dp_dz2(gs,grid,val);
    else if (m1== 1)
      return get_dp_dxz(gs,grid,val);
    else if (m1== 2)
      return get_dp_dx2y2(gs,grid,val);
  }
  else if (l1==3)
  {
    if (m1==-3)
      return get_dp_fm3(gs,grid,val);
    else if (m1==-2)
      return get_dp_fm2(gs,grid,val);
    else if (m1==-1)
      return get_dp_fm1(gs,grid,val);
    else if (m1== 0)
      return get_dp_f0(gs,grid,val);
    else if (m1== 1)
      return get_dp_fp1(gs,grid,val);
    else if (m1== 2)
      return get_dp_fp2(gs,grid,val);
    else if (m1== 3)
      return get_dp_fp3(gs,grid,val);
  }
  else if (l1==4)
  {
    if (m1==-4)
      return get_dp_gm4(gs,grid,val);
    else if (m1==-3)
      return get_dp_gm3(gs,grid,val);
    else if (m1==-2)
      return get_dp_gm2(gs,grid,val);
    else if (m1==-1)
      return get_dp_gm1(gs,grid,val);
    else if (m1== 0)
      return get_dp_g0(gs,grid,val);
    else if (m1== 1)
      return get_dp_gp1(gs,grid,val);
    else if (m1== 2)
      return get_dp_gp2(gs,grid,val);
    else if (m1== 3)
      return get_dp_gp3(gs,grid,val);
    else if (m1== 4)
      return get_dp_gp4(gs,grid,val);
  }
  else if (l1==5)
  { 
    if (m1==-5)
      return get_dp_hm5(gs,grid,val);
    else if (m1==-4)
      return get_dp_hm4(gs,grid,val);
    else if (m1==-3)
      return get_dp_hm3(gs,grid,val);
    else if (m1==-2)
      return get_dp_hm2(gs,grid,val);
    else if (m1==-1)
      return get_dp_hm1(gs,grid,val);
    else if (m1== 0)
      return get_dp_h0(gs,grid,val);
    else if (m1== 1)
      return get_dp_hp1(gs,grid,val);
    else if (m1== 2)
      return get_dp_hp2(gs,grid,val);
    else if (m1== 3)
      return get_dp_hp3(gs,grid,val);
    else if (m1== 4)
      return get_dp_hp4(gs,grid,val);
    else if (m1== 5)
      return get_dp_hp5(gs,grid,val);
  }
  return;
}

void eval_p(int gs, float* grid, float* val, int n1, int l1, int m1, float zeta1)
{
  if (n1==1)
    return get_p_1s(gs,grid,val,zeta1);
  else if (n1==2)
  {
    if (l1==0)
      return get_p_2s(gs,grid,val,zeta1);
    else 
    {
      if (m1==1)
        return get_p_2px(gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_2py(gs,grid,val,zeta1);
      else
        return get_p_2pz(gs,grid,val,zeta1);
    }
  }
  else if (n1==3)
  {
    if (l1==0)
    {
      return get_p_3s(gs,grid,val,zeta1);
    }
    else if (l1==1)
    {
      if (m1==1)
        return get_p_3px(gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_3py(gs,grid,val,zeta1);
      else
        return get_p_3pz(gs,grid,val,zeta1);
    }
    else if (l1==2)
    {
      if (m1==-2)
        return get_p_3dxy(gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_3dyz(gs,grid,val,zeta1);
      else if (m1== 0)
        return get_p_3dz2(gs,grid,val,zeta1);
      else if (m1== 1)
        return get_p_3dxz(gs,grid,val,zeta1);
      else if (m1== 2)
        return get_p_3dx2y2(gs,grid,val,zeta1);
    }
  }
  else if (n1==4)
  {
    if (l1==0)
    {
      return get_p_4s(gs,grid,val,zeta1);
    }
    else if (l1==1)
    {
      if (m1==1)
        return get_p_4px(gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_4py(gs,grid,val,zeta1);
      else
        return get_p_4pz(gs,grid,val,zeta1);
    }
    else if (l1==2)
    {
      //printf(" get_p_4d. m1: %2i \n",m1);
      if (m1==-2)
        return get_p_4dxy(gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_4dyz(gs,grid,val,zeta1);
      else if (m1== 0)
        return get_p_4dz2(gs,grid,val,zeta1);
      else if (m1== 1)
        return get_p_4dxz(gs,grid,val,zeta1);
      else if (m1== 2)
        return get_p_4dx2y2(gs,grid,val,zeta1);
    }
    else if (l1==3)
    {
      if (m1==-3)
        return get_p_4fm3(gs,grid,val,zeta1);
      else if (m1==-2)
        return get_p_4fm2(gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_4fm1(gs,grid,val,zeta1);
      else if (m1== 0)
        return get_p_4f0 (gs,grid,val,zeta1);
      else if (m1== 1)
        return get_p_4fp1(gs,grid,val,zeta1);
      else if (m1== 2)
        return get_p_4fp2(gs,grid,val,zeta1);
      else if (m1== 3)
        return get_p_4fp3(gs,grid,val,zeta1);
    }
  }
  else if (n1==5)
  {
    if (l1==0)
    {
      return get_p_5s(gs,grid,val,zeta1);
    }
    else if (l1==1)
    {
      if (m1==1)
        return get_p_5px(gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_5py(gs,grid,val,zeta1);
      else
        return get_p_5pz(gs,grid,val,zeta1);
    }
    else if (l1==2)
    {
      if (m1==-2)
        return get_p_5dxy(gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_5dyz(gs,grid,val,zeta1);
      else if (m1==0)
        return get_p_5dz2(gs,grid,val,zeta1);
      else if (m1==1)
        return get_p_5dxz(gs,grid,val,zeta1);
      else if (m1==2)
        return get_p_5dx2y2(gs,grid,val,zeta1);
    }
    else if (l1==3)
    {
      if (m1==-3)
        return get_p_5fm3(gs,grid,val,zeta1);
      else if (m1==-2)
        return get_p_5fm2(gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_5fm1(gs,grid,val,zeta1);
      else if (m1==0)
        return get_p_5f0 (gs,grid,val,zeta1);
      else if (m1==1)
        return get_p_5fp1(gs,grid,val,zeta1);
      else if (m1==2)
        return get_p_5fp2(gs,grid,val,zeta1);
      else if (m1==3)
        return get_p_5fp3(gs,grid,val,zeta1);
    }
    else if (l1==4)
    {
      if (m1==-4)
        return get_p_5gm4(gs,grid,val,zeta1);
      else if (m1==-3)
        return get_p_5gm3(gs,grid,val,zeta1);
      else if (m1==-2)
        return get_p_5gm2(gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_5gm1(gs,grid,val,zeta1);
      else if (m1== 0)
        return get_p_5g0(gs,grid,val,zeta1);
      else if (m1== 1)
        return get_p_5gp1(gs,grid,val,zeta1);
      else if (m1== 2)
        return get_p_5gp2(gs,grid,val,zeta1);
      else if (m1== 3)
        return get_p_5gp3(gs,grid,val,zeta1);
      else if (m1== 4)
        return get_p_5gp4(gs,grid,val,zeta1);
    }
  }
  else if (n1==6)
  {
    if (l1==0)
    {
      return get_p_6s(gs,grid,val,zeta1);
    }
    else if (l1==1)
    {
      if (m1==1)
        return get_p_6px(gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_6py(gs,grid,val,zeta1);
      else
        return get_p_6pz(gs,grid,val,zeta1);
    }
    else if (l1==2)
    {
      if (m1==-2)
        return get_p_6dxy(gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_6dyz(gs,grid,val,zeta1);
      else if (m1==0)
        return get_p_6dz2(gs,grid,val,zeta1);
      else if (m1==1)
        return get_p_6dxz(gs,grid,val,zeta1);
      else if (m1==2)
        return get_p_6dx2y2(gs,grid,val,zeta1);
    }
    else if (l1==3)
    {
      if (m1==-3)
        return get_p_6fm3(gs,grid,val,zeta1);
      else if (m1==-2)
        return get_p_6fm2(gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_6fm1(gs,grid,val,zeta1);
      else if (m1==0)
        return get_p_6f0 (gs,grid,val,zeta1);
      else if (m1==1)
        return get_p_6fp1(gs,grid,val,zeta1);
      else if (m1==2)
        return get_p_6fp2(gs,grid,val,zeta1);
      else if (m1==3)
        return get_p_6fp3(gs,grid,val,zeta1);
    }
    else if (l1==4)
    {
      printf(" WARNING: 6g derivatives not available \n");
    }
    else if (l1==5)
    {
      if (m1==-5)
        return get_p_6hm5(gs,grid,val,zeta1);
      else if (m1==-4)
        return get_p_6hm4(gs,grid,val,zeta1);
      else if (m1==-3)
        return get_p_6hm3(gs,grid,val,zeta1);
      else if (m1==-2)
        return get_p_6hm2(gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_6hm1(gs,grid,val,zeta1);
      else if (m1== 0)
        return get_p_6h0(gs,grid,val,zeta1);
      else if (m1== 1)
        return get_p_6hp1(gs,grid,val,zeta1);
      else if (m1== 2)
        return get_p_6hp2(gs,grid,val,zeta1);
      else if (m1== 3)
        return get_p_6hp3(gs,grid,val,zeta1);
      else if (m1== 4)
        return get_p_6hp4(gs,grid,val,zeta1);
      else if (m1== 5)
        return get_p_6hp5(gs,grid,val,zeta1);
    }
  }
  else if (n1==7)
  {
    if (l1==0)
    {
      return get_p_7s(gs,grid,val,zeta1);
    }
    else if (l1==1)
    {
      if (m1==1)
        return get_p_6px(gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_6py(gs,grid,val,zeta1);
      else
        return get_p_6pz(gs,grid,val,zeta1);
    }
    else if (l1==2)
    {
      if (m1==-2)
        return get_p_7dxy(gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_7dyz(gs,grid,val,zeta1);
      else if (m1==0)
        return get_p_7dz2(gs,grid,val,zeta1);
      else if (m1==1)
        return get_p_7dxz(gs,grid,val,zeta1);
      else if (m1==2)
        return get_p_7dx2y2(gs,grid,val,zeta1);
    }
    else if (l1==3)
    {
    }
  }
  else if (n1==8)
  {
    if (l1==0)
    {
      return get_p_8s(gs,grid,val,zeta1);
    }
  }


  return;
}
