#include "pVp.h"

//sign convention:
// plain derivatives (d/dx,d/dy,d/dz)

//tests through 5g look okay
//6h untested

void get_p_1sd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double ezor = ezr*zeta/r;

    val[3*i]   *= -x*ezor;
    val[3*i+1] *= -y*ezor;
    val[3*i+2] *= -z*ezor;
  }
  return;
}

void get_p_2sd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr*(1.-zr)/r;

    val[3*i]   *= x*ezor;
    val[3*i+1] *= y*ezor;
    val[3*i+2] *= z*ezor;
  }
  return;
}

void get_p_3sd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr*(2.-zr);

    val[3*i]   *= x*ezor;
    val[3*i+1] *= y*ezor;
    val[3*i+2] *= z*ezor;
  }
  return;
}

void get_p_4sd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr*(3.-zr)*r;

    val[3*i]   *= x*ezor;
    val[3*i+1] *= y*ezor;
    val[3*i+2] *= z*ezor;
  }
  return;
}

void get_p_5sd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr*(4.-zr)*r*r;

    val[3*i]   *= x*ezor;
    val[3*i+1] *= y*ezor;
    val[3*i+2] *= z*ezor;
  }
  return;
}

void get_p_6sd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr*(5.-zr)*r*r*r;

    val[3*i]   *= x*ezor;
    val[3*i+1] *= y*ezor;
    val[3*i+2] *= z*ezor;
  }
  return;
}

void get_p_7sd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double r2 = r*r;
    double ezor = ezr*(6.-zr)*r2*r2;

    val[3*i]   *= x*ezor;
    val[3*i+1] *= y*ezor;
    val[3*i+2] *= z*ezor;
  }
  return;
}

void get_p_8sd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double r2 = r*r; double r3 = r2*r;
    double ezor = ezr*(7.-zr)*r3*r2;

    val[3*i]   *= x*ezor;
    val[3*i+1] *= y*ezor;
    val[3*i+2] *= z*ezor;
  }
  return;
}

void get_p_2pxd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;
    double rmzx2 = r-zeta*x*x;

    val[3*i]   *= rmzx2*ezor;
    val[3*i+1] *= -zeta*x*y*ezor;
    val[3*i+2] *= -zeta*x*z*ezor;
  }
  return;
}

void get_p_2pyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;
    double rmzy2 = r-zeta*y*y;

    val[3*i]   *= -zeta*x*y*ezor;
    val[3*i+1] *= rmzy2*ezor;
    val[3*i+2] *= -zeta*y*z*ezor;
  }
  return;
}

void get_p_2pzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;
    double rmzz2 = r-zeta*z*z;

    val[3*i]   *= -zeta*x*z*ezor;
    val[3*i+1] *= -zeta*y*z*ezor;
    val[3*i+2] *= rmzz2*ezor;
  }
  return;
}

void get_p_3pxd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;
    double omzr = 1.-zr;

    val[3*i]   *= (y*y+z*z+x*x*(2.-zr))*ezor;
    val[3*i+1] *= x*y*omzr*ezor;
    val[3*i+2] *= x*z*omzr*ezor;
  }
  return;
}

void get_p_3pyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;
    double omzr = 1.-zr;

    val[3*i]   *= x*y*omzr*ezor;
    val[3*i+1] *= (x*x+z*z+y*y*(2.-zr))*ezor;
    val[3*i+2] *= y*z*omzr*ezor;
  }
  return;
}

void get_p_3pzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;
    double omzr = 1.-zr;

    val[3*i]   *= x*z*omzr*ezor;
    val[3*i+1] *= y*z*omzr*ezor;
    val[3*i+2] *= (x*x+y*y+z*z*(2.-zr))*ezor;
  }
  return;
}

void get_p_4pxd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double tmzr = 2.-zr;

    val[3*i]   *= (y*y+z*z+x*x*(3.-zr))*ezr;
    val[3*i+1] *= x*y*tmzr*ezr;
    val[3*i+2] *= x*z*tmzr*ezr;
  }
  return;
}

void get_p_4pyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double tmzr = 2.-zr;

    val[3*i]   *= x*y*tmzr*ezr;
    val[3*i+1] *= (x*x+z*z+y*y*(3.-zr))*ezr;
    val[3*i+2] *= y*z*tmzr*ezr;
  }
  return;
}

void get_p_4pzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double tmzr = 2.-zr;

    val[3*i]   *= x*z*tmzr*ezr;
    val[3*i+1] *= y*z*tmzr*ezr;
    val[3*i+2] *= (x*x+y*y+z*z*(3.-zr))*ezr;
  }
  return;
}

void get_p_5pxd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double x2 = x*x;
    double r2 = r*r;
    double zr2 = zeta*r2;

    val[3*i]   *= (-x2*zr2+r*(3.*x2+r2))*ezr;
    val[3*i+1] *= -x*y*(zr2-3.*r)*ezr;
    val[3*i+2] *= -x*z*(zr2-3.*r)*ezr;
  }
  return;
}

void get_p_5pyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;
    double y2 = y*y;
    double r2 = r*r;
    double zr2 = zeta*r2;

    val[3*i]   *= -x*y*(zr2-3.*r)*ezr;
    val[3*i+1] *= (-y2*zr2+r*(3.*y2+r2))*ezr;
    val[3*i+2] *= -y*z*(zr2-3.*r)*ezr;
  }
  return;
}

void get_p_5pzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double z2 = z*z;
    double r2 = r*r;
    double zr2 = zeta*r2;

    val[3*i]   *= -x*z*(zr2-3.*r)*ezr;
    val[3*i+1] *= -y*z*(zr2-3.*r)*ezr;
    val[3*i+2] *= (-z2*zr2+r*(3.*z2+r2))*ezr;
  }
  return;
}

void get_p_6pxd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double r2 = r*r;
    double zr2 = zeta*r2;

    val[3*i]   *= r2*(y2+z2+x2*(5.-zr))*ezr;
    val[3*i+1] *= x*y*r2*(4.-zr)*ezr;
    val[3*i+2] *= x*z*r2*(4.-zr)*ezr;
  }
  return;
}

void get_p_6pyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double r2 = r*r;
    double zr2 = zeta*r2;

    val[3*i]   *= x*y*r2*(4.-zr)*ezr;
    val[3*i+1] *= r2*(x2+z2+y2*(5.-zr))*ezr;
    val[3*i+2] *= y*z*r2*(4.-zr)*ezr;
  }
  return;
}

void get_p_6pzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double r2 = r*r;
    double zr2 = zeta*r2;

    val[3*i]   *= x*z*r2*(4.-zr)*ezr;
    val[3*i+1] *= y*z*r2*(4.-zr)*ezr;
    val[3*i+2] *= r2*(x2+y2+z2*(5.-zr))*ezr;
  }
  return;
}

void get_p_7pxd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double x2 = x*x;
    double r2 = r*r;
    double zr2 = zeta*r2;

    val[3*i]   *= r2*(-zeta*x2*r2 + r*(5.*x2+r2))*ezr;
    val[3*i+1] *= x*y*r2*(5.*r-zr2)*ezr;
    val[3*i+2] *= x*z*r2*(5.*r-zr2)*ezr;
  }
  return;
}

void get_p_7pyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;
    double y2 = y*y;
    double r2 = r*r;
    double zr2 = zeta*r2;

    val[3*i]   *= x*y*r2*(5.*r-zr2)*ezr;
    val[3*i+1] *= r2*(-y2*zr2+r*(5.*y2+r2))*ezr;
    val[3*i+2] *= y*z*r2*(5.*r-zr2)*ezr;
  }
  return;
}

void get_p_7pzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double z2 = z*z;
    double r2 = r*r;
    double zr2 = zeta*r2;

    val[3*i]   *= x*z*r2*(5.*r-zr2)*ezr;
    val[3*i+1] *= y*z*r2*(5.*r-zr2)*ezr;
    val[3*i+2] *= r2*(-z2*zr2+r*(5.*z2+r2))*ezr;
  }
  return;
}

void get_p_3dxyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)                 
  for (int i=0;i<gs;i++)                                                                  
  {
    double x = grid[6*i];  
    double y = grid[6*i+1]; 
    double z = grid[6*i+2]; 
    double r = grid[6*i+3]; 
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r; 

    val[3*i]   *= y*ezor*(r-zeta*x*x);
    val[3*i+1] *= x*ezor*(r-zeta*y*y);
    val[3*i+2] *= -zeta*x*y*z*ezor;
  }
  return;
}

void get_p_3dxzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;

    val[3*i]   *= z*ezor*(r-zeta*x*x);
    val[3*i+1] *= -zeta*x*y*z*ezor;
    val[3*i+2] *= x*ezor*(r-zeta*z*z);
  }
  return;
}

void get_p_3dyzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;

    val[3*i]   *= -zeta*x*y*z*ezor;
    val[3*i+1] *= z*ezor*(r-zeta*y*y);
    val[3*i+2] *= y*ezor*(r-zeta*z*z);
  }
  return;
}

void get_p_3dx2y2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;
    double zx2my2 = zeta*(x*x-y*y);

    val[3*i]   *= x*ezor*(2.*r-zx2my2);
    val[3*i+1] *= -y*ezor*(2.*r+zx2my2);
    val[3*i+2] *= -z*zx2my2*ezor;
  }
  return;
}

void get_p_3dz2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;
    double ztz2r2 = zeta*(3.*z*z-r*r);

    val[3*i]   *= -x*ezor*(ztz2r2+2.*r);
    val[3*i+1] *= -y*ezor*(ztz2r2+2.*r);
    val[3*i+2] *= z*ezor*(-ztz2r2+4.*r);
  }
  return;
}

void get_p_4dxyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)                 
  for (int i=0;i<gs;i++)                                                                  
  {
    double x = grid[6*i];  
    double y = grid[6*i+1]; 
    double z = grid[6*i+2]; 
    double r = grid[6*i+3]; 
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r; 
    double x2 = x*x; double y2 = y*y; double z2 = z*z;

    val[3*i]   *= y*ezor*(y2 + z2 + x2*(2.-zr));
    val[3*i+1] *= x*ezor*(x2 + z2 + y2*(2.-zr));
    val[3*i+2] *= x*y*z*(1.-zeta*r)*ezor;
  }
  return;
}

void get_p_4dxzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;
    double x2 = x*x; double y2 = y*y; double z2 = z*z;

    val[3*i]   *= z*(y2 + z2 + x2*(2.-zr))*ezor;
    val[3*i+1] *= x*y*z*(1.-zr)*ezor;
    val[3*i+2] *= x*(x2 + y2 + z2*(2.-zeta*r))*ezor;
  }
  return;
}

void get_p_4dyzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;
    double x2 = x*x; double y2 = y*y; double z2 = z*z;

    val[3*i]   *= x*y*z*(1.-zr)*ezor;
    val[3*i+1] *= z*(x2 + z2 + y2*(2.-zr))*ezor;
    val[3*i+2] *= y*(x2 + y2 + z2*(2.-zr))*ezor;
  }
  return;
}

void get_p_4dx2y2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;
    double x2 = x*x; double y2 = y*y; double z2 = z*z;

    val[3*i]   *=  x*(2.*z2 + x2*(3.-zr) + y2*(1.+zr))*ezor;
    val[3*i+1] *= -y*(2.*z2 + y2*(3.-zr) + x2*(1.+zr))*ezor;
    val[3*i+2] *= -z*(x-y)*(x+y)*(-1.+zr)*ezor;
  }
  return;
}

void get_p_4dz2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;
    double ztz2r2 = zeta*(3.*z*z-r*r);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double f0 = y2-2.*z2;
    double f1 = zr*f0;

    val[3*i]   *= x*(-3.*y2 + f1 + x2*(-3.+zr))*ezor;
    val[3*i+1] *= y*(-3.*y2 + f1 + x2*(-3.+zr))*ezor;
    val[3*i+2] *= z*(f1 + 3.*(y2+2.*z2) + x2*(3.+zr))*ezor;
  }
  return;
}

void get_p_5dxyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)                 
  for (int i=0;i<gs;i++)                                                                  
  {
    double x = grid[6*i];  
    double y = grid[6*i+1]; 
    double z = grid[6*i+2]; 
    double r = grid[6*i+3]; 
    double zr = zeta*r;
    double ezr = exp(-zr);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;

    val[3*i]   *= y*ezr*(y2 + z2 + x2*(3.-zr));
    val[3*i+1] *= x*ezr*(x2 + z2 + y2*(3.-zr));
    val[3*i+2] *= x*y*z*(2.-zeta*r)*ezr;
  }
  return;
}

void get_p_5dxzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;

    val[3*i]   *= z*(y2 + z2 + x2*(3.-zr))*ezr;
    val[3*i+1] *= x*y*z*(2.-zr)*ezr;
    val[3*i+2] *= x*(x2 + y2 + z2*(3.-zeta*r))*ezr;
  }
  return;
}

void get_p_5dyzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;

    val[3*i]   *= x*y*z*(2.-zr)*ezr;
    val[3*i+1] *= z*(x2 + z2 + y2*(3.-zr))*ezr;
    val[3*i+2] *= y*(x2 + y2 + z2*(3.-zr))*ezr;
  }
  return;
}

void get_p_5dx2y2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    //double ezor = ezr/r;
    double x2 = x*x; double y2 = y*y; double z2 = z*z;

    val[3*i]   *=  x*(2.*z2 + y2*zr + x2*(4.-zr))*ezr;
    val[3*i+1] *= -y*(2.*z2 + x2*zr + y2*(4.-zr))*ezr;
    val[3*i+2] *=  z*(x-y)*(x+y)*(2.-zr)*ezr;
  }
  return;
}

void get_p_5dz2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ztz2r2 = zeta*(3.*z*z-r*r);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double f0 = -4.*y2 + 2.*z2;
    double f1 = zr*(y2-2.*z2);

    val[3*i]   *= x*(f0 + f1 + x2*(-4.+zr))*ezr;
    val[3*i+1] *= y*(f0 + f1 + x2*(-4.+zr))*ezr;
    val[3*i+2] *= z*(f1 + 2.*(y2+4.*z2) + x2*(2.+zr))*ezr;
  }
  return;
}

void get_p_6dxyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)                 
  for (int i=0;i<gs;i++)                                                                  
  {
    double x = grid[6*i];  
    double y = grid[6*i+1]; 
    double z = grid[6*i+2]; 
    double r = grid[6*i+3]; 
    double zr = zeta*r;
    double ezr = exp(-zr);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double r2 = r*r;

    val[3*i]   *= y*ezr*(-zeta*x2*r2 + r*(r2+3.*x2));
    val[3*i+1] *= x*ezr*(-zeta*y2*r2 + r*(r2+3.*y2));
    val[3*i+2] *= x*y*z*(3.*r-zeta*r2)*ezr;
  }
  return;
}

void get_p_6dxzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double r2 = r*r;

    val[3*i]   *= z*(-zeta*x2*r2 + r*(r2+3.*x2))*ezr;
    val[3*i+1] *= x*y*z*(3.*r - zeta*r2)*ezr;
    val[3*i+2] *= x*(-zeta*z2*r2 + r*(r2+3.*z2))*ezr;
  }
  return;
}

void get_p_6dyzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double r2 = r*r;

    val[3*i]   *= x*y*z*(3.*r-zeta*r2)*ezr;
    val[3*i+1] *= z*(-zeta*y2*r2 + r*(r2+3.*y2))*ezr;
    val[3*i+2] *= y*(-zeta*z2*r2 + r*(r2+3.*z2))*ezr;
  }
  return;
}

void get_p_6dx2y2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double r2 = r*r;

    val[3*i]   *= x*((5.*x2-y2+2.*z2)*r - zeta*(x-y)*(x+y)*r2)*ezr;
    val[3*i+1] *= y*((x2-5.*y2-2.*z2)*r - zeta*(x-y)*(x+y)*r2)*ezr;
    val[3*i+2] *= z*(x-y)*(x+y)*(3.*r - zeta*r2)*ezr;
  }
  return;
}

void get_p_6dz2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ztz2r2 = zeta*(3.*z*z-r*r);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double r2 = r*r;
    double f0 = -5.*(x2+y2)+4.*z2;
    double f1 = zeta*(x2+y2-2.*z2)*r2;

    val[3*i]   *= x*(f1 + r*f0)*ezr;
    val[3*i+1] *= y*(f1 + r*f0)*ezr;
    val[3*i+2] *= z*(f1 + r*(x2+y2+10.*z2))*ezr;
  }
  return;
}

void get_p_7dxyd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)                 
  for (int i=0;i<gs;i++)                                                                  
  {
    double x = grid[6*i];  
    double y = grid[6*i+1]; 
    double z = grid[6*i+2]; 
    double r = grid[6*i+3]; 
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r; 
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double r2 = r*r;
    double f1 = 5.-zr;

    val[3*i]   *= y*r2*ezr*(y2 + z2 + x2*f1);
    val[3*i+1] *= x*r2*ezr*(x2 + z2 + y2*f1);
    val[3*i+2] *= x*y*z*r2*(4.*r-zeta*r2)*ezr;
  }
  return;
}

void get_p_7dxzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double r2 = r*r;

    val[3*i]   *= z*r2*(y2 + z2 + x2*(5.-zr))*ezr;
    val[3*i+1] *= x*y*z*r2*(4.*r - zeta*r2)*ezr;
    val[3*i+2] *= x*r2*(x2 + y2 + z2*(5.-zr))*ezr;
  }
  return;
}

void get_p_7dyzd(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ezor = ezr/r;
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double r2 = r*r;

    val[3*i]   *= x*y*z*r2*(4.*r-zeta*r2)*ezr;
    val[3*i+1] *= z*r2*(x2 + z2 + y2*(5.-zr))*ezr;
    val[3*i+2] *= y*r2*(x2 + y2 + z2*(5.-zr))*ezr;
  }
  return;
}

void get_p_7dx2y2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double r2 = r*r;
    double f0 = 2.-zr;
    double f1 = -6.+zr;

    val[3*i]   *= x*r2*(2.*z2 - y2*f0 - x2*f1)*ezr;
    val[3*i+1] *= y*r2*(-2.*z2 + x2*f0 + y2*f1)*ezr;
    val[3*i+2] *= z*r2*(x-y)*(x+y)*(4.-zr)*ezr;
  }
  return;
}

void get_p_7dz2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zr = zeta*r;
    double ezr = exp(-zr);
    double ztz2r2 = zeta*(3.*z*z-r*r);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double r2 = r*r;
    double f0 = -6.+zr;

    val[3*i]   *= x*r2*(6.*(z-y)*(y+z) + (y2-2.*z2)*zr + x2*f0)*ezr;
    val[3*i+1] *= y*r2*(6.*(z-y)*(y+z) + (y2-2.*z2)*zr + x2*f0)*ezr;
    val[3*i+2] *= z*r2*(zeta*(x2+y2)*r - 2.*z2*f0)*ezr;
  }
  return;
}

void get_p_4fm3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x;
    double y2 = y*y;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double zty2m3x2 = (y2-3.*x2)*zeta;

    val[3*i]   *= x*y*ezor*(zty2m3x2 + 6.*r);
    val[3*i+1] *= ezor*(zty2m3x2*y2 + 3.*(x-y)*(x+y)*r);
    val[3*i+2] *= y*z*ezor*zty2m3x2;
  }
  return;
}

void get_p_4fm2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double y2m3x2 = y2-3.*x2;

    val[3*i]   *= y*z*ezor*(r-zeta*x2);
    val[3*i+1] *= x*z*ezor*(r-zeta*y2);
    val[3*i+2] *= x*y*ezor*(r-zeta*z2);
  }
  return;
}

void get_p_4fm1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double f1 = x2+y2-4.*z2;

    val[3*i]   *= x*y*ezor*(zeta*f1-2.*r);
    val[3*i+1] *= ezor*(zeta*y2*f1 - (x2+3.*y2-4.*z2)*r);
    val[3*i+2] *= y*z*ezor*(zeta*f1+8.*r);
  }
  return;
}

void get_p_4f0d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double f1 = x2+y2;
    double f2 = 3.*zeta*f1-2.*zeta*z2;

    val[3*i]   *= x*z*ezor*(f2 - 6.*r);
    val[3*i+1] *= y*z*ezor*(3.*zeta*f1 - 2.*zeta*z2 - 6.*r);
    val[3*i+2] *= ezor*(f2*z2 - 3.*(f1-2.*z2)*r);
  }
  return;
}

void get_p_4fp1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double f1 = x2+y2-4.*z2;

    val[3*i]   *= ezor*(zeta*x2*f1 - (3.*x2+y2-4.*z2)*r);
    val[3*i+1] *= x*y*ezor*(zeta*f1 - 2.*r);
    val[3*i+2] *= x*z*ezor*(zeta*f1 + 8.*r);
  }
  return;
}

void get_p_4fp2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double f1 = zeta*(y2-x2);

    val[3*i]   *= x*z*ezor*(f1+2.*r);
    val[3*i+1] *= y*z*ezor*(f1-2.*r);
    val[3*i+2] *= (x-y)*(x+y)*ezor*(-zeta*z2+r);
  }
  return;
}

void get_p_4fp3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double f1 = x2-3.*y2;

    val[3*i]   *= ezor*(-zeta*(x2*x2-3.*x2*y2) + 3.*(x-y)*(x+y)*r);
    val[3*i+1] *= -x*y*ezor*(zeta*f1 + 6.*r);
    val[3*i+2] *= -x*z*ezor*(zeta*f1);
  }
  return;
}

void get_p_5fm3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double zty2m3x2 = (y2-3.*x2)*zeta;

    val[3*i]   *= x*y*ezor*(6.*z2 + x2*(9.-3.*zr) + y2*(5.+zr));
    val[3*i+1] *= ezor*(3.*x2*x2 - 3.*y2*z2 + y2*y2*(zr-4.) + 3.*x2*(z2+y2*(1.-zr)));
    val[3*i+2] *= y*z*ezor*zty2m3x2*(zr-1.);
  }
  return;
}

void get_p_5fm2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double y2m3x2 = y2-3.*x2;

    val[3*i]   *= y*z*ezor*(y2+z2 + x2*(2.-zr));
    val[3*i+1] *= x*z*ezor*(x2+z2 + y2*(2.-zr));
    val[3*i+2] *= x*y*ezor*(x2+y2 + z2*(2.-zr));
  }
  return;
}

void get_p_5fm1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double f1 = x2+y2-4.*z2;

    val[3*i]   *= x*y*ezor*(-3.*y2 + 2.*z2 + zeta*(y2-4.*z2)*r + x2*(-3.+zr));
    val[3*i+1] *= ezor*(-x2*x2+4.*z2*z2+y2*z2*(5.-4.*zr) + y2*y2*(zr-4.) + x2*(3.*z2+y2*(zr-5.)));
    val[3*i+2] *= y*z*ezor*(7.*y2 + 12.*z2 + zeta*(y2-4.*z2)*r + x2*(7.+zr));
  }
  return;
}

void get_p_5f0d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double f1 = x2+y2;
    double f2 = 3.*zeta*f1-2.*zeta*z2;

    val[3*i]   *= x*z*ezor*(-9.*y2 - 4.*z2 + zeta*(3.*y2-2.*z2)*r + 3.*x2*(zr-3.));
    val[3*i+1] *= y*z*ezor*(-9.*y2 - 4.*z2 + zeta*(3.*y2-2.*z2)*r + 3.*x2*(zr-3.));
    val[3*i+2] *= ezor*(-3.*x2*x2-3.*y2*y2 + 3.*zeta*y2*z2*r - 2.*z2*z2*(-4.+zr) + x2*(-6.*y2+3.*zeta*z2*r));
  }
  return;
}

void get_p_5fp1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double f1 = x2+y2-4.*z2;

    val[3*i]   *= ezor*(-y2*y2 + 3.*y2*z2 + 4.*z2*z2 + x2*x2*(-4.+zr) + x2*(z2*(5-4.*zr) + y2*(-5.+zr)));
    val[3*i+1] *= x*y*ezor*(-3.*y2 + 2.*z2 + zeta*(y2-4.*z2)*r + x2*(zr-3.));
    val[3*i+2] *= x*z*ezor*(7.*y2 + 12.*z2 + zeta*(y2-4.*z2)*r + x2*(7.+zr));
  }
  return;
}

void get_p_5fp2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double f1 = zeta*(y2-x2);

    val[3*i]   *=  x*z*ezor*(2.*z2 + x2*(3.-zr) + y2*(1.+zr));
    val[3*i+1] *= -y*z*ezor*(2.*z2 + y2*(3.-zr) + x2*(1.+zr));
    val[3*i+2] *=  (x-y)*(x+y)*ezor*(x2+y2 + z2*(2.-zr));
  }
  return;
}

void get_p_5fp3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double f1 = x2-3.*y2;

    val[3*i]   *= ezor*(-3.*y2*(y2+z2) + x2*x2*(4.-zr) + 3.*x2*(z2+y2*(zr-1.)));
    val[3*i+1] *= -x*y*ezor*(6.*z2+y2*(9.-3.*zr) + x2*(5.+zr));
    val[3*i+2] *= -x*z*ezor*(x2-3.*y2)*(zr-1.);
  }
  return;
}

void get_p_6fm3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double zr = zeta*r;
    double ezr = exp(-zr);

    val[3*i]   *= x*y*(6.*z2 - 3.*x2*(zr-4.) + y2*(4+zr))*ezr;
    val[3*i+1] *= (3.*x2*x2 - 3.*y2*z2 + y2*y2*(zr-5.) + 3.*x2*(z2+y*(2.-zr)))*ezr;
    val[3*i+2] *= y*z*(y2-3.*x2)*(zr-2.)*ezr;
  }
  return;
}

void get_p_6fm2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double zr = zeta*r;
    double ezr = exp(-zr);
    double f1 = 3.-zr;

    val[3*i]   *= y*z*(y2+z2+x2*f1)*ezr;
    val[3*i+1] *= x*z*(x2+z2+y2*f1)*ezr;
    val[3*i+2] *= x*y*(x2+y2+z2*f1)*ezr;
  }
  return;
}

void get_p_6fm1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double zr = zeta*r;
    double ezr = exp(-zr);
    double f1 = zeta*(y2-4.*z2)*r;

    val[3*i]   *= x*y*(-4.*y2+6.*z2+f1 + x2*(zr-4.))*ezr;
    val[3*i+1] *= (-x2*x2 + 4.*z2*z2 + y2*z2*(9.-4.*zr) + y2*y2*(zr-5.) + x2*(3.*z2+y2*(zr-6.)))*ezr;
    val[3*i+2] *= y*z*(6.*y2+16.*z2 + f1 + x2*(6.+zr))*ezr;
  }
  return;
}

void get_p_6f0d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double zr = zeta*r;
    double ezr = exp(-zr);
    double f1 = zeta*(3.*y2-2.*z2)*r;
    double f2 = -2.*(6.*y2+z2) + 3.*x2*(zr-4.);

    val[3*i]   *= x*y*(f1 + f2)*ezr;
    val[3*i+1] *= y*z*(f1 + f2)*ezr;
    val[3*i+2] *= (-3.*x2*x2 - 3.*y2*y2 - 2.*z2*z2*(zr-5.) + 3.*y2*z2*(zr-1.) - 3.*x2*(2.*y2+z2*(1.-zr)))*ezr;
  }
  return;
}

void get_p_6fp1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double zr = zeta*r;
    double ezr = exp(-zr);
    double f1 = zeta*(y2-4.*z2)*r;

    val[3*i]   *= (-y2*y2+3.*y2*z2+4.*z2*z2+ x2*x2*(zr-5.) + x2*(z2*(9.-4.*zr)+y2*(zr-6.)))*ezr;
    val[3*i+1] *= x*y*(-4.*y2 + 6.*z2 + f1 + x2*(zr-4.))*ezr;
    val[3*i+2] *= x*z*(6.*y2 + 16.*z2 + f1 + x2*(6.+zr))*ezr;
  }
  return;
}

void get_p_6fp2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double zr = zeta*r;
    double ezr = exp(-zr);

    val[3*i]   *= x*z*(2.*z2 + zr*y2 + x2*(4.-zr))*ezr;
    val[3*i+1] *= y*z*(-2.*z2 - zr*x2 + y2*(zr-4.))*ezr;
    val[3*i+2] *= (x-y)*(x+y)*(x2 + y2 + z2*(3.-zr))*ezr;
  }
  return;
}

void get_p_6fp3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double zr = zeta*r;
    double ezr = exp(-zr);

    val[3*i]   *= (-3.*y2*(y2+z2) + x2*x2*(5.-zr) + 3.*x2*(z2+y2*(zr-2.)))*ezr;
    val[3*i+1] *= x*y*(-6.*z2+3.*y2*(zr-4.) + x2*(4.+zr))*ezr;
    val[3*i+2] *= x*z*(x2-3.*y2)*(2.-zr)*ezr;
  }
  return;
}

void get_p_7fm3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r2 = r*r;
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2;
    double zr = zeta*r;
    double ezr = exp(-zr);

    val[3*i]   *= x*y*(-zeta*(3.*x2-y2)*r2+3.*r*(5.*x2+y2+2.*z2))*ezr;
    val[3*i+1] *= (zeta*y2*(-3.*x2+y2)*r2+3*r*(x4+3.*x2*y2-2.*y4+(x-y)*(x+y)*z2))*ezr;
    val[3*i+2] *= y*(-3*x2+y2)*z*(-3*r+zeta*r2)*ezr;
  }
  return;
}

void get_p_7fm2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r2 = r*r;
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double zr = zeta*r;
    double ezr = exp(-zr);

    val[3*i]   *= y*z*(-zeta*x2*r2+r*(4*x2+y2+z2))*ezr;
    val[3*i+1] *= x*z*(-zeta*y2*r2+r*(4*x2+y2+z2))*ezr;
    val[3*i+2] *= x*y*(-zeta*z2*r2+r*(4*x2+y2+z2))*ezr;
  }
  return;
}

void get_p_7fm1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r2 = r*r;
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;
    double zr = zeta*r;
    double ezr = exp(-zr);

    val[3*i]   *= x*y*(-5*(x2+y2-2*z2)*r+zeta*(x2+y2-4*z2)*r2)*ezr;
    val[3*i+1] *= (zeta*y2*(x2+y2-4*z2)*r2-r*(x4+6*y4-13*y2*z2-4*z4+x2*(7*y2-3*z2)))*ezr;
    val[3*i+2] *= y*z*(zeta*(x2+y2-4*z2)*r2+5*r*(x2+y2+4*z2))*ezr;
  }
  return;
}

void get_p_7f0d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r2 = r*r;
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double z4 = z2*z2;
    double zr = zeta*r;
    double ezr = exp(-zr);
    double x2y2 = x2+y2;

    val[3*i]   *= x*z*(-15*x2y2*r+zeta*(3*x2y2-2*z2)*r2)*ezr;
    val[3*i+1] *= y*z*(-15*x2y2*r+zeta*(3*x2y2-2*z2)*r2)*ezr;
    val[3*i+2] *= (zeta*z2*(3*x2y2-2*z2)*r2-3*r*(x2y2*x2y2+2*x2y2*z2-4*z4))*ezr;
  }
  return;
}

void get_p_7fp1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r2 = r*r;
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;
    double zr = zeta*r;
    double ezr = exp(-zr);

    val[3*i]   *= (zeta*x2*(x2+y2-4*z2)*r2-r*(6*x4+7*x2*y2+y4-(13*x2+3*y2)*z2-4*z4))*ezr;
    val[3*i+1] *= x*y*(-5*(x2+y2-2*z2)*r+zeta*(x2+y2-4*z2)*(x2+y2+z2))*ezr;
    val[3*i+2] *= x*z*(zeta*(x2+y2-4*z2)*r2+5*r*(x2+y2+4*z2))*ezr;
  }
  return;
}

void get_p_7fp2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r2 = r*r;
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double zr = zeta*r;
    double ezr = exp(-zr);

    val[3*i]   *= x*z*((-5*x2+y2-2*z2)*r+zeta*(x-y)*(x+y)*r2)*ezr;
    val[3*i+1] *= y*z*((x2-5*y2-2*z2)*r-zeta*(x-y)*(x+y)*r2)*ezr;
    val[3*i+2] *= (x-y)*(x+y)*(-zeta*z2*r2+r*(x2+y2+4*z2))*ezr;
  }
  return;
}

void get_p_7fp3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r2 = r*r;
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2;
    double zr = zeta*r;
    double ezr = exp(-zr);
    double x23y2 = x2-3.*y2;

    val[3*i]   *= (-zeta*x2*x23y2*r2+3*r*(2*x4+x2*(-3*y2+z2)-y2*(y2+z2)))*ezr;
    val[3*i+1] *= x*y*(zeta*x23y2*r2+3*r*(x2+5*y2+2*z2))*ezr;
    val[3*i+2] *= x*x23y2*z*(-3*r+zeta*r2)*ezr;
  }
  return;
}

void get_p_5gm4d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double f0 = zeta*(y2-x2);
    double xmy = 3.*x2-y2;
    double ymx = 3.*y2-x2;

    val[3*i]   *= y*(x2*f0 + xmy*r)*ezor;
    val[3*i+1] *= x*(y2*f0 - ymx*r)*ezor;
    val[3*i+2] *= x*y*z*f0*ezor;
  }
  return;
}

void get_p_5gm3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double f0 = (x-y)*(x+y);
    double xmy = 3.*x2-y2;

    val[3*i]   *= x*y*z*(-zeta*xmy + 6.*r)*ezor;
    val[3*i+1] *= z*(-zeta*y2*xmy + 3.*f0*r)*ezor;
    val[3*i+2] *= y*xmy*(-zeta*z2+r)*ezor;
  }
  return;
}

void get_p_5gm2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double f0 = 6.*z2-x2-y2;
    double f1 = f0-2.*x2;
    double f2 = f0-2.*y2;

    val[3*i]   *= y*(-zeta*x2*f0 + f1*r)*ezor;
    val[3*i+1] *= x*(-zeta*y2*f0 + f2*r)*ezor;
    val[3*i+2] *= x*y*z*(-zeta*f0 + 12.*r)*ezor;
  }
  return;
}

void get_p_5gm1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double f0 = zeta*(3.*(x2+y2)-4.*z2);

    val[3*i]   *= x*y*z*(f0 - 6.*r)*ezor;
    val[3*i+1] *= z*(y2*f0 + (4.*z2-3.*x2-9.*y2)*r)*ezor;
    val[3*i+2] *= y*(z2*f0 - 3.*(x2+y2-4.*z2)*r)*ezor;
  }
  return;
}

void get_p_5g0d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double x2py2 = x2+y2;
    double xy2 = x2py2*x2py2;
    double f0 = zeta*(-3.*xy2+24.*x2py2*z2-8.*z2*z2);

    //val[3*i]   *= x*(12.*x2py2*r - 48.*z2*r + f0)*ezor;
    //val[3*i+1] *= y*(12.*x2py2*r - 48.*z2*r + f0)*ezor;
    //val[3*i+2] *= z*(-48.*x2py2*r + 32.*z2*r + f0)*ezor;

    val[3*i]   *= x*(12.*(x2py2-4.*z2)*r + f0)*ezor;
    val[3*i+1] *= y*(12.*(x2py2-4.*z2)*r + f0)*ezor;
    val[3*i+2] *= z*(-16.*(3.*x2py2-2.*z2)*r + f0)*ezor;
  }
  return;
}

void get_p_5gp1d(int tid, int gs, double* grid, double* val, double zeta)
{
 //rederived this one. could use as 6g,7g, etc
  //int nlm = 0;
  //int pw = nlm-3;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    //double x4 = x2*x2; double y4 = y2*z2; double z4 = z2*z2;
    double zr = zeta*r;
    //double ezr = exp(-zr);
    //double rp = pow(r,pw);
    double ezor = exp(-zr)/r;
    double x2py2 = x2+y2;
    double f0 = 4.*z2-9.*x2-3.*y2;
    double f1 = 3.*x2py2 - 4.*z2;

    //val[3*i+0] *= z*(zeta*x2*(3.*x4 + 3.*y4 - y2*z2 - 4.*z4 + x2*(6.*y2 - z2)) + r*((-9. - 3.*nlm)*x4 - 3.*y4 + y2*z2 + 4.*z4 + x2*((-12. - 3.*nlm)*y2 + (-5. + 4.*nlm)*z2)))*rp*ezr;
    //val[3*i+1] *= x*y*z*(r*((-6. - 3.*nlm)*x2 + (-6. - 3.*nlm)*y2 + (-6. + 4.*nlm)*z2) + zeta*(3.*x4 + 3.*y4 - y2*z2 - 4.*z4 + x2*(6.*y2 - z2)))*rp*ezr;
    //val[3*i+2] *= x*(r*(-3.*x4 - 6.*x2*y2 - 3.*y4 + (9. - 3.*nlm)*(x2 + y2)*z2 + (12. + 4.*nlm)*z4) + zeta*z2*(3.*x4 + 3.*y4 - y2*z2 - 4.*z4 + x2*(6.*y2 - z2)))*rp*ezr;
    val[3*i]   *= z*(zeta*x2*f1+r*f0)*ezor;
    val[3*i+1] *= x*y*z*(zeta*f1 - 6.*r)*ezor;
    val[3*i+2] *= x*(zeta*f1*z2 - 3.*r*(x2py2-4.*z2))*ezor;
  }
  return;
}

void get_p_5gp2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double x2py2 = x2+y2;
    double xy2 = x2py2*x2py2;
    double f0 = x2py2-6.*z2;
    double f1 = (x-y)*(x+y);

    val[3*i]   *= x*(zeta*f1*f0 - 4.*(x2-3.*z2)*r)*ezor;
    val[3*i+1] *= y*(zeta*f1*f0 + 4.*(y2-3.*z2)*r)*ezor;
    val[3*i+2] *= z*f1*(zeta*f0 + 12.*r)*ezor;
  }
  return;
}

void get_p_5gp3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double f0 = x2-3.*y2;

    val[3*i]   *= z*(-zeta*x2*f0 + 3.*(x-y)*(x+y)*r)*ezor;
    val[3*i+1] *= -x*y*z*(zeta*f0 + 6.*r)*ezor;
    val[3*i+2] *= x*f0*(-zeta*z2 + r)*ezor;
  }
  return;
}

void get_p_5gp4d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;
    double f0 = x2*x2-6.*x2*y2+y2*y2;

    val[3*i]   *= x*(-zeta*f0 + 4.*(x2-3.*y2)*r)*ezor;
    val[3*i+1] *= -y*(zeta*f0 + 4.*(3.*x2-y2)*r)*ezor;
    val[3*i+2] *= z*zeta*f0*ezor;
  }
  return;
}

void get_p_6hm5d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;
    double zr = zeta*r;
    double ezr = exp(-zr);
    //double ezor = exp(-zr)/r;
    double f1 = 5.*x4-10.*x2*y+y4;

    val[3*i]   *= x*y*(20.*(y2-x2) - (10.*x2*y2-5.*x4-y4)*zeta/r) * ezr;
    val[3*i+1] *= (y2*(20.*x2-4.*y2) + (-5.*x4+10.*x2*y2-y4) - y2*(-5.*x4+10.*x2*y2-y4)*zeta/r) * ezr;
    val[3*i+2] *= y*z*(-5.*x4+10.*x2*y2-y4)*zeta/r * ezr;
    //val[3*i]   *= -x*y*(zeta*f1 + 20.*(y2-x2)*r)*ezor;
    //val[3*i+1] *= (-zeta*y2*f1 + 5.*(x4-6.*x2*y2+y4)*r)*ezor;
    //val[3*i+2] *= -y*z*f1*ezor;
  }
  return;
}

void get_p_6hm4d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;
    double zr = zeta*r;
    double ezr = exp(-zr);
    //double f0 = (x-y)*(x+y);

    val[3*i]   *= (-2.*x2*y*z + y*z*(y2-x2) - x2*y*z*(y2-x2)*zeta/r) * ezr;
    val[3*i+1] *= ( 2.*x2*y*z + x*z*(y2-x2) - x*y2*z*(y2-x2)*zeta/r) * ezr;
    val[3*i+2] *= (x*y*(y2-x2) - x*y*(y2-x2)*z2*zeta/r) * ezr;
    //val[3*i]   *= y*z*(zeta*x2*f0 + (y2-3.*x2)*r)*ezor;
    //val[3*i+1] *= x*z*(zeta*y2*f0 + (3.*y2-x2)*r)*ezor;
    //val[3*i+2] *= x*y*f0*(r-zeta*z2)*ezor;
  }
  return;
}

void get_p_6hm3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z; double r2 = r*r;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;
    double zr = zeta*r;
    double ezr = exp(-zr);
    double f0 = (9.*(y2-3.*x2)*z2 + (3.*x2-y2)*r2)*zeta/r;

    val[3*i]   *= y*(2.*x*(3.*x2-y2) - 54.*x*z2 + 6.*x*r2 - x*f0) * ezr;
    val[3*i+1] *= (y*(2.*y*(3.*x2-y2) + 18.*y*z2 - 2.*y*r2) + 9.*(y2-3.*x2)*z2 + (3.*x2-y2)*r2 - y2*f0) * ezr;
    val[3*i+2] *= z*(y*(2.*(3.*x2-y2) + 18.*(y2-3.*x2)) - y*f0) * ezr;
    //val[3*i]   *= (-zeta*x2*f0 + 5.*(x4-6.*x2*y2+y4)*r)*ezor;
    //val[3*i+1] *= x*y*(zeta*f0 + 20.*(x-y)*(x+y)*r)*ezor;
    //val[3*i+2] *= x*z*f0*ezor;
  }
  return;
}

void get_p_6hm2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;
    double zr = zeta*r;
    double ezr = exp(-zr);
    double f1 = 2.*z2-x2-y2;
    double f0 = x*y*z*f1*zeta/r;

    val[3*i]   *= (-2.*x2*y*z + y*z*f1 - x*f0) * ezr;
    val[3*i+1] *= (-2.*x*y2*z + x*z*f1 - y*f0) * ezr;
    val[3*i+2] *= ( 4.*x*y*z2 + x*y*f1 - z*f0) * ezr;
    //val[3*i]   *= x*y*(zeta*f2*f1 - 3.*f0*r)*ezor;
    //val[3*i+1] *= (zeta*y2*f2*f1 + r*(-3.*x4-6.*x2*y2+5.*y4+24.*(x-y)*(x+y)*z2))*ezor;
    //val[3*i+2] *= y*z*f2*(zeta*f1 + 16.*r)*ezor;
  }
  return;
}

void get_p_6hm1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z; double r2 = r*r;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2; double r4 = r2*r2;
    double zr = zeta*r;
    double ezr = exp(-zr);
    double f0 = 28.*z2-4.*r2;
    double f1 = (-21.*z4 + 14.*z2*r2 - r4)*zeta/r;

    val[3*i]   *= (x*f0 - x*f1) * ezr;
    val[3*i+1] *= (y*f0 - y*f1) * ezr;
    val[3*i+2] *= (-56.*z*z2 + 24.*z*r2 - z*f1) * ezr;
    //val[3*i]   *= x*x*(zeta*f0-4.*(x2-3.*y2)*r)*ezor;
    //val[3*i+1] *= y*z*(zeta*f0+4.*(3.*x2-y2)*r)*ezor;
    //val[3*i+2] *= f0*(-zeta*z2+r)*ezor;
  }
  return;
}

void get_p_6h0d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z; double r2 = r*r;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2; double r4 = r2*r2;
    double zr = zeta*r;
    double ezr = exp(-zr);
    double f0 = -140.*z2 + 60.*r2;
    double f1 = 112.*z2 - 80.*r2;
    double f2 = (63.*z4 - 70.*z2*r2 + 15.*r4)*zeta/r;

    val[3*i]   *= (x*f0 - x*f2)*ezr;
    val[3*i+1] *= (y*f0 - y*f2)*ezr;
    val[3*i+2] *= (z*f1 - z*f2)*ezr;
    //val[3*i]   *= y*z*(-zeta*x2*f0 + (2.*x2+f0)*r)*ezor;
    //val[3*i+1] *= x*z*(-zeta*y2*f0 + (2.*y2+f0)*r)*ezor;
    //val[3*i+2] *= x*y*(-zeta*z2*f0 + (-4.*z2+f0)*r)*ezor;
  }
  return;
}

void get_p_6hp1d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z; double r2 = r*r;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2; double r4 = r2*r2;
    double zr = zeta*r;
    double ezr = exp(-zr);
    double f0 = 28.*z2 - 4.*r2;
    double f1 = -56.*z2 + 24.*r2;
    double f2 = (-21.*z4 + 14.*z2*r2 - r4)*zeta/r;

    val[3*i]   *= (x*f0 - x*f2)*ezr;
    val[3*i+1] *= (y*f0 - y*f2)*ezr;
    val[3*i+2] *= (z*f1 - z*f2)*ezr;
    //val[3*i]   *= (-zeta*x2*f1*f0 + r*(5.*x4-3.*y2*(y2-8.*z2)-6.*x2*(y2+4.*z2)))*ezor;
    //val[3*i+1] *= -x*y*(zeta*f1*f0 + 4.*r*(x2+3.*(y2-4.*z2)))*ezor;
    //val[3*i+2] *= -x*z*f1*(zeta*f0+16.*r)*ezor;
  }
  return;
}

void get_p_6hp2d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z; double r2 = r*r;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;
    double zr = zeta*r;
    double ezr = exp(-zr);
    double f0 = (x+y)*(x-y);
    double f1 = (3.*z2 - r2);
    double f2 = z*f0*f1*zeta/r;

    val[3*i]   *= (-2.*x*f0 - (x-y)*z*f1 - (x+y)*z*f1 + x*f2)*ezr;
    val[3*i+1] *= (-2.*y*f0 - (x-y)*z*f1 + (x+y)*z*f1 + y*f2)*ezr;
    val[3*i+2] *= (f0*z2 - f0*f1 + z*f2)*ezr;
    //val[3*i]   *= -x*y*(-4.*(f0-6.*z2)*r + zeta*(f1-12.*f0*z2+8.*z4))*ezor;
    //val[3*i+1] *= (zeta*y2*(-f1+12.*f0*z2-8.*z4) + r*(x4+6.*x2*y2+5.*y4-12.*(x2+3.*y2)*z2+8.*z4))*ezor;
    //val[3*i+2] *= -y*z*((24.*f0-32.*z2)*r + zeta*(f1-12.*f0*z2+8.*z4))*ezor;
  }
  return;
}

void get_p_6hp3d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z; double r2 = r*r;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;
    double zr = zeta*r;
    double ezr = exp(-zr);
    double f0 = x2-3.*y2;
    double f1 = (-9.*f0*z2+f0*r2)*zeta/r;

    val[3*i]   *= (x*(2.*x*f0 - 18.*x*z2 + 2.*x*r2) - 9.*f0*z2 + f0*r2 - x2*f1) * ezr;
    val[3*i+1] *= (x*(2.*y*f0 + 54.*y*z2 - 6.*y*r2) - x*y*f1) * ezr;
    val[3*i+2] *= (x*(2.*f0*z - 18.*f0*z) - x*z*f1) * ezr;
    //val[3*i]   *= x*z*(-zeta*f0*f1 + 4.*(x-z)*(x+z)*r)*ezor;
    //val[3*i+1] *= y*z*(-zeta*f0*f1 + 4.*(z2-y2)*r)*ezor;
    //val[3*i+2] *= f0*(-zeta*z2*f1 + (f1-4.*z2)*r)*ezor;
  }
  return;
}

void get_p_6hp4d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;
    double zr = zeta*r;
    double ezr = exp(-zr);
    double f0 = (x4 - 6.*x2*y2 + y4)*z*zeta/r;

    val[3*i]   *= ((4.*x2 - 12.*y2)*x*z - x*f0) * ezr;
    val[3*i+1] *= ((-12.*x2 + 4.*y2)*y*z - y*f0) * ezr;
    val[3*i+2] *= (x4 - 6.*x2*y2 + y4 - z*f0) * ezr;
    //val[3*i]   *= x*z*(20.*f1*r + zeta*f2)*ezor;
    //val[3*i+1] *= y*z*(20.*f1*r + zeta*f2)*ezor;
    //val[3*i+2] *= (zeta*z2*f2 + 5.*r*(3.*f0*f0-24.*f0*z2+8.*z4))*ezor;
  }
  return;
}

void get_p_6hp5d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;
    double zr = zeta*r;
    double ezr = exp(-zr);
    double f0 = (-x4 + 10.*x2*y2 - 5.*y4);
    double f1 = x*f0*zeta/r;

    val[3*i]   *= (x2*(-4.*x2 + 20.*y2) + f0 - x*f1) * ezr;
    val[3*i+1] *= (x*y*(20.*x2 - 20.*y2) - y*f1) * ezr;
    val[3*i+2] *= -(z*f1) * ezr;
    //val[3*i]   *= (-x2*f1 + r*(5.*x4+6.*x2*y2+y4-12.*(3.*x2+y2)*z2 + 8.*z4))*ezor;
    //val[3*i+1] *= x*y*(4.*(f0-6.*z2)*r - f1)*ezor;
    //val[3*i+2] *= -x*z*((24.*f0-32.*z2)*r + f1)*ezor;
  }
  return;
}

#if 0
//incomplete
void get_p_7hm5d(int tid, int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;
    double zr = zeta*r;
    double ezor = exp(-zr)/r;

    val[3*i]   *= -x*y*(20.*y2*z2+5.*x2*x2*(zr-5.) + y4*(19.+zr) - 10.*x2*(2.*z2+y2*(zr-1.)))*ezor;
    val[3*i+1] *= (5.*x2*x4+ 5.*y4*z2 + y4*y2*(6.-zr) + 5.*x4*(z2-y2*(4.+zr)) + 5.*x2*(-6.y2*z2 + y4*(2.*zr-6.*y2*z2)))*ezor;
    val[3*i+2] *= -y*z*(5.*x4-10.*x2*y2+y4)*(zr-1.))*ezor;
  }
  return;
}
#endif

void get_dp_pxd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double xor3 = x/r/r/r;

    val[3*i+0] *= 1./r - x*xor3;
    val[3*i+1] *= -y*xor3;
    val[3*i+2] *= -z*xor3;
  }
  return;
}

void get_dp_pyd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double yor3 = y/r/r/r;

    val[3*i+0] *= -x*yor3;
    val[3*i+1] *= 1./r - y*yor3;
    val[3*i+2] *= -z*yor3;
  }
  return;
}

void get_dp_pzd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double zor3 = z/r/r/r;

    val[3*i+0] *= -x*zor3;
    val[3*i+1] *= -y*zor3;
    val[3*i+2] *= 1./r - z*zor3;
  }
  return;
}

void get_dp_dxyd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r2 = r*r; double txyrm4 = 2.*x*y/r2/r2;

    val[3*i+0] *= -x*txyrm4 + y/r2;
    val[3*i+1] *= -y*txyrm4 + x/r2;
    val[3*i+2] *= -z*txyrm4;
  }
  return;
}

void get_dp_dyzd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r2 = r*r; double tyzrm4 = 2.*y*z/r2/r2;

    val[3*i+0] *= -x*tyzrm4;
    val[3*i+1] *= -y*tyzrm4 + z/r2;
    val[3*i+2] *= -z*tyzrm4 + y/r2;
  }
  return;
}

void get_dp_dz2d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r2 = r*r; double szrm4 = 6.*z/r2/r2;

    val[3*i+0] *= -x*z*szrm4;
    val[3*i+1] *= -y*z*szrm4;
    val[3*i+2] *= szrm4*(x*x+y*y);
  }
  return;
}

void get_dp_dxzd(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r2 = r*r; double txzrm4 = 2.*x*z/r2/r2;

    val[3*i+0] *= -x*txzrm4 + z/r2;
    val[3*i+1] *= -y*txzrm4;
    val[3*i+2] *= -z*txzrm4 + x/r2;
  }
  return;
}

void get_dp_dx2y2d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double r2 = r*r; double tx2y2rm4 = 2.*(x*x-y*y)/r2/r2;

    val[3*i+0] *= -x*tx2y2rm4 + 2.*x/r2;
    val[3*i+1] *= -y*tx2y2rm4 - 2.*y/r2;
    val[3*i+2] *= -z*tx2y2rm4;
  }
  return;
}

void get_dp_fm3d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y;
    double rm3 = 1./r/r/r;
    double rm5 = rm3/r/r;
    double f0 = 3.*x2-y2;
    double f1 = 3.*y*f0*rm5;

    val[3*i+0] *= -x*f1 + 6.*x*y*rm3;
    val[3*i+1] *= -y*f1 - 2.*y2*rm3 + f0*rm3;
    val[3*i+2] *= -z*f1;
  }
  return;
}

void get_dp_fm2d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y;
    double rm3 = 1./r/r/r;
    double rm5 = rm3/r/r;
    double f1 = 3.*x*y*z*rm5;

    val[3*i+0] *= -x*f1 + y*z*rm3;
    val[3*i+1] *= -y*f1 + x*z*rm3;
    val[3*i+2] *= -z*f1 + x*y*rm3;
  }
  return;
}

void get_dp_fm1d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y;
    double rm3 = 1./r/r/r;
    double rm5 = rm3/r/r;
    double f0 = 4.*z*z-x*x-y*y;
    double f1 = 3.*y*f0*rm5;

    val[3*i+0] *= -x*f1 - 2.*x*y*rm3;
    val[3*i+1] *= -y*f1 - 2.*y*y*rm3 + f0*rm3;
    val[3*i+2] *= -z*f1 + 8.*y*z*rm3;
  }
  return;
}

void get_dp_f0d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y;
    double rm3 = 1./r/r/r;
    double rm5 = rm3/r/r;
    double f0 = 2.*z*z-3.*x*x-3.*y*y;
    double f1 = 3.*z*f0*rm5;

    val[3*i+0] *= -x*f1 - 6.*x*z*rm3;
    val[3*i+1] *= -y*f1 - 6.*y*z*rm3;
    val[3*i+2] *= -z*f1 + 4.*z*z*rm3 + f0*rm3;
  }
  return;
}

void get_dp_fp1d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y;
    double rm3 = 1./r/r/r;
    double rm5 = rm3/r/r;
    double f0 = 4.*z*z-x*x-y*y;
    double f1 = 3.*x*f0*rm5;

    val[3*i+0] *= -x*f1 - 2.*x*x*rm3 + f0*rm3;
    val[3*i+1] *= -y*f1 - 2.*x*y*rm3;
    val[3*i+2] *= -z*f1 + 8.*x*z*rm3;
  }
  return;
}

void get_dp_fp2d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y;
    double rm3 = 1./r/r/r;
    double rm5 = rm3/r/r;
    double f0 = x*x-y*y;
    double f1 = 3.*z*f0*rm5;

    val[3*i+0] *= -x*f1 + 2.*x*z*rm3;
    val[3*i+1] *= -y*f1 - 2.*y*z*rm3;
    val[3*i+2] *= -z*f1 + f0*rm3;
  }
  return;
}

void get_dp_fp3d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y;
    double rm3 = 1./r/r/r;
    double rm5 = rm3/r/r;
    double f0 = x*x-3.*y*y;
    double f1 = 3.*x*f0*rm5;

    val[3*i+0] *= -x*f1 + f0*rm3 + 2.*x*x*rm3;
    val[3*i+1] *= -y*f1 - 6.*x*y*rm3;
    val[3*i+2] *= -z*f1;
  }
  return;
}

void get_dp_gm4d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double rm2 = 1./r/r;
    double f0 = (x2-y2)*rm2;
    double xyr = x*y*rm2;
    double f1 = 4.*xyr*f0*rm2;

    val[3*i+0] *= -x*f1 + 2.*x*xyr*rm2 + y*f0*rm2;
    val[3*i+1] *= -y*f1 - 2.*y*xyr*rm2 + x*f0*rm2;
    val[3*i+2] *= -z*f1;
  }
  return;
}

void get_dp_gm3d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double rm2 = 1./r/r;
    double f0 = (3.*x2-y2)*rm2;
    double yzr = y*z*rm2;
    double f1 = yzr*f0*rm2;

    val[3*i+0] *= -4.*x*f1 + 6.*x*yzr*rm2;
    val[3*i+1] *= -4.*y*f1 - 2.*y*yzr*rm2 + z*f0*rm2;
    val[3*i+2] *= -4.*z*f1 + y*f0*rm2;
  }
  return;
}

void get_dp_gm2d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double rm2 = 1./r/r;
    double f0 = (6.*z2-x2-y2)*rm2;
    double xyr = x*y*rm2;
    double f1 = xyr*f0*rm2;

    val[3*i+0] *= -4.*x*f1 + 2.*x*xyr*rm2 + y*f0*rm2;
    val[3*i+1] *= -4.*y*f1 - 2.*y*xyr*rm2 + x*f0*rm2;
    val[3*i+2] *= -4.*z*f1 + 12.*xyr*z*rm2;
  }
  return;
}

void get_dp_gm1d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double r2 = r*r;
    double rm2 = 1./r/r;
    double x2r = x2*rm2; double y2r = y2*rm2; double z2r = z2*rm2;
    double x2y2r = (x2+y2)*rm2;
    //double f0 = 4.*z2-3.*x2-3.*y2;
    double f0 = 7.*z2r-3.;
    double yzr = y*z*rm2;
    double f1 = 4.*yzr*f0*rm2;

    //val[3*i+0] *= -x*f1 - 6.*x*yz*rm4;
    //val[3*i+1] *= -y*f1 - 6.*y*yz*rm4 + z*f0*rm4;
    //val[3*i+2] *= -z*f1 + 8.*z*yz*rm4 + y*f0*rm4;

    val[3*i+0] *= 2.*x*yzr*(3.*x2r+3.*y2r-11.*z2r)*rm2;
    val[3*i+1] *= z*(-3.*(-x2r*x2r+y2r*y2r) + z2r*(x2r-21.*y2r+4.*z2r))*rm2;
    val[3*i+2] *= y*(-3.*x2y2r*x2y2r + 21.*x2y2r*z2r - 4.*z2r*z2r)*rm2;
  }
  return;
}

void get_dp_g0d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double rm2 = 1./r/r;
    double z2r = z2*rm2;
    double x2py2 = x2+y2;
    double x2py2r = x2py2*rm2;
   // double f0 = 3.*x2py2*x2py2 - 24.*x2py2*z2 + 8.*z2*z2;
   // double f1 = 4.*f0*rm6;

   // val[3*i+0] *= -x*f1 + x*(12.*x2py2 - 48.*z2)*rm4;
   // val[3*i+1] *= -y*f1 + y*(12.*x2py2 - 48.*z2)*rm4;
   // val[3*i+2] *= -z*f1 + (-48.*z*x2py2 + 32.*z*z2)*rm4;

    double f2 = (3.*x2py2-4.*z2)*rm2;
    double f3 = 20.*z2r*f2*rm2;

    val[3*i+0] *= x*f3;
    val[3*i+1] *= y*f3;
    val[3*i+2] *= -20.*z*x2py2r*f2*rm2;
  }
  return;
}

void get_dp_gp1d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double rm2 = 1./r/r;
    double xzr = x*z*rm2;
    double f0 = (4.*z2-3.*x2-3.*y2)*rm2;
    double f1 = 4.*xzr*f0*rm2;

    val[3*i+0] *= -x*f1 - 6.*x*xzr*rm2 + z*f0*rm2;
    val[3*i+1] *= -y*f1 - 6.*y*xzr*rm2;
    val[3*i+2] *= -z*f1 + 8.*z*xzr*rm2 + x*f0*rm2;
  }
  return;
}

void get_dp_gp2d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double rm2 = 1./r/r;
    double x2my2r = (x2-y2)*rm2;
    double sz2x2y2r = (6.*z2-x2-y2)*rm2;
    double f0 = x2my2r*sz2x2y2r;
    double f1 = 4.*f0*rm2;

    val[3*i+0] *= -x*f1 - 2.*x*x2my2r*rm2 + 2.*x*sz2x2y2r*rm2;
    val[3*i+1] *= -y*f1 - 2.*y*x2my2r*rm2 - 2.*y*sz2x2y2r*rm2;
    val[3*i+2] *= -z*f1 + 12.*z*x2my2r*rm2;
  }
  return;
}

void get_dp_gp3d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double rm2 = 1./r/r;
    double x2r = x2*rm2;
    double xzr = x*z*rm2;
    double x2my2r = (x2-3.*y2)*rm2;
    double f1 = 4.*xzr*x2my2r*rm2;

    val[3*i+0] *= -x*f1 + 2.*x2r*z*rm2 + z*x2my2r*rm2;
    val[3*i+1] *= -y*f1 - 6.*y*xzr*rm2;
    val[3*i+2] *= -z*f1 + x*x2my2r*rm2;
  }
  return;
}

void get_dp_gp4d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double rm2 = 1./r/r;
    double x2m3y2r = (x2-3.*y2)*rm2;
    double y2m3x2r = (y2-3.*x2)*rm2;
    double f0 = x2*x2m3y2r + y2*y2m3x2r;
    double f1 = 4.*f0*rm2*rm2;

  //double check this
    val[3*i+0] *= -x*f1 + 4.*x*x2m3y2r*rm2;
    val[3*i+1] *= -y*f1 + 4.*y*y2m3x2r*rm2;
    val[3*i+2] *= -z*f1;
  }
  return;
}

void get_dp_hm5d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;

    val[3*i+0] *= 5.*x*y*(x4+5.*y4+4.*y2*z2-2.*x2*(5.*y2+2.*z2))*rm7;
    val[3*i+1] *= -5.*(x4*x2-10.*x4*y2+5.*x2*y4+(x4-6.*x2*y2+y4)*z2)*rm7;
    val[3*i+2] *= 5.*y*(5.*x4-10.*x2*y2+y4)*z*rm7;
  }
  return;
}

void get_dp_hm4d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;

    val[3*i+0] *= y*z*(2.*x4+y2*(y2+z2)-x2*(7.*y2+3.*z2))*rm7;
    val[3*i+1] *= -x*z*(x4+2.*y4-3.*y2*z2+x2*(-7.*y2+z2))*rm7;
    val[3*i+2] *= -x*(x-y)*y*(x+y)*(x2+y2-4.*z2)*rm7;
  }
  return;
}

void get_dp_hm3d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;

    val[3*i+0] *= -3.*x*y*(x4-3.*y4+28.*y2*z2+16.*z4-2*x2*(y2+14.*z2))*rm7;
    val[3*i+1] *= 3.*(x2*(x2-3.*y2)*(x2+y2)-7.*(x4-6.*x2*y2+y4)*z2+8.*(-x2+y2)*z4)*rm7;
    val[3*i+2] *= 3.*y*(-3.*x2+y2)*z*(7.*(x2+y2)-8.*z2)*rm7;
  }
  return;
}

void get_dp_hm2d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;

    val[3*i+0] *= y*z*(2*x4-y4+y2*z2+2.*z4+x2*(y2-11.*z2))*rm7;
    val[3*i+1] *= x*z*(-x4+2.*y4-11.*y2*z2+2.*z4+x2*(y2+z2))*rm7;
    val[3*i+2] *= x*y*(-(x2+y2)*(x2+y2)+10.*(x2+y2)*z2-4.*z4)*rm7;
  }
  return;
}

void get_dp_hm1d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;

    val[3*i+0] *= x*y*((x2+y2)*(x2+y2)-40*(x2+y2)*z2+64.*z4)*rm7;
    val[3*i+1] *= -x2*(x2+y2)*(x2+y2)+(11.*x2-29.*y2)*(x2+y2)*z2+4.*(x2+17.*y2)*z4-8.*z4*z2*rm7;
    val[3*i+2] *= y*z*(29.*(x2+y2)*(x2+y2)-68.*(x2+y2)*z2+8.*z4)*rm7;
  }
  return;
}

void get_dp_h0d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;

    val[3*i+0] *= 15.*x*z*(-(x2+y2)*(x2+y2)+12.*(x2+y2)*z2-8.*z4)*rm7;
    val[3*i+1] *= 15.*y*z*(-(x2+y2)*(x2+y2)+12.*(x2+y2)*z2-8.*z4)*rm7;
    val[3*i+2] *= 15.*(x2+y2)*((x2+y2)*(x2+y2)-12.*(x2+y2)*z2+8.*z4)*rm7;
  }
  return;
}

void get_dp_hp1d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;

    val[3*i+0] *= -y2*(x2+y2)*(x2+y2)+(x2+y2)*(-29.*x2+11.*y2)*z2+4.*(17.*x2+y2)*z4-8.*z4*z2*rm7;
    val[3*i+1] *= x*y*((x2+y2)*(x2+y2)-40*(x2+y2)*z2+64.*z4)*rm7;
    val[3*i+2] *= x*z*(29.*(x2+y2)*(x2+y2)-68.*(x2+y2)*z2+8.*z4)*rm7;
  }
  return;
}

void get_dp_hp2d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;

    val[3*i+0] *= x*z*(x4-5.*y4+14.*y2*z2+4.*z4-2*x2*(2*y2+5.*z2))*rm7;
    val[3*i+1] *= -y*z*(-5.*x4-4.*x2*y2+y4+2.*(7.*x2-5.*y2)*z2+4.*z4)*rm7;
    val[3*i+2] *= -((x-y)*(x+y)*((x2+y2)*(x2+y2)-10.*(x2+y2)*z2+4.*z4))*rm7;
  }
  return;
}

void get_dp_hp3d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;

    val[3*i+0] *= 9.*x4*y2+6.*x2*y4-3.*y4*y2+21*(x4-6.*x2*y2+y4)*z2+24.*(-x2+y2)*z4*rm7;
    val[3*i+1] *= 3.*x*y*(-3.*x4-2*x2*y2+y4+28.*(x-y)*(x+y)*z2+16.*z4)*rm7;
    val[3*i+2] *= -3.*x*(x2-3.*y2)*z*(7.*(x2+y2)-8.*z2)*rm7;
  }
  return;
}

void get_dp_hp4d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;

    val[3*i+0] *= -x*z*(x4-22.*x2*y2+17.*y4-4.*(x2-3.*y2)*z2)*rm7;
    val[3*i+1] *= -y*z*(17.*x4+y4-4.*y2*z2+x2*(-22.*y2+12.*z2))*rm7;
    val[3*i+2] *= (x4-6.*x2*y2+y4)*(x2+y2-4.*z2)*rm7;
  }
  return;
}

void get_dp_hp5d(int tid, int gs, double* grid, double* val)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs]) async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i];
    double y = grid[6*i+1];
    double z = grid[6*i+2];
    double r = grid[6*i+3];
    double rm7 = pow(r,-7.);
    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double x4 = x2*x2; double y4 = y2*y2; double z4 = z2*z2;

    val[3*i+0] *= -5.*(5.*x4*y2-10.*x2*y4+y4*y2+(x4-6.*x2*y2+y4)*z2)*rm7;
    val[3*i+1] *= 5.*x*y*(5.*x4-10.*x2*y2+y4+4.*(x-y)*(x+y)*z2)*rm7;
    val[3*i+2] *= 5.*x*(x4-10.*x2*y2+5.*y4)*z*rm7;
  }
  return;
}



//derivatives of the spherical part only
void eval_dp_3rd(int tid, int gs, double* grid, double* val, int n1, int l1, int m1)
{
  if (l1==0)
  {
    acc_assign(3*gs,val,0.);
  }
  else if (l1==1)
  {
    if (m1==1)
      return get_dp_pxd(tid,gs,grid,val);
    else if (m1==-1)
      return get_dp_pyd(tid,gs,grid,val);
    else
      return get_dp_pzd(tid,gs,grid,val);
  }
  else if (l1==2)
  {
    if (m1==-2)
      return get_dp_dxyd(tid,gs,grid,val);
    else if (m1==-1)
      return get_dp_dyzd(tid,gs,grid,val);
    else if (m1== 0)
      return get_dp_dz2d(tid,gs,grid,val);
    else if (m1== 1)
      return get_dp_dxzd(tid,gs,grid,val);
    else if (m1== 2)
      return get_dp_dx2y2d(tid,gs,grid,val);
  }
  else if (l1==3)
  {
    if (m1==-3)
      return get_dp_fm3d(tid,gs,grid,val);
    else if (m1==-2)
      return get_dp_fm2d(tid,gs,grid,val);
    else if (m1==-1)
      return get_dp_fm1d(tid,gs,grid,val);
    else if (m1== 0)
      return get_dp_f0d(tid,gs,grid,val);
    else if (m1== 1)
      return get_dp_fp1d(tid,gs,grid,val);
    else if (m1== 2)
      return get_dp_fp2d(tid,gs,grid,val);
    else if (m1== 3)
      return get_dp_fp3d(tid,gs,grid,val);
  }
  else if (l1==4)
  {
    if (m1==-4)
      return get_dp_gm4d(tid,gs,grid,val);
    else if (m1==-3)
      return get_dp_gm3d(tid,gs,grid,val);
    else if (m1==-2)
      return get_dp_gm2d(tid,gs,grid,val);
    else if (m1==-1)
      return get_dp_gm1d(tid,gs,grid,val);
    else if (m1== 0)
      return get_dp_g0d(tid,gs,grid,val);
    else if (m1== 1)
      return get_dp_gp1d(tid,gs,grid,val);
    else if (m1== 2)
      return get_dp_gp2d(tid,gs,grid,val);
    else if (m1== 3)
      return get_dp_gp3d(tid,gs,grid,val);
    else if (m1== 4)
      return get_dp_gp4d(tid,gs,grid,val);
  }
  else if (l1==5)
  {
    if (m1==-5)
      return get_dp_hm5d(tid,gs,grid,val);
    else if (m1==-4)
      return get_dp_hm4d(tid,gs,grid,val);
    else if (m1==-3)
      return get_dp_hm3d(tid,gs,grid,val);
    else if (m1==-2)
      return get_dp_hm2d(tid,gs,grid,val);
    else if (m1==-1)
      return get_dp_hm1d(tid,gs,grid,val);
    else if (m1== 0)
      return get_dp_h0d(tid,gs,grid,val);
    else if (m1== 1)
      return get_dp_hp1d(tid,gs,grid,val);
    else if (m1== 2)
      return get_dp_hp2d(tid,gs,grid,val);
    else if (m1== 3)
      return get_dp_hp3d(tid,gs,grid,val);
    else if (m1== 4)
      return get_dp_hp4d(tid,gs,grid,val);
    else if (m1== 5)
      return get_dp_hp5d(tid,gs,grid,val);
  }
  return;
}

//derivatives including the exp term
void eval_pd(int tid, int gs, double* grid, double* val, int n1, int l1, int m1, double zeta1)
{
  if (n1==1)
    return get_p_1sd(tid,gs,grid,val,zeta1);
  else if (n1==2)
  {
    if (l1==0)
      return get_p_2sd(tid,gs,grid,val,zeta1);
    else
    {
      if (m1==1)
        return get_p_2pxd(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_2pyd(tid,gs,grid,val,zeta1);
      else
        return get_p_2pzd(tid,gs,grid,val,zeta1);
    }
  }
  else if (n1==3)
  {
    if (l1==0)
    {
      return get_p_3sd(tid,gs,grid,val,zeta1);
    }
    else if (l1==1)
    {
      if (m1==1)
        return get_p_3pxd(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_3pyd(tid,gs,grid,val,zeta1);
      else
        return get_p_3pzd(tid,gs,grid,val,zeta1);
    }
    else if (l1==2)
    {
      if (m1==-2)
        return get_p_3dxyd(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_3dyzd(tid,gs,grid,val,zeta1);
      else if (m1== 0)
        return get_p_3dz2d(tid,gs,grid,val,zeta1);
      else if (m1== 1)
        return get_p_3dxzd(tid,gs,grid,val,zeta1);
      else if (m1== 2)
        return get_p_3dx2y2d(tid,gs,grid,val,zeta1);
    }
  }
  else if (n1==4)
  {
    if (l1==0)
    {
      return get_p_4sd(tid,gs,grid,val,zeta1);
    }
    else if (l1==1)
    {
      if (m1==1)
        return get_p_4pxd(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_4pyd(tid,gs,grid,val,zeta1);
      else
        return get_p_4pzd(tid,gs,grid,val,zeta1);
    }
    else if (l1==2)
    {
      if (m1==-2)
        return get_p_4dxyd(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_4dyzd(tid,gs,grid,val,zeta1);
      else if (m1== 0)
        return get_p_4dz2d(tid,gs,grid,val,zeta1);
      else if (m1== 1)
        return get_p_4dxzd(tid,gs,grid,val,zeta1);
      else if (m1== 2)
        return get_p_4dx2y2d(tid,gs,grid,val,zeta1);
    }
    else if (l1==3)
    {
      if (m1==-3)
        return get_p_4fm3d(tid,gs,grid,val,zeta1);
      else if (m1==-2)
        return get_p_4fm2d(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_4fm1d(tid,gs,grid,val,zeta1);
      else if (m1== 0)
        return get_p_4f0d(tid,gs,grid,val,zeta1);
      else if (m1== 1)
        return get_p_4fp1d(tid,gs,grid,val,zeta1);
      else if (m1== 2)
        return get_p_4fp2d(tid,gs,grid,val,zeta1);
      else if (m1== 3)
        return get_p_4fp3d(tid,gs,grid,val,zeta1);
    }
  }
  else if (n1==5)
  {
    if (l1==0)
    {
      return get_p_5sd(tid,gs,grid,val,zeta1);
    }
    else if (l1==1)
    {
      if (m1==1)
        return get_p_5pxd(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_5pyd(tid,gs,grid,val,zeta1);
      else
        return get_p_5pzd(tid,gs,grid,val,zeta1);
    }
    else if (l1==2)
    {
      if (m1==-2)
        return get_p_5dxyd(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_5dyzd(tid,gs,grid,val,zeta1);
      else if (m1==0)
        return get_p_5dz2d(tid,gs,grid,val,zeta1);
      else if (m1==1)
        return get_p_5dxzd(tid,gs,grid,val,zeta1);
      else if (m1==2)
        return get_p_5dx2y2d(tid,gs,grid,val,zeta1);
    }
    else if (l1==3)
    {
      if (m1==-3)
        return get_p_5fm3d(tid,gs,grid,val,zeta1);
      else if (m1==-2)
        return get_p_5fm2d(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_5fm1d(tid,gs,grid,val,zeta1);
      else if (m1==0)
        return get_p_5f0d(tid,gs,grid,val,zeta1);
      else if (m1==1)
        return get_p_5fp1d(tid,gs,grid,val,zeta1);
      else if (m1==2)
        return get_p_5fp2d(tid,gs,grid,val,zeta1);
      else if (m1==3)
        return get_p_5fp3d(tid,gs,grid,val,zeta1);
    }
    else if (l1==4)
    {
      if (m1==-4)
        return get_p_5gm4d(tid,gs,grid,val,zeta1);
      else if (m1==-3)
        return get_p_5gm3d(tid,gs,grid,val,zeta1);
      else if (m1==-2)
        return get_p_5gm2d(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_5gm1d(tid,gs,grid,val,zeta1);
      else if (m1== 0)
        return get_p_5g0d(tid,gs,grid,val,zeta1);
      else if (m1== 1)
        return get_p_5gp1d(tid,gs,grid,val,zeta1);
      else if (m1== 2)
        return get_p_5gp2d(tid,gs,grid,val,zeta1);
      else if (m1== 3)
        return get_p_5gp3d(tid,gs,grid,val,zeta1);
      else if (m1== 4)
        return get_p_5gp4d(tid,gs,grid,val,zeta1);
    }
  }
  else if (n1==6)
  {
    if (l1==0)
    {
      return get_p_6sd(tid,gs,grid,val,zeta1);
    }
    else if (l1==1)
    {
      if (m1==1)
        return get_p_6pxd(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_6pyd(tid,gs,grid,val,zeta1);
      else
        return get_p_6pzd(tid,gs,grid,val,zeta1);
    }
    else if (l1==2)
    {
      if (m1==-2)
        return get_p_6dxyd(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_6dyzd(tid,gs,grid,val,zeta1);
      else if (m1==0)
        return get_p_6dz2d(tid,gs,grid,val,zeta1);
      else if (m1==1)
        return get_p_6dxzd(tid,gs,grid,val,zeta1);
      else if (m1==2)
        return get_p_6dx2y2d(tid,gs,grid,val,zeta1);
    }
    else if (l1==3)
    {
      if (m1==-3)
        return get_p_6fm3d(tid,gs,grid,val,zeta1);
      else if (m1==-2)
        return get_p_6fm2d(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_6fm1d(tid,gs,grid,val,zeta1);
      else if (m1==0)
        return get_p_6f0d(tid,gs,grid,val,zeta1);
      else if (m1==1)
        return get_p_6fp1d(tid,gs,grid,val,zeta1);
      else if (m1==2)
        return get_p_6fp2d(tid,gs,grid,val,zeta1);
      else if (m1==3)
        return get_p_6fp3d(tid,gs,grid,val,zeta1);
    }
    else if (l1==4)
    {
      printf(" WARNING: 6g derivatives not available \n");
    }
    else if (l1==5)
    {
      if (m1==-5)
        return get_p_6hm5d(tid,gs,grid,val,zeta1);
      else if (m1==-4)
        return get_p_6hm4d(tid,gs,grid,val,zeta1);
      else if (m1==-3)
        return get_p_6hm3d(tid,gs,grid,val,zeta1);
      else if (m1==-2)
        return get_p_6hm2d(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_6hm1d(tid,gs,grid,val,zeta1);
      else if (m1== 0)
        return get_p_6h0d(tid,gs,grid,val,zeta1);
      else if (m1== 1)
        return get_p_6hp1d(tid,gs,grid,val,zeta1);
      else if (m1== 2)
        return get_p_6hp2d(tid,gs,grid,val,zeta1);
      else if (m1== 3)
        return get_p_6hp3d(tid,gs,grid,val,zeta1);
      else if (m1== 4)
        return get_p_6hp4d(tid,gs,grid,val,zeta1);
      else if (m1== 5)
        return get_p_6hp5d(tid,gs,grid,val,zeta1);
    }
  }
  else if (n1==7)
  {
    if (l1==0)
    {
      return get_p_7sd(tid,gs,grid,val,zeta1);
    }
    else if (l1==1)
    {
      if (m1==1)
        return get_p_6pxd(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_6pyd(tid,gs,grid,val,zeta1);
      else
        return get_p_6pzd(tid,gs,grid,val,zeta1);
    }
    else if (l1==2)
    {
      if (m1==-2)
        return get_p_7dxyd(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_7dyzd(tid,gs,grid,val,zeta1);
      else if (m1==0)
        return get_p_7dz2d(tid,gs,grid,val,zeta1);
      else if (m1==1)
        return get_p_7dxzd(tid,gs,grid,val,zeta1);
      else if (m1==2)
        return get_p_7dx2y2d(tid,gs,grid,val,zeta1);
    }
    else if (l1==3)
    {
      if (m1==-3)
        return get_p_7fm3d(tid,gs,grid,val,zeta1);
      else if (m1==-2)
        return get_p_7fm2d(tid,gs,grid,val,zeta1);
      else if (m1==-1)
        return get_p_7fm1d(tid,gs,grid,val,zeta1);
      else if (m1==0)
        return get_p_7f0d(tid,gs,grid,val,zeta1);
      else if (m1==1)
        return get_p_7fp1d(tid,gs,grid,val,zeta1);
      else if (m1==2)
        return get_p_7fp2d(tid,gs,grid,val,zeta1);
      else if (m1==3)
        return get_p_7fp3d(tid,gs,grid,val,zeta1);
    }
  }
  else if (n1==8)
  {
    if (l1==0)
    {
      return get_p_8sd(tid,gs,grid,val,zeta1);
    }
    else
    {
      printf(" ERROR: no n>7 l>0 derivatives available \n");
    }
  }

  if (tid<0)
  {
    #pragma acc wait
  }

  return;
}

