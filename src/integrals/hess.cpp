#include "pVp.h"

//this code is now double precision

void get_h_1s(int gs, double* grid, double* val, double zeta)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3]; double r3 = r*r*r;
    double ezr = exp(-zeta*r);
    double ezor = ezr/r3;

    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double xy = x*y; double xz = x*z; double yz = y*z;
    double zr = zeta*r; double ozr = 1.f+zr;

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

void get_h_2s(int gs, double* grid, double* val, double zeta)
{
  double zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3]; double r2 = r*r; double r3 = r2*r;
    double ezr = exp(-zeta*r);
    double ezor = ezr/r3;

    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double xy = x*y; double xz = x*z; double yz = y*z;
    double zr = zeta*r; double ozr = 1.f+zr;

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

void get_h_2px(int gs, double* grid, double* val, double zeta)
{
  double zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3]; double r3 = r*r*r;
    double ezr = exp(-zeta*r);
    double ezor = ezr/r3;

    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double xy = x*y; double xz = x*z; double yz = y*z;
    double zr = zeta*r; double ozr = 1.f+zr;

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

void get_h_2py(int gs, double* grid, double* val, double zeta)
{
  double zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3]; double r3 = r*r*r;
    double ezr = exp(-zeta*r);
    double ezor = ezr/r3;

    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double xy = x*y; double xz = x*z; double yz = y*z;
    double zr = zeta*r; double ozr = 1.f+zr;

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

void get_h_2pz(int gs, double* grid, double* val, double zeta)
{
  double zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3]; double r3 = r*r*r;
    double ezr = exp(-zeta*r);
    double ezor = ezr/r3;

    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double xy = x*y; double xz = x*z; double yz = y*z;
    double zr = zeta*r; double ozr = 1.f+zr;

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

void get_h_3dxy(int gs, double* grid, double* val, double zeta)
{
  double zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3]; double r3 = r*r*r;
    double ezr = exp(-zeta*r);
    double ezor = ezr/r3;

    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double xy = x*y; double xz = x*z; double yz = y*z;
    double zr = zeta*r;

   //xx, xy, xz, yy, yz, zz
    val[6*i]  *= ezor*x*y*zeta*(-3*(y2+z2)+x2*(-2+zr));
    val[6*i+1] *= ezor*(-x2*x2*zeta+(y2+z2)*(r-y2*zeta)+ x2*(r-z2*zeta+y2*zeta*(-1+zr)));
    val[6*i+2] *= -ezor*y*z*zeta*(y2+z2-r*x2*zeta);
    val[6*i+3] *= ezor*x*y*zeta*(-3*x2-3*z2+y2*(-2+zr));
    val[6*i+4] *= -ezor*x*z*zeta*(x2+z2-r*y2*zeta);
    val[6*i+5] *= -ezor*x*y*zeta*(x2+y2-r*z2*zeta);
  }
  return;
}

void get_h_3dyz(int gs, double* grid, double* val, double zeta)
{
  double zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3]; double r3 = r*r*r;
    double ezr = exp(-zeta*r);
    double ezor = ezr/r3;

    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double xy = x*y; double xz = x*z; double yz = y*z;
    double zr = zeta*r;

   //xx, xy, xz, yy, yz, zz
    val[6*i]   *= -ezor*y*z*zeta*(y2+z2-r*x2*zeta);
    val[6*i+1] *= -ezor*x*z*zeta*(x2+z2-r*y2*zeta);
    val[6*i+2] *= -ezor*x*y*zeta*(x2+y2-r*z2*zeta);
    val[6*i+3] *= ezor*y*z*zeta*(-3*x2-3*z2+y2*(-2+zr));
    val[6*i+4] *= ezor*(-y2*y2*zeta+z2*(r-z2*zeta)+x2*(r-(y2+z2)*zeta)+y2*(r+z2*zeta*(-1+zr)));
    val[6*i+5] *= ezor*y*z*zeta*(-3*x2-3*y2+z2*(-2+zr));
  }
  return;
}

void get_h_3dz2(int gs, double* grid, double* val, double zeta)
{
  double zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3]; double r3 = r*r*r;
    double ezr = exp(-zeta*r);
    double ezor = ezr/r3;

    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double xy = x*y; double xz = x*z; double yz = y*z;
    double zr = zeta*r;

   //xx, xy, xz, yy, yz, zz
    val[6*i]   *= ezor*(x2*x2*zeta*(4-zr)+x2*(-2*r+y2*zeta*(5-zr) + z2*zeta*(5+2*zr))+(y2+z2)*(y2*zeta-2*(r+z2*zeta)));
    val[6*i+1] *= ezor*x*y*zeta*(3*(y2+2*z2)-r*(y2-2*z2)*zeta+x2*(3-zr));
    val[6*i+2] *= -ezor*x*z*zeta*(3*y2+r*(y2-2*z2)*zeta+x2*(3+zr));
    val[6*i+3] *= ezor*(x2*x2*zeta+y2*y2*zeta*(4-zr)-2*z2*(r+z2*zeta) - x2*(2*r+z2*zeta+y2*zeta*(-5+zr))+y2*(-2*r+z2*zeta*(5+2*zr)));
    val[6*i+4] *= -ezor*y*z*zeta*(3*y2+r*(y2-2*z2)*zeta+x2*(3+zr));
    val[6*i+5] *= ezor*(4*r*z2+x2*x2*zeta+y2*y2*zeta-8*z2*z2*zeta+2*r*z2*z2*zt2 + x2*(4*r+2*y2*zeta-z2*zeta*(10+zr))-y2*(-4*r+z2*zeta*(10+zr)));
  }
  return;
}

void get_h_3dxz(int gs, double* grid, double* val, double zeta)
{
  double zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3]; double r3 = r*r*r;
    double ezr = exp(-zeta*r);
    double ezor = ezr/r3;

    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double xy = x*y; double xz = x*z; double yz = y*z;
    double zr = zeta*r;

   //xx, xy, xz, yy, yz, zz
    val[6*i]   *= ezor*x*z*zeta*(-3*(y2+z2)+x2*(-2+zr));
    val[6*i+1] *= -ezor*y*z*zeta*(y2+z2-r*x2*zeta);
    val[6*i+2] *= ezor*(-x2*x2*zeta+(y2+z2)*(r-z2*zeta)+x2*(r-y2*zeta+z2*zeta*(-1+zr)));
    val[6*i+3] *= -ezor*x*z*zeta*(x2+z2-r*y2*zeta);
    val[6*i+4] *= -ezor*x*y*zeta*(x2+y2-r*z2*zeta);
    val[6*i+5] *= ezor*x*z*zeta*(-3*x2-3*y2+z2*(-2+zr));
  }
  return;
}

void get_h_3dx2y2(int gs, double* grid, double* val, double zeta)
{
  double zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3]; double r3 = r*r*r;
    double ezr = exp(-zeta*r);
    double ezor = ezr/r3;

    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double xy = x*y; double xz = x*z; double yz = y*z;
    double zr = zeta*r;

   //xx, xy, xz, yy, yz, zz
    val[6*i]   *= ezor*(x2*x2*zeta*(-4+zr)+(y2+z2)*(2*r+y2*zeta)-x2*(-2*r+5*z2*zeta+y2*zeta*(5+zr)));
    val[6*i+1] *= ezor*x*(x-y)*y*(x+y)*zeta*(1+zr);
    val[6*i+2] *= ezor*x*z*zeta*(-2*z2+x2*(-1+zr)-y2*(3+zr));
    val[6*i+3] *= ezor*(-2*r*z2-x2*x2*zeta+y2*y2*zeta*(4-zr)+y2*(-2*r+5*z2*zeta)+x2*(-2*r-z2*zeta+y2*zeta*(5+zr)));
    val[6*i+4] *= ezor*y*z*zeta*(2*z2+y2*(1-zr)+x2*(3+zr));
    val[6*i+5] *= -ezor*(x-y)*(x+y)*zeta*(x2+y2-r*z2*zeta);
  }
  return;
}

void get_h_4fm3(int gs, double* grid, double* val, double zeta)
{
  double zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3]; double r3 = r*r*r;
    double ezr = exp(-zeta*r);
    double ezor = ezr/r3;

    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double xy = x*y; double xz = x*z; double yz = y*z;
    double zr = zeta*r;

   //xx, xy, xz, yy, yz, zz
    val[6*i]   *= ezor*y*(3*x2*x2*zeta*(-4+zr)+(y2+z2)*(6*r+y2*zeta)-x2*(-6*r+15*z2*zeta+y2*zeta*(15+zr)));
    val[6*i+1] *= ezor*x*(6*r*z2-3*x2*x2*zeta-y2*y2*zeta*(4+zr)+y2*(6*r-3*z2*zeta)+3*x2*(2*r-z2*zeta+y2*zeta*(-1+zr)));
    val[6*i+2] *= ezor*x*y*z*zeta*(-6*z2+3*x2*(-1+zr)-y2*(7+zr));
    val[6*i+3] *= -ezor*y*(6*r*z2+9*x2*x2*zeta+y2*y2*zeta*(-6+zr)+y2*(6*r-7*z2*zeta)+x2*(6*r+9*z2*zeta-y2*zeta*(1+3*zr)));
    val[6*i+4] *= ezor*z*zeta*(-3*x2*x2+3*y2*z2+y2*y2*(2-zr)+3*x2*(-z2+y2*(1+zr)));
    val[6*i+5] *= ezor*y*(-3*x2+y2)*zeta*(x2+y2-r*z2*zeta);
  }
  return;
}

void get_h_4fm2(int gs, double* grid, double* val, double zeta)
{
  double zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3]; double r3 = r*r*r;
    double ezr = exp(-zeta*r);
    double ezor = ezr/r3;

    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double xy = x*y; double xz = x*z; double yz = y*z;
    double zr = zeta*r;

   //xx, xy, xz, yy, yz, zz
    val[6*i]   *= ezor*x*y*z*zeta*(-3*(y2+z2)+x2*(-2+zr));
    val[6*i+1] *= ezor*z*(-x2*x2*zeta+(y2+z2)*(r-y2*zeta)+x2*(r-z2*zeta+y2*zeta*(-1+zr)));
    val[6*i+2] *= ezor*y*(-x2*x2*zeta+(y2+z2)*(r-z2*zeta)+x2*(r-y2*zeta+z2*zeta*(-1+zr)));
    val[6*i+3] *= ezor*x*y*z*zeta*(-3*x2-3*z2+y2*(-2+zr));
    val[6*i+4] *= ezor*x*(-y2*y2*zeta+z2*(r-z2*zeta)+x2*(r-(y2+z2)*zeta)+y2*(r+z2*zeta*(-1+zr)));
    val[6*i+5] *= ezor*x*y*z*zeta*(-3*x2-3*y2+z2*(-2+zr));
  }
  return;
}

void get_h_4fm1(int gs, double* grid, double* val, double zeta)
{
  double zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3]; double r3 = r*r*r;
    double ezr = exp(-zeta*r);
    double ezor = ezr/r3;

    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double xy = x*y; double xz = x*z; double yz = y*z;
    double zr = zeta*r;

   //xx, xy, xz, yy, yz, zz
    val[6*i]   *= ezor*y*(x2*x2*zeta*(4-zr)+x2*(-2*r+y2*zeta*(5-zr)+z2*zeta*(5+4*zr))+(y2+z2)*(y2*zeta-2*(r+2*z2*zeta)));
    val[6*i+1] *= ezor*x*(x2*x2*zeta+y2*y2*zeta*(4-zr)-2*z2*(r+2*z2*zeta)-x2*(2*r+3*z2*zeta+y2*zeta*(-5+zr))+y2*(-2*r+z2*zeta*(5+4*zr)));
    val[6*i+2] *= -ezor*x*y*z*zeta*(7*y2+2*z2+r*(y2-4*z2)*zeta+x2*(7+zr));
    val[6*i+3] *= -ezor*y*(-3*x2*x2*zeta+y2*y2*zeta*(-6+zr)+6*z2*(r+2*z2*zeta)+y2*(6*r+z2*zeta*(1-4*zr))+x2*(6*r+9*z2*zeta+y2*zeta*(-9+zr)));
    val[6*i+4] *= -ezor*z*(-8*r*z2-x2*x2*zeta+4*z2*z2*zeta+y2*y2*zeta*(6+zr)+y2*(-8*r+z2*zeta*(5-4*zr))+x2*(-8*r+3*z2*zeta+y2*zeta*(5+zr)));
    val[6*i+5] *= ezor*y*(8*r*z2+x2*x2*zeta+y2*y2*zeta+4*z2*z2*zeta*(-4+zr)+x2*(8*r+2*y2*zeta-z2*zeta*(20+zr))-y2*(-8*r+z2*zeta*(20+zr)));
  }
  return;
}

void get_h_4f0(int gs, double* grid, double* val, double zeta)
{
  double zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3]; double r3 = r*r*r;
    double ezr = exp(-zeta*r);
    double ezor = ezr/r3;

    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double xy = x*y; double xz = x*z; double yz = y*z;
    double zr = zeta*r;

   //xx, xy, xz, yy, yz, zz
    val[6*i]   *= -ezor*z*(3*x2*x2*zeta*(-4+zr)+(y2+z2)*(6*r-3*y2*zeta+2*z2*zeta)+x2*(6*r+3*y2*zeta*(-5+zr)-z2*zeta*(15+2*zr)));
    val[6*i+1] *= ezor*x*y*z*zeta*(9*y2+14*z2+r*(-3*y2+2*z2)*zeta+x2*(9-3*zr));
    val[6*i+2] *= ezor*x*(-6*r*z2+3*x2*x2*zeta+3*y2*y2*zeta+2*z2*z2*zeta*(1+zr)-3*r*y2*(2+z2*zt2)-3*x2*(-2*y2*zeta+r*(2+z2*zt2)));
    val[6*i+3] *= -ezor*z*(6*r*z2-3*x2*x2*zeta+2*z2*z2*zeta+3*y2*y2*zeta*(-4+zr)+x2*(6*r-z2*zeta+3*y2*zeta*(-5+zr))+y2*(6*r-z2*zeta*(15+2*zr)));
    val[6*i+4] *= ezor*y*(-6*r*z2+3*x2*x2*zeta+3*y2*y2*zeta+2*z2*z2*zeta*(1+zr)-3*r*y2*(2+z2*zt2)-3*x2*(-2*y2*zeta+r*(2+z2*zt2)));
    val[6*i+5] *= ezor*z*(9*x2*x2*zeta+9*y2*y2*zeta+2*(6*r*z2+z2*z2*zeta*(-6+zr))+y2*(12*r-z2*zeta*(8+3*zr))+x2*(12*r+18*y2*zeta-z2*zeta*(8+3*zr)));
  }
  return;
}

void get_h_4fp1(int gs, double* grid, double* val, double zeta)
{
  double zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3]; double r3 = r*r*r;
    double ezr = exp(-zeta*r);
    double ezor = ezr/r3;

    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double xy = x*y; double xz = x*z; double yz = y*z;
    double zr = zeta*r;

   //xx, xy, xz, yy, yz, zz
    val[6*i]   *= -ezor*x*(x2*x2*zeta*(-6+zr)+3*(y2+z2)*(2*r-y2*zeta+4*z2*zeta)+x2*(6*r+z2*zeta*(1-4*zr)+y2*zeta*(-9+zr)));
    val[6*i+1] *= ezor*y*(x2*x2*zeta*(4-zr)+x2*(-2*r+y2*zeta*(5-zr)+z2*zeta*(5+4*zr))+(y2+z2)*(y2*zeta-2*(r+2*z2*zeta)));
    val[6*i+2] *= -ezor*z*(x2*x2*zeta*(6+zr)-(y2+z2)*(8*r+y2*zeta-4*z2*zeta)+x2*(-8*r+z2*zeta*(5-4*zr)+y2*zeta*(5+zr)));
    val[6*i+3] *= ezor*x*(x2*x2*zeta+y2*y2*zeta*(4-zr)-2*z2*(r+2*z2*zeta)-x2*(2*r+3*z2*zeta+y2*zeta*(-5+zr))+y2*(-2*r+z2*zeta*(5+4*zr)));
    val[6*i+4] *= -ezor*x*y*z*zeta*(7*y2+2*z2+r*(y2-4*z2)*zeta+x2*(7+zr));
    val[6*i+5] *= ezor*x*(8*r*z2+x2*x2*zeta+y2*y2*zeta+4*z2*z2*zeta*(-4+zr)+x2*(8*r+2*y2*zeta-z2*zeta*(20+zr))-y2*(-8*r+z2*zeta*(20+zr)));
  }
  return;
}

void get_h_4fp2(int gs, double* grid, double* val, double zeta)
{
  double zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3]; double r3 = r*r*r;
    double ezr = exp(-zeta*r);
    double ezor = ezr/r3;

    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double xy = x*y; double xz = x*z; double yz = y*z;
    double zr = zeta*r;

   //xx, xy, xz, yy, yz, zz
    val[6*i]   *= ezor*z*(x2*x2*zeta*(-4+zr)+(y2+z2)*(2*r+y2*zeta)-x2*(-2*r+5*z2*zeta+y2*zeta*(5+zr)));
    val[6*i+1] *= ezor*x*(x-y)*y*(x+y)*z*zeta*(1+zr);
    val[6*i+2] *= ezor*x*(-x2*x2*zeta+y2*y2*zeta+2*z2*(r-z2*zeta)+x2*(2*r+z2*zeta*(-2+zr))-y2*(-2*r+z2*zeta*(2+zr)));
    val[6*i+3] *= -ezor*z*(2*r*z2+x2*x2*zeta+y2*y2*zeta*(-4+zr)+y2*(2*r-5*z2*zeta)+x2*(2*r+z2*zeta-y2*zeta*(5+zr)));
    val[6*i+4] *= ezor*y*(-2*r*z2-x2*x2*zeta+y2*y2*zeta+2*z2*z2*zeta-y2*(2*r+z2*zeta*(-2+zr))+x2*(-2*r+z2*zeta*(2+zr)));
    val[6*i+5] *= -ezor*(x-y)*(x+y)*z*zeta*(3*x2+3*y2+z2*(2-zr));
  }
  return;
}

void get_h_4fp3(int gs, double* grid, double* val, double zeta)
{
  double zt2 = zeta*zeta;

 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:6*gs]) ////async(tid)
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i]; double y = grid[6*i+1]; double z = grid[6*i+2];
    double r = grid[6*i+3]; double r3 = r*r*r;
    double ezr = exp(-zeta*r);
    double ezor = ezr/r3;

    double x2 = x*x; double y2 = y*y; double z2 = z*z;
    double xy = x*y; double xz = x*z; double yz = y*z;
    double zr = zeta*r;

   //xx, xy, xz, yy, yz, zz
    val[6*i]   *= ezor*(x2*x2*x*zeta*(-6+zr)+3*x*(y2+z2)*(2*r+3*y2*zeta)-x2*x*(-6*r+7*z2*zeta+y2*zeta*(1+3*zr)));
    val[6*i+1] *= ezor*y*(x2*x2*zeta*(4+zr)+3*(y2+z2)*(-2*r+y2*zeta)+3*x2*(-2*r+z2*zeta+y2*zeta*(1-zr)));
    val[6*i+2] *= ezor*z*zeta*(3*y2*(y2+z2)+x2*x2*(-2+zr)-3*x2*(z2+y2*(1+zr)));
    val[6*i+3] *= -ezor*x*(x2*x2*zeta+x2*(6*r+z2*zeta-y2*zeta*(15+zr))+3*(2*r*z2+y2*y2*zeta*(-4+zr)+y2*(2*r-5*z2*zeta)));
    val[6*i+4] *= ezor*x*y*z*zeta*(6*z2+y2*(3-3*zr)+x2*(7+zr));
    val[6*i+5] *= -ezor*x*(x2-3*y2)*zeta*(x2+y2-r*z2*zeta);
  }
  return;
}

void eval_h(int gs, double* grid, double* val, int n1, int l1, int m1, double zeta1)
{
  if (n1==1)
    return get_h_1s(gs,grid,val,zeta1);
  else if (n1==2)
  {
    if (l1==0)
      return get_h_2s(gs,grid,val,zeta1);
    else
    {
      if (m1==1)
        return get_h_2px(gs,grid,val,zeta1);
      else if (m1==-1)
        return get_h_2py(gs,grid,val,zeta1);
      else
        return get_h_2pz(gs,grid,val,zeta1);
    }
  }
  else if (n1==3)
  {
    printf("  WARNING: n=3 Hessian is being tested \n");
    if (l1==2)
    {
      if (m1==-2)
        return get_h_3dxy(gs,grid,val,zeta1);
      else if (m1==-1)
        return get_h_3dyz(gs,grid,val,zeta1);
      else if (m1==0)
        return get_h_3dz2(gs,grid,val,zeta1);
      else if (m1==1)
        return get_h_3dxz(gs,grid,val,zeta1);
      else if (m1==2)
        return get_h_3dx2y2(gs,grid,val,zeta1);
    }
    else
      printf("  ERROR: Hessian only for n=l+1 \n");
  }
  else if (n1==4)
  {
    printf("  WARNING: n=4 Hessian being tested \n");
    if (l1==3)
    {
      if (m1==-3)
        return get_h_4fm3(gs,grid,val,zeta1);
      else if (m1==-2)
        return get_h_4fm2(gs,grid,val,zeta1);
      else if (m1==-1)
        return get_h_4fm1(gs,grid,val,zeta1);
      else if (m1==0)
        return get_h_4f0(gs,grid,val,zeta1);
      else if (m1==1)
        return get_h_4fp1(gs,grid,val,zeta1);
      else if (m1==2)
        return get_h_4fp2(gs,grid,val,zeta1);
      else if (m1==3)
        return get_h_4fp3(gs,grid,val,zeta1);
    }
    else
      printf("  ERROR: Hessian only for n=l+1 \n");
  }
  else
  {
    printf("  ERROR: n>4, l>3 Hessian not available \n");
  }

  return;
}

void eval_h(int gs, float* grid, float* val, int n1, int l1, int m1, float zeta1)
{
  int gs6 = 6*gs;
  double* gridd = new double[gs6];
  double* vald = new double[gs];

  #pragma acc enter data create(gridd[0:gs6],vald[0:gs])
  eval_h(gs,gridd,vald,n1,l1,m1,zeta1);

  #pragma acc parallel loop present(val[0:gs],vald[0:gs])
  for (int j=0;j<gs;j++)
    val[j] = vald[j];

  #pragma acc exit data delete(vald[0:gs],gridd[0:gs6])

  delete [] gridd;
  delete [] vald;

  return;
}
