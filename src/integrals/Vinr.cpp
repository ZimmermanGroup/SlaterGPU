#include "Vinr.h"
#include "gamma.h"

#define USE_GAM 1
//fixed this
#define EVL64 1
#define RGLIMIT 0

//this contains norms and Inl(r) functions

//Notes:
//1. in Inr functions, exp could be expf
//2. need to carefully check double vs float
//3. need 6h derivative (untested)

// most of these ftns are *=
// but Yukawa potential is +=


void eval_ke(int gs, float* grid, float* val, int n, int l, float zeta)
{
  float f1 = n*(n-1) - l*(l+1);
  float f2 = 2.f*n*zeta;
  float f3 = zeta*zeta;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float or1 = 1.f/r;
    float or2 = or1*or1;

    float term = f1*or2 - f2*or1 + f3;
    val[i] *= term;
  }

  return;
}

void eval_ke(int gs, double* grid, double* val, int n, int l, double zeta)
{
  double f1 = n*(n-1) - l*(l+1);
  double f2 = 2.*n*zeta;
  double f3 = zeta*zeta;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double or1 = 1./r;
    double or2 = or1*or1;

    double term = f1*or2 - f2*or1 + f3;
    val[i] *= term;
  }

  return;
}

void eval_ke_erf(int gs, float* grid, float* val, int n, int l, float zeta, float ef1)
{
 //use this on 1s only
  if (n>1) return eval_ke(gs,grid,val,n,l,zeta);

  float f1 = n*(n-1) - l*(l+1);
  float f2 = -2.f*n*zeta;
  float f3 = zeta*zeta;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float or1 = 1.f/r;
    //float or2 = or1*or1;
    float erf1 = erff(ef1*r);
 
    or1 *= erf1;
    //or2 *= erf1;

    //float term = f1*or2 - f2*or1 + f3;
    float term = f2*or1 + f3;
    val[i] *= term;
  }

  return;
}

void eval_ke3(int gs, float* grid, float* val, int n, int l, float zeta)
{
  float f1 = n*(n-1) - l*(l+1);
  float f2 = 2.f*n*zeta;
  float f3 = zeta*zeta;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float or1 = 1.f/r;
    float or2 = or1*or1;

    float term = f1*or2 - f2*or1 + f3;
    val[3*i+0] *= term;
    val[3*i+1] *= term;
    val[3*i+2] *= term;
  }

  return;
}

void eval_dke(int gs, float* grid, float* val, int n, int l, float zeta)
{
  double f1 = n*(n-1) - l*(l+1);
  double f2 = 2.*n*zeta;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double or1 = 1./r;
    double or2 = or1*or1;

    double term = -2.*f1*or2*or2 + f2*or2*or1;
    val[3*i+0] *= x*term;
    val[3*i+1] *= y*term;
    val[3*i+2] *= z*term;
  }

  return;
}

void eval_ne(int gs, float* grid, float** val, int s1, int s2, int natoms, int* atno, float* coords, float A0, float B0, float C0)
{
  int ng = s2-s1;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:ng][0:gs],atno[0:natoms],coords[0:3*natoms])
#endif
  for (int i=0;i<gs;i++)
  {
    float x1 = grid[6*i]+A0;
    float y1 = grid[6*i+1]+B0;
    float z1 = grid[6*i+2]+C0;

    float ne1 = 0.;
   #pragma acc loop reduction(+:ne1)
    for (int j=0;j<natoms;j++)
    {
      float Zeff = atno[j];
      float x2 = x1-coords[3*j];
      float y2 = y1-coords[3*j+1];
      float z2 = z1-coords[3*j+2];
      float Ra = sqrtf(x2*x2+y2*y2+z2*z2);

      ne1 += Zeff/Ra;
    }
   #pragma acc loop
    for (int k=0;k<ng;k++)
      val[k][i] *= ne1;
  }

  return;
}

void eval_ne_3(int gs, float* grid, float** val, int s1, int s2, int natoms, int* atno, float* coords, float A0, float B0, float C0)
{
  int ng = s2-s1;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:ng][0:3*gs],atno[0:natoms],coords[0:3*natoms])
#endif
  for (int i=0;i<gs;i++)
  {
    float x1 = grid[6*i]+A0;
    float y1 = grid[6*i+1]+B0;
    float z1 = grid[6*i+2]+C0;

    float ne1 = 0.;
   #pragma acc loop reduction(+:ne1)
    for (int j=0;j<natoms;j++)
    {
      float Zeff = atno[j];
      float x2 = x1-coords[3*j];
      float y2 = y1-coords[3*j+1];
      float z2 = z1-coords[3*j+2];
      float Ra = sqrtf(x2*x2+y2*y2+z2*z2);

      ne1 += Zeff/Ra;
    }
   #pragma acc loop
    for (int k=0;k<ng;k++)
   #pragma acc loop
    for (int m=0;m<3;m++)
      val[k][3*i+m] *= ne1;
  }

  return;
}

void eval_d_ne_3(int gs, float* grid, float** val, int s1, int s2, int natoms, int* atno, float* coords, float A0, float B0, float C0)
{
  int ng = s2-s1;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:ng][0:3*gs],atno[0:natoms],coords[0:3*natoms])
#endif
  for (int i=0;i<gs;i++)
  {
    float x1 = grid[6*i]+A0;
    float y1 = grid[6*i+1]+B0;
    float z1 = grid[6*i+2]+C0;

    float ne1 = 0.;
   #pragma acc loop reduction(+:ne1)
    for (int j=0;j<natoms;j++)
    {
      float Zeff = atno[j];
      float x2 = x1-coords[3*j];
      float y2 = y1-coords[3*j+1];
      float z2 = z1-coords[3*j+2];
      float Ra = sqrtf(x2*x2+y2*y2+z2*z2);

      ne1 += Zeff/Ra;
    }
   #pragma acc loop
    for (int k=0;k<ng;k++)
   #pragma acc loop
    for (int m=0;m<3;m++)
      val[k][3*i+m] *= ne1;
  }

  return;
}


//derivatives of Inr

void eval_inr_6h_r1d(int gs, float* grid, float* val, float zeta)
{
 //test this
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz*oz2;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;
  double oz9 = oz4*oz5;
  double oz10 = oz5*oz5;
  double oz11 = oz6*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r; double r3 = r2*r; double r4 = r2*r2; double r5 = r4*r;
    double ezr = exp(-zeta*r);
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double f1 = 2874009600./r4/r3*oz7*oz6;
    double v1 = ezr*(f1*r*zeta + 1437004800.*oz11/r5 + 479001600.*oz10/r4 + 119750400.*oz9/r3 + 23950080.*oz8/r2 + 3991680.*oz7/r + 570240.*oz6 + 71280.*r*oz5 + 7920.*r2*oz4 + 792.*r3*oz3 + 77.*r4*oz2 + 11.*r5*oz);
    v1 += (ezr-1.)*f1;

   #if RGLIMIT
    if (r<0.1) v1 = 0.;
   #endif

    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_5g_r1d(int gs, float* grid, float* val, float zeta)
{
 //careful, test this
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz*oz2;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;
  double oz9 = oz4*oz5;
  double oz10 = oz5*oz5;
  double oz11 = oz6*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r; double r3 = r2*r; double r4 = r2*r2;
    double ezr = exp(-zeta*r);
    double xor1 = x/r; double yor = y/r; double zor = z/r;

   //exp(-x) = 1-x;
    double f1 = 18144000./r3/r3*oz11;
    double v1 = ezr*(f1*r*zeta + 9072000.*oz9/r4 + 3024000.*oz8/r3 + 756000.*oz7/r2 + 151200.*oz6/r + 25200.*oz5 + 3600.*r*oz4 + 450.*r2*oz3 + 54.*r3*oz2 + 9.*r4*oz);
    //double f1 = 18144000.*oz11;
    //double v1 = ezr*f1*r5
    //printf("  r: %12.10f  f1: %12.10f v1: %12.10f  em1: %12.10f \n",r,f1,v1,(ezr-1.)*f1);
    v1 += (ezr-1.)*f1;

   #if RGLIMIT
    if (r<0.1) v1 = 0.;
   #endif

    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_6f_r1d(int gs, float* grid, float* val, float zeta)
{
 //careful, test this
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz*oz2;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;
  double oz9 = oz4*oz5;
  double oz10 = oz5*oz5;
  double oz11 = oz5*oz6;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r; double r3 = r2*r; double r4 = r2*r2; double r5 = r3*r2;
    double ezr = exp(-zeta*r);
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double f1 = 14515200./r5*oz11;
    double v1 = ezr*(f1*r*zeta + 7257600.*oz9/r3 + 2419200.*oz8/r2 + 604800.*oz7/r + 120960.*oz6 + 20160.*r*oz5 + 2898.*r2*oz4 + 378.*r3*oz3 + 49.*r4*oz2 + 7.*r5*oz);
    v1 += (ezr-1.)*f1;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_5f_r1d(int gs, float* grid, float* val, float zeta)
{
 //careful, test this
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz*oz2;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;
  double oz9 = oz4*oz5;
  double oz10 = oz5*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r; double r3 = r2*r; double r4 = r2*r2;
    double ezr = exp(-zeta*r);
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double f1 = 1451520./r2/r3*oz10;
    double v1 = ezr*(f1*r*zeta + 725760.*oz8/r3 + 241920.*oz7/r2 + 60480.*oz6/r + 12096.*oz5 + 2016.*r*oz4 + 294.*r2*oz3 + 42.*r3*oz2 + 7.*r4*oz);
    v1 += (ezr-1.)*f1;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_4f_r1d(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz*oz2;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;
  double oz9 = oz4*oz5;

  //#pragma acc update self(grid[0:6*gs],val[0:3*gs])

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r; double r3 = r2*r;
    double ezr = exp(-zeta*r);
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double f1 = 161280./r2/r3*oz9;
    double v1 = ezr*(f1*r*zeta + 80640.*oz7/r3 + 26880.*oz6/r2 + 6720.*oz5/r + 1344.*oz4 + 224.*r*oz3 + 35.*r2*oz2 + 7.*r3*oz);
    v1 += (ezr-1.)*f1;
    //printf("  v1s: %2.1e %2.1e %2.1e %2.1e %2.1e %2.1e %2.1e %2.1e \n",f1*r*zeta,80640.*oz7/r3,26880.*oz6/r2,6720.*oz5/r,1344.*oz4,224.*r*oz3,35.*r2*oz2,7.*r3*oz);
    //printf("  xyzr: %8.5f %8.5f %8.5f %8.5f  ezr: %10.8f zeta: %5.3f  f1/v1: %2.1e %2.1e \n",x,y,z,r,ezr,zeta,f1,v1);
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;

    //val[i] *= -ezr*(720.f*oz7r3 + 720.f*oz6r2 + 360.f*oz5r + 120.f*oz4 + 30.f*r*oz3 + 5.f*r2*oz2) + 720.f*oz7r3;
  }

  return;
}

void eval_inr_7d_r1d(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz*oz2;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;
  double oz9 = oz5*oz4;
  double oz11 = oz5*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r; double r3 = r2*r; double r4 = r2*r2;
    double ezr = exp(-zeta*r);
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double or4z11 = 10886400./r4*oz11;
    double v1 = ezr*(or4z11 + or4z11*r*zeta + 5443200.*oz9/r2 + 1814400.*oz8/r + 453600.*oz7 + 90960.*r*oz6 + 15360.*r2*oz5 + 2280.*oz4*r3 + 310.*r4*oz3 + 40.*r3*r2*oz2 + 5.*r3*r3*oz) - or4z11;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_6d_r1d(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz*oz2;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;
  double oz9 = oz5*oz4;
  double oz10 = oz5*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r; double r3 = r2*r; double r4 = r2*r2;
    double ezr = exp(-zeta*r);
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double or4z10 = 1088640./r4*oz10;
    double v1 = ezr*(or4z10 + or4z10*r*zeta + 544320.*oz8/r2 + 181440.*oz7/r + 45360.*oz6 + 9120.*r*oz5 + 1560.*r2*oz4 + 240.*oz3*r3 + 35.*r4*oz2 + 5.*r3*r2*oz) - or4z10;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_5d_r1d(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz*oz2;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;
  //double oz8 = oz4*oz4;
  double oz9 = oz5*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r; double r3 = r2*r; double r4 = r2*r2;
    double ezr = exp(-zeta*r);
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double or4z9 = 120960./r4*oz9;
    double v1 = ezr*(or4z9 + or4z9*r*zeta + 60480.*oz7/r2 + 20160.*oz6/r + 5040.*oz5 + 1020.*r*oz4 + 180.*r2*oz3 + 30.*oz2*r3 + 5.*r4*oz) - or4z9;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;

    //val[i] *= -ezr*(720.f*oz7r3 + 720.f*oz6r2 + 360.f*oz5r + 120.f*oz4 + 30.f*r*oz3 + 5.f*r2*oz2) + 720.f*oz7r3;
  }

  return;
}

void eval_inr_4d_r1d(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz*oz2;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;
  //double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r; double r3 = r2*r;
    double ezr = exp(-zeta*r);
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double for4z8 = 15120.*oz8/r2/r2;
    double v1 = ezr*(for4z8 + for4z8*r*zeta + 7560.*oz6/r2 + 2520.*oz5/r + 630.*oz4 + 130.*r*oz3 + 25.*r2*oz2 + 5.*oz*r3) - for4z8;

    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;

    //val[i] *= -ezr*(720.f*oz7r3 + 720.f*oz6r2 + 360.f*oz5r + 120.f*oz4 + 30.f*r*oz3 + 5.f*r2*oz2) + 720.f*oz7r3;
  }

  return;
}

void eval_inr_3d_r1d(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz*oz2;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz7 = oz3*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r; double r3 = r2*r;
    double ezr = exp(-zeta*r);
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double tor4z7 = 2160./r2/r2*oz7;
    double v1 = ezr*(tor4z7 + tor4z7*r*zeta + 1080.*oz5/r2 + 360.*oz4/r + 90.*oz3 + 20.*r*oz2 + 5.*r2*oz) - tor4z7;
    //v1 += (ezr-1.)*tor4z7);

    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;

    //val[i] *= -ezr*(720.f*oz7r3 + 720.f*oz6r2 + 360.f*oz5r + 120.f*oz4 + 30.f*r*oz3 + 5.f*r2*oz2) + 720.f*oz7r3;
  }

  return;
}

void eval_inr_7p_r1d(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz3*oz2;
  double oz6 = oz3*oz3;
  double oz7 = oz4*oz3;
  double oz8 = oz4*oz4;
  double oz9 = oz4*oz5;
  double oz10 = oz5*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r; double r3 = r2*r;
    double ezr = exp(-zeta*r);
    double tor3z10 = 725760.*oz10/r3;
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double v1 = ezr*(tor3z10 + tor3z10*r*zeta + 362880.*oz8/r + 121680.*oz7 + 30960.*oz6*r + 6408.*oz5*r2 + 1128.*oz4*r3 + 174.*oz3*r2*r2 + 24.*oz2*r3*r2 + 3.*r3*r3*oz) - tor3z10;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_6p_r1d(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz3*oz2;
  double oz6 = oz3*oz3;
  double oz7 = oz4*oz3;
  double oz8 = oz4*oz4;
  double oz9 = oz4*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r; double r3 = r2*r;
    double ezr = exp(-zeta*r);
    double tor3z9 = 80640.*oz9/r3;
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double v1 = ezr*(tor3z9 + tor3z9*r*zeta + 40320.*oz7/r + 13560.*oz6 + 3480.*oz5*r + 732.*oz4*r2 + 132.*oz3*r3 + 21.*oz2*r2*r2 + 3.*oz*r3*r2) - tor3z9;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_5p_r1d(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz3*oz2;
  double oz6 = oz3*oz3;
  double oz7 = oz4*oz3;
  double oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r; double r3 = r2*r;
    double ezr = exp(-zeta*r);
    double tor3z8 = 10080.*oz8/r3;
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double v1 = ezr*(tor3z8 + tor3z8*r*zeta + 5040.*oz6/r + 1704.*oz5 + 444.*oz4*r + 96.*oz3*r2 + 18.*oz2*r3 + 3.*oz*r2*r2) - tor3z8;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;

    //val[i] *= -ezr*(120.f*oz6r2 + 120.f*oz5r + 60.f*oz4 + 18.f*r*oz3 + 3.f*r2*oz2) + 120.f*oz6r2;
  }

  return;
}

void eval_inr_4p_r1d(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz3*oz2;
  double oz6 = oz3*oz3;
  double oz7 = oz4*oz3;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r; double r3 = r2*r;
    double ezr = exp(-zeta*r);
    double for3z7 = 1440.*oz7/r3;
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double v1 = (ezr*(for3z7 + for3z7*r*zeta + 720.*oz5/r + 246.*oz4 + 66.*oz3*r + 15.*oz2*r2 + 3.*oz*r3) - for3z7);
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;

    //val[i] *= -ezr*(120.f*oz6r2 + 120.f*oz5r + 60.f*oz4 + 18.f*r*oz3 + 3.f*r2*oz2) + 120.f*oz6r2;
  }

  return;
}

void eval_inr_3p_r1d(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz3*oz2;
  double oz6 = oz3*oz3;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r;
    double ezr = exp(-zeta*r);
    double tor3z6 = 240.*oz6/r2/r;
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double v1 = (ezr*(tor3z6 + tor3z6*r*zeta + 120.*oz4/r + 42.*oz3 + 12.*oz2*r + 3.*oz*r2) - tor3z6);
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;

    //printf("  3p_r1d.  r: %9.6f  v1(b): %12.10f  add: %12.10f  v1: %12.10f \n",r,v1-(ezr-1.)*tor3z6,(ezr-1.)*tor3z6,v1);

    //val[i] *= -ezr*(120.f*oz6r2 + 120.f*oz5r + 60.f*oz4 + 18.f*r*oz3 + 3.f*r2*oz2) + 120.f*oz6r2;
  }

  return;
}

void eval_inr_2p_r1d(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz3*oz2;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    double x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double r2 = r*r;
    double ezr = exp(-zeta*r);
    double foz5r3 = 48.*oz5/r2/r;
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double v1 = (ezr*(foz5r3 + foz5r3*r*zeta + 24.*oz3/r + 9.*oz2 + 3.*r*oz) - foz5r3);
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;

    //val[i] *= -ezr*(toz5r2+24.f*oz4/r+12.f*oz3+3.f*r*oz2) + toz5r2;
  }

  return;
}


void eval_inr_8s_r1d(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;
  double oz9 = oz4*oz5;
  double oz10 = oz5*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double r2 = r*r; double r3 = r2*r; double r4 = r2*r2; double r5 = r3*r2;
    double for2z10 = 362880.*oz10/r2;
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double v1 = ezr*(for2z10 + for2z10*zeta*r + 181440.*oz8 + 60480.*r*oz7 + 15120.*r2*oz6 + 3024.*r3*oz5 + 504.*r4*oz4 + 72.*r5*oz3 + 9.*r3*r3*oz2 + r3*r4*oz) - for2z10;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_7s_r1d(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;
  double oz9 = oz4*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double r2 = r*r; double r3 = r2*r; double r4 = r2*r2; double r5 = r3*r2;
    double for2z9 = 40320.*oz9/r2;
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double v1 = ezr*(for2z9 + for2z9*zeta*r + 20160.*oz7 + 6720.*r*oz6 + 1680.*r2*oz5 + 336.*r3*oz4 + 56.*r4*oz3 + 8.*r5*oz2 + r3*r3*oz) - for2z9;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_6s_r1d(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double r2 = r*r; double r3 = r2*r; double r4 = r2*r2; double r5 = r3*r2;
    double for2z8 = 5040.*oz8/r2;
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double v1 = ezr*(for2z8 + for2z8*zeta*r + 2520.*oz6 + 840.*r*oz5 + 210.*r2*oz4 + 42.*r3*oz3 + 7.*r4*oz2 + r5*oz) - for2z8;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_5s_r1d(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double r2 = r*r; double r3 = r2*r; double r4 = r2*r2;
    double sor2z7 = 720.*oz7/r2;
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double v1 = (ezr*(sor2z7 + sor2z7*zeta*r + 360.*oz5 + 120.*r*oz4 + 30.*r2*oz3 + 6.*r3*oz2 + r4*oz) - sor2z7);
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_4s_r1d(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double r2 = r*r; double r3 = r2*r;
    double or2z6 = 120.*oz6/r2;
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double v1 = (ezr*(or2z6 + or2z6*zeta*r + 60.*oz4 + 20.*r*oz3 + 5.*r2*oz2 + r3*oz) - or2z6);
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_3s_r1d(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double r2 = r*r;
    double tor2z5 = 24.*oz5/r2;
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double v1 = (ezr*(tor2z5 + tor2z5*zeta*r + 12.*oz3 + 4.*r*oz2 + r2*oz) - tor2z5);
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;

    //val[i] *= -ezr*(24.f*oz5r + 18.f*oz4 + 6.f*roz3 + r2oz2) + 24.f*oz5r;
  }

  return;
}

void eval_inr_2s_r1d(int gs, float* grid, float* val, float zeta)
{
  double oz = 1.f/zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double ezr = expf(-zeta*r);
    double soz4r2 = 6.*oz4/r/r;
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double v1 = (ezr*(soz4r2 + 6.*oz3/r + 3.*oz2 + r*oz) - soz4r2);
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_1s_r1d(int gs, float* grid, float* val, float zeta)
{
  //float oz = 1.f/zeta;
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x = grid[6*i+0]; float y = grid[6*i+1]; float z = grid[6*i+2];
    float r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double tz3r2 = 2.*oz3/r/r;
    double xor1 = x/r; double yor = y/r; double zor = z/r;

    double v1 = (ezr*(tz3r2 + 2.*oz2/r + oz) - tz3r2);
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}


////r1 functions////

void eval_inr_6h_r1(int gs, float* grid, float* val, float zeta)
{
  //float oz = 1.f/zeta;
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz3*oz;
  double oz5 = oz4*oz;
  double oz6 = oz5*oz;
  double oz7 = oz6*oz;
  double oz8 = oz7*oz;
  double oz9 = oz8*oz;
  double oz10 = oz9*oz;
  double oz11 = oz10*oz;
  double oz12 = oz11*oz;
  double oz13 = oz12*oz;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];

    double r2 = r*r; double r3 = r2*r; double r4 = r2*r2; double r5 = r2*r3; double r6 = r3*r3;

#if 1
    double ezr = exp(-zeta*r);
    double foz13r6 = 479001600.*oz13/r6;
  
    val[i] *= -ezr*(foz13r6 + 479001600.*oz12/r5 + 239500800.*oz11/r4 + 79833600.*oz10/r3 + 19958400.*oz9/r2
                    + 3991680.*oz8/r + 665280.*oz7 + 95040.*r*oz6 + 11880.*r2*oz5 + 1320.*r3*oz4 + 132.*r4*oz3 + 11.*r5*oz2)
              + foz13r6;
  #if RGLIMIT
    if (r<0.1) val[i] = 0.;
  #endif
#else
    float ezr = expf(-zeta*r);
    float foz13r6 = 479001600.f*oz13/r6;
  
    val[i] *= -ezr*(foz13r6 + 479001600.f*oz12/r5 + 239500800.f*oz11/r4 + 79833600.f*oz10/r3 + 19958400.f*oz9/r2
                    + 3991680.f*oz8/r + 665280.f*oz7 + 95040.f*r*oz6 + 11880.f*r2*oz5 + 1320.f*r3*oz4 + 132.f*r4*oz3 + 11.f*r5*oz2)
              + foz13r6;
#endif
  }

  return;
}

void eval_inr_5g_r1(int gs, float* grid, float* val, float zeta)
{
  //float oz = 1.f/zeta;
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz3*oz;
  double oz5 = oz4*oz;
  double oz6 = oz5*oz;
  double oz7 = oz6*oz;
  double oz8 = oz7*oz;
  double oz9 = oz8*oz;
  double oz10 = oz9*oz;
  double oz11 = oz10*oz;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];

    double r2 = r*r; double r3 = r2*r; double r4 = r2*r2; double r5 = r2*r3;

#if 1
    double ezr = exp(-zeta*r);
    double toz11r5 = 3628800.*oz11/r5;

    val[i] *= -ezr*(toz11r5 + 3628800.*oz10/r4 + 1814400.*oz9/r3 + 604800.*oz8/r2 + 151200.*oz7/r
                    + 30240.*oz6 + 5040.*r*oz5 + 720.*r2*oz4 + 90.*r3*oz3 + 9.*r4*oz2) + toz11r5;
   #if RGLIMIT
    if (r<0.1) val[i] = 0.;
   #endif
#else
    float ezr = expf(-zeta*r);
    float toz11r5 = 3628800.f*oz11/r5;

    val[i] *= -ezr*(toz11r5 + 3628800.f*oz10/r4 + 1814400.f*oz9/r3 + 604800.f*oz8/r2 + 151200.f*oz7/r
                    + 30240.f*oz6 + 5040.f*r*oz5 + 720.f*r2*oz4 + 90.f*r3*oz3 + 9.f*r4*oz2) + toz11r5;
#endif
  }

  return;
}

void eval_inr_6f_r1(int gs, float* grid, float* val, float zeta)
{
 //untested
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz3*oz;
  double oz5 = oz4*oz;
  double oz6 = oz5*oz;
  double oz7 = oz6*oz;
  double oz8 = oz7*oz;
  double oz9 = oz8*oz;
  double oz10 = oz9*oz;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];

    double r2 = r*r; double r3 = r2*r; double r4 = r2*r2;

    double ezr = exp(-zeta*r);
    double toz11r4 = 3628800.*oz9/r4;

    val[i] *= oz2*(-ezr*(toz11r4 + toz11r4*r*zeta + 1814400.*oz7/r2 + 604800.*oz6/r + 151200.*oz5
                    + 30240.*oz4*r + 5040.*r2*oz3 + 714.*r3*oz2 + 84.*r4*oz + 7.*r3*r2) + toz11r4);
  }

  return;
}

void eval_inr_5f_r1(int gs, float* grid, float* val, float zeta)
{
 //untested
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz3*oz;
  double oz5 = oz4*oz;
  double oz6 = oz5*oz;
  double oz7 = oz6*oz;
  double oz8 = oz7*oz;
  double oz9 = oz8*oz;
  double oz10 = oz9*oz;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];

    double r2 = r*r; double r3 = r2*r; double r4 = r2*r2;

    double ezr = exp(-zeta*r);
    double toz10r4 = 362880.*oz8/r4;

    val[i] *= oz2*(-ezr*(toz10r4 + toz10r4*r*zeta + 181440.*oz6/r2 + 60480.*oz5/r + 15120.*oz4
                    + 3024.*oz3*r + 504.*r2*oz2 + 70.*r3*oz + 7.*r4) + toz10r4);
  }

  return;
}

void eval_inr_4f_r1(int gs, float* grid, float* val, float zeta)
{
  //float oz = 1./zeta;
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz*oz2;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;
  double oz9 = oz4*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];

    double r2 = r*r; double r3 = r*r2; double r4 = r2*r2;

    double oz9r4 = oz9/r4;
    double oz8r3 = oz8/r3;
    double oz7r2 = oz7/r2;
    double oz6r = oz6/r;
    double roz4 = r*oz4;
    double r2oz3 = r2*oz3;
    double r3oz2 = r3*oz2;

#if 1
    double ezr = exp(-zeta*r);
    double foz9r4 = 40320.*oz9r4;

    val[i] *= -ezr*(foz9r4 + 40320.*oz8r3 + 20160.*oz7r2 + 6720.*oz6r + 1680.*oz5 + 336.*roz4 + 56.*r2oz3 + 7.*r3oz2)+foz9r4;
#else  
    float ezr = expf(-zeta*r);
    float foz9r4 = 40320.f*oz9r4;

    val[i] *= -ezr*(foz9r4 + 40320.f*oz8r3 + 20160.f*oz7r2 + 6720.f*oz6r + 1680.f*oz5 + 336.f*roz4 + 56.f*r2oz3 + 7.f*r3oz2)+foz9r4;
#endif
  }

  return;
}

void eval_inr_7d_r1(int gs, float* grid, float* val, float zeta)
{
 //untested
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz*oz2;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz2*oz4;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;
  double oz9 = oz5*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double r2 = r*r; double r3 = r2*r; double r4 = r2*r2;

    double foz11r3 = 3628800.*oz9/r3;
    val[i] *= oz2*(foz11r3 - ezr*(foz11r3 + foz11r3*zeta*r + 1814400.*oz7/r + 604800.*oz6 + 151200.*r*oz5 + 30120.*r2*oz4 + 4920.*r3*oz3 + 660.*r4*oz2 + + 70.*r3*r2*oz + 5.*r3*r3));
  }
 
  return;
}

void eval_inr_6d_r1(int gs, float* grid, float* val, float zeta)
{
 //untested
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz*oz2;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz2*oz4;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;
  double oz9 = oz5*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double r2 = r*r; double r3 = r2*r; double r4 = r2*r2;

    double foz10r3 = 362880.*oz8/r3;
    val[i] *= oz2*(foz10r3 - ezr*(foz10r3 + foz10r3*zeta*r + 181440.*oz6/r + 60480.*oz5 + 15120.*r*oz4 + 3000.*r2*oz3 + 480.*r3*oz2 + 60.*r4*oz + 5.*r3*r2));
  }
 
  return;
}

void eval_inr_5d_r1(int gs, float* grid, float* val, float zeta)
{
 //untested
  //float oz = 1.f/zeta;
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz*oz2;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz2*oz4;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;
  double oz9 = oz5*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double r2 = r*r; double r3 = r2*r; double r4 = r2*r2;

    double foz9r3 = 40320.*oz7/r3;

    val[i] *= oz2*(foz9r3 - ezr*(foz9r3 + 40320.*oz6/r2 + 20160.*oz5/r + 6720.*oz4 + 1680.*r*oz3 + 330.*r2*oz2 + 50.*r3*oz + 5.*r4));
  }
 
  return;
}

void eval_inr_4d_r1(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz*oz2;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz2*oz4;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double r2 = r*r; double r3 = r2*r;

    double foz8r3 = 5040.*oz8/r3;

    val[i] *= foz8r3 - ezr*(foz8r3 + 5040.*oz7/r2 + 2520.*oz6/r + 840.*oz5 + 210.*r*oz4 + 40.*r2*oz3 + 5.*r3*oz2);
  }
 
  return;
}

void eval_inr_3d_r1(int gs, float* grid, float* val, float zeta)
{
#if 1
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz*oz2;
  double oz4 = oz2*oz2;
  double oz7 = oz3*oz4;
#else
  float oz = 1.f/zeta;
  float oz2 = oz*oz;
  float oz3 = oz*oz2;
  float oz4 = oz2*oz2;
  float oz7 = oz3*oz4;
#endif

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float r2 = r*r; float r3 = r2*r;

   #if 1
    double oz6r2 = oz2*oz4/r2;
    double oz5r = oz2*oz3/r;

    double ezr = exp(-zeta*r);
    double oz7r3 = oz7/r3;

    val[i] *= -ezr*(720.*oz7r3 + 720.*oz6r2 + 360.*oz5r + 120.*oz4 + 30.*r*oz3 + 5.*r2*oz2) + 720.*oz7r3;
   #else
    float oz6r2 = oz2*oz4/r2;
    float oz5r = oz2*oz3/r;

    float ezr = expf(-zeta*r);
    float oz7r3 = oz7/r3;

    val[i] *= -ezr*(720.f*oz7r3 + 720.f*oz6r2 + 360.f*oz5r + 120.f*oz4 + 30.f*r*oz3 + 5.f*r2*oz2) + 720.f*oz7r3;
   #endif
  }

  return;
}

void eval_inr_7p_r1(int gs, float* grid, float* val, float zeta)
{
 //untested

  //float oz = 1.f/zeta;
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz3*oz2;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double r2 = r*r;
    double r3 = r2*r;

    double fr2z10 = 362880.*oz8/r2;
    val[i] *= oz2*(fr2z10 - ezr*(fr2z10 + fr2z10*zeta*r + 181440.*oz6 + 59760.*oz5*r + 14400.*oz4*r2 + 2664.*oz3*r3 + 384.*r2*r2*oz2 + 42.*oz*r3*r2 + 3.*r3*r3));
  }

  return;
}

void eval_inr_6p_r1(int gs, float* grid, float* val, float zeta)
{
 //untested

  //float oz = 1.f/zeta;
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz3*oz2;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double r2 = r*r;
    double r3 = r2*r;

    double fr2z9 = 40320.*oz7/r2;
    val[i] *= oz2*(fr2z9 - ezr*(fr2z9 + fr2z9*zeta*r + 20160.*oz5 + 6600.*oz4*r + 1560.*oz3*r2 + 276.*oz2*r3 + 36.*oz*r2*r2 + 3.*r3*r2));
  }

  return;
}

void eval_inr_5p_r1(int gs, float* grid, float* val, float zeta)
{
 //untested

  //float oz = 1.f/zeta;
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz3*oz2;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double r2 = r*r;
    double r3 = r2*r;

#if 1
    double fr2z8 = 5040.*oz6/r2;
    val[i] *= oz2*(fr2z8 - ezr*(fr2z8 + fr2z8*zeta*r + 2520.*oz4 + 816.*oz3*r + 186.*oz2*r2 + 30.*oz*r3 + 3.*r2*r2));
#else
    float fr2z8 = 5040.f*oz8/r2;
    val[i] *= fr2z8 - ezr*(fr2z8 + 5040.f*oz7/r + 2520.f*oz6 + 816.f*oz5*r + 186.f*oz4*r2 + 30.f*oz3*r3 + 3.f*oz2*r2*r2);
#endif
  }

  return;
}

void eval_inr_4p_r1(int gs, float* grid, float* val, float zeta)
{
 #if 1
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz3*oz2;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4; 
 #else
  float oz = 1.f/zeta;
  float oz2 = oz*oz;
  float oz3 = oz2*oz;
  float oz4 = oz2*oz2;
  float oz5 = oz3*oz2;
  float oz6 = oz3*oz3;
  float oz7 = oz3*oz4; 
 #endif

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float r2 = r*r;
    float r3 = r2*r;

   #if 1
    double ezr = expf(-zeta*r);
    double soz7r2 = 720.*oz7/r2;

    double v1 = soz7r2 - ezr*(soz7r2 + 720.*oz6/r + 360.*oz5 + 114.*r*oz4 + 24.*r2*oz3 + 3.*r3*oz2);
    val[i] *= v1;
   #else
    float ezr = expf(-zeta*r);
    float soz7r2 = 720.f*oz7/r2;

    float v1 = soz7r2 - ezr*(soz7r2 + 720.f*oz6/r + 360.f*oz5 + 114.f*r*oz4 + 24.f*r2*oz3 + 3.f*r3*oz2);
    val[i] *= v1;
   #endif
  }
 
  return;
}

void eval_inr_3p_r1(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz3*oz2;
  double oz6 = oz3*oz3;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float r2 = r*r;

   #if 1
    double ezr = exp(-zeta*r);
    double oz5r = oz5/r;
    double oz6r2 = 120.*oz6/r2;

    val[i] *= -ezr*(oz6r2 + 120.*oz5r + 60.*oz4 + 18.*r*oz3 + 3.*r2*oz2) + oz6r2;
   #else
    float ezr = expf(-zeta*r);
    float oz5r = oz5/r;
    float oz6r2 = 120.f*oz6/r2;

    val[i] *= -ezr*(oz6r2 + 120.f*oz5r + 60.f*oz4 + 18.f*r*oz3 + 3.f*r2*oz2) + oz6r2;
   #endif
  }

  return;
}

void eval_inr_2p_r1(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz3*oz2;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float r2 = r*r;
    double ezr = exp(-zeta*r);
    double toz5r2 = 24.f*oz5/r2;

    val[i] *= -ezr*(toz5r2+24.f*oz4/r+12.f*oz3+3.f*r*oz2) + toz5r2;
  }

  return;
}

void eval_inr_8s_r1(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double r2 = r*r; double r3 = r2*r; double r4 = r2*r2;

    double forz10 = 362880./r*oz8;
    val[i] *= oz2*(forz10 - ezr*(forz10 + 322560.*oz7 + 141120.*r*oz6 + 40320.*r2*oz5 + 8400.*r3*oz4 + 1344.*r4*oz3 + 168.*r3*r2*oz2 + 16.*r3*r3*oz + r4*r3));
  }

  return;
}

void eval_inr_7s_r1(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double r2 = r*r; double r3 = r2*r;

    double forz9 = 40320./r*oz7;
    val[i] *= oz2*(forz9 - ezr*(forz9 + 35280.*oz6 + 15120.*r*oz5 + 4200.*r2*oz4 + 840.*r3*oz3 + 126.*r2*r2*oz2 + 14.*r3*r2*oz + r3*r3));
  }

  return;
}

void eval_inr_6s_r1(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;
  double oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double r2 = r*r; double r3 = r2*r;

    double forz8 = 5040./r*oz8;
    val[i] *= forz8 - ezr*(forz8 + 4320.*oz7 + 1800.*r*oz6 + 480.*r2*oz5 + 90.*r3*oz4 + 12.*r2*r2*oz3 + r2*r3*oz2);
  }

  return;
}

void eval_inr_5s_r1(int gs, float* grid, float* val, float zeta)
{
 //untested

  //float oz = 1.f/zeta;
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;
  double oz7 = oz3*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double r2 = r*r;

#if 1
    double sorz7 = 720./r*oz7;
    val[i] *= sorz7 - ezr*(sorz7 + 600.*oz6 + 240.*r*oz5 + 60.*r2*oz4 + 10.*r2*r*oz3 + r2*r2*oz2);
#else
    float sorz7 = 720.f/r*oz7;
    val[i] *= sorz7 - ezr*(sorz7 + 600.f*oz6 + 240.f*r*oz5 + 60.f*r2*oz4 + 10.f*r2*r*oz3 + r2*r2*oz2);
#endif
  }

  return;
}

void eval_inr_4s_r1(int gs, float* grid, float* val, float zeta)
{
  double oz = 1./zeta;
  double oz2 = oz*oz;
  double oz3 = oz2*oz;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;
  double oz6 = oz3*oz3;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    double ezr = exp(-zeta*r);

    double r2 = r*r;
    double r3 = r2*r;
    double ooz6r = 120.*oz6/r;

    val[i] *= ooz6r - ezr*(ooz6r + 96.*oz5 + 36.*r*oz4 + 8.*r2*oz3 + r3*oz2);
  }

  return;
}

void eval_inr_3s_r1(int gs, float* grid, float* val, float zeta)
{
  double oz2 = 1./zeta/zeta;
  double oz3 = oz2/zeta;
  double oz4 = oz2*oz2;
  double oz5 = oz2*oz3;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    double ezr = exp(-zeta*r);

    double oz5r = oz5/r;
    double r2oz2 = r*r*oz2;
    double roz3 = oz3*r;
    double toz5r = 24.*oz5r;

    val[i] *= -ezr*(toz5r + 18.*oz4 + 6.*roz3 + r2oz2) + toz5r;
  }

  return;
}

void eval_inr_2s_r1(int gs, float* grid, float* val, float zeta)
{
  double oz2 = 1./zeta/zeta;
  double oz3 = oz2/zeta;
  double oz4 = oz2*oz2;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
  #if 1
    double ezr = exp(-zeta*r);
    double roz2 = r*oz2;
    double soz4r = 6.*oz4/r;

    val[i] *= -ezr*(4.*oz3+roz2+soz4r)+soz4r;
  #else
    double ezr = exp(-zeta*r);
    double roz2 = r*oz2;
    double soz4r = 6.f*oz4/r;

    val[i] *= -ezr*(4.f*oz3+roz2+soz4r)+soz4r;
  #endif
  }

  return;
}

void eval_inr_1s_r1(int gs, float* grid, float* val, float zeta)
{
  double oz2 = 1./zeta/zeta;
  double toz3 = 2.*oz2/zeta;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    double ezr = exp(-zeta*r);
    double tz3r = toz3/r;
    val[i] *= -ezr*(oz2+tz3r)+tz3r;
  }

  return;
}

void eval_inr_d(int gs, double* grid, double* val, int n1, int l1, double zeta1)
{
 #if EVL64
  #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
  for (int i=0;i<gs;i++)
  {
    double r = grid[6*i+3];
    double v = dvinr_gam(n1,l1,r,zeta1);
    val[3*i+0] *= v*grid[6*i+0]/r;
    val[3*i+1] *= v*grid[6*i+1]/r;
    val[3*i+2] *= v*grid[6*i+2]/r;
  }
  return;
 #else
  #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float v = dvinr_gamf(n1,l1,r,zeta1);
    val[3*i+0] *= v*grid[6*i+0]/r;
    val[3*i+1] *= v*grid[6*i+1]/r;
    val[3*i+2] *= v*grid[6*i+2]/r;
  }
  return;
 #endif
}


void eval_inr_d(int gs, float* grid, float* val, int n1, int l1, float zeta1)
{
#if USE_GAM
 #if EVL64
  #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    double v = dvinr_gam(n1,l1,r,zeta1);
    val[3*i+0] *= v*grid[6*i+0]/r;
    val[3*i+1] *= v*grid[6*i+1]/r;
    val[3*i+2] *= v*grid[6*i+2]/r;
  }
  return;
 #else
  #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float v = dvinr_gamf(n1,l1,r,zeta1);
    val[3*i+0] *= v*grid[6*i+0]/r;
    val[3*i+1] *= v*grid[6*i+1]/r;
    val[3*i+2] *= v*grid[6*i+2]/r;
  }
  return;
 #endif
#endif

  if (n1==1)
    eval_inr_1s_r1d(gs,grid,val,zeta1);
  else if (n1==2)
  {
    if (l1==0)
      eval_inr_2s_r1d(gs,grid,val,zeta1);
    else if (l1==1)
      eval_inr_2p_r1d(gs,grid,val,zeta1);
  }
  else if (n1==3)
  {
    if (l1==0)
      eval_inr_3s_r1d(gs,grid,val,zeta1);
    else if (l1==1)
      eval_inr_3p_r1d(gs,grid,val,zeta1);
    else if (l1==2)
      eval_inr_3d_r1d(gs,grid,val,zeta1);
  }
  else if (n1==4)
  {
    if (l1==0)
      eval_inr_4s_r1d(gs,grid,val,zeta1);
    else if (l1==1)
      eval_inr_4p_r1d(gs,grid,val,zeta1);
    else if (l1==2)
      eval_inr_4d_r1d(gs,grid,val,zeta1);
    else if (l1==3)
      eval_inr_4f_r1d(gs,grid,val,zeta1);
  }
  else if (n1==5)
  {
    if (l1==0)
      eval_inr_5s_r1d(gs,grid,val,zeta1);
    else if (l1==1)
      eval_inr_5p_r1d(gs,grid,val,zeta1);
    else if (l1==2)
      eval_inr_5d_r1d(gs,grid,val,zeta1);
    else if (l1==3)
      eval_inr_5f_r1d(gs,grid,val,zeta1);
    else if (l1==4)
      eval_inr_5g_r1d(gs,grid,val,zeta1);
  }
  else if (n1==6)
  {
    if (l1==0)
      eval_inr_6s_r1d(gs,grid,val,zeta1);
    else if (l1==1)
      eval_inr_6p_r1d(gs,grid,val,zeta1);
    else if (l1==2)
      eval_inr_6d_r1d(gs,grid,val,zeta1);
    else if (l1==3)
      eval_inr_6f_r1d(gs,grid,val,zeta1);
    //else if (l1==4)
    //
    else if (l1==5)
      eval_inr_6h_r1d(gs,grid,val,zeta1);
  }
  else if (n1==7)
  {
    if (l1==0)
      eval_inr_7s_r1d(gs,grid,val,zeta1);
    else if (l1==1)
      eval_inr_7p_r1d(gs,grid,val,zeta1);
    else if (l1==2)
      eval_inr_7d_r1d(gs,grid,val,zeta1);
    //else if (l1==3)
    //else if (l1==4)
    //else if (l1==5)
  }
  else if (n1==7)
  {
    if (l1==0)
      eval_inr_8s_r1d(gs,grid,val,zeta1);
  }
}

void eval_inr_r12(int gs, double* grid, double* val, int n1, int l1, double zeta1, int index)
{
 //double precision evaluation, gamma ftn only

  //printf("  eval_inr_r12 for n/l/zt: %i %i %8.5f \n",n1,l1,zeta1);
  #if USE_ACC
  #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
  #endif
  for (int i=0;i<gs;i++)
  {
    double v1d = vinr_gam(n1,l1,grid[6*i+3],zeta1);
    val[i] *= v1d;
  }
}

void eval_inr_r12(int gs, double* grid, double* val, int n1, int l1, double zeta1, int index, int tid)
{
 //double precision evaluation, gamma ftn only

  //printf("  eval_inr_r12 for n/l/zt: %i %i %8.5f \n",n1,l1,zeta1);
  #if USE_ACC
  #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) async(tid+1)
  #endif
  for (int i=0;i<gs;i++)
  {
    double v1d = vinr_gam(n1,l1,grid[6*i+3],zeta1);
    val[i] *= v1d; 
  }
  #pragma acc wait(tid+1)
}


void eval_inr_r12(int gs, float* grid, float* val, int n1, int l1, float zeta1, int index)
{
#if USE_GAM
 #if EVL64
  #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
  for (int i=0;i<gs;i++)
  {
    double v1d = vinr_gam(n1,l1,grid[6*i+3],zeta1);
    val[i] *= v1d;
  }
 #else
  #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
  for (int i=0;i<gs;i++)
  {
    float v1 = vinr_gamf(n1,l1,grid[6*i+3],zeta1);
    val[i] *= v1;
  }
 #endif
  return;
#endif

 #if !USE_GAM
  if (index==3)
  {
    if (n1==1)
      eval_inr_1s_r1(gs,grid,val,zeta1);
    else if (n1==2)
    {
      if (l1==0)
        eval_inr_2s_r1(gs,grid,val,zeta1);
      else if (l1==1)
        eval_inr_2p_r1(gs,grid,val,zeta1);
    }
    else if (n1==3)
    {
      if (l1==0)
        eval_inr_3s_r1(gs,grid,val,zeta1);
      else if (l1==1)
        eval_inr_3p_r1(gs,grid,val,zeta1);
      else if (l1==2)
        eval_inr_3d_r1(gs,grid,val,zeta1);
    }
    else if (n1==4)
    {
      //printf(" inr_d n1: 4 l1: %i \n",l1);
      if (l1==0)
        eval_inr_4s_r1(gs,grid,val,zeta1);
      else if (l1==1)
        eval_inr_4p_r1(gs,grid,val,zeta1);
      else if (l1==2)
        eval_inr_4d_r1(gs,grid,val,zeta1);
      else if (l1==3)
        eval_inr_4f_r1(gs,grid,val,zeta1);
    }
    else if (n1==5)
    {
      if (l1==0)
        eval_inr_5s_r1(gs,grid,val,zeta1);
      else if (l1==1)
        eval_inr_5p_r1(gs,grid,val,zeta1);
      else if (l1==2)
        eval_inr_5d_r1(gs,grid,val,zeta1);
      else if (l1==3)
        eval_inr_5f_r1(gs,grid,val,zeta1);
      else if (l1==4)
        eval_inr_5g_r1(gs,grid,val,zeta1);
    }
    else if (n1==6)
    {
      if (l1==0)
        eval_inr_6s_r1(gs,grid,val,zeta1);
      else if (l1==1)
        eval_inr_6p_r1(gs,grid,val,zeta1);
      else if (l1==2)
        eval_inr_6d_r1(gs,grid,val,zeta1);
      else if (l1==3)
        eval_inr_6f_r1(gs,grid,val,zeta1);
      //else if (l1==4)
      //
      else if (l1==5)
        eval_inr_6h_r1(gs,grid,val,zeta1);
    }
    else if (n1==7)
    {
      if (l1==0)
        eval_inr_7s_r1(gs,grid,val,zeta1);
      else if (l1==1)
        eval_inr_7p_r1(gs,grid,val,zeta1);
      else if (l1==2)
        eval_inr_7d_r1(gs,grid,val,zeta1);
      //else if (l1==3)
      //else if (l1==4)
      //else if (l1==5)
    }
    else if (n1==8)
    {
      if (l1==0)
        eval_inr_8s_r1(gs,grid,val,zeta1);
    }
  } 
  if (index==4)
  {
    printf(" ERROR: shouldn't be here in index==4 \n");
    exit(1);
  }
 #endif

  return;
}

void get_inr_2s_f(int nrad, float zeta, float* r, float* inr)
{
  float oz2 = 1./zeta/zeta;
  float oz3 = oz2/zeta;

#if USE_ACC
 #pragma acc parallel loop independent present(r[0:nrad],inr[0:nrad])
#endif
  for (int i=0;i<nrad;i++)
  {
    float r1 = r[i];
    float zr = zeta*r1;
    float roz2 = r1*oz2;
    float soz4r = 6.f*oz2*oz2/r1;
    float ezr = expf(-zeta*r1);

    float t1 = -ezr*(4.f*oz3+roz2+soz4r);
    float t2 = soz4r;
    inr[i] = t1+t2;
  }

  return;
}

void get_inr_2p_f(int nrad, float zeta, float* r, float* inr)
{
  float oz = 1./zeta;
  float oz2 = oz*oz;
  float oz3 = oz2*oz;
  float oz4 = oz2*oz2;
  float oz5 = oz3*oz2;

#if USE_ACC
 #pragma acc parallel loop independent present(r[0:nrad],inr[0:nrad])
#endif
  for (int i=0;i<nrad;i++)
  {
    float r1 = r[i];
    float zr = zeta*r1;
    float r2 = r1*r1;
    float ezr = expf(-zeta*r1);
    float toz5r2 = 24.*oz5/r2;

    float t1 = -ezr*(toz5r2+24.*oz4/r1+12.*oz3+3.*r1*oz2);
    float t2 = toz5r2;
    inr[i] = t1+t2;
  }

  return;
}

void get_inr_1s_f(int nrad, float zeta, float* r, float* inr)
{
  float oz2 = 1./zeta/zeta;
  float oz3 = oz2/zeta;

#if USE_ACC
 #pragma acc parallel loop independent present(r[0:nrad],inr[0:nrad])
#endif
  for (int i=0;i<nrad;i++)
  {
    float r1 = r[i];
    float zr = zeta*r1;
    float ezr = expf(-zeta*r1);
    float tz3r = 2.*oz3/r1;
  
    float t1 = -ezr*(oz2 + tz3r);
    float t2 = tz3r;

    inr[i] = t1+t2;
  }

  return;
}

void get_inr_1s(int nrad, double zeta, double* r, double* inr)
{
  double oz2 = 1./zeta/zeta;
  double oz3 = oz2/zeta;

  for (int i=0;i<nrad;i++)
  {
    double r1 = r[i];
    double zr = zeta*r1;
    double ezr = exp(-zeta*r1);
    double tz3r = 2.*oz3/r1;
  
    double t1 = -ezr*(oz2 + tz3r);
    double t2 = tz3r;

    inr[i] = t1+t2;
  }

  return;
}

//inr over radial grid only
void get_inr(int n1, int l1, double zeta1, int nrad, float* r, float* inr)
{
  if (n1==1 && l1==0)
    get_inr_1s_f(nrad,zeta1,r,inr);
  else if (n1==2)
  {
    if (l1==0)
      get_inr_2s_f(nrad,zeta1,r,inr);
    else if (l1==1)
      get_inr_2p_f(nrad,zeta1,r,inr);
  }
  else if (n1==3)
  {
    //if (l1==0)
    //  get_inr_3d_f(nrad,zeta1,r,inr);
  }   
  else
  {
    //printf(" INR not available for n1>2 \n"); exit(1);
  }
  return;
}

//#pragma acc routine seq
size_t mchn(int n, int k)
{
  if (k > n) return 0;
  if (k * 2 > n) k = n-k;
  if (k == 0) return 1;

  size_t i1 = n;
  for(int i=2;i<=k;++i)
  {
    i1 *= (n-i+1);
    i1 /= i;
  }
  return i1;
}

void eval_inr_yukawa_mt(int gs, double* grid, double* val, int n1, int l1, double zt1, double gam1)
{
  //Mathematica-derived expressions
  double zt = zt1;
  double g = gam1;
  int gs6 = 6*gs;

  if (n1==1)
  {
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
     //included a 0.5
      double r = grid[6*j+3];
      val[j] = 0.5*(4*exp(2*r*zt)*zt + 2*exp(r*(g + zt))*(pow(g,2)*r - zt*(2 + r*zt)))/
              (exp(r*(g + 2*zt))*r*pow(pow(g,2) - pow(zt,2),2));
    }
  }
  else if (n1==2 && l1==0)
  {
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
     //included a 0.5
      double r = grid[6*j+3];
      val[j] = 0.5*(-4*exp(2*r*zt)*(pow(g,2) + 3*pow(zt,2)) + 2*exp(r*(g + zt))*(pow(g,4)*pow(r,2) - 2*pow(g,2)*(-1 + r*zt*(2 + r*zt)) + pow(zt,2)*(6 + r*zt*(4 + r*zt))))/
               (exp(r*(g + 2*zt))*r*pow(pow(g,2) - pow(zt,2),3));
    }
  }
  else if (n1==2 && l1==1)
  {
   #if 0
    double zt2 = zt1*zt1;
    double zt3 = zt1*zt2;
    double zt4 = zt2*zt2;
    double g2 = gam1*gam1;
    double gpz = gam1+zt1;
    double gmz = gam1-zt1;
    double gpz3 = pow(gpz,3);
    double gmz3 = pow(gmz,3);
   #endif

   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
     //included a 0.75
      double r = grid[6*j+3];
      val[j] = 0.75*(-32*exp(r*zt)*(1 + g*r)*zt + 4*exp(g*r)* (pow(g,4)*pow(r,3) - 2*pow(g,2)*pow(r,2)*zt*(2 + r*zt) + zt*(2 + r*zt)*(4 + r*zt*(2 + r*zt))))/
      (exp(r*(g + zt))*pow(r,2)*pow(pow(g,2) - pow(zt,2),3));
    }
  }
  else if (n1==3 && l1==0)
  {
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
     //included 0.5
      double r = grid[6*j+3];
      val[j] = 0.5*(48*exp(2*r*zt)*zt*(pow(g,2) + pow(zt,2)) + 
     2*exp(r*(g + zt))*(pow(g,6)*pow(r,3) - 3*pow(g,4)*r*(-2 + r*zt*(2 + r*zt)) + 
        3*pow(g,2)*zt*(-8 + r*zt*pow(2 + r*zt,2)) - 
        pow(zt,3)*(24 + r*zt*(18 + r*zt*(6 + r*zt)))))/
   (exp(r*(g + 2*zt))*r*pow(pow(g,2) - pow(zt,2),4));
    }
  }
  else if (n1==3 && l1==1)
  {
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
      double r = grid[6*j+3];
      val[j] = 0.75*(2*(((1 + g*r)*(16*exp(r*(g + zt))*pow(g,5) + 
            80*exp(r*(g + zt))*pow(g,3)*pow(zt,2) + 
            exp(2*g*r)*pow(g + zt,4)*
             (2*zt + (-2 + g*r - r*zt)*
                (g*(4 + g*r*(-2 + g*r)) + g*r*(3 - 2*g*r)*zt + r*(-1 + g*r)*pow(zt,2))) - 
            pow(g - zt,4)*(2*zt + (g*(4 + g*r*(2 + g*r)) + g*r*(3 + 2*g*r)*zt + 
                  r*(1 + g*r)*pow(zt,2))*(2 + r*(g + zt)))))/exp(r*(2*g + zt)) + 
       (2*pow(g - zt,4)*(2*zt + (g*(4 + g*r*(2 + g*r)) + g*r*(3 + 2*g*r)*zt + 
               r*(1 + g*r)*pow(zt,2))*(2 + r*(g + zt)))*(g*r*cosh(g*r) - sinh(g*r)))/
        exp(r*(g + zt))))/(pow(g,3)*pow(r,2)*pow(g - zt,4)*pow(g + zt,4));
    }
  }
  else if (n1==3 && l1==2)
  {
    double fs = 5./6.;
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
      double r = grid[6*j+3];
      val[j] = fs*(288*exp(r*zt)*(3 + g*r*(3 + g*r))*zt + 
     6*exp(g*r)*(pow(g,6)*pow(r,5) - 3*pow(g,4)*pow(r,4)*zt*(2 + r*zt) + 
        3*pow(g,2)*pow(r,2)*zt*(2 + r*zt)*(4 + r*zt*(2 + r*zt)) - 
        zt*(144 + r*zt*(12 + pow(r,2)*pow(zt,2))*(12 + r*zt*(6 + r*zt)))))/
   (exp(r*(g + zt))*pow(r,3)*pow(pow(g,2) - pow(zt,2),4));
    }
  }
  else if (n1==4 && l1==0)
  {
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
     //included 0.5
      double r = grid[6*j+3];
      val[j] = 0.5*(-48*exp(2*r*zt)*(pow(g,4) + 10*pow(g,2)*pow(zt,2) + 5*pow(zt,4)) + 
     2*exp(r*(g + zt))*(pow(g,8)*pow(r,4) - 
        4*pow(g,6)*pow(r,2)*(-1 + r*zt)*(3 + r*zt) - 
        4*pow(g,2)*pow(zt,2)*(-60 + pow(r,2)*pow(zt,2)*(15 + r*zt*(6 + r*zt))) + 
        6*pow(g,4)*(4 + r*zt*(-16 + r*zt*(2 + r*zt*(4 + r*zt)))) + 
        pow(zt,4)*(120 + r*zt*(96 + r*zt*(36 + r*zt*(8 + r*zt))))))/
   (exp(r*(g + 2*zt))*r*pow(pow(g,2) - pow(zt,2),5));
    }
  }
  else if (n1==4 && l1==1)
  {
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
     //included 0.75
      double r = grid[6*j+3];
      val[j] = 0.75*(4*(-48*exp(r*zt)*(1 + g*r)*zt*(3*pow(g,2) + 5*pow(zt,2)) + 
       exp(g*r)*(pow(g,8)*pow(r,5) - 
          2*pow(g,6)*pow(r,3)*(-5 + 2*r*zt*(2 + r*zt)) + 
          6*pow(g,4)*pow(r,2)*zt*(-12 + r*zt*(1 + r*zt)*(3 + r*zt)) - 
          2*pow(g,2)*zt*(-72 + r*zt*(-72 + r*zt*(24 + r*zt*(33 + 2*r*zt*(6 + r*zt))))) + 
          pow(zt,3)*(240 + r*zt*(240 + r*zt*(120 + r*zt*(38 + r*zt*(8 + r*zt))))))))/
   (exp(r*(g + zt))*pow(r,2)*pow(pow(g,2) - pow(zt,2),5));
    }
  }
  else if (n1==4 && l1==2)
  {
    double fs = 5./6.;
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
      double r = grid[6*j+3];
      val[j] = fs*(6*(-48*exp(r*zt)*(3 + g*r*(3 + g*r))*(pow(g,2) + 7*pow(zt,2)) + 
       exp(g*r)*(pow(g,8)*pow(r,6) - 
          2*pow(g,6)*pow(r,4)*(-3 + 2*r*zt*(2 + r*zt)) + 
          6*pow(g,4)*pow(r,2)*(-4 + r*zt*(-4 + r*zt*(5 + r*zt*(4 + r*zt)))) - 
          2*pow(g,2)*(-72 + r*zt*(-72 + 
                r*zt*(48 + r*zt*(72 + r*zt*(39 + 2*r*zt*(6 + r*zt)))))) + 
          pow(zt,2)*(1008 + r*zt*(1008 + 
                r*zt*(504 + r*zt*(168 + r*zt*(42 + r*zt*(8 + r*zt)))))))))/
   (exp(r*(g + zt))*pow(r,3)*pow(g - zt,5)*pow(g + zt,5));
    }
  }
  else if (n1==4 && l1==3)
  {
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
     //included 0.875
      double r = grid[6*j+3];
      val[j] = 0.875*(8*(-384*exp(r*zt)*(15 + g*r*(15 + g*r*(6 + g*r)))*zt + 
       exp(g*r)*(pow(g,8)*pow(r,7) - 4*pow(g,6)*pow(r,6)*zt*(2 + r*zt) + 
          6*pow(g,4)*pow(r,4)*zt*(2 + r*zt)*(4 + r*zt*(2 + r*zt)) - 
          4*pow(g,2)*pow(r,2)*zt*
           (144 + r*zt*(12 + pow(r,2)*pow(zt,2))*(12 + r*zt*(6 + r*zt))) + 
          zt*(5760 + r*zt*(5760 + r*zt*
                 (2880 + r*zt*(960 + r*zt*(240 + r*zt*(48 + r*zt*(8 + r*zt))))))))))/
      (exp(r*(g + zt))*pow(r,4)*pow(pow(g,2) - pow(zt,2),5));
    }
  }
  else if (n1==5 && l1==0)
  {
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
     //included 0.5
      double r = grid[6*j+3];
      val[j] = 0.5*(480*exp(2*r*zt)*zt*(3*pow(g,2) + pow(zt,2))*(pow(g,2) + 3*pow(zt,2)) + 
     2*exp(r*(g + zt))*(pow(g,10)*pow(r,5) + 
        10*pow(g,6)*r*(12 - 24*r*zt + 4*pow(r,3)*pow(zt,3) + 
           pow(r,4)*pow(zt,4)) - 5*pow(g,8)*pow(r,3)*(-4 + r*zt*(2 + r*zt)) - 
        10*pow(g,4)*zt*(72 + r*zt*
            (-108 + r*zt*(-24 + r*zt*(12 + r*zt*(6 + r*zt))))) + 
        5*pow(g,2)*pow(zt,3)*(-480 + 
           r*zt*(-120 + r*zt*(48 + r*zt*(32 + r*zt*(8 + r*zt))))) - 
        pow(zt,5)*(720 + r*zt*(600 + r*zt*(240 + r*zt*(60 + r*zt*(10 + r*zt)))))))/
       (exp(r*(g + 2*zt))*r*pow(pow(g,2) - pow(zt,2),6));
    }
  }
  else if (n1==5 && l1==1)
  {
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
     //included 0.75
      double r = grid[6*j+3];
      val[j] = 0.75*(2*((g*(1 + g*r)*(pow(g + zt,7)*(120*r*(g - zt) + 
               exp(r*(g - zt))*r*(-g + zt)*
                (120 - r*(g - zt)*(120 + 
                     r*(g - zt)*(-60 + r*(g - zt)*(20 + r*(g - zt)*(-5 + g*r - r*zt))))
                  )) + r*pow(g - zt,7)*(g + zt)*
             (120 - (120 + r*(g + zt)*
                   (120 + r*(g + zt)*
                      (60 + r*(g + zt)*(20 + r*(g + zt)*(5 + r*(g + zt))))))/
                exp(r*(g + zt)))))/exp(g*r) + 
       (g - zt)*(g + zt)*(((1 + g*r)*
             (pow(g + zt,6)*(24*r*(g - zt) + 
                  exp(r*(g - zt))*r*(-g + zt)*
                   (24 + r*(g - zt)*
                      (-24 + r*(g - zt)*(12 + r*(g - zt)*(-4 + g*r - r*zt))))) + 
               r*pow(g - zt,6)*(g + zt)*
                (24 - (24 + r*(g + zt)*
                      (24 + r*(g + zt)*(12 + r*(g + zt)*(4 + r*(g + zt)))))/
                   exp(r*(g + zt)))))/exp(g*r) + 
          (2*r*pow(g - zt,6)*(g + zt)*
             (24 + r*(g + zt)*(24 + r*(g + zt)*(12 + r*(g + zt)*(4 + r*(g + zt)))))*
             (g*r*cosh(g*r) - sinh(g*r)))/exp(r*(g + zt))) + 
       (2*g*r*pow(g - zt,7)*(g + zt)*
          (120 + r*(g + zt)*(120 + 
               r*(g + zt)*(60 + r*(g + zt)*(20 + r*(g + zt)*(5 + r*(g + zt))))))*
          (g*r*cosh(g*r) - sinh(g*r)))/exp(r*(g + zt))))/
   (pow(g,3)*pow(r,3)*pow(g - zt,7)*pow(g + zt,7));
    }
  }
  else if (n1==5 && l1==2)
  {
    double fs = 5./6.;
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
     //included fs
      double r = grid[6*j+3];
      val[j] = fs*(6*(384*exp(r*zt)*(3 + g*r*(3 + g*r))*zt*(3*pow(g,2) + 7*pow(zt,2)) + 
       exp(g*r)*(pow(g,10)*pow(r,7) + 
          pow(g,8)*pow(r,5)*(14 - 5*r*zt*(2 + r*zt)) + 
          2*pow(g,6)*pow(r,4)*zt*(-72 + r*zt*(12 + 5*r*zt*(4 + r*zt))) - 
          2*pow(g,4)*pow(r,2)*zt*
           (-288 + r*zt*(-288 + r*zt*(24 + r*zt*(78 + 5*r*zt*(6 + r*zt))))) + 
          pow(g,2)*zt*(-3456 + r*zt*
              (-3456 + r*zt*(-384 + 
                   r*zt*(768 + r*zt*(528 + r*zt*(184 + 5*r*zt*(8 + r*zt))))))) - 
          pow(zt,3)*(8064 + r*zt*
              (8064 + r*zt*(4032 + 
                   r*zt*(1344 + r*zt*(336 + r*zt*(66 + r*zt*(10 + r*zt))))))))))/
         (exp(r*(g + zt))*pow(r,3)*pow(g - zt,6)*pow(g + zt,6));
    }
  }
  else if (n1==5 && l1==3)
  {
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
     //included 0.??
      double r = grid[6*j+3];
      val[j] = 0.875*(8*(384*exp(r*zt)*(15 + g*r*(15 + g*r*(6 + g*r)))*(pow(g,2) + 9*pow(zt,2)) + 
       exp(g*r)*(pow(g,10)*pow(r,8) + 
          pow(g,8)*pow(r,6)*(8 - 5*r*zt*(2 + r*zt)) + 
          2*pow(g,6)*pow(r,4)*(-24 + r*zt*(-24 + r*zt*(24 + 5*r*zt*(4 + r*zt)))) - 
          2*pow(g,4)*pow(r,2)*(-288 + 
             r*zt*(-288 + r*zt*(72 + r*zt*(168 + r*zt*(96 + 5*r*zt*(6 + r*zt)))))) + 
          pow(g,2)*(-5760 + r*zt*
              (-5760 + r*zt*(2304 + 
                   r*zt*(4224 + r*zt*
                       (2352 + r*zt*(816 + r*zt*(208 + 5*r*zt*(8 + r*zt)))))))) - 
          pow(zt,2)*(51840 + r*zt*
              (51840 + r*zt*(25920 + 
                   r*zt*(8640 + r*zt*
                       (2160 + r*zt*(432 + r*zt*(72 + r*zt*(10 + r*zt)))))))))))/
         (exp(r*(g + zt))*pow(r,4)*pow(g - zt,6)*pow(g + zt,6));
    }
  }
  else if (n1==5 && l1==4)
  {
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
     //included 0.9
      double r = grid[6*j+3];
      val[j] = 0.9*(10*(3840*exp(r*zt)*(105 + g*r*(105 + g*r*(45 + g*r*(10 + g*r))))*zt + 
       exp(g*r)*(pow(g,10)*pow(r,9) - 5*pow(g,8)*pow(r,8)*zt*(2 + r*zt) + 
          10*pow(g,6)*pow(r,6)*zt*(2 + r*zt)*(4 + r*zt*(2 + r*zt)) - 
          10*pow(g,4)*pow(r,4)*zt*
           (144 + r*zt*(12 + pow(r,2)*pow(zt,2))*(12 + r*zt*(6 + r*zt))) + 
          5*pow(g,2)*pow(r,2)*zt*
           (5760 + r*zt*(5760 + r*zt*
                 (2880 + r*zt*(960 + r*zt*(240 + r*zt*(48 + r*zt*(8 + r*zt))))))) - 
          zt*(403200 + r*zt*(403200 + 
                r*zt*(201600 + r*zt*(67200 + 
                      r*zt*(16800 + r*zt*(3360 + r*zt*(560 + r*zt*(80 + r*zt*(10 + r*zt)))))
                      )))))))/
         (exp(r*(g + zt))*pow(r,5)*pow(pow(g,2) - pow(zt,2),6));
    }
  }
  else if (n1==6 && l1==0)
  {
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
     //included 0.5
      double r = grid[6*j+3];
      val[j] = 0.5*(-1440*exp(2*r*zt)*(pow(g,6) + 21*pow(g,4)*pow(zt,2) + 
        35*pow(g,2)*pow(zt,4) + 7*pow(zt,6)) + 
     2*exp(r*(g + zt))*(pow(g,12)*pow(r,6) - 
        6*pow(g,10)*pow(r,4)*(-5 + r*zt*(2 + r*zt)) + 
        15*pow(g,8)*pow(r,2)*(-2 + r*zt)*(-12 + r*zt*(10 + r*zt*(6 + r*zt))) - 
        20*pow(g,6)*(-36 + r*zt*(216 + 
              r*zt*(3 + r*zt)*(-48 + pow(r,2)*pow(zt,2)*(3 + r*zt)))) + 
        15*pow(g,4)*pow(zt,2)*(1008 + 
           r*zt*(-672 - 336*r*zt + pow(r,3)*pow(zt,3)*(28 + r*zt*(8 + r*zt)))) - 
        6*pow(g,2)*pow(zt,4)*(-4200 - 1680*r*zt + 
           pow(r,3)*pow(zt,3)*(160 + r*zt*(55 + r*zt*(10 + r*zt)))) + 
        pow(zt,6)*(5040 + r*zt*(4320 + 
              r*zt*(1800 + r*zt*(480 + r*zt*(90 + r*zt*(12 + r*zt))))))))/
       (exp(r*(g + 2*zt))*r*pow(pow(g,2) - pow(zt,2),7));
    }
  }
  else if (n1==6 && l1==1)
  {
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
     //included 0.??
      double r = grid[6*j+3];
      val[j] = 0.75*(-2*((g*(1 + g*r)*(pow(g + zt,8)*
             (720*r*(g - zt) + exp(r*(g - zt))*r*(-g + zt)*
                (720 + r*(g - zt)*(-720 + 
                     r*(g - zt)*(360 + 
                        r*(g - zt)*
                         (-120 + r*(g - zt)*(30 + r*(g - zt)*(-6 + g*r - r*zt))))))) - 
            r*pow(g - zt,8)*(g + zt)*
             (720 - (720 + r*(g + zt)*
                   (720 + r*(g + zt)*
                      (360 + r*(g + zt)*
                         (120 + r*(g + zt)*(30 + r*(g + zt)*(6 + r*(g + zt)))))))/
                exp(r*(g + zt)))))/exp(g*r) + 
       (g - zt)*(g + zt)*(((1 + g*r)*
             (pow(g + zt,7)*(120*r*(g - zt) + 
                  exp(r*(g - zt))*r*(-g + zt)*
                   (120 - r*(g - zt)*
                      (120 + r*(g - zt)*
                         (-60 + r*(g - zt)*(20 + r*(g - zt)*(-5 + g*r - r*zt)))))) - 
               r*pow(g - zt,7)*(g + zt)*
                (120 - (120 + r*(g + zt)*
                      (120 + r*(g + zt)*
                         (60 + r*(g + zt)*(20 + r*(g + zt)*(5 + r*(g + zt))))))/
                   exp(r*(g + zt)))))/exp(g*r) - 
          (2*r*pow(g - zt,7)*(g + zt)*
             (120 + r*(g + zt)*(120 + 
                  r*(g + zt)*(60 + r*(g + zt)*(20 + r*(g + zt)*(5 + r*(g + zt))))))*
             (g*r*cosh(g*r) - sinh(g*r)))/exp(r*(g + zt))) - 
       (2*g*r*pow(g - zt,8)*(g + zt)*
          (720 + r*(g + zt)*(720 + 
               r*(g + zt)*(360 + r*(g + zt)*
                   (120 + r*(g + zt)*(30 + r*(g + zt)*(6 + r*(g + zt)))))))*
          (g*r*cosh(g*r) - sinh(g*r)))/exp(r*(g + zt))))/
       (pow(g,3)*pow(r,3)*pow(g - zt,8)*pow(g + zt,8));
    }
  }
  else if (n1==6 && l1==2)
  {
    double fs = 5./6.;
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
     //included fs
      double r = grid[6*j+3];
      val[j] = fs*(-3*((pow(g,2)*(3 + g*r*(3 + g*r))*
          (pow(g + zt,8)*(720*r*(g - zt) + 
               exp(r*(g - zt))*r*(-g + zt)*
                (720 + r*(g - zt)*(-720 + 
                     r*(g - zt)*(360 + 
                        r*(g - zt)*
                         (-120 + r*(g - zt)*(30 + r*(g - zt)*(-6 + g*r - r*zt))))))) + 
            r*pow(g - zt,8)*(g + zt)*
             (720 - (720 + r*(g + zt)*
                   (720 + r*(g + zt)*
                      (360 + r*(g + zt)*
                         (120 + r*(g + zt)*(30 + r*(g + zt)*(6 + r*(g + zt)))))))/
                exp(r*(g + zt)))))/exp(g*r) + 
       (2*pow(g,2)*r*pow(g - zt,8)*(g + zt)*
          (720 + r*(g + zt)*(720 + 
               r*(g + zt)*(360 + r*(g + zt)*
                   (120 + r*(g + zt)*(30 + r*(g + zt)*(6 + r*(g + zt)))))))*
          (3*g*r*cosh(g*r) - (3 + pow(g,2)*pow(r,2))*sinh(g*r)))/
        exp(r*(g + zt)) + 3*g*(g - zt)*(g + zt)*
        (((3 + g*r*(3 + g*r))*(pow(g + zt,7)*
                (120*r*(g - zt) + exp(r*(g - zt))*r*(-g + zt)*
                   (120 - r*(g - zt)*
                      (120 + r*(g - zt)*
                         (-60 + r*(g - zt)*(20 + r*(g - zt)*(-5 + g*r - r*zt)))))) + 
               r*pow(g - zt,7)*(g + zt)*
                (120 - (120 + r*(g + zt)*
                      (120 + r*(g + zt)*
                         (60 + r*(g + zt)*(20 + r*(g + zt)*(5 + r*(g + zt))))))/
                   exp(r*(g + zt)))))/exp(g*r) + 
          (2*r*pow(g - zt,7)*(g + zt)*
             (120 + r*(g + zt)*(120 + 
                  r*(g + zt)*(60 + r*(g + zt)*(20 + r*(g + zt)*(5 + r*(g + zt))))))*
             (3*g*r*cosh(g*r) - (3 + pow(g,2)*pow(r,2))*sinh(g*r)))/
           exp(r*(g + zt))) + 3*pow(g - zt,2)*pow(g + zt,2)*
        (((3 + g*r*(3 + g*r))*(pow(g + zt,6)*
                (24*r*(g - zt) + exp(r*(g - zt))*r*(-g + zt)*
                   (24 + r*(g - zt)*
                      (-24 + r*(g - zt)*(12 + r*(g - zt)*(-4 + g*r - r*zt))))) + 
               r*pow(g - zt,6)*(g + zt)*
                (24 - (24 + r*(g + zt)*
                      (24 + r*(g + zt)*(12 + r*(g + zt)*(4 + r*(g + zt)))))/
                   exp(r*(g + zt)))))/exp(g*r) - 
          (2*r*pow(g - zt,6)*(g + zt)*
             (24 + r*(g + zt)*(24 + r*(g + zt)*(12 + r*(g + zt)*(4 + r*(g + zt)))))*
             (-3*g*r*cosh(g*r) + (3 + pow(g,2)*pow(r,2))*sinh(g*r)))/
           exp(r*(g + zt)))))/
       (pow(g,5)*pow(r,4)*pow(g - zt,8)*pow(g + zt,8));
    }
  }
  else if (n1==6 && l1==3)
  {
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
     //included 0.875
      double r = grid[6*j+3];
      val[j] = 0.875*(8*(-11520*exp(r*zt)*(15 + g*r*(15 + g*r*(6 + g*r)))*zt*
        (pow(g,2) + 3*pow(zt,2)) + 
       exp(g*r)*(pow(g,12)*pow(r,9) - 
          6*pow(g,10)*pow(r,7)*(-1 + r*zt)*(3 + r*zt) + 
          15*pow(g,8)*pow(r,6)*zt*(-16 + r*zt*(2 + r*zt*(4 + r*zt))) + 
          15*pow(g,4)*pow(r,2)*zt*(2 + r*zt)*
           (-576 + r*zt*(-24 + 12*r*zt + pow(r,3)*pow(zt,3))*
              (12 + r*zt*(6 + r*zt))) - 
          20*pow(g,6)*pow(r,4)*zt*
           (-72 - 72*r*zt + pow(r,3)*pow(zt,3)*(15 + r*zt*(6 + r*zt))) - 
          6*pow(g,2)*zt*(-28800 + 
             r*zt*(-28800 + r*zt*(-5760 + 
                   r*zt*(3840 + r*zt*
                       (3120 + r*zt*(1200 + r*zt*(320 + r*zt*(65 + r*zt*(10 + r*zt)))))
                      )))) + pow(zt,3)*
           (518400 + r*zt*(518400 + 
                r*zt*(259200 + r*zt*
                    (86400 + r*zt*(21600 + 
                         r*zt*(4320 + r*zt*(720 + r*zt*(102 + r*zt*(12 + r*zt))))))))))
       ))/(exp(r*(g + zt))*pow(r,4)*pow(g - zt,7)*pow(g + zt,7));
    }
  }
  else if (n1==6 && l1==4)
  {
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
     //included 0.9
      double r = grid[6*j+3];
      val[j] = 0.9*(10*(-3840*exp(r*zt)*(105 + g*r*(105 + g*r*(45 + g*r*(10 + g*r))))*
        (pow(g,2) + 11*pow(zt,2)) + 
       exp(g*r)*(pow(g,12)*pow(r,10) - 
          2*pow(g,10)*pow(r,8)*(-5 + 3*r*zt*(2 + r*zt)) + 
          5*pow(g,8)*pow(r,6)*(-16 + r*zt*(-16 + r*zt*(14 + 3*r*zt*(4 + r*zt)))) - 
          20*pow(g,6)*pow(r,4)*
           (-72 + r*zt*(-72 + r*zt*(8 + r*zt*(32 + r*zt*(19 + r*zt*(6 + r*zt)))))) + 
          5*pow(g,4)*pow(r,2)*(-5760 + 
             r*zt*(-5760 + r*zt*(288 + 
                   r*zt*(2208 + r*zt*
                       (1344 + r*zt*(480 + r*zt*(124 + 3*r*zt*(8 + r*zt)))))))) - 
          2*pow(g,2)*(-201600 + r*zt*
              (-201600 + r*zt*(57600 + 
                   r*zt*(124800 + r*zt*
                       (70800 + r*zt*
                          (24720 + 
                            r*zt*(6320 + r*zt*(1280 + r*zt*(215 + 3*r*zt*(10 + r*zt))))
                            )))))) + 
          pow(zt,2)*(4435200 + r*zt*
              (4435200 + r*zt*(2217600 + 
                   r*zt*(739200 + r*zt*
                       (184800 + r*zt*
                          (36960 + 
                            r*zt*(6160 + r*zt*(880 + r*zt*(110 + r*zt*(12 + r*zt)))))))
                   ))))))/
         (exp(r*(g + zt))*pow(r,5)*pow(g - zt,7)*pow(g + zt,7));
    }
  }
  else if (n1==6 && l1==5)
  {
    double et = 11./12.;
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
     //included 11/12
      double r = grid[6*j+3];
      val[j] = et*(12*(-46080*exp(r*zt)*(945 + g*r*(945 + g*r*(420 + g*r*(105 + g*r*(15 + g*r)))))*
        zt + exp(g*r)*(pow(g,12)*pow(r,11) - 
          6*pow(g,10)*pow(r,10)*zt*(2 + r*zt) + 
          15*pow(g,8)*pow(r,8)*zt*(2 + r*zt)*(4 + r*zt*(2 + r*zt)) - 
          20*pow(g,6)*pow(r,6)*zt*
           (144 + r*zt*(12 + pow(r,2)*pow(zt,2))*(12 + r*zt*(6 + r*zt))) + 
          15*pow(g,4)*pow(r,4)*zt*
           (5760 + r*zt*(5760 + r*zt*
                 (2880 + r*zt*(960 + r*zt*(240 + r*zt*(48 + r*zt*(8 + r*zt))))))) - 
          6*pow(g,2)*pow(r,2)*zt*
           (403200 + r*zt*(403200 + 
                r*zt*(201600 + r*zt*
                    (67200 + r*zt*(16800 + 
                         r*zt*(3360 + r*zt*(560 + r*zt*(80 + r*zt*(10 + r*zt)))))))))\
           + zt*(43545600 + r*zt*(43545600 + 
                r*zt*(21772800 + r*zt*
                    (7257600 + r*zt*
                       (1814400 + r*zt*
                          (362880 + 
                            r*zt*(60480 + 
                               r*zt*
                               (8640 + r*zt*(1080 + r*zt*(120 + r*zt*(12 + r*zt))))))))
                   ))))))/
         (exp(r*(g + zt))*pow(r,6)*pow(pow(g,2) - pow(zt,2),7));
    }
  }
  else
  {
    printf("\n ERROR: Yukawa potential for nl: %i %i not available \n",n1,l1);
  }

 //check for incorrect values
 #pragma acc parallel loop present(val[0:gs])
  for (int j=0;j<gs;j++)
  {
    if (val[j]<0.)
      val[j] = 0.;
  }
 #pragma acc parallel loop present(val[0:gs])
  for (int j=0;j<gs;j++)
  {
    if (val[j]>1000.)
      val[j] = 0.;
  }

  return;
}

size_t mchn_v2(int l, int m)
{
  size_t v = 1;
  for (int i=2;i<l+m;i++)
    v *= i;
  for (int i=1;i<m;i++)
    v /= i;
  for (int i=1;i<l-m;i++)
    v /= i;
  return v;
}

void eval_inr_yukawa(int gs, double* grid, double* val, int n1, int l1, double zeta1, double gam1)
{
  if (gam1<0.01) printf(" WARNING: Yukawa potential not accurate for small gamma \n");
  return eval_inr_yukawa_mt(gs,grid,val,n1,l1,zeta1,gam1);

 //this ftn has numerical issues at low r or l>1
  int gs6 = 6*gs;

  double tlp1 = 2.*l1+1.;
  int lm = l1%2;
  double s1 = 1;
  if (lm==0) s1 = -1;
  const double f0 = s1*tlp1;

  bool zmg_zero = 0;
  if (zeta1-gam1<=0.) zmg_zero = 1;

  const double zpg = zeta1+gam1;
  const double zmg = zeta1-gam1;

  printf("  inr_yukawa. gamma: %8.5f  zeta1: %8.5f l1: %i  f0: %8.5f  zmg_zero: %i \n",gam1,zeta1,l1,f0,zmg_zero);

  if (!zmg_zero)
  for (int m=0;m<=l1;m++)
  {
    double lcm = mchn_v2(l1,m);
    double tm1 = pow(2.,-(m+1));
    double gm = pow(gam1,-m);
    double zgpnm = pow(zpg,-(n1-m+1));
    double zgmnm = pow(zmg,-(n1-m+1));
    int nm1 = n1-m+1;
    double s2 = pow(-1.,l1-m+1);

    printf("  m: %i lcm: %8.5f tm1: %8.5f gm: %8.5f zgnm: %8.5f %8.5f s2: %8.5f \n",m,lcm,tm1,gm,zgpnm,zgmnm,s2);

  #if USE_ACC
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
  #endif
    for (int j=0;j<gs;j++)
    {
      double r = grid[6*j+3];
      double zpgr = zpg*r;
      double zmgr = zmg*r;
      double gr = gam1*r;

     //i^l+2 terms accounted for in f0, above
      double jligr = spherical_bessel_1(l1,gr);
      double higr = spherical_bessel_3(l1,gr);

      double Gnmzgr = igamc(nm1,zpgr); //Upper incomplete gamma ftn of (n-m+1,(zeta+gamma)*r)
      double gnmzpgr = igam(nm1,zpgr); //lower incomplete gamma ftn of (n-m+1,(zeta+gamma)*r)
      double gnmzmgr = igam(nm1,zmgr); //lower

      double f1 = lcm*tm1*gm;
      double v1a = 2.*jligr*Gnmzgr*zgpnm;
      double v1b = gnmzpgr*zgpnm;
      double v1c = s2*gnmzmgr*zgmnm;
      double v1 = f1*( v1a + higr*(v1b + v1c) );

      printf("   b's: %8.5f %8.5f   g's:  %8.5f %8.5f %8.5f   v's: %8.5f %8.5f %8.5f %8.5f  \n", jligr,higr,Gnmzgr,gnmzpgr,gnmzmgr,v1a,v1b,v1c,v1);

      val[j] += f0*v1;
    }
  }

  if (zmg_zero)
  for (int m=0;m<=l1;m++)
  {
    double lcm = mchn_v2(l1,m);
    double tm1 = pow(2.,-(m+1));
    double gm = pow(gam1,-m);
    double zgpnm = pow(zpg,-(n1-m+1));
    double zgmnm = pow(zmg,-(n1-m+1));
    int nm1 = n1-m+1;
    double s2 = pow(-1.,l1-m+1);

    printf("  m: %i lcm: %8.5f tm1: %8.5f gm: %8.5f zgnm: %8.5f %8.5f s2: %8.5f \n",m,lcm,tm1,gm,zgpnm,zgmnm,s2);

  #if USE_ACC
   #pragma acc parallel loop present(val[0:gs],grid[0:gs6])
  #endif
    for (int j=0;j<gs;j++)
    {
      double r = grid[6*j+3];
      double zpgr = zpg*r;
      double zmgr = zmg*r;
      double gr = gam1*r;

     //i^l+2 terms accounted for in f0, above
      double jligr = spherical_bessel_1(l1,gr);
      double higr = spherical_bessel_3(l1,gr);

      double Gnmzgr = igamc(nm1,zpgr); //Upper incomplete gamma ftn of (n-m+1,(zeta+gamma)*r)
      double gnmzpgr = igam(nm1,zpgr); //lower incomplete gamma ftn of (n-m+1,(zeta+gamma)*r)
      double gnmzmgr = igamn(nm1,zmgr); //lower incomplete gamma, new version for negative x

      double f1 = lcm*tm1*gm;
      double v1a = 2.*jligr*Gnmzgr*zgpnm;
      double v1b = gnmzpgr*zgpnm;
      double v1c = s2*gnmzmgr*zgmnm;
      double v1 = f1*( v1a + higr*(v1b + v1c) );

      printf(" r: %8.5f   b's: %8.5f %8.5f   g's:  %8.5f %8.5f %8.5f   v's: %8.5f %8.5f %8.5f %8.5f  \n",r,jligr,higr,Gnmzgr,gnmzpgr,gnmzmgr,v1a,v1b,v1c,v1);

      val[j] += f0*v1;
    }
  }

  return;
}

double norm_sh_theta_phi(int l, int m)
{
  double num = (2.*l+1)*fact(l-abs(m));
  double den = 2.*PI*fact(l+abs(m));
  double v1 = sqrt(num/den);
  if (m==0) v1/sqrt(2.);
  return v1;
}

double norm_sh(int l, int m)
{
 #if CART_D
  if (l==2)
  {
    if (m<3) 
      return 1.092548430592080*0.577350269189626;
    else
      return 1.092548430592080;
  }
 #endif

 //Note: rounded the last digit up for s,p,d
 // l==5 are numerical
  if (l==0)
    return 0.2820947917738782;
  if (l==1) //p
    return 0.4886025119029200;
  if (l==2) //d
  {
    if (m==0)
      return 0.3153915652525202;
    if (m==2)
      return 0.5462742152960397;
    return 1.092548430592080;
  }
  if (l==3) //f
  {
    if (m==-3) return 0.5900435899266435;
    if (m==-2) return 2.890611442640554;
    if (m==-1) return 0.4570457994644658;
    if (m== 0) return 0.3731763325901154;
    if (m== 1) return 0.4570457994644658;
    if (m== 2) return 1.445305721320277;
    if (m== 3) return 0.5900435899266435;
  }
  if (l==4) //g
  {
    if (m==-4) return 2.503342941796705;
    if (m==-3 || m==3) return 1.770130769779931;
    if (m==-2) return 0.946174695757560;
    if (m==-1) return 0.669046543557289;
    if (m== 0) return 0.105785546915204;
    if (m== 1) return 0.669046543557289;
    if (m== 2) return 0.473087347878780;
    if (m== 4) return 0.625835735449176;
  }
  if (l==5) //h
  {
    if (m==-5) return 0.656382056840170;
    if (m==-4) return 8.302649259524165;
    if (m==-3) return 0.489238299435250;
    if (m==-2) return 4.793536784973324;
    if (m==-1) return 0.452946651195697;
    if (m== 0) return 0.116950322453424;
    if (m== 1) return 0.452946651195697;
    if (m== 2) return 2.396768392486662;
    if (m== 3) return 0.489238299435250;
    if (m== 4) return 2.075662314881041;
    if (m== 5) return 0.656382056840170;
  }
  if (l==6) //i
  {
    if (m==-6) return 1.366368210383829;
    if (m==-5) return 2.366619162231752;
    if (m==-4) return 2.018259602914897;
    if (m==-3) return 0.921205259514923;
    if (m==-2) return 0.921205259514923;
    if (m==-1) return 0.582621362518731;
    if (m== 0) return 0.063569202267628;
    if (m== 1) return 0.582621362518731;
    if (m== 2) return 0.460602629757462;
    if (m== 3) return 0.921205259514923;
    if (m== 4) return 0.504564900728724;
    if (m== 5) return 2.366619162231752;
    if (m== 6) return 0.683184105191914;
  }
  if (l==7) //j
  {
    if (m==-7) return 0.707162732524596;
    if (m==-6) return 5.2919213236038;
    if (m==-5) return 0.51891557872026;
    if (m==-4) return 4.151324629762083;
    if (m==-3) return 0.156458933862294;
    if (m==-2) return 0.442532692444983;
    if (m==-1) return 0.090331607582517;
    if (m== 0) return 0.068284276912005;
    if (m== 1) return 0.090331607582517;
    if (m== 2) return 0.221266346222491;
    if (m== 3) return 0.156458933862294;
    if (m== 4) return 1.037831157440521;
    if (m== 5) return 0.51891557872026;
    if (m== 6) return 2.6459606618019;
    if (m== 7) return 0.707162732524596;
  }

  printf(" ERROR: norm not implemented for l=%i m=%i \n",l,m);
  exit(1);
  return 1.;
}

double norm_sv(int n, int l, int m, double zeta)
{
  double num = 4.*PI*pow(2.*zeta,n+0.5);
  double den = sqrt(fact(2*n))*(2.*l+1.);
  double val = num/den;
  val *= norm_sh(l,m);
  return val;
}

double norm(int n, int l, int m, double zeta)
{
  double num = pow(2.*zeta,n+0.5);
  double den = sqrt(fact(2*n));

  double val = num/den;
  val *= norm_sh(l,m);
  return val;
}

size_t fact(size_t N) 
{
  if (N==0) return 1;
  return N*fact(N-1);
}

