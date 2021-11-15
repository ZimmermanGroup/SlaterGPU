#include "Vinr.h"

// for getting cutoff
// #include "get_rcutoff.h"
// sorting and search
#if TEST_SORT
#include "lib/thrust_util.h"
#endif
// taylor series approx
// #include "taylor.h"
#include "gamma.h"

#define RGLIMIT 0
#include "fp_def.h"

//this contains norms and Inl(r) functions

//Notes:
//1. in Inr functions, exp could be expf
//2. need to carefully check FP2 vs FP1
//3. need 6h derivative (untested)


void eval_ke(int gs, FP1* grid, FP1* val, int n, int l, FP1 zeta)
{
  FP1 f1 = n*(n-1) - l*(l+1);
  FP1 f2 = 2.f*n*zeta;
  FP1 f3 = zeta*zeta;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r = grid[6*i+3];
    FP1 or1 = 1.f/r;
    FP1 or2 = or1*or1;

    FP1 term = f1*or2 - f2*or1 + f3;
    val[i] *= term;
  }

  return;
}

void eval_ke3(int gs, FP1* grid, FP1* val, int n, int l, FP1 zeta)
{
  FP1 f1 = n*(n-1) - l*(l+1);
  FP1 f2 = 2.f*n*zeta;
  FP1 f3 = zeta*zeta;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r = grid[6*i+3];
    FP1 or1 = 1.f/r;
    FP1 or2 = or1*or1;

    FP1 term = f1*or2 - f2*or1 + f3;
    val[3*i+0] *= term;
    val[3*i+1] *= term;
    val[3*i+2] *= term;
  }

  return;
}

void eval_dke(int gs, FP1* grid, FP1* val, int n, int l, FP1 zeta)
{
  FP2 f1 = n*(n-1) - l*(l+1);
  FP2 f2 = 2.*n*zeta;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 or1 = 1./r;
    FP2 or2 = or1*or1;

    FP2 term = -2.*f1*or2*or2 + f2*or2*or1;
    val[3*i+0] *= x*term;
    val[3*i+1] *= y*term;
    val[3*i+2] *= z*term;
  }

  return;
}

void eval_ne(int gs, FP1* grid, FP1** val, int s1, int s2, int natoms, int* atno, FP1* coords, FP1 A0, FP1 B0, FP1 C0)
{
  int ng = s2-s1;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:ng][0:gs],atno[0:natoms],coords[0:3*natoms])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x1 = grid[6*i]+A0;
    FP1 y1 = grid[6*i+1]+B0;
    FP1 z1 = grid[6*i+2]+C0;

    FP1 ne1 = 0.;
   #pragma acc loop reduction(+:ne1)
    for (int j=0;j<natoms;j++)
    {
      FP1 Zeff = atno[j];
      FP1 x2 = x1-coords[3*j];
      FP1 y2 = y1-coords[3*j+1];
      FP1 z2 = z1-coords[3*j+2];
      FP1 Ra = sqrtf(x2*x2+y2*y2+z2*z2);

      ne1 += Zeff/Ra;
    }
   #pragma acc loop
    for (int k=0;k<ng;k++)
      val[k][i] *= ne1;
  }

  return;
}

void eval_ne_3(int gs, FP1* grid, FP1** val, int s1, int s2, int natoms, int* atno, FP1* coords, FP1 A0, FP1 B0, FP1 C0)
{
  int ng = s2-s1;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:ng][0:3*gs],atno[0:natoms],coords[0:3*natoms])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x1 = grid[6*i]+A0;
    FP1 y1 = grid[6*i+1]+B0;
    FP1 z1 = grid[6*i+2]+C0;

    FP1 ne1 = 0.;
   #pragma acc loop reduction(+:ne1)
    for (int j=0;j<natoms;j++)
    {
      FP1 Zeff = atno[j];
      FP1 x2 = x1-coords[3*j];
      FP1 y2 = y1-coords[3*j+1];
      FP1 z2 = z1-coords[3*j+2];
      FP1 Ra = sqrtf(x2*x2+y2*y2+z2*z2);

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

void eval_d_ne_3(int gs, FP1* grid, FP1** val, int s1, int s2, int natoms, int* atno, FP1* coords, FP1 A0, FP1 B0, FP1 C0)
{
  int ng = s2-s1;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:ng][0:3*gs],atno[0:natoms],coords[0:3*natoms])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x1 = grid[6*i]+A0;
    FP1 y1 = grid[6*i+1]+B0;
    FP1 z1 = grid[6*i+2]+C0;

    FP1 ne1 = 0.;
   #pragma acc loop reduction(+:ne1)
    for (int j=0;j<natoms;j++)
    {
      FP1 Zeff = atno[j];
      FP1 x2 = x1-coords[3*j];
      FP1 y2 = y1-coords[3*j+1];
      FP1 z2 = z1-coords[3*j+2];
      FP1 Ra = sqrtf(x2*x2+y2*y2+z2*z2);

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

void eval_inr_6h_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
 //test this
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz*oz2;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;
  FP2 oz9 = oz4*oz5;
  FP2 oz10 = oz5*oz5;
  FP2 oz11 = oz6*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 r2 = r*r; FP2 r3 = r2*r; FP2 r4 = r2*r2; FP2 r5 = r4*r;
    FP2 ezr = exp(-zeta*r);
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 f1 = 2874009600./r4/r3*oz7*oz6;
    FP2 v1 = ezr*(f1*r*zeta + 1437004800.*oz11/r5 + 479001600.*oz10/r4 + 119750400.*oz9/r3 + 23950080.*oz8/r2 + 3991680.*oz7/r + 570240.*oz6 + 71280.*r*oz5 + 7920.*r2*oz4 + 792.*r3*oz3 + 77.*r4*oz2 + 11.*r5*oz);
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

void eval_inr_5g_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
 //careful, test this
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz*oz2;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;
  FP2 oz9 = oz4*oz5;
  FP2 oz10 = oz5*oz5;
  FP2 oz11 = oz6*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 r2 = r*r; FP2 r3 = r2*r; FP2 r4 = r2*r2;
    FP2 ezr = exp(-zeta*r);
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

   //exp(-x) = 1-x;
    FP2 f1 = 18144000./r3/r3*oz11;
    FP2 v1 = ezr*(f1*r*zeta + 9072000.*oz9/r4 + 3024000.*oz8/r3 + 756000.*oz7/r2 + 151200.*oz6/r + 25200.*oz5 + 3600.*r*oz4 + 450.*r2*oz3 + 54.*r3*oz2 + 9.*r4*oz);
    //FP2 f1 = 18144000.*oz11;
    //FP2 v1 = ezr*f1*r5
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

void eval_inr_6f_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
 //careful, test this
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz*oz2;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;
  FP2 oz9 = oz4*oz5;
  FP2 oz10 = oz5*oz5;
  FP2 oz11 = oz5*oz6;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 r2 = r*r; FP2 r3 = r2*r; FP2 r4 = r2*r2; FP2 r5 = r3*r2;
    FP2 ezr = exp(-zeta*r);
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 f1 = 14515200./r5*oz11;
    FP2 v1 = ezr*(f1*r*zeta + 7257600.*oz9/r3 + 2419200.*oz8/r2 + 604800.*oz7/r + 120960.*oz6 + 20160.*r*oz5 + 2898.*r2*oz4 + 378.*r3*oz3 + 49.*r4*oz2 + 7.*r5*oz);
    v1 += (ezr-1.)*f1;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_5f_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
 //careful, test this
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz*oz2;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;
  FP2 oz9 = oz4*oz5;
  FP2 oz10 = oz5*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 r2 = r*r; FP2 r3 = r2*r; FP2 r4 = r2*r2;
    FP2 ezr = exp(-zeta*r);
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 f1 = 1451520./r2/r3*oz10;
    FP2 v1 = ezr*(f1*r*zeta + 725760.*oz8/r3 + 241920.*oz7/r2 + 60480.*oz6/r + 12096.*oz5 + 2016.*r*oz4 + 294.*r2*oz3 + 42.*r3*oz2 + 7.*r4*oz);
    v1 += (ezr-1.)*f1;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_4f_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz*oz2;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;
  FP2 oz9 = oz4*oz5;

  //#pragma acc update self(grid[0:6*gs],val[0:3*gs])

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 r2 = r*r; FP2 r3 = r2*r;
    FP2 ezr = exp(-zeta*r);
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 f1 = 161280./r2/r3*oz9;
    FP2 v1 = ezr*(f1*r*zeta + 80640.*oz7/r3 + 26880.*oz6/r2 + 6720.*oz5/r + 1344.*oz4 + 224.*r*oz3 + 35.*r2*oz2 + 7.*r3*oz);
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

void eval_inr_7d_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz*oz2;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;
  FP2 oz9 = oz5*oz4;
  FP2 oz11 = oz5*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 r2 = r*r; FP2 r3 = r2*r; FP2 r4 = r2*r2;
    FP2 ezr = exp(-zeta*r);
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 or4z11 = 10886400./r4*oz11;
    FP2 v1 = ezr*(or4z11 + or4z11*r*zeta + 5443200.*oz9/r2 + 1814400.*oz8/r + 453600.*oz7 + 90960.*r*oz6 + 15360.*r2*oz5 + 2280.*oz4*r3 + 310.*r4*oz3 + 40.*r3*r2*oz2 + 5.*r3*r3*oz) - or4z11;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_6d_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz*oz2;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;
  FP2 oz9 = oz5*oz4;
  FP2 oz10 = oz5*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 r2 = r*r; FP2 r3 = r2*r; FP2 r4 = r2*r2;
    FP2 ezr = exp(-zeta*r);
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 or4z10 = 1088640./r4*oz10;
    FP2 v1 = ezr*(or4z10 + or4z10*r*zeta + 544320.*oz8/r2 + 181440.*oz7/r + 45360.*oz6 + 9120.*r*oz5 + 1560.*r2*oz4 + 240.*oz3*r3 + 35.*r4*oz2 + 5.*r3*r2*oz) - or4z10;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_5d_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz*oz2;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;
  //FP2 oz8 = oz4*oz4;
  FP2 oz9 = oz5*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 r2 = r*r; FP2 r3 = r2*r; FP2 r4 = r2*r2;
    FP2 ezr = exp(-zeta*r);
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 or4z9 = 120960./r4*oz9;
    FP2 v1 = ezr*(or4z9 + or4z9*r*zeta + 60480.*oz7/r2 + 20160.*oz6/r + 5040.*oz5 + 1020.*r*oz4 + 180.*r2*oz3 + 30.*oz2*r3 + 5.*r4*oz) - or4z9;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;

    //val[i] *= -ezr*(720.f*oz7r3 + 720.f*oz6r2 + 360.f*oz5r + 120.f*oz4 + 30.f*r*oz3 + 5.f*r2*oz2) + 720.f*oz7r3;
  }

  return;
}

void eval_inr_4d_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz*oz2;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;
  //FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 r2 = r*r; FP2 r3 = r2*r;
    FP2 ezr = exp(-zeta*r);
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 for4z8 = 15120.*oz8/r2/r2;
    FP2 v1 = ezr*(for4z8 + for4z8*r*zeta + 7560.*oz6/r2 + 2520.*oz5/r + 630.*oz4 + 130.*r*oz3 + 25.*r2*oz2 + 5.*oz*r3) - for4z8;

    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;

    //val[i] *= -ezr*(720.f*oz7r3 + 720.f*oz6r2 + 360.f*oz5r + 120.f*oz4 + 30.f*r*oz3 + 5.f*r2*oz2) + 720.f*oz7r3;
  }

  return;
}

void eval_inr_3d_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz*oz2;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz7 = oz3*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 r2 = r*r; FP2 r3 = r2*r;
    FP2 ezr = exp(-zeta*r);
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 tor4z7 = 2160./r2/r2*oz7;
    FP2 v1 = ezr*(tor4z7 + tor4z7*r*zeta + 1080.*oz5/r2 + 360.*oz4/r + 90.*oz3 + 20.*r*oz2 + 5.*r2*oz) - tor4z7;
    //v1 += (ezr-1.)*tor4z7);

    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;

    //val[i] *= -ezr*(720.f*oz7r3 + 720.f*oz6r2 + 360.f*oz5r + 120.f*oz4 + 30.f*r*oz3 + 5.f*r2*oz2) + 720.f*oz7r3;
  }

  return;
}

void eval_inr_7p_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz3*oz2;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz4*oz3;
  FP2 oz8 = oz4*oz4;
  FP2 oz9 = oz4*oz5;
  FP2 oz10 = oz5*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 r2 = r*r; FP2 r3 = r2*r;
    FP2 ezr = exp(-zeta*r);
    FP2 tor3z10 = 725760.*oz10/r3;
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 v1 = ezr*(tor3z10 + tor3z10*r*zeta + 362880.*oz8/r + 121680.*oz7 + 30960.*oz6*r + 6408.*oz5*r2 + 1128.*oz4*r3 + 174.*oz3*r2*r2 + 24.*oz2*r3*r2 + 3.*r3*r3*oz) - tor3z10;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_6p_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz3*oz2;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz4*oz3;
  FP2 oz8 = oz4*oz4;
  FP2 oz9 = oz4*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 r2 = r*r; FP2 r3 = r2*r;
    FP2 ezr = exp(-zeta*r);
    FP2 tor3z9 = 80640.*oz9/r3;
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 v1 = ezr*(tor3z9 + tor3z9*r*zeta + 40320.*oz7/r + 13560.*oz6 + 3480.*oz5*r + 732.*oz4*r2 + 132.*oz3*r3 + 21.*oz2*r2*r2 + 3.*oz*r3*r2) - tor3z9;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_5p_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz3*oz2;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz4*oz3;
  FP2 oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 r2 = r*r; FP2 r3 = r2*r;
    FP2 ezr = exp(-zeta*r);
    FP2 tor3z8 = 10080.*oz8/r3;
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 v1 = ezr*(tor3z8 + tor3z8*r*zeta + 5040.*oz6/r + 1704.*oz5 + 444.*oz4*r + 96.*oz3*r2 + 18.*oz2*r3 + 3.*oz*r2*r2) - tor3z8;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;

    //val[i] *= -ezr*(120.f*oz6r2 + 120.f*oz5r + 60.f*oz4 + 18.f*r*oz3 + 3.f*r2*oz2) + 120.f*oz6r2;
  }

  return;
}

void eval_inr_4p_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz3*oz2;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz4*oz3;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 r2 = r*r; FP2 r3 = r2*r;
    FP2 ezr = exp(-zeta*r);
    FP2 for3z7 = 1440.*oz7/r3;
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 v1 = (ezr*(for3z7 + for3z7*r*zeta + 720.*oz5/r + 246.*oz4 + 66.*oz3*r + 15.*oz2*r2 + 3.*oz*r3) - for3z7);
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;

    //val[i] *= -ezr*(120.f*oz6r2 + 120.f*oz5r + 60.f*oz4 + 18.f*r*oz3 + 3.f*r2*oz2) + 120.f*oz6r2;
  }

  return;
}

void eval_inr_3p_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz3*oz2;
  FP2 oz6 = oz3*oz3;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 r2 = r*r;
    FP2 ezr = exp(-zeta*r);
    FP2 tor3z6 = 240.*oz6/r2/r;
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 v1 = (ezr*(tor3z6 + tor3z6*r*zeta + 120.*oz4/r + 42.*oz3 + 12.*oz2*r + 3.*oz*r2) - tor3z6);
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;

    //printf("  3p_r1d.  r: %9.6f  v1(b): %12.10f  add: %12.10f  v1: %12.10f \n",r,v1-(ezr-1.)*tor3z6,(ezr-1.)*tor3z6,v1);

    //val[i] *= -ezr*(120.f*oz6r2 + 120.f*oz5r + 60.f*oz4 + 18.f*r*oz3 + 3.f*r2*oz2) + 120.f*oz6r2;
  }

  return;
}

void eval_inr_2p_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz3*oz2;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP2 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 r2 = r*r;
    FP2 ezr = exp(-zeta*r);
    FP2 foz5r3 = 48.*oz5/r2/r;
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 v1 = (ezr*(foz5r3 + foz5r3*r*zeta + 24.*oz3/r + 9.*oz2 + 3.*r*oz) - foz5r3);
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;

    //val[i] *= -ezr*(toz5r2+24.f*oz4/r+12.f*oz3+3.f*r*oz2) + toz5r2;
  }

  return;
}


void eval_inr_8s_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;
  FP2 oz9 = oz4*oz5;
  FP2 oz10 = oz5*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);
    FP2 r2 = r*r; FP2 r3 = r2*r; FP2 r4 = r2*r2; FP2 r5 = r3*r2;
    FP2 for2z10 = 362880.*oz10/r2;
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 v1 = ezr*(for2z10 + for2z10*zeta*r + 181440.*oz8 + 60480.*r*oz7 + 15120.*r2*oz6 + 3024.*r3*oz5 + 504.*r4*oz4 + 72.*r5*oz3 + 9.*r3*r3*oz2 + r3*r4*oz) - for2z10;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_7s_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;
  FP2 oz9 = oz4*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);
    FP2 r2 = r*r; FP2 r3 = r2*r; FP2 r4 = r2*r2; FP2 r5 = r3*r2;
    FP2 for2z9 = 40320.*oz9/r2;
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 v1 = ezr*(for2z9 + for2z9*zeta*r + 20160.*oz7 + 6720.*r*oz6 + 1680.*r2*oz5 + 336.*r3*oz4 + 56.*r4*oz3 + 8.*r5*oz2 + r3*r3*oz) - for2z9;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_6s_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);
    FP2 r2 = r*r; FP2 r3 = r2*r; FP2 r4 = r2*r2; FP2 r5 = r3*r2;
    FP2 for2z8 = 5040.*oz8/r2;
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 v1 = ezr*(for2z8 + for2z8*zeta*r + 2520.*oz6 + 840.*r*oz5 + 210.*r2*oz4 + 42.*r3*oz3 + 7.*r4*oz2 + r5*oz) - for2z8;
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_5s_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);
    FP2 r2 = r*r; FP2 r3 = r2*r; FP2 r4 = r2*r2;
    FP2 sor2z7 = 720.*oz7/r2;
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 v1 = (ezr*(sor2z7 + sor2z7*zeta*r + 360.*oz5 + 120.*r*oz4 + 30.*r2*oz3 + 6.*r3*oz2 + r4*oz) - sor2z7);
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_4s_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);
    FP2 r2 = r*r; FP2 r3 = r2*r;
    FP2 or2z6 = 120.*oz6/r2;
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 v1 = (ezr*(or2z6 + or2z6*zeta*r + 60.*oz4 + 20.*r*oz3 + 5.*r2*oz2 + r3*oz) - or2z6);
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_3s_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);
    FP2 r2 = r*r;
    FP2 tor2z5 = 24.*oz5/r2;
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 v1 = (ezr*(tor2z5 + tor2z5*zeta*r + 12.*oz3 + 4.*r*oz2 + r2*oz) - tor2z5);
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;

    //val[i] *= -ezr*(24.f*oz5r + 18.f*oz4 + 6.f*roz3 + r2oz2) + 24.f*oz5r;
  }

  return;
}

void eval_inr_2s_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1.f/zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 ezr = expf(-zeta*r);
    FP2 soz4r2 = 6.*oz4/r/r;
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 v1 = (ezr*(soz4r2 + 6.*oz3/r + 3.*oz2 + r*oz) - soz4r2);
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}

void eval_inr_1s_r1d(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  //FP1 oz = 1.f/zeta;
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:3*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid[6*i+0]; FP1 y = grid[6*i+1]; FP1 z = grid[6*i+2];
    FP1 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);
    FP2 tz3r2 = 2.*oz3/r/r;
    FP2 xor1 = x/r; FP2 yor = y/r; FP2 zor = z/r;

    FP2 v1 = (ezr*(tz3r2 + 2.*oz2/r + oz) - tz3r2);
    val[3*i+0] *= v1*xor1;
    val[3*i+1] *= v1*yor;
    val[3*i+2] *= v1*zor;
  }

  return;
}


////r1 functions////

void eval_inr_6h_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  //FP1 oz = 1.f/zeta;
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz3*oz;
  FP2 oz5 = oz4*oz;
  FP2 oz6 = oz5*oz;
  FP2 oz7 = oz6*oz;
  FP2 oz8 = oz7*oz;
  FP2 oz9 = oz8*oz;
  FP2 oz10 = oz9*oz;
  FP2 oz11 = oz10*oz;
  FP2 oz12 = oz11*oz;
  FP2 oz13 = oz12*oz;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP2 r = grid[6*i+3];

    FP2 r2 = r*r; FP2 r3 = r2*r; FP2 r4 = r2*r2; FP2 r5 = r2*r3; FP2 r6 = r3*r3;

#if 1
    FP2 ezr = exp(-zeta*r);
    FP2 foz13r6 = 479001600.*oz13/r6;
  
    val[i] *= -ezr*(foz13r6 + 479001600.*oz12/r5 + 239500800.*oz11/r4 + 79833600.*oz10/r3 + 19958400.*oz9/r2
                    + 3991680.*oz8/r + 665280.*oz7 + 95040.*r*oz6 + 11880.*r2*oz5 + 1320.*r3*oz4 + 132.*r4*oz3 + 11.*r5*oz2)
              + foz13r6;
  #if RGLIMIT
    if (r<0.1) val[i] = 0.;
  #endif
#else
    FP1 ezr = expf(-zeta*r);
    FP1 foz13r6 = 479001600.f*oz13/r6;
  
    val[i] *= -ezr*(foz13r6 + 479001600.f*oz12/r5 + 239500800.f*oz11/r4 + 79833600.f*oz10/r3 + 19958400.f*oz9/r2
                    + 3991680.f*oz8/r + 665280.f*oz7 + 95040.f*r*oz6 + 11880.f*r2*oz5 + 1320.f*r3*oz4 + 132.f*r4*oz3 + 11.f*r5*oz2)
              + foz13r6;
#endif
  }

  return;
}

void eval_inr_5g_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  //FP1 oz = 1.f/zeta;
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz3*oz;
  FP2 oz5 = oz4*oz;
  FP2 oz6 = oz5*oz;
  FP2 oz7 = oz6*oz;
  FP2 oz8 = oz7*oz;
  FP2 oz9 = oz8*oz;
  FP2 oz10 = oz9*oz;
  FP2 oz11 = oz10*oz;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP2 r = grid[6*i+3];

    FP2 r2 = r*r; FP2 r3 = r2*r; FP2 r4 = r2*r2; FP2 r5 = r2*r3;

#if 1
    FP2 ezr = exp(-zeta*r);
    FP2 toz11r5 = 3628800.*oz11/r5;

    val[i] *= -ezr*(toz11r5 + 3628800.*oz10/r4 + 1814400.*oz9/r3 + 604800.*oz8/r2 + 151200.*oz7/r
                    + 30240.*oz6 + 5040.*r*oz5 + 720.*r2*oz4 + 90.*r3*oz3 + 9.*r4*oz2) + toz11r5;
   #if RGLIMIT
    if (r<0.1) val[i] = 0.;
   #endif
#else
    FP1 ezr = expf(-zeta*r);
    FP1 toz11r5 = 3628800.f*oz11/r5;

    val[i] *= -ezr*(toz11r5 + 3628800.f*oz10/r4 + 1814400.f*oz9/r3 + 604800.f*oz8/r2 + 151200.f*oz7/r
                    + 30240.f*oz6 + 5040.f*r*oz5 + 720.f*r2*oz4 + 90.f*r3*oz3 + 9.f*r4*oz2) + toz11r5;
#endif
  }

  return;
}

void eval_inr_6f_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
 //untested
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz3*oz;
  FP2 oz5 = oz4*oz;
  FP2 oz6 = oz5*oz;
  FP2 oz7 = oz6*oz;
  FP2 oz8 = oz7*oz;
  FP2 oz9 = oz8*oz;
  FP2 oz10 = oz9*oz;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP2 r = grid[6*i+3];

    FP2 r2 = r*r; FP2 r3 = r2*r; FP2 r4 = r2*r2;

    FP2 ezr = exp(-zeta*r);
    FP2 toz11r4 = 3628800.*oz9/r4;

    val[i] *= oz2*(-ezr*(toz11r4 + toz11r4*r*zeta + 1814400.*oz7/r2 + 604800.*oz6/r + 151200.*oz5
                    + 30240.*oz4*r + 5040.*r2*oz3 + 714.*r3*oz2 + 84.*r4*oz + 7.*r3*r2) + toz11r4);
  }

  return;
}

void eval_inr_5f_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
 //untested
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz3*oz;
  FP2 oz5 = oz4*oz;
  FP2 oz6 = oz5*oz;
  FP2 oz7 = oz6*oz;
  FP2 oz8 = oz7*oz;
  FP2 oz9 = oz8*oz;
  FP2 oz10 = oz9*oz;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP2 r = grid[6*i+3];

    FP2 r2 = r*r; FP2 r3 = r2*r; FP2 r4 = r2*r2;

    FP2 ezr = exp(-zeta*r);
    FP2 toz10r4 = 362880.*oz8/r4;

    val[i] *= oz2*(-ezr*(toz10r4 + toz10r4*r*zeta + 181440.*oz6/r2 + 60480.*oz5/r + 15120.*oz4
                    + 3024.*oz3*r + 504.*r2*oz2 + 70.*r3*oz + 7.*r4) + toz10r4);
  }

  return;
}

void eval_inr_4f_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  //FP1 oz = 1./zeta;
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz*oz2;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;
  FP2 oz9 = oz4*oz5;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP2 r = grid[6*i+3];

    FP2 r2 = r*r; FP2 r3 = r*r2; FP2 r4 = r2*r2;

    FP2 oz9r4 = oz9/r4;
    FP2 oz8r3 = oz8/r3;
    FP2 oz7r2 = oz7/r2;
    FP2 oz6r = oz6/r;
    FP2 roz4 = r*oz4;
    FP2 r2oz3 = r2*oz3;
    FP2 r3oz2 = r3*oz2;

#if 1
    FP2 ezr = exp(-zeta*r);
    FP2 foz9r4 = 40320.*oz9r4;

    val[i] *= -ezr*(foz9r4 + 40320.*oz8r3 + 20160.*oz7r2 + 6720.*oz6r + 1680.*oz5 + 336.*roz4 + 56.*r2oz3 + 7.*r3oz2)+foz9r4;
#else  
    FP1 ezr = expf(-zeta*r);
    FP1 foz9r4 = 40320.f*oz9r4;

    val[i] *= -ezr*(foz9r4 + 40320.f*oz8r3 + 20160.f*oz7r2 + 6720.f*oz6r + 1680.f*oz5 + 336.f*roz4 + 56.f*r2oz3 + 7.f*r3oz2)+foz9r4;
#endif
  }

  return;
}

void eval_inr_7d_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
 //untested
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz*oz2;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz2*oz4;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;
  FP2 oz9 = oz5*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP2 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);
    FP2 r2 = r*r; FP2 r3 = r2*r; FP2 r4 = r2*r2;

    FP2 foz11r3 = 3628800.*oz9/r3;
    val[i] *= oz2*(foz11r3 - ezr*(foz11r3 + foz11r3*zeta*r + 1814400.*oz7/r + 604800.*oz6 + 151200.*r*oz5 + 30120.*r2*oz4 + 4920.*r3*oz3 + 660.*r4*oz2 + + 70.*r3*r2*oz + 5.*r3*r3));
  }
 
  return;
}

void eval_inr_6d_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
 //untested
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz*oz2;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz2*oz4;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;
  FP2 oz9 = oz5*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP2 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);
    FP2 r2 = r*r; FP2 r3 = r2*r; FP2 r4 = r2*r2;

    FP2 foz10r3 = 362880.*oz8/r3;
    val[i] *= oz2*(foz10r3 - ezr*(foz10r3 + foz10r3*zeta*r + 181440.*oz6/r + 60480.*oz5 + 15120.*r*oz4 + 3000.*r2*oz3 + 480.*r3*oz2 + 60.*r4*oz + 5.*r3*r2));
  }
 
  return;
}

void eval_inr_5d_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
 //untested
  //FP1 oz = 1.f/zeta;
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz*oz2;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz2*oz4;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;
  FP2 oz9 = oz5*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP2 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);
    FP2 r2 = r*r; FP2 r3 = r2*r; FP2 r4 = r2*r2;

    FP2 foz9r3 = 40320.*oz7/r3;

    val[i] *= oz2*(foz9r3 - ezr*(foz9r3 + 40320.*oz6/r2 + 20160.*oz5/r + 6720.*oz4 + 1680.*r*oz3 + 330.*r2*oz2 + 50.*r3*oz + 5.*r4));
  }
 
  return;
}

void eval_inr_4d_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz*oz2;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz2*oz4;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);
    FP2 r2 = r*r; FP2 r3 = r2*r;

    FP2 foz8r3 = 5040.*oz8/r3;

    val[i] *= foz8r3 - ezr*(foz8r3 + 5040.*oz7/r2 + 2520.*oz6/r + 840.*oz5 + 210.*r*oz4 + 40.*r2*oz3 + 5.*r3*oz2);
  }
 
  return;
}

void eval_inr_3d_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
#if 1
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz*oz2;
  FP2 oz4 = oz2*oz2;
  FP2 oz7 = oz3*oz4;
#else
  FP1 oz = 1.f/zeta;
  FP1 oz2 = oz*oz;
  FP1 oz3 = oz*oz2;
  FP1 oz4 = oz2*oz2;
  FP1 oz7 = oz3*oz4;
#endif

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r = grid[6*i+3];
    FP1 r2 = r*r; FP1 r3 = r2*r;

   #if 1
    FP2 oz6r2 = oz2*oz4/r2;
    FP2 oz5r = oz2*oz3/r;

    FP2 ezr = exp(-zeta*r);
    FP2 oz7r3 = oz7/r3;

    val[i] *= -ezr*(720.*oz7r3 + 720.*oz6r2 + 360.*oz5r + 120.*oz4 + 30.*r*oz3 + 5.*r2*oz2) + 720.*oz7r3;
   #else
    FP1 oz6r2 = oz2*oz4/r2;
    FP1 oz5r = oz2*oz3/r;

    FP1 ezr = expf(-zeta*r);
    FP1 oz7r3 = oz7/r3;

    val[i] *= -ezr*(720.f*oz7r3 + 720.f*oz6r2 + 360.f*oz5r + 120.f*oz4 + 30.f*r*oz3 + 5.f*r2*oz2) + 720.f*oz7r3;
   #endif
  }

  return;
}

void eval_inr_7p_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
 //untested

  //FP1 oz = 1.f/zeta;
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz3*oz2;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP2 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);
    FP2 r2 = r*r;
    FP2 r3 = r2*r;

    FP2 fr2z10 = 362880.*oz8/r2;
    val[i] *= oz2*(fr2z10 - ezr*(fr2z10 + fr2z10*zeta*r + 181440.*oz6 + 59760.*oz5*r + 14400.*oz4*r2 + 2664.*oz3*r3 + 384.*r2*r2*oz2 + 42.*oz*r3*r2 + 3.*r3*r3));
  }

  return;
}

void eval_inr_6p_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
 //untested

  //FP1 oz = 1.f/zeta;
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz3*oz2;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP2 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);
    FP2 r2 = r*r;
    FP2 r3 = r2*r;

    FP2 fr2z9 = 40320.*oz7/r2;
    val[i] *= oz2*(fr2z9 - ezr*(fr2z9 + fr2z9*zeta*r + 20160.*oz5 + 6600.*oz4*r + 1560.*oz3*r2 + 276.*oz2*r3 + 36.*oz*r2*r2 + 3.*r3*r2));
  }

  return;
}

void eval_inr_5p_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
 //untested

  //FP1 oz = 1.f/zeta;
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz3*oz2;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP2 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);
    FP2 r2 = r*r;
    FP2 r3 = r2*r;

#if 1
    FP2 fr2z8 = 5040.*oz6/r2;
    val[i] *= oz2*(fr2z8 - ezr*(fr2z8 + fr2z8*zeta*r + 2520.*oz4 + 816.*oz3*r + 186.*oz2*r2 + 30.*oz*r3 + 3.*r2*r2));
#else
    FP1 fr2z8 = 5040.f*oz8/r2;
    val[i] *= fr2z8 - ezr*(fr2z8 + 5040.f*oz7/r + 2520.f*oz6 + 816.f*oz5*r + 186.f*oz4*r2 + 30.f*oz3*r3 + 3.f*oz2*r2*r2);
#endif
  }

  return;
}

void eval_inr_4p_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
 #if 1
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz3*oz2;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4; 
 #else
  FP1 oz = 1.f/zeta;
  FP1 oz2 = oz*oz;
  FP1 oz3 = oz2*oz;
  FP1 oz4 = oz2*oz2;
  FP1 oz5 = oz3*oz2;
  FP1 oz6 = oz3*oz3;
  FP1 oz7 = oz3*oz4; 
 #endif

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r = grid[6*i+3];
    FP1 r2 = r*r;
    FP1 r3 = r2*r;

   #if 1
    FP2 ezr = expf(-zeta*r);
    FP2 soz7r2 = 720.*oz7/r2;

    FP2 v1 = soz7r2 - ezr*(soz7r2 + 720.*oz6/r + 360.*oz5 + 114.*r*oz4 + 24.*r2*oz3 + 3.*r3*oz2);
    val[i] *= v1;
   #else
    FP1 ezr = expf(-zeta*r);
    FP1 soz7r2 = 720.f*oz7/r2;

    FP1 v1 = soz7r2 - ezr*(soz7r2 + 720.f*oz6/r + 360.f*oz5 + 114.f*r*oz4 + 24.f*r2*oz3 + 3.f*r3*oz2);
    val[i] *= v1;
   #endif
  }
 
  return;
}

void eval_inr_3p_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz3*oz2;
  FP2 oz6 = oz3*oz3;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r = grid[6*i+3];
    FP1 r2 = r*r;

   #if 1
    FP2 ezr = exp(-zeta*r);
    FP2 oz5r = oz5/r;
    FP2 oz6r2 = 120.*oz6/r2;

    val[i] *= -ezr*(oz6r2 + 120.*oz5r + 60.*oz4 + 18.*r*oz3 + 3.*r2*oz2) + oz6r2;
   #else
    FP1 ezr = expf(-zeta*r);
    FP1 oz5r = oz5/r;
    FP1 oz6r2 = 120.f*oz6/r2;

    val[i] *= -ezr*(oz6r2 + 120.f*oz5r + 60.f*oz4 + 18.f*r*oz3 + 3.f*r2*oz2) + oz6r2;
   #endif
  }

  return;
}

void eval_inr_2p_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz3*oz2;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r = grid[6*i+3];
    FP1 r2 = r*r;
    FP2 ezr = exp(-zeta*r);
    FP2 toz5r2 = 24.f*oz5/r2;

    val[i] *= -ezr*(toz5r2+24.f*oz4/r+12.f*oz3+3.f*r*oz2) + toz5r2;
  }

  return;
}

void eval_inr_8s_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);
    FP2 r2 = r*r; FP2 r3 = r2*r; FP2 r4 = r2*r2;

    FP2 forz10 = 362880./r*oz8;
    val[i] *= oz2*(forz10 - ezr*(forz10 + 322560.*oz7 + 141120.*r*oz6 + 40320.*r2*oz5 + 8400.*r3*oz4 + 1344.*r4*oz3 + 168.*r3*r2*oz2 + 16.*r3*r3*oz + r4*r3));
  }

  return;
}

void eval_inr_7s_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);
    FP2 r2 = r*r; FP2 r3 = r2*r;

    FP2 forz9 = 40320./r*oz7;
    val[i] *= oz2*(forz9 - ezr*(forz9 + 35280.*oz6 + 15120.*r*oz5 + 4200.*r2*oz4 + 840.*r3*oz3 + 126.*r2*r2*oz2 + 14.*r3*r2*oz + r3*r3));
  }

  return;
}

void eval_inr_6s_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;
  FP2 oz8 = oz4*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);
    FP2 r2 = r*r; FP2 r3 = r2*r;

    FP2 forz8 = 5040./r*oz8;
    val[i] *= forz8 - ezr*(forz8 + 4320.*oz7 + 1800.*r*oz6 + 480.*r2*oz5 + 90.*r3*oz4 + 12.*r2*r2*oz3 + r2*r3*oz2);
  }

  return;
}

void eval_inr_5s_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
 //untested

  //FP1 oz = 1.f/zeta;
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;
  FP2 oz7 = oz3*oz4;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP2 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);
    FP2 r2 = r*r;

#if 1
    FP2 sorz7 = 720./r*oz7;
    val[i] *= sorz7 - ezr*(sorz7 + 600.*oz6 + 240.*r*oz5 + 60.*r2*oz4 + 10.*r2*r*oz3 + r2*r2*oz2);
#else
    FP1 sorz7 = 720.f/r*oz7;
    val[i] *= sorz7 - ezr*(sorz7 + 600.f*oz6 + 240.f*r*oz5 + 60.f*r2*oz4 + 10.f*r2*r*oz3 + r2*r2*oz2);
#endif
  }

  return;
}

void eval_inr_4s_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz = 1./zeta;
  FP2 oz2 = oz*oz;
  FP2 oz3 = oz2*oz;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;
  FP2 oz6 = oz3*oz3;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);

    FP2 r2 = r*r;
    FP2 r3 = r2*r;
    FP2 ooz6r = 120.*oz6/r;

    val[i] *= ooz6r - ezr*(ooz6r + 96.*oz5 + 36.*r*oz4 + 8.*r2*oz3 + r3*oz2);
  }

  return;
}

void eval_inr_3s_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz2 = 1./zeta/zeta;
  FP2 oz3 = oz2/zeta;
  FP2 oz4 = oz2*oz2;
  FP2 oz5 = oz2*oz3;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);

    FP2 oz5r = oz5/r;
    FP2 r2oz2 = r*r*oz2;
    FP2 roz3 = oz3*r;
    FP2 toz5r = 24.*oz5r;

    val[i] *= -ezr*(toz5r + 18.*oz4 + 6.*roz3 + r2oz2) + toz5r;
  }

  return;
}

void eval_inr_2s_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz2 = 1./zeta/zeta;
  FP2 oz3 = oz2/zeta;
  FP2 oz4 = oz2*oz2;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r = grid[6*i+3];
  #if 1
    FP2 ezr = exp(-zeta*r);
    FP2 roz2 = r*oz2;
    FP2 soz4r = 6.*oz4/r;

    val[i] *= -ezr*(4.*oz3+roz2+soz4r)+soz4r;
  #else
    FP2 ezr = exp(-zeta*r);
    FP2 roz2 = r*oz2;
    FP2 soz4r = 6.f*oz4/r;

    val[i] *= -ezr*(4.f*oz3+roz2+soz4r)+soz4r;
  #endif
  }

  return;
}

void eval_inr_1s_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  FP2 oz2 = 1./zeta/zeta;
  FP2 toz3 = 2.*oz2/zeta;

#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r = grid[6*i+3];
    FP2 ezr = exp(-zeta*r);
    FP2 tz3r = toz3/r;
    val[i] *= -ezr*(oz2+tz3r)+tz3r;
  }

  return;
}

void eval_inr_d(int gsin, FP1* gridin, FP1* valin, int n1, int l1, FP1 zeta1)
{

  #if 1
  #if EVL64
  #pragma acc parallel loop independent present(gridin[0:6*gsin],valin[0:3*gsin])
  for (int i = 0; i < gsin; i++) {
      FP2 r = gridin[6*i+3];
      FP2 v = dvinr_gam(n1,l1,r,zeta1);
      valin[3*i+0] *= v*gridin[6*i+0]/r;
      valin[3*i+1] *= v*gridin[6*i+1]/r;
      valin[3*i+2] *= v*gridin[6*i+2]/r;
    }
  return;
  #else
  #pragma acc parallel loop independent present(gridin[0:6*gsin],valin[0:3*gsin])
  for (int i = 0; i < gsin; i++) {
      FP1 r = gridin[6*i+3];
      FP1 v = dvinr_gamf(n1,l1,r,zeta1);
      valin[3*i+0] *= v*gridin[6*i+0]/r;
      valin[3*i+1] *= v*gridin[6*i+1]/r;
      valin[3*i+2] *= v*gridin[6*i+2]/r;
    }
  return;
  #endif
  #else
  if (1) {
    // FP1 rc = (FP1)get_rcutoff(n1,l1,zeta1);
    int idx_rc = 0;
    // #pragma acc host_data use_device(gridin)
    // {
    //   idx_rc = lower_bound(gridin,gsin,rc);
    // }
    // // printf("n:%d l:%d z: % .3e rc);
    // dvinr_tay(n1,l1,idx_rc,gridin,valin,zeta1);
    FP1 * grid = gridin + (idx_rc*6);
    FP1 * val = valin + idx_rc;
    int gs = gsin - idx_rc;

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
  #endif
}

void eval_inr_r12(int gsin, FP1* gridin, FP1* valin, int n1, int l1, FP1 zeta1)
{

  #if EVL64
  #pragma acc parallel loop independent present(gridin[0:6*gsin],valin[0:gsin])
  for (int i = 0; i < gsin; i++) {
    valin[i] *= vinr_gam(n1,l1,gridin[6*i+3],zeta1);
  }
  #else
  #pragma acc parallel loop independent present(gridin[0:6*gsin],valin[0:gsin])
  for (int i = 0; i < gsin; i++) {
    valin[i] *= vinr_gamf(n1,l1,gridin[6*i+3],zeta1);
  }
  return;
  #endif

  #if 0
  if (1)
  {
    // FP1 rc = (FP1)get_rcutoff(n1,l1,zeta1);
    int idx_rc = 0;
    // #pragma acc host_data use_device(gridin)
    // {
    //  idx_rc = lower_bound(gridin,gsin,rc);
    // }
    // printf("n:%d l:%d z: % .3e rc);
    // vinr_tay(n1,l1,idx_rc,gridin,valin,zeta1);
    FP1 * grid = gridin + (idx_rc*6);
    FP1 * val = valin + idx_rc;
    int gs = gsin - idx_rc;
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

  return;
  #endif
}

void eval_exp_r2(int tid, int gs, FP1* grid, FP1* val, FP1 zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r = grid[6*i+4];
    FP1 ezr = expf(-zeta*r);
    val[i] *= ezr;
  }
  return;
}

void eval_exp_r1(int tid, int gs, FP1* grid, FP1* val, FP1 zeta)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs]) //async(tid)
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r = grid[6*i+3];
    FP1 ezr = expf(-zeta*r);
    val[i] *= ezr;
  }
  return;
}

void eval_exp_r2(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  eval_exp_r2(0,gs,grid,val,zeta);
 #if USE_ACC
  //#pragma acc wait
 #endif
}

void eval_exp_r1(int gs, FP1* grid, FP1* val, FP1 zeta)
{
  eval_exp_r1(0,gs,grid,val,zeta);
 #if USE_ACC
  //#pragma acc wait
 #endif
}

void get_inr_2s_f(int nrad, FP1 zeta, FP1* r, FP1* inr)
{
  FP1 oz2 = 1./zeta/zeta;
  FP1 oz3 = oz2/zeta;

#if USE_ACC
 #pragma acc parallel loop independent present(r[0:nrad],inr[0:nrad])
#endif
  for (int i=0;i<nrad;i++)
  {
    FP1 r1 = r[i];
    FP1 zr = zeta*r1;
    FP1 roz2 = r1*oz2;
    FP1 soz4r = 6.f*oz2*oz2/r1;
    FP1 ezr = expf(-zeta*r1);

    FP1 t1 = -ezr*(4.f*oz3+roz2+soz4r);
    FP1 t2 = soz4r;
    inr[i] = t1+t2;
  }

  return;
}

void get_inr_2p_f(int nrad, FP1 zeta, FP1* r, FP1* inr)
{
  FP1 oz = 1./zeta;
  FP1 oz2 = oz*oz;
  FP1 oz3 = oz2*oz;
  FP1 oz4 = oz2*oz2;
  FP1 oz5 = oz3*oz2;

#if USE_ACC
 #pragma acc parallel loop independent present(r[0:nrad],inr[0:nrad])
#endif
  for (int i=0;i<nrad;i++)
  {
    FP1 r1 = r[i];
    FP1 zr = zeta*r1;
    FP1 r2 = r1*r1;
    FP1 ezr = expf(-zeta*r1);
    FP1 toz5r2 = 24.*oz5/r2;

    FP1 t1 = -ezr*(toz5r2+24.*oz4/r1+12.*oz3+3.*r1*oz2);
    FP1 t2 = toz5r2;
    inr[i] = t1+t2;
  }

  return;
}

void get_inr_1s_f(int nrad, FP1 zeta, FP1* r, FP1* inr)
{
  FP1 oz2 = 1./zeta/zeta;
  FP1 oz3 = oz2/zeta;

#if USE_ACC
 #pragma acc parallel loop independent present(r[0:nrad],inr[0:nrad])
#endif
  for (int i=0;i<nrad;i++)
  {
    FP1 r1 = r[i];
    FP1 zr = zeta*r1;
    FP1 ezr = expf(-zeta*r1);
    FP1 tz3r = 2.*oz3/r1;
  
    FP1 t1 = -ezr*(oz2 + tz3r);
    FP1 t2 = tz3r;

    inr[i] = t1+t2;
  }

  return;
}

void get_inr_1s(int nrad, FP2 zeta, FP2* r, FP2* inr)
{
  FP2 oz2 = 1./zeta/zeta;
  FP2 oz3 = oz2/zeta;

  for (int i=0;i<nrad;i++)
  {
    FP2 r1 = r[i];
    FP2 zr = zeta*r1;
    FP2 ezr = exp(-zeta*r1);
    FP2 tz3r = 2.*oz3/r1;
  
    FP2 t1 = -ezr*(oz2 + tz3r);
    FP2 t2 = tz3r;

    inr[i] = t1+t2;
  }

  return;
}

//inr over radial grid only
void get_inr(int n1, int l1, FP2 zeta1, int nrad, FP1* r, FP1* inr)
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

FP2 norm_sh(int l, int m)
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
  if (l==1)
    return 0.4886025119029200;
  if (l==2)
  {
    if (m==0)
      return 0.3153915652525202;
    if (m==2)
      return 0.5462742152960397;
    return 1.092548430592080;
  }
  if (l==3)
  {
    if (m==-3) return 0.5900435899266435;
    if (m==-2) return 2.890611442640554;
    if (m==-1) return 0.4570457994644658;
    if (m== 0) return 0.3731763325901154;
    if (m== 1) return 0.4570457994644658;
    if (m== 2) return 1.445305721320277;
    if (m== 3) return 0.5900435899266435;
  }
  if (l==4)
  {
    if (m==-4) return 2.503342941797;
    if (m==-3 || m==3) return 1.770130769780;
    if (m==-2) return 0.946174695758;
    if (m==-1) return 0.669046543557;
    //if (m==-1) return 0.423142187661;
    if (m== 0) return 0.105785546915;
    //if (m== 0) return 0.042963579755;
    if (m== 1) return 0.669046543557;
    if (m== 2) return 0.473087347879;
    if (m== 4) return 0.625835735449;
    //if (m==-4 || m==-3 || m==-2) return 0.62583574;
    //if (m==-1 || m== 0 || m== 1) return 2.50334294;
    //if (m== 2 || m== 3 || m== 4) return 1.77013077;
  }
  if (l==5)
  {
    if (m==-5) return 0.6563820568;
    if (m==-4) return 8.3026492595;
    if (m==-3) return 0.6563820568;
    if (m==-2) return 0.4892382994;
    if (m==-1) return 2.0756623149;
    if (m== 0) return 4.793536785;
    if (m== 1) return 0.4892382994;
    if (m== 2) return 0.4529466512;
    if (m== 3) return 2.3967683925;
    if (m== 4) return 0.1169503225;
    if (m== 5) return 0.4529466512;
  }
  if (l==6)
  {

  }

  printf(" ERROR: norm not implemented for l=%i m=%i \n",l,m);
  exit(1);
  return 1.;
}

FP2 norm_sv(int n, int l, int m, FP2 zeta)
{
  FP2 num = 4.*PI*pow(2.*zeta,n+0.5);
  FP2 den = sqrt(fact(2*n))*(2.*l+1.);
  FP2 val = num/den;
  val *= norm_sh(l,m);
  return val;
}

FP2 norm(int n, int l, int m, FP2 zeta)
{
  FP2 num = pow(2.*zeta,n+0.5);
  FP2 den = sqrt(fact(2*n));

  FP2 val = num/den;
  val *= norm_sh(l,m);
  return val;
}

int fact(int N) 
{
  if (N==0) return 1;
  return N*fact(N-1);
}

