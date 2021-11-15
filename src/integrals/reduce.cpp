#include "integrals.h"
#include "fp_def.h"

#if RED_DOUBLE
void reduce_2c1(int s1, int s2, int gs, FP1** val1, FP1** val3, int iN, int N, FP2* An)
#else
void reduce_2c1(int s1, int s2, int gs, FP1** val1, FP1** val3, int iN, int N, FP1* An)
#endif
{
  int N2 = N*N;

#if USE_ACC
 #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs],val3[0:iN][0:gs],An[0:N2])
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s1;i2<s2;i2++)
    {
      int ii1 = i1-s1;
      int ii2 = i2-s1;
     #if RED_DOUBLE
      FP2 val = 0.;
     #else
      FP1 val = 0.f;
     #endif

     #pragma acc loop reduction(+:val)
      for (int j=0;j<gs;j++)
        val += val1[ii1][j]*val3[ii2][j];
   
      An[i1*N+i2] = val;
    }
  } //MM reduction
  #pragma acc wait

#else
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s1;i2<s2;i2++)
  {
    int ii1 = i1-s1; int ii2 = i2-s1;
    FP1 val = 0.f;
    for (int j=0;j<gs;j++)
      val += val1[ii1][j]*val3[ii2][j];
    An[i1*N+i2] = val;
  }
#endif

  return;
}

#if RED_DOUBLE
void reduce_2c2(int s1, int s2, int s3, int s4, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, int iN, int N, FP2* An)
#else
void reduce_2c2(int s1, int s2, int s3, int s4, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, int iN, int N, FP1* An)
#endif
{
  int N2 = N*N;

#if USE_ACC
 #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs],val4[0:iN][0:gs],An[0:N2])
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s3;i2<s4;i2++)
    {
      int ii1 = i1-s1;
      int ii2 = i2-s3;
    
     #if RED_DOUBLE
      FP2 val = 0.;
     #else
      FP1 val = 0.f;
     #endif

     #pragma acc loop reduction(+:val)
      for (int j=0;j<gs;j++)
        val += val1[ii1][j]*val3[ii2][j];
     #pragma acc loop reduction(+:val)
      for (int j=0;j<gs;j++) 
        val += val2[ii1][j]*val4[ii2][j];
        
      An[i1*N+i2] = val;
    }
  } //MM reduction
  #pragma acc wait

#else
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3;
    FP1 val = 0.f;
    for (int j=0;j<gs;j++)
      val += val1[ii1][j]*val3[ii2][j];
    for (int j=0;j<gs;j++)
      val += val2[ii1][j]*val4[ii2][j];
    An[i1*N+i2] = val;
  }
#endif

  return;
}

#if RED_DOUBLE
void reduce_2c3(int s1, int s2, int s3, int s4, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, int iN, int N, FP2* An)
#else
void reduce_2c3(int s1, int s2, int s3, int s4, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, int iN, int N, FP1* An)
#endif
{
  int N2 = N*N;
  //printf("   reduce 2c3: %i-%i %i-%i \n",s1,s2,s3,s4);

#if USE_ACC
 #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs],val4[0:iN][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs],An[0:N2])
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s3;i2<s4;i2++)
    {
      int ii1 = i1-s1;
      int ii2 = i2-s3;
    
     #if RED_DOUBLE
      FP2 val = 0.;
     #else
      FP1 val = 0.f;
     #endif

     #pragma acc loop reduction(+:val)
      for (int j=0;j<gs;j++)
        val += val1[ii1][j]*val3[ii2][j];
     #pragma acc loop reduction(+:val)
      for (int j=0;j<gs;j++) 
        val += val2[ii1][j]*val4[ii2][j];
     #pragma acc loop reduction(+:val)
      for (int j=0;j<gs;j++) 
        val += val5[ii1][j]*val6[ii2][j];
        
      An[i1*N+i2] = val;
    }
  } //MM reduction
  #pragma acc wait

#else
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3;
    FP1 val = 0.f;
    for (int j=0;j<gs;j++)
      val += val1[ii1][j]*val3[ii2][j];
    for (int j=0;j<gs;j++)
      val += val2[ii1][j]*val4[ii2][j];
    for (int j=0;j<gs;j++)
      val += val5[ii1][j]*val6[ii2][j];
    An[i1*N+i2] = val;
  }
#endif

  return;
}

#if RED_DOUBLE
void reduce_2c1d(int m1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val3, FP1** val1x, FP1** val3x, int iN, int N, int natoms, FP2 scalar, FP2* xyz_grad)
#else
void reduce_2c1d(int m1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val3, FP1** val1x, FP1** val3x, int iN, int N, int natoms, FP2 scalar, FP1* xyz_g)
#endif
{
  int N2 = N*N;
  int N3 = 3*natoms;
  int gs3 = 3*gs;

  FP2 valxt = 0.;
  FP2 valyt = 0.;
  FP2 valzt = 0.;

#if USE_ACC
 #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs],val3[0:iN][0:gs],val1x[0:iN][0:gs3],val3x[0:iN][0:gs3],dpq[0:N2],norms[0:N2],xyz_grad[0:3*natoms]) reduction(+:valxt,valyt,valzt)
#endif
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s3;i2<s4;i2++)
    {
      int ii1 = i1-s1;
      int ii2 = i2-s3;
     #if RED_DOUBLE
      FP2 valx, valy, valz;
      valx = valy = valz = 0.;
     #else
      FP1 valx, valy, valz;
      valx = valy = valz = 0.f;
     #endif

     #pragma acc loop reduction(+:valx,valy,valz)
      for (int j=0;j<gs;j++)
      {
        FP1 v1 = val1[ii1][j];
        valx += v1*val3x[ii2][3*j+0]; valy += v1*val3x[ii2][3*j+1]; valz += v1*val3x[ii2][3*j+2];

        FP1 v3 = val3[ii2][j];
        valx += v3*val1x[ii1][3*j+0]; valy += v3*val1x[ii1][3*j+1]; valz += v3*val1x[ii1][3*j+2];
      }

      FP2 nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
      valx *= nd1; valy *= nd1; valz *= nd1;
      //printf("  i1/2(m):  %2i %2i (%i)  valxyz: %9.5f %9.5f %9.5f  norm: %8.3f dpq: %8.5f  nd1: %8.5f \n",i1,i2,m1,valx,valy,valz,norms[i1*N+i2],dpq[i1*N+i2],nd1);
      //printf("  i1/2(m):  %2i %2i (%i)  valxyz: %9.5f %9.5f %9.5f \n",i1,i2,m1,valx,valy,valz);

      valxt += valx; valyt += valy; valzt += valz;
    }

  } //MM reduction
  #pragma acc wait

 #pragma acc serial present(xyz_grad[0:N3])
  {
    xyz_grad[3*m1+0] += scalar*valxt;
    xyz_grad[3*m1+1] += scalar*valyt;
    xyz_grad[3*m1+2] += scalar*valzt;
  }

  return;
}

#if RED_DOUBLE
void reduce_2c1ds(int m1, int n1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3x, FP1** val4x, int iN, int N, int natoms, FP2 scalar, FP2* xyz_grad)
#else
void reduce_2c1ds(int m1, int n1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3x, FP1** val4x, int iN, int N, int natoms, FP2 scalar, FP1* xyz_g)
#endif
{
  int N2 = N*N;
  int N3 = 3*natoms;
  int gs3 = 3*gs;

  FP2 valxt = 0.;
  FP2 valyt = 0.;
  FP2 valzt = 0.;

#if USE_ACC
 #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs],val2[0:iN][0:gs],val3x[0:iN][0:gs3],val4x[0:iN][0:gs3],dpq[0:N2],norms[0:N2],xyz_grad[0:3*natoms]) reduction(+:valxt,valyt,valzt)
#endif
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s3;i2<s4;i2++)
    {
      int ii1 = i1-s1;
      int ii2 = i2-s3;
     #if RED_DOUBLE
      FP2 valx, valy, valz;
      valx = valy = valz = 0.;
     #else
      FP1 valx, valy, valz;
      valx = valy = valz = 0.f;
     #endif

     #pragma acc loop reduction(+:valx,valy,valz)
      for (int j=0;j<gs;j++)
      {
        FP1 v1 = val1[ii1][j];
        valx += v1*val3x[ii2][3*j+0];
        valy += v1*val3x[ii2][3*j+1];
        valz += v1*val3x[ii2][3*j+2];

        FP1 v2 = val2[ii1][j];
        valx += v2*val4x[ii2][3*j+0];
        valy += v2*val4x[ii2][3*j+1];
        valz += v2*val4x[ii2][3*j+2];
      }

      FP2 nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
      valx *= nd1; valy *= nd1; valz *= nd1;
      //printf("  i1/2(m):  %2i %2i (%i)  valxyz: %9.5f %9.5f %9.5f  norm: %8.3f dpq: %8.5f  nd1: %8.5f \n",i1,i2,m1,valx,valy,valz,norms[i1*N+i2],dpq[i1*N+i2],nd1);
      //if (nd1>0.) printf("  i1/2(m):  %2i %2i (%i)  valxyz: %9.5f %9.5f %9.5f \n",i1,i2,m1,valx,valy,valz);

      valxt += valx;
      valyt += valy;
      valzt += valz;
    }

  } //MM reduction
  #pragma acc wait

 #pragma acc serial present(xyz_grad[0:N3])
  {
    xyz_grad[3*m1+0] += scalar*valxt;
    xyz_grad[3*m1+1] += scalar*valyt;
    xyz_grad[3*m1+2] += scalar*valzt;
    xyz_grad[3*n1+0] -= scalar*valxt;
    xyz_grad[3*n1+1] -= scalar*valyt;
    xyz_grad[3*n1+2] -= scalar*valzt;
  }

  return;
}

#if RED_DOUBLE
void reduce_2c2d(int m1, int n1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, int iN, int N, int natoms, FP2 scalar, FP2* xyz_grad)
#else
void reduce_2c2d(int m1, int n1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, int iN, int N, int natoms, FP2 scalar, FP1* xyz_grad)
#endif
{
  int N2 = N*N;
  int gs3 = 3*gs;

  FP2 valxt1 = 0.; FP2 valyt1 = 0.; FP2 valzt1 = 0.;
  FP2 valxt2 = 0.; FP2 valyt2 = 0.; FP2 valzt2 = 0.;

 #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs],val4[0:iN][0:gs],norms[0:N2],val1x[0:iN][0:gs3],val2x[0:iN][0:gs3],val3x[0:iN][0:gs3],val4x[0:iN][0:gs3],dpq[0:N2]) reduction(+:valxt1,valxt2,valyt1,valyt2,valzt1,valzt2)
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s3;i2<s4;i2++)
    {
      FP2 valx1 = 0.; FP2 valy1 = 0.; FP2 valz1 = 0.;
      FP2 valx2 = 0.; FP2 valy2 = 0.; FP2 valz2 = 0.;
     #pragma acc loop reduction(+:valx1,valx2,valy1,valy2,valz1,valz2)
      for (int j=0;j<gs;j++)
      {
        int ii1 = i1-s1;
        int ii2 = i2-s3;

        FP2 v1 = val1[ii1][j]; FP2 v2 = val2[ii1][j];
        valx2 += v1*val3x[ii2][3*j+0]; valy2 += v1*val3x[ii2][3*j+1]; valz2 += v1*val3x[ii2][3*j+2];
        valx2 += v2*val4x[ii2][3*j+0]; valy2 += v2*val4x[ii2][3*j+1]; valz2 += v2*val4x[ii2][3*j+2];

        FP2 v3 = val3[ii2][j]; FP2 v4 = val4[ii2][j];
        valx1 += v3*val1x[ii1][3*j+0]; valy1 += v3*val1x[ii1][3*j+1]; valz1 += v3*val1x[ii1][3*j+2];
        valx1 += v4*val2x[ii1][3*j+0]; valy1 += v4*val2x[ii1][3*j+1]; valz1 += v4*val2x[ii1][3*j+2];
      }
      FP2 nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
      valx1 *= nd1; valy1 *= nd1; valz1 *= nd1;
      valx2 *= nd1; valy2 *= nd1; valz2 *= nd1;

      valxt1 += valx1; valyt1 += valy1; valzt1 += valz1;
      valxt2 += valx2; valyt2 += valy2; valzt2 += valz2;

      //printf("  i1/2(mn): %2i %2i (%i%i) valxyz: %9.5f %9.5f %9.5f  %9.5f %9.5f %9.5f  norm: %8.3f dpq: %8.5f \n",i1,i2,m1,n1,valx1,valy1,valz1,valx2,valy2,valz2,norms[i1*N+i2],dpq[i1*N+i2]);
      //if (nd1>0.)
      //  printf("  i1/2(mn): %2i %2i (%i%i) valxyz: %9.5f %9.5f %9.5f  %9.5f %9.5f %9.5f \n",i1,i2,m1,n1,valx1,valy1,valz1,valx2,valy2,valz2);
      //else if (nd1>0.)
      //  printf("  i1/2(mn): %2i %2i (%i%i) valxyz: %9.5f %9.5f %9.5f  \n",i1,i2,m1,n1,valx1+valx2,valy1+valy2,valz1+valz2);
    }
  }

 #if 0
  FP2 vxa = (valxt1 + valxt2)*0.5;
  FP2 vya = (valyt1 + valyt2)*0.5;
  FP2 vza = (valzt1 + valzt2)*0.5;
  valxt1 -= vxa; valyt1 -= vya; valzt1 -= vza;
  valxt2 -= vxa; valyt2 -= vya; valzt2 -= vza;
 #endif

 #pragma acc serial present(xyz_grad[0:3*natoms])
  {
    xyz_grad[3*m1+0] += scalar*valxt1;
    xyz_grad[3*m1+1] += scalar*valyt1;
    xyz_grad[3*m1+2] += scalar*valzt1;
    xyz_grad[3*n1+0] += scalar*valxt2;
    xyz_grad[3*n1+1] += scalar*valyt2;
    xyz_grad[3*n1+2] += scalar*valzt2;
  }

  return;
}

#if RED_DOUBLE
void reduce_2c2dh(int m1, int n1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3x, FP1** val4x, int iN, int N, int natoms, FP2 scalar, FP2* xyz_grad)
#else
void reduce_2c2dh(int m1, int n1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3x, FP1** val4x, int iN, int N, int natoms, FP2 scalar, FP1* xyz_grad)
#endif
{
  int N2 = N*N;
  int gs3 = 3*gs;

  FP2 valxt = 0.; FP2 valyt = 0.; FP2 valzt = 0.;

 #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs],val2[0:iN][0:gs],val3x[0:iN][0:gs3],val4x[0:iN][0:gs3],dpq[0:N2],norms[0:N2]) reduction(+:valxt,valyt,valzt)
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s3;i2<s4;i2++)
    {
      FP2 valx = 0.; FP2 valy = 0.; FP2 valz = 0.;
     #pragma acc loop reduction(+:valx,valy,valz)
      for (int j=0;j<gs;j++)
      {
        int ii1 = i1-s1;
        int ii2 = i2-s3;

        FP2 v1 = val1[ii1][j]; FP2 v2 = val2[ii1][j];
        valx += v1*val3x[ii2][3*j+0]; valy += v1*val3x[ii2][3*j+1]; valz += v1*val3x[ii2][3*j+2];
        valx += v2*val4x[ii2][3*j+0]; valy += v2*val4x[ii2][3*j+1]; valz += v2*val4x[ii2][3*j+2];
      }
      FP2 nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
      valx *= nd1; valy *= nd1; valz *= nd1;

      valxt += valx; valyt += valy; valzt += valz;

      //printf("  i1/2(mn): %2i %2i (%i%i) valxyz: %9.5f %9.5f %9.5f  %9.5f %9.5f %9.5f  norm: %8.3f dpq: %8.5f \n",i1,i2,m1,n1,valx1,valy1,valz1,valx2,valy2,valz2,norms[i1*N+i2],dpq[i1*N+i2]);
      //if (m1!=n1 && nd1>0.)
      //  printf("  i1/2(mn): %2i %2i (%i%i) valxyz: %9.5f %9.5f %9.5f  %9.5f %9.5f %9.5f \n",i1,i2,m1,n1,valx1,valy1,valz1,valx2,valy2,valz2);
      //else if (nd1>0.)
      //  printf("  i1/2(mn): %2i %2i (%i%i) valxyz: %9.5f %9.5f %9.5f  \n",i1,i2,m1,n1,valx1+valx2,valy1+valy2,valz1+valz2);
    }
  }

 #pragma acc serial present(xyz_grad[0:3*natoms])
  {
    xyz_grad[3*m1+0] += scalar*valxt;
    xyz_grad[3*m1+1] += scalar*valyt;
    xyz_grad[3*m1+2] += scalar*valzt;
    xyz_grad[3*n1+0] -= scalar*valxt;
    xyz_grad[3*n1+1] -= scalar*valyt;
    xyz_grad[3*n1+2] -= scalar*valzt;
  }

  return;
}

#if RED_DOUBLE
void reduce_2c2ds(int m1, int n1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, int iN, int N, int natoms, FP2 scalar, FP2* xyz_grad)
#else
void reduce_2c2ds(int m1, int n1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, int iN, int N, int natoms, FP2 scalar, FP1* xyz_grad)
#endif
{
  int N2 = N*N;
  int gs3 = 3*gs;

  FP2 valxt1 = 0.; FP2 valyt1 = 0.; FP2 valzt1 = 0.;
  FP2 valxt2 = 0.; FP2 valyt2 = 0.; FP2 valzt2 = 0.;

 #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs],val4[0:iN][0:gs],norms[0:N2],val1x[0:iN][0:gs3],val2x[0:iN][0:gs3],val3x[0:iN][0:gs3],val4x[0:iN][0:gs3],dpq[0:N2]) reduction(+:valxt1,valxt2,valyt1,valyt2,valzt1,valzt2)
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s3;i2<s4;i2++)
    {
      FP2 valx1 = 0.; FP2 valy1 = 0.; FP2 valz1 = 0.;
      FP2 valx2 = 0.; FP2 valy2 = 0.; FP2 valz2 = 0.;
     #pragma acc loop reduction(+:valx1,valx2,valy1,valy2,valz1,valz2)
      for (int j=0;j<gs;j++)
      {
        int ii1 = i1-s1;
        int ii2 = i2-s3;

        FP2 v1 = val1[ii1][j]; FP2 v2 = val2[ii1][j];
        valx2 += v1*val3x[ii2][3*j+0]; valy2 += v1*val3x[ii2][3*j+1]; valz2 += v1*val3x[ii2][3*j+2];
        valx2 += v2*val4x[ii2][3*j+0]; valy2 += v2*val4x[ii2][3*j+1]; valz2 += v2*val4x[ii2][3*j+2];

        FP2 v3 = val3[ii2][j]; FP2 v4 = val4[ii2][j];
        valx1 += v3*val1x[ii1][3*j+0]; valy1 += v3*val1x[ii1][3*j+1]; valz1 += v3*val1x[ii1][3*j+2];
        valx1 += v4*val2x[ii1][3*j+0]; valy1 += v4*val2x[ii1][3*j+1]; valz1 += v4*val2x[ii1][3*j+2];
      }
      FP2 nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
      valx1 *= nd1; valy1 *= nd1; valz1 *= nd1;
      valx2 *= nd1; valy2 *= nd1; valz2 *= nd1;

      valxt1 += valx1; valyt1 += valy1; valzt1 += valz1;
      valxt2 += valx2; valyt2 += valy2; valzt2 += valz2;

      //printf("  i1/2(mn): %2i %2i (%i%i) valxyz: %9.5f %9.5f %9.5f  %9.5f %9.5f %9.5f  norm: %8.3f dpq: %8.5f \n",i1,i2,m1,n1,valx1,valy1,valz1,valx2,valy2,valz2,norms[i1*N+i2],dpq[i1*N+i2]);
      //if (m1!=n1 && nd1>0.)
      //  printf("  i1/2(mn): %2i %2i (%i%i) valxyz: %9.5f %9.5f %9.5f  %9.5f %9.5f %9.5f \n",i1,i2,m1,n1,valx1,valy1,valz1,valx2,valy2,valz2);
      //else if (nd1>0.)
      //  printf("  i1/2(mn): %2i %2i (%i%i) valxyz: %9.5f %9.5f %9.5f  \n",i1,i2,m1,n1,valx1+valx2,valy1+valy2,valz1+valz2);
    }
  }

 #if 0
  FP2 vxa = (valxt1 + valxt2)*0.5;
  FP2 vya = (valyt1 + valyt2)*0.5;
  FP2 vza = (valzt1 + valzt2)*0.5;
  valxt1 -= vxa; valyt1 -= vya; valzt1 -= vza;
  valxt2 -= vxa; valyt2 -= vya; valzt2 -= vza;
 #endif

 #pragma acc serial present(xyz_grad[0:3*natoms])
  {
    xyz_grad[3*m1+0] += scalar*valxt1;
    xyz_grad[3*m1+1] += scalar*valyt1;
    xyz_grad[3*m1+2] += scalar*valzt1;
    xyz_grad[3*n1+0] -= scalar*valxt2;
    xyz_grad[3*n1+1] -= scalar*valyt2;
    xyz_grad[3*n1+2] -= scalar*valzt2;
  }

  return;
}

#if RED_DOUBLE
void reduce_2c3d(int m1, int n1, int p1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, FP1** val5x, FP1** val6x, int iN, int N, int natoms, FP2 scalar, FP2* xyz_grad)
#else
void reduce_2c3d(int m1, int n1, int p1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, int iN, int N, int natoms, FP2 scalar, FP1* xyz_grad)
#endif
{
  int N2 = N*N;
  int gs3 = 3*gs;

  FP2 valxt1 = 0.; FP2 valyt1 = 0.; FP2 valzt1 = 0.;
  FP2 valxt2 = 0.; FP2 valyt2 = 0.; FP2 valzt2 = 0.;

#if USE_ACC
 #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs],val4[0:iN][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs],val1x[0:iN][0:gs3],val2x[0:iN][0:gs3],val3x[0:iN][0:gs3],val4x[0:iN][0:gs3],val5x[0:iN][0:gs3],val6x[0:iN][0:gs3],norms[0:N2],dpq[0:N2],xyz_grad[0:3*natoms]) reduction(+:valxt1,valyt1,valzt1,valxt2,valyt2,valzt2)
#endif
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s3;i2<s4;i2++)
    {
      int ii1 = i1-s1;
      int ii2 = i2-s3;
    
      FP2 valx1 = 0.; FP2 valy1 = 0.; FP2 valz1 = 0.;
      FP2 valx2 = 0.; FP2 valy2 = 0.; FP2 valz2 = 0.;
     #pragma acc loop reduction(+:valx1,valx2,valy1,valy2,valz1,valz2)
      for (int j=0;j<gs;j++)
      {
        FP1 v1 = val1[ii1][j]; FP1 v2 = val2[ii1][j]; FP1 v5 = val5[ii1][j];
        valx2 += v1*val3x[ii2][3*j+0]; valx2 += v2*val4x[ii2][3*j+0]; valx2 += v5*val6x[ii2][3*j+0];
        valy2 += v1*val3x[ii2][3*j+1]; valy2 += v2*val4x[ii2][3*j+1]; valy2 += v5*val6x[ii2][3*j+1];
        valz2 += v1*val3x[ii2][3*j+2]; valz2 += v2*val4x[ii2][3*j+2]; valz2 += v5*val6x[ii2][3*j+2];
        FP1 v3 = val3[ii2][j]; FP1 v4 = val4[ii2][j]; FP1 v6 = val6[ii2][j];
        valx1 += v3*val1x[ii1][3*j+0]; valx1 += v4*val2x[ii1][3*j+0]; valx1 += v6*val5x[ii1][3*j+0];
        valy1 += v3*val1x[ii1][3*j+1]; valy1 += v4*val2x[ii1][3*j+1]; valy1 += v6*val5x[ii1][3*j+1];
        valz1 += v3*val1x[ii1][3*j+2]; valz1 += v4*val2x[ii1][3*j+2]; valz1 += v6*val5x[ii1][3*j+2];
      }

      FP2 nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
      valx1 *= nd1; valy1 *= nd1; valz1 *= nd1;
      valx2 *= nd1; valy2 *= nd1; valz2 *= nd1;
      //if (nd1>0.) printf(" 2c3d i1/2: %2i %2i valxyz: %9.5f %9.5f %9.5f  %9.5f %9.5f %9.5f  norm: %8.3f dpq: %8.5f \n",i1,i2,valx1,valy1,valz1,valx2,valy2,valz2,norms[i1*N+i2],dpq[i1*N+i2]);
      valxt1 += valx1; valyt1 += valy1; valzt1 += valz1;
      valxt2 += valx2; valyt2 += valy2; valzt2 += valz2;

    }
  } //MM reduction

 #pragma acc serial present(xyz_grad[0:3*natoms])
  {
    xyz_grad[3*m1+0] += scalar*valxt1; xyz_grad[3*m1+1] += scalar*valyt1; xyz_grad[3*m1+2] += scalar*valzt1;
    xyz_grad[3*n1+0] += scalar*valxt2; xyz_grad[3*n1+1] += scalar*valyt2; xyz_grad[3*n1+2] += scalar*valzt2;

    xyz_grad[3*p1+0] -= scalar*valxt1; xyz_grad[3*p1+1] -= scalar*valyt1; xyz_grad[3*p1+2] -= scalar*valzt1;
    xyz_grad[3*p1+0] -= scalar*valxt2; xyz_grad[3*p1+1] -= scalar*valyt2; xyz_grad[3*p1+2] -= scalar*valzt2;
  }

  return;
}

#if RED_DOUBLE
void reduce_2c3de(int m1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* Pao, FP1** val1, FP1** val2, FP1** val3, FP1** val1x, FP1** val2x, FP1** val3x, int N, int imaxN, int natoms, FP2 scalar, FP2* xyz_grad)
#else
void reduce_2c3de(int m1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* Pao, FP1** val1, FP1** val2, FP1** val3, FP1** val1x, FP1** val2x, FP1** val3x, int N, int imaxN, FP1 scalar, FP1* xyz_grad)
#endif
{
  int N2 = N*N;
  int N3 = 3*natoms;
  int gs3 = 3*gs;

  FP2 valxt = 0.; FP2 valyt = 0.; FP2 valzt = 0.;

#if USE_ACC
 #pragma acc parallel loop collapse(2) present(val1[0:imaxN][0:gs],val2[0:imaxN][0:gs],val3[0:imaxN][0:gs],val1x[0:imaxN][0:gs3],val2x[0:imaxN][0:gs3],val3x[0:imaxN][0:gs3],norms[0:N2],Pao[0:N2]) reduction(+:valxt,valyt,valzt)
#endif
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3;

    FP2 valx = 0.; FP2 valy = 0.; FP2 valz = 0.;

   #pragma acc loop reduction(+:valx,valy,valz)
    for (int j=0;j<gs;j++)
    {
      FP2 v1 = val1[ii1][j]; FP2 v2 = val2[ii1][j]; FP2 v3 = val3[ii1][j];

      valx += v1*val1x[ii2][3*j+0]; valy += v1*val1x[ii2][3*j+1]; valz += v1*val1x[ii2][3*j+2];
      valx += v2*val2x[ii2][3*j+0]; valy += v2*val2x[ii2][3*j+1]; valz += v2*val2x[ii2][3*j+2];
      valx += v3*val3x[ii2][3*j+0]; valy += v3*val3x[ii2][3*j+1]; valz += v3*val3x[ii2][3*j+2];
    }

    FP2 nd1 = Pao[i1*N+i2]*norms[i1*N+i2];
    valx *= nd1; valy *= nd1; valz *= nd1;
    valxt += valx; valyt += valy; valzt += valz;

    if (nd1>0.) printf(" 2c3de i12: %2i %2i valxyz: %9.5f %9.5f %9.5f  nd1: %7.3f \n",i1,i2,valx,valy,valz,nd1);
  } //i1,i2,i3

 #pragma acc serial present(xyz_grad[0:N3])
  {
    xyz_grad[3*m1+0] += scalar*valxt;
    xyz_grad[3*m1+1] += scalar*valyt;
    xyz_grad[3*m1+2] += scalar*valzt;
  }

  return;
}

#if RED_DOUBLE
void reduce_3c2d112(int m1, int n1, int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP2* norms1, FP2* norms2, FP2* dC, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, FP1** val5x, FP1** val6x, int N, int Naux, int imaxN, int imaxNa, int natoms, FP2 scalar, FP2* xyz_grad)
#else
void reduce_3c2d112(int m1, int n1, int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP2* norms1, FP2* norms2, FP2* dC, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, FP1** val5x, FP1** val6x, int N, int Naux, int imaxN, int imaxNa, FP1 scalar, FP1* xyz_grad)
#endif
{
  int N2 = N*N;
  int N2a = N2*Naux;
  int N3 = 3*natoms;
  int gs3 = 3*gs;

  FP2 valxt1 = 0.; FP2 valyt1 = 0.; FP2 valzt1 = 0.;
  FP2 valxt2 = 0.; FP2 valyt2 = 0.; FP2 valzt2 = 0.;

#if USE_ACC
 #pragma acc parallel loop collapse(3) present(val1[0:imaxNa][0:gs],val2[0:imaxN][0:gs],val3[0:imaxN][0:gs],val4[0:imaxNa][0:gs],val5[0:imaxN][0:gs],val6[0:imaxN][0:gs],val1x[0:imaxNa][0:gs3],val2x[0:imaxN][0:gs3],val3x[0:imaxN][0:gs3],val4x[0:imaxNa][0:gs],val5x[0:imaxN][0:gs],val6x[0:imaxN][0:gs],norms1[0:Naux],norms2[0:N2],dC[0:N2a]) reduction(+:valxt1,valyt1,valzt1,valxt2,valyt2,valzt2)
#endif
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  for (int i3=s5;i3<s6;i3++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3; int ii3 = i3-s5;

    FP2 valx1 = 0.; FP2 valy1 = 0.; FP2 valz1 = 0.;
    FP2 valx2 = 0.; FP2 valy2 = 0.; FP2 valz2 = 0.;

   #pragma acc loop reduction(+:valx1,valy1,valz1,valx2,valy2,valz2)
    for (int j=0;j<gs;j++)
    {
      FP2 v1 = val1[ii1][j]; FP2 v2 = val2[ii2][j]; FP2 v3 = val3[ii3][j];
      FP2 v12 = v1*v2; FP2 v13 = v1*v3; FP2 v23 = v2*v3;

      valx2 += v12*val3x[ii3][3*j+0]; valy2 += v12*val3x[ii3][3*j+1]; valz2 += v12*val3x[ii3][3*j+2];
      valx1 += v13*val2x[ii2][3*j+0]; valy1 += v13*val2x[ii2][3*j+1]; valz1 += v13*val2x[ii2][3*j+2];
      valx1 += v23*val1x[ii1][3*j+0]; valy1 += v23*val1x[ii1][3*j+1]; valz1 += v23*val1x[ii1][3*j+2];

      FP2 v4 = val4[ii1][j]; FP2 v5 = val5[ii2][j]; FP2 v6 = val6[ii3][j];
      FP2 v45 = v4*v5; FP2 v46 = v4*v6; FP2 v56 = v5*v6;

      valx2 += v45*val6x[ii3][3*j+0]; valy2 += v45*val6x[ii3][3*j+1]; valz2 += v45*val6x[ii3][3*j+2];
      valx1 += v46*val5x[ii2][3*j+0]; valy1 += v46*val5x[ii2][3*j+1]; valz1 += v46*val5x[ii2][3*j+2];
      valx1 += v56*val4x[ii1][3*j+0]; valy1 += v56*val4x[ii1][3*j+1]; valz1 += v56*val4x[ii1][3*j+2];
    }

    FP2 nd1 = dC[i1*N2+i2*N+i3]*norms1[i1]*norms2[i2*N+i3];
    valx1 *= nd1; valy1 *= nd1; valz1 *= nd1;
    valx2 *= nd1; valy2 *= nd1; valz2 *= nd1;
    valxt1 += valx1; valyt1 += valy1; valzt1 += valz1;
    valxt2 += valx2; valyt2 += valy2; valzt2 += valz2;

    //printf("  i13: %2i %2i %2i valxyz: %10.6f %10.6f %10.6f  %10.6f %10.6f %10.6f  norm: %8.3f dC: %8.5f \n",i1,i2,i3,valx1,valy1,valz1,valx2,valy2,valz2,norms1[i1]*norms2[i2*N+i3],dC[i1*N2+i2*N+i3]);
  } //i1,i2,i3
 
 #pragma acc serial present(xyz_grad[0:N3])
  {
    xyz_grad[3*m1+0] += valxt1*scalar;
    xyz_grad[3*m1+1] += valyt1*scalar;
    xyz_grad[3*m1+2] += valzt1*scalar;
    xyz_grad[3*n1+0] += valxt2*scalar;
    xyz_grad[3*n1+1] += valyt2*scalar;
    xyz_grad[3*n1+2] += valzt2*scalar;
  }

  return;
}

#if RED_DOUBLE
void reduce_3c2d122(int m1, int n1, int s1, int s2, int s3, int s4, int gs, FP2* norms1, FP2* norms2, FP2* dC, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, FP1** val5x, FP1** val6x, int N, int Naux, int imaxN, int imaxNa, int natoms, FP2 scalar, FP2* xyz_grad)
#else
void reduce_3c2d122(int m1, int n1, int s1, int s2, int s3, int s4, int gs, FP2* norms1, FP2* norms2, FP2* dC, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, FP1** val5x, FP1** val6x, int N, int Naux, int imaxN, int imaxNa, FP1 scalar, FP1* xyz_grad)
#endif
{
  int N2 = N*N;
  int N2a = N2*Naux;
  int N3 = 3*natoms;
  int gs3 = 3*gs;

  FP2 valxt1 = 0.; FP2 valyt1 = 0.; FP2 valzt1 = 0.;
  FP2 valxt2 = 0.; FP2 valyt2 = 0.; FP2 valzt2 = 0.;

#if USE_ACC
 #pragma acc parallel loop collapse(3) present(val1[0:imaxNa][0:gs],val2[0:imaxN][0:gs],val3[0:imaxN][0:gs],val4[0:imaxNa][0:gs],val5[0:imaxN][0:gs],val6[0:imaxN][0:gs],val1x[0:imaxNa][0:gs3],val2x[0:imaxN][0:gs3],val3x[0:imaxN][0:gs3],val4x[0:imaxNa][0:gs],val5x[0:imaxN][0:gs],val6x[0:imaxN][0:gs],norms1[0:Naux],norms2[0:N2],dC[0:N2a]) reduction(+:valxt1,valyt1,valzt1,valxt2,valyt2,valzt2)
#endif
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  for (int i3=s3;i3<s4;i3++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3; int ii3 = i3-s3;

    FP2 valx1 = 0.; FP2 valy1 = 0.; FP2 valz1 = 0.;
    FP2 valx2 = 0.; FP2 valy2 = 0.; FP2 valz2 = 0.;

   #pragma acc loop reduction(+:valx1,valy1,valz1,valx2,valy2,valz2)
    for (int j=0;j<gs;j++)
    {
      FP2 v1 = val1[ii1][j]; FP2 v2 = val2[ii2][j]; FP2 v3 = val3[ii3][j];
      FP2 v12 = v1*v2; FP2 v13 = v1*v3; FP2 v23 = v2*v3;

      valx2 += v12*val3x[ii3][3*j+0]; valy2 += v12*val3x[ii3][3*j+1]; valz2 += v12*val3x[ii3][3*j+2];
      valx2 += v13*val2x[ii2][3*j+0]; valy2 += v13*val2x[ii2][3*j+1]; valz2 += v13*val2x[ii2][3*j+2];
      valx1 += v23*val1x[ii1][3*j+0]; valy1 += v23*val1x[ii1][3*j+1]; valz1 += v23*val1x[ii1][3*j+2];

      FP2 v4 = val4[ii1][j]; FP2 v5 = val5[ii2][j]; FP2 v6 = val6[ii3][j];
      FP2 v45 = v4*v5; FP2 v46 = v4*v6; FP2 v56 = v5*v6;

      valx2 += v45*val6x[ii3][3*j+0]; valy2 += v45*val6x[ii3][3*j+1]; valz2 += v45*val6x[ii3][3*j+2];
      valx2 += v46*val5x[ii2][3*j+0]; valy2 += v46*val5x[ii2][3*j+1]; valz2 += v46*val5x[ii2][3*j+2];
      valx1 += v56*val4x[ii1][3*j+0]; valy1 += v56*val4x[ii1][3*j+1]; valz1 += v56*val4x[ii1][3*j+2];
    }

    FP2 nd1 = dC[i1*N2+i2*N+i3]*norms1[i1]*norms2[i2*N+i3];
    valx1 *= nd1; valy1 *= nd1; valz1 *= nd1;
    valx2 *= nd1; valy2 *= nd1; valz2 *= nd1;
    valxt1 += valx1; valyt1 += valy1; valzt1 += valz1;
    valxt2 += valx2; valyt2 += valy2; valzt2 += valz2;

    //printf("  i13: %2i %2i %2i valxyz: %10.6f %10.6f %10.6f  %10.6f %10.6f %10.6f  norm: %8.3f dC: %8.5f \n",i1,i2,i3,valx1,valy1,valz1,valx2,valy2,valz2,norms1[i1]*norms2[i2*N+i3],dC[i1*N2+i2*N+i3]);
  } //i1,i2,i3
 
 #pragma acc serial present(xyz_grad[0:N3])
  {
    xyz_grad[3*m1+0] += valxt1*scalar;
    xyz_grad[3*m1+1] += valyt1*scalar;
    xyz_grad[3*m1+2] += valzt1*scalar;
    xyz_grad[3*n1+0] += valxt2*scalar;
    xyz_grad[3*n1+1] += valyt2*scalar;
    xyz_grad[3*n1+2] += valzt2*scalar;
  }

  return;
}

#if RED_DOUBLE
void reduce_3c3d(int m1, int n1, int p1, int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP2* norms1, FP2* norms2, FP2* dC, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val7, FP1** val8, FP1** val9, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, FP1** val5x, FP1** val6x, FP1** val7x, FP1** val8x, FP1** val9x, int N, int Naux, int imaxN, int imaxNa, int natoms, FP2 scalar, FP2* xyz_grad)
#else
void reduce_3c3d(int m1, int n1, int p1, int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP2* norms1, FP2* norms2, FP2* dC, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val7, FP1** val8, FP1** val9, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, FP1** val5x, FP1** val6x, FP1** val7x, FP1** val8x, FP1** val9x, int N, int Naux, int imaxN, int imaxNa, FP1 scalar, FP1* xyz_grad)
#endif
{
  int N2 = N*N;
  int N2a = N2*Naux;
  int N3 = 3*natoms;
  int gs3 = 3*gs;

  FP2 valxt1 = 0.; FP2 valyt1 = 0.; FP2 valzt1 = 0.;
  FP2 valxt2 = 0.; FP2 valyt2 = 0.; FP2 valzt2 = 0.;
  FP2 valxt3 = 0.; FP2 valyt3 = 0.; FP2 valzt3 = 0.;

#if USE_ACC
 #pragma acc parallel loop collapse(3) present(val1[0:imaxNa][0:gs],val2[0:imaxN][0:gs],val3[0:imaxN][0:gs],val4[0:imaxNa][0:gs],val5[0:imaxN][0:gs],val6[0:imaxN][0:gs],val7[0:imaxNa][0:gs],val8[0:imaxN][0:gs],val9[0:imaxN][0:gs],val1x[0:imaxNa][0:gs3],val2x[0:imaxN][0:gs3],val3x[0:imaxN][0:gs3],val4x[0:imaxNa][0:gs3],val5x[0:imaxN][0:gs3],val6x[0:imaxN][0:gs3],val7x[0:imaxNa][0:gs3],val8x[0:imaxN][0:gs3],val9x[0:imaxN][0:gs3],norms1[0:Naux],norms2[0:N2],dC[0:N2a]) reduction(+:valxt1,valyt1,valzt1,valxt2,valyt2,valzt2,valxt3,valyt3,valzt3)
#endif
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  for (int i3=s5;i3<s6;i3++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3; int ii3 = i3-s5;

    FP2 valx1 = 0.; FP2 valy1 = 0.; FP2 valz1 = 0.;
    FP2 valx2 = 0.; FP2 valy2 = 0.; FP2 valz2 = 0.;
    FP2 valx3 = 0.; FP2 valy3 = 0.; FP2 valz3 = 0.;

   #pragma acc loop reduction(+:valx1,valy1,valz1,valx2,valy2,valz2,valx3,valy3,valz3)
    for (int j=0;j<gs;j++)
    {
      FP2 v1 = val1[ii1][j]; FP2 v2 = val2[ii2][j]; FP2 v3 = val3[ii3][j];
      FP2 v12 = v1*v2; FP2 v13 = v1*v3; FP2 v23 = v2*v3;
      valx3 += v12*val3x[ii3][3*j+0]; valy3 += v12*val3x[ii3][3*j+1]; valz3 += v12*val3x[ii3][3*j+2];
      valx2 += v13*val2x[ii2][3*j+0]; valy2 += v13*val2x[ii2][3*j+1]; valz2 += v13*val2x[ii2][3*j+2];
      valx1 += v23*val1x[ii1][3*j+0]; valy1 += v23*val1x[ii1][3*j+1]; valz1 += v23*val1x[ii1][3*j+2];

      FP2 v4 = val4[ii1][j]; FP2 v5 = val5[ii2][j]; FP2 v6 = val6[ii3][j];
      FP2 v45 = v4*v5; FP2 v46 = v4*v6; FP2 v56 = v5*v6;
      valx3 += v45*val6x[ii3][3*j+0]; valy3 += v45*val6x[ii3][3*j+1]; valz3 += v45*val6x[ii3][3*j+2];
      valx2 += v46*val5x[ii2][3*j+0]; valy2 += v46*val5x[ii2][3*j+1]; valz2 += v46*val5x[ii2][3*j+2];
      valx1 += v56*val4x[ii1][3*j+0]; valy1 += v56*val4x[ii1][3*j+1]; valz1 += v56*val4x[ii1][3*j+2];

      FP2 v7 = val7[ii1][j]; FP2 v8 = val8[ii2][j]; FP2 v9 = val9[ii3][j];
      FP2 v78 = v7*v8; FP2 v79 = v7*v9; FP2 v89 = v8*v9;
      valx3 += v78*val9x[ii3][3*j+0]; valy3 += v78*val9x[ii3][3*j+1]; valz3 += v78*val9x[ii3][3*j+2];
      valx2 += v79*val8x[ii2][3*j+0]; valy2 += v79*val8x[ii2][3*j+1]; valz2 += v79*val8x[ii2][3*j+2];
      valx1 += v89*val7x[ii1][3*j+0]; valy1 += v89*val7x[ii1][3*j+1]; valz1 += v89*val7x[ii1][3*j+2];
    }

    FP2 nd1 = dC[i1*N2+i2*N+i3]*norms1[i1]*norms2[i2*N+i3];
    valxt1 += valx1*nd1; valyt1 += valy1*nd1; valzt1 += valz1*nd1;
    valxt2 += valx2*nd1; valyt2 += valy2*nd1; valzt2 += valz2*nd1;
    valxt3 += valx3*nd1; valyt3 += valy3*nd1; valzt3 += valz3*nd1;

    //if (nd1>0.) printf("  i13: %2i %2i %2i valxyz: %10.6f %10.6f %10.6f  %10.6f %10.6f %10.6f  %10.6f %10.6f %10.6f  nd1: %8.5f \n",i1,i2,i3,valx1,valy1,valz1,valx2,valy2,valz2,valx3,valy3,valz3,norms1[i1]*norms2[i2*N+i3]*dC[i1*N2+i2*N+i3]);
  } //i1,i2,i3
 
 #pragma acc serial present(xyz_grad[0:N3])
  {
    xyz_grad[3*m1+0] += valxt1*scalar;
    xyz_grad[3*m1+1] += valyt1*scalar;
    xyz_grad[3*m1+2] += valzt1*scalar;
    xyz_grad[3*n1+0] += valxt2*scalar;
    xyz_grad[3*n1+1] += valyt2*scalar;
    xyz_grad[3*n1+2] += valzt2*scalar;
    xyz_grad[3*p1+0] += valxt3*scalar;
    xyz_grad[3*p1+1] += valyt3*scalar;
    xyz_grad[3*p1+2] += valzt3*scalar;
  }

  return;
}






#if RED_DOUBLE
void reduce_3c1(int s1, int s2, int s3, int s4, int gs, FP1** val1, FP1** val2, FP1** val3, FP1* valt1, int N, int Naux, int imaxN, int imaxNa, FP2* C)
#else
void reduce_3c1(int s1, int s2, int s3, int s4, int gs, FP1** val1, FP1** val2, FP1** val3, FP1* valt1, int N, int Naux, int imaxN, int imaxNa, FP1* C)
#endif
{
  int N2 = N*N;
  int N2a = N2*Naux;

  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3;

   #pragma acc parallel loop present(val1[0:imaxNa][0:gs],val2[0:imaxN][0:gs],valt1[0:gs])
    for (int j=0;j<gs;j++)
      valt1[j] = val1[ii1][j] * val2[ii2][j];

   #if USE_ACC
    #pragma acc parallel loop present(valt1[0:gs],val3[0:imaxN][0:gs],C[0:N2a]) 
   #endif
    for (int i3=s3;i3<s4;i3++)
    {
      int ii3 = i3-s3;

     #if RED_DOUBLE
      FP2 val = 0.;
     #else
      FP1 val = 0.f;
     #endif

     #pragma acc loop reduction(+:val)
      for (int j=0;j<gs;j++)
        val += valt1[j] * val3[ii3][j];
 
      C[i1*N2+i2*N+i3] = val;

      //printf(" m: %i  i1/2/3: %2i %2i %2i val: %8.5f \n",m,i1,i2,i3,val);
    } //i3
  } //i1,i2
  return;
}

#if RED_DOUBLE
void reduce_3c1b(int s1, int s2, int s3, int s4, int gs, FP1** val1, FP1** val2, FP1** val3, int N, int Naux, int imaxN, int imaxNa, FP2* C)
#else
void reduce_3c1b(int s1, int s2, int s3, int s4, int gs, FP1** val1, FP1** val2, FP1** val3, int N, int Naux, int imaxN, int imaxNa, FP1* C)
#endif
{
  int N2 = N*N;
  int N2a = N2*Naux;

#if USE_ACC
 #pragma acc parallel loop collapse(3) present(val1[0:imaxNa][0:gs],val2[0:imaxN][0:gs],val3[0:imaxN][0:gs],C[0:N2a]) 
#endif
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  for (int i3=s3;i3<s4;i3++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3; int ii3 = i3-s3;

   #if RED_DOUBLE
    FP2 val = 0.;
   #else
    FP1 val = 0.f;
   #endif

   #pragma acc loop reduction(+:val)
    for (int j=0;j<gs;j++)
      val += val1[ii1][j] * val2[ii2][j] * val3[ii3][j];
 
    C[i1*N2+i2*N+i3] = val;

    //printf(" m: %i  i1/2/3: %2i %2i %2i val: %8.5f \n",m,i1,i2,i3,val);
  } //i1,i2,i3

  return;
}

#if RED_DOUBLE
void reduce_3c2(int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1* valt1, FP1* valt2, int N, int Naux, int imaxN, int imaxNa, FP2* C)
#else
void reduce_3c2(int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1* valt1, FP1* valt2, int N, int Naux, int imaxN, int imaxNa, FP1* C)
#endif
{
  int N2 = N*N;
  int N2a = N2*Naux;

  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3;
  #if USE_ACC
   #pragma acc parallel loop present(val1[0:imaxNa][0:gs],val2[0:imaxN][0:gs],valt1[0:gs]) 
  #endif
    for (int j=0;j<gs;j++)
      valt1[j] = val1[ii1][j] * val2[ii2][j];
  #if USE_ACC
   #pragma acc parallel loop present(val4[0:imaxNa][0:gs],val5[0:imaxN][0:gs],valt2[0:gs]) 
  #endif
    for (int j=0;j<gs;j++)
      valt2[j] = val4[ii1][j] * val5[ii2][j];

  #if USE_ACC
   #pragma acc parallel loop present(val3[0:imaxN][0:gs],val6[0:imaxN][0:gs],valt1[0:gs],valt2[0:gs],C[0:N2a]) 
  #endif
    for (int i3=s5;i3<s6;i3++)
    {
      int ii3 = i3-s5;

     #if RED_DOUBLE
      FP2 val = 0.;
     #else
      FP1 val = 0.f;
     #endif

    #if USE_ACC
     #pragma acc loop reduction(+:val)
    #endif
      for (int j=0;j<gs;j++)
        val += valt1[j] * val3[ii3][j];
    #if USE_ACC
     #pragma acc loop reduction(+:val)
    #endif
      for (int j=0;j<gs;j++)
        val += valt2[j] * val6[ii3][j];
 
      C[i1*N2+i3*N+i2] = val;
    }
  }

  return;
}

#if RED_DOUBLE
void reduce_3c2b(int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, int N, int Naux, int imaxN, int imaxNa, FP2* C)
#else
void reduce_3c2b(int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, int N, int Naux, int imaxN, int imaxNa, FP1* C)
#endif
{
  int N2 = N*N;
  int N2a = N2*Naux;

 #if USE_ACC
  #pragma acc parallel loop collapse(3) present(val1[0:imaxNa][0:gs],val2[0:imaxN][0:gs],val3[0:imaxN][0:gs],val4[0:imaxNa][0:gs],val5[0:imaxN][0:gs],val6[0:imaxN][0:gs],C[0:N2a])
 #endif
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  for (int i3=s5;i3<s6;i3++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3; int ii3 = i3-s5;

   #if RED_DOUBLE
    FP2 val = 0.;
   #else
    FP1 val = 0.f;
   #endif

   #pragma acc loop reduction(+:val)
    for (int j=0;j<gs;j++)
      val += val1[ii1][j] * val2[ii2][j] * val3[ii3][j];
   #pragma acc loop reduction(+:val)
    for (int j=0;j<gs;j++)
      val += val4[ii1][j] * val5[ii2][j] * val6[ii3][j];
 
    C[i1*N2+i3*N+i2] = val; 
  }

  return;
}

#if RED_DOUBLE
void reduce_3c3(int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val7, FP1** val8, FP1** val9, FP1* valt1, FP1* valt2, FP1* valt3, int N, int Naux, int imaxN, int imaxNa, FP2* C)
#else
void reduce_3c3(int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val7, FP1** val8, FP1** val9, FP1* valt1, FP1* valt2, FP1* valt3, int N, int Naux, int imaxN, int imaxNa, FP1* C)
#endif
{
  int N2 = N*N;
  int N2a = N2*Naux;

  for (int i1=s1;i1<s2;i1++)
 #pragma acc loop
  for (int i2=s3;i2<s4;i2++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3; 
   #if USE_ACC
    #pragma acc parallel loop present(val1[0:imaxNa][0:gs],val2[0:imaxN][0:gs],valt1[0:gs]) 
   #endif
    for (int j=0;j<gs;j++)
      valt1[j] = val1[ii1][j] * val2[ii2][j];
   #if USE_ACC
    #pragma acc parallel loop present(val4[0:imaxNa][0:gs],val5[0:imaxN][0:gs],valt2[0:gs]) 
   #endif
    for (int j=0;j<gs;j++)
      valt2[j] = val4[ii1][j] * val5[ii2][j];
   #if USE_ACC
    #pragma acc parallel loop present(val7[0:imaxNa][0:gs],val8[0:imaxN][0:gs],valt3[0:gs]) 
   #endif
    for (int j=0;j<gs;j++)
      valt3[j] = val7[ii1][j] * val8[ii2][j];

   #if USE_ACC
    #pragma acc parallel loop present(val3[0:imaxN][0:gs],val6[0:imaxN][0:gs],val9[0:imaxN][0:gs],valt1[0:gs],valt2[0:gs],valt3[0:gs],C[0:N2a]) 
   #endif
    for (int i3=s5;i3<s6;i3++)
    {
      int ii3 = i3-s5;

     #if RED_DOUBLE
      FP2 val = 0.;
     #else
      FP1 val = 0.f;
     #endif

     #pragma acc loop reduction(+:val)
      for (int j=0;j<gs;j++)
        val += valt1[j] * val3[ii3][j];
     #pragma acc loop reduction(+:val)
      for (int j=0;j<gs;j++)
        val += valt2[j] * val6[ii3][j];
     #pragma acc loop reduction(+:val)
      for (int j=0;j<gs;j++)
        val += valt3[j] * val9[ii3][j];

      C[i1*N2+i2*N+i3] = val;

     } //i3
   } //loop i1,i2
  return;
}

#if RED_DOUBLE
void reduce_3c3b(int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val7, FP1** val8, FP1** val9, int N, int Naux, int imaxN, int imaxNa, FP2* C)
#else
void reduce_3c3b(int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val7, FP1** val8, FP1** val9, int N, int Naux, int imaxN, int imaxNa, FP1* C)
#endif
{
  int N2 = N*N;
  int N2a = N2*Naux;

 #if USE_ACC
  #pragma acc parallel loop collapse(3) present(val1[0:imaxNa][0:gs],val2[0:imaxN][0:gs],val3[0:imaxN][0:gs],val4[0:imaxNa][0:gs],val5[0:imaxN][0:gs],val6[0:imaxN][0:gs],val7[0:imaxNa][0:gs],val8[0:imaxN][0:gs],val9[0:imaxN][0:gs],C[0:N2a])
 #endif
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  for (int i3=s5;i3<s6;i3++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3; int ii3 = i3-s5;

   #if RED_DOUBLE
    FP2 val = 0.;
   #else
    FP1 val = 0.f;
   #endif

   #pragma acc loop reduction(+:val)
    for (int j=0;j<gs;j++)
      val += val1[ii1][j] * val2[ii2][j] * val3[ii3][j];
   #pragma acc loop reduction(+:val)
    for (int j=0;j<gs;j++)
      val += val4[ii1][j] * val5[ii2][j] * val6[ii3][j];
   #pragma acc loop reduction(+:val)
    for (int j=0;j<gs;j++)
      val += val7[ii1][j] * val8[ii2][j] * val9[ii3][j];
    
    C[i1*N2+i2*N+i3] = val;
  }

  return;
}


#if RED_DOUBLE
void reduce_2c2v(int p1, int s1, int s2, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, int iN, int N, int nc, FP2 scalar, FP2* V)
#else
void reduce_2c2v(int p1, int s1, int s2, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, int iN, int N, int nc, FP2 scalar, FP1* V)
#endif
{
  int N2 = N*N;

  FP2 valt = 0.;

 #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs],val4[0:iN][0:gs],norms[0:N2],dpq[0:N2]) reduction(+:valt)
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s1;i2<s2;i2++)
    {
      FP2 val = 0.;
     #pragma acc loop reduction(+:val)
      for (int j=0;j<gs;j++)
      {
        int ii1 = i1-s1;
        int ii2 = i2-s1;

        FP2 v3 = val3[ii2][j]; FP2 v4 = val4[ii2][j];
        val += v3*val1[ii1][j];
        val += v4*val2[ii1][j]; 
      }

      FP2 nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
      valt += val*nd1;
    }
  }

 #pragma acc serial present(V[0:nc])
  {
    V[p1] += scalar*valt;
  }

  return;
}

#if RED_DOUBLE
void reduce_2c2vd(int p1, int s1, int s2, int gs, FP2* norms, FP2* dpq, FP1** val1x, FP1** val2x, FP1** val3, FP1** val4, int iN, int N, int nc, FP2 scalar, FP2* dV)
#else
void reduce_2c2vd(int p1, int s1, int s2, int gs, FP2* norms, FP2* dpq, FP1** val1x, FP1** val2x, FP1** val3, FP1** val4, int iN, int N, int nc, FP2 scalar, FP1* dV)
#endif
{
  int N2 = N*N;
  int gs3 = 3*gs;

  FP2 valxt = 0.; FP2 valyt = 0.; FP2 valzt = 0.;

 #pragma acc parallel loop collapse(2) present(val1x[0:iN][0:gs3],val2x[0:iN][0:gs3],val3[0:iN][0:gs],val4[0:iN][0:gs],norms[0:N2],dpq[0:N2]) reduction(+:valxt,valyt,valzt)
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s1;i2<s2;i2++)
    {
      FP2 valx = 0.; FP2 valy = 0.; FP2 valz = 0.;
     #pragma acc loop reduction(+:valx,valy,valz)
      for (int j=0;j<gs;j++)
      {
        int ii1 = i1-s1;
        int ii2 = i2-s1;

        FP2 v3 = val3[ii2][j]; FP2 v4 = val4[ii2][j];
        valx += v3*val1x[ii1][3*j+0]; valy += v3*val1x[ii1][3*j+1]; valz += v3*val1x[ii1][3*j+2];
        valx += v4*val2x[ii1][3*j+0]; valy += v4*val2x[ii1][3*j+1]; valz += v4*val2x[ii1][3*j+2];
      }
      FP2 nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
      valx *= nd1; valy *= nd1; valz *= nd1;

      valxt += valx; valyt += valy; valzt += valz;

      // if (nd1>0.)
      //  printf("  i1/2(mn): %2i %2i (%i%i) valxyz: %9.5f %9.5f %9.5f  \n",i1,i2,m1,n1,valx+valx2,valy+valy2,valz+valz2);
    }
  }

 #pragma acc serial present(dV[0:3*nc])
  {
    dV[3*p1+0] += scalar*valxt;
    dV[3*p1+1] += scalar*valyt;
    dV[3*p1+2] += scalar*valzt;
  }

  return;
}

#if RED_DOUBLE
void reduce_2c3v(int p1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, int iN, int N, int nc, FP2 scalar, FP2* V)
#else
void reduce_2c3v(int p1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, int iN, int N, int nc, FP2 scalar, FP1* V)
#endif
{
  int N2 = N*N;

  FP2 valt = 0.;

 #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs],val4[0:iN][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs],norms[0:N2],dpq[0:N2]) reduction(+:valt)
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s3;i2<s4;i2++)
    {
      FP2 val = 0.;
     #pragma acc loop reduction(+:val)
      for (int j=0;j<gs;j++)
      {
        int ii1 = i1-s1;
        int ii2 = i2-s3;

        FP2 v1 = val1[ii1][j]; FP2 v2 = val2[ii1][j]; FP2 v3 = val3[ii1][j];
        val += v1*val4[ii2][j];
        val += v2*val5[ii2][j]; 
        val += v3*val6[ii2][j]; 
      }

      FP2 nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
      valt += val*nd1;

      printf("  i12: %i-%i  val: %8.5f  norm/dpq: %8.5f %8.5f  \n",i1,i2,val*nd1,norms[i1*N+i2],dpq[i1*N+i2]);
    }
  }

 #pragma acc serial present(V[0:nc])
  {
    V[p1] += scalar*valt;
  }

  return;
}

#if RED_DOUBLE
void reduce_2c3vd(int p1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4, FP1** val5, FP1** val6, int iN, int N, int nc, FP2 scalar, FP2* dV)
#else
void reduce_2c3vd(int p1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4, FP1** val5, FP1** val6, int iN, int N, int nc, FP2 scalar, FP1* dV)
#endif
{
  int N2 = N*N;
  int gs3 = 3*gs;

  FP2 valxt = 0.; FP2 valyt = 0.; FP2 valzt = 0.;

 #pragma acc parallel loop collapse(2) present(val1x[0:iN][0:gs3],val2x[0:iN][0:gs3],val3x[0:iN][0:gs],val4[0:iN][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs],norms[0:N2],dpq[0:N2]) reduction(+:valxt,valyt,valzt)
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s3;i2<s4;i2++)
    {
      FP2 valx = 0.; FP2 valy = 0.; FP2 valz = 0.;
     #pragma acc loop reduction(+:valx,valy,valz)
      for (int j=0;j<gs;j++)
      {
        int ii1 = i1-s1;
        int ii2 = i2-s3;

        FP2 v4 = val4[ii2][j]; FP2 v5 = val5[ii2][j]; FP2 v6 = val6[ii2][j];
        valx += v4*val1x[ii1][3*j+0]; valy += v4*val1x[ii1][3*j+1]; valz += v4*val1x[ii1][3*j+2];
        valx += v5*val2x[ii1][3*j+0]; valy += v5*val2x[ii1][3*j+1]; valz += v5*val2x[ii1][3*j+2];
        valx += v6*val3x[ii1][3*j+0]; valy += v6*val3x[ii1][3*j+1]; valz += v6*val3x[ii1][3*j+2];
      }
      FP2 nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
      valx *= nd1; valy *= nd1; valz *= nd1;

      valxt += valx; valyt += valy; valzt += valz;

      // if (nd1>0.)
      //  printf("  i1/2(mn): %2i %2i (%i%i) valxyz: %9.5f %9.5f %9.5f  \n",i1,i2,m1,n1,valx+valx2,valy+valy2,valz+valz2);
    }
  }

 #pragma acc serial present(dV[0:3*nc])
  {
    dV[3*p1+0] += scalar*valxt;
    dV[3*p1+1] += scalar*valyt;
    dV[3*p1+2] += scalar*valzt;
  }

  return;
}
