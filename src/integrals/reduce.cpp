#include "integrals.h"


void reduce_4c_1d(int s1, int s2, int s3, int s4, int gs, int gsp, int M, float* grid1, float* grid2, float** val1a, float** val2a, float** val3a, float** val4a, float* wt1, float* wt2, float* gt)
{
 //2221
  int iN = M;
  int M2 = M*M;
  int M3 = M2*M;
  int gs6 = 6*gs;
  int gsp6 = 6*gsp;

  #pragma acc parallel loop present(gt[0:M2*M2],grid1[0:gs6],grid2[0:gsp6],val1a[0:iN][0:gs],val2a[0:iN][0:gs],val3a[0:iN][0:gsp],val4a[0:iN][0:gsp],wt1[0:gs],wt2[0:gsp])
  for (int i1=s3;i1<s4;i1++)
  #pragma acc loop independent
  for (int i2=s3;i2<=i1;i2++)
  #pragma acc loop independent
  for (int i3=s3;i3<s4;i3++)
  #pragma acc loop independent
  for (int i4=s1;i4<s2;i4++)
  {
    int ii1 = i1-s3; int ii2 = i2-s3;
    int ii3 = i3-s3; int ii4 = i4-s1;
    float* val1n = val1a[ii1]; float* val2n = val2a[ii2];
    float* val3n = val3a[ii3]; float* val4n = val4a[ii4];

    float den1 = 0.f;
    double v1 = 0.;
    #pragma acc loop reduction(+:v1)
    for (int j=0;j<gs;j++) //cannot use collapse: total index too high
    {
      float x1 = grid1[6*j+0]; float y1 = grid1[6*j+1]; float z1 = grid1[6*j+2];

      float vj = val1n[j]*val2n[j]*wt1[j];
      float v2 = 0.f;

      if (fabs(vj)>1.e-10)
      #pragma acc loop reduction(+:v2)
      for (int k=0;k<gsp;k++)
      {
        float x2 = grid2[6*k+0]; float y2 = grid2[6*k+1]; float z2 = grid2[6*k+2];
        float x12 = x1-x2; float y12 = y1-y2; float z12 = z1-z2;

        float r1 = sqrtf(x12*x12+y12*y12+z12*z12)+1.e-12;

        v2 += val3n[k]*val4n[k]*wt2[k]/r1;
      }

      v1 += v2*vj;
    }

    gt[ii1*M3+ii2*M2+ii3*M+ii4] += v1;
    //printf("  v1(%i %i %i %i): %8.5f \n",i1,i2,i3,i4,v1);
  }

  return;
}

void reduce_4c_1c(int s1, int s2, int s3, int s4, int gs, int gsp, int M, float* grid1, float* grid2, float** val1a, float** val2a, float** val3a, float** val4a, float* wt1, float* wt2, float* gt)
{
 //1221
  int iN = M;
  int M2 = M*M;
  int M3 = M2*M;
  int gs6 = 6*gs;
  int gsp6 = 6*gsp;

  #pragma acc parallel loop present(gt[0:M2*M2],grid1[0:gs6],grid2[0:gsp6],val1a[0:iN][0:gs],val2a[0:iN][0:gs],val3a[0:iN][0:gsp],val4a[0:iN][0:gsp],wt1[0:gs],wt2[0:gsp])
  for (int i1=s1;i1<s2;i1++)
  #pragma acc loop independent
  for (int i2=s3;i2<s4;i2++)
  #pragma acc loop independent
  for (int i3=s3;i3<s4;i3++)
  #pragma acc loop independent
  for (int i4=s1;i4<s2;i4++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3;
    int ii3 = i3-s3; int ii4 = i4-s1;
    float* val1n = val1a[ii1]; float* val2n = val2a[ii2];
    float* val3n = val3a[ii3]; float* val4n = val4a[ii4];

    float den1 = 0.f;
    double v1 = 0.;
    #pragma acc loop reduction(+:v1)
    for (int j=0;j<gs;j++) //cannot use collapse: total index too high
    {
      float x1 = grid1[6*j+0]; float y1 = grid1[6*j+1]; float z1 = grid1[6*j+2];

      float vj = val1n[j]*val2n[j]*wt1[j];
      float v2 = 0.f;

      if (fabs(vj)>1.e-10)
      #pragma acc loop reduction(+:v2)
      for (int k=0;k<gsp;k++)
      {
        float x2 = grid2[6*k+0]; float y2 = grid2[6*k+1]; float z2 = grid2[6*k+2];
        float x12 = x1-x2; float y12 = y1-y2; float z12 = z1-z2;

        float r1 = sqrtf(x12*x12+y12*y12+z12*z12)+1.e-12;

        v2 += val3n[k]*val4n[k]*wt2[k]/r1;
      }

      v1 += v2*vj;
    }

    gt[ii1*M3+ii2*M2+ii3*M+ii4] += v1;
    //printf("  v1(%i %i %i %i): %8.5f \n",i1,i2,i3,i4,v1);
  }

  return;
}

void reduce_4c_1b(int s1, int s2, int s3, int s4, int gs, int gsp, int M, float* grid1, float* grid2, float** val1a, float** val2a, float** val3a, float** val4a, float* wt1, float* wt2, float* gt)
{
 //1121
  int iN = M;
  int M2 = M*M;
  int M3 = M2*M;
  int gs6 = 6*gs;
  int gsp6 = 6*gsp;

  #pragma acc parallel loop present(gt[0:M2*M2],grid1[0:gs6],grid2[0:gsp6],val1a[0:iN][0:gs],val2a[0:iN][0:gs],val3a[0:iN][0:gsp],val4a[0:iN][0:gsp],wt1[0:gs],wt2[0:gsp])
  for (int i1=s1;i1<s2;i1++)
  #pragma acc loop independent
  for (int i2=s1;i2<=i1;i2++)
  #pragma acc loop independent
  for (int i3=s3;i3<s4;i3++)
  #pragma acc loop independent
  for (int i4=s1;i4<s2;i4++)
  {
    int ii1 = i1-s1; int ii2 = i2-s1;
    int ii3 = i3-s3; int ii4 = i4-s1;
    float* val1n = val1a[ii1]; float* val2n = val2a[ii2];
    float* val3n = val3a[ii3]; float* val4n = val4a[ii4];

    float den1 = 0.f;
    double v1 = 0.;
    #pragma acc loop reduction(+:v1)
    for (int j=0;j<gs;j++) //cannot use collapse: total index too high
    {
      float x1 = grid1[6*j+0]; float y1 = grid1[6*j+1]; float z1 = grid1[6*j+2];

      float vj = val1n[j]*val2n[j]*wt1[j];
      float v2 = 0.f;

      if (fabs(vj)>1.e-10)
      #pragma acc loop reduction(+:v2)
      for (int k=0;k<gsp;k++)
      {
        float x2 = grid2[6*k+0]; float y2 = grid2[6*k+1]; float z2 = grid2[6*k+2];
        float x12 = x1-x2; float y12 = y1-y2; float z12 = z1-z2;

        float r1 = sqrtf(x12*x12+y12*y12+z12*z12)+1.e-12;

        v2 += val3n[k]*val4n[k]*wt2[k]/r1;
      }

      v1 += v2*vj;
    }

    gt[ii1*M3+ii2*M2+ii3*M+ii4] += v1;
    //printf("  v1(%i %i %i %i): %8.5f \n",i1,i2,i3,i4,v1);
  }

  return;
}

void reduce_4c_1(int s1, int s2, int s3, int s4, int gs, int gsp, int M, float* grid1, float* grid2, float** val1a, float** val2a, float** val3a, float** val4a, float* wt1, float* wt2, float* gt)
{
 //1122
  int iN = M;
  int M2 = M*M;
  int M3 = M2*M;
  int gs6 = 6*gs;
  int gsp6 = 6*gsp;

 //collapse i1,i3?

   //alternatively, could compute potential first
  #pragma acc parallel loop present(gt[0:M2*M2],grid1[0:gs6],grid2[0:gsp6],val1a[0:iN][0:gs],val2a[0:iN][0:gs],val3a[0:iN][0:gsp],val4a[0:iN][0:gsp],wt1[0:gs],wt2[0:gsp])
  for (int i1=s1;i1<s2;i1++)
  #pragma acc loop independent
  for (int i2=s1;i2<=i1;i2++)
  #pragma acc loop independent
  for (int i3=s3;i3<s4;i3++)
  #pragma acc loop independent
  for (int i4=s3;i4<=i3;i4++)
  {
    int ii1 = i1-s1; int ii2 = i2-s1;
    int ii3 = i3-s3; int ii4 = i4-s3;
    float* val1n = val1a[ii1]; float* val2n = val2a[ii2];
    float* val3n = val3a[ii3]; float* val4n = val4a[ii4];

    double v1 = 0.;
    #pragma acc loop reduction(+:v1)
    for (int j=0;j<gs;j++) //cannot use collapse: total index too high
    {
      float x1 = grid1[6*j+0]; float y1 = grid1[6*j+1]; float z1 = grid1[6*j+2];

      float vj = val1n[j]*val2n[j]*wt1[j];
      float v2 = 0.f;

      if (fabs(vj)>1.e-10)
      #pragma acc loop reduction(+:v2)
      for (int k=0;k<gsp;k++)
      {
        float x2 = grid2[6*k+0]; float y2 = grid2[6*k+1]; float z2 = grid2[6*k+2];
        float x12 = x1-x2; float y12 = y1-y2; float z12 = z1-z2;

        float r1 = sqrtf(x12*x12+y12*y12+z12*z12)+1.e-12;

        v2 += val3n[k]*val4n[k]*wt2[k]/r1;
      }

      v1 += v2*vj;
    }

    gt[ii1*M3+ii2*M2+ii3*M+ii4] += v1;
    //printf("  v1(%i %i %i %i): %8.5f \n",i1,i2,i3,i4,v1);
  }

  return;
}

//fully double precision
void reduce_2c1(int s1, int s2, int gs, double** val1, double** val3, int iN, int N, double* An)
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
      double val = 0.;

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
    double val = 0.;
    for (int j=0;j<gs;j++)
      val += val1[ii1][j]*val3[ii2][j];
    An[i1*N+i2] = val;
  }
#endif

  return;
}

//fully double precision
void reduce_2c1(int s1, int s2, int s3, int s4, int gs, double** val1, double** val3, int iN, int N, double* An)
{
  int N2 = N*N;

#if USE_ACC
 #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs],val3[0:iN][0:gs],An[0:N2])
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s3;i2<s4;i2++)
    {
      int ii1 = i1-s1;
      int ii2 = i2-s3;
      double val = 0.;

     #pragma acc loop reduction(+:val)
      for (int j=0;j<gs;j++)
        val += val1[ii1][j]*val3[ii2][j];
   
      An[i1*N+i2] = val;
    }
  } //MM reduction
  #pragma acc wait

#else
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3;
    double val = 0.;
    for (int j=0;j<gs;j++)
      val += val1[ii1][j]*val3[ii2][j];
    An[i1*N+i2] = val;
  }
#endif

  return;
}

#if RED_DOUBLE
void reduce_2c1(int s1, int s2, int gs, float** val1, float** val3, int iN, int N, double* An)
#else
void reduce_2c1(int s1, int s2, int gs, float** val1, float** val3, int iN, int N, float* An)
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
      double val = 0.;
     #else
      float val = 0.f;
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
    float val = 0.f;
    for (int j=0;j<gs;j++)
      val += val1[ii1][j]*val3[ii2][j];
    An[i1*N+i2] = val;
  }
#endif

  return;
}

//fully double precision
void reduce_2c2(int s1, int s2, int s3, int s4, int gs, double** val1, double** val2, double** val3, double** val4, int iN, int N, double* An)
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
    
      double val = 0.;

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
    double val = 0.;
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
void reduce_2c2(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, float** val4, int iN, int N, double* An)
#else
void reduce_2c2(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, float** val4, int iN, int N, float* An)
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
      double val = 0.;
     #else
      float val = 0.f;
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
    float val = 0.f;
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
void reduce_2c3(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, int iN, int N, double* An)
#else
void reduce_2c3(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, int iN, int N, float* An)
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
      double val = 0.;
     #else
      float val = 0.f;
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
    float val = 0.f;
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
void reduce_2c1d(int m1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val3, float** val1x, float** val3x, int iN, int N, int natoms, double scalar, double* xyz_grad)
#else
void reduce_2c1d(int m1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val3, float** val1x, float** val3x, int iN, int N, int natoms, double scalar, float* xyz_g)
#endif
{
  int N2 = N*N;
  int N3 = 3*natoms;
  int gs3 = 3*gs;

  double valxt = 0.;
  double valyt = 0.;
  double valzt = 0.;

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
      double valx, valy, valz;
      valx = valy = valz = 0.;
     #else
      float valx, valy, valz;
      valx = valy = valz = 0.f;
     #endif

     #pragma acc loop reduction(+:valx,valy,valz)
      for (int j=0;j<gs;j++)
      {
        float v1 = val1[ii1][j];
        valx += v1*val3x[ii2][3*j+0]; valy += v1*val3x[ii2][3*j+1]; valz += v1*val3x[ii2][3*j+2];

        float v3 = val3[ii2][j];
        valx += v3*val1x[ii1][3*j+0]; valy += v3*val1x[ii1][3*j+1]; valz += v3*val1x[ii1][3*j+2];
      }

      double nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
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
void reduce_2c1ds(int m1, int n1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3x, float** val4x, int iN, int N, int natoms, double scalar, double* xyz_grad)
#else
void reduce_2c1ds(int m1, int n1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3x, float** val4x, int iN, int N, int natoms, double scalar, float* xyz_g)
#endif
{
  int N2 = N*N;
  int N3 = 3*natoms;
  int gs3 = 3*gs;

  double valxt = 0.;
  double valyt = 0.;
  double valzt = 0.;

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
      double valx, valy, valz;
      valx = valy = valz = 0.;
     #else
      float valx, valy, valz;
      valx = valy = valz = 0.f;
     #endif

     #pragma acc loop reduction(+:valx,valy,valz)
      for (int j=0;j<gs;j++)
      {
        float v1 = val1[ii1][j];
        valx += v1*val3x[ii2][3*j+0];
        valy += v1*val3x[ii2][3*j+1];
        valz += v1*val3x[ii2][3*j+2];

        float v2 = val2[ii1][j];
        valx += v2*val4x[ii2][3*j+0];
        valy += v2*val4x[ii2][3*j+1];
        valz += v2*val4x[ii2][3*j+2];
      }

      double nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
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
void reduce_2c2d(int m1, int n1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, float** val1x, float** val2x, float** val3x, float** val4x, int iN, int N, int natoms, double scalar, double* xyz_grad)
#else
void reduce_2c2d(int m1, int n1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, float** val1x, float** val2x, float** val3x, float** val4x, int iN, int N, int natoms, double scalar, float* xyz_grad)
#endif
{
  int N2 = N*N;
  int gs3 = 3*gs;

  double valxt1 = 0.; double valyt1 = 0.; double valzt1 = 0.;
  double valxt2 = 0.; double valyt2 = 0.; double valzt2 = 0.;

 #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs],val4[0:iN][0:gs],norms[0:N2],val1x[0:iN][0:gs3],val2x[0:iN][0:gs3],val3x[0:iN][0:gs3],val4x[0:iN][0:gs3],dpq[0:N2]) reduction(+:valxt1,valxt2,valyt1,valyt2,valzt1,valzt2)
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s3;i2<s4;i2++)
    {
      double valx1 = 0.; double valy1 = 0.; double valz1 = 0.;
      double valx2 = 0.; double valy2 = 0.; double valz2 = 0.;
     #pragma acc loop reduction(+:valx1,valx2,valy1,valy2,valz1,valz2)
      for (int j=0;j<gs;j++)
      {
        int ii1 = i1-s1;
        int ii2 = i2-s3;

        double v1 = val1[ii1][j]; double v2 = val2[ii1][j];
        valx2 += v1*val3x[ii2][3*j+0]; valy2 += v1*val3x[ii2][3*j+1]; valz2 += v1*val3x[ii2][3*j+2];
        valx2 += v2*val4x[ii2][3*j+0]; valy2 += v2*val4x[ii2][3*j+1]; valz2 += v2*val4x[ii2][3*j+2];

        double v3 = val3[ii2][j]; double v4 = val4[ii2][j];
        valx1 += v3*val1x[ii1][3*j+0]; valy1 += v3*val1x[ii1][3*j+1]; valz1 += v3*val1x[ii1][3*j+2];
        valx1 += v4*val2x[ii1][3*j+0]; valy1 += v4*val2x[ii1][3*j+1]; valz1 += v4*val2x[ii1][3*j+2];
      }
      double nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
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
  double vxa = (valxt1 + valxt2)*0.5;
  double vya = (valyt1 + valyt2)*0.5;
  double vza = (valzt1 + valzt2)*0.5;
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
void reduce_2c2dh(int m1, int n1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3x, float** val4x, int iN, int N, int natoms, double scalar, double* xyz_grad)
#else
void reduce_2c2dh(int m1, int n1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3x, float** val4x, int iN, int N, int natoms, double scalar, float* xyz_grad)
#endif
{
  int N2 = N*N;
  int gs3 = 3*gs;

  double valxt = 0.; double valyt = 0.; double valzt = 0.;

 #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs],val2[0:iN][0:gs],val3x[0:iN][0:gs3],val4x[0:iN][0:gs3],dpq[0:N2],norms[0:N2]) reduction(+:valxt,valyt,valzt)
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s3;i2<s4;i2++)
    {
      double valx = 0.; double valy = 0.; double valz = 0.;
     #pragma acc loop reduction(+:valx,valy,valz)
      for (int j=0;j<gs;j++)
      {
        int ii1 = i1-s1;
        int ii2 = i2-s3;

        double v1 = val1[ii1][j]; double v2 = val2[ii1][j];
        valx += v1*val3x[ii2][3*j+0]; valy += v1*val3x[ii2][3*j+1]; valz += v1*val3x[ii2][3*j+2];
        valx += v2*val4x[ii2][3*j+0]; valy += v2*val4x[ii2][3*j+1]; valz += v2*val4x[ii2][3*j+2];
      }
      double nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
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
void reduce_2c2ds(int m1, int n1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, float** val1x, float** val2x, float** val3x, float** val4x, int iN, int N, int natoms, double scalar, double* xyz_grad)
#else
void reduce_2c2ds(int m1, int n1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, float** val1x, float** val2x, float** val3x, float** val4x, int iN, int N, int natoms, double scalar, float* xyz_grad)
#endif
{
  int N2 = N*N;
  int gs3 = 3*gs;

  double valxt1 = 0.; double valyt1 = 0.; double valzt1 = 0.;
  double valxt2 = 0.; double valyt2 = 0.; double valzt2 = 0.;

 #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs],val4[0:iN][0:gs],norms[0:N2],val1x[0:iN][0:gs3],val2x[0:iN][0:gs3],val3x[0:iN][0:gs3],val4x[0:iN][0:gs3],dpq[0:N2]) reduction(+:valxt1,valxt2,valyt1,valyt2,valzt1,valzt2)
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s3;i2<s4;i2++)
    {
      double valx1 = 0.; double valy1 = 0.; double valz1 = 0.;
      double valx2 = 0.; double valy2 = 0.; double valz2 = 0.;
     #pragma acc loop reduction(+:valx1,valx2,valy1,valy2,valz1,valz2)
      for (int j=0;j<gs;j++)
      {
        int ii1 = i1-s1;
        int ii2 = i2-s3;

        double v1 = val1[ii1][j]; double v2 = val2[ii1][j];
        valx2 += v1*val3x[ii2][3*j+0]; valy2 += v1*val3x[ii2][3*j+1]; valz2 += v1*val3x[ii2][3*j+2];
        valx2 += v2*val4x[ii2][3*j+0]; valy2 += v2*val4x[ii2][3*j+1]; valz2 += v2*val4x[ii2][3*j+2];

        double v3 = val3[ii2][j]; double v4 = val4[ii2][j];
        valx1 += v3*val1x[ii1][3*j+0]; valy1 += v3*val1x[ii1][3*j+1]; valz1 += v3*val1x[ii1][3*j+2];
        valx1 += v4*val2x[ii1][3*j+0]; valy1 += v4*val2x[ii1][3*j+1]; valz1 += v4*val2x[ii1][3*j+2];
      }
      double nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
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
  double vxa = (valxt1 + valxt2)*0.5;
  double vya = (valyt1 + valyt2)*0.5;
  double vza = (valzt1 + valzt2)*0.5;
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
void reduce_2c3d(int m1, int n1, int p1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val1x, float** val2x, float** val3x, float** val4x, float** val5x, float** val6x, int iN, int N, int natoms, double scalar, double* xyz_grad)
#else
void reduce_2c3d(int m1, int n1, int p1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, int iN, int N, int natoms, double scalar, float* xyz_grad)
#endif
{
  int N2 = N*N;
  int gs3 = 3*gs;

  double valxt1 = 0.; double valyt1 = 0.; double valzt1 = 0.;
  double valxt2 = 0.; double valyt2 = 0.; double valzt2 = 0.;

#if USE_ACC
 #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs],val4[0:iN][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs],val1x[0:iN][0:gs3],val2x[0:iN][0:gs3],val3x[0:iN][0:gs3],val4x[0:iN][0:gs3],val5x[0:iN][0:gs3],val6x[0:iN][0:gs3],norms[0:N2],dpq[0:N2],xyz_grad[0:3*natoms]) reduction(+:valxt1,valyt1,valzt1,valxt2,valyt2,valzt2)
#endif
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s3;i2<s4;i2++)
    {
      int ii1 = i1-s1;
      int ii2 = i2-s3;
    
      double valx1 = 0.; double valy1 = 0.; double valz1 = 0.;
      double valx2 = 0.; double valy2 = 0.; double valz2 = 0.;
     #pragma acc loop reduction(+:valx1,valx2,valy1,valy2,valz1,valz2)
      for (int j=0;j<gs;j++)
      {
        float v1 = val1[ii1][j]; float v2 = val2[ii1][j]; float v5 = val5[ii1][j];
        valx2 += v1*val3x[ii2][3*j+0]; valx2 += v2*val4x[ii2][3*j+0]; valx2 += v5*val6x[ii2][3*j+0];
        valy2 += v1*val3x[ii2][3*j+1]; valy2 += v2*val4x[ii2][3*j+1]; valy2 += v5*val6x[ii2][3*j+1];
        valz2 += v1*val3x[ii2][3*j+2]; valz2 += v2*val4x[ii2][3*j+2]; valz2 += v5*val6x[ii2][3*j+2];
        float v3 = val3[ii2][j]; float v4 = val4[ii2][j]; float v6 = val6[ii2][j];
        valx1 += v3*val1x[ii1][3*j+0]; valx1 += v4*val2x[ii1][3*j+0]; valx1 += v6*val5x[ii1][3*j+0];
        valy1 += v3*val1x[ii1][3*j+1]; valy1 += v4*val2x[ii1][3*j+1]; valy1 += v6*val5x[ii1][3*j+1];
        valz1 += v3*val1x[ii1][3*j+2]; valz1 += v4*val2x[ii1][3*j+2]; valz1 += v6*val5x[ii1][3*j+2];
      }

      double nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
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
void reduce_2c3de(int m1, int s1, int s2, int s3, int s4, int gs, double* norms, double* Pao, float** val1, float** val2, float** val3, float** val1x, float** val2x, float** val3x, int N, int imaxN, int natoms, double scalar, double* xyz_grad)
#else
void reduce_2c3de(int m1, int s1, int s2, int s3, int s4, int gs, double* norms, double* Pao, float** val1, float** val2, float** val3, float** val1x, float** val2x, float** val3x, int N, int imaxN, float scalar, float* xyz_grad)
#endif
{
  int N2 = N*N;
  int N3 = 3*natoms;
  int gs3 = 3*gs;

  double valxt = 0.; double valyt = 0.; double valzt = 0.;

#if USE_ACC
 #pragma acc parallel loop collapse(2) present(val1[0:imaxN][0:gs],val2[0:imaxN][0:gs],val3[0:imaxN][0:gs],val1x[0:imaxN][0:gs3],val2x[0:imaxN][0:gs3],val3x[0:imaxN][0:gs3],norms[0:N2],Pao[0:N2]) reduction(+:valxt,valyt,valzt)
#endif
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3;

    double valx = 0.; double valy = 0.; double valz = 0.;

   #pragma acc loop reduction(+:valx,valy,valz)
    for (int j=0;j<gs;j++)
    {
      double v1 = val1[ii1][j]; double v2 = val2[ii1][j]; double v3 = val3[ii1][j];

      valx += v1*val1x[ii2][3*j+0]; valy += v1*val1x[ii2][3*j+1]; valz += v1*val1x[ii2][3*j+2];
      valx += v2*val2x[ii2][3*j+0]; valy += v2*val2x[ii2][3*j+1]; valz += v2*val2x[ii2][3*j+2];
      valx += v3*val3x[ii2][3*j+0]; valy += v3*val3x[ii2][3*j+1]; valz += v3*val3x[ii2][3*j+2];
    }

    double nd1 = Pao[i1*N+i2]*norms[i1*N+i2];
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
void reduce_3c2d112(int m1, int n1, int s1, int s2, int s3, int s4, int s5, int s6, int gs, double* norms1, double* norms2, double* dC, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val1x, float** val2x, float** val3x, float** val4x, float** val5x, float** val6x, int N, int Naux, int imaxN, int imaxNa, int natoms, double scalar, double* xyz_grad)
#else
void reduce_3c2d112(int m1, int n1, int s1, int s2, int s3, int s4, int s5, int s6, int gs, double* norms1, double* norms2, double* dC, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val1x, float** val2x, float** val3x, float** val4x, float** val5x, float** val6x, int N, int Naux, int imaxN, int imaxNa, float scalar, float* xyz_grad)
#endif
{
  int N2 = N*N;
  int N2a = N2*Naux;
  int N3 = 3*natoms;
  int gs3 = 3*gs;

  double valxt1 = 0.; double valyt1 = 0.; double valzt1 = 0.;
  double valxt2 = 0.; double valyt2 = 0.; double valzt2 = 0.;

#if USE_ACC
 #pragma acc parallel loop collapse(3) present(val1[0:imaxNa][0:gs],val2[0:imaxN][0:gs],val3[0:imaxN][0:gs],val4[0:imaxNa][0:gs],val5[0:imaxN][0:gs],val6[0:imaxN][0:gs],val1x[0:imaxNa][0:gs3],val2x[0:imaxN][0:gs3],val3x[0:imaxN][0:gs3],val4x[0:imaxNa][0:gs],val5x[0:imaxN][0:gs],val6x[0:imaxN][0:gs],norms1[0:Naux],norms2[0:N2],dC[0:N2a]) reduction(+:valxt1,valyt1,valzt1,valxt2,valyt2,valzt2)
#endif
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  for (int i3=s5;i3<s6;i3++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3; int ii3 = i3-s5;

    double valx1 = 0.; double valy1 = 0.; double valz1 = 0.;
    double valx2 = 0.; double valy2 = 0.; double valz2 = 0.;

   #pragma acc loop reduction(+:valx1,valy1,valz1,valx2,valy2,valz2)
    for (int j=0;j<gs;j++)
    {
      double v1 = val1[ii1][j]; double v2 = val2[ii2][j]; double v3 = val3[ii3][j];
      double v12 = v1*v2; double v13 = v1*v3; double v23 = v2*v3;

      valx2 += v12*val3x[ii3][3*j+0]; valy2 += v12*val3x[ii3][3*j+1]; valz2 += v12*val3x[ii3][3*j+2];
      valx1 += v13*val2x[ii2][3*j+0]; valy1 += v13*val2x[ii2][3*j+1]; valz1 += v13*val2x[ii2][3*j+2];
      valx1 += v23*val1x[ii1][3*j+0]; valy1 += v23*val1x[ii1][3*j+1]; valz1 += v23*val1x[ii1][3*j+2];

      double v4 = val4[ii1][j]; double v5 = val5[ii2][j]; double v6 = val6[ii3][j];
      double v45 = v4*v5; double v46 = v4*v6; double v56 = v5*v6;

      valx2 += v45*val6x[ii3][3*j+0]; valy2 += v45*val6x[ii3][3*j+1]; valz2 += v45*val6x[ii3][3*j+2];
      valx1 += v46*val5x[ii2][3*j+0]; valy1 += v46*val5x[ii2][3*j+1]; valz1 += v46*val5x[ii2][3*j+2];
      valx1 += v56*val4x[ii1][3*j+0]; valy1 += v56*val4x[ii1][3*j+1]; valz1 += v56*val4x[ii1][3*j+2];
    }

    double nd1 = dC[i1*N2+i2*N+i3]*norms1[i1]*norms2[i2*N+i3];
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
void reduce_3c2d122(int m1, int n1, int s1, int s2, int s3, int s4, int gs, double* norms1, double* norms2, double* dC, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val1x, float** val2x, float** val3x, float** val4x, float** val5x, float** val6x, int N, int Naux, int imaxN, int imaxNa, int natoms, double scalar, double* xyz_grad)
#else
void reduce_3c2d122(int m1, int n1, int s1, int s2, int s3, int s4, int gs, double* norms1, double* norms2, double* dC, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val1x, float** val2x, float** val3x, float** val4x, float** val5x, float** val6x, int N, int Naux, int imaxN, int imaxNa, float scalar, float* xyz_grad)
#endif
{
  int N2 = N*N;
  int N2a = N2*Naux;
  int N3 = 3*natoms;
  int gs3 = 3*gs;

  double valxt1 = 0.; double valyt1 = 0.; double valzt1 = 0.;
  double valxt2 = 0.; double valyt2 = 0.; double valzt2 = 0.;

#if USE_ACC
 #pragma acc parallel loop collapse(3) present(val1[0:imaxNa][0:gs],val2[0:imaxN][0:gs],val3[0:imaxN][0:gs],val4[0:imaxNa][0:gs],val5[0:imaxN][0:gs],val6[0:imaxN][0:gs],val1x[0:imaxNa][0:gs3],val2x[0:imaxN][0:gs3],val3x[0:imaxN][0:gs3],val4x[0:imaxNa][0:gs],val5x[0:imaxN][0:gs],val6x[0:imaxN][0:gs],norms1[0:Naux],norms2[0:N2],dC[0:N2a]) reduction(+:valxt1,valyt1,valzt1,valxt2,valyt2,valzt2)
#endif
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  for (int i3=s3;i3<s4;i3++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3; int ii3 = i3-s3;

    double valx1 = 0.; double valy1 = 0.; double valz1 = 0.;
    double valx2 = 0.; double valy2 = 0.; double valz2 = 0.;

   #pragma acc loop reduction(+:valx1,valy1,valz1,valx2,valy2,valz2)
    for (int j=0;j<gs;j++)
    {
      double v1 = val1[ii1][j]; double v2 = val2[ii2][j]; double v3 = val3[ii3][j];
      double v12 = v1*v2; double v13 = v1*v3; double v23 = v2*v3;

      valx2 += v12*val3x[ii3][3*j+0]; valy2 += v12*val3x[ii3][3*j+1]; valz2 += v12*val3x[ii3][3*j+2];
      valx2 += v13*val2x[ii2][3*j+0]; valy2 += v13*val2x[ii2][3*j+1]; valz2 += v13*val2x[ii2][3*j+2];
      valx1 += v23*val1x[ii1][3*j+0]; valy1 += v23*val1x[ii1][3*j+1]; valz1 += v23*val1x[ii1][3*j+2];

      double v4 = val4[ii1][j]; double v5 = val5[ii2][j]; double v6 = val6[ii3][j];
      double v45 = v4*v5; double v46 = v4*v6; double v56 = v5*v6;

      valx2 += v45*val6x[ii3][3*j+0]; valy2 += v45*val6x[ii3][3*j+1]; valz2 += v45*val6x[ii3][3*j+2];
      valx2 += v46*val5x[ii2][3*j+0]; valy2 += v46*val5x[ii2][3*j+1]; valz2 += v46*val5x[ii2][3*j+2];
      valx1 += v56*val4x[ii1][3*j+0]; valy1 += v56*val4x[ii1][3*j+1]; valz1 += v56*val4x[ii1][3*j+2];
    }

    double nd1 = dC[i1*N2+i2*N+i3]*norms1[i1]*norms2[i2*N+i3];
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
void reduce_3c3d(int m1, int n1, int p1, int s1, int s2, int s3, int s4, int s5, int s6, int gs, double* norms1, double* norms2, double* dC, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val7, float** val8, float** val9, float** val1x, float** val2x, float** val3x, float** val4x, float** val5x, float** val6x, float** val7x, float** val8x, float** val9x, int N, int Naux, int imaxN, int imaxNa, int natoms, double scalar, double* xyz_grad)
#else
void reduce_3c3d(int m1, int n1, int p1, int s1, int s2, int s3, int s4, int s5, int s6, int gs, double* norms1, double* norms2, double* dC, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val7, float** val8, float** val9, float** val1x, float** val2x, float** val3x, float** val4x, float** val5x, float** val6x, float** val7x, float** val8x, float** val9x, int N, int Naux, int imaxN, int imaxNa, float scalar, float* xyz_grad)
#endif
{
  int N2 = N*N;
  int N2a = N2*Naux;
  int N3 = 3*natoms;
  int gs3 = 3*gs;

  double valxt1 = 0.; double valyt1 = 0.; double valzt1 = 0.;
  double valxt2 = 0.; double valyt2 = 0.; double valzt2 = 0.;
  double valxt3 = 0.; double valyt3 = 0.; double valzt3 = 0.;

#if USE_ACC
 #pragma acc parallel loop collapse(3) present(val1[0:imaxNa][0:gs],val2[0:imaxN][0:gs],val3[0:imaxN][0:gs],val4[0:imaxNa][0:gs],val5[0:imaxN][0:gs],val6[0:imaxN][0:gs],val7[0:imaxNa][0:gs],val8[0:imaxN][0:gs],val9[0:imaxN][0:gs],val1x[0:imaxNa][0:gs3],val2x[0:imaxN][0:gs3],val3x[0:imaxN][0:gs3],val4x[0:imaxNa][0:gs3],val5x[0:imaxN][0:gs3],val6x[0:imaxN][0:gs3],val7x[0:imaxNa][0:gs3],val8x[0:imaxN][0:gs3],val9x[0:imaxN][0:gs3],norms1[0:Naux],norms2[0:N2],dC[0:N2a]) reduction(+:valxt1,valyt1,valzt1,valxt2,valyt2,valzt2,valxt3,valyt3,valzt3)
#endif
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  for (int i3=s5;i3<s6;i3++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3; int ii3 = i3-s5;

    double valx1 = 0.; double valy1 = 0.; double valz1 = 0.;
    double valx2 = 0.; double valy2 = 0.; double valz2 = 0.;
    double valx3 = 0.; double valy3 = 0.; double valz3 = 0.;

   #pragma acc loop reduction(+:valx1,valy1,valz1,valx2,valy2,valz2,valx3,valy3,valz3)
    for (int j=0;j<gs;j++)
    {
      double v1 = val1[ii1][j]; double v2 = val2[ii2][j]; double v3 = val3[ii3][j];
      double v12 = v1*v2; double v13 = v1*v3; double v23 = v2*v3;
      valx3 += v12*val3x[ii3][3*j+0]; valy3 += v12*val3x[ii3][3*j+1]; valz3 += v12*val3x[ii3][3*j+2];
      valx2 += v13*val2x[ii2][3*j+0]; valy2 += v13*val2x[ii2][3*j+1]; valz2 += v13*val2x[ii2][3*j+2];
      valx1 += v23*val1x[ii1][3*j+0]; valy1 += v23*val1x[ii1][3*j+1]; valz1 += v23*val1x[ii1][3*j+2];

      double v4 = val4[ii1][j]; double v5 = val5[ii2][j]; double v6 = val6[ii3][j];
      double v45 = v4*v5; double v46 = v4*v6; double v56 = v5*v6;
      valx3 += v45*val6x[ii3][3*j+0]; valy3 += v45*val6x[ii3][3*j+1]; valz3 += v45*val6x[ii3][3*j+2];
      valx2 += v46*val5x[ii2][3*j+0]; valy2 += v46*val5x[ii2][3*j+1]; valz2 += v46*val5x[ii2][3*j+2];
      valx1 += v56*val4x[ii1][3*j+0]; valy1 += v56*val4x[ii1][3*j+1]; valz1 += v56*val4x[ii1][3*j+2];

      double v7 = val7[ii1][j]; double v8 = val8[ii2][j]; double v9 = val9[ii3][j];
      double v78 = v7*v8; double v79 = v7*v9; double v89 = v8*v9;
      valx3 += v78*val9x[ii3][3*j+0]; valy3 += v78*val9x[ii3][3*j+1]; valz3 += v78*val9x[ii3][3*j+2];
      valx2 += v79*val8x[ii2][3*j+0]; valy2 += v79*val8x[ii2][3*j+1]; valz2 += v79*val8x[ii2][3*j+2];
      valx1 += v89*val7x[ii1][3*j+0]; valy1 += v89*val7x[ii1][3*j+1]; valz1 += v89*val7x[ii1][3*j+2];
    }

    double nd1 = dC[i1*N2+i2*N+i3]*norms1[i1]*norms2[i2*N+i3];
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
void reduce_3c1(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, float* valt1, int N, int Naux, int imaxN, int imaxNa, double* C)
#else
void reduce_3c1(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, float* valt1, int N, int Naux, int imaxN, int imaxNa, float* C)
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
      double val = 0.;
     #else
      float val = 0.f;
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

//fully double precision
void reduce_3c1b(int s1, int s2, int s3, int s4, int gs, double** val1, double** val2, double** val3, int N, int Naux, int imaxN, int imaxNa, double* C)
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

    double val = 0.;

   #pragma acc loop reduction(+:val)
    for (int j=0;j<gs;j++)
      val += val1[ii1][j] * val2[ii2][j] * val3[ii3][j];
 
    C[i1*N2+i2*N+i3] += val;
  } //i1,i2,i3

  return;
}

//for PS batching
void reduce_3c1br(int s1, int s2, int s3, int s4, int gs, double** val1, double** val2, double** val3, int N, int Naux, int imaxN, int imaxNa, double* C)
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

    double val = 0.;

   #pragma acc loop reduction(+:val)
    for (int j=0;j<gs;j++)
      val += val1[ii1][j] * val2[ii2][j] * val3[ii3][j];
 
    C[i1*N2+i2*N+i3] += val;
  } //i1,i2,i3

  return;
}

//fully double precision
void reduce_3c1b(int s1, int s2, int s3, int s4, int s5, int s6, int gs, double** val1, double** val2, double** val3, int N, int Naux, int imaxN, int imaxNa, double* C)
{
  int N2 = N*N;
  int N2a = N2*Naux;

#if USE_ACC
 #pragma acc parallel loop collapse(3) present(val1[0:imaxNa][0:gs],val2[0:imaxN][0:gs],val3[0:imaxN][0:gs],C[0:N2a]) 
#endif
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  for (int i3=s5;i3<s6;i3++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3; int ii3 = i3-s5;

    double val = 0.;

   #pragma acc loop reduction(+:val)
    for (int j=0;j<gs;j++)
      val += val1[ii1][j] * val2[ii2][j] * val3[ii3][j];
 
    C[i1*N2+i2*N+i3] += val;
  } //i1,i2,i3

  return;
}

#if RED_DOUBLE
void reduce_3c1b(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, int N, int Naux, int imaxN, int imaxNa, double* C)
#else
void reduce_3c1b(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, int N, int Naux, int imaxN, int imaxNa, float* C)
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
    double val = 0.;
   #else
    float val = 0.f;
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
void reduce_3c2(int s1, int s2, int s3, int s4, int s5, int s6, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float* valt1, float* valt2, int N, int Naux, int imaxN, int imaxNa, double* C)
#else
void reduce_3c2(int s1, int s2, int s3, int s4, int s5, int s6, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float* valt1, float* valt2, int N, int Naux, int imaxN, int imaxNa, float* C)
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
      double val = 0.;
     #else
      float val = 0.f;
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
void reduce_3c2b(int s1, int s2, int s3, int s4, int s5, int s6, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, int N, int Naux, int imaxN, int imaxNa, double* C)
#else
void reduce_3c2b(int s1, int s2, int s3, int s4, int s5, int s6, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, int N, int Naux, int imaxN, int imaxNa, float* C)
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
    double val = 0.;
   #else
    float val = 0.f;
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
void reduce_3c3(int s1, int s2, int s3, int s4, int s5, int s6, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val7, float** val8, float** val9, float* valt1, float* valt2, float* valt3, int N, int Naux, int imaxN, int imaxNa, double* C)
#else
void reduce_3c3(int s1, int s2, int s3, int s4, int s5, int s6, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val7, float** val8, float** val9, float* valt1, float* valt2, float* valt3, int N, int Naux, int imaxN, int imaxNa, float* C)
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
      double val = 0.;
     #else
      float val = 0.f;
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
void reduce_3c3b(int s1, int s2, int s3, int s4, int s5, int s6, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val7, float** val8, float** val9, int N, int Naux, int imaxN, int imaxNa, double* C)
#else
void reduce_3c3b(int s1, int s2, int s3, int s4, int s5, int s6, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val7, float** val8, float** val9, int N, int Naux, int imaxN, int imaxNa, float* C)
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
    double val = 0.;
   #else
    float val = 0.f;
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
void reduce_2c2v(int p1, int s1, int s2, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, int iN, int N, int nc, double scalar, double* V)
#else
void reduce_2c2v(int p1, int s1, int s2, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, int iN, int N, int nc, double scalar, float* V)
#endif
{
  int N2 = N*N;

  double valt = 0.;

 #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs],val4[0:iN][0:gs],norms[0:N2],dpq[0:N2]) reduction(+:valt)
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s1;i2<s2;i2++)
    {
      double val = 0.;
     #pragma acc loop reduction(+:val)
      for (int j=0;j<gs;j++)
      {
        int ii1 = i1-s1;
        int ii2 = i2-s1;

        double v3 = val3[ii2][j]; double v4 = val4[ii2][j];
        val += v3*val1[ii1][j];
        val += v4*val2[ii1][j]; 
      }

      double nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
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
void reduce_2c2vd(int p1, int s1, int s2, int gs, double* norms, double* dpq, float** val1x, float** val2x, float** val3, float** val4, int iN, int N, int nc, double scalar, double* dV)
#else
void reduce_2c2vd(int p1, int s1, int s2, int gs, double* norms, double* dpq, float** val1x, float** val2x, float** val3, float** val4, int iN, int N, int nc, double scalar, float* dV)
#endif
{
  int N2 = N*N;
  int gs3 = 3*gs;

  double valxt = 0.; double valyt = 0.; double valzt = 0.;

 #pragma acc parallel loop collapse(2) present(val1x[0:iN][0:gs3],val2x[0:iN][0:gs3],val3[0:iN][0:gs],val4[0:iN][0:gs],norms[0:N2],dpq[0:N2]) reduction(+:valxt,valyt,valzt)
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s1;i2<s2;i2++)
    {
      double valx = 0.; double valy = 0.; double valz = 0.;
     #pragma acc loop reduction(+:valx,valy,valz)
      for (int j=0;j<gs;j++)
      {
        int ii1 = i1-s1;
        int ii2 = i2-s1;

        double v3 = val3[ii2][j]; double v4 = val4[ii2][j];
        valx += v3*val1x[ii1][3*j+0]; valy += v3*val1x[ii1][3*j+1]; valz += v3*val1x[ii1][3*j+2];
        valx += v4*val2x[ii1][3*j+0]; valy += v4*val2x[ii1][3*j+1]; valz += v4*val2x[ii1][3*j+2];
      }
      double nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
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
void reduce_2c3v(int p1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, int iN, int N, int nc, double scalar, double* V)
#else
void reduce_2c3v(int p1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, int iN, int N, int nc, double scalar, float* V)
#endif
{
  int N2 = N*N;

  double valt = 0.;

 #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs],val4[0:iN][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs],norms[0:N2],dpq[0:N2]) reduction(+:valt)
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s3;i2<s4;i2++)
    {
      double val = 0.;
     #pragma acc loop reduction(+:val)
      for (int j=0;j<gs;j++)
      {
        int ii1 = i1-s1;
        int ii2 = i2-s3;

        double v1 = val1[ii1][j]; double v2 = val2[ii1][j]; double v3 = val3[ii1][j];
        val += v1*val4[ii2][j];
        val += v2*val5[ii2][j]; 
        val += v3*val6[ii2][j]; 
      }

      double nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
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
void reduce_2c3vd(int p1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1x, float** val2x, float** val3x, float** val4, float** val5, float** val6, int iN, int N, int nc, double scalar, double* dV)
#else
void reduce_2c3vd(int p1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1x, float** val2x, float** val3x, float** val4, float** val5, float** val6, int iN, int N, int nc, double scalar, float* dV)
#endif
{
  int N2 = N*N;
  int gs3 = 3*gs;

  double valxt = 0.; double valyt = 0.; double valzt = 0.;

 #pragma acc parallel loop collapse(2) present(val1x[0:iN][0:gs3],val2x[0:iN][0:gs3],val3x[0:iN][0:gs],val4[0:iN][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs],norms[0:N2],dpq[0:N2]) reduction(+:valxt,valyt,valzt)
  for (int i1=s1;i1<s2;i1++)
  {
    for (int i2=s3;i2<s4;i2++)
    {
      double valx = 0.; double valy = 0.; double valz = 0.;
     #pragma acc loop reduction(+:valx,valy,valz)
      for (int j=0;j<gs;j++)
      {
        int ii1 = i1-s1;
        int ii2 = i2-s3;

        double v4 = val4[ii2][j]; double v5 = val5[ii2][j]; double v6 = val6[ii2][j];
        valx += v4*val1x[ii1][3*j+0]; valy += v4*val1x[ii1][3*j+1]; valz += v4*val1x[ii1][3*j+2];
        valx += v5*val2x[ii1][3*j+0]; valy += v5*val2x[ii1][3*j+1]; valz += v5*val2x[ii1][3*j+2];
        valx += v6*val3x[ii1][3*j+0]; valy += v6*val3x[ii1][3*j+1]; valz += v6*val3x[ii1][3*j+2];
      }
      double nd1 = norms[i1*N+i2]*dpq[i1*N+i2];
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
