#include "integrals.h"













void copy_grid(int gs, double* grid1, double* grid2)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid1[0:6*gs],grid2[0:6*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    grid1[6*i+0] = grid2[6*i+0];
    grid1[6*i+1] = grid2[6*i+1];
    grid1[6*i+2] = grid2[6*i+2];
    grid1[6*i+3] = grid2[6*i+3];
  }

  return;
}




void recenter_grid_zero(int gs, double* grid, double x2, double y2, double z2)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    grid[6*i+0] += x2;
    grid[6*i+1] += y2;
    grid[6*i+2] += z2;
    double xn = grid[6*i+0];
    double yn = grid[6*i+1];
    double zn = grid[6*i+2];

    double r1 = sqrt(xn*xn+yn*yn+zn*zn);
    grid[6*i+3] = r1;
  }

  return;
}

void recenter_grid(int gs, double* grid, double x2, double y2, double z2)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    grid[6*i+0] += x2;
    grid[6*i+1] += y2;
    grid[6*i+2] += z2;
    double xn = grid[6*i+0];
    double yn = grid[6*i+1];
    double zn = grid[6*i+2];

   //distance to other center
    double r2 = sqrt(xn*xn+yn*yn+zn*zn);
    grid[6*i+4] = r2;
  }

  return;
}

void add_r1_to_grid(int gs, double* grid1, double A2, double B2, double C2)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid1[0:6*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    double x1 = grid1[6*i+0];
    double y1 = grid1[6*i+1];
    double z1 = grid1[6*i+2];
    double x12 = x1-A2;
    double y12 = y1-B2;
    double z12 = z1-C2;
    double r2 = sqrt(x12*x12+y12*y12+z12*z12);
    grid1[6*i+3] = r2;
  }

  return;
}

void add_r2_to_grid(int gs, double* grid1, double A2, double B2, double C2)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid1[0:6*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    double x1 = grid1[6*i+0];
    double y1 = grid1[6*i+1];
    double z1 = grid1[6*i+2];
    double x12 = x1-A2;
    double y12 = y1-B2;
    double z12 = z1-C2;
    double r2 = sqrt(x12*x12+y12*y12+z12*z12);
    grid1[6*i+4] = r2;
  }

  return;
}

/*
void add_r2_to_grid(int gs, double* grid1, float A2, float B2, float C2)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid1[0:6*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    double x1 = grid1[6*i+0];
    double y1 = grid1[6*i+1];
    double z1 = grid1[6*i+2];
    double x12 = x1-A2;
    double y12 = y1-B2;
    double z12 = z1-C2;
    double r2 = sqrt(x12*x12+y12*y12+z12*z12);
    grid1[6*i+4] = r2;
  }

  return;
}
*/

void generate_central_grid_2d(bool use_murak, double* grid1, double* wt1, float Z1, int nrad, int nang, double* ang_g, double* ang_w)
{
  double* r = new double[nrad];
  double* w = new double[nrad];

 #if USE_ACC
  #pragma acc enter data create(r[0:nrad],w[0:nrad])
 #endif

  if (use_murak)
    get_murak_grid(nrad,r,w,Z1,3);
  else
  {
    int m = 3;
    double zeta = Z1/10.; //exponent
    get_murak_grid_zeta(nrad,r,w,zeta,m);
  }

  int gs = nrad*nang;

#if USE_ACC
 #pragma acc parallel loop independent present(r[0:nrad],w[0:nrad],ang_g[0:3*nang],ang_w[0:nang],grid1[0:6*gs],wt1[0:gs])
#endif
  for (int i=0;i<nrad;i++)
  {
    double r1 = r[i];
    double wr1 = w[i];

  #if USE_ACC
   #pragma acc loop independent 
  #endif
    for (int j=0;j<nang;j++)
    {
      double w1 = wr1*ang_w[j];

     //grid positions
      double x1 = r1 * ang_g[3*j+0];
      double y1 = r1 * ang_g[3*j+1];
      double z1 = r1 * ang_g[3*j+2];

      int wg = i*nang+j;
      grid1[6*wg+0] = x1;
      grid1[6*wg+1] = y1;
      grid1[6*wg+2] = z1;
      grid1[6*wg+3] = r1;
      grid1[6*wg+4] = r1;
      grid1[6*wg+5] = r1;
      wt1[wg] = w1;
    }
  }

#if USE_ACC
  #pragma acc exit data delete(r[0:nrad],w[0:nrad])
#endif

  delete [] r;
  delete [] w;

  return;
}

#if 0
void generate_central_grid(float* grid1, float* wt1, float* val1, int need_inr, float Z1, int n1, int l1, float zeta1, int nrad, int nang, float* ang_g, float* ang_w)
{
  printf("\n WARNING: shouldn't be here in generate_central_grid() \n");

  float* r = new float[nrad];
  float* w = new float[nrad];
  float* er = new float[nrad];
  float* inr = new float[nrad];

 #if DEBUG
  printf("  generate_central_grid for: %i %i (zeta: %8.5f) \n",n1,l1,zeta1);
 #endif

 #if USE_ACC
  #pragma acc enter data create(r[0:nrad],w[0:nrad],er[0:nrad],inr[0:nrad])
 #endif

  get_murak_grid_f(nrad,r,w,er,Z1,zeta1,3);
  if (need_inr)
    get_inr(n1,l1,zeta1,nrad,r,inr);
  else
    acc_assign(nrad,inr,0.);

  int gs = nrad*nang;

#if USE_ACC
 #pragma acc parallel loop independent present(r[0:nrad],w[0:nrad],inr[0:nrad],ang_g[0:3*nang],ang_w[0:nang],grid1[0:6*gs],wt1[0:gs],val1[0:gs])
#endif
  for (int i=0;i<nrad;i++)
  {
    float r1 = r[i];
    float wr1 = w[i];
    float inr1 = inr[i];
 
  #if USE_ACC
   #pragma acc loop independent 
  #endif
    for (int j=0;j<nang;j++)
    {
      float w1 = wr1*ang_w[j];

     //grid positions
      float x1 = r1 * ang_g[3*j+0];
      float y1 = r1 * ang_g[3*j+1];
      float z1 = r1 * ang_g[3*j+2];

      int wg = i*nang+j;
      grid1[6*wg+0] = x1;
      grid1[6*wg+1] = y1;
      grid1[6*wg+2] = z1;
      grid1[6*wg+3] = r1;
      grid1[6*wg+4] = r1;
      grid1[6*wg+5] = r1;
      wt1[wg] = w1;
      val1[wg] = inr1;
    }
  }

#if USE_ACC
  #pragma acc exit data delete(r[0:nrad],w[0:nrad],er[0:nrad],inr[0:nrad])
#endif

  delete [] r;
  delete [] w;
  delete [] er;
  delete [] inr;

  return;
}
#endif

int get_imax_n2i(int natoms, int N, vector<vector<double> >& basis, int* n2i)
{
  for (int n=0;n<natoms;n++) n2i[n] = 0;

  int imaxN = 1;
  int wa = 0;
  int iprev = 0;
  for (int i=0;i<N;i++)
  {
    int wa2 = basis[i][9];
    if (wa2!=wa)
    {
      int cmaxN = i-iprev; 
      if (cmaxN>imaxN) imaxN = cmaxN;
      n2i[wa] = i;
      wa = wa2;
      iprev = i;
    }
  }
  if (N-iprev>imaxN) imaxN = N-iprev;
  n2i[natoms-1] = N;

  return imaxN;
}
