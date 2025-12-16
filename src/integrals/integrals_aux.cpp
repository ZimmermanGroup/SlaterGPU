#include "integrals.h"
#include "grid_util.h"

#include "../../include/cintwrapper.h"
#include "../../include/cintprep.h"

#include "write.h"

void auto_crash()
{
  //return;

 //debug: what is on devices?
  int nang = 100;
  double* ang_t = new double[nang];

 //#pragma acc exit data delete(ang_t[0:nang])
 #pragma acc parallel loop present(ang_t[0:nang])
  for (int j=0;j<2;j++)
    printf(" debug: %8.5f \n",ang_t[j]);

  delete [] ang_t;
  return;
}

void print_array(int size, float* vec)
{
  for (int i=0;i<size;i++)
    printf(" %6.3e",vec[i]);
  printf("\n");
}

void clean_small_values(int N, float* S)
{
  int N2 = N*N;
 #pragma acc parallel loop independent present(S[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i+1;j<N;j++)
  if (fabs(S[i*N+j])<CL_THRESH)
    S[i*N+j] = S[j*N+i] = 0.;
}

void clean_small_values(int N, double* S)
{
  int N2 = N*N;
 #pragma acc parallel loop independent present(S[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i+1;j<N;j++)
  if (fabs(S[i*N+j])<CL_THRESH)
    S[i*N+j] = S[j*N+i] = 0.;
}

void acc_assign(int size, float* vec, float v1)
{
 #pragma acc parallel loop independent present(vec[0:size])
  for (int j=0;j<size;j++)
    vec[j] = v1;
}

void acc_assign(int size, double* vec, double v1)
{
 #pragma acc parallel loop independent present(vec[0:size])
  for (int j=0;j<size;j++)
    vec[j] = v1;
}

void acc_assign(int tid, int size, double* vec, double v1)
{
 #pragma acc parallel loop independent present(vec[0:size]) //async(tid+1)
  for (int j=0;j<size;j++)
    vec[j] = v1;
}

void acc_assign(int size, float* vec1, float* vec2, float v1)
{
 #pragma acc parallel loop independent present(vec1[0:size],vec2[0:size])
  for (int j=0;j<size;j++)
    vec1[j] = vec2[j] = v1;
}

void acc_assign(int size, float* vec1, float* vec2, float* vec3, float v1)
{
 #pragma acc parallel loop independent present(vec1[0:size],vec2[0:size],vec3[0:size])
  for (int j=0;j<size;j++)
    vec1[j] = vec2[j] = vec3[j] = v1;
}

float acc_sum(int size, float* vec)
{
  float sum = 0.;
 #pragma acc parallel loop present(vec[0:size]) reduction(+:sum)
  for (int i=0;i<size;i++)
    sum += vec[i];

  return sum;
}

void acc_copy(int tid, int size, double* v1, double* v2)
{
  #pragma acc parallel loop independent present(v1[0:size],v2[0:size]) //async(tid+1)
  for (int i=0;i<size;i++)
    v1[i] = v2[i];

  return;
}

void acc_copy(int size, double* v1, double* v2)
{
  #pragma acc parallel loop independent present(v1[0:size],v2[0:size])
  for (int i=0;i<size;i++)
    v1[i] = v2[i];

  return;
}

void acc_copyf(int tid, int size, float* v1, float* v2)
{
  #pragma acc parallel loop independent present(v1[0:size],v2[0:size]) //async(tid+1)
  for (int i=0;i<size;i++)
    v1[i] = v2[i];

  return;
}

void acc_copyf(int size, float* v1, float* v2)
{
  #pragma acc parallel loop independent present(v1[0:size],v2[0:size])
  for (int i=0;i<size;i++)
    v1[i] = v2[i];

  return;
}

void acc_copyf(int size, float* v1, float* v2, float* v3, float* v4)
{
  #pragma acc parallel loop independent present(v1[0:size],v2[0:size],v3[0:size],v4[0:size])
  for (int i=0;i<size;i++)
  {
    v1[i] = v2[i];
    v3[i] = v4[i];
  }

  return;
}

void acc_copyf(int size, float* v1, float* v2, float* v3, float* v4, float* v5, float* v6)
{
  #pragma acc parallel loop independent present(v1[0:size],v2[0:size],v3[0:size],v4[0:size],v5[0:size],v6[0:size])
  for (int i=0;i<size;i++)
  {
    v1[i] = v2[i];
    v3[i] = v4[i];
    v5[i] = v6[i];
  }

  return;
}

void eliminate_small_wt_3(int size, float* wt1, float* wt2, float* wt3)
{
 #pragma acc parallel loop independent present(wt1[0:size],wt2[0:size],wt3[0:size])
  for (int i=0;i<size;i++)
  {
    if (wt1[i]<WT_THRESH)
      wt1[i] = 0.;
    if (wt2[i]<WT_THRESH)
      wt2[i] = 0.;
    if (wt3[i]<WT_THRESH)
      wt3[i] = 0.;
  }

}

void eliminate_small_wt_3(int s1, int size, float* wt1, float* wt2, float* wt3)
{
 #pragma acc parallel loop independent present(wt1[0:size],wt2[0:size],wt3[0:size])
  for (int i=s1;i<size;i++)
  {
    if (wt1[i]<WT_THRESH_D)
      wt1[i] = 0.;
    if (wt2[i]<WT_THRESH_D)
      wt2[i] = 0.;
    if (wt3[i]<WT_THRESH_D)
      wt3[i] = 0.;
  }

}

void eliminate_small_wt(int size, float* wt1)
{
 #pragma acc parallel loop independent present(wt1[0:size])
  for (int i=0;i<size;i++)
  if (wt1[i]<WT_THRESH)
  {
    //printf(" small wt: %12.10f \n",wt1[i]);
    wt1[i] = 0.;
  }

}

void eliminate_small_wt(int s1, int size, float* wt1)
{
 #pragma acc parallel loop independent present(wt1[0:size])
  for (int i=s1;i<size;i++)
  if (wt1[i]<WT_THRESH_D)
  {
    //printf(" small wt: %12.10f \n",wt1[i]);
    wt1[i] = 0.;
  }

}

void copy_grid(int gs, float* grid1, float* grid2)
{
 #pragma acc parallel loop independent present(grid1[0:6*gs],grid2[0:6*gs])
  for (int i=0;i<gs;i++)
  {
    grid1[6*i+0] = grid2[6*i+0];
    grid1[6*i+1] = grid2[6*i+1];
    grid1[6*i+2] = grid2[6*i+2];
    grid1[6*i+3] = grid2[6*i+3];
  }

  return;
}

void copy_grid(int gs, double* grid1, float* grid2)
{
 #pragma acc parallel loop independent present(grid1[0:6*gs],grid2[0:6*gs])
  for (int j=0;j<6*gs;j++)
    grid1[j] = grid2[j];

  return;
}

void copy_grid(int gs, double* grid1, double* grid2)
{
 #pragma acc parallel loop independent present(grid1[0:6*gs],grid2[0:6*gs])
  for (int i=0;i<gs;i++)
  {
    grid1[6*i+0] = grid2[6*i+0];
    grid1[6*i+1] = grid2[6*i+1];
    grid1[6*i+2] = grid2[6*i+2];
    grid1[6*i+3] = grid2[6*i+3];
  }

  return;
}

void copy_grid(int gs, float* grid1, float* wt1, float* grid2, float* wt2)
{
 #pragma acc parallel loop independent present(grid1[0:6*gs],wt1[0:gs],grid2[0:6*gs],wt2[0:gs])
  for (int i=0;i<gs;i++)
  {
    grid2[6*i+0] = grid1[6*i+0];
    grid2[6*i+1] = grid1[6*i+1];
    grid2[6*i+2] = grid1[6*i+2];
    grid2[6*i+3] = grid1[6*i+3];
    wt2[i] = wt1[i];
  }

  return;
}


void recenter_grid_zero(int gs, float* grid, float x2, float y2, float z2)
{
 #pragma acc parallel loop independent present(grid[0:6*gs])
  for (int i=0;i<gs;i++)
  {
    grid[6*i+0] += x2;
    grid[6*i+1] += y2;
    grid[6*i+2] += z2;
    float xn = grid[6*i+0];
    float yn = grid[6*i+1];
    float zn = grid[6*i+2];

    float r1 = sqrtf(xn*xn+yn*yn+zn*zn);
    grid[6*i+3] = r1;
  }

  return;
}

void recenter_grid_zero(int gs, double* grid, double x2, double y2, double z2)
{
 #pragma acc parallel loop independent present(grid[0:6*gs])
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

void recenter_grid_zero(int tid, int gs, double* grid, double x2, double y2, double z2)
{
 #pragma acc parallel loop independent present(grid[0:6*gs]) //async(tid+1)
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

void recenter_grid(int gs, float* grid, float x2, float y2, float z2)
{
 #pragma acc parallel loop independent present(grid[0:6*gs])
  for (int i=0;i<gs;i++)
  {
    grid[6*i+0] += x2;
    grid[6*i+1] += y2;
    grid[6*i+2] += z2;
    float xn = grid[6*i+0];
    float yn = grid[6*i+1];
    float zn = grid[6*i+2];

   //distance to other center
    float r2 = sqrtf(xn*xn+yn*yn+zn*zn);
    grid[6*i+4] = r2;
  }

  return;
}

void recenter_grid(int gs, double* grid, double x2, double y2, double z2)
{
 #pragma acc parallel loop independent present(grid[0:6*gs])
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

void recenter_grid_exp(int gs, float* grid, float* wt, float* val, float x2, float y2, float z2, float zeta2)
{
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
  for (int i=0;i<gs;i++)
  {
    grid[6*i+0] += x2;
    grid[6*i+1] += y2;
    grid[6*i+2] += z2;
    float xn = grid[6*i+0];
    float yn = grid[6*i+1];
    float zn = grid[6*i+2];
    float r2 = grid[6*i+3];

    float r1 = sqrtf(xn*xn+yn*yn+zn*zn);
    float er2 = expf(-zeta2*r2);

   //distance to first center
    grid[6*i+4] = r1;
    //wt[i] *= er2;
    val[i] = er2;
  }

  return;
}

void add_r123_to_grid(int gs, float* grid1, float A1, float B1, float C1, float A2, float B2, float C2, float A3, float B3, float C3)
{
 #pragma acc parallel loop independent present(grid1[0:6*gs])
  for (int i=0;i<gs;i++)
  {
    float x = grid1[6*i+0];
    float y = grid1[6*i+1];
    float z = grid1[6*i+2];
    float x1 = x-A1; float y1 = y-B1; float z1 = z-C1;
    float x2 = x-A2; float y2 = y-B2; float z2 = z-C2;
    float x3 = x-A3; float y3 = y-B3; float z3 = z-C3;
    float r1 = sqrtf(x1*x1+y1*y1+z1*z1);
    float r2 = sqrtf(x2*x2+y2*y2+z2*z2);
    float r3 = sqrtf(x3*x3+y3*y3+z3*z3);
    grid1[6*i+3] = r1;
    grid1[6*i+4] = r2;
    grid1[6*i+5] = r3;
  }

  return;
}

void add_r1_to_grid_6z(int gs, float* grid1, float* grid2, float* grid3, float* grid4, float* grid5, float* grid6)
{
 #pragma acc parallel loop independent present(grid1[0:6*gs],grid2[0:6*gs],grid3[0:6*gs],grid4[0:6*gs],grid5[0:6*gs],grid6[0:6*gs])
  for (int i=0;i<gs;i++)
  {
    float x1 = grid1[6*i+0]; float y1 = grid1[6*i+1]; float z1 = grid1[6*i+2];
    float r1 = sqrtf(x1*x1+y1*y1+z1*z1);
    grid1[6*i+3] = r1;

    float x2 = grid2[6*i+0]; float y2 = grid2[6*i+1]; float z2 = grid2[6*i+2];
    float r2 = sqrtf(x2*x2+y2*y2+z2*z2);
    grid2[6*i+3] = r2;

    float x3 = grid3[6*i+0]; float y3 = grid3[6*i+1]; float z3 = grid3[6*i+2];
    float r3 = sqrtf(x3*x3+y3*y3+z3*z3);
    grid3[6*i+3] = r3;

    float x4 = grid4[6*i+0]; float y4 = grid4[6*i+1]; float z4 = grid4[6*i+2];
    float r4 = sqrtf(x4*x4+y4*y4+z4*z4);
    grid4[6*i+3] = r4;

    float x5 = grid5[6*i+0]; float y5 = grid5[6*i+1]; float z5 = grid5[6*i+2];
    float r5 = sqrtf(x5*x5+y5*y5+z5*z5);
    grid5[6*i+3] = r5;

    float x6 = grid6[6*i+0]; float y6 = grid6[6*i+1]; float z6 = grid6[6*i+2];
    float r6 = sqrtf(x6*x6+y6*y6+z6*z6);
    grid6[6*i+3] = r6;
  }

  return;
}

void add_r1_to_grid(int gs, float* grid1, float A2, float B2, float C2)
{
  int gs6 = 6*gs;
 #pragma acc parallel loop present(grid1[0:gs6])
  for (int i=0;i<gs;i++)
  {
    float x1 = grid1[6*i+0];
    float y1 = grid1[6*i+1];
    float z1 = grid1[6*i+2];
    float x12 = x1-A2;
    float y12 = y1-B2;
    float z12 = z1-C2;
    float r2 = sqrtf(x12*x12+y12*y12+z12*z12);
    grid1[6*i+3] = r2;
  }

  return;
}

void add_r1_to_grid(int tid, int gs, double* grid1, double A2, double B2, double C2)
{
  int gs6 = 6*gs;

 #pragma acc parallel loop present(grid1[0:gs6]) //async(tid+1)
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

void add_r1_to_grid(int gs, double* grid1, double A2, double B2, double C2)
{
  int gs6 = 6*gs;

 #pragma acc parallel loop present(grid1[0:gs6])
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

void add_r2_to_grid(int gs, float* grid1, float A2, float B2, float C2)
{
 #pragma acc parallel loop independent present(grid1[0:6*gs])
  for (int i=0;i<gs;i++)
  {
    float x1 = grid1[6*i+0];
    float y1 = grid1[6*i+1];
    float z1 = grid1[6*i+2];
    float x12 = x1-A2;
    float y12 = y1-B2;
    float z12 = z1-C2;
    float r2 = sqrtf(x12*x12+y12*y12+z12*z12);
    grid1[6*i+4] = r2;
  }

  return;
}

void add_r2_to_grid(int gs, double* grid1, double A2, double B2, double C2)
{
 #pragma acc parallel loop independent present(grid1[0:6*gs])
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

void add_r3_to_grid(int gs, float* grid1, float A3, float B3, float C3)
{
 #pragma acc parallel loop independent present(grid1[0:6*gs])
  for (int i=0;i<gs;i++)
  {
    float x1 = grid1[6*i+0];
    float y1 = grid1[6*i+1];
    float z1 = grid1[6*i+2];
    float x13 = x1-A3;
    float y13 = y1-B3;
    float z13 = z1-C3;
    float r3 = sqrtf(x13*x13+y13*y13+z13*z13);
    grid1[6*i+5] = r3;
  }

  return;
}

int find_center_of_grid(float Z1, int nrad)
{
  float* r = new float[nrad];
  float* w = new float[nrad];

  #pragma acc enter data create(r[0:nrad],w[0:nrad])

  get_murak_grid_f(nrad,r,w,Z1,3);

  int s1 = 0;
 #pragma acc parallel loop present(w[0:nrad])
  for (int i=0;i<nrad;i++)
  if (w[i]<WT_THRESH_D)
    s1 = i;
  s1++;

  //printf(" s1: %2i \n",s1);

  #pragma acc exit data delete(r[0:nrad],w[0:nrad])

  return s1;
}

void get_angular_grid(int size_ang, double* ang_g, double* ang_w)
{
  double* xa = new double[size_ang];
  double* ya = new double[size_ang];
  double* za = new double[size_ang];
  ld_by_order(size_ang,xa,ya,za,ang_w);
  for (int i=0;i<size_ang;i++)
  {
    ang_g[3*i+0] = xa[i];
    ang_g[3*i+1] = ya[i];
    ang_g[3*i+2] = za[i];
  }
  delete [] xa;
  delete [] ya;
  delete [] za;

#if 0
  float sumawt = 0.;
  for (int j=0;j<size_ang;j++)
    sumawt += ang_w[j];
  printf(" ang_w: ");
  for (int j=0;j<size_ang;j++)
    printf(" %8.5f",ang_w[j]);
  printf("\n");
  printf(" sum(ang_w): %8.5f \n",sumawt);
#endif

  return;
}

void rgrid_one_atom(int nrad, int Z1, float* r1)
{
  float w1[nrad];
  #pragma acc enter data create(r1[0:nrad],w1[0:nrad])

  get_murak_grid_f(nrad,r1,w1,Z1,3);
  #pragma acc exit data copyout(r1[0:nrad])
  #pragma acc exit data delete(w1[0:nrad])

  return;
}

void generate_central_grid_3d(int m, double* grid1, double* wt1, float Z1, int nrad, int nang, double* ang_g, double* ang_w)
{
  double* r = new double[nrad];
  double* w = new double[nrad];

  #pragma acc enter data create(r[0:nrad],w[0:nrad])

  //printf("  in generate_central_grid_3d. m_murak: %i \n",m);

  if (m<1) { printf("\n  ERROR: m cannot be less than 1 in get_murak_grid \n"); exit(-1); }

  double zeta = Z1/10.; //exponent
  if (Z1<0.) zeta = -Z1;
  get_murak_grid_zeta(nrad,r,w,zeta,m);

  int gs = nrad*nang;

  for (int i=0;i<nrad;i++)
  {
   #pragma acc parallel loop present(r[0:nrad],w[0:nrad],ang_g[0:3*nang],ang_w[0:nang],grid1[0:6*gs],wt1[0:gs])
    for (int j=0;j<nang;j++)
    {
      double r1 = r[i];
      double wr1 = w[i];
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

  #pragma acc exit data delete(r[0:nrad],w[0:nrad])

  delete [] r;
  delete [] w;

  return;
}

void generate_central_grid_2d(int tid, int wb, int nb, bool use_murak, double* grid1, double* wt1, float Z1, int nrad, int nang, double* ang_g, double* ang_w)
{
  double* r = new double[nrad];
  double* w = new double[nrad];

  #pragma acc enter data create(r[0:nrad],w[0:nrad]) //async(tid+1)

  //printf("  in generate_central_grid_2d. use_murak: %i \n",(int)use_murak);

  //use_murak = 1;
 #if 0
  get_eumac_grid(nrad,r,w,25,1);
 #endif
 #if 1
  if (use_murak)
    get_murak_grid(tid,nrad,r,w,Z1,3);
  else
  {
    int m = 3;
    double zeta = Z1/10.; //exponent
    if (Z1<0.) zeta = -Z1;
    get_murak_grid_zeta(tid,nrad,r,w,zeta,m);
  }
 #endif

  int gs = nrad*nang/nb;

  //int ic = 0;
 #pragma acc parallel loop independent present(r[0:nrad],w[0:nrad],ang_g[0:3*nang],ang_w[0:nang],grid1[0:6*gs],wt1[0:gs]) //async(tid+1)
  for (int i=0;i<nrad;i++)
  if (i%nb==wb)
  {
    int i1 = i/nb;
    //printf(" (%i/%i) \n",i1,ic);

    #pragma acc loop independent
    for (int j=0;j<nang;j++)
    {
      double r1 = r[i];
      double wr1 = w[i];
      double w1 = wr1*ang_w[j];

     //grid positions
      double x1 = r1 * ang_g[3*j+0];
      double y1 = r1 * ang_g[3*j+1];
      double z1 = r1 * ang_g[3*j+2];

      int wg = i1*nang+j;
      grid1[6*wg+0] = x1;
      grid1[6*wg+1] = y1;
      grid1[6*wg+2] = z1;
      grid1[6*wg+3] = r1;
      grid1[6*wg+4] = r1;
      grid1[6*wg+5] = r1;
      wt1[wg] = w1;
    }
    //ic++;
  }

  #pragma acc exit data delete(r[0:nrad],w[0:nrad]) //async(tid+1)

  delete [] r;
  delete [] w;

  if (tid<0)
  {
    #pragma acc wait
  }
  acc_wait_all();

  return;
}

void generate_central_grid_2d(int tid, bool use_murak, double* grid1, double* wt1, float Z1, int nrad, int nang, double* ang_g, double* ang_w)
{
  return generate_central_grid_2d(tid,0,1,use_murak,grid1,wt1,Z1,nrad,nang,ang_g,ang_w);
}

void generate_central_grid_2d(int tid, int wb, int nb, bool use_murak, double* grid1, double* wt1, float Z1, int nrad, int nang)
{
  double* ang_g = new double[nang*3];
  double* ang_w = new double[nang];

  get_angular_grid(nang,ang_g,ang_w);
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])

  generate_central_grid_2d(tid,wb,nb,use_murak,grid1,wt1,Z1,nrad,nang,ang_g,ang_w);

  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  delete [] ang_g;
  delete [] ang_w;
}

void generate_central_grid_2d(int tid, bool use_murak, double* grid1, double* wt1, float Z1, int nrad, int nang)
{
  return generate_central_grid_2d(tid,0,1,use_murak,grid1,wt1,Z1,nrad,nang);
}

void generate_central_grid_2(float* grid1, float* wt1, float Z1, int nrad, int nang, float* ang_g, float* ang_w)
{
  float* r = new float[nrad];
  float* w = new float[nrad];

  #pragma acc enter data create(r[0:nrad],w[0:nrad])

  get_murak_grid_f(nrad,r,w,Z1,3);

  int gs = nrad*nang;

 #pragma acc parallel loop independent present(r[0:nrad],w[0:nrad],ang_g[0:3*nang],ang_w[0:nang],grid1[0:6*gs],wt1[0:gs])
  for (int i=0;i<nrad;i++)
  {
    float r1 = r[i];
    float wr1 = w[i];

   #pragma acc loop independent
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
    }
  }

  #pragma acc exit data delete(r[0:nrad],w[0:nrad])

  delete [] r;
  delete [] w;

  return;
}

void transpose_C(int Naux, int N, float* C)
{
  int N2 = N*N;
  int N2a = N2*Naux;
  int nna = N*Naux;

  float* tmp = new float[N2a];

  for (int i=0;i<Naux;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    tmp[j*nna+k*Naux+i] = C[i*N2+j*N+k];

  for (int i=0;i<N2a;i++)
    C[i] = tmp[i];

  delete [] tmp;

  return;
}

void transpose_C(int Naux, int N, double* C)
{
  int N2 = N*N;
  int N2a = N2*Naux;
  int nna = N*Naux;

  double* tmp = new double[N2a];

  for (int i=0;i<Naux;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    tmp[j*nna+k*Naux+i] = C[i*N2+j*N+k];

  for (int i=0;i<N2a;i++)
    C[i] = tmp[i];

  delete [] tmp;

  return;
}

void copy_symm(int natoms, int N, int Naux, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, float* C, int type)
{
 //copy symmetric 2-atom terms
  int N2 = N*N;

  for (int m=0;m<natoms;m++)
  {
    for (int i1=0;i1<Naux;i1++)
    if (basis_aux[i1][9]==m)
    {
      for (int n=0;n<natoms;n++)
      if (m!=n)
      {
       //2-atom
        if (type==1)
        {
          for (int i2=0;i2<N;i2++)
          if (basis[i2][9]==n)
          {
            for (int i3=0;i3<N;i3++)
            if (basis[i3][9]==m)
              C[i1*N2+i3*N+i2] = C[i1*N2+i2*N+i3];
          }
        }
        else
        {
          for (int i2=0;i2<N;i2++)
          if (basis[i2][9]==m)
          {
            for (int i3=0;i3<N;i3++)
            if (basis[i3][9]==n)
              C[i1*N2+i3*N+i2] = C[i1*N2+i2*N+i3];
          }
        }

       //3-atom
        for (int p=n+1;p<natoms;p++)
        if (p!=m)
        {
          for (int i2=0;i2<N;i2++)
          if (basis[i2][9]==n)
          {
            for (int i3=0;i3<N;i3++)
            if (basis[i3][9]==p)
              C[i1*N2+i3*N+i2] = C[i1*N2+i2*N+i3];
          }
        }

      }
    }
  }
  return;
}

void copy_symm(int natoms, int N, int Naux, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, double* C, int type)
{
 //copy symmetric 2-atom terms
  int N2 = N*N;

  for (int m=0;m<natoms;m++)
  {
    for (int i1=0;i1<Naux;i1++)
    if (basis_aux[i1][9]==m)
    {
      for (int n=0;n<natoms;n++)
      if (m!=n)
      {
       //2-atom
        if (type==1)
        {
          for (int i2=0;i2<N;i2++)
          if (basis[i2][9]==n)
          {
            for (int i3=0;i3<N;i3++)
            if (basis[i3][9]==m)
              C[i1*N2+i3*N+i2] = C[i1*N2+i2*N+i3];
          }
        }
        else
        {
          for (int i2=0;i2<N;i2++)
          if (basis[i2][9]==m)
          {
            for (int i3=0;i3<N;i3++)
            if (basis[i3][9]==n)
              C[i1*N2+i3*N+i2] = C[i1*N2+i2*N+i3];
          }
        }

       //3-atom
        for (int p=n+1;p<natoms;p++)
        if (p!=m)
        {
          for (int i2=0;i2<N;i2++)
          if (basis[i2][9]==n)
          {
            for (int i3=0;i3<N;i3++)
            if (basis[i3][9]==p)
              C[i1*N2+i3*N+i2] = C[i1*N2+i2*N+i3];
          }
        }

      }
    }
  }
  return;
}

void copy_symm_3c_ps(int natoms, int N, int Naux, int* n2i, int* na2i, double* C)
{
  int N2 = N*N;
  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];
    int s3 = s1; int s4 = s2;
    int s5 = 0; if (m>0) s5 = na2i[m-1];  int s6 = na2i[m];

    for (int n=0;n<natoms;n++)
    if (n!=m)
    {
      s3 = 0; if (n>0) s3 = n2i[n-1]; s4 = n2i[n];
      s5 = 0; if (n>0) s5 = na2i[n-1]; s6 = na2i[n];

      for (int i1=s5;i1<s6;i1++)
      for (int i2=s1;i2<s2;i2++)
      for (int i3=s3;i3<s4;i3++)
        C[i1*N2+i3*N+i2] = C[i1*N2+i2*N+i3];

     //nothing to do with second set
      s3 = s1; s4 = s2;
    }
  }

  for (int m=0;m<natoms;m++)
  {
    for (int n=m+1;n<natoms;n++)
    {
      int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      for (int p=0;p<natoms;p++)
      if (p!=m && p!=n)
      {
        int s5 = 0; if (p>0) s5 = na2i[p-1]; int s6 = na2i[p];

        for (int i1=s5;i1<s6;i1++)
        for (int i2=s1;i2<s2;i2++)
        for (int i3=s3;i3<s4;i3++)
          C[i1*N2+i3*N+i2] = C[i1*N2+i2*N+i3];
      }
    }
  }

  return;
}

void copy_symm_4c_ps(int natoms, int* n2i, int N, double* olp)
{
  int N2 = N*N;
  int N3 = N2*N;
  int N4 = N3*N;

  for (int m=0;m<natoms;m++)
  {
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];
    for (int n=m+1;n<natoms;n++)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

     #pragma acc parallel loop collapse(4) present(olp[0:N4])
      for (int i=s1;i<s2;i++)
      for (int j=s3;j<s4;j++)
      for (int k=s3;k<s4;k++)
      for (int l=s3;l<s4;l++)
      {
        double v1 = olp[i*N3+j*N2+k*N+l];
        olp[j*N3+i*N2+k*N+l] = v1;
        olp[k*N3+l*N2+i*N+j] = v1;
        olp[k*N3+l*N2+j*N+i] = v1;
      }
     #pragma acc parallel loop collapse(4) present(olp[0:N4])
      for (int i=s1;i<s2;i++)
      for (int j=s1;j<s2;j++)
      for (int k=s3;k<s4;k++)
      for (int l=s3;l<s4;l++)
      {
        double v1 = olp[i*N3+j*N2+k*N+l];
        olp[k*N3+l*N2+i*N+j] = v1;

        olp[i*N3+l*N2+k*N+j] = v1;
        olp[l*N3+i*N2+k*N+j] = v1;
        olp[i*N3+l*N2+j*N+k] = v1;
        olp[l*N3+i*N2+j*N+k] = v1;
      }
     #pragma acc parallel loop collapse(4) present(olp[0:N4])
      for (int i=s1;i<s2;i++)
      for (int j=s1;j<s2;j++)
      for (int k=s1;k<s2;k++)
      for (int l=s3;l<s4;l++)
      {
        double v1 = olp[i*N3+j*N2+k*N+l];
        olp[l*N3+i*N2+j*N+k] = v1;
        olp[i*N3+l*N2+j*N+k] = v1;
        olp[i*N3+j*N2+l*N+k] = v1;
      }
    } //2-center cases

    for (int n=m+1;n<natoms;n++)
    {
      for (int p=n+1;p<natoms;p++)
      if (p!=m)
      {
        int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];
        int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];
        int s5 = 0; if (p>0) s5 = n2i[p-1]; int s6 = n2i[p];

       #pragma acc parallel loop collapse(4) present(olp[0:N4])
        for (int i=s1;i<s2;i++)
        for (int j=s1;j<s2;j++)
        for (int k=s3;k<s4;k++)
        for (int l=s5;l<s6;l++)
        {
          double v1 = olp[i*N3+j*N2+k*N+l];
          olp[i*N3+j*N2+l*N+k] = v1;

          olp[i*N3+k*N2+j*N+l] = v1;
          olp[i*N3+l*N2+j*N+k] = v1;

          olp[i*N3+k*N2+l*N+j] = v1;
          olp[i*N3+l*N2+k*N+j] = v1;

          olp[k*N3+i*N2+j*N+l] = v1;
          olp[l*N3+i*N2+j*N+k] = v1;

          olp[k*N3+l*N2+i*N+j] = v1;
          olp[l*N3+k*N2+i*N+j] = v1;

          olp[k*N3+i*N2+l*N+j] = v1;
          olp[l*N3+i*N2+k*N+j] = v1;
        }

       #pragma acc parallel loop collapse(4) present(olp[0:N4])
        for (int i=s1;i<s2;i++)
        for (int j=s3;j<s4;j++)
        for (int k=s3;k<s4;k++)
        for (int l=s5;l<s6;l++)
        {
          double v1 = olp[i*N3+j*N2+k*N+l];
          olp[l*N3+j*N2+k*N+i] = v1;

          olp[j*N3+k*N2+i*N+l] = v1;
          olp[j*N3+k*N2+l*N+i] = v1;

          olp[j*N3+i*N2+k*N+l] = v1;
          olp[j*N3+l*N2+k*N+i] = v1;

          olp[j*N3+i*N2+l*N+k] = v1;
          olp[j*N3+l*N2+i*N+k] = v1;

          olp[i*N3+l*N2+j*N+k] = v1;
          olp[l*N3+i*N2+j*N+k] = v1;

          olp[i*N3+j*N2+l*N+k] = v1;
          olp[l*N3+j*N2+i*N+k] = v1;
        }

       #pragma acc parallel loop collapse(4) present(olp[0:N4])
        for (int i=s1;i<s2;i++)
        for (int j=s3;j<s4;j++)
        for (int k=s5;k<s6;k++)
        for (int l=s5;l<s6;l++)
        {
          double v1 = olp[i*N3+j*N2+k*N+l];
          olp[j*N3+i*N2+k*N+l] = v1;

          olp[k*N3+l*N2+i*N+j] = v1;
          olp[k*N3+l*N2+j*N+i] = v1;

          olp[k*N3+i*N2+l*N+j] = v1;
          olp[k*N3+j*N2+l*N+i] = v1;

          olp[k*N3+i*N2+j*N+l] = v1;
          olp[k*N3+j*N2+i*N+l] = v1;

          olp[i*N3+k*N2+l*N+j] = v1;
          olp[j*N3+k*N2+l*N+i] = v1;

          olp[i*N3+k*N2+j*N+l] = v1;
          olp[j*N3+k*N2+i*N+l] = v1;
        }
      }
    } //3-center cases

   //4-center
    for (int n=m+1;n<natoms;n++)
    {
      for (int p=n+1;p<natoms;p++)
      if (p!=m)
      for (int q=p+1;q<natoms;q++)
      if (q!=m && q!=n)
      {
        int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];
        int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];
        int s5 = 0; if (p>0) s5 = n2i[p-1]; int s6 = n2i[p];
        int s7 = 0; if (q>0) s7 = n2i[q-1]; int s8 = n2i[q];

       #pragma acc parallel loop collapse(4) present(olp[0:N4])
        for (int i=s1;i<s2;i++)
        for (int j=s3;j<s4;j++)
        for (int k=s5;k<s6;k++)
        for (int l=s7;l<s8;l++)
        {
          double v1 = olp[i*N3+j*N2+k*N+l];
          //printf("   D: %i%i%i%i: %8.3e \n",i,j,k,l,v1);
          olp[i*N3+j*N2+l*N+k] = v1;
          olp[i*N3+k*N2+j*N+l] = v1;
          olp[i*N3+k*N2+l*N+j] = v1;
          olp[i*N3+l*N2+j*N+k] = v1;
          olp[i*N3+l*N2+k*N+j] = v1;

          olp[j*N3+i*N2+k*N+l] = v1;
          olp[j*N3+i*N2+l*N+k] = v1;
          olp[j*N3+k*N2+i*N+l] = v1;
          olp[j*N3+k*N2+l*N+i] = v1;
          olp[j*N3+l*N2+i*N+k] = v1;
          olp[j*N3+l*N2+k*N+i] = v1;

          olp[k*N3+i*N2+j*N+l] = v1;
          olp[k*N3+i*N2+l*N+j] = v1;
          olp[k*N3+j*N2+i*N+l] = v1;
          olp[k*N3+j*N2+l*N+i] = v1;
          olp[k*N3+l*N2+i*N+j] = v1;
          olp[k*N3+l*N2+j*N+i] = v1;

          olp[l*N3+i*N2+k*N+j] = v1;
          olp[l*N3+i*N2+j*N+k] = v1;
          olp[l*N3+k*N2+i*N+j] = v1;
          olp[l*N3+k*N2+j*N+i] = v1;
          olp[l*N3+j*N2+i*N+k] = v1;
          olp[l*N3+j*N2+k*N+i] = v1;
        }
      }
    } //4-center cases
  }

  return;
}

void copy_symm_4c_ps_cpu(int natoms, int* n2i, int N, double* olp)
{
  int N2 = N*N;
  int N3 = N2*N;
  //int N4 = N3*N;

  for (int m=0;m<natoms;m++)
  {
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];
    for (int n=m+1;n<natoms;n++)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      for (int i=s1;i<s2;i++)
      for (int j=s3;j<s4;j++)
      for (int k=s3;k<s4;k++)
      for (int l=s3;l<s4;l++)
      {
        double v1 = olp[i*N3+j*N2+k*N+l];
        olp[j*N3+i*N2+k*N+l] = v1;
        olp[k*N3+l*N2+i*N+j] = v1;
        olp[k*N3+l*N2+j*N+i] = v1;
      }
      for (int i=s1;i<s2;i++)
      for (int j=s1;j<s2;j++)
      for (int k=s3;k<s4;k++)
      for (int l=s3;l<s4;l++)
      {
        double v1 = olp[i*N3+j*N2+k*N+l];
        olp[k*N3+l*N2+i*N+j] = v1;

        olp[i*N3+l*N2+k*N+j] = v1;
        olp[l*N3+i*N2+k*N+j] = v1;
        olp[i*N3+l*N2+j*N+k] = v1;
        olp[l*N3+i*N2+j*N+k] = v1;
      }
      for (int i=s1;i<s2;i++)
      for (int j=s1;j<s2;j++)
      for (int k=s1;k<s2;k++)
      for (int l=s3;l<s4;l++)
      {
        double v1 = olp[i*N3+j*N2+k*N+l];
        olp[l*N3+i*N2+j*N+k] = v1;
        olp[i*N3+l*N2+j*N+k] = v1;
        olp[i*N3+j*N2+l*N+k] = v1;
      }
    } //2-center cases

    for (int n=m+1;n<natoms;n++)
    {
      for (int p=n+1;p<natoms;p++)
      if (p!=m)
      {
        int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];
        int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];
        int s5 = 0; if (p>0) s5 = n2i[p-1]; int s6 = n2i[p];

        for (int i=s1;i<s2;i++)
        for (int j=s1;j<s2;j++)
        for (int k=s3;k<s4;k++)
        for (int l=s5;l<s6;l++)
        {
          double v1 = olp[i*N3+j*N2+k*N+l];
          olp[i*N3+j*N2+l*N+k] = v1;

          olp[i*N3+k*N2+j*N+l] = v1;
          olp[i*N3+l*N2+j*N+k] = v1;

          olp[i*N3+k*N2+l*N+j] = v1;
          olp[i*N3+l*N2+k*N+j] = v1;

          olp[k*N3+i*N2+j*N+l] = v1;
          olp[l*N3+i*N2+j*N+k] = v1;

          olp[k*N3+l*N2+i*N+j] = v1;
          olp[l*N3+k*N2+i*N+j] = v1;

          olp[k*N3+i*N2+l*N+j] = v1;
          olp[l*N3+i*N2+k*N+j] = v1;
        }

        for (int i=s1;i<s2;i++)
        for (int j=s3;j<s4;j++)
        for (int k=s3;k<s4;k++)
        for (int l=s5;l<s6;l++)
        {
          double v1 = olp[i*N3+j*N2+k*N+l];
          olp[l*N3+j*N2+k*N+i] = v1;

          olp[j*N3+k*N2+i*N+l] = v1;
          olp[j*N3+k*N2+l*N+i] = v1;

          olp[j*N3+i*N2+k*N+l] = v1;
          olp[j*N3+l*N2+k*N+i] = v1;

          olp[j*N3+i*N2+l*N+k] = v1;
          olp[j*N3+l*N2+i*N+k] = v1;

          olp[i*N3+l*N2+j*N+k] = v1;
          olp[l*N3+i*N2+j*N+k] = v1;

          olp[i*N3+j*N2+l*N+k] = v1;
          olp[l*N3+j*N2+i*N+k] = v1;
        }

        for (int i=s1;i<s2;i++)
        for (int j=s3;j<s4;j++)
        for (int k=s5;k<s6;k++)
        for (int l=s5;l<s6;l++)
        {
          double v1 = olp[i*N3+j*N2+k*N+l];
          olp[j*N3+i*N2+k*N+l] = v1;

          olp[k*N3+l*N2+i*N+j] = v1;
          olp[k*N3+l*N2+j*N+i] = v1;

          olp[k*N3+i*N2+l*N+j] = v1;
          olp[k*N3+j*N2+l*N+i] = v1;

          olp[k*N3+i*N2+j*N+l] = v1;
          olp[k*N3+j*N2+i*N+l] = v1;

          olp[i*N3+k*N2+l*N+j] = v1;
          olp[j*N3+k*N2+l*N+i] = v1;

          olp[i*N3+k*N2+j*N+l] = v1;
          olp[j*N3+k*N2+i*N+l] = v1;
        }
      }
    } //3-center cases

   //4-center
    for (int n=m+1;n<natoms;n++)
    {
      for (int p=n+1;p<natoms;p++)
      if (p!=m)
      for (int q=p+1;q<natoms;q++)
      if (q!=m && q!=n)
      {
        int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];
        int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];
        int s5 = 0; if (p>0) s5 = n2i[p-1]; int s6 = n2i[p];
        int s7 = 0; if (q>0) s7 = n2i[q-1]; int s8 = n2i[q];

        for (int i=s1;i<s2;i++)
        for (int j=s3;j<s4;j++)
        for (int k=s5;k<s6;k++)
        for (int l=s7;l<s8;l++)
        {
          double v1 = olp[i*N3+j*N2+k*N+l];
          //printf("   D: %i%i%i%i: %8.3e \n",i,j,k,l,v1);
          olp[i*N3+j*N2+l*N+k] = v1;
          olp[i*N3+k*N2+j*N+l] = v1;
          olp[i*N3+k*N2+l*N+j] = v1;
          olp[i*N3+l*N2+j*N+k] = v1;
          olp[i*N3+l*N2+k*N+j] = v1;

          olp[j*N3+i*N2+k*N+l] = v1;
          olp[j*N3+i*N2+l*N+k] = v1;
          olp[j*N3+k*N2+i*N+l] = v1;
          olp[j*N3+k*N2+l*N+i] = v1;
          olp[j*N3+l*N2+i*N+k] = v1;
          olp[j*N3+l*N2+k*N+i] = v1;

          olp[k*N3+i*N2+j*N+l] = v1;
          olp[k*N3+i*N2+l*N+j] = v1;
          olp[k*N3+j*N2+i*N+l] = v1;
          olp[k*N3+j*N2+l*N+i] = v1;
          olp[k*N3+l*N2+i*N+j] = v1;
          olp[k*N3+l*N2+j*N+i] = v1;

          olp[l*N3+i*N2+k*N+j] = v1;
          olp[l*N3+i*N2+j*N+k] = v1;
          olp[l*N3+k*N2+i*N+j] = v1;
          olp[l*N3+k*N2+j*N+i] = v1;
          olp[l*N3+j*N2+i*N+k] = v1;
          olp[l*N3+j*N2+k*N+i] = v1;
        }
      }
    } //4-center cases
  }

  return;
}

int get_imax_n2ip(int Nmax, int natoms, int N, vector<vector<double> >& basis, vector<vector<int> >& n2ip)
{
  n2ip.clear();

  int n2i[natoms];
  int iN0 = get_imax_n2i(natoms,N,basis,n2i);

  int imaxN = 0;
  for (int n=0;n<natoms;n++)
  {
    int s1 = 0; if (n>0) s1 = n2i[n-1]; int s2 = n2i[n];

    int size = s2-s1;
    imaxN = max(size,imaxN);
    //printf("  working on s12: %3i %3i  size: %3i \n",s1,s2,size);

    if (size<Nmax)
    {
      vector<int> p1;
      p1.push_back(s1);
      p1.push_back(s2);
      n2ip.push_back(p1);
    }
    else
    {
      vector<int> p1;
      p1.push_back(s1);

      while (size>Nmax)
      {
        s1 += Nmax;
        p1.push_back(s1);
        size -= Nmax;
      }
      p1.push_back(s2);

      n2ip.push_back(p1);
    }
  }
  imaxN = min(Nmax,imaxN);

 #if 0
  printf("\n     n2ip: \n");
  for (int n=0;n<natoms;n++)
  {
    printf("      ");
    for (int j=0;j<n2ip[n].size();j++)
      printf("  %2i",n2ip[n][j]);
    printf("\n");
  }
 #endif

  return imaxN;
}

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
  n2i[wa] = N;

  for (int n=wa;n<natoms;n++)
    n2i[n] = N;

  return imaxN;
}

int get_natoms_with_basis(int natoms, int* atno, vector<vector<double> >& basis)
{
  int N = basis.size();

  int natoms1 = 0;
  for (int n=0;n<natoms;n++)
  {
    if (atno[n]==0)
    {
      bool found = 0;
      for (int j=0;j<N;j++)
      if (basis[j][9]==n)
      { found = 1; break ;}
      if (!found)
        break;
    }
    else
      natoms1 = n+1;
  }
  
  return natoms1;
}

vector<vector<double> > setup_integrals_gsgpu(vector<vector<double> >& basis_aux, int natoms, int* atno, double* coords, int& nbas, int& nenv, int& N, int& Naux, int& nbas_ri, int* &atm, int* &bas, double* &env, int prl)
{
  string basfile = "basis";
  string auxfile = "aux";

  BT::DO_CART = false;
  bool do_ri = true;
  CINTPrep inp(do_ri);
  //inp.read_xyz(xyzfile);
  inp.assign_coords(natoms,atno,coords,true);
  inp.read_bas(basfile);
  bool found_ri = inp.read_bas_ri(auxfile);
  if (!found_ri) inp.do_ri = false;

  inp.prep_env();

  nbas = inp.get_nbas();
  nenv = inp.get_nenv();
  N = inp.get_var_dim();
  if (found_ri)
  {
    Naux = inp.get_var_dim_ri();
    nbas_ri = inp.get_nbas_ri();
  }
  else
    nbas_ri = Naux = 0;

  if (prl>-1) printf("\n  setup_integrals_g. natoms: %2i nbas: %2i nenv: %2i N: %2i Naux: %2i nbas_ri: %2i \n",natoms,nbas,nenv,N,Naux,nbas_ri);

  atm = inp.get_atm();
  bas = inp.get_bas();
  env = inp.get_env();

  map< int, basis_t > basmap = inp.basmap;

  vector<vector<double> > basis;
  for (int n=0;n<natoms;n++)
  {
    int i1 = atno[n];
    basis_t basis1 = basmap[i1];
    int Z = basis1.nuc;
    int nshell = basis1.shells.size();
    int nbasat = inp.anum_to_N[Z];
    if (prl>0)
      printf("  atom %2i has atomic num %2i with %2i shells and %3i basis \n",n+1,Z,nshell,nbasat);

    if (prl>2)
    for (int j=0;j<nshell;j++)
    {
      int size1 = basis1.exps[j].size();
      int l1 = basis1.shells[j];
      printf("   shell %i has %2i Gaussians and l1: %i \n",j,size1,l1);
      for (int k=0;k<size1;k++)
        printf("   zeta:  %8.5f  norm:  %8.5f / %8.5f \n",basis1.exps[j][k],basis1.coef[j][k],get_gto_norm(l1,basis1.exps[j][k]));
    }
  }
  printf("\n");

  double n0 = norm_sh(0,0);
  for (int n=0;n<natoms;n++)
  {
    int i1 = atno[n];
    basis_t basis1 = basmap[i1];
    int Z = basis1.nuc;
    int nshell = basis1.shells.size();
    //int nbasat = inp.anum_to_N[Z];

    for (int j=0;j<nshell;j++)
    {
      int size1 = basis1.exps[j].size();
      int l1 = basis1.shells[j];
      int lmin = -l1; int lmax = l1;

      for (int m=lmin;m<=lmax;m++)
      {
        vector<double> b1; for (int p=0;p<10;p++) b1.push_back(0);
        b1[0] = l1+1; b1[1] = l1; b1[2] = m;
       //using b1[3]==b1[4] to store # of gaussians
        b1[3] = b1[4] = size1;
        b1[5] = coords[3*n+0]; b1[6] = coords[3*n+1]; b1[7] = coords[3*n+2];
        b1[8] = Z; b1[9] = n;

       //zeta then normalization coeffs
        for (int k=0;k<size1;k++)
          b1.push_back(basis1.exps[j][k]);
        for (int k=0;k<size1;k++)
        {
          double norm1 = get_gto_norm(l1,basis1.exps[j][k])*basis1.coef[j][k];
          if (l1>0) norm1 *= norm_sh(l1,m)/n0;
          b1.push_back(norm1);
        }

        basis.push_back(b1);
      }
    }
  }

  if (prl>1)
  for (int i=0;i<basis.size();i++)
  {
    vector<double> basis1 = basis[i];
    int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; int ng = basis1[3];
    double x1 = basis1[5]; double y1 = basis1[6]; double z1 = basis1[7]; double Z1 = basis1[8];
    printf("  basis[%2i] nlm:  %i  %i %2i (XYZ: %8.5f %8.5f %8.5f Z: %3.2f) ng: %2i \n",i,n1,l1,m1,x1,y1,z1,Z1,ng);
    for (int j=0;j<ng;j++)
    {
      printf("   %8.5f  %8.5f \n",basis1[10+j],basis1[10+ng+j]);
    }
  }

 //just to make basis_aux's size correct
  basis_aux.clear();
  for (int i=0;i<Naux;i++)
  {
    vector<double> b1; b1.push_back(1);
    basis_aux.push_back(b1);
  }

  return basis;
}

void compute_integrals_g(int natm, int nbas, int nenv, int N, int Naux, int nbas_ri, int* atm, int* bas, double* env, double* S, double* T, double* jH1, double* A, double* C, double* pvp, int prl)
{
  get_overlap(S, N, natm, nbas, nenv, atm, bas, env);
  get_hcore(jH1, N, natm, nbas, nenv, atm, bas, env);
  if (T!=NULL)
    get_tcore(T, N, natm, nbas, nenv, atm, bas, env);
  printf("before gen_pvp\n");
  gen_pvp(pvp, N, natm, nbas, nenv, atm, bas, env);
  printf("after gen_pvp\n");

  if (prl>1)
  {
    printf("\n S: \n");
    print_square(N,S);
    printf("\n jH1: \n");
    print_square(N,jH1);
  }

  if (Naux>0)
  {
    gen_eri_2c(A, Naux, natm, nbas, nenv, nbas_ri, atm, bas, env);
    gen_eri_3c(C, N, Naux, natm, nbas, nenv, nbas_ri, atm, bas, env);
  }
  else
  {
    printf("  no auxiliary basis \n");
  }

  if (prl>1)
  {
    printf("\n A: \n");
    for (int m=0;m<Naux;m++)
    {
      for (int n=0;n<Naux;n++)
        printf(" %8.5f",A[m*Naux+n]);
      printf("\n");
    }

    int N2 = N*N;
    printf("\n C: \n");
    for (int m=0;m<N2;m++)
    {
      for (int n=0;n<Naux;n++)
        printf(" %8.5f",C[m*Naux+n]);
      printf("\n");
    }
  }

  return;
}
void compute_gaussian_integrals_to_disk(int N, int Naux, int natoms, int nbas, int nenv, int nbas_ri, int* atm, int* bas, double* env)
{
  int N2 = N*N;
  int Na2 = Naux*Naux;
  int N2a = Naux*N2;

  double* S = new double[N2]; double* jH1 = new double[N2]; double* A = new double[Na2]; double* C = new double[N2a];
  double* T = new double[N2]; double* En = new double[N2];
  double* pvp = new double[N2];
  for (int m=0;m<N2;m++) En[m] = 0.;
  
  printf("before compute_integrals_g\n");
  compute_integrals_g(natoms,nbas,nenv,N,Naux,nbas_ri,atm,bas,env,S,NULL,jH1,A,C,pvp,1);
  printf("after compute_integrals_g\n");

 //note we don't have separate En and T, so using sum of the two
  for (int m=0;m<N2;m++) T[m] = jH1[m];
  write_S_En_T(N,S,En,T);
  write_square(N,pvp,"pvp",2);
  write_square(Naux,A,"A",2);
  write_C(Naux,N2,C);

  return;
}
