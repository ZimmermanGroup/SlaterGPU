#include "integrals.h"
#include "braggslater.h"
#include <chrono>
#define TIMER_KD 1

#define TEST_SORT 0
#include "fp_def.h"

/*
 //current status of compute_all_2/3c code:
  1. add n=5 functions (need Inr)
  2. optimization of reduce functions

  4. need to check 4s,4p,4d functions for accuracy
  5. Cartesian d functions


 //completed
  1. implement d+f functions (done)
  1. b. implement g+h functions (done)
  2. need En/pVp --> 3c integrations (done)
  3. 4f functions have some kind of bug (fixed)
  4. 5g functions m==0 component has bug (high angular momentum sensitive to small weights in grid, fixed)
  5. 3c integration
  6. 2s, 3s, 3p functions + derivatives
  7. 4s, 4p, 4d functions

*/
using namespace std;

#include <string>
// void auto_crash();
// void print_duration(chrono::high_resolution_clock::time_point t1, chrono::high_resolution_clock::time_point t2, string name);

void print_grid(int gs, FP1* grid1, FP1* grid2, FP1* wt1, FP1* wt2, int prl)
{
  if (gs>300) return;
  if (prl>0)
  {
    //printf("\n printing grid1 \n");
    for (int i=0;i<gs;i++)
      //printf(" %8.5f %8.5f %8.5f  wt: %10.6f  r1/2: %5.2f %5.2f \n",grid1[6*i+0],grid1[6*i+1],grid1[6*i+2],wt1[i],grid1[6*i+3],grid1[6*i+4]);
      printf(" %9.5f %9.5f %9.5f  wt: %7.4e  r1/2: %5.2f %5.2f \n",grid1[6*i+0],grid1[6*i+1],grid1[6*i+2],wt1[i],grid1[6*i+3],grid1[6*i+4]);
    if (grid2!=NULL && wt2!=NULL)
    {
      //printf("\n printing grid2 \n");
      for (int i=0;i<gs;i++)
        //printf(" %8.5f %8.5f %8.5f  wt: %10.6f  r1/2: %5.2f %5.2f \n",grid2[6*i+0],grid2[6*i+1],grid2[6*i+2],wt2[i],grid2[6*i+3],grid2[6*i+4]);
        printf(" %9.5f %9.5f %9.5f  wt: %7.4e  r1/2: %5.2f %5.2f \n",grid2[6*i+0],grid2[6*i+1],grid2[6*i+2],wt2[i],grid2[6*i+3],grid2[6*i+4]);
    }
  }
  return;
}

void print_grid(int gs, FP1* grid1, FP1* grid2, FP1* grid3, FP1* wt1, FP1* wt2, FP1* wt3, int prl)
{
  if (gs>300) return;
  if (prl>0)
  {
    //printf("\n printing grid1 \n");
    for (int i=0;i<gs;i++)
      //printf(" %8.5f %8.5f %8.5f  wt: %10.6f  r123: %5.2f %5.2f %5.2f \n",grid1[6*i+0],grid1[6*i+1],grid1[6*i+2],wt1[i],grid1[6*i+3],grid1[6*i+4],grid1[6*i+5]);
      printf(" %9.5f %9.5f %9.5f  wt: %7.4e  r123: %5.2f %5.2f %5.2f \n",grid1[6*i+0],grid1[6*i+1],grid1[6*i+2],wt1[i],grid1[6*i+3],grid1[6*i+4],grid1[6*i+5]);
    if (grid2!=NULL && wt2!=NULL)
    {
      //printf("\n printing grid2 \n");
      for (int i=0;i<gs;i++)
        //printf(" %8.5f %8.5f %8.5f  wt: %10.6f  r123: %5.2f %5.2f %5.2f \n",grid2[6*i+0],grid2[6*i+1],grid2[6*i+2],wt2[i],grid2[6*i+3],grid2[6*i+4],grid2[6*i+5]);
        printf(" %9.5f %9.5f %9.5f  wt: %7.4e  r123: %5.2f %5.2f %5.2f \n",grid2[6*i+0],grid2[6*i+1],grid2[6*i+2],wt2[i],grid2[6*i+3],grid2[6*i+4],grid2[6*i+5]);
    }
    if (grid3!=NULL && wt3!=NULL)
    {
      //printf("\n printing grid3 \n");
      for (int i=0;i<gs;i++)
        //printf(" %8.5f %8.5f %8.5f  wt: %10.6f  r123: %5.2f %5.2f %5.2f \n",grid3[6*i+0],grid3[6*i+1],grid3[6*i+2],wt3[i],grid3[6*i+3],grid3[6*i+4],grid3[6*i+5]);
        printf(" %9.5f %9.5f %9.5f  wt: %7.4e  r123: %5.2f %5.2f %5.2f \n",grid3[6*i+0],grid3[6*i+1],grid3[6*i+2],wt3[i],grid3[6*i+3],grid3[6*i+4],grid3[6*i+5]);
    }
  }
  return;
}

void print_array(int size, FP1* vec)
{
  for (int i=0;i<size;i++)
    printf(" %6.3e",vec[i]);
  printf("\n");
}

void clean_small_values(int N, FP1* S) 
{
  int N2 = N*N;
 #pragma acc parallel loop independent present(S[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i+1;j<N;j++)
  if (fabs(S[i*N+j])<CL_THRESH)
    S[i*N+j] = S[j*N+i] = 0.;
}

// if !EVL64 commented out by Soumi
//#if !EVL64
void clean_small_values(int N, FP2* S)
{
  int N2 = N*N;
 #pragma acc parallel loop independent present(S[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i+1;j<N;j++)
  if (fabs(S[i*N+j])<CL_THRESH)
    S[i*N+j] = S[j*N+i] = 0.;
}
//#endif

void acc_assign(int size, FP1* vec, FP1 v1) 
{
#if USE_ACC
 #pragma acc parallel loop independent present(vec[0:size])
#endif
  for (int j=0;j<size;j++)
    vec[j] = v1;
}

// if !EVL64 commented out by Soumi
//#if !EVL64
void acc_assign(int size, FP2* vec, FP2 v1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(vec[0:size])
#endif
  for (int j=0;j<size;j++)
    vec[j] = v1;
}
//#endif

void acc_assign(int size, FP1* vec1, FP1* vec2, FP1 v1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(vec1[0:size],vec2[0:size])
#endif
  for (int j=0;j<size;j++)
    vec1[j] = vec2[j] = v1;
}

void acc_assign(int size, FP1* vec1, FP1* vec2, FP1* vec3, FP1 v1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(vec1[0:size],vec2[0:size],vec3[0:size])
#endif
  for (int j=0;j<size;j++)
    vec1[j] = vec2[j] = vec3[j] = v1;
}

FP1 acc_sum(int size, FP1* vec)
{
  FP1 sum = 0.;
#if USE_ACC
 #pragma acc parallel loop present(vec[0:size]) reduction(+:sum)
#endif
  for (int i=0;i<size;i++)
    sum += vec[i];

  return sum;
}

void acc_copyf(int tid, int size, FP1* v1, FP1* v2)
{
 #if USE_ACC
  #pragma acc parallel loop independent present(v1[0:size],v2[0:size]) //async(tid)
 #endif
  for (int i=0;i<size;i++)
    v1[i] = v2[i];

  return;
}

void acc_copyf(int size, FP1* v1, FP1* v2)
{
 #if USE_ACC
  #pragma acc parallel loop independent present(v1[0:size],v2[0:size])
 #endif
  for (int i=0;i<size;i++)
    v1[i] = v2[i];

  return;
}

void acc_copyf(int size, FP1* v1, FP1* v2, FP1* v3, FP1* v4)
{
 #if USE_ACC
  #pragma acc parallel loop independent present(v1[0:size],v2[0:size],v3[0:size],v4[0:size])
 #endif
  for (int i=0;i<size;i++)
  {
    v1[i] = v2[i];
    v3[i] = v4[i];
  }

  return;
}

void acc_copyf(int size, FP1* v1, FP1* v2, FP1* v3, FP1* v4, FP1* v5, FP1* v6)
{
 #if USE_ACC
  #pragma acc parallel loop independent present(v1[0:size],v2[0:size],v3[0:size],v4[0:size],v5[0:size],v6[0:size])
 #endif
  for (int i=0;i<size;i++)
  {
    v1[i] = v2[i];
    v3[i] = v4[i];
    v5[i] = v6[i];
  }

  return;
}

void eliminate_small_wt_3(int size, FP1* wt1, FP1* wt2, FP1* wt3)
{
#if USE_ACC
 #pragma acc parallel loop independent present(wt1[0:size],wt2[0:size],wt3[0:size])
#endif
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

void eliminate_small_wt_3(int s1, int size, FP1* wt1, FP1* wt2, FP1* wt3)
{
#if USE_ACC
 #pragma acc parallel loop independent present(wt1[0:size],wt2[0:size],wt3[0:size])
#endif
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

void eliminate_small_wt(int size, FP1* wt1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(wt1[0:size])
#endif
  for (int i=0;i<size;i++) 
  if (wt1[i]<WT_THRESH)
  {
    //printf(" small wt: %12.10f \n",wt1[i]);
    wt1[i] = 0.;
  }

}

void eliminate_small_wt(int s1, int size, FP1* wt1)
{
#if USE_ACC
 #pragma acc parallel loop independent present(wt1[0:size])
#endif
  for (int i=s1;i<size;i++) 
  if (wt1[i]<WT_THRESH_D)
  {
    //printf(" small wt: %12.10f \n",wt1[i]);
    wt1[i] = 0.;
  }

}

void copy_grid(int gs, FP1* grid1, FP1* grid2)
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

void copy_grid(int gs, FP1* grid1, FP1* wt1, FP1* grid2, FP1* wt2)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid1[0:6*gs],wt1[0:gs],grid2[0:6*gs],wt2[0:gs])
#endif
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


void recenter_grid_zero(int gs, FP1* grid, FP1 x2, FP1 y2, FP1 z2)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    grid[6*i+0] += x2;
    grid[6*i+1] += y2;
    grid[6*i+2] += z2;
    FP1 xn = grid[6*i+0];
    FP1 yn = grid[6*i+1];
    FP1 zn = grid[6*i+2];

    FP1 r1 = sqrtf(xn*xn+yn*yn+zn*zn);
    grid[6*i+3] = r1;
  }

  return;
}

void recenter_grid(int gs, FP1* grid, FP1 x2, FP1 y2, FP1 z2)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    grid[6*i+0] += x2;
    grid[6*i+1] += y2;
    grid[6*i+2] += z2;
    FP1 xn = grid[6*i+0];
    FP1 yn = grid[6*i+1];
    FP1 zn = grid[6*i+2];

   //distance to other center
    FP1 r2 = sqrtf(xn*xn+yn*yn+zn*zn);
    grid[6*i+4] = r2;
  }

  return;
}

void recenter_grid_exp(int gs, FP1* grid, FP1* wt, FP1* val, FP1 x2, FP1 y2, FP1 z2, FP1 zeta2)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid[0:6*gs],val[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    grid[6*i+0] += x2;
    grid[6*i+1] += y2;
    grid[6*i+2] += z2;
    FP1 xn = grid[6*i+0];
    FP1 yn = grid[6*i+1];
    FP1 zn = grid[6*i+2];
    FP1 r2 = grid[6*i+3];

    FP1 r1 = sqrtf(xn*xn+yn*yn+zn*zn);
    FP1 er2 = expf(-zeta2*r2);

   //distance to first center
    grid[6*i+4] = r1;
    //wt[i] *= er2;
    val[i] = er2;
  }

  return;
}

void add_r123_to_grid(int gs, FP1* grid1, FP1 A1, FP1 B1, FP1 C1, FP1 A2, FP1 B2, FP1 C2, FP1 A3, FP1 B3, FP1 C3)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid1[0:6*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x = grid1[6*i+0];
    FP1 y = grid1[6*i+1];
    FP1 z = grid1[6*i+2];
    FP1 x1 = x-A1; FP1 y1 = y-B1; FP1 z1 = z-C1;
    FP1 x2 = x-A2; FP1 y2 = y-B2; FP1 z2 = z-C2;
    FP1 x3 = x-A3; FP1 y3 = y-B3; FP1 z3 = z-C3;
    FP1 r1 = sqrtf(x1*x1+y1*y1+z1*z1);
    FP1 r2 = sqrtf(x2*x2+y2*y2+z2*z2);
    FP1 r3 = sqrtf(x3*x3+y3*y3+z3*z3);
    grid1[6*i+3] = r1;
    grid1[6*i+4] = r2;
    grid1[6*i+5] = r3;
  }

  return;
}

void add_r1_to_grid_6z(int gs, FP1* grid1, FP1* grid2, FP1* grid3, FP1* grid4, FP1* grid5, FP1* grid6)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid1[0:6*gs],grid2[0:6*gs],grid3[0:6*gs],grid4[0:6*gs],grid5[0:6*gs],grid6[0:6*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x1 = grid1[6*i+0]; FP1 y1 = grid1[6*i+1]; FP1 z1 = grid1[6*i+2];
    FP1 r1 = sqrtf(x1*x1+y1*y1+z1*z1);
    grid1[6*i+3] = r1;

    FP1 x2 = grid2[6*i+0]; FP1 y2 = grid2[6*i+1]; FP1 z2 = grid2[6*i+2];
    FP1 r2 = sqrtf(x2*x2+y2*y2+z2*z2);
    grid2[6*i+3] = r2;

    FP1 x3 = grid3[6*i+0]; FP1 y3 = grid3[6*i+1]; FP1 z3 = grid3[6*i+2];
    FP1 r3 = sqrtf(x3*x3+y3*y3+z3*z3);
    grid3[6*i+3] = r3;

    FP1 x4 = grid4[6*i+0]; FP1 y4 = grid4[6*i+1]; FP1 z4 = grid4[6*i+2];
    FP1 r4 = sqrtf(x4*x4+y4*y4+z4*z4);
    grid4[6*i+3] = r4;

    FP1 x5 = grid5[6*i+0]; FP1 y5 = grid5[6*i+1]; FP1 z5 = grid5[6*i+2];
    FP1 r5 = sqrtf(x5*x5+y5*y5+z5*z5);
    grid5[6*i+3] = r5;

    FP1 x6 = grid6[6*i+0]; FP1 y6 = grid6[6*i+1]; FP1 z6 = grid6[6*i+2];
    FP1 r6 = sqrtf(x6*x6+y6*y6+z6*z6);
    grid6[6*i+3] = r6;
  }

  return;
}

void add_r1_to_grid(int gs, FP1* grid1, FP1 A2, FP1 B2, FP1 C2)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid1[0:6*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x1 = grid1[6*i+0];
    FP1 y1 = grid1[6*i+1];
    FP1 z1 = grid1[6*i+2];
    FP1 x12 = x1-A2;
    FP1 y12 = y1-B2;
    FP1 z12 = z1-C2;
    FP1 r2 = sqrtf(x12*x12+y12*y12+z12*z12);
    grid1[6*i+3] = r2;
  }

  return;
}

void add_r2_to_grid(int gs, FP1* grid1, FP1 A2, FP1 B2, FP1 C2)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid1[0:6*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x1 = grid1[6*i+0];
    FP1 y1 = grid1[6*i+1];
    FP1 z1 = grid1[6*i+2];
    FP1 x12 = x1-A2;
    FP1 y12 = y1-B2;
    FP1 z12 = z1-C2;
    FP1 r2 = sqrtf(x12*x12+y12*y12+z12*z12);
    grid1[6*i+4] = r2;
  }

  return;
}

void add_r3_to_grid(int gs, FP1* grid1, FP1 A3, FP1 B3, FP1 C3)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid1[0:6*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 x1 = grid1[6*i+0];
    FP1 y1 = grid1[6*i+1];
    FP1 z1 = grid1[6*i+2];
    FP1 x13 = x1-A3;
    FP1 y13 = y1-B3;
    FP1 z13 = z1-C3;
    FP1 r3 = sqrtf(x13*x13+y13*y13+z13*z13);
    grid1[6*i+5] = r3;
  }

  return;
}

int find_center_of_grid(FP1 Z1, int nrad)
{
  FP1* r = new FP1[nrad];
  FP1* w = new FP1[nrad];

 #if USE_ACC
  #pragma acc enter data create(r[0:nrad],w[0:nrad])
 #endif

  get_murak_grid_f(nrad,r,w,Z1,3);

  int s1 = 0;
#if USE_ACC
 #pragma acc parallel loop present(w[0:nrad])
#endif
  for (int i=0;i<nrad;i++)
  if (w[i]<WT_THRESH_D)
    s1 = i; 
  s1++;

  //printf(" s1: %2i \n",s1);

  #pragma acc exit data delete(r[0:nrad],w[0:nrad])

  return s1;
}

void generate_central_grid_2(FP1* grid1, FP1* wt1, FP1 Z1, int nrad, int nang, FP1* ang_g, FP1* ang_w)
{
  FP1* r = new FP1[nrad];
  FP1* w = new FP1[nrad];

 #if USE_ACC
  #pragma acc enter data create(r[0:nrad],w[0:nrad])
 #endif

  get_murak_grid_f(nrad,r,w,Z1,3);

  int gs = nrad*nang;

#if USE_ACC
 #pragma acc parallel loop independent present(r[0:nrad],w[0:nrad],ang_g[0:3*nang],ang_w[0:nang],grid1[0:6*gs],wt1[0:gs])
#endif
  for (int i=0;i<nrad;i++)
  {
    FP1 r1 = r[i];
    FP1 wr1 = w[i];

  #if USE_ACC
   #pragma acc loop independent 
  #endif
    for (int j=0;j<nang;j++)
    {
      FP1 w1 = wr1*ang_w[j];

     //grid positions
      FP1 x1 = r1 * ang_g[3*j+0];
      FP1 y1 = r1 * ang_g[3*j+1];
      FP1 z1 = r1 * ang_g[3*j+2];

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

void generate_central_grid(FP1* grid1, FP1* wt1, FP1* val1, int need_inr, FP1 Z1, int n1, int l1, FP1 zeta1, int nrad, int nang, FP1* ang_g, FP1* ang_w)
{
  printf("\n WARNING: shouldn't be here in generate_central_grid() \n");

  FP1* r = new FP1[nrad];
  FP1* w = new FP1[nrad];
  FP1* er = new FP1[nrad];
  FP1* inr = new FP1[nrad];

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
    FP1 r1 = r[i];
    FP1 wr1 = w[i];
    FP1 inr1 = inr[i];
 
  #if USE_ACC
   #pragma acc loop independent 
  #endif
    for (int j=0;j<nang;j++)
    {
      FP1 w1 = wr1*ang_w[j];

     //grid positions
      FP1 x1 = r1 * ang_g[3*j+0];
      FP1 y1 = r1 * ang_g[3*j+1];
      FP1 z1 = r1 * ang_g[3*j+2];

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

FP1 becke_a(int Z1, int Z2)
{
  FP1 x = get_radii_2_2(Z1)/get_radii_2_2(Z2);

  FP1 u = (x-1.f)/(x+1.f);
  FP1 a = u/(u*u-1.f);

  if (a>0.5f) a = 0.5f;
  else if (a<-0.5f) a = -0.5f;
  return a;
}

#pragma acc routine seq
FP1 bf3(FP1 f1)
{
  //default is order 3
  #define ORDER4 0

  f1 = 1.5f*f1 - 0.5f*f1*f1*f1;
  f1 = 1.5f*f1 - 0.5f*f1*f1*f1;
  f1 = 1.5f*f1 - 0.5f*f1*f1*f1;
 #if ORDER4
  f1 = 1.5f*f1 - 0.5f*f1*f1*f1;
 #endif
  f1 = 0.5f*(1.f-f1);
  f1 = f1*f1;

  return f1;
}

#pragma acc routine seq
FP1 get_3c_bw(FP1 r1, FP1 r2, FP1 r3, FP1 a12, FP1 a13, FP1 a23, FP1 a21, FP1 a31, FP1 a32, FP1 oR12, FP1 oR13, FP1 oR23)
{
  FP1 mu12 = (r1-r2)*oR12;
  FP1 mu13 = (r1-r3)*oR13;
  FP1 mu23 = (r2-r3)*oR23;
  FP1 omm12 = 1.f-mu12*mu12;
  FP1 omm13 = 1.f-mu13*mu13;
  FP1 omm23 = 1.f-mu23*mu23;

 //original
  FP1 nu12 =  mu12 + a12*omm12;
  FP1 nu21 = -mu12 + a21*omm12;
  FP1 nu13 =  mu13 + a13*omm13;
  FP1 nu31 = -mu13 + a31*omm13;
  FP1 nu23 =  mu23 + a23*omm23;
  FP1 nu32 = -mu23 + a32*omm23;

  FP1 f12 = bf3(nu12);
  FP1 f21 = bf3(nu21);
  FP1 f13 = bf3(nu13);
  FP1 f31 = bf3(nu31);
  FP1 f23 = bf3(nu23);
  FP1 f32 = bf3(nu32);

  FP1 norm = f12*f13 + f21*f23 + f31*f32 + 1.e-10f;
  FP1 s1 = f12*f13/norm;

  return s1;
}

void becke_weight_3c(int gs, FP1* grid1, FP1* wt1, FP1* grid2, FP1* wt2, FP1* grid3, FP1* wt3, int Z1, int Z2, int Z3, FP1 A2, FP1 B2, FP1 C2, FP1 A3, FP1 B3, FP1 C3)
{
  FP1 A23 = A3-A2; FP1 B23 = B3-B2; FP1 C23 = C3-C2;

  const FP1 R12 = sqrtf(A2*A2+B2*B2+C2*C2);
  const FP1 R13 = sqrtf(A3*A3+B3*B3+C3*C3);
  const FP1 R23 = sqrtf(A23*A23+B23*B23+C23*C23);
  const FP1 oR12 = 1./R12;
  const FP1 oR13 = 1./R13;
  const FP1 oR23 = 1./R23;

  FP1 a12 = becke_a(Z1,Z2);
  FP1 a13 = becke_a(Z1,Z3);
  FP1 a23 = becke_a(Z2,Z3);
  FP1 a21 = becke_a(Z2,Z1);
  FP1 a31 = becke_a(Z3,Z1);
  FP1 a32 = becke_a(Z3,Z2); 

#if USE_ACC
 #pragma acc parallel loop independent present(grid1[0:6*gs],wt1[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r1 = grid1[6*i+3];
    FP1 r2 = grid1[6*i+4];
    FP1 r3 = grid1[6*i+5];
    FP1 s1 = get_3c_bw(r1,r2,r3,a12,a13,a23,a21,a31,a32,oR12,oR13,oR23);

    wt1[i] *= s1;
  }

#if USE_ACC
 #pragma acc parallel loop independent present(grid2[0:6*gs],wt2[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r1 = grid2[6*i+3];
    FP1 r2 = grid2[6*i+4];
    FP1 r3 = grid2[6*i+5];
    FP1 s1 = get_3c_bw(r1,r2,r3,a21,a23,a13,a12,a32,a31,oR12,oR23,oR13);

    wt2[i] *= s1;
  }

#if USE_ACC
 #pragma acc parallel loop independent present(grid3[0:6*gs],wt3[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r1 = grid3[6*i+3];
    FP1 r2 = grid3[6*i+4];
    FP1 r3 = grid3[6*i+5];
    FP1 s1 = get_3c_bw(r1,r2,r3,a32,a31,a21,a23,a13,a12,oR23,oR13,oR12);

    wt3[i] *= s1;
  }

  return;
}

void becke_weight_2c(int gs, FP1* grid1, FP1* wt1, FP1* grid2, FP1* wt2, int Z1, int Z2, FP1 A2, FP1 B2, FP1 C2)
{
 //first center is at 0,0,0 and second at A2,B2,C2
  FP1 R  = sqrtf(A2*A2+B2*B2+C2*C2);
  const FP1 oR = 1./R;

  const FP1 a1 = becke_a(Z1,Z2);
  const FP1 a2 = becke_a(Z2,Z1);
  //const FP1 a1 = becke_a(Z2,Z1);
  //const FP1 a2 = becke_a(Z1,Z2);

  //printf(" a1/2: %8.5f %8.5f \n",a1,a2);

#if USE_ACC
 #pragma acc parallel loop independent present(grid1[0:6*gs],wt1[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r1 = grid1[6*i+3];
    FP1 r2 = grid1[6*i+4];
    FP1 mu = (r1-r2)*oR;
    FP1 omm2 = 1.f-mu*mu;
    FP1 nu1 =  mu + a1*omm2;
    FP1 nu2 = -mu + a2*omm2;

    FP1 f1 = bf3(nu1);
    FP1 f2 = bf3(nu2);

    FP1 norm = f1 + f2 + 1.e-10f;
    FP1 s1 = f1/norm;

    //printf(" r1/2: %8.5f %8.5f s1: %8.5f \n",r1,r2,s1);
    wt1[i] *= s1;
  }

#if USE_ACC
 #pragma acc parallel loop independent present(grid2[0:6*gs],wt2[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    FP1 r1 = grid2[6*i+3];
    FP1 r2 = grid2[6*i+4];
    FP1 mu = (r1-r2)*oR;
    FP1 omm2 = 1.f-mu*mu;
    FP1 nu1 =  mu + a2*omm2;
    FP1 nu2 = -mu + a1*omm2;

    FP1 f1 = bf3(nu1);
    FP1 f2 = bf3(nu2);

    FP1 norm = f1 + f2 + 1.e-10f;
    FP1 s1 = f1/norm;

    //printf(" r1/2: %8.5f %8.5f s1: %8.5f \n",r1,r2,s1);
    wt2[i] *= s1;
  }

  return;
}


void transpose_C(int Naux, int N, FP1* C) 
{
  int N2 = N*N;
  int N2a = N2*Naux;
  int nna = N*Naux;

  FP1* tmp = new FP1[N2a];

  for (int i=0;i<Naux;i++) 
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    tmp[j*nna+k*Naux+i] = C[i*N2+j*N+k];

  for (int i=0;i<N2a;i++)
    C[i] = tmp[i];

  delete [] tmp;

  return;
}

// if !EVL64 commented out by Soumi
//#if !EVL64
void transpose_C(int Naux, int N, FP2* C)
{
  int N2 = N*N;
  int N2a = N2*Naux;
  int nna = N*Naux;

  FP2* tmp = new FP2[N2a];

  for (int i=0;i<Naux;i++) 
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    tmp[j*nna+k*Naux+i] = C[i*N2+j*N+k];

  for (int i=0;i<N2a;i++)
    C[i] = tmp[i];

  delete [] tmp;

  return;
}
//#endif

void copy_symm(int natoms, int N, int Naux, vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, FP1* C, int type)  
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

// if !EVL64 commented out by Soumi
//#if !EVL64
void copy_symm(int natoms, int N, int Naux, vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, FP2* C, int type)
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
//#endif

int get_imax_n2i(int natoms, int N, vector<vector<FP2> >& basis, int* n2i)
{
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

#if RED_DOUBLE
void compute_all_3c_para(int ngpu, int natoms, int* atno, FP1* coords, 
                         vector<vector<FP2> > &basis, 
                         vector<vector<FP2> > &basis_aux, 
                         int nrad, int nang, FP2* ang_g0, FP2* ang_w0, 
                         FP2* C, int prl)
#else
void compute_all_3c_para(int ngpu, int natoms, int* atno, FP1* coords, 
                         vector<vector<FP2> > &basis, 
                         vector<vector<FP2> > &basis_aux, 
                         int nrad, int nang, FP2* ang_g0, FP2* ang_w0, 
                         FP1* C, int prl)
#endif
{
 #if !USE_ACC
  return compute_all_3c_v2(natoms,atno,coords,basis,basis_aux,nrad,nang,ang_g0,ang_w0,C,prl);
 #endif

  int nomp = ngpu;
 //#pragma omp parallel
  //nomp = omp_get_num_threads(); 

  if (prl>1) printf(" beginning compute_all_3c_para. nthreads: %2i \n",nomp);

  int N = basis.size();
  int N2 = N*N;
  int Naux = basis_aux.size();
  int N2a = N2*Naux;
  int gs = nrad*nang;
  int gs6 = 6*gs;

  int estart = find_center_of_grid(1,nrad)*nang;
 
  int* na2i = new int[natoms];
  int iNa = get_imax_n2i(natoms,Naux,basis_aux,na2i);

  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);

  //printf("  iN/a: %i %i \n",iN,iNa);

  FP1* grid1 = new FP1[gs6];
  FP1* wt1 = new FP1[gs];

  FP1* grid2 = new FP1[gs6];
  FP1* wt2 = new FP1[gs];

  FP1* grid3 = new FP1[gs6];
  FP1* wt3 = new FP1[gs];

  FP1** val1 = new FP1*[iNa];
  FP1** val2 = new FP1*[iN];
  FP1** val3 = new FP1*[iN];
  FP1** val4 = new FP1*[iNa];
  FP1** val5 = new FP1*[iN];
  FP1** val6 = new FP1*[iN];
  FP1** val7 = new FP1*[iNa];
  FP1** val8 = new FP1*[iN];
  FP1** val9 = new FP1*[iN];
  for (int i=0;i<iNa;i++) val1[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val2[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val3[i] = new FP1[gs];
  for (int i=0;i<iNa;i++) val4[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val5[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val6[i] = new FP1[gs];
  for (int i=0;i<iNa;i++) val7[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val8[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val9[i] = new FP1[gs];
  FP1* valt1 = new FP1[gs];
  FP1* valt2 = new FP1[gs];
  FP1* valt3 = new FP1[gs];

  FP1* grid1s = new FP1[gs6];
  FP1* grid2s = new FP1[gs6];
  FP1* grid3s = new FP1[gs6];

  FP1* grid1p = new FP1[gs6];
  FP1* grid2p = new FP1[gs6];
  FP1* grid3p = new FP1[gs6];

  FP1* wtt1 = new FP1[gs];
  FP1* wtt2 = new FP1[gs];

  FP1* ang_g = new FP1[3*nang];
  FP1* ang_w = new FP1[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

  //using namespace std::chrono;
  //high_resolution_clock::time_point t0 = high_resolution_clock::now();
 
#if USE_ACC
 #pragma omp parallel for schedule(static,1) num_threads(nomp)
#endif
  for (int n=0;n<nomp;n++)
  {
    int tid = omp_get_thread_num();
    #if USE_ACC
    acc_set_device_num(tid,acc_device_nvidia);
    #endif

    #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
    #pragma acc enter data copyin(n2i[0:natoms],na2i[0:natoms])

    #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
    #pragma acc enter data create(grid2[0:gs6],wt2[0:gs])
    #pragma acc enter data create(grid3[0:gs6],wt3[0:gs])
    #pragma acc enter data create(val1[0:iNa][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs]) 
    #pragma acc enter data create(val4[0:iNa][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs]) 
    #pragma acc enter data create(val7[0:iNa][0:gs],val8[0:iN][0:gs],val9[0:iN][0:gs]) 

    #pragma acc enter data create(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6])
    #pragma acc enter data create(grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6])

    #pragma acc enter data create(wtt1[0:gs],wtt2[0:gs],valt1[0:gs],valt2[0:gs],valt3[0:gs])
    if (tid>0)
    {
      #pragma acc enter data create(C[0:N2a])
    }
    acc_assign(N2a,C,0.);
  }
  #if USE_ACC
  acc_set_device_num(0,acc_device_nvidia);
  #endif

  //high_resolution_clock::time_point t1 = high_resolution_clock::now();

#if USE_ACC
 #pragma omp parallel for num_threads(nomp)
#endif
  for (int m=0;m<natoms;m++)
  {
    int tid = omp_get_thread_num();
    #if USE_ACC
    acc_set_device_num(tid,acc_device_nvidia);
    #endif
    //printf("  launch %i/%i \n",m,tid);

    int s1 = 0; if (m>0) s1 = na2i[m-1]; int s2 = na2i[m];
    int s3 = 0; if (m>0) s3 = n2i[m-1]; int s4 = n2i[m];

    FP1 Z1 = (FP1)atno[m];
    FP1 A1 = coords[3*m+0]; FP1 B1 = coords[3*m+1]; FP1 C1 = coords[3*m+2];
    
    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);

   #pragma acc parallel loop present(val1[0:iNa][0:gs])
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val1[ii1][j] = 1.f;
    }

   #pragma acc parallel loop present(val2[0:iN][0:gs])
    for (int ii2=0;ii2<s4-s3;ii2++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val2[ii2][j] = 1.f;
    }
   #pragma acc parallel loop present(val3[0:iN][0:gs],wt1[0:gs])
    for (int ii3=0;ii3<s4-s3;ii3++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val3[ii3][j] = wt1[j];
    }

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<FP2> basis1 = basis_aux[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

      eval_inr_r12(gs,grid1,val1[ii1],n1,l1,zeta1);
      eval_sh_3r(gs,grid1,val1[ii1],n1,l1,m1);
    } //loop i1

    for (int i2=s3;i2<s4;i2++)
    {
      int ii2 = i2-s3;

      vector<FP2> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

      eval_sh(ii2,gs,grid1,val2[ii2],n2,l2,m2,zeta2);
    } //loop i2

    for (int i3=s3;i3<s4;i3++)
    {
      int ii3 = i3-s3;

      vector<FP2> basis3 = basis[i3];
      int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];

      eval_sh(ii3,gs,grid1,val3[ii3],n3,l3,m3,zeta3);
    } //loop i3
    
    reduce_3c1b(s1,s2,s3,s4,gs,val1,val2,val3,N,Naux,iN,iNa,C);
    // printf("%d %d %d\n",m,m,m);

#if 1
   //two-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      // printf("%d %d %d\n",m,n,n);
     //s12 over atom m, aux function
     //s34 over atom m, regular function
     //s56 over atom n, regular function
      int s5 = 0; if (n>0) s5 = n2i[n-1]; int s6 = n2i[n];

      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2);
      recenter_grid(gs,grid2,A12,B12,C12);

      copy_grid(gs,grid1s,grid1); 
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid2,wt2,Z1,Z2,A12,B12,C12);
      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wt2);

      add_r1_to_grid(gs,grid2,0.,0.,0.);

     #pragma acc parallel loop collapse(2) present(val4[0:iNa][0:gs])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
        for (int j=0;j<gs;j++)
          val4[ii1][j] = 1.f;
      }

     #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs],val5[0:iN][0:gs])
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
        for (int j=0;j<gs;j++)
        {
          val2[ii2][j] = 1.f;
          val5[ii2][j] = 1.f;
        }
      }

     #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],wtt1[0:gs],wt2[0:gs])
      for (int ii3=0;ii3<s6-s5;ii3++)
      {
        for (int j=0;j<gs;j++)
        {
          val3[ii3][j] = wtt1[j];
          val6[ii3][j] = wt2[j];
        }
      }

     //i1 on atom m
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis_aux[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

        eval_inr_r12(gs,grid2,val4[ii1],n1,l1,zeta1);
        eval_sh_3r(gs,grid2,val4[ii1],n1,l1,m1);
      }

     //i2 on atom m
      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;

        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1,val2[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2,val5[ii2],n2,l2,m2,zeta2);
      }

     //i3 on atom n
      for (int i3=s5;i3<s6;i3++)
      {
        int ii3 = i3-s5;

        vector<FP2> basis3 = basis[i3];
        int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];
    
        eval_sh(ii3,gs,grid1s,val3[ii3],n3,l3,m3,zeta3);
        eval_sh(ii3,gs,grid2s,val6[ii3],n3,l3,m3,zeta3);
      }

      reduce_3c2b(s1,s2,s3,s4,s5,s6,gs,val1,val2,val3,val4,val5,val6,N,Naux,iN,iNa,C);

     #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs],val5[0:iN][0:gs])
      for (int ii2=0;ii2<s6-s5;ii2++)
      {
        for (int j=0;j<gs;j++)
        {
          val2[ii2][j] = 1.f;
          val5[ii2][j] = 1.f;
        }
      }

     #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],wtt1[0:gs],wt2[0:gs])
      for (int ii3=0;ii3<s6-s5;ii3++)
      {
        for (int j=0;j<gs;j++)
        {
          val3[ii3][j] = wtt1[j];
          val6[ii3][j] = wt2[j];
        }
      }

     //now do i2+i3 on atom n 
     //i2 on atom n
      for (int i2=s5;i2<s6;i2++)
      {
        int ii2 = i2-s5;

        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1s,val2[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2s,val5[ii2],n2,l2,m2,zeta2);
      }

     //i3 on atom n
      for (int i3=s5;i3<s6;i3++)
      {
        int ii3 = i3-s5;

        vector<FP2> basis3 = basis[i3];
        int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];
    
        eval_sh(ii3,gs,grid1s,val3[ii3],n3,l3,m3,zeta3);
        eval_sh(ii3,gs,grid2s,val6[ii3],n3,l3,m3,zeta3);
      }

      reduce_3c2b(s1,s2,s5,s6,s5,s6,gs,val1,val2,val3,val4,val5,val6,N,Naux,iN,iNa,C);
 
    } //loop n over second atom
#endif 

#if 1
   //three-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1);
      recenter_grid(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

     //testing condition
      for (int p=n+1;p<natoms;p++)
      if (p!=m)
      //for (int p=0;p<natoms;p++)
      //if (p!=m && p!=n)
      {
        // printf("%d %d %d\n",m,n,p);
        int s5 = 0; if (p>0) s5 = n2i[p-1]; int s6 = n2i[p];

        FP1 Z3 = (FP1)atno[p];
        FP1 A3 = coords[3*p+0]; FP1 B3 = coords[3*p+1]; FP1 C3 = coords[3*p+2];
        FP1 A13 = A3-A1; FP1 B13 = B3-B1; FP1 C13 = C3-C1;

        generate_central_grid_2(grid3,wt3,Z3,nrad,nang,ang_g,ang_w);
        copy_grid(gs,grid3p,grid3); //grid 3 centered on atom 3
        recenter_grid(gs,grid3,A13,B13,C13); //grid 3 centered on atom 1

        copy_grid(gs,grid3s,grid3);
        recenter_grid(gs,grid3s,-A12,-B12,-C12); //grid 3 centered on atom 2

        copy_grid(gs,grid2p,grid2); //grid 2 centered on atom 3
        recenter_grid(gs,grid2p,-A13,-B13,-C13);

        copy_grid(gs,grid1p,grid1);
        recenter_grid(gs,grid1p,-A13,-B13,-C13); //grid 1 centered on atom 3

        add_r1_to_grid_6z(gs,grid1s,grid2s,grid3s,grid1p,grid2p,grid3p);

      //need to keep all of these distances in order
        add_r3_to_grid(gs,grid1,A13,B13,C13);
        add_r123_to_grid(gs,grid2,A12,B12,C12,0.,0.,0.,A13,B13,C13);
        add_r123_to_grid(gs,grid3,A13,B13,C13,A12,B12,C12,0.,0.,0.);

        acc_copyf(gs,wtt1,wt1,wtt2,wt2);

        becke_weight_3c(gs,grid1,wtt1,grid2,wtt2,grid3,wt3,Z1,Z2,Z3,A12,B12,C12,A13,B13,C13);
        eliminate_small_wt_3(estart,gs,wtt1,wtt2,wt3);

        add_r1_to_grid(gs,grid2,0.,0.,0.);
        add_r1_to_grid(gs,grid3,0.,0.,0.);


       #pragma acc parallel loop collapse(2) present(val4[0:iNa][0:gs],val7[0:iNa][0:gs])
        for (int ii1=0;ii1<s2-s1;ii1++)
        {
          for (int j=0;j<gs;j++)
          {
            val4[ii1][j] = 1.f;
            val7[ii1][j] = 1.f;
          }
        }

       #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs],val5[0:iN][0:gs],val8[0:iN][0:gs])
        for (int ii2=0;ii2<s4-s3;ii2++)
        {
          for (int j=0;j<gs;j++)
          {
            val2[ii2][j] = 1.f;
            val5[ii2][j] = 1.f;
            val8[ii2][j] = 1.f;
          }
        }

       #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],val9[0:iN][0:gs],wtt1[0:gs],wtt2[0:gs],wt3[0:gs])
        for (int ii3=0;ii3<s6-s5;ii3++)
        {
          for (int j=0;j<gs;j++)
          {
            val3[ii3][j] = wtt1[j];
            val6[ii3][j] = wtt2[j];
            val9[ii3][j] = wt3[j];
          }
        }

 
        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<FP2> basis1 = basis_aux[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

          eval_inr_r12(gs,grid2,val4[ii1],n1,l1,zeta1);
          eval_sh_3r(gs,grid2,val4[ii1],n1,l1,m1);
          eval_inr_r12(gs,grid3,val7[ii1],n1,l1,zeta1);
          eval_sh_3r(gs,grid3,val7[ii1],n1,l1,m1);
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<FP2> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

          eval_sh(ii2,gs,grid3s,val8[ii2],n2,l2,m2,zeta2);
          eval_sh(ii2,gs,grid2s,val5[ii2],n2,l2,m2,zeta2);
          eval_sh(ii2,gs,grid1s,val2[ii2],n2,l2,m2,zeta2);
        }

        for (int i3=s5;i3<s6;i3++)
        {
          int ii3 = i3-s5;
          vector<FP2> basis3 = basis[i3];
          int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];
    
          eval_sh(ii3,gs,grid3p,val9[ii3],n3,l3,m3,zeta3);
          eval_sh(ii3,gs,grid2p,val6[ii3],n3,l3,m3,zeta3);
          eval_sh(ii3,gs,grid1p,val3[ii3],n3,l3,m3,zeta3);
        }

        reduce_3c3b(s1,s2,s3,s4,s5,s6,gs,val1,val2,val3,val4,val5,val6,val7,val8,val9,N,Naux,iN,iNa,C);

      } //loop p over third atom
    } //loop n over second atom
#endif
  } //loop m over natoms

  //high_resolution_clock::time_point t2 = high_resolution_clock::now();

 //collect integrals from all GPUs
  FP2* C_all = new FP2[N2a]();

  for (int n=0;n<nomp;n++)
  {
    #if USE_ACC
    acc_set_device_num(n,acc_device_nvidia);
    #endif

    #pragma acc update self(C[0:N2a])

   #if 0
    printf("\n C(%i): \n",n);
    for (int i=0;i<Naux;i++)
    {
      printf(" i: %i \n",i);
      for (int j=0;j<N;j++)
      for (int k=0;k<N;k++)
        printf("  %12.6f",C[i*N2+j*N+k]);
      printf("\n");
    }
   #endif

    for (int i=0;i<N2a;i++)
      C_all[i] += C[i];

    #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
    if (n>0)
    {
      #pragma acc exit data delete(C[0:N2a])
    }
  }
  #if USE_ACC
  acc_set_device_num(0,acc_device_nvidia);
  #endif

  for (int i=0;i<N2a;i++)
    C[i] = C_all[i];

  delete [] C_all;

 //apply symmetry and normalization
  copy_symm(natoms,N,Naux,basis,basis_aux,C,1);

  for (int i=0;i<Naux;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    C[i*N2+j*N+k] *= norm_sv(basis_aux[i][0],basis_aux[i][1],basis_aux[i][2],basis_aux[i][3])*basis[j][4]*basis[k][4];

  if (prl>2)
  {
    printf("\n C: \n");
    for (int i=0;i<Naux;i++)
    {
      printf(" i: %i \n",i);
      for (int j=0;j<N;j++)
      for (int k=0;k<N;k++)
        printf("  %12.6f",C[i*N2+j*N+k]);
      printf("\n");
    }
  }

  transpose_C(Naux,N,C);
  #pragma acc update device(C[0:N2a])


 //CPMZ check this
#if USE_ACC
 #pragma omp parallel for schedule(static,1) num_threads(nomp)
#endif
  for (int n=0;n<nomp;n++)
  {
    int tid = omp_get_thread_num();
    #if USE_ACC
    acc_set_device_num(tid,acc_device_nvidia);
    #endif

    #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
    #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs])
    #pragma acc exit data delete(grid3[0:gs6],wt3[0:gs])

    #pragma acc exit data delete(val1[0:iNa][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs]) 
    #pragma acc exit data delete(val4[0:iNa][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs]) 
    #pragma acc exit data delete(val7[0:iNa][0:gs],val8[0:iN][0:gs],val9[0:iN][0:gs]) 

    #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6]) 
    #pragma acc exit data delete(grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6]) 

    #pragma acc exit data delete(wtt1[0:gs],wtt2[0:gs],valt1[0:gs],valt2[0:gs],valt3[0:gs])
    #pragma acc exit data delete(n2i[0:natoms],na2i[0:natoms])
  }

 //dominated by step2
  //high_resolution_clock::time_point t3 = high_resolution_clock::now();
  //print_duration(t0,t1,"v3 step1");
  //print_duration(t1,t2,"v3 step2");
  //print_duration(t2,t3,"v3 step3");
  //auto_crash();

  delete [] n2i;
  delete [] na2i;

  delete [] ang_g;
  delete [] ang_w;

  delete [] grid1s;
  delete [] grid2s;
  delete [] grid3s;
  delete [] grid1p;
  delete [] grid2p;
  delete [] grid3p;

  delete [] wtt1;
  delete [] wtt2;
  delete [] grid1;
  delete [] grid2;
  delete [] grid3;
  delete [] wt1;
  delete [] wt2;
  delete [] wt3;
  for (int i=0;i<iNa;i++) delete [] val1[i];
  for (int i=0;i<iN;i++) delete [] val2[i];
  for (int i=0;i<iN;i++) delete [] val3[i];
  for (int i=0;i<iNa;i++) delete [] val4[i];
  for (int i=0;i<iN;i++) delete [] val5[i];
  for (int i=0;i<iN;i++) delete [] val6[i];
  for (int i=0;i<iNa;i++) delete [] val7[i];
  for (int i=0;i<iN;i++) delete [] val8[i];
  for (int i=0;i<iN;i++) delete [] val9[i];
  delete [] val1;
  delete [] val2;
  delete [] val3;
  delete [] val4;
  delete [] val5;
  delete [] val6;
  delete [] val7;
  delete [] val8;
  delete [] val9;
  delete [] valt1;
  delete [] valt2;
  delete [] valt3;

  //printf("\n ENDING EARLY \n");
  //exit(0);

  return;
}

void compute_all_3c_para2(int ngpu, int natoms, int* atno, FP1* coords, 
                         vector<vector<FP2> > &basis, 
                         vector<vector<FP2> > &basis_aux, 
                         int nrad, int nang, FP2* ang_g0, FP2* ang_w0, 
                         FP2* C, int prl)
{
 #if !USE_ACC
  return compute_all_3c_v2(natoms,atno,coords,basis,basis_aux,nrad,nang,ang_g0,ang_w0,C,prl);
 #endif

  int nomp = ngpu;
 //#pragma omp parallel
  //nomp = omp_get_num_threads(); 

  if (prl>1) printf(" beginning compute_all_3c_para. nthreads: %2i \n",nomp);

  int N = basis.size();
  int N2 = N*N;
  int Naux = basis_aux.size();
  int N2a = N2*Naux;
  int gs = nrad*nang;
  int gs6 = 6*gs;

  int estart = find_center_of_grid(1,nrad)*nang;
 
  int* na2i = new int[natoms];
  int iNa = get_imax_n2i(natoms,Naux,basis_aux,na2i);

  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);

  //printf("  iN/a: %i %i \n",iN,iNa);

  FP1* grid1 = new FP1[gs6];
  FP1* wt1 = new FP1[gs];

  FP1* grid2 = new FP1[gs6];
  FP1* wt2 = new FP1[gs];

  FP1* grid3 = new FP1[gs6];
  FP1* wt3 = new FP1[gs];

  FP1** val1 = new FP1*[iNa];
  FP1** val2 = new FP1*[iN];
  FP1** val3 = new FP1*[iN];
  FP1** val4 = new FP1*[iNa];
  FP1** val5 = new FP1*[iN];
  FP1** val6 = new FP1*[iN];
  FP1** val7 = new FP1*[iNa];
  FP1** val8 = new FP1*[iN];
  FP1** val9 = new FP1*[iN];
  for (int i=0;i<iNa;i++) val1[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val2[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val3[i] = new FP1[gs];
  for (int i=0;i<iNa;i++) val4[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val5[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val6[i] = new FP1[gs];
  for (int i=0;i<iNa;i++) val7[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val8[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val9[i] = new FP1[gs];
  FP1* valt1 = new FP1[gs];
  FP1* valt2 = new FP1[gs];
  FP1* valt3 = new FP1[gs];

  FP1* grid1s = new FP1[gs6];
  FP1* grid2s = new FP1[gs6];
  FP1* grid3s = new FP1[gs6];

  FP1* grid1p = new FP1[gs6];
  FP1* grid2p = new FP1[gs6];
  FP1* grid3p = new FP1[gs6];

  FP1* wtt1 = new FP1[gs];
  FP1* wtt2 = new FP1[gs];

  FP1* ang_g = new FP1[3*nang];
  FP1* ang_w = new FP1[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

  //using namespace std::chrono;
  //high_resolution_clock::time_point t0 = high_resolution_clock::now();
 
  #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = omp_get_thread_num();
    #if USE_ACC
    acc_set_device_num(tid,acc_device_nvidia);
    #endif

    #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
    #pragma acc enter data copyin(n2i[0:natoms],na2i[0:natoms])

    #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
    #pragma acc enter data create(grid2[0:gs6],wt2[0:gs])
    #pragma acc enter data create(grid3[0:gs6],wt3[0:gs])
    #pragma acc enter data create(val1[0:iNa][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs]) 
    #pragma acc enter data create(val4[0:iNa][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs]) 
    #pragma acc enter data create(val7[0:iNa][0:gs],val8[0:iN][0:gs],val9[0:iN][0:gs]) 

    #pragma acc enter data create(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6])
    #pragma acc enter data create(grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6])

    #pragma acc enter data create(wtt1[0:gs],wtt2[0:gs],valt1[0:gs],valt2[0:gs],valt3[0:gs])
    if (tid>0)
    {
      #pragma acc enter data create(C[0:N2a])
    }
    acc_assign(N2a,C,0.);
  }
  #if USE_ACC
  acc_set_device_num(0,acc_device_nvidia);
  #endif

  //high_resolution_clock::time_point t1 = high_resolution_clock::now();
  //Setting up manager-worker scheme
  int work_size = natoms;
  if (natoms > 1) {
    work_size += natoms * (natoms - 1);
  }
  if (natoms > 2) {
    work_size += natoms * (natoms - 1) * (natoms - 2) / 2;
  }
  int * work = new int[3*work_size]();
  int count = 0;
  for (int i = 0; i < natoms; i++) {
    work[3*i+0] = i;
    work[3*i+1] = i;
    work[3*i+2] = i;
    // printf("%d %d %d\n",i,i,i);
    count++;
  }
  for (int i = 0; i < natoms; i++) {
    for (int j = 0; j < natoms; j++) {
      if (i != j) {
        work[3*(count)+0] = i;
        work[3*(count)+1] = j;
        work[3*(count)+2] = j;
        // printf("%d %d %d\n",i,j,j);
        count++;
      }
    }
  }
  for (int i = 0; i < natoms; i++) {
    for (int j = 0; j < natoms; j++) {
      if (i != j) {
        for (int k = j+1; k < natoms; k++) {
          if (k != i) {
            work[3*(count)+0] = i;
            work[3*(count)+1] = j;
            work[3*(count)+2] = k;
            // printf("%d %d %d\n",i,j,k);
            count++;
          }
        }
      }
    }
  }
  printf("\n");

  // Starting integrals
  int tasks = 0;
  #pragma omp parallel num_threads(nomp)
  {
    int tid = omp_get_thread_num();
    int num_threads = omp_get_num_threads();
    int task = tid;

    #if USE_ACC
    acc_set_device_num(tid,acc_device_nvidia);
    #endif

    #pragma omp single
    tasks = num_threads;

    while (tasks <= count) {
      int m = work[3*task+0];
      int n = work[3*task+1];
      int p = work[3*task+2];

      // Compute grid for m first
      int s1 = 0; if (m>0) s1 = na2i[m-1]; int s2 = na2i[m];
      int s3 = 0; if (m>0) s3 = n2i[m-1]; int s4 = n2i[m];

      FP1 Z1 = (FP1)atno[m];
      FP1 A1 = coords[3*m+0]; 
      FP1 B1 = coords[3*m+1]; 
      FP1 C1 = coords[3*m+2];

      generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);
      #pragma acc parallel loop present(val1[0:iNa][0:gs])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
      #pragma acc loop
        for (int j=0;j<gs;j++)
          val1[ii1][j] = 1.f;
      }

      #pragma acc parallel loop present(val2[0:iN][0:gs])
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
      #pragma acc loop
        for (int j=0;j<gs;j++)
          val2[ii2][j] = 1.f;
      }
      #pragma acc parallel loop present(val3[0:iN][0:gs],wt1[0:gs])
      for (int ii3=0;ii3<s4-s3;ii3++)
      {
        #pragma acc loop
        for (int j=0;j<gs;j++)
          val3[ii3][j] = wt1[j];
      }
      // Done with m grid
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis_aux[i1];
        int n1 = basis1[0]; 
        int l1 = basis1[1]; 
        int m1 = basis1[2]; 
        FP1 zeta1 = basis1[3];

        eval_inr_r12(gs,grid1,val1[ii1],n1,l1,zeta1);
        eval_sh_3r(gs,grid1,val1[ii1],n1,l1,m1);
      } //loop i1

      // single atom, note by transitive property n == p... >.<
      if (m == n && m == p) {
        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;

          vector<FP2> basis2 = basis[i2];
          int n2 = basis2[0]; 
          int l2 = basis2[1]; 
          int m2 = basis2[2]; 
          FP1 zeta2 = basis2[3];

          eval_sh(ii2,gs,grid1,val2[ii2],n2,l2,m2,zeta2);
        } //loop i2

        for (int i3=s3;i3<s4;i3++)
        {
          int ii3 = i3-s3;

          vector<FP2> basis3 = basis[i3];
          int n3 = basis3[0]; 
          int l3 = basis3[1]; 
          int m3 = basis3[2]; 
          FP1 zeta3 = basis3[3];

          eval_sh(ii3,gs,grid1,val3[ii3],n3,l3,m3,zeta3);
        } //loop i3
        
        reduce_3c1b(s1,s2,s3,s4,gs,val1,val2,val3,N,Naux,iN,iNa,C);
      } // end single atom
      // two atoms
      else if (m != n && n == p) {
        int s5 = 0; if (n>0) s5 = n2i[n-1]; int s6 = n2i[n];

        FP1 Z2 = (FP1)atno[n];
        FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
        FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;
        add_r2_to_grid(gs,grid1,A12,B12,C12);

        generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
        copy_grid(gs,grid2s,grid2);
        recenter_grid(gs,grid2,A12,B12,C12);

        copy_grid(gs,grid1s,grid1); 
        recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

        acc_copyf(gs,wtt1,wt1);
        becke_weight_2c(gs,grid1,wtt1,grid2,wt2,Z1,Z2,A12,B12,C12);
        eliminate_small_wt(estart,gs,wtt1);
        eliminate_small_wt(estart,gs,wt2);

        add_r1_to_grid(gs,grid2,0.,0.,0.);

        #pragma acc parallel loop collapse(2) present(val4[0:iNa][0:gs])
        for (int ii1=0;ii1<s2-s1;ii1++)
        {
          for (int j=0;j<gs;j++)
            val4[ii1][j] = 1.f;
        }

        #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs],val5[0:iN][0:gs])
        for (int ii2=0;ii2<s4-s3;ii2++)
        {
          for (int j=0;j<gs;j++)
          {
            val2[ii2][j] = 1.f;
            val5[ii2][j] = 1.f;
          }
        }

        #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],wtt1[0:gs],wt2[0:gs])
        for (int ii3=0;ii3<s6-s5;ii3++)
        {
          for (int j=0;j<gs;j++)
          {
            val3[ii3][j] = wtt1[j];
            val6[ii3][j] = wt2[j];
          }
        }

        //i1 on atom m
        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<FP2> basis1 = basis_aux[i1];
          int n1 = basis1[0]; 
          int l1 = basis1[1]; 
          int m1 = basis1[2]; 
          FP1 zeta1 = basis1[3];

          eval_inr_r12(gs,grid2,val4[ii1],n1,l1,zeta1);
          eval_sh_3r(gs,grid2,val4[ii1],n1,l1,m1);
        }

        //i2 on atom m
        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;

          vector<FP2> basis2 = basis[i2];
          int n2 = basis2[0]; 
          int l2 = basis2[1]; 
          int m2 = basis2[2]; 
          FP1 zeta2 = basis2[3];

          eval_sh(ii2,gs,grid1,val2[ii2],n2,l2,m2,zeta2);
          eval_sh(ii2,gs,grid2,val5[ii2],n2,l2,m2,zeta2);
        }

        //i3 on atom n
        for (int i3=s5;i3<s6;i3++)
        {
          int ii3 = i3-s5;

          vector<FP2> basis3 = basis[i3];
          int n3 = basis3[0]; 
          int l3 = basis3[1]; 
          int m3 = basis3[2]; 
          FP1 zeta3 = basis3[3];
      
          eval_sh(ii3,gs,grid1s,val3[ii3],n3,l3,m3,zeta3);
          eval_sh(ii3,gs,grid2s,val6[ii3],n3,l3,m3,zeta3);
        }

        reduce_3c2b(s1,s2,s3,s4,s5,s6,gs,val1,val2,val3,val4,val5,val6,N,Naux,iN,iNa,C);

        #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs],val5[0:iN][0:gs])
        for (int ii2=0;ii2<s6-s5;ii2++)
        {
          for (int j=0;j<gs;j++)
          {
            val2[ii2][j] = 1.f;
            val5[ii2][j] = 1.f;
          }
        }

        #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],wtt1[0:gs],wt2[0:gs])
        for (int ii3=0;ii3<s6-s5;ii3++)
        {
          for (int j=0;j<gs;j++)
          {
            val3[ii3][j] = wtt1[j];
            val6[ii3][j] = wt2[j];
          }
        }

        // now do i2+i3 on atom n 
        // i2 on atom n
        for (int i2=s5;i2<s6;i2++)
        {
          int ii2 = i2-s5;

          vector<FP2> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

          eval_sh(ii2,gs,grid1s,val2[ii2],n2,l2,m2,zeta2);
          eval_sh(ii2,gs,grid2s,val5[ii2],n2,l2,m2,zeta2);
        }

        //i3 on atom n
        for (int i3=s5;i3<s6;i3++)
        {
          int ii3 = i3-s5;

          vector<FP2> basis3 = basis[i3];
          int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];
      
          eval_sh(ii3,gs,grid1s,val3[ii3],n3,l3,m3,zeta3);
          eval_sh(ii3,gs,grid2s,val6[ii3],n3,l3,m3,zeta3);
        }

        reduce_3c2b(s1,s2,s5,s6,s5,s6,gs,val1,val2,val3,val4,val5,val6,N,Naux,iN,iNa,C);
      } // end 2 atom
      // 3 atoms
      else if (m != n && m != p) {
        // printf("%d %d %d\n",m,n,p);
        int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

        FP1 Z2 = (FP1)atno[n];
        FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
        FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

        //grid1 at 0,0,0 now has r1 at 3, r2 at 4
        add_r2_to_grid(gs,grid1,A12,B12,C12);

        generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
        copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
        recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

        copy_grid(gs,grid1s,grid1);
        recenter_grid(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

        int s5 = 0; if (p>0) s5 = n2i[p-1]; int s6 = n2i[p];

        FP1 Z3 = (FP1)atno[p];
        FP1 A3 = coords[3*p+0]; FP1 B3 = coords[3*p+1]; FP1 C3 = coords[3*p+2];
        FP1 A13 = A3-A1; FP1 B13 = B3-B1; FP1 C13 = C3-C1;

        generate_central_grid_2(grid3,wt3,Z3,nrad,nang,ang_g,ang_w);
        copy_grid(gs,grid3p,grid3); //grid 3 centered on atom 3
        recenter_grid(gs,grid3,A13,B13,C13); //grid 3 centered on atom 1

        copy_grid(gs,grid3s,grid3);
        recenter_grid(gs,grid3s,-A12,-B12,-C12); //grid 3 centered on atom 2

        copy_grid(gs,grid2p,grid2); //grid 2 centered on atom 3
        recenter_grid(gs,grid2p,-A13,-B13,-C13);

        copy_grid(gs,grid1p,grid1);
        recenter_grid(gs,grid1p,-A13,-B13,-C13); //grid 1 centered on atom 3

        add_r1_to_grid_6z(gs,grid1s,grid2s,grid3s,grid1p,grid2p,grid3p);

        //need to keep all of these distances in order
        add_r3_to_grid(gs,grid1,A13,B13,C13);
        add_r123_to_grid(gs,grid2,A12,B12,C12,0.,0.,0.,A13,B13,C13);
        add_r123_to_grid(gs,grid3,A13,B13,C13,A12,B12,C12,0.,0.,0.);

        acc_copyf(gs,wtt1,wt1,wtt2,wt2);

        becke_weight_3c(gs,grid1,wtt1,grid2,wtt2,grid3,wt3,Z1,Z2,Z3,A12,B12,C12,A13,B13,C13);
        eliminate_small_wt_3(estart,gs,wtt1,wtt2,wt3);

        add_r1_to_grid(gs,grid2,0.,0.,0.);
        add_r1_to_grid(gs,grid3,0.,0.,0.);

        #pragma acc parallel loop collapse(2) present(val4[0:iNa][0:gs],val7[0:iNa][0:gs])
        for (int ii1=0;ii1<s2-s1;ii1++)
        {
          for (int j=0;j<gs;j++)
          {
            val4[ii1][j] = 1.f;
            val7[ii1][j] = 1.f;
          }
        }

        #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs],val5[0:iN][0:gs],val8[0:iN][0:gs])
        for (int ii2=0;ii2<s4-s3;ii2++)
        {
          for (int j=0;j<gs;j++)
          {
            val2[ii2][j] = 1.f;
            val5[ii2][j] = 1.f;
            val8[ii2][j] = 1.f;
          }
        }

        #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],val9[0:iN][0:gs],wtt1[0:gs],wtt2[0:gs],wt3[0:gs])
        for (int ii3=0;ii3<s6-s5;ii3++)
        {
          for (int j=0;j<gs;j++)
          {
            val3[ii3][j] = wtt1[j];
            val6[ii3][j] = wtt2[j];
            val9[ii3][j] = wt3[j];
          }
        }

 
        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<FP2> basis1 = basis_aux[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

          eval_inr_r12(gs,grid2,val4[ii1],n1,l1,zeta1);
          eval_sh_3r(gs,grid2,val4[ii1],n1,l1,m1);
          eval_inr_r12(gs,grid3,val7[ii1],n1,l1,zeta1);
          eval_sh_3r(gs,grid3,val7[ii1],n1,l1,m1);
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<FP2> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

          eval_sh(ii2,gs,grid3s,val8[ii2],n2,l2,m2,zeta2);
          eval_sh(ii2,gs,grid2s,val5[ii2],n2,l2,m2,zeta2);
          eval_sh(ii2,gs,grid1s,val2[ii2],n2,l2,m2,zeta2);
        }

        for (int i3=s5;i3<s6;i3++)
        {
          int ii3 = i3-s5;
          vector<FP2> basis3 = basis[i3];
          int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];
    
          eval_sh(ii3,gs,grid3p,val9[ii3],n3,l3,m3,zeta3);
          eval_sh(ii3,gs,grid2p,val6[ii3],n3,l3,m3,zeta3);
          eval_sh(ii3,gs,grid1p,val3[ii3],n3,l3,m3,zeta3);
        }

        reduce_3c3b(s1,s2,s3,s4,s5,s6,gs,val1,val2,val3,val4,val5,val6,val7,val8,val9,N,Naux,iN,iNa,C);
      } // end 3 atoms

      #pragma omp atomic capture
      task = tasks++;
    } // while loop for manager worker
  } // omp parallel region
  delete [] work;

  //high_resolution_clock::time_point t2 = high_resolution_clock::now();

 //collect integrals from all GPUs
  FP2* C_all = new FP2[N2a]();

  for (int n=0;n<nomp;n++)
  {
    #if USE_ACC
    acc_set_device_num(n,acc_device_nvidia);
    #endif

    #pragma acc update self(C[0:N2a])

    for (int i=0;i<N2a;i++)
      C_all[i] += C[i];

    #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
    if (n>0)
    {
      #pragma acc exit data delete(C[0:N2a])
    }
  }
  #if USE_ACC
  acc_set_device_num(0,acc_device_nvidia);
  #endif

  for (int i=0;i<N2a;i++)
    C[i] = C_all[i];

  delete [] C_all;

  //apply symmetry and normalization
  copy_symm(natoms,N,Naux,basis,basis_aux,C,1);

  for (int i=0;i<Naux;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    C[i*N2+j*N+k] *= norm_sv(basis_aux[i][0],basis_aux[i][1],basis_aux[i][2],basis_aux[i][3])*basis[j][4]*basis[k][4];

  if (prl>2)
  {
    printf("\n C: \n");
    for (int i=0;i<Naux;i++)
    {
      printf(" i: %i \n",i);
      for (int j=0;j<N;j++)
      for (int k=0;k<N;k++)
        printf("  %12.6f",C[i*N2+j*N+k]);
      printf("\n");
    }
  }

  transpose_C(Naux,N,C);
  #pragma acc update device(C[0:N2a])

  #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = omp_get_thread_num();
    #if USE_ACC
    acc_set_device_num(tid,acc_device_nvidia);
    #endif

    #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
    #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs])
    #pragma acc exit data delete(grid3[0:gs6],wt3[0:gs])

    #pragma acc exit data delete(val1[0:iNa][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs]) 
    #pragma acc exit data delete(val4[0:iNa][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs]) 
    #pragma acc exit data delete(val7[0:iNa][0:gs],val8[0:iN][0:gs],val9[0:iN][0:gs]) 

    #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6]) 
    #pragma acc exit data delete(grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6]) 

    #pragma acc exit data delete(wtt1[0:gs],wtt2[0:gs],valt1[0:gs],valt2[0:gs],valt3[0:gs])
    #pragma acc exit data delete(n2i[0:natoms],na2i[0:natoms])
  }

  delete [] n2i;
  delete [] na2i;

  delete [] ang_g;
  delete [] ang_w;

  delete [] grid1s;
  delete [] grid2s;
  delete [] grid3s;
  delete [] grid1p;
  delete [] grid2p;
  delete [] grid3p;

  delete [] wtt1;
  delete [] wtt2;
  delete [] grid1;
  delete [] grid2;
  delete [] grid3;
  delete [] wt1;
  delete [] wt2;
  delete [] wt3;
  for (int i=0;i<iNa;i++) delete [] val1[i];
  for (int i=0;i<iN;i++) delete [] val2[i];
  for (int i=0;i<iN;i++) delete [] val3[i];
  for (int i=0;i<iNa;i++) delete [] val4[i];
  for (int i=0;i<iN;i++) delete [] val5[i];
  for (int i=0;i<iN;i++) delete [] val6[i];
  for (int i=0;i<iNa;i++) delete [] val7[i];
  for (int i=0;i<iN;i++) delete [] val8[i];
  for (int i=0;i<iN;i++) delete [] val9[i];
  delete [] val1;
  delete [] val2;
  delete [] val3;
  delete [] val4;
  delete [] val5;
  delete [] val6;
  delete [] val7;
  delete [] val8;
  delete [] val9;
  delete [] valt1;
  delete [] valt2;
  delete [] valt3;

  return;
}

void compute_all_3c_para3(int ngpu, int natoms, int* atno, FP1* coords, 
                         vector<vector<FP2> > &basis, 
                         vector<vector<FP2> > &basis_aux, 
                         int nrad, int nang, FP2* ang_g0, FP2* ang_w0, 
                         FP2* C, int prl)
{
 #if !USE_ACC
  return compute_all_3c_v2(natoms,atno,coords,basis,basis_aux,nrad,nang,ang_g0,ang_w0,C,prl);
 #endif

  int nomp = ngpu;
  // printf(" beginning compute_all_3c_para3. nthreads: %2i \n",nomp);

  int N = basis.size();
  int N2 = N*N;
  int Naux = basis_aux.size();
  int N2a = N2*Naux;
  int gs = nrad*nang;
  int gs6 = 6*gs;

  int estart = find_center_of_grid(1,nrad)*nang;

  FP1* grid1 = new FP1[gs6];
  FP1* wt1 = new FP1[gs];

  FP1* grid2 = new FP1[gs6];
  FP1* wt2 = new FP1[gs];

  FP1* grid3 = new FP1[gs6];
  FP1* wt3 = new FP1[gs];

  int* na2i = new int[natoms];
  int iNa = get_imax_n2i(natoms,Naux,basis_aux,na2i);

  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);

  FP1 * val1 = new FP1[gs];
  FP1 * val2 = new FP1[gs];
  FP1 * val3 = new FP1[gs];
  FP1 * val4 = new FP1[gs];
  FP1 * val5 = new FP1[gs];
  FP1 * val6 = new FP1[gs];
  FP1 * val7 = new FP1[gs];
  FP1 * val8 = new FP1[gs];
  FP1 * val9 = new FP1[gs];

  FP1* grid1s = new FP1[gs6];
  FP1* grid2s = new FP1[gs6];
  FP1* grid3s = new FP1[gs6];

  FP1* grid1p = new FP1[gs6];
  FP1* grid2p = new FP1[gs6];
  FP1* grid3p = new FP1[gs6];

  FP1* wtt1 = new FP1[gs];
  FP1* wtt2 = new FP1[gs];

  FP1* ang_g = new FP1[3*nang];
  FP1* ang_w = new FP1[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

  #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = omp_get_thread_num();
    #if USE_ACC
    acc_set_device_num(tid,acc_device_nvidia);
    #endif

    #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
    #pragma acc enter data copyin(n2i[0:natoms],na2i[0:natoms])

    #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
    #pragma acc enter data create(grid2[0:gs6],wt2[0:gs])
    #pragma acc enter data create(grid3[0:gs6],wt3[0:gs])
    #pragma acc enter data create(val1[0:gs],val2[0:gs],val3[0:gs]) 
    #pragma acc enter data create(val4[0:gs],val5[0:gs],val6[0:gs]) 
    #pragma acc enter data create(val7[0:gs],val8[0:gs],val9[0:gs]) 

    #pragma acc enter data create(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6])
    #pragma acc enter data create(grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6])

    #pragma acc enter data create(wtt1[0:gs],wtt2[0:gs])
    if (tid>0)
    {
      #pragma acc enter data create(C[0:N2a])
    }
    acc_assign(N2a,C,0.);
  }
  #if USE_ACC
  acc_set_device_num(0,acc_device_nvidia);
  #endif

  int tasks = 0;
  int count = N2a;
  #pragma omp parallel num_threads(nomp)
  {
    int tid = omp_get_thread_num();
    int num_threads = omp_get_num_threads();
    int task = tid;

    #if USE_ACC
    acc_set_device_num(tid,acc_device_nvidia);
    #endif
    
    #pragma omp single
    tasks = num_threads;

    while (task < count) {
      int mb = task / N2;
      int nb = (task - mb*N2)/N;
      int pb = task - mb * N2 - nb * N;
      
      int m = basis_aux[mb][9];
      int n = basis[nb][9];
      int p = basis[pb][9];

      if (1) {//(pb <= nb && p <= n) {
        FP1 Z1 = (FP1)atno[m];
        FP1 A1 = coords[3*m+0];
        FP1 B1 = coords[3*m+1];
        FP1 C1 = coords[3*m+2];

        FP1 Z2 = (FP1)atno[n];
        FP1 A2 = coords[3*n+0];
        FP1 B2 = coords[3*n+1];
        FP1 C2 = coords[3*n+2];

        FP1 Z3 = (FP1)atno[p];
        FP1 A3 = coords[3*p+0];
        FP1 B3 = coords[3*p+1];
        FP1 C3 = coords[3*p+2];

        generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);
        #pragma acc parallel loop present(val1[0:gs])
        for (int j=0;j<gs;j++)
          val1[j] = 1.;

        #pragma acc parallel loop present(val2[0:gs])
        for (int j=0;j<gs;j++)
          val2[j] = 1.;
        
        #pragma acc parallel loop present(val3[0:gs],wt1[0:gs])
        for (int j=0;j<gs;j++)
          val3[j] = wt1[j];
        // Done with m grid

        vector<FP2> basis1 = basis_aux[mb];
        int n1 = basis1[0]; 
        int l1 = basis1[1]; 
        int m1 = basis1[2]; 
        FP1 zeta1 = basis1[3];
        vector<FP2> basis2 = basis[nb];
        int n2 = basis2[0]; 
        int l2 = basis2[1]; 
        int m2 = basis2[2]; 
        FP1 zeta2 = basis2[3];
        vector<FP2> basis3 = basis[pb];
        int n3 = basis3[0]; 
        int l3 = basis3[1]; 
        int m3 = basis3[2]; 
        FP1 zeta3 = basis3[3];

        eval_inr_r12(gs,grid1,val1,n1,l1,zeta1);
        eval_sh_3r(gs,grid1,val1,n1,l1,m1);

        // single atom, note by transitive property n == p... >.<
        if (m == n && m == p) {  
          eval_sh(0,gs,grid1,val2,n2,l2,m2,zeta2);
          eval_sh(0,gs,grid1,val3,n3,l3,m3,zeta3);

          #pragma acc kernels present(val1[0:gs],val2[0:gs],val3[0:gs],C[0:N2a])
          {
            FP2 val = 0.;
            #pragma acc loop reduction(+:val) 
            for (int j = 0; j < gs; j++) {
              val += val1[j] * val2[j] * val3[j];
            }
            C[mb * N2 + nb * N + pb] = val;
          }
        } // end single atom
        // two atoms -- due to previous if, m == n == p is excluded
        else if (m != n && n == p) {
          FP1 A12, B12, C12;
          A12 = A2 - A1;
          B12 = B2 - B1;
          C12 = C2 - C1;
        
          add_r2_to_grid(gs,grid1,A12,B12,C12);

          generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
          copy_grid(gs,grid2s,grid2);
          recenter_grid(gs,grid2,A12,B12,C12);

          copy_grid(gs,grid1s,grid1); 
          recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

          acc_copyf(gs,wtt1,wt1);
          becke_weight_2c(gs,grid1,wtt1,grid2,wt2,Z1,Z2,A12,B12,C12);
          eliminate_small_wt(estart,gs,wtt1);
          eliminate_small_wt(estart,gs,wt2);

          add_r1_to_grid(gs,grid2,0.,0.,0.);

          #pragma acc parallel loop present(val4[0:gs])
          for (int j = 0; j < gs; j++) {
            val4[j] = 1;
          }
          #pragma acc parallel loop present(val2[0:gs],val5[0:gs])
          for (int j = 0; j < gs; j++) {
            val2[j] = 1;
            val5[j] = 1;
          }
          #pragma acc parallel loop present(val3[0:gs],val6[0:gs],wtt1[0:gs],wt2[0:gs])
          for (int j = 0; j < gs; j++) {
            val3[j] = wtt1[j];
            val6[j] = wt2[j];
          }

          // i1 on atom m
          eval_inr_r12(gs,grid2,val4,n1,l1,zeta1);
          eval_sh_3r(gs,grid2,val4,n1,l1,m1);

          eval_sh(0,gs,grid1s,val2,n2,l2,m2,zeta2);
          eval_sh(0,gs,grid2s,val5,n2,l2,m2,zeta2);
          eval_sh(0,gs,grid1s,val3,n3,l3,m3,zeta3);
          eval_sh(0,gs,grid2s,val6,n3,l3,m3,zeta3);
          #pragma acc kernels present(\
            val1[0:gs],val2[0:gs],val3[0:gs],val4[0:gs],val5[0:gs],val6[0:gs])
          {
            FP2 val = 0.;
            #pragma acc loop reduction(+:val)
            for (int j = 0; j < gs; j++) {
              val += val1[j] * val2[j] * val3[j];
            }

            #pragma acc loop reduction(+:val)
            for (int j = 0; j < gs; j++) {
              val += val4[j] * val5[j] * val6[j];
            }

            C[mb * N2 + nb * N + pb] = val;
          }
        }
        else if (m == n || m == p) {
          FP1 A12, B12, C12;
          if (m == n) {
            A12 = A3 - A1;
            B12 = B3 - B1;
            C12 = C3 - C1;
          }
          else if (m == p) {
            A12 = A2 - A1;
            B12 = B2 - B1;
            C12 = C2 - C1;
          }
          add_r2_to_grid(gs,grid1,A12,B12,C12);

          generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
          copy_grid(gs,grid2s,grid2);
          recenter_grid(gs,grid2,A12,B12,C12);

          copy_grid(gs,grid1s,grid1); 
          recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

          acc_copyf(gs,wtt1,wt1);
          becke_weight_2c(gs,grid1,wtt1,grid2,wt2,Z1,Z2,A12,B12,C12);
          eliminate_small_wt(estart,gs,wtt1);
          eliminate_small_wt(estart,gs,wt2);

          add_r1_to_grid(gs,grid2,0.,0.,0.);

          #pragma acc parallel loop present(val4[0:gs])
          for (int j = 0; j < gs; j++) {
            val4[j] = 1;
          }
          #pragma acc parallel loop present(val2[0:gs],val5[0:gs])
          for (int j = 0; j < gs; j++) {
            val2[j] = 1;
            val5[j] = 1;
          }
          #pragma acc parallel loop present(val3[0:gs],val6[0:gs],wtt1[0:gs],wt2[0:gs])
          for (int j = 0; j < gs; j++) {
            val3[j] = wtt1[j];
            val6[j] = wt2[j];
          }

          // i1
          eval_inr_r12(gs,grid2,val4,n1,l1,zeta1);
          eval_sh_3r(gs,grid2,val4,n1,l1,m1);

          if (m == n) {
            eval_sh(0,gs,grid1,val2,n2,l2,m2,zeta2);
            eval_sh(0,gs,grid2,val5,n2,l2,m2,zeta2);
            eval_sh(0,gs,grid1s,val3,n3,l3,m3,zeta3);
            eval_sh(0,gs,grid2s,val6,n3,l3,m3,zeta3);
          }
          else if (m == p) {
            eval_sh(0,gs,grid1s,val2,n2,l2,m2,zeta2);
            eval_sh(0,gs,grid2s,val5,n2,l2,m2,zeta2);
            eval_sh(0,gs,grid1,val3,n3,l3,m3,zeta3);
            eval_sh(0,gs,grid2,val6,n3,l3,m3,zeta3);
          }
          #pragma acc kernels present(\
            val1[0:gs],val2[0:gs],val3[0:gs],val4[0:gs],val5[0:gs],val6[0:gs])
          {
            FP2 val = 0.;
            #pragma acc loop reduction(+:val)
            for (int j = 0; j < gs; j++) {
              val += val1[j] * val2[j] * val3[j];
            }

            #pragma acc loop reduction(+:val)
            for (int j = 0; j < gs; j++) {
              val += val4[j] * val5[j] * val6[j];
            }

            C[mb * N2 + nb * N + pb] = val;
          }
        }
        // 3 atoms
        else if (m != n && m != p && n != p) {
          FP1 A12 = A2 - A1;
          FP1 B12 = B2 - B1;
          FP1 C12 = C2 - C1;
          FP1 A13 = A3 - A1;
          FP1 B13 = B3 - B1;
          FP1 C13 = C3 - C1;

          add_r2_to_grid(gs,grid1,A12,B12,C12);

          generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
          copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
          recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

          copy_grid(gs,grid1s,grid1);
          recenter_grid(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

          generate_central_grid_2(grid3,wt3,Z3,nrad,nang,ang_g,ang_w);
          copy_grid(gs,grid3p,grid3); //grid 3 centered on atom 3
          recenter_grid(gs,grid3,A13,B13,C13); //grid 3 centered on atom 1

          copy_grid(gs,grid3s,grid3);
          recenter_grid(gs,grid3s,-A12,-B12,-C12); //grid 3 centered on atom 2

          copy_grid(gs,grid2p,grid2); //grid 2 centered on atom 3
          recenter_grid(gs,grid2p,-A13,-B13,-C13);

          copy_grid(gs,grid1p,grid1);
          recenter_grid(gs,grid1p,-A13,-B13,-C13); //grid 1 centered on atom 3

          add_r1_to_grid_6z(gs,grid1s,grid2s,grid3s,grid1p,grid2p,grid3p);

          //need to keep all of these distances in order
          add_r3_to_grid(gs,grid1,A13,B13,C13);
          add_r123_to_grid(gs,grid2,A12,B12,C12,0.,0.,0.,A13,B13,C13);
          add_r123_to_grid(gs,grid3,A13,B13,C13,A12,B12,C12,0.,0.,0.);

          acc_copyf(gs,wtt1,wt1,wtt2,wt2);

          becke_weight_3c(gs,grid1,wtt1,grid2,wtt2,grid3,wt3,Z1,Z2,Z3,A12,B12,C12,A13,B13,C13);
          eliminate_small_wt_3(estart,gs,wtt1,wtt2,wt3);

          add_r1_to_grid(gs,grid2,0.,0.,0.);
          add_r1_to_grid(gs,grid3,0.,0.,0.);

          #pragma acc parallel loop present(val4[0:gs],val7[0:gs])
          for (int j=0;j<gs;j++)
          {
            val4[j] = 1.;
            val7[j] = 1.;
          }

          #pragma acc parallel loop present(val2[0:gs],val5[0:gs],val8[0:gs])
          for (int j=0;j<gs;j++)
          {
            val2[j] = 1.;
            val5[j] = 1.;
            val8[j] = 1.;
          }

          #pragma acc parallel loop present(\
            val3[0:gs],val6[0:gs],val9[0:gs],wtt1[0:gs],wtt2[0:gs],wt3[0:gs])
          for (int j=0;j<gs;j++)
          {
            val3[j] = wtt1[j];
            val6[j] = wtt2[j];
            val9[j] = wt3[j];
          }

          eval_inr_r12(gs,grid2,val4,n1,l1,zeta1);
          eval_sh_3r(gs,grid2,val4,n1,l1,m1);
          eval_inr_r12(gs,grid3,val7,n1,l1,zeta1);
          eval_sh_3r(gs,grid3,val7,n1,l1,m1);

          eval_sh(0,gs,grid3s,val8,n2,l2,m2,zeta2);
          eval_sh(0,gs,grid2s,val5,n2,l2,m2,zeta2);
          eval_sh(0,gs,grid1s,val2,n2,l2,m2,zeta2);
          
          eval_sh(0,gs,grid3p,val9,n3,l3,m3,zeta3);
          eval_sh(0,gs,grid2p,val6,n3,l3,m3,zeta3);
          eval_sh(0,gs,grid1p,val3,n3,l3,m3,zeta3);

          #pragma acc kernels present(\
            val1[0:gs],val2[0:gs],val3[0:gs],\
            val4[0:gs],val5[0:gs],val6[0:gs],\
            val7[0:gs],val8[0:gs],val9[0:gs],\
            C[0:N2a])
          {
            FP1 val = 0;
            #pragma acc loop reduction(+:val)
            for (int j = 0; j < gs; j++) {
              val += val1[j] * val2[j] * val3[j];
            }
            #pragma acc loop reduction(+:val)
            for (int j = 0; j < gs; j++) {
              val += val4[j] * val5[j] * val6[j];
            }
            #pragma acc loop reduction(+:val)
            for (int j = 0; j < gs; j++) {
              val += val7[j] * val8[j] * val9[j];
            }
            C[mb * N2 + nb * N + pb] = val;
            C[mb * N2 + pb * N + nb] = val;
          }
        } // end 3 atoms
      }
      #pragma omp atomic capture
      task = tasks++;
    } // while loop for manager worker
  } // omp parallel region

  FP2* C_all = new FP2[N2a]();

  for (int n=0;n<nomp;n++)
  {
    #if USE_ACC
    acc_set_device_num(n,acc_device_nvidia);
    #endif

    #pragma acc update self(C[0:N2a])

    for (int i=0;i<N2a;i++)
      C_all[i] += C[i];

    #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
    if (n>0)
    {
      #pragma acc exit data delete(C[0:N2a])
    }
  }
  #if USE_ACC
  acc_set_device_num(0,acc_device_nvidia);
  #endif

  for (int i=0;i<N2a;i++)
    C[i] = C_all[i];

  delete [] C_all;
  //apply symmetry and normalization
  copy_symm(natoms,N,Naux,basis,basis_aux,C,1);

  for (int i=0;i<Naux;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    C[i*N2+j*N+k] *= norm_sv(basis_aux[i][0],basis_aux[i][1],basis_aux[i][2],basis_aux[i][3])*basis[j][4]*basis[k][4];

  if (prl>2)
  {
    printf("\n C: \n");
    for (int i=0;i<Naux;i++)
    {
      printf(" i: %i \n",i);
      for (int j=0;j<N;j++)
      for (int k=0;k<N;k++)
        printf("  %12.6f",C[i*N2+j*N+k]);
      printf("\n");
    }
  }

  transpose_C(Naux,N,C);
  #pragma acc update device(C[0:N2a])

  #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = omp_get_thread_num();
    #if USE_ACC
    acc_set_device_num(tid,acc_device_nvidia);
    #endif

    #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
    #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs])
    #pragma acc exit data delete(grid3[0:gs6],wt3[0:gs])

    #pragma acc exit data delete(val1[0:gs],val2[0:gs],val3[0:gs]) 
    #pragma acc exit data delete(val4[0:gs],val5[0:gs],val6[0:gs]) 
    #pragma acc exit data delete(val7[0:gs],val8[0:gs],val9[0:gs]) 

    #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6]) 
    #pragma acc exit data delete(grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6]) 

    #pragma acc exit data delete(wtt1[0:gs],wtt2[0:gs])
    #pragma acc exit data delete(n2i[0:natoms],na2i[0:natoms])
  }

  delete [] n2i;
  delete [] na2i;

  delete [] ang_g;
  delete [] ang_w;

  delete [] grid1s;
  delete [] grid2s;
  delete [] grid3s;
  delete [] grid1p;
  delete [] grid2p;
  delete [] grid3p;

  delete [] wtt1;
  delete [] wtt2;
  delete [] grid1;
  delete [] grid2;
  delete [] grid3;
  delete [] wt1;
  delete [] wt2;
  delete [] wt3;
  
  delete [] val1;
  delete [] val2;
  delete [] val3;
  delete [] val4;
  delete [] val5;
  delete [] val6;
  delete [] val7;
  delete [] val8;
  delete [] val9;

  return;
}


#if RED_DOUBLE
void compute_all_3c_v2(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* C, int prl)
#else
void compute_all_3c_v2(int natoms, int* atno, FP1* coords, 
                       vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP1* C, int prl)
#endif
{
  if (prl>1) printf(" beginning compute_all_3c_v2 \n");

  int N = basis.size();
  int N2 = N*N;
  int Naux = basis_aux.size();
  int N2a = N2*Naux;
  int gs = nrad*nang;
  int gs6 = 6*gs;

  int estart = find_center_of_grid(1,nrad)*nang;
 
  int* na2i = new int[natoms];
  int iNa = get_imax_n2i(natoms,Naux,basis_aux,na2i);

  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);

  //printf("  iN/a: %i %i \n",iN,iNa);

  FP1* grid1 = new FP1[gs6];
  FP1* wt1 = new FP1[gs];

  FP1* grid2 = new FP1[gs6];
  FP1* wt2 = new FP1[gs];

  FP1* grid3 = new FP1[gs6];
  FP1* wt3 = new FP1[gs];

  FP1** val1 = new FP1*[iNa];
  FP1** val2 = new FP1*[iN];
  FP1** val3 = new FP1*[iN];
  FP1** val4 = new FP1*[iNa];
  FP1** val5 = new FP1*[iN];
  FP1** val6 = new FP1*[iN];
  FP1** val7 = new FP1*[iNa];
  FP1** val8 = new FP1*[iN];
  FP1** val9 = new FP1*[iN];
  for (int i=0;i<iNa;i++) val1[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val2[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val3[i] = new FP1[gs];
  for (int i=0;i<iNa;i++) val4[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val5[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val6[i] = new FP1[gs];
  for (int i=0;i<iNa;i++) val7[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val8[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val9[i] = new FP1[gs];
  FP1* valt1 = new FP1[gs];
  FP1* valt2 = new FP1[gs];
  FP1* valt3 = new FP1[gs];

  FP1* grid1s = new FP1[gs6];
  FP1* grid2s = new FP1[gs6];
  FP1* grid3s = new FP1[gs6];

  FP1* grid1p = new FP1[gs6];
  FP1* grid2p = new FP1[gs6];
  FP1* grid3p = new FP1[gs6];

  FP1* wtt1 = new FP1[gs];
  FP1* wtt2 = new FP1[gs];

  FP1* ang_g = new FP1[3*nang];
  FP1* ang_w = new FP1[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];


 #if USE_ACC
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data copyin(n2i[0:natoms],na2i[0:natoms])

  #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
  #pragma acc enter data create(grid2[0:gs6],wt2[0:gs])
  #pragma acc enter data create(grid3[0:gs6],wt3[0:gs])
  #pragma acc enter data create(val1[0:iNa][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs]) 
  #pragma acc enter data create(val4[0:iNa][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs]) 
  #pragma acc enter data create(val7[0:iNa][0:gs],val8[0:iN][0:gs],val9[0:iN][0:gs]) 

  #pragma acc enter data create(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6])
  #pragma acc enter data create(grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6])

  #pragma acc enter data create(wtt1[0:gs],wtt2[0:gs],valt1[0:gs],valt2[0:gs],valt3[0:gs])
  //#pragma acc enter data create(C[0:N2a])
 #endif
  acc_assign(N2a,C,0.);

#if TIMER_KD
  size_t vinr_time = 0;
  size_t red_time = 0;
  auto t1 = chrono::high_resolution_clock::now();
  auto t2 = chrono::high_resolution_clock::now();
#endif

  for (int m=0;m<natoms;m++)
  {
    
    int s1 = 0; if (m>0) s1 = na2i[m-1]; int s2 = na2i[m];
    int s3 = 0; if (m>0) s3 = n2i[m-1]; int s4 = n2i[m];
    //printf("  m: %i  s1/2->3/4: %i %i - %i %i \n",m,s1,s2,s3,s4);

    FP1 Z1 = (FP1)atno[m];
    FP1 A1 = coords[3*m+0]; FP1 B1 = coords[3*m+1]; FP1 C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);
    //#pragma acc update self(grid1[0:gs6],wt1[0:gs])
    //print_grid(gs,grid1,NULL,wt1,NULL,prl);

    //CPMZ need more of this
    #pragma acc parallel loop present(val1[0:iNa][0:gs])
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val1[ii1][j] = 1.f;
    }

    #pragma acc parallel loop present(val2[0:iN][0:gs])
    for (int ii2=0;ii2<s4-s3;ii2++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val2[ii2][j] = 1.f;
    }
    #pragma acc parallel loop present(val3[0:iN][0:gs],wt1[0:gs])
    for (int ii3=0;ii3<s4-s3;ii3++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val3[ii3][j] = wt1[j];
    }
    
    #if TIMER_KD
    t1 = chrono::high_resolution_clock::now();
    #endif
    //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<FP2> basis1 = basis_aux[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

      eval_inr_r12(gs,grid1,val1[ii1],n1,l1,zeta1);
      eval_sh_3r(gs,grid1,val1[ii1],n1,l1,m1);
    } //loop i1
    
    #if TIMER_KD
    t2 = chrono::high_resolution_clock::now();
    vinr_time += chrono::duration_cast<chrono::nanoseconds>(t2-t1).count();
    #endif

    for (int i2=s3;i2<s4;i2++)
    {
      int ii2 = i2-s3;

      vector<FP2> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

      eval_sh(ii2,gs,grid1,val2[ii2],n2,l2,m2,zeta2);
    } //loop i2

    for (int i3=s3;i3<s4;i3++)
    {
      int ii3 = i3-s3;

      vector<FP2> basis3 = basis[i3];
      int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];

      eval_sh(ii3,gs,grid1,val3[ii3],n3,l3,m3,zeta3);
    } //loop i3

    
    #if TIMER_KD
    t1 = chrono::high_resolution_clock::now();
    #endif

    #if REDUCEV==1
    reduce_3c1(s1,s2,s3,s4,gs,val1,val2,val3,valt1,N,Naux,iN,iNa,C);
    #else
    reduce_3c1b(s1,s2,s3,s4,gs,val1,val2,val3,N,Naux,iN,iNa,C);
    #endif

    #if TIMER_KD
    t2 = chrono::high_resolution_clock::now();
    red_time += chrono::duration_cast<chrono::nanoseconds>(t2-t1).count();;
    #endif


    //two-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      
     //s12 over atom m, aux function
     //s34 over atom m, regular function
     //s56 over atom n, regular function
      int s5 = 0; if (n>0) s5 = n2i[n-1]; int s6 = n2i[n];
      //printf("  n: %i  s1/2->3/4->5/6: %i %i - %i %i - %i %i \n",n,s1,s2,s3,s4,s5,s6);

      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2);
      recenter_grid(gs,grid2,A12,B12,C12);

      copy_grid(gs,grid1s,grid1); 
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid2,wt2,Z1,Z2,A12,B12,C12);
      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wt2);

      add_r1_to_grid(gs,grid2,0.,0.,0.);

      //#pragma acc update self(grid1[0:gs6],wt1[0:gs],grid2[0:gs6],wt2[0:gs])
      //print_grid(gs,grid1,grid2,wt1,wt2,prl);

      #pragma acc parallel loop collapse(2) present(val4[0:iNa][0:gs])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
        for (int j=0;j<gs;j++)
          val4[ii1][j] = 1.f;
      }

      #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs],val5[0:iN][0:gs])
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
        for (int j=0;j<gs;j++)
        {
          val2[ii2][j] = 1.f;
          val5[ii2][j] = 1.f;
        }
      }

      #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],wtt1[0:gs],wt2[0:gs])
      for (int ii3=0;ii3<s6-s5;ii3++)
      {
        for (int j=0;j<gs;j++)
        {
          val3[ii3][j] = wtt1[j];
          val6[ii3][j] = wt2[j];
        }
      }

      #if TIMER_KD
      t1 = chrono::high_resolution_clock::now();
      #endif
      //i1 on atom m
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis_aux[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];
        //printf("  m: %i i1: %i   nlm: %i %i %2i zeta: %8.5f \n",m,i1,n1,l1,m1,zeta1);

        eval_inr_r12(gs,grid2,val4[ii1],n1,l1,zeta1);
        eval_sh_3r(gs,grid2,val4[ii1],n1,l1,m1);
      }
      
      #if TIMER_KD
      t2 = chrono::high_resolution_clock::now();
      vinr_time += chrono::duration_cast<chrono::nanoseconds>(t2-t1).count();
      #endif

      //i2 on atom m
      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;

        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];
        //printf("   n: %i i2: %i   nlm: %i %i %2i zeta: %8.5f \n",n,i2,n2,l2,m2,zeta2);

        eval_sh(ii2,gs,grid1,val2[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2,val5[ii2],n2,l2,m2,zeta2);
      }

      //i3 on atom n
      for (int i3=s5;i3<s6;i3++)
      {
        int ii3 = i3-s5;

        vector<FP2> basis3 = basis[i3];
        int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];
    
        eval_sh(ii3,gs,grid1s,val3[ii3],n3,l3,m3,zeta3);
        eval_sh(ii3,gs,grid2s,val6[ii3],n3,l3,m3,zeta3);
      }

      
      #if TIMER_KD
      t1 = chrono::high_resolution_clock::now();
      #endif

      #if REDUCEV==1
      reduce_3c2(s1,s2,s3,s4,s5,s6,gs,val1,val2,val3,val4,val5,val6,valt1,valt2,N,Naux,iN,iNa,C);
      #else
      reduce_3c2b(s1,s2,s3,s4,s5,s6,gs,val1,val2,val3,val4,val5,val6,N,Naux,iN,iNa,C);
      #endif

      #if TIMER_KD
      t2 = chrono::high_resolution_clock::now();
      red_time += chrono::duration_cast<chrono::nanoseconds>(t2-t1).count();
      #endif
      

      #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs],val5[0:iN][0:gs])
      for (int ii2=0;ii2<s6-s5;ii2++)
      {
        for (int j=0;j<gs;j++)
        {
          val2[ii2][j] = 1.f;
          val5[ii2][j] = 1.f;
        }
      }

      #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],wtt1[0:gs],wt2[0:gs])
      for (int ii3=0;ii3<s6-s5;ii3++)
      {
        for (int j=0;j<gs;j++)
        {
          val3[ii3][j] = wtt1[j];
          val6[ii3][j] = wt2[j];
        }
      }

      //now do i2+i3 on atom n 
      //i2 on atom n
      for (int i2=s5;i2<s6;i2++)
      {
        int ii2 = i2-s5;

        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];
        //printf("   n: %i i2: %i   nlm: %i %i %2i zeta: %8.5f \n",n,i2,n2,l2,m2,zeta2);

        eval_sh(ii2,gs,grid1s,val2[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2s,val5[ii2],n2,l2,m2,zeta2);
      }

      //i3 on atom n
      for (int i3=s5;i3<s6;i3++)
      {
        int ii3 = i3-s5;

        vector<FP2> basis3 = basis[i3];
        int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];
    
        eval_sh(ii3,gs,grid1s,val3[ii3],n3,l3,m3,zeta3);
        eval_sh(ii3,gs,grid2s,val6[ii3],n3,l3,m3,zeta3);
      }

      #if TIMER_KD
      t1 = chrono::high_resolution_clock::now();
      #endif

      #if REDUCEV==1
      reduce_3c2(s1,s2,s5,s6,s5,s6,gs,val1,val2,val3,val4,val5,val6,valt1,valt2,N,Naux,iN,iNa,C);
      #else
      reduce_3c2b(s1,s2,s5,s6,s5,s6,gs,val1,val2,val3,val4,val5,val6,N,Naux,iN,iNa,C);
      #endif

      #if TIMER_KD
      t2 = chrono::high_resolution_clock::now();
      red_time += chrono::duration_cast<chrono::nanoseconds>(t2-t1).count();;
      #endif

    } //loop n over second atom



    //three-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

      //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1);
      recenter_grid(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      //testing condition
      for (int p=n+1;p<natoms;p++)
      if (p!=m)
      //for (int p=0;p<natoms;p++)
      //if (p!=m && p!=n)
      {
        int s5 = 0; if (p>0) s5 = n2i[p-1]; int s6 = n2i[p];
        //printf("  p: %i  s1/2->3/4->5/6: %i %i - %i %i - %i %i \n",p,s1,s2,s3,s4,s5,s6);

        FP1 Z3 = (FP1)atno[p];
        FP1 A3 = coords[3*p+0]; FP1 B3 = coords[3*p+1]; FP1 C3 = coords[3*p+2];
        FP1 A13 = A3-A1; FP1 B13 = B3-B1; FP1 C13 = C3-C1;

        generate_central_grid_2(grid3,wt3,Z3,nrad,nang,ang_g,ang_w);
        copy_grid(gs,grid3p,grid3); //grid 3 centered on atom 3
        recenter_grid(gs,grid3,A13,B13,C13); //grid 3 centered on atom 1

        copy_grid(gs,grid3s,grid3);
        recenter_grid(gs,grid3s,-A12,-B12,-C12); //grid 3 centered on atom 2

        copy_grid(gs,grid2p,grid2); //grid 2 centered on atom 3
        recenter_grid(gs,grid2p,-A13,-B13,-C13);

        copy_grid(gs,grid1p,grid1);
        recenter_grid(gs,grid1p,-A13,-B13,-C13); //grid 1 centered on atom 3

        add_r1_to_grid_6z(gs,grid1s,grid2s,grid3s,grid1p,grid2p,grid3p);

        //need to keep all of these distances in order
        add_r3_to_grid(gs,grid1,A13,B13,C13);
        add_r123_to_grid(gs,grid2,A12,B12,C12,0.,0.,0.,A13,B13,C13);
        add_r123_to_grid(gs,grid3,A13,B13,C13,A12,B12,C12,0.,0.,0.);

        acc_copyf(gs,wtt1,wt1,wtt2,wt2);

        becke_weight_3c(gs,grid1,wtt1,grid2,wtt2,grid3,wt3,Z1,Z2,Z3,A12,B12,C12,A13,B13,C13);
        eliminate_small_wt_3(estart,gs,wtt1,wtt2,wt3);

        add_r1_to_grid(gs,grid2,0.,0.,0.);
        add_r1_to_grid(gs,grid3,0.,0.,0.);


        #pragma acc parallel loop collapse(2) present(val4[0:iNa][0:gs],val7[0:iNa][0:gs])
        for (int ii1=0;ii1<s2-s1;ii1++)
        {
          for (int j=0;j<gs;j++)
          {
            val4[ii1][j] = 1.f;
            val7[ii1][j] = 1.f;
          }
        }

        #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs],val5[0:iN][0:gs],val8[0:iN][0:gs])
        for (int ii2=0;ii2<s4-s3;ii2++)
        {
          for (int j=0;j<gs;j++)
          {
            val2[ii2][j] = 1.f;
            val5[ii2][j] = 1.f;
            val8[ii2][j] = 1.f;
          }
        }

        #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],val9[0:iN][0:gs],wtt1[0:gs],wtt2[0:gs],wt3[0:gs])
        for (int ii3=0;ii3<s6-s5;ii3++)
        {
          for (int j=0;j<gs;j++)
          {
            val3[ii3][j] = wtt1[j];
            val6[ii3][j] = wtt2[j];
            val9[ii3][j] = wt3[j];
          }
        }

        #if TIMER_KD
        t1 = chrono::high_resolution_clock::now();
        #endif
        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<FP2> basis1 = basis_aux[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

          eval_inr_r12(gs,grid2,val4[ii1],n1,l1,zeta1);
          eval_sh_3r(gs,grid2,val4[ii1],n1,l1,m1);
          eval_inr_r12(gs,grid3,val7[ii1],n1,l1,zeta1);
          eval_sh_3r(gs,grid3,val7[ii1],n1,l1,m1);
        }
        #if TIMER_KD
        t2 = chrono::high_resolution_clock::now();
        vinr_time += chrono::duration_cast<chrono::nanoseconds>(t2-t1).count();
        #endif
        
        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<FP2> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

          eval_sh(ii2,gs,grid3s,val8[ii2],n2,l2,m2,zeta2);
          eval_sh(ii2,gs,grid2s,val5[ii2],n2,l2,m2,zeta2);
          eval_sh(ii2,gs,grid1s,val2[ii2],n2,l2,m2,zeta2);
        }

        for (int i3=s5;i3<s6;i3++)
        {
          int ii3 = i3-s5;
          vector<FP2> basis3 = basis[i3];
          int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];
    
          eval_sh(ii3,gs,grid3p,val9[ii3],n3,l3,m3,zeta3);
          eval_sh(ii3,gs,grid2p,val6[ii3],n3,l3,m3,zeta3);
          eval_sh(ii3,gs,grid1p,val3[ii3],n3,l3,m3,zeta3);
        }

        
        #if TIMER_KD
        t1 = chrono::high_resolution_clock::now();
        #endif

        #if REDUCEV==1
        reduce_3c3(s1,s2,s3,s4,s5,s6,gs,val1,val2,val3,val4,val5,val6,val7,val8,val9,valt1,valt2,valt3,N,Naux,iN,iNa,C);
        #else
        reduce_3c3b(s1,s2,s3,s4,s5,s6,gs,val1,val2,val3,val4,val5,val6,val7,val8,val9,N,Naux,iN,iNa,C);
        #endif

        #if TIMER_KD
        t2 = chrono::high_resolution_clock::now();
        red_time += chrono::duration_cast<chrono::nanoseconds>(t2-t1).count();;
        #endif

      } //loop p over third atom
    } //loop n over second atom

  } //loop m over natoms

  #if USE_ACC
    #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
    //#pragma acc exit data copyout(C[0:N2a])
    #pragma acc update self(C[0:N2a])
  #endif

  #if TIMER_KD
  printf("\nVINR time: %5.3e s REDUCE time: %5.3e s\n\n",(FP2)vinr_time/1.e9,(FP2)red_time/1.e9);
  #endif

  copy_symm(natoms,N,Naux,basis,basis_aux,C,1);

  for (int i=0;i<Naux;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    C[i*N2+j*N+k] *= norm_sv(basis_aux[i][0],basis_aux[i][1],basis_aux[i][2],basis_aux[i][3])*basis[j][4]*basis[k][4];

  if (prl>2)
  {
    printf("\n C: \n");
    for (int i=0;i<Naux;i++)
    {
      printf(" i: %i \n",i);
      for (int j=0;j<N;j++)
      for (int k=0;k<N;k++)
        printf("  %12.6f",C[i*N2+j*N+k]);
      printf("\n");
    }
  }

  transpose_C(Naux,N,C);
  #pragma acc update device(C[0:N2a])


  //CPMZ check this
  #if USE_ACC
    #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
    #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs])
    #pragma acc exit data delete(grid3[0:gs6],wt3[0:gs])

    #pragma acc exit data delete(val1[0:iNa][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs]) 
    #pragma acc exit data delete(val4[0:iNa][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs]) 
    #pragma acc exit data delete(val7[0:iNa][0:gs],val8[0:iN][0:gs],val9[0:iN][0:gs]) 

    #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6]) 
    #pragma acc exit data delete(grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6]) 

    #pragma acc exit data delete(wtt1[0:gs],wtt2[0:gs],valt1[0:gs],valt2[0:gs],valt3[0:gs])
    #pragma acc exit data delete(n2i[0:natoms],na2i[0:natoms])
  #endif

  delete [] n2i;
  delete [] na2i;

  delete [] ang_g;
  delete [] ang_w;

  delete [] grid1s;
  delete [] grid2s;
  delete [] grid3s;
  delete [] grid1p;
  delete [] grid2p;
  delete [] grid3p;

  delete [] wtt1;
  delete [] wtt2;
  delete [] grid1;
  delete [] grid2;
  delete [] grid3;
  delete [] wt1;
  delete [] wt2;
  delete [] wt3;
  for (int i=0;i<iNa;i++) delete [] val1[i];
  for (int i=0;i<iN;i++) delete [] val2[i];
  for (int i=0;i<iN;i++) delete [] val3[i];
  for (int i=0;i<iNa;i++) delete [] val4[i];
  for (int i=0;i<iN;i++) delete [] val5[i];
  for (int i=0;i<iN;i++) delete [] val6[i];
  for (int i=0;i<iNa;i++) delete [] val7[i];
  for (int i=0;i<iN;i++) delete [] val8[i];
  for (int i=0;i<iN;i++) delete [] val9[i];
  delete [] val1;
  delete [] val2;
  delete [] val3;
  delete [] val4;
  delete [] val5;
  delete [] val6;
  delete [] val7;
  delete [] val8;
  delete [] val9;
  delete [] valt1;
  delete [] valt2;
  delete [] valt3;

  return;
}

void compute_all_3c(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* C, int prl) //FP1* C changed to FP2* C by Soumi
{
  if (prl>1) printf(" beginning compute_all_3c \n");

  int N = basis.size();
  int N2 = N*N;
  int Naux = basis_aux.size();
  int N2a = N2*Naux;

  int gs = nrad*nang;

  FP1* grid1 = new FP1[6*gs];
  FP1* wt1 = new FP1[gs];
  FP1* val1 = new FP1[gs];

  FP1* grid2 = new FP1[6*gs];
  FP1* wt2 = new FP1[gs];
  FP1* val2 = new FP1[gs];

  FP1* grid3 = new FP1[6*gs];
  FP1* wt3 = new FP1[gs];
  FP1* val3 = new FP1[gs];

  FP1* grid1s = new FP1[6*gs];
  FP1* grid2s = new FP1[6*gs];
  FP1* grid3s = new FP1[6*gs];

  FP1* grid1p = new FP1[6*gs];
  FP1* grid2p = new FP1[6*gs];
  FP1* grid3p = new FP1[6*gs];

  FP1* valt1 = new FP1[gs];
  FP1* valt2 = new FP1[gs];
  FP1* valt3 = new FP1[gs];
  FP1* valt4 = new FP1[gs];
  FP1* valt5 = new FP1[gs];
  FP1* valt6 = new FP1[gs];
  FP1* wtt1 = new FP1[gs];
  FP1* wtt2 = new FP1[gs];

  FP1* ang_g = new FP1[3*nang];
  FP1* ang_w = new FP1[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

  int imaxN = 1;
  int* n2i = new int[natoms];
  int wa = 0;
  for (int i=0;i<N;i++)
  {
    int wa2 = basis[i][9];
    if (wa2!=wa)
    {
      int cmaxN = wa2-wa; 
      if (cmaxN>imaxN) imaxN = cmaxN;
      n2i[wa] = i;
      wa = wa2;
    }
  }
  n2i[natoms-1] = N;

 #if USE_ACC
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data copyin(n2i[0:natoms])

  #pragma acc enter data create(grid1[0:6*gs],wt1[0:gs],val1[0:gs])
  #pragma acc enter data create(grid2[0:6*gs],wt2[0:gs],val2[0:gs]) 
  #pragma acc enter data create(grid3[0:6*gs],wt3[0:gs],val3[0:gs]) 

  #pragma acc enter data create(grid1s[0:6*gs],grid2s[0:6*gs],grid3s[0:6*gs])
  #pragma acc enter data create(grid1p[0:6*gs],grid2p[0:6*gs],grid3p[0:6*gs])

  #pragma acc enter data create(valt1[0:gs],valt2[0:gs],valt3[0:gs],valt4[0:gs],valt5[0:gs],valt6[0:gs])
  #pragma acc enter data create(wtt1[0:gs],wtt2[0:gs])
  #pragma acc enter data create(C[0:N2a])
 #endif
  acc_assign(N2a,C,0.);

  for (int m=0;m<natoms;m++)
  {
    //int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    FP1 Z1 = (FP1)atno[m];
    FP1 A1 = coords[3*m+0]; FP1 B1 = coords[3*m+1]; FP1 C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);
    //generate_central_grid(grid1,wt1,val1,0,Z1,1,0,1.,nrad,nang,ang_g,ang_w);
    //#pragma acc update self(grid1[0:6*gs],wt1[0:gs])
    //print_grid(gs,grid1,NULL,wt1,NULL,prl);

   //first compute single atom ints
    for (int i1=0;i1<Naux;i1++)
    if (basis_aux[i1][9]==m)
    {
      vector<FP2> basis1 = basis_aux[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];
      //printf("\n  m: %i i1: %2i   nlm: %i %i %2i zeta: %8.5f \n",m,i1,n1,l1,m1,zeta1);

      acc_assign(gs,val1,1.);

      eval_inr_r12(gs,grid1,val1,n1,l1,zeta1);
      eval_sh_3r(gs,grid1,val1,n1,l1,m1);

      //#pragma acc update self(val1[0:gs])
      //print_array(gs,val1);

      for (int i2=0;i2<N;i2++)
      if (basis[i2][9]==m)
      {
        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];
        //printf("    n: %i i2: %2i   nlm: %i %i %2i zeta: %8.5f   (b2) \n",m,i2,n2,l2,m2,zeta2);

        acc_copyf(gs,valt1,val1);

        eval_sh(i2,gs,grid1,valt1,n2,l2,m2,zeta2);

        for (int i3=0;i3<N;i3++)
        if (basis[i3][9]==m)
        {
          vector<FP2> basis3 = basis[i3];
          int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];

          acc_copyf(gs,valt2,valt1);

          eval_sh(i2,gs,grid1,valt2,n3,l3,m3,zeta3);

          FP1 val = 0.;
        #if USE_ACC
         #pragma acc parallel loop independent present(valt2[0:gs],wt1[0:gs]) reduction(+:val)
        #endif
          for (int j=0;j<gs;j++)
            val += valt2[j] * wt1[j];
 
        #if USE_ACC
         #pragma acc serial present(C[0:N2a])
        #endif
          C[i1*N2+i2*N+i3] = val;

          //printf(" m: %i  i1/2/3: %2i %2i %2i val: %8.5f \n",m,i1,i2,i3,val);

        } //loop i3
      } //loop i2
    } //loop i1

   //two-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      //generate_central_grid(grid2,wt2,val2,0,Z2,1,0,1.,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2);
      recenter_grid(gs,grid2,A12,B12,C12);

      copy_grid(gs,grid1s,grid1); 
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid2,wt2,Z1,Z2,A12,B12,C12);
      eliminate_small_wt(gs,wtt1);
      eliminate_small_wt(gs,wt2);

      add_r1_to_grid(gs,grid2,0.,0.,0.);

      //#pragma acc update self(grid1[0:6*gs],wt1[0:gs],grid2[0:6*gs],wt2[0:gs])
      //print_grid(gs,grid1,grid2,wt1,wt2,prl);

      for (int i1=0;i1<Naux;i1++)
      if (basis_aux[i1][9]==m)
      {
        vector<FP2> basis1 = basis_aux[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];
        //printf("  m: %i i1: %i   nlm: %i %i %2i zeta: %8.5f \n",m,i1,n1,l1,m1,zeta1);

        acc_assign(gs,val1,1.);
        acc_assign(gs,val2,1.);

        eval_inr_r12(gs,grid1,val1,n1,l1,zeta1);
        eval_sh_3r(gs,grid1,val1,n1,l1,m1);
        eval_inr_r12(gs,grid2,val2,n1,l1,zeta1);
        eval_sh_3r(gs,grid2,val2,n1,m1,l1);

        for (int i2=0;i2<N;i2++)
        if (basis[i2][9]==n)
        {
          vector<FP2> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];
          //printf("   n: %i i2: %i   nlm: %i %i %2i zeta: %8.5f \n",n,i2,n2,l2,m2,zeta2);

          acc_copyf(gs,valt1,val1);
          acc_copyf(gs,valt2,val2);

          eval_sh(i2,gs,grid2s,valt2,n2,l2,m2,zeta2);
          eval_sh(i2,gs,grid1s,valt1,n2,l2,m2,zeta2);

         //third AO function on same center as aux function
          for (int i3=0;i3<N;i3++)
          if (basis[i3][9]==m) 
          {
            vector<FP2> basis3 = basis[i3];
            int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];
    
            acc_copyf(gs,valt3,valt1);
            acc_copyf(gs,valt4,valt2);

            eval_sh(i3,gs,grid2,valt4,n3,l3,m3,zeta3);
            eval_sh(i3,gs,grid1,valt3,n3,l3,m3,zeta3);

            FP1 val = 0.;
          #if USE_ACC
           #pragma acc parallel loop independent present(valt3[0:gs],wtt1[0:gs]) reduction(+:val)
          #endif
            for (int j=0;j<gs;j++)
              val += valt3[j] * wtt1[j];
          #if USE_ACC
           #pragma acc parallel loop independent present(valt4[0:gs],wt2[0:gs]) reduction(+:val)
          #endif
            for (int j=0;j<gs;j++)
              val += valt4[j] * wt2[j];

          #if USE_ACC
           #pragma acc serial present(C[0:N2a])
          #endif
            C[i1*N2+i2*N+i3] = val;

            //printf(" 1. i1/2/3: %i %i %i / %i %i %i val: %5.3f \n",i1,i2,i3,m,n,m,val);
          } //loop i3 over third basis function

         //third AO function on same center as second AO function
          for (int i3=0;i3<N;i3++)
          if (basis[i3][9]==n) 
          {
            vector<FP2> basis3 = basis[i3];
            int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];
    
            acc_copyf(gs,valt3,valt1);
            acc_copyf(gs,valt4,valt2);

           //CPMZ check these
            eval_sh(i3,gs,grid2s,valt4,n3,l3,m3,zeta3);
            eval_sh(i3,gs,grid1s,valt3,n3,l3,m3,zeta3);

            FP1 val = 0.;
          #if USE_ACC
           #pragma acc parallel loop independent present(valt3[0:gs],wtt1[0:gs]) reduction(+:val)
          #endif
            for (int j=0;j<gs;j++)
              val += valt3[j] * wtt1[j];
          #if USE_ACC
           #pragma acc parallel loop independent present(valt4[0:gs],wt2[0:gs]) reduction(+:val)
          #endif
            for (int j=0;j<gs;j++)
              val += valt4[j] * wt2[j];

          #if USE_ACC
           #pragma acc serial present(C[0:N2a])
          #endif
            C[i1*N2+i2*N+i3] = val;

            //printf(" 2. i1/2/3: %i %i %i / %i %i %i val: %5.3f \n",i1,i2,i3,m,n,n,val);

          } //loop i3 over third basis function

        } //loop i2 over second basis function
      } //loop i1 over first basis function

    } //loop n over second atom
 
   //three-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      //generate_central_grid(grid2,wt2,val2,0,Z2,1,0,1.,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1);
      recenter_grid(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      for (int p=0;p<natoms;p++)
      if (p!=m && p!=n)
      {
        FP1 Z3 = (FP1)atno[p];
        FP1 A3 = coords[3*p+0]; FP1 B3 = coords[3*p+1]; FP1 C3 = coords[3*p+2];
        FP1 A13 = A3-A1; FP1 B13 = B3-B1; FP1 C13 = C3-C1;
        //FP1 A23 = A3-A2; FP1 B23 = B3-B2; FP1 C23 = C3-C2;

        generate_central_grid_2(grid3,wt3,Z3,nrad,nang,ang_g,ang_w);
        //generate_central_grid(grid3,wt3,val3,0,Z3,1,0,1.,nrad,nang,ang_g,ang_w);
        copy_grid(gs,grid3p,grid3); //grid 3 centered on atom 3
        recenter_grid(gs,grid3,A13,B13,C13); //grid 3 centered on atom 1

        copy_grid(gs,grid3s,grid3);
        recenter_grid(gs,grid3s,-A12,-B12,-C12); //grid 3 centered on atom 2

        copy_grid(gs,grid2p,grid2); //grid 2 centered on atom 3
        recenter_grid(gs,grid2p,-A13,-B13,-C13);

        copy_grid(gs,grid1p,grid1);
        recenter_grid(gs,grid1p,-A13,-B13,-C13); //grid 1 centered on atom 3

        add_r1_to_grid_6z(gs,grid1s,grid2s,grid3s,grid1p,grid2p,grid3p);

      //need to keep all of these distances in order
        add_r3_to_grid(gs,grid1,A13,B13,C13);
        add_r123_to_grid(gs,grid2,A12,B12,C12,0.,0.,0.,A13,B13,C13);
        add_r123_to_grid(gs,grid3,A13,B13,C13,A12,B12,C12,0.,0.,0.);

        acc_copyf(gs,wtt1,wt1);
        acc_copyf(gs,wtt2,wt2);

        becke_weight_3c(gs,grid1,wtt1,grid2,wtt2,grid3,wt3,Z1,Z2,Z3,A12,B12,C12,A13,B13,C13);
        eliminate_small_wt_3(gs,wtt1,wtt2,wt3);

        add_r1_to_grid(gs,grid2,0.,0.,0.);
        add_r1_to_grid(gs,grid3,0.,0.,0.);

        //#pragma acc update self(grid1[0:6*gs],wtt1[0:gs],grid2[0:6*gs],wtt2[0:gs],grid3[0:6*gs],wt3[0:gs])
        //print_grid(gs,grid1,grid2,wtt1,wtt2,prl);
        //print_grid(gs,grid3,NULL,wt3,NULL,prl);
 
        for (int i1=0;i1<Naux;i1++)
        if (basis_aux[i1][9]==m)
        {
          vector<FP2> basis1 = basis_aux[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];
          //printf("  m: %i i1: %i   nlm: %i %i %2i zeta: %8.5f \n",m,i1,n1,l1,m1,zeta1);

          acc_assign(gs,val1,1.);
          acc_assign(gs,val2,1.);
          acc_assign(gs,val3,1.);

          eval_inr_r12(gs,grid1,val1,n1,l1,zeta1);
          eval_sh_3r(gs,grid1,val1,n1,l1,m1);
          eval_inr_r12(gs,grid2,val2,n1,l1,zeta1);
          eval_sh_3r(gs,grid2,val2,n1,l1,m1);
          eval_inr_r12(gs,grid3,val3,n1,l1,zeta1);
          eval_sh_3r(gs,grid3,val3,n1,l1,m1);

          for (int i2=0;i2<N;i2++)
          if (basis[i2][9]==n)
          {
            vector<FP2> basis2 = basis[i2];
            int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];
            //printf("   n: %i i2: %i   nlm: %i %i %2i zeta: %8.5f \n",n,i2,n2,l2,m2,zeta2);

            acc_copyf(gs,valt1,val1);
            acc_copyf(gs,valt2,val2);
            acc_copyf(gs,valt3,val3);

            eval_sh(i2,gs,grid3s,valt3,n2,l2,m2,zeta2);
            eval_sh(i2,gs,grid2s,valt2,n2,l2,m2,zeta2);
            eval_sh(i2,gs,grid1s,valt1,n2,l2,m2,zeta2);

            for (int i3=0;i3<N;i3++)
            if (basis[i3][9]==p) 
            {
              vector<FP2> basis3 = basis[i3];
              int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];
    
              acc_copyf(gs,valt4,valt1);
              acc_copyf(gs,valt5,valt2);
              acc_copyf(gs,valt6,valt3);

              eval_sh(i3,gs,grid3p,valt6,n3,l3,m3,zeta3);
              eval_sh(i3,gs,grid2p,valt5,n3,l3,m3,zeta3);
              eval_sh(i3,gs,grid1p,valt4,n3,l3,m3,zeta3);

              FP1 val = 0.;
            #if USE_ACC
             #pragma acc parallel loop independent present(valt4[0:gs],wtt1[0:gs]) reduction(+:val)
            #endif
              for (int j=0;j<gs;j++)
                val += valt4[j] * wtt1[j];
            #if USE_ACC
             #pragma acc parallel loop independent present(valt5[0:gs],wtt2[0:gs]) reduction(+:val)
            #endif
              for (int j=0;j<gs;j++)
                val += valt5[j] * wtt2[j];
            #if USE_ACC
             #pragma acc parallel loop independent present(valt6[0:gs],wt3[0:gs]) reduction(+:val)
            #endif
              for (int j=0;j<gs;j++)
                val += valt6[j] * wt3[j];

            #if USE_ACC
             #pragma acc serial present(C[0:N2a])
            #endif
              C[i1*N2+i2*N+i3] = val;

              //printf(" 3c i1/2/3: %i %i %i / %i %i %i val: %5.3f \n",i1,i2,i3,m,n,m,val);
            } //loop i3
          } //loop i2
        } //loop i1

      } //loop p over third atom
    } //loop n over second atom

  } //loop m over natoms

 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc exit data copyout(C[0:N2a])
 #endif

#if DEBUG
  if (prl>0)
  {
    printf("\n C(raw): \n");
    for (int j=0;j<N;j++)
    for (int k=0;k<N;k++)
    {
      for (int i=0;i<Naux;i++)
        printf(" %8.5f",C[i*N2+j*N+k]);
      printf("\n");
    }
  }
#endif

  copy_symm(natoms,N,Naux,basis,basis_aux,C,1);

#if 0
  if (prl>1)
  {
    printf("\n C(raw): \n");
    for (int j=0;j<N;j++)
    for (int k=0;k<N;k++)
    {
      for (int i=0;i<Naux;i++)
        printf("  %12.6f",C[i*N2+j*N+k]);
      printf("\n");
    }
  }
#endif

  for (int i=0;i<Naux;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    C[i*N2+j*N+k] *= norm_sv(basis_aux[i][0],basis_aux[i][1],basis_aux[i][2],basis_aux[i][3])*basis[j][4]*basis[k][4];

  if (prl>2)
  {
    printf("\n C: \n");
    for (int i=0;i<Naux;i++)
    {
      printf(" i: %i \n",i);
      for (int j=0;j<N;j++)
      for (int k=0;k<N;k++)
        printf("  %12.6f",C[i*N2+j*N+k]);
      printf("\n");
    }
  }

 //CPMZ check this
#if USE_ACC
  #pragma acc exit data delete(grid1[0:6*gs],wt1[0:gs],val1[0:gs])
  #pragma acc exit data delete(grid2[0:6*gs],wt2[0:gs],val2[0:gs]) 
  #pragma acc exit data delete(grid3[0:6*gs],wt3[0:gs],val3[0:gs]) 

  #pragma acc exit data delete(grid1s[0:6*gs],grid2s[0:6*gs],grid3s[0:6*gs]) 
  #pragma acc exit data delete(grid1p[0:6*gs],grid2p[0:6*gs],grid3p[0:6*gs]) 

  #pragma acc exit data delete(valt1[0:gs],valt2[0:gs],valt3[0:gs],valt4[0:gs],valt5[0:gs],valt6[0:gs])
  #pragma acc exit data delete(wtt1[0:gs],wtt2[0:gs])
  #pragma acc exit data delete(C[0:N2a])
  #pragma acc exit data delete(n2i[0:natoms])
#endif

  delete [] n2i;

  delete [] ang_g;
  delete [] ang_w;

  delete [] grid1s;
  delete [] grid2s;
  delete [] grid3s;
  delete [] grid1p;
  delete [] grid2p;
  delete [] grid3p;

  delete [] valt1;
  delete [] valt2;
  delete [] valt3;
  delete [] wtt1;
  delete [] wtt2;
  delete [] grid1;
  delete [] grid2;
  delete [] grid3;
  delete [] wt1;
  delete [] wt2;
  delete [] wt3;
  delete [] val1;
  delete [] val2;
  delete [] val3;

  return;
}

void norm_2c(int N, vector<vector<FP2> > &basis, FP1* S)
{
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    S[i*N+j] *= basis[i][4]*basis[j][4];

  for (int i=0;i<N;i++)
  for (int j=i+1;j<N;j++)
    S[j*N+i] = S[i*N+j];

  return;
}

#if !EVL64
void norm_2c(int N, vector<vector<FP2> > &basis, FP2* S)
{
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    S[i*N+j] *= basis[i][4]*basis[j][4];

  for (int i=0;i<N;i++)
  for (int j=i+1;j<N;j++)
    S[j*N+i] = S[i*N+j];

  return;
}
#endif

#if RED_DOUBLE
void compute_VdV(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, int nc, FP1* coordsc, FP2* Pao, FP2* V, FP2* dV, int prl)
#else
void compute_VdV(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, int nc, FP1* coordsc, FP2* Pao, FP1* V, FP1* dV, int prl)
#endif
{
  if (nc<1) return;

  if (prl>1) printf(" beginning compute_VdV \n");

  int N = basis.size();
  int N2 = N*N;
  int nc3 = 3*nc;

  int gs = nrad*nang;
  int gs3 = 3*gs; int gs6 = 6*gs;

  int estart = find_center_of_grid(1,nrad)*nang;

  FP1* grid1 = new FP1[gs6];
  FP1* wt1 = new FP1[gs];

  FP1* grid2 = new FP1[gs6];
  FP1* wt2 = new FP1[gs];

  FP1* grid3 = new FP1[gs6];
  FP1* wt3 = new FP1[gs];

  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);
  //printf("  iN: %i \n",imaxN);

  FP2* norms = new FP2[N2];
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    norms[i*N+j] = basis[i][4]*basis[j][4];

 //intermediate storage
  int iN = imaxN;
  FP1* grid1s = new FP1[gs6]; FP1* grid2s = new FP1[gs6]; FP1* grid3s = new FP1[gs6];
  FP1* grid1p = new FP1[gs6]; FP1* grid2p = new FP1[gs6]; FP1* grid3p = new FP1[gs6];
  FP1** valS1 = new FP1*[iN]; for (int i=0;i<iN;i++) valS1[i] = new FP1[gs];
  FP1** valS2 = new FP1*[iN]; for (int i=0;i<iN;i++) valS2[i] = new FP1[gs];
  FP1** valS3 = new FP1*[iN]; for (int i=0;i<iN;i++) valS3[i] = new FP1[gs];
  FP1** valS4 = new FP1*[iN]; for (int i=0;i<iN;i++) valS4[i] = new FP1[gs];
  FP1** valS5 = new FP1*[iN]; for (int i=0;i<iN;i++) valS5[i] = new FP1[gs];
  FP1** valS6 = new FP1*[iN]; for (int i=0;i<iN;i++) valS6[i] = new FP1[gs];

  FP1** valt1 = new FP1*[iN]; for (int i=0;i<iN;i++) valt1[i] = new FP1[gs];
  FP1** valt2 = new FP1*[iN]; for (int i=0;i<iN;i++) valt2[i] = new FP1[gs];
  FP1** valt3 = new FP1*[iN]; for (int i=0;i<iN;i++) valt3[i] = new FP1[gs];

  FP1** valtv1 = new FP1*[iN]; for (int i=0;i<iN;i++) valtv1[i] = new FP1[gs3];
  FP1** valtv2 = new FP1*[iN]; for (int i=0;i<iN;i++) valtv2[i] = new FP1[gs3];
  FP1** valtv3 = new FP1*[iN]; for (int i=0;i<iN;i++) valtv3[i] = new FP1[gs3];

  FP1* wtt1 = new FP1[gs];
  FP1* wtt2 = new FP1[gs];

 #if RED_DOUBLE
  FP2* V1 = new FP2[nc];
  FP2* dV1 = new FP2[nc3]; 
 #else
  FP1* V1 = new FP1[nc];
  FP1* dV1 = new FP1[nc3]; 
 #endif

  FP1* ang_g = new FP1[3*nang];
  FP1* ang_w = new FP1[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

 #if USE_ACC
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data copyin(n2i[0:natoms],norms[0:N2])
  #pragma acc enter data copyin(coords[0:3*natoms],atno[0:natoms])
  #pragma acc enter data copyin(Pao[0:N2])

  #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
  #pragma acc enter data create(grid2[0:gs6],wt2[0:gs]) 
  #pragma acc enter data create(grid3[0:gs6],wt3[0:gs]) 
  #pragma acc enter data create(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs],valS5[0:iN][0:gs],valS6[0:iN][0:gs])
  #pragma acc enter data create(valt1[0:iN][0:gs],valt2[0:iN][0:gs],valt3[0:iN][0:gs])
  #pragma acc enter data create(valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3],valtv3[0:iN][0:gs3])
  #pragma acc enter data create(grid1s[0:gs6],grid1p[0:gs6],grid2s[0:gs6],grid2p[0:gs6],grid3s[0:gs6],grid3p[0:gs6])
  #pragma acc enter data create(wtt1[0:gs],wtt2[0:gs])
  #pragma acc enter data create(V[0:nc],dV[0:nc3],V1[0:nc],dV1[0:nc3])
 #endif
  acc_assign(nc,V,0.);
  acc_assign(nc3,dV,0.);

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    FP1 Z1 = (FP1)atno[m];
    FP1 A1 = coords[3*m+0]; FP1 B1 = coords[3*m+1]; FP1 C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);

  #if USE_ACC
   #pragma acc parallel loop collapse(2) present(valS1[0:iN][0:gs],valS3[0:iN][0:gs])
  #endif
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
      for (int j=0;j<gs;j++)
        valS1[ii1][j] = valS3[ii1][j] = 1.f;
    }

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<FP2> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

      eval_sh(ii1,gs,grid1,valS1[ii1],n1,l1,m1,zeta1);
    }

    for (int i2=s1;i2<s2;i2++)
    {
      int ii2 = i2-s1;
      vector<FP2> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

      eval_sh(ii2,gs,grid1,valS3[ii2],n2,l2,m2,zeta2);
    } //loop i2 evaluate

   #if USE_ACC
    #pragma acc wait
   #endif


   //2 basis on same center (m), 1 charge elsewhere (p)
    for (int p=0;p<nc;p++)
    {
      FP1 Zb = 1.;
      FP1 An = coordsc[3*p+0]; FP1 Bn = coordsc[3*p+1]; FP1 Cn = coordsc[3*p+2];
      FP1 A1n = An-A1; FP1 B1n = Bn-B1; FP1 C1n = Cn-C1;

      add_r2_to_grid(gs,grid1,A1n,B1n,C1n);

      generate_central_grid_2(grid3,wt3,Zb,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid3p,grid3); //grid 3 centered on charge
      recenter_grid(gs,grid3,A1n,B1n,C1n); //grid 3 centered on atom 1

      copy_grid(gs,grid1p,grid1);
      recenter_grid(gs,grid1p,-A1n,-B1n,-C1n); //grid 1 centered on charge

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid3,wt3,Z1,Zb,A1n,B1n,C1n);

      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wt3);

      add_r1_to_grid(gs,grid3,0.,0.,0.);
      add_r2_to_grid(gs,grid3,A1n,B1n,C1n);

      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc parallel loop present(valt1[0:iN][0:gs],valS1[0:iN][0:gs],wtt1[0:gs])
        for (int j=0;j<gs;j++)
          valt1[ii1][j] = valS1[ii1][j]*wtt1[j];

       #pragma acc parallel loop present(valS2[0:iN][0:gs],wt3[0:gs])
        for (int j=0;j<gs;j++)
          valS2[ii1][j] = wt3[j];

     //CPMZ check
     //  #pragma acc parallel loop present(valS3[0:iN][0:gs])
     //   for (int j=0;j<gs;j++)
     //     valS3[ii1][j] = 1.f;

       #pragma acc parallel loop present(valS4[0:iN][0:gs])
        for (int j=0;j<gs;j++)
          valS4[ii1][j] = 1.f;
      }


      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

        eval_sh(ii1,gs,grid3,valS2[ii1],n1,l1,m1,zeta1);
      }

      for (int i2=s1;i2<s2;i2++)
      {
        int ii2 = i2-s1;
        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

        eval_sh(ii2,gs,grid3,valS4[ii2],n2,l2,m2,zeta2);
      } //loop i2 evaluate


      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc parallel loop collapse(2) present(valtv1[0:iN][0:gs3],valt1[0:iN][0:gs])
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valtv1[ii1][3*j+k] = valt1[ii1][j];

       #pragma acc parallel loop collapse(2) present(valtv2[0:iN][0:gs3],valS2[0:iN][0:gs])
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valtv2[ii1][3*j+k] = valS2[ii1][j];
      }

     #pragma acc parallel loop present(grid1[0:gs6],grid3[0:gs6],grid1p[0:gs6],grid3p[0:gs6],valt1[0:iN][0:gs],valS2[0:iN][0:gs],valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3])
      for (int j=0;j<gs;j++)
      {
        FP1 x1 = grid1p[6*j+0]; FP1 y1 = grid1p[6*j+1]; FP1 z1 = grid1p[6*j+2];
        FP1 x3 = grid3p[6*j+0]; FP1 y3 = grid3p[6*j+1]; FP1 z3 = grid3p[6*j+2];
        FP1 Rn1 = grid1[6*j+4]+1.e-20f;
        FP1 Rn3 = grid3[6*j+4]+1.e-20f;
        FP1 r12 = Rn1*Rn1; FP1 r32 = Rn3*Rn3;

        FP1 ne1 = -1.f/Rn1; FP1 ne3 = -1.f/Rn3;
        FP1 dx1 = x1*ne1/r12; FP1 dy1 = y1*ne1/r12; FP1 dz1 = z1*ne1/r12;
        FP1 dx3 = x3*ne3/r32; FP1 dy3 = y3*ne3/r32; FP1 dz3 = z3*ne3/r32;

        //ne1 = ne3 = 1.;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valt1[ii1][j] *= ne1;
       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valS2[ii1][j] *= ne3;

       #pragma acc loop 
        for (int ii1=0;ii1<s2-s1;ii1++)
        {
          valtv1[ii1][3*j+0] *= dx1;
          valtv1[ii1][3*j+1] *= dy1;
          valtv1[ii1][3*j+2] *= dz1;
        }
       #pragma acc loop 
        for (int ii1=0;ii1<s2-s1;ii1++)
        {
          valtv2[ii1][3*j+0] *= dx3;
          valtv2[ii1][3*j+1] *= dy3;
          valtv2[ii1][3*j+2] *= dz3;
        }
      }

     //collect 2-atom values
      acc_assign(nc,V1,0.);
      acc_assign(nc3,dV1,0.);

      reduce_2c2v (p,s1,s2,gs,norms,Pao,valt1,valS2,valS3,valS4,iN,N,nc,1.,V1);
      reduce_2c2vd(p,s1,s2,gs,norms,Pao,valtv1,valtv2,valS3,valS4,iN,N,nc,1.,dV1);

     #pragma acc parallel loop present(V[0:nc],V1[0:nc])
      for (int j=0;j<nc;j++)
        V[j] += V1[j];
     #pragma acc parallel loop present(dV[0:nc3],dV1[0:nc3])
      for (int j=0;j<nc3;j++)
        dV[j] += dV1[j];

    } //loop p over nuclear center


   //two-center basis
    if (0)
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1); 
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

    #if USE_ACC
     #pragma acc parallel loop present(valS1[0:iN][0:gs],valS2[0:iN][0:gs])
    #endif
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS1[ii1][j] = valS2[ii1][j] = 1.f;
      }

    #if USE_ACC
     #pragma acc parallel loop present(valS3[0:iN][0:gs],valS4[0:iN][0:gs])
    #endif
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS3[ii2][j] = valS4[ii2][j] = 1.f;
      }

      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

        eval_sh(ii1,gs,grid1,valS1[ii1],n1,l1,m1,zeta1);
        eval_sh(ii1,gs,grid2,valS2[ii1],n1,l1,m1,zeta1);
      }

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1s,valS3[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2s,valS4[ii2],n2,l2,m2,zeta2);
      }


     //2 basis on diff centers, 1 charge elsewhere
      for (int p=0;p<nc;p++)
      {
        FP1 Zb = 1.;
        FP1 An = coordsc[3*p+0]; FP1 Bn = coordsc[3*p+1]; FP1 Cn = coordsc[3*p+2];
        FP1 A1n = An-A1; FP1 B1n = Bn-B1; FP1 C1n = Cn-C1;

        generate_central_grid_2(grid3,wt3,Zb,nrad,nang,ang_g,ang_w);
        copy_grid(gs,grid3p,grid3); //grid 3 centered on charge
        recenter_grid(gs,grid3,A1n,B1n,C1n); //grid 3 centered on atom 1

        copy_grid(gs,grid1p,grid1);
        recenter_grid(gs,grid1p,-A1n,-B1n,-C1n); //grid 1 centered on charge
        copy_grid(gs,grid2p,grid2);
        recenter_grid(gs,grid2p,-A1n,-B1n,-C1n); //grid 2 centered on charge

        copy_grid(gs,grid3s,grid3);
        recenter_grid(gs,grid3s,-A12,-B12,-C12); //grid 3 centered on atom 2

        //add_r1_to_grid_6z(gs,grid1s,grid2s,grid3s,grid1p,grid2p,grid3p);
        add_r1_to_grid(gs,grid3s,0.,0.,0.);
      
      //need to keep all of these distances in order
        //add_r3_to_grid(gs,grid1,A1n,B1n,C1n);
        add_r123_to_grid(gs,grid1,0.,0.,0.,A12,B12,C12,A1n,B1n,C1n);
        add_r123_to_grid(gs,grid2,A12,B12,C12,0.,0.,0.,A1n,B1n,C1n);
      	add_r123_to_grid(gs,grid3,A1n,B1n,C1n,A12,B12,C12,0.,0.,0.);
     
        acc_copyf(gs,wtt1,wt1,wtt2,wt2);

        becke_weight_3c(gs,grid1,wtt1,grid2,wtt2,grid3,wt3,Z1,Z2,Zb,A12,B12,C12,A1n,B1n,C1n);
        eliminate_small_wt_3(estart,gs,wtt1,wtt2,wt3);

        add_r2_to_grid(gs,grid1,A1n,B1n,C1n);
        add_r1_to_grid(gs,grid2,0.,0.,0.);
        add_r1_to_grid(gs,grid3,0.,0.,0.);
        add_r2_to_grid(gs,grid2,A1n,B1n,C1n);
        add_r2_to_grid(gs,grid3,A1n,B1n,C1n);

        for (int ii1=0;ii1<s2-s1;ii1++)
        {
         #pragma acc parallel loop present(valt1[0:iN][0:gs],valS1[0:iN][0:gs],wtt1[0:gs])
          for (int j=0;j<gs;j++)
            valt1[ii1][j] = valS1[ii1][j]*wtt1[j];

         #pragma acc parallel loop present(valt2[0:iN][0:gs],valS2[0:iN][0:gs],wtt2[0:gs])
          for (int j=0;j<gs;j++)
            valt2[ii1][j] = valS2[ii1][j]*wtt2[j];

         #pragma acc parallel loop present(valS5[0:iN][0:gs],wt3[0:gs])
          for (int j=0;j<gs;j++)
            valS5[ii1][j] = wt3[j];

         #pragma acc parallel loop collapse(2) present(valtv1[0:iN][0:gs3],valS1[0:iN][0:gs],wtt1[0:gs])
          for (int j=0;j<gs;j++)
          for (int k=0;k<3;k++)
            valtv1[ii1][3*j+k] = valS1[ii1][j]*wtt1[j];

         #pragma acc parallel loop collapse(2) present(valtv2[0:iN][0:gs3],valS2[0:iN][0:gs],wtt2[0:gs])
          for (int j=0;j<gs;j++)
          for (int k=0;k<3;k++)
            valtv2[ii1][3*j+k] = valS2[ii1][j]*wtt2[j];
        }

      #if USE_ACC
       #pragma acc parallel loop collapse(2) present(valS6[0:iN][0:gs])
      #endif
        for (int ii2=0;ii2<s4-s3;ii2++)
        {
          for (int j=0;j<gs;j++)
            valS6[ii2][j] = 1.f;
        }

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<FP2> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

          eval_sh(ii1,gs,grid3,valS5[ii1],n1,l1,m1,zeta1);
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<FP2> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

          eval_sh(ii2,gs,grid3s,valS6[ii2],n2,l2,m2,zeta2);
        } //loop i2 evaluate

       #pragma acc parallel loop collapse(2) present(valtv3[0:iN][0:gs3],valS5[0:iN][0:gs])
        for (int ii1=0;ii1<s2-s1;ii1++)
        {
          for (int j=0;j<gs;j++)
          {
            FP1 v1wt = valS5[ii1][j];
            valtv3[ii1][3*j+0] = v1wt;
            valtv3[ii1][3*j+1] = v1wt;
            valtv3[ii1][3*j+2] = v1wt;
          }
        }


       #pragma acc parallel loop present(grid1[0:gs6],grid2[0:gs6],grid3[0:gs6],grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6],valt1[0:iN][0:gs],valt2[0:iN][0:gs],valS5[0:iN][0:gs],valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3],valtv3[0:iN][0:gs3])
        for (int j=0;j<gs;j++)
        {
          FP1 x1 = grid1p[6*j+0]; FP1 y1 = grid1p[6*j+1]; FP1 z1 = grid1p[6*j+2];
          FP1 x2 = grid2p[6*j+0]; FP1 y2 = grid2p[6*j+1]; FP1 z2 = grid2p[6*j+2];
          FP1 x3 = grid3p[6*j+0]; FP1 y3 = grid3p[6*j+1]; FP1 z3 = grid3p[6*j+2];
          FP1 Rn1 = grid1[6*j+4]+1.e-20f;
          FP1 Rn2 = grid2[6*j+4]+1.e-20f;
          FP1 Rn3 = grid3[6*j+4]+1.e-20f;
          FP1 r12 = Rn1*Rn1; FP1 r22 = Rn2*Rn2; FP1 r32 = Rn3*Rn3;
          FP1 ne1 = -1.f/Rn1; FP1 ne2 = -1.f/Rn2; FP1 ne3 = -1.f/Rn3;
          FP1 dx1 = x1*ne1/r12; FP1 dy1 = y1*ne1/r12; FP1 dz1 = z1*ne1/r12;
          FP1 dx2 = x2*ne2/r22; FP1 dy2 = y2*ne2/r22; FP1 dz2 = z2*ne2/r22;
          FP1 dx3 = x3*ne3/r32; FP1 dy3 = y3*ne3/r32; FP1 dz3 = z3*ne3/r32;

          ne1 = ne2 = ne3 = 1.;

         #pragma acc loop
          for (int ii1=0;ii1<s2-s1;ii1++)
            valt1[ii1][j] *= ne1;
         #pragma acc loop
          for (int ii1=0;ii1<s2-s1;ii1++)
            valt2[ii1][j] *= ne2;
         #pragma acc loop
          for (int ii1=0;ii1<s2-s1;ii1++)
            valS5[ii1][j] *= ne3;

         #pragma acc loop
          for (int ii1=0;ii1<s2-s1;ii1++)
          {
            valtv1[ii1][3*j+0] *= dx1;
            valtv1[ii1][3*j+1] *= dy1;
            valtv1[ii1][3*j+2] *= dz1;
          }
         #pragma acc loop
          for (int ii1=0;ii1<s2-s1;ii1++)
          {
            valtv2[ii1][3*j+0] *= dx2;
            valtv2[ii1][3*j+1] *= dy2;
            valtv2[ii1][3*j+2] *= dz2;
          }
         #pragma acc loop
          for (int ii1=0;ii1<s2-s1;ii1++)
          {
            valtv3[ii1][3*j+0] *= dx3;
            valtv3[ii1][3*j+1] *= dy3;
            valtv3[ii1][3*j+2] *= dz3;
          }
        }

        //collect 3-atom values
        acc_assign(nc,V1,0.);
        acc_assign(nc3,dV1,0.);

        reduce_2c3v (p,s1,s2,s3,s4,gs,norms,Pao,valt1, valt2, valS5, valS3,valS4,valS6,iN,N,nc,1.,V1);
        reduce_2c3vd(p,s1,s2,s3,s4,gs,norms,Pao,valtv1,valtv2,valtv3,valS3,valS4,valS6,iN,N,nc,1.,dV1);

       #pragma acc parallel loop present(V[0:nc],V1[0:nc])
        for (int j=0;j<nc;j++)
          V[j] += V1[j];
       #pragma acc parallel loop present(dV[0:nc3],dV1[0:nc3])
        for (int j=0;j<nc3;j++)
          dV[j] += dV1[j];

      } //loop p over nuclear center

    } //loop n over second atom


  } //loop m over natoms

 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc exit data delete(V1[0:nc],dV1[0:nc3])
  #pragma acc exit data copyout(V[0:nc],dV[0:nc3])
 #endif

  if (prl>0)
  {
    printf("\n V(no nn): \n");
    for (int i=0;i<nc;i++)
      printf(" %10.5f",V[i]);
    printf("\n");
  }

  if (prl>0)
  {
    printf("\n dV(no nn): \n");
    for (int i=0;i<nc;i++)
    {
      for (int j=0;j<3;j++)
        printf(" %12.6f",dV[3*i+j]);
      printf("\n");
    }
    printf("\n");
  }

 //include nuclear contribution to V
  if (1)
  for (int i=0;i<natoms;i++)
  {
    FP1 Z1 = atno[i];
    FP1 A1 = coords[3*i+0]; FP1 B1 = coords[3*i+1]; FP1 C1 = coords[3*i+2];
    for (int p=0;p<nc;p++)
    {
      FP1 A1n = A1 - coordsc[3*p+0];
      FP1 B1n = B1 - coordsc[3*p+1];
      FP1 C1n = C1 - coordsc[3*p+2];
      FP2 r1n = sqrt(A1n*A1n+B1n*B1n+C1n*C1n);
      FP2 r2 = r1n*r1n;

      FP1 ne1 = Z1/r1n;
      FP1 dx = A1n*ne1/r2; FP1 dy = B1n*ne1/r2; FP1 dz = C1n*ne1/r2;

      V[p] += ne1;
      dV[3*p+0] += dx;
      dV[3*p+1] += dy;
      dV[3*p+2] += dz;
    }
  }

  if (prl>1)
  {
    printf("\n V: \n");
    for (int i=0;i<nc;i++)
      printf(" %10.5f",V[i]);
    printf("\n");
  }

  if (prl>1)
  {
    printf("\n dV: \n");
    for (int i=0;i<nc;i++)
    {
      for (int j=0;j<3;j++)
        printf(" %12.6f",dV[3*i+j]);
      printf("\n");
    }
    printf("\n");
  }

 //CPMZ check this
#if USE_ACC
  #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
  #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs]) 
  #pragma acc exit data delete(grid3[0:gs6],wt3[0:gs]) 
  #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6],grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6])
  #pragma acc exit data delete(wtt1[0:gs],wtt2[0:gs])
  #pragma acc exit data delete(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs],valS5[0:iN][0:gs],valS6[0:iN][0:gs])
  #pragma acc exit data delete(valt1[0:iN][0:gs],valt2[0:iN][0:gs],valt3[0:iN][0:gs],valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3],valtv3[0:iN][0:gs3])
  #pragma acc exit data delete(n2i[0:natoms])
  #pragma acc exit data delete(coords[0:3*natoms],atno[0:natoms])
#endif

  delete [] ang_g;
  delete [] ang_w;

  delete [] V1;
  delete [] dV1;

  delete [] n2i;
  delete [] norms;

  delete [] grid1s;
  delete [] grid2s;
  delete [] grid3s;
  delete [] grid1p;
  delete [] grid2p;
  delete [] grid3p;

  for (int i=0;i<iN;i++) delete [] valS1[i];
  for (int i=0;i<iN;i++) delete [] valS2[i];
  for (int i=0;i<iN;i++) delete [] valS3[i];
  for (int i=0;i<iN;i++) delete [] valS4[i];
  for (int i=0;i<iN;i++) delete [] valS5[i];
  for (int i=0;i<iN;i++) delete [] valS6[i];

  for (int i=0;i<iN;i++) delete [] valt1[i];
  for (int i=0;i<iN;i++) delete [] valt2[i];
  for (int i=0;i<iN;i++) delete [] valt3[i];
  for (int i=0;i<iN;i++) delete [] valtv1[i];
  for (int i=0;i<iN;i++) delete [] valtv2[i];
  for (int i=0;i<iN;i++) delete [] valtv3[i];

  delete [] valS1; delete [] valS2; delete [] valS3; delete [] valS4; delete [] valS5; delete [] valS6;
  delete [] valt1; delete [] valt2; delete [] valt3; delete [] valtv1; delete [] valtv2; delete [] valtv3;
  delete [] wtt1;
  delete [] wtt2;

  delete [] grid1;
  delete [] grid2;
  delete [] grid3;
  delete [] wt1;
  delete [] wt2;
  delete [] wt3;

  return;
}

#if RED_DOUBLE
void compute_Enp_para(int ngpu, int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* En, FP2* pVp, int prl)
#else
void compute_Enp_para(int ngpu, int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP1* En, FP1* pVp, int prl)
#endif
{
 #if !USE_ACC
  printf("\n ERROR: compute_Enp_para requires OpenACC \n");
  exit(1);
 #endif

  int nomp = ngpu;
 //#pragma omp parallel
  //nomp = omp_get_num_threads(); 

  if (prl>1) printf(" beginning compute_Enp_para (%i) \n",nomp);

  int N = basis.size();
  int N2 = N*N;

  int gs = nrad*nang;
  int gs3 = 3*gs; int gs6 = 6*gs;

  int estart = find_center_of_grid(1,nrad)*nang;

  FP1* grid1 = new FP1[gs6];
  FP1* wt1 = new FP1[gs];

  FP1* grid2 = new FP1[gs6];
  FP1* wt2 = new FP1[gs];

  FP1* grid3 = new FP1[gs6];
  FP1* wt3 = new FP1[gs];

  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);
  //printf("  iN: %i \n",imaxN);

 //intermediate storage
  int iN = imaxN;
  FP1* grid1s = new FP1[gs6]; FP1* grid2s = new FP1[gs6]; FP1* grid3s = new FP1[gs6];
  FP1* grid1p = new FP1[gs6]; FP1* grid2p = new FP1[gs6]; FP1* grid3p = new FP1[gs6];
  FP1** valS1 = new FP1*[iN]; for (int i=0;i<iN;i++) valS1[i] = new FP1[gs];
  FP1** valS2 = new FP1*[iN]; for (int i=0;i<iN;i++) valS2[i] = new FP1[gs];
  FP1** valS3 = new FP1*[iN]; for (int i=0;i<iN;i++) valS3[i] = new FP1[gs];
  FP1** valS4 = new FP1*[iN]; for (int i=0;i<iN;i++) valS4[i] = new FP1[gs];
  FP1** valS5 = new FP1*[iN]; for (int i=0;i<iN;i++) valS5[i] = new FP1[gs];
  FP1** valS6 = new FP1*[iN]; for (int i=0;i<iN;i++) valS6[i] = new FP1[gs];
  FP1** valV1 = new FP1*[iN]; for (int i=0;i<iN;i++) valV1[i] = new FP1[gs3];
  FP1** valV2 = new FP1*[iN]; for (int i=0;i<iN;i++) valV2[i] = new FP1[gs3];
  FP1** valV3 = new FP1*[iN]; for (int i=0;i<iN;i++) valV3[i] = new FP1[gs3];
  FP1** valV4 = new FP1*[iN]; for (int i=0;i<iN;i++) valV4[i] = new FP1[gs3];
  FP1** valV5 = new FP1*[iN]; for (int i=0;i<iN;i++) valV5[i] = new FP1[gs3];
  FP1** valV6 = new FP1*[iN]; for (int i=0;i<iN;i++) valV6[i] = new FP1[gs3];

  FP1** valt1 = new FP1*[iN]; for (int i=0;i<iN;i++) valt1[i] = new FP1[gs];
  FP1** valt2 = new FP1*[iN]; for (int i=0;i<iN;i++) valt2[i] = new FP1[gs];
  FP1** valt3 = new FP1*[iN]; for (int i=0;i<iN;i++) valt3[i] = new FP1[gs];
  FP1** valtv1 = new FP1*[iN]; for (int i=0;i<iN;i++) valtv1[i] = new FP1[gs3];
  FP1** valtv2 = new FP1*[iN]; for (int i=0;i<iN;i++) valtv2[i] = new FP1[gs3];
  FP1** valtv3 = new FP1*[iN]; for (int i=0;i<iN;i++) valtv3[i] = new FP1[gs3];

  FP1* wtt1 = new FP1[gs];
  FP1* wtt2 = new FP1[gs];

 #if RED_DOUBLE
  FP2* En1 = new FP2[N2];
  FP2* pVp1 = new FP2[N2]; 
 #else
  FP1* En1 = new FP1[N2];
  FP1* pVp1 = new FP1[N2]; 
 #endif

  FP1* ang_g = new FP1[3*nang];
  FP1* ang_w = new FP1[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

#if USE_ACC
 #pragma omp parallel for schedule(static,1) num_threads(nomp)
#endif
  for (int n=0;n<nomp;n++)
  {
    int tid = omp_get_thread_num();
    #if USE_ACC
    acc_set_device_num(tid,acc_device_nvidia);
    #endif
    #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
    #pragma acc enter data copyin(n2i[0:natoms])
    #pragma acc enter data copyin(coords[0:3*natoms],atno[0:natoms])

    #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
    #pragma acc enter data create(grid2[0:gs6],wt2[0:gs]) 
    #pragma acc enter data create(grid3[0:gs6],wt3[0:gs]) 
    #pragma acc enter data create(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs],valS5[0:iN][0:gs],valS6[0:iN][0:gs])
    #pragma acc enter data create(valV1[0:iN][0:gs3],valV2[0:iN][0:gs3],valV3[0:iN][0:gs3],valV4[0:iN][0:gs3],valV5[0:iN][0:gs3],valV6[0:iN][0:gs3])
    #pragma acc enter data create(valt1[0:iN][0:gs],valt2[0:iN][0:gs],valt3[0:iN][0:gs])
    #pragma acc enter data create(valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3],valtv3[0:iN][0:gs3])
    #pragma acc enter data create(grid1s[0:gs6],grid1p[0:gs6],grid2s[0:gs6],grid2p[0:gs6],grid3s[0:gs6],grid3p[0:gs6])
    #pragma acc enter data create(wtt1[0:gs],wtt2[0:gs])
    #pragma acc enter data create(En1[0:N2],pVp1[0:N2])
    if (tid>0)
    {
      #pragma acc enter data create(En[0:N2],pVp[0:N2])
    }
    acc_assign(N2,En,0.);
    acc_assign(N2,pVp,0.);
  }

#if USE_ACC
 #pragma omp parallel for num_threads(nomp)
#endif
  for (int m=0;m<natoms;m++)
  {
    int tid = omp_get_thread_num();
    #if USE_ACC
    acc_set_device_num(tid,acc_device_nvidia);
    #endif
    //printf("  launch %i/%i \n",m,tid);

   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    FP1 Z1 = (FP1)atno[m];
    FP1 A1 = coords[3*m+0]; FP1 B1 = coords[3*m+1]; FP1 C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);

  #if USE_ACC
   #pragma acc parallel loop present(valS1[0:iN][0:gs],valV1[0:iN][0:gs3])
  #endif
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valS1[ii1][j] = 1.f;
     #pragma acc loop
      for (int j=0;j<gs3;j++)
        valV1[ii1][j] = 1.f;
    }

  #if USE_ACC
   #pragma acc parallel loop present(valS3[0:iN][0:gs],valV3[0:iN][0:gs3])
  #endif
    for (int ii2=0;ii2<s2-s1;ii2++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valS3[ii2][j] = 1.f;
     #pragma acc loop
      for (int j=0;j<gs;j++)
     #pragma acc loop
      for (int k=0;k<3;k++)
        valV3[ii2][3*j+k] = 1.f;
    }

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<FP2> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

      eval_sh(ii1,gs,grid1,valS1[ii1],n1,l1,m1,zeta1);
      eval_p(gs,grid1,valV1[ii1],n1,l1,m1,zeta1);
    }

    for (int i2=s1;i2<s2;i2++)
    {
      int ii2 = i2-s1;
      vector<FP2> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

      eval_sh(ii2,gs,grid1,valS3[ii2],n2,l2,m2,zeta2);
      eval_p(gs,grid1,valV3[ii2],n2,l2,m2,zeta2);
    } //loop i2 evaluate

   #if USE_ACC
    #pragma acc wait
   #endif


  //first, do single center
   #pragma acc parallel loop present(wt1[0:gs],valt1[0:iN][0:gs],valtv1[0:iN][0:gs3],valS1[0:iN][0:gs],valV1[0:iN][0:gs3])
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valt1[ii1][j] = valS1[ii1][j]*wt1[j];
     #pragma acc loop collapse(2)
      for (int j=0;j<gs;j++)
      for (int k=0;k<3;k++)
        valtv1[ii1][3*j+k] = valV1[ii1][3*j+k]*wt1[j];
    }

   #pragma acc parallel loop present(grid1[0:gs6],valt1[0:iN][0:gs],valtv1[0:iN][0:gs3])
    for (int j=0;j<gs;j++)
    {
      FP1 Rn1 = grid1[6*j+3]+1.e-20f;
      FP1 ne1 = Z1/Rn1;

     #pragma acc loop
      for (int ii1=0;ii1<s2-s1;ii1++)
        valt1[ii1][j] *= ne1;

     #pragma acc loop collapse(2)
      for (int ii1=0;ii1<s2-s1;ii1++)
      for (int k=0;k<3;k++)
        valtv1[ii1][3*j+k] *= ne1;
    }

   //collect 1-atom values
    acc_assign(N2,En1,0.);
    acc_assign(N2,pVp1,0.);

    reduce_2c1(s1,s2,gs,valt1,valS3,iN,N,En1);
    reduce_2c1(s1,s2,gs3,valtv1,valV3,iN,N,pVp1);

   #pragma acc parallel loop present(En[0:N2],En1[0:N2],pVp[0:N2],pVp1[0:N2])
    for (int j=0;j<N2;j++)
    {
      En[j] += En1[j];
      pVp[j] += pVp1[j];
    }

   //2 basis on same center, then 1 nucleus elsewhere
    for (int p=0;p<natoms;p++)
    if (p!=m)
    {
      FP1 Zn = (FP1)atno[p];
      FP1 An = coords[3*p+0]; FP1 Bn = coords[3*p+1]; FP1 Cn = coords[3*p+2];
      FP1 A1n = An-A1; FP1 B1n = Bn-B1; FP1 C1n = Cn-C1;

      add_r2_to_grid(gs,grid1,A1n,B1n,C1n);

      generate_central_grid_2(grid3,wt3,Zn,nrad,nang,ang_g,ang_w);
      //copy_grid(gs,grid3p,grid3); //grid 3 centered on atom 3
      recenter_grid(gs,grid3,A1n,B1n,C1n); //grid 3 centered on atom 1

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid3,wt3,Z1,Zn,A1n,B1n,C1n);

      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wt3);

      add_r1_to_grid(gs,grid3,0.,0.,0.);
      add_r2_to_grid(gs,grid3,A1n,B1n,C1n);

      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc parallel loop present(valt1[0:iN][0:gs],valS1[0:iN][0:gs],wtt1[0:gs])
        for (int j=0;j<gs;j++)
          valt1[ii1][j] = valS1[ii1][j]*wtt1[j];

       #pragma acc parallel loop collapse(2) present(valtv1[0:iN][0:gs3],valV1[0:iN][0:gs3],wtt1[0:gs])
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valtv1[ii1][3*j+k] = valV1[ii1][3*j+k]*wtt1[j];

       #pragma acc parallel loop present(valS2[0:iN][0:gs],wt3[0:gs])
        for (int j=0;j<gs;j++)
          valS2[ii1][j] = wt3[j];

       #pragma acc parallel loop collapse(2) present(valV2[0:iN][0:gs3],wt3[0:gs])
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valV2[ii1][3*j+k] = wt3[j];
      }


      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

        eval_sh(ii1,gs,grid3,valS2[ii1],n1,l1,m1,zeta1);
        eval_p(gs,grid3,valV2[ii1],n1,l1,m1,zeta1);
      }

      for (int i2=s1;i2<s2;i2++)
      {
        int ii2 = i2-s1;
        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

        acc_assign(gs,valS4[ii2],1.);
        acc_assign(gs3,valV4[ii2],1.);

        eval_sh(ii2,gs,grid3,valS4[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid3,valV4[ii2],n2,l2,m2,zeta2);
      } //loop i2 evaluate

     #pragma acc parallel loop present(grid1[0:gs6],grid3[0:gs6],valt1[0:iN][0:gs],valS2[0:iN][0:gs],valtv1[0:iN][0:gs3],valV2[0:iN][0:gs3])
      for (int j=0;j<gs;j++)
      {
        FP1 Rn1 = grid1[6*j+4]+1.e-20f;
        FP1 Rn3 = grid3[6*j+4]+1.e-20f;
        FP1 ne1 = Zn/Rn1;
        FP1 ne3 = Zn/Rn3;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valt1[ii1][j] *= ne1;
       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valS2[ii1][j] *= ne3;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
       #pragma acc loop
        for (int k=0;k<3;k++)
          valtv1[ii1][3*j+k] *= ne1;
       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
       #pragma acc loop
        for (int k=0;k<3;k++)
          valV2[ii1][3*j+k] *= ne3;
      }

     //collect 2-atom values
      acc_assign(N2,En1,0.);
      acc_assign(N2,pVp1,0.);

      reduce_2c2(s1,s2,s1,s2,gs, valt1, valS2,valS3,valS4,iN,N,En1);
      reduce_2c2(s1,s2,s1,s2,gs3,valtv1,valV2,valV3,valV4,iN,N,pVp1);

     #pragma acc parallel loop present(En[0:N2],En1[0:N2],pVp[0:N2],pVp1[0:N2])
      for (int j=0;j<N2;j++)
      {
        En[j] += En1[j];
        pVp[j] += pVp1[j];
      }

    } //loop p over nuclear center




   //two-center basis
    for (int n=m+1;n<natoms;n++)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1); 
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      acc_copyf(gs,wtt1,wt1);
      acc_copyf(gs,wtt2,wt2);
      becke_weight_2c(gs,grid1,wtt1,grid2,wtt2,Z1,Z2,A12,B12,C12);

      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wtt2);

     //needs to happen after Becke weighting
      add_r1_to_grid(gs,grid2,0.,0.,0.);
      add_r2_to_grid(gs,grid2,A12,B12,C12);


    #if USE_ACC
     #pragma acc parallel loop present(valS1[0:iN][0:gs],valV1[0:iN][0:gs3],valS2[0:iN][0:gs],valV2[0:iN][0:gs3])
    #endif
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS1[ii1][j] = valS2[ii1][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          valV1[ii1][j] = valV2[ii1][j] = 1.f;
      }

    #if USE_ACC
     #pragma acc parallel loop present(valS3[0:iN][0:gs],valV3[0:iN][0:gs3],valS4[0:iN][0:gs],valV4[0:iN][0:gs3])
    #endif
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS3[ii2][j] = valS4[ii2][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          valV3[ii2][j] = valV4[ii2][j] = 1.f;
      }

      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

        eval_sh(ii1,gs,grid1,valS1[ii1],n1,l1,m1,zeta1);
        eval_sh(ii1,gs,grid2,valS2[ii1],n1,l1,m1,zeta1);
        eval_p(gs,grid1,valV1[ii1],n1,l1,m1,zeta1);
        eval_p(gs,grid2,valV2[ii1],n1,l1,m1,zeta1);
      }

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1s,valS3[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2s,valS4[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid1s,valV3[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid2s,valV4[ii2],n2,l2,m2,zeta2);
      }

     //2 basis on diff centers, 1 nucleus matches
     #pragma acc parallel loop present(wtt1[0:gs],valt1[0:iN][0:gs],valtv1[0:iN][0:gs3],valS1[0:iN][0:gs],valV1[0:iN][0:gs3])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valt1[ii1][j] = valS1[ii1][j]*wtt1[j];

       #pragma acc loop collapse(2)
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valtv1[ii1][3*j+k] = valV1[ii1][3*j+k]*wtt1[j];
      }

     #pragma acc parallel loop present(wtt2[0:gs],valt2[0:iN][0:gs],valtv2[0:iN][0:gs3],valS2[0:iN][0:gs],valV2[0:iN][0:gs3])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valt2[ii1][j] = valS2[ii1][j]*wtt2[j];

       #pragma acc loop collapse(2)
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valtv2[ii1][3*j+k] = valV2[ii1][3*j+k]*wtt2[j];
      }

     #pragma acc parallel loop present(grid1[0:gs6],grid2[0:gs6],valt1[0:iN][0:gs],valt2[0:iN][0:gs],valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3])
      for (int j=0;j<gs;j++)
      {
        FP1 Rn1a = grid1[6*j+3]+1.e-20f;
        FP1 Rn2a = grid1[6*j+4]+1.e-20f;
        FP1 Rn1b = grid2[6*j+3]+1.e-20f;
        FP1 Rn2b = grid2[6*j+4]+1.e-20f;
        FP1 ne1 = Z1/Rn1a+Z2/Rn2a;
        FP1 ne2 = Z1/Rn1b+Z2/Rn2b;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valt1[ii1][j] *= ne1;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valt2[ii1][j] *= ne2;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
       #pragma acc loop
        for (int k=0;k<3;k++)
          valtv1[ii1][3*j+k] *= ne1;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
       #pragma acc loop
        for (int k=0;k<3;k++)
          valtv2[ii1][3*j+k] *= ne2;
      }

     //collect 2c-atom same-atom values
      acc_assign(N2,En1,0.);
      acc_assign(N2,pVp1,0.);

      reduce_2c2(s1,s2,s3,s4,gs,valt1,valt2,valS3,valS4,iN,N,En1);
      reduce_2c2(s1,s2,s3,s4,gs3,valtv1,valtv2,valV3,valV4,iN,N,pVp1);

     #pragma acc parallel loop present(En[0:N2],En1[0:N2])
      for (int j=0;j<N2;j++)
        En[j] += En1[j];
     #pragma acc parallel loop present(pVp[0:N2],pVp1[0:N2])
      for (int j=0;j<N2;j++)
        pVp[j] += pVp1[j];


     //2 basis on diff centers, then 1 nucleus elsewhere
      for (int p=0;p<natoms;p++)
      if (p!=m && p!=n)
      {
        FP1 Zn = (FP1)atno[p];
        FP1 An = coords[3*p+0]; FP1 Bn = coords[3*p+1]; FP1 Cn = coords[3*p+2];
        FP1 A1n = An-A1; FP1 B1n = Bn-B1; FP1 C1n = Cn-C1;

        generate_central_grid_2(grid3,wt3,Zn,nrad,nang,ang_g,ang_w);
        //copy_grid(gs,grid3p,grid3); //grid 3 centered on atom 3
        recenter_grid(gs,grid3,A1n,B1n,C1n); //grid 3 centered on atom 1

        copy_grid(gs,grid3s,grid3);
        recenter_grid(gs,grid3s,-A12,-B12,-C12); //grid 3 centered on atom 2

        //copy_grid(gs,grid2p,grid2); //grid 2 centered on atom 3
        //recenter_grid(gs,grid2p,-A1n,-B1n,-C1n);
     
        //copy_grid(gs,grid1p,grid1);
        //recenter_grid(gs,grid1p,-A1n,-B1n,-C1n); //grid 1 centered on atom 3

        //add_r1_to_grid_6z(gs,grid1s,grid2s,grid3s,grid1p,grid2p,grid3p);
        add_r1_to_grid(gs,grid3s,0.,0.,0.);
      
      //need to keep all of these distances in order
        //add_r3_to_grid(gs,grid1,A1n,B1n,C1n);
        add_r123_to_grid(gs,grid1,0.,0.,0.,A12,B12,C12,A1n,B1n,C1n);
        add_r123_to_grid(gs,grid2,A12,B12,C12,0.,0.,0.,A1n,B1n,C1n);
      	add_r123_to_grid(gs,grid3,A1n,B1n,C1n,A12,B12,C12,0.,0.,0.);
     
        acc_copyf(gs,wtt1,wt1,wtt2,wt2);

        becke_weight_3c(gs,grid1,wtt1,grid2,wtt2,grid3,wt3,Z1,Z2,Zn,A12,B12,C12,A1n,B1n,C1n);
        eliminate_small_wt_3(estart,gs,wtt1,wtt2,wt3);

        add_r2_to_grid(gs,grid1,A1n,B1n,C1n);
        add_r1_to_grid(gs,grid2,0.,0.,0.);
        add_r1_to_grid(gs,grid3,0.,0.,0.);
        add_r2_to_grid(gs,grid2,A1n,B1n,C1n);
        add_r2_to_grid(gs,grid3,A1n,B1n,C1n);

        for (int ii1=0;ii1<s2-s1;ii1++)
        {
         #pragma acc parallel loop present(valt1[0:iN][0:gs],valS1[0:iN][0:gs],wtt1[0:gs])
          for (int j=0;j<gs;j++)
            valt1[ii1][j] = valS1[ii1][j]*wtt1[j];

         #pragma acc parallel loop present(valt2[0:iN][0:gs],valS2[0:iN][0:gs],wtt2[0:gs])
          for (int j=0;j<gs;j++)
            valt2[ii1][j] = valS2[ii1][j]*wtt2[j];

         #pragma acc parallel loop present(valS5[0:iN][0:gs],wt3[0:gs])
          for (int j=0;j<gs;j++)
            valS5[ii1][j] = wt3[j];

         #pragma acc parallel loop collapse(2) present(valtv1[0:iN][0:gs3],valV1[0:iN][0:gs3],wtt1[0:gs])
          for (int j=0;j<gs;j++)
          for (int k=0;k<3;k++)
            valtv1[ii1][3*j+k] = valV1[ii1][3*j+k]*wtt1[j];

         #pragma acc parallel loop collapse(2) present(valtv2[0:iN][0:gs3],valV2[0:iN][0:gs3],wtt2[0:gs])
          for (int j=0;j<gs;j++)
          for (int k=0;k<3;k++)
            valtv2[ii1][3*j+k] = valV2[ii1][3*j+k]*wtt2[j];

         #pragma acc parallel loop collapse(2) present(valV5[0:iN][0:gs3],wt3[0:gs])
          for (int j=0;j<gs;j++)
          for (int k=0;k<3;k++)
            valV5[ii1][3*j+k] = wt3[j];
        }

      #if USE_ACC
       #pragma acc parallel loop present(valS6[0:iN][0:gs],valV6[0:iN][0:gs3])
      #endif
        for (int ii2=0;ii2<s4-s3;ii2++)
        {
         #pragma acc loop
          for (int j=0;j<gs;j++)
            valS6[ii2][j] = 1.f;

         #pragma acc loop
          for (int j=0;j<gs3;j++)
            valV6[ii2][j] = 1.f;
        }

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<FP2> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

          eval_sh(ii1,gs,grid3,valS5[ii1],n1,l1,m1,zeta1);
          eval_p(gs,grid3,valV5[ii1],n1,l1,m1,zeta1);
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<FP2> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

          eval_sh(ii2,gs,grid3s,valS6[ii2],n2,l2,m2,zeta2);
          eval_p(gs,grid3s,valV6[ii2],n2,l2,m2,zeta2);
        } //loop i2 evaluate

       #pragma acc parallel loop present(grid1[0:gs6],grid2[0:gs6],grid3[0:gs6],valt1[0:iN][0:gs],valt2[0:iN][0:gs],valS5[0:iN][0:gs],valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3],valV5[0:iN][0:gs3])
        for (int j=0;j<gs;j++)
        {
          FP1 Rn1 = grid1[6*j+4]+1.e-20f; FP1 Rn2 = grid2[6*j+4]+1.e-20f; FP1 Rn3 = grid3[6*j+4]+1.e-20f;
          FP1 ne1 = Zn/Rn1; FP1 ne2 = Zn/Rn2; FP1 ne3 = Zn/Rn3;

         #pragma acc loop
          for (int ii1=0;ii1<s2-s1;ii1++)
          {
            valt1[ii1][j] *= ne1;
            valt2[ii1][j] *= ne2;
            valS5[ii1][j] *= ne3;
          }

         #pragma acc loop collapse(2)
          for (int ii1=0;ii1<s2-s1;ii1++)
          for (int k=0;k<3;k++)
          {
            valtv1[ii1][3*j+k] *= ne1;
            valtv2[ii1][3*j+k] *= ne2;
            valV5[ii1][3*j+k] *= ne3;
          }
        }

        //collect 3-atom values
        acc_assign(N2,En1,0.);
        acc_assign(N2,pVp1,0.);

        reduce_2c3(s1,s2,s3,s4,gs, valt1, valt2, valS3,valS4,valS5,valS6,iN,N,En1);
        reduce_2c3(s1,s2,s3,s4,gs3,valtv1,valtv2,valV3,valV4,valV5,valV6,iN,N,pVp1);

       #pragma acc parallel loop present(En[0:N2],En1[0:N2],pVp[0:N2],pVp1[0:N2])
        for (int j=0;j<N2;j++) 
        {
          En[j] += En1[j];
          pVp[j] += pVp1[j];
        }

      } //loop p over nuclear center

    } //loop n over second atom


  } //loop m over natoms


 //gather from all gpus
  FP2* En_all = new FP2[N2]();
  FP2* pVp_all = new FP2[N2]();

  for (int n=0;n<nomp;n++)
  {
    #if USE_ACC
    acc_set_device_num(n,acc_device_nvidia);
    #endif

    #pragma acc update self(En[0:N2],pVp[0:N2])

    for (int i=0;i<N2;i++)
      En_all[i] += En[i];
    for (int i=0;i<N2;i++)
      pVp_all[i] += pVp[i];

    #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
    #pragma acc exit data delete(En1[0:N2],pVp1[0:N2])
  }
  #if USE_ACC
  acc_set_device_num(0,acc_device_nvidia);
  #endif

  for (int i=0;i<N2;i++)
    En[i] = En_all[i];
  for (int i=0;i<N2;i++)
    pVp[i] = pVp_all[i];

  delete [] En_all;
  delete [] pVp_all;

  #pragma acc update device(En[0:N2],pVp[0:N2])

 //done gathering


  FP2* norm = new FP2[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];
  #pragma acc enter data copyin(norm[0:N])

 #pragma acc parallel loop independent present(En[0:N2],pVp[0:N2],norm[0:N])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i;j<N;j++)
  {
    FP2 n12 = norm[i]*norm[j];
    En[i*N+j]  *= -n12;
    pVp[i*N+j] *= -n12;
  }

 #pragma acc parallel loop independent present(En[0:N2],pVp[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i+1;j<N;j++)
  {
    En[j*N+i] = En[i*N+j];
    pVp[j*N+i] = pVp[i*N+j];
  }

  #pragma acc exit data delete(norm[0:N])
  delete [] norm;

  clean_small_values(N,En);
  clean_small_values(N,pVp);


  if (prl>1)
  {
    printf("\n En: \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %10.5f",En[i*N+j]);
      printf("\n");
    }
  }

  if (prl>2)
  {
    printf("\n pVp: \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %12.6f",pVp[i*N+j]);
      printf("\n");
    }
    printf("\n");
  }

 //CPMZ check this
#if USE_ACC
 #pragma omp parallel for schedule(static,1) num_threads(nomp)
#endif
  for (int n=0;n<nomp;n++)
  {
    int tid = omp_get_thread_num();
    #if USE_ACC
    acc_set_device_num(tid,acc_device_nvidia);
    #endif

    #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
    #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs]) 
    #pragma acc exit data delete(grid3[0:gs6],wt3[0:gs]) 
    #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6],grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6])
    #pragma acc exit data delete(wtt1[0:gs],wtt2[0:gs])
    #pragma acc exit data delete(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs],valS5[0:iN][0:gs],valS6[0:iN][0:gs])
    #pragma acc exit data delete(valV1[0:iN][0:gs3],valV2[0:iN][0:gs3],valV3[0:iN][0:gs3],valV4[0:iN][0:gs3],valV5[0:iN][0:gs3],valV6[0:iN][0:gs3])
    #pragma acc exit data delete(valt1[0:iN][0:gs],valt2[0:iN][0:gs],valt3[0:iN][0:gs],valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3],valtv3[0:iN][0:gs3])
    #pragma acc exit data delete(n2i[0:natoms])
    #pragma acc exit data delete(coords[0:3*natoms],atno[0:natoms])

    if (tid>0)
    {
      #pragma acc exit data delete(pVp[0:N2],En[0:N2])
    }
  }
  #if USE_ACC
  acc_set_device_num(0,acc_device_nvidia);
  #endif

  delete [] ang_g;
  delete [] ang_w;

  delete [] pVp1;
  delete [] En1;

  delete [] n2i;

  delete [] grid1s;
  delete [] grid2s;
  delete [] grid3s;
  delete [] grid1p;
  delete [] grid2p;
  delete [] grid3p;

  for (int i=0;i<iN;i++) delete [] valS1[i];
  for (int i=0;i<iN;i++) delete [] valS2[i];
  for (int i=0;i<iN;i++) delete [] valS3[i];
  for (int i=0;i<iN;i++) delete [] valS4[i];
  for (int i=0;i<iN;i++) delete [] valS5[i];
  for (int i=0;i<iN;i++) delete [] valS6[i];
  for (int i=0;i<iN;i++) delete [] valV1[i];
  for (int i=0;i<iN;i++) delete [] valV2[i];
  for (int i=0;i<iN;i++) delete [] valV3[i];
  for (int i=0;i<iN;i++) delete [] valV4[i];
  for (int i=0;i<iN;i++) delete [] valV5[i];
  for (int i=0;i<iN;i++) delete [] valV6[i];

  for (int i=0;i<iN;i++) delete [] valt1[i];
  for (int i=0;i<iN;i++) delete [] valt2[i];
  for (int i=0;i<iN;i++) delete [] valt3[i];
  for (int i=0;i<iN;i++) delete [] valtv1[i];
  for (int i=0;i<iN;i++) delete [] valtv2[i];
  for (int i=0;i<iN;i++) delete [] valtv3[i];

  delete [] valS1; delete [] valS2; delete [] valS3; delete [] valS4; delete [] valS5; delete [] valS6;
  delete [] valV1; delete [] valV2; delete [] valV3; delete [] valV4; delete [] valV5; delete [] valV6;
  delete [] valt1; delete [] valt2; delete [] valt3; delete [] valtv1; delete [] valtv2; delete [] valtv3;
  delete [] wtt1;
  delete [] wtt2;

  delete [] grid1;
  delete [] grid2;
  delete [] grid3;
  delete [] wt1;
  delete [] wt2;
  delete [] wt3;

  return;
}

#if RED_DOUBLE
void compute_Enp(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* En, FP2* pVp, int prl)
#else
void compute_Enp(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP1* En, FP1* pVp, int prl)
#endif
{
  if (prl>1) printf(" beginning compute_Enp \n");

  int N = basis.size();
  int N2 = N*N;

  int gs = nrad*nang;
  int gs3 = 3*gs; int gs6 = 6*gs;

  int estart = find_center_of_grid(1,nrad)*nang;

  FP1* grid1 = new FP1[gs6];
  FP1* wt1 = new FP1[gs];

  FP1* grid2 = new FP1[gs6];
  FP1* wt2 = new FP1[gs];

  FP1* grid3 = new FP1[gs6];
  FP1* wt3 = new FP1[gs];

  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);
  //printf("  iN: %i \n",imaxN);

 //intermediate storage
  int iN = imaxN;
  FP1* grid1s = new FP1[gs6]; FP1* grid2s = new FP1[gs6]; FP1* grid3s = new FP1[gs6];
  FP1* grid1p = new FP1[gs6]; FP1* grid2p = new FP1[gs6]; FP1* grid3p = new FP1[gs6];
  FP1** valS1 = new FP1*[iN]; for (int i=0;i<iN;i++) valS1[i] = new FP1[gs];
  FP1** valS2 = new FP1*[iN]; for (int i=0;i<iN;i++) valS2[i] = new FP1[gs];
  FP1** valS3 = new FP1*[iN]; for (int i=0;i<iN;i++) valS3[i] = new FP1[gs];
  FP1** valS4 = new FP1*[iN]; for (int i=0;i<iN;i++) valS4[i] = new FP1[gs];
  FP1** valS5 = new FP1*[iN]; for (int i=0;i<iN;i++) valS5[i] = new FP1[gs];
  FP1** valS6 = new FP1*[iN]; for (int i=0;i<iN;i++) valS6[i] = new FP1[gs];
  FP1** valV1 = new FP1*[iN]; for (int i=0;i<iN;i++) valV1[i] = new FP1[gs3];
  FP1** valV2 = new FP1*[iN]; for (int i=0;i<iN;i++) valV2[i] = new FP1[gs3];
  FP1** valV3 = new FP1*[iN]; for (int i=0;i<iN;i++) valV3[i] = new FP1[gs3];
  FP1** valV4 = new FP1*[iN]; for (int i=0;i<iN;i++) valV4[i] = new FP1[gs3];
  FP1** valV5 = new FP1*[iN]; for (int i=0;i<iN;i++) valV5[i] = new FP1[gs3];
  FP1** valV6 = new FP1*[iN]; for (int i=0;i<iN;i++) valV6[i] = new FP1[gs3];

  FP1** valt1 = new FP1*[iN]; for (int i=0;i<iN;i++) valt1[i] = new FP1[gs];
  FP1** valt2 = new FP1*[iN]; for (int i=0;i<iN;i++) valt2[i] = new FP1[gs];
  FP1** valt3 = new FP1*[iN]; for (int i=0;i<iN;i++) valt3[i] = new FP1[gs];
  FP1** valtv1 = new FP1*[iN]; for (int i=0;i<iN;i++) valtv1[i] = new FP1[gs3];
  FP1** valtv2 = new FP1*[iN]; for (int i=0;i<iN;i++) valtv2[i] = new FP1[gs3];
  FP1** valtv3 = new FP1*[iN]; for (int i=0;i<iN;i++) valtv3[i] = new FP1[gs3];

  FP1* wtt1 = new FP1[gs];
  FP1* wtt2 = new FP1[gs];

 #if RED_DOUBLE
  FP2* En1 = new FP2[N2];
  FP2* pVp1 = new FP2[N2]; 
 #else
  FP1* En1 = new FP1[N2];
  FP1* pVp1 = new FP1[N2]; 
 #endif

  FP1* ang_g = new FP1[3*nang];
  FP1* ang_w = new FP1[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

 #if USE_ACC
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data copyin(n2i[0:natoms])
  #pragma acc enter data copyin(coords[0:3*natoms],atno[0:natoms])

  #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
  #pragma acc enter data create(grid2[0:gs6],wt2[0:gs]) 
  #pragma acc enter data create(grid3[0:gs6],wt3[0:gs]) 
  #pragma acc enter data create(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs],valS5[0:iN][0:gs],valS6[0:iN][0:gs])
  #pragma acc enter data create(valV1[0:iN][0:gs3],valV2[0:iN][0:gs3],valV3[0:iN][0:gs3],valV4[0:iN][0:gs3],valV5[0:iN][0:gs3],valV6[0:iN][0:gs3])
  #pragma acc enter data create(valt1[0:iN][0:gs],valt2[0:iN][0:gs],valt3[0:iN][0:gs])
  #pragma acc enter data create(valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3],valtv3[0:iN][0:gs3])
  #pragma acc enter data create(grid1s[0:gs6],grid1p[0:gs6],grid2s[0:gs6],grid2p[0:gs6],grid3s[0:gs6],grid3p[0:gs6])
  #pragma acc enter data create(wtt1[0:gs],wtt2[0:gs])
  #pragma acc enter data create(En1[0:N2],pVp1[0:N2])
  //#pragma acc enter data create(En[0:N2],pVp[0:N2])
 #endif
  acc_assign(N2,En,0.);
  acc_assign(N2,pVp,0.);

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    FP1 Z1 = (FP1)atno[m];
    FP1 A1 = coords[3*m+0]; FP1 B1 = coords[3*m+1]; FP1 C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);

  #if USE_ACC
   #pragma acc parallel loop present(valS1[0:iN][0:gs],valV1[0:iN][0:gs3])
  #endif
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valS1[ii1][j] = 1.f;
     #pragma acc loop
      for (int j=0;j<gs3;j++)
        valV1[ii1][j] = 1.f;
    }

  #if USE_ACC
   #pragma acc parallel loop present(valS3[0:iN][0:gs],valV3[0:iN][0:gs3])
  #endif
    for (int ii2=0;ii2<s2-s1;ii2++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valS3[ii2][j] = 1.f;
     #pragma acc loop
      for (int j=0;j<gs;j++)
     #pragma acc loop
      for (int k=0;k<3;k++)
        valV3[ii2][3*j+k] = 1.f;
    }

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<FP2> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

      eval_sh(ii1,gs,grid1,valS1[ii1],n1,l1,m1,zeta1);
      eval_p(gs,grid1,valV1[ii1],n1,l1,m1,zeta1);
    }

    for (int i2=s1;i2<s2;i2++)
    {
      int ii2 = i2-s1;
      vector<FP2> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

      eval_sh(ii2,gs,grid1,valS3[ii2],n2,l2,m2,zeta2);
      eval_p(gs,grid1,valV3[ii2],n2,l2,m2,zeta2);
    } //loop i2 evaluate

   #if USE_ACC
    #pragma acc wait
   #endif


  //first, do single center
   #pragma acc parallel loop present(wt1[0:gs],valt1[0:iN][0:gs],valtv1[0:iN][0:gs3],valS1[0:iN][0:gs],valV1[0:iN][0:gs3])
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valt1[ii1][j] = valS1[ii1][j]*wt1[j];
     #pragma acc loop collapse(2)
      for (int j=0;j<gs;j++)
      for (int k=0;k<3;k++)
        valtv1[ii1][3*j+k] = valV1[ii1][3*j+k]*wt1[j];
    }

   #pragma acc parallel loop present(grid1[0:gs6],valt1[0:iN][0:gs],valtv1[0:iN][0:gs3])
    for (int j=0;j<gs;j++)
    {
      FP1 Rn1 = grid1[6*j+3]+1.e-20f;
      FP1 ne1 = Z1/Rn1;

     #pragma acc loop
      for (int ii1=0;ii1<s2-s1;ii1++)
        valt1[ii1][j] *= ne1;

     #pragma acc loop collapse(2)
      for (int ii1=0;ii1<s2-s1;ii1++)
      for (int k=0;k<3;k++)
        valtv1[ii1][3*j+k] *= ne1;
    }

   //collect 1-atom values
    acc_assign(N2,En1,0.);
    acc_assign(N2,pVp1,0.);

    reduce_2c1(s1,s2,gs,valt1,valS3,iN,N,En1);
    reduce_2c1(s1,s2,gs3,valtv1,valV3,iN,N,pVp1);

   #pragma acc parallel loop present(En[0:N2],En1[0:N2],pVp[0:N2],pVp1[0:N2])
    for (int j=0;j<N2;j++)
    {
      En[j] += En1[j];
      pVp[j] += pVp1[j];
    }

   //2 basis on same center, then 1 nucleus elsewhere
    for (int p=0;p<natoms;p++)
    if (p!=m)
    {
      FP1 Zn = (FP1)atno[p];
      FP1 An = coords[3*p+0]; FP1 Bn = coords[3*p+1]; FP1 Cn = coords[3*p+2];
      FP1 A1n = An-A1; FP1 B1n = Bn-B1; FP1 C1n = Cn-C1;

      add_r2_to_grid(gs,grid1,A1n,B1n,C1n);

      generate_central_grid_2(grid3,wt3,Zn,nrad,nang,ang_g,ang_w);
      //copy_grid(gs,grid3p,grid3); //grid 3 centered on atom 3
      recenter_grid(gs,grid3,A1n,B1n,C1n); //grid 3 centered on atom 1

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid3,wt3,Z1,Zn,A1n,B1n,C1n);

      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wt3);

      add_r1_to_grid(gs,grid3,0.,0.,0.);
      add_r2_to_grid(gs,grid3,A1n,B1n,C1n);

      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc parallel loop present(valt1[0:iN][0:gs],valS1[0:iN][0:gs],wtt1[0:gs])
        for (int j=0;j<gs;j++)
          valt1[ii1][j] = valS1[ii1][j]*wtt1[j];

       #pragma acc parallel loop collapse(2) present(valtv1[0:iN][0:gs3],valV1[0:iN][0:gs3],wtt1[0:gs])
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valtv1[ii1][3*j+k] = valV1[ii1][3*j+k]*wtt1[j];

       #pragma acc parallel loop present(valS2[0:iN][0:gs],wt3[0:gs])
        for (int j=0;j<gs;j++)
          valS2[ii1][j] = wt3[j];

       #pragma acc parallel loop collapse(2) present(valV2[0:iN][0:gs3],wt3[0:gs])
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valV2[ii1][3*j+k] = wt3[j];
      }


      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

        eval_sh(ii1,gs,grid3,valS2[ii1],n1,l1,m1,zeta1);
        eval_p(gs,grid3,valV2[ii1],n1,l1,m1,zeta1);
      }

      for (int i2=s1;i2<s2;i2++)
      {
        int ii2 = i2-s1;
        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

        acc_assign(gs,valS4[ii2],1.);
        acc_assign(gs3,valV4[ii2],1.);

        eval_sh(ii2,gs,grid3,valS4[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid3,valV4[ii2],n2,l2,m2,zeta2);
      } //loop i2 evaluate

     #pragma acc parallel loop present(grid1[0:gs6],grid3[0:gs6],valt1[0:iN][0:gs],valS2[0:iN][0:gs],valtv1[0:iN][0:gs3],valV2[0:iN][0:gs3])
      for (int j=0;j<gs;j++)
      {
        FP1 Rn1 = grid1[6*j+4]+1.e-20f;
        FP1 Rn3 = grid3[6*j+4]+1.e-20f;
        FP1 ne1 = Zn/Rn1;
        FP1 ne3 = Zn/Rn3;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valt1[ii1][j] *= ne1;
       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valS2[ii1][j] *= ne3;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
       #pragma acc loop
        for (int k=0;k<3;k++)
          valtv1[ii1][3*j+k] *= ne1;
       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
       #pragma acc loop
        for (int k=0;k<3;k++)
          valV2[ii1][3*j+k] *= ne3;
      }

     //collect 2-atom values
      acc_assign(N2,En1,0.);
      acc_assign(N2,pVp1,0.);

      reduce_2c2(s1,s2,s1,s2,gs, valt1, valS2,valS3,valS4,iN,N,En1);
      reduce_2c2(s1,s2,s1,s2,gs3,valtv1,valV2,valV3,valV4,iN,N,pVp1);

     #pragma acc parallel loop present(En[0:N2],En1[0:N2],pVp[0:N2],pVp1[0:N2])
      for (int j=0;j<N2;j++)
      {
        En[j] += En1[j];
        pVp[j] += pVp1[j];
      }

    } //loop p over nuclear center




   //two-center basis
    for (int n=m+1;n<natoms;n++)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1); 
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      acc_copyf(gs,wtt1,wt1);
      acc_copyf(gs,wtt2,wt2);
      becke_weight_2c(gs,grid1,wtt1,grid2,wtt2,Z1,Z2,A12,B12,C12);

      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wtt2);

     //needs to happen after Becke weighting
      add_r1_to_grid(gs,grid2,0.,0.,0.);
      add_r2_to_grid(gs,grid2,A12,B12,C12);


    #if USE_ACC
     #pragma acc parallel loop present(valS1[0:iN][0:gs],valV1[0:iN][0:gs3],valS2[0:iN][0:gs],valV2[0:iN][0:gs3])
    #endif
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS1[ii1][j] = valS2[ii1][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          valV1[ii1][j] = valV2[ii1][j] = 1.f;
      }

    #if USE_ACC
     #pragma acc parallel loop present(valS3[0:iN][0:gs],valV3[0:iN][0:gs3],valS4[0:iN][0:gs],valV4[0:iN][0:gs3])
    #endif
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS3[ii2][j] = valS4[ii2][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          valV3[ii2][j] = valV4[ii2][j] = 1.f;
      }

      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

        eval_sh(ii1,gs,grid1,valS1[ii1],n1,l1,m1,zeta1);
        eval_sh(ii1,gs,grid2,valS2[ii1],n1,l1,m1,zeta1);
        eval_p(gs,grid1,valV1[ii1],n1,l1,m1,zeta1);
        eval_p(gs,grid2,valV2[ii1],n1,l1,m1,zeta1);
      }

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1s,valS3[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2s,valS4[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid1s,valV3[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid2s,valV4[ii2],n2,l2,m2,zeta2);
      }

     //2 basis on diff centers, 1 nucleus matches
     #pragma acc parallel loop present(wtt1[0:gs],valt1[0:iN][0:gs],valtv1[0:iN][0:gs3],valS1[0:iN][0:gs],valV1[0:iN][0:gs3])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valt1[ii1][j] = valS1[ii1][j]*wtt1[j];

       #pragma acc loop collapse(2)
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valtv1[ii1][3*j+k] = valV1[ii1][3*j+k]*wtt1[j];
      }

     #pragma acc parallel loop present(wtt2[0:gs],valt2[0:iN][0:gs],valtv2[0:iN][0:gs3],valS2[0:iN][0:gs],valV2[0:iN][0:gs3])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valt2[ii1][j] = valS2[ii1][j]*wtt2[j];

       #pragma acc loop collapse(2)
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valtv2[ii1][3*j+k] = valV2[ii1][3*j+k]*wtt2[j];
      }

     #pragma acc parallel loop present(grid1[0:gs6],grid2[0:gs6],valt1[0:iN][0:gs],valt2[0:iN][0:gs],valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3])
      for (int j=0;j<gs;j++)
      {
        FP1 Rn1a = grid1[6*j+3]+1.e-20f;
        FP1 Rn2a = grid1[6*j+4]+1.e-20f;
        FP1 Rn1b = grid2[6*j+3]+1.e-20f;
        FP1 Rn2b = grid2[6*j+4]+1.e-20f;
        FP1 ne1 = Z1/Rn1a+Z2/Rn2a;
        FP1 ne2 = Z1/Rn1b+Z2/Rn2b;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valt1[ii1][j] *= ne1;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valt2[ii1][j] *= ne2;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
       #pragma acc loop
        for (int k=0;k<3;k++)
          valtv1[ii1][3*j+k] *= ne1;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
       #pragma acc loop
        for (int k=0;k<3;k++)
          valtv2[ii1][3*j+k] *= ne2;
      }

     //collect 2c-atom same-atom values
      acc_assign(N2,En1,0.);
      acc_assign(N2,pVp1,0.);

      reduce_2c2(s1,s2,s3,s4,gs,valt1,valt2,valS3,valS4,iN,N,En1);
      reduce_2c2(s1,s2,s3,s4,gs3,valtv1,valtv2,valV3,valV4,iN,N,pVp1);

     #pragma acc parallel loop present(En[0:N2],En1[0:N2])
      for (int j=0;j<N2;j++)
        En[j] += En1[j];
     #pragma acc parallel loop present(pVp[0:N2],pVp1[0:N2])
      for (int j=0;j<N2;j++)
        pVp[j] += pVp1[j];


     //2 basis on diff centers, then 1 nucleus elsewhere
      for (int p=0;p<natoms;p++)
      if (p!=m && p!=n)
      {
        FP1 Zn = (FP1)atno[p];
        FP1 An = coords[3*p+0]; FP1 Bn = coords[3*p+1]; FP1 Cn = coords[3*p+2];
        FP1 A1n = An-A1; FP1 B1n = Bn-B1; FP1 C1n = Cn-C1;

        generate_central_grid_2(grid3,wt3,Zn,nrad,nang,ang_g,ang_w);
        //copy_grid(gs,grid3p,grid3); //grid 3 centered on atom 3
        recenter_grid(gs,grid3,A1n,B1n,C1n); //grid 3 centered on atom 1

        copy_grid(gs,grid3s,grid3);
        recenter_grid(gs,grid3s,-A12,-B12,-C12); //grid 3 centered on atom 2

        //copy_grid(gs,grid2p,grid2); //grid 2 centered on atom 3
        //recenter_grid(gs,grid2p,-A1n,-B1n,-C1n);
     
        //copy_grid(gs,grid1p,grid1);
        //recenter_grid(gs,grid1p,-A1n,-B1n,-C1n); //grid 1 centered on atom 3

        //add_r1_to_grid_6z(gs,grid1s,grid2s,grid3s,grid1p,grid2p,grid3p);
        add_r1_to_grid(gs,grid3s,0.,0.,0.);
      
      //need to keep all of these distances in order
        //add_r3_to_grid(gs,grid1,A1n,B1n,C1n);
        add_r123_to_grid(gs,grid1,0.,0.,0.,A12,B12,C12,A1n,B1n,C1n);
        add_r123_to_grid(gs,grid2,A12,B12,C12,0.,0.,0.,A1n,B1n,C1n);
      	add_r123_to_grid(gs,grid3,A1n,B1n,C1n,A12,B12,C12,0.,0.,0.);
     
        acc_copyf(gs,wtt1,wt1,wtt2,wt2);

        becke_weight_3c(gs,grid1,wtt1,grid2,wtt2,grid3,wt3,Z1,Z2,Zn,A12,B12,C12,A1n,B1n,C1n);
        eliminate_small_wt_3(estart,gs,wtt1,wtt2,wt3);

        add_r2_to_grid(gs,grid1,A1n,B1n,C1n);
        add_r1_to_grid(gs,grid2,0.,0.,0.);
        add_r1_to_grid(gs,grid3,0.,0.,0.);
        add_r2_to_grid(gs,grid2,A1n,B1n,C1n);
        add_r2_to_grid(gs,grid3,A1n,B1n,C1n);

        for (int ii1=0;ii1<s2-s1;ii1++)
        {
         #pragma acc parallel loop present(valt1[0:iN][0:gs],valS1[0:iN][0:gs],wtt1[0:gs])
          for (int j=0;j<gs;j++)
            valt1[ii1][j] = valS1[ii1][j]*wtt1[j];

         #pragma acc parallel loop present(valt2[0:iN][0:gs],valS2[0:iN][0:gs],wtt2[0:gs])
          for (int j=0;j<gs;j++)
            valt2[ii1][j] = valS2[ii1][j]*wtt2[j];

         #pragma acc parallel loop present(valS5[0:iN][0:gs],wt3[0:gs])
          for (int j=0;j<gs;j++)
            valS5[ii1][j] = wt3[j];

         #pragma acc parallel loop collapse(2) present(valtv1[0:iN][0:gs3],valV1[0:iN][0:gs3],wtt1[0:gs])
          for (int j=0;j<gs;j++)
          for (int k=0;k<3;k++)
            valtv1[ii1][3*j+k] = valV1[ii1][3*j+k]*wtt1[j];

         #pragma acc parallel loop collapse(2) present(valtv2[0:iN][0:gs3],valV2[0:iN][0:gs3],wtt2[0:gs])
          for (int j=0;j<gs;j++)
          for (int k=0;k<3;k++)
            valtv2[ii1][3*j+k] = valV2[ii1][3*j+k]*wtt2[j];

         #pragma acc parallel loop collapse(2) present(valV5[0:iN][0:gs3],wt3[0:gs])
          for (int j=0;j<gs;j++)
          for (int k=0;k<3;k++)
            valV5[ii1][3*j+k] = wt3[j];
        }

      #if USE_ACC
       #pragma acc parallel loop present(valS6[0:iN][0:gs],valV6[0:iN][0:gs3])
      #endif
        for (int ii2=0;ii2<s4-s3;ii2++)
        {
         #pragma acc loop
          for (int j=0;j<gs;j++)
            valS6[ii2][j] = 1.f;

         #pragma acc loop
          for (int j=0;j<gs3;j++)
            valV6[ii2][j] = 1.f;
        }

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<FP2> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

          eval_sh(ii1,gs,grid3,valS5[ii1],n1,l1,m1,zeta1);
          eval_p(gs,grid3,valV5[ii1],n1,l1,m1,zeta1);
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<FP2> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

          eval_sh(ii2,gs,grid3s,valS6[ii2],n2,l2,m2,zeta2);
          eval_p(gs,grid3s,valV6[ii2],n2,l2,m2,zeta2);
        } //loop i2 evaluate

       #pragma acc parallel loop present(grid1[0:gs6],grid2[0:gs6],grid3[0:gs6],valt1[0:iN][0:gs],valt2[0:iN][0:gs],valS5[0:iN][0:gs],valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3],valV5[0:iN][0:gs3])
        for (int j=0;j<gs;j++)
        {
          FP1 Rn1 = grid1[6*j+4]+1.e-20f; FP1 Rn2 = grid2[6*j+4]+1.e-20f; FP1 Rn3 = grid3[6*j+4]+1.e-20f;
          FP1 ne1 = Zn/Rn1; FP1 ne2 = Zn/Rn2; FP1 ne3 = Zn/Rn3;

         #pragma acc loop
          for (int ii1=0;ii1<s2-s1;ii1++)
          {
            valt1[ii1][j] *= ne1;
            valt2[ii1][j] *= ne2;
            valS5[ii1][j] *= ne3;
          }

         #pragma acc loop collapse(2)
          for (int ii1=0;ii1<s2-s1;ii1++)
          for (int k=0;k<3;k++)
          {
            valtv1[ii1][3*j+k] *= ne1;
            valtv2[ii1][3*j+k] *= ne2;
            valV5[ii1][3*j+k] *= ne3;
          }
        }

        //collect 3-atom values
        acc_assign(N2,En1,0.);
        acc_assign(N2,pVp1,0.);

        reduce_2c3(s1,s2,s3,s4,gs, valt1, valt2, valS3,valS4,valS5,valS6,iN,N,En1);
        reduce_2c3(s1,s2,s3,s4,gs3,valtv1,valtv2,valV3,valV4,valV5,valV6,iN,N,pVp1);

       #pragma acc parallel loop present(En[0:N2],En1[0:N2],pVp[0:N2],pVp1[0:N2])
        for (int j=0;j<N2;j++) 
        {
          En[j] += En1[j];
          pVp[j] += pVp1[j];
        }

      } //loop p over nuclear center

    } //loop n over second atom


  } //loop m over natoms

  FP2* norm = new FP2[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];
  #pragma acc enter data copyin(norm[0:N])

 #pragma acc parallel loop independent present(En[0:N2],pVp[0:N2],norm[0:N])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i;j<N;j++)
  {
    FP2 n12 = norm[i]*norm[j];
    En[i*N+j]  *= -n12;
    pVp[i*N+j] *= -n12;
  }

 #pragma acc parallel loop independent present(En[0:N2],pVp[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i+1;j<N;j++)
  {
    En[j*N+i] = En[i*N+j];
    pVp[j*N+i] = pVp[i*N+j];
  }

  #pragma acc exit data delete(norm[0:N])
  delete [] norm;

  clean_small_values(N,En);
  clean_small_values(N,pVp);

 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc exit data delete(En1[0:N2],pVp1[0:N2])
  //#pragma acc exit data copyout(En[0:N2],pVp[0:N2])
  #pragma acc update self(En[0:N2],pVp[0:N2])
 #endif


  if (prl>1)
  {
    printf("\n En: \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %10.5f",En[i*N+j]);
      printf("\n");
    }
  }

  if (prl>2)
  {
    printf("\n pVp: \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %12.6f",pVp[i*N+j]);
      printf("\n");
    }
    printf("\n");
  }

 //CPMZ check this
#if USE_ACC
  #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
  #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs]) 
  #pragma acc exit data delete(grid3[0:gs6],wt3[0:gs]) 
  #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6],grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6])
  #pragma acc exit data delete(wtt1[0:gs],wtt2[0:gs])
  #pragma acc exit data delete(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs],valS5[0:iN][0:gs],valS6[0:iN][0:gs])
  #pragma acc exit data delete(valV1[0:iN][0:gs3],valV2[0:iN][0:gs3],valV3[0:iN][0:gs3],valV4[0:iN][0:gs3],valV5[0:iN][0:gs3],valV6[0:iN][0:gs3])
  #pragma acc exit data delete(valt1[0:iN][0:gs],valt2[0:iN][0:gs],valt3[0:iN][0:gs],valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3],valtv3[0:iN][0:gs3])
  #pragma acc exit data delete(n2i[0:natoms])
  #pragma acc exit data delete(coords[0:3*natoms],atno[0:natoms])
#endif

  delete [] ang_g;
  delete [] ang_w;

  delete [] pVp1;
  delete [] En1;

  delete [] n2i;

  delete [] grid1s;
  delete [] grid2s;
  delete [] grid3s;
  delete [] grid1p;
  delete [] grid2p;
  delete [] grid3p;

  for (int i=0;i<iN;i++) delete [] valS1[i];
  for (int i=0;i<iN;i++) delete [] valS2[i];
  for (int i=0;i<iN;i++) delete [] valS3[i];
  for (int i=0;i<iN;i++) delete [] valS4[i];
  for (int i=0;i<iN;i++) delete [] valS5[i];
  for (int i=0;i<iN;i++) delete [] valS6[i];
  for (int i=0;i<iN;i++) delete [] valV1[i];
  for (int i=0;i<iN;i++) delete [] valV2[i];
  for (int i=0;i<iN;i++) delete [] valV3[i];
  for (int i=0;i<iN;i++) delete [] valV4[i];
  for (int i=0;i<iN;i++) delete [] valV5[i];
  for (int i=0;i<iN;i++) delete [] valV6[i];

  for (int i=0;i<iN;i++) delete [] valt1[i];
  for (int i=0;i<iN;i++) delete [] valt2[i];
  for (int i=0;i<iN;i++) delete [] valt3[i];
  for (int i=0;i<iN;i++) delete [] valtv1[i];
  for (int i=0;i<iN;i++) delete [] valtv2[i];
  for (int i=0;i<iN;i++) delete [] valtv3[i];

  delete [] valS1; delete [] valS2; delete [] valS3; delete [] valS4; delete [] valS5; delete [] valS6;
  delete [] valV1; delete [] valV2; delete [] valV3; delete [] valV4; delete [] valV5; delete [] valV6;
  delete [] valt1; delete [] valt2; delete [] valt3; delete [] valtv1; delete [] valtv2; delete [] valtv3;
  delete [] wtt1;
  delete [] wtt2;

  delete [] grid1;
  delete [] grid2;
  delete [] grid3;
  delete [] wt1;
  delete [] wt2;
  delete [] wt3;

  return;
}

#if RED_DOUBLE
void compute_ST(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* S, FP2* T, int prl)
#else
void compute_ST(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP1* S, FP1* T, int prl)
#endif
{
  if (prl>1) printf(" beginning compute_ST \n");

  int N = basis.size();
  int N2 = N*N;

  int gs = nrad*nang;
  int gs6 = 6*gs;

  int estart = find_center_of_grid(1,nrad)*nang;

  FP1* grid1 = new FP1[gs6];
  FP1* wt1 = new FP1[gs];

  FP1* grid2 = new FP1[gs6];
  FP1* wt2 = new FP1[gs];

  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);
  //printf("  iN: %i \n",imaxN);

 //intermediate storage
  int iN = imaxN;
  FP1* grid1s = new FP1[gs6];
  FP1* grid2s = new FP1[gs6];
  FP1** valS1 = new FP1*[iN]; for (int i=0;i<iN;i++) valS1[i] = new FP1[gs];
  FP1** valS2 = new FP1*[iN]; for (int i=0;i<iN;i++) valS2[i] = new FP1[gs];
  FP1** valS3 = new FP1*[iN]; for (int i=0;i<iN;i++) valS3[i] = new FP1[gs];
  FP1** valS4 = new FP1*[iN]; for (int i=0;i<iN;i++) valS4[i] = new FP1[gs];
  FP1** valT1 = new FP1*[iN]; for (int i=0;i<iN;i++) valT1[i] = new FP1[gs];
  FP1** valT2 = new FP1*[iN]; for (int i=0;i<iN;i++) valT2[i] = new FP1[gs];
  FP1* wtt1 = new FP1[gs];

  FP1* ang_g = new FP1[3*nang];
  FP1* ang_w = new FP1[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

 #if USE_ACC
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data copyin(n2i[0:natoms])
  #pragma acc enter data copyin(coords[0:3*natoms],atno[0:natoms])

  #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
  #pragma acc enter data create(grid2[0:gs6],wt2[0:gs]) 
  #pragma acc enter data create(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs])
  #pragma acc enter data create(valT1[0:iN][0:gs],valT2[0:iN][0:gs])
  #pragma acc enter data create(grid1s[0:gs6],grid2s[0:gs6],wtt1[0:gs])
  //#pragma acc enter data create(S[0:N2],T[0:N2])
 #endif
  acc_assign(N2,S,0.);
  acc_assign(N2,T,0.);

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    FP1 Z1 = (FP1)atno[m];
    FP1 A1 = coords[3*m+0]; FP1 B1 = coords[3*m+1]; FP1 C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);

  #if USE_ACC
   #pragma acc parallel loop present(valS1[0:iN][0:gs])
  #endif
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valS1[ii1][j] = 1.f;
    }

  #if USE_ACC
   #pragma acc parallel loop present(valS3[0:iN][0:gs],wt1[0:gs])
  #endif
    for (int ii2=0;ii2<s2-s1;ii2++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valS3[ii2][j] = wt1[j];
    }

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<FP2> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

     //S
      eval_sh(ii1,gs,grid1,valS1[ii1],n1,l1,m1,zeta1);
    }

    for (int i2=s1;i2<s2;i2++)
    {
      int ii2 = i2-s1;
      vector<FP2> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

     //S
      eval_sh(ii2,gs,grid1,valS3[ii2],n2,l2,m2,zeta2);
    } //loop i2 evaluate

  #if USE_ACC
   #pragma acc parallel loop present(valS1[0:iN][0:gs],valT1[0:iN][0:gs])
  #endif
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valT1[ii1][j] = valS1[ii1][j];
    }
 
   //KE terms
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<FP2> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

      eval_ke(gs,grid1,valT1[ii1],n1,l1,zeta1);
    }

   #if USE_ACC
    #pragma acc wait
   #endif

    reduce_2c1(s1,s2,gs,valS1,valS3,iN,N,S);
    reduce_2c1(s1,s2,gs,valT1,valS3,iN,N,T);


   //two-atom ints
    for (int n=m+1;n<natoms;n++)
    //if (m!=n)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];
      //printf(" mn: %i %i s1-4: %i %i - %i %i \n",m,n,s1,s2,s3,s4);

      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      acc_copyf(gs,wtt1,wt1);
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      becke_weight_2c(gs,grid1,wtt1,grid2,wt2,Z1,Z2,A12,B12,C12);

    #if TEST_SORT
     #pragma acc host_data use_device(wtt1,grid1)
      sortByKeys(wtt1,(xyz<fp_t>*)grid1,gs);
     #pragma acc host_data use_device(wtt1,grid1)
      sortByKeys(wt2,(xyz<fp_t>*)grid2,gs);
    #else
      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wt2);
    #endif

      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2s,-A12,-B12,-C12);

      copy_grid(gs,grid1s,grid1); 
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

     //needs to happen after Becke weighting
      add_r1_to_grid(gs,grid2,0.,0.,0.);

    #if USE_ACC
     #pragma acc parallel loop present(valS1[0:iN][0:gs],valS2[0:iN][0:gs])
    #endif
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS1[ii1][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS2[ii1][j] = 1.f;
      }

    #if USE_ACC
     #pragma acc parallel loop present(valS3[0:iN][0:gs],valS4[0:iN][0:gs],wtt1[0:gs],wt2[0:gs])
    #endif
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS3[ii2][j] = wtt1[j];
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS4[ii2][j] = wt2[j];
      }

      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

       //S
        eval_sh(ii1,gs,grid1,valS1[ii1],n1,l1,m1,zeta1);
        eval_sh(ii1,gs,grid2,valS2[ii1],n1,l1,m1,zeta1);
      }

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

       //S
        eval_sh(ii2,gs,grid2s,valS4[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid1s,valS3[ii2],n2,l2,m2,zeta2);
      }

    #if USE_ACC
     #pragma acc parallel loop present(valS1[0:iN][0:gs],valT1[0:iN][0:gs])
    #endif
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valT1[ii1][j] = valS1[ii1][j];
      }
    #if USE_ACC
     #pragma acc parallel loop present(valS2[0:iN][0:gs],valT2[0:iN][0:gs])
    #endif
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valT2[ii1][j] = valS2[ii1][j];
      }

     //KE terms
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

        eval_ke(gs,grid1,valT1[ii1],n1,l1,zeta1);
        eval_ke(gs,grid2,valT2[ii1],n1,l1,zeta1);
      }

      reduce_2c2(s1,s2,s3,s4,gs,valS1,valS2,valS3,valS4,iN,N,S);
      reduce_2c2(s1,s2,s3,s4,gs,valT1,valT2,valS3,valS4,iN,N,T);

    } //loop n over second atom

  } //loop m over natoms


  FP2* norm = new FP2[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];
  #pragma acc enter data copyin(norm[0:N])

  #pragma acc parallel loop independent present(S[0:N2],T[0:N2],norm[0:N])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i;j<N;j++)
  {
    FP2 n12 = norm[i]*norm[j];
    S[i*N+j] *= n12;
    T[i*N+j] *= -0.5*n12;
  }

 #pragma acc parallel loop independent present(S[0:N2],T[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i;j<N;j++)
  {
    S[j*N+i] = S[i*N+j];
    T[j*N+i] = T[i*N+j];
  }

  #pragma acc exit data delete(norm[0:N])
  delete [] norm;

  clean_small_values(N,S);
  clean_small_values(N,T);


 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  //#pragma acc exit data copyout(S[0:N2],T[0:N2])
  #pragma acc update self(S[0:N2],T[0:N2])
 #endif


  if (prl>1 || prl==-1)
  {
    printf("\n S: \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %12.6f",S[i*N+j]);
      printf("\n");
    }
  }

  if (prl>1)
  {
    printf("\n T: \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %12.6f",T[i*N+j]);
      printf("\n");
    }
  }

 //CPMZ check this
#if USE_ACC
  #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
  #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs]) 
  #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6],wtt1[0:gs])
  #pragma acc exit data delete(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs])
  #pragma acc exit data delete(valT1[0:iN][0:gs],valT2[0:iN][0:gs])
  #pragma acc exit data delete(n2i[0:natoms])
  #pragma acc exit data delete(coords[0:3*natoms],atno[0:natoms])
#endif

  delete [] ang_g;
  delete [] ang_w;

  delete [] n2i;

  delete [] grid1s;
  delete [] grid2s;

  for (int i=0;i<iN;i++) delete [] valS1[i];
  for (int i=0;i<iN;i++) delete [] valS2[i];
  for (int i=0;i<iN;i++) delete [] valS3[i];
  for (int i=0;i<iN;i++) delete [] valS4[i];
  for (int i=0;i<iN;i++) delete [] valT1[i];
  for (int i=0;i<iN;i++) delete [] valT2[i];
  delete [] valS1; delete [] valS2; delete [] valS3; delete [] valS4;
  delete [] valT1; delete [] valT2;
  delete [] wtt1;

  delete [] grid1;
  delete [] grid2;
  delete [] wt1;
  delete [] wt2;

  return;
}

#if RED_DOUBLE
void compute_all_2c_v2(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* An, int prl)
#else
void compute_all_2c_v2(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP1* An, int prl)
#endif
{
  if (prl>1) printf(" beginning compute_all_2c_v2 \n");

 //2c integrals are all in auxiliary basis
  int N = basis.size();
  int N2 = N*N;

  int gs = nrad*nang;
  int gs6 = 6*gs;

  int estart = find_center_of_grid(1,nrad)*nang;

  FP1* grid1 = new FP1[gs6];
  FP1* wt1 = new FP1[gs];

  FP1* grid2 = new FP1[gs6];
  FP1* wt2 = new FP1[gs];

  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);
  //printf("  iN: %i \n",imaxN);

 //intermediate storage
  int iN = imaxN;
  FP1* grid1s = new FP1[gs6];
  FP1* grid2s = new FP1[gs6];
  FP1** val1 = new FP1*[iN];
  FP1** val2 = new FP1*[iN];
  FP1** val3 = new FP1*[iN];
  FP1** val4 = new FP1*[iN];
  for (int i=0;i<iN;i++)
    val1[i] = new FP1[gs];
  for (int i=0;i<iN;i++)
    val2[i] = new FP1[gs];
  for (int i=0;i<iN;i++)
    val3[i] = new FP1[gs];
  for (int i=0;i<iN;i++)
    val4[i] = new FP1[gs];
  FP1* wtt1 = new FP1[gs];

  FP1* ang_g = new FP1[3*nang];
  FP1* ang_w = new FP1[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

 #if USE_ACC
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data copyin(n2i[0:natoms])

  #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
  #pragma acc enter data create(grid2[0:gs6],wt2[0:gs]) 
  #pragma acc enter data create(val1[0:iN][0:gs],val2[0:iN][0:gs],wtt1[0:gs])
  #pragma acc enter data create(val3[0:iN][0:gs],val4[0:iN][0:gs])
  #pragma acc enter data create(grid1s[0:gs6],grid2s[0:gs6])
  //#pragma acc enter data create(An[0:N2])
 #endif
  acc_assign(N2,An,0.);

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    FP1 Z1 = (FP1)atno[m];
    FP1 A1 = coords[3*m+0]; FP1 B1 = coords[3*m+1]; FP1 C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);

  #if USE_ACC
   #pragma acc parallel loop present(val1[0:iN][0:gs])
  #endif
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val1[ii1][j] = 1.f;
    }

  #if USE_ACC
   #pragma acc parallel loop present(val3[0:iN][0:gs],wt1[0:gs])
  #endif
    for (int ii2=0;ii2<s2-s1;ii2++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val3[ii2][j] = wt1[j];
    }

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<FP2> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

      //acc_assign(gs,val1[ii1],1.);

      eval_inr_r12(gs,grid1,val1[ii1],n1,l1,zeta1);
      eval_sh_3r(gs,grid1,val1[ii1],n1,l1,m1);
    }

    for (int i2=s1;i2<s2;i2++)
    {
      int ii2 = i2-s1;
      vector<FP2> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

      eval_sh(ii2,gs,grid1,val3[ii2],n2,l2,m2,zeta2);

    } //loop i2 evaluate
   #if USE_ACC
    #pragma acc wait
   #endif

    reduce_2c1(s1,s2,gs,val1,val3,iN,N,An);

   //two-atom ints
    for (int n=m+1;n<natoms;n++)
    //if (m!=n)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];
      //printf(" mn: %i %i s1-4: %i %i - %i %i \n",m,n,s1,s2,s3,s4);

      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1); 
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid2,wt2,Z1,Z2,A12,B12,C12);

      //#pragma acc update self(grid1[0:gs6],wtt1[0:gs],grid2[0:gs6],wt2[0:gs])
      //print_grid(gs,grid1,grid2,wtt1,wt2,prl);

      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wt2);

     //needs to happen after Becke weighting
      add_r1_to_grid(gs,grid2,0.,0.,0.);

    #if USE_ACC
     #pragma acc parallel loop present(val1[0:iN][0:gs],val2[0:iN][0:gs])
    #endif
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val1[ii1][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val2[ii1][j] = 1.f;
      }

    #if USE_ACC
     #pragma acc parallel loop present(val3[0:iN][0:gs],val4[0:iN][0:gs],wtt1[0:gs],wt2[0:gs])
    #endif
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val3[ii2][j] = wtt1[j];
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val4[ii2][j] = wt2[j];
      }

      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];
        //printf("  m: %i i1: %i   nlm: %i %i %2i zeta: %8.5f \n",m,i1,n1,l1,m1,zeta1);

        //acc_assign(gs,val1[ii1],1.);
        //acc_assign(gs,val2[ii1],1.);

        eval_inr_r12(gs,grid1,val1[ii1],n1,l1,zeta1);
        eval_sh_3r(gs,grid1,val1[ii1],n1,l1,m1);
        eval_inr_r12(gs,grid2,val2[ii1],n1,l1,zeta1);
        eval_sh_3r(gs,grid2,val2[ii1],n1,l1,m1);
      }

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];
        //printf("   n: %i i2: %i   nlm: %i %i %2i zeta: %8.5f \n",n,i2,n2,l2,m2,zeta2);

        eval_sh(ii2,gs,grid2s,val4[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid1s,val3[ii2],n2,l2,m2,zeta2);
      }

      reduce_2c2(s1,s2,s3,s4,gs,val1,val2,val3,val4,iN,N,An);

    } //loop n over second atom

  } //loop m over natoms

  FP2* norm1 = new FP2[N];
  FP2* norm2 = new FP2[N];
  for (int i=0;i<N;i++)
  {
    norm1[i] = norm_sv(basis[i][0],basis[i][1],basis[i][2],basis[i][3]);
    norm2[i] = basis[i][4];
  }
  #pragma acc enter data copyin(norm1[0:N],norm2[0:N])

 #pragma acc parallel loop independent present(An[0:N2],norm1[0:N],norm2[0:N])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i;j<N;j++)
  {
    FP2 n12 = norm1[i]*norm2[j];
    An[i*N+j] *= n12;
  }

  #pragma acc exit data delete(norm1[0:N],norm2[0:N])
  delete [] norm1;
  delete [] norm2;

 #pragma acc parallel loop independent present(An[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i+1;j<N;j++)
    An[j*N+i] = An[i*N+j];

  //clean_small_values(N,An);

 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  //#pragma acc exit data copyout(An[0:N2])
  #pragma acc update self(An[0:N2])
 #endif


 #if 0
  for (int m=0;m<N;m++)
  {
    FP2 val = An[m*N+m];
    if (fabs(val)<1.e-2)
    {
      printf(" WARNING: small diagonal element in A: %3i \n",m);
    }
    FP2 v2 = 1./sqrt(val);
    //Anorm[m] = v2;
    for (int n=0;n<N;n++)
    {
      An[m*N+n] *= v2;
      An[n*N+m] *= v2;
    }
  }
 #endif


  if (prl>2)
  {
    printf("\n A: \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %12.6f",An[i*N+j]);
      printf("\n");
    }
  }

 //CPMZ check this
#if USE_ACC
  #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
  #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs]) 
  #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6])
  #pragma acc exit data delete(val1[0:iN][0:gs],val2[0:iN][0:gs],wtt1[0:gs])
  #pragma acc exit data delete(val3[0:iN][0:gs],val4[0:iN][0:gs])
  #pragma acc exit data delete(n2i[0:natoms])
#endif

  delete [] ang_g;
  delete [] ang_w;

  delete [] n2i;

  delete [] grid1s;
  delete [] grid2s;

  for (int i=0;i<iN;i++) delete [] val1[i];
  for (int i=0;i<iN;i++) delete [] val2[i];
  for (int i=0;i<iN;i++) delete [] val3[i];
  for (int i=0;i<iN;i++) delete [] val4[i];
  delete [] val1;
  delete [] val2;
  delete [] val3;
  delete [] val4;
  delete [] wtt1;

  delete [] grid1;
  delete [] grid2;
  delete [] wt1;
  delete [] wt2;

  return;
}

void compute_all_2c(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP1* A, int prl)
{
  if (prl>1) printf(" beginning compute_all_2c \n");

 //2c integrals are all in auxiliary basis
  int N = basis.size();
  int N2 = N*N;

  int gs = nrad*nang;
  FP1* grid1 = new FP1[6*gs];
  FP1* wt1 = new FP1[gs];
  FP1* val1 = new FP1[gs];

  FP1* grid2 = new FP1[6*gs];
  FP1* wt2 = new FP1[gs];
  FP1* val2 = new FP1[gs];

  int* i2m = new int[N];

 #if 1
  FP1* grid1s = new FP1[6*gs];
  FP1* grid2s = new FP1[6*gs];
  FP1** valt1 = new FP1*[N];
  FP1** valt2 = new FP1*[N];
  for (int i=0;i<N;i++)
    valt1[i] = new FP1[gs];
  for (int i=0;i<N;i++)
    valt2[i] = new FP1[gs];
 #else
  FP1* valt1 = new FP1[gs];
  FP1* valt2 = new FP1[gs];
 #endif
  FP1* wtt1 = new FP1[gs];

  FP1* ang_g = new FP1[3*nang];
  FP1* ang_w = new FP1[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

 #if USE_ACC
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data create(grid1[0:6*gs],wt1[0:gs],val1[0:gs])
  #pragma acc enter data create(grid2[0:6*gs],wt2[0:gs],val2[0:gs]) 
  #pragma acc enter data create(valt1[0:N][0:gs],valt2[0:N][0:gs],wtt1[0:gs])
  #pragma acc enter data create(grid1s[0:6*gs],grid2s[0:6*gs])
  //#pragma acc enter data create(valt1[0:gs],valt2[0:gs],wtt1[0:gs])
  #pragma acc enter data create(A[0:N2]) 
  #pragma acc enter data create(i2m[0:N])
 #endif

  for (int m=0;m<natoms;m++)
  {
    FP1 Z1 = (FP1)atno[m];
    FP1 A1 = coords[3*m+0]; FP1 B1 = coords[3*m+1]; FP1 C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);
    //#pragma acc update self(grid1[0:6*gs],wt1[0:gs])
    //print_grid(gs,grid1,NULL,wt1,NULL,prl);

    for (int i1=0;i1<N;i1++) i2m[i1] = 0;
    for (int i1=0;i1<N;i1++) if (basis[i1][9]==m) i2m[i1] = 1;
   #if USE_ACC
    #pragma acc update device(i2m[0:N])
   #endif

   //first compute single atom ints
    for (int i1=0;i1<N;i1++)
    if (basis[i1][9]==m)
    {
      vector<FP2> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];
      //printf("\n  m: %i i1: %2i   nlm: %i %i %2i zeta: %8.5f \n",m,i1,n1,l1,m1,zeta1);

      acc_assign(gs,val1,1.);

      eval_inr_r12(gs,grid1,val1,n1,l1,zeta1);
      eval_sh_3r(gs,grid1,val1,n1,l1,m1);

      //#pragma acc update self(val1)
      //printf(" val1[%i] \n",i1);
      //print_array(gs,val1);

    #if USE_ACC
     #pragma acc parallel loop present(val1[0:gs],valt1[0:N][0:gs],i2m[0:N])
    #endif
      for (int i2=0;i2<N;i2++)
      if (i2m[i2])
      {
      #if USE_ACC
       #pragma acc loop
      #endif
        for (int j=0;j<gs;j++)
          valt1[i2][j] = val1[j];
      }

      for (int i2=0;i2<N;i2++)
      if (basis[i2][9]==m)
      {
        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];
        //printf("    n: %i i2: %2i   nlm: %i %i %2i zeta: %8.5f   (b2) \n",m,i2,n2,l2,m2,zeta2);

        //acc_copyf(i2,gs,valt1[i2],val1);

        eval_sh(i2,gs,grid1,valt1[i2],n2,l2,m2,zeta2);

        //#pragma acc update self(valt1[i2])
        //printf(" val2[%i] \n",i2);
        //print_array(gs,valt1[i2]);

      } //loop i2 evaluate
     #if USE_ACC
      #pragma acc wait
     #endif

    #if USE_ACC
     #pragma acc parallel loop present(valt1[0:N][0:gs],wt1[0:gs],A[0:N2],i2m[0:N])
    #endif
      for (int i2=0;i2<N;i2++)
      if (i2m[i2])
      {
        FP1 val = 0.;
      #if USE_ACC
       #pragma acc loop reduction(+:val)
      #endif
        for (int j=0;j<gs;j++)
          val += valt1[i2][j] * wt1[j];
        A[i1*N+i2] = val;

      } //loop i2 reduce

    } //loop i1


   //two-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1); 
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid2,wt2,Z1,Z2,A12,B12,C12);
      eliminate_small_wt(gs,wtt1);
      eliminate_small_wt(gs,wt2);

      add_r1_to_grid(gs,grid2,0.,0.,0.);

      //#pragma acc update self(grid1[0:6*gs],wtt1[0:gs],grid2[0:6*gs],wt2[0:gs])
      //print_grid(gs,grid1,grid2,wtt1,wt2,prl);

      for (int i1=0;i1<N;i1++) i2m[i1] = 0;
      for (int i1=0;i1<N;i1++) if (basis[i1][9]==n) i2m[i1] = 1;
     #if USE_ACC
      #pragma acc update device(i2m[0:N])
     #endif

      for (int i1=0;i1<N;i1++)
      if (basis[i1][9]==m)
      {
        vector<FP2> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];
        //printf("  m: %i i1: %i   nlm: %i %i %2i zeta: %8.5f \n",m,i1,n1,l1,m1,zeta1);

        acc_assign(gs,val1,1.);
        acc_assign(gs,val2,1.);

        eval_inr_r12(gs,grid1,val1,n1,l1,zeta1);
        eval_sh_3r(gs,grid1,val1,n1,l1,m1);
        eval_inr_r12(gs,grid2,val2,n1,l1,zeta1);
        eval_sh_3r(gs,grid2,val2,n1,l1,m1);

      #if USE_ACC
       #pragma acc parallel loop present(val1[0:gs],valt1[0:N][0:gs],val2[0:gs],valt2[0:N][0:gs],i2m[0:N])
      #endif
        for (int i2=0;i2<N;i2++)
        if (i2m[i2])
        {
        #if USE_ACC
         #pragma acc loop
        #endif
          for (int j=0;j<gs;j++)
            valt1[i2][j] = val1[j];
        #if USE_ACC
         #pragma acc loop
        #endif
          for (int j=0;j<gs;j++)
            valt2[i2][j] = val2[j];
        }

       //launch async processes over i2
        for (int i2=i1+1;i2<N;i2++)
        if (basis[i2][9]==n)
        {
          vector<FP2> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];
          //printf("   n: %i i2: %i   nlm: %i %i %2i zeta: %8.5f \n",n,i2,n2,l2,m2,zeta2);

          //acc_copyf(i2,gs,valt2[i2],val2);
          //acc_copyf(i2,gs,valt1[i2],val1);

          eval_sh(i2,gs,grid2s,valt2[i2],n2,l2,m2,zeta2);
          eval_sh(i2,gs,grid1s,valt1[i2],n2,l2,m2,zeta2);
        }
       #if USE_ACC
        #pragma acc wait
       #endif

      #if USE_ACC
       #pragma acc parallel loop present(valt1[0:N][0:gs],wtt1[0:gs],valt2[0:N][0:gs],wt2[0:gs],A[0:N2],i2m[0:N])
      #endif
        for (int i2=i1+1;i2<N;i2++)
        if (i2m[i2])
        {
          FP1 val = 0.;
        #if USE_ACC
         #pragma acc loop reduction(+:val)
        #endif
          for (int j=0;j<gs;j++)
            val += valt1[i2][j] * wtt1[j];
        #if USE_ACC
         #pragma acc loop reduction(+:val)
        #endif
          for (int j=0;j<gs;j++)
            val += valt2[i2][j] * wt2[j];

          A[i1*N+i2] = val;

        } //loop i2 reduce

      } //loop i1 over first basis function

    } //loop n over second atom

  } //loop m over natoms

 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc exit data copyout(A[0:N2])
 #endif

#if DEBUG
  if (prl>0)
  {
    printf("\n A(raw): \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %8.5f",A[i*N+j]);
      printf("\n");
    }
  }
#endif

  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    A[i*N+j] *= norm_sv(basis[i][0],basis[i][1],basis[i][2],basis[i][3])*basis[j][4];

  for (int i=0;i<N;i++)
  for (int j=i+1;j<N;j++)
    A[j*N+i] = A[i*N+j];

  if (prl>2)
  {
    printf("\n A: \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %12.6f",A[i*N+j]);
      printf("\n");
    }
  }

#if USE_ACC
  #pragma acc exit data delete(grid1[0:6*gs],wt1[0:gs],val1[0:gs])
  #pragma acc exit data delete(grid2[0:6*gs],wt2[0:gs],val2[0:gs]) 
  #pragma acc exit data delete(grid1s[0:6*gs],grid2s[0:6*gs])
  #pragma acc exit data delete(valt1[0:N][0:gs],valt2[0:N][0:gs],wtt1[0:gs])
  #pragma acc exit data delete(A[0:N2])
  #pragma acc exit data delete(i2m[0:N])
#endif

  delete [] ang_g;
  delete [] ang_w;

  delete [] i2m;

 #if 1
  for (int i=0;i<N;i++) delete [] valt1[i];
  for (int i=0;i<N;i++) delete [] valt2[i];
 #endif
  delete [] valt1;
  delete [] valt2;
  delete [] wtt1;

  delete [] grid1;
  delete [] grid2;
  delete [] wt1;
  delete [] wt2;
  delete [] val1;
  delete [] val2;

  return;
}

