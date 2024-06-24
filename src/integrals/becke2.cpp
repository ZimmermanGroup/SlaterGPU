#include "becke.h"
#include "braggslater.h"
#include "murak.h"
#include "integrals.h"
#include "read.h"
//#include "vxc.h"
//#include "gauss.h"
//#include <xc.h>

#define pr_inc 5

void generate_central_grid_2(float* grid1, float* wt1, float Z1, int nrad, int nang, float* ang_g, float* ang_w);
void copy_grid(int gs, float* grid1, float* grid2);
void recenter_grid(int gs, float* grid, float x2, float y2, float z2);
void recenter_grid_zero(int gs, float* grid, float x2, float y2, float z2);
void add_r1_to_grid(int gs, float* grid1, float A2, float B2, float C2);
void add_r2_to_grid(int gs, float* grid1, float A2, float B2, float C2);
void print_square(int N, double* A);
//void print_vxc(int nrad, int nang, int natoms, float* grid, double* vxc);
void print_vec(int gsa, float* grid, float* A);
void print_vec(int gsa, float* grid, float* A, float* B);

using namespace std;

//#pragma acc routine seq
int get_atom_grid(int n, int natoms, int gsa)
{
  int f1 = gsa/natoms;
  int f2 = n%f1;
  int wa = (n-f2)/f1;

  //printf("  n: %4i wa: %2i \n",n,wa);

  return wa;
}

float sigmoidf(float v1)
{
  float val = 1.f/(1.f+expf(-v1));
  return 1.f - val;
  //return val;
}

float tanhh(float a)
{
  return (tanh(a)+1.)*0.5;
}


//n-center grid wts, stretched boundaries
// ta>-1 --> just get the partition for atom ta
// alpha<0. gives rescaled aij, alpha>0. gives sigmoid dropoff

void atomic_domain_cell_wt(const int ta, const float alpha, const float beta, const int natoms, const int gc, const int gs, float* grid1, float* wt1, int* atno, float* coords0)
{
  int gsa = gs*natoms;
  if (natoms==1)
  {
    if (ta>=0)
    for (int j=0;j<gsa;j++)
      wt1[j] = 1.;
    #pragma acc update device(wt1[0:gsa])

    return;
  }

  //printf("  atomic domain cells (alpha: %6.3f beta: %6.3f ta: %i) \n",alpha,beta,ta);
  //for (int n=0;n<natoms;n++)
  //  printf("  a(%i): %8.5f \n",n+1,get_radii_2(atno[n]));

  if (natoms<2) return;
  if (ta<0) { printf("\n ERROR: atomic weights with ta<0 not available \n"); exit(-1); }

  float cx0 = 0.f; float cy0 = 0.f; float cz0 = 0.f;
  if (ta>=0)
  {
    cx0 = coords0[3*ta+0];
    cy0 = coords0[3*ta+1];
    cz0 = coords0[3*ta+2];
  }
  float coords[3*natoms];
  for (int j=0;j<natoms;j++)
  {
    coords[3*j+0] = coords0[3*j+0] - cx0;
    coords[3*j+1] = coords0[3*j+1] - cy0;
    coords[3*j+2] = coords0[3*j+2] - cz0;
  }

 #if 0
  printf("  coords (shifted): \n");
  for (int n=0;n<natoms;n++)
    printf("   %5.3f %5.3f %5.3f \n",coords[3*n+0],coords[3*n+1],coords[3*n+2]);

  printf("  first few gridpts \n");
  for (int n=0;n<natoms;n++)
  {
    printf("   atom %i \n",n+1);
    int i1 = gc*gs*n;
    for (int j=0;j<3;j++)
      printf("  %8.5f %8.5f %8.5f \n",grid1[i1+gc*j+0],grid1[i1+gc*j+1],grid1[i1+gc*j+2]);
  }
 #endif


  int nat2 = natoms*natoms;

  float Ra[nat2];
 //relative atomic positions
  for (int i=0;i<natoms;i++)
  for (int j=0;j<i;j++)
  {
    float xa1 = coords[3*i+0]; float ya1 = coords[3*i+1]; float za1 = coords[3*i+2];
    float xa2 = coords[3*j+0]; float ya2 = coords[3*j+1]; float za2 = coords[3*j+2];
    float dax = xa2-xa1; float day = ya2-ya1; float daz = za2-za1;
    float ra1 = sqrtf(dax*dax+day*day+daz*daz);

    Ra[i*natoms+j] = ra1;
    Ra[j*natoms+i] = ra1;
  }

  float aij[nat2];
  for (int i=0;i<natoms;i++)
  for (int j=0;j<i;j++)
  {
    int Z1 = atno[i]; int Z2 = atno[j];
    float a12 = becke_a(Z1,Z2);
    float a21 = becke_a(Z2,Z1);
    aij[i*natoms+j] = a12;
    aij[j*natoms+i] = a21;
  }

 //rescaled partition parameters
  float bij[nat2];
  if (alpha>0.)
  {
    for (int i=0;i<nat2;i++)
      bij[i] = aij[i];
  }
  else
  for (int i=0;i<natoms;i++)
  for (int j=0;j<i;j++)
  {
    int Z1 = atno[i]; int Z2 = atno[j];
    float b12 = becke_a(-alpha,Z1,Z2);
    float b21 = becke_a(-alpha,Z2,Z1);
    bij[i*natoms+j] = b12;
    bij[j*natoms+i] = b21;
  }

 #if 0
  printf("  b: \n");
  print_square(natoms,bij);
 #endif

  //#pragma acc enter data copyin(aij[0:nat2],Ra[0:nat2])

  float fp[nat2];
  float fpb[nat2];
  //#pragma acc enter data create(fp[0:nat2],fpb[0:nat2])
  //#pragma acc parallel loop present(fp[0:nat2])
  for (int i=0;i<natoms;i++)
    fp[i*natoms+i] = fpb[i*natoms+i] = 1.;

  float tansc = alpha;
  if (alpha<=0.) tansc = 10.f;

  //#pragma acc parallel loop private(fp)
  for (int n=0;n<gsa;n++)
  {
    float x = grid1[gc*n+0]; float y = grid1[gc*n+1]; float z = grid1[gc*n+2];

    int wa = get_atom_grid(n,natoms,gsa);

    int i = ta;
    for (int j=0;j<natoms;j++)
    if (j!=i)
    //#pragma acc parallel loop collapse(2) present(coords[0:3*natoms],aij[0:nat2],Ra[0:nat2])
    //for (int i=0;i<natoms;i++)
    //for (int j=0;j<i;j++)
    {
      float a12 = aij[i*natoms+j];
      float a21 = aij[j*natoms+i];
      float b12 = bij[i*natoms+j];
      float b21 = bij[j*natoms+i];
      float oR12 = 1.f/Ra[i*natoms+j];

      float dx1 = x-coords[3*i+0]; float dy1 = y-coords[3*i+1]; float dz1 = z-coords[3*i+2];
      float dx2 = x-coords[3*j+0]; float dy2 = y-coords[3*j+1]; float dz2 = z-coords[3*j+2];
      float r1 = sqrtf(dx1*dx1+dy1*dy1+dz1*dz1);
      float r2 = sqrtf(dx2*dx2+dy2*dy2+dz2*dz2);

      float mu12 = (r1-r2)*oR12;
      float omm12 = 1.f-mu12*mu12;

      float nu12 =  mu12 + a12*omm12;
      float nu21 = -mu12 + a21*omm12;
      float f12 = bf3(nu12);
      float f21 = bf3(nu21);

      if (alpha>0.)
      {
        //f12 = tanhh(tansc*nu21); //note 12/21 swap
        //f21 = tanhh(tansc*nu12);
        f12 = sigmoidf(tansc*nu12);
        f21 = sigmoidf(tansc*nu21);
        double ch12 = cosh(fabs(nu12));
        double och = 1./ch12;
        f12 = och*och;
        f21 = och*och;
      }

      fp[i*natoms+j] = f12;
      fp[j*natoms+i] = f21;

     //"b" terms

      //const float nu0 = beta*oR12;
      const float nu0 = beta;
      if (alpha>0.)
      {
        //mu12 = r1-r2; //sigmoid not dependent on R12
        omm12 = 1.f-mu12*mu12;

        nu12 =  mu12 + b12*omm12 - nu0;
        nu21 = -mu12 + b21*omm12 + nu0;

        f12 = sigmoidf(tansc*nu12);
        f21 = sigmoidf(tansc*nu21);
        //f12 = tanhh(tansc*nu21); //note swap
        //f21 = tanhh(tansc*nu12);

        //f12 = bf3(nu12);
        //f21 = bf3(nu21);
      }
      else
      {
       //standard partition, except b term
        //omm12 = 1.f-mu12*mu12;
        nu12 =  mu12 + b12*omm12;
        nu21 = -mu12 + b21*omm12;

        f12 = bf3(nu12);
        f21 = bf3(nu21);
      }

      fpb[i*natoms+j] = f12;
      fpb[j*natoms+i] = f21;

    }

   //doing the Becke grid weight
    if (ta==-1)
    {
     //total weight does not normalize to unity 
     //  when alpha != 1.0
      float s1 = 1.;
      for (int j=0;j<natoms;j++)
        s1 *= fpb[wa*natoms+j];

      float norm = 1.e-15;
      //#pragma acc loop collapse(2)
      for (int i=0;i<natoms;i++)
      {
        float s0 = 1.;
        for (int j=0;j<natoms;j++)
          s0 *= fp[i*natoms+j];
        norm += s0;
      }

     //#pragma acc serial present(wt1[0:gsa])
      wt1[n] *= s1/norm;
    }
    else
    {
     //save the cell partition
      float s1 = 1.;
      for (int j=0;j<natoms;j++)
        s1 *= fpb[ta*natoms+j]; //ta here, wa above

      float norm = 1.e-15;
      if (0)
      for (int i=0;i<natoms;i++)
      {
        float s0 = 1.;
        for (int j=0;j<natoms;j++)
          s0 *= fp[i*natoms+j];
        norm += s0;
      }

      //wt1[n] = s1/norm;
      wt1[n] = s1;
    }

  }

  //#pragma acc exit data delete(fp[0:nat2],fpb[0:nat2])
  //#pragma acc exit data delete(aij[0:nat2],Ra[0:nat2])
  #pragma acc update device(wt1[0:gsa])

  return;
}

//n-center grid wts (float)
void becke_weight_nc(const int natoms, const int gc, const int gs, float* grid1, float* wt1, int* atno, float* coords)
{
  //printf("  testing becke_weight \n");

  if (natoms<2) return;

  int gsa = gs*natoms;
  int nat2 = natoms*natoms;

  float Ra[nat2];
 //atomic positions
  for (int i=0;i<natoms;i++)
  for (int j=0;j<i;j++)
  {
    float xa1 = coords[3*i+0]; float ya1 = coords[3*i+1]; float za1 = coords[3*i+2];
    float xa2 = coords[3*j+0]; float ya2 = coords[3*j+1]; float za2 = coords[3*j+2];
    float dax = xa2-xa1; float day = ya2-ya1; float daz = za2-za1;
    float ra1 = sqrtf(dax*dax+day*day+daz*daz);

    Ra[i*natoms+j] = ra1;
    Ra[j*natoms+i] = ra1;
    //printf(" Ra(%i-%i): %8.5f \n",i,j,ra1);
  }


  float aij[nat2];
  for (int i=0;i<natoms;i++)
  for (int j=0;j<i;j++)
  {
    int Z1 = atno[i]; int Z2 = atno[j];
    float a12 = becke_a(Z1,Z2);
    float a21 = becke_a(Z2,Z1);
    aij[i*natoms+j] = a12;
    aij[j*natoms+i] = a21;
  }

  //#pragma acc enter data copyin(aij[0:nat2],Ra[0:nat2])

  float fp[nat2];
  //#pragma acc enter data create(fp[0:nat2])

  //#pragma acc parallel loop private(fp)
  for (int n=0;n<gsa;n++)
  {
    float x = grid1[gc*n+0]; float y = grid1[gc*n+1]; float z = grid1[gc*n+2];

    int wa = get_atom_grid(n,natoms,gsa);

    //#pragma acc parallel loop present(fp[0:nat2])
    for (int i=0;i<natoms;i++)
      fp[i*natoms+i] = 1.;

    //#pragma acc parallel loop collapse(2) present(coords[0:3*natoms],aij[0:nat2],Ra[0:nat2])
    for (int i=0;i<natoms;i++)
    for (int j=0;j<i;j++)
    {
      float a12 = aij[i*natoms+j];
      float a21 = aij[j*natoms+i];
      float oR12 = 1.f/Ra[i*natoms+j];

      float dx1 = x-coords[3*i+0]; float dy1 = y-coords[3*i+1]; float dz1 = z-coords[3*i+2];
      float dx2 = x-coords[3*j+0]; float dy2 = y-coords[3*j+1]; float dz2 = z-coords[3*j+2];
      float r1 = sqrtf(dx1*dx1+dy1*dy1+dz1*dz1);
      float r2 = sqrtf(dx2*dx2+dy2*dy2+dz2*dz2);

      float mu12 = (r1-r2)*oR12;
      float omm12 = 1.f-mu12*mu12;
      float nu12 =  mu12 + a12*omm12;
      float nu21 = -mu12 + a21*omm12;
      float f12 = bf3(nu12);
      float f21 = bf3(nu21);

      fp[i*natoms+j] = f12;
      fp[j*natoms+i] = f21;

      //printf(" r1/2: %8.5f %8.5f s1: %8.5f \n",r1,r2,f12/(f12+f21));
    }

    float s1;
    float norm = 1.e-15;
    //#pragma acc loop collapse(2)
    for (int i=0;i<natoms;i++)
    {
      float s0 = 1.;
      for (int j=0;j<natoms;j++)
        s0 *= fp[i*natoms+j];
      norm += s0;

      if (i==wa)
        s1 = s0;
    }

   //#pragma acc serial present(wt1[0:gsa])
    wt1[n] *= s1/norm;
  }

  //#pragma acc exit data delete(fp[0:nat2])
  //#pragma acc exit data delete(aij[0:nat2],Ra[0:nat2])
  #pragma acc update device(wt1[0:gsa])

  return;
}

//n-center grid wts (double)
void becke_weight_nc(const int natoms, const int gc, const int gs, double* grid1, double* wt1, int* atno, float* coords)
{
  //printf("  testing becke_weight \n");

  if (natoms<2) return;

  int gsa = gs*natoms;
  int nat2 = natoms*natoms;

  double Ra[nat2];
 //atomic positions
  for (int i=0;i<natoms;i++)
  for (int j=0;j<i;j++)
  {
    double xa1 = coords[3*i+0]; double ya1 = coords[3*i+1]; double za1 = coords[3*i+2];
    double xa2 = coords[3*j+0]; double ya2 = coords[3*j+1]; double za2 = coords[3*j+2];
    double dax = xa2-xa1; double day = ya2-ya1; double daz = za2-za1;
    double ra1 = sqrtf(dax*dax+day*day+daz*daz);

    Ra[i*natoms+j] = ra1;
    Ra[j*natoms+i] = ra1;
    //printf(" Ra(%i-%i): %8.5f \n",i,j,ra1);
  }


  double aij[nat2];
  for (int i=0;i<natoms;i++)
  for (int j=0;j<i;j++)
  {
    int Z1 = atno[i]; int Z2 = atno[j];
    double a12 = becke_a(Z1,Z2);
    double a21 = becke_a(Z2,Z1);
    aij[i*natoms+j] = a12;
    aij[j*natoms+i] = a21;
  }

  //#pragma acc enter data copyin(aij[0:nat2],Ra[0:nat2])

  double fp[nat2];
  //#pragma acc enter data create(fp[0:nat2])

  //#pragma acc parallel loop private(fp)
  for (int n=0;n<gsa;n++)
  {
    double x = grid1[gc*n+0]; double y = grid1[gc*n+1]; double z = grid1[gc*n+2];

    int wa = get_atom_grid(n,natoms,gsa);

    //#pragma acc parallel loop present(fp[0:nat2])
    for (int i=0;i<natoms;i++)
      fp[i*natoms+i] = 1.;

    //#pragma acc parallel loop collapse(2) present(coords[0:3*natoms],aij[0:nat2],Ra[0:nat2])
    for (int i=0;i<natoms;i++)
    for (int j=0;j<i;j++)
    {
      double a12 = aij[i*natoms+j];
      double a21 = aij[j*natoms+i];
      double oR12 = 1./Ra[i*natoms+j];

      double dx1 = x-coords[3*i+0]; double dy1 = y-coords[3*i+1]; double dz1 = z-coords[3*i+2];
      double dx2 = x-coords[3*j+0]; double dy2 = y-coords[3*j+1]; double dz2 = z-coords[3*j+2];
      double r1 = sqrt(dx1*dx1+dy1*dy1+dz1*dz1);
      double r2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2);

      double mu12 = (r1-r2)*oR12;
      double omm12 = 1.-mu12*mu12;
      double nu12 =  mu12 + a12*omm12;
      double nu21 = -mu12 + a21*omm12;
      double f12 = bf3d(nu12);
      double f21 = bf3d(nu21);

      fp[i*natoms+j] = f12;
      fp[j*natoms+i] = f21;

      //printf(" r1/2: %8.5f %8.5f s1: %8.5f \n",r1,r2,f12/(f12+f21));
    }

    double s1;
    double norm = 1.e-20;
    //#pragma acc loop collapse(2)
    for (int i=0;i<natoms;i++)
    {
      double s0 = 1.;
      for (int j=0;j<natoms;j++)
        s0 *= fp[i*natoms+j];
      norm += s0;

      if (i==wa)
        s1 = s0;
    }

   //#pragma acc serial present(wt1[0:gsa])
    wt1[n] *= s1/norm;
  }

  //#pragma acc exit data delete(aij[0:nat2],Ra[0:nat2])
  #pragma acc update device(wt1[0:gsa])

  return;
}

#pragma acc routine seq
float get_3c_bw(float r1, float r2, float r3, float a12, float a13, float a23, float a21, float a31, float a32, float oR12, float oR13, float oR23)
{
  float mu12 = (r1-r2)*oR12;
  float mu13 = (r1-r3)*oR13;
  float mu23 = (r2-r3)*oR23;
  float omm12 = 1.f-mu12*mu12;
  float omm13 = 1.f-mu13*mu13;
  float omm23 = 1.f-mu23*mu23;

 //original
  float nu12 =  mu12 + a12*omm12;
  float nu21 = -mu12 + a21*omm12;
  float nu13 =  mu13 + a13*omm13;
  float nu31 = -mu13 + a31*omm13;
  float nu23 =  mu23 + a23*omm23;
  float nu32 = -mu23 + a32*omm23;

  float f12 = bf3(nu12);
  float f21 = bf3(nu21);
  float f13 = bf3(nu13);
  float f31 = bf3(nu31);
  float f23 = bf3(nu23);
  float f32 = bf3(nu32);

  float norm = f12*f13 + f21*f23 + f31*f32 + 1.e-10f;
  float s1 = f12*f13/norm;

  return s1;
}

void becke_weight_3c(int gs, float* grid1, float* wt1, float* grid2, float* wt2, float* grid3, float* wt3, 
                     int Z1, int Z2, int Z3, float A2, float B2, float C2, float A3, float B3, float C3)
{
  float A23 = A3-A2; float B23 = B3-B2; float C23 = C3-C2;

  const float R12 = sqrtf(A2*A2+B2*B2+C2*C2);
  const float R13 = sqrtf(A3*A3+B3*B3+C3*C3);
  const float R23 = sqrtf(A23*A23+B23*B23+C23*C23);
  const float oR12 = 1./R12;
  const float oR13 = 1./R13;
  const float oR23 = 1./R23;

  float a12 = becke_a(Z1,Z2);
  float a13 = becke_a(Z1,Z3);
  float a23 = becke_a(Z2,Z3);
  float a21 = becke_a(Z2,Z1);
  float a31 = becke_a(Z3,Z1);
  float a32 = becke_a(Z3,Z2); 

#if USE_ACC
 #pragma acc parallel loop independent present(grid1[0:6*gs],wt1[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r1 = grid1[6*i+3];
    float r2 = grid1[6*i+4];
    float r3 = grid1[6*i+5];
    float s1 = get_3c_bw(r1,r2,r3,a12,a13,a23,a21,a31,a32,oR12,oR13,oR23);

    wt1[i] *= s1;
  }

#if USE_ACC
 #pragma acc parallel loop independent present(grid2[0:6*gs],wt2[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r1 = grid2[6*i+3];
    float r2 = grid2[6*i+4];
    float r3 = grid2[6*i+5];
    float s1 = get_3c_bw(r1,r2,r3,a21,a23,a13,a12,a32,a31,oR12,oR23,oR13);

    wt2[i] *= s1;
  }

#if USE_ACC
 #pragma acc parallel loop independent present(grid3[0:6*gs],wt3[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r1 = grid3[6*i+3];
    float r2 = grid3[6*i+4];
    float r3 = grid3[6*i+5];
    float s1 = get_3c_bw(r1,r2,r3,a32,a31,a21,a23,a13,a12,oR23,oR13,oR12);

    wt3[i] *= s1;
  }

  return;
}

void becke_weight_2c(int gs, float* grid1, float* wt1, float* grid2, float* wt2, 
                     int Z1, int Z2, float A2, float B2, float C2)
{
 //first center is at 0,0,0 and second at A2,B2,C2
  float R  = sqrtf(A2*A2+B2*B2+C2*C2);
  const float oR = 1./R;

  const float a1 = becke_a(Z1,Z2);
  const float a2 = becke_a(Z2,Z1);

  //printf(" a1/2: %8.5f %8.5f \n",a1,a2);

#if USE_ACC
 #pragma acc parallel loop independent present(grid1[0:6*gs],wt1[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r1 = grid1[6*i+3];
    float r2 = grid1[6*i+4];
    float mu = (r1-r2)*oR;
    float omm2 = 1.f-mu*mu;
    float nu1 =  mu + a1*omm2;
    float nu2 = -mu + a2*omm2;

    float f1 = bf3(nu1);
    float f2 = bf3(nu2);

    float norm = f1 + f2 + 1.e-10f;
    float s1 = f1/norm;

    //printf(" r1/2: %8.5f %8.5f s1: %8.5f \n",r1,r2,s1);
    wt1[i] *= s1;
  }

#if USE_ACC
 #pragma acc parallel loop independent present(grid2[0:6*gs],wt2[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float r1 = grid2[6*i+3];
    float r2 = grid2[6*i+4];
    float mu = (r1-r2)*oR;
    float omm2 = 1.f-mu*mu;
    float nu1 =  mu + a2*omm2;
    float nu2 = -mu + a1*omm2;

    float f1 = bf3(nu1);
    float f2 = bf3(nu2);

    float norm = f1 + f2 + 1.e-10f;
    float s1 = f1/norm;

    //printf(" r1/2: %8.5f %8.5f s1: %8.5f \n",r1,r2,s1);
    wt2[i] *= s1;
  }

  return;
}

void becke_weight_2d(int gs, double* grid1, double* wt1, double* grid2, double* wt2, 
                     double zeta1, double zeta2, double A2, double B2, double C2)
{
 //first center is at 0,0,0 and second at A2,B2,C2
  double R  = sqrt(A2*A2+B2*B2+C2*C2);
  const double oR = 1./R;

  const double a1 = becke_a_zeta(zeta1,zeta2);
  const double a2 = becke_a_zeta(zeta2,zeta1);

  //printf(" a1/2: %8.5f %8.5f \n",a1,a2);

#if USE_ACC
 #pragma acc parallel loop independent present(grid1[0:6*gs],wt1[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    double r1 = grid1[6*i+3];
    double r2 = grid1[6*i+4];
    double mu = (r1-r2)*oR;
    double omm2 = 1.f-mu*mu;
    double nu1 =  mu + a1*omm2;
    double nu2 = -mu + a2*omm2;

    double f1 = bf3d(nu1);
    double f2 = bf3d(nu2);

    double norm = f1 + f2 + 1.e-10f;
    double s1 = f1/norm;

    //printf(" r1/2: %8.5f %8.5f s1: %8.5f \n",r1,r2,s1);
    wt1[i] *= s1;
  }

#if USE_ACC
 #pragma acc parallel loop independent present(grid2[0:6*gs],wt2[0:gs])
#endif
  for (int i=0;i<gs;i++)
  {
    double r1 = grid2[6*i+3];
    double r2 = grid2[6*i+4];
    double mu = (r1-r2)*oR;
    double omm2 = 1.f-mu*mu;
    double nu1 =  mu + a2*omm2;
    double nu2 = -mu + a1*omm2;

    double f1 = bf3d(nu1);
    double f2 = bf3d(nu2);

    double norm = f1 + f2 + 1.e-10f;
    double s1 = f1/norm;

    //printf(" r1/2: %8.5f %8.5f s1: %8.5f \n",r1,r2,s1);
    wt2[i] *= s1;
  }

  return;
}

void test_2c_becke()
{
 #if 0
  if (natoms==2)
  {
    int m1 = 0; int m2 = 1;
    float Z1 = (float)atno[m1]; float Z2 = (float)atno[m2];
    float A1 = coords[3*m1+0]; float B1 = coords[3*m1+1]; float C1 = coords[3*m1+2];
    float A2 = coords[3*m2+0]; float B2 = coords[3*m2+1]; float C2 = coords[3*m2+2];
    float A12 = A2-A1; float B12 = B2-B1; float C12 = C2-C1;

    float* grid1 = new float[gs6];
    float* wt1 = new float[gs];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);
    add_r2_to_grid(gs,grid1,A12,B12,C12);

    float* grid2 = new float[gs6];
    float* wt2 = new float[gs];
    generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);

    recenter_grid(gs,grid2,A12,B12,C12);

    double* wtall = new double[natoms*gs];
    double* gridall = new double[natoms*gsc];
    for (int m=0;m<gs;m++) wtall[m] = wt1[m];
    for (int m=0;m<gs;m++) wtall[gs+m] = wt2[m];

    for (int m=0;m<gs;m++)
    {
      gridall[gc*m+0] = grid1[6*m+0];
      gridall[gc*m+1] = grid1[6*m+1];
      gridall[gc*m+2] = grid1[6*m+2];
    }
    for (int m=0;m<gs;m++)
    {
      gridall[gs3+gc*m+0] = grid2[6*m+0];
      gridall[gs3+gc*m+1] = grid2[6*m+1];
      gridall[gs3+gc*m+2] = grid2[6*m+2];
    }

    becke_weight_2c(gs,grid1,wt1,grid2,wt2,Z1,Z2,A12,B12,C12);
    

    if (0)
    {
      printf("  grid1: \n");
      for (int m=0;m<gs;m++)
        printf(" %8.5f %8.5f %8.5f \n",grid1[6*m+0],grid1[6*m+1],grid1[6*m+2]);
      printf("\n");
    }
    if (0)
    {
      printf("  grid2: \n");
      for (int m=0;m<gs;m++)
        printf(" %8.5f %8.5f %8.5f \n",grid2[6*m+0],grid2[6*m+1],grid2[6*m+2]);
      printf("\n");
    }

    printf("  wt1: ");
    for (int m=0;m<gs;m++)
      printf(" %4.2e",wt1[m]);
    printf("\n");
    printf("  wt2: ");
    for (int m=0;m<gs;m++)
      printf(" %4.2e",wt2[m]);
    printf("\n");

    becke_weight_nc(natoms,gc,gs,gridall,wtall,atno,coords);

    printf("  wt1: ");
    for (int m=0;m<gs;m++)
      printf(" %4.2e",wtall[m]);
    printf("\n");
    printf("  wt2: ");
    for (int m=0;m<gs;m++)
      printf(" %4.2e",wtall[gs+m]);
    printf("\n");

  }
 #endif
  return;
}

double becke_a_zeta(double zeta1, double zeta2)
{
  double x = zeta2/zeta1;

  float u = (x-1.)/(x+1.);
  float a = u/(u*u-1.);

  if (a>0.5) a = 0.5;
  else if (a<-0.5) a = -0.5;
  return a;
}

float becke_a(float alpha, int Z1, int Z2)
{
  float x = alpha*get_radii_2(Z1)/get_radii_2(Z2);

  float u = (x-1.f)/(x+1.f);
  float a = u/(u*u-1.f);

  if (a>0.5f) a = 0.5f;
  else if (a<-0.5f) a = -0.5f;
  return a;
}

float becke_a(int Z1, int Z2)
{
  float x = get_radii_2(Z1)/get_radii_2(Z2);

  float u = (x-1.f)/(x+1.f);
  float a = u/(u*u-1.f);

  if (a>0.5f) a = 0.5f;
  else if (a<-0.5f) a = -0.5f;
  return a;
}

#pragma acc routine seq
float bf3(float f1)
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
double bf3d(double f1)
{
  //default is order 3
  #define ORDER4 0

  f1 = 1.5*f1 - 0.5*f1*f1*f1;
  f1 = 1.5*f1 - 0.5*f1*f1*f1;
  f1 = 1.5*f1 - 0.5*f1*f1*f1;
 #if ORDER4
  f1 = 1.5*f1 - 0.5*f1*f1*f1;
 #endif
  f1 = 0.5*(1.-f1);
  f1 = f1*f1;

  return f1;
}

void get_becke_grid_full(float alpha, int natoms, int* atno, float* coords, int nrad, int nang, float* ang_g, float* ang_w, const int gc, float* grid, float* wt)
{
  int gs = nrad*nang;
  int gsa = gs*natoms;
  int gsc = gc*gs;
  int gsac = gsc*natoms;
  //int N3 = 3*natoms;

  for (int n=0;n<natoms;n++)
  {
    float r1[nrad]; float w1[nrad];
    int Z1 = atno[n];
    #pragma acc enter data create(r1[0:nrad],w1[0:nrad])
    get_murak_grid_f(nrad,r1,w1,Z1,3);

    if (0)
    {
      printf(" r: ");
      for (int j=0;j<nrad;j++)
        printf(" %8.5f",r1[j]);
      printf("\n");
    }

    #pragma acc parallel loop present(r1[0:nrad],w1[0:nrad],ang_g[0:3*nang],ang_w[0:nang],grid[0:gsac],wt[0:gsa])
    for (int m=0;m<nrad;m++)
    {
      double r0 = r1[m];
      double wr0 = w1[m];

      int wg0 = n*gs+m*nang;
      #pragma acc loop
      for (int p=0;p<nang;p++)
      {
        int wg = wg0 + p;
        grid[gc*wg+0] = r0*ang_g[3*p+0];
        grid[gc*wg+1] = r0*ang_g[3*p+1];
        grid[gc*wg+2] = r0*ang_g[3*p+2];
        //grid[gc*wg+3] = r0;
        wt[wg] = wr0*ang_w[p];
        //wg++;
      }
    } //loop r over radial points

    double A0 = coords[3*n+0]; double B0 = coords[3*n+1]; double C0 = coords[3*n+2];
    float* grid1 = &grid[gsc*n];

    #pragma acc parallel loop present(grid1[0:gsc])
    for (int m=0;m<gs;m++)
    {
      grid1[gc*m+0] += A0;
      grid1[gc*m+1] += B0;
      grid1[gc*m+2] += C0;
    }

    #pragma acc exit data delete(r1[0:nrad],w1[0:nrad])
  } //loop n over natoms

  #pragma acc update self(grid[0:gsac],wt[0:gsa])

  float beta = 0.; //not implemented here
  if (alpha!=0.)
    atomic_domain_cell_wt(-1,alpha,beta,natoms,gc,gs,grid,wt,atno,coords);
  else
    becke_weight_nc(natoms,gc,gs,grid,wt,atno,coords);

  //distance from center (arbitrary center)
  #pragma acc parallel loop present(grid[0:gsac])
  for (int m=0;m<gsa;m++)
  {
    float x1 = grid[6*m+0]; float y1 = grid[6*m+1]; float z1 = grid[6*m+2];
    float r1 = sqrtf(x1*x1+y1*y1+z1*z1);
    grid[gc*m+3] = r1;
  }

  #pragma acc update self(grid[0:gsac])

  return;
}

void get_becke_grid_full(int natoms, int* atno, float* coords, int nrad, int nang, float* ang_g, float* ang_w, const int gc, float* grid, float* wt)
{
  return get_becke_grid_full(0.,natoms,atno,coords,nrad,nang,ang_g,ang_w,gc,grid,wt);
}

void get_becke_grid_full(int natoms, int* atno, float* coords, int nrad, int nang, double* ang_g, double* ang_w, const int gc, double* grid, double* wt)
{
  int gs = nrad*nang;
  int gsa = gs*natoms;
  int gsc = gc*gs;
  int gsac = gsc*natoms;

  for (int n=0;n<natoms;n++)
  {
    double r1[nrad]; double w1[nrad];
    int Z1 = atno[n];
    #pragma acc enter data create(r1[0:nrad],w1[0:nrad])
    get_murak_grid(nrad,r1,w1,Z1,3);

    #pragma acc parallel loop present(r1[0:nrad],w1[0:nrad],ang_g[0:3*nang],ang_w[0:nang],grid[0:gsac],wt[0:gsa])
    for (int m=0;m<nrad;m++)
    {
      double r0 = r1[m];
      double wr0 = w1[m];

      int wg0 = n*gs+m*nang;
      #pragma acc loop
      for (int p=0;p<nang;p++)
      {
        int wg = wg0 + p;
        grid[gc*wg+0] = r0*ang_g[3*p+0];
        grid[gc*wg+1] = r0*ang_g[3*p+1];
        grid[gc*wg+2] = r0*ang_g[3*p+2];
        wt[wg] = wr0*ang_w[p];
      }
    } //loop r over radial points

    double A0 = coords[3*n+0]; double B0 = coords[3*n+1]; double C0 = coords[3*n+2];
    double* grid1 = &grid[gsc*n];

    #pragma acc parallel loop present(grid1[0:gsc])
    for (int m=0;m<gs;m++)
    {
      grid1[gc*m+0] += A0;
      grid1[gc*m+1] += B0;
      grid1[gc*m+2] += C0;
    }

    #pragma acc exit data delete(r1[0:nrad],w1[0:nrad])
  } //loop n over natoms

  #pragma acc update self(grid[0:gsac],wt[0:gsa])
  becke_weight_nc(natoms,gc,gs,grid,wt,atno,coords);

  //distance from center (arbitrary center)
  #pragma acc parallel loop present(grid[0:gsac])
  for (int m=0;m<gsa;m++)
  {
    double x1 = grid[6*m+0]; double y1 = grid[6*m+1]; double z1 = grid[6*m+2];
    double r1 = sqrt(x1*x1+y1*y1+z1*z1);
    grid[gc*m+3] = r1;
  }

  #pragma acc update self(grid[0:gsac])

  return;
}

void get_becke_grid_full(int natoms, int* atno, float* coords, int nrad, int nang, double* ang_g, double* ang_w, const int gc, float* grid, float* wt)
{
  float* ang_gf = new float[3*nang];
  float* ang_wf = new float[nang];

  #pragma acc enter data create(ang_gf[0:3*nang],ang_wf[0:nang])

  #pragma acc parallel loop present(ang_gf[0:3*nang],ang_g[0:3*nang])
  for (int m=0;m<3*nang;m++)
    ang_gf[m] = ang_g[m];
  #pragma acc parallel loop present(ang_wf[0:nang],ang_w[0:nang])
  for (int m=0;m<nang;m++)
    ang_wf[m] = ang_w[m];

  get_becke_grid_full(natoms,atno,coords,nrad,nang,ang_gf,ang_wf,gc,grid,wt);

  #pragma acc exit data delete(ang_gf[0:3*nang],ang_wf[0:nang])

  delete [] ang_gf;
  delete [] ang_wf;
}

void get_becke_grid_full(int natoms, int* atno, double* coords, int nrad, int nang, double* ang_g, double* ang_w, const int gc, double* grid, double* wt)
{
  float coordsf[3*natoms];
  for (int j=0;j<3*natoms;j++)
    coordsf[j] = coords[j];
  return get_becke_grid_full(natoms,atno,coordsf,nrad,nang,ang_g,ang_w,gc,grid,wt);
}

void add_r1_to_grid4(int gs, float* grid1, float A2, float B2, float C2)
{
#if USE_ACC
 #pragma acc parallel loop independent present(grid1[0:4*gs])
#endif
  for (int i=0;i<gs;i++)
  {
    float x1 = grid1[4*i+0]; 
    float y1 = grid1[4*i+1];
    float z1 = grid1[4*i+2];
    float x12 = x1-A2;
    float y12 = y1-B2;
    float z12 = z1-C2;
    float r1 = sqrtf(x12*x12+y12*y12+z12*z12);
    grid1[4*i+3] = r1;
  }
    
  return;
}

void density_in_basis2(int natoms, int* atno, float* coords, vector<vector<double> > &basis1, vector<vector<double> > &basis2, int No, double* jCA, int gsa, float* grid, float* wt, double* Paom, int prl)
{
  if (basis2.size()<1) return;

 //basis1 --> basis for No MOs
 //basis2 --> new AO basis

  if (prl>-1) printf("\n density_in_basis2 (No: %i) \n",No);

  //int gs = gsa/natoms;
  int gsa6 = 6*gsa;

  int N1 = basis1.size();
  //int N12 = N1*N1;

  if (No>N1) { printf("\n ERROR: No cannot exceed basis1.size() in compute_delta \n"); return; }

  int* n2i1 = new int[natoms];
  int imaxN1 = get_imax_n2i(natoms,N1,basis1,n2i1);

  int N2 = basis2.size();
  int N22 = N2*N2;

  for (int i=0;i<N22;i++) 
    Paom[i] = 0.;

  int* n2i2 = new int[natoms];
  int imaxN2 = get_imax_n2i(natoms,N2,basis2,n2i2);

  float* grid1 = new float[gsa6];

  int iN1 = N1; //working with MOs, so need all AOs at once
  int iN2 = N2;
  float** val1 = new float*[iN1];
  float** val2 = new float*[iN2];
  for (int i=0;i<iN1;i++)
    val1[i] = new float[gsa];
  for (int i=0;i<iN2;i++)
    val2[i] = new float[gsa];

  float** val3 = new float*[No];
  for (int i=0;i<No;i++)
    val3[i] = new float[gsa];

  #pragma acc enter data create(grid1[0:gsa6])
  #pragma acc enter data create(val1[0:iN1][0:gsa],val2[0:iN2][0:gsa],val3[0:No][0:gsa])

 //construct AO quantities for basis1
  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i1[m-1]; int s2 = n2i1[m];
  
    float Z1 = (float)atno[m];
    float A1 = coords[3*m+0]; float B1 = coords[3*m+1]; float C1 = coords[3*m+2];

    copy_grid(gsa,grid1,grid);
    recenter_grid_zero(gsa,grid1,-A1,-B1,-C1);

    #pragma acc parallel loop collapse(2) present(val1[0:iN1][0:gsa])
    for (int i1=s1;i1<s2;i1++)
    for (int j=0;j<gsa;j++)
      val1[i1][j] = 1.f;

    for (int i1=s1;i1<s2;i1++)
    {
      vector<double> basisn = basis1[i1];
      int n1 = basisn[0]; int l1 = basisn[1]; int m1 = basisn[2]; double zeta1 = basisn[3];
      float nm1 = basisn[4];

      eval_sh(i1,gsa,grid1,val1[i1],n1,l1,m1,zeta1);
      #pragma acc parallel loop present(val1[0:iN1][0:gsa])
      for (int j=0;j<gsa;j++)
        val1[i1][j] *= nm1;
    } //loop i1

  } //loop m over natoms


 //construct AO quantities for basis2
  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i2[m-1]; int s2 = n2i2[m];
  
    float Z1 = (float)atno[m];
    float A1 = coords[3*m+0]; float B1 = coords[3*m+1]; float C1 = coords[3*m+2];

    copy_grid(gsa,grid1,grid);
    recenter_grid_zero(gsa,grid1,-A1,-B1,-C1);

    #pragma acc parallel loop collapse(2) present(val2[0:iN2][0:gsa])
    for (int i2=s1;i2<s2;i2++)
    for (int j=0;j<gsa;j++)
      val2[i2][j] = 1.f;

    for (int i2=s1;i2<s2;i2++)
    {
      vector<double> basisn = basis2[i2];
      int n1 = basisn[0]; int l1 = basisn[1]; int m1 = basisn[2]; double zeta1 = basisn[3];
      float nm1 = basisn[4];

      eval_sh(i2,gsa,grid1,val2[i2],n1,l1,m1,zeta1);
      #pragma acc parallel loop present(val2[0:iN2][0:gsa])
      for (int j=0;j<gsa;j++)
        val2[i2][j] *= nm1;

    } //loop i2

  } //loop m over natoms


 //get the MOs (basis1)
  #pragma acc parallel loop collapse(2) present(val1[0:iN1][0:gsa],val3[0:No][0:gsa],jCA[0:N2])
  for (int i3=0;i3<No;i3++)
  for (int j=0;j<gsa;j++)
  {
    float v2 = 0.f;
    #pragma acc loop reduction(+:v2)
    for (int k=0;k<iN1;k++)
      v2 += jCA[k*N1+i3]*val1[k][j];
    val3[i3][j] = v2;
  }

  int No1 = No;
  double b12[No1*N2];

 //projection of orbital-wise density into basis 2
  for (int i3=0;i3<No1;i3++)
  {
    float* valm = val3[i3]; //MO

   //inner loop is over pairs of basis functions
    for (int i2=0;i2<N2;i2++)
    {
      float* valn = val2[i2]; //AOn

      double val = 0.;
      #pragma acc parallel loop present(valm[0:gsa],valn[0:gsa],wt[0:gsa]) reduction(+:val)
      for (int j=0;j<gsa;j++)
        val += valm[j]*valn[j]*wt[j];

      b12[i3*N2+i2] = val;
    }
  }

  for (int i1=0;i1<N2;i1++)
  for (int i2=0;i2<N2;i2++)
  {
    double valt = 0.;
    for (int i3=0;i3<No1;i3++)
      valt += b12[i3*N2+i1]*b12[i3*N2+i2];

    Paom[i1*N2+i2] = valt;
  }


 //cleanup
  #pragma acc exit data delete(grid1[0:gsa6])
  #pragma acc exit data delete(val1[0:iN1][0:gsa],val2[0:iN2][0:gsa],val3[0:No][0:gsa])

  delete [] n2i1;
  delete [] n2i2;

  delete [] grid1;

  for (int i=0;i<iN1;i++)
    delete [] val1[i];
  for (int i=0;i<iN2;i++)
    delete [] val2[i];
  for (int i=0;i<No;i++)
    delete [] val3[i];
  delete [] val1;
  delete [] val2;
  delete [] val3;

  return;
}

void batomic_charges(float alpha, int natoms, int* atno, float* coords, int nrad, int nang, float* ang_g, float* ang_w, float* chg, double* rho, int prl)
{
  int gc = 6;
  int gs = nrad*nang;
  int gsa = natoms*gs;
  float* grid = new float[gc*gsa];
  float* wt = new float[gsa];

  #pragma acc enter data create(grid[0:6*gsa],wt[0:gsa])

  get_becke_grid_full(alpha,natoms,atno,coords,nrad,nang,ang_g,ang_w,gc,grid,wt);

  for (int n=0;n<natoms;n++)
  {
    double* rho1 = &rho[n*gs];
    float* wt1 = &wt[n*gs];

    double dt = 0.;
    #pragma acc parallel loop present(rho1[0:gs],wt1[0:gs]) reduction(+:dt)
    for (int j=0;j<gs;j++)
      dt += rho1[j]*wt1[j];

    chg[n] = dt;
  
    if (prl>0)
      printf("    atom %2i batomic charge: %8.6f \n",n+1,dt);
  }

  #pragma acc exit data delete(grid[0:6*gsa],wt[0:gsa])

  delete [] grid;
  delete [] wt;

  return;
}

void atomic_charges(int natoms, int gs, float* chg, double* rho, float* zta, float** gridall, float** wta, int prl)
{
 //integrate rho(r) exp(-zeta r) on a per atom basis
  int gsa = gs*natoms;

  float* val1 = new float[gsa];
  #pragma acc enter data create(val1[0:gsa])

 //wta --> atomic wts without renormalization
  for (int n=0;n<natoms;n++)
  {
    double zeta1 = zta[n];
    double* rho1 = rho;
    //double* rho1 = &rho[n*gs];
    float* grid1 = gridall[n];
    float* wt1 = wta[n];

    double norm1 = norm(1,0,0,zeta1/2.);

    #pragma acc parallel loop present(val1[0:gsa],wt1[0:gsa])
    for (int j=0;j<gsa;j++)
      val1[j] = 1.f;
    eval_sh(0,gsa,grid1,val1,1,0,0,zeta1);

   //need to think through this function
    double dt = 0.;
    #pragma acc parallel loop present(rho1[0:gsa],val1[0:gsa],wt1[0:gsa]) reduction(+:dt)
    for (int j=0;j<gsa;j++)
      dt += sqrtf(rho1[j]*val1[j])*wt1[j];
    dt *= norm1;

    chg[n] = dt;
  
    if (prl>0)   
      printf("    atom %2i atomic charge: %8.6f  (zeta: %8.5f norm: %8.5f) \n",n+1,dt,zeta1,norm1);
  }

  #pragma acc exit data delete(val1[0:gsa])
  delete [] val1;

  return;
}

void becke_charges(int natoms, int gs, double* chg, double* rho, double* wt, int prl)
{
  if (prl>0)
    printf("    Becke charges:  ");

  for (int n=0;n<natoms;n++)
  {
    double* rho1 = &rho[n*gs];
    double* wt1 = &wt[n*gs];

    double dt = 0.;
    #pragma acc parallel loop present(rho1[0:gs],wt1[0:gs]) reduction(+:dt)
    for (int j=0;j<gs;j++)
      dt += rho1[j]*wt1[j];

    chg[n] = dt;

    if (prl>0)
      printf(" %8.6f",dt);
  }
  if (prl>0)
    printf("\n");

  return;
}

void becke_charges(int natoms, int gs, float* chg, double* rho, float* wt, int prl)
{
  if (prl>0)
    printf("    Becke charges:  ");

  for (int n=0;n<natoms;n++)
  {
    double* rho1 = &rho[n*gs];
    float* wt1 = &wt[n*gs];

    double dt = 0.;
    #pragma acc parallel loop present(rho1[0:gs],wt1[0:gs]) reduction(+:dt)
    for (int j=0;j<gs;j++)
      dt += rho1[j]*wt1[j];

    chg[n] = dt;

    if (prl>0)
      printf(" %8.6f",dt);
  }
  if (prl>0)
    printf("\n");

  return;
}
