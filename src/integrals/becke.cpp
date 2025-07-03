#include "becke.h"
#include "braggslater.h"
#include "murak.h"
#include "integrals.h"
#include "read.h"
#include "jellium_ints.h"
#include "print.h"
//#include "vxc.h"
#include "gauss.h"

#define pr_inc 5

void generate_central_grid_2(float* grid1, float* wt1, float Z1, int nrad, int nang, float* ang_g, float* ang_w);
void copy_grid(int gs, float* grid1, float* grid2);
void recenter_grid(int gs, float* grid, float x2, float y2, float z2);
void recenter_grid_zero(int gs, float* grid, float x2, float y2, float z2);
void add_r1_to_grid(int gs, float* grid1, float A2, float B2, float C2);
void add_r2_to_grid(int gs, float* grid1, float A2, float B2, float C2);
void print_square(int N, double* A);
void print_vxc(int nrad, int nang, int natoms, float* grid, double* vxc);
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

double sigmoid(double v1)
{
  double val = 1./(1.f+exp(-v1));
  return 1. - val;
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

void atomic_domain_cell_wt(const int ta, const float alpha, const float beta, const int natoms, const int gc, const int gs, double* grid1, double* wt1, int* atno, double* coords0)
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

  double cx0 = 0.; double cy0 = 0.; double cz0 = 0.;
  if (ta>=0)
  {
    cx0 = coords0[3*ta+0];
    cy0 = coords0[3*ta+1];
    cz0 = coords0[3*ta+2];
  }
  double coords[3*natoms];
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

  double Ra[nat2];
 //relative atomic positions
  for (int i=0;i<natoms;i++)
  for (int j=0;j<i;j++)
  {
    double xa1 = coords[3*i+0]; double ya1 = coords[3*i+1]; double za1 = coords[3*i+2];
    double xa2 = coords[3*j+0]; double ya2 = coords[3*j+1]; double za2 = coords[3*j+2];
    double dax = xa2-xa1; double day = ya2-ya1; double daz = za2-za1;
    double ra1 = sqrt(dax*dax+day*day+daz*daz);

    Ra[i*natoms+j] = ra1;
    Ra[j*natoms+i] = ra1;
  }

  double aij[nat2];
  for (int i=0;i<natoms;i++)
  for (int j=0;j<i;j++)
  {
    int Z1 = atno[i]; int Z2 = atno[j];
    double a12 = becke_ad(Z1,Z2);
    double a21 = becke_ad(Z2,Z1);
    aij[i*natoms+j] = a12;
    aij[j*natoms+i] = a21;
  }

 //rescaled partition parameters
  double bij[nat2];
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
    double b12 = becke_a(-alpha,Z1,Z2);
    double b21 = becke_a(-alpha,Z2,Z1);
    bij[i*natoms+j] = b12;
    bij[j*natoms+i] = b21;
  }

 #if 0
  printf("  b: \n");
  print_square(natoms,bij);
 #endif

  //#pragma acc enter data copyin(aij[0:nat2],Ra[0:nat2])

  double fp[nat2];
  double fpb[nat2];
  //#pragma acc enter data create(fp[0:nat2],fpb[0:nat2])
  //#pragma acc parallel loop present(fp[0:nat2])
  for (int i=0;i<natoms;i++)
    fp[i*natoms+i] = fpb[i*natoms+i] = 1.;

  double tansc = alpha;
  if (alpha<=0.) tansc = 10.;

  //#pragma acc parallel loop private(fp)
  for (int n=0;n<gsa;n++)
  {
    double x = grid1[gc*n+0]; double y = grid1[gc*n+1]; double z = grid1[gc*n+2];

    int wa = get_atom_grid(n,natoms,gsa);

    int i = ta;
    for (int j=0;j<natoms;j++)
    if (j!=i)
    //#pragma acc parallel loop collapse(2) present(coords[0:3*natoms],aij[0:nat2],Ra[0:nat2])
    //for (int i=0;i<natoms;i++)
    //for (int j=0;j<i;j++)
    {
      double a12 = aij[i*natoms+j];
      double a21 = aij[j*natoms+i];
      double b12 = bij[i*natoms+j];
      double b21 = bij[j*natoms+i];
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

      if (alpha>0.)
      {
        //f12 = tanhh(tansc*nu21); //note 12/21 swap
        //f21 = tanhh(tansc*nu12);
        f12 = sigmoid(tansc*nu12);
        f21 = sigmoid(tansc*nu21);
        double ch12 = cosh(fabs(nu12));
        double och = 1./ch12;
        f12 = och*och;
        f21 = och*och;
      }

      fp[i*natoms+j] = f12;
      fp[j*natoms+i] = f21;

     //"b" terms

      //const double nu0 = beta*oR12;
      const double nu0 = beta;
      if (alpha>0.)
      {
        //mu12 = r1-r2; //sigmoid not dependent on R12
        omm12 = 1.-mu12*mu12;

        nu12 =  mu12 + b12*omm12 - nu0;
        nu21 = -mu12 + b21*omm12 + nu0;

        f12 = sigmoid(tansc*nu12);
        f21 = sigmoid(tansc*nu21);
        //f12 = tanhh(tansc*nu21); //note swap
        //f21 = tanhh(tansc*nu12);

        //f12 = bf3d(nu12);
        //f21 = bf3d(nu21);
      }
      else
      {
       //standard partition, except b term
        //omm12 = 1.-mu12*mu12;
        nu12 =  mu12 + b12*omm12;
        nu21 = -mu12 + b21*omm12;

        f12 = bf3d(nu12);
        f21 = bf3d(nu21);
      }

      fpb[i*natoms+j] = f12;
      fpb[j*natoms+i] = f21;

    }

   //doing the Becke grid weight
    if (ta==-1)
    {
     //total weight does not normalize to unity
     //  when alpha != 1.0
      double s1 = 1.;
      for (int j=0;j<natoms;j++)
        s1 *= fpb[wa*natoms+j];

      double norm = 1.e-15;
      //#pragma acc loop collapse(2)
      for (int i=0;i<natoms;i++)
      {
        double s0 = 1.;
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
      double s1 = 1.;
      for (int j=0;j<natoms;j++)
        s1 *= fpb[ta*natoms+j]; //ta here, wa above

      double norm = 1.e-15;
      if (0)
      for (int i=0;i<natoms;i++)
      {
        double s0 = 1.;
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

void atomic_domain_cell_wt(const int ta, const float alpha, const float beta, const int natoms, const int gc, const int gs, float* grid1, float* wt1, int* atno, float* coords0)
{
  double* coordsd = new double[3*natoms];
  for (int i=0;i<3*natoms;i++)
    coordsd[i] = coords0[i];

  int gs6 = gc*gs;
  double* grid1d = new double[gs6];
  double* wt1d = new double[gs];

  for (int j=0;j<gs6;j++)
    grid1d[j] = grid1[j];
  for (int j=0;j<gs;j++)
    wt1d[j] = wt1[j];

  atomic_domain_cell_wt(ta,alpha,beta,natoms,gc,gs,grid1d,wt1d,atno,coordsd);

  for (int j=0;j<gs;j++)
    wt1[j] = wt1d[j];

  delete [] coordsd;
  delete [] grid1d;
  delete [] wt1d;

  return;
}

void atomic_domain_cell_wt(const int ta, const float alpha, const float beta, const int natoms, const int gc, const int gs, double* grid1, double* wt1, int* atno, float* coords0)
{
  double* coordsd = new double[3*natoms];
  for (int i=0;i<3*natoms;i++)
    coordsd[i] = coords0[i];

  atomic_domain_cell_wt(ta,alpha,beta,natoms,gc,gs,grid1,wt1,atno,coordsd);

  delete [] coordsd;

  return;
}

//n-center grid wts (float)
void becke_weight_nc(const int natoms, const int gc, const int gs, float* grid1, float* wt1, int* atno, float* coords)
{
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
      float oR12 = 1./Ra[i*natoms+j];

      float dx1 = x-coords[3*i+0]; float dy1 = y-coords[3*i+1]; float dz1 = z-coords[3*i+2];
      float dx2 = x-coords[3*j+0]; float dy2 = y-coords[3*j+1]; float dz2 = z-coords[3*j+2];
      float r1 = sqrtf(dx1*dx1+dy1*dy1+dz1*dz1);
      float r2 = sqrtf(dx2*dx2+dy2*dy2+dz2*dz2);

      float mu12 = (r1-r2)*oR12;
      float omm12 = 1.-mu12*mu12;
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

void becke_weight_nc(const int natoms, const int gc, const int gs, float* grid1, float* wt1, int* atno, double* coords)
{
  float coordsf[3*natoms];
  for (int n=0;n<3*natoms;n++)
    coordsf[n] = coords[n];

  becke_weight_nc(natoms,gc,gs,grid1,wt1,atno,coordsf);

  return;
}

//n-center grid wts (double)
void becke_weight_nc(const int natoms, const int gc, const int gs, double* grid1, double* wt1, int* atno, double* coords)
{
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
    double ra1 = sqrt(dax*dax+day*day+daz*daz);

    Ra[i*natoms+j] = ra1;
    Ra[j*natoms+i] = ra1;
    //printf(" Ra(%i-%i): %8.5f \n",i,j,ra1);
  }


  double aij[nat2];
  for (int i=0;i<natoms;i++)
  for (int j=0;j<i;j++)
  {
    int Z1 = atno[i]; int Z2 = atno[j];
    double a12 = becke_ad(Z1,Z2);
    double a21 = becke_ad(Z2,Z1);
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

 #pragma acc parallel loop independent present(grid1[0:6*gs],wt1[0:gs])
  for (int i=0;i<gs;i++)
  {
    float r1 = grid1[6*i+3];
    float r2 = grid1[6*i+4];
    float r3 = grid1[6*i+5];
    float s1 = get_3c_bw(r1,r2,r3,a12,a13,a23,a21,a31,a32,oR12,oR13,oR23);

    wt1[i] *= s1;
  }

 #pragma acc parallel loop independent present(grid2[0:6*gs],wt2[0:gs])
  for (int i=0;i<gs;i++)
  {
    float r1 = grid2[6*i+3];
    float r2 = grid2[6*i+4];
    float r3 = grid2[6*i+5];
    float s1 = get_3c_bw(r1,r2,r3,a21,a23,a13,a12,a32,a31,oR12,oR23,oR13);

    wt2[i] *= s1;
  }

 #pragma acc parallel loop independent present(grid3[0:6*gs],wt3[0:gs])
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

  //printf("Z1/2: %i %i \n",Z1,Z2);
  //printf("1/R: %8.5f \n",oR);
  //printf(" a1/2: %8.5f %8.5f \n",a1,a2);

 #pragma acc parallel loop independent present(grid1[0:6*gs],wt1[0:gs])
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

  //printf("done with g1 \n");

 #pragma acc parallel loop independent present(grid2[0:6*gs],wt2[0:gs])
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

 #pragma acc parallel loop independent present(grid1[0:6*gs],wt1[0:gs])
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

 #pragma acc parallel loop independent present(grid2[0:6*gs],wt2[0:gs])
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

double becke_ad(int Z1, int Z2)
{
  double x = get_radii_2(Z1)/get_radii_2(Z2);

  double u = (x-1.)/(x+1.);
  double a = u/(u*u-1.);

  if (a>0.5) a = 0.5;
  else if (a<-0.5) a = -0.5;
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

//all-float version
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

    float A0 = coords[3*n+0]; float B0 = coords[3*n+1]; float C0 = coords[3*n+2];
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

void get_becke_grid_full(int natoms, int* atno, double* coords, int nrad, int nang, float* ang_g, float* ang_w, const int gc, float* grid, float* wt)
{
  float coordsf[3*natoms];
  for (int n=0;n<3*natoms;n++)
    coordsf[n] = coords[n];
  get_becke_grid_full(0.,natoms,atno,coordsf,nrad,nang,ang_g,ang_w,gc,grid,wt);

  return;
}

void get_becke_grid_full(int natoms, int* atno, double* coords, int nrad, int nang, double* ang_g, double* ang_w, const int gc, double* grid, double* wt)
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

void get_becke_grid_full(int natoms, int* atno, double* coords, int nrad, int nang, double* ang_g, double* ang_w, const int gc, float* grid, float* wt)
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


void add_r1_to_grid4(int gs, float* grid1, float A2, float B2, float C2)
{
 #pragma acc parallel loop independent present(grid1[0:4*gs])
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

void compute_fxc(bool gbasis, int natoms, int* atno, double* coords, vector<vector<double> > &basis, bool need_wt, int gsa, float* grid, float* wt, double* vxc, double* fxc, int prl)
{
  if (!gbasis)
  { printf("\n ERROR shouldn't be here in compute_fxc (Gaussian basis version) \n"); exit(-1); }

  if (prl>1) printf("  compute_fxc: Gaussian basis \n");

  //int gsa = gs*natoms;
  int gsa6 = gsa*6;

  int N = basis.size();
  int N2 = N*N;
  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);

  int iN = imaxN;
  float** val1 = new float*[iN];
  float** val2 = new float*[iN];
  for (int i=0;i<iN;i++)
    val1[i] = new float[gsa];
  for (int i=0;i<iN;i++)
    val2[i] = new float[gsa];

  float* val0 = new float[gsa];

  float* grid1 = new float[gsa6];
  float* grid2 = new float[gsa6];

  #pragma acc enter data create(grid1[0:gsa6],grid2[0:gsa6])
  #pragma acc enter data create(val0[0:gsa],val1[0:iN][0:gsa],val2[0:iN][0:gsa])

  for (int j=0;j<N2;j++)
    fxc[j] = 0.;

  const int ig = 10; //first index to exp in Gaussian basis

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    float Z1 = (float)atno[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];

    copy_grid(gsa,grid1,grid);
    recenter_grid_zero(gsa,grid1,-A1,-B1,-C1);

    //#pragma acc update self(grid1[0:gsa6])

    #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gsa])
    for (int i1=0;i1<s2-s1;i1++)
    for (int j=0;j<gsa;j++)
      val1[i1][j] = 0.f;

    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int l1 = basis1[1]; int m1 = basis1[2]; int ng = basis1[3];
      int in = ig + ng; //index to find norm

      float* valm = val1[ii1];
      for (int j=0;j<ng;j++)
      {
        double zeta1 = basis1[ig+j]; double norm1 = basis1[in+j];

        eval_gh(gsa,grid1,val0,l1,m1,norm1,zeta1);

       #pragma acc parallel loop present(val0[0:gsa],valm[0:gsa])
        for (int k=0;k<gsa;k++)
          valm[k] += val0[k];
      }
    }

   //compute all
    for (int i1=s1;i1<s2;i1++)
    {
      float* valm = val1[i1-s1];

      for (int i2=s1;i2<s2;i2++)
      {
        float* valn = val1[i2-s1];

        double f1 = 0.;
        if (need_wt)
       #pragma acc parallel loop present(vxc[0:gsa],valm[0:gsa],valn[0:gsa],wt[0:gsa]) reduction(+:f1)
        for (int j=0;j<gsa;j++)
          f1 += vxc[j]*valm[j]*valn[j]*wt[j];
        else
       #pragma acc parallel loop present(vxc[0:gsa],valm[0:gsa],valn[0:gsa]) reduction(+:f1)
        for (int j=0;j<gsa;j++)
          f1 += vxc[j]*valm[j]*valn[j];

        fxc[i1*N+i2] = f1;
      }
    }

   //second atom
    for (int n=m+1;n<natoms;n++)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      float Z2 = (float)atno[m];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];

      #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gsa])
      for (int i2=0;i2<s4-s3;i2++)
      for (int j=0;j<gsa;j++)
        val2[i2][j] = 0.f;

      copy_grid(gsa,grid2,grid);
      recenter_grid_zero(gsa,grid2,-A2,-B2,-C2);

      //#pragma acc update self(grid2[0:gsa6])

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;

        vector<double> basis2 = basis[i2];
        int l2 = basis2[1]; int m2 = basis2[2]; int ng = basis2[3];
        int in = ig + ng; //index to find norm

        float* valn = val2[ii2];
        for (int j=0;j<ng;j++)
        {
          double zeta2 = basis2[ig+j]; double norm2 = basis2[in+j];

          eval_gh(gsa,grid2,val0,l2,m2,norm2,zeta2);

         #pragma acc parallel loop present(val0[0:gsa],valn[0:gsa])
          for (int k=0;k<gsa;k++)
            valn[k] += val0[k];
        }
      }

     //two-atom Pao elements added to grid
      for (int i1=s1;i1<s2;i1++)
      {
        float* valn = val1[i1-s1];

        for (int i2=s3;i2<s4;i2++)
        {
          float* valm = val2[i2-s3];

          double f1 = 0.;
          if (need_wt)
         #pragma acc parallel loop present(vxc[0:gsa],valm[0:gsa],valn[0:gsa],wt[0:gsa]) reduction(+:f1)
          for (int j=0;j<gsa;j++)
            f1 += vxc[j]*valm[j]*valn[j]*wt[j];
          else
         #pragma acc parallel loop present(vxc[0:gsa],valm[0:gsa],valn[0:gsa]) reduction(+:f1)
          for (int j=0;j<gsa;j++)
            f1 += vxc[j]*valm[j]*valn[j];

          fxc[i1*N+i2] = fxc[i2*N+i1] = f1;
        }
      }
    } //loop n over second atom
  }

  #pragma acc update device(fxc[0:N2])

  if (prl>2)
  {
    printf("\n fxc: \n");
    print_square(N,fxc);
  }

  #pragma acc exit data delete(grid1[0:gsa6],grid2[0:gsa6])
  #pragma acc exit data delete(val0[0:gsa],val1[0:iN][0:gsa],val2[0:iN][0:gsa])

  delete [] grid1;
  delete [] grid2;
  delete [] n2i;

  delete [] val0;
  for (int i=0;i<iN;i++)
    delete [] val1[i];
  for (int i=0;i<iN;i++)
    delete [] val2[i];
  delete [] val1;
  delete [] val2;

  return;
}

void compute_fxc(bool gbasis, int natoms, int* atno, double* coords, vector<vector<double> > &basis, bool need_wt, int gsa, float* grid, float* wt, float* vxc, double* fxc, int prl)
{
  double* vxcd = new double[gsa];
  #pragma acc enter data create(vxcd[0:gsa])

  #pragma acc parallel loop present(vxcd[0:gsa],vxc[0:gsa])
  for (int j=0;j<gsa;j++)
    vxcd[j] = vxc[j];

  compute_fxc(gbasis,natoms,atno,coords,basis,need_wt,gsa,grid,wt,vxcd,fxc,prl);

  #pragma acc exit data delete(vxcd[0:gsa])
  delete [] vxcd;

  return;
}

//Gaussian basis
void compute_rho(bool gbasis, int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* Pao, int nrad, int gsa, float* grid, double* rho, double* drho, int prl)
{
  if (!gbasis) return compute_rho(natoms,atno,coords,basis,Pao,nrad,gsa,grid,rho,drho,prl);

  if (prl>1) printf("  compute rho: Gaussian basis \n");

  //int gs = gs*natoms;
  int gsa6 = gsa*6;

  int N = basis.size();
  //int N2 = N*N;
  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);

  int iN = imaxN;
  float** val1 = new float*[iN];
  float** val2 = new float*[iN];
  for (int i=0;i<iN;i++)
    val1[i] = new float[gsa];
  for (int i=0;i<iN;i++)
    val2[i] = new float[gsa];

  float* val0 = new float[gsa];

  float* grid1 = new float[gsa6];
  float* grid2 = new float[gsa6];

  #pragma acc enter data create(grid1[0:gsa6],grid2[0:gsa6])
  #pragma acc enter data create(val0[0:gsa],val1[0:iN][0:gsa],val2[0:iN][0:gsa])

  if (prl>1)
  {
    printf("\n Pao: \n");
    print_square(N,Pao);
  }

  #pragma acc parallel loop present(rho[0:gsa])
  for (int j=0;j<gsa;j++)
    rho[j] = 0.;

  const int ig = 10; //first index to exp in Gaussian basis

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    float Z1 = (float)atno[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];

    copy_grid(gsa,grid1,grid);
    recenter_grid_zero(gsa,grid1,-A1,-B1,-C1);

    //#pragma acc update self(grid1[0:gsa6])

    #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gsa])
    for (int i1=0;i1<s2-s1;i1++)
    for (int j=0;j<gsa;j++)
      val1[i1][j] = 0.f;

    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;
    
      vector<double> basis1 = basis[i1];
      int l1 = basis1[1]; int m1 = basis1[2]; int ng = basis1[3];
      int in = ig + ng; //index to find norm

      //printf("  (1) basis %i has %2i terms \n",i1,ng);
      float* valm = val1[ii1];
      for (int j=0;j<ng;j++)
      {
        double zeta1 = basis1[ig+j]; double norm1 = basis1[in+j];

        //printf("   (1)  evaluating zeta/norm: %8.5f %8.5f  l/m: %i %i \n",zeta1,norm1,l1,m1);
        eval_gh(gsa,grid1,val0,l1,m1,norm1,zeta1);

       #pragma acc parallel loop present(val0[0:gsa],valm[0:gsa])
        for (int k=0;k<gsa;k++)
          valm[k] += val0[k];
      }
    }

    //printf("\n val1[0]: \n");
    //print_vec(gsa,grid1,val1[0]);

   //compute all
    for (int i1=s1;i1<s2;i1++)
    {
      float* valm = val1[i1-s1];

      for (int i2=s1;i2<s2;i2++)
      {
        float* valn = val1[i2-s1];

        float d12 = Pao[i1*N+i2];
       #pragma acc parallel loop present(rho[0:gsa],valm[0:gsa],valn[0:gsa])
        for (int j=0;j<gsa;j++)
          rho[j] += d12*valm[j]*valn[j];
      }
    }

    for (int n=m+1;n<natoms;n++)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      float Z2 = (float)atno[m];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];

      #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gsa])
      for (int i2=0;i2<s4-s3;i2++)
      for (int j=0;j<gsa;j++)
        val2[i2][j] = 0.f;

      copy_grid(gsa,grid2,grid);
      recenter_grid_zero(gsa,grid2,-A2,-B2,-C2);

      //#pragma acc update self(grid2[0:gsa6])

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;

        vector<double> basis2 = basis[i2];
        int l2 = basis2[1]; int m2 = basis2[2]; int ng = basis2[3];
        int in = ig + ng; //index to find norm

        float* valn = val2[ii2];
        for (int j=0;j<ng;j++)
        {
          double zeta2 = basis2[ig+j]; double norm2 = basis2[in+j];

          eval_gh(gsa,grid2,val0,l2,m2,norm2,zeta2);

         #pragma acc parallel loop present(val0[0:gsa],valn[0:gsa])
          for (int k=0;k<gsa;k++)
            valn[k] += val0[k];
        }
      }

      //printf("\n val2[0]: \n");
      //print_vec(gsa,grid2,val2[0]);

     //two-atom Pao elements added to grid
      for (int i1=s1;i1<s2;i1++)
      {
        float* valn = val1[i1-s1];

        for (int i2=s3;i2<s4;i2++)
        {
          float* valm = val2[i2-s3];

          float d1 = 2.f*Pao[i1*N+i2];
          #pragma acc parallel loop present(rho[0:gsa],valn[0:gsa],valm[0:gsa])
          for (int j=0;j<gsa;j++)
            rho[j] += d1*valn[j]*valm[j];
        }
      }
    } //loop n over second atom
  }

  #pragma acc parallel loop present(rho[0:gsa])
  for (int j=0;j<gsa;j++)
    rho[j] *= 2.;

  if (prl>2)
  {
    #pragma acc update self(rho[0:gsa])
    printf("\n density: \n");
    print_vec(gsa,grid,rho);
  }

  #pragma acc exit data delete(grid1[0:gsa6],grid2[0:gsa6])
  #pragma acc exit data delete(val0[0:gsa],val1[0:iN][0:gsa],val2[0:iN][0:gsa])

  delete [] grid1;
  delete [] grid2;
  delete [] n2i;

  delete [] val0;
  for (int i=0;i<iN;i++)
    delete [] val1[i];
  for (int i=0;i<iN;i++)
    delete [] val2[i];
  delete [] val1;
  delete [] val2;

  return;
}

void compute_rho(int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* Pao, int nrad, int gsa, float* grid, double* rho, double* drho, int prl)
{
  //int gsa = natoms*gs;
  int gsa3 = 3*gsa;
  int gsa6 = 6*gsa;

  //printf(" in compute_rho \n");
  int nang = gsa/nrad;
  #include "jsetup.cpp"
  int gs2 = gsa;

  if (sgs_basis && natoms>1)
  { printf("\n ERROR: compute_rho doesn't support natoms>1 SGS basis \n"); exit(-1); }
  if (sgs_basis && drho!=NULL)
  { printf("\n WARNING: compute_rho doesn't support drho for SGS basis \n"); drho = NULL; }


  int N = basis.size();
  int N2 = N*N;
  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);

  double norm[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];
  double* Paon = new double[N2];
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    Paon[i*N+j] = Pao[i*N+j]*norm[i]*norm[j];

  if (prl>1)
  {
    printf("\n Paon: \n");
    print_square(N,Paon);
  }

  int iN = imaxN;
  float** val1 = new float*[iN];
  float** val2 = new float*[iN];
  for (int i=0;i<iN;i++)
    val1[i] = new float[gsa];
  for (int i=0;i<iN;i++)
    val2[i] = new float[gsa];

  float** val1p = new float*[iN];
  float** val2p = new float*[iN];
  for (int i=0;i<iN;i++)
    val1p[i] = new float[gsa3];
  for (int i=0;i<iN;i++)
    val2p[i] = new float[gsa3];

  float* grid1 = new float[gsa6];
  float* grid2 = new float[gsa6];

  float* delt = new float[gsa3];


  #pragma acc enter data create(grid1[0:gsa6],grid2[0:gsa6],delt[0:gsa3])
  #pragma acc enter data create(val1[0:iN][0:gsa],val2[0:iN][0:gsa],val1p[0:iN][0:gsa3],val2p[0:iN][0:gsa3])

  #pragma acc parallel loop present(delt[0:gsa3])
  for (int m=0;m<gsa3;m++)
    delt[m] = 0.f;

  #pragma acc parallel loop present(rho[0:gsa])
  for (int j=0;j<gsa;j++)
    rho[j] = 0.f;
  if (drho!=NULL)
  #pragma acc parallel loop present(drho[0:gsa])
  for (int j=0;j<gsa;j++)
    drho[j] = 0.f;

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    float Z1 = (float)atno[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];

    copy_grid(gsa,grid1,grid);
    recenter_grid_zero(gsa,grid1,-A1,-B1,-C1);

    #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gsa])
    for (int i1=0;i1<s2-s1;i1++)
    for (int j=0;j<gsa;j++)
      val1[i1][j] = 1.f;

    if (drho!=NULL)
    #pragma acc parallel loop collapse(2) present(val1p[0:iN][0:gsa3])
    for (int i1=0;i1<s2-s1;i1++)
    for (int j=0;j<gsa3;j++)
      val1p[i1][j] = 1.f;

    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

      if (ss_basis)
        eval_ss(ii1,gs2,grid1,val1[ii1],n1,l1,m1,zeta1,Rc);
      else if (sgs_basis)
        eval_sgs(ii1,gs1,gs2,grid1,val1[ii1],n1,l1,m1,zeta1,Rc);
      else
        eval_sh(ii1,gsa,grid1,val1[ii1],n1,l1,m1,zeta1);
      if (drho!=NULL)
        eval_p(gsa,grid1,val1p[ii1],n1,l1,m1,zeta1);
    }

   //single-atom Pao elements added to grid
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;
      float* valn = val1[ii1];
      float* valpn = val1p[ii1];

      for (int i2=s1;i2<s2;i2++)
      {
        int ii2 = i2-s1;
        float d1 = Paon[i1*N+i2];
        float* valm = val1[ii2];
        float* valpm = val1p[ii2];

       #pragma acc parallel loop present(rho[0:gsa],valn[0:gsa],valm[0:gsa])
        for (int j=0;j<gsa;j++)
          rho[j] += d1*valn[j]*valm[j];

        if (drho!=NULL)
       #pragma acc parallel loop present(delt[0:gsa3],valn[0:gsa],valpn[0:gsa3],valm[0:gsa],valpm[0:gsa3])
        for (int j=0;j<gsa;j++)
        {
          delt[3*j+0] += d1*(valn[j]*valpm[3*j+0]+valpn[3*j+0]*valm[j]);
          delt[3*j+1] += d1*(valn[j]*valpm[3*j+1]+valpn[3*j+1]*valm[j]);
          delt[3*j+2] += d1*(valn[j]*valpm[3*j+2]+valpn[3*j+2]*valm[j]);
        }
      }
    }

    for (int n=m+1;n<natoms;n++)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      float Z2 = (float)atno[m];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];

      #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gsa])
      for (int i2=0;i2<s4-s3;i2++)
      for (int j=0;j<gsa;j++)
        val2[i2][j] = 1.f;
      if (drho!=NULL)
      #pragma acc parallel loop collapse(2) present(val2p[0:iN][0:gsa3])
      for (int i2=0;i2<s4-s3;i2++)
      for (int j=0;j<gsa3;j++)
        val2p[i2][j] = 1.f;

      copy_grid(gsa,grid2,grid);
      recenter_grid_zero(gsa,grid2,-A2,-B2,-C2);

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;

        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        if (sgs_basis)
        { printf("\n ERROR: shouldn't be here in compute rho for SGS basis \n"); exit(-1); }
        else
          eval_sh(ii2,gsa,grid2,val2[ii2],n2,l2,m2,zeta2);
        if (drho!=NULL)
          eval_p(gsa,grid2,val2p[ii2],n2,l2,m2,zeta2);
      }

     //two-atom Pao elements added to grid
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;
        float* valn = val1[ii1];
        float* valpn = val1p[ii1];

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          float d1 = 2.f*Paon[i1*N+i2];
          float* valm = val2[ii2];
          float* valpm = val2p[ii2];

          #pragma acc parallel loop present(rho[0:gsa],valn[0:gsa],valm[0:gsa])
          for (int j=0;j<gsa;j++)
            rho[j] += d1*valn[j]*valm[j];

          if (drho!=NULL)
          #pragma acc parallel loop present(delt[0:gsa3],valn[0:gsa],valpn[0:gsa3],valm[0:gsa],valpm[0:gsa3])
          for (int j=0;j<gsa;j++)
          {
            delt[3*j+0] += d1*(valn[j]*valpm[3*j+0]+valpn[3*j+0]*valm[j]);
            delt[3*j+1] += d1*(valn[j]*valpm[3*j+1]+valpn[3*j+1]*valm[j]);
            delt[3*j+2] += d1*(valn[j]*valpm[3*j+2]+valpn[3*j+2]*valm[j]);
          }
        }
      }

    } //loop n over unique atoms

  } //loop m over natoms


  if (drho!=NULL)
  #pragma acc parallel loop present(drho[0:gsa],delt[0:gsa3])
  for (int m=0;m<gsa;m++)
    drho[m] = delt[3*m+0]*delt[3*m+0]+delt[3*m+1]*delt[3*m+1]+delt[3*m+2]*delt[3*m+2];


 //cleanup
  #pragma acc exit data delete(grid1[0:gsa6],grid2[0:gsa6],delt[0:gsa3])
  #pragma acc exit data delete(val1[0:iN][0:gsa],val2[0:iN][0:gsa],val1p[0:iN][0:gsa3],val2p[0:iN][0:gsa3])

  delete [] Paon;
  delete [] n2i;

  delete [] delt;

  delete [] grid1;
  delete [] grid2;

  for (int i=0;i<iN;i++)
    delete [] val1[i];
  for (int i=0;i<iN;i++)
    delete [] val2[i];
  delete [] val1;
  delete [] val2;

  for (int i=0;i<iN;i++)
    delete [] val1p[i];
  for (int i=0;i<iN;i++)
    delete [] val2p[i];
  delete [] val1p;
  delete [] val2p;

  return;
}

//evaluates in double precision
void compute_rhod(int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* Pao, int nrad, int gsa, double* grid, double* rho, double* drho, double* Td, int prl)
{
  int gsa3 = 3*gsa;
  int gsa6 = 6*gsa;

  int tid = -1; //no multi-gpu parallel

  int nang = gsa/nrad;
  #include "jsetup.cpp"
  int gs2 = gsa;
  //if (sgs_basis) printf("  compute_rhod gs12: %4i %4i \n",gs1,gs2);

  if (sgs_basis && natoms>1){ printf("\n ERROR: compute_rho doesn't support natoms>1 for SGS basis \n"); exit(-1); }
  if (sgs_basis && drho!=NULL) printf("\n WARNING: compute_rho for SGS basis does not support drho \n");

  //if (Td!=NULL && natoms>1) printf("\n WARNING: testing Td in compute_rhod for natoms>1 \n");

  int N = basis.size();
  int N2 = N*N;
  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);

  double norm[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];
  double* Paon = new double[N2];
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    Paon[i*N+j] = Pao[i*N+j]*norm[i]*norm[j];

  if (prl>1)
  {
    printf("\n Paon: \n");
    print_square(N,Paon);
  }

  int iN = imaxN;
  double** val1 = new double*[iN];
  double** val2 = new double*[iN];
  for (int i=0;i<iN;i++)
    val1[i] = new double[gsa];
  for (int i=0;i<iN;i++)
    val2[i] = new double[gsa];

  double** val1p = new double*[iN];
  double** val2p = new double*[iN];
  for (int i=0;i<iN;i++)
    val1p[i] = new double[gsa3];
  for (int i=0;i<iN;i++)
    val2p[i] = new double[gsa3];

  double* grid1 = new double[gsa6];
  double* grid2 = new double[gsa6];

  double* tmp = new double[gsa3];
  double* delt = new double[gsa3];


  #pragma acc enter data create(grid1[0:gsa6],grid2[0:gsa6],tmp[0:gsa3],delt[0:gsa3])
  #pragma acc enter data create(val1[0:iN][0:gsa],val2[0:iN][0:gsa],val1p[0:iN][0:gsa3],val2p[0:iN][0:gsa3])

  #pragma acc parallel loop present(delt[0:gsa3])
  for (int m=0;m<gsa3;m++)
    delt[m] = 0.;

  #pragma acc parallel loop present(rho[0:gsa])
  for (int j=0;j<gsa;j++)
    rho[j] = 0.;
  if (drho!=NULL)
  #pragma acc parallel loop present(drho[0:gsa])
  for (int j=0;j<gsa;j++)
    drho[j] = 0.;

  if (Td!=NULL)
  #pragma acc parallel loop present(Td[0:gsa])
  for (int j=0;j<gsa;j++)
    Td[j] = 0.;

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    double Z1 = (double)atno[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];

    copy_grid(gsa,grid1,grid);
    recenter_grid_zero(gsa,grid1,-A1,-B1,-C1);

    #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gsa])
    for (int i1=0;i1<s2-s1;i1++)
    for (int j=0;j<gsa;j++)
      val1[i1][j] = 1.;

    if (drho!=NULL)
    #pragma acc parallel loop collapse(2) present(val1p[0:iN][0:gsa3])
    for (int i1=0;i1<s2-s1;i1++)
    for (int j=0;j<gsa3;j++)
      val1p[i1][j] = 1.;

    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

      if (ss_basis)
        eval_ssd(ii1,gs2,grid1,val1[ii1],n1,l1,m1,zeta1,Rc);
      else if (sgs_basis)
        eval_sgsd(ii1,gs1,gs2,grid1,val1[ii1],n1,l1,m1,zeta1,Rc);
      else
        eval_shd(ii1,gsa,grid1,val1[ii1],n1,l1,m1,zeta1);
      if (drho!=NULL)
      {
        if (ss_basis)
          eval_ss_pd(gs2,grid1,val1p[ii1],tmp,n1,l1,m1,zeta1,Rc);
        else
          eval_pd(tid,gsa,grid1,val1p[ii1],n1,l1,m1,zeta1);
      }
    }

   //single-atom Pao elements added to grid
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;
      double* valn = val1[ii1];
      double* valpn = val1p[ii1];

      for (int i2=s1;i2<s2;i2++)
      {
        int ii2 = i2-s1;
        double d1 = Paon[i1*N+i2];
        double* valm = val1[ii2];
        double* valpm = val1p[ii2];

       #pragma acc parallel loop present(rho[0:gsa],valn[0:gsa],valm[0:gsa])
        for (int j=0;j<gsa;j++)
          rho[j] += d1*valn[j]*valm[j];

        if (drho!=NULL)
       #pragma acc parallel loop present(delt[0:gsa3],valn[0:gsa],valpn[0:gsa3],valm[0:gsa],valpm[0:gsa3])
        for (int j=0;j<gsa;j++)
        {
          delt[3*j+0] += d1*(valn[j]*valpm[3*j+0]+valpn[3*j+0]*valm[j]);
          delt[3*j+1] += d1*(valn[j]*valpm[3*j+1]+valpn[3*j+1]*valm[j]);
          delt[3*j+2] += d1*(valn[j]*valpm[3*j+2]+valpn[3*j+2]*valm[j]);
        }

        d1 *= 0.5;
        if (Td!=NULL)
       #pragma acc parallel loop present(Td[0:gsa],valpn[0:gsa3],valpm[0:gsa3])
        for (int j=0;j<gsa;j++)
        {
          double vx = valpn[3*j+0]*valpm[3*j+0]; double vy = valpn[3*j+1]*valpm[3*j+1]; double vz = valpn[3*j+2]*valpm[3*j+2];
          Td[j] += d1*(vx + vy + vz);
        }
      }
    }

    for (int n=m+1;n<natoms;n++)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      double Z2 = (double)atno[m];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];

      #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gsa])
      for (int i2=0;i2<s4-s3;i2++)
      for (int j=0;j<gsa;j++)
        val2[i2][j] = 1.;
      #pragma acc parallel loop collapse(2) present(val2p[0:iN][0:gsa3])
      for (int i2=0;i2<s4-s3;i2++)
      for (int j=0;j<gsa3;j++)
        val2p[i2][j] = 1.;

      copy_grid(gsa,grid2,grid);
      recenter_grid_zero(gsa,grid2,-A2,-B2,-C2);

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;

        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        eval_shd(ii2,gsa,grid2,val2[ii2],n2,l2,m2,zeta2);
        if (drho!=NULL)
          eval_pd(tid,gsa,grid2,val2p[ii2],n2,l2,m2,zeta2);
      }

     //two-atom Pao elements added to grid
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;
        double* valn = val1[ii1];
        double* valpn = val1p[ii1];

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          double d1 = 2.*Paon[i1*N+i2];
          double* valm = val2[ii2];
          double* valpm = val2p[ii2];

          #pragma acc parallel loop present(rho[0:gsa],valn[0:gsa],valm[0:gsa])
          for (int j=0;j<gsa;j++)
            rho[j] += d1*valn[j]*valm[j];

          if (drho!=NULL)
          #pragma acc parallel loop present(delt[0:gsa3],valn[0:gsa],valpn[0:gsa3],valm[0:gsa],valpm[0:gsa3])
          for (int j=0;j<gsa;j++)
          {
            delt[3*j+0] += d1*(valn[j]*valpm[3*j+0]+valpn[3*j+0]*valm[j]);
            delt[3*j+1] += d1*(valn[j]*valpm[3*j+1]+valpn[3*j+1]*valm[j]);
            delt[3*j+2] += d1*(valn[j]*valpm[3*j+2]+valpn[3*j+2]*valm[j]);
          }

          d1 *= 0.5; //CPMZ not sure
          if (0)
          if (Td!=NULL)
         #pragma acc parallel loop present(Td[0:gsa],valpn[0:gsa3],valpm[0:gsa3])
          for (int j=0;j<gsa;j++)
          {
            double vx = valpn[3*j+0]*valpm[3*j+0]; double vy = valpn[3*j+1]*valpm[3*j+1]; double vz = valpn[3*j+2]*valpm[3*j+2];
            Td[j] += d1*(vx + vy + vz);
          }

        }
      }

    } //loop n over unique atoms

  } //loop m over natoms


  if (drho!=NULL)
  #pragma acc parallel loop present(drho[0:gsa],delt[0:gsa3])
  for (int m=0;m<gsa;m++)
    drho[m] = delt[3*m+0]*delt[3*m+0]+delt[3*m+1]*delt[3*m+1]+delt[3*m+2]*delt[3*m+2];


 //cleanup
  #pragma acc exit data delete(grid1[0:gsa6],grid2[0:gsa6],tmp[0:gsa3],delt[0:gsa3])
  #pragma acc exit data delete(val1[0:iN][0:gsa],val2[0:iN][0:gsa],val1p[0:iN][0:gsa3],val2p[0:iN][0:gsa3])

  delete [] Paon;
  delete [] n2i;

  delete [] tmp;
  delete [] delt;

  delete [] grid1;
  delete [] grid2;

  for (int i=0;i<iN;i++)
    delete [] val1[i];
  for (int i=0;i<iN;i++)
    delete [] val2[i];
  delete [] val1;
  delete [] val2;

  for (int i=0;i<iN;i++)
    delete [] val1p[i];
  for (int i=0;i<iN;i++)
    delete [] val2p[i];
  delete [] val1p;
  delete [] val2p;

  return;
}

void compute_rhod(int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* Pao, int nrad, int gsa, float* gridf, double* rho, double* drho, int prl)
{
  int gsa6 = 6*gsa;

  double* grid = new double[gsa6];
  #pragma acc enter data create(grid[0:gsa6])

  #pragma acc parallel loop present(grid[0:gsa6],gridf[0:gsa6])
  for (int j=0;j<gsa6;j++)
    grid[j] = gridf[j];

  #pragma acc update self(grid[0:gsa6])
  compute_rhod(natoms,atno,coords,basis,Pao,nrad,gsa,grid,rho,drho,NULL,prl);

  #pragma acc exit data delete(grid[0:gsa6])
  delete [] grid;

  return;
}


void compute_delta(int natoms, int* atno, double* coords, vector<vector<double> > &basis1, vector<vector<double> > &basis2, int No, double* jCA, bool gga, bool tau, float* rho, int gsa, float* grid, float* wt, double* diff, int prl)
{
 //basis1 --> new AO basis
 //basis2 --> MO basis

  if (prl>1) printf("\n compute_delta \n");

 //expects rho to be normalized to one
  int gs = gsa/natoms;
  int gsa6 = 6*gsa;

  int N1 = basis1.size();
  int N12 = N1*N1;

  if (No>N1) { printf("\n ERROR: No cannot exceed basis1.size() in compute_delta \n"); return; }

  int* n2i1 = new int[natoms];
  int imaxN1 = get_imax_n2i(natoms,N1,basis1,n2i1);

  int N2 = basis2.size();

  int* n2i2 = new int[natoms];
  int imaxN2 = get_imax_n2i(natoms,N2,basis2,n2i2);

  float* grid1 = new float[gsa6];

  int iN1 = N1;
  int iN2 = N2; //working with MOs, so need all AOs at once
  float** val1 = new float*[iN1];
  float** val2 = new float*[iN2];
  for (int i=0;i<iN1;i++)
    val1[i] = new float[gsa];
  for (int i=0;i<iN2;i++)
    val2[i] = new float[gsa];

  float** val3 = new float*[No];
  for (int i=0;i<No;i++)
    val3[i] = new float[gsa];

  float* rhoc = new float[gsa];
  float* rhow = new float[gsa];
  float* phichi = new float[gsa];

  #pragma acc enter data create(grid1[0:gsa6])
  #pragma acc enter data create(val1[0:iN1][0:gsa],val2[0:iN2][0:gsa])
  #pragma acc enter data create(val3[0:No][0:gsa],phichi[0:gsa],rhow[0:gsa],rhoc[0:gsa])

  for (int j=0;j<N12;j++)
    diff[j] = 0.f;

  #pragma acc parallel loop present(rhow[0:gsa],rho[0:gsa],wt[0:gsa])
  for (int j=0;j<gsa;j++)
    rhow[j] = rho[j]*wt[j];

 //construct AO quantities for basis1
  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i1[m-1]; int s2 = n2i1[m];

    float Z1 = (float)atno[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];

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

      if (tau)
      {
        //#pragma acc parallel loop present(val1[0:iN1][0:gsa],val1T[0:iN1][0:gsa])
        //for (int j=0;j<gsa;j++)
        //  val1T[i1][j] = val1[i1][j];
        eval_ke(gsa,grid1,val1[i1],n1,l1,zeta1);
      }

    } //loop i1

  } //loop m over natoms


 //construct AO quantities for basis2
  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i2[m-1]; int s2 = n2i2[m];

    float Z1 = (float)atno[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];

    copy_grid(gsa,grid1,grid);
    recenter_grid_zero(gsa,grid1,-A1,-B1,-C1);

    #pragma acc parallel loop collapse(2) present(val2[0:iN2][0:gsa])
    for (int i1=s1;i1<s2;i1++)
    for (int j=0;j<gsa;j++)
      val2[i1][j] = 1.f;

    for (int i1=s1;i1<s2;i1++)
    {
      vector<double> basisn = basis2[i1];
      int n1 = basisn[0]; int l1 = basisn[1]; int m1 = basisn[2]; double zeta1 = basisn[3];
      float nm1 = basisn[4];

      eval_sh(i1,gsa,grid1,val2[i1],n1,l1,m1,zeta1);
      #pragma acc parallel loop present(val2[0:iN2][0:gsa])
      for (int j=0;j<gsa;j++)
        val2[i1][j] *= nm1;

    } //loop i1

  } //loop m over natoms


 //get the MOs
  #pragma acc parallel loop collapse(2) present(val2[0:iN2][0:gsa],val3[0:No][0:gsa],jCA[0:N2])
  for (int i3=0;i3<No;i3++)
  for (int j=0;j<gsa;j++)
  {
    float v2 = 0.f;
    #pragma acc loop reduction(+:v2)
    for (int k=0;k<iN2;k++)
      v2 += jCA[k*N2+i3]*val2[k][j];
    val3[i3][j] = v2;
  }


 //here compute diff
  int No1 = No;
  for (int i1=0;i1<iN1;i1++)
  for (int i3=0;i3<No1;i3++)
  {
    float* valm = val1[i1]; //AO 
    float* valn = val3[i3]; //MO

    int at1 = basis1[i1][9];
    int sh = gs*at1;

    #pragma acc parallel loop present(phichi[0:gsa])
    for (int j=0;j<gsa;j++)
      phichi[j] = 0.f;

   //overlap of AO (or Laplacian of AO) with MO
    float nm1 = 0.f;
    #pragma acc parallel loop present(phichi[0:gsa],valm[0:gsa],valn[0:gsa],wt[0:gsa]) reduction(+:nm1)
    for (int j=0;j<gs;j++)
    {
      float v1 = valm[sh+j]*valn[sh+j];
      phichi[sh+j] = v1;
      nm1 += v1*wt[sh+j];
    }
    if (prl>1)
      printf("   AO[%i] - MO[%i] norm: %8.5f \n",i1,i3,nm1);

    if (fabs(nm1)<1.e-3f)
      nm1 = 1.f;
    else
      nm1 = 1.f/nm1;

    #pragma acc parallel loop present(phichi[0:gsa])
    for (int j=0;j<gsa;j++)
      phichi[j] *= nm1;

    if (nm1!=1.f && 0)
    {
      vector<double> basisn = basis1[i1];
      int n1 = basisn[0]; int l1 = basisn[1]; int m1 = basisn[2]; double zeta1 = basisn[3];
      //debug
      #pragma acc update self(phichi[0:gsa])
      printf("   phi_chi[%2i / %i %i %i %8.5f]: \n",i1,n1,l1,m1,zeta1);
      print_vec(gsa,grid,rho,phichi);
    }

    #pragma acc parallel loop present(phichi[0:gsa],wt[0:gsa])
    for (int j=0;j<gsa;j++)
      phichi[j] *= wt[j];

    float valt = 0.f;
    if (nm1!=1.f)
    {
     //deal with the atomic component
      #pragma acc parallel loop present(rhoc[0:gsa])
      for (int j=0;j<gsa;j++)
        rhoc[j] = 0.f;
      #pragma acc parallel loop present(rhoc[0:gsa],rhow[0:gsa])
      for (int j=0;j<gs;j++)
        rhoc[sh+j] = rhow[sh+j];

      float valt1 = 0.f;
      #pragma acc parallel loop present(rhoc[0:gsa],phichi[0:gsa]) reduction(+:valt1)
      for (int j=0;j<gsa;j++)
        valt1 += fabs(phichi[j]-rhoc[j]);
      float valt2 = 0.f;
      #pragma acc parallel loop present(rhoc[0:gsa],phichi[0:gsa]) reduction(+:valt2)
      for (int j=0;j<gsa;j++)
        valt2 += fabs(phichi[j]+rhoc[j]);
      valt = min(valt1,valt2);
    }

    diff[i1*N1+i3] = valt;
  }


  //#pragma acc update device(diff[0:N2])

 //cleanup
  #pragma acc exit data delete(grid1[0:gsa6])
  #pragma acc exit data delete(val1[0:iN1][0:gsa],val2[0:iN2][0:gsa])
  #pragma acc exit data delete(val3[0:No][0:gsa],phichi[0:gsa],rhow[0:gsa],rhoc[0:gsa])

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
  delete [] phichi;
  delete [] rhow;
  delete [] rhoc;

  return;
}

void density_in_basis2(int natoms, int* atno, double* coords, vector<vector<double> > &basis1, vector<vector<double> > &basis2, int No, double* jCA, int gsa, float* grid, float* wt, double* Paom, int prl)
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
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];

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
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];

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

void batomic_charges(float alpha, int natoms, int* atno, double* coords, int nrad, int nang, float* ang_g, float* ang_w, float* chg, double* rho, int prl)
{
  int gc = 6;
  int gs = nrad*nang;
  int gsa = natoms*gs;

  float coordsf[3*natoms];
  for (int n=0;n<3*natoms;n++)
    coordsf[n] = coords[n];

  float* grid = new float[gc*gsa];
  float* wt = new float[gsa];

  #pragma acc enter data create(grid[0:6*gsa],wt[0:gsa])

  get_becke_grid_full(alpha,natoms,atno,coordsf,nrad,nang,ang_g,ang_w,gc,grid,wt);

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

void compute_fxc(int natoms, int* atno, double* coords, vector<vector<double> > &basis, bool gga, bool tau, bool need_wt, double* Pao, double* vxc, double* vxcs, int nrad, int gs, float* grid, float* wt, double* fxc, int prl)
{
 //need_wt==0 --> wt vxc
 //need_wt==1 --> expects vxc to be wt'd already

  if (gga && Pao==NULL)
  {
    printf("\n ERROR: gga functionals require Pao in compute_fxc \n");
    exit(1);
  }
  if ((gga || tau) && need_wt)
  {
    printf(" ERROR: compute_fxc mismatch. gga/tau not compatible with need_wt \n");
    exit(1);
  }

  int nang = gs/nrad;
  #include "jsetup.cpp"
  int gs2 = gs;
  if (!sgs_basis) gs1 = 0;
  //printf("  compute_fxc gs12: %4i %4i  ss_basis: %i \n",gs1,gs2,(int)ss_basis);

  if (sgs_basis && natoms>1)
  { printf("\n ERROR: compute_fxc doesn't support natoms>1 for SGS basis \n"); exit(-1); }

  int gs3 = 3*gs;
  int gs6 = 6*gs;

  int N = basis.size();
  int N2 = N*N;
  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);

  float* grid1 = new float[gs6];
  float* grid2 = new float[gs6];

  double norm[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];

  double* Paon = NULL;
  if (gga)
  {
   //norm to 1e-?
    Paon = new double[N2];
    for (int i=0;i<N;i++)
    for (int j=0;j<N;j++)
      Paon[i*N+j] = Pao[i*N+j]*norm[i]*norm[j];
  }

  int iN = imaxN;
  float** val1 = new float*[iN];
  float** val2 = new float*[iN];
  for (int i=0;i<iN;i++)
    val1[i] = new float[gs];
  for (int i=0;i<iN;i++)
    val2[i] = new float[gs];


  double* grho = NULL;
  float** val1p = NULL;
  float** val2p = NULL;
  if (gga)
  {
    grho = new double[gs3];

    val1p = new float*[iN];
    val2p = new float*[iN];
    for (int i=0;i<iN;i++)
      val1p[i] = new float[gs3];
    for (int i=0;i<iN;i++)
      val2p[i] = new float[gs3];
  }

  float** val1T = NULL;
  if (tau)
  {
    val1T = new float*[iN];
    for (int i=0;i<iN;i++)
      val1T[i] = new float[gs];
  }

  #pragma acc enter data create(grid1[0:gs6],grid2[0:gs6])
  #pragma acc enter data create(val1[0:iN][0:gs],val2[0:iN][0:gs])
  if (gga)
  {
    #pragma acc enter data create(grho[0:gs3],val1p[0:iN][0:gs3],val2p[0:iN][0:gs3])
  }
  if (tau)
  {
    #pragma acc enter data create(val1T[0:iN][0:gs])
  }

  for (int j=0;j<N2;j++)
    fxc[j] = 0.f;

  if (gga)
  #pragma acc parallel loop present(grho[0:gs3])
  for (int j=0;j<gs3;j++)
    grho[j] = 0.;

  //first assemble grho
  if (gga)
  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    float Z1 = (float)atno[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];

    copy_grid(gs,grid1,grid);
    recenter_grid_zero(gs,grid1,-A1,-B1,-C1);

    #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs])
    for (int i1=0;i1<s2-s1;i1++)
    for (int j=0;j<gs;j++)
      val1[i1][j] = 1.f;

    if (gga)
    #pragma acc parallel loop collapse(2) present(val1p[0:iN][0:gs3])
    for (int i1=0;i1<s2-s1;i1++)
    for (int j=0;j<gs3;j++)
      val1p[i1][j] = 1.f;

    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

      eval_sh(ii1,gs,grid1,val1[ii1],n1,l1,m1,zeta1);
      eval_p(gs,grid1,val1p[ii1],n1,l1,m1,zeta1);
    }

    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;
      float* valn = val1[ii1];
      float* valpn = val1p[ii1];

      for (int i2=s1;i2<s2;i2++)
      {
        int ii2 = i2-s1;
        float* valm = val1[ii2];
        float* valpm = val1p[ii2];

        double d1 = Paon[i1*N+i2];

       #pragma acc parallel loop present(grho[0:gs3],valn[0:gs],valpn[0:gs3],valm[0:gs],valpm[0:gs3])
        for (int j=0;j<gs;j++)
        {
          grho[3*j+0] += d1*(valn[j]*valpm[3*j+0]+valpn[3*j+0]*valm[j]);
          grho[3*j+1] += d1*(valn[j]*valpm[3*j+1]+valpn[3*j+1]*valm[j]);
          grho[3*j+2] += d1*(valn[j]*valpm[3*j+2]+valpn[3*j+2]*valm[j]);
        }
      }
    }

    for (int n=m+1;n<natoms;n++)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      float Z2 = (float)atno[m];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];

      #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs])
      for (int i2=0;i2<s4-s3;i2++)
      for (int j=0;j<gs;j++)
        val2[i2][j] = 1.f;

      #pragma acc parallel loop collapse(2) present(val2p[0:iN][0:gs3])
      for (int i2=0;i2<s4-s3;i2++)
      for (int j=0;j<gs3;j++)
        val2p[i2][j] = 1.f;

      copy_grid(gs,grid2,grid);
      recenter_grid_zero(gs,grid2,-A2,-B2,-C2);

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;

        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        eval_sh(ii2,gs,grid2,val2[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid2,val2p[ii2],n2,l2,m2,zeta2);
      }

      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;
        float* valn = val1[ii1];
        float* valpn = val1p[ii1];

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          float* valm = val2[ii2];
          float* valpm = val2p[ii2];

          double d1 = Paon[i1*N+i2];

         #pragma acc parallel loop present(grho[0:gs3],valn[0:gs],valpn[0:gs3],valm[0:gs],valpm[0:gs3])
          for (int j=0;j<gs;j++)
          {
            grho[3*j+0] += d1*(valn[j]*valpm[3*j+0]+valpn[3*j+0]*valm[j]);
            grho[3*j+1] += d1*(valn[j]*valpm[3*j+1]+valpn[3*j+1]*valm[j]);
            grho[3*j+2] += d1*(valn[j]*valpm[3*j+2]+valpn[3*j+2]*valm[j]);
          }
        } //loop i2
      } //loop i1

    } //loop n
  } //loop m


  if (gga && prl>-1 && gs<1000)
  {
    #pragma acc update self(grho[0:gs3])
    printf(" grho: \n");
    for (int j=0;j<gs;j++)
      printf("   %8.5f %8.5f %8.5f \n",grho[3*j+0],grho[3*j+1],grho[3*j+2]);
    printf("\n");
  }

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    float Z1 = (float)atno[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];

    copy_grid(gs,grid1,grid);
    recenter_grid_zero(gs,grid1,-A1,-B1,-C1);

    #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs])
    for (int i1=0;i1<s2-s1;i1++)
    for (int j=0;j<gs;j++)
      val1[i1][j] = 1.f;

    if (gga)
    #pragma acc parallel loop collapse(2) present(val1p[0:iN][0:gs3])
    for (int i1=0;i1<s2-s1;i1++)
    for (int j=0;j<gs3;j++)
      val1p[i1][j] = 1.f;

    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

      if (ss_basis)
        eval_ss(ii1,gs2,grid1,val1[ii1],n1,l1,m1,zeta1,Rc);
      else if (sgs_basis)
        eval_sgs(ii1,gs1,gs2,grid1,val1[ii1],n1,l1,m1,zeta1,Rc);
      else
        eval_sh(ii1,gs,grid1,val1[ii1],n1,l1,m1,zeta1);
      if (gga)
        eval_p(gs,grid1,val1p[ii1],n1,l1,m1,zeta1);
    }

   //single-atom elements over grid
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;
      float* valn = val1[ii1];
      float* valpn = NULL; if (gga) valpn = val1p[ii1];
      //float* valtn = NULL; if (tau) valtn = val1T[ii1];

      for (int i2=s1;i2<s2;i2++)
      {
        int ii2 = i2-s1;
        float* valm = val1[ii2];
        float* valpm = NULL; if (gga) valpm = val1p[ii2];

        double valt = 0.;
        if (need_wt)
        {
          #pragma acc parallel loop present(vxc[0:gs],valm[0:gs],valn[0:gs],wt[0:gs]) reduction(+:valt)
          for (int j=0;j<gs;j++)
            valt += valn[j]*valm[j]*vxc[j]*wt[j];
        }
        else
        {
          #pragma acc parallel loop present(vxc[0:gs],valm[0:gs],valn[0:gs]) reduction(+:valt)
          for (int j=0;j<gs;j++)
            valt += valn[j]*valm[j]*vxc[j];
        }

        if (gga)
        #pragma acc parallel loop present(grho[0:gs3],vxcs[0:gs],valm[0:gs],valpm[0:gs3],valn[0:gs],valpn[0:gs3]) reduction(+:valt)
        for (int j=0;j<gs;j++)
        {
          double grx = grho[3*j+0]; double gry = grho[3*j+1]; double grz = grho[3*j+2];
          double valx = valn[j]*valpm[3*j+0]+valpn[3*j+0]*valm[j];
          double valy = valn[j]*valpm[3*j+1]+valpn[3*j+1]*valm[j];
          double valz = valn[j]*valpm[3*j+2]+valpn[3*j+2]*valm[j];
          valt += vxcs[j]*(valx*grx+valy*gry+valz*grz);
        }

       #if 0
        if (tau) //incomplete
        {
          valt = 0.;
          #pragma acc parallel loop present(vxc[0:gs],valtn[0:gs],valm[0:gs]) reduction(+:valt)
          for (int j=0;j<gs;j++)
            valt += valtn[j]*valm[j]*vxc[j];
        }
       #endif

        fxc[i1*N+i2] = valt;
      }
    }

    for (int n=m+1;n<natoms;n++)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      float Z2 = (float)atno[m];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];

      #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs])
      for (int i2=0;i2<s4-s3;i2++)
      for (int j=0;j<gs;j++)
        val2[i2][j] = 1.f;

      if (gga)
      #pragma acc parallel loop collapse(2) present(val2p[0:iN][0:gs3])
      for (int i2=0;i2<s4-s3;i2++)
      for (int j=0;j<gs3;j++)
        val2p[i2][j] = 1.f;

      copy_grid(gs,grid2,grid);
      recenter_grid_zero(gs,grid2,-A2,-B2,-C2);

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;

        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        eval_sh(ii2,gs,grid2,val2[ii2],n2,l2,m2,zeta2);
        if (gga)
          eval_p(gs,grid2,val2p[ii2],n2,l2,m2,zeta2);
      }

     //two-atom elements over grid
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;
        float* valn = val1[ii1];
        float* valpn = NULL; if (gga) valpn = val1p[ii1];
        //float* valtn = NULL; if (tau) valtn = val1T[ii1];

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          float* valm = val2[ii2];
          float* valpm = NULL; if (gga) valpm = val2p[ii2];
          //float* valtm = NULL; if (tau) valtm = val2T[ii2];

          double valt = 0.;
          if (need_wt)
          {
            #pragma acc parallel loop present(vxc[0:gs],valm[0:gs],valn[0:gs],wt[0:gs]) reduction(+:valt)
            for (int j=0;j<gs;j++)
              valt += valn[j]*valm[j]*vxc[j]*wt[j];
          }
          else
          {
            #pragma acc parallel loop present(vxc[0:gs],valm[0:gs],valn[0:gs]) reduction(+:valt)
            for (int j=0;j<gs;j++)
              valt += valn[j]*valm[j]*vxc[j];
          }

          if (gga)
          #pragma acc parallel loop present(grho[0:gs3],vxcs[0:gs],valm[0:gs],valpm[0:gs3],valn[0:gs],valpn[0:gs3]) reduction(+:valt)
          for (int j=0;j<gs;j++)
          {
            double grx = grho[3*j+0]; double gry = grho[3*j+1]; double grz = grho[3*j+2];
            double valx = valn[j]*valpm[3*j+0]+valpn[3*j+0]*valm[j];
            double valy = valn[j]*valpm[3*j+1]+valpn[3*j+1]*valm[j];
            double valz = valn[j]*valpm[3*j+2]+valpn[3*j+2]*valm[j];
            valt += vxcs[j]*(valx*grx+valy*gry+valz*grz);
          }

         #if 0
          if (tau) //incomplete
          {
            valt = 0.;
            #pragma acc parallel loop present(vxc[0:gs],valtn[0:gs],valm[0:gs]) reduction(+:valt)
            for (int j=0;j<gs;j++)
              valt += valtn[j]*valm[j]*vxc[j];
          }
         #endif

          fxc[i1*N+i2] = fxc[i2*N+i1] = valt;

        } //loop i2
      } //loop i1

    } //loop n over unique atoms
  } //loop m over natoms

  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    fxc[i*N+j] *= norm[i]*norm[j];

  #pragma acc update device(fxc[0:N2])


 //cleanup
  #pragma acc exit data delete(grid1[0:gs6],grid2[0:gs6])
  #pragma acc exit data delete(val1[0:iN][0:gs],val2[0:iN][0:gs])
  if (gga)
  {
    #pragma acc exit data delete(grho[0:gs3],val1p[0:iN][0:gs3],val2p[0:iN][0:gs3])
  }
  if (tau)
  {
    #pragma acc exit data delete(val1T[0:iN][0:gs])
  }

  delete [] n2i;

  delete [] grid1;
  delete [] grid2;

  for (int i=0;i<iN;i++)
    delete [] val1[i];
  for (int i=0;i<iN;i++)
    delete [] val2[i];
  delete [] val1;
  delete [] val2;

  if (gga)
  {
    delete [] Paon;
    delete [] grho;
    for (int i=0;i<iN;i++)
      delete [] val1p[i];
    for (int i=0;i<iN;i++)
      delete [] val2p[i];
    delete [] val1p;
    delete [] val2p;
  }

  if (tau)
  {
    for (int i=0;i<iN;i++)
      delete [] val1T[i];
    delete [] val1T;
  }

  return;
}

void compute_fxcd(int natoms, int* atno, double* coords, vector<vector<double> > &basis, bool gga, bool tau, bool need_wt, double* Pao, double* vxc, double* vxcs, int nrad, int gs, double* grid, double* wt, double* fxc, int prl)
{
 //need_wt==0 --> wt vxc
 //need_wt==1 --> expects vxc to be wt'd already

  if (gga && Pao==NULL)
  {
    printf("\n ERROR: gga functionals require Pao in compute_fxc \n");
    exit(1);
  }
  if ((gga || tau) && need_wt)
  {
    printf(" ERROR: compute_fxc mismatch. gga/tau not compatible with need_wt \n");
    exit(1);
  }

  int tid = -1; //no multi-gpu parallel

  int nang = gs/nrad;
  #include "jsetup.cpp"
  int gs2 = gs;
  //printf("  compute_fxcd gs12: %4i %4i  ss_basis: %i \n",gs1,gs2,(int)ss_basis);

  if (sgs_basis && natoms>1)
  { printf("\n ERROR: compute_fxcd doesn't support natoms>1 or drho for sgs_basis \n"); exit(-1); }

  int gs3 = 3*gs;
  int gs6 = 6*gs;

  int N = basis.size();
  int N2 = N*N;
  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);

  double* grid1 = new double[gs6];
  double* grid2 = new double[gs6];

  double norm[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];

  double* Paon = NULL;
  if (gga)
  {
   //norm to 1e-?
    Paon = new double[N2];
    for (int i=0;i<N;i++)
    for (int j=0;j<N;j++)
      Paon[i*N+j] = Pao[i*N+j]*norm[i]*norm[j];
  }

  int iN = imaxN;
  double** val1 = new double*[iN];
  double** val2 = new double*[iN];
  for (int i=0;i<iN;i++)
    val1[i] = new double[gs];
  for (int i=0;i<iN;i++)
    val2[i] = new double[gs];


  double* grho = NULL;
  double** val1p = NULL;
  double** val2p = NULL;
  if (gga)
  {
    grho = new double[gs3];

    val1p = new double*[iN];
    val2p = new double*[iN];
    for (int i=0;i<iN;i++)
      val1p[i] = new double[gs3];
    for (int i=0;i<iN;i++)
      val2p[i] = new double[gs3];
  }

  double** val1T = NULL;
  if (tau)
  {
    val1T = new double*[iN];
    for (int i=0;i<iN;i++)
      val1T[i] = new double[gs];
  }

  #pragma acc enter data create(grid1[0:gs6],grid2[0:gs6])
  #pragma acc enter data create(val1[0:iN][0:gs],val2[0:iN][0:gs])
  if (gga)
  {
    #pragma acc enter data create(grho[0:gs3],val1p[0:iN][0:gs3],val2p[0:iN][0:gs3])
  }
  if (tau)
  {
    #pragma acc enter data create(val1T[0:iN][0:gs])
  }

  for (int j=0;j<N2;j++)
    fxc[j] = 0.;

  if (gga)
  #pragma acc parallel loop present(grho[0:gs3])
  for (int j=0;j<gs3;j++)
    grho[j] = 0.;

  //first assemble grho
  if (gga)
  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    double Z1 = (double)atno[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];

    copy_grid(gs,grid1,grid);
    recenter_grid_zero(gs,grid1,-A1,-B1,-C1);

    #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs])
    for (int i1=0;i1<s2-s1;i1++)
    for (int j=0;j<gs;j++)
      val1[i1][j] = 1.;

    #pragma acc parallel loop collapse(2) present(val1p[0:iN][0:gs3])
    for (int i1=0;i1<s2-s1;i1++)
    for (int j=0;j<gs3;j++)
      val1p[i1][j] = 1.;

    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

      eval_shd(ii1,gs,grid1,val1[ii1],n1,l1,m1,zeta1);
      eval_pd(tid,gs,grid1,val1p[ii1],n1,l1,m1,zeta1);
    }

    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;
      double* valn = val1[ii1];
      double* valpn = val1p[ii1];

      for (int i2=s1;i2<s2;i2++)
      {
        int ii2 = i2-s1;
        double* valm = val1[ii2];
        double* valpm = val1p[ii2];

        double d1 = Paon[i1*N+i2];

       #pragma acc parallel loop present(grho[0:gs3],valn[0:gs],valpn[0:gs3],valm[0:gs],valpm[0:gs3])
        for (int j=0;j<gs;j++)
        {
          grho[3*j+0] += d1*(valn[j]*valpm[3*j+0]+valpn[3*j+0]*valm[j]);
          grho[3*j+1] += d1*(valn[j]*valpm[3*j+1]+valpn[3*j+1]*valm[j]);
          grho[3*j+2] += d1*(valn[j]*valpm[3*j+2]+valpn[3*j+2]*valm[j]);
        }
      }
    }

    for (int n=m+1;n<natoms;n++)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      double Z2 = (double)atno[m];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];

      #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs])
      for (int i2=0;i2<s4-s3;i2++)
      for (int j=0;j<gs;j++)
        val2[i2][j] = 1.;

      #pragma acc parallel loop collapse(2) present(val2p[0:iN][0:gs3])
      for (int i2=0;i2<s4-s3;i2++)
      for (int j=0;j<gs3;j++)
        val2p[i2][j] = 1.;

      copy_grid(gs,grid2,grid);
      recenter_grid_zero(gs,grid2,-A2,-B2,-C2);

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;

        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        eval_shd(ii2,gs,grid2,val2[ii2],n2,l2,m2,zeta2);
        eval_pd(tid,gs,grid2,val2p[ii2],n2,l2,m2,zeta2);
      }

      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;
        double* valn = val1[ii1];
        double* valpn = val1p[ii1];

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          double* valm = val2[ii2];
          double* valpm = val2p[ii2];

          double d1 = Paon[i1*N+i2];

         #pragma acc parallel loop present(grho[0:gs3],valn[0:gs],valpn[0:gs3],valm[0:gs],valpm[0:gs3])
          for (int j=0;j<gs;j++)
          {
            grho[3*j+0] += d1*(valn[j]*valpm[3*j+0]+valpn[3*j+0]*valm[j]);
            grho[3*j+1] += d1*(valn[j]*valpm[3*j+1]+valpn[3*j+1]*valm[j]);
            grho[3*j+2] += d1*(valn[j]*valpm[3*j+2]+valpn[3*j+2]*valm[j]);
          }
        } //loop i2
      } //loop i1

    } //loop n
  } //loop m


  if (prl>1 && gs<1000)
  {
    #pragma acc update self(grho[0:gs3])
    printf(" grho: \n");
    for (int j=0;j<gs;j++)
      printf("   %8.5f %8.5f %8.5f \n",grho[3*j+0],grho[3*j+1],grho[3*j+2]);
    printf("\n");
  }

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    double Z1 = (double)atno[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];

    copy_grid(gs,grid1,grid);
    recenter_grid_zero(gs,grid1,-A1,-B1,-C1);

    #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs])
    for (int i1=0;i1<s2-s1;i1++)
    for (int j=0;j<gs;j++)
      val1[i1][j] = 1.;

    if (gga)
    #pragma acc parallel loop collapse(2) present(val1p[0:iN][0:gs3])
    for (int i1=0;i1<s2-s1;i1++)
    for (int j=0;j<gs3;j++)
      val1p[i1][j] = 1.;

    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

      if (ss_basis)
        eval_ssd(ii1,gs,grid1,val1[ii1],n1,l1,m1,zeta1,Rc);
      else if (sgs_basis)
        eval_sgsd(ii1,gs1,gs2,grid1,val1[ii1],n1,l1,m1,zeta1,Rc);
      else
        eval_shd(ii1,gs,grid1,val1[ii1],n1,l1,m1,zeta1);
      if (gga)
        eval_pd(tid,gs,grid1,val1p[ii1],n1,l1,m1,zeta1);
    }

   //single-atom elements over grid
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;
      double* valn = val1[ii1];
      double* valpn = NULL; if (gga) valpn = val1p[ii1];
      //double* valtn = NULL; if (tau) valtn = val1T[ii1];

      for (int i2=s1;i2<s2;i2++)
      {
        int ii2 = i2-s1;
        double* valm = val1[ii2];
        double* valpm = NULL; if (gga) valpm = val1p[ii2];

        double valt = 0.;
        if (need_wt)
        {
          #pragma acc parallel loop present(vxc[0:gs],valm[0:gs],valn[0:gs],wt[0:gs]) reduction(+:valt)
          for (int j=0;j<gs;j++)
            valt += valn[j]*valm[j]*vxc[j]*wt[j];
        }
        else
        {
          #pragma acc parallel loop present(vxc[0:gs],valm[0:gs],valn[0:gs]) reduction(+:valt)
          for (int j=0;j<gs;j++)
            valt += valn[j]*valm[j]*vxc[j];
        }

        if (gga)
        #pragma acc parallel loop present(grho[0:gs3],vxcs[0:gs],valm[0:gs],valpm[0:gs3],valn[0:gs],valpn[0:gs3]) reduction(+:valt)
        for (int j=0;j<gs;j++)
        {
          double grx = grho[3*j+0]; double gry = grho[3*j+1]; double grz = grho[3*j+2];
          double valx = valn[j]*valpm[3*j+0]+valpn[3*j+0]*valm[j];
          double valy = valn[j]*valpm[3*j+1]+valpn[3*j+1]*valm[j];
          double valz = valn[j]*valpm[3*j+2]+valpn[3*j+2]*valm[j];
          valt += vxcs[j]*(valx*grx+valy*gry+valz*grz);
        }

       #if 0
        if (tau) //incomplete
        {
          valt = 0.;
          #pragma acc parallel loop present(vxc[0:gs],valtn[0:gs],valm[0:gs]) reduction(+:valt)
          for (int j=0;j<gs;j++)
            valt += valtn[j]*valm[j]*vxc[j];
        }
       #endif

        fxc[i1*N+i2] = valt;
      }
    }

    for (int n=m+1;n<natoms;n++)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      double Z2 = (double)atno[m];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];

      #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs])
      for (int i2=0;i2<s4-s3;i2++)
      for (int j=0;j<gs;j++)
        val2[i2][j] = 1.;

      if (gga)
      #pragma acc parallel loop collapse(2) present(val2p[0:iN][0:gs3])
      for (int i2=0;i2<s4-s3;i2++)
      for (int j=0;j<gs3;j++)
        val2p[i2][j] = 1.;

      copy_grid(gs,grid2,grid);
      recenter_grid_zero(gs,grid2,-A2,-B2,-C2);

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;

        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        eval_shd(ii2,gs,grid2,val2[ii2],n2,l2,m2,zeta2);
        if (gga)
          eval_pd(tid,gs,grid2,val2p[ii2],n2,l2,m2,zeta2);
      }

     //two-atom elements over grid
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;
        double* valn = val1[ii1];
        double* valpn = NULL; if (gga) valpn = val1p[ii1];
        //double* valtn = NULL; if (tau) valtn = val1T[ii1];

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          double* valm = val2[ii2];
          double* valpm = NULL; if (gga) valpm = val2p[ii2];
          //double* valtm = NULL; if (tau) valtm = val2T[ii2];

          double valt = 0.;
          if (need_wt)
          {
            #pragma acc parallel loop present(vxc[0:gs],valm[0:gs],valn[0:gs],wt[0:gs]) reduction(+:valt)
            for (int j=0;j<gs;j++)
              valt += valn[j]*valm[j]*vxc[j]*wt[j];
          }
          else
          {
            #pragma acc parallel loop present(vxc[0:gs],valm[0:gs],valn[0:gs]) reduction(+:valt)
            for (int j=0;j<gs;j++)
              valt += valn[j]*valm[j]*vxc[j];
          }

          if (gga)
          #pragma acc parallel loop present(grho[0:gs3],vxcs[0:gs],valm[0:gs],valpm[0:gs3],valn[0:gs],valpn[0:gs3]) reduction(+:valt)
          for (int j=0;j<gs;j++)
          {
            double grx = grho[3*j+0]; double gry = grho[3*j+1]; double grz = grho[3*j+2];
            double valx = valn[j]*valpm[3*j+0]+valpn[3*j+0]*valm[j];
            double valy = valn[j]*valpm[3*j+1]+valpn[3*j+1]*valm[j];
            double valz = valn[j]*valpm[3*j+2]+valpn[3*j+2]*valm[j];
            valt += vxcs[j]*(valx*grx+valy*gry+valz*grz);
          }

         #if 0
          if (tau) //incomplete
          {
            valt = 0.;
            #pragma acc parallel loop present(vxc[0:gs],valtn[0:gs],valm[0:gs]) reduction(+:valt)
            for (int j=0;j<gs;j++)
              valt += valtn[j]*valm[j]*vxc[j];
          }
         #endif

          fxc[i1*N+i2] = fxc[i2*N+i1] = valt;

        } //loop i2
      } //loop i1

    } //loop n over unique atoms
  } //loop m over natoms

  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    fxc[i*N+j] *= norm[i]*norm[j];

  #pragma acc update device(fxc[0:N2])


 //cleanup
  #pragma acc exit data delete(grid1[0:gs6],grid2[0:gs6])
  #pragma acc exit data delete(val1[0:iN][0:gs],val2[0:iN][0:gs])
  if (gga)
  {
    #pragma acc exit data delete(grho[0:gs3],val1p[0:iN][0:gs3],val2p[0:iN][0:gs3])
  }
  if (tau)
  {
    #pragma acc exit data delete(val1T[0:iN][0:gs])
  }

  delete [] n2i;

  delete [] grid1;
  delete [] grid2;

  for (int i=0;i<iN;i++)
    delete [] val1[i];
  for (int i=0;i<iN;i++)
    delete [] val2[i];
  delete [] val1;
  delete [] val2;

  if (gga)
  {
    delete [] Paon;
    delete [] grho;
    for (int i=0;i<iN;i++)
      delete [] val1p[i];
    for (int i=0;i<iN;i++)
      delete [] val2p[i];
    delete [] val1p;
    delete [] val2p;
  }

  if (tau)
  {
    for (int i=0;i<iN;i++)
      delete [] val1T[i];
    delete [] val1T;
  }

  return;
}

void compute_fxc(int natoms, int* atno, double* coords, vector<vector<double> > &basis, bool gga, bool tau, bool need_wt, double* Pao, float* vxc, float* vxcs, int nrad, int gsa, float* grid, float* wt, double* fxc, int prl)
{
  if (vxcs!=NULL) { printf("\n ERROR in compute_fxc (float version) vxcs must be NULL \n"); exit(1); }

  double* vxcd = new double[gsa];
  #pragma acc enter data create(vxcd[0:gsa])

  #pragma acc parallel loop present(vxcd[0:gsa],vxc[0:gsa])
  for (int j=0;j<gsa;j++)
    vxcd[j] = vxc[j];

  compute_fxc(natoms,atno,coords,basis,gga,tau,need_wt,Pao,vxcd,NULL,nrad,gsa,grid,wt,fxc,prl);

  #pragma acc exit data delete(vxcd[0:gsa])
  delete [] vxcd;

  return;
}


void compute_dft_grad(int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* Pao, bool is_gga, double* vxc, double* vxcs, int gs, float* grid, double* grad, int prl)
{
  printf("\n WARNING: compute_dft_grad is in testing \n");

  int gsa = natoms*gs;
  int gsa3 = 3*gsa;
  int gsa6 = 6*gsa;

  int N = basis.size();
  int N2 = N*N;
  int N3 = 3*natoms;
  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);

  float* grid1 = new float[gsa6];
  float* grid2 = new float[gsa6];

  double norm[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];

  double* Paon = new double[N2];
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    Paon[i*N+j] = Pao[i*N+j]*norm[i]*norm[j];

  //double* norm12 = new double[N2];
  //for (int i=0;i<N;i++)
  //for (int j=0;j<N;j++)
  //  norm12[i*N+j] = norm[i]*norm[j];

  int iN = imaxN;
  float** val1 = new float*[iN];
  float** val2 = new float*[iN];
  for (int i=0;i<iN;i++)
    val1[i] = new float[gsa];
  for (int i=0;i<iN;i++)
    val2[i] = new float[gsa];

  float** val1p = new float*[iN];
  float** val2p = new float*[iN];
  for (int i=0;i<iN;i++)
    val1p[i] = new float[gsa3];
  for (int i=0;i<iN;i++)
    val2p[i] = new float[gsa3];

  float** val1h = new float*[iN];
  float** val2h = new float*[iN];
  if (is_gga)
  for (int i=0;i<iN;i++)
    val1h[i] = new float[gsa6];
  if (is_gga)
  for (int i=0;i<iN;i++)
    val2h[i] = new float[gsa6];

  #pragma acc enter data create(grid1[0:gsa6],grid2[0:gsa6])
  #pragma acc enter data create(val1[0:iN][0:gsa],val2[0:iN][0:gsa],val1p[0:iN][0:gsa3],val2p[0:iN][0:gsa3])
  if (is_gga)
  {
    #pragma acc enter data create(val1h[0:iN][0:gsa6],val2h[0:iN][0:gsa6])
  }

  for (int j=0;j<N3;j++)
    grad[j] = 0.;

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    float Z1 = (float)atno[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];

    copy_grid(gsa,grid1,grid);
    recenter_grid_zero(gsa,grid1,-A1,-B1,-C1);

    #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gsa])
    for (int i1=0;i1<s2-s1;i1++)
    for (int j=0;j<gsa;j++)
      val1[i1][j] = 1.f;

    #pragma acc parallel loop collapse(2) present(val1p[0:iN][0:gsa3])
    for (int i1=0;i1<s2-s1;i1++)
    for (int j=0;j<gsa3;j++)
      val1p[i1][j] = 1.f;

    if (is_gga)
    #pragma acc parallel loop collapse(2) present(val1h[0:iN][0:gsa6])
    for (int i1=0;i1<s2-s1;i1++)
    for (int j=0;j<gsa6;j++)
      val1h[i1][j] = 1.f;

    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

      eval_sh(ii1,gsa,grid1,val1[ii1],n1,l1,m1,zeta1);
      eval_p(gsa,grid1,val1p[ii1],n1,l1,m1,zeta1);
      if (is_gga)
        eval_h(gsa,grid1,val1h[ii1],n1,l1,m1,zeta1);
    }

   #if 0
   //single-atom elements should be zero
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;
      float* valm = val1[ii1];
      float* valpm = val1p[ii1];
      //float* valhm = val1h[ii1];

      for (int i2=s1;i2<s2;i2++)
      {
        int ii2 = i2-s1;
        float* valn = val1[ii2];
        float* valpn = val1p[ii2];
        //float* valhn = val1h[ii2];

        //double valt = 0.;
        //for (int j=0;j<gsa;j++)
        //  valt += valm[j]*valn[j]*vxc[j];

        double d12 = Paon[i1*N+i2];
       #pragma acc parallel loop present(valm[0:gsa],valn[0:gsa],vxc[0:gsa]) reduction(+:gradx,grady,gradz)
        for (int j=0;j<gsa;j++)
        {
          double valx = valm[j]*valpn[3*j+0] + valpm[3*j+0]*valn[j];
          double valy = valm[j]*valpn[3*j+1] + valpm[3*j+1]*valn[j];
          double valz = valm[j]*valpn[3*j+2] + valpm[3*j+2]*valn[j];
          grad[3*m+0] += vxc[j]*valx*d12;
          grad[3*m+1] += vxc[j]*valy*d12;
          grad[3*m+2] += vxc[j]*valz*d12;
          //valt += vxcs[j]*(valx*valx+valy*valy+valz*valz);
        }
      }
    } //single-atom i1,i2
   #endif

    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      float Z2 = (float)atno[m];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];

      #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gsa])
      for (int i2=0;i2<s4-s3;i2++)
      for (int j=0;j<gsa;j++)
        val2[i2][j] = 1.f;

      #pragma acc parallel loop collapse(2) present(val2p[0:iN][0:gsa3])
      for (int i2=0;i2<s4-s3;i2++)
      for (int j=0;j<gsa3;j++)
        val2p[i2][j] = 1.f;

      if (is_gga)
      #pragma acc parallel loop collapse(2) present(val2h[0:iN][0:gsa6])
      for (int i2=0;i2<s4-s3;i2++)
      for (int j=0;j<gsa6;j++)
        val2h[i2][j] = 1.f;

      copy_grid(gsa,grid2,grid);
      recenter_grid_zero(gsa,grid2,-A2,-B2,-C2);

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;

        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        eval_sh(ii2,gsa,grid2,val2[ii2],n2,l2,m2,zeta2);
        eval_p(gsa,grid2,val2p[ii2],n2,l2,m2,zeta2);
        if (is_gga)
          eval_h(gsa,grid2,val2h[ii2],n2,l2,m2,zeta2);
      }

     //two-atom Pao elements added to grid
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;
        float* valm = val1[ii1];
        float* valpm = val1p[ii1];
        //float* valhm = val1h[ii1];

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          float* valn = val2[ii2];
          float* valpn = val2p[ii2];
          float* valhn = val2h[ii2];

          float gradx = 0.; float grady = 0.; float gradz = 0.;

          float d12 = Paon[i1*N+i2];
          if (fabs(d12)>1.e-12)
          {
           #pragma acc parallel loop present(vxc[0:gsa],valm[0:gsa],valpn[0:gsa3]) reduction(+:gradx,grady,gradz)
            for (int j=0;j<gsa;j++)
            {
              float vxc1 = vxc[j];
              float valm1 = valm[j];
              float valxn = valm1*valpn[3*j+0]; float valyn = valm1*valpn[3*j+1]; float valzn = valm1*valpn[3*j+2];
              gradx += vxc1*valxn;
              grady += vxc1*valyn;
              gradz += vxc1*valzn;
            }
            gradx *= d12; grady *= d12; gradz *= d12;

            //printf("  i1/2: %2i %2i  gradxyz: %8.5f %8.5f %8.5f \n",i1,i2,gradx,grady,gradz);
            grad[3*n+0] += gradx;
            grad[3*n+1] += grady;
            grad[3*n+2] += gradz;
          }

          if (is_gga)
          {
            printf("  TESTING: gga part of gradient \n");
            gradx = grady = gradz = 0.;

           #pragma acc parallel loop present(vxc[0:gsa],vxcs[0:gsa],valm[0:gsa],valn[0:gsa],valpm[0:gsa3],valpn[0:gsa3],valhn[0:gsa6]) reduction(+:gradx,grady,gradz)
            for (int j=0;j<gsa;j++)
            {
              float vmj = valm[j]; float vnj = valn[j]; 
              float valx1 = valpm[3*j+0]*vnj; float valy1 = valpm[3*j+1]*vnj; float valz1 = valpm[3*j+2]*vnj;
              float valx2 = vmj*valpn[3*j+0]; float valy2 = vmj*valpn[3*j+1]; float valz2 = vmj*valpn[3*j+2];
              float valx = valx1 + valx2; float valy = valy1 + valy2; float valz = valz1 + valz2;

             //xx, xy, xz, yy, yz, zz
              float vx_x = valpm[3*j+0]*valpn[3*j+0] + vmj*valhn[6*j+0];
              float vx_y = valpm[3*j+0]*valpn[3*j+1] + vmj*valhn[6*j+1];
              float vx_z = valpm[3*j+0]*valpn[3*j+2] + vmj*valhn[6*j+2];

              float vy_x = valpm[3*j+1]*valpn[3*j+0] + vmj*valhn[6*j+1];
              float vy_y = valpm[3*j+1]*valpn[3*j+1] + vmj*valhn[6*j+3];
              float vy_z = valpm[3*j+1]*valpn[3*j+2] + vmj*valhn[6*j+4];

              float vz_x = valpm[3*j+2]*valpn[3*j+0] + vmj*valhn[6*j+2];
              float vz_y = valpm[3*j+2]*valpn[3*j+1] + vmj*valhn[6*j+4];
              float vz_z = valpm[3*j+2]*valpn[3*j+2] + vmj*valhn[6*j+5];

              float vxcsd = vxcs[j];
              gradx += vxcsd*(vx_x*valx + vy_x*valy + vz_x*valz);
              grady += vxcsd*(vx_y*valx + vy_y*valy + vz_y*valz);
              gradz += vxcsd*(vx_z*valx + vy_z*valy + vz_z*valz);
            }
            grad[3*n+0] += gradx*d12;
            grad[3*n+1] += grady*d12;
            grad[3*n+2] += gradz*d12;
          }

        } //loop i2
      } //loop i1

    } //loop n over unique atoms

  } //loop m over natoms

  for (int j=0;j<N3;j++) 
    grad[j] *= -2.;

  printf("  xc_grad: \n");
  print_gradient(natoms,grad);


 //cleanup
  #pragma acc exit data delete(grid1[0:gsa6],grid2[0:gsa6])
  #pragma acc exit data delete(val1[0:iN][0:gsa],val2[0:iN][0:gsa],val1p[0:iN][0:gsa3],val2p[0:iN][0:gsa3])

  if (is_gga)
  {
    #pragma acc exit data delete(val1h[0:iN][0:gsa6],val2h[0:iN][0:gsa6])
  }

  delete [] n2i;

  delete [] Paon;
  //delete [] norm12;

  delete [] grid1;
  delete [] grid2;

  for (int i=0;i<iN;i++)
    delete [] val1[i];
  for (int i=0;i<iN;i++)
    delete [] val2[i];
  delete [] val1;
  delete [] val2;

  for (int i=0;i<iN;i++)
    delete [] val1p[i];
  for (int i=0;i<iN;i++)
    delete [] val2p[i];
  delete [] val1p;
  delete [] val2p;

  if (is_gga)
  for (int i=0;i<iN;i++)
    delete [] val1h[i];
  if (is_gga)
  for (int i=0;i<iN;i++)
    delete [] val2h[i];
  delete [] val1h;
  delete [] val2h;

  return;
}

void test_becke_weight(int natoms, int* atno, double* coords, vector<vector<double> >& basis, int nrad, int nang, double* ang_g0, double* ang_w0)
{
  printf("\n\n Becke weight test function (not enabled) \n");
#if 0

  int prl = 1;

  int gc = 6;
  int gs = nrad*nang;
  int gsc = gc*gs;
  int gsa = natoms*gs;

  int N = basis.size();
  int N2 = N*N;
  int N3 = 3*natoms;

  float* ang_g = new float[3*nang];
  float* ang_w = new float[nang];
  for (int i=0;i<3*nang;i++) ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++) ang_w[i] = ang_w0[i];

  //printf(" ang_g: \n");
  //for (int m=0;m<nang;m++)
  //  printf(" %8.5f %8.5f %8.5f \n",ang_g[3*m],ang_g[3*m+1],ang_g[3*m+2]);
  //printf("\n");

  float* grid = new float[natoms*gsc];
  float* wt = new float[gsa];

  #pragma acc enter data create(grid[0:natoms*gsc],wt[0:gsa])
  #pragma acc enter data copyin(coords[0:N3],ang_g[0:3*nang],ang_w[0:nang])

  get_becke_grid_full(natoms,atno,coords,nrad,nang,ang_g,ang_w,gc,grid,wt);

  if (0)
  {
    printf("  grid: \n");
    for (int n=0;n<natoms;n++)
    {
      for (int m=0;m<gs;m++)
        printf(" %8.5f %8.5f %8.5f \n",grid[n*gsc+gc*m+0],grid[n*gsc+gc*m+1],grid[n*gsc+gc*m+2]);
      printf("\n");
    }
    printf("\n");
  }

  if (0)
  {
    printf(" wt: \n");
    for (int n=0;n<natoms;n++)
    {
      for (int m=0;m<gs;m++)
        printf(" %4.2e",wt[n*gs+m]);
      printf("\n");
    }
  }

  double* Pao = new double[N2]();
  read_square(N,Pao,"Pao");

  double* rho = new double[gsa]();
  double* drho = new double[gsa]();
  #pragma acc enter data create(rho[0:gsa],drho[0:gsa])

  compute_rho(natoms,atno,coords,basis,Pao,nrad,gsa,grid,rho,drho,1);
  #pragma acc update self(rho[0:gsa],drho[0:gsa])

  printf(" rho: \n");
  for (int n=0;n<natoms;n++)
  {
    for (int m=0;m<gs;m+=pr_inc)
      printf(" %8.5f",rho[n*gs+m]);
    printf("\n");
  }

  printf(" drho: \n");
  for (int n=0;n<natoms;n++)
  {
    for (int m=0;m<gs;m+=pr_inc)
      printf(" %8.5f",drho[n*gs+m]);
    printf("\n");
  }

 // 1 + 7 --> LDA
 // 101 + 130 --> PBE

  int x_type, c_type;
  if (!read_dft(x_type,c_type,1))
  {
    printf("  ERROR: please provide DFT functional \n");
    return;
  }

 //exchange component
  int func_id = x_type;
  xc_func_type func1;
  xc_func_init(&func1, func_id, XC_UNPOLARIZED);

  bool is_gga = 0;
  if (func1.info->family==XC_FAMILY_GGA)
  {
    printf("\n found a GGA functional \n");
    is_gga = 1;
  }
  else if (func1.info->family==XC_FAMILY_LDA)
    printf("\n found an LDA functional \n");

  double* ex = new double[gsa];
  double* vx = new double[gsa];
  double* vxs = new double[gsa]();

  if (func1.info->family==XC_FAMILY_LDA)
    xc_lda_exc_vxc(&func1,gsa,rho,ex,vx);
  else if (func1.info->family==XC_FAMILY_GGA)
    xc_gga_exc_vxc(&func1,gsa,rho,drho,ex,vx,vxs);
  else
  {
    printf(" ERROR: unrecognized DFT type \n");
    exit(1);
  }

 //correlation component
  func_id = c_type;
  xc_func_type func2;
  xc_func_init(&func2, func_id, XC_UNPOLARIZED);

  double* ec = new double[gsa];
  double* vc = new double[gsa];
  double* vcs = new double[gsa]();

  if (func1.info->family==XC_FAMILY_LDA)
    xc_lda_exc_vxc(&func2,gsa,rho,ec,vc);
  else if (func1.info->family==XC_FAMILY_GGA)
    xc_gga_exc_vxc(&func2,gsa,rho,drho,ec,vc,vcs);
  else
  {
    printf(" ERROR: unrecognized DFT type \n");
    exit(1);
  }

  double* exc = new double[gsa];
  double* vxc = new double[gsa];
  double* vxcs = new double[gsa];
  for (int m=0;m<gsa;m++)
    exc[m] = ex[m]+ec[m];
  for (int m=0;m<gsa;m++)
    vxc[m] = vx[m]+vc[m];
  for (int m=0;m<gsa;m++)
    vxcs[m] = vxs[m]+vcs[m];

  if (prl>0)
  {
    printf(" exc: \n");
    for (int n=0;n<natoms;n++)
    {
      for (int m=0;m<gs;m+=pr_inc)
        printf(" %8.5f",exc[n*gs+m]);
      printf("\n");
    }
    printf("\n");
    printf(" vxc: \n");
    for (int n=0;n<natoms;n++)
    {
    for (int m=0;m<gs;m+=pr_inc)
        printf(" %8.5f",vxc[n*gs+m]);
      printf("\n");
    }
    printf("\n");
  }
  if (func1.info->family==XC_FAMILY_GGA)
  {
    printf(" vxcs: \n");
    for (int n=0;n<natoms;n++)
    {
      for (int m=0;m<gs;m+=pr_inc)
        printf(" %8.5f",vxcs[n*gs+m]);
      printf("\n");
    }
    printf("\n");
  }

  printf("\n  vxc along z axis \n");
  print_vxc(nrad,nang,natoms,grid,vxc);

  double den = 0.;
  for (int j=0;j<gsa;j++)
    den += rho[j]*wt[j];
  printf(" density: %12.8f \n",den);

  for (int j=0;j<gsa;j++)
    exc[j] *= wt[j];
  for (int j=0;j<gsa;j++)
    vxc[j] *= wt[j];
  for (int j=0;j<gsa;j++)
    vxcs[j] *= wt[j];

  double Exc = 0.;
  for (int j=0;j<gsa;j++)
    Exc += rho[j]*exc[j];
  printf(" Exc: %12.8f \n",Exc);

  if (!is_gga)
  for (int j=0;j<gsa;j++)
    vxcs[j] = 0.;

  #pragma acc enter data copyin(vxc[0:gsa],vxcs[0:gsa])

  double* fxc = new double[N2];
  for (int m=0;m<N2;m++) fxc[m] = 0.;
  compute_fxc(natoms,atno,coords,basis,is_gga,0,0,Pao,vxc,vxcs,nrad,gsa,grid,wt,fxc,1);

  printf("\n fxc: \n");
  print_square(N,fxc);


  double* grad = new double[N3];
  compute_dft_grad(natoms,atno,coords,basis,Pao,is_gga,vxc,vxcs,gs,grid,grad,1);

  delete [] grad;


 //cleanup
  #pragma acc exit data delete(vxc[0:gsa],vxcs[0:gsa])
  #pragma acc exit data delete(rho[0:gsa],drho[0:gsa])

  #pragma acc exit data delete(grid[0:natoms*gsc],wt[0:gsa])
  #pragma acc exit data delete(coords[0:N3],ang_g[0:3*nang],ang_w[0:nang])

  delete [] grid;
  delete [] wt;

  delete [] Pao;
  delete [] rho;
  delete [] drho;

  delete [] exc;
  delete [] vxc;

  delete [] ex;
  delete [] ec;
  delete [] vx;
  delete [] vc;
  delete [] vxs;
  delete [] vcs;

  delete [] fxc;
#endif

  return;
}

