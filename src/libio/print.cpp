#include "print.h"

//this file includes functions that print to screen
// and don't require any particular calculation to do so

void print_coords(int natoms, double* coords)
{
  float B2A = 1.f/A2B;
  for (int n=0;n<natoms;n++)
    printf(" %10.5f %10.5f %10.5f \n",coords[3*n+0]*B2A,coords[3*n+1]*B2A,coords[3*n+2]*B2A);
}

void print_gradient(int natoms, double* grad)
{
  for (int n=0;n<natoms;n++)
    printf(" %10.5f %10.5f %10.5f \n",grad[3*n+0],grad[3*n+1],grad[3*n+2]);
}

void print_square_diff(int N, double* S1, double* S2)
{
  for (int n=0;n<N;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %10.5f",S1[n*N+m],S1[n*N+m]);
    printf("\n");
  }
}

void print_square_fine(int N, float* S)
{
  for (int n=0;n<N;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %12.8f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_fine(int M, int N, double* S)
{
  for (int n=0;n<M;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %12.8f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_fine(int N, double* S)
{
  for (int n=0;n<N;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %12.8f",S[n*N+m]);
    printf("\n");
  }
}

void print_square(int N, double* S)
{
  for (int n=0;n<N;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}

void print_square(int N, float* S)
{
  for (int n=0;n<N;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}

void print_square(int M, int N, float* S)
{
  for (int n=0;n<M;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}

void print_square(int M, int N, double* S)
{
  for (int n=0;n<M;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}

void print_mos_col(int M, int N, vector<vector<double> > basis, double* jCA)
{
  if (M<1) return;
  for (int j=0;j<N;j++)
  {
    int n1 = basis[j][0]; int l1 = basis[j][1]; int m1 = basis[j][2];

    printf("  %i%i%2i ",n1,l1,m1);
    for (int k=0;k<M;k++)
      printf(" %10.5f",jCA[j*N+k]);
    printf("\n");
  }
}

void print_square_col(int M, int N, double* S)
{
  for (int n=0;n<N;n++)
  {
    for (int m=0;m<M;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_col_sm(int M, int N, double* S)
{
  for (int n=0;n<N;n++)
  {
    for (int m=0;m<M;m++)
      printf(" %6.3f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_ss(int N, double* S)
{
  for (int n=0;n<N;n++)
  {
    printf("  ");
    for (int m=0;m<N;m++)
      printf(" %8.5f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_sm(int N, double* S)
{
  for (int n=0;n<N;n++)
  {
    printf("  ");
    for (int m=0;m<N;m++)
      printf(" %5.2f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_sm(int N, float* S)
{
  for (int n=0;n<N;n++)
  {
    printf("  ");
    for (int m=0;m<N;m++)
      printf(" %5.2f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_ss_sm(int N, double* S)
{
  for (int n=0;n<N;n++)
  {
    printf("  ");
    for (int m=0;m<N;m++)
      printf(" %6.3f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_nxn(int No, int N, float* S)
{
  for (int n=0;n<No;n++)
  {
    for (int m=0;m<No;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_nxn(int No, int N, double* S)
{
  for (int n=0;n<No;n++)
  {
    for (int m=0;m<No;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}

void print_rectangle(int N1, int N2, double* S)
{
  for (int n=0;n<N1;n++)
  {
    printf("  ");
    for (int m=0;m<N2;m++)
      printf(" %10.5f",S[n*N2+m]);
    printf("\n");
  }
}

void print_rectangle_e(int N1, int N2, double* S)
{
  for (int n=0;n<N1;n++)
  {
    printf("  ");
    for (int m=0;m<N2;m++)
      printf(" %5.1e",S[n*N2+m]);
    printf("\n");
  }
}

void print_rectangle_sm(int N1, int N2, double* S)
{
  for (int n=0;n<N1;n++)
  {
    printf("  ");
    for (int m=0;m<N2;m++)
      printf(" %7.4f",S[n*N2+m]);
    printf("\n");
  }
}

void print_rdm(int M, float* rdm)
{
  int M2 = M*M;
  int M3 = M*M2;
  printf("\n rdm: \n");
  for (int p=0;p<M;p++)
  {
    for (int q=0;q<M;q++)
      print_square_fine(M,&rdm[p*M3+q*M2]);
    printf("\n");
  }
  printf("\n");

  return;
}

void print_rdm(int M, double* rdm)
{
  int M2 = M*M;
  int M3 = M*M2;
  printf("\n rdm: \n");
  for (int p=0;p<M;p++)
  {
    for (int q=0;q<M;q++)
      print_square_fine(M,&rdm[p*M3+q*M2]);
    printf("\n");
  }
  printf("\n");

  return;
}

void print_vec(int gsa, float* grid, float* vxc)
{
  printf("	 z        vector \n");
  for (int m=0;m<gsa;m++)
  {
    float x1 = grid[6*m+0];
    float y1 = grid[6*m+1];
    float z1 = grid[6*m+2];
    //float r1 = grid[6*m+3];
    if (x1==0.f && y1==0.f)
      printf("  %9.5f  %8.5f \n",z1,vxc[m]);
  }
  return;
}

void print_vec(int gsa, float* grid, double* vxc)
{
  printf("	 z        vector \n");
  for (int m=0;m<gsa;m++)
  {
    float x1 = grid[6*m+0];
    float y1 = grid[6*m+1];
    float z1 = grid[6*m+2];
    //float r1 = grid[6*m+3];
    if (x1==0.f && y1==0.f)
      printf("  %9.5f  %8.5f \n",z1,vxc[m]);
  }
  return;
}

void print_vec(int gsa, float* grid, float* A, float* B)
{
  printf("       z        vector     vector \n");
  for (int m=0;m<gsa;m++)
  {
    float x1 = grid[6*m+0];
    float y1 = grid[6*m+1];
    float z1 = grid[6*m+2];
    //float r1 = grid[6*m+3];
    if (x1==0.f && y1==0.f)
      printf("  %9.5f  %8.5f %8.5f \n",z1,A[m],B[m]);
  }
  return;
}

void print_vec(int gsa, float* grid, double* A, double* B)
{
  printf("       z        vector     vector \n");
  for (int m=0;m<gsa;m++)
  {
    float x1 = grid[6*m+0];
    float y1 = grid[6*m+1];
    float z1 = grid[6*m+2];
    //float r1 = grid[6*m+3];
    if (x1==0.f && y1==0.f)
      printf("  %9.5f  %8.5f %8.5f \n",z1,A[m],B[m]);
  }
  return;
}


void print_dft_vals(int natoms, int gs, double* grid, double* rho, double* drho, double* Td, double* vc, int zpos, bool scirep)
{
  double minz = 1.e-6;
  if (scirep) minz = 1.e-12;

  printf("       z           rho         drho        tau       vc \n");

  if (scirep)
  for (int n=0;n<natoms;n++)
  for (int m=0;m<gs;m++)
  {
    int mn = n*gs+m;
    float x1 = grid[6*mn+0]; float y1 = grid[6*mn+1]; float z1 = grid[6*mn+2];
    if (!zpos || z1>0.)
    if (fabs(z1)>minz)
    if (fabs(x1)<1.e-8 && fabs(y1)<1.e-8)
      printf("  %13.8e  %13.8e  %13.8e  %13.8e  %13.8e \n",z1,rho[mn],drho[mn],Td[mn],vc[mn]);
  }
  else
  for (int n=0;n<natoms;n++)
  for (int m=0;m<gs;m++)
  {
    int mn = n*gs+m;
    float x1 = grid[6*mn+0]; float y1 = grid[6*mn+1]; float z1 = grid[6*mn+2];
    if (!zpos || z1>0.)
    if (fabs(z1)>minz)
    if (fabs(x1)<1.e-8 && fabs(y1)<1.e-8)
      printf("  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f \n",z1,rho[mn],drho[mn],Td[mn],vc[mn]);
  }

  return;
}


void print_vc(int natoms, int gs, float* grid, float* rho, float* vc, int mode)
{
  printf("        z        rho       vc \n");
  for (int n=0;n<natoms;n++)
  for (int m=0;m<gs;m++)
  {
    int mn = n*gs+m;
    float x1 = grid[6*mn+0]; float y1 = grid[6*mn+1]; float z1 = grid[6*mn+2];
    if (fabs(z1)>1.e-5)
    if (fabs(z1)<8. || mode>1)
    if (fabs(x1)<1.e-8 && fabs(y1)<1.e-8)
    {
      printf("  %9.5f  %8.5f  %8.5f \n",z1,rho[mn],vc[mn]);
      //printf("  %9.5f %9.5f %9.5f  %8.5f  %8.5f \n",x1,y1,z1,rho[mn],vc[mn]);
    }
  }

  return;
}

void print_vc(int natoms, int gs, float* grid, double* rho, double* vc, int mode)
{
  printf("        z        rho       vc \n");
  for (int n=0;n<natoms;n++)
  for (int m=0;m<gs;m++)
  {
    int mn = n*gs+m;
    float x1 = grid[6*mn+0]; float y1 = grid[6*mn+1]; float z1 = grid[6*mn+2];
    if (fabs(z1)>1.e-5)
    if (fabs(z1)<8. || mode>1)
    if (fabs(x1)<1.e-8 && fabs(y1)<1.e-8)
    {
      printf("  %9.5f  %8.5f  %8.5f \n",z1,rho[mn],vc[mn]);
    }
  }

  return;
}

void print_vc_shift(int natoms, int gs, double* grid, float* rho, double* vc, int mode)
{
  printf("        z        rho       vc \n");

  if (natoms==1)
  for (int n=0;n<natoms;n++)
  for (int m=0;m<gs;m++)
  {
    int mn = n*gs+m;
    float x1 = grid[6*mn+0]; float y1 = grid[6*mn+1]; float z1 = grid[6*mn+2];
    if (fabs(z1)>1.e-5)
    if (fabs(x1)<1.e-8 && fabs(y1)<1.e-8)
    {
      if (z1> 80.) z1 -= 95.;
      if (z1>0.)
        printf("  %9.5f  %8.5f  %8.5f \n",z1,rho[mn],vc[mn]);
    }
  }
  else
  for (int n=0;n<natoms;n++)
  for (int m=0;m<gs;m++)
  {
    int mn = n*gs+m;
    float x1 = grid[6*mn+0]; float y1 = grid[6*mn+1]; float z1 = grid[6*mn+2];
    if (fabs(z1)>1.e-5)
    if (fabs(z1)<8. || mode>1)
    if (fabs(x1)<1.e-8 && fabs(y1)<1.e-8)
    {
      if (z1> 80.) z1 -= 95.;
      if (z1<-80.) z1 += 95;
      printf("  %9.5f  %8.5f  %8.5f \n",z1,rho[mn],vc[mn]);
    }
  }

  return;
}

void print_vc_shift(int natoms, int gs, float* grid, double* rho, double* vc, int mode)
{
  printf("        z        rho       vc \n");

  if (natoms==1)
  for (int n=0;n<natoms;n++)
  for (int m=0;m<gs;m++)
  {
    int mn = n*gs+m;
    float x1 = grid[6*mn+0]; float y1 = grid[6*mn+1]; float z1 = grid[6*mn+2];
    if (fabs(z1)>1.e-5)
    if (fabs(x1)<1.e-8 && fabs(y1)<1.e-8)
    {
      if (z1> 80.) z1 -= 95.;
      if (z1>0.)
        printf("  %9.5f  %8.5f  %8.5f \n",z1,rho[mn],vc[mn]);
    }
  }
  else
  for (int n=0;n<natoms;n++)
  for (int m=0;m<gs;m++)
  {
    int mn = n*gs+m;
    float x1 = grid[6*mn+0]; float y1 = grid[6*mn+1]; float z1 = grid[6*mn+2];
    if (fabs(z1)>1.e-5)
    if (fabs(z1)<8. || mode>1)
    if (fabs(x1)<1.e-8 && fabs(y1)<1.e-8)
    {
      if (z1> 80.) z1 -= 95.;
      if (z1<-80.) z1 += 95;
      printf("  %9.5f  %8.5f  %8.5f \n",z1,rho[mn],vc[mn]);
    }
  }

  return;
}

void print_vc_shift(int natoms, int gs, double* grid, double* rho, double* vc, int mode)
{
  printf("        z        rho       vc \n");

  if (natoms==1)
  for (int n=0;n<natoms;n++)
  for (int m=0;m<gs;m++)
  {
    int mn = n*gs+m;
    float x1 = grid[6*mn+0]; float y1 = grid[6*mn+1]; float z1 = grid[6*mn+2];
    if (fabs(z1)>1.e-5)
    if (fabs(x1)<1.e-8 && fabs(y1)<1.e-8)
    {
      if (z1> 80.) z1 -= 95.;
      if (z1>0.)
        printf("  %9.5f  %8.5f  %8.5f \n",z1,rho[mn],vc[mn]);
    }
  }
  else
  for (int n=0;n<natoms;n++)
  for (int m=0;m<gs;m++)
  {
    int mn = n*gs+m;
    float x1 = grid[6*mn+0]; float y1 = grid[6*mn+1]; float z1 = grid[6*mn+2];
    if (fabs(z1)>1.e-5)
    if (fabs(z1)<8. || mode>1)
    if (fabs(x1)<1.e-8 && fabs(y1)<1.e-8)
    {
      if (z1> 80.) z1 -= 95.;
      if (z1<-80.) z1 += 95;
      printf("  %9.5f  %8.5f  %8.5f \n",z1,rho[mn],vc[mn]);
    }
  }

  return;
}

void print_wt(int natoms, int gs, float* grid, float* rho, float* wt, int mode)
{
  printf("        z        rho       wt \n");
  for (int n=0;n<natoms;n++)
  for (int m=0;m<gs;m++)
  {
    int mn = n*gs+m;
    float x1 = grid[6*mn+0]; float y1 = grid[6*mn+1]; float z1 = grid[6*mn+2];
    if (fabs(z1)<8. || mode>1)
    if (fabs(x1)<1.e-8 && fabs(y1)<1.e-8)
    {
      printf("  %9.5f  %8.5f  %8.5e \n",z1,rho[mn],wt[mn]);
    }
  }

  return;
}

void print_vc(int natoms, int gs, float* grid, float* rho, double* vc, int mode)
{
  printf("        z        rho       vc \n");
  for (int n=0;n<natoms;n++)
  for (int m=0;m<gs;m++)
  {
    int mn = n*gs+m;
    float x1 = grid[6*mn+0]; float y1 = grid[6*mn+1]; float z1 = grid[6*mn+2];
    if (fabs(z1)>1.e-5)
    if (fabs(z1)<8. || mode>1)
    if (fabs(x1)<1.e-8 && fabs(y1)<1.e-8)
    {
      printf("  %9.5f  %8.5f  %8.5f \n",z1,rho[mn],vc[mn]);
    }
  }

  return;
}

void print_vc_shift(int natoms, int gs, float* grid, float* rho, float* vc, int mode)
{
  printf("        z        rho       vc \n");
  if (natoms==1)
  for (int n=0;n<natoms;n++)
  for (int m=0;m<gs;m++)
  {
    int mn = n*gs+m;
    float x1 = grid[6*mn+0]; float y1 = grid[6*mn+1]; float z1 = grid[6*mn+2];
    if (fabs(z1)>1.e-5)
    if (fabs(z1)<8. || mode>1)
    if (fabs(x1)<1.e-8 && fabs(y1)<1.e-8)
    {
      if (z1> 80.) z1 -= 95.;

      if (z1>0.)
        printf("  %9.5f  %8.5f  %8.5f \n",z1,rho[mn],vc[mn]);
    }
  }
  else
  for (int n=0;n<natoms;n++)
  for (int m=0;m<gs;m++)
  {
    int mn = n*gs+m;
    float x1 = grid[6*mn+0]; float y1 = grid[6*mn+1]; float z1 = grid[6*mn+2];
    if (fabs(z1)>1.e-5)
    if (fabs(z1)<8. || mode>1)
    if (fabs(x1)<1.e-8 && fabs(y1)<1.e-8)
    {
      if (z1> 80.) z1 -= 95.;
      if (z1<-80.) z1 += 95;

      printf("  %9.5f  %8.5f  %8.5f \n",z1,rho[mn],vc[mn]);
    }
  }

  return;
}

void print_vc_shift(int natoms, int gs, float* grid, float* rho, double* vc, int mode)
{
  printf("        z        rho       vc \n");
  for (int n=0;n<natoms;n++)
  for (int m=0;m<gs;m++)
  {
    int mn = n*gs+m;
    float x1 = grid[6*mn+0]; float y1 = grid[6*mn+1]; float z1 = grid[6*mn+2];
    if (fabs(z1)>1.e-5)
    if (fabs(z1)<8. || mode>1)
    if (fabs(x1)<1.e-8 && fabs(y1)<1.e-8)
    {
      if (z1> 80.) z1 -= 95.;
      if (z1<-80.) z1 += 95;

      printf("  %9.5f  %8.5f  %8.5f \n",z1,rho[mn],vc[mn]);
    }
  }

  return;
}

void print_vc_shift(int natoms, int gs, float* grid, double* rho, float* vc, int mode)
{
  printf("        z        rho       vc \n");
  for (int n=0;n<natoms;n++)
  for (int m=0;m<gs;m++)
  {
    int mn = n*gs+m;
    float x1 = grid[6*mn+0]; float y1 = grid[6*mn+1]; float z1 = grid[6*mn+2];
    if (fabs(z1)>1.e-5)
    if (fabs(z1)<8. || mode>1)
    if (fabs(x1)<1.e-8 && fabs(y1)<1.e-8)
    {
      if (z1> 80.) z1 -= 95.;
      if (z1<-80.) z1 += 95;

      printf("  %9.5f  %8.5f  %8.5f \n",z1,rho[mn],vc[mn]);
    }
  }

  return;
}


void print_vxc(int nrad, int nang, int natoms, float* grid, double* vxc)
{
  int gs = nrad*nang;

  int wg = 0;
  if (natoms<=2)
  {
    printf("       z        vxc \n");
    for (int n=0;n<natoms;n++)
    for (int m=0;m<gs;m++)
    {
      if (wg%2==0)
      {
        int mn = n*gs+m;
        float x1 = grid[6*mn+0]; float y1 = grid[6*mn+1]; float z1 = grid[6*mn+2];
        if (x1==0. && y1==0.)
        {
          printf("  %9.5f   %8.5f \n",z1,vxc[mn]);
        }
      }
      wg++;
    }
  }
  else
  {
    printf("       z        vxch \n");
    for (int n=0;n<natoms;n++)
    for (int m=0;m<gs;m++)
    {
      int mn = n*gs+m;
      float x1 = grid[6*mn+0];
      float y1 = grid[6*mn+1];
      float z1 = grid[6*mn+2];
      //float r1 = grid[6*mn+3];
      if (x1==0.f && y1==0.f)
        printf("  %9.5f  %8.5f \n",z1,vxc[mn]);
    }
  }
  return;
}

void print_vxc(int nrad, int nang, int natoms, double* grid, double* vxc)
{
  int gs = nrad*nang;
  int gsa = natoms*gs;
  int gs6 = 6*gsa;

  float gridf[gs6];
  for (int m=0;m<gs6;m++)
    gridf[m] = grid[m];
  print_vxc(nrad,nang,natoms,gridf,vxc);

  return;
}

void print_vxc(int nrad, int nang, int natoms, float* grid, float* vxcf)
{
  int gs = nrad*nang;
  int gsa = natoms*gs;

  double vxc[gsa];
  for (int m=0;m<gsa;m++)
    vxc[m] = vxcf[m];
  print_vxc(nrad,nang,natoms,grid,vxc);

  return;
}

void print_sphere(int nrad, int nang, int natoms, float* grid, double* vxc)
{
  return;
  if (natoms>1) return;

  printf(" vxc at each radial value \n");

  //int gs = nrad*nang;

  for (int m=0;m<nrad;m++)
  for (int n=0;n<nang;n++)
  {
    int i1 = m*nang+n;
    float r1 = grid[6*i1+3];
    printf(" %9.5f %8.5f \n",r1,vxc[i1]);
  }
}

void print_axes(int nrad, int nang, int natoms, float* grid, float* r1, float* rho, float* T, float* gf, float* vxch)
{
  int gs = nrad*nang;

  if (natoms<=2)
  {
    printf("       z        rho        T         gf      vxc(h) \n");
    for (int n=0;n<natoms;n++)
    for (int m=0;m<gs;m++)
    {
      int mn = n*gs+m;
      float x1 = grid[6*mn+0]; float y1 = grid[6*mn+1]; float z1 = grid[6*mn+2];
      if (fabs(z1)<8. || natoms>1)
      if (fabs(x1)<1.e-8 && fabs(y1)<1.e-8)
      {
        printf("  %9.5f  %8.5f  %8.5f  %8.5f  %8.5f \n",z1,rho[mn],T[mn],gf[mn],vxch[mn]);
      }
    }
  }
  else
  {
    printf("       r        rho        T         gf      vxch \n");
    for (int n=0;n<natoms;n++)
    for (int m=0;m<gs;m++)
    {
      int mn = n*gs+m;
      float x1 = grid[6*mn+0];
      float y1 = grid[6*mn+1];
      float z1 = grid[6*mn+2];
      //float r1 = grid[6*mn+3];
      float rho1 = fabs(rho[mn]);
      //if (fabs(z1)<8.)
      if (fabs(z1)<16.)
      if (rho1>1.e-5 && x1==0.f && y1==0.f)
        printf("  %9.5f  %8.5f  %8.5f  %8.5f  %8.5f \n",z1,rho[mn],T[mn],gf[mn],vxch[mn]);
    }
  }

 #if 0
  int oat = gs*0; //which atom's grid
  int offset = 4;
  for (int m=0;m<nrad;m++)
    printf(" %8.5f",r1[m]);
  printf("\n");
  for (int m=0;m<nrad;m++)
    printf(" %8.5f",rho_wf[oat+m*nang+offset]);
  printf("\n");
  for (int m=0;m<nrad;m++)
    printf(" %8.5f",T_wf[oat+m*nang+offset]);
  printf("\n");
  for (int m=0;m<nrad;m++)
    printf(" %8.5f",vxch[oat+m*nang+offset]);
  printf("\n");

  printf("  grid points: \n");
  for (int m=0;m<nrad;m++)
  {
    int wg = oat + m*nang+offset;
    printf(" %8.5f %8.5f %8.5f \n",grid[6*wg],grid[6*wg+1],grid[6*wg+2]);
  }
 #endif

  return;
}

void print_axes(int nrad, int nang, int natoms, float* grid, float* r1, float* rho, float* T, float* gf, double* vxc)
{
  int gs = nrad*nang;
  int gsa = gs*natoms;
  float* vxcf = new float[gsa];
  for (int m=0;m<gsa;m++)
    vxcf[m] = vxc[m];

  print_axes(nrad,nang,natoms,grid,r1,rho,T,gf,vxcf);

  delete [] vxcf;

  return;
}
