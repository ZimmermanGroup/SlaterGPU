#include "grid_util.h"
#include "lebedev2.h"
#include <cstdio>
#include <algorithm>

using namespace std;

void compare_pao_12(int ngpu, bool gbasis, int natoms, int* atno, double* coords, int nrad, int nang, double* ang_g, double* ang_w, vector<vector<double> >& basis, cusolverDnHandle_t cu_hdl)
{
  int gs = nrad*nang;
  int gsa = natoms*gs;
  int gsa6 = gsa*6;
  int dim = gsa;

  int N = basis.size();
  int N2 = N*N;

  double* Pao1 = new double[N2];
  double* Pao2 = new double[N2];
  bool found1 = read_square(N,Pao1,"Pao_1");
  bool found2 = read_square(N,Pao2,"Pao_2");
  if (!found1 || !found2) { printf("\n\n  couldn't find Pao_1 or Pao_2 file \n"); return; }

  double* rho1 = new double[gsa];
  double* rho2 = new double[gsa];
  double* grid = new double[gsa6];
  float* gridf = new float[gsa6];
  double* wt = new double[gsa];

  printf("\n gbasis: %i natoms: %i nrad: %3i nang: %4i gsa: %4i \n",(int)gbasis,natoms,nrad,nang,gsa);

  #pragma acc enter data create(rho1[0:gsa],rho2[0:gsa])
  #pragma acc enter data create(grid[0:gsa6],wt[0:gsa],gridf[0:gsa6])
  #pragma acc enter data copyin(coords[0:3*natoms],ang_g[0:3*nang],ang_w[0:nang])

  printf("\n comparing Pao_1 to Pao_2 \n");

  if (read_int("GRID_SQUARE"))
  {
    int nsq = 1;
    for (int nq=0;nq<1000;nq++)
    {
      int gsize = nq*nq*nq;
      if (gsize>gsa) { nsq = nq-1; break; }
    }
    int gsize = nsq*nsq*nsq;
    printf("\n using square grid (size: %4i) \n",gsize);
    double a = 3;
    double volume = a*a*a;
    double sp = a/(nsq+1);
    double wt1 = volume/gsize;
    for (int j=0;j<gsize;j++)
      wt[j] = wt1;

    int ng = 0;
    int nsq2 = nsq/2;
    for (int i=0;i<nsq;i++)
    {
      double x1 = sp*(i-nsq2);
      for (int j=0;j<nsq;j++)
      {
        double y1 = sp*(j-nsq2);
        for (int k=0;k<nsq;k++)
        {
          double z1 = sp*(k-nsq2);

          grid[6*ng+0] = x1;
          grid[6*ng+1] = y1;
          grid[6*ng+2] = z1;
          ng++;
        }
      }
    }
    for (int i=6*ng;i<gsa6;i++)
      grid[i] = 0.;
    for (int i=ng;i<gsa;i++)
      wt[i] = 0.;
    dim = gsize;

    #pragma acc update device(grid[0:gsa6],wt[0:gs])
  }
  else
    get_becke_grid_full(natoms,atno,coords,nrad,nang,ang_g,ang_w,6,grid,wt);

  if (gbasis)
  {
    for (int j=0;j<gsa6;j++) gridf[j] = grid[j];
    #pragma acc update device(gridf[0:gsa6])

    compute_rho(1,natoms,atno,coords,basis,Pao1,nrad,dim,gridf,rho1,NULL,1);
    compute_rho(1,natoms,atno,coords,basis,Pao2,nrad,dim,gridf,rho2,NULL,1);
  }
  else
  {
    compute_rhod(natoms,atno,coords,basis,Pao1,nrad,dim,grid,rho1,NULL,NULL,1);
    compute_rhod(natoms,atno,coords,basis,Pao2,nrad,dim,grid,rho2,NULL,NULL,1);
  }

  #pragma acc update self(rho1[0:gsa],rho2[0:gsa],grid[0:gsa6],wt[0:gsa])

  if (natoms==1)
  {
    printf("       r        rho1      rho2 \n");
    for (int j=0;j<nrad;j++)
      printf("  %8.5f  %8.5f \n",grid[6*j*nang+3],rho1[j*nang],rho2[j*nang]);
    printf("\n");
  }

  double diff1 = 0.;
  double diff2 = 0.;
 #pragma acc parallel loop present(rho1[0:gsa],rho2[0:gsa],wt[0:gsa]) reduction(+:diff1,diff2)
  for (int j=0;j<gsa;j++)
  {
    double d12 = rho1[j]-rho2[j];
    diff1 += abs(d12)*wt[j];
    diff2 += d12*d12*wt[j];
  }

  diff2 = sqrt(diff2);
  printf("\n L1/L2: %10.6f %10.6f \n",diff1,diff2);

  #pragma acc exit data delete(rho1[0:gsa],rho2[0:gsa],grid[0:gsa6],wt[0:gsa],gridf[0:gsa],coords[0:3*natoms],ang_g[0:3*nang],ang_w[0:nang])

  delete [] rho1;
  delete [] rho2;
  delete [] grid;
  delete [] gridf;
  delete [] wt;
  delete [] Pao1;
  delete [] Pao2;

  return;
}

void save_grid_rho(bool gbasis, int natoms, int* atno, double* coords, int nrad, int nang, double* ang_g, double* ang_w, vector<vector<double> >& basis)
{
  int gs = nrad*nang;
  int gsa = natoms*gs;
  int gsa6 = gsa*6;

  int N = basis.size();
  int N2 = N*N;

  double* Pao = new double[N2];
  bool found = read_square(N,Pao,"Pao_ci");
  if (!found) read_square(N,Pao,"Pao");
  if (!found) { printf("  couldn't find Pao or Pao_ci file \n"); return; }

  double* rho = new double[gsa];
  double* grid = new double[gsa6];
  float* gridf = new float[gsa6];
  double* wt = new double[gsa];

  printf("\n gbasis: %i natoms: %i nrad: %3i nang: %4i gsa: %4i \n",(int)gbasis,natoms,nrad,nang,gsa);

  #pragma acc enter data create(rho[0:gsa])
  #pragma acc enter data create(grid[0:gsa6],wt[0:gsa],gridf[0:gsa6])
  #pragma acc enter data copyin(coords[0:3*natoms],ang_g[0:3*nang],ang_w[0:nang])

  get_becke_grid_full(natoms,atno,coords,nrad,nang,ang_g,ang_w,6,grid,wt);
  if (gbasis)
  {
    //for (int j=0;j<gsa6;j++) gridf[j] = grid[j];
    //#pragma acc update device(gridf[0:gsa6])

    compute_rhodg(1,natoms,atno,coords,basis,Pao,nrad,gsa,grid,rho,NULL,1);
  }
  else
    compute_rhod(natoms,atno,coords,basis,Pao,nrad,gsa,grid,rho,NULL,NULL,1);

  #pragma acc update self(rho[0:gsa],grid[0:gsa6],wt[0:gsa])

  if (natoms==1)
  {
    printf("       r        rho \n");
    for (int j=0;j<nrad;j++)
      printf("  %8.5f  %8.5f \n",grid[6*j*nang+3],rho[j*nang]);
    printf("\n");
  }

  double dent = 0.;
 #pragma acc parallel loop present(rho[0:gsa],wt[0:gsa]) reduction(+:dent)
  for (int j=0;j<gsa;j++)
    dent += rho[j]*wt[j];
  printf("\n total e-s: %8.5f \n",dent);

  printf("\n  writing RHO_WF to disk \n");
  write_vector(gsa,rho,"RHO_WF");
  printf("\n  writing GRID_WTS to disk \n");
  write_grid(natoms,nrad,nang,grid,wt);

  #pragma acc exit data delete(rho[0:gsa],grid[0:gsa6],wt[0:gsa],gridf[0:gsa],coords[0:3*natoms],ang_g[0:3*nang],ang_w[0:nang])

  delete [] rho;
  delete [] grid;
  delete [] gridf;
  delete [] wt;
  delete [] Pao;

  return;
}

void save_grid_ao_basis(bool gbasis, int natoms, int* atno, double* coords, int nrad, int nang, double* ang_g, double* ang_w, vector<vector<double> >& basis)
{
  int gs = nrad*nang;
  int gsa = natoms*gs;
  int gsa6 = gsa*6;

  int N = basis.size();
  int N2 = N*N;

  double* Pao = new double[N2];
  double* rho = new double[gsa];
  double* grid = new double[gsa6];
  double* wt = new double[gsa];

  printf("\n gbasis: %i natoms: %i nrad: %3i nang: %4i gsa: %4i \n",(int)gbasis,natoms,nrad,nang,gsa);

  #pragma acc enter data create(rho[0:gsa])
  #pragma acc enter data create(grid[0:gsa6],wt[0:gsa])
  #pragma acc enter data copyin(coords[0:3*natoms],ang_g[0:3*nang],ang_w[0:nang])

  get_becke_grid_full(natoms,atno,coords,nrad,nang,ang_g,ang_w,6,grid,wt);

  for (int i=0;i<N;i++)
  for (int j=0;j<=i;j++)
  {
    for (int k=0;k<N2;k++) Pao[k] = 0.;
    Pao[i*N+j] = Pao[j*N+i] = 0.5;
    if (i==j) Pao[i*N+j] = 1.;

    if (gbasis)
      compute_rhodg(1,natoms,atno,coords,basis,Pao,nrad,gsa,grid,rho,NULL,1);
    else
      compute_rhod(natoms,atno,coords,basis,Pao,nrad,gsa,grid,rho,NULL,NULL,1);
    #pragma acc update self(rho[0:gsa])

    double dent = 0.;
   #pragma acc parallel loop present(rho[0:gsa],wt[0:gsa]) reduction(+:dent)
    for (int k=0;k<gsa;k++)
      dent += rho[k]*wt[k];
    printf("  %2i-%2i total e-s: %8.5f \n",i,j,dent);

    string filename = "RHO_AO_"+to_string(i+1)+"_"+to_string(j+1);
    write_vector(gsa,rho,filename);
  }

  #pragma acc exit data delete(rho[0:gsa],grid[0:gsa6],wt[0:gsa],coords[0:3*natoms],ang_g[0:3*nang],ang_w[0:nang])

  delete [] rho;
  delete [] grid;
  delete [] wt;
  delete [] Pao;

  return;
}

