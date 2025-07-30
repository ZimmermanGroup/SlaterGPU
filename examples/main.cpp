#if USE_ACC
#include "accel.h"
#endif
#include <omp.h>
#include <cstdio>
#include <vector>
#include <string>
#include <chrono>

#include "read.h"
#include "write.h"
#include "lebedev2.h"
#include "integrals.h"
#include "grid_util.h"

#include "prosph.h"
#include "quad.h"

#include "../include/cintwrapper.h"
#include "../include/cintprep.h"

using namespace std;



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

void compute_integrals_g(int natm, int nbas, int nenv, int N, int Naux, int nbas_ri, int* atm, int* bas, double* env, double* S, double* T, double* jH1, double* A, double* C, int prl)
{
  get_overlap(S, N, natm, nbas, nenv, atm, bas, env);
  get_hcore(jH1, N, natm, nbas, nenv, atm, bas, env);
  if (T!=NULL)
    get_tcore(T, N, natm, nbas, nenv, atm, bas, env);

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
  for (int m=0;m<N2;m++) En[m] = 0.;

  compute_integrals_g(natoms,nbas,nenv,N,Naux,nbas_ri,atm,bas,env,S,NULL,jH1,A,C,1);

 //note we don't have separate En and T, so using sum of the two
  for (int m=0;m<N2;m++) T[m] = jH1[m];
  write_S_En_T(N,S,En,T);
  write_square(Naux,A,"A",2);
  write_C(Naux,N2,C);

  return;
}


void compute_ps_integrals_to_disk(int natoms, int* atno, double* coords, vector<vector<double> > basis, vector<vector<double> > basis_aux, int prl)
{
 //prolate spheroidal coordinates
  if (prl > 0) printf("Computing prolate spheroidal integrals:\n");

  int N = basis.size();
  int N2 = N*N;
  int Naux = basis_aux.size();
  int Naux2 = Naux*Naux;
  int N2a = N2*Naux;

  int nmu, nnu, nphi;
  read_gridps(nmu,nnu,nphi,1);
  int nquad = 8; int nquad2 = 8;
  read_quad(nquad,nquad2);

  int nsplit = 1;
  int nsplit_read = read_int("NSPLIT");
  if (nsplit_read>0) nsplit = nsplit_read;

  double gamma = 0.5;
  double gm_read = read_float("GAMMA");
  if (gm_read>0.) gamma = gm_read;

  int nbatch = 1;
  int nbatch_read = read_int("NBATCH");
  if (nbatch_read>0) nbatch = nbatch_read;

 #if 0
  nmu = 8;
  nnu = 4;
  nphi = 4;
 #endif

  printf("  PS grid size: %2i %2i %2i quad_order: %2i %2i  w/quad: %6i \n",nmu,nnu,nphi,nquad,nquad2,nmu*nnu*nphi*nquad*nquad*nquad);

  double Ssp[N2];
  double Tsp[N2];
  double Ensp[N2];
  double pVpsp[N2];
  double* Asp = new double[Naux2]();
  double* Aysp = new double[Naux2]();
  double* Csp = new double[N2a]();
  double* Cysp = new double[N2a]();
  for (int j=0;j<N2;j++) Ssp[j] = Tsp[j] = Ensp[j] = pVpsp[j] = 0.;

  compute_STEn_ps(natoms,atno,coords,basis,nquad,nmu,nnu,nphi,Ssp,Tsp,Ensp,prl);
  compute_pVp_ps(natoms,atno,coords,basis,nquad,nmu,nnu,nphi,pVpsp,prl);
  compute_pVp_3c_ps(natoms,atno,coords,basis,nquad,nquad2,nsplit,nmu,nnu,nphi,pVpsp,prl);
  compute_2c_ps(0,0,gamma,natoms,atno,coords,basis_aux,nquad,nmu,nnu,nphi,Asp,prl);
  compute_3c_ps(0,0,gamma,nbatch,natoms,atno,coords,basis,basis_aux,nquad,nquad2,nsplit,nmu,nnu,nphi,Ensp,Csp,prl);

  if (prl > 0) printf("Printing PS Integral Files:\n");
  write_S_En_T(N,Ssp,Ensp,Tsp);
  write_square(N,pVpsp,"pVp",2);
  write_square(Naux,Asp,"A",2);
  write_C(Naux,N2,Csp);

  delete [] Asp;
  delete [] Csp;

  return;
}


int main(int argc, char* argv[]) {
  printf("Running test code for SlaterGPU\n");
  int nomp = 1;
  #pragma omp parallel
  nomp = omp_get_num_threads();

  int ngpu = 0;
  #if USE_ACC
  ngpu = acc_get_num_devices(acc_device_nvidia);
  if (ngpu != nomp) {
    ngpu = min(nomp,ngpu);
    nomp = ngpu;
  }
  acc_init(acc_device_nvidia);
  printf("This is a GPU run with %d GPUs found\n",ngpu);
  #else
  printf(" This is a CPU run\n");
  #endif

  int nrad = 8;
  int nang = 4;
  read_nrad_nang(nrad,nang,1);

  int charge = 0;
  int nup = 0;
  int * atno = new int[100]();
  double * coords;
  vector< vector< double > > basis;
  vector< vector< double > > basis_aux;
  double Enn = 0;

  bool gbasis = read_basis();

  // seting up
  int natoms = initialize(gbasis,basis,basis_aux,atno,coords,charge,nup,Enn,1);
  float * coordsf = new float[3*natoms];
  for (int i = 0; i < 3*natoms; i++) {
    coordsf[i] = coords[i];
  }

  int Na, Nb, Nc;
  int No = electron_count(charge,natoms,atno,Nc,Na,Nb);

  printf(" There are %d occupied orbitals (%d core)\n",No,Nc);
  printf(" Na: %d Nb: %d\n",Na,Nb);

  int N = basis.size();
  int Naux = basis_aux.size();

  int nbas,nenv,nbas_ri;
  int* atm; int* bas; double* env;
  if (gbasis)
  {
    basis = setup_integrals_gsgpu(basis_aux,natoms,atno,coords,nbas,nenv,N,Naux,nbas_ri,atm,bas,env,1);

    printf("  Gaussian basis sizes, N: %2i Naux: %3i \n",N,Naux);
    printf("  basis size: %2i \n",basis.size());
  }

  int N2 = N*N;
  int Naux2 = Naux*Naux;
  int N2a = N2*Naux;

  int size_ang = order_table(nang);
  size_t gs = nrad * size_ang;
  printf(" nrad: %d nang: %d size_ang: %d\n",nrad,nang,size_ang);
  printf(" grid size: %u\n",gs);

  double * ang_g = new double[3*size_ang];
  double * ang_w = new double[size_ang];
  get_angular_grid(size_ang,ang_g,ang_w);

  // Doing work
  if (true) {
    double * A = new double[Naux2]();
    double * Anorm = new double[Naux];
    double * C = new double[N2a]();
    double * S = new double[N2]();
    double * En = new double[N2];
    double * T = new double[N2];
    double * pVp = new double[N2];

    //Setup for VdV and 4c integral tests

    double* Pao = new double[N2]();
    //need Pao file!!
    Pao[0] = 1.;
    read_square(N,Pao,"Pao");

    int nc = 1;
    float* coordsc = new float[3*nc];
    coordsc[0] = 0.; coordsc[1] = 0.; coordsc[2] = coordsf[2]+20.;

    double* V = new double[nc];
    double* dV = new double[3*nc];

    double* g = new double[N2*N2]();

    int prl = 1;
    int c4 = 0;

    #pragma acc enter data create(A[0:Naux2],C[0:N2a],S[0:N2],En[0:N2],T[0:N2],pVp[0:N2])
    #pragma acc enter data create(V[0:nc],dV[0:3*nc],Pao[0:N2],coordsc[0:3*nc])
    #pragma acc enter data create(g[0:N2*N2])
    printf("1e ints: %d\n2c2e ints: %d\n3c3e ints: %d\n4c4e ints: %d\n",N2, Naux2, N2a, N2*N2);

    if (gbasis)
    {
      printf("\n TESTING gaussian integrals \n");
      compute_gaussian_integrals_to_disk(N,Naux,natoms,nbas,nenv,nbas_ri,atm,bas,env);
    }
    else if (check_PS() > 0)
      compute_ps_integrals_to_disk(natoms,atno,coords,basis,basis_aux,prl);
    else
    {
      printf("Computing Standard Integrals:\n");

      auto t1 = chrono::high_resolution_clock::now();
      compute_ST(natoms, atno, coordsf, basis, nrad, size_ang, ang_g, ang_w, S, T, prl);

      auto t2 = chrono::high_resolution_clock::now();

      bool do_2c_v1 = read_int("DO_2c_V1");
      if (do_2c_v1)
      {
        if (prl > 0) printf("  using compute_all_2c instead of compute_all_2c_v2 \n");
        //compute_all_2c(natoms,atno,coordsf,basis_aux,nrad,size_ang,ang_g,ang_w,A,prl);
      }
      else
      {
        if (prl > 0) printf("  using compute_all_2c_v2 \n");
        compute_all_2c_v2(0,natoms,atno,coordsf,basis_aux,nrad,size_ang,ang_g,ang_w,A,prl);
      }

      auto t3 = chrono::high_resolution_clock::now();
      auto t4 = chrono::high_resolution_clock::now();
      if (ngpu == 1)
      {
        compute_Enp(natoms,atno,coordsf,basis,nrad,size_ang,ang_g,ang_w,En,pVp,prl);
        auto t4 = chrono::high_resolution_clock::now();

        bool do_3c_v1 = read_int("DO_3C_V1");
        if (do_3c_v1)
        {
          if (prl > 0) printf("  using compute_all_3c instead of compute_all_3c_v2 \n");
          //compute_all_3c(natoms,atno,coordsf,basis,basis_aux,nrad,size_ang,ang_g,ang_w,C,prl);
        }
        else
          if (prl > 0) printf("  using compute_all_3c_v2 \n");
          compute_all_3c_v2(0,natoms,atno,coordsf,basis,basis_aux,nrad,size_ang,ang_g,ang_w,C,prl);
      }
      else
      {
        compute_Enp_para(ngpu,natoms,atno,coordsf,basis,nrad,size_ang,ang_g,ang_w,En,pVp,prl);
        auto t4 = chrono::high_resolution_clock::now();

        compute_all_3c_para(ngpu,0,natoms,atno,coordsf,basis,basis_aux,nrad,size_ang,ang_g,ang_w,C,prl);
      }

      auto t5 = chrono::high_resolution_clock::now();
      compute_VdV(natoms,atno,coordsf,basis,nrad,size_ang,ang_g,ang_w,nc,coordsc,Pao,V,dV,prl);

      auto t6 = chrono::high_resolution_clock::now();
      if (c4 > 0)
        compute_all_4c_v2(natoms,atno,coordsf,basis,nrad,size_ang,ang_g,ang_w,g,prl);

      auto t7 = chrono::high_resolution_clock::now();

      auto elapsed12 = chrono::duration_cast<chrono::nanoseconds>(t2-t1).count();
      auto elapsed23 = chrono::duration_cast<chrono::nanoseconds>(t3-t2).count();
      auto elapsed34 = chrono::duration_cast<chrono::nanoseconds>(t4-t3).count();
      auto elapsed45 = chrono::duration_cast<chrono::nanoseconds>(t5-t4).count();
      auto elapsed56 = chrono::duration_cast<chrono::nanoseconds>(t6-t5).count();
      auto elapsed67 = chrono::duration_cast<chrono::nanoseconds>(t7-t6).count();
      auto elapsed17 = chrono::duration_cast<chrono::nanoseconds>(t7-t1).count();

      printf("-------------------------------\n");
      printf("Integral ST   time: %5.3e s\n",(double)elapsed12/1.e9);
      printf("Integral 2c2e time: %5.3e s\n",(double)elapsed23/1.e9);

      if (nomp == 1)
      {
        printf("Integral Vne  time: %5.3e s\n",(double)elapsed34/1.e9);
        printf("Integral 3c2e time: %5.3e s\n",(double)elapsed45/1.e9);
      }
      else
      {
        printf("Integral Vne  (para)  time: %5.3e s\n",(double)elapsed34/1.e9);
        printf("Integral 3c2e (para)  time: %5.3e s\n",(double)elapsed45/1.e9);
      }

      printf("Integral VdV  time: %5.3e s\n",(double)elapsed56/1.e9);
      if (c4 > 0)
        printf("Integral 4c (v2)    time: %5.3e s\n",(double)elapsed67/1.e9);

      printf("-------------------------------\n");
      printf("Integral tot  time: %5.3e s\n",(double)elapsed17/1.e9);
      printf("-------------------------------\n");

      #pragma acc exit data copyout(A[0:Naux2],C[0:N2a],S[0:N2],En[0:N2],T[0:N2],pVp[0:N2])
      #pragma acc exit data copyout(V[0:nc],dV[0:3*nc],Pao[0:N2],coordsc[0:3*nc],g[0:N2*N2])

      if (prl > 0) printf("Printing Standard Integral Files:\n");
      write_S_En_T(N,S,En,T);
      write_square(Naux,A,"A",2);
      write_square(N,pVp,"pVp",2);
      write_C(Naux, N2, C);
      write_square(N2,g,"g",2);
    }

    delete [] A, Anorm, C, S, En, T, pVp;
    delete [] V, dV, Pao, coordsc;

    delete [] g;
  }

  delete [] ang_g, ang_w;

  printf("\n  Done!\n");
  return 0;

}
