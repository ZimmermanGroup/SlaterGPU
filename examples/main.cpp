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
#include "scf_util.h"
#include "prosph.h"
#include "quad.h"
#include "jellium_ints.h"


using namespace std;


void save_mog(int natoms, int* atno, double* coords, int nrad, int gsa, vector<vector<double> > basis,
              int No, double* grid, double* wt, int prl)
{
  if (nrad<1) nrad = 1; 
  if (natoms>3) return;

  //assumes jCA on gpu

  //int gs = gsa/natoms;
  //int gsa3 = 3*gsa;
  int gsa6 = 6*gsa;

  printf("\n getting MOs on the grid (size: %6i) \n",gsa);
  if (gsa<1) { exit(-1); }

  int nang = gsa/nrad/natoms;
  //#include "slatergpu/jsetup.cpp"
  int gs2 = gsa; 

  int N = basis.size();
  int N2 = N*N; 
  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);

  int the_num = read_int("TNUM");

  float norm[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];

  int iN = N; 
  float** val1 = new float*[iN];
  for (int i=0;i<iN;i++)
    val1[i] = new float[gsa];

  float* grid1 = new float[gsa6];

  #pragma acc enter data copyin(norm[0:N])
  #pragma acc enter data create(val1[0:iN][0:gsa])
  #pragma acc enter data create(grid1[0:gsa6])

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    float Z1 = (float)atno[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2]; 

    //copy_grid(gsa,grid1,grid);
    #pragma acc parallel loop independent present(grid1[0:gsa6],grid[0:gsa6])  //need to check if this should be 6*gs or gsa6 gsa seems like it would be for Becke
    for (int j=0;j<gsa6;j++)
      grid1[j] = grid[j];
    recenter_grid_zero(gsa,grid1,-A1,-B1,-C1);
    
    #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gsa])
    for (int i1=s1;i1<s2;i1++)
    for (int j=0;j<gsa;j++)
      val1[i1][j] = 1.f; 

    for (int i1=s1;i1<s2;i1++)
    {    
      int ii1 = i1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; float zeta1 = basis1[3];

      #if 0
      if (ss_basis)
        eval_ss(ii1,gs2,grid,val1[ii1],n1,l1,m1,zeta1,Rc);
      else if (sgs_basis)
        eval_sgs(ii1,gs1,gs2,grid,val1[ii1],n1,l1,m1,zeta1,Rc);
      else 
      #endif
        eval_sh(ii1,gsa,grid1,val1[ii1],n1,l1,m1,zeta1);
    }    

  } //loop m over natoms

  //double* tm1 = new double[gsa];
  //#pragma acc enter data create(tm1[0:gsa])
  #if 0
  printf("  z:");
  if (natoms>1)
  for (int k=0;k<gsa;k++)
  {
    float x1 = grid[6*k+0]; float y1 = grid[6*k+1]; float z1 = grid[6*k+2];
    if (x1==0. && y1==0.)
      printf(" %9.5f",z1);
  }
  else 
  for (int k=0;k<gsa;k++)
  {
    float x1 = grid[6*k+0]; float y1 = grid[6*k+1]; float z1 = grid[6*k+2];
    if (x1==0. && y1==0. && z1>0.)
      printf(" %9.5f",z1);
  }
  printf("\n");
  #endif

  #pragma acc update self(wt[0:gsa],val1[0:iN][0:gsa])

  printf("\n grid: \n");
  for (int i=0;i<gsa;i++)
    printf("%12.12f ",grid[i]);

  printf("\n weights: \n");
  //for (int i=0;i<gsa;i++)
  //  printf("%12.12f ",wt[i]);
  printf("%12.12f \n",wt[100000]);

  double bnorm = 0.;
  printf("\n bnorm: \n");
  for (int i=0;i<gsa;i++)
  {
    printf("val1: %12.12f \n",val1[0][i]);
    bnorm += val1[0][i]*val1[0][i]*wt[i];
  }
  printf("%12.12f \n",bnorm);
  
//added loop ADS
  for (int k1=0;k1<the_num+1;k1++)
  {
    double* mogrid = new double[N*gsa];
    string str_it = to_string(k1);
    double* jCA = new double[N2];
    read_square(N,jCA,("CAoptL2_v"+str_it).c_str());
    //printf("CAoptL2_v%s", str_it.c_str());
    //print_square(N,jCA);
    double* tm1 = new double[gsa];
    #pragma acc enter data copyin(jCA[0:N2])
    #pragma acc enter data create(tm1[0:gsa])
    double normt = 0.;
    //MOs
    for (int i1=0;i1<N;i1++)
    {
      #pragma acc parallel loop present(tm1[0:gsa])
      for (int j=0;j<gsa;j++)
        tm1[j] = 0.f; 

      #pragma acc parallel loop present(tm1[0:gsa],val1[0:iN][0:gsa],norm[0:N],jCA[0:N2])
      for (int j=0;j<gsa;j++)
      {    
        float v1 = 0.f; 
        #pragma acc loop reduction(+:v1)
        for (int k=0;k<N;k++)
          v1 += jCA[k*N+i1]*norm[k]*val1[k][j];
        tm1[j] = v1;
      }    

      float sign1 = 1.;
      #pragma acc serial present(tm1[0:gsa])
      if (tm1[0]<0.)
        sign1 *= -1.; 

      #pragma acc update self(tm1[0:gsa])

      normt = 0.;

      for (int i3=0;i3<gsa;i3++)
      {
        normt += tm1[i3]*tm1[i3]*wt[i3];
      }

      printf("norm: %12.12f \n", normt);

      for (int i2=0;i2<gsa;i2++)
      {
        mogrid[i1*gsa+i2] = tm1[i2];
      }
    }
    //printf("norm: %12.12f \n", normt);
    #pragma acc exit data delete(tm1[0:gsa],jCA[0:N2])
    delete [] tm1; 
    delete [] jCA;

    if(k1 == 0)
    {
      printf("MOg_0: \n");
      for (int l=0;l<N;l++)
      {
        for (int l1=0;l1<gsa;l1++)
        {
          //printf("%12.12f ", mogrid[l*gsa+l1]);
        }
        //printf("\n");
      }
    }
  
    write_rect(No,gsa,mogrid,"MOg_"+str_it);

    delete [] mogrid;
  } //ADS

  //cleanup
   #pragma acc exit data delete(norm[0:N])
   #pragma acc exit data delete(val1[0:iN][0:gsa])
   #pragma acc exit data delete(grid1[0:gsa6])

   delete [] n2i; 

   delete [] grid1;

   for (int i=0;i<iN;i++)
     delete [] val1[i];
   delete [] val1;

   return;
}

void do_jellium(int natoms, int Ne, vector<vector<double> > basis, vector<vector<double> > basis_aux,
             int nrad, int nang, double* ang_g, double* ang_w, int prl)
{
  //if (natoms>1)
  //  return spheroidal_jellium(natoms,z0,Ne,basis,basis_aux,nrad,nang,ang_g,ang_w,prl,cu_hdl);

  double Rc = 1.;
  double Rc_read = read_float("RC");
  if (Rc_read>0.) Rc = Rc_read;

  bool slater = 0;
  int order = read_int("JELLIUM");
  if (order<2) slater = 1;

  printf("\n TESTING jellium functions.   order: %i \n",order);


  double Zg = read_float("ZG");
  double ztg = read_float("ZTG");
  if (ztg<0.) { printf("\n ERROR: Gaussian impurity zeta cannot be less than 0. \n"); exit(-1); }

  int No = Ne/2;
  //int Ne = 2*No;
  double rs = Rc*pow(Ne,-1./3.);

  printf("  Ne: %2i  rs: %8.5f  Rc: %8.5f  Zg: %6.3f  ztg: %5.3f \n",Ne,rs,Rc,Zg,ztg);
  print_murak_rmax(Rc,nrad,basis);

  int N = basis.size();
  int N2 = N*N;
  int Naux = basis_aux.size();
  int Na2 = Naux*Naux;
  int N2a = Naux*N2;


  //update norms in basis
  for (int i=0;i<N;i++)
  {
    int n1 = basis[i][0]; int l1 = basis[i][1]; int m1 = basis[i][2]; double zt1 = basis[i][3];
    double norm0 = basis[i][4];
    double norm1 = norm0; //norm_sgs(n1,l1,m1,zt1,Rc);
    basis[i][4] = norm1;
    //printf("  basis %2i  n1: %i  norm: %8.5f -> %8.5f \n",i,n1,norm0,norm1);
  }

  int atno[natoms];
  atno[0] = 1;
  //MD float coordsf[3]; coordsf[0] = coordsf[1] = coordsf[2] = 0.;
  double coords[3]; coords[0] = coords[1] = coords[2] = 0.;


 //need Jellium integrals in SGS basis
  double* En = new double[N2]();
  double* S = new double[N2];
  double* T = new double[N2];
  double* A = new double[Na2];
  double* C = new double[N2a];

  #pragma acc enter data create(S[0:N2],En[0:N2],T[0:N2],A[0:Na2],C[0:N2a])

  bool found_ints = read_integrals("0",N,Naux,S,T,En,A,C,prl);
  //bool found_a = read_square(Naux,A,"A");
  //if (!found_a) { printf("\n ERROR: couldn't find 2c integrals \n"); exit(-1); }
  #pragma acc update device(S[0:N2],En[0:N2],T[0:N2],A[0:Na2],C[0:N2a])

  bool have_ints = read_int("JINTS");
  if (have_ints)
  {
    bool have_A = check_file("A");
    bool have_C = check_file("Ciap");
    bool have_1 = check_file("SENT");

    if (!have_A || !have_C || !have_1) have_ints = 0;
  }
  if (!have_ints)
  {
    compute_STEn_jellium(slater,order,Zg,ztg,rs,Rc,Ne,1,basis,nrad,nang,ang_g,ang_w,S,T,En,prl);
    compute_C_jellium(slater,Rc,basis,basis_aux,nrad,nang,ang_g,ang_w,C);

    write_S_En_T(N,S,En,T);
    write_C(Naux,N2,C);
    write_int(1,"JINTS");
  }
  bool need_olintsd = !check_file("olintsd");
  int coulomb = read_int("COULOMB");
  if (need_olintsd && coulomb<2)
  {
    double* ol = new double[N2*N2];

    compute_4c_ol_jellium(slater,Rc,basis,nrad,nang,ang_g,ang_w,ol,prl);
    write_square(N2,ol,"olintsd",1);

    delete [] ol;
  }
  //#pragma acc update device(S[0:N2],En[0:N2],T[0:N2],A[0:Na2],C[0:N2a])

  printf("\n done with Jellium integrals \n\n");

// ADD THIS SECTION HERE to test compute_Exyz_slaussian
  printf("\n Testing compute_Exyz_slaussian \n");
  
  double* E = new double[3*N2]();  // Initialize to zeros
  #pragma acc enter data create(E[0:3*N2])
  
  compute_Exyz_slaussian(Rc, basis, nrad, nang, ang_g, ang_w, E, prl);
  
  #pragma acc exit data delete(E[0:3*N2])
  delete [] E;
  
  printf("\n done testing compute_Exyz_slaussian \n\n");
  // END OF ADDED SECTION

   //cleanup
  #pragma acc exit data delete(En[0:N2])
  #pragma acc enter data create(S[0:N2],T[0:N2],A[0:Na2],C[0:N2a])
  delete [] En;
  delete [] S;
  delete [] T;
  delete [] A;
  delete [] C;

  return;
}

void compute_ps_integrals_to_disk(int natoms, int* atno, double* coords, vector<vector<double> > basis, vector<vector<double> > basis_aux, int prl, int jellium)
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

  if (jellium>0)
  {
    compute_2c_ps(0,0,gamma,nbatch,natoms,atno,coords,basis_aux,nquad,nmu,nnu,nphi,Asp,prl);
    if (prl > 0) printf("Writing A to disk to be used for jellium\n");
    write_square(Naux,Asp,"A",2);
  }
  else
  {
    compute_STEn_ps(natoms,atno,coords,basis,nquad,nmu,nnu,nphi,Ssp,Tsp,Ensp,prl);
    bool x2c = read_int("X2C");
    if (x2c)
    {
      compute_pVp_ps(natoms,atno,coords,basis,nquad,nmu,nnu,nphi,pVpsp,prl);
      compute_pVp_3c_ps(natoms,atno,coords,basis,nquad,nquad2,nsplit,nmu,nnu,nphi,pVpsp,prl);
    }
    compute_2c_ps(0,0,gamma,nbatch,natoms,atno,coords,basis_aux,nquad,nmu,nnu,nphi,Asp,prl);
    compute_3c_ps(0,0,gamma,nbatch,natoms,atno,coords,basis,basis_aux,nquad,nquad2,nsplit,nmu,nnu,nphi,Ensp,Csp,prl);

    if (prl > 0) printf("Printing PS Integral Files:\n");
    write_S_En_T(N,Ssp,Ensp,Tsp);
    if (x2c)
      write_square(N,pVpsp,"pVp",2);
    write_square(Naux,Asp,"A",2);
    write_C(Naux,N2,Csp);
  }

  delete [] Asp;
  delete [] Csp;
  delete [] Cysp;
  delete [] Aysp;

  return;
}


int main(int argc, char* argv[]) {
  printf("Running test code for SlaterGPU\n");
  int nomp = 1;
  #pragma omp parallel
  nomp = omp_get_num_threads();

  int save_mog_ads = 1;
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

  int Ne = Na + Nb;

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

  int size_ang = get_npts_from_order(nang);
  size_t gs = nrad * size_ang;
  printf(" nrad: %d nang: %d size_ang: %d\n",nrad,nang,size_ang);
  printf(" grid size: %u\n",gs);

  double * ang_g = new double[3*size_ang];
  double * ang_w = new double[size_ang];
  get_angular_grid(size_ang,ang_g,ang_w);

  int jellium = read_int("JELLIUM");
  bool do_ps_integrals = check_PS();

  int prl = 1;

  if (jellium>0)
  {
    //get A matrix
    compute_ps_integrals_to_disk(natoms,atno,coords,basis,basis_aux,prl,jellium);
    do_ps_integrals = 0;
  }

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
   
    //int prl = 1;
    int c4 = 0;
    double* g = nullptr;
     
    if(c4 > 0)
      double* g = new double[N2*N2]();

    #pragma acc enter data create(A[0:Naux2],C[0:N2a],S[0:N2],En[0:N2],T[0:N2],pVp[0:N2])
    #pragma acc enter data create(V[0:nc],dV[0:3*nc],Pao[0:N2],coordsc[0:3*nc])
    if(c4 > 0)
      #pragma acc enter data create(g[0:N2*N2])
    printf("1e ints: %d\n2c2e ints: %d\n3c3e ints: %d\n4c4e ints: %d\n",N2, Naux2, N2a, N2*N2);

    if (save_mog_ads)
    {
      int gsa = gs*natoms;
      int gsa6 = gsa*6;
      int Nos = No + (Na - Nb);
      double* grid =  new double[gsa6];
      double* wt = new double[gsa];
      for (int i=0;i<gsa6;i++)
      {
        if(i<gsa) wt[i] = 0.;
        grid[i] = 0.;
      }
      #pragma acc enter data copyin(grid[0:gsa6],wt[0:gsa],ang_g[0:3*size_ang],ang_w[0:size_ang])
      get_becke_grid_full(natoms,atno,coords,nrad,nang,ang_g,ang_w,6,grid,wt);
      save_mog(natoms,atno,coords,nrad,gsa,basis,Nos,grid,wt,prl);
      #pragma acc exit data delete(grid[0:gsa6],wt[0:gsa],ang_g[0:3*size_ang],ang_w[0:size_ang])
      delete [] grid;
      delete [] wt;
      return 0;
    }

    if (gbasis)
    {
      printf("\n TESTING gaussian integrals \n");
      compute_gaussian_integrals_to_disk(N,Naux,natoms,nbas,nenv,nbas_ri,atm,bas,env);
    }
    else if (jellium > 0)
    {
      do_jellium(natoms,Ne,basis,basis_aux,nrad,nang,ang_g,ang_w,prl);
    }
    else if (do_ps_integrals)
      compute_ps_integrals_to_disk(natoms,atno,coords,basis,basis_aux,prl,jellium);
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
      //if (c4 > 0)
      //  compute_all_4c_v2(natoms,atno,coordsf,basis,nrad,size_ang,ang_g,ang_w,g,prl);

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
      //if (c4 > 0)
      //  printf("Integral 4c (v2)    time: %5.3e s\n",(double)elapsed67/1.e9);

      printf("-------------------------------\n");
      printf("Integral tot  time: %5.3e s\n",(double)elapsed17/1.e9);
      printf("-------------------------------\n");

      #pragma acc exit data copyout(A[0:Naux2],C[0:N2a],S[0:N2],En[0:N2],T[0:N2],pVp[0:N2])
      #pragma acc exit data copyout(V[0:nc],dV[0:3*nc],Pao[0:N2],coordsc[0:3*nc])

      if(c4 > 0)
        #pragma acc exit data copyout(g[0:N2*N2])

      if (prl > 0) printf("Printing Standard Integral Files:\n");
      write_S_En_T(N,S,En,T);
      write_square(Naux,A,"A",2);
      bool x2c = read_int("X2C");
      if (x2c)
        write_square(N,pVp,"pVp",2);
      write_C(Naux, N2, C);

      if(c4 > 0)
        write_square(N2,g,"g",2);
    }

    delete [] A, Anorm, C, S, En, T, pVp;
    delete [] V, dV, Pao, coordsc;
    if(c4 > 0)  
      delete [] g;
  }

  delete [] ang_g, ang_w;

  printf("\n  Done!\n");
  return 0;

}
