#include "integrals.h"
#include "prosph.h"
#include "quad.h"

#define SWITCH_23 0

//now contains 4c_ol_fast, when memory isn't an issue
//original 4c_ol function now partitioned at significant increased cost

//add collapse to initializations (done)
//Vk edit 1
void auto_crash();


void get_ztm_lm(int s1, int s2, vector<vector<double> >& basis, double& ztm, int& lm)
{
  ztm = 10; lm = 0;

 //minimum value of zeta
 //maximum value of l over basis
  for (int i1=s1;i1<s2;i1++)
  if (basis[i1][3]<ztm)
    ztm = basis[i1][3];
  for (int i1=s1;i1<s2;i1++)
  if (basis[i1][1]>lm)
    lm = basis[i1][1];

 #if 0
  lm = 0;
 //maximum values of zeta/l over basis
  for (int i1=s1;i1<s2;i1++)
  if (basis[i1][3]>ztm)
    ztm = basis[i1][3];
  for (int i1=s1;i1<s2;i1++)
  if (basis[i1][1]>lm)
    lm = basis[i1][1];
 #endif

  return;
}

void reduce_3cenp(double Z3, int s1, int s2, int s3, int s4, int N, int iN, int gs, double* grid, double** valS1, double** valS2, double* valt, double* wt, double* pVp1)
{
  int gs3 = 3*gs;
  int gs6 = 6*gs;
  int N2 = N*N;

 #pragma acc parallel loop present(grid[0:gs6],wt[0:gs],valt[0:gs])
  for (int j=0;j<gs;j++)
  {
    double R = grid[6*j+3];
    valt[j] = -Z3/R*wt[j];
  }

 #pragma acc parallel loop collapse(2) present(valS1[0:iN][0:gs3],valS2[0:iN][0:gs3],valt[0:gs],pVp1[0:N2])
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3;
    double* valm = valS1[ii1]; 
    double* valn = valS2[ii2];

    double val = 0.;
   #pragma acc loop reduction(+:val)
    for (int j=0;j<gs;j++)
    {
      double dp = valm[3*j+0]*valn[3*j+0] + valm[3*j+1]*valn[3*j+1] + valm[3*j+2]*valn[3*j+2];
      val += dp*valt[j];
    }

    pVp1[i1*N+i2] += val;
  }

  return;
}

//flattened reduce_3cenp overloaded to accept double *

void reduce_3cenp(double Z3, int s1, int s2, int s3, int s4, int N, int iN, int gs,
                  double* grid, double* valS1, double* valS2, double* valt, double* wt, double* pVp1)
{
  int gs3 = 3 * gs;
  int gs6 = 6 * gs;
  int N2 = N * N;

  // Compute valt
  #pragma acc parallel loop present(grid[0:gs6], wt[0:gs], valt[0:gs])
  for (int j = 0; j < gs; j++) {
    double R = grid[6 * j + 3];
    valt[j] = -Z3 / R * wt[j];
  }

  // Contract with flattened valS1/valS2
  #pragma acc parallel loop collapse(2) present(valS1[0:iN*gs3], valS2[0:iN*gs3], valt[0:gs], pVp1[0:N2])
  for (int i1 = s1; i1 < s2; i1++)
  for (int i2 = s3; i2 < s4; i2++) {
    int ii1 = i1 - s1;
    int ii2 = i2 - s3;

    double* valm = &valS1[ii1 * gs3];
    double* valn = &valS2[ii2 * gs3];

    double val = 0.0;
    #pragma acc loop reduction(+:val)
    for (int j = 0; j < gs; j++) {
      double dp = valm[3 * j + 0] * valn[3 * j + 0] +
                  valm[3 * j + 1] * valn[3 * j + 1] +
                  valm[3 * j + 2] * valn[3 * j + 2];
      val += dp * valt[j];
    }

    pVp1[i1 * N + i2] += val;
  }
}



void reduce_3cenp_v2(double Z3, int s1, int s2, int s3, int s4, int N, int iN, int gs0, int qos, int gs, double* grid, double** valS1, double** valS2, double* valt, double* wt, double* pVp1)
{
  return reduce_3cenp(Z3,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,pVp1);

 //slower version of integration
 // integrates within cells first, then adds it up
 // could batch over multiple cells
 // tested against first version, no real difference
  int gs3 = 3*gs;
  int gs6 = 6*gs;
  int N2 = N*N;

  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3;
    double* valm = valS1[ii1]; 
    double* valn = valS2[ii2];

   #pragma acc parallel loop present(valt[0:gs],grid[0:gs6],wt[0:gs])
    for (int j=0;j<gs;j++)
    {
      double R = grid[6*j+3];
      valt[j] = -Z3/R*wt[j];
    }

   #pragma acc parallel loop present(valm[0:gs3],valn[0:gs3],valt[0:gs]) 
    for (int j=0;j<gs0;j++)
    {
      double val1 = 0.;
     #pragma acc loop reduction(+:val1)
      for (int k=0;k<qos;k++)
      {
        int jk = j*qos+k;
        double dp = valm[3*jk+0]*valn[3*jk+0] + valm[3*jk+1]*valn[3*jk+1] + valm[3*jk+2]*valn[3*jk+2];
        val1 += dp*valt[jk];
      }
      valt[j] = val1;
    }

    double val = 0.;
   #pragma acc parallel loop present(valt[0:gs]) reduction(+:val)
    for (int j=0;j<gs0;j++)
      val += valt[j];

   #pragma acc serial present(pVp1[0:N2])
    pVp1[i1*N+i2] += val;
  }

  return;
}

//overloaded and flattened
void reduce_3cen(double Z3, int s1, int s2, int s3, int s4,
                 int N, int iN, int gs, double* grid,
                 double* valS1, double* valS2, double* valt,
                 double* wt, double* En1)
{
  int gs6 = 6 * gs;
  int N2 = N * N;

  if (wt == nullptr) {
    #pragma acc parallel loop present(grid[0:gs6], valt[0:gs])
    for (int j = 0; j < gs; j++) {
      double R = grid[6 * j + 3];
      valt[j] = -Z3 / R;
    }
  } else {
    #pragma acc parallel loop present(grid[0:gs6], wt[0:gs], valt[0:gs])
    for (int j = 0; j < gs; j++) {
      double R = grid[6 * j + 3];
      valt[j] = -Z3 / R * wt[j];
    }
  }

  #pragma acc parallel loop collapse(2) present(valS1[0:iN * gs], valS2[0:iN * gs], valt[0:gs], En1[0:N2])
  for (int i1 = s1; i1 < s2; i1++) {
    for (int i2 = s3; i2 < s4; i2++) {
      int ii1 = i1 - s1;
      int ii2 = i2 - s3;

      double val = 0.0;

      #pragma acc loop reduction(+:val)
      for (int j = 0; j < gs; j++) {
        val += valS1[ii1 * gs + j] * valS2[ii2 * gs + j] * valt[j];
      }

      En1[i1 * N + i2] += val;
    }
  }

  return;
}


void reduce_3cen(double Z3, int s1, int s2, int s3, int s4, int N, int iN, int gs, double* grid, double** valS1, double** valS2, double* valt, double* wt, double* En1)
{
  int gs6 = 6*gs;
  int N2 = N*N;

  if (wt==NULL)
 #pragma acc parallel loop present(grid[0:gs6],valt[0:gs])
  for (int j=0;j<gs;j++)
  {
    double R = grid[6*j+3];
    valt[j] = -Z3/R;
  }

  if (wt!=NULL)
 #pragma acc parallel loop present(grid[0:gs6],wt[0:gs],valt[0:gs])
  for (int j=0;j<gs;j++)
  {
    double R = grid[6*j+3];
    valt[j] = -Z3/R*wt[j];

    //printf(" b13: %i %i   R: %8.5f \n",s1,s3,R);
  }

 #pragma acc parallel loop collapse(2) present(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valt[0:gs],En1[0:N2])
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3;
    double* valm = valS1[ii1]; 
    double* valn = valS2[ii2];

    double val = 0.;
   #pragma acc loop reduction(+:val)
    for (int j=0;j<gs;j++)
      val += valm[j]*valn[j]*valt[j];

    En1[i1*N+i2] += val;
  }

  return;
}

void init_s12nw(int s1, int s2, int s3, int s4, int iN, int gs, double** valS1, double** valS2)
{
 #pragma acc parallel loop collapse(2) present(valS1[0:iN][0:gs])
  for (int ii1=0;ii1<s2-s1;ii1++)
  {
    for (int j=0;j<gs;j++)
      valS1[ii1][j] = 1.;
  }

 #pragma acc parallel loop collapse(2) present(valS2[0:iN][0:gs])
  for (int ii1=0;ii1<s4-s3;ii1++)
  {
    for (int j=0;j<gs;j++)
      valS2[ii1][j] = 1.;
  }
}

// Overloaded init_s12nw for flat double* arrays
void init_s12nw(int s1, int s2, int s3, int s4, int iN, int gs, double* valS1, double* valS2)
{
  // Initialize valS1 block (flattened)
  #pragma acc parallel loop collapse(2) present(valS1[0:iN * gs])
  for (int ii1 = 0; ii1 < s2 - s1; ii1++)
  {
    for (int j = 0; j < gs; j++)
    {
      valS1[ii1 * gs + j] = 1.0;
    }
  }

  // Initialize valS2 block (flattened)
  #pragma acc parallel loop collapse(2) present(valS2[0:iN * gs])
  for (int ii2 = 0; ii2 < s4 - s3; ii2++)
  {
    for (int j = 0; j < gs; j++)
    {
      valS2[ii2 * gs + j] = 1.0;
    }
  }
} // END init_s12nw (flat)


void eval_s12(int s1, int s2, int s3, int s4, vector<vector<double> >& basis, int iN, int gs, double* grid, double** valS1, double** valS2)
{
  for (int i1=s1;i1<s2;i1++)
  {
    int ii1 = i1-s1;

    vector<double> basis1 = basis[i1];
    int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];
  
   //S
    eval_shd(ii1,gs,grid,valS1[ii1],n1,l1,m1,zeta1);
  }

  for (int i1=s3;i1<s4;i1++)
  {
    int ii1 = i1-s3;

    vector<double> basis1 = basis[i1];
    int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];
  
   //S
    eval_shd(ii1,gs,grid,valS2[ii1],n1,l1,m1,zeta1);
  }

  return;
}

void get_atnod(int natoms, vector<vector<double> > basis, double* atnod)
{
  for (int n=0;n<basis.size();n++)
  {
    int wa = basis[n][9];
    atnod[wa] = basis[n][8];
  }

  //printf(" atnod:");
  //for (int n=0;n<natoms;n++)
  //  printf(" %4.1f",atnod[n]);
  //printf("\n");
}

void compute_STEn_ps(int natoms, int* atno, double* coords, vector<vector<double> > &basis, int quad_order, int nmu, int nnu, int nphi, double* S, double* T, double* En, int prl)
{
  if (prl>-1) printf("  beginning compute_STEn_ps \n");

  int nomp_max = 1;
 #pragma omp parallel
  nomp_max = omp_get_num_threads();
  
  int ngpu = 0;
 #if USE_ACC
  ngpu = acc_get_num_devices(acc_device_nvidia);
  int nomp = ngpu;
 #else
  int nomp = nomp_max;
 #endif

  int N = basis.size();
  int N2 = N*N;

  double atnod[natoms];
  get_atnod(natoms,basis,atnod);

  int qos = quad_order*quad_order*quad_order;
  int gs = nmu*nnu*nphi*qos;
  int gs6 = 6*gs;

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

  double* grid = new double[gs6];
  double* wt = new double[gs];

  double gpumem = (double)acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  double togb = 1./1024./1024./1024.;

 //this calculation not accurate yet
  int Nmax = 150;
  double mem0 = gs*7.+1.*nmu*nnu*nphi;
  while (Nmax>0)
  {
    //double mem1 = 8.*(gs*3.*Nmax + 1.*mem0);
    double mem1 = 8.*(3.*gs*Nmax + mem0 + 3.*N2 + 4.*gs6 + 1.*gs);
    if (mem1<gpumem)
    {
      if (prl>1) printf("    mem0: %5.3f mem1: %5.3f \n",mem0*togb,mem1*togb);
      break;
    }
    Nmax--;
  }
  //if (Nmax>8)
  //  Nmax -= 8;
  //Nmax = 50;
  //int NNm = read_int("NMAX");
  //Nmax = NNm;
  if (Nmax<=0) { printf("\n ERROR: couldn't calculate gpu memory requirements (mem0: %5.3f mem1: %5.3f) \n",mem0*togb,8.*(3.*gs*Nmax + mem0 + 3.*N2 + 4.*gs6 + 1.*gs)*togb); exit(-1); }

  vector<vector<int> > n2ip;
  int imaxN = get_imax_n2ip(Nmax,natoms,N,basis,n2ip);
  if (prl>1) printf("   imaxN: %2i \n",imaxN);

  //double gsxvalsv = 8.*(imaxN*gs*2. + 1.*mem0); //vals+grid/wt+gridm
  //double gsxvalsv_gb = gsxvalsv/1024./1024./1024.;
  //printf("   estimated memory needed (2c): %6.3f GB \n",gsxvalsv_gb);

  double gpumem_gb = gpumem/1024./1024./1024.;
  printf("   gpu memory available: %6.3f GB \n",gpumem_gb);

  //int* n2i = new int[natoms];
  //int imaxN = get_imax_n2i(natoms,N,basis,n2i);
  //printf(" imaxN: %2i \n",imaxN);

 //intermediate storage
  int iN = imaxN;
  double** valS1 = new double*[iN]; for (int i=0;i<iN;i++) valS1[i] = new double[gs];
  double** valS2 = new double*[iN]; for (int i=0;i<iN;i++) valS2[i] = new double[gs];
  double** valT1 = new double*[iN]; for (int i=0;i<iN;i++) valT1[i] = new double[gs];
  double* valt = new double[gs];

 #if USE_ACC
 //#pragma omp parallel for schedule(static) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = n;
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc enter data create(S[0:N2],T[0:N2],En[0:N2])

    #pragma acc enter data create(grid[0:gs6],wt[0:gs])
    #pragma acc enter data create(valS1[0:iN][0:gs],valS2[0:iN][0:gs])
    #pragma acc enter data create(valT1[0:iN][0:gs])
    #pragma acc enter data create(valt[0:gs])

    acc_assign(N2,S,0.);
    acc_assign(N2,T,0.);
    acc_assign(N2,En,0.);
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  double gpumem_2 = 1.*acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  printf("   after alloc, gpu memory available: %6.3f GB \n",gpumem_2*togb); 

 #pragma omp parallel for schedule(dynamic) num_threads(nomp)
  for (int m=0;m<natoms;m++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);

    //double Z1 = (double)atno[m];
    double Z1 = atnod[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];
    double coordn[6];
    coordn[0] = coordn[1] = coordn[2] = 0.;

    generate_ps_quad_grid(Z1,1,coordn,quad_order,quad_order,nmu,nnu,nphi,grid,wt);
    //generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);
    add_r1_to_grid(gs,grid,0.,0.,0.);

   //working on this block of the matrix
    //int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    //j>i
    for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
    for (int sp2=sp1;sp2<n2ip[m].size()-1;sp2++)
    {
      int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
      int s3 = n2ip[m][sp2]; int s4 = n2ip[m][sp2+1];

     #pragma acc parallel loop present(valS1[0:iN][0:gs])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS1[ii1][j] = 1.;
      }

     #pragma acc parallel loop present(valS2[0:iN][0:gs],wt[0:gs])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS2[ii1][j] = wt[j];
      }

     //first compute single atom ints
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

       //S
        eval_shd(ii1,gs,grid,valS1[ii1],n1,l1,m1,zeta1);
      }

      for (int i1=s3;i1<s4;i1++)
      {
        int ii1 = i1-s3;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

       //S
        eval_shd(ii1,gs,grid,valS2[ii1],n1,l1,m1,zeta1);
      }

     #pragma acc parallel loop present(valS1[0:iN][0:gs],valT1[0:iN][0:gs])
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

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_ke(gs,grid,valT1[ii1],n1,l1,zeta1);
      }

      #pragma acc wait

      reduce_2c1(s1,s2,s3,s4,gs,valS1,valS2,iN,N,S);
      reduce_2c1(s1,s2,s3,s4,gs,valT1,valS2,iN,N,T);

     /////////////////////////////////////////////////////////////////////
     //electron-nuclear attraction
      reduce_3cen(Z1,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,NULL,En);
     /////////////////////////////////////////////////////////////////////
    } //loop sp over s12

   //two-atom ints
    for (int n=m+1;n<natoms;n++)
    {
      //double Z2 = (double)atno[n];
      double Z2 = atnod[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;
      coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;

      generate_ps_quad_grid(0.,2,coordn,quad_order,quad_order,nmu,nnu,nphi,grid,wt);
      add_r1_to_grid(gs,grid,0.,0.,0.);

      //int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

     //basis ftns on each center
      for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
      for (int sp2=0;sp2<n2ip[n].size()-1;sp2++)
      {
        int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
        int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2+1];

       #pragma acc parallel loop present(valS1[0:iN][0:gs])
        for (int ii1=0;ii1<s2-s1;ii1++)
        {
         #pragma acc loop
          for (int j=0;j<gs;j++)
            valS1[ii1][j] = 1.;
        }

       #pragma acc parallel loop present(valS2[0:iN][0:gs],wt[0:gs])
        for (int ii1=0;ii1<s4-s3;ii1++)
        {
         #pragma acc loop
          for (int j=0;j<gs;j++)
            valS2[ii1][j] = wt[j];
        }

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<double> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

         //S
          eval_shd(ii1,gs,grid,valS1[ii1],n1,l1,m1,zeta1);
        }

       #pragma acc parallel loop present(valS1[0:iN][0:gs],valT1[0:iN][0:gs])
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

          vector<double> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

          eval_ke(gs,grid,valT1[ii1],n1,l1,zeta1);
        }

       //second center
        recenter_grid_zero(gs,grid,-A12,-B12,-C12);

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

         //S
          eval_shd(ii2,gs,grid,valS2[ii2],n2,l2,m2,zeta2);
        }

        reduce_2c1(s1,s2,s3,s4,gs,valS1,valS2,iN,N,S);
        reduce_2c1(s1,s2,s3,s4,gs,valT1,valS2,iN,N,T);

       ///////////////////////////////////////////////////////////////////
       //electron-nuclear attraction centers 1+2
       #pragma acc parallel loop present(valS2[0:iN][0:gs],wt[0:gs])
        for (int ii1=0;ii1<s4-s3;ii1++)
        {
         #pragma acc loop
          for (int j=0;j<gs;j++)
            valS2[ii1][j] /= wt[j];
        }

        reduce_3cen(Z2,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,En);
        recenter_grid_zero(gs,grid,A12,B12,C12);
        reduce_3cen(Z1,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,En);
        ////////////////////////////////////////////////////////////////////

      } //loop sp over s14

     //En for basis ftns on same center (1)
      for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
      for (int sp2=sp1;sp2<n2ip[m].size()-1;sp2++)
      {
        int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
        int s3 = n2ip[m][sp2]; int s4 = n2ip[m][sp2+1];

        //printf("  En 2c(1). s12: %2i %2i s34: %2i %2i \n",s1,s2,s3,s4);

        init_s12nw(s1,s2,s3,s4,iN,gs,valS1,valS2);
        eval_s12(s1,s2,s3,s4,basis,iN,gs,grid,valS1,valS2);

       //second center
        recenter_grid_zero(gs,grid,-A12,-B12,-C12);

       ///////////////////////////////////////////////////////////////////
       //electron-nuclear attraction center 2
        reduce_3cen(Z2,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,En);
        ////////////////////////////////////////////////////////////////////

        recenter_grid_zero(gs,grid,A12,B12,C12);

      } //loop sp over s14

     //En for basis ftns on same center (2)
      for (int sp1=0;sp1<n2ip[n].size()-1;sp1++)
      for (int sp2=sp1;sp2<n2ip[n].size()-1;sp2++)
      {
        int s1 = n2ip[n][sp1]; int s2 = n2ip[n][sp1+1];
        int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2+1];

        //printf("  En 2c(2). s12: %2i %2i s34: %2i %2i \n",s1,s2,s3,s4);

        init_s12nw(s1,s2,s3,s4,iN,gs,valS1,valS2);

       //second center
        recenter_grid_zero(gs,grid,-A12,-B12,-C12);

        eval_s12(s1,s2,s3,s4,basis,iN,gs,grid,valS1,valS2);

        recenter_grid_zero(gs,grid,A12,B12,C12);

       ///////////////////////////////////////////////////////////////////
       //electron-nuclear attraction center 1
        reduce_3cen(Z1,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,En);
        ////////////////////////////////////////////////////////////////////

      } //loop sp over s14

    } //loop n over second atom

  } //loop m over natoms
  acc_set_device_num(0,acc_device_nvidia);

  if (nomp>1)
  {
   //gather parallel fragments
    double St[N2]; double Tt[N2]; double Ent[N2];
    for (int j=0;j<N2;j++)
      St[j] = Tt[j] = Ent[j] = 0.;

    for (int n=0;n<nomp;n++)
    {
      int tid = n;
      acc_set_device_num(tid,acc_device_nvidia);

      #pragma acc update self(S[0:N2],T[0:N2],En[0:N2])

      for (int j=0;j<N2;j++)
      {
        St[j] += S[j];
        Tt[j] += T[j];
        Ent[j] += En[j];
      }
    }

    for (int j=0;j<N2;j++)
    {
      S[j] = St[j];
      En[j] = Ent[j];
      T[j] = Tt[j];
    }

    acc_set_device_num(0,acc_device_nvidia);
    #pragma acc update device(S[0:N2],T[0:N2],En[0:N2])
  }

  double norm[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];
  #pragma acc enter data copyin(norm[0:N])

  #pragma acc parallel loop collapse(2) independent present(S[0:N2],T[0:N2],En[0:N2],norm[0:N])
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
  {
    double n12 = norm[i]*norm[j];
    S[i*N+j] *= n12;
    T[i*N+j] *= -0.5*n12;
    En[i*N+j] *= n12;
  }
  #pragma acc exit data delete(norm[0:N])

 //symmetrize wrt ij
 #pragma acc parallel loop independent present(S[0:N2],T[0:N2],En[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i+1;j<N;j++)
  {
    S[j*N+i] = S[i*N+j];
    T[j*N+i] = T[i*N+j];
    En[j*N+i] = En[i*N+j];
  }

  #pragma acc update self(S[0:N2],T[0:N2],En[0:N2])

  if (prl>0)
  {
    printf("  S diagonals accuracy: \n");
    int nlow = 0;
    for (int j=0;j<N;j++)
    {
      double v1 = log10(fabs(1.-S[j*N+j])+1.e-16);
      printf(" %4.1f",v1);

      if (v1>-9.)
        nlow++;
    }
    printf("\n");
    if (nlow)
      printf("   WARNING: found %2i low-accuracy diagonal elements \n",nlow);
  }

 //might as well eliminate errors on diagonal
  for (int i=0;i<N;i++)
    S[i*N+i] = 1.;

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

 #if USE_ACC
 //#pragma omp parallel for schedule(static) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = n;
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc exit data delete(S[0:N2],T[0:N2],En[0:N2])
    #pragma acc exit data delete(grid[0:gs6],wt[0:gs])
    #pragma acc exit data delete(valS1[0:iN][0:gs],valS2[0:iN][0:gs])
    #pragma acc exit data delete(valT1[0:iN][0:gs])
    #pragma acc exit data delete(valt[0:gs])
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  //printf(" done with dealloc in STEn integrals \n"); fflush(stdout);
  //auto_crash();

  //delete [] n2i;

  for (int i=0;i<iN;i++) delete [] valS1[i];
  for (int i=0;i<iN;i++) delete [] valS2[i];
  for (int i=0;i<iN;i++) delete [] valT1[i];
  delete [] valS1; delete [] valS2;
  delete [] valT1;
  delete [] valt;

  delete [] grid;
  delete [] wt;

  return;
}

/////TESTING NEEDED

void compute_pVp_ps(int natoms, int* atno, double* coords, vector<vector<double> > &basis, int quad_order, int nmu, int nnu, int nphi, double* pVp, int prl)
{
 //barely tested
  if (prl>-1) printf("  beginning compute_pVp_ps \n");

  int nomp_max = 1;
 #pragma omp parallel
  nomp_max = omp_get_num_threads();
  
  int ngpu = 0;
 #if USE_ACC
  ngpu = acc_get_num_devices(acc_device_nvidia);
  int nomp = ngpu;
 #else
  int nomp = nomp_max;
 #endif

  int N = basis.size();
  int N2 = N*N;

  double atnod[natoms];
  get_atnod(natoms,basis,atnod);

  int qos = quad_order*quad_order*quad_order;
  int gs0 = nmu*nnu*nphi;
  int gs = gs0*qos;
  int gs3 = 3*gs;
  int gs6 = 6*gs;

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

  double* grid = new double[gs6];
  double* wt = new double[gs];

  double gpumem = (double)acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  double togb = 1./1024./1024./1024.;

 //this calculation not accurate yet
  int Nmax = 150;
  double mem0 = gs*7.+1.*nmu*nnu*nphi;
  while (Nmax>0)
  {
    //double mem1 = 8.*(gs*6.*Nmax + 1.*mem0);
    double mem1 = 8.*(mem0 + 1.*N2 + 3.*gs6 + 2.*gs + 2.*Nmax*gs3);
    if (mem1<gpumem)
    {
      if (prl>1) printf("    mem0: %5.3f mem1: %5.3f \n",mem0*togb,mem1*togb);
      break;
    }
    Nmax--;
  }


  if (Nmax<=0) { printf("\n ERROR: couldn't calculate gpu memory requirements \n"); exit(-1); }

/*
The following tries flattening 2D arrays.
*/

//Double * starts here VK final check


  vector<vector<int> > n2ip;
  int imaxN = get_imax_n2ip(Nmax, natoms, N, basis, n2ip);
  if (prl > 1) printf("   imaxN: %2i \n", imaxN);

  double gpumem_gb = gpumem / 1024. / 1024. / 1024.;
  printf("   gpu memory available: %6.3f GB \n", gpumem_gb);

  // intermediate storage
  int iN = imaxN;
  double* valS1 = new double[iN * gs3];
  double* valS2 = new double[iN * gs3];
  double* valt  = new double[gs];

  #if USE_ACC
  for (int n = 0; n < nomp; n++)
  {
    int tid = n;
    acc_set_device_num(tid, acc_device_nvidia);

    #pragma acc enter data create(pVp[0:N2])
    #pragma acc enter data create(grid[0:gs6], wt[0:gs])
    #pragma acc enter data create(valS1[0:iN * gs3], valS2[0:iN * gs3])
    #pragma acc enter data create(valt[0:gs])

    acc_assign(N2, pVp, 0.);
  }
  acc_set_device_num(0, acc_device_nvidia);
  #endif

  double gpumem_2 = 1. * acc_get_property(0, acc_device_nvidia, acc_property_free_memory);
  printf("   after alloc, gpu memory available: %6.3f GB \n", gpumem_2 * togb);
  printf("\n Double * Problematic zone starts\n");

  #pragma omp parallel for schedule(dynamic) num_threads(nomp)
  for (int m = 0; m < natoms; m++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid, acc_device_nvidia);

    double Z1 = atnod[m];
    double A1 = coords[3 * m + 0]; double B1 = coords[3 * m + 1]; double C1 = coords[3 * m + 2];
    double coordn[6] = {0.};

    generate_ps_quad_grid(Z1, 1, coordn, quad_order, quad_order, nmu, nnu, nphi, grid, wt);
    add_r1_to_grid(gs, grid, 0., 0., 0.);

    for (int sp1 = 0; sp1 < n2ip[m].size() - 1; sp1++)
    for (int sp2 = sp1; sp2 < n2ip[m].size() - 1; sp2++)
    {
      int s1 = n2ip[m][sp1];     int s2 = n2ip[m][sp1 + 1];
      int s3 = n2ip[m][sp2];     int s4 = n2ip[m][sp2 + 1];

      #pragma acc parallel loop present(valS1[0:iN * gs3])
      for (int ii1 = 0; ii1 < s2 - s1; ii1++)
      {
        #pragma acc loop
        for (int j = 0; j < gs3; j++)
          valS1[ii1 * gs3 + j] = 1.;
      }

      #pragma acc parallel loop present(valS2[0:iN * gs3])
      for (int ii1 = 0; ii1 < s4 - s3; ii1++)
      {
        #pragma acc loop
        for (int j = 0; j < gs3; j++)
          valS2[ii1 * gs3 + j] = 1.;
      }

      // evaluate valS1
      for (int i1 = s1; i1 < s2; i1++)
      {
        int ii1 = i1 - s1;
        vector<double> basis1 = basis[i1];
        int n1 = basis1[0], l1 = basis1[1], m1 = basis1[2];
        double zeta1 = basis1[3];

        eval_pd(gs, grid, &valS1[ii1 * gs3], n1, l1, m1, zeta1);
      }

      // evaluate valS2
      for (int i2 = s3; i2 < s4; i2++)
      {
        int ii2 = i2 - s3;
        vector<double> basis1 = basis[i2];
        int n1 = basis1[0], l1 = basis1[1], m1 = basis1[2];
        double zeta1 = basis1[3];

        eval_pd(gs, grid, &valS2[ii2 * gs3], n1, l1, m1, zeta1);
      }

      #pragma acc wait

      reduce_3cenp(Z1, s1, s2, s3, s4, N, iN, gs, grid, valS1, valS2, valt, wt, pVp);
    }

//NEWER needs fixing to support double *

//two-atom ints
    for (int n=m+1;n<natoms;n++)
    {
      double Z2 = atnod[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;
      coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;

     //pull this out of loop later on
      generate_ps_quad_grid(0.,2,coordn,quad_order,quad_order,nmu,nnu,nphi,grid,wt);
      add_r1_to_grid(gs,grid,0.,0.,0.);

     //basis functions on each center
      for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
      for (int sp2=0;sp2<n2ip[n].size()-1;sp2++)
      {
        int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
        int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2+1];

        init_s12nw(s1,s2,s3,s4,iN,gs3,valS1,valS2);

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<double> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];
  
         //p
          //eval_pd(gs,grid,valS1[ii1],n1,l1,m1,zeta1);
	        eval_pd(gs, grid, &valS1[ii1 * gs3], n1, l1, m1, zeta1);
        }

       //second center
        recenter_grid_zero(gs,grid,-A12,-B12,-C12);

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

         //p
          //eval_pd(gs,grid,valS2[ii2],n2,l2,m2,zeta2);
          eval_pd(gs, grid, &valS2[ii2 * gs3], n2, l2, m2, zeta2);
        }

       ///////////////////////////////////////////////////////////////////
       //electron-nuclear attraction centers 1+2

        reduce_3cenp(Z2,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,pVp);
        recenter_grid_zero(gs,grid,A12,B12,C12);
        reduce_3cenp(Z1,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,pVp);
        ////////////////////////////////////////////////////////////////////

      } //loop sp over s14

     //basis functions on center 1
      for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
      for (int sp2=sp1;sp2<n2ip[m].size()-1;sp2++)
      {
        int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
        int s3 = n2ip[m][sp2]; int s4 = n2ip[m][sp2+1];

        init_s12nw(s1,s2,s3,s4,iN,gs3,valS1,valS2);

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<double> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];
  
         //p
          //eval_pd(gs,grid,valS1[ii1],n1,l1,m1,zeta1);
          eval_pd(gs, grid, &valS1[ii1 * gs3], n1, l1, m1, zeta1);
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

         //p
          //eval_pd(gs,grid,valS2[ii2],n2,l2,m2,zeta2);
          eval_pd(gs, grid, &valS2[ii2 * gs3], n2, l2, m2, zeta2);
        }

       //second center
        recenter_grid_zero(gs,grid,-A12,-B12,-C12);

       ///////////////////////////////////////////////////////////////////
       //electron-nuclear attraction center 2
        reduce_3cenp(Z2,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,pVp);
        ////////////////////////////////////////////////////////////////////

        recenter_grid_zero(gs,grid,A12,B12,C12);

      } //loop sp over s14

     //basis functions on center 2
      for (int sp1=0;sp1<n2ip[n].size()-1;sp1++)
      for (int sp2=0;sp2<n2ip[n].size()-1;sp2++)
      {
        int s1 = n2ip[n][sp1]; int s2 = n2ip[n][sp1+1];
        int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2+1];

        init_s12nw(s1,s2,s3,s4,iN,gs3,valS1,valS2);

       //second center
        recenter_grid_zero(gs,grid,-A12,-B12,-C12);

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<double> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];
  
         //p
          //eval_pd(gs,grid,valS1[ii1],n1,l1,m1,zeta1);
          eval_pd(gs, grid, &valS1[ii1 * gs3], n1, l1, m1, zeta1);
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

         //p
          //eval_pd(gs,grid,valS2[ii2],n2,l2,m2,zeta2);
          eval_pd(gs, grid, &valS2[ii2 * gs3], n2, l2, m2, zeta2);
        }

        recenter_grid_zero(gs,grid,A12,B12,C12);

       ///////////////////////////////////////////////////////////////////
       //electron-nuclear attraction center 1
        reduce_3cenp(Z1,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,pVp);
        ////////////////////////////////////////////////////////////////////

      } //loop sp over s14

    } //loop n over second atom

  } //loop m over natoms
  acc_set_device_num(0,acc_device_nvidia);

  if (nomp>1)
  {
   //gather parallel fragments
    double pVpt[N2];
    for (int j=0;j<N2;j++)
      pVpt[j] = 0.;

    for (int n=0;n<nomp;n++)
    {
      int tid = n;
      acc_set_device_num(tid,acc_device_nvidia);

      #pragma acc update self(pVp[0:N2])

      for (int j=0;j<N2;j++)
        pVpt[j] += pVp[j];
    }

    for (int j=0;j<N2;j++)
      pVp[j] = pVpt[j];

    acc_set_device_num(0,acc_device_nvidia);
    #pragma acc update device(pVp[0:N2])
  }

  double norm[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];
  #pragma acc enter data copyin(norm[0:N])

  #pragma acc parallel loop collapse(2) independent present(pVp[0:N2],norm[0:N])
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
  {
    double n12 = norm[i]*norm[j];
    pVp[i*N+j] *= n12;
  }
  #pragma acc exit data delete(norm[0:N])

 //symmetrize wrt ij
 #pragma acc parallel loop independent present(pVp[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i+1;j<N;j++)
  {
    pVp[j*N+i] = pVp[i*N+j];
  }

  #pragma acc update self(pVp[0:N2])

  if (prl>1)
  {
    printf("\n pVp (2c): \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %12.6f",pVp[i*N+j]);
      printf("\n");
    }
  }

 #if USE_ACC
 //#pragma omp parallel for schedule(static) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = n;
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc exit data delete(pVp[0:N2])
    #pragma acc exit data delete(grid[0:gs6], wt[0:gs])
    #pragma acc exit data delete(valS1[0:iN * gs3], valS2[0:iN * gs3])
    #pragma acc exit data delete(valt[0:gs])
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif


  delete [] valS1; delete [] valS2;
  delete [] valt;

  delete [] grid;
  delete [] wt;

  printf("\nYippiii!\n");

  return;

}

////END

/*
void compute_pVp_ps(int natoms, int* atno, double* coords, vector<vector<double> > &basis, int quad_order, int nmu, int nnu, int nphi, double* pVp, int prl)
{
 //barely tested
  if (prl>-1) printf("  beginning compute_pVp_ps \n");

  int nomp_max = 1;
 #pragma omp parallel
  nomp_max = omp_get_num_threads();
  
  int ngpu = 0;
 #if USE_ACC
  ngpu = acc_get_num_devices(acc_device_nvidia);
  int nomp = ngpu;
 #else
  int nomp = nomp_max;
 #endif

  int N = basis.size();
  int N2 = N*N;

  double atnod[natoms];
  get_atnod(natoms,basis,atnod);

  int qos = quad_order*quad_order*quad_order;
  int gs0 = nmu*nnu*nphi;
  int gs = gs0*qos;
  int gs3 = 3*gs;
  int gs6 = 6*gs;

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

  double* grid = new double[gs6];
  double* wt = new double[gs];

  double gpumem = (double)acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  double togb = 1./1024./1024./1024.;

 //this calculation not accurate yet
  int Nmax = 150;
  double mem0 = gs*7.+1.*nmu*nnu*nphi;
  while (Nmax>0)
  {
    //double mem1 = 8.*(gs*6.*Nmax + 1.*mem0);
    double mem1 = 8.*(mem0 + 1.*N2 + 3.*gs6 + 2.*gs + 2.*Nmax*gs3);
    if (mem1<gpumem)
    {
      if (prl>1) printf("    mem0: %5.3f mem1: %5.3f \n",mem0*togb,mem1*togb);
      break;
    }
    Nmax--;
  }


  if (Nmax<=0) { printf("\n ERROR: couldn't calculate gpu memory requirements \n"); exit(-1); }



//Double * starts here VK final check



#if 1

  vector<vector<int> > n2ip;
  int imaxN = get_imax_n2ip(Nmax, natoms, N, basis, n2ip);
  if (prl > 1) printf("   imaxN: %2i \n", imaxN);

  double gpumem_gb = gpumem / 1024. / 1024. / 1024.;
  printf("   gpu memory available: %6.3f GB \n", gpumem_gb);

  // intermediate storage
  int iN = imaxN;
  double* valS1 = new double[iN * gs3];
  double* valS2 = new double[iN * gs3];
  double* valt  = new double[gs];

  #if USE_ACC
  for (int n = 0; n < nomp; n++)
  {
    int tid = n;
    acc_set_device_num(tid, acc_device_nvidia);

    #pragma acc enter data create(pVp[0:N2])
    #pragma acc enter data create(grid[0:gs6], wt[0:gs])
    #pragma acc enter data create(valS1[0:iN * gs3], valS2[0:iN * gs3])
    #pragma acc enter data create(valt[0:gs])

    acc_assign(N2, pVp, 0.);
  }
  acc_set_device_num(0, acc_device_nvidia);
  #endif

  double gpumem_2 = 1. * acc_get_property(0, acc_device_nvidia, acc_property_free_memory);
  printf("   after alloc, gpu memory available: %6.3f GB \n", gpumem_2 * togb);
  printf("\n Double * Problematic zone starts\n");

  #pragma omp parallel for schedule(dynamic) num_threads(nomp)
  for (int m = 0; m < natoms; m++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid, acc_device_nvidia);

    double Z1 = atnod[m];
    double A1 = coords[3 * m + 0]; double B1 = coords[3 * m + 1]; double C1 = coords[3 * m + 2];
    double coordn[6] = {0.};

    generate_ps_quad_grid(Z1, 1, coordn, quad_order, quad_order, nmu, nnu, nphi, grid, wt);
    add_r1_to_grid(gs, grid, 0., 0., 0.);

    for (int sp1 = 0; sp1 < n2ip[m].size() - 1; sp1++)
    for (int sp2 = sp1; sp2 < n2ip[m].size() - 1; sp2++)
    {
      int s1 = n2ip[m][sp1];     int s2 = n2ip[m][sp1 + 1];
      int s3 = n2ip[m][sp2];     int s4 = n2ip[m][sp2 + 1];

      #pragma acc parallel loop present(valS1[0:iN * gs3])
      for (int ii1 = 0; ii1 < s2 - s1; ii1++)
      {
        #pragma acc loop
        for (int j = 0; j < gs3; j++)
          valS1[ii1 * gs3 + j] = 1.;
      }

      #pragma acc parallel loop present(valS2[0:iN * gs3])
      for (int ii1 = 0; ii1 < s4 - s3; ii1++)
      {
        #pragma acc loop
        for (int j = 0; j < gs3; j++)
          valS2[ii1 * gs3 + j] = 1.;
      }

      // evaluate valS1
      for (int i1 = s1; i1 < s2; i1++)
      {
        int ii1 = i1 - s1;
        vector<double> basis1 = basis[i1];
        int n1 = basis1[0], l1 = basis1[1], m1 = basis1[2];
        double zeta1 = basis1[3];

        eval_pd(gs, grid, &valS1[ii1 * gs3], n1, l1, m1, zeta1);
      }

      // evaluate valS2
      for (int i2 = s3; i2 < s4; i2++)
      {
        int ii2 = i2 - s3;
        vector<double> basis1 = basis[i2];
        int n1 = basis1[0], l1 = basis1[1], m1 = basis1[2];
        double zeta1 = basis1[3];

        eval_pd(gs, grid, &valS2[ii2 * gs3], n1, l1, m1, zeta1);
      }

      #pragma acc wait

      reduce_3cenp(Z1, s1, s2, s3, s4, N, iN, gs, grid, valS1, valS2, valt, wt, pVp);
    }
  }

  printf("\n Double * Problematic zone ends\n");

  // GPU cleanup
  #if USE_ACC
  #pragma acc exit data delete(pVp[0:N2])
  #pragma acc exit data delete(grid[0:gs6], wt[0:gs])
  #pragma acc exit data delete(valS1[0:iN * gs3], valS2[0:iN * gs3])
  #pragma acc exit data delete(valt[0:gs])
  #endif

  // CPU cleanup
  delete[] valS1;
  delete[] valS2;
  delete[] valt;

  exit(0);

#endif

}

*/

//Double * ends here VK final check

//Double ** starts here VK final check

/*

#if 1

  vector<vector<int> > n2ip;
  int imaxN = get_imax_n2ip(Nmax, natoms, N, basis, n2ip);
  if (prl > 1) printf("   imaxN: %2i \n", imaxN);

  double gpumem_gb = gpumem / 1024. / 1024. / 1024.;
  printf("   gpu memory available: %6.3f GB \n", gpumem_gb);

  // intermediate storage
  int iN = imaxN;
  double** valS1 = new double*[iN]; for (int i = 0; i < iN; i++) valS1[i] = new double[gs3];
  double** valS2 = new double*[iN]; for (int i = 0; i < iN; i++) valS2[i] = new double[gs3];
  double* valt  = new double[gs];

  #if USE_ACC
  for (int n = 0; n < nomp; n++)
  {
    int tid = n;
    acc_set_device_num(tid, acc_device_nvidia);

    #pragma acc enter data create(pVp[0:N2])
    #pragma acc enter data create(grid[0:gs6], wt[0:gs])
    #pragma acc enter data create(valS1[0:iN][0:gs3], valS2[0:iN][0:gs3])
    #pragma acc enter data create(valt[0:gs])

    acc_assign(N2, pVp, 0.);
  }
  acc_set_device_num(0, acc_device_nvidia);
  #endif

  double gpumem_2 = 1. * acc_get_property(0, acc_device_nvidia, acc_property_free_memory);
  printf("   after alloc, gpu memory available: %6.3f GB \n", gpumem_2 * togb);
  printf("\n Double ** Problematic zone starts\n");

  #pragma omp parallel for schedule(dynamic) num_threads(nomp)
  for (int m = 0; m < natoms; m++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid, acc_device_nvidia);

    double Z1 = atnod[m];
    double A1 = coords[3 * m + 0]; double B1 = coords[3 * m + 1]; double C1 = coords[3 * m + 2];
    double coordn[6] = {0.};

    generate_ps_quad_grid(Z1, 1, coordn, quad_order, quad_order, nmu, nnu, nphi, grid, wt);
    add_r1_to_grid(gs, grid, 0., 0., 0.);

    for (int sp1 = 0; sp1 < n2ip[m].size() - 1; sp1++)
    for (int sp2 = sp1; sp2 < n2ip[m].size() - 1; sp2++)
    {
      int s1 = n2ip[m][sp1];     int s2 = n2ip[m][sp1 + 1];
      int s3 = n2ip[m][sp2];     int s4 = n2ip[m][sp2 + 1];

      #pragma acc parallel loop present(valS1[0:iN][0:gs3])
      for (int ii1 = 0; ii1 < s2 - s1; ii1++)
      {
        #pragma acc loop
        for (int j = 0; j < gs3; j++)
          valS1[ii1][j] = 1.;
      }

      #pragma acc parallel loop present(valS2[0:iN][0:gs3])
      for (int ii1 = 0; ii1 < s4 - s3; ii1++)
      {
        #pragma acc loop
        for (int j = 0; j < gs3; j++)
          valS2[ii1][j] = 1.;
      }

      // evaluate valS1
      for (int i1 = s1; i1 < s2; i1++)
      {
        int ii1 = i1 - s1;
        vector<double> basis1 = basis[i1];
        int n1 = basis1[0], l1 = basis1[1], m1 = basis1[2];
        double zeta1 = basis1[3];

        eval_pd(gs, grid, valS1[ii1], n1, l1, m1, zeta1);
      }

      // evaluate valS2
      for (int i2 = s3; i2 < s4; i2++)
      {
        int ii2 = i2 - s3;
        vector<double> basis1 = basis[i2];
        int n1 = basis1[0], l1 = basis1[1], m1 = basis1[2];
        double zeta1 = basis1[3];

        eval_pd(gs, grid, valS2[ii2], n1, l1, m1, zeta1);
      }

      #pragma acc wait

      // reduce_3cenp(Z1, s1, s2, s3, s4, N, iN, gs, grid, valS1, valS2, valt, wt, pVp);
    }
  }

  printf("\n Double ** Problematic zone ends\n");

  // GPU cleanup
  #if USE_ACC
  #pragma acc exit data delete(pVp[0:N2])
  #pragma acc exit data delete(grid[0:gs6], wt[0:gs])
  #pragma acc exit data delete(valS1[0:iN][0:gs3], valS2[0:iN][0:gs3])
  #pragma acc exit data delete(valt[0:gs])
  #endif

  // CPU cleanup
  for (int i = 0; i < iN; i++) {
    delete[] valS1[i];
    delete[] valS2[i];
  }
  delete[] valS1;
  delete[] valS2;
  delete[] valt;

  exit(0);  // clean exit

#endif

}

*/

//Double ** ends here VK final check




//OLD CHECKS


//Double * starts VK edit

/*

#if 1

  vector<vector<int> > n2ip;
  int imaxN = get_imax_n2ip(Nmax, natoms, N, basis, n2ip);
  if (prl > 1) printf("   imaxN: %2i \n", imaxN);

  double gpumem_gb = gpumem / 1024. / 1024. / 1024.;
  printf("   gpu memory available: %6.3f GB \n", gpumem_gb);

  // intermediate storage
  int iN = imaxN;
  double* valS1 = new double[iN * gs3];
  double* valS2 = new double[iN * gs3];
  double* valt  = new double[gs];

  #if USE_ACC
  for (int n = 0; n < nomp; n++)
  {
    int tid = n;
    acc_set_device_num(tid, acc_device_nvidia);

    #pragma acc enter data create(pVp[0:N2])
    #pragma acc enter data create(grid[0:gs6], wt[0:gs])
    #pragma acc enter data create(valS1[0:iN * gs3], valS2[0:iN * gs3])
    #pragma acc enter data create(valt[0:gs])

    acc_assign(N2, pVp, 0.);
  }
  acc_set_device_num(0, acc_device_nvidia);
  #endif

  double gpumem_2 = 1. * acc_get_property(0, acc_device_nvidia, acc_property_free_memory);
  printf("   after alloc, gpu memory available: %6.3f GB \n", gpumem_2 * togb);

  #pragma omp parallel for schedule(dynamic) num_threads(nomp)
  for (int m = 0; m < natoms; m++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid, acc_device_nvidia);

    double Z1 = atnod[m];
    double A1 = coords[3 * m + 0]; double B1 = coords[3 * m + 1]; double C1 = coords[3 * m + 2];
    double coordn[6] = {0.};

    generate_ps_quad_grid(Z1, 1, coordn, quad_order, quad_order, nmu, nnu, nphi, grid, wt);
    add_r1_to_grid(gs, grid, 0., 0., 0.);

    printf("\n Double * Problematic zone starts");

    for (int sp1 = 0; sp1 < n2ip[m].size() - 1; sp1++)
    for (int sp2 = sp1; sp2 < n2ip[m].size() - 1; sp2++)
    {
      int s1 = n2ip[m][sp1];     int s2 = n2ip[m][sp1 + 1];
      int s3 = n2ip[m][sp2];     int s4 = n2ip[m][sp2 + 1];

      #pragma acc parallel loop present(valS1[0:iN * gs3])
      for (int ii1 = 0; ii1 < s2 - s1; ii1++)
      {
        #pragma acc loop
        for (int j = 0; j < gs3; j++)
          valS1[ii1 * gs3 + j] = 1.;
      }

      #pragma acc parallel loop present(valS2[0:iN * gs3])
      for (int ii1 = 0; ii1 < s4 - s3; ii1++)
      {
        #pragma acc loop
        for (int j = 0; j < gs3; j++)
          valS2[ii1 * gs3 + j] = 1.;
      }

      // evaluate valS1
      for (int i1 = s1; i1 < s2; i1++)
      {
        int ii1 = i1 - s1;
        vector<double> basis1 = basis[i1];
        int n1 = basis1[0], l1 = basis1[1], m1 = basis1[2];
        double zeta1 = basis1[3];

        eval_pd(gs, grid, &valS1[ii1 * gs3], n1, l1, m1, zeta1);
      }

      // evaluate valS2
      for (int i2 = s3; i2 < s4; i2++)
      {
        int ii2 = i2 - s3;
        vector<double> basis1 = basis[i2];
        int n1 = basis1[0], l1 = basis1[1], m1 = basis1[2];
        double zeta1 = basis1[3];

        eval_pd(gs, grid, &valS2[ii2 * gs3], n1, l1, m1, zeta1);
      }

      #pragma acc wait

      //electron-nuclear attraction 
      // reduce_3cenp(Z1, s1, s2, s3, s4, N, iN, gs, grid, valS1, valS2, valt, wt, pVp);
    }

    printf("\n Double * Problematic zone ends");
    exit(-1);
  }
}
#endif

*/

//Double * ends VK edit

//Double ** starts VK edit

/*

#if 1

  vector<vector<int> > n2ip;
  int imaxN = get_imax_n2ip(Nmax, natoms, N, basis, n2ip);
  if (prl > 1) printf("   imaxN: %2i \n", imaxN);

  double gpumem_gb = gpumem / 1024. / 1024. / 1024.;
  printf("   gpu memory available: %6.3f GB \n", gpumem_gb);

  // intermediate storage
  int iN = imaxN;
  double** valS1 = new double*[iN]; for (int i = 0; i < iN; i++) valS1[i] = new double[gs3];
  double** valS2 = new double*[iN]; for (int i = 0; i < iN; i++) valS2[i] = new double[gs3];
  double* valt  = new double[gs];

  #if USE_ACC
  for (int n = 0; n < nomp; n++)
  {
    int tid = n;
    acc_set_device_num(tid, acc_device_nvidia);

    #pragma acc enter data create(pVp[0:N2])
    #pragma acc enter data create(grid[0:gs6], wt[0:gs])
    #pragma acc enter data create(valS1[0:iN][0:gs3], valS2[0:iN][0:gs3])
    #pragma acc enter data create(valt[0:gs])

    acc_assign(N2, pVp, 0.);
  }
  acc_set_device_num(0, acc_device_nvidia);
  #endif

  double gpumem_2 = 1. * acc_get_property(0, acc_device_nvidia, acc_property_free_memory);
  printf("   after alloc, gpu memory available: %6.3f GB \n", gpumem_2 * togb);

  #pragma omp parallel for schedule(dynamic) num_threads(nomp)
  for (int m = 0; m < natoms; m++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid, acc_device_nvidia);

    double Z1 = atnod[m];
    double A1 = coords[3 * m + 0]; double B1 = coords[3 * m + 1]; double C1 = coords[3 * m + 2];
    double coordn[6] = {0.};

    generate_ps_quad_grid(Z1, 1, coordn, quad_order, quad_order, nmu, nnu, nphi, grid, wt);
    add_r1_to_grid(gs, grid, 0., 0., 0.);

    printf("\n Double ** Problematic zone starts");

    for (int sp1 = 0; sp1 < n2ip[m].size() - 1; sp1++)
    for (int sp2 = sp1; sp2 < n2ip[m].size() - 1; sp2++)
    {
      int s1 = n2ip[m][sp1];     int s2 = n2ip[m][sp1 + 1];
      int s3 = n2ip[m][sp2];     int s4 = n2ip[m][sp2 + 1];

      #pragma acc parallel loop present(valS1[0:iN][0:gs3])
      for (int ii1 = 0; ii1 < s2 - s1; ii1++)
      {
        #pragma acc loop
        for (int j = 0; j < gs3; j++)
          valS1[ii1][j] = 1.;
      }

      #pragma acc parallel loop present(valS2[0:iN][0:gs3])
      for (int ii1 = 0; ii1 < s4 - s3; ii1++)
      {
        #pragma acc loop
        for (int j = 0; j < gs3; j++)
          valS2[ii1][j] = 1.;
      }

      // evaluate valS1
      for (int i1 = s1; i1 < s2; i1++)
      {
        int ii1 = i1 - s1;
        vector<double> basis1 = basis[i1];
        int n1 = basis1[0], l1 = basis1[1], m1 = basis1[2];
        double zeta1 = basis1[3];

        eval_pd(gs, grid, valS1[ii1], n1, l1, m1, zeta1);
      }

      // evaluate valS2
      for (int i2 = s3; i2 < s4; i2++)
      {
        int ii2 = i2 - s3;
        vector<double> basis1 = basis[i2];
        int n1 = basis1[0], l1 = basis1[1], m1 = basis1[2];
        double zeta1 = basis1[3];

        eval_pd(gs, grid, valS2[ii2], n1, l1, m1, zeta1);
      }

      #pragma acc wait

      //electron-nuclear attraction 
      // reduce_3cenp(Z1, s1, s2, s3, s4, N, iN, gs, grid, valS1, valS2, valt, wt, pVp);
    }

    printf("\n Double ** Problematic zone ends");
    exit(-1);
  }
}

#endif

*/

//Double ** ends VK edit


/////////END/////////


/*
  vector<vector<int> > n2ip;
  int imaxN = get_imax_n2ip(Nmax,natoms,N,basis,n2ip);
  if (prl>1)printf("   imaxN: %2i \n",imaxN);

  double gpumem_gb = gpumem/1024./1024./1024.;
  printf("   gpu memory available: %6.3f GB \n",gpumem_gb);

 //intermediate storage
  int iN = imaxN;
  double** valS1 = new double*[iN]; for (int i=0;i<iN;i++) valS1[i] = new double[gs3];
  double** valS2 = new double*[iN]; for (int i=0;i<iN;i++) valS2[i] = new double[gs3];
  double* valt = new double[gs];

 #if USE_ACC
 //#pragma omp parallel for schedule(static) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = n;
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc enter data create(pVp[0:N2])

    #pragma acc enter data create(grid[0:gs6],wt[0:gs])
    #pragma acc enter data create(valS1[0:iN][0:gs3],valS2[0:iN][0:gs3])
    #pragma acc enter data create(valt[0:gs])

    acc_assign(N2,pVp,0.);
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  double gpumem_2 = 1.*acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  printf("   after alloc, gpu memory available: %6.3f GB \n",gpumem_2*togb); 

 #pragma omp parallel for schedule(dynamic) num_threads(nomp)
  for (int m=0;m<natoms;m++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);

    double Z1 = atnod[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];
    double coordn[6];
    coordn[0] = coordn[1] = coordn[2] = 0.;

    generate_ps_quad_grid(Z1,1,coordn,quad_order,quad_order,nmu,nnu,nphi,grid,wt);
    add_r1_to_grid(gs,grid,0.,0.,0.);
    printf("\n Problematic zone starts");
    //j>i
    for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
    for (int sp2=sp1;sp2<n2ip[m].size()-1;sp2++)
    {
      int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
      int s3 = n2ip[m][sp2]; int s4 = n2ip[m][sp2+1];

      //printf(" pVp_2c s12: %2i %2i s34: %2i %2i \n",s1,s2,s3,s4);

     #pragma acc parallel loop present(valS1[0:iN][0:gs3])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          valS1[ii1][j] = 1.;
      }

     #pragma acc parallel loop present(valS2[0:iN][0:gs3])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          valS2[ii1][j] = 1.;
      }

     //first compute single atom ints
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_pd(gs,grid,valS1[ii1],n1,l1,m1,zeta1);
      }

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;

        vector<double> basis1 = basis[i2];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_pd(gs,grid,valS2[ii2],n1,l1,m1,zeta1);
      }

      #pragma acc wait

     /////////////////////////////////////////////////////////////////////
     //electron-nuclear attraction
      //reduce_3cenp(Z1,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,pVp);   VK tests
     /////////////////////////////////////////////////////////////////////

    } //loop sp over s12
    printf("\n Problematic zone ends");
    exit(-1);

   //two-atom ints
    for (int n=m+1;n<natoms;n++)
    {
      double Z2 = atnod[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;
      coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;

     //pull this out of loop later on
      generate_ps_quad_grid(0.,2,coordn,quad_order,quad_order,nmu,nnu,nphi,grid,wt);
      add_r1_to_grid(gs,grid,0.,0.,0.);

     //basis functions on each center
      for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
      for (int sp2=0;sp2<n2ip[n].size()-1;sp2++)
      {
        int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
        int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2+1];

        init_s12nw(s1,s2,s3,s4,iN,gs3,valS1,valS2);

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<double> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];
  
         //p
          eval_pd(gs,grid,valS1[ii1],n1,l1,m1,zeta1);
        }

       //second center
        recenter_grid_zero(gs,grid,-A12,-B12,-C12);

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

         //p
          eval_pd(gs,grid,valS2[ii2],n2,l2,m2,zeta2);
        }

       ///////////////////////////////////////////////////////////////////
       //electron-nuclear attraction centers 1+2

        reduce_3cenp(Z2,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,pVp);
        recenter_grid_zero(gs,grid,A12,B12,C12);
        reduce_3cenp(Z1,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,pVp);
        ////////////////////////////////////////////////////////////////////

      } //loop sp over s14

     //basis functions on center 1
      for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
      for (int sp2=sp1;sp2<n2ip[m].size()-1;sp2++)
      {
        int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
        int s3 = n2ip[m][sp2]; int s4 = n2ip[m][sp2+1];

        init_s12nw(s1,s2,s3,s4,iN,gs3,valS1,valS2);

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<double> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];
  
         //p
          eval_pd(gs,grid,valS1[ii1],n1,l1,m1,zeta1);
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

         //p
          eval_pd(gs,grid,valS2[ii2],n2,l2,m2,zeta2);
        }

       //second center
        recenter_grid_zero(gs,grid,-A12,-B12,-C12);

       ///////////////////////////////////////////////////////////////////
       //electron-nuclear attraction center 2
        reduce_3cenp(Z2,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,pVp);
        ////////////////////////////////////////////////////////////////////

        recenter_grid_zero(gs,grid,A12,B12,C12);

      } //loop sp over s14

     //basis functions on center 2
      for (int sp1=0;sp1<n2ip[n].size()-1;sp1++)
      for (int sp2=0;sp2<n2ip[n].size()-1;sp2++)
      {
        int s1 = n2ip[n][sp1]; int s2 = n2ip[n][sp1+1];
        int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2+1];

        init_s12nw(s1,s2,s3,s4,iN,gs3,valS1,valS2);

       //second center
        recenter_grid_zero(gs,grid,-A12,-B12,-C12);

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<double> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];
  
         //p
          eval_pd(gs,grid,valS1[ii1],n1,l1,m1,zeta1);
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

         //p
          eval_pd(gs,grid,valS2[ii2],n2,l2,m2,zeta2);
        }

        recenter_grid_zero(gs,grid,A12,B12,C12);

       ///////////////////////////////////////////////////////////////////
       //electron-nuclear attraction center 1
        reduce_3cenp(Z1,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,pVp);
        ////////////////////////////////////////////////////////////////////

      } //loop sp over s14

    } //loop n over second atom

  } //loop m over natoms
  acc_set_device_num(0,acc_device_nvidia);

  if (nomp>1)
  {
   //gather parallel fragments
    double pVpt[N2];
    for (int j=0;j<N2;j++)
      pVpt[j] = 0.;

    for (int n=0;n<nomp;n++)
    {
      int tid = n;
      acc_set_device_num(tid,acc_device_nvidia);

      #pragma acc update self(pVp[0:N2])

      for (int j=0;j<N2;j++)
        pVpt[j] += pVp[j];
    }

    for (int j=0;j<N2;j++)
      pVp[j] = pVpt[j];

    acc_set_device_num(0,acc_device_nvidia);
    #pragma acc update device(pVp[0:N2])
  }

  double norm[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];
  #pragma acc enter data copyin(norm[0:N])

  #pragma acc parallel loop collapse(2) independent present(pVp[0:N2],norm[0:N])
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
  {
    double n12 = norm[i]*norm[j];
    pVp[i*N+j] *= n12;
  }
  #pragma acc exit data delete(norm[0:N])

 //symmetrize wrt ij
 #pragma acc parallel loop independent present(pVp[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i+1;j<N;j++)
  {
    pVp[j*N+i] = pVp[i*N+j];
  }

  #pragma acc update self(pVp[0:N2])

  if (prl>1)
  {
    printf("\n pVp (2c): \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %12.6f",pVp[i*N+j]);
      printf("\n");
    }
  }

 #if USE_ACC
 //#pragma omp parallel for schedule(static) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = n;
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc exit data delete(pVp[0:N2])
    #pragma acc exit data delete(grid[0:gs6],wt[0:gs])
    #pragma acc exit data delete(valS1[0:iN][0:gs3],valS2[0:iN][0:gs3])
    #pragma acc exit data delete(valt[0:gs])
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  for (int i=0;i<iN;i++) delete [] valS1[i];
  for (int i=0;i<iN;i++) delete [] valS2[i];
  delete [] valS1; delete [] valS2;
  delete [] valt;

  delete [] grid;
  delete [] wt;

  return;
}

*/

void compute_2c_ps(bool do_overlap, bool do_yukawa, double gamma, int natoms, int* atno, double* coords, vector<vector<double> > &basis, int quad_order, int nmu, int nnu, int nphi, double* A, int prl)
{
  if (do_overlap && prl>1) { printf("\n WARNING: testing do_overlap in compute_2c_ps \n"); }

  if (prl>-1) { if (do_yukawa) printf("  beginning compute_2c_ps (Yukawa. gamma: %5.3f) \n",gamma); else printf("  beginning compute_2c_ps \n"); }

  int nomp_max = 1;
 #pragma omp parallel
  nomp_max = omp_get_num_threads();

  int ngpu = 0;
 #if USE_ACC
  ngpu = acc_get_num_devices(acc_device_nvidia);
  int nomp = ngpu;
 #else
  int nomp = nomp_max;
 #endif

  bool dy = do_yukawa;

  int N = basis.size();
  int N2 = N*N;

  if (N<1) { printf(" ERROR: cannot compute 2c integrals, no RI basis functions \n"); exit(-1); }

  //double atnod[natoms];
  //get_atnod(natoms,basis,atnod);

  int qos = quad_order*quad_order*quad_order;
  int gs = nmu*nnu*nphi*qos;
  int gs6 = 6*gs;

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

 //expand volume?
  double cfn = 1.1;
  int lmax = 0;
  double ztminl[10];
  for (int l=0;l<10;l++) ztminl[l] = 10000.;
  for (int l=0;l<10;l++)
  {
    for (int j=0;j<N;j++)
    {
      int l1 = basis[j][1];
      if (l1>lmax) lmax = l1;
      if (l1==l && basis[j][3]<ztminl[l])
        ztminl[l] = basis[j][3];
    }
    if (lmax==l && prl>1)
      printf("   l: %i ztmin: %8.5f \n",l,ztminl[l]);
  }

  cfn = 1.;
  //if (maxl>3) cfn = 1.1;

  //printf("   cfn: %5.2f \n",cfn);

  double* grid = new double[gs6];
  double* wt = new double[gs];

  double gpumem = (double)acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  double togb = 1./1024./1024./1024.;

 //this calculation not accurate yet
  int Nmax = 300;
  double mem0 = gs*7.+1.*nmu*nnu*nphi;
  while (Nmax>0)
  {
    //double mem1 = 8.*(gs*2.*Nmax + 1.*mem0);
    double mem1 = 8.*(gs*2.*Nmax + 1.*mem0 + 1.*N2 + 2.*gs6 + 2.*gs);
    if (mem1<gpumem)
    {
      if (prl>1) printf("    mem0: %5.3f mem1: %5.3f \n",mem0*togb,mem1*togb);
      break;
    }
    Nmax--;
  }
  //if (Nmax>5)
  //  Nmax -= 5;
  if (Nmax<=0) { printf("\n ERROR: couldn't calculate gpu memory requirements \n"); exit(-1); }

  vector<vector<int> > n2ip;
  int imaxN = get_imax_n2ip(Nmax,natoms,N,basis,n2ip);
  if (prl>1)printf("   imaxN: %2i \n",imaxN);

  double gsxvalsv = 8.*(imaxN*gs*2. + 1.*mem0); //vals+grid/wt+gridm
  double gsxvalsv_gb = gsxvalsv/1024./1024./1024.;
  printf("   estimated memory needed (2c): %6.3f GB \n",gsxvalsv_gb);

  double gpumem_gb = gpumem/1024./1024./1024.;
  printf("   gpu memory available: %6.3f GB \n",gpumem_gb);

//New code
 //intermediate storage
  int iN = imaxN;
  double* valS1 = new double[iN * gs];
  double* valV2 = new double[iN * gs];


 #if USE_ACC
 #pragma omp parallel for schedule(static) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = n;
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc enter data create(A[0:N2])
    #pragma acc enter data copyin(coords[0:3*natoms],atno[0:natoms])

    #pragma acc enter data create(grid[0:gs6],wt[0:gs])
    #pragma acc enter data create(valS1[0:iN * gs], valV2[0:iN * gs])


    acc_assign(N2,A,0.);
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  double gpumem_2 = 1.*acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  printf("   after alloc, gpu memory available: %6.3f GB \n",gpumem_2*togb); 

 #pragma omp parallel for schedule(dynamic) num_threads(nomp)
for (int m = 0; m < natoms; m++) {
  int tid = omp_get_thread_num();
  acc_set_device_num(tid, acc_device_nvidia);

  double Z1 = (double)atno[m];
  double A1 = coords[3 * m + 0];
  double B1 = coords[3 * m + 1];
  double C1 = coords[3 * m + 2];
  double coordn[6] = {0.};

  generate_ps_quad_grid(cfn, Z1, 1, coordn, quad_order, quad_order, nmu, nnu, nphi, grid, wt);
  add_r1_to_grid(gs, grid, 0., 0., 0.);

  for (int sp1 = 0; sp1 < n2ip[m].size() - 1; sp1++)
  for (int sp2 = sp1; sp2 < n2ip[m].size() - 1; sp2++) {
    int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1 + 1];
    int s3 = n2ip[m][sp2]; int s4 = n2ip[m][sp2 + 1];

    if (dy) {
      #pragma acc parallel loop collapse(2) present(valS1[0:iN * gs], wt[0:gs])
      for (int ii1 = 0; ii1 < s2 - s1; ii1++)
      for (int j = 0; j < gs; j++)
        valS1[ii1 * gs + j] = wt[j];

      #pragma acc parallel loop collapse(2) present(valV2[0:iN * gs])
      for (int ii2 = 0; ii2 < s4 - s3; ii2++)
      for (int j = 0; j < gs; j++)
        valV2[ii2 * gs + j] = 0.;
    } else {
      #pragma acc parallel loop collapse(2) present(valS1[0:iN * gs])
      for (int ii1 = 0; ii1 < s2 - s1; ii1++)
      for (int j = 0; j < gs; j++)
        valS1[ii1 * gs + j] = 1.;

      #pragma acc parallel loop collapse(2) present(valV2[0:iN * gs], wt[0:gs])
      for (int ii2 = 0; ii2 < s4 - s3; ii2++)
      for (int j = 0; j < gs; j++)
        valV2[ii2 * gs + j] = wt[j];
    }

    for (int i1 = s1; i1 < s2; i1++) {
      int ii1 = i1 - s1;
      vector<double> basis1 = basis[i1];
      int n1 = basis1[0], l1 = basis1[1], m1 = basis1[2]; double zeta1 = basis1[3];
      eval_shd(ii1, gs, grid, &valS1[ii1 * gs], n1, l1, m1, zeta1);
    }

    for (int i2 = s3; i2 < s4; i2++) {
      int ii2 = i2 - s3;
      vector<double> basis2 = basis[i2];
      int n2 = basis2[0], l2 = basis2[1], m2 = basis2[2]; double zeta2 = basis2[3];

      if (do_overlap)
        eval_shd(ii2, gs, grid, &valV2[ii2 * gs], n2, l2, m2, zeta2);
      else {
        if (dy)
          eval_inr_yukawa(gs, grid, &valV2[ii2 * gs], n2, l2, zeta2, gamma);
        else
          eval_inr_r12(gs, grid, &valV2[ii2 * gs], n2, l2, zeta2, ii2);
        eval_sh_3rd(gs, grid, &valV2[ii2 * gs], n2, l2, m2);
      }
    }

    #pragma acc wait
    reduce_2c1(s1, s2, s3, s4, gs, valS1, valV2, iN, N, A);
  }

  for (int n = m + 1; n < natoms; n++) {
    double Z2 = (double)atno[n];
    double A2 = coords[3 * n + 0]; double B2 = coords[3 * n + 1]; double C2 = coords[3 * n + 2];
    double A12 = A2 - A1; double B12 = B2 - B1; double C12 = C2 - C1;
    coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;

    generate_ps_quad_grid(cfn, 0., 2, coordn, quad_order, quad_order, nmu, nnu, nphi, grid, wt);
    add_r1_to_grid(gs, grid, 0., 0., 0.);

    for (int sp1 = 0; sp1 < n2ip[m].size() - 1; sp1++)
    for (int sp2 = 0; sp2 < n2ip[n].size() - 1; sp2++) {
      int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1 + 1];
      int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2 + 1];

      if (dy) {
        #pragma acc parallel loop collapse(2) present(valS1[0:iN * gs], wt[0:gs])
        for (int ii1 = 0; ii1 < s2 - s1; ii1++)
        for (int j = 0; j < gs; j++)
          valS1[ii1 * gs + j] = wt[j];

        #pragma acc parallel loop collapse(2) present(valV2[0:iN * gs])
        for (int ii2 = 0; ii2 < s4 - s3; ii2++)
        for (int j = 0; j < gs; j++)
          valV2[ii2 * gs + j] = 0.;
      } else {
        #pragma acc parallel loop collapse(2) present(valS1[0:iN * gs])
        for (int ii1 = 0; ii1 < s2 - s1; ii1++)
        for (int j = 0; j < gs; j++)
          valS1[ii1 * gs + j] = 1.;

        #pragma acc parallel loop collapse(2) present(valV2[0:iN * gs], wt[0:gs])
        for (int ii2 = 0; ii2 < s4 - s3; ii2++)
        for (int j = 0; j < gs; j++)
          valV2[ii2 * gs + j] = wt[j];
      }

      for (int i1 = s1; i1 < s2; i1++) {
        int ii1 = i1 - s1;
        vector<double> basis1 = basis[i1];
        int n1 = basis1[0], l1 = basis1[1], m1 = basis1[2]; double zeta1 = basis1[3];
        eval_shd(ii1, gs, grid, &valS1[ii1 * gs], n1, l1, m1, zeta1);
      }

      recenter_grid_zero(gs, grid, -A12, -B12, -C12);

      for (int i2 = s3; i2 < s4; i2++) {
        int ii2 = i2 - s3;
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0], l2 = basis2[1], m2 = basis2[2]; double zeta2 = basis2[3];

        if (do_overlap)
          eval_shd(ii2, gs, grid, &valV2[ii2 * gs], n2, l2, m2, zeta2);
        else {
          if (dy)
            eval_inr_yukawa(gs, grid, &valV2[ii2 * gs], n2, l2, zeta2, gamma);
          else
            eval_inr_r12(gs, grid, &valV2[ii2 * gs], n2, l2, zeta2, ii2);
          eval_sh_3rd(gs, grid, &valV2[ii2 * gs], n2, l2, m2);
        }
      }

      reduce_2c1(s1, s2, s3, s4, gs, valS1, valV2, iN, N, A);
      recenter_grid_zero(gs, grid, A12, B12, C12);
    }
  }
}
  acc_set_device_num(0, acc_device_nvidia);

if (nomp > 1) {
  double* At = new double[N2]();

  for (int n = 0; n < nomp; n++) {
    acc_set_device_num(n, acc_device_nvidia);
    #pragma acc update self(A[0:N2])
    for (int j = 0; j < N2; j++) At[j] += A[j];
  }

  for (int j = 0; j < N2; j++) A[j] = At[j];

  delete[] At;
  acc_set_device_num(0, acc_device_nvidia);
  #pragma acc update device(A[0:N2])
}

double* norm1 = new double[N];
double* norm2 = new double[N];
if (do_overlap) {
  for (int i = 0; i < N; i++)
    norm1[i] = norm(basis[i][0], basis[i][1], basis[i][2], basis[i][3]);
} else {
  for (int i = 0; i < N; i++)
    norm1[i] = norm_sv(basis[i][0], basis[i][1], basis[i][2], basis[i][3]);
}
for (int i = 0; i < N; i++)
  norm2[i] = basis[i][4];
#pragma acc enter data copyin(norm1[0:N], norm2[0:N])

#pragma acc parallel loop independent present(A[0:N2], norm1[0:N], norm2[0:N])
for (int i = 0; i < N; i++)
#pragma acc loop independent
for (int j = 0; j < N; j++) {
  double n12 = norm2[i] * norm1[j];
  A[i * N + j] *= n12;
}

#if 0
printf(" checking asymm \n");
double ct = 1.e-3;
#pragma acc update self(A[0:N2])
for (int i = 0; i < N; i++)
for (int j = 0; j < i; j++)
if (fabs(A[i * N + j] - A[j * N + i]) > ct)
  printf("  %8.3e %8.3e  diff: %8.3e \n", A[i * N + j], A[j * N + i], A[i * N + j] - A[j * N + i]);
#endif

// symmetrize wrt ij
#pragma acc parallel loop independent present(A[0:N2])
for (int i = 0; i < N; i++)
#pragma acc loop independent
for (int j = i + 1; j < N; j++)
  A[j * N + i] = A[i * N + j];

double lt = 1.e-15;
#pragma acc parallel loop present(A[0:N2])
for (int j = 0; j < N2; j++)
  if (fabs(A[j]) < lt)
    A[j] = 0.;

#pragma acc exit data delete(norm1[0:N], norm2[0:N])
#pragma acc update self(A[0:N2])

if (prl > 1 || prl == -1) {
  printf("\n A: \n");
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++)
      printf(" %12.6f", A[i * N + j]);
    printf("\n");
  }
}

#if USE_ACC
#pragma omp parallel for schedule(static) num_threads(nomp)
for (int n = 0; n < nomp; n++) {
  int tid = n;
  acc_set_device_num(tid, acc_device_nvidia);

  #pragma acc exit data delete(A[0:N2])
  #pragma acc exit data delete(grid[0:gs6], wt[0:gs])
  #pragma acc exit data delete(valS1[0:iN * gs], valV2[0:iN * gs])
  #pragma acc exit data delete(coords[0:3 * natoms], atno[0:natoms])
}
acc_set_device_num(0, acc_device_nvidia);
#endif

delete[] valS1;
delete[] valV2;
delete[] grid;
delete[] wt;
delete[] norm1;
delete[] norm2;

return;
}
//ends

/*
 //intermediate storage
  int iN = imaxN;
  double** valS1 = new double*[iN]; for (int i=0;i<iN;i++) valS1[i] = new double[gs];
  double** valV2 = new double*[iN]; for (int i=0;i<iN;i++) valV2[i] = new double[gs];

 #if USE_ACC
 #pragma omp parallel for schedule(static) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = n;
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc enter data create(A[0:N2])
    #pragma acc enter data copyin(coords[0:3*natoms],atno[0:natoms])

    #pragma acc enter data create(grid[0:gs6],wt[0:gs])
    #pragma acc enter data create(valS1[0:iN][0:gs],valV2[0:iN][0:gs])

    acc_assign(N2,A,0.);
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  double gpumem_2 = 1.*acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  printf("   after alloc, gpu memory available: %6.3f GB \n",gpumem_2*togb); 

 #pragma omp parallel for schedule(dynamic) num_threads(nomp)
  for (int m=0;m<natoms;m++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);

    double Z1 = (double)atno[m];
    //double Z1 = atnod[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];
    double coordn[6];
    coordn[0] = 0.; coordn[1] = 0.; coordn[2] = 0.;

    generate_ps_quad_grid(cfn,Z1,1,coordn,quad_order,quad_order,nmu,nnu,nphi,grid,wt);
    add_r1_to_grid(gs,grid,0.,0.,0.);

    for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
    for (int sp2=sp1;sp2<n2ip[m].size()-1;sp2++)
    {
      int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
      int s3 = n2ip[m][sp2]; int s4 = n2ip[m][sp2+1];
      //printf("    m: %i  s1/2: %3i %3i s3/4: %3i %3i \n",m,s1,s2,s3,s4);

      if (dy)
      {
       #pragma acc parallel loop collapse(2) present(valS1[0:iN][0:gs],wt[0:gs])
        for (int ii1=0;ii1<s2-s1;ii1++)
        for (int j=0;j<gs;j++)
          valS1[ii1][j] = wt[j];

       #pragma acc parallel loop collapse(2) present(valV2[0:iN][0:gs])
        for (int ii2=0;ii2<s4-s3;ii2++)
        for (int j=0;j<gs;j++)
          valV2[ii2][j] = 0.;
      }
      else
      {
       #pragma acc parallel loop collapse(2) present(valS1[0:iN][0:gs])
        for (int ii1=0;ii1<s2-s1;ii1++)
        for (int j=0;j<gs;j++)
          valS1[ii1][j] = 1.;

       #pragma acc parallel loop collapse(2) present(valV2[0:iN][0:gs],wt[0:gs])
        for (int ii2=0;ii2<s4-s3;ii2++)
        for (int j=0;j<gs;j++)
          valV2[ii2][j] = wt[j];
      }

     //first compute single atom ints
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_shd(ii1,gs,grid,valS1[ii1],n1,l1,m1,zeta1);
      }

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        if (do_overlap)
        {
          eval_shd(ii2,gs,grid,valV2[ii2],n2,l2,m2,zeta2);
        }
        else
        {
         //V
          if (dy)
            eval_inr_yukawa(gs,grid,valV2[ii2],n2,l2,zeta2,gamma);
          else
            eval_inr_r12(gs,grid,valV2[ii2],n2,l2,zeta2,ii2);
          eval_sh_3rd(gs,grid,valV2[ii2],n2,l2,m2);
        }
      } //loop i2 evaluate

      #pragma acc wait

      reduce_2c1(s1,s2,s3,s4,gs,valS1,valV2,iN,N,A);
    } //sp partition over m


   //two-atom ints
    for (int n=m+1;n<natoms;n++)
    {
      double Z2 = (double)atno[n];
      //double Z2 = atnod[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;
      coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;

      generate_ps_quad_grid(cfn,0.,2,coordn,quad_order,quad_order,nmu,nnu,nphi,grid,wt);
      add_r1_to_grid(gs,grid,0.,0.,0.);

      for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
      for (int sp2=0;sp2<n2ip[n].size()-1;sp2++)
      {  
        int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
        int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2+1];
        //printf("    mn: %i %i s1/2: %3i %3i s3/4: %3i %3i \n",m,n,s1,s2,s3,s4);

        if (dy)
        {
         #pragma acc parallel loop collapse(2) present(valS1[0:iN][0:gs],wt[0:gs])
          for (int ii1=0;ii1<s2-s1;ii1++)
          for (int j=0;j<gs;j++)
            valS1[ii1][j] = wt[j];
         #pragma acc parallel loop collapse(2) present(valV2[0:iN][0:gs])
          for (int ii1=0;ii1<s4-s3;ii1++)
          for (int j=0;j<gs;j++)
            valV2[ii1][j] = 0.;
        }
        else
        {
         #pragma acc parallel loop collapse(2) present(valS1[0:iN][0:gs])
          for (int ii1=0;ii1<s2-s1;ii1++)
          for (int j=0;j<gs;j++)
            valS1[ii1][j] = 1.;
         #pragma acc parallel loop collapse(2) present(valV2[0:iN][0:gs],wt[0:gs])
          for (int ii1=0;ii1<s4-s3;ii1++)
          for (int j=0;j<gs;j++)
            valV2[ii1][j] = wt[j];
        }

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<double> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

         //S
          eval_shd(ii1,gs,grid,valS1[ii1],n1,l1,m1,zeta1);
        }

       //second center
        recenter_grid_zero(gs,grid,-A12,-B12,-C12);

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

          if (do_overlap)
          {
            eval_shd(ii2,gs,grid,valV2[ii2],n2,l2,m2,zeta2);
          }
          else
          {
           //V
            if (dy)
              eval_inr_yukawa(gs,grid,valV2[ii2],n2,l2,zeta2,gamma);
            else
              eval_inr_r12(gs,grid,valV2[ii2],n2,l2,zeta2,ii2);
            eval_sh_3rd(gs,grid,valV2[ii2],n2,l2,m2);
          }
        }

        reduce_2c1(s1,s2,s3,s4,gs,valS1,valV2,iN,N,A);

       //put grid back in place?
        recenter_grid_zero(gs,grid,A12,B12,C12);

      } //sp partition over n
    } //loop n over second atom

  } //loop m over natoms
  acc_set_device_num(0,acc_device_nvidia);


  if (nomp>1)
  {
   //gather parallel
    double At[N2];
    for (int j=0;j<N2;j++) At[j] = 0.;

    for (int n=0;n<nomp;n++)
    {
      int tid = n;
      acc_set_device_num(tid,acc_device_nvidia);

      #pragma acc update self(A[0:N2])

      for (int j=0;j<N2;j++)
        At[j] += A[j];
    }

    for (int j=0;j<N2;j++)
      A[j] = At[j];

    acc_set_device_num(0,acc_device_nvidia);
    #pragma acc update device(A[0:N2])
  }

  double norm1[N];
  double norm2[N];
  if (do_overlap)
  for (int i=0;i<N;i++)
    norm1[i] = norm(basis[i][0],basis[i][1],basis[i][2],basis[i][3]);
  else
  for (int i=0;i<N;i++)
    norm1[i] = norm_sv(basis[i][0],basis[i][1],basis[i][2],basis[i][3]);
  for (int i=0;i<N;i++)
    norm2[i] = basis[i][4];
  #pragma acc enter data copyin(norm1[0:N],norm2[0:N])

  #pragma acc parallel loop independent present(A[0:N2],norm1[0:N],norm2[0:N])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=0;j<N;j++)
  {
    double n12 = norm2[i]*norm1[j];
    A[i*N+j] *= n12;
  }

 #if 0
  printf(" checking asymm \n");
  double ct = 1.e-3;
  #pragma acc update self(A[0:N2])
  for (int i=0;i<N;i++)
  for (int j=0;j<i;j++)
  if (fabs(A[i*N+j]-A[j*N+i])>ct)
    printf("  %8.3e %8.3e  diff: %8.3e \n",A[i*N+j],A[j*N+i],A[i*N+j]-A[j*N+i]);
 #endif

 //symmetrize wrt ij
 #pragma acc parallel loop independent present(A[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i+1;j<N;j++)
    A[j*N+i] = A[i*N+j];

  double lt = 1.e-15;
 #pragma acc parallel loop present(A[0:N2])
  for (int j=0;j<N2;j++)
  if (fabs(A[j])<lt)
    A[j] = 0.;

  #pragma acc exit data delete(norm1[0:N],norm2[0:N])

  #pragma acc update self(A[0:N2])


  if (prl>1 || prl==-1)
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
 #pragma omp parallel for schedule(static) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = n;
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc exit data delete(A[0:N2])
    #pragma acc exit data delete(grid[0:gs6],wt[0:gs])
    #pragma acc exit data delete(valS1[0:iN][0:gs],valV2[0:iN][0:gs])
    #pragma acc exit data delete(coords[0:3*natoms],atno[0:natoms])
  }
  acc_set_device_num(0,acc_device_nvidia);
#endif

  //printf(" done with dealloc in 2c integrals \n"); fflush(stdout);
  //auto_crash();

  for (int i=0;i<iN;i++) delete [] valS1[i];
  for (int i=0;i<iN;i++) delete [] valV2[i];
  delete [] valS1; delete [] valV2;

  delete [] grid;
  delete [] wt;

  return;
}
*/

//overloaded with double *
void init_s12v3(bool dy, int s1, int s2, int s3, int s4, int s5, int s6, int iN, int iNa, int gs,
                double* val1, double* val2, double* val3, double* wt)
{
  // Initialize val2 to 1.0
  #pragma acc parallel loop collapse(2) present(val2[0:iN * gs])
  for (int ii2 = 0; ii2 < s4 - s3; ii2++) {
    for (int j = 0; j < gs; j++) {
      val2[ii2 * gs + j] = 1.0;
    }
  }

  if (dy) {
    // Set val1 = wt
    #pragma acc parallel loop collapse(2) present(val1[0:iN * gs], wt[0:gs])
    for (int ii1 = 0; ii1 < s2 - s1; ii1++) {
      for (int j = 0; j < gs; j++) {
        val1[ii1 * gs + j] = wt[j];
      }
    }

    // Set val3 = 0
    if (val3 != nullptr) {
      #pragma acc parallel loop collapse(2) present(val3[0:iNa * gs])
      for (int ii3 = 0; ii3 < s6 - s5; ii3++) {
        for (int j = 0; j < gs; j++) {
          val3[ii3 * gs + j] = 0.0;
        }
      }
    }
  } else {
    // Set val1 = 1.0
    #pragma acc parallel loop collapse(2) present(val1[0:iN * gs])
    for (int ii1 = 0; ii1 < s2 - s1; ii1++) {
      for (int j = 0; j < gs; j++) {
        val1[ii1 * gs + j] = 1.0;
      }
    }

    // Set val3 = wt
    if (val3 != nullptr) {
      #pragma acc parallel loop collapse(2) present(val3[0:iNa * gs], wt[0:gs])
      for (int ii3 = 0; ii3 < s6 - s5; ii3++) {
        for (int j = 0; j < gs; j++) {
          val3[ii3 * gs + j] = wt[j];
        }
      }
    }
  }

  return;
}


void init_s12v3(bool dy, int s1, int s2, int s3, int s4, int s5, int s6, int iN, int iNa, int gs, double** val1, double** val2, double** val3, double* wt)
{
  #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs])
  for (int ii2=0;ii2<s4-s3;ii2++)
  {
    for (int j=0;j<gs;j++)
      val2[ii2][j] = 1.;
  }

  if (dy)
  {
    #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs],wt[0:gs])
    for (int ii1=0;ii1<s2-s1;ii1++)
    for (int j=0;j<gs;j++)
      val1[ii1][j] = wt[j];

    if (val3!=NULL)
    #pragma acc parallel loop collapse(2) present(val3[0:iNa][0:gs])
    for (int ii3=0;ii3<s6-s5;ii3++)
    for (int j=0;j<gs;j++)
      val3[ii3][j] = 0.;
  }
  else
  {
    #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs])
    for (int ii1=0;ii1<s2-s1;ii1++)
    for (int j=0;j<gs;j++)
      val1[ii1][j] = 1.;

    if (val3!=NULL)
    #pragma acc parallel loop collapse(2) present(val3[0:iNa][0:gs],wt[0:gs])
    for (int ii3=0;ii3<s6-s5;ii3++)
    for (int j=0;j<gs;j++)
      val3[ii3][j] = wt[j];
  }

  return;
}

//overloaded and flattened
void eval_s12v3(bool dol, bool dy, double gamma,
                int s1, int s2, int s3, int s4, int s5, int s6,
                int gs, double* grid,
                const std::vector<std::vector<double>>& basis,
                const std::vector<std::vector<double>>& basis_aux,
                double* val1, double* val2, double* val3)
{
  // val1: (s2 - s1) * gs
  for (int i1 = s1; i1 < s2; i1++) {
    int ii1 = i1 - s1;
    const std::vector<double>& b1 = basis[i1];
    int n1 = b1[0], l1 = b1[1], m1 = b1[2]; double zeta1 = b1[3];
    eval_shd(ii1, gs, grid, &val1[ii1 * gs], n1, l1, m1, zeta1);
  }

  // val2: (s4 - s3) * gs
  for (int i2 = s3; i2 < s4; i2++) {
    int ii2 = i2 - s3;
    const std::vector<double>& b2 = basis[i2];
    int n2 = b2[0], l2 = b2[1], m2 = b2[2]; double zeta2 = b2[3];
    eval_shd(ii2, gs, grid, &val2[ii2 * gs], n2, l2, m2, zeta2);
  }

  // val3: (s6 - s5) * gs
  for (int i3 = s5; i3 < s6; i3++) {
    int ii3 = i3 - s5;
    const std::vector<double>& b3 = basis_aux[i3];
    int n3 = b3[0], l3 = b3[1], m3 = b3[2]; double zeta3 = b3[3];

    if (dol) {
      eval_shd(ii3, gs, grid, &val3[ii3 * gs], n3, l3, m3, zeta3);
    } else {
      if (dy)
        eval_inr_yukawa(gs, grid, &val3[ii3 * gs], n3, l3, zeta3, gamma);
      else
        eval_inr_r12(gs, grid, &val3[ii3 * gs], n3, l3, zeta3, ii3);
      eval_sh_3rd(gs, grid, &val3[ii3 * gs], n3, l3, m3);
    }
  }

  #pragma acc wait
  return;
}


void eval_s12v3(bool dol, bool dy, double gamma, int s1, int s2, int s3, int s4, int s5, int s6, int gs, double* grid, vector<vector<double> >& basis, vector<vector<double> >& basis_aux, double** val1, double** val2, double** val3)
{
 //single-center evaluations

  for (int i1=s1;i1<s2;i1++)
  {
    int ii1 = i1-s1;

    vector<double> basis1 = basis[i1];
    int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

    eval_shd(ii1,gs,grid,val1[ii1],n1,l1,m1,zeta1);
  }

  for (int i2=s3;i2<s4;i2++)
  {
    int ii2 = i2-s3;

    vector<double> basis2 = basis[i2];
    int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

    eval_shd(ii2,gs,grid,val2[ii2],n2,l2,m2,zeta2);
  }

  for (int i3=s5;i3<s6;i3++)
  {
    int ii3 = i3-s5;
    vector<double> basis3 = basis_aux[i3];
    int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

    if (dol)
      eval_shd(ii3,gs,grid,val3[ii3],n3,l3,m3,zeta3);
    else
    {
     //V
      if (dy)
        eval_inr_yukawa(gs,grid,val3[ii3],n3,l3,zeta3,gamma);
      else
        eval_inr_r12(gs,grid,val3[ii3],n3,l3,zeta3,ii3);
      eval_sh_3rd (gs,grid,val3[ii3],n3,l3,m3);
    }
  }

  #pragma acc wait

  return;
}

//overloaded and flattened
void eval_s12v3_2(bool dol, bool dy, double gamma,
                  int s1, int s2, int s3, int s4, int s5, int s6,
                  int gs, double* grid,
                  const std::vector<std::vector<double>>& basis,
                  const std::vector<std::vector<double>>& basis_aux,
                  double* val1, double* val2, double* val3,
                  int type, double A12, double B12, double C12,
                  double A13, double B13, double C13)
{
  // val1: (s2 - s1) * gs
  for (int i1 = s1; i1 < s2; i1++) {
    int ii1 = i1 - s1;
    const std::vector<double>& b1 = basis[i1];
    int n1 = b1[0], l1 = b1[1], m1 = b1[2]; double zeta1 = b1[3];
    eval_shd(ii1, gs, grid, &val1[ii1 * gs], n1, l1, m1, zeta1);
  }

  if (type == 1 || type == 3)
    recenter_grid_zero(gs, grid, -A12, -B12, -C12);
  if (type == 4)
    recenter_grid_zero(gs, grid, -A13, -B13, -C13);

  // val2: (s4 - s3) * gs
  for (int i2 = s3; i2 < s4; i2++) {
    int ii2 = i2 - s3;
    const std::vector<double>& b2 = basis[i2];
    int n2 = b2[0], l2 = b2[1], m2 = b2[2]; double zeta2 = b2[3];
    eval_shd(ii2, gs, grid, &val2[ii2 * gs], n2, l2, m2, zeta2);
  }

  if (type == 2)
    recenter_grid_zero(gs, grid, -A12, -B12, -C12);
  if (type == 3)
    recenter_grid_zero(gs, grid, A12 - A13, B12 - B13, C12 - C13);

  // val3: (s6 - s5) * gs
  if (val3 != nullptr) {
    for (int i3 = s5; i3 < s6; i3++) {
      int ii3 = i3 - s5;
      const std::vector<double>& b3 = basis_aux[i3];
      int n3 = b3[0], l3 = b3[1], m3 = b3[2]; double zeta3 = b3[3];

      if (dol) {
        eval_shd(ii3, gs, grid, &val3[ii3 * gs], n3, l3, m3, zeta3);
      } else {
        if (dy)
          eval_inr_yukawa(gs, grid, &val3[ii3 * gs], n3, l3, zeta3, gamma);
        else
          eval_inr_r12(gs, grid, &val3[ii3 * gs], n3, l3, zeta3, ii3);
        eval_sh_3rd(gs, grid, &val3[ii3 * gs], n3, l3, m3);
      }
    }
  }

  // reset grid
  if (type == 1 || type == 2)
    recenter_grid_zero(gs, grid, A12, B12, C12);
  if (type == 3)
    recenter_grid_zero(gs, grid, A13, B13, C13);

  #pragma acc wait
  return;
}


void eval_s12v3_2(bool dol, bool dy, double gamma, int s1, int s2, int s3, int s4, int s5, int s6, int gs, double* grid, vector<vector<double> >& basis, vector<vector<double> >& basis_aux, double** val1, double** val2, double** val3, int type, double A12, double B12, double C12, double A13, double B13, double C13)
{
 //multi-center evaluations

  for (int i1=s1;i1<s2;i1++)
  {
    int ii1 = i1-s1;

    vector<double> basis1 = basis[i1];
    int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

    eval_shd(ii1,gs,grid,val1[ii1],n1,l1,m1,zeta1);
  }

  if (type==1 || type==3) //i2/i3 on second atom
    recenter_grid_zero(gs,grid,-A12,-B12,-C12);
  if (type==4) //i2/i3 on atom 3
    recenter_grid_zero(gs,grid,-A13,-B13,-C13);

  for (int i2=s3;i2<s4;i2++)
  {
    int ii2 = i2-s3;

    vector<double> basis2 = basis[i2];
    int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

    eval_shd(ii2,gs,grid,val2[ii2],n2,l2,m2,zeta2);
  }

  if (type==2) //i3 only on second atom
    recenter_grid_zero(gs,grid,-A12,-B12,-C12);

  if (type==3)
    recenter_grid_zero(gs,grid,A12-A13,B12-B13,C12-C13);

  if (val3!=NULL)
  for (int i3=s5;i3<s6;i3++)
  {
    int ii3 = i3-s5;
    vector<double> basis3 = basis_aux[i3];
    int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

    if (dol)
      eval_shd(ii3,gs,grid,val3[ii3],n3,l3,m3,zeta3);
    else
    {
     //V
      if (dy)
        eval_inr_yukawa(gs,grid,val3[ii3],n3,l3,zeta3,gamma);
      else
        eval_inr_r12(gs,grid,val3[ii3],n3,l3,zeta3,ii3);
      eval_sh_3rd (gs,grid,val3[ii3],n3,l3,m3);
    }
  }

 //return grid to original position
  if (type==1 || type==2)
    recenter_grid_zero(gs,grid,A12,B12,C12);
  if (type==3) //type==4 does not reset
    recenter_grid_zero(gs,grid,A13,B13,C13);

  #pragma acc wait

  return;
}

//overloaded and flattened

void eval_p12(int s1, int s2, int s3, int s4, int gs, double* grid,
              std::vector<std::vector<double>>& basis,
              double* val1, double* val2,
              double A12, double B12, double C12,
              double A13, double B13, double C13)
{
  for (int i1 = s1; i1 < s2; i1++)
  {
    int ii1 = i1 - s1;

    std::vector<double> basis1 = basis[i1];
    int n1 = basis1[0];
    int l1 = basis1[1];
    int m1 = basis1[2];
    double zeta1 = basis1[3];

    eval_pd(gs, grid, &val1[ii1 * gs], n1, l1, m1, zeta1);
  }

  recenter_grid_zero(gs, grid, -A12, -B12, -C12);

  for (int i2 = s3; i2 < s4; i2++)
  {
    int ii2 = i2 - s3;

    std::vector<double> basis2 = basis[i2];
    int n2 = basis2[0];
    int l2 = basis2[1];
    int m2 = basis2[2];
    double zeta2 = basis2[3];

    eval_pd(gs, grid, &val2[ii2 * gs], n2, l2, m2, zeta2);
  }

  recenter_grid_zero(gs, grid, A12, B12, C12);

  #pragma acc wait

  return;
}


void eval_p12(int s1, int s2, int s3, int s4, int gs, double* grid, vector<vector<double> >& basis, double** val1, double** val2, double A12, double B12, double C12, double A13, double B13, double C13)
{
  for (int i1=s1;i1<s2;i1++)
  {
    int ii1 = i1-s1;

    vector<double> basis1 = basis[i1];
    int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

    eval_pd(gs,grid,val1[ii1],n1,l1,m1,zeta1);
  }

  recenter_grid_zero(gs,grid,-A12,-B12,-C12);

  for (int i2=s3;i2<s4;i2++)
  {
    int ii2 = i2-s3;

    vector<double> basis2 = basis[i2];
    int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

    eval_pd(gs,grid,val2[ii2],n2,l2,m2,zeta2);
  }

  recenter_grid_zero(gs,grid,A12,B12,C12);

  #pragma acc wait

  return;
}

void compute_pVp_3c_ps(int natoms, int* atno, double* coords, vector<vector<double> > &basis, int quad_order, int quad_r_order, int nsplit, int nmu, int nnu, int nphi, double* pVp, int prl)
{
  if (prl>-1) printf("  beginning compute_pVp_3c_ps \n");

  int nomp_max = 1;
 #pragma omp parallel
  nomp_max = omp_get_num_threads();
  
  int ngpu = 0;
 #if USE_ACC
  ngpu = acc_get_num_devices(acc_device_nvidia);
  int nomp = ngpu;
 #else
  int nomp = nomp_max;
 #endif

  //printf("  nomp: %2i \n",nomp);

  int N = basis.size();
  int N2 = N*N;

  double atnod[natoms];
  get_atnod(natoms,basis,atnod);

  int qoh = quad_r_order; //refined grid region

  if (quad_order!=quad_r_order)
  {
    printf("\n WARNING: 3c pVp requires quad_orders to be equal \n");
    qoh = quad_r_order = quad_order;
  }

  int qos = quad_order*quad_order*quad_order;
  int qosh = qoh*qoh*qoh;
  if (natoms<3) nsplit = 1;
  int nsg = nsplit*nsplit*nsplit;
  int gs0 = nmu*nnu*nphi;
  int gs = gs0*qos; //2-atom grid size
  int gsh = (gs0-8)*qos+8*nsg*qosh; //3-atom grid size
  gs0 = nmu*nnu*nphi + 8*(nsg-1); //okay since qos==qoh
  int gs3 = 3*gsh;
  int gs6 = 6*gsh;

  //printf("   qos/h: %3i %3i \n",qos,qosh);
  if (prl>1) printf("   gs(h): %8i  %8i \n",gs,gsh);

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

  double* grid = new double[gs6];
  double* wt = new double[gsh];

  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);

  double gpumem = 1.*acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  double togb = 1./1024./1024./1024.;

 //this calculation not accurate yet
  int Nmax = 150;
  double mem0 = gs6+2*gsh;
  while (Nmax>0)
  {
    //double mem1 = 8.*(2.*Nmax*gs3 + 1.*mem0);
    double mem1 = 8.*(2.*Nmax*gs3 + 1.*mem0 + 1.*N2 + 3.*gs6 + 2.*gsh);
    if (mem1<gpumem)
    {
      if (prl>1) printf("    mem0: %5.3f mem1: %5.3f \n",mem0*togb,mem1*togb);
      break;
    }
    Nmax--;
  }
  //if (Nmax>4)
  //  Nmax -= 4;
  if (Nmax<=0) { printf("\n ERROR: couldn't calculate gpu memory requirements \n"); exit(-1); }
    
  vector<vector<int> > n2ip;
  int imaxN = get_imax_n2ip(Nmax,natoms,N,basis,n2ip);
  if (prl<1) printf("   imaxN: %2i \n",imaxN);

  double gpumem_gb = gpumem/1024./1024./1024.;
  printf("   gpu memory available: %6.3f GB \n",gpumem_gb);

 //intermediate storage
  iN = imaxN;
  double* valS1 = new double[iN * gs3];
  double* valS2 = new double[iN * gs3];

  double* valt = new double[gsh];
  double* pVpp = new double[N2];

 #if USE_ACC
  #pragma omp parallel for schedule(static) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = n;
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc enter data create(pVpp[0:N2])

    #pragma acc enter data create(grid[0:gs6],wt[0:gsh])
    #pragma acc enter data create(valS1[0:iN*gs3], valS2[0:iN*gs3])
    #pragma acc enter data create(valt[0:gsh])

    acc_assign(N2,pVpp,0.);
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  double gpumem_2 = 1.*acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  printf("   after alloc, gpu memory available: %6.3f GB \n",gpumem_2*togb); 

 //3c part of pVp integrals
 #pragma omp parallel for schedule(dynamic) num_threads(nomp)
  for (int m=0;m<natoms;m++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);

    //double Z1 = atnod[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];
    double coordn[9];
    coordn[0] = coordn[1] = coordn[2] = 0.;

   //three-atom case
    for (int n=m+1;n<natoms;n++)
    {
      //double Z2 = atnod[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;

      for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
      for (int sp2=0;sp2<n2ip[n].size()-1;sp2++)
      {
        int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
        int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2+1];

        double ztm1 = 0.; int lm1 = 0; double ztm2 = 0.; int lm2 = 0;
        get_ztm_lm(s1,s2,basis,ztm1,lm1);
        get_ztm_lm(s3,s4,basis,ztm2,lm2);
        //printf("    pVp_3c ztm/lm: %8.5f %i - %8.5f %i \n",ztm1,lm1,ztm2,lm2);

        for (int p=0;p<natoms;p++)
        if (p!=m && p!=n)
        {
          //double Z3 = (double)atno[p];
          double Z3 = atnod[p];
          double A3 = coords[3*p+0]; double B3 = coords[3*p+1]; double C3 = coords[3*p+2];
          double A13 = A3-A1; double B13 = B3-B1; double C13 = C3-C1;
          coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;
          coordn[6] = A13; coordn[7] = B13; coordn[8] = C13;

         //////////////////////////////////////////////////////////////////////////////////////////////////////
         // reduction over E-n attraction
         //  s1 on atom 1, s on atom 2, v on atom 3

        #if SWITCH_23
         //   but atoms 2 and 3 are switched in grid generation
          coordn[3] = A13; coordn[4] = B13; coordn[5] = C13;
          coordn[6] = A12; coordn[7] = B12; coordn[8] = C12;
        #endif

          //generate_ps_quad_grid(3,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
          generate_ps_quad_grid_3c_refine(ztm1,ztm2,nsplit,3,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
          add_r1_to_grid(gsh,grid,0.,0.,0.);

          init_s12v3(0,s1,s2,s3,s4,0,0,iN,0,gs3,  valS1,valS2,NULL,NULL);
          eval_p12  (s1,s2,s3,s4,gsh,grid,basis,valS1,valS2,A12,B12,C12,A13,B13,C13);

          recenter_grid_zero(gsh,grid,-A13,-B13,-C13);
          reduce_3cenp(Z3,s1,s2,s3,s4,N,iN,gsh,grid,valS1,valS2,valt,wt,pVpp);
         //////////////////////////////////////////////////////////////////////////////////////////////////////

        } //loop p over third atom

      }//loop sp12

    } //loop n over second atom

  } //loop m over natoms
  acc_set_device_num(0,acc_device_nvidia);


  double norm1[N];
  for (int i=0;i<N;i++)
    norm1[i] = basis[i][4];
  #pragma acc enter data copyin(norm1[0:N])

  if (nomp>1)
  {
    double pVpt[N2];
    for (int j=0;j<N2;j++) pVpt[j] = 0.;

    for (int n=0;n<nomp;n++)
    {
      int tid = n;
      acc_set_device_num(tid,acc_device_nvidia);

      #pragma acc update self(pVpp[0:N2])

      for (int j=0;j<N2;j++)
        pVpt[j] += pVpp[j];
    }
    acc_set_device_num(0,acc_device_nvidia);

    for (int i=0;i<N;i++)
    for (int j=0;j<N;j++)
      pVpt[i*N+j] *= norm1[i]*norm1[j];

    for (int i=0;i<N;i++)
    for (int j=0;j<i;j++)
      pVpt[i*N+j] = pVpt[j*N+i];

    //add to prior pVp term
    for (int j=0;j<N2;j++)
      pVp[j] += pVpt[j];
  }
  else
  {
   #pragma acc parallel loop collapse(2) present(pVpp[0:N2],norm1[0:N])
    for (int i=0;i<N;i++)
    for (int j=0;j<N;j++)
      pVpp[i*N+j] *= norm1[i]*norm1[j];

    #pragma acc parallel loop independent present(pVpp[0:N2])
    for (int i=0;i<N;i++)
    #pragma acc loop independent
    for (int j=0;j<i;j++)
      pVpp[i*N+j] = pVpp[j*N+i];

    #pragma acc update self(pVpp[0:N2])

    for (int i=0;i<N2;i++)
      pVp[i] += pVpp[i];
  }

  const double lt = 1.e-15;
  for (int j=0;j<N2;j++)
  if (fabs(pVp[j])<lt)
    pVp[j] = 0.;

  if (prl>1)
  {
    printf("\n pVp (3c): \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %12.6f",pVp[i*N+j]);
      printf("\n");
    }
  }

 #if USE_ACC
 #pragma omp parallel for schedule(static) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = n;
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc exit data delete(pVpp[0:N2])
    #pragma acc exit data delete(norm1[0:N])
    #pragma acc exit data delete(grid[0:gs6],wt[0:gsh])
    #pragma acc exit data delete(valS1[0:iN*gs3], valS2[0:iN*gs3])
    #pragma acc exit data delete(valt[0:gsh])
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  //printf(" done with dealloc in 3c integrals \n"); fflush(stdout);
  //auto_crash();

  delete [] n2i;


  delete [] valS1; delete [] valS2;
  delete [] valt;
  delete [] pVpp;

  delete [] grid;
  delete [] wt;

  return;
}

void compute_3c_ps(bool do_overlap, bool do_yukawa, double gamma, int nbatch, int natoms, int* atno, double* coords, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, int quad_order, int quad_r_order, int nsplit, int nmu, int nnu, int nphi, double* En, double* C, int prl)
{
  if (prl>-1) { if (do_yukawa) printf("  beginning compute_3c_ps (Yukawa. gamma: %5.3f) \n",gamma); else printf("  beginning compute_3c_ps \n"); }
  if (do_overlap && prl>1) { printf("\n WARNING: testing do_overlap in compute_3c_ps \n"); }

  int nomp_max = 1;
 #pragma omp parallel
  nomp_max = omp_get_num_threads();
  
  int ngpu = 0;
 #if USE_ACC
  ngpu = acc_get_num_devices(acc_device_nvidia);
  int nomp = ngpu;
 #else
  int nomp = nomp_max;
 #endif

  double cfn_en = 1.0;

  bool dy = do_yukawa;
  bool dol = do_overlap;
  if (En!=NULL && dy)
  { printf("\n ERROR: cannot run En with Yukawa in compute_3c_ps \n"); exit(-1); }

  int N = basis.size();
  int N2 = N*N;
  int Naux = basis_aux.size();
  int N2a = N2*Naux;

  if (N<1) { printf(" ERROR: cannot compute 3c integrals, no primary basis functions \n"); exit(-1); }
  if (Naux<1) { printf(" ERROR: cannot compute 3c integrals, no RI basis functions \n"); exit(-1); }

  double atnod[natoms];
  get_atnod(natoms,basis,atnod);

  int qoh = quad_r_order; //refined grid region

  int qos = quad_order*quad_order*quad_order;
  int qosh = qoh*qoh*qoh;
  if (natoms<3) nsplit = 1;
  int nsg = nsplit*nsplit*nsplit;

  if (0)
  if (nbatch>1 && natoms>2)
  {
    printf("\n WARNING: cannot use nbatch with natoms>2 yet \n");
    nbatch = 1;
  }

 //batching both grids
  int gs = (nmu*nnu*nphi)*qos/nbatch; //2-atom grid size
  int gsh = ((nmu*nnu*nphi-8)*qos+8*nsg*qosh)/nbatch; //3-atom grid size
  int gs6 = 6*gsh;

  if (prl>1) printf("   gs(h): %8i  %8i  nb: %2i \n",gs,gsh,nbatch);

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

  double* grid = new double[gs6];
  double* wt = new double[gsh];

  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);

  double gpumem = 1.*acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  double togb = 1./1024./1024./1024.;

 //need to improve this estimate
  int Nmax = 100;
  double mem0 = gsh*iN*2. + gsh*7. + 1.*nmu*nnu*nphi + 1.*gs6 + 1.*N2 + 2.*N2a;
  while (Nmax>0)
  {
    double mem1 = 8.*(2.*gsh*iN + mem0 + 3.*gs6 + 2.*gsh + 1.*Nmax*gsh);
    if (mem1<gpumem)
    {
      if (prl>1) printf("    mem0: %6.1f mem1: %6.1f \n",mem0*togb,mem1*togb);
      break;
    }
    Nmax--;
  }
  //if (Nmax>5)
  //  Nmax -= 5;
  if (Nmax<=0) { printf("\n ERROR: couldn't calculate gpu memory requirements (mem0: %6.1f) \n",mem0*togb); exit(-1); }

  vector<vector<int> > n2aip;
  int iNa = get_imax_n2ip(Nmax,natoms,Naux,basis_aux,n2aip);
  if (prl>1) printf("   iNa: %2i \n",iNa);

  int* na2i = new int[natoms]; //needed for copy_symm
  get_imax_n2i(natoms,Naux,basis_aux,na2i);

  double gsxvalsv = 8.*(gsh*iN*2. + gsh*iNa + gsh*7. + 1.*nmu*nnu*nphi); //vals+grid/wt+gridm
  double gsxvalsv_gb = gsxvalsv*togb;
  printf("   estimated memory needed (3c): %6.3f GB \n",gsxvalsv_gb);

  double gpumem_gb = gpumem/1024./1024./1024.;
  printf("   gpu memory available: %6.3f GB \n",gpumem_gb);

  if (gsxvalsv>gpumem) { printf("\n WARNING: probably not enough memory to do 3c integrals \n"); }

 //intermediate storage
  double* valS1 = new double[iN * gsh];
  double* valS2 = new double[iN * gsh];
  double* valV3 = new double[iNa * gsh];

  double* valt = new double[gsh];
  double* Enp = new double[N2];
  double* Cp = new double[N2a];

 #if USE_ACC
  #pragma omp parallel for schedule(static) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = n;
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc enter data create(Enp[0:N2],C[0:N2a],Cp[0:N2a])

    #pragma acc enter data create(grid[0:gs6],wt[0:gsh])
    #pragma acc enter data create(valS1[0:iN * gsh], valS2[0:iN * gsh], valV3[0:iNa * gsh])
    #pragma acc enter data create(valt[0:gsh])

    acc_assign(gs6,grid,0.);
    acc_assign(gsh,grid,0.);
    acc_assign(N2,Enp,0.);
    acc_assign(N2a,Cp,0.);
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  double gpumem_2 = 1.*acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  printf("   after alloc, gpu memory available: %6.3f GB \n",gpumem_2*togb); 

 //Coulomb 3c integrals and 3c part of En integrals
 #pragma omp parallel for schedule(dynamic) num_threads(nomp)
  for (int m=0;m<natoms;m++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);

    double Z1 = atnod[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];
    double coordn[9];
    coordn[0] = coordn[1] = coordn[2] = 0.;

   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];
    int s3 = s1; int s4 = s2;

    for (int wb=0;wb<nbatch;wb++)
    {
      generate_ps_quad_grid(1.,wb,nbatch,Z1,1,coordn,quad_order,quad_order,nmu,nnu,nphi,grid,wt);
      add_r1_to_grid(gs,grid,0.,0.,0.);

      for (int sp=0;sp<n2aip[m].size()-1;sp++)
      {
        int s5 = n2aip[m][sp]; int s6 = n2aip[m][sp+1];

       //all basis on one atom
        init_s12v3(dy,s1,s2,s3,s4,s5,s6,iN,iNa,gs,valS1,valS2,valV3,wt);
        eval_s12v3(dol,dy,gamma,s1,s2,s3,s4,s5,s6,gs,grid,basis,basis_aux,valS1,valS2,valV3);

        //printf("    m: %i   s12: %2i %2i s56: %2i %2i \n",m,s1,s2,s5,s6);
        reduce_3c1b(s5,s6,s1,s2,gs,valV3,valS1,valS2,N,Naux,iN,iNa,Cp);

      } //loop sp

    } //loop wb over batches

   //two-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      //double Z2 = atnod[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;
      coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;

      s3 = 0; if (n>0) s3 = n2i[n-1]; s4 = n2i[n];

      //printf("     mn: %i %i   s12: %2i %2i s34: %2i %2i s56: %2i %2i \n",m,n,s1,s2,s3,s4,s5,s6);
      for (int wb=0;wb<nbatch;wb++)
      {
        generate_ps_quad_grid(1.,wb,nbatch,0.,2,coordn,quad_order,quad_order,nmu,nnu,nphi,grid,wt);
        add_r1_to_grid(gs,grid,0.,0.,0.);

        for (int sp=0;sp<n2aip[n].size()-1;sp++)
        {
          int s5 = n2aip[n][sp]; int s6 = n2aip[n][sp+1];
          //printf(" n: %i  s5/6: %3i %3i \n",n,s5,s6);

         //s1 on atom 1, s2v3 on atom 2
          init_s12v3  (dy,s1,s2,s3,s4,s5,s6,iN,iNa,gs,              valS1,valS2,valV3,wt);
          eval_s12v3_2(dol,dy,gamma,s1,s2,s3,s4,s5,s6,gs,grid,basis,basis_aux,valS1,valS2,valV3,1,A12,B12,C12,0.,0.,0.);

          reduce_3c1b(s5,s6,s1,s2,s3,s4,gs,valV3,valS1,valS2,N,Naux,iN,iNa,Cp);

         //s12 on atom 1, v3 on atom 2
          int s3b = s1; int s4b = s2;
          //printf("     mn: %i %i   s12: %2i %2i s34: %2i %2i s56: %2i %2i \n",m,n,s1,s2,s3,s4,s5,s6);

          init_s12v3  (dy,s1,s2,s3b,s4b,s5,s6,iN,iNa,gs,              valS1,valS2,valV3,wt);
          eval_s12v3_2(dol,dy,gamma,s1,s2,s3b,s4b,s5,s6,gs,grid,basis,basis_aux,valS1,valS2,valV3,2,A12,B12,C12,0.,0.,0.);

          reduce_3c1b(s5,s6,s1,s2,s3b,s4b,gs,valV3,valS1,valS2,N,Naux,iN,iNa,Cp);
        }
      } //loop wb over nbatch

    } //loop n over second atom

   //three-atom case
    for (int n=m+1;n<natoms;n++)
    {
      //double Z2 = atnod[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;

      for (int p=0;p<natoms;p++)
      if (p!=m && p!=n)
      {
        s1 = 0; if (m>0) s1 = n2i[m-1]; s2 = n2i[m];
        s3 = 0; if (n>0) s3 = n2i[n-1]; s4 = n2i[n];

        double Z3 = atnod[p];
        double A3 = coords[3*p+0]; double B3 = coords[3*p+1]; double C3 = coords[3*p+2];
        double A13 = A3-A1; double B13 = B3-B1; double C13 = C3-C1;
        coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;
        coordn[6] = A13; coordn[7] = B13; coordn[8] = C13;

        for (int sp=0;sp<n2aip[p].size()-1;sp++)
        {
          int s5 = n2aip[p][sp]; int s6 = n2aip[p][sp+1];
          //printf(" p: %i  s5/6: %3i %3i \n",p,s5,s6);

          //if (natoms>3) printf("     mnp: %i %i %i   s12: %3i %3i  s34: %3i %3i  s56: %3i %3i \n",m,n,p,s1,s2,s3,s4,s5,s6);

          double ztm1 = 0.; int lm1 = 0; double ztm2 = 0.; int lm2 = 0;
          get_ztm_lm(s1,s2,basis,ztm1,lm1);
          get_ztm_lm(s3,s4,basis,ztm2,lm2);
          //printf("    3c ztm/lm: %8.5f %i - %8.5f %i \n",ztm1,lm1,ztm2,lm2);

          for (int wb=0;wb<nbatch;wb++)
          {
            generate_ps_quad_grid_3c_refine(wb,nbatch,ztm1,ztm2,nsplit,3,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
            add_r1_to_grid(gsh,grid,0.,0.,0.);

           //s1 on atom 1, s2 on atom 2, v3 on atom 3
            init_s12v3  (dy,s1,s2,s3,s4,s5,s6,iN,iNa,gsh,              valS1,valS2,valV3,wt);
            eval_s12v3_2(dol,dy,gamma,s1,s2,s3,s4,s5,s6,gsh,grid,basis,basis_aux,valS1,valS2,valV3,3,A12,B12,C12,A13,B13,C13);

            reduce_3c1b(s5,s6,s1,s2,s3,s4,gsh,valV3,valS1,valS2,N,Naux,iN,iNa,Cp);
          }
        }

        if (!dol && !dy)
        {
         //////////////////////////////////////////////////////////////////////////////////////////////////////
         // reduction over E-n attraction
         //  s1 on atom 1, s on atom 2, v on atom 3

        #if SWITCH_23
         //   but atoms 2 and 3 are switched in grid generation
          coordn[3] = A13; coordn[4] = B13; coordn[5] = C13;
          coordn[6] = A12; coordn[7] = B12; coordn[8] = C12;
        #endif

          //printf("  mnp: %i %i %i  xyz: %8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f \n",m,n,p,coordn[0],coordn[1],coordn[2],coordn[3],coordn[4],coordn[5],coordn[6],coordn[7],coordn[8]);

          double ztm1 = 0.; int lm1 = 0; double ztm2 = 0.; int lm2 = 0;
          get_ztm_lm(s1,s2,basis,ztm1,lm1);
          get_ztm_lm(s3,s4,basis,ztm2,lm2);
          //printf("    3c ztm/lm: %8.5f %i - %8.5f %i \n",ztm1,lm1,ztm2,lm2);

          for (int wb=0;wb<nbatch;wb++)
          {
            //generate_ps_quad_grid(3,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
            generate_ps_quad_grid_3c_refine(cfn_en,wb,nbatch,ztm1,ztm2,nsplit,3,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
            add_r1_to_grid(gsh,grid,0.,0.,0.);

            init_s12v3  (0,s1,s2,s3,s4,0,0,iN,iNa,gsh,              valS1,valS2,NULL,wt);
            eval_s12v3_2(0,0,0.,s1,s2,s3,s4,0,0,gsh,grid,basis,basis_aux,valS1,valS2,NULL,3,A12,B12,C12,A13,B13,C13);

            recenter_grid_zero(gsh,grid,-A13,-B13,-C13);
            reduce_3cen(Z3,s1,s2,s3,s4,N,iN,gsh,grid,valS1,valS2,valt,wt,Enp);
          }
         //////////////////////////////////////////////////////////////////////////////////////////////////////
        }

      } //loop p over third atom

    } //loop n over second atom

  } //loop m over natoms
  acc_set_device_num(0,acc_device_nvidia);


  double norm1[Naux];
  double norm2[N];
  if (!dol)
  for (int i=0;i<Naux;i++)
    norm1[i] = norm_sv(basis_aux[i][0],basis_aux[i][1],basis_aux[i][2],basis_aux[i][3]);
  else
  for (int i=0;i<Naux;i++)
    norm1[i] = norm(basis_aux[i][0],basis_aux[i][1],basis_aux[i][2],basis_aux[i][3]);
  for (int i=0;i<N;i++)
    norm2[i] = basis[i][4];
  #pragma acc enter data copyin(norm1[0:Naux],norm2[0:N])

  if (nomp>1)
  {
    double Ent[N2];
    for (int j=0;j<N2;j++) Ent[j] = 0.;
    for (int j=0;j<N2a;j++) C[j] = 0.;

    for (int n=0;n<nomp;n++)
    {
      int tid = n;
      acc_set_device_num(tid,acc_device_nvidia);

      #pragma acc update self(Cp[0:N2a],Enp[0:N2])

      for (int j=0;j<N2;j++)
        Ent[j] += Enp[j];

      for (int j=0;j<N2a;j++)
        C[j] += Cp[j];
    }
    acc_set_device_num(0,acc_device_nvidia);
    #pragma acc update device(C[0:N2a])

    if (!dol && !dy)
    {
      for (int i=0;i<N;i++)
      for (int j=0;j<N;j++)
        Ent[i*N+j] *= norm2[i]*norm2[j];

      for (int i=0;i<N;i++)
      for (int j=0;j<i;j++)
        Ent[i*N+j] = Ent[j*N+i];

      //add to prior STEn term
      for (int j=0;j<N2;j++)
        En[j] += Ent[j];
    }
  }
  else
  {
    if (!dol && !dy)
    {
     #pragma acc parallel loop independent present(Enp[0:N2],norm2[0:N])
      for (int i=0;i<N;i++)
      #pragma acc loop
      for (int j=0;j<N;j++)
        Enp[i*N+j] *= norm2[i]*norm2[j];

      #pragma acc parallel loop independent present(Enp[0:N2])
      for (int i=0;i<N;i++)
      #pragma acc loop independent
      for (int j=0;j<i;j++)
        Enp[i*N+j] = Enp[j*N+i];

      #pragma acc update self(Enp[0:N2])

      for (int i=0;i<N2;i++)
        En[i] += Enp[i];
    }

   #pragma acc parallel loop present(C[0:N2a],Cp[0:N2a])
    for (int j=0;j<N2a;j++)
      C[j] = Cp[j];
  }

 #pragma acc parallel loop collapse(3) present(C[0:N2a],norm1[0:Naux],norm2[0:N])
  for (int i=0;i<Naux;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
  {
    double n123 = norm1[i]*norm2[j]*norm2[k];
    C[i*N2+j*N+k] *= n123;
  }
  #pragma acc update self(C[0:N2a])

  copy_symm_3c_ps(natoms,N,Naux,n2i,na2i,C);

  //if (!do_overlap)
    transpose_C(Naux,N,C);

  if (prl>2 || (prl>1 && N<10))
  {
    int nna = N*Naux;
    printf("\n C: \n");
    for (int i=0;i<Naux;i++)
    {
      //printf(" i: %i \n",i);
      for (int j=0;j<N;j++)
      for (int k=0;k<N;k++)
        printf("  %12.6f",C[j*nna+k*Naux+i]);
      printf("\n");
    }
  }

 #if USE_ACC
 #pragma omp parallel for schedule(static) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = n;
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc exit data delete(Enp[0:N2],C[0:N2a],Cp[0:N2a])
    #pragma acc exit data delete(norm1[0:Naux],norm2[0:N])
    #pragma acc exit data delete(grid[0:gs6],wt[0:gsh])
    #pragma acc exit data delete(valS1[0:iN * gsh], valS2[0:iN * gsh], valV3[0:iNa * gsh])
    #pragma acc exit data delete(valt[0:gsh])
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  //printf(" done with dealloc in 3c integrals \n"); fflush(stdout);
  //auto_crash();

  delete [] n2i;
  delete [] na2i;


  delete [] valS1; delete [] valS2; delete [] valV3;
  delete [] valt;
  delete [] Enp;
  delete [] Cp;

  delete [] grid;
  delete [] wt;

  return;
}

void init_s14(int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8, int iN, int gs, double** val1, double** val2, double** val3, double** val4, double* wt)
{
  #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gs])
  for (int ii1=0;ii1<s2-s1;ii1++)
  {
    for (int j=0;j<gs;j++)
      val1[ii1][j] = 1.;
  }

  #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs])
  for (int ii2=0;ii2<s4-s3;ii2++)
  {
    for (int j=0;j<gs;j++)
      val2[ii2][j] = 1.;
  }

  #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs])
  for (int ii3=0;ii3<s6-s5;ii3++)
  {
    for (int j=0;j<gs;j++)
      val3[ii3][j] = 1.;
  }

  #pragma acc parallel loop collapse(2) present(val4[0:iN][0:gs],wt[0:gs])
  for (int ii4=0;ii4<s8-s7;ii4++)
  {
    for (int j=0;j<gs;j++)
      val4[ii4][j] = wt[j];
  }

  return;
}

void eval_s14(int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8, int gs, double* grid, vector<vector<double> >& basis, double** val1, double** val2, double** val3, double** val4, int type, double A12, double B12, double C12)
{
  for (int i1=s1;i1<s2;i1++)
  {
    int ii1 = i1-s1;

    vector<double> basis1 = basis[i1];
    int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

    eval_shd(ii1,gs,grid,val1[ii1],n1,l1,m1,zeta1);
  }

  if (type==1)
    recenter_grid_zero(gs,grid,-A12,-B12,-C12);

  for (int i2=s3;i2<s4;i2++)
  {
    int ii2 = i2-s3;

    vector<double> basis2 = basis[i2];
    int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

    eval_shd(ii2,gs,grid,val2[ii2],n2,l2,m2,zeta2);
  }

  if (type==2)
    recenter_grid_zero(gs,grid,-A12,-B12,-C12);

  for (int i3=s5;i3<s6;i3++)
  {
    int ii3 = i3-s5;

    vector<double> basis3 = basis[i3];
    int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

    eval_shd(ii3,gs,grid,val3[ii3],n3,l3,m3,zeta3);
  }

  if (type==3)
    recenter_grid_zero(gs,grid,-A12,-B12,-C12);

  for (int i4=s7;i4<s8;i4++)
  {
    int ii4 = i4-s7;

    vector<double> basis4 = basis[i4];
    int n4 = basis4[0]; int l4 = basis4[1]; int m4 = basis4[2]; double zeta4 = basis4[3];

    eval_shd(ii4,gs,grid,val4[ii4],n4,l4,m4,zeta4);
  }

  if (type>0)
    recenter_grid_zero(gs,grid,A12,B12,C12);

  return;
}

void eval_s14b(int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8, int gs, double* grid, vector<vector<double> >& basis, double** val1, double** val2, double** val3, double** val4, int type, double A12, double B12, double C12, double A13, double B13, double C13)
{
  //s12 on atom 1, s3 on atom 2, s4 on atom 3 (type==1)
  //s1 on atom 1, s23 on atom 2, s4 on atom 3 (type==2)
  //s1 on atom 1, s2 on atom 2, s34 on atom 3 (type==3)
  //s1 on atom 1, s2 on atom 3, s34 on atom 2 (type==4)

  for (int i1=s1;i1<s2;i1++)
  {
    int ii1 = i1-s1;

    vector<double> basis1 = basis[i1];
    int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

    eval_shd(ii1,gs,grid,val1[ii1],n1,l1,m1,zeta1);
  }

  if (type==2 || type==3)
    recenter_grid_zero(gs,grid,-A12,-B12,-C12);
  if (type==4)
    recenter_grid_zero(gs,grid,-A13,-B13,-C13);

  for (int i2=s3;i2<s4;i2++)
  {
    int ii2 = i2-s3;

    vector<double> basis2 = basis[i2];
    int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

    eval_shd(ii2,gs,grid,val2[ii2],n2,l2,m2,zeta2);
  }

  if (type==1)
    recenter_grid_zero(gs,grid,-A12,-B12,-C12);
  if (type==3)
    recenter_grid_zero(gs,grid,A12-A13,B12-B13,C12-C13);
  if (type==4)
    recenter_grid_zero(gs,grid,A13-A12,B13-B12,C13-C12);

  for (int i3=s5;i3<s6;i3++)
  {
    int ii3 = i3-s5;

    vector<double> basis3 = basis[i3];
    int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

    eval_shd(ii3,gs,grid,val3[ii3],n3,l3,m3,zeta3);
  }

  if (type==1 || type==2)
    recenter_grid_zero(gs,grid,A12-A13,B12-B13,C12-C13);

  for (int i4=s7;i4<s8;i4++)
  {
    int ii4 = i4-s7;

    vector<double> basis4 = basis[i4];
    int n4 = basis4[0]; int l4 = basis4[1]; int m4 = basis4[2]; double zeta4 = basis4[3];

    eval_shd(ii4,gs,grid,val4[ii4],n4,l4,m4,zeta4);
  }

  if (type==1 || type==2 || type==3)
    recenter_grid_zero(gs,grid,A13,B13,C13);
  if (type==4)
    recenter_grid_zero(gs,grid,A12,B12,C12);

  return;
}

void eval_s14c(int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8, int gs, double* grid, vector<vector<double> >& basis, double** val1, double** val2, double** val3, double** val4, double A12, double B12, double C12, double A13, double B13, double C13, double A14, double B14, double C14)
{
  for (int i1=s1;i1<s2;i1++)
  {
    int ii1 = i1-s1;

    vector<double> basis1 = basis[i1];
    int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

    eval_shd(ii1,gs,grid,val1[ii1],n1,l1,m1,zeta1);
  }

  recenter_grid_zero(gs,grid,-A12,-B12,-C12);

  for (int i2=s3;i2<s4;i2++)
  {
    int ii2 = i2-s3;

    vector<double> basis2 = basis[i2];
    int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

    eval_shd(ii2,gs,grid,val2[ii2],n2,l2,m2,zeta2);
  }

  recenter_grid_zero(gs,grid,A12-A13,B12-B13,C12-C13);

  for (int i3=s5;i3<s6;i3++)
  {
    int ii3 = i3-s5;

    vector<double> basis3 = basis[i3];
    int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

    eval_shd(ii3,gs,grid,val3[ii3],n3,l3,m3,zeta3);
  }

  recenter_grid_zero(gs,grid,A13-A14,B13-B14,C13-C14);

  for (int i4=s7;i4<s8;i4++)
  {
    int ii4 = i4-s7;

    vector<double> basis4 = basis[i4];
    int n4 = basis4[0]; int l4 = basis4[1]; int m4 = basis4[2]; double zeta4 = basis4[3];

    eval_shd(ii4,gs,grid,val4[ii4],n4,l4,m4,zeta4);
  }

  recenter_grid_zero(gs,grid,A14,B14,C14);

  return;
}

void reduce_4c_ol1(int ii1, int ii2, int s5, int s6, int s7, int s8, int gs, double* valm, double* valn, double** val3, double** val4, int iN, int M, int M2, int M3, int M4, double* tmp)
{
 //integrate second pair of ftns

 #pragma acc parallel loop collapse(2) present(valm[0:gs],valn[0:gs],val3[0:iN][0:gs],val4[0:iN][0:gs],tmp[0:M4])
  for (int i3=s5;i3<s6;i3++)
  for (int i4=s7;i4<s8;i4++)
  {
    int ii3 = i3-s5;
    int ii4 = i4-s7;

    double* valp = val3[ii3];
    double* valq = val4[ii4];

    double v1 = 0.;
   #pragma acc loop reduction(+:v1)
    for (int j=0;j<gs;j++)
      v1 += valm[j]*valn[j]*valp[j]*valq[j];

    tmp[ii1*M3+ii2*M2+ii3*M+ii4] = v1;
  }

  return;
}

//worth exploring acc optimization
void reduce_4c_ol1(int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8, int gs, double** val1, double** val2, double** val3, double** val4, int N, int iN, double* tmp, double* olp)
{
  int N2 = N*N;
  int N3 = N2*N;

 //save GPU memory for N4 array
  int s12 = s2-s1;
  int s34 = s4-s3;
  int s56 = s6-s5;
  int s78 = s8-s7;
  int M = s78;
  int M2 = s56*M;
  int M3 = s34*M2;
  int M4 = s12*M3;

 //integrals have 24-fold symmetry.
 // currently using only 2-fold
  bool symm_one = 0;
  bool symm_two = 0;
  if (s1==s3)
    symm_one = 1;
  if (s5==s7)
    symm_two = 1;

  //if (symm_one)
  //  printf("    using symm_one \n");

  if (symm_one)
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<=i1;i2++)
  {
    int ii1 = i1-s1;
    int ii2 = i2-s3;

    double* valm = val1[ii1];
    double* valn = val2[ii2];

    reduce_4c_ol1(ii1,ii2,s5,s6,s7,s8,gs,valm,valn,val3,val4,iN,M,M2,M3,M4,tmp);

   #pragma acc serial present(tmp[0:M4])
    for (int i3=s5;i3<s6;i3++)
    for (int i4=s7;i4<s8;i4++)
    {
      int ii3 = i3-s5;
      int ii4 = i4-s7;

      double v1 = tmp[ii1*M3+ii2*M2+ii3*M+ii4];
      tmp[ii2*M3+ii1*M2+ii3*M+ii4] = v1;
    }
  }

  if (!symm_one)
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  {
    int ii1 = i1-s1;
    int ii2 = i2-s3;

    double* valm = val1[ii1];
    double* valn = val2[ii2];

    reduce_4c_ol1(ii1,ii2,s5,s6,s7,s8,gs,valm,valn,val3,val4,iN,M,M2,M3,M4,tmp);
  }

  #pragma acc update self(tmp[0:M4])
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  for (int i3=s5;i3<s6;i3++)
  for (int i4=s7;i4<s8;i4++)
    olp[i1*N3+i2*N2+i3*N+i4] = tmp[(i1-s1)*M3+(i2-s3)*M2+(i3-s5)*M+(i4-s7)];

  return;
}


void compute_4c_ol_ps_fast(int natoms, int* atno, double* coords, vector<vector<double> > &basis, int quad_order, int quad_r_order, int nmu, int nnu, int nphi, double* ol, int prl)
{
  if (prl>-1) printf("  beginning compute_4c_ol_ps_fast \n");

  int nomp_max = 1;
 #pragma omp parallel
  nomp_max = omp_get_num_threads();
  
  int ngpu = 0;
 #if USE_ACC
  ngpu = acc_get_num_devices(acc_device_nvidia);
  int nomp = ngpu;
 #else
  int nomp = nomp_max;
 #endif

  int N = basis.size();
  int N2 = N*N;
  int N3 = N2*N;
  int N4 = N3*N;

//  quad_r_order = quad_order;
  int qos = quad_order*quad_order*quad_order;
  int qoh = quad_r_order;
  int qosh = qoh*qoh*qoh;
  int gs = nmu*nnu*nphi*qos;
  int gsh = (nmu*nnu*nphi-8)*qos+8*qosh;
  int gshh = (nmu*nnu*nphi-16)*qos+16*qosh;
  int gs6 = 6*gshh;

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

  double* grid = new double[gs6];
  double* wt = new double[gshh];

  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);

 //intermediate storage
  double** valS1 = new double*[iN]; for (int i=0;i<iN;i++) valS1[i] = new double[gshh];
  double** valS2 = new double*[iN]; for (int i=0;i<iN;i++) valS2[i] = new double[gshh];
  double** valS3 = new double*[iN]; for (int i=0;i<iN;i++) valS3[i] = new double[gshh];
  double** valS4 = new double*[iN]; for (int i=0;i<iN;i++) valS4[i] = new double[gshh];
  double* valt = new double[gshh];
 //CPMZ reconfiguring this
  int M4 = iN*iN*iN*iN;
  double* tmp = new double[M4];
  double* olp = new double[N4]();

 #if USE_ACC
 #pragma omp parallel for schedule(dynamic) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = n;
    acc_set_device_num(tid,acc_device_nvidia);

    //#pragma acc enter data create(olp[0:N4])
    #pragma acc enter data create(tmp[0:M4])
    #pragma acc enter data create(grid[0:gs6],wt[0:gshh])
    #pragma acc enter data create(valS1[0:iN][0:gshh],valS2[0:iN][0:gshh],valS3[0:iN][0:gshh],valS4[0:iN][0:gshh])
    #pragma acc enter data create(valt[0:gshh])

    acc_assign(M4,tmp,0.);
    //acc_assign(N4,olp,0.);
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

 #pragma omp parallel for schedule(dynamic) num_threads(nomp)
  for (int m=0;m<natoms;m++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);

   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    double Z1 = (double)atno[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];
    double coordn[12];
    coordn[0] = coordn[1] = coordn[2] = 0.;

    generate_ps_quad_grid(Z1,1,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
    add_r1_to_grid(gs,grid,0.,0.,0.);

   //all basis on one atom
    init_s14(s1,s2,s1,s2,s1,s2,s1,s2,iN,gs,valS1,valS2,valS3,valS4,wt);
    eval_s14(s1,s2,s1,s2,s1,s2,s1,s2,gs,grid,basis,valS1,valS2,valS3,valS4,0,0.,0.,0.);

    //printf("    m: %i   s12: %2i %2i s56: %2i %2i \n",m,s1,s2,s5,s6);
    reduce_4c_ol1(s1,s2,s1,s2,s1,s2,s1,s2,gs,valS1,valS2,valS3,valS4,N,iN,tmp,olp);

   //two-atom ints
    for (int n=m+1;n<natoms;n++)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      double Z2 = (double)atno[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;
      coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;

      //printf("     mn: %i %i   s12: %2i %2i s34: %2i %2i \n",m,n,s1,s2,s3,s4);
      generate_ps_quad_grid(0.,2,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
      add_r1_to_grid(gs,grid,0.,0.,0.);

     //s1 on atom 1, s234 on atom 2
      init_s14(s1,s2,s3,s4,s3,s4,s3,s4,iN,gs,valS1,valS2,valS3,valS4,wt);
      eval_s14(s1,s2,s3,s4,s3,s4,s3,s4,gs,grid,basis,valS1,valS2,valS3,valS4,1,A12,B12,C12);

      reduce_4c_ol1(s1,s2,s3,s4,s3,s4,s3,s4,gs,valS1,valS2,valS3,valS4,N,iN,tmp,olp);

     //s12 on atom 1, s34 on atom 2
      init_s14(0,0,s1,s2,0,0,0,0,iN,gs,valS1,valS2,valS3,valS4,wt);
      eval_s14(0,0,s1,s2,0,0,0,0,gs,grid,basis,valS1,valS2,valS3,valS4,0,A12,B12,C12);

      reduce_4c_ol1(s1,s2,s1,s2,s3,s4,s3,s4,gs,valS1,valS2,valS3,valS4,N,iN,tmp,olp);

     //s123 on atom 1, s4 on atom 2
      init_s14(0,0,0,0,s1,s2,0,0,iN,gs,valS1,valS2,valS3,valS4,wt);
      eval_s14(0,0,0,0,s1,s2,0,0,gs,grid,basis,valS1,valS2,valS3,valS4,0,A12,B12,C12);

      reduce_4c_ol1(s1,s2,s1,s2,s1,s2,s3,s4,gs,valS1,valS2,valS3,valS4,N,iN,tmp,olp);

    } //loop n over second atom

   //three-atom case
    for (int n=m+1;n<natoms;n++)
    {
      double Z2 = (double)atno[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;
      coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;

      for (int p=n+1;p<natoms;p++)
      if (p!=m)
      {
        int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];
        int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];
        int s5 = 0; if (p>0) s5 = n2i[p-1]; int s6 = n2i[p];

        double Z3 = (double)atno[p];
        double A3 = coords[3*p+0]; double B3 = coords[3*p+1]; double C3 = coords[3*p+2];
        double A13 = A3-A1; double B13 = B3-B1; double C13 = C3-C1;
        coordn[6] = A13; coordn[7] = B13; coordn[8] = C13;

        //printf("  mnp: %i %i %i   s12: %2i %2i  s34: %2i %2i  s56: %2i %2i \n",m,n,p,s1,s2,s3,s4,s5,s6);

        double ztm = 0.; int lm = 0;
        get_ztm_lm(s5,s6,basis,ztm,lm);

        generate_ps_quad_grid(0.,3,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
        add_r1_to_grid(gsh,grid,0.,0.,0.);

       //need all unique combinations, as each atom is unique

       //s12 on atom 1, s3 on atom 2, s4 on atom 3 (type==1)
        init_s14(s1,s2,s1,s2,s3,s4,s5,s6,iN,gsh,valS1,valS2,valS3,valS4,wt);
        eval_s14b(s1,s2,s1,s2,s3,s4,s5,s6,gsh,grid,basis,valS1,valS2,valS3,valS4,1,A12,B12,C12,A13,B13,C13);

        reduce_4c_ol1(s1,s2,s1,s2,s3,s4,s5,s6,gsh,valS1,valS2,valS3,valS4,N,iN,tmp,olp);

       //s1 on atom 1, s23 on atom 2, s4 on atom 3 (type==2)
        init_s14(0,0,s3,s4,0,0,0,0,iN,gsh,valS1,valS2,valS3,valS4,wt);
        eval_s14b(0,0,s3,s4,0,0,0,0,gsh,grid,basis,valS1,valS2,valS3,valS4,2,A12,B12,C12,A13,B13,C13);

        reduce_4c_ol1(s1,s2,s3,s4,s3,s4,s5,s6,gsh,valS1,valS2,valS3,valS4,N,iN,tmp,olp);

       //s1 on atom 1, s2 on atom 2, s34 on atom 3 (type==3)
       // works but most basis ftns should be on atoms 1+2
        init_s14(0,0,0,0,s5,s6,0,0,iN,gsh,valS1,valS2,valS3,valS4,wt);
        eval_s14b(0,0,0,0,s5,s6,0,0,gsh,grid,basis,valS1,valS2,valS3,valS4,3,A12,B12,C12,A13,B13,C13);

       //s1 on atom 1, s2 on atom 3, s34 on atom 2 (type==4)
       //buggy
        //init_s14(0,0,s3,s4,s5,s6,s5,s6,iN,gsh,valS1,valS2,valS3,valS4,wt);
        //eval_s14b(0,0,s3,s4,s5,s6,s5,s6,gsh,grid,basis,valS1,valS2,valS3,valS4,4,A12,B12,C12,A13,B13,C13);

        reduce_4c_ol1(s1,s2,s3,s4,s5,s6,s5,s6,gsh,valS1,valS2,valS3,valS4,N,iN,tmp,olp);

      } //loop p over third atom

    } //loop n over second atom

   //four-atom case
    for (int n=m+1;n<natoms;n++)
    {
      double Z2 = (double)atno[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;
      coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;

      for (int p=n+1;p<natoms;p++)
      if (p!=m)
      for (int q=p+1;q<natoms;q++)
      if (q!=m && q!=n)
      {
        double Z3 = (double)atno[p];
        double A3 = coords[3*p+0]; double B3 = coords[3*p+1]; double C3 = coords[3*p+2];
        double A13 = A3-A1; double B13 = B3-B1; double C13 = C3-C1;
        coordn[6] = A13; coordn[7] = B13; coordn[8] = C13;

        double Z4 = (double)atno[q];
        double A4 = coords[3*q+0]; double B4 = coords[3*q+1]; double C4 = coords[3*q+2];
        double A14 = A4-A1; double B14 = B4-B1; double C14 = C4-C1;
        coordn[9] = A14; coordn[10] = B14; coordn[11] = C14;

       //CPMZ need to fix this up
        generate_ps_quad_grid(0.,4,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
        add_r1_to_grid(gshh,grid,0.,0.,0.);

        int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];
        int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];
        int s5 = 0; if (p>0) s5 = n2i[p-1]; int s6 = n2i[p];
        int s7 = 0; if (q>0) s7 = n2i[q-1]; int s8 = n2i[q];

        //printf("  mnpq: %i %i %i %i   s12: %2i %2i  s34: %2i %2i  s56: %2i %2i  s78: %2i %2i \n",m,n,p,q,s1,s2,s3,s4,s5,s6,s7,s8);

       //s1 on atom 1, s2 on atom 2, s3 on atom 3, s4 on atom 4
        init_s14(s1,s2,s3,s4,s5,s6,s7,s8,iN,gshh,valS1,valS2,valS3,valS4,wt);
        eval_s14c(s1,s2,s3,s4,s5,s6,s7,s8,gshh,grid,basis,valS1,valS2,valS3,valS4,A12,B12,C12,A13,B13,C13,A14,B14,C14);

        reduce_4c_ol1(s1,s2,s3,s4,s5,s6,s7,s8,gshh,valS1,valS2,valS3,valS4,N,iN,tmp,olp);

      } //loop pq over third and fourth atoms
    } //loop n over second atom
  } //loop m over natoms
  acc_set_device_num(0,acc_device_nvidia);


  double norm[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];
  //#pragma acc enter data copyin(norm[0:N])

 #if 0
 //olp no longer stored on GPU
  if (nomp>1)
  {
    double* olt = new double[N4]();
    for (int j=0;j<N4;j++) olt[j] = 0.;

    for (int n=0;n<nomp;n++)
    {
      int tid = n;
      acc_set_device_num(tid,acc_device_nvidia);

      #pragma acc update self(olp[0:N4])
      for (int j=0;j<N4;j++)
        olt[j] += olp[j];
    }
    for (int j=0;j<N4;j++)
      olp[j] = olt[j];

    acc_set_device_num(0,acc_device_nvidia);
    #pragma acc update device(olp[0:N4])

    delete [] olt;
  }
 #endif

  copy_symm_4c_ps_cpu(natoms,n2i,N,olp);

 //#pragma acc parallel loop collapse(4) present(olp[0:N4],norm[0:N])
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
  for (int l=0;l<N;l++)
  {
    double n1234 = norm[i]*norm[j]*norm[k]*norm[l];
    olp[i*N3+j*N2+k*N+l] *= n1234;
  }
  //#pragma acc update self(olp[0:N4])

  for (int j=0;j<N4;j++)
    ol[j] = olp[j];

  const double lt = 1.e-15;
  for (int j=0;j<N4;j++)
  if (fabs(ol[j])<lt)
    ol[j] = 0.;

  if (prl>2 || (prl>0 && N<4))
  {
    printf("\n olp: \n");
    for (int i=0;i<N2;i++)
    {
      print_square(N,&olp[i*N2]);
    }
  }

  #pragma acc exit data delete(norm[0:N])

 #if USE_ACC
 #pragma omp parallel for schedule(dynamic) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = n;
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc exit data delete(olp[0:N4])
    #pragma acc exit data delete(grid[0:gs6],wt[0:gshh])
    #pragma acc exit data delete(valS1[0:iN][0:gshh],valS2[0:iN][0:gshh],valS3[0:iN][0:gshh],valS4[0:iN][0:gshh])
    #pragma acc exit data delete(valt[0:gshh])
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  delete [] n2i;

  for (int i=0;i<iN;i++) delete [] valS1[i];
  for (int i=0;i<iN;i++) delete [] valS2[i];
  for (int i=0;i<iN;i++) delete [] valS3[i];
  for (int i=0;i<iN;i++) delete [] valS4[i];
  delete [] valS1; delete [] valS2; delete [] valS3; delete [] valS4;
  delete [] valt;
  delete [] olp;

  delete [] grid;
  delete [] wt;

  return;
}

void compute_4c_ol_ps(int natoms, int* atno, double* coords, vector<vector<double> > &basis, int quad_order, int quad_r_order, int nmu, int nnu, int nphi, double* ol, int prl)
{
  int Nmax = 50;
  if (basis.size()<Nmax && quad_order<10)
    return compute_4c_ol_ps_fast(natoms,atno,coords,basis,quad_order,quad_r_order,nmu,nnu,nphi,ol,prl);

  if (prl>-1) printf("  beginning compute_4c_ol_ps \n");

  int nomp_max = 1;
 #pragma omp parallel
  nomp_max = omp_get_num_threads();
  
  int ngpu = 0;
 #if USE_ACC
  ngpu = acc_get_num_devices(acc_device_nvidia);
  int nomp = ngpu;
 #else
  int nomp = nomp_max;
 #endif

  int N = basis.size();
  int N2 = N*N;
  int N3 = N2*N;
  int N4 = N3*N;

//  quad_r_order = quad_order;
  int qos = quad_order*quad_order*quad_order;
  int qoh = quad_r_order;
  int qosh = qoh*qoh*qoh;
  int gs = nmu*nnu*nphi*qos;
  int gsh = (nmu*nnu*nphi-8)*qos+8*qosh;
  int gshh = (nmu*nnu*nphi-16)*qos+16*qosh;
  int gs6 = 6*gshh;

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

  double* grid = new double[gs6];
  double* wt = new double[gshh];

  double gpumem = (double)acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  double togb = 1./1024./1024./1024.;

  printf("  gs: %8i  gpu mem total: %6.3f GB ngpu: %2i \n",gshh,gpumem*togb,ngpu); 

  vector<vector<int> > n2ip;
  int imaxN = get_imax_n2ip(Nmax,natoms,N,basis,n2ip);
  if (prl>1) printf("   imaxN: %2i \n",imaxN);

 //needed for copy_symm ftn
  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);
 //enables memory savings
  iN = imaxN;

 //intermediate storage
  double** valS1 = new double*[iN]; for (int i=0;i<iN;i++) valS1[i] = new double[gshh];
  double** valS2 = new double*[iN]; for (int i=0;i<iN;i++) valS2[i] = new double[gshh];
  double** valS3 = new double*[iN]; for (int i=0;i<iN;i++) valS3[i] = new double[gshh];
  double** valS4 = new double*[iN]; for (int i=0;i<iN;i++) valS4[i] = new double[gshh];
  double* valt = new double[gshh];
 //CPMZ reconfiguring this
  int M4 = iN*iN*iN*iN;
  double* tmp = new double[M4];
  double* olp = new double[N4]();

 #if USE_ACC
 #pragma omp parallel for schedule(dynamic) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = n;
    acc_set_device_num(tid,acc_device_nvidia);

    //#pragma acc enter data create(olp[0:N4])
    #pragma acc enter data create(tmp[0:M4])

    #pragma acc enter data create(grid[0:gs6],wt[0:gshh])
    #pragma acc enter data create(valS1[0:iN][0:gshh],valS2[0:iN][0:gshh],valS3[0:iN][0:gshh],valS4[0:iN][0:gshh])
    #pragma acc enter data create(valt[0:gshh])

    acc_assign(M4,tmp,0.);
    //acc_assign(N4,olp,0.);
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  double gpumem_2 = 1.*acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  printf("   after alloc, gpu memory available: %6.3f GB \n",gpumem_2*togb); 

 #pragma omp parallel for schedule(dynamic) num_threads(nomp)
  for (int m=0;m<natoms;m++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);

    double Z1 = (double)atno[m];
    //double Z1 = atnod[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];
    double coordn[12];
    coordn[0] = coordn[1] = coordn[2] = 0.;

    generate_ps_quad_grid(Z1,1,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
    add_r1_to_grid(gs,grid,0.,0.,0.);

   //all basis on one atom
    for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
    for (int sp2=0;sp2<n2ip[m].size()-1;sp2++)
    for (int sp3=0;sp3<n2ip[m].size()-1;sp3++)
    for (int sp4=0;sp4<n2ip[m].size()-1;sp4++)
    {
      int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
      int s3 = n2ip[m][sp2]; int s4 = n2ip[m][sp2+1];
      int s5 = n2ip[m][sp3]; int s6 = n2ip[m][sp3+1];
      int s7 = n2ip[m][sp4]; int s8 = n2ip[m][sp4+1];

      //printf("    m: %i   s12: %2i %2i s34: %2i %2i s56: %2i %2i s78: %2i %2i \n",m,s1,s2,s3,s4,s5,s6,s7,s8);

      init_s14(s1,s2,s3,s4,s5,s6,s7,s8,iN,gs,valS1,valS2,valS3,valS4,wt);
      eval_s14(s1,s2,s3,s4,s5,s6,s7,s8,gs,grid,basis,valS1,valS2,valS3,valS4,0,0.,0.,0.);

      reduce_4c_ol1(s1,s2,s3,s4,s5,s6,s7,s8,gs,valS1,valS2,valS3,valS4,N,iN,tmp,olp);
    }

   //two-atom ints
    for (int n=m+1;n<natoms;n++)
    {
      double Z2 = (double)atno[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;
      coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;

      generate_ps_quad_grid(0.,2,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
      add_r1_to_grid(gs,grid,0.,0.,0.);

      //int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];
      for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
      for (int sp2=0;sp2<n2ip[n].size()-1;sp2++)
      for (int sp3=0;sp3<n2ip[n].size()-1;sp3++)
      for (int sp4=0;sp4<n2ip[n].size()-1;sp4++)
      {
        int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
        int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2+1];
        int s5 = n2ip[n][sp3]; int s6 = n2ip[n][sp3+1];
        int s7 = n2ip[n][sp4]; int s8 = n2ip[n][sp4+1];

        //printf("     mn: %i %i   s12: %2i %2i s34: %2i %2i s56: %2i %2i s78: %2i %2i \n",m,n,s1,s2,s3,s4,s5,s6,s7,s8);

       //s1 on atom 1, s234 on atom 2
        init_s14(s1,s2,s3,s4,s5,s6,s7,s8,iN,gs,valS1,valS2,valS3,valS4,wt);
        eval_s14(s1,s2,s3,s4,s5,s6,s7,s8,gs,grid,basis,valS1,valS2,valS3,valS4,1,A12,B12,C12);

        reduce_4c_ol1(s1,s2,s3,s4,s5,s6,s7,s8,gs,valS1,valS2,valS3,valS4,N,iN,tmp,olp);
      }

      for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
      for (int sp2=0;sp2<n2ip[m].size()-1;sp2++)
      for (int sp3=0;sp3<n2ip[n].size()-1;sp3++)
      for (int sp4=0;sp4<n2ip[n].size()-1;sp4++)
      {
        int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
        int s3 = n2ip[m][sp2]; int s4 = n2ip[m][sp2+1];
        int s5 = n2ip[n][sp3]; int s6 = n2ip[n][sp3+1];
        int s7 = n2ip[n][sp4]; int s8 = n2ip[n][sp4+1];

        //printf("     mn: %i   s12: %2i %2i s34: %2i %2i s56: %2i %2i s78: %2i %2i \n",m,n,s1,s2,s3,s4,s5,s6,s7,s8);

       //s12 on atom 1, s34 on atom 2
        init_s14(s1,s2,s3,s4,s5,s6,s7,s8,iN,gs,valS1,valS2,valS3,valS4,wt);
        eval_s14(s1,s2,s3,s4,s5,s6,s7,s8,gs,grid,basis,valS1,valS2,valS3,valS4,2,A12,B12,C12);

        reduce_4c_ol1(s1,s2,s3,s4,s5,s6,s7,s8,gs,valS1,valS2,valS3,valS4,N,iN,tmp,olp);
      }

      for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
      for (int sp2=0;sp2<n2ip[m].size()-1;sp2++)
      for (int sp3=0;sp3<n2ip[m].size()-1;sp3++)
      for (int sp4=0;sp4<n2ip[n].size()-1;sp4++)
      {
        int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
        int s3 = n2ip[m][sp2]; int s4 = n2ip[m][sp2+1];
        int s5 = n2ip[m][sp3]; int s6 = n2ip[m][sp3+1];
        int s7 = n2ip[n][sp4]; int s8 = n2ip[n][sp4+1];

        //printf("     mn: %i %i   s12: %2i %2i s34: %2i %2i s56: %2i %2i s78: %2i %2i \n",m,n,s1,s2,s3,s4,s5,s6,s7,s8);

       //s123 on atom 1, s4 on atom 2
        init_s14(s1,s2,s3,s4,s5,s6,s7,s8,iN,gs,valS1,valS2,valS3,valS4,wt);
        eval_s14(s1,s2,s3,s4,s5,s6,s7,s8,gs,grid,basis,valS1,valS2,valS3,valS4,3,A12,B12,C12);

        reduce_4c_ol1(s1,s2,s3,s4,s5,s6,s7,s8,gs,valS1,valS2,valS3,valS4,N,iN,tmp,olp);
      } //loop sp12

    } //loop n over second atom

   //three-atom case
    for (int n=m+1;n<natoms;n++)
    {
      double Z2 = (double)atno[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;
      coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;

      for (int p=n+1;p<natoms;p++)
      if (p!=m)
      {
        double Z3 = (double)atno[p];
        double A3 = coords[3*p+0]; double B3 = coords[3*p+1]; double C3 = coords[3*p+2];
        double A13 = A3-A1; double B13 = B3-B1; double C13 = C3-C1;
        coordn[6] = A13; coordn[7] = B13; coordn[8] = C13;

        //printf("  mnp: %i %i %i   s12: %2i %2i  s34: %2i %2i  s56: %2i %2i \n",m,n,p,s1,s2,s3,s4,s5,s6);

        double ztm = 0.; int lm = 0;
        for (int sp4=0;sp4<n2ip[p].size()-1;sp4++)
        {
          int s7 = n2ip[p][sp4]; int s8 = n2ip[p][sp4+1];
          get_ztm_lm(s7,s8,basis,ztm,lm);
        }

        generate_ps_quad_grid(0.,3,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
        add_r1_to_grid(gsh,grid,0.,0.,0.);

       //need all unique combinations, as each atom is unique
        for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
        for (int sp2=0;sp2<n2ip[m].size()-1;sp2++)
        for (int sp3=0;sp3<n2ip[n].size()-1;sp3++)
        for (int sp4=0;sp4<n2ip[p].size()-1;sp4++)
        {
          int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
          int s3 = n2ip[m][sp2]; int s4 = n2ip[m][sp2+1];
          int s5 = n2ip[n][sp3]; int s6 = n2ip[n][sp3+1];
          int s7 = n2ip[p][sp4]; int s8 = n2ip[p][sp4+1];

         //s12 on atom 1, s3 on atom 2, s4 on atom 3 (type==1)
          init_s14(s1,s2,s3,s4,s5,s6,s7,s8,iN,gsh,valS1,valS2,valS3,valS4,wt);
          eval_s14b(s1,s2,s3,s4,s5,s6,s7,s8,gsh,grid,basis,valS1,valS2,valS3,valS4,1,A12,B12,C12,A13,B13,C13);

          reduce_4c_ol1(s1,s2,s3,s4,s5,s6,s7,s8,gsh,valS1,valS2,valS3,valS4,N,iN,tmp,olp);
        }

        for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
        for (int sp2=0;sp2<n2ip[n].size()-1;sp2++)
        for (int sp3=0;sp3<n2ip[n].size()-1;sp3++)
        for (int sp4=0;sp4<n2ip[p].size()-1;sp4++)
        {
          int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
          int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2+1];
          int s5 = n2ip[n][sp3]; int s6 = n2ip[n][sp3+1];
          int s7 = n2ip[p][sp4]; int s8 = n2ip[p][sp4+1];

         //s1 on atom 1, s23 on atom 2, s4 on atom 3 (type==2)
          init_s14(s1,s2,s3,s4,s5,s6,s7,s8,iN,gsh,valS1,valS2,valS3,valS4,wt);
          eval_s14b(s1,s2,s3,s4,s5,s6,s7,s8,gsh,grid,basis,valS1,valS2,valS3,valS4,2,A12,B12,C12,A13,B13,C13);

          reduce_4c_ol1(s1,s2,s3,s4,s5,s6,s7,s8,gsh,valS1,valS2,valS3,valS4,N,iN,tmp,olp);
        }

        for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
        for (int sp2=0;sp2<n2ip[n].size()-1;sp2++)
        for (int sp3=0;sp3<n2ip[p].size()-1;sp3++)
        for (int sp4=0;sp4<n2ip[p].size()-1;sp4++)
        {
          int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
          int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2+1];
          int s5 = n2ip[p][sp3]; int s6 = n2ip[p][sp3+1];
          int s7 = n2ip[p][sp4]; int s8 = n2ip[p][sp4+1];

         //s1 on atom 1, s2 on atom 2, s34 on atom 3 (type==3)
         // works but most basis ftns should be on atoms 1+2
          init_s14(s1,s2,s3,s4,s5,s6,s7,s8,iN,gsh,valS1,valS2,valS3,valS4,wt);
          eval_s14b(s1,s2,s3,s4,s5,s6,s7,s8,gsh,grid,basis,valS1,valS2,valS3,valS4,3,A12,B12,C12,A13,B13,C13);

          reduce_4c_ol1(s1,s2,s3,s4,s5,s6,s7,s8,gsh,valS1,valS2,valS3,valS4,N,iN,tmp,olp);
        }

      } //loop p over third atom

    } //loop n over second atom

   //four-atom case
    for (int n=m+1;n<natoms;n++)
    {
      double Z2 = (double)atno[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;
      coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;

      for (int p=n+1;p<natoms;p++)
      if (p!=m)
      for (int q=p+1;q<natoms;q++)
      if (q!=m && q!=n)
      {
        double Z3 = (double)atno[p];
        double A3 = coords[3*p+0]; double B3 = coords[3*p+1]; double C3 = coords[3*p+2];
        double A13 = A3-A1; double B13 = B3-B1; double C13 = C3-C1;
        coordn[6] = A13; coordn[7] = B13; coordn[8] = C13;

        double Z4 = (double)atno[q];
        double A4 = coords[3*q+0]; double B4 = coords[3*q+1]; double C4 = coords[3*q+2];
        double A14 = A4-A1; double B14 = B4-B1; double C14 = C4-C1;
        coordn[9] = A14; coordn[10] = B14; coordn[11] = C14;

       //CPMZ need to fix this up
        generate_ps_quad_grid(0.,4,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
        add_r1_to_grid(gshh,grid,0.,0.,0.);

        //printf("  mnpq: %i %i %i %i   s12: %2i %2i  s34: %2i %2i  s56: %2i %2i  s78: %2i %2i \n",m,n,p,q,s1,s2,s3,s4,s5,s6,s7,s8);

        for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
        for (int sp2=0;sp2<n2ip[n].size()-1;sp2++)
        for (int sp3=0;sp3<n2ip[p].size()-1;sp3++)
        for (int sp4=0;sp4<n2ip[q].size()-1;sp4++)
        {
          int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
          int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2+1];
          int s5 = n2ip[p][sp3]; int s6 = n2ip[p][sp3+1];
          int s7 = n2ip[q][sp4]; int s8 = n2ip[q][sp4+1];

         //s1 on atom 1, s2 on atom 2, s3 on atom 3, s4 on atom 4
          init_s14(s1,s2,s3,s4,s5,s6,s7,s8,iN,gshh,valS1,valS2,valS3,valS4,wt);
          eval_s14c(s1,s2,s3,s4,s5,s6,s7,s8,gshh,grid,basis,valS1,valS2,valS3,valS4,A12,B12,C12,A13,B13,C13,A14,B14,C14);

          reduce_4c_ol1(s1,s2,s3,s4,s5,s6,s7,s8,gshh,valS1,valS2,valS3,valS4,N,iN,tmp,olp);

        } //loop sp

      } //loop pq over third and fourth atoms
    } //loop n over second atom
  } //loop m over natoms
  acc_set_device_num(0,acc_device_nvidia);


  double norm[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];
  //#pragma acc enter data copyin(norm[0:N])

 #if 0
 //olp no longer stored on GPU
  if (nomp>1)
  {
    double* olt = new double[N4]();
    for (int j=0;j<N4;j++) olt[j] = 0.;

    for (int n=0;n<nomp;n++)
    {
      int tid = n;
      acc_set_device_num(tid,acc_device_nvidia);

      #pragma acc update self(olp[0:N4])
      for (int j=0;j<N4;j++)
        olt[j] += olp[j];
    }
    for (int j=0;j<N4;j++)
      olp[j] = olt[j];

    acc_set_device_num(0,acc_device_nvidia);
    #pragma acc update device(olp[0:N4])

    delete [] olt;
  }
 #endif

  copy_symm_4c_ps_cpu(natoms,n2i,N,olp);

 //#pragma acc parallel loop collapse(4) present(olp[0:N4],norm[0:N])
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
  for (int l=0;l<N;l++)
  {
    double n1234 = norm[i]*norm[j]*norm[k]*norm[l];
    olp[i*N3+j*N2+k*N+l] *= n1234;
  }
  //#pragma acc update self(olp[0:N4])

  for (int j=0;j<N4;j++)
    ol[j] = olp[j];

  const double lt = 1.e-15;
  for (int j=0;j<N4;j++)
  if (fabs(ol[j])<lt)
    ol[j] = 0.;

  if (prl>2 || (prl>0 && N<4))
  {
    printf("\n olp: \n");
    for (int i=0;i<N2;i++)
    {
      print_square(N,&olp[i*N2]);
    }
  }

  #pragma acc exit data delete(norm[0:N])

 #if USE_ACC
 #pragma omp parallel for schedule(dynamic) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = n;
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc exit data delete(olp[0:N4])
    #pragma acc exit data delete(tmp[0:M4])
    #pragma acc exit data delete(grid[0:gs6],wt[0:gshh])
    #pragma acc exit data delete(valS1[0:iN][0:gshh],valS2[0:iN][0:gshh],valS3[0:iN][0:gshh],valS4[0:iN][0:gshh])
    #pragma acc exit data delete(valt[0:gshh])
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  delete [] n2i;

  for (int i=0;i<iN;i++) delete [] valS1[i];
  for (int i=0;i<iN;i++) delete [] valS2[i];
  for (int i=0;i<iN;i++) delete [] valS3[i];
  for (int i=0;i<iN;i++) delete [] valS4[i];
  delete [] valS1; delete [] valS2; delete [] valS3; delete [] valS4;
  delete [] valt;
  delete [] tmp;
  delete [] olp;

  delete [] grid;
  delete [] wt;

  return;
}
