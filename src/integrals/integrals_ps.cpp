#include "integrals.h"
#include "prosph.h"
#include "quad.h"

#define SWITCH_23 0
#define OMP_PARA 1

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

void reduce_3cenp(int tid, double Z3, int s1, int s2, int s3, int s4, int N, int iN, int gs, double* grid, double* valS1, double* valS2, double* valt, double* wt, double* pVp1)
{
  int gs3 = 3*gs;
  int gs6 = 6*gs;
  int N2 = N*N;
  int igs3 = iN*gs3;

 #pragma acc parallel loop present(grid[0:gs6],wt[0:gs],valt[0:gs]) async(tid+1)
  for (int j=0;j<gs;j++)
  {
    double R = grid[6*j+3];
    valt[j] = -Z3/R*wt[j];
  }

 #pragma acc parallel loop collapse(2) present(valS1[0:igs3],valS2[0:igs3],valt[0:gs],pVp1[0:N2]) async(tid+1)
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3;
    double* valm = &valS1[ii1*gs3];
    double* valn = &valS2[ii2*gs3];

    double val = 0.;
   #pragma acc loop reduction(+:val)
    for (int j=0;j<gs;j++)
    {
      double dp = valm[3*j+0]*valn[3*j+0] + valm[3*j+1]*valn[3*j+1] + valm[3*j+2]*valn[3*j+2];
      val += dp*valt[j];
    }

    pVp1[i1*N+i2] += val;
  }

 /*
  if (tid<0)
  {
    acc_wait_all();
    //#pragma acc wait
  }
 */

  return;
}

void reduce_3cen(int tid, double Z3, int s1, int s2, int s3, int s4, int N, int iN, int gs, double* grid, double* valS1, double* valS2, double* valt, double* wt, double* En1)
{
  int gs6 = 6*gs;
  int N2 = N*N;
  int igs = iN*gs;

  if (wt==NULL)
 #pragma acc parallel loop present(grid[0:gs6],valt[0:gs]) async(tid+1)
  for (int j=0;j<gs;j++)
  {
    double R = grid[6*j+3];
    valt[j] = -Z3/R;
  }

  if (wt!=NULL)
 #pragma acc parallel loop present(grid[0:gs6],wt[0:gs],valt[0:gs]) async(tid+1)
  for (int j=0;j<gs;j++)
  {
    double R = grid[6*j+3];
    valt[j] = -Z3/R*wt[j];

    //printf(" b13: %i %i   R: %8.5f \n",s1,s3,R);
  }

 #pragma acc parallel loop collapse(2) present(valS1[0:igs],valS2[0:igs],valt[0:gs],En1[0:N2]) async(tid+1)
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3;
    double* valm = &valS1[ii1*gs];
    double* valn = &valS2[ii2*gs];

    double val = 0.;
   #pragma acc loop reduction(+:val)
    for (int j=0;j<gs;j++)
      val += valm[j]*valn[j]*valt[j];

    En1[i1*N+i2] += val;
  }

 /*
  if (tid<0)
  {
    acc_wait_all();
    //#pragma acc wait
  }
 */

  return;
}

void init_s12nw(int tid, int s1, int s2, int s3, int s4, int iN, int gs, double* valS1, double* valS2)
{
  int igs = iN*gs;
 #pragma acc parallel loop collapse(2) present(valS1[0:igs]) async(tid+1)
  for (int ii1=0;ii1<s2-s1;ii1++)
  {
    for (int j=0;j<gs;j++)
      valS1[ii1*gs+j] = 1.;
  }

 #pragma acc parallel loop collapse(2) present(valS2[0:igs]) async(tid+1)
  for (int ii1=0;ii1<s4-s3;ii1++)
  {
    for (int j=0;j<gs;j++)
      valS2[ii1*gs+j] = 1.;
  }
}

void eval_s12(int tid, int s1, int s2, int s3, int s4, vector<vector<double> >& basis, int iN, int gs, double* grid, double* valS1, double* valS2)
{
  for (int i1=s1;i1<s2;i1++)
  {
    int ii1 = i1-s1;

    vector<double> basis1 = basis[i1];
    int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

   //S
    eval_shd(tid,gs,grid,&valS1[ii1*gs],n1,l1,m1,zeta1);
  }

  for (int i1=s3;i1<s4;i1++)
  {
    int ii1 = i1-s3;

    vector<double> basis1 = basis[i1];
    int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

   //S
    eval_shd(tid,gs,grid,&valS2[ii1*gs],n1,l1,m1,zeta1);
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

  int ngpu = 0;
 #if USE_ACC
  ngpu = acc_get_num_devices(acc_device_nvidia);
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
  double togb = pow(1024.,-3.);

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

  if (Nmax<=0) { printf("\n ERROR: couldn't calculate gpu memory requirements (mem0: %5.3f mem1: %5.3f) \n",mem0*togb,8.*(3.*gs*Nmax + mem0 + 3.*N2 + 4.*gs6 + 1.*gs)*togb); exit(-1); }

  vector<vector<int> > n2ip;
  int imaxN = get_imax_n2ip(Nmax,natoms,N,basis,n2ip);
  if (prl>1) printf("   imaxN: %2i \n",imaxN);

  double gpumem_gb = gpumem/1024./1024./1024.;
  printf("   gpu memory available: %6.3f GB \n",gpumem_gb);

  //int* n2i = new int[natoms];
  //int imaxN = get_imax_n2i(natoms,N,basis,n2i);
  //printf(" imaxN: %2i \n",imaxN);

 //intermediate storage
  int iN = imaxN;
  int igs = iN*gs;
  double* valS1 = new double[igs];
  double* valS2 = new double[igs];
  double* valT1 = new double[igs];
  //double** valS1 = new double*[iN]; for (int i=0;i<iN;i++) valS1[i] = new double[gs];
  //double** valS2 = new double*[iN]; for (int i=0;i<iN;i++) valS2[i] = new double[gs];
  //double** valT1 = new double*[iN]; for (int i=0;i<iN;i++) valT1[i] = new double[gs];
  double* valt = new double[gs];

 #if USE_ACC
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc enter data create(S[0:N2],T[0:N2],En[0:N2]) //async(n+1)

    #pragma acc enter data create(grid[0:gs6],wt[0:gs]) //async(n+1)
    //#pragma acc enter data create(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valT1[0:iN][0:gs])
    #pragma acc enter data create(valS1[0:igs],valS2[0:igs],valT1[0:igs]) //async(n+1)
    #pragma acc enter data create(valt[0:gs]) //async(n+1)

    acc_assign(N2,S,0.);
    acc_assign(N2,T,0.);
    acc_assign(N2,En,0.);

    //acc_wait_all();
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  double gpumem_2 = 1.*acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  printf("   after alloc, gpu memory available: %6.3f GB \n",gpumem_2*togb);

 #if OMP_PARA
 #pragma omp parallel for schedule(dynamic,1) num_threads(ngpu)
 #endif
  for (int m=0;m<natoms;m++)
  {
   #if OMP_PARA
    int tid = omp_get_thread_num();
   #else
    int tid = m%ngpu;
   #endif
    acc_set_device_num(tid,acc_device_nvidia);

    double Z1 = atnod[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];
    double coordn[6];
    coordn[0] = coordn[1] = coordn[2] = 0.;

    generate_ps_quad_grid(tid,Z1,1,coordn,quad_order,quad_order,nmu,nnu,nphi,grid,wt);
    //add_r1_to_grid(tid,gs,grid,0.,0.,0.);

   //working on this block of the matrix
    //int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    //j>i
    for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
    for (int sp2=sp1;sp2<n2ip[m].size()-1;sp2++)
    {
      int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
      int s3 = n2ip[m][sp2]; int s4 = n2ip[m][sp2+1];

     #pragma acc parallel loop collapse(2) present(valS1[0:igs]) async(tid+1)
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
        for (int j=0;j<gs;j++)
          valS1[ii1*gs+j] = 1.;
      }

     #pragma acc parallel loop collapse(2) present(valS2[0:igs],wt[0:gs]) async(tid+1)
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
        for (int j=0;j<gs;j++)
          valS2[ii1*gs+j] = wt[j];
      }

     //first compute single atom ints
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

       //S
        eval_shd(tid,gs,grid,&valS1[ii1*gs],n1,l1,m1,zeta1);
      }

      for (int i1=s3;i1<s4;i1++)
      {
        int ii1 = i1-s3;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

       //S
        eval_shd(tid,gs,grid,&valS2[ii1*gs],n1,l1,m1,zeta1);
      }

     #pragma acc parallel loop collapse(2) present(valS1[0:igs],valT1[0:igs]) async(tid+1)
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
        for (int j=0;j<gs;j++)
          valT1[ii1*gs+j] = valS1[ii1*gs+j];
      }

     //KE terms
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_ked(tid,gs,grid,&valT1[ii1*gs],n1,l1,zeta1);
      }

      reduce_2c1(tid,s1,s2,s3,s4,gs,valS1,valS2,iN,N,S);
      reduce_2c1(tid,s1,s2,s3,s4,gs,valT1,valS2,iN,N,T);

     /////////////////////////////////////////////////////////////////////
     //electron-nuclear attraction
      reduce_3cen(tid,Z1,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,NULL,En);
     /////////////////////////////////////////////////////////////////////
    } //loop sp over s12

   //two-atom ints
    for (int n=m+1;n<natoms;n++)
    {
      double Z2 = atnod[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;
      coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;

      generate_ps_quad_grid(tid,0.,2,coordn,quad_order,quad_order,nmu,nnu,nphi,grid,wt);
      add_r1_to_grid(tid,gs,grid,0.,0.,0.);

      //int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

     //basis ftns on each center
      for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
      for (int sp2=0;sp2<n2ip[n].size()-1;sp2++)
      {
        int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
        int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2+1];

       #pragma acc parallel loop collapse(2) present(valS1[0:igs]) async(tid+1)
        for (int ii1=0;ii1<s2-s1;ii1++)
        {
          for (int j=0;j<gs;j++)
            valS1[ii1*gs+j] = 1.;
        }

       #pragma acc parallel loop collapse(2) present(valS2[0:igs],wt[0:gs]) async(tid+1)
        for (int ii1=0;ii1<s4-s3;ii1++)
        {
          for (int j=0;j<gs;j++)
            valS2[ii1*gs+j] = wt[j];
        }

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<double> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

         //S
          eval_shd(tid,gs,grid,&valS1[ii1*gs],n1,l1,m1,zeta1);
        }

       #pragma acc parallel loop collapse(2) present(valS1[0:igs],valT1[0:igs]) async(tid+1)
        for (int ii1=0;ii1<s2-s1;ii1++)
        {
          for (int j=0;j<gs;j++)
            valT1[ii1*gs+j] = valS1[ii1*gs+j];
        }

       //KE terms
        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<double> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

          eval_ked(tid,gs,grid,&valT1[ii1*gs],n1,l1,zeta1);
        }

       //second center
        recenter_grid_zero(tid,gs,grid,-A12,-B12,-C12);

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

         //S
          eval_shd(tid,gs,grid,&valS2[ii2*gs],n2,l2,m2,zeta2);
        }

        reduce_2c1(tid,s1,s2,s3,s4,gs,valS1,valS2,iN,N,S);
        reduce_2c1(tid,s1,s2,s3,s4,gs,valT1,valS2,iN,N,T);

       ///////////////////////////////////////////////////////////////////
       //electron-nuclear attraction centers 1+2
       #pragma acc parallel loop collapse(2) present(valS2[0:igs],wt[0:gs]) async(tid+1)
        for (int ii1=0;ii1<s4-s3;ii1++)
        {
          for (int j=0;j<gs;j++)
            valS2[ii1*gs+j] /= wt[j];
        }

        reduce_3cen(tid,Z2,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,En);
        recenter_grid_zero(tid,gs,grid,A12,B12,C12);
        reduce_3cen(tid,Z1,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,En);
        ////////////////////////////////////////////////////////////////////

      } //loop sp over s14

     //En for basis ftns on same center (1)
      for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
      for (int sp2=sp1;sp2<n2ip[m].size()-1;sp2++)
      {
        int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
        int s3 = n2ip[m][sp2]; int s4 = n2ip[m][sp2+1];

        //printf("  En 2c(1). s12: %2i %2i s34: %2i %2i \n",s1,s2,s3,s4);

        init_s12nw(tid,s1,s2,s3,s4,iN,gs,valS1,valS2);
        eval_s12(tid,s1,s2,s3,s4,basis,iN,gs,grid,valS1,valS2);

       //second center
        recenter_grid_zero(tid,gs,grid,-A12,-B12,-C12);

       ///////////////////////////////////////////////////////////////////
       //electron-nuclear attraction center 2
        reduce_3cen(tid,Z2,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,En);
        ////////////////////////////////////////////////////////////////////

        recenter_grid_zero(tid,gs,grid,A12,B12,C12);

      } //loop sp over s14

     //En for basis ftns on same center (2)
      for (int sp1=0;sp1<n2ip[n].size()-1;sp1++)
      for (int sp2=sp1;sp2<n2ip[n].size()-1;sp2++)
      {
        int s1 = n2ip[n][sp1]; int s2 = n2ip[n][sp1+1];
        int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2+1];

        //printf("  En 2c(2). s12: %2i %2i s34: %2i %2i \n",s1,s2,s3,s4);

        init_s12nw(tid,s1,s2,s3,s4,iN,gs,valS1,valS2);

       //second center
        recenter_grid_zero(tid,gs,grid,-A12,-B12,-C12);

        eval_s12(tid,s1,s2,s3,s4,basis,iN,gs,grid,valS1,valS2);

        recenter_grid_zero(tid,gs,grid,A12,B12,C12);

       ///////////////////////////////////////////////////////////////////
       //electron-nuclear attraction center 1
        reduce_3cen(tid,Z1,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,En);
        ////////////////////////////////////////////////////////////////////

      } //loop sp over s14

    } //loop n over second atom

   #if OMP_PARA
    acc_wait_all();
    #pragma acc wait
   #endif

  } //loop m over natoms

  #pragma omp barrier
 #if 0
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);
    //#pragma acc wait
    acc_wait_all();
  }
 #endif
  acc_set_device_num(0,acc_device_nvidia);

  if (ngpu>1)
  {
   //gather parallel fragments
    double St[N2]; double Tt[N2]; double Ent[N2];
    for (int j=0;j<N2;j++)
      St[j] = Tt[j] = Ent[j] = 0.;

    for (int n=0;n<ngpu;n++)
    {
      acc_set_device_num(n,acc_device_nvidia);

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
    printf("  S diagonals accuracy: ");
    double vmin = 1000.;
    double vmax = -1000.;
    for (int j=0;j<N;j++)
    {
      double v1 = log10(fabs(1.-S[j*N+j])+1.e-16);
      if (v1<vmin) vmin = v1;
      if (v1>vmax) vmax = v1;
    }
    printf("  %5.2f to %5.2f \n",vmin,vmax);
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
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc exit data delete(S[0:N2],T[0:N2],En[0:N2]) //async(n+1)
    #pragma acc exit data delete(grid[0:gs6],wt[0:gs]) //async(n+1)
    //#pragma acc exit data delete(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valT1[0:iN][0:gs])
    #pragma acc exit data delete(valS1[0:igs],valS2[0:igs],valT1[0:igs]) //async(n+1)
    #pragma acc exit data delete(valt[0:gs]) //async(n+1)

    //acc_wait_all();
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  //auto_crash();

  //delete [] n2i;

  //for (int i=0;i<iN;i++) delete [] valS1[i];
  //for (int i=0;i<iN;i++) delete [] valS2[i];
  //for (int i=0;i<iN;i++) delete [] valT1[i];
  delete [] valS1;
  delete [] valS2;
  delete [] valT1;
  delete [] valt;

  delete [] grid;
  delete [] wt;

  return;
}

void compute_pVp_ps(int natoms, int* atno, double* coords, vector<vector<double> > &basis, int quad_order, int nmu, int nnu, int nphi, double* pVp, int prl)
{
 //barely tested
  if (prl>-1) printf("  beginning compute_pVp_ps \n");

  int ngpu = 0;
 #if USE_ACC
  ngpu = acc_get_num_devices(acc_device_nvidia);
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

  vector<vector<int> > n2ip;
  int imaxN = get_imax_n2ip(Nmax,natoms,N,basis,n2ip);
  if (prl>1)printf("   imaxN: %2i \n",imaxN);

  double gpumem_gb = gpumem/1024./1024./1024.;
  printf("   gpu memory available: %6.3f GB \n",gpumem_gb);

 //intermediate storage
  int iN = imaxN;
  int igs3 = iN*gs3;
  double* valS1 = new double[igs3];
  double* valS2 = new double[igs3];
  //double** valS1 = new double*[iN]; for (int i=0;i<iN;i++) valS1[i] = new double[gs3];
  //double** valS2 = new double*[iN]; for (int i=0;i<iN;i++) valS2[i] = new double[gs3];
  double* valt = new double[gs];

 #if USE_ACC
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc enter data create(pVp[0:N2])

    #pragma acc enter data create(grid[0:gs6],wt[0:gs])
    //#pragma acc enter data create(valS1[0:iN][0:gs3],valS2[0:iN][0:gs3])
    #pragma acc enter data create(valS1[0:igs3],valS2[0:igs3])
    #pragma acc enter data create(valt[0:gs])

    acc_assign(N2,pVp,0.);
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  double gpumem_2 = 1.*acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  printf("   after alloc, gpu memory available: %6.3f GB \n",gpumem_2*togb);

 #if OMP_PARA
 #pragma omp parallel for schedule(dynamic,1) num_threads(ngpu)
 #endif
  for (int m=0;m<natoms;m++)
  {
   #if OMP_PARA
    int tid = omp_get_thread_num();
   #else
    int tid = m%ngpu;
   #endif
    acc_set_device_num(tid,acc_device_nvidia);

    double Z1 = atnod[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];
    double coordn[6];
    coordn[0] = coordn[1] = coordn[2] = 0.;

    generate_ps_quad_grid(tid,Z1,1,coordn,quad_order,quad_order,nmu,nnu,nphi,grid,wt);
    add_r1_to_grid(tid,gs,grid,0.,0.,0.);

    //j>i
    for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
    for (int sp2=sp1;sp2<n2ip[m].size()-1;sp2++)
    {
      int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
      int s3 = n2ip[m][sp2]; int s4 = n2ip[m][sp2+1];

      //printf(" pVp_2c s12: %2i %2i s34: %2i %2i \n",s1,s2,s3,s4);

     #pragma acc parallel loop collapse(2) present(valS1[0:igs3]) async(tid+1)
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
        for (int j=0;j<gs3;j++)
          valS1[ii1*gs3+j] = 1.;
      }

     #pragma acc parallel loop collapse(2) present(valS2[0:gs3]) async(tid+1)
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
        for (int j=0;j<gs3;j++)
          valS2[ii1*gs3+j] = 1.;
      }

     //first compute single atom ints
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_pd(tid,gs,grid,&valS1[ii1*gs3],n1,l1,m1,zeta1);
      }

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;

        vector<double> basis1 = basis[i2];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_pd(tid,gs,grid,&valS2[ii2*gs3],n1,l1,m1,zeta1);
      }

     /////////////////////////////////////////////////////////////////////
     //electron-nuclear attraction
      reduce_3cenp(tid,Z1,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,pVp);
     /////////////////////////////////////////////////////////////////////

    } //loop sp over s12

   //two-atom ints
    for (int n=m+1;n<natoms;n++)
    {
      double Z2 = atnod[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;
      coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;

      generate_ps_quad_grid(tid,0.,2,coordn,quad_order,quad_order,nmu,nnu,nphi,grid,wt);
      add_r1_to_grid(tid,gs,grid,0.,0.,0.);

     //basis functions on each center
      for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
      for (int sp2=0;sp2<n2ip[n].size()-1;sp2++)
      {
        int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
        int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2+1];

        init_s12nw(tid,s1,s2,s3,s4,iN,gs3,valS1,valS2);

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<double> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

         //p
          eval_pd(tid,gs,grid,&valS1[ii1*gs3],n1,l1,m1,zeta1);
        }

       //second center
        recenter_grid_zero(tid,gs,grid,-A12,-B12,-C12);

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

         //p
          eval_pd(tid,gs,grid,&valS2[ii2*gs3],n2,l2,m2,zeta2);
        }

       ///////////////////////////////////////////////////////////////////
       //electron-nuclear attraction centers 1+2

        reduce_3cenp(tid,Z2,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,pVp);
        recenter_grid_zero(tid,gs,grid,A12,B12,C12);
        reduce_3cenp(tid,Z1,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,pVp);
        ////////////////////////////////////////////////////////////////////

      } //loop sp over s14

     //basis functions on center 1
      for (int sp1=0;sp1<n2ip[m].size()-1;sp1++)
      for (int sp2=sp1;sp2<n2ip[m].size()-1;sp2++)
      {
        int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
        int s3 = n2ip[m][sp2]; int s4 = n2ip[m][sp2+1];

        init_s12nw(tid,s1,s2,s3,s4,iN,gs3,valS1,valS2);

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<double> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

         //p
          eval_pd(tid,gs,grid,&valS1[ii1*gs3],n1,l1,m1,zeta1);
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

         //p
          eval_pd(tid,gs,grid,&valS2[ii2*gs3],n2,l2,m2,zeta2);
        }

       //second center
        recenter_grid_zero(tid,gs,grid,-A12,-B12,-C12);

       ///////////////////////////////////////////////////////////////////
       //electron-nuclear attraction center 2
        reduce_3cenp(tid,Z2,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,pVp);
        ////////////////////////////////////////////////////////////////////

        recenter_grid_zero(tid,gs,grid,A12,B12,C12);

      } //loop sp over s14

     //basis functions on center 2
      for (int sp1=0;sp1<n2ip[n].size()-1;sp1++)
      for (int sp2=0;sp2<n2ip[n].size()-1;sp2++)
      {
        int s1 = n2ip[n][sp1]; int s2 = n2ip[n][sp1+1];
        int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2+1];

        init_s12nw(tid,s1,s2,s3,s4,iN,gs3,valS1,valS2);

       //second center
        recenter_grid_zero(tid,gs,grid,-A12,-B12,-C12);

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<double> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

         //p
          eval_pd(tid,gs,grid,&valS1[ii1*gs3],n1,l1,m1,zeta1);
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

         //p
          eval_pd(tid,gs,grid,&valS2[ii2*gs3],n2,l2,m2,zeta2);
        }

        recenter_grid_zero(tid,gs,grid,A12,B12,C12);

       ///////////////////////////////////////////////////////////////////
       //electron-nuclear attraction center 1
        reduce_3cenp(tid,Z1,s1,s2,s3,s4,N,iN,gs,grid,valS1,valS2,valt,wt,pVp);
        ////////////////////////////////////////////////////////////////////

      } //loop sp over s14

    } //loop n over second atom

   #if OMP_PARA
    #pragma acc wait
   #endif

  } //loop m over natoms

  #pragma omp barrier
 #if 0
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);
    //#pragma acc wait
    acc_wait_all();
  }
 #endif
  acc_set_device_num(0,acc_device_nvidia);

  if (ngpu>1)
  {
   //gather parallel fragments
    double pVpt[N2];
    for (int j=0;j<N2;j++)
      pVpt[j] = 0.;

    for (int n=0;n<ngpu;n++)
    {
      acc_set_device_num(n,acc_device_nvidia);

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
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc exit data delete(pVp[0:N2])
    #pragma acc exit data delete(grid[0:gs6],wt[0:gs])
    //#pragma acc exit data delete(valS1[0:iN][0:gs3],valS2[0:iN][0:gs3])
    #pragma acc exit data delete(valS1[0:igs3],valS2[0:igs3])
    #pragma acc exit data delete(valt[0:gs])
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  //for (int i=0;i<iN;i++) delete [] valS1[i];
  //for (int i=0;i<iN;i++) delete [] valS2[i];
  delete [] valS1; delete [] valS2;
  delete [] valt;

  delete [] grid;
  delete [] wt;

  return;
}

 //CPMZ some debug code
void give_me_an_error(bool do_overlap, bool do_yukawa, double gamma, int natoms, int* atno, double* coords, vector<vector<double> > &basis, int quad_order, int nmu, int nnu, int nphi, double* A, int prl)
{
  return;

  int N = basis.size();
  int N2 = N*N;

  int gs = 100000;
  int gs6 = 6*gs;
  int igs = N*gs;

  int ngpu = 0;
 #if USE_ACC
  ngpu = acc_get_num_devices(acc_device_nvidia);
 #endif

  double* grid = new double[gs6];
  double* wt = new double[gs];

  double* valS1 = new double[igs];
  double* valV2 = new double[igs];

 #if USE_ACC
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc enter data create(A[0:N2])

    #pragma acc enter data create(grid[0:gs6],wt[0:gs])
    #pragma acc enter data create(valS1[0:igs],valV2[0:igs])

    acc_assign(N2,A,0.);
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  printf(" we will do nothing \n");

#if USE_ACC
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc exit data delete(A[0:N2])
    #pragma acc exit data delete(grid[0:gs6],wt[0:gs])
    //#pragma acc exit data delete(valS1[0:iN][0:gs],valV2[0:iN][0:gs])
    #pragma acc exit data delete(valS1[0:igs],valV2[0:igs])
  }
  acc_set_device_num(0,acc_device_nvidia);
#endif

  //printf(" done with dealloc in 2c integrals \n"); fflush(stdout);
  //auto_crash();

  //for (int i=0;i<iN;i++) delete [] valS1[i];
  //for (int i=0;i<iN;i++) delete [] valV2[i];
  delete [] valS1; delete [] valV2;

  delete [] grid;
  delete [] wt;

  return;
}

void compute_2c_ps(bool do_overlap, bool do_yukawa, double gamma, int natoms, int* atno, double* coords, vector<vector<double> > &basis, int quad_order, int nmu, int nnu, int nphi, double* A, int prl)
{
  if (do_overlap && prl>1) { printf("\n WARNING: testing do_overlap in compute_2c_ps \n"); }

  int ngpu = 0;
 #if USE_ACC
  ngpu = acc_get_num_devices(acc_device_nvidia);
 #endif
  //if (nomp>1) { printf(" WARNING: cannot run compute_2c_ps in parallel \n"); nomp = 1; }

  bool dy = do_yukawa;

  int N = basis.size();
  int N2 = N*N;

  if (N<1) { printf(" ERROR: cannot compute 2c integrals, no RI basis functions \n"); exit(-1); }

  int qos = quad_order*quad_order*quad_order;
  int gs = 1;
  int gs6 = 1;

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

  double gpumem = (double)acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  double togb = pow(1024.,-3.);

  vector<vector<int> > n2ip;
  int imaxN = get_imax_n2ip(10000,natoms,N,basis,n2ip);
  if (prl>1) printf("   imaxN: %2i \n",imaxN);

  int nbatch = 1;
  int nbatch_max = 24;
  bool passed_mem_check = 0;
  for (int nb=1;nb<nbatch_max;nb++)
  if (nmu%nb==0)
  {
    gs = nmu*nnu*nphi*qos/nb;
    gs6 = 6*gs;

    int Nmax = imaxN;
   //unclear how this total is arrived at (i.e., the extra 2.*gs6 term)
    double mem0 = gs*7. + 2.*gs6 + 1.*nmu*nnu*nphi/nb;
    double mem1 = 8.*(gs*2.*Nmax + 1.*mem0 + 1.*N2 + 2.*gs6 + 2.*gs);

    if (mem1<gpumem)
    {
      nbatch = nb;
      passed_mem_check = 1;
      break;
    }
  }

  if (!passed_mem_check) { printf("\n ERROR: could not find a batch size where nb divides nmu \n"); exit(-1); }

  if (prl>-1) { if (do_yukawa) printf("  beginning compute_2c_ps (Yukawa. gamma: %5.3f  nbatch: %2i) \n",gamma,nbatch); else printf("  beginning compute_2c_ps (ngpu: %i  nbatch: %2i) \n",ngpu,nbatch); }

  double gpumem_gb = gpumem*togb;
  printf("   gpu memory available: %6.3f GB \n",gpumem_gb);

  double* grid = new double[gs6];
  double* wt = new double[gs];

 //intermediate storage
  int iN = imaxN;
  int igs = iN*gs;
  double* valS1 = new double[igs];
  double* valV2 = new double[igs];
  //double** valS1 = new double*[iN]; for (int i=0;i<iN;i++) valS1[i] = new double[gs];
  //double** valV2 = new double*[iN]; for (int i=0;i<iN;i++) valV2[i] = new double[gs];

 #if USE_ACC
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc enter data create(A[0:N2]) //async(n+1)

    #pragma acc enter data create(grid[0:gs6],wt[0:gs]) //async(n+1)
    //#pragma acc enter data create(valS1[0:iN][0:gs],valV2[0:iN][0:gs])
    #pragma acc enter data create(valS1[0:igs],valV2[0:igs]) //async(n+1)

    acc_assign(N2,A,0.);

    //acc_wait_all();
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

 //okay up to here
 //allocate dealloc isn't the problem
  //return;

  double gpumem_2 = 1.*acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  printf("   after alloc, gpu memory available: %6.3f GB \n",gpumem_2*togb);

 #if OMP_PARA
 //if run in serial, this is okay
  //#pragma omp parallel for schedule(dynamic,1) num_threads(ngpu)
 #endif
  for (int m=0;m<natoms;m++)
  {
   #if OMP_PARA
    int tid = omp_get_thread_num();
   #else
    int tid = m%ngpu;
   #endif
    acc_set_device_num(tid,acc_device_nvidia);

    double Z1 = (double)atno[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];
    double coordn[6];
    coordn[0] = 0.; coordn[1] = 0.; coordn[2] = 0.;

    for (int wb=0;wb<nbatch;wb++)
    {
      generate_ps_quad_grid(tid,cfn,wb,nbatch,Z1,1,coordn,quad_order,quad_order,nmu,nnu,nphi,grid,wt);

      int sp1 = 0; int sp2 = 0;
      int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
      int s3 = n2ip[m][sp2]; int s4 = n2ip[m][sp2+1];
      //printf("    m: %i  s1/2: %3i %3i s3/4: %3i %3i \n",m,s1,s2,s3,s4);

      if (dy)
      {
       #pragma acc parallel loop collapse(2) present(valS1[0:igs],wt[0:gs]) async(tid+1)
        for (int ii1=0;ii1<s2-s1;ii1++)
        for (int j=0;j<gs;j++)
          valS1[ii1*gs+j] = wt[j];

       #pragma acc parallel loop collapse(2) present(valV2[0:igs]) async(tid+1)
        for (int ii2=0;ii2<s4-s3;ii2++)
        for (int j=0;j<gs;j++)
          valV2[ii2*gs+j] = 0.;
      }
      else
      {
       #pragma acc parallel loop collapse(2) present(valS1[0:igs]) async(tid+1)
        for (int ii1=0;ii1<s2-s1;ii1++)
        for (int j=0;j<gs;j++)
          valS1[ii1*gs+j] = 1.;

       #pragma acc parallel loop collapse(2) present(valV2[0:gs],wt[0:gs]) async(tid+1)
        for (int ii2=0;ii2<s4-s3;ii2++)
        for (int j=0;j<gs;j++)
          valV2[ii2*gs+j] = wt[j];
      }

     //first compute single atom ints
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_shd(tid,gs,grid,&valS1[ii1*gs],n1,l1,m1,zeta1);
      }

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        if (do_overlap)
        {
          eval_shd(tid,gs,grid,&valV2[ii2*gs],n2,l2,m2,zeta2);
        }
        else
        {
         //V
          if (dy)
            eval_inr_yukawa(gs,grid,&valV2[ii2*gs],n2,l2,zeta2,gamma);
          else
            eval_inr_r12(tid,gs,grid,&valV2[ii2*gs],n2,l2,zeta2);
          eval_sh_3rd(tid,gs,grid,&valV2[ii2*gs],n2,l2,m2);
        }
      } //loop i2 evaluate

      reduce_2c1(tid,s1,s2,s3,s4,gs,valS1,valV2,iN,N,A);
    } //sp partition over m


   //two-atom ints
    for (int n=m+1;n<natoms;n++)
    {
      double Z2 = (double)atno[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;
      coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;

      for (int wb=0;wb<nbatch;wb++)
      {
        generate_ps_quad_grid(tid,cfn,wb,nbatch,0.,2,coordn,quad_order,quad_order,nmu,nnu,nphi,grid,wt);
        add_r1_to_grid(tid,gs,grid,0.,0.,0.);

        int sp1 = 0; int sp2 = 0;
        int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
        int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2+1];
        //printf("    mn: %i %i s1/2: %3i %3i s3/4: %3i %3i \n",m,n,s1,s2,s3,s4);

        if (dy)
        {
         #pragma acc parallel loop collapse(2) present(valS1[0:igs],wt[0:gs]) async(tid+1)
          for (int ii1=0;ii1<s2-s1;ii1++)
          for (int j=0;j<gs;j++)
            valS1[ii1*gs+j] = wt[j];
         #pragma acc parallel loop collapse(2) present(valV2[0:igs]) async(tid+1)
          for (int ii1=0;ii1<s4-s3;ii1++)
          for (int j=0;j<gs;j++)
            valV2[ii1*gs+j] = 0.;
        }
        else
        {
         #pragma acc parallel loop collapse(2) present(valS1[0:igs]) async(tid+1)
          for (int ii1=0;ii1<s2-s1;ii1++)
          for (int j=0;j<gs;j++)
            valS1[ii1*gs+j] = 1.;
         #pragma acc parallel loop collapse(2) present(valV2[0:igs],wt[0:gs]) async(tid+1)
          for (int ii1=0;ii1<s4-s3;ii1++)
          for (int j=0;j<gs;j++)
            valV2[ii1*gs+j] = wt[j];
        }

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<double> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

         //S
          eval_shd(tid,gs,grid,&valS1[ii1*gs],n1,l1,m1,zeta1);
        }

       //second center
        recenter_grid_zero(tid,gs,grid,-A12,-B12,-C12);

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

          if (do_overlap)
          {
            eval_shd(tid,gs,grid,&valV2[ii2*gs],n2,l2,m2,zeta2);
          }
          else
          {
           //V
            if (dy)
              eval_inr_yukawa(gs,grid,&valV2[ii2*gs],n2,l2,zeta2,gamma);
            else
              eval_inr_r12(tid,gs,grid,&valV2[ii2*gs],n2,l2,zeta2);
            eval_sh_3rd(tid,gs,grid,&valV2[ii2*gs],n2,l2,m2);
          }
        }

        reduce_2c1(tid,s1,s2,s3,s4,gs,valS1,valV2,iN,N,A);

      } //sp partition over n
    } //loop n over second atom

   #if OMP_PARA
    #pragma acc wait
   #endif

  } //loop m over natoms

  #pragma omp barrier
 #if 0
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);
    //#pragma acc wait
    acc_wait_all();
  }
 #endif
  acc_set_device_num(0,acc_device_nvidia);

  if (ngpu>1)
  {
   //gather parallel
    double At[N2];
    for (int j=0;j<N2;j++) At[j] = 0.;

    for (int n=0;n<ngpu;n++)
    {
      acc_set_device_num(n,acc_device_nvidia);

      #pragma acc update self(A[0:N2]) //async(n+1)

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
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc exit data delete(A[0:N2]) //async(n+1)
    #pragma acc exit data delete(grid[0:gs6],wt[0:gs]) //async(n+1)
    //#pragma acc exit data delete(valS1[0:iN][0:gs],valV2[0:iN][0:gs])
    #pragma acc exit data delete(valS1[0:igs],valV2[0:igs]) //async(n+1)

    //acc_wait_all();
  }
  acc_set_device_num(0,acc_device_nvidia);
#endif

  //printf(" done with dealloc in 2c integrals \n"); fflush(stdout);
  //auto_crash();

  //for (int i=0;i<iN;i++) delete [] valS1[i];
  //for (int i=0;i<iN;i++) delete [] valV2[i];
  delete [] valS1; delete [] valV2;

  delete [] grid;
  delete [] wt;

  return;
}

void init_s12v3(int tid, bool dy, int s1, int s2, int s3, int s4, int s5, int s6, int iN, int iNa, int gs, double* val1, double* val2, double* val3, double* wt)
{
  int igs = iN*gs;
  int iags = iNa*gs;
  #pragma acc parallel loop collapse(2) present(val2[0:igs]) async(tid+1)
  for (int ii2=0;ii2<s4-s3;ii2++)
  {
    for (int j=0;j<gs;j++)
      val2[ii2*gs+j] = 1.;
  }

  if (dy)
  {
    #pragma acc parallel loop collapse(2) present(val1[0:igs],wt[0:gs]) async(tid+1)
    for (int ii1=0;ii1<s2-s1;ii1++)
    for (int j=0;j<gs;j++)
      val1[ii1*gs+j] = wt[j];

    if (val3!=NULL)
    #pragma acc parallel loop collapse(2) present(val3[0:iags]) async(tid+1)
    for (int ii3=0;ii3<s6-s5;ii3++)
    for (int j=0;j<gs;j++)
      val3[ii3*gs+j] = 0.;
  }
  else
  {
    #pragma acc parallel loop collapse(2) present(val1[0:igs]) async(tid+1)
    for (int ii1=0;ii1<s2-s1;ii1++)
    for (int j=0;j<gs;j++)
      val1[ii1*gs+j] = 1.;

    if (val3!=NULL)
    #pragma acc parallel loop collapse(2) present(val3[0:iags],wt[0:gs]) async(tid+1)
    for (int ii3=0;ii3<s6-s5;ii3++)
    for (int j=0;j<gs;j++)
      val3[ii3*gs+j] = wt[j];
  }

  return;
}

void eval_s12v3(int tid, bool dol, bool dy, double gamma, int s1, int s2, int s3, int s4, int s5, int s6, int gs, double* grid, vector<vector<double> >& basis, vector<vector<double> >& basis_aux, double* val1, double* val2, double* val3)
{
 //single-center evaluations

  for (int i1=s1;i1<s2;i1++)
  {
    int ii1 = i1-s1;

    vector<double> basis1 = basis[i1];
    int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

    eval_shd(tid,gs,grid,&val1[ii1*gs],n1,l1,m1,zeta1);
  }

  for (int i2=s3;i2<s4;i2++)
  {
    int ii2 = i2-s3;

    vector<double> basis2 = basis[i2];
    int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

    eval_shd(tid,gs,grid,&val2[ii2*gs],n2,l2,m2,zeta2);
  }

  for (int i3=s5;i3<s6;i3++)
  {
    int ii3 = i3-s5;
    vector<double> basis3 = basis_aux[i3];
    int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

    if (dol)
      eval_shd(tid,gs,grid,&val3[ii3*gs],n3,l3,m3,zeta3);
    else
    {
     //V
      if (dy)
        eval_inr_yukawa(gs,grid,&val3[ii3*gs],n3,l3,zeta3,gamma);
      else
        eval_inr_r12(tid,gs,grid,&val3[ii3*gs],n3,l3,zeta3);
      eval_sh_3rd(tid,gs,grid,&val3[ii3*gs],n3,l3,m3);
    }
  }

  return;
}

void eval_s12v3_2(int tid, bool dol, bool dy, double gamma, int s1, int s2, int s3, int s4, int s5, int s6, int gs, double* grid, vector<vector<double> >& basis, vector<vector<double> >& basis_aux, double* val1, double* val2, double* val3, int type, double A12, double B12, double C12, double A13, double B13, double C13)
{
 //multi-center evaluations

  for (int i1=s1;i1<s2;i1++)
  {
    int ii1 = i1-s1;

    vector<double> basis1 = basis[i1];
    int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

    eval_shd(tid,gs,grid,&val1[ii1*gs],n1,l1,m1,zeta1);
  }

  if (type==1 || type==3) //i2/i3 on second atom
    recenter_grid_zero(tid,gs,grid,-A12,-B12,-C12);
  if (type==4) //i2/i3 on atom 3
    recenter_grid_zero(tid,gs,grid,-A13,-B13,-C13);

  for (int i2=s3;i2<s4;i2++)
  {
    int ii2 = i2-s3;

    vector<double> basis2 = basis[i2];
    int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

    eval_shd(tid,gs,grid,&val2[ii2*gs],n2,l2,m2,zeta2);
  }

  if (type==2) //i3 only on second atom
    recenter_grid_zero(tid,gs,grid,-A12,-B12,-C12);

  if (type==3)
    recenter_grid_zero(tid,gs,grid,A12-A13,B12-B13,C12-C13);

  if (val3!=NULL)
  for (int i3=s5;i3<s6;i3++)
  {
    int ii3 = i3-s5;
    vector<double> basis3 = basis_aux[i3];
    int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

    if (dol)
      eval_shd(tid,gs,grid,&val3[ii3*gs],n3,l3,m3,zeta3);
    else
    {
     //V
      if (dy)
        eval_inr_yukawa(gs,grid,&val3[ii3*gs],n3,l3,zeta3,gamma);
      else
        eval_inr_r12(tid,gs,grid,&val3[ii3*gs],n3,l3,zeta3);
      eval_sh_3rd(tid,gs,grid,&val3[ii3*gs],n3,l3,m3);
    }
  }

 //return grid to original position
  if (type==1 || type==2)
    recenter_grid_zero(tid,gs,grid,A12,B12,C12);
  if (type==3) //type==4 does not reset
    recenter_grid_zero(tid,gs,grid,A13,B13,C13);

  return;
}

void eval_p12(int tid, int s1, int s2, int s3, int s4, int gs, double* grid, vector<vector<double> >& basis, double* val1, double* val2, double A12, double B12, double C12, double A13, double B13, double C13)
{
  //if (tid>-1) { printf("\n ERROR eval_p12 not ready for tid \n"); exit(-1); }

  int gs3 = gs*3;
  for (int i1=s1;i1<s2;i1++)
  {
    int ii1 = i1-s1;

    vector<double> basis1 = basis[i1];
    int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

    eval_pd(tid,gs,grid,&val1[ii1*gs3],n1,l1,m1,zeta1);
  }

  recenter_grid_zero(tid,gs,grid,-A12,-B12,-C12);

  for (int i2=s3;i2<s4;i2++)
  {
    int ii2 = i2-s3;

    vector<double> basis2 = basis[i2];
    int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

    eval_pd(tid,gs,grid,&val2[ii2*gs3],n2,l2,m2,zeta2);
  }

  recenter_grid_zero(tid,gs,grid,A12,B12,C12);

  return;
}

void compute_pVp_3c_ps(int natoms, int* atno, double* coords, vector<vector<double> > &basis, int quad_order, int quad_r_order, int nsplit, int nmu, int nnu, int nphi, double* pVp, int prl)
{
  if (prl>-1) printf("  beginning compute_pVp_3c_ps \n");

  int ngpu = 0;
 #if USE_ACC
  ngpu = acc_get_num_devices(acc_device_nvidia);
 #endif

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
  if (Nmax<=0) { printf("\n ERROR: couldn't calculate gpu memory requirements \n"); exit(-1); }

  vector<vector<int> > n2ip;
  int imaxN = get_imax_n2ip(Nmax,natoms,N,basis,n2ip);
  if (prl<1) printf("   imaxN: %2i \n",imaxN);

  double gpumem_gb = gpumem/1024./1024./1024.;
  printf("   gpu memory available: %6.3f GB \n",gpumem_gb);

 //intermediate storage
  iN = imaxN;
  int igs3 = iN*gs3;
  double* valS1 = new double[igs3];
  double* valS2 = new double[igs3];
  //double** valS1 = new double*[iN];  for (int i=0;i<iN;i++)  valS1[i] = new double[gs3];
  //double** valS2 = new double*[iN];  for (int i=0;i<iN;i++)  valS2[i] = new double[gs3];
  double* valt = new double[gsh];
  double* pVpp = new double[N2];

 #if USE_ACC
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc enter data create(pVpp[0:N2])

    #pragma acc enter data create(grid[0:gs6],wt[0:gsh])
    //#pragma acc enter data create(valS1[0:iN][0:gs3],valS2[0:iN][0:gs3])
    #pragma acc enter data create(valS1[0:igs3],valS2[0:igs3])
    #pragma acc enter data create(valt[0:gsh])

    acc_assign(N2,pVpp,0.);
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  double gpumem_2 = 1.*acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  printf("   after alloc, gpu memory available: %6.3f GB \n",gpumem_2*togb);

 //3c part of pVp integrals
 #if OMP_PARA
 #pragma omp parallel for schedule(dynamic,1) num_threads(ngpu)
 #endif
  for (int m=0;m<natoms;m++)
  {
   #if OMP_PARA
    int tid = omp_get_thread_num();
   #else
    int tid = m%ngpu;
   #endif
    acc_set_device_num(tid,acc_device_nvidia);

    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];
    double coordn[9];
    coordn[0] = coordn[1] = coordn[2] = 0.;

   //three-atom case
    for (int n=m+1;n<natoms;n++)
    {
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

          generate_ps_quad_grid_3c_refine(tid,ztm1,ztm2,nsplit,3,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
          add_r1_to_grid(tid,gsh,grid,0.,0.,0.);

          init_s12v3(tid,0,s1,s2,s3,s4,0,0,iN,0,gs3,  valS1,valS2,NULL,NULL);
          eval_p12  (tid,s1,s2,s3,s4,gsh,grid,basis,valS1,valS2,A12,B12,C12,A13,B13,C13);

          recenter_grid_zero(tid,gsh,grid,-A13,-B13,-C13);
          reduce_3cenp(tid,Z3,s1,s2,s3,s4,N,iN,gsh,grid,valS1,valS2,valt,wt,pVpp);
         //////////////////////////////////////////////////////////////////////////////////////////////////////

          acc_wait_all();
        } //loop p over third atom

      }//loop sp12

    } //loop n over second atom

   #if OMP_PARA
    acc_wait_all();
    #pragma acc wait
   #endif

  } //loop m over natoms

  #pragma omp barrier
 #if 0
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);
    //#pragma acc wait
    acc_wait_all();
  }
 #endif
  acc_set_device_num(0,acc_device_nvidia);


  double norm1[N];
  for (int i=0;i<N;i++)
    norm1[i] = basis[i][4];
  #pragma acc enter data copyin(norm1[0:N])

  if (ngpu>1)
  {
    double pVpt[N2];
    for (int j=0;j<N2;j++) pVpt[j] = 0.;

    for (int n=0;n<ngpu;n++)
    {
      acc_set_device_num(n,acc_device_nvidia);

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
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc exit data delete(pVpp[0:N2])
    #pragma acc exit data delete(norm1[0:N])
    #pragma acc exit data delete(grid[0:gs6],wt[0:gsh])
    //#pragma acc exit data delete(valS1[0:iN][0:gs3],valS2[0:iN][0:gs3])
    #pragma acc exit data delete(valS1[0:igs3],valS2[0:igs3])
    #pragma acc exit data delete(valt[0:gsh])
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  //printf(" done with dealloc in 3c integrals \n"); fflush(stdout);
  //auto_crash();

  delete [] n2i;

  //for (int i=0;i<iN;i++) delete [] valS1[i];
  //for (int i=0;i<iN;i++) delete [] valS2[i];
  delete [] valS1; delete [] valS2;
  delete [] valt;
  delete [] pVpp;

  delete [] grid;
  delete [] wt;

  return;
}

void compute_3c_ps(bool do_overlap, bool do_yukawa, double gamma, int nbatch_min, int natoms, int* atno, double* coords, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, int quad_order, int quad_r_order, int nsplit, int nmu, int nnu, int nphi, double* En, double* C, int prl)
{
  if (do_overlap && prl>1) { printf("\n WARNING: testing do_overlap in compute_3c_ps \n"); }

  int ngpu = 0;
 #if USE_ACC
  ngpu = acc_get_num_devices(acc_device_nvidia);
 #endif

  double cfn_en = 1.0;

  bool dy = do_yukawa;
  bool dol = do_overlap;
  if (En!=NULL && dy)
  { printf("\n ERROR: cannot run En with Yukawa in compute_3c_ps \n"); exit(-1); }

  int N = basis.size();
  int N2 = N*N;
  int Naux = basis_aux.size();

  if (N<1) { printf(" ERROR: cannot compute 3c integrals, no primary basis functions \n"); exit(-1); }
  if (Naux<1) { printf(" ERROR: cannot compute 3c integrals, no RI basis functions \n"); exit(-1); }

  double atnod[natoms];
  get_atnod(natoms,basis,atnod);

  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);

  vector<vector<int> > n2aip;
  int iNa = get_imax_n2ip(10000,natoms,Naux,basis_aux,n2aip);
  if (prl>1) printf("   iNa: %2i \n",iNa);

  int* na2i = new int[natoms]; //needed for copy_symm
  get_imax_n2i(natoms,Naux,basis_aux,na2i);

  int N2a = N2*Naux;
  int N2b = N2*iNa; //reduce C storage on GPU

  double gpumem = 1.*acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  double togb = pow(1024.,-3.);

  int qoh = quad_r_order; //refined grid region
  if (quad_order!=quad_r_order)
  {
    printf("  WARNING: quad_orders mismatached in compute_3c_ps. %2i %2i \n",quad_order,quad_r_order);
  }

  int qos = quad_order*quad_order*quad_order;
  int qosh = qoh*qoh*qoh;
  if (natoms<3) nsplit = 1;
  int nsg = nsplit*nsplit*nsplit;

 //batching both grids
  int gs = 1;
  int gsh = 1;
  int gs6 = 1;

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

  int nbatch = 1;

  //sets a minimum amount of batching
  //int nbatch_read = read_int("NBATCH");
  if (nbatch_min>1) nbatch = nbatch_min;

  int nbatch_max = 24;
  bool passed_mem_check = 0;
  for (int nb=nbatch;nb<nbatch_max;nb++)
 //maybe this should try dividing qos
  if (nmu%nb==0)
  {
    //printf("   nb loop: %2i \n",nb);
    gs =(nmu*nnu*nphi)*qos/nb; //2-atom grid size
    gsh =((nmu*nnu*nphi-8)*qos+8*nsg*qosh)/nb; //3-atom grid size
    gs6 = 6*gsh;

    int Nmax = iNa;
    double mem0 = gsh*iN*2. + gsh*7. + 1.*nmu*nnu*nphi + 1.*gs6 + 1.*N2 + N2b;
    double mem1 = 8.*(2.*gsh*iN + mem0 + 3.*gs6 + 2.*gsh + 1.*Nmax*gsh);

    //bool check1 = nmu%nb==0;
    //bool check2 = ((nmu*nnu*nphi-8)*qos + 8*nsg*qosh)%nb==0;
    if (mem1<gpumem)
    {
      nbatch = nb;
      passed_mem_check = 1;
      break;
    }
  }

  if (!passed_mem_check) { printf("\n ERROR: could not find a batch size where nb divides nmu \n"); exit(-1); }
  if (prl>-1) { if (do_yukawa) printf("  beginning compute_3c_ps (Yukawa. gamma: %5.3f  nbatch: %2i) \n",gamma,nbatch); else printf("  beginning compute_3c_ps (ngpu: %i  nbatch: %2i) \n",ngpu,nbatch); }

  double* grid = new double[gs6];
  double* wt = new double[gsh];

  double gsxvalsv = 8.*(gsh*iN*2. + gsh*iNa + gsh*7. + 1.*nmu*nnu*nphi); //vals+grid/wt+gridm
  double gsxvalsv_gb = gsxvalsv*togb;
  //printf("   estimated memory needed (3c): %6.3f GB \n",gsxvalsv_gb);

  double gpumem_gb = gpumem*togb;
  printf("   gpu memory available: %6.3f GB \n",gpumem_gb);

  if (gsxvalsv>gpumem) { printf("\n WARNING: probably not enough memory to do 3c integrals \n"); }

 //intermediate storage
  int igsh = iN*gsh;
  int iagsh = iNa*gsh;
  double* valS1 = new double[igsh];
  double* valS2 = new double[igsh];
  double* valV3 = new double[iagsh];
  //double** valS1 = new double*[iN];  for (int i=0;i<iN;i++)  valS1[i] = new double[gsh];
  //double** valS2 = new double*[iN];  for (int i=0;i<iN;i++)  valS2[i] = new double[gsh];
  //double** valV3 = new double*[iNa]; for (int i=0;i<iNa;i++) valV3[i] = new double[gsh];
  double* valt = new double[gsh];
  double* Enp = new double[N2];
 //tmp storage on GPU
  double** Ctp = new double*[ngpu];
  for (int n=0;n<ngpu;n++)
    Ctp[n] = new double[N2b];

 //will increment this as calculation proceeds
  for (int j=0;j<N2a;j++) C[j] = 0.;

 #if USE_ACC
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    double* Cp = Ctp[n];
    #pragma acc enter data create(Enp[0:N2],Cp[0:N2b])

    #pragma acc enter data create(grid[0:gs6],wt[0:gsh])
    //#pragma acc enter data create(valS1[0:iN][0:gsh],valS2[0:iN][0:gsh],valV3[0:iNa][0:gsh])
    #pragma acc enter data create(valS1[0:igsh],valS2[0:igsh],valV3[0:iagsh])
    #pragma acc enter data create(valt[0:gsh])

    acc_assign(gs6,grid,0.);
    //acc_assign(gsh,grid,0.);
    acc_assign(N2,Enp,0.);
    acc_assign(N2b,Cp,0.);
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  double gpumem_2 = 1.*acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  printf("   after alloc, gpu memory available: %6.3f GB \n",gpumem_2*togb);

 //Coulomb 3c integrals and 3c part of En integrals
 #if OMP_PARA
  #pragma omp parallel for schedule(dynamic,1) num_threads(ngpu)
 #endif
  for (int m=0;m<natoms;m++)
  {
   #if OMP_PARA
    int tid = omp_get_thread_num();
   #else
    int tid = m%ngpu;
   #endif
    acc_set_device_num(tid,acc_device_nvidia);
    double* Cp = Ctp[tid];

    double Z1 = atnod[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];
    double coordn[9];
    coordn[0] = coordn[1] = coordn[2] = 0.;

   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];
    int s3 = s1; int s4 = s2;
    int sp = 0;
    int s5 = n2aip[m][sp]; int s6 = n2aip[m][sp+1];

    #pragma acc parallel loop present(Cp[0:N2b])
    for (int j=0;j<N2b;j++)
      Cp[j] = 0.;

    for (int wb=0;wb<nbatch;wb++)
    {
      generate_ps_quad_grid(tid,1.,wb,nbatch,Z1,1,coordn,quad_order,quad_order,nmu,nnu,nphi,grid,wt);
      add_r1_to_grid(tid,gs,grid,0.,0.,0.);

     //all basis on one atom
      init_s12v3(tid,dy,s1,s2,s3,s4,s5,s6,iN,iNa,gs,valS1,valS2,valV3,wt);
      eval_s12v3(tid,dol,dy,gamma,s1,s2,s3,s4,s5,s6,gs,grid,basis,basis_aux,valS1,valS2,valV3);

      //printf("    m: %i   s12: %2i %2i s56: %2i %2i \n",m,s1,s2,s5,s6);
      reduce_3c1b(tid,s5,s6,s1,s2,gs,valV3,valS1,valS2,N,Naux,iN,iNa,NULL,Cp);
      acc_wait_all();

    } //loop wb over batches

    #pragma acc update self(Cp[0:N2b])
    //acc_wait_all();

    for (int i1=s5;i1<s6;i1++)
    for (int i2=s1;i2<s2;i2++)
    for (int i3=s1;i3<s2;i3++)
    {
      int ii1 = i1-s5; int ii2 = i2-s1; int ii3 = i3-s1;
     #pragma omp atomic
      C[i1*N2+i2*N+i3] += Cp[ii1*N2+i2*N+i3];
    }

   //two-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      //double Z2 = atnod[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;
      coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;

      s3 = 0; if (n>0) s3 = n2i[n-1]; s4 = n2i[n];
      int sp = 0;
      int s5 = n2aip[n][sp]; int s6 = n2aip[n][sp+1];

      #pragma acc parallel loop present(Cp[0:N2b])
      for (int j=0;j<N2b;j++)
        Cp[j] = 0.;

      //printf("     mn: %i %i   s12: %2i %2i s34: %2i %2i s56: %2i %2i \n",m,n,s1,s2,s3,s4,s5,s6);
      for (int wb=0;wb<nbatch;wb++)
      {
        generate_ps_quad_grid(tid,1.,wb,nbatch,0.,2,coordn,quad_order,quad_order,nmu,nnu,nphi,grid,wt);
        add_r1_to_grid(tid,gs,grid,0.,0.,0.);

        //printf(" n: %i  s5/6: %3i %3i \n",n,s5,s6);

       //s1 on atom 1, s2v3 on atom 2
        init_s12v3  (tid,dy,s1,s2,s3,s4,s5,s6,iN,iNa,gs,              valS1,valS2,valV3,wt);
        eval_s12v3_2(tid,dol,dy,gamma,s1,s2,s3,s4,s5,s6,gs,grid,basis,basis_aux,valS1,valS2,valV3,1,A12,B12,C12,0.,0.,0.);

        reduce_3c1b(tid,s5,s6,s1,s2,s3,s4,gs,valV3,valS1,valS2,N,Naux,iN,iNa,NULL,Cp);
        acc_wait_all();

       //s12 on atom 1, v3 on atom 2
        int s3b = s1; int s4b = s2;
        //printf("     mn: %i %i   s12: %2i %2i s34: %2i %2i s56: %2i %2i \n",m,n,s1,s2,s3,s4,s5,s6);

        init_s12v3  (tid,dy,s1,s2,s3b,s4b,s5,s6,iN,iNa,gs,              valS1,valS2,valV3,wt);
        eval_s12v3_2(tid,dol,dy,gamma,s1,s2,s3b,s4b,s5,s6,gs,grid,basis,basis_aux,valS1,valS2,valV3,2,A12,B12,C12,0.,0.,0.);

        reduce_3c1b(tid,s5,s6,s1,s2,s3b,s4b,gs,valV3,valS1,valS2,N,Naux,iN,iNa,NULL,Cp);
        acc_wait_all();
      } //loop wb over nbatch

      #pragma acc update self(Cp[0:N2b])
      //acc_wait_all();

      for (int i1=s5;i1<s6;i1++)
      for (int i2=s1;i2<s2;i2++)
      for (int i3=s3;i3<s4;i3++)
      {
        int ii1 = i1-s5; int ii2 = i2-s1; int ii3 = i3-s3;
       #pragma omp atomic
        C[i1*N2+i2*N+i3] += Cp[ii1*N2+i2*N+i3];
      }

      for (int i1=s5;i1<s6;i1++)
      for (int i2=s1;i2<s2;i2++)
      for (int i3=s1;i3<s2;i3++)
      {
        int ii1 = i1-s5; int ii2 = i2-s1; int ii3 = i3-s3;
       #pragma omp atomic
        C[i1*N2+i2*N+i3] += Cp[ii1*N2+i2*N+i3];
      }

    } //loop n over second atom

   #if OMP_PARA
    #pragma acc wait
   #endif

  } //loop m over first atom

  //acc_wait_all();
  //#pragma omp barrier
  //#pragma acc wait

 #if OMP_PARA
 #pragma omp parallel for schedule(dynamic,1) num_threads(ngpu)
 #endif
  for (int m=0;m<natoms;m++)
  {
   #if OMP_PARA
    int tid = omp_get_thread_num();
   #else
    int tid = m%ngpu;
   #endif
    acc_set_device_num(tid,acc_device_nvidia);
    double* Cp = Ctp[tid];

    double Z1 = atnod[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];
    double coordn[9];
    coordn[0] = coordn[1] = coordn[2] = 0.;

   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];
    int s3 = s1; int s4 = s2;

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

        int sp = 0;
        int s5 = n2aip[p][sp]; int s6 = n2aip[p][sp+1];

        //if (natoms>3) printf("     mnp: %i %i %i   s12: %3i %3i  s34: %3i %3i  s56: %3i %3i \n",m,n,p,s1,s2,s3,s4,s5,s6);

        double ztm1 = 0.; int lm1 = 0; double ztm2 = 0.; int lm2 = 0;
        get_ztm_lm(s1,s2,basis,ztm1,lm1);
        get_ztm_lm(s3,s4,basis,ztm2,lm2);
        //printf("    3c ztm/lm: %8.5f %i - %8.5f %i \n",ztm1,lm1,ztm2,lm2);

        #pragma acc parallel loop present(Cp[0:N2b])
        for (int j=0;j<N2b;j++)
          Cp[j] = 0.;

        for (int wb=0;wb<nbatch;wb++)
        {
          generate_ps_quad_grid_3c_refine(tid,wb,nbatch,ztm1,ztm2,nsplit,3,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
          add_r1_to_grid(tid,gsh,grid,0.,0.,0.);

         //s1 on atom 1, s2 on atom 2, v3 on atom 3
          init_s12v3  (tid,dy,s1,s2,s3,s4,s5,s6,iN,iNa,gsh,              valS1,valS2,valV3,wt);
          eval_s12v3_2(tid,dol,dy,gamma,s1,s2,s3,s4,s5,s6,gsh,grid,basis,basis_aux,valS1,valS2,valV3,3,A12,B12,C12,A13,B13,C13);

          reduce_3c1b(tid,s5,s6,s1,s2,s3,s4,gsh,valV3,valS1,valS2,N,Naux,iN,iNa,NULL,Cp);
          acc_wait_all();
        }

        #pragma acc update self(Cp[0:N2b])

        for (int i1=s5;i1<s6;i1++)
        for (int i2=s1;i2<s2;i2++)
        for (int i3=s3;i3<s4;i3++)
        {
          int ii1 = i1-s5; int ii2 = i2-s1; int ii3 = i3-s3;
         #pragma omp atomic
          C[i1*N2+i2*N+i3] += Cp[ii1*N2+i2*N+i3];
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
            generate_ps_quad_grid_3c_refine(tid,cfn_en,wb,nbatch,ztm1,ztm2,nsplit,3,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
            add_r1_to_grid(tid,gsh,grid,0.,0.,0.);

            init_s12v3  (tid,0,s1,s2,s3,s4,0,0,iN,iNa,gsh,              valS1,valS2,NULL,wt);
            eval_s12v3_2(tid,0,0,0.,s1,s2,s3,s4,0,0,gsh,grid,basis,basis_aux,valS1,valS2,NULL,3,A12,B12,C12,A13,B13,C13);

            recenter_grid_zero(tid,gsh,grid,-A13,-B13,-C13);
            reduce_3cen(tid,Z3,s1,s2,s3,s4,N,iN,gsh,grid,valS1,valS2,valt,wt,Enp);
            acc_wait_all();
          }
         //////////////////////////////////////////////////////////////////////////////////////////////////////
        }

      } //loop p over third atom

    } //loop n over second atom

  } //loop m over natoms

  acc_wait_all();
  #pragma omp barrier
 #if 0
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);
    //#pragma acc wait
    acc_wait_all();
  }
 #endif
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
  //#pragma acc enter data copyin(norm1[0:Naux],norm2[0:N])

 #if 0
 //aggregate C elements
  for (int n=0;n<ngpu;n++)
  for (int j=0;j<N2a;j++)
    C[j] += Ct[n][j];
 #endif

  if (ngpu>1)
  {
    double Ent[N2];
    for (int j=0;j<N2;j++) Ent[j] = 0.;
    //for (int j=0;j<N2a;j++) C[j] = 0.;

    for (int n=0;n<ngpu;n++)
    {
      acc_set_device_num(n,acc_device_nvidia);

      //#pragma acc update self(Cp[0:N2a],Enp[0:N2])
      #pragma acc update self(Enp[0:N2])

      for (int j=0;j<N2;j++)
        Ent[j] += Enp[j];

      //for (int j=0;j<N2a;j++)
      //  C[j] += Cp[j];
    }
    acc_set_device_num(0,acc_device_nvidia);
    //#pragma acc update device(C[0:N2a])

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
      #pragma acc update self(Enp[0:N2])
     //#pragma acc parallel loop collapse(2) present(Enp[0:N2],norm2[0:N])
      for (int i=0;i<N;i++)
      for (int j=0;j<N;j++)
        Enp[i*N+j] *= norm2[i]*norm2[j];

      //#pragma acc parallel loop independent present(Enp[0:N2])
      for (int i=0;i<N;i++)
      //#pragma acc loop independent
      for (int j=0;j<i;j++)
        Enp[i*N+j] = Enp[j*N+i];

      //#pragma acc update self(Enp[0:N2])

      for (int i=0;i<N2;i++)
        En[i] += Enp[i];
    }

   //#pragma acc parallel loop present(C[0:N2a],Cp[0:N2a])
   // for (int j=0;j<N2a;j++)
   //   C[j] = Cp[j];
  }

 //#pragma acc parallel loop collapse(3) present(C[0:N2a],norm1[0:Naux],norm2[0:N])
  for (int i=0;i<Naux;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
  {
    double n123 = norm1[i]*norm2[j]*norm2[k];
    C[i*N2+j*N+k] *= n123;
  }
  //#pragma acc update self(C[0:N2a])

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
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    double* Cp = Ctp[n];
    #pragma acc exit data delete(Enp[0:N2],Cp[0:N2b])
    #pragma acc exit data delete(grid[0:gs6],wt[0:gsh])
    //#pragma acc exit data delete(valS1[0:iN][0:gsh],valS2[0:iN][0:gsh],valV3[0:iNa][0:gsh])
    #pragma acc exit data delete(valS1[0:igsh],valS2[0:igsh],valV3[0:iagsh])
    #pragma acc exit data delete(valt[0:gsh])
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  //#pragma acc exit data delete(norm1[0:Naux],norm2[0:N])

  //printf(" done with dealloc in 3c integrals \n"); fflush(stdout);
  //auto_crash();

  delete [] n2i;
  delete [] na2i;

  //for (int i=0;i<iN;i++) delete [] valS1[i];
  //for (int i=0;i<iN;i++) delete [] valS2[i];
  //for (int i=0;i<iNa;i++) delete [] valV3[i];
  delete [] valS1; delete [] valS2; delete [] valV3;
  delete [] valt;
  delete [] Enp;
  for (int n=0;n<ngpu;n++)
    delete [] Ctp[n];
  delete [] Ctp;
  //delete [] Cp;

  delete [] grid;
  delete [] wt;

  return;
}

void init_s14(int tid, int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8, int iN, int gs, double* val1, double* val2, double* val3, double* val4, double* wt)
{
  int igs = iN*gs;
  #pragma acc parallel loop collapse(2) present(val1[0:igs]) async(tid+1)
  for (int ii1=0;ii1<s2-s1;ii1++)
  {
    for (int j=0;j<gs;j++)
      val1[ii1*gs+j] = 1.;
  }

  #pragma acc parallel loop collapse(2) present(val2[0:igs]) async(tid+1)
  for (int ii2=0;ii2<s4-s3;ii2++)
  {
    for (int j=0;j<gs;j++)
      val2[ii2*gs+j] = 1.;
  }

  #pragma acc parallel loop collapse(2) present(val3[0:igs]) async(tid+1)
  for (int ii3=0;ii3<s6-s5;ii3++)
  {
    for (int j=0;j<gs;j++)
      val3[ii3*gs+j] = 1.;
  }

  #pragma acc parallel loop collapse(2) present(val4[0:igs],wt[0:gs]) async(tid+1)
  for (int ii4=0;ii4<s8-s7;ii4++)
  {
    for (int j=0;j<gs;j++)
      val4[ii4*gs+j] = wt[j];
  }

  return;
}

void eval_s14(int tid, int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8, int gs, double* grid, vector<vector<double> >& basis, double* val1, double* val2, double* val3, double* val4, int type, double A12, double B12, double C12)
{
  for (int i1=s1;i1<s2;i1++)
  {
    int ii1 = i1-s1;

    vector<double> basis1 = basis[i1];
    int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

    eval_shd(tid,gs,grid,&val1[ii1*gs],n1,l1,m1,zeta1);
  }

  if (type==1)
    recenter_grid_zero(tid,gs,grid,-A12,-B12,-C12);

  for (int i2=s3;i2<s4;i2++)
  {
    int ii2 = i2-s3;

    vector<double> basis2 = basis[i2];
    int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

    eval_shd(tid,gs,grid,&val2[ii2*gs],n2,l2,m2,zeta2);
  }

  if (type==2)
    recenter_grid_zero(tid,gs,grid,-A12,-B12,-C12);

  for (int i3=s5;i3<s6;i3++)
  {
    int ii3 = i3-s5;

    vector<double> basis3 = basis[i3];
    int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

    eval_shd(tid,gs,grid,&val3[ii3*gs],n3,l3,m3,zeta3);
  }

  if (type==3)
    recenter_grid_zero(tid,gs,grid,-A12,-B12,-C12);

  for (int i4=s7;i4<s8;i4++)
  {
    int ii4 = i4-s7;

    vector<double> basis4 = basis[i4];
    int n4 = basis4[0]; int l4 = basis4[1]; int m4 = basis4[2]; double zeta4 = basis4[3];

    eval_shd(tid,gs,grid,&val4[ii4*gs],n4,l4,m4,zeta4);
  }

  if (type>0)
    recenter_grid_zero(tid,gs,grid,A12,B12,C12);

  return;
}

void eval_s14b(int tid, int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8, int gs, double* grid, vector<vector<double> >& basis, double* val1, double* val2, double* val3, double* val4, int type, double A12, double B12, double C12, double A13, double B13, double C13)
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

    eval_shd(tid,gs,grid,&val1[ii1*gs],n1,l1,m1,zeta1);
  }

  if (type==2 || type==3)
    recenter_grid_zero(tid,gs,grid,-A12,-B12,-C12);
  if (type==4)
    recenter_grid_zero(tid,gs,grid,-A13,-B13,-C13);

  for (int i2=s3;i2<s4;i2++)
  {
    int ii2 = i2-s3;

    vector<double> basis2 = basis[i2];
    int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

    eval_shd(tid,gs,grid,&val2[ii2*gs],n2,l2,m2,zeta2);
  }

  if (type==1)
    recenter_grid_zero(tid,gs,grid,-A12,-B12,-C12);
  if (type==3)
    recenter_grid_zero(tid,gs,grid,A12-A13,B12-B13,C12-C13);
  if (type==4)
    recenter_grid_zero(tid,gs,grid,A13-A12,B13-B12,C13-C12);

  for (int i3=s5;i3<s6;i3++)
  {
    int ii3 = i3-s5;

    vector<double> basis3 = basis[i3];
    int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

    eval_shd(tid,gs,grid,&val3[ii3*gs],n3,l3,m3,zeta3);
  }

  if (type==1 || type==2)
    recenter_grid_zero(tid,gs,grid,A12-A13,B12-B13,C12-C13);

  for (int i4=s7;i4<s8;i4++)
  {
    int ii4 = i4-s7;

    vector<double> basis4 = basis[i4];
    int n4 = basis4[0]; int l4 = basis4[1]; int m4 = basis4[2]; double zeta4 = basis4[3];

    eval_shd(tid,gs,grid,&val4[ii4*gs],n4,l4,m4,zeta4);
  }

  if (type==1 || type==2 || type==3)
    recenter_grid_zero(tid,gs,grid,A13,B13,C13);
  if (type==4)
    recenter_grid_zero(tid,gs,grid,A12,B12,C12);

  return;
}

void eval_s14c(int tid, int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8, int gs, double* grid, vector<vector<double> >& basis, double* val1, double* val2, double* val3, double* val4, double A12, double B12, double C12, double A13, double B13, double C13, double A14, double B14, double C14)
{
  for (int i1=s1;i1<s2;i1++)
  {
    int ii1 = i1-s1;

    vector<double> basis1 = basis[i1];
    int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

    eval_shd(tid,gs,grid,&val1[ii1*gs],n1,l1,m1,zeta1);
  }

  recenter_grid_zero(tid,gs,grid,-A12,-B12,-C12);

  for (int i2=s3;i2<s4;i2++)
  {
    int ii2 = i2-s3;

    vector<double> basis2 = basis[i2];
    int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

    eval_shd(tid,gs,grid,&val2[ii2*gs],n2,l2,m2,zeta2);
  }

  recenter_grid_zero(tid,gs,grid,A12-A13,B12-B13,C12-C13);

  for (int i3=s5;i3<s6;i3++)
  {
    int ii3 = i3-s5;

    vector<double> basis3 = basis[i3];
    int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

    eval_shd(tid,gs,grid,&val3[ii3*gs],n3,l3,m3,zeta3);
  }

  recenter_grid_zero(tid,gs,grid,A13-A14,B13-B14,C13-C14);

  for (int i4=s7;i4<s8;i4++)
  {
    int ii4 = i4-s7;

    vector<double> basis4 = basis[i4];
    int n4 = basis4[0]; int l4 = basis4[1]; int m4 = basis4[2]; double zeta4 = basis4[3];

    eval_shd(tid,gs,grid,&val4[ii4*gs],n4,l4,m4,zeta4);
  }

  recenter_grid_zero(tid,gs,grid,A14,B14,C14);


  return;
}

void reduce_4c_ol1(int tid, int ii1, int ii2, int s5, int s6, int s7, int s8, int gs, double* valm, double* valn, double* val3, double* val4, int iN, int M, int M2, int M3, int M4, double* tmp)
{
 //integrate second pair of ftns
  int igs = iN*gs;

  if (s5!=s7)
  {
   #pragma acc parallel loop collapse(2) present(valm[0:gs],valn[0:gs],val3[0:igs],val4[0:igs],tmp[0:M4]) async(tid+1)
    for (int i3=s5;i3<s6;i3++)
    for (int i4=s7;i4<s8;i4++)
    {
      int ii3 = i3-s5;
      int ii4 = i4-s7;

      double* valp = &val3[ii3*gs];
      double* valq = &val4[ii4*gs];

      double v1 = 0.;
     #pragma acc loop reduction(+:v1)
      for (int j=0;j<gs;j++)
        v1 += valm[j]*valn[j]*valp[j]*valq[j];

      tmp[ii1*M3+ii2*M2+ii3*M+ii4] = v1;
    }
  }
  else
  {
   #pragma acc parallel loop present(valm[0:gs],valn[0:gs],val3[0:igs],val4[0:igs],tmp[0:M4]) async(tid+1)
    for (int i3=s5;i3<s6;i3++)
   #pragma acc loop
    for (int i4=s5;i4<=i3;i4++)
    {
      int ii3 = i3-s5;
      int ii4 = i4-s5;

      double* valp = &val3[ii3*gs];
      double* valq = &val4[ii4*gs];

      double v1 = 0.;
     #pragma acc loop reduction(+:v1)
      for (int j=0;j<gs;j++)
        v1 += valm[j]*valn[j]*valp[j]*valq[j];

      tmp[ii1*M3+ii2*M2+ii3*M+ii4] = v1;
      tmp[ii1*M3+ii2*M2+ii4*M+ii3] = v1;
    }
  }

  return;
}

//worth exploring acc optimization
void reduce_4c_ol1(int tid, int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8, int gs, double* val1, double* val2, double* val3, double* val4, int N, int iN, double* tmp, double* olp)
{
  int N2 = N*N;
  int N3 = N2*N;

  int igs = iN*gs;

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

    double* valm = &val1[ii1*gs];
    double* valn = &val2[ii2*gs];

    reduce_4c_ol1(tid,ii1,ii2,s5,s6,s7,s8,gs,valm,valn,val3,val4,iN,M,M2,M3,M4,tmp);

   #pragma acc serial present(tmp[0:M4]) async(tid+1)
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

    double* valm = &val1[ii1*gs];
    double* valn = &val2[ii2*gs];

    reduce_4c_ol1(tid,ii1,ii2,s5,s6,s7,s8,gs,valm,valn,val3,val4,iN,M,M2,M3,M4,tmp);
  }

  #pragma acc update self(tmp[0:M4]) async(tid+1)

 //CPMZ will limit parallelization?
  //#pragma acc wait(tid)
  acc_wait_all();

 //this last assignment may be problematic in parallel.
  //testing omp atomic statement
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  for (int i3=s5;i3<s6;i3++)
  for (int i4=s7;i4<s8;i4++)
  {
   #pragma omp atomic
    olp[i1*N3+i2*N2+i3*N+i4] += tmp[(i1-s1)*M3+(i2-s3)*M2+(i3-s5)*M+(i4-s7)];
  }

  return;
}

void reweight_core(int tid, const double beta, int natoms, int* atno, double* coords, int gs, double* grid, double* wt)
{
  int gs6 = 6*gs;

  printf("   reweight_core for %i atoms. Z: ",natoms);
  for (int n=0;n<natoms;n++)
    printf(" %i",atno[n]);
  printf("\n");

  for (int n=0;n<natoms;n++)
  {
    double alpha = 2.*atno[n];
    double A1 = coords[3*n+0]; double B1 = coords[3*n+1]; double C1 = coords[3*n+2];

   #pragma acc parallel loop present(grid[0:gs6],wt[0:gs]) async(tid+1)
    for (int j=0;j<gs;j++)
    {
      double x12 = grid[6*j+0]-A1;
      double y12 = grid[6*j+1]-B1;
      double z12 = grid[6*j+2]-C1;
      double r12 = x12*x12+y12*y12+z12*z12;

      wt[j] *= 1.+beta*exp(-alpha*r12);
    }
  }

  return;
}

void compute_4c_ol_ps(int nbatch_min, int natoms, int* atno, double* coords, vector<vector<double> > &basis, int quad_order, int quad_r_order, int nmu, int nnu, int nphi, double* ol, int prl)
{
  int ngpu = 0;
 #if USE_ACC
  ngpu = acc_get_num_devices(acc_device_nvidia);
 #endif

  double rw_core = read_float("CORE4C");
  if (rw_core<0.) rw_core = 0.;

  int N = basis.size();
  size_t N2 = N*N;
  size_t N3 = N2*N;
  size_t N4 = N3*N;

  vector<vector<int> > n2ip;
  int imaxN = get_imax_n2ip(10000,natoms,N,basis,n2ip);
  if (prl>1 || imaxN>80) printf("   imaxN: %2i \n",imaxN);

 //needed for copy_symm ftn
  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);

  double gpumem = (double)acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  double togb = pow(1024.,-3.);

  int qos = quad_order*quad_order*quad_order;
  int qoh = quad_r_order;
  int qosh = qoh*qoh*qoh;

  int gs = 1;
  int gsh = 1;
  int gshh = 1;
  int gs6 = 1;

  int M4 = iN*iN*iN*iN;

  int nbatch = 1;
  //sets a minimum amount of batching
  //int nbatch_read = read_int("NBATCH");
  if (nbatch_min>1) nbatch = nbatch_min;

  int nbatch_max = 24;
  bool passed_mem_check = 0;
  for (int nb=nbatch;nb<nbatch_max;nb++)
  if (nmu%nb==0)
  {
    gs = nmu*nnu*nphi*qos/nb;
    gsh = ((nmu*nnu*nphi-8)*qos+8*qosh)/nb;
    gshh = ((nmu*nnu*nphi-16)*qos+16*qosh)/nb;
    gs6 = 6*gshh;

    double mem1 = 8*(4.*iN*gshh + 1.*gs6 + 2.*gshh + nmu*nnu*nphi/nb + 1.*M4);

    if (mem1<gpumem)
    {
      nbatch = nb;
      passed_mem_check = 1;
      break;
    }
  }

  if (!passed_mem_check) { printf("\n ERROR: could not find a batch size where nb divides nmu \n"); exit(-1); }

  if (prl>-1) printf("  beginning compute_4c_ol_ps.  ngpu: %i  nbatch: %2i \n",ngpu,nbatch);
  //if (ngpu>1) printf("   TESTING parallel 4c_ol_ps \n");


 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

  double* grid = new double[gs6];
  double* wt = new double[gshh];

  printf("  gs: %8i  gpu mem total: %6.3f GB ngpu: %2i \n",gshh,gpumem*togb,ngpu);

 //intermediate storage
  int igshh = iN*gshh;
  double* valS1 = new double[igshh];
  double* valS2 = new double[igshh];
  double* valS3 = new double[igshh];
  double* valS4 = new double[igshh];
  //double** valS1 = new double*[iN]; for (int i=0;i<iN;i++) valS1[i] = new double[gshh];
  //double** valS2 = new double*[iN]; for (int i=0;i<iN;i++) valS2[i] = new double[gshh];
  //double** valS3 = new double*[iN]; for (int i=0;i<iN;i++) valS3[i] = new double[gshh];
  //double** valS4 = new double*[iN]; for (int i=0;i<iN;i++) valS4[i] = new double[gshh];
  double* valt = new double[gshh];
  double** tmpp = new double*[ngpu];
  for (int n=0;n<ngpu;n++)
    tmpp[n] = new double[M4]; //on gpu
  double* olp = new double[N4](); //not on gpu

 #if USE_ACC
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    double* tmp = tmpp[n];
    //#pragma acc enter data create(olp[0:N4])
    #pragma acc enter data create(tmp[0:M4])

    #pragma acc enter data create(grid[0:gs6],wt[0:gshh])
    //#pragma acc enter data create(valS1[0:iN][0:gshh],valS2[0:iN][0:gshh],valS3[0:iN][0:gshh],valS4[0:iN][0:gshh])
    #pragma acc enter data create(valS1[0:igshh],valS2[0:igshh],valS3[0:igshh],valS4[0:igshh])
    #pragma acc enter data create(valt[0:gshh])

    acc_assign(M4,tmp,0.);
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  double gpumem_2 = 1.*acc_get_property(0,acc_device_nvidia,acc_property_free_memory);
  printf("   after alloc, gpu memory available: %6.3f GB \n",gpumem_2*togb);

 #if OMP_PARA
 #pragma omp parallel for schedule(dynamic,1) num_threads(ngpu)
 #endif
  for (int m=0;m<natoms;m++)
  {
   #if OMP_PARA
    int tid = omp_get_thread_num();
   #else
    int tid = m%ngpu;
   #endif
    acc_set_device_num(tid,acc_device_nvidia);
    double* tmp = tmpp[tid];

    double Z1 = (double)atno[m];
    //double Z1 = atnod[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];
    double coordn[12];
    coordn[0] = coordn[1] = coordn[2] = 0.;
    int atnon[4]; atnon[0] = atno[m];

    for (int wb=0;wb<nbatch;wb++)
    {
      generate_ps_quad_grid(tid,1.,wb,nbatch,Z1,1,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
      add_r1_to_grid(tid,gs,grid,0.,0.,0.);
      if (rw_core>0.) reweight_core(tid,rw_core,1,atnon,coordn,gs,grid,wt);

      int sp1 = 0; int sp2 = 0; int sp3 = 0; int sp4 = 0;
      int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
      int s3 = n2ip[m][sp2]; int s4 = n2ip[m][sp2+1];
      int s5 = n2ip[m][sp3]; int s6 = n2ip[m][sp3+1];
      int s7 = n2ip[m][sp4]; int s8 = n2ip[m][sp4+1];

      //printf("    m: %i   s12: %2i %2i s34: %2i %2i s56: %2i %2i s78: %2i %2i \n",m,s1,s2,s3,s4,s5,s6,s7,s8);

      init_s14(tid,s1,s2,s3,s4,s5,s6,s7,s8,iN,gs,valS1,valS2,valS3,valS4,wt);
      eval_s14(tid,s1,s2,s3,s4,s5,s6,s7,s8,gs,grid,basis,valS1,valS2,valS3,valS4,0,0.,0.,0.);

      reduce_4c_ol1(tid,s1,s2,s3,s4,s5,s6,s7,s8,gs,valS1,valS2,valS3,valS4,N,iN,tmp,olp);
    }

   //two-atom ints
    for (int n=m+1;n<natoms;n++)
    {
      double Z2 = (double)atno[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;
      coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;
      atnon[1] = atno[n];

      for (int wb=0;wb<nbatch;wb++)
      {
        generate_ps_quad_grid(tid,1.,wb,nbatch,0.,2,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
        add_r1_to_grid(tid,gs,grid,0.,0.,0.);
        if (rw_core>0.) reweight_core(tid,rw_core,2,atnon,coordn,gs,grid,wt);

        int sp1 = 0; int sp2 = 0; int sp3 = 0; int sp4 = 0;
        int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
        int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2+1];
        int s5 = n2ip[n][sp3]; int s6 = n2ip[n][sp3+1];
        int s7 = n2ip[n][sp4]; int s8 = n2ip[n][sp4+1];

        //printf("     mn: %i %i   s12: %2i %2i s34: %2i %2i s56: %2i %2i s78: %2i %2i \n",m,n,s1,s2,s3,s4,s5,s6,s7,s8);

       //s1 on atom 1, s234 on atom 2
        init_s14(tid,s1,s2,s3,s4,s5,s6,s7,s8,iN,gs,valS1,valS2,valS3,valS4,wt);
        eval_s14(tid,s1,s2,s3,s4,s5,s6,s7,s8,gs,grid,basis,valS1,valS2,valS3,valS4,1,A12,B12,C12);

        reduce_4c_ol1(tid,s1,s2,s3,s4,s5,s6,s7,s8,gs,valS1,valS2,valS3,valS4,N,iN,tmp,olp);

        s1 = n2ip[m][sp1]; s2 = n2ip[m][sp1+1];
        s3 = n2ip[m][sp2]; s4 = n2ip[m][sp2+1];
        s5 = n2ip[n][sp3]; s6 = n2ip[n][sp3+1];
        s7 = n2ip[n][sp4]; s8 = n2ip[n][sp4+1];

        //printf("     mn: %i   s12: %2i %2i s34: %2i %2i s56: %2i %2i s78: %2i %2i \n",m,n,s1,s2,s3,s4,s5,s6,s7,s8);

       //s12 on atom 1, s34 on atom 2
        init_s14(tid,s1,s2,s3,s4,s5,s6,s7,s8,iN,gs,valS1,valS2,valS3,valS4,wt);
        eval_s14(tid,s1,s2,s3,s4,s5,s6,s7,s8,gs,grid,basis,valS1,valS2,valS3,valS4,2,A12,B12,C12);

        reduce_4c_ol1(tid,s1,s2,s3,s4,s5,s6,s7,s8,gs,valS1,valS2,valS3,valS4,N,iN,tmp,olp);

        s1 = n2ip[m][sp1]; s2 = n2ip[m][sp1+1];
        s3 = n2ip[m][sp2]; s4 = n2ip[m][sp2+1];
        s5 = n2ip[m][sp3]; s6 = n2ip[m][sp3+1];
        s7 = n2ip[n][sp4]; s8 = n2ip[n][sp4+1];

        //printf("     mn: %i %i   s12: %2i %2i s34: %2i %2i s56: %2i %2i s78: %2i %2i \n",m,n,s1,s2,s3,s4,s5,s6,s7,s8);

       //s123 on atom 1, s4 on atom 2
        init_s14(tid,s1,s2,s3,s4,s5,s6,s7,s8,iN,gs,valS1,valS2,valS3,valS4,wt);
        eval_s14(tid,s1,s2,s3,s4,s5,s6,s7,s8,gs,grid,basis,valS1,valS2,valS3,valS4,3,A12,B12,C12);

        reduce_4c_ol1(tid,s1,s2,s3,s4,s5,s6,s7,s8,gs,valS1,valS2,valS3,valS4,N,iN,tmp,olp);

      } //batches loop

    } //loop n over second atom

   //three-atom case
    for (int n=m+1;n<natoms;n++)
    {
      double Z2 = (double)atno[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;
      coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;
      atnon[1] = atno[n];

      for (int p=n+1;p<natoms;p++)
      if (p!=m)
      {
        double Z3 = (double)atno[p];
        double A3 = coords[3*p+0]; double B3 = coords[3*p+1]; double C3 = coords[3*p+2];
        double A13 = A3-A1; double B13 = B3-B1; double C13 = C3-C1;
        coordn[6] = A13; coordn[7] = B13; coordn[8] = C13;
        atnon[2] = atno[p];

        //printf("  mnp: %i %i %i   s12: %2i %2i  s34: %2i %2i  s56: %2i %2i \n",m,n,p,s1,s2,s3,s4,s5,s6);

        double ztm = 0.; int lm = 0;
        for (int sp4=0;sp4<n2ip[p].size()-1;sp4++)
        {
          int s7 = n2ip[p][sp4]; int s8 = n2ip[p][sp4+1];
          get_ztm_lm(s7,s8,basis,ztm,lm);
        }

        for (int wb=0;wb<nbatch;wb++)
        {
          generate_ps_quad_grid(tid,1.,wb,nbatch,0.,3,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
          add_r1_to_grid(tid,gsh,grid,0.,0.,0.);
          if (rw_core>0.) reweight_core(tid,rw_core,3,atnon,coordn,gs,grid,wt);

          int sp1 = 0; int sp2 = 0; int sp3 = 0; int sp4 = 0;

          int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
          int s3 = n2ip[m][sp2]; int s4 = n2ip[m][sp2+1];
          int s5 = n2ip[n][sp3]; int s6 = n2ip[n][sp3+1];
          int s7 = n2ip[p][sp4]; int s8 = n2ip[p][sp4+1];

         //s12 on atom 1, s3 on atom 2, s4 on atom 3 (type==1)
          init_s14(tid,s1,s2,s3,s4,s5,s6,s7,s8,iN,gsh,valS1,valS2,valS3,valS4,wt);
          eval_s14b(tid,s1,s2,s3,s4,s5,s6,s7,s8,gsh,grid,basis,valS1,valS2,valS3,valS4,1,A12,B12,C12,A13,B13,C13);

          reduce_4c_ol1(tid,s1,s2,s3,s4,s5,s6,s7,s8,gsh,valS1,valS2,valS3,valS4,N,iN,tmp,olp);

          s1 = n2ip[m][sp1]; s2 = n2ip[m][sp1+1];
          s3 = n2ip[n][sp2]; s4 = n2ip[n][sp2+1];
          s5 = n2ip[n][sp3]; s6 = n2ip[n][sp3+1];
          s7 = n2ip[p][sp4]; s8 = n2ip[p][sp4+1];

         //s1 on atom 1, s23 on atom 2, s4 on atom 3 (type==2)
          init_s14(tid,s1,s2,s3,s4,s5,s6,s7,s8,iN,gsh,valS1,valS2,valS3,valS4,wt);
          eval_s14b(tid,s1,s2,s3,s4,s5,s6,s7,s8,gsh,grid,basis,valS1,valS2,valS3,valS4,2,A12,B12,C12,A13,B13,C13);

          reduce_4c_ol1(tid,s1,s2,s3,s4,s5,s6,s7,s8,gsh,valS1,valS2,valS3,valS4,N,iN,tmp,olp);

          s1 = n2ip[m][sp1]; s2 = n2ip[m][sp1+1];
          s3 = n2ip[n][sp2]; s4 = n2ip[n][sp2+1];
          s5 = n2ip[p][sp3]; s6 = n2ip[p][sp3+1];
          s7 = n2ip[p][sp4]; s8 = n2ip[p][sp4+1];

         //s1 on atom 1, s2 on atom 2, s34 on atom 3 (type==3)
         // works but most basis ftns should be on atoms 1+2
          init_s14(tid,s1,s2,s3,s4,s5,s6,s7,s8,iN,gsh,valS1,valS2,valS3,valS4,wt);
          eval_s14b(tid,s1,s2,s3,s4,s5,s6,s7,s8,gsh,grid,basis,valS1,valS2,valS3,valS4,3,A12,B12,C12,A13,B13,C13);

          reduce_4c_ol1(tid,s1,s2,s3,s4,s5,s6,s7,s8,gsh,valS1,valS2,valS3,valS4,N,iN,tmp,olp);

        } //batches loop

      } //loop p over third atom

    } //loop n over second atom

   //four-atom case
    for (int n=m+1;n<natoms;n++)
    {
      double Z2 = (double)atno[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];
      double A12 = A2-A1; double B12 = B2-B1; double C12 = C2-C1;
      coordn[3] = A12; coordn[4] = B12; coordn[5] = C12;
      atnon[1] = atno[n];

      for (int p=n+1;p<natoms;p++)
      if (p!=m)
      for (int q=p+1;q<natoms;q++)
      if (q!=m && q!=n)
      {
        double Z3 = (double)atno[p];
        double A3 = coords[3*p+0]; double B3 = coords[3*p+1]; double C3 = coords[3*p+2];
        double A13 = A3-A1; double B13 = B3-B1; double C13 = C3-C1;
        coordn[6] = A13; coordn[7] = B13; coordn[8] = C13;
        atnon[2] = atno[p];

        double Z4 = (double)atno[q];
        double A4 = coords[3*q+0]; double B4 = coords[3*q+1]; double C4 = coords[3*q+2];
        double A14 = A4-A1; double B14 = B4-B1; double C14 = C4-C1;
        coordn[9] = A14; coordn[10] = B14; coordn[11] = C14;
        atnon[3] = atno[q];

       //CPMZ need to fix this up
        for (int wb=0;wb<nbatch;wb++)
        {
          generate_ps_quad_grid(tid,1.,wb,nbatch,0.,4,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
          add_r1_to_grid(tid,gshh,grid,0.,0.,0.);
          if (rw_core>0.) reweight_core(tid,rw_core,4,atnon,coordn,gs,grid,wt);

        //printf("  mnpq: %i %i %i %i   s12: %2i %2i  s34: %2i %2i  s56: %2i %2i  s78: %2i %2i \n",m,n,p,q,s1,s2,s3,s4,s5,s6,s7,s8);

          int sp1 = 0; int sp2 = 0; int sp3 = 0; int sp4 = 0;

          int s1 = n2ip[m][sp1]; int s2 = n2ip[m][sp1+1];
          int s3 = n2ip[n][sp2]; int s4 = n2ip[n][sp2+1];
          int s5 = n2ip[p][sp3]; int s6 = n2ip[p][sp3+1];
          int s7 = n2ip[q][sp4]; int s8 = n2ip[q][sp4+1];

         //s1 on atom 1, s2 on atom 2, s3 on atom 3, s4 on atom 4
          init_s14(tid,s1,s2,s3,s4,s5,s6,s7,s8,iN,gshh,valS1,valS2,valS3,valS4,wt);
          eval_s14c(tid,s1,s2,s3,s4,s5,s6,s7,s8,gshh,grid,basis,valS1,valS2,valS3,valS4,A12,B12,C12,A13,B13,C13,A14,B14,C14);

          reduce_4c_ol1(tid,s1,s2,s3,s4,s5,s6,s7,s8,gshh,valS1,valS2,valS3,valS4,N,iN,tmp,olp);

        } //batches loop

      } //loop pq over third and fourth atoms
    } //loop n over second atom

   #if OMP_PARA
    #pragma acc wait
   #endif

  } //loop m over natoms

  #pragma omp barrier
 #if 0
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);
    //#pragma acc wait
    acc_wait_all();
  }
 #endif
  acc_set_device_num(0,acc_device_nvidia);


  double norm[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];

  copy_symm_4c_ps_cpu(natoms,n2i,N,olp);

  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
  for (int l=0;l<N;l++)
  {
    double n1234 = norm[i]*norm[j]*norm[k]*norm[l];
    olp[i*N3+j*N2+k*N+l] *= n1234;
  }

  for (size_t j=0;j<N4;j++)
    ol[j] = olp[j];

  const double lt = 1.e-15;
  for (size_t j=0;j<N4;j++)
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
 //#pragma omp parallel for schedule(static,1) num_threads(ngpu)
  for (int n=0;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    double* tmp = tmpp[n];
    //#pragma acc exit data delete(olp[0:N4])
    #pragma acc exit data delete(tmp[0:M4])
    #pragma acc exit data delete(grid[0:gs6],wt[0:gshh])
    //#pragma acc exit data delete(valS1[0:iN][0:gshh],valS2[0:iN][0:gshh],valS3[0:iN][0:gshh],valS4[0:iN][0:gshh])
    #pragma acc exit data delete(valS1[0:igshh],valS2[0:igshh],valS3[0:igshh],valS4[0:igshh])
    #pragma acc exit data delete(valt[0:gshh])
  }
  acc_set_device_num(0,acc_device_nvidia);
 #endif

  delete [] n2i;

  //for (int i=0;i<iN;i++) delete [] valS1[i];
  //for (int i=0;i<iN;i++) delete [] valS2[i];
  //for (int i=0;i<iN;i++) delete [] valS3[i];
  //for (int i=0;i<iN;i++) delete [] valS4[i];
  delete [] valS1; delete [] valS2; delete [] valS3; delete [] valS4;
  delete [] valt;
  for (int n=0;n<ngpu;n++)
    delete [] tmpp[n];
  delete [] tmpp;
  delete [] olp;

  delete [] grid;
  delete [] wt;

  return;
}
