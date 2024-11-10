#include "cusp.h"

#define PSP_THRESH_WARNING 1.e-8
#define PSP_THRESH 1.e-10
#define DOUBLE_CUSP 1

using namespace std;

void print_rectangle_sm(int N1, int N2, double* S);


void compute_jellium_cusp(int order, int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* pb, int prl)
{
  if (natoms>1) { printf(" ERROR: cannot use jellium cusp with natoms>1 \n"); exit(-1); }

  double Rc = read_float("RC");
 //#include "jsetup.cpp"
  double Zg = read_float("ZG");
  double ztg = read_float("ZTG");
  int Ne = read_int("NEJ");
  bool ss_basis = 0; if (basis[0].size()>10) ss_basis = 1;
  if (!ss_basis) { printf("  WARNING: no cusp condition applied (jellium) \n"); return; }

 //must use either n=l+1 or n=l+2 for all basis ftns
  int cusp_type = 1;
  if (basis[0][1]>0 || basis[0][0]>2) { printf("\n ERROR: cusp requires first basis to be 1s or 2s \n"); exit(-1); }
  if (basis[0][0]==2)
    cusp_type = 2;

  //if (cusp_type>1) { printf("\n ERROR: cusp not ready beyond 1s \n"); exit(-1); }

  printf("  cusp w/Jellium (type: %i) \n",cusp_type);

  double f1 = -0.5 * Ne * pow(Rc,-3);
 //CPMZ testing
  f1 -= Zg*ztg;
  if (order>2)
  {
    if (cusp_type==1)
      f1 = 0.; //r2 term zero, suitable for 1s basis
    else
      f1 = -(1.*Ne)/order*pow(Rc,-(order+1)); //r^order term
  }
  else if (cusp_type>1)
    printf("\n WARNING: jellium order %i incompatible w/ cusp_type %i \n",order,cusp_type);

  int N = basis.size();

 //only one condition
  double* pb1 = new double[N]();
  {
    for (int i1=0;i1<N;i1++)
    {
      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; double zeta1 = basis1[3];

     //testing
      if (l1==0)
      {
        if (n1==1)
          pb1[i1] = f1 + pow(zeta1,4)/12.; //SS_V3 basis
        else if (n1==2)
          pb1[i1] = f1 + pow(zeta1,6)/240.; //SS_V3 basis
        else
        {
          printf("\n ERROR: cusp doesn't support jellium 3s basis in \n"); exit(-1);
        }
      }
      else
        pb1[i1] = 0.;
    }
  }

  for (int i=0;i<N;i++)
    pb[i] = pb1[i];

  delete [] pb1;

  return;
}

void compute_cusp(int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* pb, int prl)
{
  int jellium = read_int("JELLIUM");
  if (jellium>0)
    return compute_jellium_cusp(jellium,natoms,atno,coords,basis,pb,prl);

  int N = basis.size();
  int nN = natoms*N;

  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);

  double* pb1 = new double[nN]();
  double* pb2 = new double[nN]();

 #if DOUBLE_CUSP
  double* grid1 = new double[6]();
  double* tmp1 = new double[1]();
 #else
  float* grid1 = new float[6]();
  float* tmp1 = new float[1]();
 #endif

  #pragma acc enter data create(grid1[0:6],tmp1[0:1])

  int not_X[natoms];
  for (int n=0;n<natoms;n++) not_X[n] = 0;
  for (int n=0;n<natoms;n++)
  if (atno[n]>0)
    not_X[n] = 1;

  double* norm = new double[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];

  for (int n=0;n<natoms;n++)
  if (not_X[n])
  {
    int s1 = 0; if (n>0) s1 = n2i[n-1]; int s2 = n2i[n];

    double Z1 = (double)atno[n];
    double A1 = coords[3*n+0]; double B1 = coords[3*n+1]; double C1 = coords[3*n+2];

   //basis functions on all atoms, evaluated at atom n
    for (int i1=0;i1<N;i1++)
    {
      int ind1 = n*N+i1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

      double A2 = basis1[5]; double B2 = basis1[6]; double C2 = basis1[7];
      double A12 = A1-A2; double B12 = B1-B2; double C12 = C1-C2;

      double r1 = sqrt(A12*A12+B12*B12+C12*C12);
     #pragma acc serial present(grid1[0:6],tmp1[0:1])
      {
        grid1[0] = A12; grid1[1] = B12; grid1[2] = C12;
        grid1[3] = grid1[4] = grid1[5] = r1;

        tmp1[0] = 1.;
      }
      //printf("  i1: %i nlm: %i %i %i zeta: %5.3f  r: %5.3f \n",i1,n1,l1,m1,zeta1,r1);

     //not much to compute
     #if DOUBLE_CUSP
      eval_shd(i1,1,grid1,tmp1,n1,l1,m1,zeta1);
     #else
      eval_sh(i1,1,grid1,tmp1,n1,l1,m1,zeta1);
     #endif
      #pragma acc update self(tmp1[0:1])

      pb2[ind1] = Z1*norm[i1]*tmp1[0];

      //printf("    Z1: %5.3f norm: %5.3f tmp1: %10.6f v1: %10.6f \n",Z1,norm[i1],tmp1[0],pb2[ind1]);
    }

   //basis functions on atom n
    for (int i1=s1;i1<s2;i1++)
    {
      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; double zeta1 = basis1[3];

      if (n1<3 && l1==0)
      {
        int ind1 = n*N+i1;

        double v1 = norm[i1];
        if (n1==1) v1 *= -zeta1;
        pb1[ind1] = v1;
        //printf("    zeta: %5.3f norm: %5.3f val1: %10.6f \n",zeta1,norm[i1],v1);
      }
    }


  }

 //accumulate individual terms
  for (int i=0;i<nN;i++)
    pb[i] = pb1[i] + pb2[i];

  #pragma acc exit data delete(grid1[0:6],tmp1[0:1])

  delete [] n2i;
  delete [] norm;

  delete [] grid1;
  delete [] tmp1;
  delete [] pb1;
  delete [] pb2;

  if (prl>1)
  {
    printf("  pB: \n");
    print_rectangle_sm(natoms,N,pb);
  }

  return;
}

void compute_Pc(int natoms, int N, double* pB, double* Pc)
{
  int prl = 1;

  int N2 = N*N;
  int n2 = natoms*natoms;
  double W[n2];

  for (int i=0;i<natoms;i++)
  for (int j=0;j<natoms;j++)
  {
    double w1 = 0.;
    for (int m=0;m<N;m++)
      w1 += pB[i*N+m]*pB[j*N+m];
    W[i*natoms+j] = w1;
  }
  invert_stable_cpu(W,natoms,0.);
  //invert_stable_cpu(W,natoms,1.e-10,1);

  for (int i=0;i<N2;i++) Pc[i] = 0.;
  for (int i=0;i<N;i++)
    Pc[i*N+i] = 1.;

  for (int m=0;m<natoms;m++)
  for (int n=0;n<natoms;n++)
  {
    for (int i=0;i<N;i++)
    for (int j=0;j<N;j++)
      Pc[i*N+j] -= pB[m*N+i]*W[m*natoms+n]*pB[n*N+j];
  }

  if (prl>1)
  {
    printf("  pB: ");
    for (int j=0;j<N*natoms;j++)
      printf(" %8.5f",pB[j]);
    printf("\n");
  }
  if (prl>1)
  {
    printf("  Pc: \n");
    print_square(N,Pc);
  }

 #if 0
  double tmp1[N2];
  double tmp1e[N];
  for (int i=0;i<N2;i++)
    tmp1[i] = Pc[i];
  la_diag(N,N,tmp1,tmp1e);

  printf(" eigenvalues:");
  for (int i=0;i<N;i++)
    printf(" %8.5f",tmp1e[i]);
  printf("\n");
 #endif

  return;
}

void evaluate_cusps(int nrad, int nang, int natoms, int* atno, double* coordsf, vector<vector<double> > basis, double* Pao, float* grid)
{
  printf("\n checking cusp condition \n");
  printf("  WARNING: no angular wts included (yet) \n");

  int gs = nrad*nang;
  int gsa = gs*natoms;

  double* rho = new double[gsa]();
  #pragma acc enter data create(rho[0:gsa])

  compute_rho(natoms,atno,coordsf,basis,Pao,nrad,gsa,grid,rho,NULL,1);
  #pragma acc exit data copyout(rho[0:gsa])

  float rn[nrad];
  for (int n=0;n<natoms;n++)
  {
    double Zn = atno[n];
    rgrid_one_atom(nrad,atno[n],rn);
    int i0 = n*gs;

    for (int m=0;m<5;m++)
    {
      double ra = rn[m]; double rb = rn[m+1]; double dr = rb-ra;
      printf(" [%i]  r: %10.8f r+1: %10.8f ",n+1,ra,rb);

      int i1a = i0+m*nang;
      int i1b = i0+(m+1)*nang;
      double dt = 0.;
      double drhodr = 0.;
      for (int k=0;k<nang;k++)
      {
        int i2a = i1a+k; int i2b = i1b+k;
        double d1 = rho[i2a]; double d2 = rho[i2b];
  
 	dt += (d1+d2);
        drhodr += (d2-d1)/dr;
      }
 
      double rho_avg = dt/(2*nang);
      drhodr /= nang;
      double tzr = -2.*Zn*rho_avg;
      double del = (tzr-drhodr)/tzr;
      printf("  -2Zrho:  %11.6f  drho/dr:  %11.6f  (fraction: %8.5f) \n",tzr,drhodr,del);
    }
  }
  
  delete [] rho;
  return;
}

void check_cusp(int No, int natoms, int* atno, vector<vector<double> >& basis, double* pB, double* jCA)
{
  int N = basis.size();

  printf("  checking cusps \n");
  for (int a1=0;a1<natoms;a1++)
  {
    printf("  atom %2i: ",a1+1);
    for (int i=0;i<No;i++)
    {
      double v1 = 0.;
      for (int j=0;j<N;j++)
        v1 += pB[a1*N+j]*jCA[j*N+i];
      printf("  %2i - %8.5f",i+1,v1);
    }
    printf("\n");
  }

  return;
}

int prepare_PSP(int natoms, int N, double* S, double* Pc, double* Xp, cusolverDnHandle_t cu_hdl, int prl)
{
 //CPMZ need to clean this up
  int N2 = N*N;
  double thresh = PSP_THRESH;
  double threshw = PSP_THRESH_WARNING;
  double INF = 1.e300; //numeric_limits<double>::infinity();

  double* X = new double[N2];
  double* Xs = new double[N2];
  double* Xe = new double[N];
  #pragma acc enter data create(X[0:N2],Xs[0:N2],Xe[0:N])
  double tmp1[N2];

 //build PSP
  double PS[N2];
  for (int i=0;i<N2;i++) PS[i] = 0.;
  mat_times_mat_cpu(PS,Pc,S,N);

  for (int i=0;i<N2;i++) X[i] = 0.;
  mat_times_mat_cpu(X,PS,Pc,N);
  for (int i=0;i<N2;i++)
    X[i] *= -1;

  //double PSP[N2];
  //for (int i=0;i<N2;i++)
  //  PSP[i] = X[i];

  #pragma acc update device(X[0:N2])

 //get X, which diagonalizes PSP
  for (int i=0;i<N;i++) Xe[i] = 0.;
 #if USE_ACC
  #pragma acc update device(Xe[0:N])
  diagonalize_cusolver(N,N,X,Xe,cu_hdl);
  #pragma acc update self(Xe[0:N])
 #else
  la_diag(N,N,X,Xe);
 #endif

  for (int i=0;i<N;i++)
  if (std::isnan(Xe[i]))
  {
    printf("\n ERROR: NaN found in PSP eigenvalues \n");
    exit(1);
  }

  if (prl>-3)
  {
    int nprint = natoms+3; if (prl>1) nprint = N;
    printf("  PSP eigenvalues:");
    int N1 = N; if (N1>nprint) N1 = nprint;
    for (int i=0;i<N1;i++)
      printf(" %4.1e",Xe[i]);
    if (N1!=N) printf(" ...");
    for (int i=max(N1,N-nprint);i<N;i++)
      printf(" %4.1e",Xe[i]);
    printf("\n");

    int warning = -natoms;
    for (int i=0;i<N;i++)
    if (fabs(Xe[i])<threshw)
      warning++;
    if (warning)
      printf("  WARNING: found too many relatively low eigenvalues (%i) \n",warning);
  }

  int nsmall1 = 0;
  for (int i=0;i<N;i++)
  if (fabs(Xe[i])<thresh || Xe[i]>0.)
    nsmall1++;
  if (natoms!=nsmall1 && prl>1)
    printf(" WARNING: incorrect number of cusp constraints: %2i found vs %2i expected \n",nsmall1,natoms);
  if (nsmall1<natoms)
  {
    printf(" ERROR: insufficient constraints in cusp \n");
    exit(1);
  }

  #pragma acc update self(X[0:N2])

  double Xt[N2];
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    Xt[i*N+j] = X[j*N+i];

 #if 0
  ao_to_mo_cpu(N,PSP,Pc,Xt,tmp1);
  printf("  diagonalized PSP: \n");
  print_square(N,Pc);
 #endif

  double Sp[N2];
  ao_to_mo_cpu(N,S,Sp,Xt,tmp1);
  for (int i=0;i<N2;i++)
    Xs[i] = Sp[i];

 //leaving the order alone
  for (int i=0;i<N2;i++)
    Xs[i] *= -1;

  //nsmall = natoms;
  for (int i=N-nsmall1;i<N;i++)
  for (int j=0;j<N;j++)
    Xs[i*N+j] = Xs[j*N+i] = 0.;
  for (int i=N-nsmall1;i<N;i++)
    Xs[i*N+i] = INF;

  if (prl>1)
  {
    printf("\n Xs: \n");
    print_square(N,Xs);
  }


  if (prl>0) printf("  orthogonalizer \n");

 #if USE_ACC && 0
 //possible bug: diagonalize doesn't work for diagonal matrices?
  #pragma acc update device(Xs[0:N2])
  diagonalize_cusolver(N,N,Xs,Xe,cu_hdl);
  #pragma acc update self(Xs[0:N2],Xe[0:N])
 #else
  la_diag(N,N,Xs,Xe);
 #endif

  if (prl>2)
  {
    printf("\n Xs eigenvectors: \n");
    print_square(N,Xs);
  }

  if (prl>2)
  { //same as PSP eigenvalues
    printf("   Xs ev:");
    for (int n=0;n<N;n++)
      printf(" %8.5f",Xe[n]);
    printf("\n");
  }

  int nsmall = 0;
  for (int i=0;i<N;i++)
  //if (fabs(Xe[i])<thresh)
  if (fabs(Xe[i])<thresh || Xe[i]>0.)
  {
    Xe[i] = INF;
    nsmall++;
  }
  //printf("   %2i small ev \n",nsmall);

 #if 1
  //nsmall = natoms;
  for (int i=0;i<N;i++)
    Xe[i] = fabs(Xe[i]);
 #endif

  if (prl>1)
  {
    printf("   ev:");
    for (int n=0;n<N;n++)
      printf(" %8.5f",Xe[n]);
    printf("\n");
  }

  for (int i=0;i<N2;i++) tmp1[i] = 0.;
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    tmp1[i*N+j] += Xs[j*N+i]/sqrt(Xe[j]);
  for (int i=0;i<N2;i++)
    Xs[i] = tmp1[i];

#if 0
 #if USE_ACC
  #pragma acc update device(Xs[0:N2])
  nsmall = mat_root_inv_stable_cusolver(Xs,N,thresh2,cu_hdl);
  #pragma acc update self(Xs[0:N2])
 #else
  nsmall = mat_root_inv_stable_cpu(Xs,N,thresh_inv,prl-1);
 #endif
#endif

  if (prl>2)
  {
    printf("\n Xs': \n");
    print_square(N,Xs);
  }

 //this is the object we are looking for
 #if 1
  for (int i=0;i<N2;i++) Xp[i] = 0.;
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    Xp[i*N+j] += Xt[i*N+k] * Xs[k*N+j]; //for Xs = symm ortho
 #else
  if (0) //this only works if S is diagonalized by Xt
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    Xp[i*N+j] = Xt[i*N+j] / sqrt(Xe[j]);
 #endif

 #if 0
  printf("  eliminating %i rows of X \n",natoms);
  for (int n=N-natoms;n<N;n++)
  for (int j=0;j<N;j++)
    Xp[n*N+j] = 0.;
 #endif

  #pragma acc update device(Xp[0:N2])

  if (prl>-1)
  {
    printf(" Xp: \n");
    print_square(N,Xp);
  }

  #pragma acc exit data delete(X[0:N2],Xs[0:N2],Xe[0:N])
  delete [] X;
  delete [] Xs;
  delete [] Xe;

  return nsmall1;
}
