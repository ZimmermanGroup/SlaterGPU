#include "cusp.h"

#define PSP_THRESH_WARNING 1.e-8
#define PSP_THRESH 1.e-10
#define DOUBLE_CUSP 1

using namespace std;

void print_rectangle_sm(int N1, int N2, double* S);


void compute_diatomic_symm(int natoms, int* atno, vector<vector<double> > basis, vector<double*>& pB_all, int prl)
{
  if (natoms!=2) return;

  if (atno[0]!=atno[1])
  {
    printf("  compute_diatomic_symm requires Z1=Z2 \n");
    return;
  }

  printf("\n TESTING compute_diatomic_symm \n");

  int nb0 = pB_all.size();
  int N = basis.size();
  int nN = natoms*N;
  int M = N/2;
  double ort2 = 1./sqrt(2.);

  bool is_dfghi[N];
  for (int j=0;j<N;j++) is_dfghi[j] = 0;

 //select l>1, m!=0 basis
  int nc = 0;
  for (int j=0;j<N;j++)
  if (basis[j][1]>1)
  {
    //printf("   basis %i is l>1 \n",j);
    bool is_z = !basis[j][2];
    double zt1 = basis[j][3];

    if (!is_z) //&& zt1 < 2.1
    {
      is_dfghi[j] = 1;
      nc++;
    }
  }
  nc /= 2;

  for (int j=0;j<nc;j++)
  {
    double* pB1 = new double[nN]();
    pB_all.push_back(pB1);
  }

  printf("  nb0: %i  nc: %i  size pB_all: %i \n",nb0,nc,pB_all.size());

  int wc = 0;
  for (int j=0;j<M;j++)
  if (is_dfghi[j])
  {
    double* pB1 = pB_all[nb0+wc];
   #if 0
    if (basis[j][2]==0)
    {
     //z symmetry
      pB1[j] = pB1[N+j] = ort2;
      pB1[M+j] = pB1[N+M+j] = ort2;
    }
    else
   #endif
    {
     //xy symmetry
      pB1[j] = pB1[N+j] = ort2;
      pB1[M+j] = pB1[N+M+j] = -ort2;
    }
    wc++;
  }

  return;
}

void compute_jellium_cusp(int order, int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* pB1, double* pB2, int prl)
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

  double f0 = 0.;
  double f1 = -0.5 * Ne * pow(Rc,-3);
 //CPMZ testing
  //f1 -= Zg*ztg;

  double c = 1./Rc; //characteristic exponent of jellium potential
  double c2 = c*c;
  double c3 = c2*c;
  double c4 = c2*c2;
  double c5 = c4*c;
  double c6 = c4*c2;
  double c7 = c5*c2;
  double c8 = c4*c4;
  if (order>2)
  {
    double d1 = 3.*sqrt(3.)*Ne / (2.*c*c*c*PI);
    f0 = -d1*c*c; //zero order term
    //f0 = -1.1;
    f1 =  d1*pow(c,4)/3.; //2nd order term
    //f1 = d1*c*c;

   #if 0
    if (cusp_type==1)
      f1 = 0.; //r2 term zero, suitable for 1s basis
    else
      f1 = -(1.*Ne)/order*pow(Rc,-(order+1)); //r^order term
   #endif
  }
  else if (cusp_type>1)
    printf("\n WARNING: jellium order %i incompatible w/ cusp_type %i \n",order,cusp_type);

  int N = basis.size();
  double* norm = new double[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];

  double SHIFTER = read_float("SHIFTER");
 // conditions
  double* pb1 = new double[N]();
  double* pb2 = new double[N]();
  {
    for (int i1=0;i1<N;i1++)
    {
      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; double zeta1 = basis1[3];

     //testing
      if (l1==0)
      {
        double a = zeta1;

        pb1[i1] = (f1 - pow(a,4)/3.)*norm[i1]; //SS_V4 basis
       #if 0
        if (n1==1)
          pb1[i1] = f1 + pow(zeta1,4)/12.; //SS_V3 basis
        else if (n1==2)
          pb1[i1] = f1 + pow(zeta1,6)/240.; //SS_V3 basis
        else
        {
          //printf("\n ERROR: cusp doesn't support jellium 3s basis in \n"); exit(-1);
        }
       #endif
      }
      else
        pb1[i1] = 0.;
    }
  }

  for (int i=0;i<N;i++)
    pB1[i] = pb1[i];
  if (pB2!=NULL)
  for (int i=0;i<N;i++)
    pB2[i] = pb2[i];

  delete [] pb1;
  delete [] pb2;
  delete [] norm;

  return;
}

void compute_jellium_v1_cusp(int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* pB1, double* pB2, int prl)
{
  double Rc = read_float("RC");
  printf("  compute_jellium_v1_cusp for jellium order: %i  Rc: %8.5f \n",1,Rc);

  double c = 1./Rc;
  double c2 = c*c;

  int N = basis.size();

  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);

  double* norm = new double[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];

  for (int i=0;i<N;i++)
    pB1[i] = 0.;
  if (pB2!=NULL)
  for (int i=0;i<N;i++)
    pB2[i] = 0.;

  for (int n=0;n<natoms;n++)
  //if (not_X[n])
  {
    int s1 = 0; if (n>0) s1 = n2i[n-1]; int s2 = n2i[n];

    double Z1 = (double)atno[n];
    double A1 = coords[3*n+0]; double B1 = coords[3*n+1]; double C1 = coords[3*n+2];
    for (int i1=s1;i1<s2;i1++)
    {
      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

      double t1 = 1.;
      if (n1==2) t1 = 0.;

      if (n1<3 && l1==0)
      {
        int ind1 = n*N+i1;

        double v1 = norm[i1];
        if (n1==1) v1 *= -zeta1;
        pB1[ind1] = v1;

        if (pB2!=NULL)
          pB2[ind1] = zeta1*zeta1/2.*t1*norm[i1] - c2/2.*t1*norm[i1];
      }

    } //loop i1 over basis
  } //loop n over atoms

  delete [] n2i;
  delete [] norm;

  return;
}

void compute_cusp(int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* pB1, double* pB2, int prl)
{
  int jellium = read_int("JELLIUM");
  if (jellium!=3 && jellium>1)
    return compute_jellium_cusp(jellium,natoms,atno,coords,basis,pB1,pB2,prl);
  if (jellium==3 && natoms>1) { printf("\n ERROR: cannot do jellium cusp (type 3) \n"); exit(-1); }

 //Slater basis set, jellium has no 1/r term:
  if (jellium==1)
    return compute_jellium_v1_cusp(natoms,atno,coords,basis,pB1,pB2,prl);

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

    } //nuclear charge Z at atom

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
    pB1[i] = pb1[i] + pb2[i];
  if (pB2!=NULL)
  for (int i=0;i<nN;i++)
    pB2[i] = 0.;

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
    print_rectangle_sm(natoms,N,pB1);
    if (pB2!=NULL)
    {
      printf("  pB2: \n");
      print_rectangle_sm(natoms,N,pB2);
    }
  }

  return;
}

void compute_cusp(int natoms, int* atno, double* coords, vector<vector<double> > &basis, vector<double*> pB, int prl)
{
  double* pB1 = pB[0];
  double* pB2 = NULL; if (pB.size()>1) pB2 = pB[1];
  return compute_cusp(natoms,atno,coords,basis,pB1,pB2,prl);
}

void project_Pc(int natoms, int N, double* pB, double* Pc)
{
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

  //printf(" W: \n");
  //print_square(natoms,W);

  bool is_zero = 1;
  for (int i=0;i<natoms;i++)
  for (int j=0;j<natoms;j++)
  if (W[i*natoms+j]!=0.)
  { is_zero = 0; break; }

  if (is_zero)
  {
    printf("   WARNING: project_Pc cannot invert constraint matrix \n");
    return;
  }

  invert_stable_cpu(W,natoms,1.e-12);
  //invert_stable_cpu(W,natoms,0.);
  //invert_stable_cpu(W,natoms,1.e-10,1);

  //printf(" Winv: \n");
  //print_square(natoms,W);

  for (int m=0;m<natoms;m++)
  for (int n=0;n<natoms;n++)
  {
    for (int i=0;i<N;i++)
    for (int j=0;j<N;j++)
      Pc[i*N+j] -= pB[m*N+i]*W[m*natoms+n]*pB[n*N+j];
  }

  return;
}


void compute_Pc(int natoms, int N, vector<double*> pB_all, vector<double*> Pc_all)
{
  int prl = 1;

  int N2 = N*N;

  for (int n=0;n<Pc_all.size();n++)
  {
    double* Pc = Pc_all[n];
    double* pB = pB_all[n];

    for (int i=0;i<N2;i++) Pc[i] = 0.;
    for (int i=0;i<N;i++)
      Pc[i*N+i] = 1.;

    project_Pc(natoms,N,pB,Pc);

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

void compute_Pc(int natoms, int N, double* pB, double* Pc)
{
  vector<double*> pB_all;
  vector<double*> Pc_all;
  pB_all.push_back(pB);
  Pc_all.push_back(Pc);
  return compute_Pc(natoms,N,pB_all,Pc_all);
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

void check_cusp(int No, int natoms, int* atno, vector<vector<double> >& basis, vector<double*> pB, double* jCA)
{
  double* pB1 = pB[0];
  return check_cusp(No,natoms,atno,basis,pB1,jCA);
}

void project_S(int N, double* S, vector<double*> Pc_all, double* X)
{
  int N2 = N*N;
  double PS[N2];

  for (int i=0;i<N2;i++)
    X[i] = S[i];

  int nPc = Pc_all.size();
  for (int n=0;n<nPc;n++)
  {
    double* Pc = Pc_all[n];

   //build PSP
    for (int i=0;i<N2;i++) PS[i] = 0.;
    mat_times_mat_cpu(PS,Pc,X,N);

    for (int i=0;i<N2;i++) X[i] = 0.;
    mat_times_mat_cpu(X,PS,Pc,N);
  }

  return;
}

int prepare_PSP(int natoms, int N, double* S, vector<double*> Pc_all, double* Xp, cusolverDnHandle_t cu_hdl, int prl)
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

  project_S(N,S,Pc_all,X);

  int nPc = Pc_all.size();
  int nlowexpected = nPc + natoms-1;
  //printf(" X: \n");
  //print_square(N,X);

  for (int i=0;i<N2;i++)
    X[i] *= -1;

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
    int nprint1 = 5; if (nprint1>N) nprint1 = N;
    int nprint = nlowexpected+2; if (prl>1) nprint = N;
    printf("  PSP eigenvalues:");
    int N1 = N; if (N1>nprint1) N1 = nprint1;
    for (int i=0;i<N1;i++)
      printf(" %4.1e",Xe[i]);
    if (N1!=N) printf(" ...");
    for (int i=max(N1,N-nprint);i<N;i++)
      printf(" %4.1e",Xe[i]);
    printf("\n");

    int nlow = 0;
    //if (Pc2!=NULL) warning *= 2;
    for (int i=0;i<N;i++)
    if (fabs(Xe[i])<threshw)
      nlow++;
    if (nlow>nlowexpected)
    {
      printf("  WARNING: found too many relatively low eigenvalues (%i) \n",nlow-nlowexpected);
      for (int i=0;i<nlow;i++)
      {
        int i1 = N-i-1;
        printf("    %i:",i1);
        for (int j=0;j<N;j++)
        if (fabs(X[i1*N+j])>0.2)
          printf("   %i %5.2f",j,X[i1*N+j]);
        printf("\n");
      }
    }
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

int prepare_PSP(int natoms, int N, double* S, double* Pc, double* Xp, cusolverDnHandle_t cu_hdl, int prl)
{
  vector<double*> Pc_all;
  Pc_all.push_back(Pc);
  return prepare_PSP(natoms,N,S,Pc_all,Xp,cu_hdl,prl);
}
