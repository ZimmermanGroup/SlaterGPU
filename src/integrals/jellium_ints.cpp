#include "jellium_ints.h"

//debug: use only Slaters or Gaussians
//#define USE_SLATER 0
#define USE_GAUSS 0
#define ATOMIC_VR 0

#define SHIFT_EN 0

//will be divided by 10 in murak_zeta
#define Z1J 0.2
#define MMURAK 2

//to do
// 1. add Gaussian impurity (done)
// 2. test Gaussian basis set (done)
// 3. RI basis for Jellium? (Slater)


int get_gs1(int nrad, int nang, double* grid, double Rc)
{
  int gs1 = 0;
  for (int j=0;j<nrad;j++)
  if (grid[6*nang*j+3]>Rc)
  {
    gs1 = j;
    break;
  }
  gs1 *= nang;

  return gs1;
}

int get_gs1(int nrad, int nang, float* gridf, double Rc)
{
  int gs1 = 0;
  for (int j=0;j<nrad;j++)
  if (gridf[6*nang*j+3]>Rc)
  {
    gs1 = j;
    break;
  }
  gs1 *= nang;

  return gs1;
}

void symmetrize_diatomic(int N1, int N, double* jCA, int* aol, int* aom)
{
  //separates xy DOF from z

 //find largest component in jCA
  for (int j=0;j<N1;j++)
  {
    double cmax = 0.;
    int cl = -1;
    int cm = -111;

    for (int k=0;k<N;k++)
    {
      double cn = fabs(jCA[k*N+j]);
      if (cn>cmax)
      {
        cl = aol[k];
        cm = aom[k];
        cmax = cn;
      }
    }

    //printf("  orb %2i lm: %i %2i \n",j,cl,cm);
    if (cl<0 || cm<-110)
      printf(" ERROR: couldn't find lm (%i %i) in symmetrize_mos (%2i) \n",cl,cm,j);

   //eliminate nonmatching terms
    if (cm==0) //z-axis basis
    {
      for (int k=0;k<N;k++)
      if (aom[k]!=0)
        jCA[k*N+j] = 0.;
    }
    else //everything else
    {
      for (int k=0;k<N;k++)
      if (aom[k]!=cm)
        jCA[k*N+j] = 0.;
    }
  }

  return;
}

void symmetrize_atom(int N1, int N, double* jCA, int* aol, int* aom)
{
 //separates xyz DOF

 //find largest component in jCA
  for (int j=0;j<N1;j++)
  {
    double cmax = 0.;
    int cl = -1;
    int cm = -111;

    for (int k=0;k<N;k++)
    {
      double cn = fabs(jCA[k*N+j]);
      if (cn>cmax)
      {
        cl = aol[k];
        cm = aom[k];
        cmax = cn;
      }
    }

    //printf("  orb %2i lm: %i %2i \n",j,cl,cm);
    if (cl<0 || cm<-110)
      printf(" ERROR: couldn't find lm (%i %i) in symmetrize_mos (%2i) \n",cl,cm,j);

   //eliminate nonmatching terms
    for (int k=0;k<N;k++)
    if (aol[k]!=cl || aom[k]!=cm)
      jCA[k*N+j] = 0.;
  }

  return;
}

void symmetrize_MOs(int natoms, int N1, vector<vector<double> > basis, double* jCA, double* S)
{
  if (natoms>2) return;

  printf("\n symmetrize_MOs N1: %2i \n",N1);

  int N = basis.size();
  int N2 = N*N;

  int prl = 1;
  if (prl>1)
  {
    int Npr = max(N1,4);
    printf("\n before symmetrize: \n");
    print_mos_col(Npr,N,basis,jCA);
  }

  bool equal = 0;
  if (natoms==2)
  {
    equal = 1;
    double Z1 = basis[0][8];
    for (int j=1;j<N;j++)
    if (basis[j][8]!=Z1)
    {
      equal = 0; break;
    }

   //found homodiatomic
    if (equal)
    {
      printf("  equalizing homodiatomic \n");
      int M = N/2;
      for (int i=0;i<N1;i++)
      for (int j=0;j<M;j++)
      {
        double s1 = sign(jCA[j*N+i]);
        double v1 = fabs(jCA[j*N+i]);
        double s2 = sign(jCA[(M+j)*N+i]);
        double v2 = fabs(jCA[(M+j)*N+i]);

        //printf("  sv1: %3.1f %8.5f   sv2: %3.1f %8.5f \n",s1,v1,s2,v2);
        double v12 = 0.5*(v1+v2);
        jCA[j*N+i] = v12*s1;
        jCA[(M+j)*N+i] = v12*s2;
      }
    }
  }

  int aol[N]; //angular momentum
  int aom[N]; //m quantum #
  for (int j=0;j<N;j++)
  {
    aol[j] = basis[j][1];
    aom[j] = basis[j][2];
  }

  if (natoms==1)
    symmetrize_atom(N1,N,jCA,aol,aom);
  else if (natoms==2) // && equal)
    symmetrize_diatomic(N1,N,jCA,aol,aom);

  double* O = new double[N2];
  compute_CSC(N,jCA,S,jCA,O);

  double norms[N1];
  for (int j=0;j<N1;j++)
    norms[j] = 1./sqrt(O[j*N+j]);

  printf(" norms: ");
  for (int j=0;j<N1;j++)
    printf(" %12.8e",norms[j]);
  printf("\n");

  for (int j=0;j<N1;j++)
  for (int k=0;k<N;k++)
    jCA[k*N+j] *= norms[j];

  if (prl>1)
  {
    printf("\n after symmetrize: \n");
    print_mos_col(4,N,basis,jCA);
  }

  delete [] O;

  return;
}

//auxiliary basis potential
// this doesn't work since there is no preconditioner
void compute_V_jellium(vector<vector<double> > basis, int nrad, double* rgrid, double** V)
{
  int gs = nrad;
  //int gs6 = 6*gs;
  int N = basis.size();

  double fp = 4.*PI;

 //step size for SD solver
  double alpha = 0.0001;
  int maxsteps = read_int("VSTEPS");
  if (maxsteps<1) maxsteps = 30;
  double rmax = 0.1;

  double* U1 = new double[gs];
  double* U2 = new double[gs];
  #pragma acc enter data create(U1[0:gs],U2[0:gs])

  double* grad = new double[gs];
  double* step = new double[gs];
  double* valt = new double[gs];
  #pragma acc enter data create(grad[0:gs],step[0:gs],valt[0:gs])

  for (int i1=0;i1<N;i1++)
  {
    int n1 = basis[i1][0]; int l1 = basis[i1][1]; int m1 = basis[i1][2];
    double zt1 = basis[i1][3]; double norm1 = basis[i1][4];
    int llp = l1*(l1+1);

    if (l1>0) { printf("\n WARNING: code not ready for l>0 \n"); exit(-1); }

   //initial potential
    double normv = norm_sv(n1,l1,m1,zt1);
    double oz2 = 1./zt1/zt1;
    double toz3 = 2.*oz2/zt1;
   #pragma acc parallel loop present(U1[0:gs],rgrid[0:gs])
    for (int j=0;j<gs;j++)
    {
      double r = rgrid[j];
      double ezr = exp(-zt1*r);
      double tz3r = toz3/r;

      U1[j] = -ezr*(oz2+tz3r)+tz3r;
      U1[j] *= normv*r;
      //U1[j] = 0.;
    }

   #pragma acc parallel loop present(valt[0:gs])
    for (int j=0;j<gs;j++)
      valt[j] = norm1;

  //evaluate ns ftns
    int n1m = n1-1;
   #pragma acc parallel loop present(valt[0:gs],rgrid[0:gs])
    for (int j=0;j<gs;j++)
    {
      double r = rgrid[j];
      double rp = pow(r,n1m);
      double ezr = exp(-zt1*r);
      valt[j] *= rp*ezr;
    }

    #pragma acc update self(valt[0:gs],U1[0:gs])
    printf("  r          rho         V0 \n");
    for (int j=0;j<gs;j++)
      printf("  %8.5f  %8.5f  %8.5f \n",rgrid[j],valt[j],U1[j]/rgrid[j]);
    printf("\n");

   //check this for l1>0
    //eval_shd(i1,gs,grid,valt,n1,l1,m1,zt1);

    double gp = 1.;
    for (int ns=0;ns<maxsteps;ns++)
    {
     #pragma acc parallel loop present(grad[0:gs])
      for (int j=0;j<gs;j++)
        grad[j] = 0.;

     #pragma acc parallel loop present(rgrid[0:gs],grad[0:gs],valt[0:gs],U1[0:gs],U2[0:gs])
      for (int j=1;j<gs-1;j++)
      {
        int j0 = j-1;
        int j1 = j;
        int j2 = j+1;

        double rm = rgrid[j0];
        double r  = rgrid[j1];
        double rp = rgrid[j2];
        double r2 = r*r;

        double h1 = r-rm;
        double h2 = rp-r;
        double h12 = (h1+h2)*0.5;

        double um = U1[j0];
        double u  = U1[j1];
        double up = U1[j2];

        double v1 = valt[j];
        double d2u = ((up-u)/h2 - (u-um)/h1)/h12;

        double rhs = fp*r*v1;
        double res = d2u - llp/r2 + rhs;

        if (j<10 || (j>gs/3 && j<gs/2))
        {
          printf("  r: %6.1e %6.1e %6.1e  h: %6.1e %6.1e %6.1e \n",rm,r,rp,h1,h2,h12);
          printf("  r: %6.1e  du12:  %9.2e %9.2e \n",r,(up-u)/h2,(u-um)/h1);
          printf("  r: %6.1e  d2u: %9.6f  u: %9.6f  rrho: %9.6f  res: %9.6f \n",r,d2u,u,rhs,res);
        }

        if (res>rmax) res = rmax; if (res<-rmax) res = -rmax;
        if (res>rhs) res = rhs; if (res<-rhs) res = -rhs;
        grad[j] = res;
      }

     //Hessian?
     #pragma acc parallel loop present(grad[0:gs],rgrid[0:gs])
      for (int j=1;j<gs-1;j++)
      {
        double dr = rgrid[j]-rgrid[j-1];
        grad[j] *= dr;
      }

    //forward finite difference from pt 0 equals that of pt 1
     #pragma acc serial present(grad[0:gs])
      {
        grad[0] = grad[1];
        grad[gs-1] = grad[gs-2];
        grad[0] = 0.;
      }

      double g1 = 0.;
     #pragma acc parallel loop present(grad[0:gs]) reduction(+:g1)
      for (int j=0;j<gs;j++)
        g1 += grad[j]*grad[j];

      double beta = g1/gp;
      if (beta>1.) beta = 1.;
      if (ns==0) beta = 0.;
      //beta = 0.;
      gp = g1;

      printf("  iter %4i  g1: %9.2e  beta: %8.5f \n",ns,g1,beta);

     #pragma acc parallel loop present(grad[0:gs],step[0:gs])
      for (int j=0;j<gs;j++)
        step[j] = alpha*grad[j] + beta*step[j];

     #pragma acc parallel loop present(U1[0:gs],step[0:gs])
      for (int j=0;j<gs;j++)
        U1[j] += step[j];

     #pragma acc parallel loop present(U1[0:gs])
      for (int j=0;j<gs;j++)
      if (U1[j]<0.)
        U1[j] = 0.;

    } //loop ns over opt steps

    printf("\n V(final) \n");
    #pragma acc update self(U1[0:gs])
    printf("  r          U1 \n");
    for (int j=0;j<gs;j++)
      printf("  %8.5f  %8.5f \n",rgrid[j],U1[j]/rgrid[j]);
    printf("\n");

    double* V1 = V[i1];
   #pragma acc parallel loop present(U1[0:gs],V1[0:gs])
    for (int j=0;j<gs;j++)
      V1[j] = U1[j];
  }

  #pragma acc exit data delete(step[0:gs],valt[0:gs])
  #pragma acc exit data delete(U1[0:gs],U2[0:gs])
  delete [] step;
  delete [] valt;
  delete [] U1;
  delete [] U2;

  return;
}

void compute_C_jellium(bool use_slater, double Rc, vector<vector<double> > basis, vector<vector<double> > basis_aux, int nrad, int nang,
                       double* ang_g, double* ang_w, double* C)
{
  bool sgs_basis = 0;

  int N = basis.size();
  int Naux = basis_aux.size();
  int N2 = N*N;
  int N2a = N2*Naux;

  int gs = nrad*nang;
  int gs6 = 6*gs;

  double norm1[Naux];
  double norm2[N];
  for (int i=0;i<Naux;i++)
  {
   //aux basis is Slater only
    norm1[i] = norm_sv(basis_aux[i][0],basis_aux[i][1],basis_aux[i][2],basis_aux[i][3]);
  }
  for (int i=0;i<N;i++)
  {
   //primary basis is renormalized by integration
    norm2[i] = basis[i][4];
  }

  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])

  double* grid = new double[gs6];
  double* wt = new double[gs];
  #pragma acc enter data create(grid[0:gs6],wt[0:gs])

  int m_murak = MMURAK;
  double Z1 = Z1J;
  //bool murak_zeta = 1;
  //generate_central_grid_2d(-1,!murak_zeta,grid,wt,Z1,nrad,nang,ang_g,ang_w);
  generate_central_grid_3d(m_murak,grid,wt,Z1,nrad,nang,ang_g,ang_w);

  #pragma acc update self(grid[0:gs6])
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])

 #if 0
  printf("\n grid: \n");
  for (int j=0;j<gs;j++)
    printf("  %8.5f %8.5f %8.5f  r: %8.5f \n",grid[6*j+0],grid[6*j+1],grid[6*j+2],grid[6*j+3]);
 #endif

  int gs1 = get_gs1(nrad,nang,grid,Rc);
  int gs2 = gs;

  if (use_slater)
  {
    printf("  WARNING: Slater only \n");
    gs1 = 0;
  }
 #if USE_GAUSS
  printf("  DEBUG: Gaussian only \n");
  gs1 = nrad*nang;
 #endif

  int igs = N*gs;
  int iags = Naux*gs;
  double* valS = new double[igs];
  double* valV = new double[iags];
  //double** valS = new double*[N];
  //double** valV = new double*[Naux];
  //for (int i=0;i<N;i++)
  //  valS[i] = new double[gs];
  //for (int i=0;i<Naux;i++)
  //  valV[i] = new double[gs];
  //#pragma acc enter data create(valS[0:N][0:gs],valV[0:Naux][0:gs])
  #pragma acc enter data create(valS[0:igs],valV[0:iags])

  //instead of this, using Slater basis for potential
  //compute_V_jellium(basis_aux,gsr,rgrid,V);

 //Slater potentials
  for (int i1=0;i1<Naux;i1++)
  {
    vector<double> basis1 = basis_aux[i1];
    int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zt1 = basis1[3];
    //double zt1s = 2.*Rc*zt1;

    #pragma acc parallel loop present(valV[0:iags],wt[0:gs])
    for (int j=0;j<gs;j++)
      valV[i1*gs+j] = wt[j];

    eval_inr_r12(-1,gs,grid,&valV[i1*gs],n1,l1,zt1);
    eval_sh_3rd(gs,grid,&valV[i1*gs],n1,l1,m1);
  }

 //SGS or SS basis
  for (int i2=0;i2<N;i2++)
  {
    vector<double> basis2 = basis[i2];
    int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zt2 = basis2[3];

    #pragma acc parallel loop present(valS[0:igs])
    for (int j=0;j<gs;j++)
      valS[i2*gs+j] = 1.;

    if (use_slater)
      eval_shd(i2,gs2,grid,&valS[i2*gs],n2,l2,m2,zt2);
    else if (sgs_basis)
      eval_sgsd(i2,gs1,gs2,grid,&valS[i2*gs],n2,l2,m2,zt2,Rc);
    else
      eval_ssd(i2,gs,grid,&valS[i2*gs],n2,l2,m2,zt2,Rc);
  }

 //need to rearrange and eliminate redundancy
  for (int i1=0;i1<Naux;i1++)
  {
    double* valv = &valV[i1*gs];

   #pragma acc parallel loop collapse(2) present(valv[0:gs],valS[0:igs],C[0:N2a])
    for (int i2=0;i2<N;i2++)
    for (int i3=0;i3<N;i3++)
    {
      double* valm = &valS[i2*gs];
      double* valn = &valS[i3*gs];

      double v1 = 0.;
     #pragma acc loop reduction(+:v1)
      for (int j=0;j<gs;j++)
        v1 += valv[j]*valm[j]*valn[j];

      C[i1*N2+i2*N+i3] = v1;
    }
  }

  #pragma acc update self(C[0:N2a])

  for (int i1=0;i1<Naux;i1++)
  for (int i2=0;i2<N;i2++)
  for (int i3=0;i3<N;i3++)
    C[i1*N2+i2*N+i3] *= norm1[i1]*norm2[i2]*norm2[i3];

  double thresh = 1.e-12;
  for (int i=0;i<N2a;i++)
  if (fabs(C[i])<thresh)
    C[i] = 0.;

  transpose_C(Naux,N,C);

  bool normalize_A = check_file("An");
  if (normalize_A)
  {
    printf("  normalizing C matrix \n");
    double Anorm[Naux];
    read_array(Naux,Anorm,"An");

    int nna = N*Naux;
    for (int p=0;p<Naux;p++)
    {
      double v2 = Anorm[p];

      for (int i=0;i<N;i++)
      for (int j=0;j<N;j++)
        C[i*nna+j*Naux+p] *= v2;
    }
  }

  #pragma acc update device(C[0:N2a])

 //cleanup
  #pragma acc exit data delete(valS[0:igs],valV[0:iags])
  //#pragma acc exit data delete(valS[0:N][0:gs],valV[0:Naux][0:gs])
  //for (int i=0;i<N;i++)
  //  delete [] valS[i];
  delete [] valS;
  //for (int i=0;i<Naux;i++)
  //  delete [] valV[i];
  delete [] valV;

  #pragma acc exit data delete(grid[0:gs6],wt[0:gs])
  delete [] grid;
  delete [] wt;

  return;
}

void compute_Vr_jellium(int order, double Zg, double ztg, double Rc, int Ne, int gs1, int gs2, double* grid, double* Vr)
{
  int gs = gs2;
  int gs6 = 6*gs;

 #if ATOMIC_VR
  gs1 = 0;
 #endif

  if (order==2)
  {
    if (Rc<1.e-9) { printf(" ERROR: invalid Rc in compute_Vr_jellium \n"); Rc = 1.; }
    printf("  compute_Vr_jellium  order: %i  Ne: %2i  Zg/ztg: %8.5f %8.5f  Rc: %8.5f  gs12: %4i %4i \n",order,Ne,Zg,ztg,Rc,gs1,gs2);

  //using standard Jellium, though these two possibilities have merit
   #if 0
   //which potential type to use
    #define VR2 1

    double a = 10./Rc;
    double c1 = -Ne;

   //strength of attraction potential
   #if VR2
    double c2 = c1*a/6.;
   #else
    double c2 = c1*a/8.;
   #endif

    printf("  WARNING: testing new jellium potential a/c1/c2: %8.5f %8.5f %8.5f \n",a,c1,c2);

   #pragma acc parallel loop present(Vr[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
      double r = grid[6*j+3];
      double ar = a*r;
     #if VR2
     //(1+ar+a2r2/3)
      double num = a*(3.+ar*(3.-ar));
      double den = 18. + 6.*ar*(3.+ar);
     #else
     //(1+ar+2a2r2/5+a3r3/15)
      double num = a*(9.+ar*(9.+ar*(2.-ar)));
      double den = 4.*(30.+2.*ar*(15.+ar*(6.+ar)));
     #endif
      Vr[j] = c1*num/den + c2;
    }
   #endif

   #if 1
   //standard Jellium (quadratic then 1/r)
    //double c1 = -0.5*Ne/Rc;
    //double c2 = Rc*Rc;
    double c1 = -(1.*Ne)/order/Rc;
    double c2 = 1.+order;
    double c3 = pow(Rc,order);
    double c4 = -Ne;

   #pragma acc parallel loop present(Vr[0:gs],grid[0:gs6])
    for (int j=0;j<gs1;j++)
    {
      double r = grid[6*j+3];
      double rn = pow(r,order);
      Vr[j] = c1*(c2-rn/c3);
    }

   #pragma acc parallel loop present(Vr[0:gs],grid[0:gs6])
    for (int j=gs1;j<gs2;j++)
    {
      double r = grid[6*j+3];
      Vr[j] = c4/r;
    }
   #endif

  }
  else if (order==1 || order==4) //SS_V4 and SS_V5
  {
    double a = 1./Rc;
    double a2 = a*a;
    double a3 = a2*a;
    double c1 = 3.*sqrt(3.)*Ne / (2.*a3*PI);
    double c2 = c1*a2/2.;
    //c2 = 0.;

    printf("  compute_Vr_jellium. order: %i  Rc: %8.5f  a: %8.5f \n",order,Rc,a);

   #pragma acc parallel loop present(Vr[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
      double r = grid[6*j+3];
      //double r2 = r*r;
      double ar = a*r;

      Vr[j] = c1*a2*(-3. + ar*(-3.+ar))/(6.+ 2.*ar*(3.+ar)) - c2;
    }
  }
  else if (order==5) //not used
  {
    double a = 1./Rc;

    printf(" INCOMPLETE:  compute_Vr_jellium. order: %i  Rc: %8.5f  a: %8.5f \n",order,Rc,a);

   #pragma acc parallel loop present(Vr[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
      double r = grid[6*j+3];
      double ar = a*r;

      //Vr[j] = 1.;
    }
  }
  else if (order==6)
  {
    printf("\n ERROR: order 5 jellium to be implemented \n");
    exit(-1);
  }
  else
  {
    printf("\n ERROR: no jellium potential available at order %i \n",order);
    exit(-1);
  }

  if (Zg!=0.)
  {
   #pragma acc parallel loop present(Vr[0:gs],grid[0:gs6])
    for (int j=0;j<gs;j++)
    {
      double r = grid[6*j+3];
      Vr[j] -= Zg*exp(-ztg*r*r);
    }
  }

  #pragma acc update self(Vr[0:gs])

  return;
}

void compute_Vr_jellium(double Zg, double ztg, double Rc, int Ne, int gs1, int gs2, double* grid, double* Vr)
{
  int order = read_int("JELLIUM");
  if (order<2) order = 2;

  return compute_Vr_jellium(order,Zg,ztg,Rc,Ne,gs1,gs2,grid,Vr);
}

void compute_Vr_jellium(double Zg, double ztg, double Rc, int Ne, int gs1, int gs2, float* gridf, double* Vr)
{
  int gs = gs2;
  int gs6 = 6*gs;

  if (Rc<1.e-9) { printf(" ERROR: invalid Rc in compute_Vr_jellium \n"); Rc = 1.; }

  double* grid = new double[gs6];
  #pragma acc enter data create(grid[0:gs6])

 #pragma acc parallel loop present(grid[0:gs6],gridf[0:gs6])
  for (int j=0;j<gs6;j++)
    grid[j] = gridf[j];

  compute_Vr_jellium(Zg,ztg,Rc,Ne,gs1,gs2,grid,Vr);

  #pragma acc exit data delete(grid[0:gs6])
  delete [] grid;

  return;
}

//integrates Vr, use the previous ftn instead
void compute_Vr_jellium(double rs, double Rc, int nrad, int nang, double* ang_g, double* ang_w, double* grid, double* Vr)
{
  int gs = nrad*nang;
  int gs6 = 6*gs;
  int gsr = (nrad+1)*nang;
  int gsr6 = 6*gsr;

  double* gridr = new double[gsr6];
  double* wtr = new double[gsr];
  #pragma acc enter data create(gridr[0:gsr6],wtr[0:gsr])

  int m_murak = MMURAK;
  double Z1 = Z1J;
  //bool murak_zeta = 1;
  //generate_central_grid_2d(-1,1,gridr,wtr,1.,nrad+1,nang,ang_g,ang_w);
  generate_central_grid_3d(m_murak,gridr,wtr,Z1,nrad,nang,ang_g,ang_w);


 //value of density at each pt with r<Rc
  double tofp = 3./(4.*PI);
  double tofprs = tofp*pow(rs,-3.);

  #pragma acc parallel loop present(Vr[0:gs])
  for (int j=0;j<gs;j++)
    Vr[j] = 0.;

  #pragma acc parallel loop present(Vr[0:gs],grid[0:gs6],gridr[0:gsr6],wtr[0:gsr])
  for (int j=0;j<gs;j++)
  {
    double x1 = grid [6*j+0]; double y1 = grid [6*j+1]; double z1 = grid [6*j+2];

    double v1 = 0.;
   #pragma acc loop reduction(+:v1)
    for (int k=0;k<gsr;k++)
    if (gridr[6*k+3]<Rc)
    {
      double x2 = gridr[6*k+0]; double y2 = gridr[6*k+1]; double z2 = gridr[6*k+2];
      double wt2 = wtr[k];

      double x12 = x1-x2; double y12 = y1-y2; double z12 = z1-z2;
      double r12 = sqrt(x12*x12+y12*y12+z12*z12);

      //printf("  xyz: %8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f  r12: %8.5f  wt2: %8.5f \n",x1,y1,z1,x2,y2,z2,r12,wt2);
      v1 -= wt2/r12;
    }

    Vr[j] = v1*tofprs;
  }

  #pragma acc exit data delete(gridr[0:gsr6],wtr[0:gsr])
  delete [] gridr;
  delete [] wtr;

  return;
}

void fill_permute_ol4(int N, double* ol)
{
  int N2 = N*N;
  int N3 = N2*N;

  const double thresh = 1.e-12;

  for (int i1=0;i1<N;i1++)
  for (int i2=0;i2<=i1;i2++)
  for (int i3=0;i3<=i2;i3++)
  for (int i4=0;i4<=i3;i4++)
  {
    //1234
    //1243
    //1324
    //1342
    //1423
    //1432

    double v1 = ol[i1*N3+i2*N2+i3*N+i4];
    if (fabs(v1)<thresh) v1 = 0.;

    ol[i1*N3+i2*N2+i3*N+i4] = v1;
    ol[i1*N3+i2*N2+i4*N+i3] = v1;
    ol[i1*N3+i3*N2+i2*N+i4] = v1;
    ol[i1*N3+i3*N2+i4*N+i2] = v1;
    ol[i1*N3+i4*N2+i2*N+i3] = v1;
    ol[i1*N3+i4*N2+i3*N+i2] = v1;

    //2134
    //2143
    //2314
    //2341
    //2413
    //2431

    ol[i2*N3+i1*N2+i3*N+i4] = v1;
    ol[i2*N3+i1*N2+i4*N+i3] = v1;
    ol[i2*N3+i3*N2+i1*N+i4] = v1;
    ol[i2*N3+i3*N2+i4*N+i1] = v1;
    ol[i2*N3+i4*N2+i1*N+i3] = v1;
    ol[i2*N3+i4*N2+i3*N+i1] = v1;

    //3124
    //3142
    //3214
    //3241
    //3412
    //3421

    ol[i3*N3+i1*N2+i2*N+i4] = v1;
    ol[i3*N3+i1*N2+i4*N+i2] = v1;
    ol[i3*N3+i2*N2+i1*N+i4] = v1;
    ol[i3*N3+i2*N2+i4*N+i1] = v1;
    ol[i3*N3+i4*N2+i1*N+i2] = v1;
    ol[i3*N3+i4*N2+i2*N+i1] = v1;

    //4123
    //4132
    //4213
    //4231
    //4312
    //4321

    ol[i4*N3+i1*N2+i2*N+i3] = v1;
    ol[i4*N3+i1*N2+i3*N+i2] = v1;
    ol[i4*N3+i2*N2+i1*N+i3] = v1;
    ol[i4*N3+i2*N2+i3*N+i1] = v1;
    ol[i4*N3+i3*N2+i1*N+i2] = v1;
    ol[i4*N3+i3*N2+i2*N+i1] = v1;
  }
}

void compute_4c_ol_jellium(bool use_slater, double Rc, vector<vector<double> > basis,
   	 int nrad, int nang, double* ang_g, double* ang_w, double* ol, int prl)
{
  bool sgs_basis = 0;

  if (!sgs_basis)
    printf(" compute_4c_ol_jellium (ss basis) \n");
  else
    printf(" compute_4c_ol_jellium (sgs basis) \n");

  int N = basis.size();
  int N2 = N*N;
  int N3 = N2*N;

  int gs = nrad*nang;
  int gs6 = 6*gs;

  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])

  double* grid = new double[gs6];
  double* wt = new double[gs];
  #pragma acc enter data create(grid[0:gs6],wt[0:gs])

 //grid probably should be a bit larger
  int m_murak = MMURAK;
  double Z1 = Z1J;
  //bool murak_zeta = 1;
  //generate_central_grid_2d(-1,!murak_zeta,grid,wt,Z1,nrad,nang,ang_g,ang_w);
  generate_central_grid_3d(m_murak,grid,wt,Z1,nrad,nang,ang_g,ang_w);
  #pragma acc update self(grid[0:gs6])

  int gs1 = get_gs1(nrad,nang,grid,Rc);
  int gs2 = gs;

  if (use_slater)
  {
    printf("  DEBUG: Slater only \n");
    gs1 = 0;
  }
 #if USE_GAUSS
  printf("  DEBUG: Gaussian only \n");
  gs1 = nrad*nang;
 #endif

  int igs = N*gs;
  double* valS1 = new double[igs];
  double* valS2 = new double[igs];
  //double** valS1 = new double*[N];
  //double** valS2 = new double*[N];
  //for (int i=0;i<N;i++) valS1[i] = new double[gs];
  //for (int i=0;i<N;i++) valS2[i] = new double[gs];
  //#pragma acc enter data create(valS1[0:N][0:gs],valS2[0:N][0:gs])
  #pragma acc enter data create(valS1[0:igs],valS2[0:igs])

 #pragma acc parallel loop collapse(2) present(valS1[0:igs])
  for (int i1=0;i1<N;i1++)
  for (int j=0;j<gs;j++)
    valS1[i1*gs+j] = 1.;

 #pragma acc parallel loop collapse(2) present(valS2[0:igs],wt[0:gs])
  for (int i2=0;i2<N;i2++)
  for (int j=0;j<gs;j++)
    valS2[i2*gs+j] = wt[j];

  for (int i1=0;i1<N;i1++)
  {
    vector<double> basis1 = basis[i1];
    int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3]; //double norm1 = basis1[4];
    double zeta2 = 2.*zeta1*Rc;
    double norm2 = exp(zeta1*Rc*Rc);

    if (use_slater)
    {
      eval_shd(i1,gs2,grid,&valS1[i1*gs],n1,l1,m1,zeta1);
      eval_shd(i1,gs2,grid,&valS2[i1*gs],n1,l1,m1,zeta1);
    }
    else if (sgs_basis)
    {
      eval_sgsd(i1,gs1,gs2,grid,&valS1[i1*gs],n1,l1,m1,zeta1,Rc);
      eval_sgsd(i1,gs1,gs2,grid,&valS2[i1*gs],n1,l1,m1,zeta1,Rc);
    }
    else
    {
      eval_ssd(i1,gs,grid,&valS1[i1*gs],n1,l1,m1,zeta1,Rc);
      eval_ssd(i1,gs,grid,&valS2[i1*gs],n1,l1,m1,zeta1,Rc);
    }
  }

 //slow loops. need collapse
  for (int i1=0;i1<N;i1++)
  for (int i2=0;i2<=i1;i2++)
  for (int i3=0;i3<=i2;i3++)
  for (int i4=0;i4<=i3;i4++)
  {
    double* valm = &valS1[i1*gs];
    double* valn = &valS1[i2*gs];
    double* valp = &valS1[i3*gs];
    double* valq = &valS2[i4*gs];

    double v1 = 0.;
   #pragma acc parallel loop present(valm[0:gs],valn[0:gs],valp[0:gs],valq[0:gs]) reduction(+:v1)
    for (int j=0;j<gs;j++)
      v1 += valm[j]*valn[j]*valp[j]*valq[j];

   //need to symmetrize
    double n1234 = basis[i1][4]*basis[i2][4]*basis[i3][4]*basis[i4][4];
    ol[i1*N3+i2*N2+i3*N+i4] = n1234*v1;
  }

 //symmetrize
  fill_permute_ol4(N,ol);

  #pragma acc exit data delete(valS1[0:igs],valS2[0:igs])
  //#pragma acc exit data delete(valS1[0:N][0:gs],valS2[0:N][0:gs])
  //for (int i=0;i<N;i++)
  //  delete [] valS1[i];
  //for (int i=0;i<N;i++)
  //  delete [] valS2[i];
  delete [] valS1;
  delete [] valS2;

  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc exit data delete(grid[0:gs6],wt[0:gs])
  delete [] grid;
  delete [] wt;

  return;
}

void compute_STEn_jellium(bool use_slater, int order, double Zg, double ztg, double rs, double Rc, int No, bool update_norm, vector<vector<double> >& basis,
         int nrad, int nang, double* ang_g, double* ang_w, double* S, double* T, double* En, int prl)
{
 //in this ftn, rs isn't used. instead, compute_Vr is already analytically evaluated.

  bool sgs_basis = 0;

  int N = basis.size();
  int N2 = N*N;

  int Ne = 2*No;

  int gs = nrad*nang;
  int gs6 = 6*gs;

  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])

  double* grid = new double[gs6];
  double* wt = new double[gs];
  #pragma acc enter data create(grid[0:gs6],wt[0:gs])


 //grid probably should be a bit larger
  int m_murak = MMURAK;
  double Z1 = Z1J;
  bool murak_zeta = 1;
  //generate_central_grid_2d(-1,!murak_zeta,grid,wt,Z1,nrad,nang,ang_g,ang_w);
  generate_central_grid_3d(m_murak,grid,wt,Z1,nrad,nang,ang_g,ang_w);
  #pragma acc update self(grid[0:gs6])

  int gs1 = get_gs1(nrad,nang,grid,Rc);
  int gs2 = gs;
  printf("  gs1: %4i  gs2: %4i \n",gs1,gs2);

 #if 0
  printf("\n grid: \n");
  for (int j=0;j<gs;j++)
    printf(" (%3i)   %8.5f %8.5f %8.5f  r: %8.5f \n",j,grid[6*j+0],grid[6*j+1],grid[6*j+2],grid[6*j+3]);
 #endif

  #pragma acc parallel loop present(S[0:N2],En[0:N2],T[0:N2])
  for (int j=0;j<N2;j++)
    S[j] = En[j] = T[j] = 0.;

  {
   //multi-GPU not ready
    int tid = -1;

    int gs3 = 3*gs;

   //compute the potential
    double* Vr = new double[gs];
    #pragma acc enter data create(Vr[0:gs])

    compute_Vr_jellium(order,Zg,ztg,Rc,Ne,gs1,gs2,grid,Vr);
    //compute_Vr_jellium(rs,Rc,nrad,nang,ang_g,ang_w,grid,Vr);

    if (use_slater)
    {
      printf("  DEBUG: Slater only \n");
      gs1 = 0;
    }
   #if USE_GAUSS
    printf("  DEBUG: Gaussian only \n");
    gs1 = nrad*nang;
   #endif

    if (prl>1)
    {
      #pragma acc update self(Vr[0:gs])
      printf("\n       r      Vr: \n");
      for (int j=0;j<nrad;j++)
        printf("   %8.5f  %9.6f \n",grid[6*nang*j+3],Vr[j*nang]);
    }

   //evaluate basis ftns
    int igs = N*gs;
    int igs3 = N*gs3;
    double* valS1 = new double[igs];
    double* valS2 = new double[igs];
    double* valp1 = new double[igs3];
    double* valL1 = new double[igs];
    //double** valS1 = new double*[N];
    //double** valS2 = new double*[N];
    //double** valp1 = new double*[N];
    //double** valL1 = new double*[N];
    //for (int i=0;i<N;i++) valS1[i] = new double[gs];
    //for (int i=0;i<N;i++) valS2[i] = new double[gs];
    //for (int i=0;i<N;i++) valp1[i] = new double[gs3];
    //for (int i=0;i<N;i++) valL1[i] = new double[gs];
    double* tmp = new double[gs3];
    //#pragma acc enter data create(valS1[0:N][0:gs],valS2[0:N][0:gs],valp1[0:N][0:gs3],valL1[0:N][0:gs],tmp[0:gs3])
    #pragma acc enter data create(valS1[0:igs],valS2[0:igs],valp1[0:igs3],valL1[0:igs],tmp[0:gs3])

   #pragma acc parallel loop collapse(2) present(valS1[0:igs])
    for (int i1=0;i1<N;i1++)
    for (int j=0;j<gs;j++)
      valS1[i1*gs+j] = 1.;

   #pragma acc parallel loop collapse(2) present(valS2[0:igs],wt[0:gs])
    for (int i2=0;i2<N;i2++)
    for (int j=0;j<gs;j++)
      valS2[i2*gs+j] = wt[j];

   #pragma acc parallel loop collapse(2) present(valL1[0:igs])
    for (int i1=0;i1<N;i1++)
    for (int j=0;j<gs;j++)
      valL1[i1*gs+j] = 1.;

    for (int i1=0;i1<N;i1++)
    {
      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3]; //double norm1 = basis1[4];
      double zeta2 = 2.*zeta1*Rc;
      double norm2 = exp(zeta1*Rc*Rc);
      if (use_slater)
      {
        zeta2 = zeta1; norm2 = 1.;
      }

      if (use_slater)
      {
        eval_shd(i1,gs,grid,&valS1[i1*gs],n1,l1,m1,zeta1);
        eval_shd(i1,gs,grid,&valS2[i1*gs],n1,l1,m1,zeta1);
      }
      else if (sgs_basis)
      {
        eval_sgsd(i1,gs1,gs2,grid,&valS1[i1*gs],n1,l1,m1,zeta1,Rc);
        eval_sgsd(i1,gs1,gs2,grid,&valS2[i1*gs],n1,l1,m1,zeta1,Rc);
      }
      else
      {
        eval_ssd(i1,gs,grid,&valS1[i1*gs],n1,l1,m1,zeta1,Rc);
        eval_ssd(i1,gs,grid,&valS2[i1*gs],n1,l1,m1,zeta1,Rc);
      }

     //KE
      if (sgs_basis)
      {
       #pragma acc parallel loop present(tmp[0:gs3])
        for (int j=0;j<gs3;j++)
          tmp[j] = 1.;

       //Gaussian portion of KE
        //eval_gh_ke(gs1,grid,valKg[i1],n1,l1,1.,zeta1);
        eval_pd_gh(gs1,grid,tmp,n1,l1,m1,1.,zeta1);

        #pragma acc parallel loop present(valp1[0:igs3],tmp[0:gs3])
        for (int j=0;j<3*gs1;j++)
          valp1[i1*3*gs+j] = tmp[j];

        eval_pd(tid,gs2,grid,tmp,n1,l1,m1,zeta2);

        #pragma acc parallel loop present(valp1[0:igs3],tmp[0:gs3])
        for (int j=3*gs1;j<3*gs2;j++)
          valp1[i1*gs3+j] = norm2*tmp[j];
      }
      else
      {
       //evaluates ftn and L
        eval_ss_ked(i1,gs,grid,&valL1[i1*gs],n1,l1,m1,zeta1,Rc);
      }
    }

   //contraction
    //printf("\n TESTING: ** -> * in jellium_ints! \n");
    reduce_2c1(tid,0,N,gs,valS1,valS2,N,N,S);

    if (sgs_basis)
   //PD KE contraction
    for (int i1=0;i1<N;i1++)
    for (int i2=0;i2<N;i2++)
    {
      double* valm = &valp1[i1*gs3];
      double* valn = &valp1[i2*gs3];

      double vt = 0.;
     #pragma acc parallel loop present(valm[0:gs3],valn[0:gs3],wt[0:gs]) reduction(+:vt)
      for (int j=0;j<gs;j++)
      {
        double vx = valm[3*j+0]; double wx = valn[3*j+0];
        double vy = valm[3*j+1]; double wy = valn[3*j+1];
        double vz = valm[3*j+2]; double wz = valn[3*j+2];

        vt += (vx*wx+vy*wy+vz*wz)*wt[j];
      }
      //printf("  i12: %2i %2i  vt: %8.5f \n",i1,i2,vt);

     #pragma acc serial present(T[0:N2])
      T[i1*N+i2] = -vt;
    }
    else
    {
      //printf("\n WARNING: testing fix for ** in jellium_ints! \n");
     //Laplacian KE contraction
      reduce_2c1(tid,0,N,gs,valL1,valS2,N,N,T);
    }

   #pragma acc parallel loop collapse(2) present(valS1[0:igs],Vr[0:gs])
    for (int i1=0;i1<N;i1++)
    for (int j=0;j<gs;j++)
      valS1[i1*gs+j] *= Vr[j];

   //contraction
    //printf("\n TESTING: ** -> * in jellium_ints! \n");
    reduce_2c1(tid,0,N,gs,valS1,valS2,N,N,En);

    if (prl>0) printf("  done integrating <B|O|B> \n\n");
    #pragma acc update self(S[0:N2],En[0:N2],T[0:N2])

   //cleanup
    #pragma acc exit data delete(Vr[0:gs])
    delete [] Vr;

    //#pragma acc exit data delete(valS1[0:N][0:gs],valS2[0:N][0:gs],valp1[0:N][0:gs3],valL1[0:N][0:gs],tmp[0:gs3])
    #pragma acc exit data delete(valS1[0:igs],valS2[0:igs],valp1[0:igs3],valL1[0:igs],tmp[0:gs3])
    //for (int i=0;i<N;i++)
    //  delete [] valS1[i];
    //for (int i=0;i<N;i++)
    //  delete [] valS2[i];
    //for (int i=0;i<N;i++)
    //  delete [] valp1[i];
    //for (int i=0;i<N;i++)
    //  delete [] valL1[i];
    delete [] valS1;
    delete [] valS2;
    delete [] valp1;
    delete [] valL1;
    delete [] tmp;

  }

 #if 0
  printf("\n S (before norm) \n");
  print_square(N,S);
 #endif

  {
    double norm[N];
    for (int i=0;i<N;i++)
    {
      double n1 = sqrt(S[i*N+i]);
      norm[i] = 1./n1;
    }
    if (prl>-1)
    {
      printf("  norms: ");
      for (int i=0;i<N;i++)
        printf(" %8.5f",norm[i]);
      printf("\n");
    }

    if (update_norm)
    {
      for (int i=0;i<N;i++)
        basis[i][4] = norm[i];

     //keep the norms for later
      write_vector(N,norm,"NGS");
    }

    for (int i=0;i<N;i++)
    for (int j=0;j<N;j++)
    {
      double n12 = norm[i]*norm[j];
      S[i*N+j] *= n12;
      T[i*N+j] *= -0.5*n12;
      En[i*N+j] *= n12;
    }

   #if SHIFT_EN
    printf("\n shifting diagonals of En to break degeneracy \n");
    const double sf = 1.e-6;
    //int (l=0;l<10;l++)
    for (int i=0;i<N;i++)
    //if (basis[i][1]==l)
    {
      int m = basis[i][2];
      En[i*N+i] += m*sf;
    }
   #endif

    double thresh = 1.e-12;
    for (int i=0;i<N2;i++)
    if (fabs(S[i])<thresh)
      S[i] = 0.;
    for (int i=0;i<N2;i++)
    if (fabs(T[i])<thresh)
      T[i] = 0.;
    for (int i=0;i<N2;i++)
    if (fabs(En[i])<thresh)
      En[i] = 0.;

    #pragma acc update device(S[0:N2],En[0:N2],T[0:N2])
  }

  if (prl>1 || (prl>0 && N<10))
  {
    printf("\n S: \n");
    print_square(N,S);
    printf("\n En: \n");
    print_square(N,En);
    printf("\n T: \n");
    print_square(N,T);
  }

  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc exit data delete(grid[0:gs6],wt[0:gs])
  delete [] grid;
  delete [] wt;

  return;
}

double norm_sgs(int n1, int l1, int m1, double zt1, double Rc)
{
  int nml = n1-l1-1;
  double norm1 = 1.;
  double a = zt1;
  //printf("  n1: %i l1: %i nml: %i \n",n1,l1,nml);

  if (nml==0)
    norm1 = (sqrt(PI/2.)* sqrt(((1 + 4*a*pow(Rc,2))/(exp(2*a*pow(Rc,2))*
             pow(Rc,3)) +  2*pow(a,1.5)*sqrt(2*PI)*erf(sqrt(2)*sqrt(a)*Rc))/pow(a,3)))/2.;

  else if (nml==1)
    norm1 = (sqrt(PI/2.)*sqrt(((3 + 4*a*pow(Rc,2)*(3 + 2*a*pow(Rc,2)*(3 + a*pow(Rc,2))))/(exp(2*a*pow(Rc,2))*pow(Rc,5)) +
          6*pow(a,2.5)*sqrt(2*PI)*erf(sqrt(2)*sqrt(a)*Rc))/pow(a,5)))/4.;

  else if (nml==2)
    norm1 = (sqrt(PI)*sqrt((45 + 4*a*pow(Rc,2)*(45 + 2*a*pow(Rc,2)*(45 + 2*a*pow(Rc,2)*
                  (30 + a*pow(Rc,2)*(15 + 4*a*pow(Rc,2))))))/(pow(a,7)*exp(2*a*pow(Rc,2))*pow(Rc,7)) + (60*sqrt(2*PI)*erf(sqrt(2)*sqrt(a)*Rc))/pow(a,3.5)))/16.;

  else
  { printf("\n ERROR: norm_sgs not available for n-l-1 = %i \n",nml); exit(-1); }

  norm1 = 1./norm1;

  norm1 *= norm_sh(l1,m1)/0.2820947917738782;

  return norm1;
}

void print_murak_rmax(double Rc, int nrad, vector<vector<double> > basis)
{
  double zetamin = 100.;
  for (int j=0;j<basis.size();j++)
  if (basis[j][3]<zetamin)
    zetamin = basis[j][3];
  zetamin *= Rc;

  int m = MMURAK;
  double Z1 = Z1J/10.;
  //double Z1 = zetamin;

  int gsr = nrad;
  double* rgrid = new double[gsr];
  double* rwt = new double[gsr];
  #pragma acc enter data create(rgrid[0:gsr],rwt[0:gsr])
  get_murak_grid_zeta(nrad,rgrid,rwt,Z1,m);
  #pragma acc update self(rgrid[0:gsr])

  printf("  rmax(murak): %8.5f %8.5f %8.5f %8.5f  \n",rgrid[gsr-4],rgrid[gsr-3],rgrid[gsr-2],rgrid[gsr-1]);

  #pragma acc exit data delete(rgrid[0:gsr],rwt[0:gsr])

  delete [] rgrid;
  delete [] rwt;
  return;
}

void compute_Exyz_slaussian(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g, double* ang_w, double* E, int prl, double Rc)
{
  if (prl>1) printf(" beginning compute_E_slaussian (double precision) \n");

  int N = basis.size();
  int N2 = N*N;

  printf("  compute_E. nrad/ang: %3i %3i \n",nrad,nang);

  int gs = nrad*nang;
  int gs6 = 6*gs;

  int estart = find_center_of_grid(1,nrad)*nang;

  double* grid1m = new double[gs6];
  double* grid1n = new double[gs6];
  double* wt1 = new double[gs];

  double* grid2m = new double[gs6];
  double* grid2n = new double[gs6];
  double* wt2 = new double[gs];

  double* val1m = new double[gs];
  double* val1n = new double[gs];
  double* val2m = new double[gs];
  double* val2n = new double[gs];

  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);
  //printf("  iN: %i \n",imaxN);

 #if USE_ACC
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data copyin(n2i[0:natoms])
  #pragma acc enter data copyin(coords[0:3*natoms],atno[0:natoms])

  #pragma acc enter data create(grid1m[0:gs6],grid1n[0:gs6],wt1[0:gs])
  #pragma acc enter data create(grid2m[0:gs6],grid2n[0:gs6],wt2[0:gs])
  #pragma acc enter data create(val1m[0:gs],val1n[0:gs],val2m[0:gs],val2n[0:gs])
 #endif
  acc_assign(3*N2,E,0.);

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    float Z1 = (float)atno[m];
    float A1 = coords[3*m+0]; float B1 = coords[3*m+1]; float C1 = coords[3*m+2];

    for (int i1=s1;i1<s2;i1++)
    for (int i2=s1;i2<=i1;i2++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

      vector<double> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

      float z12 = zeta1 + zeta2;
     //new grid with zeta dependence
      generate_central_grid_2d(-1,0,grid1m,wt1,z12,nrad,nang,ang_g,ang_w);

      #pragma acc parallel loop present(val1m[0:gs],wt1[0:gs])
      for (int j=0;j<gs;j++)
        val1m[j] = wt1[j];
      #pragma acc parallel loop present(val2m[0:gs])
      for (int j=0;j<gs;j++)
        val2m[j] = 1.;

     //S
      eval_ssd(ii1,gs,grid1m,val1m,n1,l1,m1,zeta1,Rc); //basis 1
      eval_ssd(ii1,gs,grid1m,val2m,n2,l2,m2,zeta2,Rc); //basis 2

      double valx = 0.; double valy = 0.; double valz = 0.;
     #pragma acc parallel loop present(val1m[0:gs],val2m[0:gs],grid1m[0:gs6]) reduction(+:valx,valy,valz)
      for (int j=0;j<gs;j++)
      {
       //assumes common 0,0,0 origin
        double x = grid1m[6*j+0]+A1;
        double y = grid1m[6*j+1]+B1;
        double z = grid1m[6*j+2]+C1;

        valx += val1m[j]*val2m[j]*x;
        valy += val1m[j]*val2m[j]*y;
        valz += val1m[j]*val2m[j]*z;
      }

     #pragma acc serial present(E[0:3*N2])
      {
        E[i1*N+i2]      = E[i2*N+i1] = valx;
        E[N2+i1*N+i2]   = E[N2+i2*N+i1] = valy;
        E[2*N2+i1*N+i2] = E[2*N2+i2*N+i1] = valz;
      }

    } //pairs of basis on single atoms


   if (natoms>1) { printf("\n WARNING: can only do 1 atom in compute_Exyz \n"); exit(-1); }
   //complete but needs testing
   #if 0
    for (int n=m+1;n<natoms;n++)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];
      //printf(" mn: %i %i s1-4: %i %i - %i %i \n",m,n,s1,s2,s3,s4);

      float Z2 = (float)atno[n];
      float A2 = coords[3*n+0]; float B2 = coords[3*n+1]; float C2 = coords[3*n+2];
      float A12 = A2-A1; float B12 = B2-B1; float C12 = C2-C1;

      for (int i1=s1;i1<s2;i1++)
      for (int i2=s3;i2<s4;i2++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

       //new grid with zeta dependence
        generate_central_grid_2d(-1,0,grid1m,wt1,zeta1,nrad,nang,ang_g,ang_w);
        generate_central_grid_2d(-1,0,grid2m,wt2,zeta2,nrad,nang,ang_g,ang_w);

       //grid1 at 0,0,0 now has r1 at 3, r2 at 4
        add_r2_to_grid(gs,grid1m,A12,B12,C12);
        recenter_grid(gs,grid2m,A12,B12,C12);

       //optimize this
        becke_weight_2d(gs,grid1m,wt1,grid2m,wt2,zeta1,zeta2,A12,B12,C12);
        //becke_weight_2d(gs,grid1m,wt1,grid2m,wt2,Z1,Z2,A12,B12,C12);

        copy_grid(gs,grid2n,grid2m);                                 ed on atom 1
        recenter_grid(gs,grid2n,-A12,-B12,-C12);      //grid 2 center
        copy_grid(gs,grid1n,grid1m);
        recenter_grid_zero(gs,grid1n,-A12,-B12,-C12); //grid 1 centered on atom 2

       //needs to happen after becke weighting
        add_r1_to_grid(gs,grid2m,0.,0.,0.);

        #pragma acc parallel loop present(val1n[0:gs],val2m[0:gs])
        for (int j=0;j<gs;j++)
          val2m[j] = val2n[j] = 1.;
        #pragma acc parallel loop present(val1m[0:gs],val2n[0:gs],wt1[0:gs],wt2[0:gs])
        for (int j=0;j<gs;j++)
        {
          val1m[j] = wt1[j];
          val1n[j] = wt2[j];
        }

       //S
        eval_ssd(ii1,gs,grid1m,val1m,n1,l1,m1,zeta1,Rc); //basis 1 on center 1
        eval_ssd(ii1,gs,grid2m,val1n,n1,l1,m1,zeta1,Rc); //basis 1 on center 2
        eval_ssd(ii1,gs,grid1n,val2m,n2,l2,m2,zeta2,Rc); //basis 2 on center 1
        eval_ssd(ii1,gs,grid2n,val2n,n2,l2,m2,zeta2,Rc); //basis 2 on center 2

        double valx = 0.; double valy = 0.; double valz = 0.;
       #pragma acc parallel loop present(val1m[0:gs],val1n[0:gs],val2m[0:gs],val2n[0:gs],grid1m[0:gs6],grid1n[0:gs6]) reduction(+:valx,valy,valz)
        for (int j=0;j<gs;j++)
        {
          double x1 = grid1m[6*j+0]+A1;
          double y1 = grid1m[6*j+1]+B1;
          double z1 = grid1m[6*j+2]+C1;
          double x2 = grid2n[6*j+0]+A2;
          double y2 = grid2n[6*j+1]+B2;
          double z2 = grid2n[6*j+2]+C2;
          valx += val1m[j]*val2m[j]*x1 + val1n[j]*val2n[j]*x2;
          valy += val1m[j]*val2m[j]*y1 + val1n[j]*val2n[j]*y2;
          valz += val1m[j]*val2m[j]*z1 + val1n[j]*val2n[j]*z2;
        }

       #pragma acc serial present(E[0:N2])
        {
          E[i1*N+i2]      = E[i2*N+i1] = valx;
          E[N2+i1*N+i2]   = E[N2+i2*N+i1] = valy;
          E[2*N2+i1*N+i2] = E[2*N2+i2*N+i1] = valz;
        }
      }

    } //loop n>m
   #endif

  } //loop m over natoms


  double* norm = new double[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];
  #pragma acc enter data copyin(norm[0:N])

  #pragma acc parallel loop independent present(E[0:3*N2],norm[0:N])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=0;j<N;j++)
  {
    double n12 = norm[i]*norm[j];
    E[i*N+j]      *= n12;
    E[N2+i*N+j]   *= n12;
    E[2*N2+i*N+j] *= n12;
  }

 #if 0
 #pragma acc parallel loop independent present(S[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=0;j<i;j++)
  {
    E[j*N+i] = S[i*N+j];
  }
 #endif

  #pragma acc exit data delete(norm[0:N])
  delete [] norm;

  //clean_small_values(N,E);


 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc update self(E[0:3*N2])
 #endif


  if (prl>1 || prl==-1)
  {
    printf("\n Exyz: \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %15.12f",E[i*N+j]);
      printf("\n");
    }
  }

#if USE_ACC
  #pragma acc exit data delete(grid1m[0:gs6],grid1n[0:gs6],wt1[0:gs])
  #pragma acc exit data delete(grid2m[0:gs6],grid2n[0:gs6],wt2[0:gs])
  #pragma acc exit data delete(val1m[0:gs],val1n[0:gs],val2m[0:gs],val2n[0:gs])
  #pragma acc exit data delete(n2i[0:natoms])
  #pragma acc exit data delete(coords[0:3*natoms],atno[0:natoms])
#endif

  delete [] n2i;

  delete [] val1m;
  delete [] val1n;
  delete [] val2m;
  delete [] val2n;
  delete [] grid1m;
  delete [] grid1n;
  delete [] grid2m;
  delete [] grid2n;
  delete [] wt1;
  delete [] wt2;

  return;
}
