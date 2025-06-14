#include "gauss.h"

#define RHO_ZERO 1.e-12
#define RDEN 1.e-20
#define pr_inc 65



//eval_gh: check factor of sqrt(2)?




void integrate_hole_para_gh(double* rdm, bool full_rdm, bool hfx_on, int Nc, int No, int M, int natoms, int* atno, double* coords, int gs, int gsb, vector<vector<double> >& basis, 
                    double* Pao, double* Pmo, double* jCA, float* grid, float* gridb, float* wt, float* wtb, float* rho, float* vxch, int prl)
{
  #define OUTER_DIVIDE 1

  int ngpu = acc_get_num_devices(acc_device_nvidia);
  if (ngpu<1) ngpu = 1;

  printf("\n creating xc hole potential in Gaussian basis (ngpu: %i) \n",ngpu);
  fflush(stdout);

  if (Nc>0)
  {
    printf("\n ERROR: Nc must be zero \n");
    exit(1);
  }

  int N = basis.size();
  int N2 = N*N;

  int M2 = M*M;
  int M3 = M*M2;
  int M4 = M2*M2;

  float* Pmof = new float[N2];
  for (int m=0;m<N2;m++)
    Pmof[m] = Pmo[m];

  if (M<8 && prl>-2)
  {
    printf("\n jCA: \n");
    print_square(N,jCA);
    printf("\n Pao: \n");
    print_square(N,Pao);
    printf("\n Pmo: \n");
    print_square(N,Pmo);
    print_rdm(M,rdm);
  }

  if (!full_rdm)
  {
   //removes one operation from within loop
    for (int p=0;p<M;p++)
    for (int q=0;q<M;q++)
    for (int r=0;r<M;r++)
    for (int s=0;s<M;s++)
      rdm[p*M3+q*M2+r*M+s] -= Pmo[p*N+q]*Pmo[r*N+s];
  }

  float* rdmf = new float[M4];
  for (int m=0;m<M4;m++)
    rdmf[m] = rdm[m];

  int gsa = natoms*gs;
  int gsa6 = 6*gsa;
  int gsba = natoms*gsb;
  int gsba6 = 6*gsba;

  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);

  int iN = N;

 //tmp arrays for summing over contracted Gaussians
  float* val0 = new float[gsa];
  float* val0b = new float[gsba];

 // 1 --> AO
 // 2 --> MO
  float** val1 = new float*[iN];
  for (int i=0;i<iN;i++)
    val1[i] = new float[gsa];
  float** val2 = new float*[iN];
  for (int i=0;i<iN;i++)
    val2[i] = new float[gsa];

  float** val1b = new float*[iN];
  for (int i=0;i<iN;i++)
    val1b[i] = new float[gsba];
  float** val2b = new float*[iN];
  for (int i=0;i<iN;i++)
    val2b[i] = new float[gsba];

  float* grid1 = new float[gsa6];
  float* grid1b = new float[gsba6];

  float* vxcht = new float[gsa];
  float* rdm2 = new float[M2];

  #pragma omp parallel for schedule(static) num_threads(ngpu)
  for (int g=0;g<ngpu;g++)
  {
    acc_set_device_num(g,acc_device_nvidia);

    if (g>0)
    {
      #pragma acc enter data create(grid[0:gsa6])
    }
    #pragma acc enter data create(vxcht[0:gsa])

    #pragma acc enter data copyin(rdmf[0:M4],Pmof[0:N2])
    #pragma acc enter data create(rdm2[0:M2])

    #pragma acc enter data copyin(jCA[0:N2],Pmo[0:N2])
    #pragma acc enter data create(val0[0:gsa],val0b[0:gsba],val1[0:iN][0:gsa],val2[0:iN][0:gsa],val1b[0:iN][0:gsba],val2b[0:iN][0:gsba])

    if (g==0)
    {
      #pragma acc enter data create(grid1[0:gsa6],grid1b[0:gsba6])
    }

    #pragma acc parallel loop present(vxcht[0:gsa])
    for (int j=0;j<gsa;j++)
      vxcht[j] = 0.f;
  }
  acc_set_device_num(0,acc_device_nvidia);


  const int ig = 10; //first index to exp in Gaussian basis

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    float Z1 = (float)atno[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];

    copy_grid(gsa,grid1,grid);
    recenter_grid_zero(gsa,grid1,-A1,-B1,-C1);
    copy_grid(gsba,grid1b,gridb);
    recenter_grid_zero(gsba,grid1b,-A1,-B1,-C1);

    #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gsa])
    for (int i1=s1;i1<s2;i1++)
    for (int j=0;j<gsa;j++)
      val1[i1][j] = 0.f;

    #pragma acc parallel loop collapse(2) present(val1b[0:iN][0:gsba])
    for (int i1=s1;i1<s2;i1++)
    for (int j=0;j<gsba;j++)
      val1b[i1][j] = 0.f;

    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1;

      vector<double> basis1 = basis[i1];
      int l1 = basis1[1]; int m1 = basis1[2]; int ng = basis1[3];
      int in = ig + ng; //index to find norm

      float* valm = val1[ii1];
      float* valn = val1b[ii1];
      for (int j=0;j<ng;j++)
      {
        double zeta1 = basis1[ig+j]; double norm1 = basis1[in+j];

        eval_gh(gsa,grid1,val0,l1,m1,norm1,zeta1);

        #pragma acc parallel loop present(val0[0:gsa],valm[0:gsa])
        for (int k=0;k<gsa;k++)
          valm[k] += val0[k];

        eval_gh(gsba,grid1b,val0b,l1,m1,norm1,zeta1);

        #pragma acc parallel loop present(val0b[0:gsba],valn[0:gsba])
        for (int k=0;k<gsba;k++)
          valn[k] += val0b[k];
      }
    }
  }


 #if OUTER_DIVIDE
 //need rho
  #pragma acc parallel loop present(rho[0:gsa])
  for (int j=0;j<gsa;j++)
    rho[j] = RHO_ZERO;

  for (int i1=0;i1<N;i1++)
  {
    int ii1 = i1;
    float* valn = val1[ii1];

    for (int i2=0;i2<N;i2++)
    {
      float d1 = Pao[i1*N+i2];

      int ii2 = i2;
      float* valm = val1[ii2];

      #pragma acc parallel loop present(rho[0:gsa],valm[0:gsa],valn[0:gsa])
      for (int j=0;j<gsa;j++)
        rho[j] += d1*valn[j]*valm[j];
    }
  }
 #endif

 //AO-->MO
  #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gsa])
  for (int i1=0;i1<N;i1++)
  for (int j=0;j<gsa;j++)
    val2[i1][j] = 0.f;
  #pragma acc parallel loop collapse(2) present(val2b[0:iN][0:gsba])
  for (int i1=0;i1<N;i1++)
  for (int j=0;j<gsba;j++)
    val2b[i1][j] = 0.f;

  #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gsa],val2[0:iN][0:gsa],jCA[0:N2])
  for (int i1=0;i1<N;i1++)
  for (int j=0;j<gsa;j++)
  {
    float v2 = 0.f;
    #pragma acc loop reduction(+:v2)
    for (int k=0;k<N;k++)
      v2 += jCA[k*N+i1]*val1[k][j];
    val2[i1][j] = v2;
  }
  #pragma acc parallel loop collapse(2) present(val1b[0:iN][0:gsba],val2b[0:iN][0:gsba],jCA[0:N2])
  for (int i1=0;i1<N;i1++)
  for (int j=0;j<gsba;j++)
  {
    float v2 = 0.f;
    #pragma acc loop reduction(+:v2)
    for (int k=0;k<N;k++)
      v2 += jCA[k*N+i1]*val1b[k][j];
    val2b[i1][j] = v2;
  }

  copy_to_all_gpu(ngpu,gsa6,grid,0);
  for (int i1=0;i1<N;i1++)
    copy_to_all_gpu(ngpu,gsa,val2[i1],0);
  for (int i1=0;i1<N;i1++)
    copy_to_all_gpu(ngpu,gsba,val2b[i1],0);


  int gdiv = gsba/ngpu+1;
  double den = 0.;
  #pragma omp parallel for schedule(static) num_threads(ngpu) reduction(+:den)
  for (int g=0;g<ngpu;g++)
  {
    acc_set_device_num(g,acc_device_nvidia);
    int s1 = g*gdiv;
    int s2 = (g+1)*gdiv; if (g==ngpu-1) s2 = gsba;

    if (prl>1) printf("  s1/2: %6i %6i \n",s1,s2);
    for (int j=s1;j<s2;j++)
    {
      float x1 = gridb[6*j+0]; float y1 = gridb[6*j+1]; float z1 = gridb[6*j+2];
      float wt1 = wtb[j];

      #pragma acc parallel loop present(rdm2[0:M2])
      for (int m=0;m<M2;m++)
        rdm2[m] = 0.f;

      float dtt = 0.f;
      #pragma acc parallel loop collapse(2) present(rdm2[0:M2],rdmf[0:M4],Pmof[0:N2],val2b[0:iN][0:gsba]) reduction(+:dtt)
      for (int r=0;r<M;r++)
      for (int s=0;s<M;s++)
      {
        float dt = 0.f;
        //float p1 = Pmof[r*N+s];

        #pragma acc loop collapse(2) reduction(+:dt)
        for (int p=0;p<M;p++)
        for (int q=0;q<M;q++)
        {
          float d0 = rdmf[p*M3+q*M2+r*M+s];
          d0 *= val2b[p][j]*val2b[q][j];
          dt += d0;
        }
        rdm2[r*M+s] = dt;
        dtt += dt;
      }

      double d1b = 0.;
      #pragma acc parallel loop collapse(2) present(Pmof[0:N2],val2b[0:iN][0:gsba]) reduction(+:d1b)
      for (int p=0;p<N;p++)
      for (int q=0;q<N;q++)
      {
        float d0 = Pmof[p*N+q];
        d0 *= val2b[p][j]*val2b[q][j];

        d1b += d0;
      }
      den += d1b*wt1;

      //if (fabs(dtt)>1.e-8 && fabs(d1b)>1.e-8)
      #pragma acc parallel loop present(vxcht[0:gsa],Pmo[0:N2],Pmof[0:N2],rdm2[0:M2],grid[0:gsa6],val2[0:iN][0:gsa],val2b[0:iN][0:gsba])
      for (int k=0;k<gsa;k++)
      {
        float x12 = grid[6*k+0]-x1; float y12 = grid[6*k+1]-y1; float z12 = grid[6*k+2]-z1;

        float rn = sqrtf(x12*x12+y12*y12+z12*z12)+RDEN;
        float or1 = 1.f/rn;

       #if !OUTER_DIVIDE
        double d1 = RHO_ZERO;
        #pragma acc loop collapse(2) reduction(+:d1)
        for (int p=0;p<N;p++)
        for (int q=0;q<N;q++)
        {
          float d0 = Pmof[p*N+q];
          d0 *= val2[p][k]*val2[q][k];

          d1 += d0;
        }
       #endif

        double d2 = 0.;
        //if (fabs(d1)>1.e-8)
        #pragma acc loop collapse(2) reduction(+:d2)
        for (int r=0;r<M;r++)
        for (int s=0;s<M;s++)
        {
          float d0 = rdm2[r*M+s];
          d0 *= val2[r][k]*val2[s][k];
          d2 += d0;
        }
       #if !OUTER_DIVIDE
        d2 /= d1;
       #endif

        vxcht[k] += d2*wt1*or1;

      } //loop k (on gpu)
    } //loop j

  } //loop g (omp)

  for (int j=0;j<gsa;j++)
    vxch[j] = 0.;
  //int gdiv2 = gsa/ngpu+1;
  for (int g=0;g<ngpu;g++)
  {
    acc_set_device_num(g,acc_device_nvidia);
    #pragma acc update self(vxcht[0:gsa])

    //#pragma omp parallel for schedule(static,gdiv2) num_threads(ngpu)
    for (int k=0;k<gsa;k++)
      vxch[k] += vxcht[k];
  }
  acc_set_device_num(0,acc_device_nvidia);
  #pragma acc update device(vxch[0:gsa])

  double Pao_den = 0.;
 #if OUTER_DIVIDE
  #pragma acc parallel loop present(rho[0:gsa],wt[0:gsa]) reduction(+:Pao_den)
  for (int k=0;k<gsa;k++)
    Pao_den += rho[k]*wt[k];

  #pragma acc parallel loop present(vxch[0:gsa],rho[0:gsa])
  for (int k=0;k<gsa;k++)
    vxch[k] /= rho[k];

  #pragma acc update self(vxch[0:gsa])
 #endif


  printf("\n");
  printf("  total den (Pmo/Pao): %9.6f %9.6f \n",den,Pao_den);


  if (full_rdm)
    printf("\n vxch(full): \n");
  else
    printf("\n vxch: \n");
  for (int n=0;n<natoms;n++)
  {
    for (int m=0;m<gs;m+=pr_inc)
      printf(" %8.5f",vxch[n*gs+m]);
    printf("\n");
  }
  printf("\n");

 //cleanup
  #pragma omp parallel for schedule(static) num_threads(ngpu)
  for (int g=0;g<ngpu;g++)
  {
    acc_set_device_num(g,acc_device_nvidia);

    if (g>0)
    {
      #pragma acc exit data delete(grid[0:gsa6])
    }
    #pragma acc exit data delete(vxcht[0:gsa])
    #pragma acc exit data delete(rdmf[0:M4])
    #pragma acc exit data delete(rdm2[0:M2],Pmof[0:N2])
    #pragma acc exit data delete(jCA[0:N2],Pmo[0:N2])
    #pragma acc exit data delete(val0[0:gsa],val0b[0:gsba],val1[0:iN][0:gsa],val2[0:iN][0:gsa],val1b[0:iN][0:gsba],val2b[0:iN][0:gsba])
    if (g==0)
    {
      #pragma acc exit data delete(grid1[0:gsa6],grid1b[0:gsba6])
    }
  }

  delete [] n2i;

  delete [] grid1;
  delete [] grid1b;

  delete [] val0;
  delete [] val0b;

  for (int i=0;i<iN;i++)
    delete [] val1[i];
  delete [] val1;

  for (int i=0;i<iN;i++)
    delete [] val2[i];
  delete [] val2;

  for (int i=0;i<iN;i++)
    delete [] val1b[i];
  delete [] val1b;

  for (int i=0;i<iN;i++)
    delete [] val2b[i];
  delete [] val2b;

  delete [] Pmof;
  delete [] vxcht;

  delete [] rdmf;

  return;
}

void wf_to_grid_gh_ke_2(int natoms, int* atno, double* coords, vector<vector<double> > basis, int nbas, int nenv, int N, int* atm, int* bas, double* env,
                        double* jCA, int gs, float* grid, float* wt, double* Td, int prl)
{
  if (prl>-1) printf("  compute_ke: Gaussian basis (testing) \n");
  printf("\n WARNING: updated eval_gh_ke ftns to new procedure \n");

  int gsa = gs*natoms;
  int gsa3 = gsa*3;
  int gsa6 = gsa*6;

  int N0 = basis.size();
  if (N0!=N) { printf("\n ERROR: mismatch in basis sizes \n"); exit(-1); }
  //int N2 = N*N;
  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);

  printf("  imaxN: %2i \n",imaxN);

  int iN = imaxN;
  float** val1  = new float*[iN];
  float** val1d = new float*[iN];
  for (int i=0;i<iN;i++)
    val1[i]  = new float[gsa];
  for (int i=0;i<iN;i++)
    val1d[i] = new float[gsa3];

  float** val2 = new float*[iN];
  for (int i=0;i<iN;i++)
    val2[i] = new float[gsa];
  float** val2d = new float*[iN];
  for (int i=0;i<iN;i++)
    val2d[i] = new float[gsa3];

  float* valmo = new float[gsa];

  float* grid1 = new float[gsa6];
  float* grid2 = new float[gsa6];

  #pragma acc enter data create(grid1[0:gsa6],grid2[0:gsa6])
  #pragma acc enter data create(valmo[0:gsa],val1[0:iN][0:gsa],val1d[0:iN][0:gsa3],val2[0:iN][0:gsa],val2d[0:iN][0:gsa3])

  if (prl>1)
  {
    printf("\n jCA: \n");
    print_square(N,jCA);
  }

  #pragma acc parallel loop present(valmo[0:gsa])
  for (int j=0;j<gsa;j++)
    valmo[j] = 0.f;

  #pragma acc parallel loop present(Td[0:gsa])
  for (int j=0;j<gsa;j++)
    Td[j] = 0.;

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    float Z1 = (float)atno[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];

    copy_grid(gsa,grid1,grid);
    recenter_grid_zero(gsa,grid1,-A1,-B1,-C1);

   //tricky due to non-constant jumps in loop
   //#pragma omp parallel for
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      int nshl = eval_gh_full(gs,grid,&val1[ii1],i1,natoms,nbas,nenv,N,atm,bas,env);

      //eval_gh_ke_full(gs,grid,&val1d[ii1],i1,natoms,nbas,nenv,N,atm,bas,env);
      printf("  i1: %2i  nshl: %2i \n",i1,nshl);

      i1 += nshl-1;
    }
    #pragma acc update device(val1[0:iN][0:gsa])



   #if 0
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; int ng = basis1[3];
      int in = ig + ng; //index to find norm

      float* valm = val1[ii1];
      float* valL = val1L[ii1];
      for (int j=0;j<ng;j++)
      {
        double zeta1 = basis1[ig+j]; double norm1 = basis1[in+j];

        eval_gh(gsa,grid1,val0,l1,m1,norm1,zeta1);

       #pragma acc parallel loop present(val0[0:gsa],valm[0:gsa])
        for (int k=0;k<gsa;k++)
          valm[k] += val0[k];

        eval_gh_ke(gsa,grid1,val0,n1,l1,norm1,zeta1);

       #pragma acc parallel loop present(val0[0:gsa],valL[0:gsa])
        for (int k=0;k<gsa;k++)
          valL[k] += val0[k];
      }
    }
   #endif

   //compute all
    for (int i1=s1;i1<s2;i1++)
    {
      float* valm = val1[i1-s1];
      //float* vald = val1d[i1-s1];

      float c1 = jCA[i1*N+0];

     #pragma acc parallel loop present(valmo[0:gsa],valm[0:gsa])
      for (int j=0;j<gsa;j++)
        valmo[j] += c1*valm[j];

    // #pragma acc parallel loop present(Td[0:gsa],vald[0:gsa])
    //  for (int j=0;j<gsa;j++)
    //    TL[j] += c1*vald[j];
    }

    if (0)
    for (int n=0;n<natoms;n++)
    {
     //working on this block of the matrix
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      float Z2 = (float)atno[n];
      double A2 = coords[3*n+0]; double B2 = coords[3*n+1]; double C2 = coords[3*n+2];

      copy_grid(gsa,grid2,grid);
      recenter_grid_zero(gsa,grid1,-A2,-B2,-C2);

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;

        int nshl = eval_gh_full(gs,grid,&val2[ii2],i2,natoms,nbas,nenv,N,atm,bas,env);
        printf("  i2: %2i  nshl: %2i \n",i2,nshl);

        i2 += nshl-1;
      }
      #pragma acc update device(val2[0:iN][0:gsa])
    }

  }


  {
   //check total density
    float* rho = new float[gsa];
    #pragma acc enter data create(rho[0:gsa])

    #pragma acc parallel loop present(rho[0:gsa],valmo[0:gsa])
    for (int j=0;j<gsa;j++)
      rho[j] = valmo[j]*valmo[j];

    double d1 = 0.;
    #pragma acc parallel loop present(rho[0:gsa],wt[0:gsa]) reduction(+:d1)
    for (int j=0;j<gsa;j++)
      d1 += rho[j]*wt[j];

    printf("   total density (v2): %9.6f \n",d1);

    #pragma acc exit data delete(rho[0:gsa])
    delete [] rho;
  }

  float ket = 0.;
  #pragma acc parallel loop present(Td[0:gsa],valmo[0:gsa],wt[0:gsa]) reduction(+:ket)
  for (int j=0;j<gsa;j++)
    ket += valmo[j]*Td[j]*wt[j];
  printf("\n total KE: %9.6f \n",ket);


  float thresh = 1.e-12;
  #pragma acc parallel loop present(valmo[0:gsa])
  for (int j=0;j<gsa;j++)
  if (fabs(valmo[j])<thresh)
    valmo[j] = thresh;

  #pragma acc parallel loop present(Td[0:gsa],valmo[0:gsa])
  for (int j=0;j<gsa;j++)
    Td[j] *= -0.5/valmo[j];

  if (prl>-2)
  {
    #pragma acc update self(valmo[0:gsa])
    printf("\n MO: \n");
    print_vec(gsa,grid,valmo);
  }
  if (0)
  {
    #pragma acc update self(Td[0:gsa])
    printf("\n ke: \n");
    print_vec(gsa,grid,Td);
  }

  #pragma acc exit data delete(grid1[0:gsa6],grid2[0:gsa6])
  #pragma acc exit data delete(valmo[0:gsa],val1[0:iN][0:gsa],val1d[0:iN][0:gsa3],val2[0:iN][0:gsa],val2d[0:iN][0:gsa3])

  delete [] grid1;
  delete [] grid2;
  delete [] n2i;

  for (int i=0;i<iN;i++)
    delete [] val1[i];
  for (int i=0;i<iN;i++)
    delete [] val1d[i];
  delete [] val1;
  delete [] val1d;

  for (int i=0;i<iN;i++)
    delete [] val2[i];
  delete [] val2;
  for (int i=0;i<iN;i++)
    delete [] val2d[i];
  delete [] val2d;


  return;
}

void wf_to_grid_gh_ke(int natoms, int* atno, double* coords, vector<vector<double> > basis, double* jCA, int gs, float* grid, float* wt, double* TL, int prl)
{
  if (prl>-1) printf("  compute_ke: Gaussian basis (single orbital only) \n");
  printf("\n WARNING: updated eval_gh_ke ftns to new procedure \n");

  int gsa = gs*natoms;
  int gsa6 = gsa*6;

  int N = basis.size();
  //int N2 = N*N;
  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);

  int iN = imaxN;
  float** val1  = new float*[iN];
  float** val1L = new float*[iN];
  for (int i=0;i<iN;i++)
    val1[i]  = new float[gsa];
  for (int i=0;i<iN;i++)
    val1L[i] = new float[gsa];

  float* val0 = new float[gsa];
  float* valmo = new float[gsa];

  float* grid1 = new float[gsa6];
  float* grid2 = new float[gsa6];

  #pragma acc enter data create(grid1[0:gsa6],grid2[0:gsa6])
  #pragma acc enter data create(val0[0:gsa],valmo[0:gsa],val1[0:iN][0:gsa],val1L[0:iN][0:gsa])

  if (prl>1)
  {
    printf("\n jCA: \n");
    print_square(N,jCA);
  }

  #pragma acc parallel loop present(valmo[0:gsa])
  for (int j=0;j<gsa;j++)
    valmo[j] = 0.f;

  #pragma acc parallel loop present(TL[0:gsa])
  for (int j=0;j<gsa;j++)
    TL[j] = 0.;

  const int ig = 10; //first index to exp in Gaussian basis

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    float Z1 = (float)atno[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];

    copy_grid(gsa,grid1,grid);
    recenter_grid_zero(gsa,grid1,-A1,-B1,-C1);

    #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gsa],val1L[0:iN][0:gsa])
    for (int i1=0;i1<s2-s1;i1++)
    for (int j=0;j<gsa;j++)
      val1[i1][j] = val1L[i1][j] = 0.f;

    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; int ng = basis1[3];
      int in = ig + ng; //index to find norm

      float* valm = val1[ii1];
      float* valL = val1L[ii1];
      for (int j=0;j<ng;j++)
      {
        double zeta1 = basis1[ig+j]; double norm1 = basis1[in+j];

        eval_gh(gsa,grid1,val0,l1,m1,norm1,zeta1);

       #pragma acc parallel loop present(val0[0:gsa],valm[0:gsa])
        for (int k=0;k<gsa;k++)
          valm[k] += val0[k];

        eval_gh_ke(gsa,grid1,val0,n1,l1,norm1,zeta1);

       #pragma acc parallel loop present(val0[0:gsa],valL[0:gsa])
        for (int k=0;k<gsa;k++)
          valL[k] += val0[k];
      }
    }

   //compute all
    for (int i1=s1;i1<s2;i1++)
    {
      float* valm = val1[i1-s1];
      float* valL = val1L[i1-s1];

      float c1 = jCA[i1*N+0];

     #pragma acc parallel loop present(valmo[0:gsa],valm[0:gsa])
      for (int j=0;j<gsa;j++)
        valmo[j] += c1*valm[j];

     #pragma acc parallel loop present(TL[0:gsa],valL[0:gsa])
      for (int j=0;j<gsa;j++)
        TL[j] += c1*valL[j];
    }
  }


  {
   //check total density
    float* rho = new float[gsa];
    #pragma acc enter data create(rho[0:gsa])

    #pragma acc parallel loop present(rho[0:gsa],valmo[0:gsa])
    for (int j=0;j<gsa;j++)
      rho[j] = valmo[j]*valmo[j];

    double d1 = 0.;
    #pragma acc parallel loop present(rho[0:gsa],wt[0:gsa]) reduction(+:d1)
    for (int j=0;j<gsa;j++)
      d1 += rho[j]*wt[j];

    printf("   total density (v1): %9.6f \n",d1);

    #pragma acc exit data delete(rho[0:gsa])
    delete [] rho;
  }


  float ket = 0.;
  #pragma acc parallel loop present(TL[0:gsa],valmo[0:gsa],wt[0:gsa]) reduction(+:ket)
  for (int j=0;j<gsa;j++)
    ket += valmo[j]*TL[j]*wt[j];
  printf("\n total KE: %9.6f \n",ket);


  float thresh = 1.e-12;
  #pragma acc parallel loop present(valmo[0:gsa])
  for (int j=0;j<gsa;j++)
  if (fabs(valmo[j])<thresh)
    valmo[j] = thresh;

  #pragma acc parallel loop present(TL[0:gsa],valmo[0:gsa])
  for (int j=0;j<gsa;j++)
    TL[j] *= -0.5/valmo[j];

  if (prl>-2)
  {
    #pragma acc update self(valmo[0:gsa])
    printf("\n MO: \n");
    print_vec(gsa,grid,valmo);
  }
  if (prl>2)
  {
    #pragma acc update self(TL[0:gsa])
    printf("\n ke: \n");
    print_vec(gsa,grid,TL);
  }

  #pragma acc exit data delete(grid1[0:gsa6],grid2[0:gsa6])
  #pragma acc exit data delete(val0[0:gsa],valmo[0:gsa],val1[0:iN][0:gsa],val1L[0:iN][0:gsa])

  delete [] grid1;
  delete [] grid2;
  delete [] n2i;

  delete [] val0;
  for (int i=0;i<iN;i++)
    delete [] val1[i];
  for (int i=0;i<iN;i++)
    delete [] val1L[i];
  delete [] val1;
  delete [] val1L;


  return;
}

void eval_gh_ke(int gs, double* grid, double* val1, int n1, int l1, double norm1, const double zeta1)
{
  int gs6 = 6*gs;

  int nnm = n1*(n1-1);
  int llp = l1*(l1+1);
  double fz2 = 4.*zeta1*zeta1;
  double ta = 2.*zeta1*(1.+2.*n1);

  #pragma acc parallel loop present(val1[0:gs],grid[0:gs6])
  for (int j=0;j<gs;j++)
  {
    double r = grid[6*j+3];
    double r2 = r*r;
    val1[j] = -ta + (nnm-llp)/r2 + fz2*r2;
  }

  if (norm1!=1.)
  #pragma acc parallel loop present(val1[0:gs])
  for (int j=0;j<gs;j++)
    val1[j] *= norm1;

  return;
}

void eval_gh_ke(int gs, float* grid, float* val1, int n1, int l1, float norm1, const float zeta1)
{
  int gs6 = 6*gs;

  int nnm = n1*(n1-1);
  int llp = l1*(l1+1);
  float fz2 = 4.*zeta1*zeta1;
  float ta = 2.*zeta1*(1+2*n1);
 #pragma acc parallel loop present(val1[0:gs],grid[0:gs6])
  for (int j=0;j<gs;j++)
  {
    float r = grid[6*j+3];
    float r2 = r*r;
    val1[j] = -ta + (nnm-llp)/r2 + fz2;
  }

  #pragma acc parallel loop present(val1[0:gs])
  for (int j=0;j<gs;j++)
    val1[j] *= norm1;

  return;
}

int eval_gh_full(int gs, float* grid, float** val1, int i1, int natoms, int nbas, int nenv, int N, int* atm, int* bas, double* env)
{
  int shl_size = 0;
  printf("\n ERROR: eval_gh_full not available \n");
  //gen_gto_on_grid<float>(nbas,natoms,nenv,bas,atm,env,
  //    grid,gs,N,i1,val1,shl_size);

  if (shl_size<1) { printf("\n ERROR: shl_size must be >0 (eval_gh_full) \n"); exit(-1); }

  return shl_size;
}

void eval_pd_gh(int gs, double* grid, double* val, int n1, int l1, int m1, double norm1, double zeta1)
{
  int gs3 = 3*gs;
  int nlm = n1-l1-1;
  double ntm = 0.5*nlm-1.;
  double tz = 2.*zeta1;
  double hz = 3.*zeta1;
  double fz = 4.*zeta1;
  double iz = 5.*zeta1;
  double sz = 6.*zeta1;
  double ez = 8.*zeta1;

  if (l1==0)
  {
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);
      double v1 = nlm-2.*zeta1*r2;

      val[3*i+0] = g1*rp*v1*x;
      val[3*i+1] = g1*rp*v1*y;
      val[3*i+2] = g1*rp*v1*z;
    }
  }
  else if (l1==1)
  {
    if (m1==1) //x
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);
      double v0 = nlm-2.*zeta1*r2;
      double v1 = 1.+v0;

      val[3*i+0] = g1*rp*(y2+z2+x2*v1);
      val[3*i+1] = g1*rp*v0*x*y;
      val[3*i+2] = g1*rp*v0*x*z;
    }
    else if (m1==-1) //y
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);
      double v0 = nlm-2.*zeta1*r2;
      double v1 = 1.+v0;

      val[3*i+0] = g1*rp*v0*x*y;
      val[3*i+1] = g1*rp*(x2+z2+y2*v1);
      val[3*i+2] = g1*rp*v0*y*z;
    }
    else if (m1==0) //z
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);
      double v0 = nlm-2.*zeta1*r2;
      double v1 = 1.+v0;

      val[3*i+0] = g1*rp*v0*x*z;
      val[3*i+1] = g1*rp*v0*y*z;
      val[3*i+2] = g1*rp*(x2+y2+z2*v1);
    }
  }
  else if (l1==2)
  {
    if (m1==-2) //xy
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);
      double v0 = nlm-2.*zeta1*r2;
      double v1 = 1.+v0;

      val[3*i+0] = g1*rp*(y2+z2+x2*v1)*y;
      val[3*i+1] = g1*rp*(x2+z2+y2*v1)*x;
      val[3*i+2] = g1*rp*v0*x*y*z;
    }
    else if (m1==-1) //yz
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);
      double v0 = nlm-2.*zeta1*r2;
      double v1 = 1.+v0;

      val[3*i+0] = g1*rp*v0*x*y*z;
      val[3*i+1] = g1*rp*(x2+z2+y2*v1)*z;
      val[3*i+2] = g1*rp*(x2+y2+z2*v1)*y;
    }
    else if (m1==0) //z2
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;
      double z4 = z2*z2;
      double x2y2 = x2+y2;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);
      double v1 = -2.-nlm+2.*zeta1*x2y2;
      double v2 = 1.-nlm+zeta1*x2y2;
      double v3 = 4.-nlm+2.*zeta1*x2y2;
      double v4 = 2.+nlm-zeta1*x2y2;

      val[3*i+0] = g1*rp*(x2y2*v1 - 2.*v2*z2 - fz*z4)*x;
      val[3*i+1] = g1*rp*(x2y2*v1 - 2.*v2*z2 - fz*z4)*y;
      val[3*i+2] = g1*rp*(x2y2*v3 + 2.*v4*z2 - fz*z4)*z;
    }
    else if (m1==1) //xz
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);
      double v0 = nlm-2.*zeta1*r2;
      double v1 = 1.+v0;

      val[3*i+0] = g1*rp*(y2+z2+x2*v1)*z;
      val[3*i+1] = g1*rp*v0*x*y*z;
      val[3*i+2] = g1*rp*(x2+y2+z2*v1)*x;
    }
    else if (m1==2) //x2-y2
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);
      double y2z2 = y2+z2;
      double v0 = 2.+nlm-tz*z2;
      double v1 = 2.-nlm+tz*y2z2;

      val[3*i+0] = g1*rp*(-tz*x2*x2 + 2.*z2 + x2*v0 + y2*v1)*x;
      val[3*i+1] = g1*rp*(-tz*x2*x2 - 2.*z2 + x2*(v0-4.) + y2*(v1-4.))*y;
      val[3*i+2] = g1*rp*(nlm-tz*r2)*(x-y)*(x+y)*z;
    }
  }
  else if (l1==3)
  {
    if (m1==-3) //(3x2-y2)y
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;
      double x4 = x2*x2;
      double y4 = y2*y2;
      double z4 = z2*z2;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);

      val[3*i+0] = g1*rp*x*y*(-sz*x4 + 6.*z2 + x2*(6. + 3.*nlm - fz*y2 - sz*z2) +  y2*(6. - nlm + tz*(y2 + z2)));
      val[3*i+1] = g1*rp*(tz*y2*y4 + x4*(3. - sz*y2) - 3.*y2*z2 - y4*(3. + nlm - tz*z2) + x2*(3.*z2 + y2*(3.*nlm - fz*y2 - sz*z2)));
      val[3*i+2] = g1*rp*y*(-3.*x2 + y2)*z*(-nlm + tz*r2);
      //val[3*i+0] = g1*rp*(-sz*x4 + 6.*z2 + x2*(6. + 3.*nlm - fz*y2 - sz*z2) + y2*(6. - nlm + tz*(y2+z2)))*x*y;
      //val[3*i+1] = g1*rp*(tz*y2*y4 + x4*(3. - sz*y2) - 3.*y2*z2 - y4*(3. + nlm - tz*z2) + x2*(3.*z2 + y2*(3.*nlm - fz*y2 - sz*z2)));
      //val[3*i+2] = g1*rp*y*(-3.*x2 + y2)*z*(-nlm + tz*r2)*y*z;
    }
    else if (m1==-2) //xyz
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;
      double x4 = x2*x2;
      double y4 = y2*y2;
      double z4 = z2*z2;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);
      double x2y2 = x2+y2;
      double y2z2 = y2+z2;

      val[3*i+0] = g1*rp*(y2 + z2 + x2*(1. + nlm - tz*r2))*y*z;
      val[3*i+1] = g1*rp*(x2*(1. - tz*y2) + z2 + y2*(1. + nlm - tz*y2z2))*x*z;
      val[3*i+2] = g1*rp*(x2 + y2 + (1. + nlm - tz*x2y2)*z2 - tz*z4)*x*y;
    }
    else if (m1==-1) //(4z2-x2-y2)y
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;
      double x4 = x2*x2;
      double y4 = y2*y2;
      double z4 = z2*z2;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);
      double x2y2 = x2+y2;

      val[3*i+0] = g1*rp*(x2y2*(-2. - nlm + tz*x2y2) - 2.*(1. - 2.*nlm + hz*x2y2)*z2 - ez*z4)*x*y;
      val[3*i+1] = g1*rp*(x2y2*(-x2 - (3. + nlm - tz*x2)*y2 + tz*y4) + (3.*x2 + (1. + 4.*nlm - sz*x2)*y2 - sz*y4)*z2 +  4.*(1. - tz*y2)*z4);
      val[3*i+2] = g1*rp*(x2y2*(8. - nlm + tz*x2y2) + 2.*(4. + 2.*nlm - hz*x2y2)*z2 - ez*z4)*y*z;
    }
    else if (m1==0) //(2z2-3x2-3y2)z
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;
      double x4 = x2*x2;
      double y4 = y2*y2;
      double z4 = z2*z2;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);
      double x2y2 = x2+y2;

      val[3*i+0] = g1*rp*(3.*x2y2*(-2. - nlm + tz*x2y2) + 2.*(-3. + nlm + zeta1*x2y2)*z2 - fz*z4)*x*z;
      val[3*i+1] = g1*rp*(3.*x2y2*(-2. - nlm + tz*x2y2) + 2.*(-3. + nlm + zeta1*x2y2)*z2 - fz*z4)*y*z;
      val[3*i+2] = g1*rp*(-3.*x2y2*x2y2 + 3.*x2y2*(1. - nlm + tz*x2y2)*z2 + 2.*(3. + nlm + zeta1*x2y2)*z4 - fz*z4*z2);
    }
    else if (m1==1) //(4z2-x2-y2)x
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;
      double x4 = x2*x2;
      double y4 = y2*y2;
      double z4 = z2*z2;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);
      double x2y2 = x2+y2;

      val[3*i+0] = g1*rp*(x2y2*(-y2 + x2*(-3. - nlm + tz*x2y2)) + (3.*y2 + x2*(1. + 4.*nlm - sz*x2y2))*z2 + 4.*(1. - tz*x2)*z4);
      val[3*i+1] = g1*rp*(x2y2*(-2. - nlm + tz*x2y2) - 2.*(1. - 2.*nlm + hz*x2y2)*z2 - ez*z4)*x*y;
      val[3*i+2] = g1*rp*(x2y2*(8. - nlm + tz*x2y2) + 2.*(4. + 2.*nlm - hz*x2y2)*z2 - ez*z4)*x*z;
    }
    else if (m1==2) //(x2-y2)z
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;
      double x4 = x2*x2;
      double y4 = y2*y2;
      double z4 = z2*z2;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);
      double x2y2 = x2+y2;
      double y2z2 = y2+z2;

      val[3*i+0] = g1*rp*(-tz*x4 + 2.*z2 + x2*(2. + nlm - tz*z2) +  y2*(2. - nlm + tz*y2z2))*x*z;
      val[3*i+1] = g1*rp*(-tz*x4 - 2.*z2 + x2*(-2. + nlm - tz*z2) +  y2*(-2. - nlm + tz*y2z2))*y*z;
      val[3*i+2] = g1*rp*(x-y)*(x+y)*(-x2 - y2 - (1. + nlm - tz*x2y2)*z2 + tz*z4);
    }
    else if (m1==3) //(x2-3y2)x
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;
      double x4 = x2*x2;
      double y4 = y2*y2;
      double z4 = z2*z2;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);
      double y2z2 = y2+z2;

      val[3*i+0] = g1*rp*(-tz*x4*x2 - 3.*y2*y2z2 + x4*(3. + nlm + fz*y2 - tz*z2) + 3.*x2*(z2 + y2*(-nlm + tz*y2z2)));
      val[3*i+1] = g1*rp*(-tz*x4 + sz*y4 - 6.*z2 - 3.*y2*(2. + nlm - tz*z2) + x2*(-6. + nlm + fz*y2 - tz*z2))*x*y;
      val[3*i+2] = g1*rp*x*(x2 - 3.*y2)*z*(nlm - tz*r2);
    }
  }
  else if (l1==4)
  {
    if (m1==-4) //x*y * (x*x - y*y)
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;
      double x4 = x2*x2;
      double y4 = y2*y2;
      double z4 = z2*z2;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);

      val[3*i+0] = g1*rp*y*(-tz*x4*x2 - y2*(y2 + z2) + x4*(3. + nlm - tz*z2) +  x2*(3.*z2 + y2*(2. - nlm + tz*(y2 + z2))));
      val[3*i+1] = g1*rp*x*(tz*y4*y2 + x4*(1. - tz*y2) - 3.*y2*z2 - y4*(3. + nlm - tz*z2) + x2*(z2 + y2*(-2. + nlm - tz*z2)));
      val[3*i+2] = g1*rp*x*(x-y)*y*(x+y)*z*(nlm - tz*r2);
    }
    else if (m1==-3) //y*z * (3.*x*x - y*y)
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;
      double x4 = x2*x2;
      double y4 = y2*y2;
      double z4 = z2*z2;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);

      val[3*i+0] = g1*rp*x*y*z*(-sz*x4 + tz*y4 + 6.*z2 + x2*(6. + 3.*nlm - fz*y2 - sz*z2) + y2*(6. - nlm + tz*z2));
      val[3*i+1] = g1*rp*z*(tz*y4*y2 + x4*(3. - sz*y2) - 3.*y2*z2 +  y4*(-3. - nlm + tz*z2) +  x2*(3.*z2 + y2*(3.*nlm - fz*y2 - sz*z2)));
      val[3*i+2] = g1*rp*y*(x4*(3. - sz*z2) + x2*(2.*y2 + (3. + 3.*nlm - fz*y2)*z2 - sz*z4) + y2*(-y2 + (-1. - nlm + tz*y2)*z2 + tz*z4));
    }
    else if (m1==-2) //x*y * (6.*z*z - x*x - y*y)
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;
      double x4 = x2*x2;
      double y4 = y2*y2;
      double z4 = z2*z2;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);

      val[3*i+0] = g1*rp*2.*y*(zeta1*x4*x2 - 0.5*y4 + 2.5*y2*z2 + 3.*z4 + x4*(-1.5 - 0.5*nlm + tz*y2 - iz*z2) + x2*(zeta1*y4 + z2*(1.5 + 3.*nlm - sz*z2) +  y2*(-2. - 0.5*nlm - iz*z2)));
      val[3*i+1] = g1*rp*2.*x*(zeta1*y4*y2 + x4*(-0.5 + zeta1*y2) + 3.*z4 + y2*z2*(1.5 + 3.*nlm - sz*z2) + y4*(-1.5 - 0.5*nlm - iz*z2) +  x2*(tz*y4 + 2.5*z2 + y2*(-2. - 0.5*nlm - iz*z2)));
      val[3*i+2] = g1*rp*2.*x*y*z*(zeta1*x4 + zeta1*y4 + z2*(6. + 3.*nlm - sz*z2) +  y2*(6. - 0.5*nlm - iz*z2) + x2*(6. - 0.5*nlm + tz*y2 - iz*z2));
    }
    else if (m1==-1) //y*z * (4.*z*z - 3.*x*x - 3.*y*y)
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;
      double x4 = x2*x2;
      double y4 = y2*y2;
      double z4 = z2*z2;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);

      val[3*i+0] = g1*rp*x*y*z*(sz*x4 + sz*y4 + z2*(-6. + 4.*nlm - ez*z2) + y2*(-6. - 3.*nlm - tz*z2) + x2*(-6. - 3.*nlm + 12.*zeta1*y2 - tz*z2));
      val[3*i+1] = g1*rp*z*(sz*y4*y2 + x4*(-3. + sz*y2) + 4.*z4 + y2*z2*(-5. + 4.*nlm - ez*z2) + y4*(-9. - 3.*nlm - tz*z2) + x2*(12.*zeta1*y4 + 1.*z2 + y2*(-12. - 3.*nlm - tz*z2)));
      val[3*i+2] = g1*rp*y*(-3.*y4 + y2*(9. - 3.*nlm + sz*y2)*z2 + (12. + 4.*nlm - tz*y2)*z4 - ez*z4*z2 + x4*(-3. + sz*z2) + x2*(-6.*y2 + (9. - 3.*nlm + 12.*zeta1*y2)*z2 - tz*z4));
    }
    else if (m1==0) //(35.*z2*z2 - 30.*z2*r2 + 3.*r2*r2)
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;
      double x4 = x2*x2;
      double y4 = y2*y2;
      double z4 = z2*z2;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);

      val[3*i+0] = g1*rp*x*(-sz*x4*x2 - sz*y4*y2 + z4*(-48. + 8.*nlm - 16.*zeta1*z2) + y2*z2*(-36. - 24.*nlm + 32.*zeta1*z2) + y4*(12. + 3.*nlm + 42.*zeta1*z2) + x4*(12. + 3.*nlm - 18.*zeta1*y2 + 42.*zeta1*z2) + x2*(-18.*zeta1*y4 + z2*(-36. - 24.*nlm + 32.*zeta1*z2) +  y2*(24. + 6.*nlm + 84.*zeta1*z2)));
      val[3*i+1] = g1*rp*y*(-sz*x4*x2 - sz*y4*y2 + z4*(-48. + 8.*nlm - 16.*zeta1*z2) + y2*z2*(-36. - 24.*nlm + 32.*zeta1*z2) + y4*(12. + 3.*nlm + 42.*zeta1*z2) + x4*(12. + 3.*nlm - 18.*zeta1*y2 + 42.*zeta1*z2) + x2*(-18.*zeta1*y4 + z2*(-36. - 24.*nlm + 32.*zeta1*z2) + y2*(24. + 6.*nlm + 84.*zeta1*z2)));
      val[3*i+2] = g1*rp*z*(-sz*x4*x2 - sz*y4*y2 + z4*(32. + 8.*nlm - 16.*zeta1*z2) + y2*z2*(-16. - 24.*nlm + 32.*zeta1*z2) + y4*(-48. + 3.*nlm + 42.*zeta1*z2) + x4*(-48. + 3.*nlm - 18.*zeta1*y2 + 42.*zeta1*z2) +  x2*(-18.*zeta1*y4 + z2*(-16. - 24.*nlm + 32.*zeta1*z2) +  y2*(-96. + 6.*nlm + 84.*zeta1*z2)));
    }
    else if (m1==1) //x*z * (4.*z*z - 3.*x*x - 3.*y*y)
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;
      double x4 = x2*x2;
      double y4 = y2*y2;
      double z4 = z2*z2;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);

      val[3*i+0] = g1*rp*z*(sz*x4*x2 - 3.*y4 + y2*z2 + 4.*z4 + x4*(-9. - 3.*nlm + 12.*zeta1*y2 - tz*z2) + x2*(sz*y4 + z2*(-5. + 4.*nlm - ez*z2) + y2*(-12. - 3.*nlm - tz*z2)));
      val[3*i+1] = g1*rp*x*y*z*(sz*x4 + sz*y4 + z2*(-6. + 4.*nlm - ez*z2) + y2*(-6. - 3.*nlm - tz*z2) + x2*(-6. - 3.*nlm + 12.*zeta1*y2 - tz*z2));
      val[3*i+2] = g1*rp*x*(-3.*y4 + y2*(9. - 3.*nlm + sz*y2)*z2 + (12. + 4.*nlm - tz*y2)*z4 - ez*z4*z2 + x4*(-3. + sz*z2) + x2*(-6.*y2 + (9. - 3.*nlm + 12.*zeta1*y2)*z2 - tz*z4));
      //val[3*i+0] = g1*rp*z*(sz*x4*x2 - 3.*y4 + y2*z2 + 4.*z4 + x4*(-9. - 3.*nlm + 12.*zeta1*y2 - tz*z2) + x2*(sz*y4 + z2*(-5. + 4.*nlm - ez*z2) + y2*(-12. - 3.*nlm - tz*z2)));
      //val[3*i+1] = g1*rp*x*y*z*(sz*x4 + sz*y4 + z2*(-6. + 4.*nlm - ez*z2) + y2*(-6. - 3.*nlm - tz*z2) + x2*(-6. - 3.*nlm + 12.*zeta1*y2 - tz*z2));
      //val[3*i+2] = g1*rp*x*(-3.*y4 + y2*(9. - 3.*nlm + sz*y2)*z2 + (12. + 4.*nlm - tz*y2)*z4 - ez*z4*z2 + x4*(-3. + sz*z2) + x2*(-6.*y2 + (9. - 3.*nlm + 12.*zeta1*y2)*z2 - tz*z4));
    }
    else if (m1==2) //(x2 - y2) * (6.*z*z - x2 - y2)
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;
      double x4 = x2*x2;
      double y4 = y2*y2;
      double z4 = z2*z2;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);

      val[3*i+0] = g1*rp*2.*x*(zeta1*x4*x2 + 0.5*nlm*y4 - zeta1*y4*y2 + y2*(6. - 3.*nlm + iz*y2)*z2 + (6. + sz*y2)*z4 + x4*(-2. - 0.5*nlm + zeta1*y2 - iz*z2)+ x2*(-2.*y2 - zeta1*y4 + z2*(4. + 3.*nlm - sz*z2)));
      val[3*i+1] = g1*rp*2.*y*(2.*x2*y2 + 2.*y4 - 6.*x2*z2 - 4.*y2*z2 - 6.*z4 + nlm*(-0.5*x4 + 0.5*y4 + 3.*x2*z2 - 3.*y2*z2) + zeta1*(x4*x2 - y4*y2 + 5.*y4*z2 + 6.*y2*z4 + x4*(y2 - 5.*z2) + x2*(-1.*y4 - 6.*z4)));
      val[3*i+2] = g1*rp*2.*z*(zeta1*x4*x2 - zeta1*y4*y2 + x4*(6. - 0.5*nlm + zeta1*y2 - iz*z2) + y4*(-6. + 0.5*nlm + iz*z2) + y2*z2*(-6. - 3.*nlm + sz*z2) + x2*((6. + 3.*nlm)*z2 + zeta1*(-y4 - 6.*z4)));
    }
    else if (m1==3) //x*z * (x*x - 3.*y*y)
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;
      double x4 = x2*x2;
      double y4 = y2*y2;
      double z4 = z2*z2;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);

      val[3*i+0] = -g1*rp*2.*z*(zeta1*x4*x2 + 1.5*y4 + 1.5*y2*z2 + x4*(-1.5 - 0.5*nlm - tz*y2 + zeta1*z2) + x2*(-1.5*z2 + y2*(1.5*nlm - hz*y2 - hz*z2)));
      val[3*i+1] = -g1*rp*2.*x*y*z*(zeta1*x4 - hz*y4 + 3.*z2 + y2*(3. + 1.5*nlm - hz*z2) + x2*(3. - 0.5*nlm - tz*y2 + zeta1*z2));
      val[3*i+2] = -g1*rp*2.*x*(x4*(-0.5 + zeta1*z2) + y2*(1.5*y2 + (1.5 + 1.5*nlm - hz*y2)*z2 - hz*z4) + x2*(y2 + (-0.5 - 0.5*nlm - tz*y2)*z2 + zeta1*z4));
    }
    else if (m1==4) //(x*x * (x*x - 3.*y*y) - y*y * (3.*x*x - y*y))
    #pragma acc parallel loop present(grid[0:6*gs],val[0:gs3])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;
      double x4 = x2*x2;
      double y4 = y2*y2;
      double z4 = z2*z2;

      double g1 = exp(-zeta1*r2);
      double rp = pow(r2,ntm);

      val[3*i+0] = -g1*rp*2.*x*(zeta1*x4*x2 + zeta1*y4*y2 + 6.*y2*z2 + y4*(6. - 0.5*nlm + zeta1*z2) + x4*(-2. - 0.5*nlm - iz*y2 + zeta1*z2) + x2*(-iz*y4 - 2.*z2 + y2*(4. + 3.*nlm - sz*z2)));
      val[3*i+1] = -g1*rp*2.*y*(zeta1*x4*x2 + zeta1*y4*y2 - 2.*y2*z2 + y4*(-2. - 0.5*nlm + zeta1*z2) + x4*(6. - 0.5*nlm - iz*y2 + zeta1*z2) + x2*(-iz*y4 + 6.*z2 + y2*(4. + 3.*nlm - sz*z2)));
      val[3*i+2] = -g1*rp*2.*z*(nlm*(-0.5*x4 + 3.*x2*y2 - 0.5*y4) + zeta1*(x4*x2 + x4*(-5.*y2 + z2) + y4*(y2 + z2) + x2*(-5.*y4 - 6.*y2*z2)));
    }
  }
  else if (l1>4)
  {
    printf("\n ERROR: cannot eval_pd_gh for l1>3 \n");
  }

  return;
}

#if 0
//positive definite KEO
void eval_pdke_gh(int gs, double* grid, double* val1, int n1, int l1, int m1, double norm1, double zeta1)
{
  int nlm = n1-l1-1;
  int nm = n1-1;

 //s ftns
  if (l1==0)
  {
    #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
    for (int i=0;i<gs;i++)
    {
      double r = grid[6*i+3];
      double r2 = r*r;

      double g2 = exp(-2.*zeta1*r2);
      double rp = pow(r2,nm);
      double v1 = n1-2.*zeta1*r2;

      val1[i] = rp*g2*v1*v1;
    }
  }

 //p ftns
  if (l1==1)
  {
    if (m1==1) //x
    #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;

      double g2 = exp(-2.*zeta1*r2);
      double rp = pow(r2,nm);
      double v1 = 1.+n1-2.*zeta1*r2;

      val1[i] = rp*g2*(y2+z2+x2*v1*v1);
    }

    else if (m1==-1)//y
    #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;

      double g2 = exp(-2.*zeta1*r2);
      double rp = pow(r2,nm);
      double v1 = 1.+n1-2.*zeta1*r2;

      val1[i] = rp*g2*(x2+z2+y2*v1*v1);
    }

    else if (m1==0) //z
    #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
    for (int i=0;i<gs;i++)
    {
      double x = grid[6*i+0];
      double y = grid[6*i+1];
      double z = grid[6*i+2];
      double r = grid[6*i+3];
      double r2 = r*r;
      double x2 = x*x;
      double y2 = y*y;
      double z2 = z*z;

      double g2 = exp(-2.*zeta1*r2);
      double rp = pow(r2,nm);
      double v1 = 1.+n1-2.*zeta1*r2;

      val1[i] = rp*g2*(x2+y2+z2*v1*v1);
    }
  }

  #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
  for (int i=0;i<gs;i++)
    val1[i] *= norm1;

  return;
}
#endif

void eval_gh(int gs, float* grid, float* val1, int l1, int m1, const float norm1, const float zeta1)
{
  #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
  for (int i=0;i<gs;i++)
  {
    float r = grid[6*i+3];
    float g1 = expf(-r*r*zeta1);
    val1[i] = norm1*g1;
  }

 //most l1 values are untested
 //especially l1>=3

  if (l1==0) return;

  if (l1==1)
  {
   //this label differs from the Slater basis,
   // though the ordering is the same
    if (m1==-1)      //px
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        val1[i] *= x;
      }
    }
    else if (m1==0) //py
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float y = grid[6*i+1];
        val1[i] *= y;
      }
    }
    else           //pz
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float z = grid[6*i+2];
        val1[i] *= z;
      }
    }
  }

  else if (l1==2)
  {
    if (m1==-2)      //dxy
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        val1[i] *= x*y;
      }
    }
    else if (m1==-1)      //dyz
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        val1[i] *= y*z;
      }
    }
    else if (m1==0)      //dz2
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        val1[i] *= 3.f*z*z-r*r;
      }
    }
    else if (m1==1)      //dxz
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float z = grid[6*i+2];
        val1[i] *= x*z;
      }
    }
    else if (m1==2)      //dx2-y2
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        val1[i] *= x*x-y*y;
      }
    }
  }

  else if (l1==3)
  {
    if (m1==-3)      //fy(3x2-y2)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        val1[i] *= y*(3.f*x*x-y*y);
      }
    }
    else if (m1==-2)      //fxyz
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+1];
        val1[i] *= x*y*z;
      }
    }
    else if (m1==-1)      //fy(5z2-r2)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        val1[i] *= y*(5.f*z*z-r*r);
      }
    }
    else if (m1==0)      //f5z3-3zr2
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        val1[i] *= z*(5.f*z*z-3.f*r*r);
      }
    }
    else if (m1==1)      //fx(5z2-r2)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        val1[i] *= x*(5.f*z*z-r*r);
      }
    }
    else if (m1==2)      //f(x2-y2)z
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        val1[i] *= (x*x-y*y)*z;
      }
    }
    else if (m1==3)      //fx(x2-3y2)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        val1[i] *= x*(x*x-3.f*y*y);
      }
    }
  }

 //not tested
  else if (l1==4)
  {
    if (m1==-4)      //g x*y * (x*x - y*y)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        val1[i] *= x*y * (x*x - y*y);
      }
    }
    else if (m1==-3)      //g y*z * (3.f*x*x - y*y)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        val1[i] *= y*z * (3.f*x*x - y*y);
      }
    }
    else if (m1==-2)      //g x*y * (6.f*z*z - x*x - y*y)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        val1[i] *= x*y * (6.f*z*z - x*x - y*y);
      }
    }
    else if (m1==-1)      //g y*z * (4.f*z*z - 3.f*x*x - 3.f*y*y)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        val1[i] *= y*z * (4.f*z*z - 3.f*x*x - 3.f*y*y);
      }
    }
    else if (m1==0)       //g (35.*z2*z2 - 30.*z2*r2 + 3.*r2*r2)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        float z2 = z*z;
        float r2 = r*r;
        val1[i] *= (35.*z2*z2 - 30.*z2*r2 + 3.*r2*r2);
      }
    }
    else if (m1==1)       //g x*z * (4.*z*z - 3.f*x*x - 3.f*y*y)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        val1[i] *= x*z * (4.*z*z - 3.f*x*x - 3.f*y*y);
      }
    }
    else if (m1==2)       //g (x2 - y2) * (6.f*z*z - x2 - y2)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float x2 = x*x;
        float y2 = y*y;
        val1[i] *= (x2 - y2) * (6.f*z*z - x2 - y2);
      }
    }
    else if (m1==3)       //g x*z * (x*x - 3.f*y*y)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        val1[i] *= x*z * (x*x - 3.f*y*y);
      }
    }
    else if (m1==4)       //g (x2 * (x2 - 3.f*y2) - y2 * (3.f*x2 - y2))
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float x2 = x*x;
        float y2 = y*y;
        val1[i] *= (x2 * (x2 - 3.f*y2) - y2 * (3.f*x2 - y2));
      }
    }
  }

 //should have same order as libcint now
  else if (l1==5)
  {
    if (m1==-5)      //h (5.f*x2*x2 - 10.f*x2*y2 + y2*y2)*y
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float x2 = x*x;
        float y2 = y*y;
        val1[i] *= -(5.f*x2*x2 - 10.f*x2*y2 + y2*y2)*y;
      }
    }
    else if (m1==-4)      //h (y2 - x2)*x*y*z
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float x2 = x*x;
        float y2 = y*y;
        val1[i] *= (y2 - x2)*x*y*z;
      }
    }
   //ordering differed (lower should be +5)
    else if (m1==-3)      //h (x2*x2 - 10.f*x2*y2 + 5.f*y2*y2)*x
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        float x2 = x*x;
        float y2 = y*y;
        float z2 = z*z;
        float r2 = r*r;
        val1[i] *= y * (r2*(3.f*x2-y2) + 9.f*z2*(y2-3.f*x2));
        //val1[i] *= (x2*x2 - 10.f*x2*y2 + 5.f*y2*y2)*x;
      }
    }
   //ordering differed (lower should be +2)
    else if (m1==-2)      //h y * (y2-3.f*x2)*(x2 + y2 - 8.f*z2)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        float x2 = x*x;
        float y2 = y*y;
        float z2 = z*z;
        float r2 = r*r;
        val1[i] *= x*y*z*(3.f*z2-r2);
        //val1[i] *= y * (y2-3.f*x2)*(x2 + y2 - 8.f*z2);
      }
    }
   //ordering differed (lower should be +4)
    else if (m1==-1)      //h z * (x2*x2 - 6.f*x2*y2 + y2*y2)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        float z2 = z*z;
        float r2 = r*r;
        val1[i] *= y * (14.f*r2*z2 - 21.f*z2*z2 - r2*r2);
        //val1[i] *= z * (x2*x2 - 6.f*x2*y2 + y2*y2);
      }
    }
   //ordering differed (lower should be -2)
    else if (m1==0)      //h (x2 + y2 - 2.f*z2)*x*y*z
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        float z2 = z*z;
        float r2 = r*r;
        val1[i] *= z*(15.f*r2*r2 - 70.f*r2*z2 + 63.f*z2*z2);
        //val1[i] *= (x2 + y2 - 2.f*z2)*x*y*z;
      }
    }
   //ordering differed (lower should be +3??)
    else if (m1==1)      //h x * (x2-3.f*y2)*(x2 + y2 - 8.f*z2)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        float z2 = z*z;
        float r2 = r*r;
        val1[i] *= x*(-r2*r2 + 14.f*r2*z2 - 21.f*z2*z2);
        //val1[i] *= x * (x2-3.f*y2)*(x2 + y2 - 8.f*z2);
      }
    }
   //ordering differed (lower should be ??)
    else if (m1==2)      //h y * (x2*x2 + y2*y2 - 12.f*y2*z2 + 8.f*z2*z2 + 2.f*x2 * (y2-6.f*z2))
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        float z2 = z*z;
        float r2 = r*r;
        val1[i] *= z * (r2 - 3.f*z2)*(x-y)*(x+y);
        //val1[i] *= y * (x2*x2 + y2*y2 - 12.f*y2*z2 + 8.f*z2*z2 + 2.f*x2 * (y2-6.f*z2));
      }
    }
   //ordering differed (lower should be +2)
    else if (m1==3)      //h (x2 - y2) * (x2 + y2 - 2.f*z2) * z
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        float x2 = x*x;
        float y2 = y*y;
        float z2 = z*z;
        float r2 = r*r;
        val1[i] *= x * (r2*(x2-3.f*y2) + 9.f*z2*(3.f*y2-x2));
        //val1[i] *= (x2 - y2) * (x2 + y2 - 2.f*z2) * z;
      }
    }
   //ordering differ (lower should be ??)
    else if (m1==4)      //h z * (15.f*x2*x2 + 15.f*y2*y2 - 40.f*y2*z2 + 8.f*z2*z2 + 10.f*x2 * (3.f*y2 - 4.f*z2))
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float x2 = x*x;
        float y2 = y*y;
        val1[i] *= z * (x2*x2 - 6.f*x2*y2 + y2*y2);
        //val1[i] *= z * (15.f*x2*x2 + 15.f*y2*y2 - 40.f*y2*z2 + 8.f*z2*z2 + 10.f*x2 * (3.f*y2 - 4.f*z2));
      }
    }
   //ordering differed (lower should be ??)
    else if (m1==5)      //h x * (x2*x2 + y2*y2 - 12.f*y2*z2 + 8.f*z2*z2 + 2.f*x2* (y2 - 6.f*z2))
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float x2 = x*x;
        float y2 = y*y;
        val1[i] *= x * (10.f*x2*y2 - x2*x2 - 5.f*y2*y2);
        //val1[i] *= x * (x2*x2 + y2*y2 - 12.f*y2*z2 + 8.f*z2*z2 + 2.f*x2* (y2 - 6.f*z2));
      }
    }
  }

  else if (l1==6) //i functions
  {
    if (m1==-6)      //i xy (-3x4 + 10x2y2 - 3y4)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float x2 = x*x;
        float y2 = y*y;

        val1[i] *= x*y*(10.f*x2*y2 - 3.f*x2*x2 - 3.f*y2*y2);
      }
    }
    else if (m1==-5)      //i yz (-5x4 + 10x2y2-y4)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float x2 = x*x;
        float y2 = y*y;

        val1[i] *= y*z * (10.f*x2*y2 - 5.f*x2*x2 - y2*y2);
      }
    }
    else if (m1==-4)      //i xy (r2(x2-y2) + 11z2(y2-x2))
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+3];
        float r = grid[6*i+3];
        float x2 = x*x;
        float y2 = y*y;
        float z2 = z*z;
        float r2 = r*r;

        val1[i] *= x*y * (r2-11.f*z2)*(x2-y2);
      }
    }
    else if (m1==-3)      //i yz (3r2(3x2-y2) + 11z2(y2-3x2))
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        float x2 = x*x;
        float y2 = y*y;
        float z2 = z*z;
        float r2 = r*r;

        val1[i] *= y*z * (3.f*r2-11.f*z2)*(3.f*x2-y2);
      }
    }
    else if (m1==-2)      //i y (-r4 + 18r2z2 - 33z4)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        float z2 = z*z;
        float r2 = r*r;

        val1[i] *= y * (18.f*r2*z2 - r2*r2 - 33.f*z2*z2);
      }
    }
    else if (m1==-1)      //i yz (-5r4 + 30r2z2 - 33z4)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        float z2 = z*z;
        float r2 = r*r;

        val1[i] *= y*z * (30.f*r2*z2 - 5.f*r2*r2 - 33.f*z2*z2);
      }
    }
    else if (m1==0)      //i (-5r6 + 105r4z2 - 315r2z4 + 231z6)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        float x2 = x*x;
        float y2 = y*y;
        float z2 = z*z; float z4 = z2*z2;
        float r2 = r*r; float r4 = r2*r2;

        val1[i] *= 105.f*r4*z2 - 5.f*r2*r4 - 315.f*r2*z4 + 231.f*z4*z2;
      }
    }
    else if (m1==1)      //i xz (-5r4 + 30r2z2 - 33z4)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        float x2 = x*x;
        float y2 = y*y;
        float z2 = z*z;
        float r2 = r*r;

        val1[i] *= x*z * (30.f*r2*z2 - 5.f*r2*r2 - 33.f*z2*z2);
      }
    }
    else if (m1==2)      //i x (r4 - 18r2z2 + 33z4)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        float x2 = x*x;
        float y2 = y*y;
        float z2 = z*z;
        float r2 = r*r;

        val1[i] *= x * (r2*r2 - 18.f*r2*z2 + 33.f*z2*z2);
      }
    }
    else if (m1==3)      //i xz (3r2(x2-3y2) + 11z2(-x2+3y2))
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        float x2 = x*x;
        float y2 = y*y;
        float z2 = z*z;
        float r2 = r*r;

        val1[i] *= x*z * (3.f*r2-11.f*z2)*(x2-3.f*y2);
      }
    }
    else if (m1==4)      //i (r2(6x2y2-x4-y4) + 11z2(x4-6x2y2+y4))
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        float x2 = x*x;
        float y2 = y*y;
        float z2 = z*z;
        float r2 = r*r;

        val1[i] *= (r2-11.f*z2)*(6.f*x2*y2-x2*x2-y2*y2);
      }
    }
    else if (m1==5)      //i xz (10x2y2-x4-5y4)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        float x2 = x*x;
        float y2 = y*y;
        float z2 = z*z;
        float r2 = r*r;

        val1[i] *= x*z * (10.f*x2*y2-x2*x2-5.f*y2*y2);
      }
    }
    else if (m1==6)      //i (x6 -15x4y2 + 15x2y4 - y6)
    {
     #pragma acc parallel loop present(grid[0:6*gs],val1[0:gs])
      for (int i=0;i<gs;i++)
      {
        float x = grid[6*i+0];
        float y = grid[6*i+1];
        float z = grid[6*i+2];
        float r = grid[6*i+3];
        float x2 = x*x; float x4 = x2*x2;
        float y2 = y*y; float y4 = y2*y2;
        float z2 = z*z;
        float r2 = r*r;

        val1[i] *= x4*x2 + 15.f*x2*(y4-x2*y2) - y4*y2;
      }
    }
  }

  return;
}

