#include "integrals.h"

#include "fp_def.h"

//1. check estart
//2. g function grad fix?


#include <string>
// void auto_crash();
// void print_duration(chrono::high_resolution_clock::time_point t1, chrono::high_resolution_clock::time_point t2, string name);

#if USE_ACC
void copy_to_all_gpu2(int ngpu, int s1, FP2* A, int include_first)
{ 
  if (ngpu<1) return;
  if (include_first==0)
  {
    acc_set_device_num(0,acc_device_nvidia);
    #pragma acc update self(A[0:s1])
  }

  int start = 1; 
  if (include_first>0) start = 0;

 //#pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=start;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);
    #pragma acc update device(A[0:s1])
  }

  acc_set_device_num(0,acc_device_nvidia);
  return;
}
#endif

void gather_12_d_En_0(int s1, int s2, int gs, int iN, FP1** valtx1, FP1** valtx2, FP1** valS1x, FP1** valS2x, FP1* wtt1, FP1* wtt2)
{
  int gs3 = 3*gs;

 #pragma acc parallel loop present(wtt1[0:gs],valtx1[0:iN][0:gs3],valS1x[0:iN][0:gs3])
  for (int ii1=0;ii1<s2-s1;ii1++)
  {
   #pragma acc loop collapse(2)
    for (int j=0;j<gs;j++)
    for (int k=0;k<3;k++)
      valtx1[ii1][3*j+k] = valS1x[ii1][3*j+k]*wtt1[j];
  }

 #pragma acc parallel loop present(wtt2[0:gs],valtx2[0:iN][0:gs3],valS2x[0:iN][0:gs3])
  for (int ii1=0;ii1<s2-s1;ii1++)
  {
   #pragma acc loop collapse(2)
    for (int j=0;j<gs;j++)
    for (int k=0;k<3;k++)
      valtx2[ii1][3*j+k] = valS2x[ii1][3*j+k]*wtt2[j];
  } 
  
}

void gather_12_d_En(int s1, int s2, int gs, int iN, FP1** valS1, FP1** valS2, FP1** valtx1, FP1** valtx2, FP1* wtt1, FP1* wtt2)
{
  int gs3 = 3*gs;

 #pragma acc parallel loop present(wtt1[0:gs],valtx1[0:iN][0:gs3],valS1[0:iN][0:gs])
  for (int ii1=0;ii1<s2-s1;ii1++)
  {
   #pragma acc loop collapse(2)
    for (int j=0;j<gs;j++)
    for (int k=0;k<3;k++)
      valtx1[ii1][3*j+k] = valS1[ii1][j]*wtt1[j];
  }
   
 #pragma acc parallel loop present(wtt2[0:gs],valtx2[0:iN][0:gs3],valS2[0:iN][0:gs])
  for (int ii1=0;ii1<s2-s1;ii1++)
  {
   #pragma acc loop collapse(2)
    for (int j=0;j<gs;j++)
    for (int k=0;k<3;k++)
      valtx2[ii1][3*j+k] = valS2[ii1][j]*wtt2[j];
  }
   
  return;
}


void gather_125_d_En(int s1, int s2, int gs, int iN, FP1** valt1, FP1** valt2, FP1** valS1, FP1** valS2, FP1** valS5, FP1** valtx1, FP1** valtx2, FP1** valS1x, FP1** valS2x, FP1** valS5x, FP1* wtt1, FP1* wtt2, FP1* wt3)
{
  int gs3 = 3*gs;

  for (int ii1=0;ii1<s2-s1;ii1++)
  {
   #pragma acc parallel loop present(valt1[0:iN][0:gs],valS1[0:iN][0:gs],wtt1[0:gs])
    for (int j=0;j<gs;j++)
      valt1[ii1][j] = valS1[ii1][j]*wtt1[j];

   #pragma acc parallel loop present(valt2[0:iN][0:gs],valS2[0:iN][0:gs],wtt2[0:gs])
    for (int j=0;j<gs;j++)
      valt2[ii1][j] = valS2[ii1][j]*wtt2[j];

   #pragma acc parallel loop present(valS5[0:iN][0:gs],wt3[0:gs])
    for (int j=0;j<gs;j++)
      valS5[ii1][j] = wt3[j];

   #pragma acc parallel loop collapse(2) present(valtx1[0:iN][0:gs3],valS1x[0:iN][0:gs3],wtt1[0:gs])
    for (int j=0;j<gs;j++)
    for (int k=0;k<3;k++)
      valtx1[ii1][3*j+k] = valS1x[ii1][3*j+k]*wtt1[j];

   #pragma acc parallel loop collapse(2) present(valtx2[0:iN][0:gs3],valS2x[0:iN][0:gs3],wtt2[0:gs])
    for (int j=0;j<gs;j++)
    for (int k=0;k<3;k++)
      valtx2[ii1][3*j+k] = valS2x[ii1][3*j+k]*wtt2[j];

   #pragma acc parallel loop collapse(2) present(valS5x[0:iN][0:gs3],wt3[0:gs])
    for (int j=0;j<gs;j++)
    for (int k=0;k<3;k++)
      valS5x[ii1][3*j+k] = wt3[j];
  }
  return;
}

void gather_125_d_En_2(int s1, int s2, int gs, int iN, FP1** valtx1, FP1** valtx2, FP1** valtx3, FP1** valS1, FP1** valS2, FP1** valS5, FP1* wtt1, FP1* wtt2, FP1* wt3)
{
  int gs3 = 3*gs;

  for (int ii1=0;ii1<s2-s1;ii1++)
  {
   #pragma acc parallel loop present(valtx1[0:iN][0:gs3],valS1[0:iN][0:gs],wtt1[0:gs])
    for (int j=0;j<gs;j++)
      valtx1[ii1][3*j+0] = valtx1[ii1][3*j+1] = valtx1[ii1][3*j+2] = valS1[ii1][j]*wtt1[j];

   #pragma acc parallel loop present(valtx2[0:iN][0:gs3],valS2[0:iN][0:gs],wtt2[0:gs])
    for (int j=0;j<gs;j++)
      valtx2[ii1][3*j+0] = valtx2[ii1][3*j+1] = valtx2[ii1][3*j+2] = valS2[ii1][j]*wtt2[j];

   #pragma acc parallel loop present(valtx3[0:iN][0:gs3],valS5[0:iN][0:gs],wt3[0:gs])
    for (int j=0;j<gs;j++)
      valtx3[ii1][3*j+0] = valtx3[ii1][3*j+1] = valtx3[ii1][3*j+2] = valS5[ii1][j]*wt3[j];
  }
  return;
}


#if RED_DOUBLE
void compute_d_3c_para(int ngpu, int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* dC, FP2* xyz_grad, int prl)
#else
void compute_d_3c_para(int ngpu, int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* dC, FP1* xyz_grad, int prl)
#endif
{
 #if !USE_ACC
  printf("\n ERROR: compute_all_3c_para requires OpenACC \n");
  exit(1);
 #endif
  
  int nomp = ngpu;
 //#pragma omp parallel
 // nomp = omp_get_num_threads();

 //CPMZ this function not done
  if (prl>0) printf(" compute_d_3c_para (ngpu: %i) \n",nomp);

  int N = basis.size();
  int N2 = N*N;
  int N3 = 3*natoms;
  int Naux = basis_aux.size();
  int N2a = N2*Naux;
  int gs = nrad*nang;
  int gs3 = 3*gs;
  int gs6 = 6*gs;

  int estart = find_center_of_grid(1,nrad)*nang;

  int* na2i = new int[natoms];
  int iNa = get_imax_n2i(natoms,Naux,basis_aux,na2i);
  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);

  FP2* norms2 = new FP2[N2];
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    norms2[i*N+j] = basis[i][4]*basis[j][4];
  FP2* norms1 = new FP2[Naux];
  for (int i=0;i<Naux;i++)
    norms1[i] = norm_sv(basis_aux[i][0],basis_aux[i][1],basis_aux[i][2],basis_aux[i][3]);

  //printf("  iN/a: %i %i \n",iN,iNa);

  FP1* grid1 = new FP1[gs6];
  FP1* wt1 = new FP1[gs];

  FP1* grid2 = new FP1[gs6];
  FP1* wt2 = new FP1[gs];

  FP1* grid3 = new FP1[gs6];
  FP1* wt3 = new FP1[gs];

  FP1** val1 = new FP1*[iNa];
  FP1** val2 = new FP1*[iN];
  FP1** val3 = new FP1*[iN];
  FP1** val4 = new FP1*[iNa];
  FP1** val5 = new FP1*[iN];
  FP1** val6 = new FP1*[iN];
  FP1** val7 = new FP1*[iNa];
  FP1** val8 = new FP1*[iN];
  FP1** val9 = new FP1*[iN];

  for (int i=0;i<iNa;i++) val1[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val2[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val3[i] = new FP1[gs];
  for (int i=0;i<iNa;i++) val4[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val5[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val6[i] = new FP1[gs];
  for (int i=0;i<iNa;i++) val7[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val8[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val9[i] = new FP1[gs];

  FP1** val1x = new FP1*[iNa];
  FP1** val2x = new FP1*[iN];
  FP1** val3x = new FP1*[iN];
  FP1** val4x = new FP1*[iNa];
  FP1** val5x = new FP1*[iN];
  FP1** val6x = new FP1*[iN];
  FP1** val7x = new FP1*[iNa];
  FP1** val8x = new FP1*[iN];
  FP1** val9x = new FP1*[iN];

  for (int i=0;i<iNa;i++) val1x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++)  val2x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++)  val3x[i] = new FP1[gs3];
  for (int i=0;i<iNa;i++) val4x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++)  val5x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++)  val6x[i] = new FP1[gs3];
  for (int i=0;i<iNa;i++) val7x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++)  val8x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++)  val9x[i] = new FP1[gs3];

  FP1* valt1 = new FP1[gs3];
  FP1* valt2 = new FP1[gs3];

  FP1* grid1s = new FP1[gs6];
  FP1* grid2s = new FP1[gs6];
  FP1* grid3s = new FP1[gs6];

  FP1* grid1p = new FP1[gs6];
  FP1* grid2p = new FP1[gs6];
  FP1* grid3p = new FP1[gs6];

  FP1* wtt1 = new FP1[gs];
  FP1* wtt2 = new FP1[gs];

  FP1* ang_g = new FP1[3*nang];
  FP1* ang_w = new FP1[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

  using namespace std::chrono;
  high_resolution_clock::time_point t0 = high_resolution_clock::now();

#if USE_ACC
 #pragma omp parallel for schedule(static,1) num_threads(nomp)
#endif
  for (int n=0;n<nomp;n++)
  {
    int tid = omp_get_thread_num();
    //printf("\n about to set device: %i \n",tid); fflush(stdout);
    #if USE_ACC
    acc_set_device_num(tid,acc_device_nvidia);
    #endif

    #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
    #pragma acc enter data copyin(n2i[0:natoms],na2i[0:natoms])

    #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
    #pragma acc enter data create(grid2[0:gs6],wt2[0:gs])
    #pragma acc enter data create(grid3[0:gs6],wt3[0:gs])
    #pragma acc enter data create(val1[0:iNa][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs]) 
    #pragma acc enter data create(val4[0:iNa][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs]) 
    #pragma acc enter data create(val7[0:iNa][0:gs],val8[0:iN][0:gs],val9[0:iN][0:gs]) 

    #pragma acc enter data create(val1x[0:iNa][0:gs3],val2x[0:iN][0:gs3],val3x[0:iN][0:gs3]) 
    #pragma acc enter data create(val4x[0:iNa][0:gs3],val5x[0:iN][0:gs3],val6x[0:iN][0:gs3]) 
    #pragma acc enter data create(val7x[0:iNa][0:gs3],val8x[0:iN][0:gs3],val9x[0:iN][0:gs3]) 

    #pragma acc enter data create(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6])
    #pragma acc enter data create(grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6])

    #pragma acc enter data create(wtt1[0:gs],wtt2[0:gs])
    #pragma acc enter data create(valt1[0:gs3],valt2[0:gs3])
    #pragma acc enter data copyin(norms1[0:Naux],norms2[0:N2])
    #pragma acc enter data create(xyz_grad[0:N3])

    if (tid>0)
    {
      #pragma acc enter data create(dC[0:N2a])
    }
    acc_assign(N3,xyz_grad,0.);
  }
  #if USE_ACC
  acc_set_device_num(0,acc_device_nvidia);
  #endif

 //fixes undercounting, corrupts dC
 #pragma acc parallel loop independent present(dC[0:N2a],n2i[0:natoms],na2i[0:natoms])
  for (int m=0;m<natoms;m++)
 #pragma acc loop independent
  for (int n=0;n<natoms;n++)
  if (m!=n)
  {
    int s1 = 0; if (m>0) s1 = na2i[m-1]; int s2 = na2i[m];
    int s3 = 0; if (m>0) s3 = n2i[m-1]; int s4 = n2i[m];
    int s5 = 0; if (n>0) s5 = n2i[n-1]; int s6 = n2i[n];
   #pragma acc loop collapse(3) 
    for (int i1=s1;i1<s2;i1++)
    for (int i2=s3;i2<s4;i2++)
    for (int i3=s5;i3<s6;i3++)
      dC[i1*N2+i2*N+i3] *= 2.;
  }

 //distribute dC to all gpus
  #if USE_ACC
  copy_to_all_gpu2(nomp,N2a,dC,0);
  #endif

#if USE_ACC
 #pragma omp parallel for num_threads(nomp)
#endif
  for (int m=0;m<natoms;m++)
  {
    int tid = omp_get_thread_num();
    //printf("\n about to set device(compute): %i \n",tid); fflush(stdout);
    #if USE_ACC
    acc_set_device_num(tid,acc_device_nvidia);
    #endif

    int s1 = 0; if (m>0) s1 = na2i[m-1]; int s2 = na2i[m];
    int s3 = 0; if (m>0) s3 = n2i[m-1]; int s4 = n2i[m];
    //printf("  m: %i  s1/2->3/4: %i %i - %i %i \n",m,s1,s2,s3,s4);

    FP1 Z1 = (FP1)atno[m];
    FP1 A1 = coords[3*m+0]; FP1 B1 = coords[3*m+1]; FP1 C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);

   #pragma acc parallel loop present(val1[0:iNa][0:gs],val1x[0:iNa][0:gs3])
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val1[ii1][j] = 1.f;
     #pragma acc loop
      for (int j=0;j<gs3;j++)
        val1x[ii1][j] = 1.f;
    }

   #pragma acc parallel loop present(val2[0:iN][0:gs],val2x[0:iN][0:gs3])
    for (int ii2=0;ii2<s4-s3;ii2++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val2[ii2][j] = 1.f;
     #pragma acc loop
      for (int j=0;j<gs3;j++)
        val2x[ii2][j] = 1.f;
    }
   #pragma acc parallel loop present(val3[0:iN][0:gs],val3x[0:iN][0:gs3],wt1[0:gs])
    for (int ii3=0;ii3<s4-s3;ii3++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val3[ii3][j] = wt1[j];
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val3x[ii3][3*j+0] = val3x[ii3][3*j+1] = val3x[ii3][3*j+2] = wt1[j];
    }

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<FP2> basis1 = basis_aux[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

      acc_assign(gs,valt1,1.);
      eval_sh_3r(gs,grid1,valt1,n1,l1,m1);

      eval_inr_d(gs,grid1,val1x[ii1],n1,l1,zeta1);
     #pragma acc parallel loop present(valt1[0:gs],val1x[0:iNa][0:gs3])
      for (int j=0;j<gs;j++)
      {
        FP1 v1 = valt1[j];
        val1x[ii1][3*j+0] *= v1;
        val1x[ii1][3*j+1] *= v1;
        val1x[ii1][3*j+2] *= v1;
      }

      eval_inr_r12(gs,grid1,val1[ii1],n1,l1,zeta1);

      if (l1>0)
      {
       #pragma acc parallel loop present(val1[0:iNa][0:gs],valt1[0:gs3])
        for (int j=0;j<gs;j++)
          valt1[3*j+0] = valt1[3*j+1] = valt1[3*j+2] = val1[ii1][j];

        eval_dp_3r(gs,grid1,valt1,n1,l1,m1);
        
       #pragma acc parallel loop present(valt1[0:gs3],val1x[0:iNa][0:gs3])
        for (int j=0;j<gs3;j++)
          val1x[ii1][j] += valt1[j];
      }

      eval_sh_3r(gs,grid1,val1[ii1],n1,l1,m1);
    } //loop i1

    for (int i2=s3;i2<s4;i2++)
    {
      int ii2 = i2-s3;

      vector<FP2> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

      eval_sh(ii2,gs,grid1,val2[ii2],n2,l2,m2,zeta2);
      eval_p(gs,grid1,val2x[ii2],n2,l2,m2,zeta2);
    } //loop i2

    for (int i3=s3;i3<s4;i3++)
    {
      int ii3 = i3-s3;

      vector<FP2> basis3 = basis[i3];
      int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];

      eval_sh(ii3,gs,grid1,val3[ii3],n3,l3,m3,zeta3);
      eval_p(gs,grid1,val3x[ii3],n3,l3,m3,zeta3);
    } //loop i3

    //reduce_3c1d(m,s1,s2,s3,s4,gs,norms1,norms2,dC,val1,val2,val3,val1x,val2x,val3x,N,Naux,iN,iNa,natoms,1.,xyz_grad);


   //two-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
     //s12 over atom m, aux function
     //s34 over atom m, regular function
     //s56 over atom n, regular function
      int s5 = 0; if (n>0) s5 = n2i[n-1]; int s6 = n2i[n];
      //printf("  mn: %i %i  s: %i-%i  %i-%i %i-%i \n",m,n,s1,s2,s3,s4,s5,s6);

      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2);
      recenter_grid(gs,grid2,A12,B12,C12);

      copy_grid(gs,grid1s,grid1); 
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid2,wt2,Z1,Z2,A12,B12,C12);
      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wt2);

      add_r1_to_grid(gs,grid2,0.,0.,0.);

     #pragma acc parallel loop present(val4[0:iNa][0:gs],val4x[0:iNa][0:gs])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val4[ii1][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          val4x[ii1][j] = 1.f;
      }

     #pragma acc parallel loop present(val2[0:iN][0:gs],val5[0:iN][0:gs],val2x[0:iN][0:gs3],val5x[0:iN][0:gs3])
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val2[ii2][j] = val5[ii2][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          val2x[ii2][j] = val5x[ii2][j] = 1.f;
      }

     #pragma acc parallel loop present(val3[0:iN][0:gs],val6[0:iN][0:gs],val3x[0:iN][0:gs3],val6x[0:iN][0:gs3],wtt1[0:gs],wt2[0:gs])
      for (int ii3=0;ii3<s6-s5;ii3++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val3[ii3][j] = wtt1[j];
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val6[ii3][j] = wt2[j];
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val3x[ii3][3*j+0] = val3x[ii3][3*j+1] = val3x[ii3][3*j+2] = wtt1[j];
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val6x[ii3][3*j+0] = val6x[ii3][3*j+1] = val6x[ii3][3*j+2] = wt2[j];
      }

     //i1 on atom m
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis_aux[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];
        //printf("  m: %i i1: %i   nlm: %i %i %2i zeta: %8.5f \n",m,i1,n1,l1,m1,zeta1);

        acc_assign(gs,valt1,1.);
        eval_sh_3r(gs,grid2,valt1,n1,l1,m1);

        eval_inr_d(gs,grid2,val4x[ii1],n1,l1,zeta1);
       #pragma acc parallel loop present(valt1[0:gs],val4x[0:iNa][0:gs3])
        for (int j=0;j<gs;j++)
        {
          FP1 v1 = valt1[j];
          val4x[ii1][3*j+0] *= v1;
          val4x[ii1][3*j+1] *= v1;
          val4x[ii1][3*j+2] *= v1;
        }

        eval_inr_r12(gs,grid2,val4[ii1],n1,l1,zeta1);

        if (l1>0)
        {
         #pragma acc parallel loop present(val4[0:iNa][0:gs],valt1[0:gs3])
          for (int j=0;j<gs;j++)
            valt1[3*j+0] = valt1[3*j+1] = valt1[3*j+2] = val4[ii1][j];

          eval_dp_3r(gs,grid2,valt1,n1,l1,m1);
        
         #pragma acc parallel loop present(valt1[0:gs3],val4x[0:iNa][0:gs3])
          for (int j=0;j<gs3;j++)
            val4x[ii1][j] += valt1[j];
        }

        eval_sh_3r(gs,grid2,val4[ii1],n1,l1,m1);
      }

     //i2 on atom m
      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;

        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];
        //printf("   n: %i i2: %i   nlm: %i %i %2i zeta: %8.5f \n",n,i2,n2,l2,m2,zeta2);

        eval_sh(ii2,gs,grid1,val2[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2,val5[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid1,val2x[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid2,val5x[ii2],n2,l2,m2,zeta2);
      }

     //i3 on atom n
      for (int i3=s5;i3<s6;i3++)
      {
        int ii3 = i3-s5;

        vector<FP2> basis3 = basis[i3];
        int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];

        eval_sh(ii3,gs,grid1s,val3[ii3],n3,l3,m3,zeta3);
        eval_sh(ii3,gs,grid2s,val6[ii3],n3,l3,m3,zeta3);
        eval_p(gs,grid1s,val3x[ii3],n3,l3,m3,zeta3);
        eval_p(gs,grid2s,val6x[ii3],n3,l3,m3,zeta3);
      }

      reduce_3c2d112(m,n,s1,s2,s3,s4,s5,s6,gs,norms1,norms2,dC,val1,val2,val3,val4,val5,val6,val1x,val2x,val3x,val4x,val5x,val6x,N,Naux,iN,iNa,natoms,1.,xyz_grad);


     #pragma acc parallel loop present(val2[0:iN][0:gs],val5[0:iN][0:gs],val2x[0:iN][0:gs3],val5x[0:iN][0:gs3])
      for (int ii2=0;ii2<s6-s5;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val2[ii2][j] = val5[ii2][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          val2x[ii2][j] = val5x[ii2][j] = 1.f;
      }

     #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],val3x[0:iN][0:gs3],val6x[0:iN][0:gs3],wt1[0:gs],wt2[0:gs])
      for (int ii3=0;ii3<s6-s5;ii3++)
      {
        for (int j=0;j<gs;j++)
        {
          val3[ii3][j] = wtt1[j];
          val6[ii3][j] = wt2[j];
          val3x[ii3][3*j+0] = val3x[ii3][3*j+1] = val3x[ii3][3*j+2] = wtt1[j];
          val6x[ii3][3*j+0] = val6x[ii3][3*j+1] = val6x[ii3][3*j+2] = wt2[j];
        }
      }

     //now do i2+i3 on atom n 
     //i2 on atom n
      for (int i2=s5;i2<s6;i2++)
      {
        int ii2 = i2-s5;

        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1s,val2[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2s,val5[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid1s,val2x[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid2s,val5x[ii2],n2,l2,m2,zeta2);
      }

     //i3 on atom n
      for (int i3=s5;i3<s6;i3++)
      {
        int ii3 = i3-s5;

        vector<FP2> basis3 = basis[i3];
        int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];

        eval_sh(ii3,gs,grid1s,val3[ii3],n3,l3,m3,zeta3);
        eval_sh(ii3,gs,grid2s,val6[ii3],n3,l3,m3,zeta3);
        eval_p(gs,grid1s,val3x[ii3],n3,l3,m3,zeta3);
        eval_p(gs,grid2s,val6x[ii3],n3,l3,m3,zeta3);
      }

      reduce_3c2d122(m,n,s1,s2,s5,s6,gs,norms1,norms2,dC,val1,val2,val3,val4,val5,val6,val1x,val2x,val3x,val4x,val5x,val6x,N,Naux,iN,iNa,natoms,1.,xyz_grad);
 
    } //loop n over second atom



   //three-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1);
      recenter_grid(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      for (int p=0;p<natoms;p++)
      if (p!=m && p!=n)
      {
        int s5 = 0; if (p>0) s5 = n2i[p-1]; int s6 = n2i[p];

        FP1 Z3 = (FP1)atno[p];
        FP1 A3 = coords[3*p+0]; FP1 B3 = coords[3*p+1]; FP1 C3 = coords[3*p+2];
        FP1 A13 = A3-A1; FP1 B13 = B3-B1; FP1 C13 = C3-C1;

        generate_central_grid_2(grid3,wt3,Z3,nrad,nang,ang_g,ang_w);
        copy_grid(gs,grid3p,grid3); //grid 3 centered on atom 3
        recenter_grid(gs,grid3,A13,B13,C13); //grid 3 centered on atom 1

        copy_grid(gs,grid3s,grid3);
        recenter_grid(gs,grid3s,-A12,-B12,-C12); //grid 3 centered on atom 2

        copy_grid(gs,grid2p,grid2); //grid 2 centered on atom 3
        recenter_grid(gs,grid2p,-A13,-B13,-C13);

        copy_grid(gs,grid1p,grid1);
        recenter_grid(gs,grid1p,-A13,-B13,-C13); //grid 1 centered on atom 3

        add_r1_to_grid_6z(gs,grid1s,grid2s,grid3s,grid1p,grid2p,grid3p);

      //need to keep all of these distances in order
        add_r3_to_grid(gs,grid1,A13,B13,C13);
        add_r123_to_grid(gs,grid2,A12,B12,C12,0.,0.,0.,A13,B13,C13);
        add_r123_to_grid(gs,grid3,A13,B13,C13,A12,B12,C12,0.,0.,0.);

        acc_copyf(gs,wtt1,wt1,wtt2,wt2);

        becke_weight_3c(gs,grid1,wtt1,grid2,wtt2,grid3,wt3,Z1,Z2,Z3,A12,B12,C12,A13,B13,C13);
        eliminate_small_wt_3(estart,gs,wtt1,wtt2,wt3);

        //#pragma acc update self(grid1[0:gs6],wtt1[0:gs],grid2[0:gs6],wtt2[0:gs],grid3[0:gs6],wt3[0:gs])
        //print_grid(gs,grid1,grid2,grid3,wtt1,wtt2,wt3,prl);

        add_r1_to_grid(gs,grid2,0.,0.,0.);
        add_r1_to_grid(gs,grid3,0.,0.,0.);


       #pragma acc parallel loop present(val4[0:iNa][0:gs],val7[0:iNa][0:gs],val4x[0:iNa][0:gs3],val7x[0:iNa][0:gs3])
        for (int ii1=0;ii1<s2-s1;ii1++)
        {
         #pragma acc loop
          for (int j=0;j<gs;j++)
            val4[ii1][j] = val7[ii1][j] = 1.f;
         #pragma acc loop
          for (int j=0;j<gs3;j++)
            val4x[ii1][j] = val7x[ii1][j] = 1.f;
        }

       #pragma acc parallel loop present(val2[0:iN][0:gs],val5[0:iN][0:gs],val8[0:iN][0:gs],val2x[0:iN][0:gs3],val5x[0:iN][0:gs3],val8x[0:iN][0:gs3])
        for (int ii2=0;ii2<s4-s3;ii2++)
        {
         #pragma acc loop
          for (int j=0;j<gs;j++)
            val2[ii2][j] = val5[ii2][j] = val8[ii2][j] = 1.f;
         #pragma acc loop
          for (int j=0;j<gs3;j++)
            val2x[ii2][j] = val5x[ii2][j] = val8x[ii2][j] = 1.f;
        }

       #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],val9[0:iN][0:gs],val3x[0:iN][0:gs3],val6x[0:iN][0:gs3],val9x[0:iN][0:gs3],wtt1[0:gs],wtt2[0:gs],wt3[0:gs])
        for (int ii3=0;ii3<s6-s5;ii3++)
        {
          for (int j=0;j<gs;j++)
          {
            val3[ii3][j] = wtt1[j];
            val6[ii3][j] = wtt2[j];
            val9[ii3][j] = wt3[j];
            val3x[ii3][3*j+0] = val3x[ii3][3*j+1] = val3x[ii3][3*j+2] = wtt1[j];
            val6x[ii3][3*j+0] = val6x[ii3][3*j+1] = val6x[ii3][3*j+2] = wtt2[j];
            val9x[ii3][3*j+0] = val9x[ii3][3*j+1] = val9x[ii3][3*j+2] = wt3[j];
          }
        }

 
        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<FP2> basis1 = basis_aux[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

          acc_assign(gs,valt1,1.);
          acc_assign(gs,valt2,1.);

          eval_sh_3r(gs,grid2,valt1,n1,l1,m1);
          eval_sh_3r(gs,grid3,valt2,n1,l1,m1);

          eval_inr_d(gs,grid2,val4x[ii1],n1,l1,zeta1);
          eval_inr_d(gs,grid3,val7x[ii1],n1,l1,zeta1);

         #pragma acc parallel loop present(valt1[0:gs],val4x[0:iNa][0:gs3],valt2[0:gs],val7x[0:iNa][0:gs3])
          for (int j=0;j<gs;j++)
          {
            FP1 v1 = valt1[j];
            val4x[ii1][3*j+0] *= v1; val4x[ii1][3*j+1] *= v1; val4x[ii1][3*j+2] *= v1;
            FP1 v2 = valt2[j];
            val7x[ii1][3*j+0] *= v2; val7x[ii1][3*j+1] *= v2; val7x[ii1][3*j+2] *= v2;
          }

          eval_inr_r12(gs,grid2,val4[ii1],n1,l1,zeta1);
          eval_inr_r12(gs,grid3,val7[ii1],n1,l1,zeta1);
          if (l1>0)
          {
           #pragma acc parallel loop present(val4[0:iNa][0:gs],val7[0:iNa][0:gs],valt1[0:gs3],valt2[0:gs3])
            for (int j=0;j<gs;j++)
            {
              valt1[3*j+0] = valt1[3*j+1] = valt1[3*j+2] = val4[ii1][j];
              valt2[3*j+0] = valt2[3*j+1] = valt2[3*j+2] = val7[ii1][j];
            }

            eval_dp_3r(gs,grid2,valt1,n1,l1,m1);
            eval_dp_3r(gs,grid3,valt2,n1,l1,m1);

           #pragma acc parallel loop present(valt1[0:gs3],valt2[0:gs3],val4x[0:iNa][0:gs3],val7x[0:iNa][0:gs3])
            for (int j=0;j<gs3;j++)
            {
              val4x[ii1][j] += valt1[j];
              val7x[ii1][j] += valt2[j];
            }
          }

          eval_sh_3r(gs,grid2,val4[ii1],n1,l1,m1);
          eval_sh_3r(gs,grid3,val7[ii1],n1,l1,m1);
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<FP2> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

          eval_sh(ii2,gs,grid3s,val8[ii2],n2,l2,m2,zeta2);
          eval_sh(ii2,gs,grid2s,val5[ii2],n2,l2,m2,zeta2);
          eval_sh(ii2,gs,grid1s,val2[ii2],n2,l2,m2,zeta2);
          eval_p(gs,grid3s,val8x[ii2],n2,l2,m2,zeta2);
          eval_p(gs,grid2s,val5x[ii2],n2,l2,m2,zeta2);
          eval_p(gs,grid1s,val2x[ii2],n2,l2,m2,zeta2);
        }

        for (int i3=s5;i3<s6;i3++)
        {
          int ii3 = i3-s5;
          vector<FP2> basis3 = basis[i3];
          int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];
    
          eval_sh(ii3,gs,grid3p,val9[ii3],n3,l3,m3,zeta3);
          eval_sh(ii3,gs,grid2p,val6[ii3],n3,l3,m3,zeta3);
          eval_sh(ii3,gs,grid1p,val3[ii3],n3,l3,m3,zeta3);
          eval_p(gs,grid3p,val9x[ii3],n3,l3,m3,zeta3);
          eval_p(gs,grid2p,val6x[ii3],n3,l3,m3,zeta3);
          eval_p(gs,grid1p,val3x[ii3],n3,l3,m3,zeta3);
        }

        reduce_3c3d(m,n,p,s1,s2,s3,s4,s5,s6,gs,norms1,norms2,dC,val1,val2,val3,val4,val5,val6,val7,val8,val9,val1x,val2x,val3x,val4x,val5x,val6x,val7x,val8x,val9x,N,Naux,iN,iNa,natoms,1.,xyz_grad);

      } //loop p over third atom
    } //loop n over second atom

  } //loop m over natoms

  FP2* xyz_all = new FP2[N3]();
  for (int n=0;n<nomp;n++)
  {
    //printf("\n about to set device(collect): %i \n",n); fflush(stdout);
    #if USE_ACC
    acc_set_device_num(n,acc_device_nvidia);
    #endif

    #pragma acc update self(xyz_grad[0:N3])
    for (int i=0;i<N3;i++)
      xyz_all[i] += xyz_grad[i];

    #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
    #pragma acc exit data delete(xyz_grad[0:N3])
  }
  #if USE_ACC
  acc_set_device_num(0,acc_device_nvidia);
  #endif

  for (int i=0;i<N3;i++)
    xyz_grad[i] = xyz_all[i];

  delete [] xyz_all;

#if DEBUG
  if (prl>1)
  {
    printf("\n xyz_grad: \n");
    for (int i=0;i<natoms;i++)
    {
      for (int j=0;j<3;j++)
        printf(" %8.5f",xyz_grad[3*i+j]);
      printf("\n");
    }
  }
#endif

#if USE_ACC
 #pragma omp parallel for schedule(static,1) num_threads(nomp)
#endif
  for (int n=0;n<nomp;n++)
  {
    int tid = omp_get_thread_num();
    //printf("\n about to set device(delete): %i \n",tid); fflush(stdout);
    #if USE_ACC
    acc_set_device_num(tid,acc_device_nvidia);
    #endif

    #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
    #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs])
    #pragma acc exit data delete(grid3[0:gs6],wt3[0:gs])

    #pragma acc exit data delete(val1[0:iNa][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs]) 
    #pragma acc exit data delete(val4[0:iNa][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs]) 
    #pragma acc exit data delete(val7[0:iNa][0:gs],val8[0:iN][0:gs],val9[0:iN][0:gs]) 

    #pragma acc exit data delete(val1x[0:iNa][0:gs3],val2x[0:iN][0:gs3],val3x[0:iN][0:gs3]) 
    #pragma acc exit data delete(val4x[0:iNa][0:gs3],val5x[0:iN][0:gs3],val6x[0:iN][0:gs3]) 
    #pragma acc exit data delete(val7x[0:iNa][0:gs3],val8x[0:iN][0:gs3],val9x[0:iN][0:gs3]) 

    #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6]) 
    #pragma acc exit data delete(grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6]) 

    #pragma acc exit data delete(wtt1[0:gs],wtt2[0:gs])
    #pragma acc exit data delete(valt1[0:gs],valt2[0:gs])
    #pragma acc exit data delete(norms1[0:Naux],norms2[0:N2])
    #pragma acc exit data delete(n2i[0:natoms],na2i[0:natoms])

    if (tid>0)
    {
      #pragma acc exit data delete(dC[0:N2a])
    }
  }
  #if USE_ACC
  acc_set_device_num(0,acc_device_nvidia);
  #endif

  delete [] norms1;
  delete [] norms2;

  delete [] n2i;
  delete [] na2i;

  delete [] ang_g;
  delete [] ang_w;

  delete [] grid1s;
  delete [] grid2s;
  delete [] grid3s;
  delete [] grid1p;
  delete [] grid2p;
  delete [] grid3p;

  delete [] wtt1;
  delete [] wtt2;
  delete [] grid1;
  delete [] grid2;
  delete [] grid3;
  delete [] wt1;
  delete [] wt2;
  delete [] wt3;
  for (int i=0;i<iNa;i++) delete [] val1[i];
  for (int i=0;i<iN;i++) delete [] val2[i];
  for (int i=0;i<iN;i++) delete [] val3[i];
  for (int i=0;i<iNa;i++) delete [] val4[i];
  for (int i=0;i<iN;i++) delete [] val5[i];
  for (int i=0;i<iN;i++) delete [] val6[i];
  for (int i=0;i<iNa;i++) delete [] val7[i];
  for (int i=0;i<iN;i++) delete [] val8[i];
  for (int i=0;i<iN;i++) delete [] val9[i];
  delete [] val1;
  delete [] val2;
  delete [] val3;
  delete [] val4;
  delete [] val5;
  delete [] val6;
  delete [] val7;
  delete [] val8;
  delete [] val9;

  for (int i=0;i<iNa;i++) delete [] val1x[i];
  for (int i=0;i<iN;i++) delete [] val2x[i];
  for (int i=0;i<iN;i++) delete [] val3x[i];
  for (int i=0;i<iNa;i++) delete [] val4x[i];
  for (int i=0;i<iN;i++) delete [] val5x[i];
  for (int i=0;i<iN;i++) delete [] val6x[i];
  for (int i=0;i<iNa;i++) delete [] val7x[i];
  for (int i=0;i<iN;i++) delete [] val8x[i];
  for (int i=0;i<iN;i++) delete [] val9x[i];
  delete [] val1x;
  delete [] val2x;
  delete [] val3x;
  delete [] val4x;
  delete [] val5x;
  delete [] val6x;
  delete [] val7x;
  delete [] val8x;
  delete [] val9x;
  delete [] valt1;
  delete [] valt2;

  return;
}

void compute_d_3c_para2(
  int ngpu, int natoms, int* atno, FP1* coords,
  vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, 
  int nrad, int nang, FP2* ang_g0, FP2* ang_w0, 
  FP2* dC, FP2* xyz_grad, int prl) {
 #if !USE_ACC
  printf("\n ERROR: compute_all_3c_para requires OpenACC \n");
  exit(1);
 #endif
  
  int nomp = ngpu;
 //#pragma omp parallel
 // nomp = omp_get_num_threads();

 //CPMZ this function not done
  if (prl>0) printf(" compute_d_3c_para (ngpu: %i) \n",nomp);

  int N = basis.size();
  int N2 = N*N;
  int N3 = 3*natoms;
  int Naux = basis_aux.size();
  int N2a = N2*Naux;
  int gs = nrad*nang;
  int gs3 = 3*gs;
  int gs6 = 6*gs;

  int estart = find_center_of_grid(1,nrad)*nang;

  int* na2i = new int[natoms];
  int iNa = get_imax_n2i(natoms,Naux,basis_aux,na2i);
  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);

  FP2* norms2 = new FP2[N2];
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    norms2[i*N+j] = basis[i][4]*basis[j][4];
  FP2* norms1 = new FP2[Naux];
  for (int i=0;i<Naux;i++)
    norms1[i] = norm_sv(basis_aux[i][0],basis_aux[i][1],basis_aux[i][2],basis_aux[i][3]);

  //printf("  iN/a: %i %i \n",iN,iNa);

  FP1* grid1 = new FP1[gs6];
  FP1* wt1 = new FP1[gs];

  FP1* grid2 = new FP1[gs6];
  FP1* wt2 = new FP1[gs];

  FP1* grid3 = new FP1[gs6];
  FP1* wt3 = new FP1[gs];

  FP1** val1 = new FP1*[iNa];
  FP1** val2 = new FP1*[iN];
  FP1** val3 = new FP1*[iN];
  FP1** val4 = new FP1*[iNa];
  FP1** val5 = new FP1*[iN];
  FP1** val6 = new FP1*[iN];
  FP1** val7 = new FP1*[iNa];
  FP1** val8 = new FP1*[iN];
  FP1** val9 = new FP1*[iN];

  for (int i=0;i<iNa;i++) val1[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val2[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val3[i] = new FP1[gs];
  for (int i=0;i<iNa;i++) val4[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val5[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val6[i] = new FP1[gs];
  for (int i=0;i<iNa;i++) val7[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val8[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val9[i] = new FP1[gs];

  FP1** val1x = new FP1*[iNa];
  FP1** val2x = new FP1*[iN];
  FP1** val3x = new FP1*[iN];
  FP1** val4x = new FP1*[iNa];
  FP1** val5x = new FP1*[iN];
  FP1** val6x = new FP1*[iN];
  FP1** val7x = new FP1*[iNa];
  FP1** val8x = new FP1*[iN];
  FP1** val9x = new FP1*[iN];

  for (int i=0;i<iNa;i++) val1x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++)  val2x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++)  val3x[i] = new FP1[gs3];
  for (int i=0;i<iNa;i++) val4x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++)  val5x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++)  val6x[i] = new FP1[gs3];
  for (int i=0;i<iNa;i++) val7x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++)  val8x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++)  val9x[i] = new FP1[gs3];

  FP1* valt1 = new FP1[gs3];
  FP1* valt2 = new FP1[gs3];

  FP1* grid1s = new FP1[gs6];
  FP1* grid2s = new FP1[gs6];
  FP1* grid3s = new FP1[gs6];

  FP1* grid1p = new FP1[gs6];
  FP1* grid2p = new FP1[gs6];
  FP1* grid3p = new FP1[gs6];

  FP1* wtt1 = new FP1[gs];
  FP1* wtt2 = new FP1[gs];

  FP1* ang_g = new FP1[3*nang];
  FP1* ang_w = new FP1[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

  using namespace std::chrono;
  high_resolution_clock::time_point t0 = high_resolution_clock::now();

#if USE_ACC
 #pragma omp parallel for schedule(static,1) num_threads(nomp)
#endif
  for (int n=0;n<nomp;n++)
  {
    int tid = omp_get_thread_num();
    //printf("\n about to set device: %i \n",tid); fflush(stdout);
    #if USE_ACC
    acc_set_device_num(tid,acc_device_nvidia);
    #endif

    #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
    #pragma acc enter data copyin(n2i[0:natoms],na2i[0:natoms])

    #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
    #pragma acc enter data create(grid2[0:gs6],wt2[0:gs])
    #pragma acc enter data create(grid3[0:gs6],wt3[0:gs])
    #pragma acc enter data create(val1[0:iNa][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs]) 
    #pragma acc enter data create(val4[0:iNa][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs]) 
    #pragma acc enter data create(val7[0:iNa][0:gs],val8[0:iN][0:gs],val9[0:iN][0:gs]) 

    #pragma acc enter data create(val1x[0:iNa][0:gs3],val2x[0:iN][0:gs3],val3x[0:iN][0:gs3]) 
    #pragma acc enter data create(val4x[0:iNa][0:gs3],val5x[0:iN][0:gs3],val6x[0:iN][0:gs3]) 
    #pragma acc enter data create(val7x[0:iNa][0:gs3],val8x[0:iN][0:gs3],val9x[0:iN][0:gs3]) 

    #pragma acc enter data create(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6])
    #pragma acc enter data create(grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6])

    #pragma acc enter data create(wtt1[0:gs],wtt2[0:gs])
    #pragma acc enter data create(valt1[0:gs3],valt2[0:gs3])
    #pragma acc enter data copyin(norms1[0:Naux],norms2[0:N2])
    #pragma acc enter data create(xyz_grad[0:N3])

    if (tid>0)
    {
      #pragma acc enter data create(dC[0:N2a])
    }
    acc_assign(N3,xyz_grad,0.);
  }
  #if USE_ACC
  acc_set_device_num(0,acc_device_nvidia);
  #endif

 //fixes undercounting, corrupts dC
 #pragma acc parallel loop independent present(dC[0:N2a],n2i[0:natoms],na2i[0:natoms])
  for (int m=0;m<natoms;m++)
 #pragma acc loop independent
  for (int n=0;n<natoms;n++)
  if (m!=n)
  {
    int s1 = 0; if (m>0) s1 = na2i[m-1]; int s2 = na2i[m];
    int s3 = 0; if (m>0) s3 = n2i[m-1]; int s4 = n2i[m];
    int s5 = 0; if (n>0) s5 = n2i[n-1]; int s6 = n2i[n];
   #pragma acc loop collapse(3) 
    for (int i1=s1;i1<s2;i1++)
    for (int i2=s3;i2<s4;i2++)
    for (int i3=s5;i3<s6;i3++)
      dC[i1*N2+i2*N+i3] *= 2.;
  }

 //distribute dC to all gpus
  #if USE_ACC
  copy_to_all_gpu2(nomp,N2a,dC,0);
  #endif

#if USE_ACC
 #pragma omp parallel for num_threads(nomp)
#endif
  for (int m=0;m<natoms;m++)
  {
    int tid = omp_get_thread_num();
    //printf("\n about to set device(compute): %i \n",tid); fflush(stdout);
    #if USE_ACC
    acc_set_device_num(tid,acc_device_nvidia);
    #endif

    int s1 = 0; if (m>0) s1 = na2i[m-1]; int s2 = na2i[m];
    int s3 = 0; if (m>0) s3 = n2i[m-1]; int s4 = n2i[m];
    //printf("  m: %i  s1/2->3/4: %i %i - %i %i \n",m,s1,s2,s3,s4);

    FP1 Z1 = (FP1)atno[m];
    FP1 A1 = coords[3*m+0]; FP1 B1 = coords[3*m+1]; FP1 C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);

   #pragma acc parallel loop present(val1[0:iNa][0:gs],val1x[0:iNa][0:gs3])
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val1[ii1][j] = 1.f;
     #pragma acc loop
      for (int j=0;j<gs3;j++)
        val1x[ii1][j] = 1.f;
    }

   #pragma acc parallel loop present(val2[0:iN][0:gs],val2x[0:iN][0:gs3])
    for (int ii2=0;ii2<s4-s3;ii2++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val2[ii2][j] = 1.f;
     #pragma acc loop
      for (int j=0;j<gs3;j++)
        val2x[ii2][j] = 1.f;
    }
   #pragma acc parallel loop present(val3[0:iN][0:gs],val3x[0:iN][0:gs3],wt1[0:gs])
    for (int ii3=0;ii3<s4-s3;ii3++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val3[ii3][j] = wt1[j];
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val3x[ii3][3*j+0] = val3x[ii3][3*j+1] = val3x[ii3][3*j+2] = wt1[j];
    }

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<FP2> basis1 = basis_aux[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

      acc_assign(gs,valt1,1.);
      eval_sh_3r(gs,grid1,valt1,n1,l1,m1);

      eval_inr_d(gs,grid1,val1x[ii1],n1,l1,zeta1);
     #pragma acc parallel loop present(valt1[0:gs],val1x[0:iNa][0:gs3])
      for (int j=0;j<gs;j++)
      {
        FP1 v1 = valt1[j];
        val1x[ii1][3*j+0] *= v1;
        val1x[ii1][3*j+1] *= v1;
        val1x[ii1][3*j+2] *= v1;
      }

      eval_inr_r12(gs,grid1,val1[ii1],n1,l1,zeta1);

      if (l1>0)
      {
       #pragma acc parallel loop present(val1[0:iNa][0:gs],valt1[0:gs3])
        for (int j=0;j<gs;j++)
          valt1[3*j+0] = valt1[3*j+1] = valt1[3*j+2] = val1[ii1][j];

        eval_dp_3r(gs,grid1,valt1,n1,l1,m1);
        
       #pragma acc parallel loop present(valt1[0:gs3],val1x[0:iNa][0:gs3])
        for (int j=0;j<gs3;j++)
          val1x[ii1][j] += valt1[j];
      }

      eval_sh_3r(gs,grid1,val1[ii1],n1,l1,m1);
    } //loop i1

    for (int i2=s3;i2<s4;i2++)
    {
      int ii2 = i2-s3;

      vector<FP2> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

      eval_sh(ii2,gs,grid1,val2[ii2],n2,l2,m2,zeta2);
      eval_p(gs,grid1,val2x[ii2],n2,l2,m2,zeta2);
    } //loop i2

    for (int i3=s3;i3<s4;i3++)
    {
      int ii3 = i3-s3;

      vector<FP2> basis3 = basis[i3];
      int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];

      eval_sh(ii3,gs,grid1,val3[ii3],n3,l3,m3,zeta3);
      eval_p(gs,grid1,val3x[ii3],n3,l3,m3,zeta3);
    } //loop i3

    //reduce_3c1d(m,s1,s2,s3,s4,gs,norms1,norms2,dC,val1,val2,val3,val1x,val2x,val3x,N,Naux,iN,iNa,natoms,1.,xyz_grad);


   //two-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
     //s12 over atom m, aux function
     //s34 over atom m, regular function
     //s56 over atom n, regular function
      int s5 = 0; if (n>0) s5 = n2i[n-1]; int s6 = n2i[n];
      //printf("  mn: %i %i  s: %i-%i  %i-%i %i-%i \n",m,n,s1,s2,s3,s4,s5,s6);

      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2);
      recenter_grid(gs,grid2,A12,B12,C12);

      copy_grid(gs,grid1s,grid1); 
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid2,wt2,Z1,Z2,A12,B12,C12);
      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wt2);

      add_r1_to_grid(gs,grid2,0.,0.,0.);

     #pragma acc parallel loop present(val4[0:iNa][0:gs],val4x[0:iNa][0:gs])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val4[ii1][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          val4x[ii1][j] = 1.f;
      }

     #pragma acc parallel loop present(val2[0:iN][0:gs],val5[0:iN][0:gs],val2x[0:iN][0:gs3],val5x[0:iN][0:gs3])
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val2[ii2][j] = val5[ii2][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          val2x[ii2][j] = val5x[ii2][j] = 1.f;
      }

     #pragma acc parallel loop present(val3[0:iN][0:gs],val6[0:iN][0:gs],val3x[0:iN][0:gs3],val6x[0:iN][0:gs3],wtt1[0:gs],wt2[0:gs])
      for (int ii3=0;ii3<s6-s5;ii3++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val3[ii3][j] = wtt1[j];
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val6[ii3][j] = wt2[j];
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val3x[ii3][3*j+0] = val3x[ii3][3*j+1] = val3x[ii3][3*j+2] = wtt1[j];
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val6x[ii3][3*j+0] = val6x[ii3][3*j+1] = val6x[ii3][3*j+2] = wt2[j];
      }

     //i1 on atom m
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis_aux[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];
        //printf("  m: %i i1: %i   nlm: %i %i %2i zeta: %8.5f \n",m,i1,n1,l1,m1,zeta1);

        acc_assign(gs,valt1,1.);
        eval_sh_3r(gs,grid2,valt1,n1,l1,m1);

        eval_inr_d(gs,grid2,val4x[ii1],n1,l1,zeta1);
       #pragma acc parallel loop present(valt1[0:gs],val4x[0:iNa][0:gs3])
        for (int j=0;j<gs;j++)
        {
          FP1 v1 = valt1[j];
          val4x[ii1][3*j+0] *= v1;
          val4x[ii1][3*j+1] *= v1;
          val4x[ii1][3*j+2] *= v1;
        }

        eval_inr_r12(gs,grid2,val4[ii1],n1,l1,zeta1);

        if (l1>0)
        {
         #pragma acc parallel loop present(val4[0:iNa][0:gs],valt1[0:gs3])
          for (int j=0;j<gs;j++)
            valt1[3*j+0] = valt1[3*j+1] = valt1[3*j+2] = val4[ii1][j];

          eval_dp_3r(gs,grid2,valt1,n1,l1,m1);
        
         #pragma acc parallel loop present(valt1[0:gs3],val4x[0:iNa][0:gs3])
          for (int j=0;j<gs3;j++)
            val4x[ii1][j] += valt1[j];
        }

        eval_sh_3r(gs,grid2,val4[ii1],n1,l1,m1);
      }

     //i2 on atom m
      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;

        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];
        //printf("   n: %i i2: %i   nlm: %i %i %2i zeta: %8.5f \n",n,i2,n2,l2,m2,zeta2);

        eval_sh(ii2,gs,grid1,val2[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2,val5[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid1,val2x[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid2,val5x[ii2],n2,l2,m2,zeta2);
      }

     //i3 on atom n
      for (int i3=s5;i3<s6;i3++)
      {
        int ii3 = i3-s5;

        vector<FP2> basis3 = basis[i3];
        int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];

        eval_sh(ii3,gs,grid1s,val3[ii3],n3,l3,m3,zeta3);
        eval_sh(ii3,gs,grid2s,val6[ii3],n3,l3,m3,zeta3);
        eval_p(gs,grid1s,val3x[ii3],n3,l3,m3,zeta3);
        eval_p(gs,grid2s,val6x[ii3],n3,l3,m3,zeta3);
      }

      reduce_3c2d112(m,n,s1,s2,s3,s4,s5,s6,gs,norms1,norms2,dC,val1,val2,val3,val4,val5,val6,val1x,val2x,val3x,val4x,val5x,val6x,N,Naux,iN,iNa,natoms,1.,xyz_grad);


     #pragma acc parallel loop present(val2[0:iN][0:gs],val5[0:iN][0:gs],val2x[0:iN][0:gs3],val5x[0:iN][0:gs3])
      for (int ii2=0;ii2<s6-s5;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val2[ii2][j] = val5[ii2][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          val2x[ii2][j] = val5x[ii2][j] = 1.f;
      }

     #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],val3x[0:iN][0:gs3],val6x[0:iN][0:gs3],wt1[0:gs],wt2[0:gs])
      for (int ii3=0;ii3<s6-s5;ii3++)
      {
        for (int j=0;j<gs;j++)
        {
          val3[ii3][j] = wtt1[j];
          val6[ii3][j] = wt2[j];
          val3x[ii3][3*j+0] = val3x[ii3][3*j+1] = val3x[ii3][3*j+2] = wtt1[j];
          val6x[ii3][3*j+0] = val6x[ii3][3*j+1] = val6x[ii3][3*j+2] = wt2[j];
        }
      }

     //now do i2+i3 on atom n 
     //i2 on atom n
      for (int i2=s5;i2<s6;i2++)
      {
        int ii2 = i2-s5;

        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1s,val2[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2s,val5[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid1s,val2x[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid2s,val5x[ii2],n2,l2,m2,zeta2);
      }

     //i3 on atom n
      for (int i3=s5;i3<s6;i3++)
      {
        int ii3 = i3-s5;

        vector<FP2> basis3 = basis[i3];
        int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];

        eval_sh(ii3,gs,grid1s,val3[ii3],n3,l3,m3,zeta3);
        eval_sh(ii3,gs,grid2s,val6[ii3],n3,l3,m3,zeta3);
        eval_p(gs,grid1s,val3x[ii3],n3,l3,m3,zeta3);
        eval_p(gs,grid2s,val6x[ii3],n3,l3,m3,zeta3);
      }

      reduce_3c2d122(m,n,s1,s2,s5,s6,gs,norms1,norms2,dC,val1,val2,val3,val4,val5,val6,val1x,val2x,val3x,val4x,val5x,val6x,N,Naux,iN,iNa,natoms,1.,xyz_grad);
 
    } //loop n over second atom



   //three-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1);
      recenter_grid(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      for (int p=0;p<natoms;p++)
      if (p!=m && p!=n)
      {
        int s5 = 0; if (p>0) s5 = n2i[p-1]; int s6 = n2i[p];

        FP1 Z3 = (FP1)atno[p];
        FP1 A3 = coords[3*p+0]; FP1 B3 = coords[3*p+1]; FP1 C3 = coords[3*p+2];
        FP1 A13 = A3-A1; FP1 B13 = B3-B1; FP1 C13 = C3-C1;

        generate_central_grid_2(grid3,wt3,Z3,nrad,nang,ang_g,ang_w);
        copy_grid(gs,grid3p,grid3); //grid 3 centered on atom 3
        recenter_grid(gs,grid3,A13,B13,C13); //grid 3 centered on atom 1

        copy_grid(gs,grid3s,grid3);
        recenter_grid(gs,grid3s,-A12,-B12,-C12); //grid 3 centered on atom 2

        copy_grid(gs,grid2p,grid2); //grid 2 centered on atom 3
        recenter_grid(gs,grid2p,-A13,-B13,-C13);

        copy_grid(gs,grid1p,grid1);
        recenter_grid(gs,grid1p,-A13,-B13,-C13); //grid 1 centered on atom 3

        add_r1_to_grid_6z(gs,grid1s,grid2s,grid3s,grid1p,grid2p,grid3p);

      //need to keep all of these distances in order
        add_r3_to_grid(gs,grid1,A13,B13,C13);
        add_r123_to_grid(gs,grid2,A12,B12,C12,0.,0.,0.,A13,B13,C13);
        add_r123_to_grid(gs,grid3,A13,B13,C13,A12,B12,C12,0.,0.,0.);

        acc_copyf(gs,wtt1,wt1,wtt2,wt2);

        becke_weight_3c(gs,grid1,wtt1,grid2,wtt2,grid3,wt3,Z1,Z2,Z3,A12,B12,C12,A13,B13,C13);
        eliminate_small_wt_3(estart,gs,wtt1,wtt2,wt3);

        //#pragma acc update self(grid1[0:gs6],wtt1[0:gs],grid2[0:gs6],wtt2[0:gs],grid3[0:gs6],wt3[0:gs])
        //print_grid(gs,grid1,grid2,grid3,wtt1,wtt2,wt3,prl);

        add_r1_to_grid(gs,grid2,0.,0.,0.);
        add_r1_to_grid(gs,grid3,0.,0.,0.);


       #pragma acc parallel loop present(val4[0:iNa][0:gs],val7[0:iNa][0:gs],val4x[0:iNa][0:gs3],val7x[0:iNa][0:gs3])
        for (int ii1=0;ii1<s2-s1;ii1++)
        {
         #pragma acc loop
          for (int j=0;j<gs;j++)
            val4[ii1][j] = val7[ii1][j] = 1.f;
         #pragma acc loop
          for (int j=0;j<gs3;j++)
            val4x[ii1][j] = val7x[ii1][j] = 1.f;
        }

       #pragma acc parallel loop present(val2[0:iN][0:gs],val5[0:iN][0:gs],val8[0:iN][0:gs],val2x[0:iN][0:gs3],val5x[0:iN][0:gs3],val8x[0:iN][0:gs3])
        for (int ii2=0;ii2<s4-s3;ii2++)
        {
         #pragma acc loop
          for (int j=0;j<gs;j++)
            val2[ii2][j] = val5[ii2][j] = val8[ii2][j] = 1.f;
         #pragma acc loop
          for (int j=0;j<gs3;j++)
            val2x[ii2][j] = val5x[ii2][j] = val8x[ii2][j] = 1.f;
        }

       #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],val9[0:iN][0:gs],val3x[0:iN][0:gs3],val6x[0:iN][0:gs3],val9x[0:iN][0:gs3],wtt1[0:gs],wtt2[0:gs],wt3[0:gs])
        for (int ii3=0;ii3<s6-s5;ii3++)
        {
          for (int j=0;j<gs;j++)
          {
            val3[ii3][j] = wtt1[j];
            val6[ii3][j] = wtt2[j];
            val9[ii3][j] = wt3[j];
            val3x[ii3][3*j+0] = val3x[ii3][3*j+1] = val3x[ii3][3*j+2] = wtt1[j];
            val6x[ii3][3*j+0] = val6x[ii3][3*j+1] = val6x[ii3][3*j+2] = wtt2[j];
            val9x[ii3][3*j+0] = val9x[ii3][3*j+1] = val9x[ii3][3*j+2] = wt3[j];
          }
        }

 
        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<FP2> basis1 = basis_aux[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

          acc_assign(gs,valt1,1.);
          acc_assign(gs,valt2,1.);

          eval_sh_3r(gs,grid2,valt1,n1,l1,m1);
          eval_sh_3r(gs,grid3,valt2,n1,l1,m1);

          eval_inr_d(gs,grid2,val4x[ii1],n1,l1,zeta1);
          eval_inr_d(gs,grid3,val7x[ii1],n1,l1,zeta1);

         #pragma acc parallel loop present(valt1[0:gs],val4x[0:iNa][0:gs3],valt2[0:gs],val7x[0:iNa][0:gs3])
          for (int j=0;j<gs;j++)
          {
            FP1 v1 = valt1[j];
            val4x[ii1][3*j+0] *= v1; val4x[ii1][3*j+1] *= v1; val4x[ii1][3*j+2] *= v1;
            FP1 v2 = valt2[j];
            val7x[ii1][3*j+0] *= v2; val7x[ii1][3*j+1] *= v2; val7x[ii1][3*j+2] *= v2;
          }

          eval_inr_r12(gs,grid2,val4[ii1],n1,l1,zeta1);
          eval_inr_r12(gs,grid3,val7[ii1],n1,l1,zeta1);
          if (l1>0)
          {
           #pragma acc parallel loop present(val4[0:iNa][0:gs],val7[0:iNa][0:gs],valt1[0:gs3],valt2[0:gs3])
            for (int j=0;j<gs;j++)
            {
              valt1[3*j+0] = valt1[3*j+1] = valt1[3*j+2] = val4[ii1][j];
              valt2[3*j+0] = valt2[3*j+1] = valt2[3*j+2] = val7[ii1][j];
            }

            eval_dp_3r(gs,grid2,valt1,n1,l1,m1);
            eval_dp_3r(gs,grid3,valt2,n1,l1,m1);

           #pragma acc parallel loop present(valt1[0:gs3],valt2[0:gs3],val4x[0:iNa][0:gs3],val7x[0:iNa][0:gs3])
            for (int j=0;j<gs3;j++)
            {
              val4x[ii1][j] += valt1[j];
              val7x[ii1][j] += valt2[j];
            }
          }

          eval_sh_3r(gs,grid2,val4[ii1],n1,l1,m1);
          eval_sh_3r(gs,grid3,val7[ii1],n1,l1,m1);
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<FP2> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

          eval_sh(ii2,gs,grid3s,val8[ii2],n2,l2,m2,zeta2);
          eval_sh(ii2,gs,grid2s,val5[ii2],n2,l2,m2,zeta2);
          eval_sh(ii2,gs,grid1s,val2[ii2],n2,l2,m2,zeta2);
          eval_p(gs,grid3s,val8x[ii2],n2,l2,m2,zeta2);
          eval_p(gs,grid2s,val5x[ii2],n2,l2,m2,zeta2);
          eval_p(gs,grid1s,val2x[ii2],n2,l2,m2,zeta2);
        }

        for (int i3=s5;i3<s6;i3++)
        {
          int ii3 = i3-s5;
          vector<FP2> basis3 = basis[i3];
          int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];
    
          eval_sh(ii3,gs,grid3p,val9[ii3],n3,l3,m3,zeta3);
          eval_sh(ii3,gs,grid2p,val6[ii3],n3,l3,m3,zeta3);
          eval_sh(ii3,gs,grid1p,val3[ii3],n3,l3,m3,zeta3);
          eval_p(gs,grid3p,val9x[ii3],n3,l3,m3,zeta3);
          eval_p(gs,grid2p,val6x[ii3],n3,l3,m3,zeta3);
          eval_p(gs,grid1p,val3x[ii3],n3,l3,m3,zeta3);
        }

        reduce_3c3d(m,n,p,s1,s2,s3,s4,s5,s6,gs,norms1,norms2,dC,val1,val2,val3,val4,val5,val6,val7,val8,val9,val1x,val2x,val3x,val4x,val5x,val6x,val7x,val8x,val9x,N,Naux,iN,iNa,natoms,1.,xyz_grad);

      } //loop p over third atom
    } //loop n over second atom

  } //loop m over natoms

  FP2* xyz_all = new FP2[N3]();
  for (int n=0;n<nomp;n++)
  {
    //printf("\n about to set device(collect): %i \n",n); fflush(stdout);
    #if USE_ACC
    acc_set_device_num(n,acc_device_nvidia);
    #endif

    #pragma acc update self(xyz_grad[0:N3])
    for (int i=0;i<N3;i++)
      xyz_all[i] += xyz_grad[i];

    #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
    #pragma acc exit data delete(xyz_grad[0:N3])
  }
  #if USE_ACC
  acc_set_device_num(0,acc_device_nvidia);
  #endif

  for (int i=0;i<N3;i++)
    xyz_grad[i] = xyz_all[i];

  delete [] xyz_all;

#if DEBUG
  if (prl>1)
  {
    printf("\n xyz_grad: \n");
    for (int i=0;i<natoms;i++)
    {
      for (int j=0;j<3;j++)
        printf(" %8.5f",xyz_grad[3*i+j]);
      printf("\n");
    }
  }
#endif

#if USE_ACC
 #pragma omp parallel for schedule(static,1) num_threads(nomp)
#endif
  for (int n=0;n<nomp;n++)
  {
    int tid = omp_get_thread_num();
    //printf("\n about to set device(delete): %i \n",tid); fflush(stdout);
    #if USE_ACC
    acc_set_device_num(tid,acc_device_nvidia);
    #endif

    #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
    #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs])
    #pragma acc exit data delete(grid3[0:gs6],wt3[0:gs])

    #pragma acc exit data delete(val1[0:iNa][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs]) 
    #pragma acc exit data delete(val4[0:iNa][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs]) 
    #pragma acc exit data delete(val7[0:iNa][0:gs],val8[0:iN][0:gs],val9[0:iN][0:gs]) 

    #pragma acc exit data delete(val1x[0:iNa][0:gs3],val2x[0:iN][0:gs3],val3x[0:iN][0:gs3]) 
    #pragma acc exit data delete(val4x[0:iNa][0:gs3],val5x[0:iN][0:gs3],val6x[0:iN][0:gs3]) 
    #pragma acc exit data delete(val7x[0:iNa][0:gs3],val8x[0:iN][0:gs3],val9x[0:iN][0:gs3]) 

    #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6]) 
    #pragma acc exit data delete(grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6]) 

    #pragma acc exit data delete(wtt1[0:gs],wtt2[0:gs])
    #pragma acc exit data delete(valt1[0:gs],valt2[0:gs])
    #pragma acc exit data delete(norms1[0:Naux],norms2[0:N2])
    #pragma acc exit data delete(n2i[0:natoms],na2i[0:natoms])

    if (tid>0)
    {
      #pragma acc exit data delete(dC[0:N2a])
    }
  }
  #if USE_ACC
  acc_set_device_num(0,acc_device_nvidia);
  #endif

  delete [] norms1;
  delete [] norms2;

  delete [] n2i;
  delete [] na2i;

  delete [] ang_g;
  delete [] ang_w;

  delete [] grid1s;
  delete [] grid2s;
  delete [] grid3s;
  delete [] grid1p;
  delete [] grid2p;
  delete [] grid3p;

  delete [] wtt1;
  delete [] wtt2;
  delete [] grid1;
  delete [] grid2;
  delete [] grid3;
  delete [] wt1;
  delete [] wt2;
  delete [] wt3;
  for (int i=0;i<iNa;i++) delete [] val1[i];
  for (int i=0;i<iN;i++) delete [] val2[i];
  for (int i=0;i<iN;i++) delete [] val3[i];
  for (int i=0;i<iNa;i++) delete [] val4[i];
  for (int i=0;i<iN;i++) delete [] val5[i];
  for (int i=0;i<iN;i++) delete [] val6[i];
  for (int i=0;i<iNa;i++) delete [] val7[i];
  for (int i=0;i<iN;i++) delete [] val8[i];
  for (int i=0;i<iN;i++) delete [] val9[i];
  delete [] val1;
  delete [] val2;
  delete [] val3;
  delete [] val4;
  delete [] val5;
  delete [] val6;
  delete [] val7;
  delete [] val8;
  delete [] val9;

  for (int i=0;i<iNa;i++) delete [] val1x[i];
  for (int i=0;i<iN;i++) delete [] val2x[i];
  for (int i=0;i<iN;i++) delete [] val3x[i];
  for (int i=0;i<iNa;i++) delete [] val4x[i];
  for (int i=0;i<iN;i++) delete [] val5x[i];
  for (int i=0;i<iN;i++) delete [] val6x[i];
  for (int i=0;i<iNa;i++) delete [] val7x[i];
  for (int i=0;i<iN;i++) delete [] val8x[i];
  for (int i=0;i<iN;i++) delete [] val9x[i];
  delete [] val1x;
  delete [] val2x;
  delete [] val3x;
  delete [] val4x;
  delete [] val5x;
  delete [] val6x;
  delete [] val7x;
  delete [] val8x;
  delete [] val9x;
  delete [] valt1;
  delete [] valt2;

  return;
}


#if RED_DOUBLE
void compute_d_3c(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* dC, FP2* xyz_grad, int prl)
#else
void compute_d_3c(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, vector<vector<FP2> > &basis_aux, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* dC, FP1* xyz_grad, int prl)
#endif
{
  if (prl>1) printf(" compute_d_3c \n");

  int N = basis.size();
  int N2 = N*N;
  int N3 = 3*natoms;
  int Naux = basis_aux.size();
  int N2a = N2*Naux;
  int gs = nrad*nang;
  int gs3 = 3*gs;
  int gs6 = 6*gs;

  int estart = find_center_of_grid(1,nrad)*nang;

  int* na2i = new int[natoms];
  int iNa = get_imax_n2i(natoms,Naux,basis_aux,na2i);
  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);

  FP2* norms2 = new FP2[N2];
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    norms2[i*N+j] = basis[i][4]*basis[j][4];
  FP2* norms1 = new FP2[Naux];
  for (int i=0;i<Naux;i++)
    norms1[i] = norm_sv(basis_aux[i][0],basis_aux[i][1],basis_aux[i][2],basis_aux[i][3]);

  //printf("  iN/a: %i %i \n",iN,iNa);

  FP1* grid1 = new FP1[gs6];
  FP1* wt1 = new FP1[gs];

  FP1* grid2 = new FP1[gs6];
  FP1* wt2 = new FP1[gs];

  FP1* grid3 = new FP1[gs6];
  FP1* wt3 = new FP1[gs];

  FP1** val1 = new FP1*[iNa];
  FP1** val2 = new FP1*[iN];
  FP1** val3 = new FP1*[iN];
  FP1** val4 = new FP1*[iNa];
  FP1** val5 = new FP1*[iN];
  FP1** val6 = new FP1*[iN];
  FP1** val7 = new FP1*[iNa];
  FP1** val8 = new FP1*[iN];
  FP1** val9 = new FP1*[iN];

  for (int i=0;i<iNa;i++) val1[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val2[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val3[i] = new FP1[gs];
  for (int i=0;i<iNa;i++) val4[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val5[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val6[i] = new FP1[gs];
  for (int i=0;i<iNa;i++) val7[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val8[i] = new FP1[gs];
  for (int i=0;i<iN;i++)  val9[i] = new FP1[gs];

  FP1** val1x = new FP1*[iNa];
  FP1** val2x = new FP1*[iN];
  FP1** val3x = new FP1*[iN];
  FP1** val4x = new FP1*[iNa];
  FP1** val5x = new FP1*[iN];
  FP1** val6x = new FP1*[iN];
  FP1** val7x = new FP1*[iNa];
  FP1** val8x = new FP1*[iN];
  FP1** val9x = new FP1*[iN];

  for (int i=0;i<iNa;i++) val1x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++)  val2x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++)  val3x[i] = new FP1[gs3];
  for (int i=0;i<iNa;i++) val4x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++)  val5x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++)  val6x[i] = new FP1[gs3];
  for (int i=0;i<iNa;i++) val7x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++)  val8x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++)  val9x[i] = new FP1[gs3];

  FP1* valt1 = new FP1[gs3];
  FP1* valt2 = new FP1[gs3];

  FP1* grid1s = new FP1[gs6];
  FP1* grid2s = new FP1[gs6];
  FP1* grid3s = new FP1[gs6];

  FP1* grid1p = new FP1[gs6];
  FP1* grid2p = new FP1[gs6];
  FP1* grid3p = new FP1[gs6];

  FP1* wtt1 = new FP1[gs];
  FP1* wtt2 = new FP1[gs];

  FP1* ang_g = new FP1[3*nang];
  FP1* ang_w = new FP1[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];


 #if USE_ACC
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data copyin(n2i[0:natoms],na2i[0:natoms])

  #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
  #pragma acc enter data create(grid2[0:gs6],wt2[0:gs])
  #pragma acc enter data create(grid3[0:gs6],wt3[0:gs])
  #pragma acc enter data create(val1[0:iNa][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs]) 
  #pragma acc enter data create(val4[0:iNa][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs]) 
  #pragma acc enter data create(val7[0:iNa][0:gs],val8[0:iN][0:gs],val9[0:iN][0:gs]) 

  #pragma acc enter data create(val1x[0:iNa][0:gs3],val2x[0:iN][0:gs3],val3x[0:iN][0:gs3]) 
  #pragma acc enter data create(val4x[0:iNa][0:gs3],val5x[0:iN][0:gs3],val6x[0:iN][0:gs3]) 
  #pragma acc enter data create(val7x[0:iNa][0:gs3],val8x[0:iN][0:gs3],val9x[0:iN][0:gs3]) 

  #pragma acc enter data create(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6])
  #pragma acc enter data create(grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6])

  #pragma acc enter data create(wtt1[0:gs],wtt2[0:gs])
  #pragma acc enter data create(valt1[0:gs3],valt2[0:gs3])
  #pragma acc enter data copyin(norms1[0:Naux],norms2[0:N2])
 //dC[0:N2a],
  #pragma acc enter data create(xyz_grad[0:N3])
 #endif
  acc_assign(N3,xyz_grad,0.);


 //fixes undercounting, corrupts dC
 #pragma acc parallel loop independent present(dC[0:N2a],n2i[0:natoms],na2i[0:natoms])
  for (int m=0;m<natoms;m++)
 #pragma acc loop independent
  for (int n=0;n<natoms;n++)
  if (m!=n)
  {
    int s1 = 0; if (m>0) s1 = na2i[m-1]; int s2 = na2i[m];
    int s3 = 0; if (m>0) s3 = n2i[m-1]; int s4 = n2i[m];
    int s5 = 0; if (n>0) s5 = n2i[n-1]; int s6 = n2i[n];
   #pragma acc loop collapse(3) 
    for (int i1=s1;i1<s2;i1++)
    for (int i2=s3;i2<s4;i2++)
    for (int i3=s5;i3<s6;i3++)
      dC[i1*N2+i2*N+i3] *= 2.;
  }


  for (int m=0;m<natoms;m++)
  {
    int s1 = 0; if (m>0) s1 = na2i[m-1]; int s2 = na2i[m];
    int s3 = 0; if (m>0) s3 = n2i[m-1]; int s4 = n2i[m];
    //printf("  m: %i  s1/2->3/4: %i %i - %i %i \n",m,s1,s2,s3,s4);

    FP1 Z1 = (FP1)atno[m];
    FP1 A1 = coords[3*m+0]; FP1 B1 = coords[3*m+1]; FP1 C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);

   #pragma acc parallel loop present(val1[0:iNa][0:gs],val1x[0:iNa][0:gs3])
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val1[ii1][j] = 1.f;
     #pragma acc loop
      for (int j=0;j<gs3;j++)
        val1x[ii1][j] = 1.f;
    }

   #pragma acc parallel loop present(val2[0:iN][0:gs],val2x[0:iN][0:gs3])
    for (int ii2=0;ii2<s4-s3;ii2++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val2[ii2][j] = 1.f;
     #pragma acc loop
      for (int j=0;j<gs3;j++)
        val2x[ii2][j] = 1.f;
    }
   #pragma acc parallel loop present(val3[0:iN][0:gs],val3x[0:iN][0:gs3],wt1[0:gs])
    for (int ii3=0;ii3<s4-s3;ii3++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val3[ii3][j] = wt1[j];
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val3x[ii3][3*j+0] = val3x[ii3][3*j+1] = val3x[ii3][3*j+2] = wt1[j];
    }

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<FP2> basis1 = basis_aux[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

      acc_assign(gs,valt1,1.);
      eval_sh_3r(gs,grid1,valt1,n1,l1,m1);

      eval_inr_d(gs,grid1,val1x[ii1],n1,l1,zeta1);
     #pragma acc parallel loop present(valt1[0:gs],val1x[0:iNa][0:gs3])
      for (int j=0;j<gs;j++)
      {
        FP1 v1 = valt1[j];
        val1x[ii1][3*j+0] *= v1;
        val1x[ii1][3*j+1] *= v1;
        val1x[ii1][3*j+2] *= v1;
      }

      eval_inr_r12(gs,grid1,val1[ii1],n1,l1,zeta1);

      if (l1>0)
      {
       #pragma acc parallel loop present(val1[0:iNa][0:gs],valt1[0:gs3])
        for (int j=0;j<gs;j++)
          valt1[3*j+0] = valt1[3*j+1] = valt1[3*j+2] = val1[ii1][j];

        eval_dp_3r(gs,grid1,valt1,n1,l1,m1);
        
       #pragma acc parallel loop present(valt1[0:gs3],val1x[0:iNa][0:gs3])
        for (int j=0;j<gs3;j++)
          val1x[ii1][j] += valt1[j];
      }

      eval_sh_3r(gs,grid1,val1[ii1],n1,l1,m1);
    } //loop i1

    for (int i2=s3;i2<s4;i2++)
    {
      int ii2 = i2-s3;

      vector<FP2> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

      eval_sh(ii2,gs,grid1,val2[ii2],n2,l2,m2,zeta2);
      eval_p(gs,grid1,val2x[ii2],n2,l2,m2,zeta2);
    } //loop i2

    for (int i3=s3;i3<s4;i3++)
    {
      int ii3 = i3-s3;

      vector<FP2> basis3 = basis[i3];
      int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];

      eval_sh(ii3,gs,grid1,val3[ii3],n3,l3,m3,zeta3);
      eval_p(gs,grid1,val3x[ii3],n3,l3,m3,zeta3);
    } //loop i3

    //reduce_3c1d(m,s1,s2,s3,s4,gs,norms1,norms2,dC,val1,val2,val3,val1x,val2x,val3x,N,Naux,iN,iNa,natoms,1.,xyz_grad);


   //two-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
     //s12 over atom m, aux function
     //s34 over atom m, regular function
     //s56 over atom n, regular function
      int s5 = 0; if (n>0) s5 = n2i[n-1]; int s6 = n2i[n];
      //printf("  mn: %i %i  s: %i-%i  %i-%i %i-%i \n",m,n,s1,s2,s3,s4,s5,s6);

      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2);
      recenter_grid(gs,grid2,A12,B12,C12);

      copy_grid(gs,grid1s,grid1); 
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid2,wt2,Z1,Z2,A12,B12,C12);
      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wt2);

      add_r1_to_grid(gs,grid2,0.,0.,0.);

     #pragma acc parallel loop present(val4[0:iNa][0:gs],val4x[0:iNa][0:gs])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val4[ii1][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          val4x[ii1][j] = 1.f;
      }

     #pragma acc parallel loop present(val2[0:iN][0:gs],val5[0:iN][0:gs],val2x[0:iN][0:gs3],val5x[0:iN][0:gs3])
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val2[ii2][j] = val5[ii2][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          val2x[ii2][j] = val5x[ii2][j] = 1.f;
      }

     #pragma acc parallel loop present(val3[0:iN][0:gs],val6[0:iN][0:gs],val3x[0:iN][0:gs3],val6x[0:iN][0:gs3],wtt1[0:gs],wt2[0:gs])
      for (int ii3=0;ii3<s6-s5;ii3++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val3[ii3][j] = wtt1[j];
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val6[ii3][j] = wt2[j];
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val3x[ii3][3*j+0] = val3x[ii3][3*j+1] = val3x[ii3][3*j+2] = wtt1[j];
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val6x[ii3][3*j+0] = val6x[ii3][3*j+1] = val6x[ii3][3*j+2] = wt2[j];
      }

     //i1 on atom m
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis_aux[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];
        //printf("  m: %i i1: %i   nlm: %i %i %2i zeta: %8.5f \n",m,i1,n1,l1,m1,zeta1);

        acc_assign(gs,valt1,1.);
        eval_sh_3r(gs,grid2,valt1,n1,l1,m1);

        eval_inr_d(gs,grid2,val4x[ii1],n1,l1,zeta1);
       #pragma acc parallel loop present(valt1[0:gs],val4x[0:iNa][0:gs3])
        for (int j=0;j<gs;j++)
        {
          FP1 v1 = valt1[j];
          val4x[ii1][3*j+0] *= v1;
          val4x[ii1][3*j+1] *= v1;
          val4x[ii1][3*j+2] *= v1;
        }

        eval_inr_r12(gs,grid2,val4[ii1],n1,l1,zeta1);

        if (l1>0)
        {
         #pragma acc parallel loop present(val4[0:iNa][0:gs],valt1[0:gs3])
          for (int j=0;j<gs;j++)
            valt1[3*j+0] = valt1[3*j+1] = valt1[3*j+2] = val4[ii1][j];

          eval_dp_3r(gs,grid2,valt1,n1,l1,m1);
        
         #pragma acc parallel loop present(valt1[0:gs3],val4x[0:iNa][0:gs3])
          for (int j=0;j<gs3;j++)
            val4x[ii1][j] += valt1[j];
        }

        eval_sh_3r(gs,grid2,val4[ii1],n1,l1,m1);
      }

     //i2 on atom m
      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;

        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];
        //printf("   n: %i i2: %i   nlm: %i %i %2i zeta: %8.5f \n",n,i2,n2,l2,m2,zeta2);

        eval_sh(ii2,gs,grid1,val2[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2,val5[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid1,val2x[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid2,val5x[ii2],n2,l2,m2,zeta2);
      }

     //i3 on atom n
      for (int i3=s5;i3<s6;i3++)
      {
        int ii3 = i3-s5;

        vector<FP2> basis3 = basis[i3];
        int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];

        eval_sh(ii3,gs,grid1s,val3[ii3],n3,l3,m3,zeta3);
        eval_sh(ii3,gs,grid2s,val6[ii3],n3,l3,m3,zeta3);
        eval_p(gs,grid1s,val3x[ii3],n3,l3,m3,zeta3);
        eval_p(gs,grid2s,val6x[ii3],n3,l3,m3,zeta3);
      }

      reduce_3c2d112(m,n,s1,s2,s3,s4,s5,s6,gs,norms1,norms2,dC,val1,val2,val3,val4,val5,val6,val1x,val2x,val3x,val4x,val5x,val6x,N,Naux,iN,iNa,natoms,1.,xyz_grad);


     #pragma acc parallel loop present(val2[0:iN][0:gs],val5[0:iN][0:gs],val2x[0:iN][0:gs3],val5x[0:iN][0:gs3])
      for (int ii2=0;ii2<s6-s5;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val2[ii2][j] = val5[ii2][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          val2x[ii2][j] = val5x[ii2][j] = 1.f;
      }

     #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],val3x[0:iN][0:gs3],val6x[0:iN][0:gs3],wt1[0:gs],wt2[0:gs])
      for (int ii3=0;ii3<s6-s5;ii3++)
      {
        for (int j=0;j<gs;j++)
        {
          val3[ii3][j] = wtt1[j];
          val6[ii3][j] = wt2[j];
          val3x[ii3][3*j+0] = val3x[ii3][3*j+1] = val3x[ii3][3*j+2] = wtt1[j];
          val6x[ii3][3*j+0] = val6x[ii3][3*j+1] = val6x[ii3][3*j+2] = wt2[j];
        }
      }

     //now do i2+i3 on atom n 
     //i2 on atom n
      for (int i2=s5;i2<s6;i2++)
      {
        int ii2 = i2-s5;

        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1s,val2[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2s,val5[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid1s,val2x[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid2s,val5x[ii2],n2,l2,m2,zeta2);
      }

     //i3 on atom n
      for (int i3=s5;i3<s6;i3++)
      {
        int ii3 = i3-s5;

        vector<FP2> basis3 = basis[i3];
        int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];

        eval_sh(ii3,gs,grid1s,val3[ii3],n3,l3,m3,zeta3);
        eval_sh(ii3,gs,grid2s,val6[ii3],n3,l3,m3,zeta3);
        eval_p(gs,grid1s,val3x[ii3],n3,l3,m3,zeta3);
        eval_p(gs,grid2s,val6x[ii3],n3,l3,m3,zeta3);
      }

      reduce_3c2d122(m,n,s1,s2,s5,s6,gs,norms1,norms2,dC,val1,val2,val3,val4,val5,val6,val1x,val2x,val3x,val4x,val5x,val6x,N,Naux,iN,iNa,natoms,1.,xyz_grad);
 
    } //loop n over second atom



   //three-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1);
      recenter_grid(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      for (int p=0;p<natoms;p++)
      if (p!=m && p!=n)
      {
        int s5 = 0; if (p>0) s5 = n2i[p-1]; int s6 = n2i[p];

        FP1 Z3 = (FP1)atno[p];
        FP1 A3 = coords[3*p+0]; FP1 B3 = coords[3*p+1]; FP1 C3 = coords[3*p+2];
        FP1 A13 = A3-A1; FP1 B13 = B3-B1; FP1 C13 = C3-C1;

        generate_central_grid_2(grid3,wt3,Z3,nrad,nang,ang_g,ang_w);
        copy_grid(gs,grid3p,grid3); //grid 3 centered on atom 3
        recenter_grid(gs,grid3,A13,B13,C13); //grid 3 centered on atom 1

        copy_grid(gs,grid3s,grid3);
        recenter_grid(gs,grid3s,-A12,-B12,-C12); //grid 3 centered on atom 2

        copy_grid(gs,grid2p,grid2); //grid 2 centered on atom 3
        recenter_grid(gs,grid2p,-A13,-B13,-C13);

        copy_grid(gs,grid1p,grid1);
        recenter_grid(gs,grid1p,-A13,-B13,-C13); //grid 1 centered on atom 3

        add_r1_to_grid_6z(gs,grid1s,grid2s,grid3s,grid1p,grid2p,grid3p);

      //need to keep all of these distances in order
        add_r3_to_grid(gs,grid1,A13,B13,C13);
        add_r123_to_grid(gs,grid2,A12,B12,C12,0.,0.,0.,A13,B13,C13);
        add_r123_to_grid(gs,grid3,A13,B13,C13,A12,B12,C12,0.,0.,0.);

        acc_copyf(gs,wtt1,wt1,wtt2,wt2);

        becke_weight_3c(gs,grid1,wtt1,grid2,wtt2,grid3,wt3,Z1,Z2,Z3,A12,B12,C12,A13,B13,C13);
        eliminate_small_wt_3(estart,gs,wtt1,wtt2,wt3);

        //#pragma acc update self(grid1[0:gs6],wtt1[0:gs],grid2[0:gs6],wtt2[0:gs],grid3[0:gs6],wt3[0:gs])
        //print_grid(gs,grid1,grid2,grid3,wtt1,wtt2,wt3,prl);

        add_r1_to_grid(gs,grid2,0.,0.,0.);
        add_r1_to_grid(gs,grid3,0.,0.,0.);


       #pragma acc parallel loop present(val4[0:iNa][0:gs],val7[0:iNa][0:gs],val4x[0:iNa][0:gs3],val7x[0:iNa][0:gs3])
        for (int ii1=0;ii1<s2-s1;ii1++)
        {
         #pragma acc loop
          for (int j=0;j<gs;j++)
            val4[ii1][j] = val7[ii1][j] = 1.f;
         #pragma acc loop
          for (int j=0;j<gs3;j++)
            val4x[ii1][j] = val7x[ii1][j] = 1.f;
        }

       #pragma acc parallel loop present(val2[0:iN][0:gs],val5[0:iN][0:gs],val8[0:iN][0:gs],val2x[0:iN][0:gs3],val5x[0:iN][0:gs3],val8x[0:iN][0:gs3])
        for (int ii2=0;ii2<s4-s3;ii2++)
        {
         #pragma acc loop
          for (int j=0;j<gs;j++)
            val2[ii2][j] = val5[ii2][j] = val8[ii2][j] = 1.f;
         #pragma acc loop
          for (int j=0;j<gs3;j++)
            val2x[ii2][j] = val5x[ii2][j] = val8x[ii2][j] = 1.f;
        }

       #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],val9[0:iN][0:gs],val3x[0:iN][0:gs3],val6x[0:iN][0:gs3],val9x[0:iN][0:gs3],wtt1[0:gs],wtt2[0:gs],wt3[0:gs])
        for (int ii3=0;ii3<s6-s5;ii3++)
        {
          for (int j=0;j<gs;j++)
          {
            val3[ii3][j] = wtt1[j];
            val6[ii3][j] = wtt2[j];
            val9[ii3][j] = wt3[j];
            val3x[ii3][3*j+0] = val3x[ii3][3*j+1] = val3x[ii3][3*j+2] = wtt1[j];
            val6x[ii3][3*j+0] = val6x[ii3][3*j+1] = val6x[ii3][3*j+2] = wtt2[j];
            val9x[ii3][3*j+0] = val9x[ii3][3*j+1] = val9x[ii3][3*j+2] = wt3[j];
          }
        }

 
        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<FP2> basis1 = basis_aux[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

          acc_assign(gs,valt1,1.);
          acc_assign(gs,valt2,1.);

          eval_sh_3r(gs,grid2,valt1,n1,l1,m1);
          eval_sh_3r(gs,grid3,valt2,n1,l1,m1);

          eval_inr_d(gs,grid2,val4x[ii1],n1,l1,zeta1);
          eval_inr_d(gs,grid3,val7x[ii1],n1,l1,zeta1);

         #pragma acc parallel loop present(valt1[0:gs],val4x[0:iNa][0:gs3],valt2[0:gs],val7x[0:iNa][0:gs3])
          for (int j=0;j<gs;j++)
          {
            FP1 v1 = valt1[j];
            val4x[ii1][3*j+0] *= v1; val4x[ii1][3*j+1] *= v1; val4x[ii1][3*j+2] *= v1;
            FP1 v2 = valt2[j];
            val7x[ii1][3*j+0] *= v2; val7x[ii1][3*j+1] *= v2; val7x[ii1][3*j+2] *= v2;
          }

          eval_inr_r12(gs,grid2,val4[ii1],n1,l1,zeta1);
          eval_inr_r12(gs,grid3,val7[ii1],n1,l1,zeta1);
          if (l1>0)
          {
           #pragma acc parallel loop present(val4[0:iNa][0:gs],val7[0:iNa][0:gs],valt1[0:gs3],valt2[0:gs3])
            for (int j=0;j<gs;j++)
            {
              valt1[3*j+0] = valt1[3*j+1] = valt1[3*j+2] = val4[ii1][j];
              valt2[3*j+0] = valt2[3*j+1] = valt2[3*j+2] = val7[ii1][j];
            }

            eval_dp_3r(gs,grid2,valt1,n1,l1,m1);
            eval_dp_3r(gs,grid3,valt2,n1,l1,m1);

           #pragma acc parallel loop present(valt1[0:gs3],valt2[0:gs3],val4x[0:iNa][0:gs3],val7x[0:iNa][0:gs3])
            for (int j=0;j<gs3;j++)
            {
              val4x[ii1][j] += valt1[j];
              val7x[ii1][j] += valt2[j];
            }
          }

          eval_sh_3r(gs,grid2,val4[ii1],n1,l1,m1);
          eval_sh_3r(gs,grid3,val7[ii1],n1,l1,m1);
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<FP2> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

          eval_sh(ii2,gs,grid3s,val8[ii2],n2,l2,m2,zeta2);
          eval_sh(ii2,gs,grid2s,val5[ii2],n2,l2,m2,zeta2);
          eval_sh(ii2,gs,grid1s,val2[ii2],n2,l2,m2,zeta2);
          eval_p(gs,grid3s,val8x[ii2],n2,l2,m2,zeta2);
          eval_p(gs,grid2s,val5x[ii2],n2,l2,m2,zeta2);
          eval_p(gs,grid1s,val2x[ii2],n2,l2,m2,zeta2);
        }

        for (int i3=s5;i3<s6;i3++)
        {
          int ii3 = i3-s5;
          vector<FP2> basis3 = basis[i3];
          int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; FP1 zeta3 = basis3[3];
    
          eval_sh(ii3,gs,grid3p,val9[ii3],n3,l3,m3,zeta3);
          eval_sh(ii3,gs,grid2p,val6[ii3],n3,l3,m3,zeta3);
          eval_sh(ii3,gs,grid1p,val3[ii3],n3,l3,m3,zeta3);
          eval_p(gs,grid3p,val9x[ii3],n3,l3,m3,zeta3);
          eval_p(gs,grid2p,val6x[ii3],n3,l3,m3,zeta3);
          eval_p(gs,grid1p,val3x[ii3],n3,l3,m3,zeta3);
        }

        reduce_3c3d(m,n,p,s1,s2,s3,s4,s5,s6,gs,norms1,norms2,dC,val1,val2,val3,val4,val5,val6,val7,val8,val9,val1x,val2x,val3x,val4x,val5x,val6x,val7x,val8x,val9x,N,Naux,iN,iNa,natoms,1.,xyz_grad);

      } //loop p over third atom
    } //loop n over second atom

  } //loop m over natoms

 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc exit data copyout(xyz_grad[0:N3])
  //#pragma acc exit data delete(dC[0:N2a])
 #endif

#if DEBUG
  if (prl>1)
  {
    printf("\n xyz_grad: \n");
    for (int i=0;i<natoms;i++)
    {
      for (int j=0;j<3;j++)
        printf(" %8.5f",xyz_grad[3*i+j]);
      printf("\n");
    }
  }
#endif

  //printf("\n WARNING: transpose correct or not? \n");
  //transpose_C(Naux,N,C);

#if USE_ACC
  #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
  #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs])
  #pragma acc exit data delete(grid3[0:gs6],wt3[0:gs])

  #pragma acc exit data delete(val1[0:iNa][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs]) 
  #pragma acc exit data delete(val4[0:iNa][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs]) 
  #pragma acc exit data delete(val7[0:iNa][0:gs],val8[0:iN][0:gs],val9[0:iN][0:gs]) 

  #pragma acc exit data delete(val1x[0:iNa][0:gs3],val2x[0:iN][0:gs3],val3x[0:iN][0:gs3]) 
  #pragma acc exit data delete(val4x[0:iNa][0:gs3],val5x[0:iN][0:gs3],val6x[0:iN][0:gs3]) 
  #pragma acc exit data delete(val7x[0:iNa][0:gs3],val8x[0:iN][0:gs3],val9x[0:iN][0:gs3]) 

  #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6]) 
  #pragma acc exit data delete(grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6]) 

  #pragma acc exit data delete(wtt1[0:gs],wtt2[0:gs])
  #pragma acc exit data delete(valt1[0:gs],valt2[0:gs])
  #pragma acc exit data delete(norms1[0:Naux],norms2[0:N2])
  #pragma acc exit data delete(n2i[0:natoms],na2i[0:natoms])
#endif

  delete [] norms1;
  delete [] norms2;

  delete [] n2i;
  delete [] na2i;

  delete [] ang_g;
  delete [] ang_w;

  delete [] grid1s;
  delete [] grid2s;
  delete [] grid3s;
  delete [] grid1p;
  delete [] grid2p;
  delete [] grid3p;

  delete [] wtt1;
  delete [] wtt2;
  delete [] grid1;
  delete [] grid2;
  delete [] grid3;
  delete [] wt1;
  delete [] wt2;
  delete [] wt3;
  for (int i=0;i<iNa;i++) delete [] val1[i];
  for (int i=0;i<iN;i++) delete [] val2[i];
  for (int i=0;i<iN;i++) delete [] val3[i];
  for (int i=0;i<iNa;i++) delete [] val4[i];
  for (int i=0;i<iN;i++) delete [] val5[i];
  for (int i=0;i<iN;i++) delete [] val6[i];
  for (int i=0;i<iNa;i++) delete [] val7[i];
  for (int i=0;i<iN;i++) delete [] val8[i];
  for (int i=0;i<iN;i++) delete [] val9[i];
  delete [] val1;
  delete [] val2;
  delete [] val3;
  delete [] val4;
  delete [] val5;
  delete [] val6;
  delete [] val7;
  delete [] val8;
  delete [] val9;

  for (int i=0;i<iNa;i++) delete [] val1x[i];
  for (int i=0;i<iN;i++) delete [] val2x[i];
  for (int i=0;i<iN;i++) delete [] val3x[i];
  for (int i=0;i<iNa;i++) delete [] val4x[i];
  for (int i=0;i<iN;i++) delete [] val5x[i];
  for (int i=0;i<iN;i++) delete [] val6x[i];
  for (int i=0;i<iNa;i++) delete [] val7x[i];
  for (int i=0;i<iN;i++) delete [] val8x[i];
  for (int i=0;i<iN;i++) delete [] val9x[i];
  delete [] val1x;
  delete [] val2x;
  delete [] val3x;
  delete [] val4x;
  delete [] val5x;
  delete [] val6x;
  delete [] val7x;
  delete [] val8x;
  delete [] val9x;
  delete [] valt1;
  delete [] valt2;

  return;
}

#if RED_DOUBLE
void compute_d_En(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* Pao, FP2* xyz_grad, int prl)
#else
void compute_d_En(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* Pao, FP1* xyz_grad, int prl)
#endif
{
  if (prl>1) printf(" compute_d_En \n");
 //no <pVp> terms

  int N = basis.size();
  int N2 = N*N;
  int N3 = 3*natoms;

  int gs = nrad*nang;
  int gs3 = 3*gs; int gs6 = 6*gs;

  FP1* grid1 = new FP1[gs6];
  FP1* wt1 = new FP1[gs];

  FP1* grid2 = new FP1[gs6];
  FP1* wt2 = new FP1[gs];

  FP1* grid3 = new FP1[gs6];
  FP1* wt3 = new FP1[gs];

  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);
  //printf("  iN: %i \n",imaxN);

  int iN = imaxN;
  FP1* grid1s = new FP1[gs6]; FP1* grid2s = new FP1[gs6]; FP1* grid3s = new FP1[gs6];
  FP1* grid1p = new FP1[gs6]; FP1* grid2p = new FP1[gs6]; FP1* grid3p = new FP1[gs6];

  FP1** valS1 = new FP1*[iN]; for (int i=0;i<iN;i++) valS1[i] = new FP1[gs];
  FP1** valS2 = new FP1*[iN]; for (int i=0;i<iN;i++) valS2[i] = new FP1[gs];
  FP1** valS3 = new FP1*[iN]; for (int i=0;i<iN;i++) valS3[i] = new FP1[gs];
  FP1** valS4 = new FP1*[iN]; for (int i=0;i<iN;i++) valS4[i] = new FP1[gs];
  FP1** valS5 = new FP1*[iN]; for (int i=0;i<iN;i++) valS5[i] = new FP1[gs];
  FP1** valS6 = new FP1*[iN]; for (int i=0;i<iN;i++) valS6[i] = new FP1[gs];

  FP1** valS1x = new FP1*[iN]; for (int i=0;i<iN;i++) valS1x[i] = new FP1[gs3];
  FP1** valS2x = new FP1*[iN]; for (int i=0;i<iN;i++) valS2x[i] = new FP1[gs3];
  FP1** valS3x = new FP1*[iN]; for (int i=0;i<iN;i++) valS3x[i] = new FP1[gs3];
  FP1** valS4x = new FP1*[iN]; for (int i=0;i<iN;i++) valS4x[i] = new FP1[gs3];
  FP1** valS5x = new FP1*[iN]; for (int i=0;i<iN;i++) valS5x[i] = new FP1[gs3];
  FP1** valS6x = new FP1*[iN]; for (int i=0;i<iN;i++) valS6x[i] = new FP1[gs3];

  FP1** valt1 = new FP1*[iN]; for (int i=0;i<iN;i++) valt1[i] = new FP1[gs];
  FP1** valt2 = new FP1*[iN]; for (int i=0;i<iN;i++) valt2[i] = new FP1[gs];
  FP1** valt3 = new FP1*[iN]; for (int i=0;i<iN;i++) valt3[i] = new FP1[gs];
  FP1** valtx1 = new FP1*[iN]; for (int i=0;i<iN;i++) valtx1[i] = new FP1[gs3];
  FP1** valtx2 = new FP1*[iN]; for (int i=0;i<iN;i++) valtx2[i] = new FP1[gs3];
  FP1** valtx3 = new FP1*[iN]; for (int i=0;i<iN;i++) valtx3[i] = new FP1[gs3];

  FP1* wtt1 = new FP1[gs];
  FP1* wtt2 = new FP1[gs];

  FP1* ang_g = new FP1[3*nang];
  FP1* ang_w = new FP1[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

  FP2* norms = new FP2[N2];
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    norms[i*N+j] = basis[i][4]*basis[j][4];

 #if USE_ACC
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data copyin(n2i[0:natoms])
  #pragma acc enter data copyin(coords[0:3*natoms],atno[0:natoms])

  #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
  #pragma acc enter data create(grid2[0:gs6],wt2[0:gs]) 
  #pragma acc enter data create(grid3[0:gs6],wt3[0:gs]) 
  #pragma acc enter data create(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs],valS5[0:iN][0:gs],valS6[0:iN][0:gs])
  #pragma acc enter data create(valS1x[0:iN][0:gs3],valS2x[0:iN][0:gs3],valS3x[0:iN][0:gs3],valS4x[0:iN][0:gs3],valS5x[0:iN][0:gs3],valS6x[0:iN][0:gs3])
  #pragma acc enter data create(valt1[0:iN][0:gs],valt2[0:iN][0:gs],valt3[0:iN][0:gs])
  #pragma acc enter data create(valtx1[0:iN][0:gs3],valtx2[0:iN][0:gs3],valtx3[0:iN][0:gs3])
  #pragma acc enter data create(grid1s[0:gs6],grid1p[0:gs6],grid2s[0:gs6],grid2p[0:gs6],grid3s[0:gs6],grid3p[0:gs6])
  #pragma acc enter data create(wtt1[0:gs],wtt2[0:gs])
  #pragma acc enter data copyin(norms[0:N2])
 //,Pao[0:N2])
  #pragma acc enter data create(xyz_grad[0:N3])
 #endif
  acc_assign(N3,xyz_grad,0.);

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    FP1 Z1 = (FP1)atno[m];
    FP1 A1 = coords[3*m+0]; FP1 B1 = coords[3*m+1]; FP1 C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);

  #if USE_ACC
   #pragma acc parallel loop present(valS1[0:iN][0:gs],valS1x[0:iN][0:gs3],valS3[0:iN][0:gs],valS3x[0:iN][0:gs3])
  #endif
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valS1[ii1][j] = valS3[ii1][j] = 1.f;
     #pragma acc loop
      for (int j=0;j<gs3;j++)
        valS1x[ii1][j] = valS3x[ii1][j] = 1.f;
    }

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<FP2> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

      eval_sh(ii1,gs,grid1,valS1[ii1],n1,l1,m1,zeta1);
      eval_p(gs,grid1,valS1x[ii1],n1,l1,m1,zeta1);
    }

    for (int i2=s1;i2<s2;i2++)
    {
      int ii2 = i2-s1;
      vector<FP2> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

      eval_sh(ii2,gs,grid1,valS3[ii2],n2,l2,m2,zeta2);
      eval_p(gs,grid1,valS3x[ii2],n2,l2,m2,zeta2);
    } //loop i2 evaluate

   #if USE_ACC
    #pragma acc wait
   #endif


   //2 basis on same center, then 1 nucleus elsewhere
    for (int p=0;p<natoms;p++)
    if (p!=m)
    {
      FP1 Zn = (FP1)atno[p];
      FP1 An = coords[3*p+0]; FP1 Bn = coords[3*p+1]; FP1 Cn = coords[3*p+2];
      FP1 A1n = An-A1; FP1 B1n = Bn-B1; FP1 C1n = Cn-C1;

      add_r2_to_grid(gs,grid1,A1n,B1n,C1n);
      copy_grid(gs,grid1p,grid1);
      recenter_grid_zero(gs,grid1p,-A1n,-B1n,-C1n); //grid 1 centered on atom 3 (r1 set)

      generate_central_grid_2(grid3,wt3,Zn,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid3p,grid3); //grid 3 centered on atom 3
      recenter_grid(gs,grid3,A1n,B1n,C1n); //grid 3 centered on atom 1

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid3,wt3,Z1,Zn,A1n,B1n,C1n);

      eliminate_small_wt(gs,wtt1);
      eliminate_small_wt(gs,wt3);

      add_r1_to_grid(gs,grid3,0.,0.,0.);
      add_r2_to_grid(gs,grid3,A1n,B1n,C1n);

      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc parallel loop present(valt1[0:iN][0:gs],valS1[0:iN][0:gs],wtt1[0:gs])
        for (int j=0;j<gs;j++)
          valt1[ii1][j] = valS1[ii1][j]*wtt1[j];

       #pragma acc parallel loop collapse(2) present(valtx1[0:iN][0:gs3],valS1x[0:iN][0:gs3],wtt1[0:gs])
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valtx1[ii1][3*j+k] = valS1x[ii1][3*j+k]*wtt1[j];

       #pragma acc parallel loop present(valS2[0:iN][0:gs],wt3[0:gs])
        for (int j=0;j<gs;j++)
          valS2[ii1][j] = wt3[j];

       #pragma acc parallel loop collapse(2) present(valS2x[0:iN][0:gs3],wt3[0:gs])
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valS2x[ii1][3*j+k] = wt3[j];
      }


      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

        eval_sh(ii1,gs,grid3,valS2[ii1],n1,l1,m1,zeta1);
        eval_p(gs,grid3,valS2x[ii1],n1,l1,m1,zeta1);
      }

      for (int i2=s1;i2<s2;i2++)
      {
        int ii2 = i2-s1;
        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

        acc_assign(gs,valS4[ii2],1.);
        acc_assign(gs3,valS4x[ii2],1.);

        eval_sh(ii2,gs,grid3,valS4[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid3,valS4x[ii2],n2,l2,m2,zeta2);
      } //loop i2 evaluate

     #pragma acc parallel loop present(grid1[0:gs6],grid3[0:gs6],valt1[0:iN][0:gs],valS2[0:iN][0:gs],valtx1[0:iN][0:gs3],valS2x[0:iN][0:gs3])
      for (int j=0;j<gs;j++)
      {
        FP1 Rn1 = grid1[6*j+4]+1.e-20f;
        FP1 Rn3 = grid3[6*j+4]+1.e-20f;
        FP1 ne1 = Zn/Rn1;
        FP1 ne3 = Zn/Rn3;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valt1[ii1][j] *= ne1;
       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valS2[ii1][j] *= ne3;

       #pragma acc loop collapse(2)
        for (int ii1=0;ii1<s2-s1;ii1++)
        for (int k=0;k<3;k++)
          valtx1[ii1][3*j+k] *= ne1;
       #pragma acc loop collapse(2)
        for (int ii1=0;ii1<s2-s1;ii1++)
        for (int k=0;k<3;k++)
          valS2x[ii1][3*j+k] *= ne3;
      }

     //del basis on center m 
      reduce_2c2ds(m,p,s1,s2,s1,s2,gs,norms,Pao,valt1,valS2,valS3,valS4,valtx1,valS2x,valS3x,valS4x,iN,N,natoms,1.,xyz_grad);

    } //loop p over nuclear center



   //two-center basis
    //for (int n=m+1;n<natoms;n++)
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1); 
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      acc_copyf(gs,wtt1,wt1);
      acc_copyf(gs,wtt2,wt2);
      becke_weight_2c(gs,grid1,wtt1,grid2,wtt2,Z1,Z2,A12,B12,C12);

      eliminate_small_wt(gs,wtt1);
      eliminate_small_wt(gs,wtt2);

     //needs to happen after Becke weighting
      add_r1_to_grid(gs,grid2,0.,0.,0.);
      add_r2_to_grid(gs,grid2,A12,B12,C12);


    #if USE_ACC
     #pragma acc parallel loop present(valS2[0:iN][0:gs],valS2x[0:iN][0:gs3])
    #endif
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS2[ii1][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          valS2x[ii1][j] = 1.f;
      }

    #if USE_ACC
     #pragma acc parallel loop present(valS3[0:iN][0:gs],valS3x[0:iN][0:gs3],valS4[0:iN][0:gs],valS4x[0:iN][0:gs3])
    #endif
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS3[ii2][j] = valS4[ii2][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          valS3x[ii2][j] = valS4x[ii2][j] = 1.f;
      }

      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

        eval_sh(ii1,gs,grid2,valS2[ii1],n1,l1,m1,zeta1);
        eval_p(gs,grid2,valS2x[ii1],n1,l1,m1,zeta1);
      }

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1s,valS3[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2s,valS4[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid1s,valS3x[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid2s,valS4x[ii2],n2,l2,m2,zeta2);
      }


     //2 basis on diff centers, 1 nucleus matches
      gather_12_d_En_0(s1,s2,gs,iN,valtx1,valtx2,valS1x,valS2x,wtt1,wtt2);

     #pragma acc parallel loop present(grid1[0:gs6],grid2[0:gs6],valtx1[0:iN][0:gs3],valtx2[0:iN][0:gs3])
      for (int j=0;j<gs;j++)
      {
        //FP1 Rn2a = grid1s[6*j+3]; FP1 Rn2b = grid2s[6*j+3];
        FP1 Rn2a = grid1[6*j+4]+1.e-20f;
        FP1 Rn2b = grid2[6*j+4]+1.e-20f;
        FP1 ne1 = Z2/Rn2a;
        FP1 ne2 = Z2/Rn2b;

       #pragma acc loop collapse(2)
        for (int ii1=0;ii1<s2-s1;ii1++)
        for (int k=0;k<3;k++)
          valtx1[ii1][3*j+k] *= ne1;

       #pragma acc loop collapse(2)
        for (int ii1=0;ii1<s2-s1;ii1++)
        for (int k=0;k<3;k++)
          valtx2[ii1][3*j+k] *= ne2;
      }

     //del basis only, first center
      reduce_2c1ds(m,n,s3,s4,s1,s2,gs,norms,Pao,valS3,valS4,valtx1,valtx2,iN,N,natoms,0.5,xyz_grad);

      gather_12_d_En_0(s3,s4,gs,iN,valtx1,valtx2,valS3x,valS4x,wtt1,wtt2);

     #pragma acc parallel loop present(grid1[0:gs6],grid2[0:gs6],valt1[0:iN][0:gs],valt2[0:iN][0:gs],valtx1[0:iN][0:gs3],valtx2[0:iN][0:gs3])
      for (int j=0;j<gs;j++)
      {
        FP1 Rn1a = grid1[6*j+3]+1.e-20f; FP1 Rn1b = grid2[6*j+3]+1.e-20f;
        FP1 ne1 = Z1/Rn1a;
        FP1 ne2 = Z1/Rn1b;

       #pragma acc loop collapse(2)
        for (int ii2=0;ii2<s4-s3;ii2++)
        for (int k=0;k<3;k++)
          valtx1[ii2][3*j+k] *= ne1;

       #pragma acc loop collapse(2)
        for (int ii2=0;ii2<s4-s3;ii2++)
        for (int k=0;k<3;k++)
          valtx2[ii2][3*j+k] *= ne2;
      }

     //del basis only, second center
      reduce_2c1ds(n,m,s1,s2,s3,s4,gs,norms,Pao,valS1,valS2,valtx1,valtx2,iN,N,natoms,0.5,xyz_grad);

     //2 basis on diff centers, then 1 nucleus elsewhere
      for (int p=0;p<natoms;p++)
      if (p!=m && p!=n)
      {
        FP1 Zn = (FP1)atno[p];
        FP1 An = coords[3*p+0]; FP1 Bn = coords[3*p+1]; FP1 Cn = coords[3*p+2];
        FP1 A1n = An-A1; FP1 B1n = Bn-B1; FP1 C1n = Cn-C1;

        generate_central_grid_2(grid3,wt3,Zn,nrad,nang,ang_g,ang_w);
        copy_grid(gs,grid3p,grid3); //grid 3 centered on atom 3
        recenter_grid(gs,grid3,A1n,B1n,C1n); //grid 3 centered on atom 1

        copy_grid(gs,grid3s,grid3);
        recenter_grid_zero(gs,grid3s,-A12,-B12,-C12); //grid 3 centered on atom 2

        copy_grid(gs,grid1p,grid1); //grid 1 centered on atom 3
        recenter_grid_zero(gs,grid1p,-A1n,-B1n,-C1n);

        copy_grid(gs,grid2p,grid2); //grid 2 centered on atom 3
        recenter_grid_zero(gs,grid2p,-A1n,-B1n,-C1n);

       //need to keep all of these distances in order
        add_r123_to_grid(gs,grid1,0.,0.,0.,A12,B12,C12,A1n,B1n,C1n);
        add_r123_to_grid(gs,grid2,A12,B12,C12,0.,0.,0.,A1n,B1n,C1n);
      	add_r123_to_grid(gs,grid3,A1n,B1n,C1n,A12,B12,C12,0.,0.,0.);
     
        acc_copyf(gs,wtt1,wt1,wtt2,wt2);

        becke_weight_3c(gs,grid1,wtt1,grid2,wtt2,grid3,wt3,Z1,Z2,Zn,A12,B12,C12,A1n,B1n,C1n);
        eliminate_small_wt_3(gs,wtt1,wtt2,wt3);

        add_r2_to_grid(gs,grid1,A1n,B1n,C1n);
        add_r1_to_grid(gs,grid2,0.,0.,0.);
        add_r1_to_grid(gs,grid3,0.,0.,0.);
        add_r2_to_grid(gs,grid2,A1n,B1n,C1n);
        add_r2_to_grid(gs,grid3,A1n,B1n,C1n);


      #if USE_ACC
       #pragma acc parallel loop collapse(2) present(valS6[0:iN][0:gs],valS6x[0:iN][0:gs3])
      #endif
        for (int ii2=0;ii2<s4-s3;ii2++)
        {
          for (int j=0;j<gs;j++)
            valS6[ii2][j] = valS6x[ii2][3*j+0] = valS6x[ii2][3*j+1] = valS6x[ii2][3*j+2] = 1.f;
        }

        gather_125_d_En(s1,s2,gs,iN,valt1,valt2,valS1,valS2,valS5,valtx1,valtx2,valS1x,valS2x,valS5x,wtt1,wtt2,wt3);

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<FP2> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

          eval_sh(ii1,gs,grid3,valS5[ii1],n1,l1,m1,zeta1);
          eval_p(gs,grid3,valS5x[ii1],n1,l1,m1,zeta1);
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<FP2> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

          eval_sh(ii2,gs,grid3s,valS6[ii2],n2,l2,m2,zeta2);
          eval_p(gs,grid3s,valS6x[ii2],n2,l2,m2,zeta2);
        } //loop i2 evaluate

       #pragma acc parallel loop present(grid1[0:gs6],grid2[0:gs6],grid3[0:gs6],valt1[0:iN][0:gs],valt2[0:iN][0:gs],valS5[0:iN][0:gs],valtx1[0:iN][0:gs3],valtx2[0:iN][0:gs3],valS5x[0:iN][0:gs3])
        for (int j=0;j<gs;j++)
        {
          //FP1 Rn1 = grid1p[6*j+3]+1.e-20f; FP1 Rn2 = grid2p[6*j+3]+1.e-20f; FP1 Rn3 = grid3p[6*j+3]+1.e-20f;
          FP1 Rn1 = grid1[6*j+4]+1.e-20f; FP1 Rn2 = grid2[6*j+4]+1.e-20f; FP1 Rn3 = grid3[6*j+4]+1.e-20f;
          FP1 ne1 = Zn/Rn1; FP1 ne2 = Zn/Rn2; FP1 ne3 = Zn/Rn3;

         #pragma acc loop
          for (int ii1=0;ii1<s2-s1;ii1++)
          {
            valt1[ii1][j] *= ne1;
            valt2[ii1][j] *= ne2;
            valS5[ii1][j] *= ne3;
          }

         #pragma acc loop collapse(2)
          for (int ii1=0;ii1<s2-s1;ii1++)
          for (int k=0;k<3;k++)
          {
            valtx1[ii1][3*j+k] *= ne1;
            valtx2[ii1][3*j+k] *= ne2;
            valS5x[ii1][3*j+k] *= ne3;
          }
        }

        //collect 3-atom values
        int m1 = m; int n1 = n;
        reduce_2c3d(m1,n1,p,s1,s2,s3,s4,gs,norms,Pao,valt1,valt2,valS3,valS4,valS5,valS6,valtx1,valtx2,valS3x,valS4x,valS5x,valS6x,iN,N,natoms,0.5,xyz_grad);

      } //loop p over nuclear center

    } //loop n over second atom


  } //loop m over natoms

 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  //#pragma acc exit data delete(Pao[0:N2])
  #pragma acc exit data copyout(xyz_grad[0:N3])
 #endif

#if DEBUG
  if (prl>1)
  {
    printf("\n xyz_grad (En): \n");
    for (int i=0;i<natoms;i++)
    {
      for (int j=0;j<3;j++)
        printf(" %8.5f",xyz_grad[i*3+j]);
      printf("\n");
    }
  }
#endif

#if USE_ACC
  #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
  #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs]) 
  #pragma acc exit data delete(grid3[0:gs6],wt3[0:gs]) 
  #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6],grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6])
  #pragma acc exit data delete(wtt1[0:gs],wtt2[0:gs])
  #pragma acc exit data delete(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs],valS5[0:iN][0:gs],valS6[0:iN][0:gs])
  #pragma acc exit data delete(valS1x[0:iN][0:gs3],valS2x[0:iN][0:gs3],valS3x[0:iN][0:gs3],valS4x[0:iN][0:gs3],valS5x[0:iN][0:gs3],valS6x[0:iN][0:gs3])
  #pragma acc exit data delete(valt1[0:iN][0:gs],valt2[0:iN][0:gs],valt3[0:iN][0:gs],valtx1[0:iN][0:gs3],valtx2[0:iN][0:gs3],valtx3[0:iN][0:gs3])
  #pragma acc exit data delete(n2i[0:natoms])
  #pragma acc exit data delete(coords[0:3*natoms],atno[0:natoms])
  #pragma acc exit data delete(norms[0:N2])
#endif

  delete [] ang_g;
  delete [] ang_w;

  delete [] n2i;
  delete [] norms;

  delete [] grid1s;
  delete [] grid2s;
  delete [] grid3s;
  delete [] grid1p;
  delete [] grid2p;
  delete [] grid3p;

  for (int i=0;i<iN;i++) delete [] valS1[i];
  for (int i=0;i<iN;i++) delete [] valS2[i];
  for (int i=0;i<iN;i++) delete [] valS3[i];
  for (int i=0;i<iN;i++) delete [] valS4[i];
  for (int i=0;i<iN;i++) delete [] valS5[i];
  for (int i=0;i<iN;i++) delete [] valS6[i];
  for (int i=0;i<iN;i++) delete [] valS1x[i];
  for (int i=0;i<iN;i++) delete [] valS2x[i];
  for (int i=0;i<iN;i++) delete [] valS3x[i];
  for (int i=0;i<iN;i++) delete [] valS4x[i];
  for (int i=0;i<iN;i++) delete [] valS5x[i];
  for (int i=0;i<iN;i++) delete [] valS6x[i];

  for (int i=0;i<iN;i++) delete [] valt1[i];
  for (int i=0;i<iN;i++) delete [] valt2[i];
  for (int i=0;i<iN;i++) delete [] valt3[i];
  for (int i=0;i<iN;i++) delete [] valtx1[i];
  for (int i=0;i<iN;i++) delete [] valtx2[i];
  for (int i=0;i<iN;i++) delete [] valtx3[i];

  delete [] valS1; delete [] valS2; delete [] valS3; delete [] valS4; delete [] valS5; delete [] valS6;
  delete [] valS1x; delete [] valS2x; delete [] valS3x; delete [] valS4x; delete [] valS5x; delete [] valS6x;
  delete [] valt1; delete [] valt2; delete [] valt3; delete [] valtx1; delete [] valtx2; delete [] valtx3;
  delete [] wtt1;
  delete [] wtt2;

  delete [] grid1;
  delete [] grid2;
  delete [] grid3;
  delete [] wt1;
  delete [] wt2;
  delete [] wt3;

  return;
}

#if RED_DOUBLE
void compute_d_ST(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* GFao, FP2* Pao, FP2* xyz_grad, int prl)
#else
void compute_d_ST(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* GFao, FP2* Pao, FP1* xyz_grad, int prl)
#endif
{
  if (prl>1) printf(" compute_d_ST \n");

  int N = basis.size();
  int N2 = N*N;
  int N3 = 3*natoms;

  int gs = nrad*nang; 
  int gs3 = 3*gs;
  int gs6 = 6*gs;
  FP1* grid1 = new FP1[gs6];
  FP1* wt1 = new FP1[gs];

  FP1* grid2 = new FP1[gs6];
  FP1* wt2 = new FP1[gs];

  FP1* grid1s = new FP1[gs6];
  FP1* grid2s = new FP1[gs6];
  FP1* wtt1 = new FP1[gs];

  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);
  //printf("  iN: %i \n",imaxN);

  FP2* norms = new FP2[N2];
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    norms[i*N+j] = basis[i][4]*basis[j][4];

 //intermediate storage
  int iN = imaxN;
  FP1** valS1 = new FP1*[iN]; for (int i=0;i<iN;i++) valS1[i] = new FP1[gs];
  FP1** valS2 = new FP1*[iN]; for (int i=0;i<iN;i++) valS2[i] = new FP1[gs];
  FP1** valS3 = new FP1*[iN]; for (int i=0;i<iN;i++) valS3[i] = new FP1[gs];
  FP1** valS4 = new FP1*[iN]; for (int i=0;i<iN;i++) valS4[i] = new FP1[gs];

  FP1** valS1x = new FP1*[iN]; for (int i=0;i<iN;i++) valS1x[i] = new FP1[gs3];
  FP1** valS2x = new FP1*[iN]; for (int i=0;i<iN;i++) valS2x[i] = new FP1[gs3];
  FP1** valS3x = new FP1*[iN]; for (int i=0;i<iN;i++) valS3x[i] = new FP1[gs3];
  FP1** valS4x = new FP1*[iN]; for (int i=0;i<iN;i++) valS4x[i] = new FP1[gs3];

  FP1** valT1 = new FP1*[iN]; for (int i=0;i<iN;i++) valT1[i] = new FP1[gs];
  FP1** valT2 = new FP1*[iN]; for (int i=0;i<iN;i++) valT2[i] = new FP1[gs];
  FP1** valT1x = new FP1*[iN]; for (int i=0;i<iN;i++) valT1x[i] = new FP1[gs3];
  FP1** valT2x = new FP1*[iN]; for (int i=0;i<iN;i++) valT2x[i] = new FP1[gs3];

  FP1* valt = new FP1[gs3];

  FP1* ang_g = new FP1[3*nang];
  FP1* ang_w = new FP1[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

 #if USE_ACC
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data copyin(n2i[0:natoms])
  #pragma acc enter data copyin(coords[0:3*natoms],atno[0:natoms])
  //#pragma acc enter data copyin(Pao[0:N2],GFao[0:N2])
  #pragma acc enter data copyin(norms[0:N2])

  #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
  #pragma acc enter data create(grid2[0:gs6],wt2[0:gs]) 
  #pragma acc enter data create(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs])
  #pragma acc enter data create(valT1[0:iN][0:gs],valT2[0:iN][0:gs])
  #pragma acc enter data create(valS1x[0:iN][0:gs3],valS2x[0:iN][0:gs3],valS3x[0:iN][0:gs3],valS4x[0:iN][0:gs3])
  #pragma acc enter data create(valT1x[0:iN][0:gs3],valT2x[0:iN][0:gs3])
  #pragma acc enter data create(valt[0:gs3])
  #pragma acc enter data create(grid1s[0:gs6],grid2s[0:gs6],wtt1[0:gs])
  #pragma acc enter data create(xyz_grad[0:N3])
 #endif
  acc_assign(N3,xyz_grad,0.);

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    FP1 Z1 = (FP1)atno[m];
    FP1 A1 = coords[3*m+0]; FP1 B1 = coords[3*m+1]; FP1 C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);

  #if USE_ACC
   #pragma acc parallel loop present(valS1[0:iN][0:gs],valS1x[0:iN][0:gs3])
  #endif
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valS1[ii1][j] = 1.f;
     #pragma acc loop
      for (int j=0;j<gs3;j++)
        valS1x[ii1][j] = 1.f;
    }

  #if USE_ACC
   #pragma acc parallel loop present(valS3[0:iN][0:gs],valS3x[0:iN][0:gs3],wt1[0:gs])
  #endif
    for (int ii2=0;ii2<s2-s1;ii2++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valS3[ii2][j] = wt1[j];
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valS3x[ii2][3*j+0] = valS3x[ii2][3*j+1] = valS3x[ii2][3*j+2] = wt1[j];
    }

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<FP2> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

      eval_sh(ii1,gs,grid1,valS1[ii1],n1,l1,m1,zeta1);
      eval_p(gs,grid1,valS1x[ii1],n1,l1,m1,zeta1);
    }

    for (int i2=s1;i2<s2;i2++)
    {
      int ii2 = i2-s1;
      vector<FP2> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

      eval_sh(ii2,gs,grid1,valS3[ii2],n2,l2,m2,zeta2);
      eval_p(gs,grid1,valS3x[ii2],n2,l2,m2,zeta2);
    } //loop i2 evaluate

  #if USE_ACC
   #pragma acc parallel loop present(valS1[0:iN][0:gs],valT1[0:iN][0:gs],valT1x[0:iN][0:gs3])
  #endif
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valT1[ii1][j] = valS1[ii1][j];
     #pragma acc loop
      for (int j=0;j<gs3;j++)
        valT1x[ii1][j] = 1.f;
    }
 
   //KE terms
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<FP2> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

      eval_p(gs,grid1,valT1x[ii1],n1,l1,m1,zeta1);
      eval_ke3(gs,grid1,valT1x[ii1],n1,l1,zeta1);

     #pragma acc parallel loop present(valS1[0:iN][0:gs],valt[0:gs3])
      for (int j=0;j<gs;j++)
        valt[3*j+0] = valt[3*j+1] = valt[3*j+2] = valS1[ii1][j];
      eval_dke(gs,grid1,valt,n1,l1,zeta1);

     #pragma acc parallel loop present(valt[0:gs3],valT1x[0:iN][0:gs3])
      for (int j=0;j<gs3;j++)
        valT1x[ii1][j] += valt[j];

     //regular term
      eval_ke(gs,grid1,valT1[ii1],n1,l1,zeta1);
    }

   #if USE_ACC
    #pragma acc wait
   #endif

    //reduce_2c1d(m,s1,s2,s1,s2,gs,norms,GFao,valS1,valS3,valS1x,valS3x,iN,N,natoms,1.,xyz_grad);
    //reduce_2c1d(m,s1,s2,s1,s2,gs,norms,Pao, valT1,valS3,valT1x,valS3x,iN,N,natoms,-0.5,xyz_grad);


   //two-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];
      //printf(" mn: %i %i s1-4: %i %i - %i %i \n",m,n,s1,s2,s3,s4);

      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1); 
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid2,wt2,Z1,Z2,A12,B12,C12);

      //#pragma acc update self(grid1[0:6*gs],wtt1[0:gs],grid2[0:6*gs],wt2[0:gs])
      //print_grid(gs,grid1,grid2,wtt1,wt2,prl);

      eliminate_small_wt(gs,wtt1);
      eliminate_small_wt(gs,wt2);

     //needs to happen after Becke weighting
      add_r1_to_grid(gs,grid2,0.,0.,0.);

    #if USE_ACC
     #pragma acc parallel loop present(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS1x[0:iN][0:gs3],valS2x[0:iN][0:gs3])
    #endif
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS1[ii1][j] = valS2[ii1][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          valS1x[ii1][j] = valS2x[ii1][j] = 1.f;
      }

    #if USE_ACC
     #pragma acc parallel loop present(valS3[0:iN][0:gs],valS4[0:iN][0:gs],valS3x[0:iN][0:gs3],valS4x[0:iN][0:gs3],wtt1[0:gs],wt2[0:gs])
    #endif
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
        {
          valS3[ii2][j] = wtt1[j];
          valS4[ii2][j] = wt2[j];
          valS3x[ii2][3*j+0] = valS3x[ii2][3*j+1] = valS3x[ii2][3*j+2] = wtt1[j];
          valS4x[ii2][3*j+0] = valS4x[ii2][3*j+1] = valS4x[ii2][3*j+2] = wt2[j];
        }
      }

      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

        eval_sh(ii1,gs,grid1,valS1[ii1],n1,l1,m1,zeta1);
        eval_sh(ii1,gs,grid2,valS2[ii1],n1,l1,m1,zeta1);
        eval_p(gs,grid1,valS1x[ii1],n1,l1,m1,zeta1); //center 1 moves
        eval_p(gs,grid2,valS2x[ii1],n1,l1,m1,zeta1); //center 1 moves
      }

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1s,valS3[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2s,valS4[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid1s,valS3x[ii2],n2,l2,m2,zeta2); //center 2 moves
        eval_p(gs,grid2s,valS4x[ii2],n2,l2,m2,zeta2); //center 2 moves
      }

    #if USE_ACC
     #pragma acc parallel loop present(valS1[0:iN][0:gs],valT1[0:iN][0:gs],valT1x[0:iN][0:gs3])
    #endif
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valT1[ii1][j] = valS1[ii1][j];
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          valT1x[ii1][j] = 1.f;
      }
    #if USE_ACC
     #pragma acc parallel loop present(valS2[0:iN][0:gs],valT2[0:iN][0:gs],valT2x[0:iN][0:gs3])
    #endif
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valT2[ii1][j] = valS2[ii1][j];
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          valT2x[ii1][j] = 1.f;
      }

     //KE terms
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

        eval_ke(gs,grid1,valT1[ii1],n1,l1,zeta1);
        eval_ke(gs,grid2,valT2[ii1],n1,l1,zeta1);

        eval_p(gs,grid1,valT1x[ii1],n1,l1,m1,zeta1);
        eval_ke3(gs,grid1,valT1x[ii1],n1,l1,zeta1);
        eval_p(gs,grid2,valT2x[ii1],n1,l1,m1,zeta1);
        eval_ke3(gs,grid2,valT2x[ii1],n1,l1,zeta1);

       #pragma acc parallel loop present(valS1[0:iN][0:gs],valt[0:gs3])
        for (int j=0;j<gs;j++)
          valt[3*j+0] = valt[3*j+1] = valt[3*j+2] = valS1[ii1][j];
        eval_dke(gs,grid1,valt,n1,l1,zeta1);

       #pragma acc parallel loop present(valt[0:gs3],valT1x[0:iN][0:gs3])
        for (int j=0;j<gs3;j++)
          valT1x[ii1][j] += valt[j];

       #pragma acc parallel loop present(valS2[0:iN][0:gs],valt[0:gs3])
        for (int j=0;j<gs;j++)
          valt[3*j+0] = valt[3*j+1] = valt[3*j+2] = valS2[ii1][j];
        eval_dke(gs,grid2,valt,n1,l1,zeta1);

       #pragma acc parallel loop present(valt[0:gs3],valT2x[0:iN][0:gs3])
        for (int j=0;j<gs3;j++)
          valT2x[ii1][j] += valt[j];
      }

      reduce_2c2d(m,n,s1,s2,s3,s4,gs,norms,GFao,valS1,valS2,valS3,valS4,valS1x,valS2x,valS3x,valS4x,iN,N,natoms,-1.,xyz_grad);
      reduce_2c2d(m,n,s1,s2,s3,s4,gs,norms,Pao, valT1,valT2,valS3,valS4,valT1x,valT2x,valS3x,valS4x,iN,N,natoms,0.5,xyz_grad);

    } //loop n over second atom

  } //loop m over natoms

 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  //#pragma acc exit data delete(Pao[0:N2],GFao[0:N2])
  #pragma acc exit data copyout(xyz_grad[0:N3])
 #endif

#if DEBUG
  if (prl>1)
  {
    printf("\n xyz_grad: \n");
    for (int i=0;i<natoms;i++)
    {
      for (int j=0;j<3;j++)
        printf(" %8.5f",xyz_grad[3*i+j]);
      printf("\n");
    }
  }
#endif

#if USE_ACC
  #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
  #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs]) 
  #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6],wtt1[0:gs])
  #pragma acc exit data delete(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs])
  #pragma acc exit data delete(valT1[0:iN][0:gs],valT2[0:iN][0:gs])
  #pragma acc exit data delete(valS1x[0:iN][0:gs3],valS2x[0:iN][0:gs3],valS3x[0:iN][0:gs3],valS4x[0:iN][0:gs3])
  #pragma acc exit data delete(valT1x[0:iN][0:gs3],valT2x[0:iN][0:gs3])
  #pragma acc exit data delete(valt[0:gs3])
  #pragma acc exit data delete(norms[0:N2],n2i[0:natoms])
  #pragma acc exit data delete(coords[0:3*natoms],atno[0:natoms])
#endif

  delete [] ang_g;
  delete [] ang_w;

  delete [] n2i;
  delete [] norms;

  for (int i=0;i<iN;i++) delete [] valS1[i];
  for (int i=0;i<iN;i++) delete [] valS2[i];
  for (int i=0;i<iN;i++) delete [] valS3[i];
  for (int i=0;i<iN;i++) delete [] valS4[i];
  for (int i=0;i<iN;i++) delete [] valT1[i];
  for (int i=0;i<iN;i++) delete [] valT2[i];
  delete [] valS1; delete [] valS2; delete [] valS3; delete [] valS4;
  delete [] valT1; delete [] valT2;
  for (int i=0;i<iN;i++) delete [] valS1x[i];
  for (int i=0;i<iN;i++) delete [] valS2x[i];
  for (int i=0;i<iN;i++) delete [] valS3x[i];
  for (int i=0;i<iN;i++) delete [] valS4x[i];
  for (int i=0;i<iN;i++) delete [] valT1x[i];
  for (int i=0;i<iN;i++) delete [] valT2x[i];
  delete [] valS1x; delete [] valS2x; delete [] valS3x; delete [] valS4x;
  delete [] valT1x; delete [] valT2x;
  delete [] valt;

  delete [] grid1s;
  delete [] grid2s;
  delete [] wtt1;


  delete [] grid1;
  delete [] grid2;
  delete [] wt1;
  delete [] wt2;

  return;
}

#if RED_DOUBLE
void compute_d_2c(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* dpq, FP2* xyz_grad, int prl)
#else
void compute_d_2c(int natoms, int* atno, FP1* coords, vector<vector<FP2> > &basis, int nrad, int nang, FP2* ang_g0, FP2* ang_w0, FP2* dpq, FP1* xyz_grad, int prl)
#endif
{
 //nuclear derivatives of 2c Coulomb integrals
  if (prl>1) printf(" compute_d_2c \n");

 //2c integrals are all in auxiliary basis
  int N = basis.size();
  int N2 = N*N;
  int N3 = 3*natoms;

  int estart = find_center_of_grid(1,nrad)*nang;

  int gs = nrad*nang;
  int gs3 = 3*gs;
  int gs6 = 2*gs3;
  FP1* grid1 = new FP1[gs6];
  FP1* wt1 = new FP1[gs];

  FP1* grid2 = new FP1[gs6];
  FP1* wt2 = new FP1[gs];

  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);
  //printf("  iN: %i \n",imaxN);

 //intermediate storage
  int iN = imaxN;
  FP1* grid1s = new FP1[gs6];
  FP1* grid2s = new FP1[gs6];
  FP1** val1 = new FP1*[iN];
  FP1** val2 = new FP1*[iN];
  FP1** val3 = new FP1*[iN];
  FP1** val4 = new FP1*[iN];
  for (int i=0;i<iN;i++) val1[i] = new FP1[gs];
  for (int i=0;i<iN;i++) val2[i] = new FP1[gs];
  for (int i=0;i<iN;i++) val3[i] = new FP1[gs];
  for (int i=0;i<iN;i++) val4[i] = new FP1[gs];
  FP1** val1x = new FP1*[iN];
  FP1** val2x = new FP1*[iN];
  FP1** val3x = new FP1*[iN];
  FP1** val4x = new FP1*[iN];
  for (int i=0;i<iN;i++) val1x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++) val2x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++) val3x[i] = new FP1[gs3];
  for (int i=0;i<iN;i++) val4x[i] = new FP1[gs3];
  FP1* valt = new FP1[gs3];
  FP1* wtt1 = new FP1[gs];

  FP2* norms = new FP2[N2];
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    norms[i*N+j] = norm_sv(basis[i][0],basis[i][1],basis[i][2],basis[i][3])*basis[j][4];

  FP1* ang_g = new FP1[3*nang];
  FP1* ang_w = new FP1[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

 #if USE_ACC
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data copyin(n2i[0:natoms])

  #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
  #pragma acc enter data create(grid2[0:gs6],wt2[0:gs]) 
  #pragma acc enter data create(wtt1[0:gs])
  #pragma acc enter data create(val1[0:iN][0:gs],val2[0:iN][0:gs])
  #pragma acc enter data create(val3[0:iN][0:gs],val4[0:iN][0:gs])
  #pragma acc enter data create(val1x[0:iN][0:gs3],val2x[0:iN][0:gs3])
  #pragma acc enter data create(val3x[0:iN][0:gs3],val4x[0:iN][0:gs3])
  #pragma acc enter data create(valt[0:gs3])
  #pragma acc enter data create(grid1s[0:gs6],grid2s[0:gs6])
  #pragma acc enter data create(xyz_grad[0:N3])
  #pragma acc enter data copyin(norms[0:N2]) 
  //,dpq[0:N2])
 #endif
  acc_assign(N3,xyz_grad,0.);

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    FP1 Z1 = (FP1)atno[m];
    FP1 A1 = coords[3*m+0]; FP1 B1 = coords[3*m+1]; FP1 C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);

  #if USE_ACC
   #pragma acc parallel loop present(val1[0:iN][0:gs],val1x[0:iN][0:gs3])
  #endif
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val1[ii1][j] = 1.f;
     #pragma acc loop
      for (int j=0;j<gs3;j++)
        val1x[ii1][j] = 1.f;
    }

  #if USE_ACC
   #pragma acc parallel loop present(val3[0:iN][0:gs],wt1[0:gs],val3x[0:iN][0:gs3])
  #endif
    for (int ii2=0;ii2<s2-s1;ii2++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val3[ii2][j] = wt1[j];
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val3x[ii2][3*j+0] = val3x[ii2][3*j+1] = val3x[ii2][3*j+2] = wt1[j];
    }
 

   //two-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      FP1 Z2 = (FP1)atno[n];
      FP1 A2 = coords[3*n+0]; FP1 B2 = coords[3*n+1]; FP1 C2 = coords[3*n+2];
      FP1 A12 = A2-A1; FP1 B12 = B2-B1; FP1 C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1); 
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid2,wt2,Z1,Z2,A12,B12,C12);

      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wt2);

     //needs to happen after Becke weighting
      add_r1_to_grid(gs,grid2,0.,0.,0.);

    #if USE_ACC
     #pragma acc parallel loop present(val1[0:iN][0:gs],val2[0:iN][0:gs],val1x[0:iN][0:gs3],val2x[0:iN][0:gs3])
    #endif
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val1[ii1][j] = val2[ii1][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          val1x[ii1][j] = val2x[ii1][j] = 1.f;
      }

    #if USE_ACC
     #pragma acc parallel loop present(val3[0:iN][0:gs],val4[0:iN][0:gs],val3x[0:iN][0:gs3],val4x[0:iN][0:gs3],wtt1[0:gs],wt2[0:gs])
    #endif
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val3[ii2][j] = wtt1[j]; 
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val3x[ii2][3*j+0] = val3x[ii2][3*j+1] = val3x[ii2][3*j+2] = wtt1[j];
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val4[ii2][j] = wt2[j];
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val4x[ii2][3*j+0] = val4x[ii2][3*j+1] = val4x[ii2][3*j+2] = wt2[j];
      }
 
   #if 1
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

        eval_inr_r12(gs,grid1,val1[ii1],n1,l1,zeta1);
        eval_inr_r12(gs,grid2,val2[ii1],n1,l1,zeta1);

        eval_inr_d(gs,grid1,val1x[ii1],n1,l1,zeta1);
        eval_inr_d(gs,grid2,val2x[ii1],n1,l1,zeta1);

        acc_assign(gs,valt,1.);
        eval_sh_3r(gs,grid1,valt,n1,l1,m1);
       #pragma acc parallel loop present(val1x[0:iN][0:gs3],valt[0:gs])
        for (int j=0;j<gs;j++)
        {
          FP1 v1 = valt[j];
          val1x[ii1][3*j+0] *= v1; val1x[ii1][3*j+1] *= v1; val1x[ii1][3*j+2] *= v1;
        }

        acc_assign(gs,valt,1.);
        eval_sh_3r(gs,grid2,valt,n1,l1,m1);
       #pragma acc parallel loop present(val2x[0:iN][0:gs3],valt[0:gs])
        for (int j=0;j<gs;j++)
        {
          FP1 v1 = valt[j];
          val2x[ii1][3*j+0] *= v1; val2x[ii1][3*j+1] *= v1; val2x[ii1][3*j+2] *= v1;
        }

        if (l1>0)
        {
         #pragma acc parallel loop present(val1[0:iN][0:gs],valt[0:gs3])
          for (int j=0;j<gs;j++)
            valt[3*j+0] = valt[3*j+1] = valt[3*j+2] = val1[ii1][j];
          eval_dp_3r(gs,grid1,valt,n1,l1,m1);
         #pragma acc parallel loop present(valt[0:gs3],val1x[0:iN][0:gs3])
          for (int j=0;j<gs3;j++)
            val1x[ii1][j] += valt[j];

         #pragma acc parallel loop present(val2[0:iN][0:gs],valt[0:gs3])
          for (int j=0;j<gs;j++)
            valt[3*j+0] = valt[3*j+1] = valt[3*j+2] = val2[ii1][j];
          eval_dp_3r(gs,grid2,valt,n1,l1,m1);
         #pragma acc parallel loop present(valt[0:gs3],val2x[0:iN][0:gs3])
          for (int j=0;j<gs3;j++)
            val2x[ii1][j] += valt[j];
        }

        eval_sh_3r(gs,grid1,val1[ii1],n1,l1,m1);
        eval_sh_3r(gs,grid2,val2[ii1],n1,l1,m1);
      }

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];
        //printf("   n: %i i2: %i   nlm: %i %i %2i zeta: %8.5f \n",n,i2,n2,l2,m2,zeta2);

        eval_sh(ii2,gs,grid1s,val3[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2s,val4[ii2],n2,l2,m2,zeta2);

        eval_p(gs,grid1s,val3x[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid2s,val4x[ii2],n2,l2,m2,zeta2);
      }
 
    #if 0
      #pragma acc update self(wtt1[0:gs],wt2[0:gs])
      printf("\n");
      printf(" debug info: \n");
      printf(" wtt1:");
      for (int j=0;j<gs;j++)
        printf(" %6.3f",wtt1[j]);
      printf("\n");
      printf(" wt2:");
      for (int j=0;j<gs;j++)
        printf(" %6.3f",wt2[j]);
      printf("\n");

      #pragma acc update self(val1[0:iN][0:gs],val2[0:iN][0:gs])
      printf("\n");
      printf(" val1:");
      for (int j=0;j<gs;j++)
        printf(" %6.3f",val1[0][j]);
      printf("\n");
      printf(" val2:");
      for (int j=0;j<gs;j++)
        printf(" %6.3f",val2[0][j]);
      printf("\n");

      #pragma acc update self(val3[0:iN][0:gs],val4[0:iN][0:gs])
      printf(" val3:");
      for (int j=0;j<gs;j++)
        printf(" %6.3f",val3[0][j]); // /wtt1[j]);
      printf("\n");
      printf(" val4:");
      for (int j=0;j<gs;j++)
        printf(" %6.3f",val4[0][j]); // /wt2[j]);
      printf("\n\n");

      #pragma acc update self(val1x[0:iN][0:gs3],val2x[0:iN][0:gs3])
      printf(" val1x:");
      for (int i=0;i<3;i++)
      {
        for (int j=0;j<gs;j++)
          printf(" %6.3f",val1x[0][j*3+i]);
        printf("\n");
      }
      printf("\n");
      printf(" val2x:");
      for (int i=0;i<3;i++)
      {
        for (int j=0;j<gs;j++)
          printf(" %6.3f",val2x[0][j*3+i]);
        printf("\n");
      }
      printf("\n");

      #pragma acc update self(val3x[0:iN][0:gs3],val4x[0:iN][0:gs3])
      printf(" val3x:");
      for (int i=0;i<3;i++)
      {
        for (int j=0;j<gs;j++)
          printf(" %6.3f",val3x[0][j*3+i]); // /wtt1[j]);
        printf("\n");
      }
      printf("\n");
      printf(" val4x:");
      for (int i=0;i<3;i++)
      {
        for (int j=0;j<gs;j++)
          printf(" %6.3f",val4x[0][j*3+i]); // /wt2[j]);
        printf("\n");
      }
      printf("\n");
    #endif

      reduce_2c2d(m,n,s1,s2,s3,s4,gs,norms,dpq,val1,val2,val3,val4,val1x,val2x,val3x,val4x,iN,N,natoms,0.5,xyz_grad);

    #else

      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

        eval_inr_d(gs,grid1,val1x[ii1],n1,l1,zeta1);
        eval_inr_d(gs,grid2,val2x[ii1],n1,l1,zeta1);

        acc_assign(gs,valt,1.);
        eval_sh_3r(gs,grid1,valt,n1,l1,m1);
       #pragma acc parallel loop present(val1x[0:iN][0:gs3],valt[0:gs])
        for (int j=0;j<gs;j++)
        {
          FP1 v1 = valt[j];
          val1x[ii1][3*j+0] *= v1; val1x[ii1][3*j+1] *= v1; val1x[ii1][3*j+2] *= v1;
        }

        acc_assign(gs,valt,1.);
        eval_sh_3r(gs,grid2,valt,n1,l1,m1);
       #pragma acc parallel loop present(val2x[0:iN][0:gs3],valt[0:gs])
        for (int j=0;j<gs;j++)
        {
          FP1 v1 = valt[j];
          val2x[ii1][3*j+0] *= v1; val2x[ii1][3*j+1] *= v1; val2x[ii1][3*j+2] *= v1;
        }

        eval_inr_r12(gs,grid1,val1[ii1],n1,l1,zeta1);
        eval_inr_r12(gs,grid2,val2[ii1],n1,l1,zeta1);
        eval_sh_3r(gs,grid1,val1[ii1],n1,l1,m1);
        eval_sh_3r(gs,grid2,val2[ii1],n1,l1,m1);
      } //loop i1

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<FP2> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; FP1 zeta2 = basis2[3];
        //printf("   n: %i i2: %i   nlm: %i %i %2i zeta: %8.5f \n",n,i2,n2,l2,m2,zeta2);

        eval_sh(ii2,gs,grid1s,val3[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2s,val4[ii2],n2,l2,m2,zeta2);

        eval_p(gs,grid1s,val3x[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid2s,val4x[ii2],n2,l2,m2,zeta2);
      }

      reduce_2c2d(m,n,s1,s2,s3,s4,gs,norms,dpq,val1,val2,val3,val4,val1x,val2x,val3x,val4x,iN,N,natoms,0.5,xyz_grad);

     //CPMZ
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<FP2> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; FP1 zeta1 = basis1[3];

        if (l1>0)
        {
          acc_assign(gs,val1[ii1],1.);
          acc_assign(gs,val2[ii1],1.);
          eval_inr_r12(gs,grid1,val1[ii1],n1,l1,zeta1);
          eval_inr_r12(gs,grid2,val2[ii1],n1,l1,zeta1);

          acc_assign(gs3,val1x[ii1],1.);
          acc_assign(gs3,val2x[ii1],1.);

         #pragma acc parallel loop present(val1[0:iN][0:gs],valt[0:gs3])
          for (int j=0;j<gs;j++)
            valt[3*j+0] = valt[3*j+1] = valt[3*j+2] = val1[ii1][j];
          eval_dp_3r(gs,grid1,valt,n1,l1,m1);
         #pragma acc parallel loop present(valt[0:gs3],val1x[0:iN][0:gs3])
          for (int j=0;j<gs3;j++)
            val1x[ii1][j] = valt[j];

         #pragma acc parallel loop present(val2[0:iN][0:gs],valt[0:gs3])
          for (int j=0;j<gs;j++)
            valt[3*j+0] = valt[3*j+1] = valt[3*j+2] = val2[ii1][j];
          eval_dp_3r(gs,grid2,valt,n1,l1,m1);
         #pragma acc parallel loop present(valt[0:gs3],val2x[0:iN][0:gs3])
          for (int j=0;j<gs3;j++)
            val2x[ii1][j] = valt[j];

          eval_sh_3r(gs,grid1,val1[ii1],n1,l1,m1);
          eval_sh_3r(gs,grid2,val2[ii1],n1,l1,m1);
        }
        else
        {
          acc_assign(gs3,val1x[ii1],0.);
          acc_assign(gs3,val2x[ii1],0.);
        }
      } //loop i1

      reduce_2c2dh(m,n,s3,s4,s1,s2,gs,norms,dpq,val3,val4,val1x,val2x,iN,N,natoms,0.5,xyz_grad);
    #endif


    } //loop n over second atom

  } //loop m over natoms

 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc exit data copyout(xyz_grad[0:N3])
 #endif

#if DEBUG
  if (prl>1)
  {
    printf("\n xyz_grad: \n");
    for (int i=0;i<natoms;i++)
    {
      for (int j=0;j<3;j++)
        printf(" %8.5f",xyz_grad[3*i+j]);
      printf("\n");
    }
  }
#endif

  //clean_small_values(N,xyz_grad);

#if USE_ACC
  #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
  #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs]) 
  #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6])
  #pragma acc exit data delete(val1[0:iN][0:gs],val2[0:iN][0:gs],wtt1[0:gs])
  #pragma acc exit data delete(val3[0:iN][0:gs],val4[0:iN][0:gs])
  #pragma acc exit data delete(val1x[0:iN][0:gs3],val2x[0:iN][0:gs3])
  #pragma acc exit data delete(val3x[0:iN][0:gs3],val4x[0:iN][0:gs3])
  #pragma acc exit data delete(valt[0:gs3])
  #pragma acc exit data delete(n2i[0:natoms])
  #pragma acc exit data delete(norms[0:N2])
#endif

  delete [] norms;
  delete [] ang_g;
  delete [] ang_w;

  delete [] n2i;

  delete [] grid1s;
  delete [] grid2s;

  for (int i=0;i<iN;i++) delete [] val1[i];
  for (int i=0;i<iN;i++) delete [] val2[i];
  for (int i=0;i<iN;i++) delete [] val3[i];
  for (int i=0;i<iN;i++) delete [] val4[i];
  for (int i=0;i<iN;i++) delete [] val1x[i];
  for (int i=0;i<iN;i++) delete [] val2x[i];
  for (int i=0;i<iN;i++) delete [] val3x[i];
  for (int i=0;i<iN;i++) delete [] val4x[i];
  delete [] val1; delete [] val2; delete [] val3; delete [] val4;
  delete [] val1x; delete [] val2x; delete [] val3x; delete [] val4x;
  delete [] valt;
  delete [] wtt1;

  delete [] grid1;
  delete [] grid2;
  delete [] wt1;
  delete [] wt2;

  return;
}
