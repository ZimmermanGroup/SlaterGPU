#include "integrals.h"

#define TEST_SORT 0
//symmetrize wrt atom swap
#define SYMM_ST 1

#define RPAD 1.e-10

/*
 //current status of compute_all_2/3c code:
  1. add n=5 functions (need Inr)
  2. optimization of reduce functions

  4. need to check 4s,4p,4d functions for accuracy
  5. Cartesian d functions


 //completed
  1. implement d+f functions (done)
  1. b. implement g+h functions (done)
  2. need En/pVp --> 3c integrations (done)
  3. 4f functions have some kind of bug (fixed)
  4. 5g functions m==0 component has bug (high angular momentum sensitive to small weights in grid, fixed)
  5. 3c integration
  6. 2s, 3s, 3p functions + derivatives
  7. 4s, 4p, 4d functions

*/

#include <string>
void auto_crash();
void print_duration(chrono::high_resolution_clock::time_point t1, chrono::high_resolution_clock::time_point t2, string name);
void print_grid(int gs, float* grid1, float* grid2, float* wt1, float* wt2, int prl);
bool close_val(double v1, double v2);
void print_square_fine(int size, double* A);

void collect_4c_1d(int s1, int s2, int s3, int s4, int gs, int gsp, int M, int N, float* gt, double* g)
{
 //2221
  int M2 = M*M;
  int M3 = M2*M;
  int N2 = N*N;
  int N3 = N2*N;

  //printf(" collect_4c_1d: %i %i %i %i \n",s1,s2,s3,s4);

  #pragma acc update self(gt[0:M2*M2])
  for (int i1=s3;i1<s4;i1++)
  for (int i2=s3;i2<=i1;i2++)
  for (int i3=s3;i3<s4;i3++)
  for (int i4=s1;i4<s2;i4++)
  {
    int ii1 = i1-s3; int ii2 = i2-s3;
    int ii3 = i3-s3; int ii4 = i4-s1;

    float v1 = gt[ii1*M3+ii2*M2+ii3*M+ii4];
    gt[ii1*M3+ii2*M2+ii3*M+ii4] = 0.f;

    g[i1*N3+i2*N2+i3*N+i4] = v1;
    g[i2*N3+i1*N2+i3*N+i4] = v1;

    g[i1*N3+i2*N2+i4*N+i3] = v1;
    g[i2*N3+i1*N2+i4*N+i3] = v1;

    g[i3*N3+i4*N2+i1*N+i2] = v1;
    g[i3*N3+i4*N2+i2*N+i1] = v1;

    g[i4*N3+i3*N2+i1*N+i2] = v1;
    g[i4*N3+i3*N2+i2*N+i1] = v1;
  }

  #pragma acc parallel loop present(gt[0:M2*M2])
  for (int j=0;j<M2*M2;j++)
    gt[j] = 0.f;

  return;
}


void collect_4c_1c(int s1, int s2, int s3, int s4, int gs, int gsp, int M, int N, float* gt, double* g)
{
 //1221
  int M2 = M*M;
  int M3 = M2*M;
  int N2 = N*N;
  int N3 = N2*N;

  //printf(" collect_4c_1c: %i %i %i %i \n",s1,s2,s3,s4);

  #pragma acc update self(gt[0:M2*M2])
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s3;i2<s4;i2++)
  for (int i3=s3;i3<s4;i3++)
  for (int i4=s1;i4<s2;i4++)
  {
    int ii1 = i1-s1; int ii2 = i2-s3;
    int ii3 = i3-s3; int ii4 = i4-s1;

    float v1 = gt[ii1*M3+ii2*M2+ii3*M+ii4];
    gt[ii1*M3+ii2*M2+ii3*M+ii4] = 0.f;

    g[i1*N3+i2*N2+i3*N+i4] = v1;
    g[i1*N3+i2*N2+i4*N+i3] = v1;
    g[i2*N3+i1*N2+i3*N+i4] = v1;
    g[i2*N3+i1*N2+i4*N+i3] = v1;
  }

  #pragma acc parallel loop present(gt[0:M2*M2])
  for (int j=0;j<M2*M2;j++)
    gt[j] = 0.f;

  return;
}


void collect_4c_1b(int s1, int s2, int s3, int s4, int gs, int gsp, int M, int N, float* gt, double* g)
{
 //1121
  int M2 = M*M;
  int M3 = M2*M;
  int N2 = N*N;
  int N3 = N2*N;

  //printf(" collect_4c_1b: %i %i %i %i \n",s1,s2,s3,s4);

  #pragma acc update self(gt[0:M2*M2])
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s1;i2<=i1;i2++)
  for (int i3=s3;i3<s4;i3++)
  for (int i4=s1;i4<s2;i4++)
  {
    int ii1 = i1-s1; int ii2 = i2-s1;
    int ii3 = i3-s3; int ii4 = i4-s1;

    float v1 = gt[ii1*M3+ii2*M2+ii3*M+ii4];
    gt[ii1*M3+ii2*M2+ii3*M+ii4] = 0.f;

    g[i1*N3+i2*N2+i3*N+i4] = v1;
    g[i2*N3+i1*N2+i3*N+i4] = v1;
    g[i3*N3+i4*N2+i1*N+i2] = v1;
    g[i3*N3+i4*N2+i2*N+i1] = v1;

    g[i1*N3+i2*N2+i4*N+i3] = v1;
    g[i2*N3+i1*N2+i4*N+i3] = v1;
    g[i4*N3+i3*N2+i1*N+i2] = v1;
    g[i4*N3+i3*N2+i2*N+i1] = v1;
  }

  #pragma acc parallel loop present(gt[0:M2*M2])
  for (int j=0;j<M2*M2;j++)
    gt[j] = 0.f;

  return;
}

void collect_4c_1(int s1, int s2, int s3, int s4, bool lr_copy, int gs, int gsp, int M, int N, float* gt, double* g)
{
  int M2 = M*M;
  int M3 = M2*M;
  int N2 = N*N;
  int N3 = N2*N;

  //printf(" collect_4c_1: %i %i %i %i \n",s1,s2,s3,s4);

  #pragma acc update self(gt[0:M2*M2])
  for (int i1=s1;i1<s2;i1++)
  for (int i2=s1;i2<=i1;i2++)
  for (int i3=s3;i3<s4;i3++)
  for (int i4=s3;i4<=i3;i4++)
  {
    int ii1 = i1-s1; int ii2 = i2-s1;
    int ii3 = i3-s3; int ii4 = i4-s3;

    float v1 = gt[ii1*M3+ii2*M2+ii3*M+ii4];
    gt[ii1*M3+ii2*M2+ii3*M+ii4] = 0.f;

    g[i1*N3+i2*N2+i3*N+i4] = v1;
    g[i2*N3+i1*N2+i3*N+i4] = v1;
    g[i1*N3+i2*N2+i4*N+i3] = v1;
    g[i2*N3+i1*N2+i4*N+i3] = v1;
    if (lr_copy)
    {
      g[i3*N3+i4*N2+i1*N+i2] = v1;
      g[i3*N3+i4*N2+i2*N+i1] = v1;
      g[i4*N3+i3*N2+i1*N+i2] = v1;
      g[i4*N3+i3*N2+i2*N+i1] = v1;
    }
  }

  #pragma acc parallel loop present(gt[0:M2*M2])
  for (int j=0;j<M2*M2;j++)
    gt[j] = 0.f;

  return;
}


void print_grid(int gs, float* grid1, float* grid2, float* wt1, float* wt2, int prl)
{
  if (gs>300 && prl<3) return;
  if (prl>0)
  {
    //printf("\n printing grid1 \n");
    for (int i=0;i<gs;i++)
      //printf(" %8.5f %8.5f %8.5f  wt: %10.6f  r1/2: %5.2f %5.2f \n",grid1[6*i+0],grid1[6*i+1],grid1[6*i+2],wt1[i],grid1[6*i+3],grid1[6*i+4]);
      printf(" %9.5f %9.5f %9.5f  wt: %7.4e  r1/2: %5.2f %5.2f \n",grid1[6*i+0],grid1[6*i+1],grid1[6*i+2],wt1[i],grid1[6*i+3],grid1[6*i+4]);
    if (grid2!=NULL && wt2!=NULL)
    {
      //printf("\n printing grid2 \n");
      for (int i=0;i<gs;i++)
        //printf(" %8.5f %8.5f %8.5f  wt: %10.6f  r1/2: %5.2f %5.2f \n",grid2[6*i+0],grid2[6*i+1],grid2[6*i+2],wt2[i],grid2[6*i+3],grid2[6*i+4]);
        printf(" %9.5f %9.5f %9.5f  wt: %7.4e  r1/2: %5.2f %5.2f \n",grid2[6*i+0],grid2[6*i+1],grid2[6*i+2],wt2[i],grid2[6*i+3],grid2[6*i+4]);
    }
  }
  return;
}

void print_grid(int gs, float* grid1, float* grid2, float* grid3, float* wt1, float* wt2, float* wt3, int prl)
{
  if (gs>300) return;
  if (prl>0)
  {
    //printf("\n printing grid1 \n");
    for (int i=0;i<gs;i++)
      //printf(" %8.5f %8.5f %8.5f  wt: %10.6f  r123: %5.2f %5.2f %5.2f \n",grid1[6*i+0],grid1[6*i+1],grid1[6*i+2],wt1[i],grid1[6*i+3],grid1[6*i+4],grid1[6*i+5]);
      printf(" %9.5f %9.5f %9.5f  wt: %7.4e  r123: %5.2f %5.2f %5.2f \n",grid1[6*i+0],grid1[6*i+1],grid1[6*i+2],wt1[i],grid1[6*i+3],grid1[6*i+4],grid1[6*i+5]);
    if (grid2!=NULL && wt2!=NULL)
    {
      //printf("\n printing grid2 \n");
      for (int i=0;i<gs;i++)
        //printf(" %8.5f %8.5f %8.5f  wt: %10.6f  r123: %5.2f %5.2f %5.2f \n",grid2[6*i+0],grid2[6*i+1],grid2[6*i+2],wt2[i],grid2[6*i+3],grid2[6*i+4],grid2[6*i+5]);
        printf(" %9.5f %9.5f %9.5f  wt: %7.4e  r123: %5.2f %5.2f %5.2f \n",grid2[6*i+0],grid2[6*i+1],grid2[6*i+2],wt2[i],grid2[6*i+3],grid2[6*i+4],grid2[6*i+5]);
    }
    if (grid3!=NULL && wt3!=NULL)
    {
      //printf("\n printing grid3 \n");
      for (int i=0;i<gs;i++)
        //printf(" %8.5f %8.5f %8.5f  wt: %10.6f  r123: %5.2f %5.2f %5.2f \n",grid3[6*i+0],grid3[6*i+1],grid3[6*i+2],wt3[i],grid3[6*i+3],grid3[6*i+4],grid3[6*i+5]);
        printf(" %9.5f %9.5f %9.5f  wt: %7.4e  r123: %5.2f %5.2f %5.2f \n",grid3[6*i+0],grid3[6*i+1],grid3[6*i+2],wt3[i],grid3[6*i+3],grid3[6*i+4],grid3[6*i+5]);
    }
  }
  return;
}

#if RED_DOUBLE
void compute_all_3c_para(int ngpu, bool do_overlap, int natoms, int* atno, float* coords, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, int nrad, int nang, double* ang_g0, double* ang_w0, double* C, int prl)
#else
void compute_all_3c_para(int ngpu, bool do_overlap, int natoms, int* atno, float* coords, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, int nrad, int nang, double* ang_g0, double* ang_w0, float* C, int prl)
#endif
{
 #if !USE_ACC
  return compute_all_3c_v2(do_overlap,natoms,atno,coords,basis,basis_aux,nrad,nang,ang_g0,ang_w0,C,prl);
 #endif

  if (ngpu<1) { printf("\n ERROR: cannot run compute_all_3c_para, not enough gpus (%i) \n",ngpu); exit(-1); }

  int nomp = ngpu;
 //#pragma omp parallel
  //nomp = omp_get_num_threads();

  if (prl>1) printf(" beginning compute_all_3c_para. nthreads: %2i \n",nomp);

  int N = basis.size();
  int N2 = N*N;
  int Naux = basis_aux.size();
  int N2a = N2*Naux;
  int gs = nrad*nang;
  int gs6 = 6*gs;

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

  int estart = find_center_of_grid(1,nrad)*nang;

  int* na2i = new int[natoms];
  int iNa = get_imax_n2i(natoms,Naux,basis_aux,na2i);

  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);

  //printf("  iN/a: %i %i \n",iN,iNa);

  float* grid1 = new float[gs6];
  float* wt1 = new float[gs];

  float* grid2 = new float[gs6];
  float* wt2 = new float[gs];

  float* grid3 = new float[gs6];
  float* wt3 = new float[gs];

  float** val1 = new float*[iNa];
  float** val2 = new float*[iN];
  float** val3 = new float*[iN];
  float** val4 = new float*[iNa];
  float** val5 = new float*[iN];
  float** val6 = new float*[iN];
  float** val7 = new float*[iNa];
  float** val8 = new float*[iN];
  float** val9 = new float*[iN];
  for (int i=0;i<iNa;i++) val1[i] = new float[gs];
  for (int i=0;i<iN;i++)  val2[i] = new float[gs];
  for (int i=0;i<iN;i++)  val3[i] = new float[gs];
  for (int i=0;i<iNa;i++) val4[i] = new float[gs];
  for (int i=0;i<iN;i++)  val5[i] = new float[gs];
  for (int i=0;i<iN;i++)  val6[i] = new float[gs];
  for (int i=0;i<iNa;i++) val7[i] = new float[gs];
  for (int i=0;i<iN;i++)  val8[i] = new float[gs];
  for (int i=0;i<iN;i++)  val9[i] = new float[gs];
  float* valt1 = new float[gs];
  float* valt2 = new float[gs];
  float* valt3 = new float[gs];

  float* grid1s = new float[gs6];
  float* grid2s = new float[gs6];
  float* grid3s = new float[gs6];

  float* grid1p = new float[gs6];
  float* grid2p = new float[gs6];
  float* grid3p = new float[gs6];

  float* wtt1 = new float[gs];
  float* wtt2 = new float[gs];

  float* ang_g = new float[3*nang];
  float* ang_w = new float[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

  //using namespace std::chrono;
  //high_resolution_clock::time_point t0 = high_resolution_clock::now();

 #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
    #pragma acc enter data copyin(n2i[0:natoms],na2i[0:natoms])

    #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
    #pragma acc enter data create(grid2[0:gs6],wt2[0:gs])
    #pragma acc enter data create(grid3[0:gs6],wt3[0:gs])
    #pragma acc enter data create(val1[0:iNa][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs])
    #pragma acc enter data create(val4[0:iNa][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs])
    #pragma acc enter data create(val7[0:iNa][0:gs],val8[0:iN][0:gs],val9[0:iN][0:gs])

    #pragma acc enter data create(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6])
    #pragma acc enter data create(grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6])

    #pragma acc enter data create(wtt1[0:gs],wtt2[0:gs],valt1[0:gs],valt2[0:gs],valt3[0:gs])
    if (tid>0)
    {
      #pragma acc enter data create(C[0:N2a])
    }
    acc_assign(N2a,C,0.);
  }
  acc_set_device_num(0,acc_device_nvidia);

  //high_resolution_clock::time_point t1 = high_resolution_clock::now();

 #pragma omp parallel for schedule(dynamic) num_threads(nomp)
  for (int m=0;m<natoms;m++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);
    //printf("  launch %i/%i \n",m,tid);

    int s1 = 0; if (m>0) s1 = na2i[m-1]; int s2 = na2i[m];
    int s3 = 0; if (m>0) s3 = n2i[m-1]; int s4 = n2i[m];

    float Z1 = (float)atno[m];
    float A1 = coords[3*m+0]; float B1 = coords[3*m+1]; float C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);

   #pragma acc parallel loop present(val1[0:iNa][0:gs])
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val1[ii1][j] = 1.f;
    }

   #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs])
    for (int ii2=0;ii2<s4-s3;ii2++)
    {
      for (int j=0;j<gs;j++)
        val2[ii2][j] = 1.f;
    }
   #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],wt1[0:gs])
    for (int ii3=0;ii3<s4-s3;ii3++)
    {
      for (int j=0;j<gs;j++)
        val3[ii3][j] = wt1[j];
    }

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis_aux[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

      if (do_overlap)
        eval_sh(ii1,gs,grid1,val1[ii1],n1,l1,m1,zeta1);
      else
      {
        eval_inr_r12(gs,grid1,val1[ii1],n1,l1,zeta1);
        eval_sh_3r(gs,grid1,val1[ii1],n1,l1,m1);
      }
    } //loop i1

    for (int i2=s3;i2<s4;i2++)
    {
      int ii2 = i2-s3;

      vector<double> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

      eval_sh(ii2,gs,grid1,val2[ii2],n2,l2,m2,zeta2);
    } //loop i2

    for (int i3=s3;i3<s4;i3++)
    {
      int ii3 = i3-s3;

      vector<double> basis3 = basis[i3];
      int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

      eval_sh(ii3,gs,grid1,val3[ii3],n3,l3,m3,zeta3);
    } //loop i3

    reduce_3c1b(s1,s2,s3,s4,gs,val1,val2,val3,N,Naux,iN,iNa,C);


   //two-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
     //s12 over atom m, aux function
     //s34 over atom m, regular function
     //s56 over atom n, regular function
      int s5 = 0; if (n>0) s5 = n2i[n-1]; int s6 = n2i[n];

      float Z2 = (float)atno[n];
      float A2 = coords[3*n+0]; float B2 = coords[3*n+1]; float C2 = coords[3*n+2];
      float A12 = A2-A1; float B12 = B2-B1; float C12 = C2-C1;

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

     #pragma acc parallel loop collapse(2) present(val4[0:iNa][0:gs])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
        for (int j=0;j<gs;j++)
          val4[ii1][j] = 1.f;
      }

     #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs],val5[0:iN][0:gs])
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
        for (int j=0;j<gs;j++)
        {
          val2[ii2][j] = 1.f;
          val5[ii2][j] = 1.f;
        }
      }

     #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],wtt1[0:gs],wt2[0:gs])
      for (int ii3=0;ii3<s6-s5;ii3++)
      {
        for (int j=0;j<gs;j++)
        {
          val3[ii3][j] = wtt1[j];
          val6[ii3][j] = wt2[j];
        }
      }

     //i1 on atom m
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis_aux[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        if (do_overlap)
          eval_sh(ii1,gs,grid2,val4[ii1],n1,l1,m1,zeta1);
        else
        {
          eval_inr_r12(gs,grid2,val4[ii1],n1,l1,zeta1);
          eval_sh_3r(gs,grid2,val4[ii1],n1,l1,m1);
        }
      }

     //i2 on atom m
      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;

        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1,val2[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2,val5[ii2],n2,l2,m2,zeta2);
      }

     //i3 on atom n
      for (int i3=s5;i3<s6;i3++)
      {
        int ii3 = i3-s5;

        vector<double> basis3 = basis[i3];
        int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

        eval_sh(ii3,gs,grid1s,val3[ii3],n3,l3,m3,zeta3);
        eval_sh(ii3,gs,grid2s,val6[ii3],n3,l3,m3,zeta3);
      }

      reduce_3c2b(s1,s2,s3,s4,s5,s6,gs,val1,val2,val3,val4,val5,val6,N,Naux,iN,iNa,C);

     #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs],val5[0:iN][0:gs])
      for (int ii2=0;ii2<s6-s5;ii2++)
      {
        for (int j=0;j<gs;j++)
        {
          val2[ii2][j] = 1.f;
          val5[ii2][j] = 1.f;
        }
      }

     #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],wtt1[0:gs],wt2[0:gs])
      for (int ii3=0;ii3<s6-s5;ii3++)
      {
        for (int j=0;j<gs;j++)
        {
          val3[ii3][j] = wtt1[j];
          val6[ii3][j] = wt2[j];
        }
      }

     //now do i2+i3 on atom n
     //i2 on atom n
      for (int i2=s5;i2<s6;i2++)
      {
        int ii2 = i2-s5;

        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1s,val2[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2s,val5[ii2],n2,l2,m2,zeta2);
      }

     //i3 on atom n
      for (int i3=s5;i3<s6;i3++)
      {
        int ii3 = i3-s5;

        vector<double> basis3 = basis[i3];
        int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

        eval_sh(ii3,gs,grid1s,val3[ii3],n3,l3,m3,zeta3);
        eval_sh(ii3,gs,grid2s,val6[ii3],n3,l3,m3,zeta3);
      }

      reduce_3c2b(s1,s2,s5,s6,s5,s6,gs,val1,val2,val3,val4,val5,val6,N,Naux,iN,iNa,C);

    } //loop n over second atom



   //three-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      float Z2 = (float)atno[n];
      float A2 = coords[3*n+0]; float B2 = coords[3*n+1]; float C2 = coords[3*n+2];
      float A12 = A2-A1; float B12 = B2-B1; float C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1);
      recenter_grid(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

     //testing condition
      for (int p=n+1;p<natoms;p++)
      if (p!=m)
      //for (int p=0;p<natoms;p++)
      //if (p!=m && p!=n)
      {
        int s5 = 0; if (p>0) s5 = n2i[p-1]; int s6 = n2i[p];

        float Z3 = (float)atno[p];
        float A3 = coords[3*p+0]; float B3 = coords[3*p+1]; float C3 = coords[3*p+2];
        float A13 = A3-A1; float B13 = B3-B1; float C13 = C3-C1;

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

        add_r1_to_grid(gs,grid2,0.,0.,0.);
        add_r1_to_grid(gs,grid3,0.,0.,0.);


       #pragma acc parallel loop collapse(2) present(val4[0:iNa][0:gs],val7[0:iNa][0:gs])
        for (int ii1=0;ii1<s2-s1;ii1++)
        {
          for (int j=0;j<gs;j++)
          {
            val4[ii1][j] = 1.f;
            val7[ii1][j] = 1.f;
          }
        }

       #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs],val5[0:iN][0:gs],val8[0:iN][0:gs])
        for (int ii2=0;ii2<s4-s3;ii2++)
        {
          for (int j=0;j<gs;j++)
          {
            val2[ii2][j] = 1.f;
            val5[ii2][j] = 1.f;
            val8[ii2][j] = 1.f;
          }
        }

       #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],val9[0:iN][0:gs],wtt1[0:gs],wtt2[0:gs],wt3[0:gs])
        for (int ii3=0;ii3<s6-s5;ii3++)
        {
          for (int j=0;j<gs;j++)
          {
            val3[ii3][j] = wtt1[j];
            val6[ii3][j] = wtt2[j];
            val9[ii3][j] = wt3[j];
          }
        }


        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<double> basis1 = basis_aux[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

          if (do_overlap)
            eval_sh(ii1,gs,grid2,val4[ii1],n1,l1,m1,zeta1);
          else
          {
            eval_inr_r12(gs,grid2,val4[ii1],n1,l1,zeta1);
            eval_sh_3r(gs,grid2,val4[ii1],n1,l1,m1);
          }

          if (do_overlap)
            eval_sh(ii1,gs,grid3,val7[ii1],n1,l1,m1,zeta1);
          else
          {
            eval_inr_r12(gs,grid3,val7[ii1],n1,l1,zeta1);
            eval_sh_3r(gs,grid3,val7[ii1],n1,l1,m1);
          }
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

          eval_sh(ii2,gs,grid3s,val8[ii2],n2,l2,m2,zeta2);
          eval_sh(ii2,gs,grid2s,val5[ii2],n2,l2,m2,zeta2);
          eval_sh(ii2,gs,grid1s,val2[ii2],n2,l2,m2,zeta2);
        }

        for (int i3=s5;i3<s6;i3++)
        {
          int ii3 = i3-s5;
          vector<double> basis3 = basis[i3];
          int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

          eval_sh(ii3,gs,grid3p,val9[ii3],n3,l3,m3,zeta3);
          eval_sh(ii3,gs,grid2p,val6[ii3],n3,l3,m3,zeta3);
          eval_sh(ii3,gs,grid1p,val3[ii3],n3,l3,m3,zeta3);
        }

        reduce_3c3b(s1,s2,s3,s4,s5,s6,gs,val1,val2,val3,val4,val5,val6,val7,val8,val9,N,Naux,iN,iNa,C);

      } //loop p over third atom
    } //loop n over second atom

  } //loop m over natoms

  //high_resolution_clock::time_point t2 = high_resolution_clock::now();

 //collect integrals from all GPUs
  double* C_all = new double[N2a]();

  for (int n=0;n<nomp;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc update self(C[0:N2a])

    for (int i=0;i<N2a;i++)
      C_all[i] += C[i];

    #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
    if (n>0)
    {
      #pragma acc exit data delete(C[0:N2a])
    }
  }
  acc_set_device_num(0,acc_device_nvidia);

  for (int i=0;i<N2a;i++)
    C[i] = C_all[i];

  delete [] C_all;

 //apply symmetry and normalization
  copy_symm(natoms,N,Naux,basis,basis_aux,C,1);

  if (do_overlap)
  {
    for (int i=0;i<Naux;i++)
    for (int j=0;j<N;j++)
    for (int k=0;k<N;k++)
      C[i*N2+j*N+k] *= basis_aux[i][4]*basis[j][4]*basis[k][4];
  }
  else
  {
    for (int i=0;i<Naux;i++)
    for (int j=0;j<N;j++)
    for (int k=0;k<N;k++)
      C[i*N2+j*N+k] *= norm_sv(basis_aux[i][0],basis_aux[i][1],basis_aux[i][2],basis_aux[i][3])*basis[j][4]*basis[k][4];
  }

  if (prl>2)
  {
    printf("\n C: \n");
    for (int i=0;i<Naux;i++)
    {
      printf(" i: %i \n",i);
      for (int j=0;j<N;j++)
      for (int k=0;k<N;k++)
        printf("  %12.6f",C[i*N2+j*N+k]);
      printf("\n");
    }
  }

  if (!do_overlap)
    transpose_C(Naux,N,C);
  #pragma acc update device(C[0:N2a])


 //CPMZ check this
 #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
    #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs])
    #pragma acc exit data delete(grid3[0:gs6],wt3[0:gs])

    #pragma acc exit data delete(val1[0:iNa][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs])
    #pragma acc exit data delete(val4[0:iNa][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs])
    #pragma acc exit data delete(val7[0:iNa][0:gs],val8[0:iN][0:gs],val9[0:iN][0:gs])

    #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6])
    #pragma acc exit data delete(grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6])

    #pragma acc exit data delete(wtt1[0:gs],wtt2[0:gs],valt1[0:gs],valt2[0:gs],valt3[0:gs])
    #pragma acc exit data delete(n2i[0:natoms],na2i[0:natoms])
  }

 //dominated by step2
  //high_resolution_clock::time_point t3 = high_resolution_clock::now();
  //print_duration(t0,t1,"v3 step1");
  //print_duration(t1,t2,"v3 step2");
  //print_duration(t2,t3,"v3 step3");
  //auto_crash();

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
  delete [] valt1;
  delete [] valt2;
  delete [] valt3;

  return;
}

#if RED_DOUBLE
void compute_all_3c_v2(bool do_overlap, int natoms, int* atno, float* coords, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, int nrad, int nang, double* ang_g0, double* ang_w0, double* C, int prl)
#else
void compute_all_3c_v2(bool do_overlap, int natoms, int* atno, float* coords, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, int nrad, int nang, double* ang_g0, double* ang_w0, float* C, int prl)
#endif
{
  if (prl>1)
  {
    if (do_overlap)
      printf(" beginning compute_all_3c_v2 (overlap) \n");
    else
      printf(" beginning compute_all_3c_v2 \n");
  }

  int N = basis.size();
  int N2 = N*N;
  int Naux = basis_aux.size();
  int N2a = N2*Naux;
  int gs = nrad*nang;
  int gs6 = 6*gs;

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

  int estart = find_center_of_grid(1,nrad)*nang;

  int* na2i = new int[natoms];
  int iNa = get_imax_n2i(natoms,Naux,basis_aux,na2i);

  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);

  //printf("  iN/a: %i %i \n",iN,iNa);

  float* grid1 = new float[gs6];
  float* wt1 = new float[gs];

  float* grid2 = new float[gs6];
  float* wt2 = new float[gs];

  float* grid3 = new float[gs6];
  float* wt3 = new float[gs];

  float** val1 = new float*[iNa];
  float** val2 = new float*[iN];
  float** val3 = new float*[iN];
  float** val4 = new float*[iNa];
  float** val5 = new float*[iN];
  float** val6 = new float*[iN];
  float** val7 = new float*[iNa];
  float** val8 = new float*[iN];
  float** val9 = new float*[iN];
  for (int i=0;i<iNa;i++) val1[i] = new float[gs];
  for (int i=0;i<iN;i++)  val2[i] = new float[gs];
  for (int i=0;i<iN;i++)  val3[i] = new float[gs];
  for (int i=0;i<iNa;i++) val4[i] = new float[gs];
  for (int i=0;i<iN;i++)  val5[i] = new float[gs];
  for (int i=0;i<iN;i++)  val6[i] = new float[gs];
  for (int i=0;i<iNa;i++) val7[i] = new float[gs];
  for (int i=0;i<iN;i++)  val8[i] = new float[gs];
  for (int i=0;i<iN;i++)  val9[i] = new float[gs];
  float* valt1 = new float[gs];
  float* valt2 = new float[gs];
  float* valt3 = new float[gs];

  float* grid1s = new float[gs6];
  float* grid2s = new float[gs6];
  float* grid3s = new float[gs6];

  float* grid1p = new float[gs6];
  float* grid2p = new float[gs6];
  float* grid3p = new float[gs6];

  float* wtt1 = new float[gs];
  float* wtt2 = new float[gs];

  float* ang_g = new float[3*nang];
  float* ang_w = new float[nang];
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

  #pragma acc enter data create(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6])
  #pragma acc enter data create(grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6])

  #pragma acc enter data create(wtt1[0:gs],wtt2[0:gs],valt1[0:gs],valt2[0:gs],valt3[0:gs])
  //#pragma acc enter data create(C[0:N2a])
 #endif
  acc_assign(N2a,C,0.);

  for (int m=0;m<natoms;m++)
  {
    int s1 = 0; if (m>0) s1 = na2i[m-1]; int s2 = na2i[m];
    int s3 = 0; if (m>0) s3 = n2i[m-1]; int s4 = n2i[m];
    //printf("  m: %i  s1/2->3/4: %i %i - %i %i \n",m,s1,s2,s3,s4);

    float Z1 = (float)atno[m];
    float A1 = coords[3*m+0]; float B1 = coords[3*m+1]; float C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);
    //#pragma acc update self(grid1[0:gs6],wt1[0:gs])
    //print_grid(gs,grid1,NULL,wt1,NULL,prl);

 //CPMZ need more of this
   #pragma acc parallel loop present(val1[0:iNa][0:gs])
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val1[ii1][j] = 1.f;
    }

   #pragma acc parallel loop present(val2[0:iN][0:gs])
    for (int ii2=0;ii2<s4-s3;ii2++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val2[ii2][j] = 1.f;
    }
   #pragma acc parallel loop present(val3[0:iN][0:gs],wt1[0:gs])
    for (int ii3=0;ii3<s4-s3;ii3++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val3[ii3][j] = wt1[j];
    }

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis_aux[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

      if (do_overlap)
        eval_sh(ii1,gs,grid1,val1[ii1],n1,l1,m1,zeta1);
      else
      {
        eval_inr_r12(gs,grid1,val1[ii1],n1,l1,zeta1);
        eval_sh_3r(gs,grid1,val1[ii1],n1,l1,m1);
      }
    } //loop i1

    for (int i2=s3;i2<s4;i2++)
    {
      int ii2 = i2-s3;

      vector<double> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

      eval_sh(ii2,gs,grid1,val2[ii2],n2,l2,m2,zeta2);
    } //loop i2

    for (int i3=s3;i3<s4;i3++)
    {
      int ii3 = i3-s3;

      vector<double> basis3 = basis[i3];
      int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

      eval_sh(ii3,gs,grid1,val3[ii3],n3,l3,m3,zeta3);
    } //loop i3

   #if REDUCEV==1
    reduce_3c1(s1,s2,s3,s4,gs,val1,val2,val3,valt1,N,Naux,iN,iNa,C);
   #else
    reduce_3c1b(s1,s2,s3,s4,gs,val1,val2,val3,N,Naux,iN,iNa,C);
   #endif


   //two-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
     //s12 over atom m, aux function
     //s34 over atom m, regular function
     //s56 over atom n, regular function
      int s5 = 0; if (n>0) s5 = n2i[n-1]; int s6 = n2i[n];
      //printf("  n: %i  s1/2->3/4->5/6: %i %i - %i %i - %i %i \n",n,s1,s2,s3,s4,s5,s6);

      float Z2 = (float)atno[n];
      float A2 = coords[3*n+0]; float B2 = coords[3*n+1]; float C2 = coords[3*n+2];
      float A12 = A2-A1; float B12 = B2-B1; float C12 = C2-C1;

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

      //#pragma acc update self(grid1[0:gs6],wt1[0:gs],grid2[0:gs6],wt2[0:gs])
      //print_grid(gs,grid1,grid2,wt1,wt2,prl);

     #pragma acc parallel loop collapse(2) present(val4[0:iNa][0:gs])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
        for (int j=0;j<gs;j++)
          val4[ii1][j] = 1.f;
      }

     #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs],val5[0:iN][0:gs])
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
        for (int j=0;j<gs;j++)
        {
          val2[ii2][j] = 1.f;
          val5[ii2][j] = 1.f;
        }
      }

     #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],wtt1[0:gs],wt2[0:gs])
      for (int ii3=0;ii3<s6-s5;ii3++)
      {
        for (int j=0;j<gs;j++)
        {
          val3[ii3][j] = wtt1[j];
          val6[ii3][j] = wt2[j];
        }
      }

     //i1 on atom m
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis_aux[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];
        //printf("  m: %i i1: %i   nlm: %i %i %2i zeta: %8.5f \n",m,i1,n1,l1,m1,zeta1);

        if (do_overlap)
          eval_sh(ii1,gs,grid2,val4[ii1],n1,l1,m1,zeta1);
        else
        {
          eval_inr_r12(gs,grid2,val4[ii1],n1,l1,zeta1);
          eval_sh_3r(gs,grid2,val4[ii1],n1,l1,m1);
        }
      }

     //i2 on atom m
      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;

        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];
        //printf("   n: %i i2: %i   nlm: %i %i %2i zeta: %8.5f \n",n,i2,n2,l2,m2,zeta2);

        eval_sh(ii2,gs,grid1,val2[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2,val5[ii2],n2,l2,m2,zeta2);
      }

     //i3 on atom n
      for (int i3=s5;i3<s6;i3++)
      {
        int ii3 = i3-s5;

        vector<double> basis3 = basis[i3];
        int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

        eval_sh(ii3,gs,grid1s,val3[ii3],n3,l3,m3,zeta3);
        eval_sh(ii3,gs,grid2s,val6[ii3],n3,l3,m3,zeta3);
      }

     #if REDUCEV==1
      reduce_3c2(s1,s2,s3,s4,s5,s6,gs,val1,val2,val3,val4,val5,val6,valt1,valt2,N,Naux,iN,iNa,C);
     #else
      reduce_3c2b(s1,s2,s3,s4,s5,s6,gs,val1,val2,val3,val4,val5,val6,N,Naux,iN,iNa,C);
     #endif

     #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs],val5[0:iN][0:gs])
      for (int ii2=0;ii2<s6-s5;ii2++)
      {
        for (int j=0;j<gs;j++)
        {
          val2[ii2][j] = 1.f;
          val5[ii2][j] = 1.f;
        }
      }

     #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],wtt1[0:gs],wt2[0:gs])
      for (int ii3=0;ii3<s6-s5;ii3++)
      {
        for (int j=0;j<gs;j++)
        {
          val3[ii3][j] = wtt1[j];
          val6[ii3][j] = wt2[j];
        }
      }

     //now do i2+i3 on atom n
     //i2 on atom n
      for (int i2=s5;i2<s6;i2++)
      {
        int ii2 = i2-s5;

        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];
        //printf("   n: %i i2: %i   nlm: %i %i %2i zeta: %8.5f \n",n,i2,n2,l2,m2,zeta2);

        eval_sh(ii2,gs,grid1s,val2[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2s,val5[ii2],n2,l2,m2,zeta2);
      }

     //i3 on atom n
      for (int i3=s5;i3<s6;i3++)
      {
        int ii3 = i3-s5;

        vector<double> basis3 = basis[i3];
        int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

        eval_sh(ii3,gs,grid1s,val3[ii3],n3,l3,m3,zeta3);
        eval_sh(ii3,gs,grid2s,val6[ii3],n3,l3,m3,zeta3);
      }

     #if REDUCEV==1
      reduce_3c2(s1,s2,s5,s6,s5,s6,gs,val1,val2,val3,val4,val5,val6,valt1,valt2,N,Naux,iN,iNa,C);
     #else
      reduce_3c2b(s1,s2,s5,s6,s5,s6,gs,val1,val2,val3,val4,val5,val6,N,Naux,iN,iNa,C);
     #endif

    } //loop n over second atom



   //three-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      float Z2 = (float)atno[n];
      float A2 = coords[3*n+0]; float B2 = coords[3*n+1]; float C2 = coords[3*n+2];
      float A12 = A2-A1; float B12 = B2-B1; float C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1);
      recenter_grid(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

     //testing condition
      for (int p=n+1;p<natoms;p++)
      if (p!=m)
      //for (int p=0;p<natoms;p++)
      //if (p!=m && p!=n)
      {
        int s5 = 0; if (p>0) s5 = n2i[p-1]; int s6 = n2i[p];
        //printf("  p: %i  s1/2->3/4->5/6: %i %i - %i %i - %i %i \n",p,s1,s2,s3,s4,s5,s6);

        float Z3 = (float)atno[p];
        float A3 = coords[3*p+0]; float B3 = coords[3*p+1]; float C3 = coords[3*p+2];
        float A13 = A3-A1; float B13 = B3-B1; float C13 = C3-C1;

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

        add_r1_to_grid(gs,grid2,0.,0.,0.);
        add_r1_to_grid(gs,grid3,0.,0.,0.);


       #pragma acc parallel loop collapse(2) present(val4[0:iNa][0:gs],val7[0:iNa][0:gs])
        for (int ii1=0;ii1<s2-s1;ii1++)
        {
          for (int j=0;j<gs;j++)
          {
            val4[ii1][j] = 1.f;
            val7[ii1][j] = 1.f;
          }
        }

       #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gs],val5[0:iN][0:gs],val8[0:iN][0:gs])
        for (int ii2=0;ii2<s4-s3;ii2++)
        {
          for (int j=0;j<gs;j++)
          {
            val2[ii2][j] = 1.f;
            val5[ii2][j] = 1.f;
            val8[ii2][j] = 1.f;
          }
        }

       #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gs],val6[0:iN][0:gs],val9[0:iN][0:gs],wtt1[0:gs],wtt2[0:gs],wt3[0:gs])
        for (int ii3=0;ii3<s6-s5;ii3++)
        {
          for (int j=0;j<gs;j++)
          {
            val3[ii3][j] = wtt1[j];
            val6[ii3][j] = wtt2[j];
            val9[ii3][j] = wt3[j];
          }
        }


        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<double> basis1 = basis_aux[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

          if (do_overlap)
            eval_sh(ii1,gs,grid2,val4[ii1],n1,l1,m1,zeta1);
          else
          {
            eval_inr_r12(gs,grid2,val4[ii1],n1,l1,zeta1);
            eval_sh_3r(gs,grid2,val4[ii1],n1,l1,m1);
          }
          if (do_overlap)
            eval_sh(ii1,gs,grid3,val7[ii1],n1,l1,m1,zeta1);
          else
          {
            eval_inr_r12(gs,grid3,val7[ii1],n1,l1,zeta1);
            eval_sh_3r(gs,grid3,val7[ii1],n1,l1,m1);
          }
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

          eval_sh(ii2,gs,grid3s,val8[ii2],n2,l2,m2,zeta2);
          eval_sh(ii2,gs,grid2s,val5[ii2],n2,l2,m2,zeta2);
          eval_sh(ii2,gs,grid1s,val2[ii2],n2,l2,m2,zeta2);
        }

        for (int i3=s5;i3<s6;i3++)
        {
          int ii3 = i3-s5;
          vector<double> basis3 = basis[i3];
          int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

          eval_sh(ii3,gs,grid3p,val9[ii3],n3,l3,m3,zeta3);
          eval_sh(ii3,gs,grid2p,val6[ii3],n3,l3,m3,zeta3);
          eval_sh(ii3,gs,grid1p,val3[ii3],n3,l3,m3,zeta3);
        }

       #if REDUCEV==1
        reduce_3c3(s1,s2,s3,s4,s5,s6,gs,val1,val2,val3,val4,val5,val6,val7,val8,val9,valt1,valt2,valt3,N,Naux,iN,iNa,C);
       #else
        reduce_3c3b(s1,s2,s3,s4,s5,s6,gs,val1,val2,val3,val4,val5,val6,val7,val8,val9,N,Naux,iN,iNa,C);
       #endif

      } //loop p over third atom
    } //loop n over second atom

  } //loop m over natoms

 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  //#pragma acc exit data copyout(C[0:N2a])
  #pragma acc update self(C[0:N2a])
 #endif

  copy_symm(natoms,N,Naux,basis,basis_aux,C,1);

  if (do_overlap)
  {
    for (int i=0;i<Naux;i++)
    for (int j=0;j<N;j++)
    for (int k=0;k<N;k++)
      C[i*N2+j*N+k] *= basis_aux[i][4]*basis[j][4]*basis[k][4];
  }
  else
  {
    for (int i=0;i<Naux;i++)
    for (int j=0;j<N;j++)
    for (int k=0;k<N;k++)
      C[i*N2+j*N+k] *= norm_sv(basis_aux[i][0],basis_aux[i][1],basis_aux[i][2],basis_aux[i][3])*basis[j][4]*basis[k][4];
  }

  if (prl>2)
  {
    printf("\n C: \n");
    for (int i=0;i<Naux;i++)
    {
      printf(" i: %i \n",i);
      for (int j=0;j<N;j++)
      for (int k=0;k<N;k++)
        printf("  %12.6f",C[i*N2+j*N+k]);
      printf("\n");
    }
  }

  if (!do_overlap)
    transpose_C(Naux,N,C);
  #pragma acc update device(C[0:N2a])

 //CPMZ check this
#if USE_ACC
  #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
  #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs])
  #pragma acc exit data delete(grid3[0:gs6],wt3[0:gs])

  #pragma acc exit data delete(val1[0:iNa][0:gs],val2[0:iN][0:gs],val3[0:iN][0:gs])
  #pragma acc exit data delete(val4[0:iNa][0:gs],val5[0:iN][0:gs],val6[0:iN][0:gs])
  #pragma acc exit data delete(val7[0:iNa][0:gs],val8[0:iN][0:gs],val9[0:iN][0:gs])

  #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6])
  #pragma acc exit data delete(grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6])

  #pragma acc exit data delete(wtt1[0:gs],wtt2[0:gs],valt1[0:gs],valt2[0:gs],valt3[0:gs])
  #pragma acc exit data delete(n2i[0:natoms],na2i[0:natoms])
#endif

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
  delete [] valt1;
  delete [] valt2;
  delete [] valt3;

  return;
}

void compute_all_3c(int natoms, int* atno, float* coords, vector<vector<double> > &basis, vector<vector<double> > &basis_aux, int nrad, int nang, double* ang_g0, double* ang_w0, float* C, int prl)
{
  if (prl>1) printf(" beginning compute_all_3c \n");

  int N = basis.size();
  int N2 = N*N;
  int Naux = basis_aux.size();
  int N2a = N2*Naux;

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

  int gs = nrad*nang;

  float* grid1 = new float[6*gs];
  float* wt1 = new float[gs];
  float* val1 = new float[gs];

  float* grid2 = new float[6*gs];
  float* wt2 = new float[gs];
  float* val2 = new float[gs];

  float* grid3 = new float[6*gs];
  float* wt3 = new float[gs];
  float* val3 = new float[gs];

  float* grid1s = new float[6*gs];
  float* grid2s = new float[6*gs];
  float* grid3s = new float[6*gs];

  float* grid1p = new float[6*gs];
  float* grid2p = new float[6*gs];
  float* grid3p = new float[6*gs];

  float* valt1 = new float[gs];
  float* valt2 = new float[gs];
  float* valt3 = new float[gs];
  float* valt4 = new float[gs];
  float* valt5 = new float[gs];
  float* valt6 = new float[gs];
  float* wtt1 = new float[gs];
  float* wtt2 = new float[gs];

  float* ang_g = new float[3*nang];
  float* ang_w = new float[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

  int imaxN = 1;
  int* n2i = new int[natoms];
  int wa = 0;
  for (int i=0;i<N;i++)
  {
    int wa2 = basis[i][9];
    if (wa2!=wa)
    {
      int cmaxN = wa2-wa;
      if (cmaxN>imaxN) imaxN = cmaxN;
      n2i[wa] = i;
      wa = wa2;
    }
  }
  n2i[natoms-1] = N;

 #if USE_ACC
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data copyin(n2i[0:natoms])

  #pragma acc enter data create(grid1[0:6*gs],wt1[0:gs],val1[0:gs])
  #pragma acc enter data create(grid2[0:6*gs],wt2[0:gs],val2[0:gs])
  #pragma acc enter data create(grid3[0:6*gs],wt3[0:gs],val3[0:gs])

  #pragma acc enter data create(grid1s[0:6*gs],grid2s[0:6*gs],grid3s[0:6*gs])
  #pragma acc enter data create(grid1p[0:6*gs],grid2p[0:6*gs],grid3p[0:6*gs])

  #pragma acc enter data create(valt1[0:gs],valt2[0:gs],valt3[0:gs],valt4[0:gs],valt5[0:gs],valt6[0:gs])
  #pragma acc enter data create(wtt1[0:gs],wtt2[0:gs])
  #pragma acc enter data create(C[0:N2a])
 #endif
  acc_assign(N2a,C,0.);

  for (int m=0;m<natoms;m++)
  {
    //int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    float Z1 = (float)atno[m];
    float A1 = coords[3*m+0]; float B1 = coords[3*m+1]; float C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);
    //generate_central_grid(grid1,wt1,val1,0,Z1,1,0,1.,nrad,nang,ang_g,ang_w);
    //#pragma acc update self(grid1[0:6*gs],wt1[0:gs])
    //print_grid(gs,grid1,NULL,wt1,NULL,prl);

   //first compute single atom ints
    for (int i1=0;i1<Naux;i1++)
    if (basis_aux[i1][9]==m)
    {
      vector<double> basis1 = basis_aux[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];
      //printf("\n  m: %i i1: %2i   nlm: %i %i %2i zeta: %8.5f \n",m,i1,n1,l1,m1,zeta1);

      acc_assign(gs,val1,1.);

      eval_inr_r12(gs,grid1,val1,n1,l1,zeta1);
      eval_sh_3r(gs,grid1,val1,n1,l1,m1);

      //#pragma acc update self(val1[0:gs])
      //print_array(gs,val1);

      for (int i2=0;i2<N;i2++)
      if (basis[i2][9]==m)
      {
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];
        //printf("    n: %i i2: %2i   nlm: %i %i %2i zeta: %8.5f   (b2) \n",m,i2,n2,l2,m2,zeta2);

        acc_copyf(gs,valt1,val1);

        eval_sh(i2,gs,grid1,valt1,n2,l2,m2,zeta2);

        for (int i3=0;i3<N;i3++)
        if (basis[i3][9]==m)
        {
          vector<double> basis3 = basis[i3];
          int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

          acc_copyf(gs,valt2,valt1);

          eval_sh(i2,gs,grid1,valt2,n3,l3,m3,zeta3);

          float val = 0.;
         #pragma acc parallel loop independent present(valt2[0:gs],wt1[0:gs]) reduction(+:val)
          for (int j=0;j<gs;j++)
            val += valt2[j] * wt1[j];

         #pragma acc serial present(C[0:N2a])
          C[i1*N2+i2*N+i3] = val;

          //printf(" m: %i  i1/2/3: %2i %2i %2i val: %8.5f \n",m,i1,i2,i3,val);

        } //loop i3
      } //loop i2
    } //loop i1

   //two-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      float Z2 = (float)atno[n];
      float A2 = coords[3*n+0]; float B2 = coords[3*n+1]; float C2 = coords[3*n+2];
      float A12 = A2-A1; float B12 = B2-B1; float C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      //generate_central_grid(grid2,wt2,val2,0,Z2,1,0,1.,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2);
      recenter_grid(gs,grid2,A12,B12,C12);

      copy_grid(gs,grid1s,grid1);
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid2,wt2,Z1,Z2,A12,B12,C12);
      eliminate_small_wt(gs,wtt1);
      eliminate_small_wt(gs,wt2);

      add_r1_to_grid(gs,grid2,0.,0.,0.);

      //#pragma acc update self(grid1[0:6*gs],wt1[0:gs],grid2[0:6*gs],wt2[0:gs])
      //print_grid(gs,grid1,grid2,wt1,wt2,prl);

      for (int i1=0;i1<Naux;i1++)
      if (basis_aux[i1][9]==m)
      {
        vector<double> basis1 = basis_aux[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];
        //printf("  m: %i i1: %i   nlm: %i %i %2i zeta: %8.5f \n",m,i1,n1,l1,m1,zeta1);

        acc_assign(gs,val1,1.);
        acc_assign(gs,val2,1.);

        eval_inr_r12(gs,grid1,val1,n1,l1,zeta1);
        eval_sh_3r(gs,grid1,val1,n1,l1,m1);
        eval_inr_r12(gs,grid2,val2,n1,l1,zeta1);
        eval_sh_3r(gs,grid2,val2,n1,m1,l1);

        for (int i2=0;i2<N;i2++)
        if (basis[i2][9]==n)
        {
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];
          //printf("   n: %i i2: %i   nlm: %i %i %2i zeta: %8.5f \n",n,i2,n2,l2,m2,zeta2);

          acc_copyf(gs,valt1,val1);
          acc_copyf(gs,valt2,val2);

          eval_sh(i2,gs,grid2s,valt2,n2,l2,m2,zeta2);
          eval_sh(i2,gs,grid1s,valt1,n2,l2,m2,zeta2);

         //third AO function on same center as aux function
          for (int i3=0;i3<N;i3++)
          if (basis[i3][9]==m)
          {
            vector<double> basis3 = basis[i3];
            int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

            acc_copyf(gs,valt3,valt1);
            acc_copyf(gs,valt4,valt2);

            eval_sh(i3,gs,grid2,valt4,n3,l3,m3,zeta3);
            eval_sh(i3,gs,grid1,valt3,n3,l3,m3,zeta3);

            float val = 0.;
           #pragma acc parallel loop independent present(valt3[0:gs],wtt1[0:gs]) reduction(+:val)
            for (int j=0;j<gs;j++)
              val += valt3[j] * wtt1[j];
           #pragma acc parallel loop independent present(valt4[0:gs],wt2[0:gs]) reduction(+:val)
            for (int j=0;j<gs;j++)
              val += valt4[j] * wt2[j];

           #pragma acc serial present(C[0:N2a])
            C[i1*N2+i2*N+i3] = val;

            //printf(" 1. i1/2/3: %i %i %i / %i %i %i val: %5.3f \n",i1,i2,i3,m,n,m,val);
          } //loop i3 over third basis function

         //third AO function on same center as second AO function
          for (int i3=0;i3<N;i3++)
          if (basis[i3][9]==n)
          {
            vector<double> basis3 = basis[i3];
            int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

            acc_copyf(gs,valt3,valt1);
            acc_copyf(gs,valt4,valt2);

           //CPMZ check these
            eval_sh(i3,gs,grid2s,valt4,n3,l3,m3,zeta3);
            eval_sh(i3,gs,grid1s,valt3,n3,l3,m3,zeta3);

            float val = 0.;
           #pragma acc parallel loop independent present(valt3[0:gs],wtt1[0:gs]) reduction(+:val)
            for (int j=0;j<gs;j++)
              val += valt3[j] * wtt1[j];
           #pragma acc parallel loop independent present(valt4[0:gs],wt2[0:gs]) reduction(+:val)
            for (int j=0;j<gs;j++)
              val += valt4[j] * wt2[j];

           #pragma acc serial present(C[0:N2a])
            C[i1*N2+i2*N+i3] = val;

            //printf(" 2. i1/2/3: %i %i %i / %i %i %i val: %5.3f \n",i1,i2,i3,m,n,n,val);

          } //loop i3 over third basis function

        } //loop i2 over second basis function
      } //loop i1 over first basis function

    } //loop n over second atom

   //three-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      float Z2 = (float)atno[n];
      float A2 = coords[3*n+0]; float B2 = coords[3*n+1]; float C2 = coords[3*n+2];
      float A12 = A2-A1; float B12 = B2-B1; float C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      //generate_central_grid(grid2,wt2,val2,0,Z2,1,0,1.,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1);
      recenter_grid(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      for (int p=0;p<natoms;p++)
      if (p!=m && p!=n)
      {
        float Z3 = (float)atno[p];
        float A3 = coords[3*p+0]; float B3 = coords[3*p+1]; float C3 = coords[3*p+2];
        float A13 = A3-A1; float B13 = B3-B1; float C13 = C3-C1;
        //float A23 = A3-A2; float B23 = B3-B2; float C23 = C3-C2;

        generate_central_grid_2(grid3,wt3,Z3,nrad,nang,ang_g,ang_w);
        //generate_central_grid(grid3,wt3,val3,0,Z3,1,0,1.,nrad,nang,ang_g,ang_w);
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

        acc_copyf(gs,wtt1,wt1);
        acc_copyf(gs,wtt2,wt2);

        becke_weight_3c(gs,grid1,wtt1,grid2,wtt2,grid3,wt3,Z1,Z2,Z3,A12,B12,C12,A13,B13,C13);
        eliminate_small_wt_3(gs,wtt1,wtt2,wt3);

        add_r1_to_grid(gs,grid2,0.,0.,0.);
        add_r1_to_grid(gs,grid3,0.,0.,0.);

        //#pragma acc update self(grid1[0:6*gs],wtt1[0:gs],grid2[0:6*gs],wtt2[0:gs],grid3[0:6*gs],wt3[0:gs])
        //print_grid(gs,grid1,grid2,wtt1,wtt2,prl);
        //print_grid(gs,grid3,NULL,wt3,NULL,prl);

        for (int i1=0;i1<Naux;i1++)
        if (basis_aux[i1][9]==m)
        {
          vector<double> basis1 = basis_aux[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];
          //printf("  m: %i i1: %i   nlm: %i %i %2i zeta: %8.5f \n",m,i1,n1,l1,m1,zeta1);

          acc_assign(gs,val1,1.);
          acc_assign(gs,val2,1.);
          acc_assign(gs,val3,1.);

          eval_inr_r12(gs,grid1,val1,n1,l1,zeta1);
          eval_sh_3r(gs,grid1,val1,n1,l1,m1);
          eval_inr_r12(gs,grid2,val2,n1,l1,zeta1);
          eval_sh_3r(gs,grid2,val2,n1,l1,m1);
          eval_inr_r12(gs,grid3,val3,n1,l1,zeta1);
          eval_sh_3r(gs,grid3,val3,n1,l1,m1);

          for (int i2=0;i2<N;i2++)
          if (basis[i2][9]==n)
          {
            vector<double> basis2 = basis[i2];
            int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];
            //printf("   n: %i i2: %i   nlm: %i %i %2i zeta: %8.5f \n",n,i2,n2,l2,m2,zeta2);

            acc_copyf(gs,valt1,val1);
            acc_copyf(gs,valt2,val2);
            acc_copyf(gs,valt3,val3);

            eval_sh(i2,gs,grid3s,valt3,n2,l2,m2,zeta2);
            eval_sh(i2,gs,grid2s,valt2,n2,l2,m2,zeta2);
            eval_sh(i2,gs,grid1s,valt1,n2,l2,m2,zeta2);

            for (int i3=0;i3<N;i3++)
            if (basis[i3][9]==p)
            {
              vector<double> basis3 = basis[i3];
              int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

              acc_copyf(gs,valt4,valt1);
              acc_copyf(gs,valt5,valt2);
              acc_copyf(gs,valt6,valt3);

              eval_sh(i3,gs,grid3p,valt6,n3,l3,m3,zeta3);
              eval_sh(i3,gs,grid2p,valt5,n3,l3,m3,zeta3);
              eval_sh(i3,gs,grid1p,valt4,n3,l3,m3,zeta3);

              float val = 0.;
             #pragma acc parallel loop independent present(valt4[0:gs],wtt1[0:gs]) reduction(+:val)
              for (int j=0;j<gs;j++)
                val += valt4[j] * wtt1[j];
             #pragma acc parallel loop independent present(valt5[0:gs],wtt2[0:gs]) reduction(+:val)
              for (int j=0;j<gs;j++)
                val += valt5[j] * wtt2[j];
             #pragma acc parallel loop independent present(valt6[0:gs],wt3[0:gs]) reduction(+:val)
              for (int j=0;j<gs;j++)
                val += valt6[j] * wt3[j];

             #pragma acc serial present(C[0:N2a])
              C[i1*N2+i2*N+i3] = val;

              //printf(" 3c i1/2/3: %i %i %i / %i %i %i val: %5.3f \n",i1,i2,i3,m,n,m,val);
            } //loop i3
          } //loop i2
        } //loop i1

      } //loop p over third atom
    } //loop n over second atom

  } //loop m over natoms

 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc exit data copyout(C[0:N2a])
 #endif

#if DEBUG
  if (prl>0)
  {
    printf("\n C(raw): \n");
    for (int j=0;j<N;j++)
    for (int k=0;k<N;k++)
    {
      for (int i=0;i<Naux;i++)
        printf(" %8.5f",C[i*N2+j*N+k]);
      printf("\n");
    }
  }
#endif

  copy_symm(natoms,N,Naux,basis,basis_aux,C,1);

#if 0
  if (prl>1)
  {
    printf("\n C(raw): \n");
    for (int j=0;j<N;j++)
    for (int k=0;k<N;k++)
    {
      for (int i=0;i<Naux;i++)
        printf("  %12.6f",C[i*N2+j*N+k]);
      printf("\n");
    }
  }
#endif

  for (int i=0;i<Naux;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    C[i*N2+j*N+k] *= norm_sv(basis_aux[i][0],basis_aux[i][1],basis_aux[i][2],basis_aux[i][3])*basis[j][4]*basis[k][4];

  if (prl>2)
  {
    printf("\n C: \n");
    for (int i=0;i<Naux;i++)
    {
      printf(" i: %i \n",i);
      for (int j=0;j<N;j++)
      for (int k=0;k<N;k++)
        printf("  %12.6f",C[i*N2+j*N+k]);
      printf("\n");
    }
  }

 //CPMZ check this
#if USE_ACC
  #pragma acc exit data delete(grid1[0:6*gs],wt1[0:gs],val1[0:gs])
  #pragma acc exit data delete(grid2[0:6*gs],wt2[0:gs],val2[0:gs])
  #pragma acc exit data delete(grid3[0:6*gs],wt3[0:gs],val3[0:gs])

  #pragma acc exit data delete(grid1s[0:6*gs],grid2s[0:6*gs],grid3s[0:6*gs])
  #pragma acc exit data delete(grid1p[0:6*gs],grid2p[0:6*gs],grid3p[0:6*gs])

  #pragma acc exit data delete(valt1[0:gs],valt2[0:gs],valt3[0:gs],valt4[0:gs],valt5[0:gs],valt6[0:gs])
  #pragma acc exit data delete(wtt1[0:gs],wtt2[0:gs])
  #pragma acc exit data delete(C[0:N2a])
  #pragma acc exit data delete(n2i[0:natoms])
#endif

  delete [] n2i;

  delete [] ang_g;
  delete [] ang_w;

  delete [] grid1s;
  delete [] grid2s;
  delete [] grid3s;
  delete [] grid1p;
  delete [] grid2p;
  delete [] grid3p;

  delete [] valt1;
  delete [] valt2;
  delete [] valt3;
  delete [] wtt1;
  delete [] wtt2;
  delete [] grid1;
  delete [] grid2;
  delete [] grid3;
  delete [] wt1;
  delete [] wt2;
  delete [] wt3;
  delete [] val1;
  delete [] val2;
  delete [] val3;

  return;
}

void norm_2c(int N, vector<vector<double> > &basis, float* S)
{
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    S[i*N+j] *= basis[i][4]*basis[j][4];

  for (int i=0;i<N;i++)
  for (int j=i+1;j<N;j++)
    S[j*N+i] = S[i*N+j];

  return;
}

void norm_2c(int N, vector<vector<double> > &basis, double* S)
{
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    S[i*N+j] *= basis[i][4]*basis[j][4];

  for (int i=0;i<N;i++)
  for (int j=i+1;j<N;j++)
    S[j*N+i] = S[i*N+j];

  return;
}

#if RED_DOUBLE
void compute_VdV(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, int nc, float* coordsc, double* Pao, double* V, double* dV, int prl)
#else
void compute_VdV(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, int nc, float* coordsc, double* Pao, float* V, float* dV, int prl)
#endif
{
  if (nc<1) return;

  if (prl>1) printf(" beginning compute_VdV \n");

  int N = basis.size();
  int N2 = N*N;
  int nc3 = 3*nc;

  int gs = nrad*nang;
  int gs3 = 3*gs; int gs6 = 6*gs;

  int estart = find_center_of_grid(1,nrad)*nang;

  float* grid1 = new float[gs6];
  float* wt1 = new float[gs];

  float* grid2 = new float[gs6];
  float* wt2 = new float[gs];

  float* grid3 = new float[gs6];
  float* wt3 = new float[gs];

  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);
  //printf("  iN: %i \n",imaxN);

  double* norms = new double[N2];
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    norms[i*N+j] = basis[i][4]*basis[j][4];

 //intermediate storage
  int iN = imaxN;
  float* grid1s = new float[gs6]; float* grid2s = new float[gs6]; float* grid3s = new float[gs6];
  float* grid1p = new float[gs6]; float* grid2p = new float[gs6]; float* grid3p = new float[gs6];
  float** valS1 = new float*[iN]; for (int i=0;i<iN;i++) valS1[i] = new float[gs];
  float** valS2 = new float*[iN]; for (int i=0;i<iN;i++) valS2[i] = new float[gs];
  float** valS3 = new float*[iN]; for (int i=0;i<iN;i++) valS3[i] = new float[gs];
  float** valS4 = new float*[iN]; for (int i=0;i<iN;i++) valS4[i] = new float[gs];
  float** valS5 = new float*[iN]; for (int i=0;i<iN;i++) valS5[i] = new float[gs];
  float** valS6 = new float*[iN]; for (int i=0;i<iN;i++) valS6[i] = new float[gs];

  float** valt1 = new float*[iN]; for (int i=0;i<iN;i++) valt1[i] = new float[gs];
  float** valt2 = new float*[iN]; for (int i=0;i<iN;i++) valt2[i] = new float[gs];
  float** valt3 = new float*[iN]; for (int i=0;i<iN;i++) valt3[i] = new float[gs];

  float** valtv1 = new float*[iN]; for (int i=0;i<iN;i++) valtv1[i] = new float[gs3];
  float** valtv2 = new float*[iN]; for (int i=0;i<iN;i++) valtv2[i] = new float[gs3];
  float** valtv3 = new float*[iN]; for (int i=0;i<iN;i++) valtv3[i] = new float[gs3];

  float* wtt1 = new float[gs];
  float* wtt2 = new float[gs];

 #if RED_DOUBLE
  double* V1 = new double[nc];
  double* dV1 = new double[nc3];
 #else
  float* V1 = new float[nc];
  float* dV1 = new float[nc3];
 #endif

  float* ang_g = new float[3*nang];
  float* ang_w = new float[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

 #if USE_ACC
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data copyin(n2i[0:natoms],norms[0:N2])
  #pragma acc enter data copyin(coords[0:3*natoms],atno[0:natoms])
  #pragma acc enter data copyin(Pao[0:N2])

  #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
  #pragma acc enter data create(grid2[0:gs6],wt2[0:gs])
  #pragma acc enter data create(grid3[0:gs6],wt3[0:gs])
  #pragma acc enter data create(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs],valS5[0:iN][0:gs],valS6[0:iN][0:gs])
  #pragma acc enter data create(valt1[0:iN][0:gs],valt2[0:iN][0:gs],valt3[0:iN][0:gs])
  #pragma acc enter data create(valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3],valtv3[0:iN][0:gs3])
  #pragma acc enter data create(grid1s[0:gs6],grid1p[0:gs6],grid2s[0:gs6],grid2p[0:gs6],grid3s[0:gs6],grid3p[0:gs6])
  #pragma acc enter data create(wtt1[0:gs],wtt2[0:gs])
  #pragma acc enter data create(V[0:nc],dV[0:nc3],V1[0:nc],dV1[0:nc3])
 #endif
  acc_assign(nc,V,0.);
  acc_assign(nc3,dV,0.);

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    float Z1 = (float)atno[m];
    float A1 = coords[3*m+0]; float B1 = coords[3*m+1]; float C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);

   #pragma acc parallel loop collapse(2) present(valS1[0:iN][0:gs],valS3[0:iN][0:gs])
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
      for (int j=0;j<gs;j++)
        valS1[ii1][j] = valS3[ii1][j] = 1.f;
    }

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

      eval_sh(ii1,gs,grid1,valS1[ii1],n1,l1,m1,zeta1);
    }

    for (int i2=s1;i2<s2;i2++)
    {
      int ii2 = i2-s1;
      vector<double> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

      eval_sh(ii2,gs,grid1,valS3[ii2],n2,l2,m2,zeta2);
    } //loop i2 evaluate

   #if USE_ACC
    #pragma acc wait
   #endif


   //2 basis on same center (m), 1 charge elsewhere (p)
    for (int p=0;p<nc;p++)
    {
      float Zb = 1.;
      float An = coordsc[3*p+0]; float Bn = coordsc[3*p+1]; float Cn = coordsc[3*p+2];
      float A1n = An-A1; float B1n = Bn-B1; float C1n = Cn-C1;

      add_r2_to_grid(gs,grid1,A1n,B1n,C1n);

      generate_central_grid_2(grid3,wt3,Zb,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid3p,grid3); //grid 3 centered on charge
      recenter_grid(gs,grid3,A1n,B1n,C1n); //grid 3 centered on atom 1

      copy_grid(gs,grid1p,grid1);
      recenter_grid(gs,grid1p,-A1n,-B1n,-C1n); //grid 1 centered on charge

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid3,wt3,Z1,Zb,A1n,B1n,C1n);

      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wt3);

      add_r1_to_grid(gs,grid3,0.,0.,0.);
      add_r2_to_grid(gs,grid3,A1n,B1n,C1n);

      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc parallel loop present(valt1[0:iN][0:gs],valS1[0:iN][0:gs],wtt1[0:gs])
        for (int j=0;j<gs;j++)
          valt1[ii1][j] = valS1[ii1][j]*wtt1[j];

       #pragma acc parallel loop present(valS2[0:iN][0:gs],wt3[0:gs])
        for (int j=0;j<gs;j++)
          valS2[ii1][j] = wt3[j];

     //CPMZ check
     //  #pragma acc parallel loop present(valS3[0:iN][0:gs])
     //   for (int j=0;j<gs;j++)
     //     valS3[ii1][j] = 1.f;

       #pragma acc parallel loop present(valS4[0:iN][0:gs])
        for (int j=0;j<gs;j++)
          valS4[ii1][j] = 1.f;
      }


      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_sh(ii1,gs,grid3,valS2[ii1],n1,l1,m1,zeta1);
      }

      for (int i2=s1;i2<s2;i2++)
      {
        int ii2 = i2-s1;
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        eval_sh(ii2,gs,grid3,valS4[ii2],n2,l2,m2,zeta2);
      } //loop i2 evaluate


      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc parallel loop collapse(2) present(valtv1[0:iN][0:gs3],valt1[0:iN][0:gs])
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valtv1[ii1][3*j+k] = valt1[ii1][j];

       #pragma acc parallel loop collapse(2) present(valtv2[0:iN][0:gs3],valS2[0:iN][0:gs])
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valtv2[ii1][3*j+k] = valS2[ii1][j];
      }

     #pragma acc parallel loop present(grid1[0:gs6],grid3[0:gs6],grid1p[0:gs6],grid3p[0:gs6],valt1[0:iN][0:gs],valS2[0:iN][0:gs],valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3])
      for (int j=0;j<gs;j++)
      {
        float x1 = grid1p[6*j+0]; float y1 = grid1p[6*j+1]; float z1 = grid1p[6*j+2];
        float x3 = grid3p[6*j+0]; float y3 = grid3p[6*j+1]; float z3 = grid3p[6*j+2];
        float Rn1 = grid1[6*j+4]+1.e-20f;
        float Rn3 = grid3[6*j+4]+1.e-20f;
        float r12 = Rn1*Rn1; float r32 = Rn3*Rn3;

        float ne1 = -1.f/Rn1; float ne3 = -1.f/Rn3;
        float dx1 = x1*ne1/r12; float dy1 = y1*ne1/r12; float dz1 = z1*ne1/r12;
        float dx3 = x3*ne3/r32; float dy3 = y3*ne3/r32; float dz3 = z3*ne3/r32;

        //ne1 = ne3 = 1.;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valt1[ii1][j] *= ne1;
       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valS2[ii1][j] *= ne3;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
        {
          valtv1[ii1][3*j+0] *= dx1;
          valtv1[ii1][3*j+1] *= dy1;
          valtv1[ii1][3*j+2] *= dz1;
        }
       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
        {
          valtv2[ii1][3*j+0] *= dx3;
          valtv2[ii1][3*j+1] *= dy3;
          valtv2[ii1][3*j+2] *= dz3;
        }
      }

     //collect 2-atom values
      acc_assign(nc,V1,0.);
      acc_assign(nc3,dV1,0.);

      reduce_2c2v (p,s1,s2,gs,norms,Pao,valt1,valS2,valS3,valS4,iN,N,nc,1.,V1);
      reduce_2c2vd(p,s1,s2,gs,norms,Pao,valtv1,valtv2,valS3,valS4,iN,N,nc,1.,dV1);

     #pragma acc parallel loop present(V[0:nc],V1[0:nc])
      for (int j=0;j<nc;j++)
        V[j] += V1[j];
     #pragma acc parallel loop present(dV[0:nc3],dV1[0:nc3])
      for (int j=0;j<nc3;j++)
        dV[j] += dV1[j];

    } //loop p over nuclear center


   //two-center basis
    if (0)
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      float Z2 = (float)atno[n];
      float A2 = coords[3*n+0]; float B2 = coords[3*n+1]; float C2 = coords[3*n+2];
      float A12 = A2-A1; float B12 = B2-B1; float C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1);
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

     #pragma acc parallel loop present(valS1[0:iN][0:gs],valS2[0:iN][0:gs])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS1[ii1][j] = valS2[ii1][j] = 1.f;
      }

     #pragma acc parallel loop present(valS3[0:iN][0:gs],valS4[0:iN][0:gs])
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS3[ii2][j] = valS4[ii2][j] = 1.f;
      }

      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_sh(ii1,gs,grid1,valS1[ii1],n1,l1,m1,zeta1);
        eval_sh(ii1,gs,grid2,valS2[ii1],n1,l1,m1,zeta1);
      }

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1s,valS3[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2s,valS4[ii2],n2,l2,m2,zeta2);
      }


     //2 basis on diff centers, 1 charge elsewhere
      for (int p=0;p<nc;p++)
      {
        float Zb = 1.;
        float An = coordsc[3*p+0]; float Bn = coordsc[3*p+1]; float Cn = coordsc[3*p+2];
        float A1n = An-A1; float B1n = Bn-B1; float C1n = Cn-C1;

        generate_central_grid_2(grid3,wt3,Zb,nrad,nang,ang_g,ang_w);
        copy_grid(gs,grid3p,grid3); //grid 3 centered on charge
        recenter_grid(gs,grid3,A1n,B1n,C1n); //grid 3 centered on atom 1

        copy_grid(gs,grid1p,grid1);
        recenter_grid(gs,grid1p,-A1n,-B1n,-C1n); //grid 1 centered on charge
        copy_grid(gs,grid2p,grid2);
        recenter_grid(gs,grid2p,-A1n,-B1n,-C1n); //grid 2 centered on charge

        copy_grid(gs,grid3s,grid3);
        recenter_grid(gs,grid3s,-A12,-B12,-C12); //grid 3 centered on atom 2

        //add_r1_to_grid_6z(gs,grid1s,grid2s,grid3s,grid1p,grid2p,grid3p);
        add_r1_to_grid(gs,grid3s,0.,0.,0.);

      //need to keep all of these distances in order
        //add_r3_to_grid(gs,grid1,A1n,B1n,C1n);
        add_r123_to_grid(gs,grid1,0.,0.,0.,A12,B12,C12,A1n,B1n,C1n);
        add_r123_to_grid(gs,grid2,A12,B12,C12,0.,0.,0.,A1n,B1n,C1n);
      	add_r123_to_grid(gs,grid3,A1n,B1n,C1n,A12,B12,C12,0.,0.,0.);

        acc_copyf(gs,wtt1,wt1,wtt2,wt2);

        becke_weight_3c(gs,grid1,wtt1,grid2,wtt2,grid3,wt3,Z1,Z2,Zb,A12,B12,C12,A1n,B1n,C1n);
        eliminate_small_wt_3(estart,gs,wtt1,wtt2,wt3);

        add_r2_to_grid(gs,grid1,A1n,B1n,C1n);
        add_r1_to_grid(gs,grid2,0.,0.,0.);
        add_r1_to_grid(gs,grid3,0.,0.,0.);
        add_r2_to_grid(gs,grid2,A1n,B1n,C1n);
        add_r2_to_grid(gs,grid3,A1n,B1n,C1n);

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

         #pragma acc parallel loop collapse(2) present(valtv1[0:iN][0:gs3],valS1[0:iN][0:gs],wtt1[0:gs])
          for (int j=0;j<gs;j++)
          for (int k=0;k<3;k++)
            valtv1[ii1][3*j+k] = valS1[ii1][j]*wtt1[j];

         #pragma acc parallel loop collapse(2) present(valtv2[0:iN][0:gs3],valS2[0:iN][0:gs],wtt2[0:gs])
          for (int j=0;j<gs;j++)
          for (int k=0;k<3;k++)
            valtv2[ii1][3*j+k] = valS2[ii1][j]*wtt2[j];
        }

       #pragma acc parallel loop collapse(2) present(valS6[0:iN][0:gs])
        for (int ii2=0;ii2<s4-s3;ii2++)
        {
          for (int j=0;j<gs;j++)
            valS6[ii2][j] = 1.f;
        }

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<double> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

          eval_sh(ii1,gs,grid3,valS5[ii1],n1,l1,m1,zeta1);
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

          eval_sh(ii2,gs,grid3s,valS6[ii2],n2,l2,m2,zeta2);
        } //loop i2 evaluate

       #pragma acc parallel loop collapse(2) present(valtv3[0:iN][0:gs3],valS5[0:iN][0:gs])
        for (int ii1=0;ii1<s2-s1;ii1++)
        {
          for (int j=0;j<gs;j++)
          {
            float v1wt = valS5[ii1][j];
            valtv3[ii1][3*j+0] = v1wt;
            valtv3[ii1][3*j+1] = v1wt;
            valtv3[ii1][3*j+2] = v1wt;
          }
        }


       #pragma acc parallel loop present(grid1[0:gs6],grid2[0:gs6],grid3[0:gs6],grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6],valt1[0:iN][0:gs],valt2[0:iN][0:gs],valS5[0:iN][0:gs],valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3],valtv3[0:iN][0:gs3])
        for (int j=0;j<gs;j++)
        {
          float x1 = grid1p[6*j+0]; float y1 = grid1p[6*j+1]; float z1 = grid1p[6*j+2];
          float x2 = grid2p[6*j+0]; float y2 = grid2p[6*j+1]; float z2 = grid2p[6*j+2];
          float x3 = grid3p[6*j+0]; float y3 = grid3p[6*j+1]; float z3 = grid3p[6*j+2];
          float Rn1 = grid1[6*j+4]+1.e-20f;
          float Rn2 = grid2[6*j+4]+1.e-20f;
          float Rn3 = grid3[6*j+4]+1.e-20f;
          float r12 = Rn1*Rn1; float r22 = Rn2*Rn2; float r32 = Rn3*Rn3;
          float ne1 = -1.f/Rn1; float ne2 = -1.f/Rn2; float ne3 = -1.f/Rn3;
          float dx1 = x1*ne1/r12; float dy1 = y1*ne1/r12; float dz1 = z1*ne1/r12;
          float dx2 = x2*ne2/r22; float dy2 = y2*ne2/r22; float dz2 = z2*ne2/r22;
          float dx3 = x3*ne3/r32; float dy3 = y3*ne3/r32; float dz3 = z3*ne3/r32;

          ne1 = ne2 = ne3 = 1.;

         #pragma acc loop
          for (int ii1=0;ii1<s2-s1;ii1++)
            valt1[ii1][j] *= ne1;
         #pragma acc loop
          for (int ii1=0;ii1<s2-s1;ii1++)
            valt2[ii1][j] *= ne2;
         #pragma acc loop
          for (int ii1=0;ii1<s2-s1;ii1++)
            valS5[ii1][j] *= ne3;

         #pragma acc loop
          for (int ii1=0;ii1<s2-s1;ii1++)
          {
            valtv1[ii1][3*j+0] *= dx1;
            valtv1[ii1][3*j+1] *= dy1;
            valtv1[ii1][3*j+2] *= dz1;
          }
         #pragma acc loop
          for (int ii1=0;ii1<s2-s1;ii1++)
          {
            valtv2[ii1][3*j+0] *= dx2;
            valtv2[ii1][3*j+1] *= dy2;
            valtv2[ii1][3*j+2] *= dz2;
          }
         #pragma acc loop
          for (int ii1=0;ii1<s2-s1;ii1++)
          {
            valtv3[ii1][3*j+0] *= dx3;
            valtv3[ii1][3*j+1] *= dy3;
            valtv3[ii1][3*j+2] *= dz3;
          }
        }

        //collect 3-atom values
        acc_assign(nc,V1,0.);
        acc_assign(nc3,dV1,0.);

        reduce_2c3v (p,s1,s2,s3,s4,gs,norms,Pao,valt1, valt2, valS5, valS3,valS4,valS6,iN,N,nc,1.,V1);
        reduce_2c3vd(p,s1,s2,s3,s4,gs,norms,Pao,valtv1,valtv2,valtv3,valS3,valS4,valS6,iN,N,nc,1.,dV1);

       #pragma acc parallel loop present(V[0:nc],V1[0:nc])
        for (int j=0;j<nc;j++)
          V[j] += V1[j];
       #pragma acc parallel loop present(dV[0:nc3],dV1[0:nc3])
        for (int j=0;j<nc3;j++)
          dV[j] += dV1[j];

      } //loop p over nuclear center

    } //loop n over second atom


  } //loop m over natoms

 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc exit data delete(V1[0:nc],dV1[0:nc3])
  #pragma acc exit data copyout(V[0:nc],dV[0:nc3])
 #endif

  if (prl>1)
  {
    printf("\n V(no nn): \n");
    for (int i=0;i<nc;i++)
      printf(" %10.5f",V[i]);
    printf("\n");
  }

  if (prl>1)
  {
    printf("\n dV(no nn): \n");
    for (int i=0;i<nc;i++)
    {
      for (int j=0;j<3;j++)
        printf(" %12.6f",dV[3*i+j]);
      printf("\n");
    }
    printf("\n");
  }

 //include nuclear contribution to V
  if (1)
  for (int i=0;i<natoms;i++)
  {
    float Z1 = atno[i];
    float A1 = coords[3*i+0]; float B1 = coords[3*i+1]; float C1 = coords[3*i+2];
    for (int p=0;p<nc;p++)
    {
      float A1n = A1 - coordsc[3*p+0];
      float B1n = B1 - coordsc[3*p+1];
      float C1n = C1 - coordsc[3*p+2];
      double r1n = sqrt(A1n*A1n+B1n*B1n+C1n*C1n);
      double r2 = r1n*r1n;

      float ne1 = Z1/r1n;
      float dx = A1n*ne1/r2; float dy = B1n*ne1/r2; float dz = C1n*ne1/r2;

      V[p] += ne1;
      dV[3*p+0] += dx;
      dV[3*p+1] += dy;
      dV[3*p+2] += dz;
    }
  }

  if (prl>1)
  {
    printf("\n V: \n");
    for (int i=0;i<nc;i++)
      printf(" %10.5f",V[i]);
    printf("\n");
  }

  if (prl>1)
  {
    printf("\n dV: \n");
    for (int i=0;i<nc;i++)
    {
      for (int j=0;j<3;j++)
        printf(" %12.6f",dV[3*i+j]);
      printf("\n");
    }
    printf("\n");
  }

 //CPMZ check this
#if USE_ACC
  #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
  #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs])
  #pragma acc exit data delete(grid3[0:gs6],wt3[0:gs])
  #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6],grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6])
  #pragma acc exit data delete(wtt1[0:gs],wtt2[0:gs])
  #pragma acc exit data delete(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs],valS5[0:iN][0:gs],valS6[0:iN][0:gs])
  #pragma acc exit data delete(valt1[0:iN][0:gs],valt2[0:iN][0:gs],valt3[0:iN][0:gs],valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3],valtv3[0:iN][0:gs3])
  #pragma acc exit data delete(n2i[0:natoms])
  #pragma acc exit data delete(coords[0:3*natoms],atno[0:natoms])
#endif

  delete [] ang_g;
  delete [] ang_w;

  delete [] V1;
  delete [] dV1;

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

  for (int i=0;i<iN;i++) delete [] valt1[i];
  for (int i=0;i<iN;i++) delete [] valt2[i];
  for (int i=0;i<iN;i++) delete [] valt3[i];
  for (int i=0;i<iN;i++) delete [] valtv1[i];
  for (int i=0;i<iN;i++) delete [] valtv2[i];
  for (int i=0;i<iN;i++) delete [] valtv3[i];

  delete [] valS1; delete [] valS2; delete [] valS3; delete [] valS4; delete [] valS5; delete [] valS6;
  delete [] valt1; delete [] valt2; delete [] valt3; delete [] valtv1; delete [] valtv2; delete [] valtv3;
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
void compute_Enp_para(int ngpu, int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* En, double* pVp, int prl)
#else
void compute_Enp_para(int ngpu, int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, float* En, float* pVp, int prl)
#endif
{
 #if !USE_ACC
  printf("\n ERROR: compute_Enp_para requires OpenACC \n");
  exit(1);
 #endif

  int nomp = ngpu;
 //#pragma omp parallel
  //nomp = omp_get_num_threads();

  if (prl>1) printf(" beginning compute_Enp_para (%i) \n",nomp);

  int N = basis.size();
  int N2 = N*N;

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

 //compute_Enp now allows non-integer Zeff values
 // elements must still be identified
  float atz[natoms];
  for (int n=0;n<natoms;n++)
  {
    for (int m=0;m<N;m++)
    if (close_val(coords[3*n+0],basis[m][5]) && close_val(coords[3*n+1],basis[m][6]) && close_val(coords[3*n+2],basis[m][7]))
    {
      atz[n] = basis[m][8];
      break;
    }
  }
  if (prl>2)
  {
    printf("  atz: ");
    for (int n=0;n<natoms;n++)
      printf(" %5.3f",atz[n]);
    printf("\n");
  }

  int gs = nrad*nang;
  int gs3 = 3*gs; int gs6 = 6*gs;

  int estart = find_center_of_grid(1,nrad)*nang;

  float* grid1 = new float[gs6];
  float* wt1 = new float[gs];

  float* grid2 = new float[gs6];
  float* wt2 = new float[gs];

  float* grid3 = new float[gs6];
  float* wt3 = new float[gs];

  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);
  //printf("  iN: %i \n",imaxN);

 //intermediate storage
  int iN = imaxN;
  float* grid1s = new float[gs6]; float* grid2s = new float[gs6]; float* grid3s = new float[gs6];
  float* grid1p = new float[gs6]; float* grid2p = new float[gs6]; float* grid3p = new float[gs6];
  float** valS1 = new float*[iN]; for (int i=0;i<iN;i++) valS1[i] = new float[gs];
  float** valS2 = new float*[iN]; for (int i=0;i<iN;i++) valS2[i] = new float[gs];
  float** valS3 = new float*[iN]; for (int i=0;i<iN;i++) valS3[i] = new float[gs];
  float** valS4 = new float*[iN]; for (int i=0;i<iN;i++) valS4[i] = new float[gs];
  float** valS5 = new float*[iN]; for (int i=0;i<iN;i++) valS5[i] = new float[gs];
  float** valS6 = new float*[iN]; for (int i=0;i<iN;i++) valS6[i] = new float[gs];
  float** valV1 = new float*[iN]; for (int i=0;i<iN;i++) valV1[i] = new float[gs3];
  float** valV2 = new float*[iN]; for (int i=0;i<iN;i++) valV2[i] = new float[gs3];
  float** valV3 = new float*[iN]; for (int i=0;i<iN;i++) valV3[i] = new float[gs3];
  float** valV4 = new float*[iN]; for (int i=0;i<iN;i++) valV4[i] = new float[gs3];
  float** valV5 = new float*[iN]; for (int i=0;i<iN;i++) valV5[i] = new float[gs3];
  float** valV6 = new float*[iN]; for (int i=0;i<iN;i++) valV6[i] = new float[gs3];

  float** valt1 = new float*[iN]; for (int i=0;i<iN;i++) valt1[i] = new float[gs];
  float** valt2 = new float*[iN]; for (int i=0;i<iN;i++) valt2[i] = new float[gs];
  float** valt3 = new float*[iN]; for (int i=0;i<iN;i++) valt3[i] = new float[gs];
  float** valtv1 = new float*[iN]; for (int i=0;i<iN;i++) valtv1[i] = new float[gs3];
  float** valtv2 = new float*[iN]; for (int i=0;i<iN;i++) valtv2[i] = new float[gs3];
  float** valtv3 = new float*[iN]; for (int i=0;i<iN;i++) valtv3[i] = new float[gs3];

  float* wtt1 = new float[gs];
  float* wtt2 = new float[gs];

 #if RED_DOUBLE
  double* En1 = new double[N2];
  double* pVp1 = new double[N2];
 #else
  float* En1 = new float[N2];
  float* pVp1 = new float[N2];
 #endif

  float* ang_g = new float[3*nang];
  float* ang_w = new float[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

 #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
    #pragma acc enter data copyin(n2i[0:natoms])
    #pragma acc enter data copyin(coords[0:3*natoms],atno[0:natoms])

    #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
    #pragma acc enter data create(grid2[0:gs6],wt2[0:gs])
    #pragma acc enter data create(grid3[0:gs6],wt3[0:gs])
    #pragma acc enter data create(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs],valS5[0:iN][0:gs],valS6[0:iN][0:gs])
    #pragma acc enter data create(valV1[0:iN][0:gs3],valV2[0:iN][0:gs3],valV3[0:iN][0:gs3],valV4[0:iN][0:gs3],valV5[0:iN][0:gs3],valV6[0:iN][0:gs3])
    #pragma acc enter data create(valt1[0:iN][0:gs],valt2[0:iN][0:gs],valt3[0:iN][0:gs])
    #pragma acc enter data create(valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3],valtv3[0:iN][0:gs3])
    #pragma acc enter data create(grid1s[0:gs6],grid1p[0:gs6],grid2s[0:gs6],grid2p[0:gs6],grid3s[0:gs6],grid3p[0:gs6])
    #pragma acc enter data create(wtt1[0:gs],wtt2[0:gs])
    #pragma acc enter data create(En1[0:N2],pVp1[0:N2])
    if (tid>0)
    {
      #pragma acc enter data create(En[0:N2],pVp[0:N2])
    }
    acc_assign(N2,En,0.);
    acc_assign(N2,pVp,0.);
  }

 #pragma omp parallel for num_threads(nomp)
  for (int m=0;m<natoms;m++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);
    //printf("  launch %i/%i \n",m,tid);

   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    float Z1 = atz[m];
    int Z1g = atno[m];
    float A1 = coords[3*m+0]; float B1 = coords[3*m+1]; float C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1g,nrad,nang,ang_g,ang_w);

   #pragma acc parallel loop present(valS1[0:iN][0:gs],valV1[0:iN][0:gs3])
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valS1[ii1][j] = 1.f;
     #pragma acc loop
      for (int j=0;j<gs3;j++)
        valV1[ii1][j] = 1.f;
    }

   #pragma acc parallel loop present(valS3[0:iN][0:gs],valV3[0:iN][0:gs3])
    for (int ii2=0;ii2<s2-s1;ii2++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valS3[ii2][j] = 1.f;
     #pragma acc loop
      for (int j=0;j<gs;j++)
     #pragma acc loop
      for (int k=0;k<3;k++)
        valV3[ii2][3*j+k] = 1.f;
    }

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

      eval_sh(ii1,gs,grid1,valS1[ii1],n1,l1,m1,zeta1);
      eval_p(gs,grid1,valV1[ii1],n1,l1,m1,zeta1);
    }

    for (int i2=s1;i2<s2;i2++)
    {
      int ii2 = i2-s1;
      vector<double> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

      eval_sh(ii2,gs,grid1,valS3[ii2],n2,l2,m2,zeta2);
      eval_p(gs,grid1,valV3[ii2],n2,l2,m2,zeta2);
    } //loop i2 evaluate

   #if USE_ACC
    #pragma acc wait
   #endif


  //first, do single center
   #pragma acc parallel loop present(wt1[0:gs],valt1[0:iN][0:gs],valtv1[0:iN][0:gs3],valS1[0:iN][0:gs],valV1[0:iN][0:gs3])
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valt1[ii1][j] = valS1[ii1][j]*wt1[j];
     #pragma acc loop collapse(2)
      for (int j=0;j<gs;j++)
      for (int k=0;k<3;k++)
        valtv1[ii1][3*j+k] = valV1[ii1][3*j+k]*wt1[j];
    }

   #pragma acc parallel loop present(grid1[0:gs6],valt1[0:iN][0:gs],valtv1[0:iN][0:gs3])
    for (int j=0;j<gs;j++)
    {
      float Rn1 = grid1[6*j+3]+RPAD;
      float ne1 = Z1/Rn1;

     #pragma acc loop
      for (int ii1=0;ii1<s2-s1;ii1++)
        valt1[ii1][j] *= ne1;

     #pragma acc loop collapse(2)
      for (int ii1=0;ii1<s2-s1;ii1++)
      for (int k=0;k<3;k++)
        valtv1[ii1][3*j+k] *= ne1;
    }

   //collect 1-atom values
    if (Z1!=0.)
    {
      acc_assign(N2,En1,0.);
      acc_assign(N2,pVp1,0.);

      reduce_2c1(s1,s2,gs,valt1,valS3,iN,N,En1);
      reduce_2c1(s1,s2,gs3,valtv1,valV3,iN,N,pVp1);

     #pragma acc parallel loop present(En[0:N2],En1[0:N2],pVp[0:N2],pVp1[0:N2])
      for (int j=0;j<N2;j++)
      {
        En[j] += En1[j];
        pVp[j] += pVp1[j];
      }
    }

   //2 basis on same center, then 1 nucleus elsewhere
    for (int p=0;p<natoms;p++)
    if (p!=m)
    {
      float Zn = atz[p];
      int Zng = atno[p];
      float An = coords[3*p+0]; float Bn = coords[3*p+1]; float Cn = coords[3*p+2];
      float A1n = An-A1; float B1n = Bn-B1; float C1n = Cn-C1;

      add_r2_to_grid(gs,grid1,A1n,B1n,C1n);

      generate_central_grid_2(grid3,wt3,Zng,nrad,nang,ang_g,ang_w);
      //copy_grid(gs,grid3p,grid3); //grid 3 centered on atom 3
      recenter_grid(gs,grid3,A1n,B1n,C1n); //grid 3 centered on atom 1

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid3,wt3,Z1,Zn,A1n,B1n,C1n);

      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wt3);

      add_r1_to_grid(gs,grid3,0.,0.,0.);
      add_r2_to_grid(gs,grid3,A1n,B1n,C1n);

      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc parallel loop present(valt1[0:iN][0:gs],valS1[0:iN][0:gs],wtt1[0:gs])
        for (int j=0;j<gs;j++)
          valt1[ii1][j] = valS1[ii1][j]*wtt1[j];

       #pragma acc parallel loop collapse(2) present(valtv1[0:iN][0:gs3],valV1[0:iN][0:gs3],wtt1[0:gs])
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valtv1[ii1][3*j+k] = valV1[ii1][3*j+k]*wtt1[j];

       #pragma acc parallel loop present(valS2[0:iN][0:gs],wt3[0:gs])
        for (int j=0;j<gs;j++)
          valS2[ii1][j] = wt3[j];

       #pragma acc parallel loop collapse(2) present(valV2[0:iN][0:gs3],wt3[0:gs])
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valV2[ii1][3*j+k] = wt3[j];
      }


      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_sh(ii1,gs,grid3,valS2[ii1],n1,l1,m1,zeta1);
        eval_p(gs,grid3,valV2[ii1],n1,l1,m1,zeta1);
      }

      for (int i2=s1;i2<s2;i2++)
      {
        int ii2 = i2-s1;
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        acc_assign(gs,valS4[ii2],1.);
        acc_assign(gs3,valV4[ii2],1.);

        eval_sh(ii2,gs,grid3,valS4[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid3,valV4[ii2],n2,l2,m2,zeta2);
      } //loop i2 evaluate

     #pragma acc parallel loop present(grid1[0:gs6],grid3[0:gs6],valt1[0:iN][0:gs],valS2[0:iN][0:gs],valtv1[0:iN][0:gs3],valV2[0:iN][0:gs3])
      for (int j=0;j<gs;j++)
      {
        float Rn1 = grid1[6*j+4]+RPAD;
        float Rn3 = grid3[6*j+4]+RPAD;
        float ne1 = Zn/Rn1;
        float ne3 = Zn/Rn3;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valt1[ii1][j] *= ne1;
       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valS2[ii1][j] *= ne3;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
       #pragma acc loop
        for (int k=0;k<3;k++)
          valtv1[ii1][3*j+k] *= ne1;
       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
       #pragma acc loop
        for (int k=0;k<3;k++)
          valV2[ii1][3*j+k] *= ne3;
      }

     //collect 2-atom values
      if (Zn!=0.)
      {
        acc_assign(N2,En1,0.);
        acc_assign(N2,pVp1,0.);

        reduce_2c2(s1,s2,s1,s2,gs, valt1, valS2,valS3,valS4,iN,N,En1);
        reduce_2c2(s1,s2,s1,s2,gs3,valtv1,valV2,valV3,valV4,iN,N,pVp1);

       #pragma acc parallel loop present(En[0:N2],En1[0:N2],pVp[0:N2],pVp1[0:N2])
        for (int j=0;j<N2;j++)
        {
          En[j] += En1[j];
          pVp[j] += pVp1[j];
        }
      }

    } //loop p over nuclear center




   //two-center basis
    for (int n=m+1;n<natoms;n++)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      float Z2 = atz[n];
      int Z2g = atno[n];
      float A2 = coords[3*n+0]; float B2 = coords[3*n+1]; float C2 = coords[3*n+2];
      float A12 = A2-A1; float B12 = B2-B1; float C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2g,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1);
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      acc_copyf(gs,wtt1,wt1);
      acc_copyf(gs,wtt2,wt2);
      becke_weight_2c(gs,grid1,wtt1,grid2,wtt2,Z1,Z2,A12,B12,C12);

      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wtt2);

     //needs to happen after Becke weighting
      add_r1_to_grid(gs,grid2,0.,0.,0.);
      add_r2_to_grid(gs,grid2,A12,B12,C12);


     #pragma acc parallel loop present(valS1[0:iN][0:gs],valV1[0:iN][0:gs3],valS2[0:iN][0:gs],valV2[0:iN][0:gs3])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS1[ii1][j] = valS2[ii1][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          valV1[ii1][j] = valV2[ii1][j] = 1.f;
      }

     #pragma acc parallel loop present(valS3[0:iN][0:gs],valV3[0:iN][0:gs3],valS4[0:iN][0:gs],valV4[0:iN][0:gs3])
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS3[ii2][j] = valS4[ii2][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          valV3[ii2][j] = valV4[ii2][j] = 1.f;
      }

      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_sh(ii1,gs,grid1,valS1[ii1],n1,l1,m1,zeta1);
        eval_sh(ii1,gs,grid2,valS2[ii1],n1,l1,m1,zeta1);
        eval_p(gs,grid1,valV1[ii1],n1,l1,m1,zeta1);
        eval_p(gs,grid2,valV2[ii1],n1,l1,m1,zeta1);
      }

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1s,valS3[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2s,valS4[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid1s,valV3[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid2s,valV4[ii2],n2,l2,m2,zeta2);
      }

     //2 basis on diff centers, 1 nucleus matches
     #pragma acc parallel loop present(wtt1[0:gs],valt1[0:iN][0:gs],valtv1[0:iN][0:gs3],valS1[0:iN][0:gs],valV1[0:iN][0:gs3])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valt1[ii1][j] = valS1[ii1][j]*wtt1[j];

       #pragma acc loop collapse(2)
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valtv1[ii1][3*j+k] = valV1[ii1][3*j+k]*wtt1[j];
      }

     #pragma acc parallel loop present(wtt2[0:gs],valt2[0:iN][0:gs],valtv2[0:iN][0:gs3],valS2[0:iN][0:gs],valV2[0:iN][0:gs3])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valt2[ii1][j] = valS2[ii1][j]*wtt2[j];

       #pragma acc loop collapse(2)
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valtv2[ii1][3*j+k] = valV2[ii1][3*j+k]*wtt2[j];
      }

     #pragma acc parallel loop present(grid1[0:gs6],grid2[0:gs6],valt1[0:iN][0:gs],valt2[0:iN][0:gs],valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3])
      for (int j=0;j<gs;j++)
      {
        float Rn1a = grid1[6*j+3]+RPAD;
        float Rn2a = grid1[6*j+4]+RPAD;
        float Rn1b = grid2[6*j+3]+RPAD;
        float Rn2b = grid2[6*j+4]+RPAD;
        float ne1 = Z1/Rn1a+Z2/Rn2a;
        float ne2 = Z1/Rn1b+Z2/Rn2b;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valt1[ii1][j] *= ne1;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valt2[ii1][j] *= ne2;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
       #pragma acc loop
        for (int k=0;k<3;k++)
          valtv1[ii1][3*j+k] *= ne1;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
       #pragma acc loop
        for (int k=0;k<3;k++)
          valtv2[ii1][3*j+k] *= ne2;
      }

     //collect 2c-atom same-atom values
      if (Z1!=0. || Z2!=0.)
      {
        acc_assign(N2,En1,0.);
        acc_assign(N2,pVp1,0.);

        reduce_2c2(s1,s2,s3,s4,gs,valt1,valt2,valS3,valS4,iN,N,En1);
        reduce_2c2(s1,s2,s3,s4,gs3,valtv1,valtv2,valV3,valV4,iN,N,pVp1);

       #pragma acc parallel loop present(En[0:N2],En1[0:N2])
        for (int j=0;j<N2;j++)
          En[j] += En1[j];
       #pragma acc parallel loop present(pVp[0:N2],pVp1[0:N2])
        for (int j=0;j<N2;j++)
          pVp[j] += pVp1[j];
      }


     //2 basis on diff centers, then 1 nucleus elsewhere
      for (int p=0;p<natoms;p++)
      if (p!=m && p!=n)
      {
        float Zn = atz[p];
        int Zng = atno[p];
        float An = coords[3*p+0]; float Bn = coords[3*p+1]; float Cn = coords[3*p+2];
        float A1n = An-A1; float B1n = Bn-B1; float C1n = Cn-C1;

        generate_central_grid_2(grid3,wt3,Zng,nrad,nang,ang_g,ang_w);
        //copy_grid(gs,grid3p,grid3); //grid 3 centered on atom 3
        recenter_grid(gs,grid3,A1n,B1n,C1n); //grid 3 centered on atom 1

        copy_grid(gs,grid3s,grid3);
        recenter_grid(gs,grid3s,-A12,-B12,-C12); //grid 3 centered on atom 2

        //copy_grid(gs,grid2p,grid2); //grid 2 centered on atom 3
        //recenter_grid(gs,grid2p,-A1n,-B1n,-C1n);

        //copy_grid(gs,grid1p,grid1);
        //recenter_grid(gs,grid1p,-A1n,-B1n,-C1n); //grid 1 centered on atom 3

        //add_r1_to_grid_6z(gs,grid1s,grid2s,grid3s,grid1p,grid2p,grid3p);
        add_r1_to_grid(gs,grid3s,0.,0.,0.);

      //need to keep all of these distances in order
        //add_r3_to_grid(gs,grid1,A1n,B1n,C1n);
        add_r123_to_grid(gs,grid1,0.,0.,0.,A12,B12,C12,A1n,B1n,C1n);
        add_r123_to_grid(gs,grid2,A12,B12,C12,0.,0.,0.,A1n,B1n,C1n);
      	add_r123_to_grid(gs,grid3,A1n,B1n,C1n,A12,B12,C12,0.,0.,0.);

        acc_copyf(gs,wtt1,wt1,wtt2,wt2);

        becke_weight_3c(gs,grid1,wtt1,grid2,wtt2,grid3,wt3,Z1,Z2,Zn,A12,B12,C12,A1n,B1n,C1n);
        eliminate_small_wt_3(estart,gs,wtt1,wtt2,wt3);

        add_r2_to_grid(gs,grid1,A1n,B1n,C1n);
        add_r1_to_grid(gs,grid2,0.,0.,0.);
        add_r1_to_grid(gs,grid3,0.,0.,0.);
        add_r2_to_grid(gs,grid2,A1n,B1n,C1n);
        add_r2_to_grid(gs,grid3,A1n,B1n,C1n);

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

         #pragma acc parallel loop collapse(2) present(valtv1[0:iN][0:gs3],valV1[0:iN][0:gs3],wtt1[0:gs])
          for (int j=0;j<gs;j++)
          for (int k=0;k<3;k++)
            valtv1[ii1][3*j+k] = valV1[ii1][3*j+k]*wtt1[j];

         #pragma acc parallel loop collapse(2) present(valtv2[0:iN][0:gs3],valV2[0:iN][0:gs3],wtt2[0:gs])
          for (int j=0;j<gs;j++)
          for (int k=0;k<3;k++)
            valtv2[ii1][3*j+k] = valV2[ii1][3*j+k]*wtt2[j];

         #pragma acc parallel loop collapse(2) present(valV5[0:iN][0:gs3],wt3[0:gs])
          for (int j=0;j<gs;j++)
          for (int k=0;k<3;k++)
            valV5[ii1][3*j+k] = wt3[j];
        }

       #pragma acc parallel loop present(valS6[0:iN][0:gs],valV6[0:iN][0:gs3])
        for (int ii2=0;ii2<s4-s3;ii2++)
        {
         #pragma acc loop
          for (int j=0;j<gs;j++)
            valS6[ii2][j] = 1.f;

         #pragma acc loop
          for (int j=0;j<gs3;j++)
            valV6[ii2][j] = 1.f;
        }

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<double> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

          eval_sh(ii1,gs,grid3,valS5[ii1],n1,l1,m1,zeta1);
          eval_p(gs,grid3,valV5[ii1],n1,l1,m1,zeta1);
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

          eval_sh(ii2,gs,grid3s,valS6[ii2],n2,l2,m2,zeta2);
          eval_p(gs,grid3s,valV6[ii2],n2,l2,m2,zeta2);
        } //loop i2 evaluate

       #pragma acc parallel loop present(grid1[0:gs6],grid2[0:gs6],grid3[0:gs6],valt1[0:iN][0:gs],valt2[0:iN][0:gs],valS5[0:iN][0:gs],valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3],valV5[0:iN][0:gs3])
        for (int j=0;j<gs;j++)
        {
          float Rn1 = grid1[6*j+4]+RPAD; float Rn2 = grid2[6*j+4]+RPAD; float Rn3 = grid3[6*j+4]+RPAD;
          float ne1 = Zn/Rn1; float ne2 = Zn/Rn2; float ne3 = Zn/Rn3;

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
            valtv1[ii1][3*j+k] *= ne1;
            valtv2[ii1][3*j+k] *= ne2;
            valV5[ii1][3*j+k] *= ne3;
          }
        }

        //collect 3-atom values
        if (Zn!=0.)
        {
          acc_assign(N2,En1,0.);
          acc_assign(N2,pVp1,0.);

          reduce_2c3(s1,s2,s3,s4,gs, valt1, valt2, valS3,valS4,valS5,valS6,iN,N,En1);
          reduce_2c3(s1,s2,s3,s4,gs3,valtv1,valtv2,valV3,valV4,valV5,valV6,iN,N,pVp1);

         #pragma acc parallel loop present(En[0:N2],En1[0:N2],pVp[0:N2],pVp1[0:N2])
          for (int j=0;j<N2;j++)
          {
            En[j] += En1[j];
            pVp[j] += pVp1[j];
          }
        }

      } //loop p over nuclear center

    } //loop n over second atom


  } //loop m over natoms


 //gather from all gpus
  double* En_all = new double[N2]();
  double* pVp_all = new double[N2]();

  for (int n=0;n<nomp;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc update self(En[0:N2],pVp[0:N2])

    for (int i=0;i<N2;i++)
      En_all[i] += En[i];
    for (int i=0;i<N2;i++)
      pVp_all[i] += pVp[i];

    #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
    #pragma acc exit data delete(En1[0:N2],pVp1[0:N2])
  }
  acc_set_device_num(0,acc_device_nvidia);

  for (int i=0;i<N2;i++)
    En[i] = En_all[i];
  for (int i=0;i<N2;i++)
    pVp[i] = pVp_all[i];

  delete [] En_all;
  delete [] pVp_all;

  #pragma acc update device(En[0:N2],pVp[0:N2])

 //done gathering


  double* norm = new double[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];
  #pragma acc enter data copyin(norm[0:N])

 #pragma acc parallel loop independent present(En[0:N2],pVp[0:N2],norm[0:N])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i;j<N;j++)
  {
    double n12 = norm[i]*norm[j];
    En[i*N+j]  *= -n12;
    pVp[i*N+j] *= -n12;
  }

 #pragma acc parallel loop independent present(En[0:N2],pVp[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i+1;j<N;j++)
  {
    En[j*N+i] = En[i*N+j];
    pVp[j*N+i] = pVp[i*N+j];
  }

  #pragma acc exit data delete(norm[0:N])
  delete [] norm;

  clean_small_values(N,En);
  clean_small_values(N,pVp);


  if (prl>1)
  {
    printf("\n En: \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %10.5f",En[i*N+j]);
      printf("\n");
    }
  }

  if (prl>2)
  {
    printf("\n pVp: \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %12.6f",pVp[i*N+j]);
      printf("\n");
    }
    printf("\n");
  }

 //CPMZ check this
 #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);

    #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
    #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs])
    #pragma acc exit data delete(grid3[0:gs6],wt3[0:gs])
    #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6],grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6])
    #pragma acc exit data delete(wtt1[0:gs],wtt2[0:gs])
    #pragma acc exit data delete(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs],valS5[0:iN][0:gs],valS6[0:iN][0:gs])
    #pragma acc exit data delete(valV1[0:iN][0:gs3],valV2[0:iN][0:gs3],valV3[0:iN][0:gs3],valV4[0:iN][0:gs3],valV5[0:iN][0:gs3],valV6[0:iN][0:gs3])
    #pragma acc exit data delete(valt1[0:iN][0:gs],valt2[0:iN][0:gs],valt3[0:iN][0:gs],valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3],valtv3[0:iN][0:gs3])
    #pragma acc exit data delete(n2i[0:natoms])
    #pragma acc exit data delete(coords[0:3*natoms],atno[0:natoms])

    if (tid>0)
    {
      #pragma acc exit data delete(pVp[0:N2],En[0:N2])
    }
  }
  acc_set_device_num(0,acc_device_nvidia);

  delete [] ang_g;
  delete [] ang_w;

  delete [] pVp1;
  delete [] En1;

  delete [] n2i;

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
  for (int i=0;i<iN;i++) delete [] valV1[i];
  for (int i=0;i<iN;i++) delete [] valV2[i];
  for (int i=0;i<iN;i++) delete [] valV3[i];
  for (int i=0;i<iN;i++) delete [] valV4[i];
  for (int i=0;i<iN;i++) delete [] valV5[i];
  for (int i=0;i<iN;i++) delete [] valV6[i];

  for (int i=0;i<iN;i++) delete [] valt1[i];
  for (int i=0;i<iN;i++) delete [] valt2[i];
  for (int i=0;i<iN;i++) delete [] valt3[i];
  for (int i=0;i<iN;i++) delete [] valtv1[i];
  for (int i=0;i<iN;i++) delete [] valtv2[i];
  for (int i=0;i<iN;i++) delete [] valtv3[i];

  delete [] valS1; delete [] valS2; delete [] valS3; delete [] valS4; delete [] valS5; delete [] valS6;
  delete [] valV1; delete [] valV2; delete [] valV3; delete [] valV4; delete [] valV5; delete [] valV6;
  delete [] valt1; delete [] valt2; delete [] valt3; delete [] valtv1; delete [] valtv2; delete [] valtv3;
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
void compute_Enp(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* En, double* pVp, int prl)
#else
void compute_Enp(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, float* En, float* pVp, int prl)
#endif
{
  if (prl>1) printf(" beginning compute_Enp \n");

  int N = basis.size();
  int N2 = N*N;

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

  float atz[natoms];
  for (int n=0;n<natoms;n++)
  {
    for (int m=0;m<N;m++)
    if (close_val(coords[3*n+0],basis[m][5]) && close_val(coords[3*n+1],basis[m][6]) && close_val(coords[3*n+2],basis[m][7]))
    {
      atz[n] = basis[m][8];
      break;
    }
  }
  if (prl>2)
  {
    printf("  atz: ");
    for (int n=0;n<natoms;n++)
      printf(" %5.3f",atz[n]);
    printf("\n");
  }

  int gs = nrad*nang;
  int gs3 = 3*gs; int gs6 = 6*gs;

  int estart = find_center_of_grid(1,nrad)*nang;

  float* grid1 = new float[gs6];
  float* wt1 = new float[gs];

  float* grid2 = new float[gs6];
  float* wt2 = new float[gs];

  float* grid3 = new float[gs6];
  float* wt3 = new float[gs];

  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);
  //printf("  iN: %i \n",imaxN);

 //intermediate storage
  int iN = imaxN;
  float* grid1s = new float[gs6]; float* grid2s = new float[gs6]; float* grid3s = new float[gs6];
  float* grid1p = new float[gs6]; float* grid2p = new float[gs6]; float* grid3p = new float[gs6];
  float** valS1 = new float*[iN]; for (int i=0;i<iN;i++) valS1[i] = new float[gs];
  float** valS2 = new float*[iN]; for (int i=0;i<iN;i++) valS2[i] = new float[gs];
  float** valS3 = new float*[iN]; for (int i=0;i<iN;i++) valS3[i] = new float[gs];
  float** valS4 = new float*[iN]; for (int i=0;i<iN;i++) valS4[i] = new float[gs];
  float** valS5 = new float*[iN]; for (int i=0;i<iN;i++) valS5[i] = new float[gs];
  float** valS6 = new float*[iN]; for (int i=0;i<iN;i++) valS6[i] = new float[gs];
  float** valV1 = new float*[iN]; for (int i=0;i<iN;i++) valV1[i] = new float[gs3];
  float** valV2 = new float*[iN]; for (int i=0;i<iN;i++) valV2[i] = new float[gs3];
  float** valV3 = new float*[iN]; for (int i=0;i<iN;i++) valV3[i] = new float[gs3];
  float** valV4 = new float*[iN]; for (int i=0;i<iN;i++) valV4[i] = new float[gs3];
  float** valV5 = new float*[iN]; for (int i=0;i<iN;i++) valV5[i] = new float[gs3];
  float** valV6 = new float*[iN]; for (int i=0;i<iN;i++) valV6[i] = new float[gs3];

  float** valt1 = new float*[iN]; for (int i=0;i<iN;i++) valt1[i] = new float[gs];
  float** valt2 = new float*[iN]; for (int i=0;i<iN;i++) valt2[i] = new float[gs];
  float** valt3 = new float*[iN]; for (int i=0;i<iN;i++) valt3[i] = new float[gs];
  float** valtv1 = new float*[iN]; for (int i=0;i<iN;i++) valtv1[i] = new float[gs3];
  float** valtv2 = new float*[iN]; for (int i=0;i<iN;i++) valtv2[i] = new float[gs3];
  float** valtv3 = new float*[iN]; for (int i=0;i<iN;i++) valtv3[i] = new float[gs3];

  float* wtt1 = new float[gs];
  float* wtt2 = new float[gs];

 #if RED_DOUBLE
  double* En1 = new double[N2];
  double* pVp1 = new double[N2]; 
 #else
  float* En1 = new float[N2];
  float* pVp1 = new float[N2]; 
 #endif

  float* ang_g = new float[3*nang];
  float* ang_w = new float[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

 #if USE_ACC
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data copyin(n2i[0:natoms])
  #pragma acc enter data copyin(coords[0:3*natoms],atno[0:natoms])

  #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
  #pragma acc enter data create(grid2[0:gs6],wt2[0:gs])
  #pragma acc enter data create(grid3[0:gs6],wt3[0:gs])
  #pragma acc enter data create(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs],valS5[0:iN][0:gs],valS6[0:iN][0:gs])
  #pragma acc enter data create(valV1[0:iN][0:gs3],valV2[0:iN][0:gs3],valV3[0:iN][0:gs3],valV4[0:iN][0:gs3],valV5[0:iN][0:gs3],valV6[0:iN][0:gs3])
  #pragma acc enter data create(valt1[0:iN][0:gs],valt2[0:iN][0:gs],valt3[0:iN][0:gs])
  #pragma acc enter data create(valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3],valtv3[0:iN][0:gs3])
  #pragma acc enter data create(grid1s[0:gs6],grid1p[0:gs6],grid2s[0:gs6],grid2p[0:gs6],grid3s[0:gs6],grid3p[0:gs6])
  #pragma acc enter data create(wtt1[0:gs],wtt2[0:gs])
  #pragma acc enter data create(En1[0:N2],pVp1[0:N2])
  //#pragma acc enter data create(En[0:N2],pVp[0:N2])
 #endif
  acc_assign(N2,En,0.);
  acc_assign(N2,pVp,0.);

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    float Z1 = atz[m];
    int Z1g = atno[m];
    float A1 = coords[3*m+0]; float B1 = coords[3*m+1]; float C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1g,nrad,nang,ang_g,ang_w);

   #pragma acc parallel loop present(valS1[0:iN][0:gs],valV1[0:iN][0:gs3])
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valS1[ii1][j] = 1.f;
     #pragma acc loop
      for (int j=0;j<gs3;j++)
        valV1[ii1][j] = 1.f;
    }

   #pragma acc parallel loop present(valS3[0:iN][0:gs],valV3[0:iN][0:gs3])
    for (int ii2=0;ii2<s2-s1;ii2++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valS3[ii2][j] = 1.f;
     #pragma acc loop
      for (int j=0;j<gs;j++)
     #pragma acc loop
      for (int k=0;k<3;k++)
        valV3[ii2][3*j+k] = 1.f;
    }

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

      eval_sh(ii1,gs,grid1,valS1[ii1],n1,l1,m1,zeta1);
      eval_p(gs,grid1,valV1[ii1],n1,l1,m1,zeta1);
    }

    for (int i2=s1;i2<s2;i2++)
    {
      int ii2 = i2-s1;
      vector<double> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

      eval_sh(ii2,gs,grid1,valS3[ii2],n2,l2,m2,zeta2);
      eval_p(gs,grid1,valV3[ii2],n2,l2,m2,zeta2);
    } //loop i2 evaluate

   #if USE_ACC
    #pragma acc wait
   #endif


  //first, do single center
   #pragma acc parallel loop present(wt1[0:gs],valt1[0:iN][0:gs],valtv1[0:iN][0:gs3],valS1[0:iN][0:gs],valV1[0:iN][0:gs3])
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valt1[ii1][j] = valS1[ii1][j]*wt1[j];
     #pragma acc loop collapse(2)
      for (int j=0;j<gs;j++)
      for (int k=0;k<3;k++)
        valtv1[ii1][3*j+k] = valV1[ii1][3*j+k]*wt1[j];
    }

   #pragma acc parallel loop present(grid1[0:gs6],valt1[0:iN][0:gs],valtv1[0:iN][0:gs3])
    for (int j=0;j<gs;j++)
    {
      float Rn1 = grid1[6*j+3]+RPAD;
      float ne1 = Z1/Rn1;

     #pragma acc loop
      for (int ii1=0;ii1<s2-s1;ii1++)
        valt1[ii1][j] *= ne1;

     #pragma acc loop collapse(2)
      for (int ii1=0;ii1<s2-s1;ii1++)
      for (int k=0;k<3;k++)
        valtv1[ii1][3*j+k] *= ne1;
    }

   //collect 1-atom values
    if (Z1!=0.)
    {
      acc_assign(N2,En1,0.);
      acc_assign(N2,pVp1,0.);

      reduce_2c1(s1,s2,gs,valt1,valS3,iN,N,En1);
      reduce_2c1(s1,s2,gs3,valtv1,valV3,iN,N,pVp1);

     #pragma acc parallel loop present(En[0:N2],En1[0:N2],pVp[0:N2],pVp1[0:N2])
      for (int j=0;j<N2;j++)
      {
        En[j] += En1[j];
        pVp[j] += pVp1[j];
      }
    }

   //2 basis on same center, then 1 nucleus elsewhere
    for (int p=0;p<natoms;p++)
    if (p!=m)
    {
      float Zn = atno[p];
      int Zng = atno[p];
      float An = coords[3*p+0]; float Bn = coords[3*p+1]; float Cn = coords[3*p+2];
      float A1n = An-A1; float B1n = Bn-B1; float C1n = Cn-C1;

      add_r2_to_grid(gs,grid1,A1n,B1n,C1n);

      generate_central_grid_2(grid3,wt3,Zng,nrad,nang,ang_g,ang_w);
      //copy_grid(gs,grid3p,grid3); //grid 3 centered on atom 3
      recenter_grid(gs,grid3,A1n,B1n,C1n); //grid 3 centered on atom 1

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid3,wt3,Z1,Zn,A1n,B1n,C1n);

      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wt3);

      add_r1_to_grid(gs,grid3,0.,0.,0.);
      add_r2_to_grid(gs,grid3,A1n,B1n,C1n);

      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc parallel loop present(valt1[0:iN][0:gs],valS1[0:iN][0:gs],wtt1[0:gs])
        for (int j=0;j<gs;j++)
          valt1[ii1][j] = valS1[ii1][j]*wtt1[j];

       #pragma acc parallel loop collapse(2) present(valtv1[0:iN][0:gs3],valV1[0:iN][0:gs3],wtt1[0:gs])
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valtv1[ii1][3*j+k] = valV1[ii1][3*j+k]*wtt1[j];

       #pragma acc parallel loop present(valS2[0:iN][0:gs],wt3[0:gs])
        for (int j=0;j<gs;j++)
          valS2[ii1][j] = wt3[j];

       #pragma acc parallel loop collapse(2) present(valV2[0:iN][0:gs3],wt3[0:gs])
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valV2[ii1][3*j+k] = wt3[j];
      }


      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_sh(ii1,gs,grid3,valS2[ii1],n1,l1,m1,zeta1);
        eval_p(gs,grid3,valV2[ii1],n1,l1,m1,zeta1);
      }

      for (int i2=s1;i2<s2;i2++)
      {
        int ii2 = i2-s1;
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        acc_assign(gs,valS4[ii2],1.);
        acc_assign(gs3,valV4[ii2],1.);

        eval_sh(ii2,gs,grid3,valS4[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid3,valV4[ii2],n2,l2,m2,zeta2);
      } //loop i2 evaluate

     #pragma acc parallel loop present(grid1[0:gs6],grid3[0:gs6],valt1[0:iN][0:gs],valS2[0:iN][0:gs],valtv1[0:iN][0:gs3],valV2[0:iN][0:gs3])
      for (int j=0;j<gs;j++)
      {
        float Rn1 = grid1[6*j+4]+RPAD;
        float Rn3 = grid3[6*j+4]+RPAD;
        float ne1 = Zn/Rn1;
        float ne3 = Zn/Rn3;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valt1[ii1][j] *= ne1;
       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valS2[ii1][j] *= ne3;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
       #pragma acc loop
        for (int k=0;k<3;k++)
          valtv1[ii1][3*j+k] *= ne1;
       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
       #pragma acc loop
        for (int k=0;k<3;k++)
          valV2[ii1][3*j+k] *= ne3;
      }

     //collect 2-atom values
      if (Zn!=0.)
      {
        acc_assign(N2,En1,0.);
        acc_assign(N2,pVp1,0.);

        reduce_2c2(s1,s2,s1,s2,gs, valt1, valS2,valS3,valS4,iN,N,En1);
        reduce_2c2(s1,s2,s1,s2,gs3,valtv1,valV2,valV3,valV4,iN,N,pVp1);

       #pragma acc parallel loop present(En[0:N2],En1[0:N2],pVp[0:N2],pVp1[0:N2])
        for (int j=0;j<N2;j++)
        {
          En[j] += En1[j];
          pVp[j] += pVp1[j];
        }
      }

    } //loop p over nuclear center




   //two-center basis
    for (int n=m+1;n<natoms;n++)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      float Z2 = atz[n];
      int Z2g = atno[n];
      float A2 = coords[3*n+0]; float B2 = coords[3*n+1]; float C2 = coords[3*n+2];
      float A12 = A2-A1; float B12 = B2-B1; float C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2g,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1); 
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      acc_copyf(gs,wtt1,wt1);
      acc_copyf(gs,wtt2,wt2);
      becke_weight_2c(gs,grid1,wtt1,grid2,wtt2,Z1,Z2,A12,B12,C12);

      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wtt2);

     //needs to happen after Becke weighting
      add_r1_to_grid(gs,grid2,0.,0.,0.);
      add_r2_to_grid(gs,grid2,A12,B12,C12);


     #pragma acc parallel loop present(valS1[0:iN][0:gs],valV1[0:iN][0:gs3],valS2[0:iN][0:gs],valV2[0:iN][0:gs3])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS1[ii1][j] = valS2[ii1][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          valV1[ii1][j] = valV2[ii1][j] = 1.f;
      }

     #pragma acc parallel loop present(valS3[0:iN][0:gs],valV3[0:iN][0:gs3],valS4[0:iN][0:gs],valV4[0:iN][0:gs3])
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS3[ii2][j] = valS4[ii2][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs3;j++)
          valV3[ii2][j] = valV4[ii2][j] = 1.f;
      }

      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_sh(ii1,gs,grid1,valS1[ii1],n1,l1,m1,zeta1);
        eval_sh(ii1,gs,grid2,valS2[ii1],n1,l1,m1,zeta1);
        eval_p(gs,grid1,valV1[ii1],n1,l1,m1,zeta1);
        eval_p(gs,grid2,valV2[ii1],n1,l1,m1,zeta1);
      }

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1s,valS3[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid2s,valS4[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid1s,valV3[ii2],n2,l2,m2,zeta2);
        eval_p(gs,grid2s,valV4[ii2],n2,l2,m2,zeta2);
      }

     //2 basis on diff centers, 1 nucleus matches
     #pragma acc parallel loop present(wtt1[0:gs],valt1[0:iN][0:gs],valtv1[0:iN][0:gs3],valS1[0:iN][0:gs],valV1[0:iN][0:gs3])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valt1[ii1][j] = valS1[ii1][j]*wtt1[j];

       #pragma acc loop collapse(2)
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valtv1[ii1][3*j+k] = valV1[ii1][3*j+k]*wtt1[j];
      }

     #pragma acc parallel loop present(wtt2[0:gs],valt2[0:iN][0:gs],valtv2[0:iN][0:gs3],valS2[0:iN][0:gs],valV2[0:iN][0:gs3])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valt2[ii1][j] = valS2[ii1][j]*wtt2[j];

       #pragma acc loop collapse(2)
        for (int j=0;j<gs;j++)
        for (int k=0;k<3;k++)
          valtv2[ii1][3*j+k] = valV2[ii1][3*j+k]*wtt2[j];
      }

     #pragma acc parallel loop present(grid1[0:gs6],grid2[0:gs6],valt1[0:iN][0:gs],valt2[0:iN][0:gs],valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3])
      for (int j=0;j<gs;j++)
      {
        float Rn1a = grid1[6*j+3]+RPAD;
        float Rn2a = grid1[6*j+4]+RPAD;
        float Rn1b = grid2[6*j+3]+RPAD;
        float Rn2b = grid2[6*j+4]+RPAD;
        float ne1 = Z1/Rn1a+Z2/Rn2a;
        float ne2 = Z1/Rn1b+Z2/Rn2b;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valt1[ii1][j] *= ne1;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
          valt2[ii1][j] *= ne2;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
       #pragma acc loop
        for (int k=0;k<3;k++)
          valtv1[ii1][3*j+k] *= ne1;

       #pragma acc loop
        for (int ii1=0;ii1<s2-s1;ii1++)
       #pragma acc loop
        for (int k=0;k<3;k++)
          valtv2[ii1][3*j+k] *= ne2;
      }

     //collect 2c-atom same-atom values
      if (Z1!=0. || Z2!=0.)
      {
        acc_assign(N2,En1,0.);
        acc_assign(N2,pVp1,0.);

        reduce_2c2(s1,s2,s3,s4,gs,valt1,valt2,valS3,valS4,iN,N,En1);
        reduce_2c2(s1,s2,s3,s4,gs3,valtv1,valtv2,valV3,valV4,iN,N,pVp1);

       #pragma acc parallel loop present(En[0:N2],En1[0:N2])
        for (int j=0;j<N2;j++)
          En[j] += En1[j];
       #pragma acc parallel loop present(pVp[0:N2],pVp1[0:N2])
        for (int j=0;j<N2;j++)
          pVp[j] += pVp1[j];
      }


     //2 basis on diff centers, then 1 nucleus elsewhere
      for (int p=0;p<natoms;p++)
      if (p!=m && p!=n)
      {
        float Zn = atz[p];
        int Zng = atno[p];
        float An = coords[3*p+0]; float Bn = coords[3*p+1]; float Cn = coords[3*p+2];
        float A1n = An-A1; float B1n = Bn-B1; float C1n = Cn-C1;

        generate_central_grid_2(grid3,wt3,Zng,nrad,nang,ang_g,ang_w);
        //copy_grid(gs,grid3p,grid3); //grid 3 centered on atom 3
        recenter_grid(gs,grid3,A1n,B1n,C1n); //grid 3 centered on atom 1

        copy_grid(gs,grid3s,grid3);
        recenter_grid(gs,grid3s,-A12,-B12,-C12); //grid 3 centered on atom 2

        //copy_grid(gs,grid2p,grid2); //grid 2 centered on atom 3
        //recenter_grid(gs,grid2p,-A1n,-B1n,-C1n);
     
        //copy_grid(gs,grid1p,grid1);
        //recenter_grid(gs,grid1p,-A1n,-B1n,-C1n); //grid 1 centered on atom 3

        //add_r1_to_grid_6z(gs,grid1s,grid2s,grid3s,grid1p,grid2p,grid3p);
        add_r1_to_grid(gs,grid3s,0.,0.,0.);
      
      //need to keep all of these distances in order
        //add_r3_to_grid(gs,grid1,A1n,B1n,C1n);
        add_r123_to_grid(gs,grid1,0.,0.,0.,A12,B12,C12,A1n,B1n,C1n);
        add_r123_to_grid(gs,grid2,A12,B12,C12,0.,0.,0.,A1n,B1n,C1n);
      	add_r123_to_grid(gs,grid3,A1n,B1n,C1n,A12,B12,C12,0.,0.,0.);
     
        acc_copyf(gs,wtt1,wt1,wtt2,wt2);

        becke_weight_3c(gs,grid1,wtt1,grid2,wtt2,grid3,wt3,Z1,Z2,Zn,A12,B12,C12,A1n,B1n,C1n);
        eliminate_small_wt_3(estart,gs,wtt1,wtt2,wt3);

        add_r2_to_grid(gs,grid1,A1n,B1n,C1n);
        add_r1_to_grid(gs,grid2,0.,0.,0.);
        add_r1_to_grid(gs,grid3,0.,0.,0.);
        add_r2_to_grid(gs,grid2,A1n,B1n,C1n);
        add_r2_to_grid(gs,grid3,A1n,B1n,C1n);

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

         #pragma acc parallel loop collapse(2) present(valtv1[0:iN][0:gs3],valV1[0:iN][0:gs3],wtt1[0:gs])
          for (int j=0;j<gs;j++)
          for (int k=0;k<3;k++)
            valtv1[ii1][3*j+k] = valV1[ii1][3*j+k]*wtt1[j];

         #pragma acc parallel loop collapse(2) present(valtv2[0:iN][0:gs3],valV2[0:iN][0:gs3],wtt2[0:gs])
          for (int j=0;j<gs;j++)
          for (int k=0;k<3;k++)
            valtv2[ii1][3*j+k] = valV2[ii1][3*j+k]*wtt2[j];

         #pragma acc parallel loop collapse(2) present(valV5[0:iN][0:gs3],wt3[0:gs])
          for (int j=0;j<gs;j++)
          for (int k=0;k<3;k++)
            valV5[ii1][3*j+k] = wt3[j];
        }

       #pragma acc parallel loop present(valS6[0:iN][0:gs],valV6[0:iN][0:gs3])
        for (int ii2=0;ii2<s4-s3;ii2++)
        {
         #pragma acc loop
          for (int j=0;j<gs;j++)
            valS6[ii2][j] = 1.f;

         #pragma acc loop
          for (int j=0;j<gs3;j++)
            valV6[ii2][j] = 1.f;
        }

        for (int i1=s1;i1<s2;i1++)
        {
          int ii1 = i1-s1;

          vector<double> basis1 = basis[i1];
          int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

          eval_sh(ii1,gs,grid3,valS5[ii1],n1,l1,m1,zeta1);
          eval_p(gs,grid3,valV5[ii1],n1,l1,m1,zeta1);
        }

        for (int i2=s3;i2<s4;i2++)
        {
          int ii2 = i2-s3;
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

          eval_sh(ii2,gs,grid3s,valS6[ii2],n2,l2,m2,zeta2);
          eval_p(gs,grid3s,valV6[ii2],n2,l2,m2,zeta2);
        } //loop i2 evaluate

       #pragma acc parallel loop present(grid1[0:gs6],grid2[0:gs6],grid3[0:gs6],valt1[0:iN][0:gs],valt2[0:iN][0:gs],valS5[0:iN][0:gs],valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3],valV5[0:iN][0:gs3])
        for (int j=0;j<gs;j++)
        {
          float Rn1 = grid1[6*j+4]+RPAD; float Rn2 = grid2[6*j+4]+RPAD; float Rn3 = grid3[6*j+4]+RPAD;
          float ne1 = Zn/Rn1; float ne2 = Zn/Rn2; float ne3 = Zn/Rn3;

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
            valtv1[ii1][3*j+k] *= ne1;
            valtv2[ii1][3*j+k] *= ne2;
            valV5[ii1][3*j+k] *= ne3;
          }
        }

        //collect 3-atom values
        if (Zn!=0.)
        {
          acc_assign(N2,En1,0.);
          acc_assign(N2,pVp1,0.);

          reduce_2c3(s1,s2,s3,s4,gs, valt1, valt2, valS3,valS4,valS5,valS6,iN,N,En1);
          reduce_2c3(s1,s2,s3,s4,gs3,valtv1,valtv2,valV3,valV4,valV5,valV6,iN,N,pVp1);

         #pragma acc parallel loop present(En[0:N2],En1[0:N2],pVp[0:N2],pVp1[0:N2])
          for (int j=0;j<N2;j++) 
          {
            En[j] += En1[j];
            pVp[j] += pVp1[j];
          }
        }

      } //loop p over nuclear center

    } //loop n over second atom


  } //loop m over natoms

  double* norm = new double[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];
  #pragma acc enter data copyin(norm[0:N])

 #pragma acc parallel loop independent present(En[0:N2],pVp[0:N2],norm[0:N])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i;j<N;j++)
  {
    double n12 = norm[i]*norm[j];
    En[i*N+j]  *= -n12;
    pVp[i*N+j] *= -n12;
  }

 #pragma acc parallel loop independent present(En[0:N2],pVp[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i+1;j<N;j++)
  {
    En[j*N+i] = En[i*N+j];
    pVp[j*N+i] = pVp[i*N+j];
  }

  #pragma acc exit data delete(norm[0:N])
  delete [] norm;

  clean_small_values(N,En);
  clean_small_values(N,pVp);

 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc exit data delete(En1[0:N2],pVp1[0:N2])
  //#pragma acc exit data copyout(En[0:N2],pVp[0:N2])
  #pragma acc update self(En[0:N2],pVp[0:N2])
 #endif


  if (prl>1)
  {
    printf("\n En: \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %10.5f",En[i*N+j]);
      printf("\n");
    }
  }

  if (prl>2)
  {
    printf("\n pVp: \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %12.6f",pVp[i*N+j]);
      printf("\n");
    }
    printf("\n");
  }

 //CPMZ check this
#if USE_ACC
  #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
  #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs])
  #pragma acc exit data delete(grid3[0:gs6],wt3[0:gs])
  #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6],grid3s[0:gs6],grid1p[0:gs6],grid2p[0:gs6],grid3p[0:gs6])
  #pragma acc exit data delete(wtt1[0:gs],wtt2[0:gs])
  #pragma acc exit data delete(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs],valS5[0:iN][0:gs],valS6[0:iN][0:gs])
  #pragma acc exit data delete(valV1[0:iN][0:gs3],valV2[0:iN][0:gs3],valV3[0:iN][0:gs3],valV4[0:iN][0:gs3],valV5[0:iN][0:gs3],valV6[0:iN][0:gs3])
  #pragma acc exit data delete(valt1[0:iN][0:gs],valt2[0:iN][0:gs],valt3[0:iN][0:gs],valtv1[0:iN][0:gs3],valtv2[0:iN][0:gs3],valtv3[0:iN][0:gs3])
  #pragma acc exit data delete(n2i[0:natoms])
  #pragma acc exit data delete(coords[0:3*natoms],atno[0:natoms])
#endif

  delete [] ang_g;
  delete [] ang_w;

  delete [] pVp1;
  delete [] En1;

  delete [] n2i;

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
  for (int i=0;i<iN;i++) delete [] valV1[i];
  for (int i=0;i<iN;i++) delete [] valV2[i];
  for (int i=0;i<iN;i++) delete [] valV3[i];
  for (int i=0;i<iN;i++) delete [] valV4[i];
  for (int i=0;i<iN;i++) delete [] valV5[i];
  for (int i=0;i<iN;i++) delete [] valV6[i];

  for (int i=0;i<iN;i++) delete [] valt1[i];
  for (int i=0;i<iN;i++) delete [] valt2[i];
  for (int i=0;i<iN;i++) delete [] valt3[i];
  for (int i=0;i<iN;i++) delete [] valtv1[i];
  for (int i=0;i<iN;i++) delete [] valtv2[i];
  for (int i=0;i<iN;i++) delete [] valtv3[i];

  delete [] valS1; delete [] valS2; delete [] valS3; delete [] valS4; delete [] valS5; delete [] valS6;
  delete [] valV1; delete [] valV2; delete [] valV3; delete [] valV4; delete [] valV5; delete [] valV6;
  delete [] valt1; delete [] valt2; delete [] valt3; delete [] valtv1; delete [] valtv2; delete [] valtv3;
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

//applied electric fields
// currently does x,y,z directions only
void compute_Exyz(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g, double* ang_w, double* E, int prl)
{
  if (prl>1) printf(" beginning compute_E (double precision) \n");

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
      eval_shd(ii1,gs,grid1m,val1m,n1,l1,m1,zeta1); //basis 1
      eval_shd(ii1,gs,grid1m,val2m,n2,l2,m2,zeta2); //basis 2

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

        copy_grid(gs,grid2n,grid2m);
        recenter_grid(gs,grid2n,-A12,-B12,-C12);      //grid 2 centered on atom 1

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
        eval_shd(ii1,gs,grid1m,val1m,n1,l1,m1,zeta1); //basis 1 on center 1
        eval_shd(ii1,gs,grid2m,val1n,n1,l1,m1,zeta1); //basis 1 on center 2
        eval_shd(ii1,gs,grid1n,val2m,n2,l2,m2,zeta2); //basis 2 on center 1
        eval_shd(ii1,gs,grid2n,val2n,n2,l2,m2,zeta2); //basis 2 on center 2

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

//high-precision overlap integrals
void compute_Sd(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g, double* ang_w, double* S, int prl)
{
  if (prl>1) printf(" beginning compute_S (double precision) \n");

  int N = basis.size();
  int N2 = N*N;

  printf("  compute_Sd. nrad/ang: %3i %3i \n",nrad,nang);

  if (prl>2)
  {
    printf("\n S(float): \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %15.12f",S[i*N+j]);
      printf("\n");
    }
    printf("\n");
  }

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

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
  acc_assign(N2,S,0.);

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

      double z12 = zeta1 + zeta2;
     //new grid with zeta dependence
      generate_central_grid_2d(-1,0,grid1m,wt1,z12,nrad,nang,ang_g,ang_w);

      #pragma acc parallel loop present(val1m[0:gs],wt1[0:gs])
      for (int j=0;j<gs;j++)
        val1m[j] = wt1[j];
      #pragma acc parallel loop present(val2m[0:gs])
      for (int j=0;j<gs;j++)
        val2m[j] = 1.;

     //S
      eval_shd(ii1,gs,grid1m,val1m,n1,l1,m1,zeta1); //basis 1
      eval_shd(ii1,gs,grid1m,val2m,n2,l2,m2,zeta2); //basis 2

     #if 0
      #pragma acc update self(val1m[0:gs],val2m[0:gs],wt1[0:gs])
      printf("   val1*val2*wt: ");
      for (int j=0;j<nrad;j++)
        printf(" %4.1e",val1m[j*nang]*val2m[j*nang]);
      printf("\n");
     #endif

      double val = 0.;
     #pragma acc parallel loop present(val1m[0:gs],val2m[0:gs]) reduction(+:val)
      for (int j=0;j<gs;j++)
        val += val1m[j]*val2m[j];

     #pragma acc serial present(S[0:N2])
      S[i1*N+i2] = S[i2*N+i1] = val;

    } //pairs of basis on single atoms


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

        copy_grid(gs,grid2n,grid2m);
        recenter_grid(gs,grid2n,-A12,-B12,-C12);      //grid 2 centered on atom 1

        copy_grid(gs,grid1n,grid1m);
        recenter_grid_zero(gs,grid1n,-A12,-B12,-C12); //grid 1 centered on atom 2

       //needs to happen after becke weighting
        add_r1_to_grid(gs,grid2m,0.,0.,0.);

        #pragma acc parallel loop present(val2n[0:gs],val2m[0:gs])
        for (int j=0;j<gs;j++)
          val2m[j] = val2n[j] = 1.;
        #pragma acc parallel loop present(val1m[0:gs],val1n[0:gs],wt1[0:gs],wt2[0:gs])
        for (int j=0;j<gs;j++)
        {
          val1m[j] = wt1[j];
          val1n[j] = wt2[j];
        }

       //S
        eval_shd(ii1,gs,grid1m,val1m,n1,l1,m1,zeta1); //basis 1 on center 1
        eval_shd(ii1,gs,grid2m,val1n,n1,l1,m1,zeta1); //basis 1 on center 2
        eval_shd(ii1,gs,grid1n,val2m,n2,l2,m2,zeta2); //basis 2 on center 1
        eval_shd(ii1,gs,grid2n,val2n,n2,l2,m2,zeta2); //basis 2 on center 2

        double val = 0.;
       #pragma acc parallel loop present(val1m[0:gs],val1n[0:gs],val2m[0:gs],val2n[0:gs]) reduction(+:val)
        for (int j=0;j<gs;j++)
          val += val1m[j]*val2m[j] + val1n[j]*val2n[j];

       #pragma acc serial present(S[0:N2])
        S[i1*N+i2] = S[i2*N+i1] = val;
      }

    } //loop n>m

  } //loop m over natoms


  double* norm = new double[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];
  #pragma acc enter data copyin(norm[0:N])

  #pragma acc parallel loop independent present(S[0:N2],norm[0:N])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=0;j<=i;j++)
  {
    double n12 = norm[i]*norm[j];
    S[i*N+j] *= n12;
  }

  if (prl>-1)
  {
    int nlow = 0;
    #pragma acc update self(S[0:N2])
    printf("  overlap diagonal accuracy (decimal points): ");
    for (int i=0;i<N;i++)
    {
      double v0 = fabs(S[i*N+i]-1.);
      double v1 = log10(v0);
      if (fabs(v1)<4) nlow++;
      printf(" %5.3f",v1);
    }
    printf("\n\n");
    if (nlow)
      printf("  WARNING: found %2i low accuracy diagonals \n",nlow);
    if (nlow>1) { printf("   therefore exiting now \n"); exit(-1); }
  }

 #pragma acc parallel loop independent present(S[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=0;j<i;j++)
  {
    S[j*N+i] = S[i*N+j];
  }

 //might as well eliminate errors on diagonal
 #pragma acc parallel loop present(S[0:N2])
  for (int i=0;i<N;i++)
    S[i*N+i] = 1.;

  #pragma acc exit data delete(norm[0:N])
  delete [] norm;

  clean_small_values(N,S);


 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc update self(S[0:N2])
 #endif


  if (prl>1 || prl==-1)
  {
    printf("\n S: \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %15.12f",S[i*N+j]);
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

  //delete [] ang_g;
  //delete [] ang_w;

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

#if RED_DOUBLE
void compute_ST(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* S, double* T, int prl)
#else
void compute_ST(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, float* S, float* T, int prl)
#endif
{
  if (prl>1) printf(" beginning compute_ST \n");

  int N = basis.size();
  int N2 = N*N;

  int gs = nrad*nang;
  int gs6 = 6*gs;

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

  int estart = find_center_of_grid(1,nrad)*nang;

  float* grid1 = new float[gs6];
  float* wt1 = new float[gs];

  float* grid2 = new float[gs6];
  float* wt2 = new float[gs];

  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);
  //printf("  iN: %i \n",imaxN);

 //intermediate storage
  int iN = imaxN;
  float* grid1s = new float[gs6];
  float* grid2s = new float[gs6];
  float** valS1 = new float*[iN]; for (int i=0;i<iN;i++) valS1[i] = new float[gs];
  float** valS2 = new float*[iN]; for (int i=0;i<iN;i++) valS2[i] = new float[gs];
  float** valS3 = new float*[iN]; for (int i=0;i<iN;i++) valS3[i] = new float[gs];
  float** valS4 = new float*[iN]; for (int i=0;i<iN;i++) valS4[i] = new float[gs];
  float** valT1 = new float*[iN]; for (int i=0;i<iN;i++) valT1[i] = new float[gs];
  float** valT2 = new float*[iN]; for (int i=0;i<iN;i++) valT2[i] = new float[gs];
  float* wtt1 = new float[gs];

  float* ang_g = new float[3*nang];
  float* ang_w = new float[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

 #if USE_ACC
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data copyin(n2i[0:natoms])
  #pragma acc enter data copyin(coords[0:3*natoms],atno[0:natoms])

  #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
  #pragma acc enter data create(grid2[0:gs6],wt2[0:gs])
  #pragma acc enter data create(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs])
  #pragma acc enter data create(valT1[0:iN][0:gs],valT2[0:iN][0:gs])
  #pragma acc enter data create(grid1s[0:gs6],grid2s[0:gs6],wtt1[0:gs])
  //#pragma acc enter data create(S[0:N2],T[0:N2])
 #endif
  acc_assign(N2,S,0.);
  acc_assign(N2,T,0.);

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    float Z1 = (float)atno[m];
    float A1 = coords[3*m+0]; float B1 = coords[3*m+1]; float C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);

   #pragma acc parallel loop present(valS1[0:iN][0:gs])
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valS1[ii1][j] = 1.f;
    }

   #pragma acc parallel loop present(valS3[0:iN][0:gs],wt1[0:gs])
    for (int ii2=0;ii2<s2-s1;ii2++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        valS3[ii2][j] = wt1[j];
    }

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

     //S
      eval_sh(ii1,gs,grid1,valS1[ii1],n1,l1,m1,zeta1);
    }

    for (int i2=s1;i2<s2;i2++)
    {
      int ii2 = i2-s1;
      vector<double> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

     //S
      eval_sh(ii2,gs,grid1,valS3[ii2],n2,l2,m2,zeta2);
    } //loop i2 evaluate

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

      eval_ke(gs,grid1,valT1[ii1],n1,l1,zeta1);
    }

   #if USE_ACC
    #pragma acc wait
   #endif

    reduce_2c1(s1,s2,gs,valS1,valS3,iN,N,S);
    reduce_2c1(s1,s2,gs,valT1,valS3,iN,N,T);


   //two-atom ints
   #if SYMM_ST
    for (int n=0;n<natoms;n++)
    if (m!=n)
   #else
    for (int n=m+1;n<natoms;n++)
   #endif
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];
      //printf(" mn: %i %i s1-4: %i %i - %i %i \n",m,n,s1,s2,s3,s4);

      float Z2 = (float)atno[n];
      float A2 = coords[3*n+0]; float B2 = coords[3*n+1]; float C2 = coords[3*n+2];
      float A12 = A2-A1; float B12 = B2-B1; float C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      acc_copyf(gs,wtt1,wt1);
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      becke_weight_2c(gs,grid1,wtt1,grid2,wt2,Z1,Z2,A12,B12,C12);

    #if TEST_SORT
     #pragma acc host_data use_device(wtt1,grid1)
      sortByKeys(wtt1,(xyz<fp_t>*)grid1,gs);
     #pragma acc host_data use_device(wtt1,grid1)
      sortByKeys(wt2,(xyz<fp_t>*)grid2,gs);
    #else
      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wt2);
    #endif

      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2s,-A12,-B12,-C12);

      copy_grid(gs,grid1s,grid1);
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

     //needs to happen after Becke weighting
      add_r1_to_grid(gs,grid2,0.,0.,0.);

     #pragma acc parallel loop present(valS1[0:iN][0:gs],valS2[0:iN][0:gs])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS1[ii1][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS2[ii1][j] = 1.f;
      }

     #pragma acc parallel loop present(valS3[0:iN][0:gs],valS4[0:iN][0:gs],wtt1[0:gs],wt2[0:gs])
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS3[ii2][j] = wtt1[j];
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valS4[ii2][j] = wt2[j];
      }

      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

       //S
        eval_sh(ii1,gs,grid1,valS1[ii1],n1,l1,m1,zeta1);
        eval_sh(ii1,gs,grid2,valS2[ii1],n1,l1,m1,zeta1);
      }

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

       //S
        eval_sh(ii2,gs,grid2s,valS4[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid1s,valS3[ii2],n2,l2,m2,zeta2);
      }

     #pragma acc parallel loop present(valS1[0:iN][0:gs],valT1[0:iN][0:gs])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valT1[ii1][j] = valS1[ii1][j];
      }
     #pragma acc parallel loop present(valS2[0:iN][0:gs],valT2[0:iN][0:gs])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valT2[ii1][j] = valS2[ii1][j];
      }

     //KE terms
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_ke(gs,grid1,valT1[ii1],n1,l1,zeta1);
        eval_ke(gs,grid2,valT2[ii1],n1,l1,zeta1);
      }

      reduce_2c2(s1,s2,s3,s4,gs,valS1,valS2,valS3,valS4,iN,N,S);
      reduce_2c2(s1,s2,s3,s4,gs,valT1,valT2,valS3,valS4,iN,N,T);

    } //loop n over second atom

  } //loop m over natoms


  double* norm = new double[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];
  #pragma acc enter data copyin(norm[0:N])

  #pragma acc parallel loop independent present(S[0:N2],T[0:N2],norm[0:N])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=0;j<N;j++)
  {
    double n12 = norm[i]*norm[j];
    S[i*N+j] *= n12;
    T[i*N+j] *= -0.5*n12;
  }

 #if !SYMM_ST
 //symmetrize wrt ij
 #pragma acc parallel loop independent present(S[0:N2],T[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i+1;j<N;j++)
  {
    S[j*N+i] = S[i*N+j];
    T[j*N+i] = T[i*N+j];
  }
 #else

 //average the two (atom m vs n)
 #pragma acc parallel loop independent present(S[0:N2],T[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i;j<N;j++)
  {
    S[j*N+i] += S[i*N+j];
    T[j*N+i] += T[i*N+j];
  }

 #pragma acc parallel loop independent present(S[0:N2],T[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i+1;j<N;j++)
  {
    S[i*N+j] = S[j*N+i];
    T[i*N+j] = T[j*N+i];
  }

 #pragma acc parallel loop independent present(S[0:N2],T[0:N2])
  for (int i=0;i<N2;i++)
  {
    S[i] *= 0.5;
    T[i] *= 0.5;
  }

 #endif

 //might as well eliminate errors on diagonal
 #pragma acc parallel loop present(S[0:N2])
  for (int i=0;i<N;i++)
    S[i*N+i] = 1.;

  #pragma acc exit data delete(norm[0:N])
  delete [] norm;

  clean_small_values(N,S);
  clean_small_values(N,T);


 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  //#pragma acc exit data copyout(S[0:N2],T[0:N2])
  #pragma acc update self(S[0:N2],T[0:N2])
 #endif


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

 //CPMZ check this
#if USE_ACC
  #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
  #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs])
  #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6],wtt1[0:gs])
  #pragma acc exit data delete(valS1[0:iN][0:gs],valS2[0:iN][0:gs],valS3[0:iN][0:gs],valS4[0:iN][0:gs])
  #pragma acc exit data delete(valT1[0:iN][0:gs],valT2[0:iN][0:gs])
  #pragma acc exit data delete(n2i[0:natoms])
  #pragma acc exit data delete(coords[0:3*natoms],atno[0:natoms])
#endif

  delete [] ang_g;
  delete [] ang_w;

  delete [] n2i;

  delete [] grid1s;
  delete [] grid2s;

  for (int i=0;i<iN;i++) delete [] valS1[i];
  for (int i=0;i<iN;i++) delete [] valS2[i];
  for (int i=0;i<iN;i++) delete [] valS3[i];
  for (int i=0;i<iN;i++) delete [] valS4[i];
  for (int i=0;i<iN;i++) delete [] valT1[i];
  for (int i=0;i<iN;i++) delete [] valT2[i];
  delete [] valS1; delete [] valS2; delete [] valS3; delete [] valS4;
  delete [] valT1; delete [] valT2;
  delete [] wtt1;

  delete [] grid1;
  delete [] grid2;
  delete [] wt1;
  delete [] wt2;

  return;
}

//high precision 2-center Coulomb
void compute_all_2c_v2d(bool do_overlap, int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* An, int prl)
{
  if (prl>1) printf(" beginning compute_all_2c_v2d \n");

 //2c integrals are all in auxiliary basis
  int N = basis.size();
  int N2 = N*N;

  int gs = nrad*nang;
  int gs6 = 6*gs;

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

  int estart = find_center_of_grid(1,nrad)*nang;

  double* grid1 = new double[gs6];
  double* wt1 = new double[gs];

  double* grid2 = new double[gs6];
  double* wt2 = new double[gs];

  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);
  //printf("  iN: %i \n",imaxN);

 //intermediate storage
  int iN = imaxN;
  double* grid1s = new double[gs6];
  double* grid2s = new double[gs6];
  double** val1 = new double*[iN];
  double** val2 = new double*[iN];
  double** val3 = new double*[iN];
  double** val4 = new double*[iN];
  for (int i=0;i<iN;i++)
    val1[i] = new double[gs];
  for (int i=0;i<iN;i++)
    val2[i] = new double[gs];
  for (int i=0;i<iN;i++)
    val3[i] = new double[gs];
  for (int i=0;i<iN;i++)
    val4[i] = new double[gs];
  double* wtt1 = new double[gs];

  double* ang_g = new double[3*nang];
  double* ang_w = new double[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

 #if USE_ACC
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data copyin(n2i[0:natoms])

  #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
  #pragma acc enter data create(grid2[0:gs6],wt2[0:gs])
  #pragma acc enter data create(val1[0:iN][0:gs],val2[0:iN][0:gs],wtt1[0:gs])
  #pragma acc enter data create(val3[0:iN][0:gs],val4[0:iN][0:gs])
  #pragma acc enter data create(grid1s[0:gs6],grid2s[0:gs6])
  //#pragma acc enter data create(An[0:N2])
 #endif
  acc_assign(N2,An,0.);

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    double Z1 = (float)atno[m];
    float A1 = coords[3*m+0]; float B1 = coords[3*m+1]; float C1 = coords[3*m+2];

    generate_central_grid_2d(-1,1,grid1,wt1,Z1,nrad,nang,ang_g,ang_w);
    //float z12 = zeta1 + zeta2;
     //new grid with zeta dependence
    //generate_central_grid_2d(-1,0,grid1m,wt1,z12,nrad,nang,ang_g,ang_w);

   #pragma acc parallel loop present(val1[0:iN][0:gs])
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val1[ii1][j] = 1.;
    }

   #pragma acc parallel loop present(val3[0:iN][0:gs],wt1[0:gs])
    for (int ii2=0;ii2<s2-s1;ii2++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val3[ii2][j] = wt1[j];
    }

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

      if (do_overlap)
        eval_shd(ii1,gs,grid1,val1[ii1],n1,l1,m1,zeta1);
      else
      {
        eval_inr_r12(-1,gs,grid1,val1[ii1],n1,l1,zeta1);
        eval_sh_3rd(gs,grid1,val1[ii1],n1,l1,m1);
      }
    }

    for (int i2=s1;i2<s2;i2++)
    {
      int ii2 = i2-s1;
      vector<double> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

      eval_shd(ii2,gs,grid1,val3[ii2],n2,l2,m2,zeta2);

    } //loop i2 evaluate
   #if USE_ACC
    #pragma acc wait
   #endif

    reduce_2c1(-1,s1,s2,gs,val1,val3,iN,N,An);

   //two-atom ints
    for (int n=m+1;n<natoms;n++)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];
      //printf(" mn: %i %i s1-4: %i %i - %i %i \n",m,n,s1,s2,s3,s4);

      float Z2 = (float)atno[n];
      float A2 = coords[3*n+0]; float B2 = coords[3*n+1]; float C2 = coords[3*n+2];
      float A12 = A2-A1; float B12 = B2-B1; float C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2d(-1,1,grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1);
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      acc_copy(gs,wtt1,wt1);
      //becke_weight_2c(gs,grid1,wtt1,grid2,wt2,zeta1,zeta2,A12,B12,C12);
      becke_weight_2d(gs,grid1,wtt1,grid2,wt2,Z1,Z2,A12,B12,C12);

      //eliminate_small_wt(estart,gs,wtt1);
      //eliminate_small_wt(estart,gs,wt2);

     //needs to happen after Becke weighting
      add_r1_to_grid(gs,grid2,0.,0.,0.);

     #pragma acc parallel loop present(val1[0:iN][0:gs],val2[0:iN][0:gs])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val1[ii1][j] = 1.;
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val2[ii1][j] = 1.;
      }

     #pragma acc parallel loop present(val3[0:iN][0:gs],val4[0:iN][0:gs],wtt1[0:gs],wt2[0:gs])
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val3[ii2][j] = wtt1[j];
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val4[ii2][j] = wt2[j];
      }

      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        if (do_overlap)
          eval_shd(ii1,gs,grid1,val1[ii1],n1,l1,m1,zeta1);
        else
        {
          eval_inr_r12(-1,gs,grid1,val1[ii1],n1,l1,zeta1);
          eval_sh_3rd(gs,grid1,val1[ii1],n1,l1,m1);
        }

        if (do_overlap)
          eval_shd(ii1,gs,grid2,val2[ii1],n1,l1,m1,zeta1);
        else
        {
          eval_inr_r12(-1,gs,grid2,val2[ii1],n1,l1,zeta1);
          eval_sh_3rd(gs,grid2,val2[ii1],n1,l1,m1);
        }
      }

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        eval_shd(ii2,gs,grid2s,val4[ii2],n2,l2,m2,zeta2);
        eval_shd(ii2,gs,grid1s,val3[ii2],n2,l2,m2,zeta2);
      }

      reduce_2c2(s1,s2,s3,s4,gs,val1,val2,val3,val4,iN,N,An);

    } //loop n over second atom

  } //loop m over natoms

  double* norm1 = new double[N];
  double* norm2 = new double[N];
  for (int i=0;i<N;i++)
  {
    if (do_overlap)
      norm1[i] = basis[i][4];
    else
      norm1[i] = norm_sv(basis[i][0],basis[i][1],basis[i][2],basis[i][3]);
    norm2[i] = basis[i][4];
  }
  #pragma acc enter data copyin(norm1[0:N],norm2[0:N])

 #pragma acc parallel loop independent present(An[0:N2],norm1[0:N],norm2[0:N])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i;j<N;j++)
  //for (int j=0;j<N;j++)
  {
    double n12 = norm1[i]*norm2[j];
    An[i*N+j] *= n12;
  }

  #pragma acc exit data delete(norm1[0:N],norm2[0:N])
  delete [] norm1;
  delete [] norm2;

 #if 1
 #pragma acc parallel loop independent present(An[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i+1;j<N;j++)
    An[j*N+i] = An[i*N+j];

 #else
 
 //average the two triangles
 #pragma acc parallel loop independent present(An[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i;j<N;j++)
    An[j*N+i] += An[i*N+j];

 #pragma acc parallel loop independent present(An[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i+1;j<N;j++)
    An[i*N+j] = An[j*N+i];

 #pragma acc parallel loop independent present(An[0:N2])
  for (int i=0;i<N2;i++)
    An[i] *= 0.5;

 #endif

  //clean_small_values(N,An);

 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  //#pragma acc exit data copyout(An[0:N2])
  #pragma acc update self(An[0:N2])
 #endif


 #if 0
  for (int m=0;m<N;m++)
  {
    double val = An[m*N+m];
    if (fabs(val)<1.e-2)
    {
      printf(" WARNING: small diagonal element in A: %3i \n",m);
    }
    double v2 = 1./sqrt(val);
    //Anorm[m] = v2;
    for (int n=0;n<N;n++)
    {
      An[m*N+n] *= v2;
      An[n*N+m] *= v2;
    }
  }
 #endif


  if (prl>2)
  {
    printf("\n A: \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %12.6f",An[i*N+j]);
      printf("\n");
    }
  }

 //CPMZ check this
#if USE_ACC
  #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
  #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs])
  #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6])
  #pragma acc exit data delete(val1[0:iN][0:gs],val2[0:iN][0:gs],wtt1[0:gs])
  #pragma acc exit data delete(val3[0:iN][0:gs],val4[0:iN][0:gs])
  #pragma acc exit data delete(n2i[0:natoms])
#endif

  delete [] ang_g;
  delete [] ang_w;

  delete [] n2i;

  delete [] grid1s;
  delete [] grid2s;

  for (int i=0;i<iN;i++) delete [] val1[i];
  for (int i=0;i<iN;i++) delete [] val2[i];
  for (int i=0;i<iN;i++) delete [] val3[i];
  for (int i=0;i<iN;i++) delete [] val4[i];
  delete [] val1;
  delete [] val2;
  delete [] val3;
  delete [] val4;
  delete [] wtt1;

  delete [] grid1;
  delete [] grid2;
  delete [] wt1;
  delete [] wt2;

  return;
}

#if RED_DOUBLE
void compute_all_2c_v2(bool do_overlap, int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* An, int prl)
#else
void compute_all_2c_v2(bool do_overlap, int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, float* An, int prl)
#endif
{
  if (prl>1) printf(" beginning compute_all_2c_v2 \n");

 //2c integrals are all in auxiliary basis
  int N = basis.size();
  int N2 = N*N;

  int gs = nrad*nang;
  int gs6 = 6*gs;

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

  int estart = find_center_of_grid(1,nrad)*nang;

  float* grid1 = new float[gs6];
  float* wt1 = new float[gs];

  float* grid2 = new float[gs6];
  float* wt2 = new float[gs];

  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);
  //printf("  iN: %i \n",imaxN);

 //intermediate storage
  int iN = imaxN;
  float* grid1s = new float[gs6];
  float* grid2s = new float[gs6];
  float** val1 = new float*[iN];
  float** val2 = new float*[iN];
  float** val3 = new float*[iN];
  float** val4 = new float*[iN];
  for (int i=0;i<iN;i++)
    val1[i] = new float[gs];
  for (int i=0;i<iN;i++)
    val2[i] = new float[gs];
  for (int i=0;i<iN;i++)
    val3[i] = new float[gs];
  for (int i=0;i<iN;i++)
    val4[i] = new float[gs];
  float* wtt1 = new float[gs];

  float* ang_g = new float[3*nang];
  float* ang_w = new float[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

 #if USE_ACC
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data copyin(n2i[0:natoms])

  #pragma acc enter data create(grid1[0:gs6],wt1[0:gs])
  #pragma acc enter data create(grid2[0:gs6],wt2[0:gs])
  #pragma acc enter data create(val1[0:iN][0:gs],val2[0:iN][0:gs],wtt1[0:gs])
  #pragma acc enter data create(val3[0:iN][0:gs],val4[0:iN][0:gs])
  #pragma acc enter data create(grid1s[0:gs6],grid2s[0:gs6])
  //#pragma acc enter data create(An[0:N2])
 #endif
  acc_assign(N2,An,0.);

  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    float Z1 = (float)atno[m];
    float A1 = coords[3*m+0]; float B1 = coords[3*m+1]; float C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);

   #pragma acc parallel loop present(val1[0:iN][0:gs])
    for (int ii1=0;ii1<s2-s1;ii1++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val1[ii1][j] = 1.f;
    }

   #pragma acc parallel loop present(val3[0:iN][0:gs],wt1[0:gs])
    for (int ii2=0;ii2<s2-s1;ii2++)
    {
     #pragma acc loop
      for (int j=0;j<gs;j++)
        val3[ii2][j] = wt1[j];
    }

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

      //acc_assign(gs,val1[ii1],1.);

      if (do_overlap)
        eval_sh(ii1,gs,grid1,val1[ii1],n1,l1,m1,zeta1);
      else
      {
        eval_inr_r12(gs,grid1,val1[ii1],n1,l1,zeta1);
        eval_sh_3r(gs,grid1,val1[ii1],n1,l1,m1);
      }
    }

    for (int i2=s1;i2<s2;i2++)
    {
      int ii2 = i2-s1;
      vector<double> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

      eval_sh(ii2,gs,grid1,val3[ii2],n2,l2,m2,zeta2);

    } //loop i2 evaluate
   #if USE_ACC
    #pragma acc wait
   #endif

    reduce_2c1(s1,s2,gs,val1,val3,iN,N,An);

   //two-atom ints
    for (int n=m+1;n<natoms;n++)
    //for (int n=0;n<natoms;n++)
    //if (m!=n)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];
      //printf(" mn: %i %i s1-4: %i %i - %i %i \n",m,n,s1,s2,s3,s4);

      float Z2 = (float)atno[n];
      float A2 = coords[3*n+0]; float B2 = coords[3*n+1]; float C2 = coords[3*n+2];
      float A12 = A2-A1; float B12 = B2-B1; float C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1); 
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid2,wt2,Z1,Z2,A12,B12,C12);

      //#pragma acc update self(grid1[0:gs6],wtt1[0:gs],grid2[0:gs6],wt2[0:gs])
      //print_grid(gs,grid1,grid2,wtt1,wt2,prl);

      eliminate_small_wt(estart,gs,wtt1);
      eliminate_small_wt(estart,gs,wt2);

     //needs to happen after Becke weighting
      add_r1_to_grid(gs,grid2,0.,0.,0.);

     #pragma acc parallel loop present(val1[0:iN][0:gs],val2[0:iN][0:gs])
      for (int ii1=0;ii1<s2-s1;ii1++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val1[ii1][j] = 1.f;
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val2[ii1][j] = 1.f;
      }

     #pragma acc parallel loop present(val3[0:iN][0:gs],val4[0:iN][0:gs],wtt1[0:gs],wt2[0:gs])
      for (int ii2=0;ii2<s4-s3;ii2++)
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val3[ii2][j] = wtt1[j];
       #pragma acc loop
        for (int j=0;j<gs;j++)
          val4[ii2][j] = wt2[j];
      }

      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;

        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];
        //printf("  m: %i i1: %i   nlm: %i %i %2i zeta: %8.5f \n",m,i1,n1,l1,m1,zeta1);

        //acc_assign(gs,val1[ii1],1.);
        //acc_assign(gs,val2[ii1],1.);

        if (do_overlap)
          eval_sh(ii1,gs,grid1,val1[ii1],n1,l1,m1,zeta1);
        else
        {
          eval_inr_r12(gs,grid1,val1[ii1],n1,l1,zeta1);
          eval_sh_3r(gs,grid1,val1[ii1],n1,l1,m1);
        }

        if (do_overlap)
          eval_sh(ii1,gs,grid2,val2[ii1],n1,l1,m1,zeta1);
        else
        {
          eval_inr_r12(gs,grid2,val2[ii1],n1,l1,zeta1);
          eval_sh_3r(gs,grid2,val2[ii1],n1,l1,m1);
        }
      }

      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];
        //printf("   n: %i i2: %i   nlm: %i %i %2i zeta: %8.5f \n",n,i2,n2,l2,m2,zeta2);

        eval_sh(ii2,gs,grid2s,val4[ii2],n2,l2,m2,zeta2);
        eval_sh(ii2,gs,grid1s,val3[ii2],n2,l2,m2,zeta2);
      }

      reduce_2c2(s1,s2,s3,s4,gs,val1,val2,val3,val4,iN,N,An);

    } //loop n over second atom

  } //loop m over natoms

  double* norm1 = new double[N];
  double* norm2 = new double[N];
  for (int i=0;i<N;i++)
  {
    if (do_overlap)
      norm1[i] = basis[i][4];
    else
      norm1[i] = norm_sv(basis[i][0],basis[i][1],basis[i][2],basis[i][3]);
    norm2[i] = basis[i][4];
  }
  #pragma acc enter data copyin(norm1[0:N],norm2[0:N])

 #pragma acc parallel loop independent present(An[0:N2],norm1[0:N],norm2[0:N])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i;j<N;j++)
  //for (int j=0;j<N;j++)
  {
    double n12 = norm1[i]*norm2[j];
    An[i*N+j] *= n12;
  }

  #pragma acc exit data delete(norm1[0:N],norm2[0:N])
  delete [] norm1;
  delete [] norm2;

 #if 1
 #pragma acc parallel loop independent present(An[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i+1;j<N;j++)
    An[j*N+i] = An[i*N+j];

 #else
 
 //average the two triangles
 #pragma acc parallel loop independent present(An[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i;j<N;j++)
    An[j*N+i] += An[i*N+j];

 #pragma acc parallel loop independent present(An[0:N2])
  for (int i=0;i<N;i++)
 #pragma acc loop independent
  for (int j=i+1;j<N;j++)
    An[i*N+j] = An[j*N+i];

 #pragma acc parallel loop independent present(An[0:N2])
  for (int i=0;i<N2;i++)
    An[i] *= 0.5;

 #endif

  //clean_small_values(N,An);

 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  //#pragma acc exit data copyout(An[0:N2])
  #pragma acc update self(An[0:N2])
 #endif


 #if 0
  for (int m=0;m<N;m++)
  {
    double val = An[m*N+m];
    if (fabs(val)<1.e-2)
    {
      printf(" WARNING: small diagonal element in A: %3i \n",m);
    }
    double v2 = 1./sqrt(val);
    //Anorm[m] = v2;
    for (int n=0;n<N;n++)
    {
      An[m*N+n] *= v2;
      An[n*N+m] *= v2;
    }
  }
 #endif


  if (prl>2)
  {
    printf("\n A: \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %12.6f",An[i*N+j]);
      printf("\n");
    }
  }

 //CPMZ check this
#if USE_ACC
  #pragma acc exit data delete(grid1[0:gs6],wt1[0:gs])
  #pragma acc exit data delete(grid2[0:gs6],wt2[0:gs])
  #pragma acc exit data delete(grid1s[0:gs6],grid2s[0:gs6])
  #pragma acc exit data delete(val1[0:iN][0:gs],val2[0:iN][0:gs],wtt1[0:gs])
  #pragma acc exit data delete(val3[0:iN][0:gs],val4[0:iN][0:gs])
  #pragma acc exit data delete(n2i[0:natoms])
#endif

  delete [] ang_g;
  delete [] ang_w;

  delete [] n2i;

  delete [] grid1s;
  delete [] grid2s;

  for (int i=0;i<iN;i++) delete [] val1[i];
  for (int i=0;i<iN;i++) delete [] val2[i];
  for (int i=0;i<iN;i++) delete [] val3[i];
  for (int i=0;i<iN;i++) delete [] val4[i];
  delete [] val1;
  delete [] val2;
  delete [] val3;
  delete [] val4;
  delete [] wtt1;

  delete [] grid1;
  delete [] grid2;
  delete [] wt1;
  delete [] wt2;

  return;
}

void compute_all_2c(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, float* A, int prl)
{
  if (prl>1) printf(" beginning compute_all_2c \n");

 //2c integrals are all in auxiliary basis
  int N = basis.size();
  int N2 = N*N;

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

  int gs = nrad*nang;
  float* grid1 = new float[6*gs];
  float* wt1 = new float[gs];
  float* val1 = new float[gs];

  float* grid2 = new float[6*gs];
  float* wt2 = new float[gs];
  float* val2 = new float[gs];

  int* i2m = new int[N];

 #if 1
  float* grid1s = new float[6*gs];
  float* grid2s = new float[6*gs];
  float** valt1 = new float*[N];
  float** valt2 = new float*[N];
  for (int i=0;i<N;i++)
    valt1[i] = new float[gs];
  for (int i=0;i<N;i++)
    valt2[i] = new float[gs];
 #else
  float* valt1 = new float[gs];
  float* valt2 = new float[gs];
 #endif
  float* wtt1 = new float[gs];

  float* ang_g = new float[3*nang];
  float* ang_w = new float[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

 #if USE_ACC
  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data create(grid1[0:6*gs],wt1[0:gs],val1[0:gs])
  #pragma acc enter data create(grid2[0:6*gs],wt2[0:gs],val2[0:gs])
  #pragma acc enter data create(valt1[0:N][0:gs],valt2[0:N][0:gs],wtt1[0:gs])
  #pragma acc enter data create(grid1s[0:6*gs],grid2s[0:6*gs])
  //#pragma acc enter data create(valt1[0:gs],valt2[0:gs],wtt1[0:gs])
  #pragma acc enter data create(A[0:N2])
  #pragma acc enter data create(i2m[0:N])
 #endif

  for (int m=0;m<natoms;m++)
  {
    float Z1 = (float)atno[m];
    float A1 = coords[3*m+0]; float B1 = coords[3*m+1]; float C1 = coords[3*m+2];

    generate_central_grid_2(grid1,wt1,Z1,nrad,nang,ang_g,ang_w);
    //#pragma acc update self(grid1[0:6*gs],wt1[0:gs])
    //print_grid(gs,grid1,NULL,wt1,NULL,prl);

    for (int i1=0;i1<N;i1++) i2m[i1] = 0;
    for (int i1=0;i1<N;i1++) if (basis[i1][9]==m) i2m[i1] = 1;
   #if USE_ACC
    #pragma acc update device(i2m[0:N])
   #endif

   //first compute single atom ints
    for (int i1=0;i1<N;i1++)
    if (basis[i1][9]==m)
    {
      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];
      //printf("\n  m: %i i1: %2i   nlm: %i %i %2i zeta: %8.5f \n",m,i1,n1,l1,m1,zeta1);

      acc_assign(gs,val1,1.);

      eval_inr_r12(gs,grid1,val1,n1,l1,zeta1);
      eval_sh_3r(gs,grid1,val1,n1,l1,m1);

      //#pragma acc update self(val1)
      //printf(" val1[%i] \n",i1);
      //print_array(gs,val1);

     #pragma acc parallel loop present(val1[0:gs],valt1[0:N][0:gs],i2m[0:N])
      for (int i2=0;i2<N;i2++)
      if (i2m[i2])
      {
       #pragma acc loop
        for (int j=0;j<gs;j++)
          valt1[i2][j] = val1[j];
      }

      for (int i2=0;i2<N;i2++)
      if (basis[i2][9]==m)
      {
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];
        //printf("    n: %i i2: %2i   nlm: %i %i %2i zeta: %8.5f   (b2) \n",m,i2,n2,l2,m2,zeta2);

        //acc_copyf(i2,gs,valt1[i2],val1);

        eval_sh(i2,gs,grid1,valt1[i2],n2,l2,m2,zeta2);

        //#pragma acc update self(valt1[i2])
        //printf(" val2[%i] \n",i2);
        //print_array(gs,valt1[i2]);

      } //loop i2 evaluate
     #if USE_ACC
      #pragma acc wait
     #endif

     #pragma acc parallel loop present(valt1[0:N][0:gs],wt1[0:gs],A[0:N2],i2m[0:N])
      for (int i2=0;i2<N;i2++)
      if (i2m[i2])
      {
        float val = 0.;
       #pragma acc loop reduction(+:val)
        for (int j=0;j<gs;j++)
          val += valt1[i2][j] * wt1[j];
        A[i1*N+i2] = val;

      } //loop i2 reduce

    } //loop i1


   //two-atom ints
    for (int n=0;n<natoms;n++)
    if (m!=n)
    {
      float Z2 = (float)atno[n];
      float A2 = coords[3*n+0]; float B2 = coords[3*n+1]; float C2 = coords[3*n+2];
      float A12 = A2-A1; float B12 = B2-B1; float C12 = C2-C1;

     //grid1 at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs,grid1,A12,B12,C12);

      generate_central_grid_2(grid2,wt2,Z2,nrad,nang,ang_g,ang_w);
      copy_grid(gs,grid2s,grid2); //grid 2 centered on atom 2
      recenter_grid(gs,grid2,A12,B12,C12); //grid 2 centered on atom 1

      copy_grid(gs,grid1s,grid1); 
      recenter_grid_zero(gs,grid1s,-A12,-B12,-C12); //grid 1 centered on atom 2

      acc_copyf(gs,wtt1,wt1);
      becke_weight_2c(gs,grid1,wtt1,grid2,wt2,Z1,Z2,A12,B12,C12);
      eliminate_small_wt(gs,wtt1);
      eliminate_small_wt(gs,wt2);

      add_r1_to_grid(gs,grid2,0.,0.,0.);

      //#pragma acc update self(grid1[0:6*gs],wtt1[0:gs],grid2[0:6*gs],wt2[0:gs])
      //print_grid(gs,grid1,grid2,wtt1,wt2,prl);

      for (int i1=0;i1<N;i1++) i2m[i1] = 0;
      for (int i1=0;i1<N;i1++) if (basis[i1][9]==n) i2m[i1] = 1;
     #if USE_ACC
      #pragma acc update device(i2m[0:N])
     #endif

      for (int i1=0;i1<N;i1++)
      if (basis[i1][9]==m)
      {
        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];
        //printf("  m: %i i1: %i   nlm: %i %i %2i zeta: %8.5f \n",m,i1,n1,l1,m1,zeta1);

        acc_assign(gs,val1,1.);
        acc_assign(gs,val2,1.);

        eval_inr_r12(gs,grid1,val1,n1,l1,zeta1);
        eval_sh_3r(gs,grid1,val1,n1,l1,m1);
        eval_inr_r12(gs,grid2,val2,n1,l1,zeta1);
        eval_sh_3r(gs,grid2,val2,n1,l1,m1);

       #pragma acc parallel loop present(val1[0:gs],valt1[0:N][0:gs],val2[0:gs],valt2[0:N][0:gs],i2m[0:N])
        for (int i2=0;i2<N;i2++)
        if (i2m[i2])
        {
         #pragma acc loop
          for (int j=0;j<gs;j++)
            valt1[i2][j] = val1[j];
         #pragma acc loop
          for (int j=0;j<gs;j++)
            valt2[i2][j] = val2[j];
        }

       //launch //async processes over i2
        for (int i2=i1+1;i2<N;i2++)
        if (basis[i2][9]==n)
        {
          vector<double> basis2 = basis[i2];
          int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];
          //printf("   n: %i i2: %i   nlm: %i %i %2i zeta: %8.5f \n",n,i2,n2,l2,m2,zeta2);

          //acc_copyf(i2,gs,valt2[i2],val2);
          //acc_copyf(i2,gs,valt1[i2],val1);

          eval_sh(i2,gs,grid2s,valt2[i2],n2,l2,m2,zeta2);
          eval_sh(i2,gs,grid1s,valt1[i2],n2,l2,m2,zeta2);
        }
       #if USE_ACC
        #pragma acc wait
       #endif

       #pragma acc parallel loop present(valt1[0:N][0:gs],wtt1[0:gs],valt2[0:N][0:gs],wt2[0:gs],A[0:N2],i2m[0:N])
        for (int i2=i1+1;i2<N;i2++)
        if (i2m[i2])
        {
          float val = 0.;
         #pragma acc loop reduction(+:val)
          for (int j=0;j<gs;j++)
            val += valt1[i2][j] * wtt1[j];
         #pragma acc loop reduction(+:val)
          for (int j=0;j<gs;j++)
            val += valt2[i2][j] * wt2[j];

          A[i1*N+i2] = val;

        } //loop i2 reduce

      } //loop i1 over first basis function

    } //loop n over second atom

  } //loop m over natoms

 #if USE_ACC
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc exit data copyout(A[0:N2])
 #endif

#if DEBUG
  if (prl>0)
  {
    printf("\n A(raw): \n");
    for (int i=0;i<N;i++)
    {
      for (int j=0;j<N;j++)
        printf(" %8.5f",A[i*N+j]);
      printf("\n");
    }
  }
#endif

  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    A[i*N+j] *= norm_sv(basis[i][0],basis[i][1],basis[i][2],basis[i][3])*basis[j][4];

  for (int i=0;i<N;i++)
  for (int j=i+1;j<N;j++)
    A[j*N+i] = A[i*N+j];

  if (prl>2)
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
  #pragma acc exit data delete(grid1[0:6*gs],wt1[0:gs],val1[0:gs])
  #pragma acc exit data delete(grid2[0:6*gs],wt2[0:gs],val2[0:gs])
  #pragma acc exit data delete(grid1s[0:6*gs],grid2s[0:6*gs])
  #pragma acc exit data delete(valt1[0:N][0:gs],valt2[0:N][0:gs],wtt1[0:gs])
  #pragma acc exit data delete(A[0:N2])
  #pragma acc exit data delete(i2m[0:N])
#endif

  delete [] ang_g;
  delete [] ang_w;

  delete [] i2m;

 #if 1
  for (int i=0;i<N;i++) delete [] valt1[i];
  for (int i=0;i<N;i++) delete [] valt2[i];
 #endif
  delete [] valt1;
  delete [] valt2;
  delete [] wtt1;

  delete [] grid1;
  delete [] grid2;
  delete [] wt1;
  delete [] wt2;
  delete [] val1;
  delete [] val2;

  return;
}


#if RED_DOUBLE
void compute_all_4c_v2(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* g, int prl)
#else
void compute_all_4c_v2(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, float* g, int prl)
#endif
{
  if (prl>-1) printf(" beginning compute_all_4c_v2 \n");
  if (natoms>2)
  {
    printf("\n ERROR: compute_all_4c for 2-atom integrals only \n");
  }

  int N = basis.size();
  int N2 = N*N;
  int N3 = N2*N;

  int gs = nrad*nang;
  int gsp = (nrad+1)*nang;
  int gs6 = 6*gs;
  int gsp6 = 6*gsp;

  if (prl>0) printf("  gs/p: %4i %4i \n",gs,gsp);

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

  int estart = find_center_of_grid(1,nrad)*nang;
  int estart2 = find_center_of_grid(1,nrad+1)*nang;

  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);

  float** grid1 = new float*[4]; for (int i=0;i<4;i++) grid1[i] = new float[gs6];
  float** wt1 = new float*[4]; for (int i=0;i<4;i++) wt1[i] = new float[gs];

  float** grid2 = new float*[4]; for (int i=0;i<4;i++) grid2[i] = new float[gsp6];
  float** wt2 = new float*[4]; for (int i=0;i<4;i++) wt2[i] = new float[gsp];

  float*** val1 = new float**[4];
  float*** val2 = new float**[4];
  float*** val3 = new float**[4];
  float*** val4 = new float**[4];
  for (int n=0;n<4;n++)
  {
    val1[n] = new float*[iN];
    for (int i=0;i<iN;i++) val1[n][i] = new float[gs];
  }
  for (int n=0;n<4;n++)
  {
    val2[n] = new float*[iN];
    for (int i=0;i<iN;i++) val2[n][i] = new float[gs];
  }
  for (int n=0;n<4;n++)
  {
    val3[n] = new float*[iN];
    for (int i=0;i<iN;i++) val3[n][i] = new float[gsp];
  }
  for (int n=0;n<4;n++)
  {
    val4[n] = new float*[iN];
    for (int i=0;i<iN;i++) val4[n][i] = new float[gsp];
  }

  float** wtt = new float*[16]; for (int i=0;i<16;i++) wtt[i] = new float[gsp];

  float* ang_g = new float[3*nang];
  float* ang_w = new float[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

  for (int i=0;i<N2*N2;i++)
    g[i] = 0.;

 //g not on gpu (too large)
  int M = iN;
  int M2 = M*M;
  //int M3 = M2*M;
  float* gt = new float[M2*M2];

 #if USE_ACC
  #pragma acc enter data create(gt[0:M2*M2])

  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data copyin(n2i[0:natoms])

  #pragma acc enter data create(grid1[0:4][0:gs6],wt1[0:4][0:gs])
  #pragma acc enter data create(grid2[0:4][0:gsp6],wt2[0:4][0:gsp])
  #pragma acc enter data create(wtt[0:16][0:gsp])

  for (int n=0;n<4;n++)
  {
    float** val1n = val1[n]; float** val2n = val2[n]; float** val3n = val3[n]; float** val4n = val4[n];
    #pragma acc enter data create(val1n[0:iN][0:gs],val2n[0:iN][0:gs],val3n[0:iN][0:gsp],val4n[0:iN][0:gsp])
  }
 #endif

  #pragma acc parallel loop present(gt[0:M2*M2])
  for (int j=0;j<M2*M2;j++)
    gt[j] = 0.f;

  for (int m=0;m<natoms;m++)
  {
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];
    float Z1 = (float)atno[m];
    float A1 = coords[3*m+0]; float B1 = coords[3*m+1]; float C1 = coords[3*m+2];

    float* grid1a = grid1[0]; float* grid2a = grid2[0];
    float* wt1a = wt1[0]; float* wt2a = wt2[0];

    generate_central_grid_2(grid1a,wt1a,Z1,nrad,  nang,ang_g,ang_w);
    generate_central_grid_2(grid2a,wt2a,Z1,nrad+1,nang,ang_g,ang_w);

    float** val1a = val1[0]; float** val2a = val2[0];
    float** val3a = val3[0]; float** val4a = val4[0];

   #pragma acc parallel loop collapse(2) present(val1a[0:iN][0:gs])
    for (int ii1=0;ii1<s2-s1;ii1++)
    for (int j=0;j<gs;j++)
      val1a[ii1][j] = 1.f;
   #pragma acc parallel loop collapse(2) present(val2a[0:iN][0:gs])
    for (int ii2=0;ii2<s2-s1;ii2++)
    for (int j=0;j<gs;j++)
      val2a[ii2][j] = 1.f;
   #pragma acc parallel loop collapse(2) present(val3a[0:iN][0:gsp])
    for (int ii3=0;ii3<s2-s1;ii3++)
    for (int j=0;j<gsp;j++)
      val3a[ii3][j] = 1.f;
   #pragma acc parallel loop collapse(2) present(val4a[0:iN][0:gsp])
    for (int ii4=0;ii4<s2-s1;ii4++)
    for (int j=0;j<gsp;j++)
      val4a[ii4][j] = 1.f;

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

      eval_sh(ii1,gs,grid1a,val1a[ii1],n1,l1,m1,zeta1);
    } //loop i1

    for (int i2=s1;i2<s2;i2++)
    {
      int ii2 = i2-s1;

      vector<double> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

      eval_sh(ii2,gs,grid1a,val2a[ii2],n2,l2,m2,zeta2);
    } //loop i2

    for (int i3=s1;i3<s2;i3++)
    {
      int ii3 = i3-s1;

      vector<double> basis3 = basis[i3];
      int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

      eval_sh(ii3,gsp,grid2a,val3a[ii3],n3,l3,m3,zeta3);
    } //loop i3

    for (int i4=s1;i4<s2;i4++)
    {
      int ii4 = i4-s1;

      vector<double> basis4 = basis[i4];
      int n4 = basis4[0]; int l4 = basis4[1]; int m4 = basis4[2]; double zeta4 = basis4[3];

      eval_sh(ii4,gsp,grid2a,val4a[ii4],n4,l4,m4,zeta4);
    } //loop i4

   //form V2a for all basis on atom m
   //contract density with potential
//    reduce_4cv_1(s1,s2,gs,gsp,grid1a,grid2a,val1a,val2a,V2a,wt1a);
//    reduce_4vv_1(s1,s2,gsp,grid2a,V2a,val3a,val4a,wt2a,gt);

    //reduce_4c_1(s1,s2,s1,s2,gs,gsp,M,grid1a,grid2a,val1a,val2a,val3a,val4a,wt1a,wt2a,gt);
    //collect_4c_1(s1,s2,s1,s2,0,gs,gsp,M,N,gt,g);


   //two-atom ints
    for (int n=0;n<m;n++)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      float Z2 = (float)atno[n];
      float A2 = coords[3*n+0]; float B2 = coords[3*n+1]; float C2 = coords[3*n+2];
      float A12 = A2-A1; float B12 = B2-B1; float C12 = C2-C1;

      float* grid1b = grid1[1]; float* grid2b = grid2[1];
      float* grid1d = grid1[3]; float* grid2d = grid2[3];
      float* wt1b = wt1[1]; float* wt2b = wt2[1];

      float** val1b = val1[1]; float** val2b = val2[1];
      float** val3b = val3[1]; float** val4b = val4[1];

     //grid1a at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs, grid1a,A12,B12,C12);
      add_r2_to_grid(gsp,grid2a,A12,B12,C12);

      generate_central_grid_2(grid1b,wt1b,Z2,nrad,  nang,ang_g,ang_w);
      generate_central_grid_2(grid2b,wt2b,Z2,nrad+1,nang,ang_g,ang_w);

      copy_grid(gs, grid1d,grid1b); //centered at atom b
      copy_grid(gsp,grid2d,grid2b); //centered at atom b

      recenter_grid(gs, grid1b,A12,B12,C12); //r1 points to self, r2 points at first atom
      recenter_grid(gsp,grid2b,A12,B12,C12);

      float* wtt1a = wtt[0]; float* wtt1b = wtt[1];
      float* wtt2a = wtt[4]; float* wtt2b = wtt[5];
      acc_copyf(gs, wtt1a,wt1a); acc_copyf(gs, wtt1b,wt1b);
      acc_copyf(gsp,wtt2a,wt2a); acc_copyf(gsp,wtt2b,wt2b);
      becke_weight_2c(gs, grid1a,wtt1a,grid1b,wtt1b,Z1,Z2,A12,B12,C12);
      becke_weight_2c(gsp,grid2a,wtt2a,grid2b,wtt2b,Z1,Z2,A12,B12,C12);

      //eliminate_small_wt(estart, gs, wtt1a);
      //eliminate_small_wt(estart, gs, wtt1b);
      //eliminate_small_wt(estart2,gsp,wtt2a);
      //eliminate_small_wt(estart2,gsp,wtt2b);

      add_r1_to_grid(gs, grid1b,0.,0.,0.); //center is on atom a
      add_r1_to_grid(gsp,grid2b,0.,0.,0.); //center is on atom a

      float* grid1c = grid1[2]; float* grid2c = grid2[2];

      copy_grid(gs, grid1c,grid1a);
      copy_grid(gsp,grid2c,grid2a);

      recenter_grid_zero(gs, grid1c,-A12,-B12,-C12);
      recenter_grid_zero(gsp,grid2c,-A12,-B12,-C12);


    //1122 (s12,s12,s34,s34)
     #pragma acc parallel loop collapse(2) present(val1b[0:iN][0:gs])
      for (int ii1=0;ii1<s2-s1;ii1++)
      for (int j=0;j<gs;j++)
        val1b[ii1][j] = 1.f;

     #pragma acc parallel loop collapse(2) present(val2b[0:iN][0:gs])
      for (int ii2=0;ii2<s2-s1;ii2++)
      for (int j=0;j<gs;j++)
        val2b[ii2][j] = 1.f;

     #pragma acc parallel loop collapse(2) present(val3b[0:iN][0:gsp])
      for (int ii3=0;ii3<s4-s3;ii3++)
      for (int j=0;j<gsp;j++)
        val3b[ii3][j] = 1.f;

     #pragma acc parallel loop collapse(2) present(val4b[0:iN][0:gsp])
      for (int ii4=0;ii4<s4-s3;ii4++)
      for (int j=0;j<gsp;j++)
        val4b[ii4][j] = 1.f;

     //i1 on atom m
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;
        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_sh(ii1,gs,grid1b,val1b[ii1],n1,l1,m1,zeta1);
      }

     //i2 on atom m
      for (int i2=s1;i2<s2;i2++)
      {
        int ii2 = i2-s1;
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1b,val2b[ii2],n2,l2,m2,zeta2);
      }

     //i3 on atom n
      for (int i3=s3;i3<s4;i3++)
      {
        int ii3 = i3-s3;
        vector<double> basis3 = basis[i3];
        int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];
    
        eval_sh(ii3,gsp,grid2d,val3b[ii3],n3,l3,m3,zeta3);
      }

     //i4 on atom n
      for (int i4=s3;i4<s4;i4++)
      {
        int ii4 = i4-s3;
        vector<double> basis4 = basis[i4];
        int n4 = basis4[0]; int l4 = basis4[1]; int m4 = basis4[2]; double zeta4 = basis4[3];
    
        eval_sh(ii4,gsp,grid2d,val4b[ii4],n4,l4,m4,zeta4);
      }

     #pragma acc parallel loop collapse(2) present(val3a[0:iN][0:gsp])
      for (int ii3=0;ii3<s4-s3;ii3++)
      for (int j=0;j<gsp;j++)
        val3a[ii3][j] = 1.f;

     #pragma acc parallel loop collapse(2) present(val4a[0:iN][0:gsp])
      for (int ii4=0;ii4<s4-s3;ii4++)
      for (int j=0;j<gsp;j++)
        val4a[ii4][j] = 1.f;

     //i3 on atom n
      for (int i3=s3;i3<s4;i3++)
      {
        int ii3 = i3-s3;
        vector<double> basis3 = basis[i3];
        int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

        eval_sh(ii3,gsp,grid2c,val3a[ii3],n3,l3,m3,zeta3);
      }

     //i4 on atom n
      for (int i4=s3;i4<s4;i4++)
      {
        int ii4 = i4-s3;
        vector<double> basis4 = basis[i4];
        int n4 = basis4[0]; int l4 = basis4[1]; int m4 = basis4[2]; double zeta4 = basis4[3];

        eval_sh(ii4,gsp,grid2c,val4a[ii4],n4,l4,m4,zeta4);
      }

     //1122
//    reduce_4cv_1(s1,s2,gs,gsp,grid1a,grid2c,val1a,val2a,V2c,wt1a);
//    reduce_4vv_4(s1,s2,gsp,grid2a,grid2c,V2a,V2c,val3a,val4a,val3b,val4b,wt2a,gt);

     //1122
      //reduce_4c_1(s1,s2,s3,s4,gs,gsp,M,grid1a,grid2a,val1a,val2a,val3a,val4a,wtt1a,wtt2a,gt);
      //reduce_4c_1(s1,s2,s3,s4,gs,gsp,M,grid1a,grid2b,val1a,val2a,val3b,val4b,wtt1a,wtt2b,gt);
      //reduce_4c_1(s1,s2,s3,s4,gs,gsp,M,grid1b,grid2a,val1b,val2b,val3a,val4a,wtt1b,wtt2a,gt);
      //reduce_4c_1(s1,s2,s3,s4,gs,gsp,M,grid1b,grid2b,val1b,val2b,val3b,val4b,wtt1b,wtt2b,gt);
      //collect_4c_1(s1,s2,s3,s4,1,gs,gsp,M,N,gt,g);



    //1121
     #pragma acc parallel loop collapse(2) present(val4a[0:iN][0:gsp])
      for (int ii4=0;ii4<s2-s1;ii4++)
      for (int j=0;j<gsp;j++)
        val4a[ii4][j] = 1.f;

     //i4 on atom m
      for (int i4=s1;i4<s2;i4++)
      {
        int ii4 = i4-s1;
        vector<double> basis4 = basis[i4];
        int n4 = basis4[0]; int l4 = basis4[1]; int m4 = basis4[2]; double zeta4 = basis4[3];

        eval_sh(ii4,gsp,grid2a,val4a[ii4],n4,l4,m4,zeta4);
      }

     #pragma acc parallel loop collapse(2) present(val4b[0:iN][0:gsp])
      for (int ii4=0;ii4<s2-s1;ii4++)
      for (int j=0;j<gsp;j++)
        val4b[ii4][j] = 1.f;

     //i2 on atom m
      for (int i4=s1;i4<s2;i4++)
      {
        int ii4 = i4-s1;
        vector<double> basis4 = basis[i4];
        int n4 = basis4[0]; int l4 = basis4[1]; int m4 = basis4[2]; double zeta4 = basis4[3];

        eval_sh(ii4,gsp,grid2b,val4b[ii4],n4,l4,m4,zeta4);
      }

     //1121
      reduce_4c_1b(s1,s2,s3,s4,gs,gsp,M,grid1a,grid2a,val1a,val2a,val3a,val4a,wtt1a,wtt2a,gt);
      reduce_4c_1b(s1,s2,s3,s4,gs,gsp,M,grid1a,grid2b,val1a,val2a,val3b,val4b,wtt1a,wtt2b,gt);
      reduce_4c_1b(s1,s2,s3,s4,gs,gsp,M,grid1b,grid2a,val1b,val2b,val3a,val4a,wtt1b,wtt2a,gt);
      reduce_4c_1b(s1,s2,s3,s4,gs,gsp,M,grid1b,grid2b,val1b,val2b,val3b,val4b,wtt1b,wtt2b,gt);
      collect_4c_1b(s1,s2,s3,s4,gs,gsp,M,N,gt,g);



    //1221
     #pragma acc parallel loop collapse(2) present(val2a[0:iN][0:gs])
      for (int ii2=0;ii2<s4-s3;ii2++)
      for (int j=0;j<gs;j++)
        val2a[ii2][j] = 1.f;

     //i2 on atom n
      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1c,val2a[ii2],n2,l2,m2,zeta2);
      }

     #pragma acc parallel loop collapse(2) present(val2b[0:iN][0:gs])
      for (int ii2=0;ii2<s4-s3;ii2++)
      for (int j=0;j<gs;j++)
        val2b[ii2][j] = 1.f;

     //i2 on atom n
      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1d,val2b[ii2],n2,l2,m2,zeta2);
      }

     //1221
      reduce_4c_1c(s1,s2,s3,s4,gs,gsp,M,grid1a,grid2a,val1a,val2a,val3a,val4a,wtt1a,wtt2a,gt);
      reduce_4c_1c(s1,s2,s3,s4,gs,gsp,M,grid1a,grid2b,val1a,val2a,val3b,val4b,wtt1a,wtt2b,gt);
      reduce_4c_1c(s1,s2,s3,s4,gs,gsp,M,grid1b,grid2a,val1b,val2b,val3a,val4a,wtt1b,wtt2a,gt);
      reduce_4c_1c(s1,s2,s3,s4,gs,gsp,M,grid1b,grid2b,val1b,val2b,val3b,val4b,wtt1b,wtt2b,gt);
      collect_4c_1c(s1,s2,s3,s4,gs,gsp,M,N,gt,g);


    //2221
     #pragma acc parallel loop collapse(2) present(val1a[0:iN][0:gs])
      for (int ii1=0;ii1<s4-s3;ii1++)
      for (int j=0;j<gs;j++)
        val1a[ii1][j] = 1.f;

     //i1 on atom n
      for (int i1=s3;i1<s4;i1++)
      {
        int ii1 = i1-s3;
        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_sh(ii1,gs,grid1c,val1a[ii1],n1,l1,m1,zeta1);
      }

     #pragma acc parallel loop collapse(2) present(val1b[0:iN][0:gs])
      for (int ii1=0;ii1<s4-s3;ii1++)
      for (int j=0;j<gs;j++)
        val1b[ii1][j] = 1.f;

     //i1 on atom n
      for (int i1=s3;i1<s4;i1++)
      {
        int ii1 = i1-s3;
        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_sh(ii1,gs,grid1d,val1b[ii1],n1,l1,m1,zeta1);
      }

      reduce_4c_1d(s1,s2,s3,s4,gs,gsp,M,grid1a,grid2a,val1a,val2a,val3a,val4a,wtt1a,wtt2a,gt);
      reduce_4c_1d(s1,s2,s3,s4,gs,gsp,M,grid1a,grid2b,val1a,val2a,val3b,val4b,wtt1a,wtt2b,gt);
      reduce_4c_1d(s1,s2,s3,s4,gs,gsp,M,grid1b,grid2a,val1b,val2b,val3a,val4a,wtt1b,wtt2a,gt);
      reduce_4c_1d(s1,s2,s3,s4,gs,gsp,M,grid1b,grid2b,val1b,val2b,val3b,val4b,wtt1b,wtt2b,gt);
      collect_4c_1d(s1,s2,s3,s4,gs,gsp,M,N,gt,g);


    } //loop n over second atom

  } //loop m over natoms

 //normalization
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
  for (int l=0;l<N;l++)
    g[i*N3+j*N2+k*N+l] *= basis[i][4]*basis[j][4]*basis[k][4]*basis[l][4];

 //show the symmetry
  if (prl>1)
  for (int p=0;p<N;p++)
  for (int q=0;q<=p;q++)
  for (int r=0;r<=p;r++)
  for (int s=0;s<=r;s++)
  {
    double v1 = g[p*N3+q*N2+r*N+s];
    double v2 = g[p*N3+q*N2+s*N+r];
    double v3 = g[q*N3+p*N2+r*N+s];
    double v4 = g[q*N3+p*N2+s*N+r];

    double v5 = g[r*N3+s*N2+p*N+q];
    double v6 = g[s*N3+r*N2+p*N+q];
    double v7 = g[r*N3+s*N2+q*N+p];
    double v8 = g[s*N3+r*N2+q*N+p];

    printf("   pqrs: %i %i %i %i: %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f \n",p,q,r,s,v1,v2,v3,v4,v5,v6,v7,v8);
  }

  if ((prl>0 && N<10) || prl>1)
  {
    printf("\n g: \n");
    for (int i=0;i<N;i++)
    for (int j=0;j<N;j++)
    {
      printf(" ij: %2i %2i \n",i,j);
      print_square(N,&g[i*N3+j*N2]);
    }
  }


#if USE_ACC
  #pragma acc exit data delete(gt[0:M2*M2])

  #pragma acc exit data delete(grid1[0:4][0:gs6],wt1[0:4][0:gs])
  #pragma acc exit data delete(grid2[0:4][0:gsp6],wt2[0:4][0:gsp])
  //#pragma acc exit data delete(grid3[0:4][0:gsp6],wt3[0:4][0:gsp])
  //#pragma acc exit data delete(grid4[0:4][0:gsp6],wt4[0:4][0:gsp])
  for (int n=0;n<4;n++)
  {
    float** val1n = val1[n]; float** val2n = val2[n]; float** val3n = val3[n]; float** val4n = val4[n];
    #pragma acc exit data delete(val1n[0:iN][0:gs],val2n[0:iN][0:gs],val3n[0:iN][0:gsp],val4n[0:iN][0:gsp])
  }

  #pragma acc exit data delete(wtt[0:16][0:gsp])
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc exit data delete(n2i[0:natoms])
#endif


  delete [] gt;
  delete [] n2i;

  delete [] ang_g;
  delete [] ang_w;

  for (int i=0;i<4;i++) 
  {
    for (int j=0;j<iN;j++)
      delete [] val1[i][j];
    delete [] val1[i];
    for (int j=0;j<iN;j++)
      delete [] val2[i][j];
    delete [] val2[i];
    for (int j=0;j<iN;j++)
      delete [] val3[i][j];
    delete [] val3[i];
    for (int j=0;j<iN;j++)
      delete [] val4[i][j];
    delete [] val4[i];
  }
 
  for (int i=0;i<4;i++)
  {
    delete [] grid1[i];
    delete [] grid2[i];
    delete [] wt1[i];
    delete [] wt2[i];
  }
  delete [] grid1;
  delete [] grid2;
  delete [] wt1;
  delete [] wt2;

  return;
}


//compares in the equality sense
bool compare_ijkl(const vector<short>& ijkl1, const vector<short>& ijkl2)
{
  short i1 = ijkl1[0]; short i2 = ijkl1[1];
  short i3 = ijkl1[2]; short i4 = ijkl1[3];

  short j1 = ijkl2[0]; short j2 = ijkl2[1];
  short j3 = ijkl2[2]; short j4 = ijkl2[3];

  bool found = false;
  if (i1==j1 && i2==j2)
  {
    if (i3==j3 && i4==j4)
      found = true;
    else if (i3==j4 && i4==j3)
      found = true;
  }
  else if (i1==j2 && i2==j1)
  {
    if (i3==j3 && i4==j4)
      found = true;
    else if (i3==j4 && i4==j3)
      found = true;
  }
  else if (i1==j3 && i2==j4)
  {
    if (i3==j1 && i4==j2)
      found = true;
    else if (i3==j2 && i4==j1)
      found = true;
  }
  else if (i1==j4 && i2==j3)
  {
    if (i3==j1 && i4==j2)
      found = true;
    else if (i3==j2 && i4==j1)
      found = true;
  }

  return found;
}

//compares in the less than sense
bool compare_ijkl_2(const vector<short>& ijkl1, const vector<short>& ijkl2)
{
  short i1 = ijkl1[0]; short i2 = ijkl1[1];
  short i3 = ijkl1[2]; short i4 = ijkl1[3];

 //put in canonical order first
  if (i1<i2) { short t1 = i1; i1 = i2; i2 = t1; } 
  if (i3<i4) { short t3 = i3; i3 = i4; i4 = t3; } 
  if (i1<i3) { short t1 = i1; short t2 = i2; i1 = i3; i2 = i4; i3 = t1; i4 = t2; }

  short j1 = ijkl2[0]; short j2 = ijkl2[1];
  short j3 = ijkl2[2]; short j4 = ijkl2[3];

  if (j1<j2) { short t1 = j1; j1 = j2; j2 = t1; } 
  if (j3<j4) { short t3 = j3; j3 = j4; j4 = t3; } 
  if (j1<j3) { short t1 = j1; short t2 = j2; j1 = j3; j2 = j4; j3 = t1; j4 = t2; }

  bool one_lt_two = false;
  if (i1==j1)
  {
    if (i3==j3)
    {
      if (i2==j2)
      {
        if (i4<j4)
          one_lt_two = true;
      }
      else if (i2<j2)
        one_lt_two = true;
    }
    else if (i3<j3)
      one_lt_two = true;
  }
  else if (i1<j1)
    one_lt_two = true;

  return one_lt_two;
}

void remove_duplicate_ijkl(vector<vector<short> >& ijkl0)
{
 #if 0
  printf("  original ijkl: \n");
  for (int i=0;i<ijkl0.size();i++)
    printf("   %2i %2i %2i %2i \n",ijkl0[i][0],ijkl0[i][1],ijkl0[i][2],ijkl0[i][3]);
 #endif

  vector<vector<short> > ijkl1 = ijkl0;

  sort(ijkl1.begin(),ijkl1.end(),compare_ijkl_2);

  int ijkl1s = ijkl1.size();
  bool list[ijkl1s];
  for (int i=0;i<ijkl1s;i++)
    list[i] = false;

  for (int i=0;i<ijkl1s-1;i++)
  {
    bool same = compare_ijkl(ijkl1[i],ijkl1[i+1]);
    if (same) list[i+1] = true;
  }

  ijkl0.clear();
  for (int i=0;i<ijkl1s;i++)
  if (!list[i])
    ijkl0.push_back(ijkl1[i]);

 #if 0
  printf("  unique ijkl: \n");
  for (int i=0;i<ijkl0.size();i++)
    printf("   %2i %2i %2i %2i \n",ijkl0[i][0],ijkl0[i][1],ijkl0[i][2],ijkl0[i][3]);
 #endif

  return;
}

vector<vector<short> > enumerate_ijkl(int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8)
{
  vector<vector<short> > ijkl;

  if (s2==s4 && s2==s6 && s2==s8)
  {
   //directly generating unique list for single-center case
    for (int i1=s1;i1<s2;i1++)
    {
      for (int i2=s1;i2<=i1;i2++)
      for (int i3=s1;i3<=i1;i3++)
      for (int i4=s1;i4<=i3;i4++)
      {
        vector<short> list1;
        list1.push_back(i1); list1.push_back(i2);
        list1.push_back(i3); list1.push_back(i4);
        ijkl.push_back(list1);
      }
    }
  }
  else
  {
    for (int i1=s1;i1<s2;i1++)
    {
      for (int i2=s3;i2<s4;i2++)
      for (int i3=s5;i3<s6;i3++)
      for (int i4=s7;i4<s8;i4++)
      {
        vector<short> list1;
        list1.push_back(i1); list1.push_back(i2);
        list1.push_back(i3); list1.push_back(i4);
        ijkl.push_back(list1);
      }
    }
    remove_duplicate_ijkl(ijkl);
  }


  return ijkl;
}

void clean_4c_grid(int gsa, int gspa, float* grid1, float* wt1, float* grid2, float* wt2)
{
  int gsa6 = 6*gsa;
  int gspa6 = 6*gspa;

  float MIN_DIST = 1.e-4;

 #pragma acc parallel loop present(grid1[0:gsa6],wt1[0:gsa],grid2[0:gspa6],wt2[0:gspa])
  for (int i=0;i<gsa;i++)
  {
    float A1 = grid1[6*i+0]; float B1 = grid1[6*i+1]; float C1 = grid1[6*i+2];
    float wt1i = wt1[i];

   #pragma acc loop
    for (int j=0;j<gspa;j++)
    {
      float A2 = grid2[6*j+0]; float B2 = grid2[6*j+1]; float C2 = grid2[6*j+2];
      float A12 = A1-A2; float B12 = B1-B2; float C12 = C1-C2;
      float r12 = sqrtf(A12*A12+B12*B12+C12*C12);
      if (r12<MIN_DIST)
      {
        if (wt1i<wt2[j])
          wt1[i] = 0.;
        else
          wt2[j] = 0.;
      }
    }
  }

  return;
}

void compute_4c_ol_oned(int ngpu, int natoms, int* atno, float* coordsf, int nrad, int nang, double* ang_g, double* ang_w, vector<vector<double> >& basis, int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8, int iN,
   int gsa, double* grid, double* wt, double** val1, double** val2, double** val3, double** val4, double* ol)
{
 //this function needs much better batching,
 // parallel utilization is poor

  int nomp = ngpu;
  int prl = 1;

  int gc = 6;
  int gsa6 = gsa*6;

  int N = basis.size();
  int N2 = N*N;
  int N3 = N2*N;

  double coords[3*natoms];
  for (int n=0;n<3*natoms;n++) coords[n] = coordsf[n];

  double minval = 1.e-16;

 //list out integrals that need to be computed
  vector<vector<short> > ijkl = enumerate_ijkl(s1,s2,s3,s4,s5,s6,s7,s8);
  int ijkl_size = ijkl.size();

  acc_set_device_num(0,acc_device_nvidia);

 //integral-specific grids
  get_becke_grid_full(natoms,atno,coords,nrad,nang,ang_g,ang_w,gc,grid,wt);

  #pragma acc update self(grid[0:gsa6],wt[0:gsa])

 #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=1;n<nomp;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc update device(grid[0:gsa6],wt[0:gsa])
  }
  acc_set_device_num(0,acc_device_nvidia);

 #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    //int tid = omp_get_thread_num();
    acc_set_device_num(n,acc_device_nvidia);

   #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gsa])
    for (int ii1=0;ii1<s2-s1;ii1++)
    for (int j=0;j<gsa;j++)
      val1[ii1][j] = 1.f;
   #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gsa])
    for (int ii2=0;ii2<s4-s3;ii2++)
    for (int j=0;j<gsa;j++)
      val2[ii2][j] = 1.f;
   #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gsa])
    for (int ii3=0;ii3<s6-s5;ii3++)
    for (int j=0;j<gsa;j++)
      val3[ii3][j] = 1.f;
   #pragma acc parallel loop collapse(2) present(val4[0:iN][0:gsa])
    for (int ii4=0;ii4<s8-s7;ii4++)
    for (int j=0;j<gsa;j++)
      val4[ii4][j] = 1.f;
  }

 #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    //int tid = omp_get_thread_num();
    acc_set_device_num(n,acc_device_nvidia);

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];
      float A1 = basis1[5]; float B1 = basis1[6]; float C1 = basis1[7];

      recenter_grid_zero(gsa,grid,-A1,-B1,-C1);
      eval_shd(ii1,gsa,grid,val1[ii1],n1,l1,m1,zeta1);
      recenter_grid_zero(gsa,grid,A1,B1,C1);
    } //loop i1

    for (int i2=s3;i2<s4;i2++)
    {
      int ii2 = i2-s3;

      vector<double> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];
      float A1 = basis2[5]; float B1 = basis2[6]; float C1 = basis2[7];

      recenter_grid_zero(gsa,grid,-A1,-B1,-C1);
      eval_shd(ii2,gsa,grid,val2[ii2],n2,l2,m2,zeta2);
      recenter_grid_zero(gsa,grid,A1,B1,C1);
    } //loop i2

    for (int i3=s5;i3<s6;i3++)
    {
      int ii3 = i3-s5;

      vector<double> basis3 = basis[i3];
      int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];
      float A1 = basis3[5]; float B1 = basis3[6]; float C1 = basis3[7];

      recenter_grid_zero(gsa,grid,-A1,-B1,-C1);
      eval_shd(ii3,gsa,grid,val3[ii3],n3,l3,m3,zeta3);
      recenter_grid_zero(gsa,grid,A1,B1,C1);
    } //loop i3

    for (int i4=s7;i4<s8;i4++)
    {
      int ii4 = i4-s7;

      vector<double> basis4 = basis[i4];
      int n4 = basis4[0]; int l4 = basis4[1]; int m4 = basis4[2]; double zeta4 = basis4[3];
      float A1 = basis4[5]; float B1 = basis4[6]; float C1 = basis4[7];

      recenter_grid_zero(gsa,grid,-A1,-B1,-C1);
      eval_shd(ii4,gsa,grid,val4[ii4],n4,l4,m4,zeta4);
      recenter_grid_zero(gsa,grid,A1,B1,C1);
    } //loop i4
  }
  //acc_set_device_num(0,acc_device_nvidia);

 #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int i0=0;i0<ijkl_size;i0++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);

    int i1 = ijkl[i0][0]; int i2 = ijkl[i0][1];
    int i3 = ijkl[i0][2]; int i4 = ijkl[i0][3];
    int ii1 = i1-s1; int ii2 = i2-s3;
    int ii3 = i3-s5; int ii4 = i4-s7;

    double v1 = 0.;
   #pragma acc parallel loop present(grid[0:gsa6],val1[0:iN][0:gsa],val2[0:iN][0:gsa],val3[0:iN][0:gsa],val4[0:iN][0:gsa],wt[0:gsa]) reduction(+:v1)
    for (int j=0;j<gsa;j++)
      v1 += val1[ii1][j]*val2[ii2][j]*val3[ii3][j]*val4[ii4][j]*wt[j];

    if (fabs(v1)<minval) v1 = 0.;
  
    ol[i1*N3+i2*N2+i3*N+i4] = v1;
    ol[i1*N3+i2*N2+i4*N+i3] = v1;
    ol[i2*N3+i1*N2+i3*N+i4] = v1;
    ol[i2*N3+i1*N2+i4*N+i3] = v1;

    ol[i3*N3+i4*N2+i1*N+i2] = v1;
    ol[i4*N3+i3*N2+i1*N+i2] = v1;
    ol[i3*N3+i4*N2+i2*N+i1] = v1;
    ol[i4*N3+i3*N2+i2*N+i1] = v1;

    if (prl>1) printf("  v1(%i %i %i %i): %8.5f \n",i1,i2,i3,i4,v1);

  } //loop i0 over ijkl

  return;
}

void compute_all_4c_ol_gend(int ngpu, int natoms, int* atno, float* coordsf, vector<vector<double> > &basis, int nrad, int nang, double* ang_g, double* ang_w, double* ol, int prl)
{
  if (prl>-1) printf(" beginning compute_all_4c_ol_gend (ngpu: %i) \n",ngpu);

  int nomp = ngpu;

  int N = basis.size();
  int N2 = N*N;
  int N3 = N2*N;

  int gs = nrad*nang;

 //handle dummy atoms with no basis ftns
  natoms = get_natoms_with_basis(natoms,atno,basis);

 //grid will contain up to 4 nuclei
  int natomsa = min(4,natoms);
  int gsa = natomsa*gs;
  int gsa6 = 6*gsa;

  if (prl>0) printf("  gs: %4i  natomsa: %i \n",gs,natomsa);

  //int estart = find_center_of_grid(1,nrad)*nang;
  //int estart2 = find_center_of_grid(1,nrad+1)*nang;

  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);

  double* grid = new double[gsa6];
  double* wt = new double[gsa];

  double** val1 = new double*[iN];
  double** val2 = new double*[iN];
  double** val3 = new double*[iN];
  double** val4 = new double*[iN];
  for (int i=0;i<iN;i++)
    val1[i] = new double[gsa];
  for (int i=0;i<iN;i++)
    val2[i] = new double[gsa];
  for (int i=0;i<iN;i++)
    val3[i] = new double[gsa];
  for (int i=0;i<iN;i++)
    val4[i] = new double[gsa];

  for (int i=0;i<N2*N2;i++)
    ol[i] = 0.;

  if (prl>0 && N>40) printf("  working on ints");

  #pragma omp parallel for schedule(static) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
    #pragma acc enter data copyin(n2i[0:natoms])

    #pragma acc enter data create(grid[0:gsa6],wt[0:gsa])

    #pragma acc enter data create(val1[0:iN][0:gsa],val2[0:iN][0:gsa],val3[0:iN][0:gsa],val4[0:iN][0:gsa])
  }
  acc_set_device_num(0,acc_device_nvidia);

  for (int m=0;m<natoms;m++)
  {
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];
    float Z1 = (float)atno[m];
    float A1 = coordsf[3*m+0]; float B1 = coordsf[3*m+1]; float C1 = coordsf[3*m+2];

   //setup for grids
    int natoms1 = 1;
    int atno1[natoms1];
    atno1[0] = atno[m];
    float coordsf1[3*natoms1];
    coordsf1[0] = A1; coordsf1[1] = B1; coordsf1[2] = C1;

    if (prl>0 && N>40) { printf("."); fflush(stdout); }
    int gsa1 = gs*natoms1;
    compute_4c_ol_oned(ngpu,natoms1,atno1,coordsf1,nrad,nang,ang_g,ang_w,basis,s1,s2,s1,s2,s1,s2,s1,s2,iN,gsa1,grid,wt,val1,val2,val3,val4,ol);

    for (int n=m+1;n<natoms;n++)
    {
      int s3 = n2i[n-1]; int s4 = n2i[n];
      float Z2 = (float)atno[n];
      float A2 = coordsf[3*n+0]; float B2 = coordsf[3*n+1]; float C2 = coordsf[3*n+2];

      int natoms2 = 2;
      int atno2[natoms2];
      atno2[0] = atno[m];
      atno2[1] = atno[n];
      float coordsf2[3*natoms2];
      coordsf2[0] = A1; coordsf2[1] = B1; coordsf2[2] = C1;
      coordsf2[3] = A2; coordsf2[4] = B2; coordsf2[5] = C2;

      int gsa2 = gs*natoms2;

      if (prl>0) printf(".");
      compute_4c_ol_oned(ngpu,natoms2,atno2,coordsf2,nrad,nang,ang_g,ang_w,basis,s1,s2,s1,s2,s3,s4,s3,s4,iN,gsa2,grid,wt,val1,val2,val3,val4,ol);
      compute_4c_ol_oned(ngpu,natoms2,atno2,coordsf2,nrad,nang,ang_g,ang_w,basis,s1,s2,s3,s4,s1,s2,s3,s4,iN,gsa2,grid,wt,val1,val2,val3,val4,ol);

      if (prl>0) printf(".");
      compute_4c_ol_oned(ngpu,natoms2,atno2,coordsf2,nrad,nang,ang_g,ang_w,basis,s1,s2,s1,s2,s1,s2,s3,s4,iN,gsa2,grid,wt,val1,val2,val3,val4,ol);
      compute_4c_ol_oned(ngpu,natoms2,atno2,coordsf2,nrad,nang,ang_g,ang_w,basis,s1,s2,s3,s4,s3,s4,s3,s4,iN,gsa2,grid,wt,val1,val2,val3,val4,ol);

      for (int p=n+1;p<natoms;p++)
      {
        int s5 = n2i[p-1]; int s6 = n2i[p];
        float Z3 = (float)atno[p];
        float A3 = coordsf[3*p+0]; float B3 = coordsf[3*p+1]; float C3 = coordsf[3*p+2];

        int natoms3 = 3;
        int atno3[natoms3];
        atno3[0] = atno[m];
        atno3[1] = atno[n];
        atno3[2] = atno[p];
        float coordsf3[3*natoms3];
        coordsf3[0] = A1; coordsf3[1] = B1; coordsf3[2] = C1;
        coordsf3[3] = A2; coordsf3[4] = B2; coordsf3[5] = C2;
        coordsf3[6] = A3; coordsf3[7] = B3; coordsf3[8] = C3;

        int gsa3 = gs*natoms3;

        if (prl>0) printf(".");
       //11 35
        compute_4c_ol_oned(ngpu,natoms3,atno3,coordsf3,nrad,nang,ang_g,ang_w,basis,s1,s2,s1,s2,s3,s4,s5,s6,iN,gsa3,grid,wt,val1,val2,val3,val4,ol);
       //13 15
        compute_4c_ol_oned(ngpu,natoms3,atno3,coordsf3,nrad,nang,ang_g,ang_w,basis,s1,s2,s3,s4,s1,s2,s5,s6,iN,gsa3,grid,wt,val1,val2,val3,val4,ol);

        if (prl>0) printf(".");
       //13 35
        compute_4c_ol_oned(ngpu,natoms3,atno3,coordsf3,nrad,nang,ang_g,ang_w,basis,s1,s2,s3,s4,s3,s4,s5,s6,iN,gsa3,grid,wt,val1,val2,val3,val4,ol);
       //15 33
        compute_4c_ol_oned(ngpu,natoms3,atno3,coordsf3,nrad,nang,ang_g,ang_w,basis,s1,s2,s5,s6,s3,s4,s3,s4,iN,gsa3,grid,wt,val1,val2,val3,val4,ol);

        if (prl>0) printf(".");
       //13 55
        compute_4c_ol_oned(ngpu,natoms3,atno3,coordsf3,nrad,nang,ang_g,ang_w,basis,s1,s2,s3,s4,s5,s6,s5,s6,iN,gsa3,grid,wt,val1,val2,val3,val4,ol);
       //15 35
        compute_4c_ol_oned(ngpu,natoms3,atno3,coordsf3,nrad,nang,ang_g,ang_w,basis,s1,s2,s5,s6,s3,s4,s5,s6,iN,gsa3,grid,wt,val1,val2,val3,val4,ol);

        for (int q=p+1;q<natoms;q++)
        {
          int s7 = n2i[q-1]; int s8 = n2i[q];
          float Z4 = (float)atno[q];
          float A4 = coordsf[3*q+0]; float B4 = coordsf[3*q+1]; float C4 = coordsf[3*q+2];

          int natoms4 = 4;
          int atno4[natoms4];
          atno4[0] = atno[m];
          atno4[1] = atno[n];
          atno4[2] = atno[p];
          atno4[3] = atno[q];
          float coordsf4[3*natoms4];
          coordsf4[0] = A1; coordsf4[1]  = B1; coordsf4[2]  = C1;
          coordsf4[3] = A2; coordsf4[4]  = B2; coordsf4[5]  = C2;
          coordsf4[6] = A3; coordsf4[7]  = B3; coordsf4[8]  = C3;
          coordsf4[9] = A4; coordsf4[10] = B4; coordsf4[11] = C4;

          int gsa4 = gs*natoms4;

         //incomplete/untested
          if (prl>0) printf(".");
          compute_4c_ol_oned(ngpu,natoms4,atno4,coordsf4,nrad,nang,ang_g,ang_w,basis,s1,s2,s3,s4,s5,s6,s7,s8,iN,gsa4,grid,wt,val1,val2,val3,val4,ol);
          compute_4c_ol_oned(ngpu,natoms4,atno4,coordsf4,nrad,nang,ang_g,ang_w,basis,s1,s2,s5,s6,s3,s4,s7,s8,iN,gsa4,grid,wt,val1,val2,val3,val4,ol);
          compute_4c_ol_oned(ngpu,natoms4,atno4,coordsf4,nrad,nang,ang_g,ang_w,basis,s1,s2,s7,s8,s3,s4,s5,s6,iN,gsa4,grid,wt,val1,val2,val3,val4,ol);

        } //loop q over natoms

      } //loop p over natoms

    } //loop n over natoms

  } //loop m over natoms

 //normalization
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
  for (int l=0;l<N;l++)
    ol[i*N3+j*N2+k*N+l] *= basis[i][4]*basis[j][4]*basis[k][4]*basis[l][4];

  if (prl>0 && N>40) { printf(" done \n"); fflush(stdout); }

 //show the symmetry
  if (prl>2)
  for (int p=0;p<N;p++)
  for (int q=0;q<=p;q++)
  for (int r=0;r<=p;r++)
  for (int s=0;s<=r;s++)
  {
    double v1 = ol[p*N3+q*N2+r*N+s];
    double v2 = ol[p*N3+q*N2+s*N+r];
    double v3 = ol[q*N3+p*N2+r*N+s];
    double v4 = ol[q*N3+p*N2+s*N+r];

    double v5 = ol[r*N3+s*N2+p*N+q];
    double v6 = ol[s*N3+r*N2+p*N+q];
    double v7 = ol[r*N3+s*N2+q*N+p];
    double v8 = ol[s*N3+r*N2+q*N+p];

    printf("   pqrs: %i %i %i %i: %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f \n",p,q,r,s,v1,v2,v3,v4,v5,v6,v7,v8);
  }

  if ((prl>1 && N<10) || prl>2)
  {
    printf("\n ol: \n");
    for (int i=0;i<N;i++)
    for (int j=0;j<N;j++)
    {
      printf(" ij: %2i %2i \n",i,j);
      print_square(N,&ol[i*N3+j*N2]);
    }
  }


  #pragma omp parallel for schedule(static) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc exit data delete(grid[0:gsa6],wt[0:gsa])
    #pragma acc exit data delete(val1[0:iN][0:gsa],val2[0:iN][0:gsa],val3[0:iN][0:gsa],val4[0:iN][0:gsa])

    #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
    #pragma acc exit data delete(n2i[0:natoms])
  }
  acc_set_device_num(0,acc_device_nvidia);

  delete [] n2i;

  for (int j=0;j<iN;j++)
    delete [] val1[j];
  delete [] val1;
  for (int j=0;j<iN;j++)
    delete [] val2[j];
  delete [] val2;
  for (int j=0;j<iN;j++)
    delete [] val3[j];
  delete [] val3;
  for (int j=0;j<iN;j++)
    delete [] val4[j];
  delete [] val4;
 
  delete [] grid;
  delete [] wt;

  return;
}

void compute_4c_ol_one(int ngpu, int natoms, int* atno, float* coordsf, int nrad, int nang, float* ang_g, float* ang_w, vector<vector<double> >& basis, int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8, int iN,
   int gsa, float* grid, float* wt, float** val1, float** val2, float** val3, float** val4, double* ol)
{
 //this function needs much better batching,
 // parallel utilization is poor

  int nomp = ngpu;
  int prl = 1;

  int gc = 6;
  int gsa6 = gsa*6;

  //double coords[3*natoms];
  //for (int n=0;n<3*natoms;n++) coords[n] = coordsf[n];

  int N = basis.size();
  int N2 = N*N;
  int N3 = N2*N;

  double minval = 1.e-16;

 //list out integrals that need to be computed
  vector<vector<short> > ijkl = enumerate_ijkl(s1,s2,s3,s4,s5,s6,s7,s8);
  int ijkl_size = ijkl.size();

  acc_set_device_num(0,acc_device_nvidia);

 //integral-specific grids
  get_becke_grid_full(0.,natoms,atno,coordsf,nrad,nang,ang_g,ang_w,gc,grid,wt);

  #pragma acc update self(grid[0:gsa6],wt[0:gsa])

 #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=1;n<nomp;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc update device(grid[0:gsa6],wt[0:gsa])
  }
  acc_set_device_num(0,acc_device_nvidia);

 #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    //int tid = omp_get_thread_num();
    acc_set_device_num(n,acc_device_nvidia);

   #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gsa])
    for (int ii1=0;ii1<s2-s1;ii1++)
    for (int j=0;j<gsa;j++)
      val1[ii1][j] = 1.f;
   #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gsa])
    for (int ii2=0;ii2<s4-s3;ii2++)
    for (int j=0;j<gsa;j++)
      val2[ii2][j] = 1.f;
   #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gsa])
    for (int ii3=0;ii3<s6-s5;ii3++)
    for (int j=0;j<gsa;j++)
      val3[ii3][j] = 1.f;
   #pragma acc parallel loop collapse(2) present(val4[0:iN][0:gsa])
    for (int ii4=0;ii4<s8-s7;ii4++)
    for (int j=0;j<gsa;j++)
      val4[ii4][j] = 1.f;
  }

 #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    //int tid = omp_get_thread_num();
    acc_set_device_num(n,acc_device_nvidia);

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];
      float A1 = basis1[5]; float B1 = basis1[6]; float C1 = basis1[7];

      recenter_grid_zero(gsa,grid,-A1,-B1,-C1);
      eval_sh(ii1,gsa,grid,val1[ii1],n1,l1,m1,zeta1);
      recenter_grid_zero(gsa,grid,A1,B1,C1);
    } //loop i1

    for (int i2=s3;i2<s4;i2++)
    {
      int ii2 = i2-s3;

      vector<double> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];
      float A1 = basis2[5]; float B1 = basis2[6]; float C1 = basis2[7];

      recenter_grid_zero(gsa,grid,-A1,-B1,-C1);
      eval_sh(ii2,gsa,grid,val2[ii2],n2,l2,m2,zeta2);
      recenter_grid_zero(gsa,grid,A1,B1,C1);
    } //loop i2

    for (int i3=s5;i3<s6;i3++)
    {
      int ii3 = i3-s5;

      vector<double> basis3 = basis[i3];
      int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];
      float A1 = basis3[5]; float B1 = basis3[6]; float C1 = basis3[7];

      recenter_grid_zero(gsa,grid,-A1,-B1,-C1);
      eval_sh(ii3,gsa,grid,val3[ii3],n3,l3,m3,zeta3);
      recenter_grid_zero(gsa,grid,A1,B1,C1);
    } //loop i3

    for (int i4=s7;i4<s8;i4++)
    {
      int ii4 = i4-s7;

      vector<double> basis4 = basis[i4];
      int n4 = basis4[0]; int l4 = basis4[1]; int m4 = basis4[2]; double zeta4 = basis4[3];
      float A1 = basis4[5]; float B1 = basis4[6]; float C1 = basis4[7];

      recenter_grid_zero(gsa,grid,-A1,-B1,-C1);
      eval_sh(ii4,gsa,grid,val4[ii4],n4,l4,m4,zeta4);
      recenter_grid_zero(gsa,grid,A1,B1,C1);
    } //loop i4
  }
  //acc_set_device_num(0,acc_device_nvidia);

 #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int i0=0;i0<ijkl_size;i0++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);

    int i1 = ijkl[i0][0]; int i2 = ijkl[i0][1];
    int i3 = ijkl[i0][2]; int i4 = ijkl[i0][3];
    int ii1 = i1-s1; int ii2 = i2-s3;
    int ii3 = i3-s5; int ii4 = i4-s7;

    //float* val1n = val1[ii1]; float* val2n = val2[ii2];
    //float* val3n = val3[ii3]; float* val4n = val4[ii4];

    //printf(" tid %i working on ijkl: %i %i %i %i / %i %i %i %i \n",tid,i1,i2,i3,i4,ii1,ii2,ii3,ii4);

    double v1 = 0.;
   #pragma acc parallel loop present(grid[0:gsa6],val1[0:iN][0:gsa],val2[0:iN][0:gsa],val3[0:iN][0:gsa],val4[0:iN][0:gsa],wt[0:gsa]) reduction(+:v1)
    for (int j=0;j<gsa;j++)
    {
      //float vj = val1n[j]*val2n[j]*wt1[j];
      //float vj = val1[ii1][j]*val2[ii2][j]*wt1[j];

      //v2 += val3n[k]*val4n[k]*wt2[k]/r1;
      //v2 += val3[ii3][k]*val4[ii4][k]*wt2[k]/r1;

      v1 += val1[ii1][j]*val2[ii2][j]*val3[ii3][j]*val4[ii4][j]*wt[j];
    }

    if (fabs(v1)<minval) v1 = 0.;
  
    ol[i1*N3+i2*N2+i3*N+i4] = v1;
    ol[i1*N3+i2*N2+i4*N+i3] = v1;
    ol[i2*N3+i1*N2+i3*N+i4] = v1;
    ol[i2*N3+i1*N2+i4*N+i3] = v1;

    ol[i3*N3+i4*N2+i1*N+i2] = v1;
    ol[i4*N3+i3*N2+i1*N+i2] = v1;
    ol[i3*N3+i4*N2+i2*N+i1] = v1;
    ol[i4*N3+i3*N2+i2*N+i1] = v1;

    if (prl>1) printf("  v1(%i %i %i %i): %8.5f \n",i1,i2,i3,i4,v1);

  } //loop i0 over ijkl

  return;
}

#if RED_DOUBLE
void compute_all_4c_ol_gen(int ngpu, int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* ol, int prl)
#else
void compute_all_4c_ol_gen(int ngpu, int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, float* ol, int prl)
#endif
{
  if (prl>-1) printf(" beginning compute_all_4c_ol_gen (ngpu: %i) \n",ngpu);

  int nomp = ngpu;

  int N = basis.size();
  int N2 = N*N;
  int N3 = N2*N;

  int gs = nrad*nang;

 //grid will contain up to 4 nuclei
  int natomsa = min(4,natoms);
  int gsa = natomsa*gs;
  int gsa6 = 6*gsa;

  if (prl>0) printf("  gs: %4i  natomsa: %i \n",gs,natomsa);

  //int estart = find_center_of_grid(1,nrad)*nang;
  //int estart2 = find_center_of_grid(1,nrad+1)*nang;

  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);

  float* grid = new float[gsa6];
  float* wt = new float[gsa];

  float** val1 = new float*[iN];
  float** val2 = new float*[iN];
  float** val3 = new float*[iN];
  float** val4 = new float*[iN];
  for (int i=0;i<iN;i++)
    val1[i] = new float[gsa];
  for (int i=0;i<iN;i++)
    val2[i] = new float[gsa];
  for (int i=0;i<iN;i++)
    val3[i] = new float[gsa];
  for (int i=0;i<iN;i++)
    val4[i] = new float[gsa];

  float* ang_g = new float[3*nang];
  float* ang_w = new float[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

  for (int i=0;i<N2*N2;i++)
    ol[i] = 0.;

  if (prl>0 && N>40) printf("  working on ints");

  #pragma omp parallel for schedule(static) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
    #pragma acc enter data copyin(n2i[0:natoms])

    #pragma acc enter data create(grid[0:gsa6],wt[0:gsa])

    #pragma acc enter data create(val1[0:iN][0:gsa],val2[0:iN][0:gsa],val3[0:iN][0:gsa],val4[0:iN][0:gsa])
  }
  acc_set_device_num(0,acc_device_nvidia);

  for (int m=0;m<natoms;m++)
  {
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];
    float Z1 = (float)atno[m];
    float A1 = coords[3*m+0]; float B1 = coords[3*m+1]; float C1 = coords[3*m+2];

   //setup for grids
    int natoms1 = 1;
    int atno1[natoms1];
    atno1[0] = atno[m];
    float coordsf1[3*natoms1];
    coordsf1[0] = A1; coordsf1[1] = B1; coordsf1[2] = C1;

    if (prl>0 && N>40) printf(".");
    int gsa1 = gs*natoms1;
    compute_4c_ol_one(ngpu,natoms1,atno1,coordsf1,nrad,nang,ang_g,ang_w,basis,s1,s2,s1,s2,s1,s2,s1,s2,iN,gsa1,grid,wt,val1,val2,val3,val4,ol);

    for (int n=m+1;n<natoms;n++)
    {
      int s3 = n2i[n-1]; int s4 = n2i[n];
      float Z2 = (float)atno[n];
      float A2 = coords[3*n+0]; float B2 = coords[3*n+1]; float C2 = coords[3*n+2];

      int natoms2 = 2;
      int atno2[natoms2];
      atno2[0] = atno[m];
      atno2[1] = atno[n];
      float coordsf2[3*natoms2];
      coordsf2[0] = A1; coordsf2[1] = B1; coordsf2[2] = C1;
      coordsf2[3] = A2; coordsf2[4] = B2; coordsf2[5] = C2;

      int gsa2 = gs*natoms2;

      if (prl>0) printf(".");
      compute_4c_ol_one(ngpu,natoms2,atno2,coordsf2,nrad,nang,ang_g,ang_w,basis,s1,s2,s1,s2,s3,s4,s3,s4,iN,gsa2,grid,wt,val1,val2,val3,val4,ol);
      compute_4c_ol_one(ngpu,natoms2,atno2,coordsf2,nrad,nang,ang_g,ang_w,basis,s1,s2,s3,s4,s1,s2,s3,s4,iN,gsa2,grid,wt,val1,val2,val3,val4,ol);

      if (prl>0) printf(".");
      compute_4c_ol_one(ngpu,natoms2,atno2,coordsf2,nrad,nang,ang_g,ang_w,basis,s1,s2,s1,s2,s1,s2,s3,s4,iN,gsa2,grid,wt,val1,val2,val3,val4,ol);
      compute_4c_ol_one(ngpu,natoms2,atno2,coordsf2,nrad,nang,ang_g,ang_w,basis,s1,s2,s3,s4,s3,s4,s3,s4,iN,gsa2,grid,wt,val1,val2,val3,val4,ol);

      for (int p=n+1;p<natoms;p++)
      {
        int s5 = n2i[p-1]; int s6 = n2i[p];
        float Z3 = (float)atno[p];
        float A3 = coords[3*p+0]; float B3 = coords[3*p+1]; float C3 = coords[3*p+2];

        int natoms3 = 3;
        int atno3[natoms3];
        atno3[0] = atno[m];
        atno3[1] = atno[n];
        atno3[2] = atno[p];
        float coordsf3[3*natoms3];
        coordsf3[0] = A1; coordsf3[1] = B1; coordsf3[2] = C1;
        coordsf3[3] = A2; coordsf3[4] = B2; coordsf3[5] = C2;
        coordsf3[6] = A3; coordsf3[7] = B3; coordsf3[8] = C3;

        int gsa3 = gs*natoms3;

        if (prl>0) printf(".");
       //11 35
        compute_4c_ol_one(ngpu,natoms3,atno3,coordsf3,nrad,nang,ang_g,ang_w,basis,s1,s2,s1,s2,s3,s4,s5,s6,iN,gsa3,grid,wt,val1,val2,val3,val4,ol);
       //13 15
        compute_4c_ol_one(ngpu,natoms3,atno3,coordsf3,nrad,nang,ang_g,ang_w,basis,s1,s2,s3,s4,s1,s2,s5,s6,iN,gsa3,grid,wt,val1,val2,val3,val4,ol);

        if (prl>0) printf(".");
       //13 35
        compute_4c_ol_one(ngpu,natoms3,atno3,coordsf3,nrad,nang,ang_g,ang_w,basis,s1,s2,s3,s4,s3,s4,s5,s6,iN,gsa3,grid,wt,val1,val2,val3,val4,ol);
       //15 33
        compute_4c_ol_one(ngpu,natoms3,atno3,coordsf3,nrad,nang,ang_g,ang_w,basis,s1,s2,s5,s6,s3,s4,s3,s4,iN,gsa3,grid,wt,val1,val2,val3,val4,ol);

        if (prl>0) printf(".");
       //13 55
        compute_4c_ol_one(ngpu,natoms3,atno3,coordsf3,nrad,nang,ang_g,ang_w,basis,s1,s2,s3,s4,s5,s6,s5,s6,iN,gsa3,grid,wt,val1,val2,val3,val4,ol);
       //15 35
        compute_4c_ol_one(ngpu,natoms3,atno3,coordsf3,nrad,nang,ang_g,ang_w,basis,s1,s2,s5,s6,s3,s4,s5,s6,iN,gsa3,grid,wt,val1,val2,val3,val4,ol);

        for (int q=p+1;q<natoms;q++)
        {
          int s7 = n2i[q-1]; int s8 = n2i[q];
          float Z4 = (float)atno[q];
          float A4 = coords[3*q+0]; float B4 = coords[3*q+1]; float C4 = coords[3*q+2];

          int natoms4 = 4;
          int atno4[natoms4];
          atno4[0] = atno[m];
          atno4[1] = atno[n];
          atno4[2] = atno[p];
          atno4[3] = atno[q];
          float coordsf4[3*natoms4];
          coordsf4[0] = A1; coordsf4[1]  = B1; coordsf4[2]  = C1;
          coordsf4[3] = A2; coordsf4[4]  = B2; coordsf4[5]  = C2;
          coordsf4[6] = A3; coordsf4[7]  = B3; coordsf4[8]  = C3;
          coordsf4[9] = A4; coordsf4[10] = B4; coordsf4[11] = C4;

          int gsa4 = gs*natoms4;

         //incomplete/untested
          if (prl>0) printf(".");
          compute_4c_ol_one(ngpu,natoms4,atno4,coordsf4,nrad,nang,ang_g,ang_w,basis,s1,s2,s3,s4,s5,s6,s7,s8,iN,gsa4,grid,wt,val1,val2,val3,val4,ol);
          compute_4c_ol_one(ngpu,natoms4,atno4,coordsf4,nrad,nang,ang_g,ang_w,basis,s1,s2,s5,s6,s3,s4,s7,s8,iN,gsa4,grid,wt,val1,val2,val3,val4,ol);
          compute_4c_ol_one(ngpu,natoms4,atno4,coordsf4,nrad,nang,ang_g,ang_w,basis,s1,s2,s7,s8,s3,s4,s5,s6,iN,gsa4,grid,wt,val1,val2,val3,val4,ol);

        } //loop q over natoms

      } //loop p over natoms

    } //loop n over natoms

  } //loop m over natoms

 //normalization
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
  for (int l=0;l<N;l++)
    ol[i*N3+j*N2+k*N+l] *= basis[i][4]*basis[j][4]*basis[k][4]*basis[l][4];

  if (prl>0 && N>40) printf(" done \n");

 //show the symmetry
  if (prl>2)
  for (int p=0;p<N;p++)
  for (int q=0;q<=p;q++)
  for (int r=0;r<=p;r++)
  for (int s=0;s<=r;s++)
  {
    double v1 = ol[p*N3+q*N2+r*N+s];
    double v2 = ol[p*N3+q*N2+s*N+r];
    double v3 = ol[q*N3+p*N2+r*N+s];
    double v4 = ol[q*N3+p*N2+s*N+r];

    double v5 = ol[r*N3+s*N2+p*N+q];
    double v6 = ol[s*N3+r*N2+p*N+q];
    double v7 = ol[r*N3+s*N2+q*N+p];
    double v8 = ol[s*N3+r*N2+q*N+p];

    printf("   pqrs: %i %i %i %i: %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f \n",p,q,r,s,v1,v2,v3,v4,v5,v6,v7,v8);
  }

  if ((prl>1 && N<10) || prl>2)
  {
    printf("\n ol: \n");
    for (int i=0;i<N;i++)
    for (int j=0;j<N;j++)
    {
      printf(" ij: %2i %2i \n",i,j);
      print_square(N,&ol[i*N3+j*N2]);
    }
  }


  #pragma omp parallel for schedule(static) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc exit data delete(grid[0:gsa6],wt[0:gsa])
    #pragma acc exit data delete(val1[0:iN][0:gsa],val2[0:iN][0:gsa],val3[0:iN][0:gsa],val4[0:iN][0:gsa])

    #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
    #pragma acc exit data delete(n2i[0:natoms])
  }
  acc_set_device_num(0,acc_device_nvidia);

  delete [] n2i;

  delete [] ang_g;
  delete [] ang_w;

  for (int j=0;j<iN;j++)
    delete [] val1[j];
  delete [] val1;
  for (int j=0;j<iN;j++)
    delete [] val2[j];
  delete [] val2;
  for (int j=0;j<iN;j++)
    delete [] val3[j];
  delete [] val3;
  for (int j=0;j<iN;j++)
    delete [] val4[j];
  delete [] val4;
 
  delete [] grid;
  delete [] wt;

  return;
}

void compute_4c_one(int ngpu, int natoms, int* atno, float* coordsf, int nrad, int nang, float* ang_g, float* ang_w, vector<vector<double> >& basis, int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8, int iN,
   int gsa, int gspa, float* grid1, float* wt1, float* grid2, float* wt2, float** val1, float** val2, float** val3, float** val4, double* g)
{
  int nomp = ngpu;

  int gc = 6;
  int gsa6 = gsa*6;
  int gspa6 = gspa*6;

  int N = basis.size();
  int N2 = N*N;
  int N3 = N2*N;
  //int M = iN;
  //int M2 = M*M;
  //int M3 = M2*M;
  //int M4 = M2*M2;

 //list out integrals that need to be computed
  vector<vector<short> > ijkl = enumerate_ijkl(s1,s2,s3,s4,s5,s6,s7,s8);
  int ijkl_size = ijkl.size();

  acc_set_device_num(0,acc_device_nvidia);

 //integral-specific grids
  get_becke_grid_full(0.,natoms,atno,coordsf,nrad,  nang,ang_g,ang_w,gc,grid1,wt1);
  get_becke_grid_full(0.,natoms,atno,coordsf,nrad+1,nang,ang_g,ang_w,gc,grid2,wt2);

 //look for close-lying grid points
  clean_4c_grid(gsa,gspa,grid1,wt1,grid2,wt2);

  #pragma acc update self(grid1[0:gsa6],wt1[0:gsa],grid2[0:gspa6],wt2[0:gspa])

 #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=1;n<nomp;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc update device(grid1[0:gsa6],wt1[0:gsa],grid2[0:gspa6],wt2[0:gspa])
  }
  acc_set_device_num(0,acc_device_nvidia);

 #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    //int tid = omp_get_thread_num();
    acc_set_device_num(n,acc_device_nvidia);

   #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gsa])
    for (int ii1=0;ii1<s2-s1;ii1++)
    for (int j=0;j<gsa;j++)
      val1[ii1][j] = 1.f;
   #pragma acc parallel loop collapse(2) present(val2[0:iN][0:gsa])
    for (int ii2=0;ii2<s4-s3;ii2++)
    for (int j=0;j<gsa;j++)
      val2[ii2][j] = 1.f;
   #pragma acc parallel loop collapse(2) present(val3[0:iN][0:gspa])
    for (int ii3=0;ii3<s6-s5;ii3++)
    for (int j=0;j<gspa;j++)
      val3[ii3][j] = 1.f;
   #pragma acc parallel loop collapse(2) present(val4[0:iN][0:gspa])
    for (int ii4=0;ii4<s8-s7;ii4++)
    for (int j=0;j<gspa;j++)
      val4[ii4][j] = 1.f;
  }

 #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    //int tid = omp_get_thread_num();
    acc_set_device_num(n,acc_device_nvidia);

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];
      float A1 = basis1[5]; float B1 = basis1[6]; float C1 = basis1[7];

      recenter_grid_zero(gsa,grid1,-A1,-B1,-C1);
      eval_sh(ii1,gsa,grid1,val1[ii1],n1,l1,m1,zeta1);
      recenter_grid_zero(gsa,grid1,A1,B1,C1);
    } //loop i1

    for (int i2=s3;i2<s4;i2++)
    {
      int ii2 = i2-s3;

      vector<double> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];
      float A1 = basis2[5]; float B1 = basis2[6]; float C1 = basis2[7];

      recenter_grid_zero(gsa,grid1,-A1,-B1,-C1);
      eval_sh(ii2,gsa,grid1,val2[ii2],n2,l2,m2,zeta2);
      recenter_grid_zero(gsa,grid1,A1,B1,C1);
    } //loop i2

    for (int i3=s5;i3<s6;i3++)
    {
      int ii3 = i3-s5;

      vector<double> basis3 = basis[i3];
      int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];
      float A1 = basis3[5]; float B1 = basis3[6]; float C1 = basis3[7];

      recenter_grid_zero(gspa,grid2,-A1,-B1,-C1);
      eval_sh(ii3,gspa,grid2,val3[ii3],n3,l3,m3,zeta3);
      recenter_grid_zero(gspa,grid2,A1,B1,C1);
    } //loop i3

    for (int i4=s7;i4<s8;i4++)
    {
      int ii4 = i4-s7;

      vector<double> basis4 = basis[i4];
      int n4 = basis4[0]; int l4 = basis4[1]; int m4 = basis4[2]; double zeta4 = basis4[3];
      float A1 = basis4[5]; float B1 = basis4[6]; float C1 = basis4[7];

      recenter_grid_zero(gspa,grid2,-A1,-B1,-C1);
      eval_sh(ii4,gspa,grid2,val4[ii4],n4,l4,m4,zeta4);
      recenter_grid_zero(gspa,grid2,A1,B1,C1);
    } //loop i4
  }
  //acc_set_device_num(0,acc_device_nvidia);

 #pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int i0=0;i0<ijkl_size;i0++)
  {
    int tid = omp_get_thread_num();
    acc_set_device_num(tid,acc_device_nvidia);

    int i1 = ijkl[i0][0]; int i2 = ijkl[i0][1];
    int i3 = ijkl[i0][2]; int i4 = ijkl[i0][3];
    //int i1 = ijkl_gpu[4*i0+0]; int i2 = ijkl_gpu[4*i0+1];
    //int i3 = ijkl_gpu[4*i0+2]; int i4 = ijkl_gpu[4*i0+3];
    int ii1 = i1-s1; int ii2 = i2-s3;
    int ii3 = i3-s5; int ii4 = i4-s7;

    //float* val1n = val1[ii1]; float* val2n = val2[ii2];
    //float* val3n = val3[ii3]; float* val4n = val4[ii4];

    //printf(" tid %i working on ijkl: %i %i %i %i / %i %i %i %i \n",tid,i1,i2,i3,i4,ii1,ii2,ii3,ii4);

    double v1 = 0.;
   #pragma acc parallel loop present(grid1[0:gsa6],grid2[0:gspa6],val1[0:iN][0:gsa],val2[0:iN][0:gsa],val3[0:iN][0:gspa],val4[0:iN][0:gspa],wt1[0:gsa],wt2[0:gspa]) reduction(+:v1)
    for (int j=0;j<gsa;j++) //cannot use collapse: total index too high
    {
      float x1 = grid1[6*j+0]; float y1 = grid1[6*j+1]; float z1 = grid1[6*j+2];

      //float vj = val1n[j]*val2n[j]*wt1[j];
      float vj = val1[ii1][j]*val2[ii2][j]*wt1[j];
      float v2 = 0.f;

      if (fabs(vj)>1.e-10)
      #pragma acc loop reduction(+:v2)
      for (int k=0;k<gspa;k++)
      {
        float x2 = grid2[6*k+0]; float y2 = grid2[6*k+1]; float z2 = grid2[6*k+2];
        float x12 = x1-x2; float y12 = y1-y2; float z12 = z1-z2;

        float r1 = sqrtf(x12*x12+y12*y12+z12*z12)+1.e-12;
 
        //v2 += val3n[k]*val4n[k]*wt2[k]/r1;
        v2 += val3[ii3][k]*val4[ii4][k]*wt2[k]/r1;
      }
  
      v1 += v2*vj;
    }
  
    g[i1*N3+i2*N2+i3*N+i4] = v1;
    g[i1*N3+i2*N2+i4*N+i3] = v1;
    g[i2*N3+i1*N2+i3*N+i4] = v1;
    g[i2*N3+i1*N2+i4*N+i3] = v1;

    g[i3*N3+i4*N2+i1*N+i2] = v1;
    g[i4*N3+i3*N2+i1*N+i2] = v1;
    g[i3*N3+i4*N2+i2*N+i1] = v1;
    g[i4*N3+i3*N2+i2*N+i1] = v1;

    //gt[ii1*M3+ii2*M2+ii3*M+ii4] = v1;
    //printf("  v1(%i %i %i %i): %8.5f \n",i1,i2,i3,i4,v1);

  } //loop i0 over ijkl

  return;
}

#if RED_DOUBLE
void compute_all_4c_gen(int ngpu, int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* g, int prl)
#else
void compute_all_4c_gen(int ngpu, int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, float* g, int prl)
#endif
{
  if (prl>-1) printf(" beginning compute_all_4c_gen (ngpu: %i) \n",ngpu);

  int nomp = ngpu;

  int N = basis.size();
  int N2 = N*N;
  int N3 = N2*N;

  int gs = nrad*nang;
  int gsp = (nrad+1)*nang;

 //grid will contain up to 4 nuclei
  int natomsa = min(4,natoms);
  int gsa = natomsa*gs;
  int gspa = natomsa*gsp;
  int gsa6 = 6*gsa;
  int gspa6 = 6*gspa;

  if (prl>0) printf("  gs/p: %4i %4i  natomsa: %i \n",gs,gsp,natomsa);


  //int estart = find_center_of_grid(1,nrad)*nang;
  //int estart2 = find_center_of_grid(1,nrad+1)*nang;

  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);

  float* grid1 = new float[gsa6];
  float* wt1 = new float[gsa];

  float* grid2 = new float[gspa6];
  float* wt2 = new float[gspa];

  float** val1 = new float*[iN];
  float** val2 = new float*[iN];
  float** val3 = new float*[iN];
  float** val4 = new float*[iN];
  for (int i=0;i<iN;i++)
    val1[i] = new float[gsa];
  for (int i=0;i<iN;i++)
    val2[i] = new float[gsa];
  for (int i=0;i<iN;i++)
    val3[i] = new float[gspa];
  for (int i=0;i<iN;i++)
    val4[i] = new float[gspa];

  float* ang_g = new float[3*nang];
  float* ang_w = new float[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

  for (int i=0;i<N2*N2;i++)
    g[i] = 0.;

 //g not on gpu (too large)
  //int M = iN;
  //int M2 = M*M;
  //int M4 = M2*M2;

  #pragma omp parallel for schedule(static) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
    #pragma acc enter data copyin(n2i[0:natoms])

    #pragma acc enter data create(grid1[0:gsa6],wt1[0:gsa])
    #pragma acc enter data create(grid2[0:gspa6],wt2[0:gspa])

    #pragma acc enter data create(val1[0:iN][0:gsa],val2[0:iN][0:gsa],val3[0:iN][0:gspa],val4[0:iN][0:gspa])
  }
  acc_set_device_num(0,acc_device_nvidia);

  if (prl>0) printf("  working on ints");
  for (int m=0;m<natoms;m++)
  {
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];
    float Z1 = (float)atno[m];
    float A1 = coords[3*m+0]; float B1 = coords[3*m+1]; float C1 = coords[3*m+2];

   //setup for grids
    int natoms1 = 1;
    int atno1[natoms1];
    atno1[0] = atno[m];
    float coordsf1[3*natoms1];
    coordsf1[0] = A1; coordsf1[1] = B1; coordsf1[2] = C1;

    if (prl>0) printf(".");
    int gsa1 = gs*natoms1; int gspa1 = gsp*natoms1;
    compute_4c_one(ngpu,natoms1,atno1,coordsf1,nrad,nang,ang_g,ang_w,basis,s1,s2,s1,s2,s1,s2,s1,s2,iN,gsa1,gspa1,grid1,wt1,grid2,wt2,val1,val2,val3,val4,g);

    for (int n=m+1;n<natoms;n++)
    {
      int s3 = n2i[n-1]; int s4 = n2i[n];
      float Z2 = (float)atno[n];
      float A2 = coords[3*n+0]; float B2 = coords[3*n+1]; float C2 = coords[3*n+2];

      int natoms2 = 2;
      int atno2[natoms2];
      atno2[0] = atno[m];
      atno2[1] = atno[n];
      float coordsf2[3*natoms2];
      coordsf2[0] = A1; coordsf2[1] = B1; coordsf2[2] = C1;
      coordsf2[3] = A2; coordsf2[4] = B2; coordsf2[5] = C2;

      int gsa2 = gs*natoms2; int gspa2 = gsp*natoms2;

      if (prl>0) printf(".");
      compute_4c_one(ngpu,natoms2,atno2,coordsf2,nrad,nang,ang_g,ang_w,basis,s1,s2,s1,s2,s3,s4,s3,s4,iN,gsa2,gspa2,grid1,wt1,grid2,wt2,val1,val2,val3,val4,g);
      compute_4c_one(ngpu,natoms2,atno2,coordsf2,nrad,nang,ang_g,ang_w,basis,s1,s2,s3,s4,s1,s2,s3,s4,iN,gsa2,gspa2,grid1,wt1,grid2,wt2,val1,val2,val3,val4,g);

      if (prl>0) printf(".");
      compute_4c_one(ngpu,natoms2,atno2,coordsf2,nrad,nang,ang_g,ang_w,basis,s1,s2,s1,s2,s1,s2,s3,s4,iN,gsa2,gspa2,grid1,wt1,grid2,wt2,val1,val2,val3,val4,g);
      compute_4c_one(ngpu,natoms2,atno2,coordsf2,nrad,nang,ang_g,ang_w,basis,s1,s2,s3,s4,s3,s4,s3,s4,iN,gsa2,gspa2,grid1,wt1,grid2,wt2,val1,val2,val3,val4,g);

      for (int p=n+1;p<natoms;p++)
      {
        int s5 = n2i[p-1]; int s6 = n2i[p];
        float Z3 = (float)atno[p];
        float A3 = coords[3*p+0]; float B3 = coords[3*p+1]; float C3 = coords[3*p+2];

        int natoms3 = 3;
        int atno3[natoms3];
        atno3[0] = atno[m];
        atno3[1] = atno[n];
        atno3[2] = atno[p];
        float coordsf3[3*natoms3];
        coordsf3[0] = A1; coordsf3[1] = B1; coordsf3[2] = C1;
        coordsf3[3] = A2; coordsf3[4] = B2; coordsf3[5] = C2;
        coordsf3[6] = A3; coordsf3[7] = B3; coordsf3[8] = C3;

        int gsa3 = gs*natoms3; int gspa3 = gsp*natoms3;

        if (prl>0) printf(".");
       //11 35
        compute_4c_one(ngpu,natoms3,atno3,coordsf3,nrad,nang,ang_g,ang_w,basis,s1,s2,s1,s2,s3,s4,s5,s6,iN,gsa3,gspa3,grid1,wt1,grid2,wt2,val1,val2,val3,val4,g);
       //13 15
        compute_4c_one(ngpu,natoms3,atno3,coordsf3,nrad,nang,ang_g,ang_w,basis,s1,s2,s3,s4,s1,s2,s5,s6,iN,gsa3,gspa3,grid1,wt1,grid2,wt2,val1,val2,val3,val4,g);

        if (prl>0) printf(".");
       //13 35
        compute_4c_one(ngpu,natoms3,atno3,coordsf3,nrad,nang,ang_g,ang_w,basis,s1,s2,s3,s4,s3,s4,s5,s6,iN,gsa3,gspa3,grid1,wt1,grid2,wt2,val1,val2,val3,val4,g);
       //15 33
        compute_4c_one(ngpu,natoms3,atno3,coordsf3,nrad,nang,ang_g,ang_w,basis,s1,s2,s5,s6,s3,s4,s3,s4,iN,gsa3,gspa3,grid1,wt1,grid2,wt2,val1,val2,val3,val4,g);

        if (prl>0) printf(".");
       //13 55
        compute_4c_one(ngpu,natoms3,atno3,coordsf3,nrad,nang,ang_g,ang_w,basis,s1,s2,s3,s4,s5,s6,s5,s6,iN,gsa3,gspa3,grid1,wt1,grid2,wt2,val1,val2,val3,val4,g);
       //15 35
        compute_4c_one(ngpu,natoms3,atno3,coordsf3,nrad,nang,ang_g,ang_w,basis,s1,s2,s5,s6,s3,s4,s5,s6,iN,gsa3,gspa3,grid1,wt1,grid2,wt2,val1,val2,val3,val4,g);

        for (int q=p+1;q<natoms;q++)
        {
          int s7 = n2i[q-1]; int s8 = n2i[q];
          float Z4 = (float)atno[q];
          float A4 = coords[3*q+0]; float B4 = coords[3*q+1]; float C4 = coords[3*q+2];

          int natoms4 = 4;
          int atno4[natoms4];
          atno4[0] = atno[m];
          atno4[1] = atno[n];
          atno4[2] = atno[p];
          atno4[3] = atno[q];
          float coordsf4[3*natoms4];
          coordsf4[0] = A1; coordsf4[1]  = B1; coordsf4[2]  = C1;
          coordsf4[3] = A2; coordsf4[4]  = B2; coordsf4[5]  = C2;
          coordsf4[6] = A3; coordsf4[7]  = B3; coordsf4[8]  = C3;
          coordsf4[9] = A4; coordsf4[10] = B4; coordsf4[11] = C4;

          int gsa4 = gs*natoms4; int gspa4 = gsp*natoms4;

         //incomplete/untested
          if (prl>0) printf(".");
          compute_4c_one(ngpu,natoms4,atno4,coordsf4,nrad,nang,ang_g,ang_w,basis,s1,s2,s3,s4,s5,s6,s7,s8,iN,gsa4,gspa4,grid1,wt1,grid2,wt2,val1,val2,val3,val4,g);
          compute_4c_one(ngpu,natoms4,atno4,coordsf4,nrad,nang,ang_g,ang_w,basis,s1,s2,s5,s6,s3,s4,s7,s8,iN,gsa4,gspa4,grid1,wt1,grid2,wt2,val1,val2,val3,val4,g);
          compute_4c_one(ngpu,natoms4,atno4,coordsf4,nrad,nang,ang_g,ang_w,basis,s1,s2,s7,s8,s3,s4,s5,s6,iN,gsa4,gspa4,grid1,wt1,grid2,wt2,val1,val2,val3,val4,g);

        } //loop q over natoms

      } //loop p over natoms

    } //loop n over natoms

  } //loop m over natoms
  if (prl>0) printf(" done \n");

 //normalization
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
  for (int l=0;l<N;l++)
    g[i*N3+j*N2+k*N+l] *= basis[i][4]*basis[j][4]*basis[k][4]*basis[l][4];

 //show the symmetry
  if (prl>1)
  for (int p=0;p<N;p++)
  for (int q=0;q<=p;q++)
  for (int r=0;r<=p;r++)
  for (int s=0;s<=r;s++)
  {
    double v1 = g[p*N3+q*N2+r*N+s];
    double v2 = g[p*N3+q*N2+s*N+r];
    double v3 = g[q*N3+p*N2+r*N+s];
    double v4 = g[q*N3+p*N2+s*N+r];

    double v5 = g[r*N3+s*N2+p*N+q];
    double v6 = g[s*N3+r*N2+p*N+q];
    double v7 = g[r*N3+s*N2+q*N+p];
    double v8 = g[s*N3+r*N2+q*N+p];

    printf("   pqrs: %i %i %i %i: %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f \n",p,q,r,s,v1,v2,v3,v4,v5,v6,v7,v8);
  }

  if ((prl>0 && N<10) || prl>1)
  {
    printf("\n g: \n");
    for (int i=0;i<N;i++)
    for (int j=0;j<N;j++)
    {
      printf(" ij: %2i %2i \n",i,j);
      print_square(N,&g[i*N3+j*N2]);
    }
  }


  #pragma omp parallel for schedule(static) num_threads(nomp)
  for (int n=0;n<nomp;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);

    #pragma acc exit data delete(grid1[0:gsa6],wt1[0:gsa])
    #pragma acc exit data delete(grid2[0:gspa6],wt2[0:gspa])
    for (int n=0;n<iN;n++)
    {
      float* val1n = val1[n]; float* val2n = val2[n]; float* val3n = val3[n]; float* val4n = val4[n];
      #pragma acc exit data delete(val1n[0:gsa],val2n[0:gsa],val3n[0:gspa],val4n[0:gspa])
    }

    #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
    #pragma acc exit data delete(n2i[0:natoms])
  }
  acc_set_device_num(0,acc_device_nvidia);

  delete [] n2i;

  delete [] ang_g;
  delete [] ang_w;

  for (int j=0;j<iN;j++)
    delete [] val1[j];
  delete [] val1;
  for (int j=0;j<iN;j++)
    delete [] val2[j];
  delete [] val2;
  for (int j=0;j<iN;j++)
    delete [] val3[j];
  delete [] val3;
  for (int j=0;j<iN;j++)
    delete [] val4[j];
  delete [] val4;
 
  delete [] grid1;
  delete [] grid2;
  delete [] wt1;
  delete [] wt2;

  return;
}

#if RED_DOUBLE
void compute_all_4c(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, double* g, int prl)
#else
void compute_all_4c(int natoms, int* atno, float* coords, vector<vector<double> > &basis, int nrad, int nang, double* ang_g0, double* ang_w0, float* g, int prl)
#endif
{
  if (prl>-1) printf(" beginning compute_all_4c \n");
  if (natoms>2)
  {
    printf("\n ERROR: compute_all_4c for 2-atom integrals only \n");
  }

  int N = basis.size();
  int N2 = N*N;
  int N3 = N2*N;

  int gs = nrad*nang;
  int gsp = (nrad+1)*nang;
  int gs6 = 6*gs;
  int gsp6 = 6*gsp;

  if (prl>0) printf("  gs/p: %4i %4i \n",gs,gsp);

  int estart = find_center_of_grid(1,nrad)*nang;
  int estart2 = find_center_of_grid(1,nrad+1)*nang;

  int* n2i = new int[natoms];
  int iN = get_imax_n2i(natoms,N,basis,n2i);

  float** grid1 = new float*[4]; for (int i=0;i<4;i++) grid1[i] = new float[gs6];
  float** wt1 = new float*[4]; for (int i=0;i<4;i++) wt1[i] = new float[gs];

  float** grid2 = new float*[4]; for (int i=0;i<4;i++) grid2[i] = new float[gsp6];
  float** wt2 = new float*[4]; for (int i=0;i<4;i++) wt2[i] = new float[gsp];

  //float** grid3 = new float*[4]; for (int i=0;i<4;i++) grid3[i] = new float[gsp6];
  //float** wt3 = new float*[4]; for (int i=0;i<4;i++) wt3[i] = new float[gsp];

  //float** grid4 = new float*[4]; for (int i=0;i<4;i++) grid4[i] = new float[gsp6];
  //float** wt4 = new float*[4]; for (int i=0;i<4;i++) wt4[i] = new float[gsp];

  float*** val1 = new float**[4];
  float*** val2 = new float**[4];
  float*** val3 = new float**[4];
  float*** val4 = new float**[4];
  for (int n=0;n<4;n++)
  {
    val1[n] = new float*[iN];
    for (int i=0;i<iN;i++) val1[n][i] = new float[gs];
  }
  for (int n=0;n<4;n++)
  {
    val2[n] = new float*[iN];
    for (int i=0;i<iN;i++) val2[n][i] = new float[gs];
  }
  for (int n=0;n<4;n++)
  {
    val3[n] = new float*[iN];
    for (int i=0;i<iN;i++) val3[n][i] = new float[gsp];
  }
  for (int n=0;n<4;n++)
  {
    val4[n] = new float*[iN];
    for (int i=0;i<iN;i++) val4[n][i] = new float[gsp];
  }

  float** wtt = new float*[16]; for (int i=0;i<16;i++) wtt[i] = new float[gsp];

  float* ang_g = new float[3*nang];
  float* ang_w = new float[nang];
  for (int i=0;i<3*nang;i++)
    ang_g[i] = ang_g0[i];
  for (int i=0;i<nang;i++)
    ang_w[i] = ang_w0[i];

  for (int i=0;i<N2*N2;i++)
    g[i] = 0.;

 //g not on gpu (too large)
  int M = iN;
  int M2 = M*M;
  //int M3 = M2*M;
  float* gt = new float[M2*M2];

 #if USE_ACC
  #pragma acc enter data create(gt[0:M2*M2])

  #pragma acc enter data copyin(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc enter data copyin(n2i[0:natoms])

  #pragma acc enter data create(grid1[0:4][0:gs6],wt1[0:4][0:gs])
  #pragma acc enter data create(grid2[0:4][0:gsp6],wt2[0:4][0:gsp])
  //#pragma acc enter data create(grid3[0:4][0:gsp6],wt3[0:4][0:gsp])
  //#pragma acc enter data create(grid4[0:4][0:gsp6],wt4[0:4][0:gsp])
  #pragma acc enter data create(wtt[0:16][0:gsp])

  for (int n=0;n<4;n++)
  {
    float** val1n = val1[n]; float** val2n = val2[n]; float** val3n = val3[n]; float** val4n = val4[n];
    #pragma acc enter data create(val1n[0:iN][0:gs],val2n[0:iN][0:gs],val3n[0:iN][0:gsp],val4n[0:iN][0:gsp])
  }
 #endif

  #pragma acc parallel loop present(gt[0:M2*M2])
  for (int j=0;j<M2*M2;j++)
    gt[j] = 0.f;

  for (int m=0;m<natoms;m++)
  {
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];
    float Z1 = (float)atno[m];
    float A1 = coords[3*m+0]; float B1 = coords[3*m+1]; float C1 = coords[3*m+2];

    float* grid1a = grid1[0]; float* grid2a = grid2[0];
    float* wt1a = wt1[0]; float* wt2a = wt2[0];

    generate_central_grid_2(grid1a,wt1a,Z1,nrad,  nang,ang_g,ang_w);
    generate_central_grid_2(grid2a,wt2a,Z1,nrad+1,nang,ang_g,ang_w);

    float** val1a = val1[0]; float** val2a = val2[0];
    float** val3a = val3[0]; float** val4a = val4[0];

   #pragma acc parallel loop collapse(2) present(val1a[0:iN][0:gs])
    for (int ii1=0;ii1<s2-s1;ii1++)
    for (int j=0;j<gs;j++)
      val1a[ii1][j] = 1.f;
   #pragma acc parallel loop collapse(2) present(val2a[0:iN][0:gs])
    for (int ii2=0;ii2<s2-s1;ii2++)
    for (int j=0;j<gs;j++)
      val2a[ii2][j] = 1.f;
   #pragma acc parallel loop collapse(2) present(val3a[0:iN][0:gsp])
    for (int ii3=0;ii3<s2-s1;ii3++)
    for (int j=0;j<gsp;j++)
      val3a[ii3][j] = 1.f;
   #pragma acc parallel loop collapse(2) present(val4a[0:iN][0:gsp])
    for (int ii4=0;ii4<s2-s1;ii4++)
    for (int j=0;j<gsp;j++)
      val4a[ii4][j] = 1.f;

   //first compute single atom ints
    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1-s1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

      eval_sh(ii1,gs,grid1a,val1a[ii1],n1,l1,m1,zeta1);
    } //loop i1

    for (int i2=s1;i2<s2;i2++)
    {
      int ii2 = i2-s1;

      vector<double> basis2 = basis[i2];
      int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

      eval_sh(ii2,gs,grid1a,val2a[ii2],n2,l2,m2,zeta2);
    } //loop i2

    for (int i3=s1;i3<s2;i3++)
    {
      int ii3 = i3-s1;

      vector<double> basis3 = basis[i3];
      int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

      eval_sh(ii3,gsp,grid2a,val3a[ii3],n3,l3,m3,zeta3);
    } //loop i3

    for (int i4=s1;i4<s2;i4++)
    {
      int ii4 = i4-s1;

      vector<double> basis4 = basis[i4];
      int n4 = basis4[0]; int l4 = basis4[1]; int m4 = basis4[2]; double zeta4 = basis4[3];

      eval_sh(ii4,gsp,grid2a,val4a[ii4],n4,l4,m4,zeta4);
    } //loop i4


    reduce_4c_1(s1,s2,s1,s2,gs,gsp,M,grid1a,grid2a,val1a,val2a,val3a,val4a,wt1a,wt2a,gt);
    collect_4c_1(s1,s2,s1,s2,0,gs,gsp,M,N,gt,g);


   //two-atom ints
    for (int n=0;n<m;n++)
    {
      int s3 = 0; if (n>0) s3 = n2i[n-1]; int s4 = n2i[n];

      float Z2 = (float)atno[n];
      float A2 = coords[3*n+0]; float B2 = coords[3*n+1]; float C2 = coords[3*n+2];
      float A12 = A2-A1; float B12 = B2-B1; float C12 = C2-C1;

      float* grid1b = grid1[1]; float* grid2b = grid2[1];
      float* grid1d = grid1[3]; float* grid2d = grid2[3];
      float* wt1b = wt1[1]; float* wt2b = wt2[1];

      float** val1b = val1[1]; float** val2b = val2[1];
      float** val3b = val3[1]; float** val4b = val4[1];

     //grid1a at 0,0,0 now has r1 at 3, r2 at 4
      add_r2_to_grid(gs, grid1a,A12,B12,C12);
      add_r2_to_grid(gsp,grid2a,A12,B12,C12);

      generate_central_grid_2(grid1b,wt1b,Z2,nrad,  nang,ang_g,ang_w);
      generate_central_grid_2(grid2b,wt2b,Z2,nrad+1,nang,ang_g,ang_w);

      copy_grid(gs, grid1d,grid1b); //centered at atom b
      copy_grid(gsp,grid2d,grid2b); //centered at atom b

      recenter_grid(gs, grid1b,A12,B12,C12); //r1 points to self, r2 points at first atom
      recenter_grid(gsp,grid2b,A12,B12,C12);

      float* wtt1a = wtt[0]; float* wtt1b = wtt[1];
      float* wtt2a = wtt[4]; float* wtt2b = wtt[5];
      acc_copyf(gs, wtt1a,wt1a); acc_copyf(gs, wtt1b,wt1b);
      acc_copyf(gsp,wtt2a,wt2a); acc_copyf(gsp,wtt2b,wt2b);
      becke_weight_2c(gs, grid1a,wtt1a,grid1b,wtt1b,Z1,Z2,A12,B12,C12);
      becke_weight_2c(gsp,grid2a,wtt2a,grid2b,wtt2b,Z1,Z2,A12,B12,C12);

      //eliminate_small_wt(estart, gs, wtt1a);
      //eliminate_small_wt(estart, gs, wtt1b);
      //eliminate_small_wt(estart2,gsp,wtt2a);
      //eliminate_small_wt(estart2,gsp,wtt2b);

      add_r1_to_grid(gs, grid1b,0.,0.,0.); //center is on atom a
      add_r1_to_grid(gsp,grid2b,0.,0.,0.); //center is on atom a

      float* grid1c = grid1[2]; float* grid2c = grid2[2];

      copy_grid(gs, grid1c,grid1a);
      copy_grid(gsp,grid2c,grid2a);

      recenter_grid_zero(gs, grid1c,-A12,-B12,-C12);
      recenter_grid_zero(gsp,grid2c,-A12,-B12,-C12);


    //1122 (s12,s12,s34,s34)
     #pragma acc parallel loop collapse(2) present(val1b[0:iN][0:gs])
      for (int ii1=0;ii1<s2-s1;ii1++)
      for (int j=0;j<gs;j++)
        val1b[ii1][j] = 1.f;

     #pragma acc parallel loop collapse(2) present(val2b[0:iN][0:gs])
      for (int ii2=0;ii2<s2-s1;ii2++)
      for (int j=0;j<gs;j++)
        val2b[ii2][j] = 1.f;

     #pragma acc parallel loop collapse(2) present(val3b[0:iN][0:gsp])
      for (int ii3=0;ii3<s4-s3;ii3++)
      for (int j=0;j<gsp;j++)
        val3b[ii3][j] = 1.f;

     #pragma acc parallel loop collapse(2) present(val4b[0:iN][0:gsp])
      for (int ii4=0;ii4<s4-s3;ii4++)
      for (int j=0;j<gsp;j++)
        val4b[ii4][j] = 1.f;

     //i1 on atom m
      for (int i1=s1;i1<s2;i1++)
      {
        int ii1 = i1-s1;
        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_sh(ii1,gs,grid1b,val1b[ii1],n1,l1,m1,zeta1);
      }

     //i2 on atom m
      for (int i2=s1;i2<s2;i2++)
      {
        int ii2 = i2-s1;
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1b,val2b[ii2],n2,l2,m2,zeta2);
      }

     //i3 on atom n
      for (int i3=s3;i3<s4;i3++)
      {
        int ii3 = i3-s3;
        vector<double> basis3 = basis[i3];
        int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];
    
        eval_sh(ii3,gsp,grid2d,val3b[ii3],n3,l3,m3,zeta3);
      }

     //i4 on atom n
      for (int i4=s3;i4<s4;i4++)
      {
        int ii4 = i4-s3;
        vector<double> basis4 = basis[i4];
        int n4 = basis4[0]; int l4 = basis4[1]; int m4 = basis4[2]; double zeta4 = basis4[3];
    
        eval_sh(ii4,gsp,grid2d,val4b[ii4],n4,l4,m4,zeta4);
      }

     #pragma acc parallel loop collapse(2) present(val3a[0:iN][0:gsp])
      for (int ii3=0;ii3<s4-s3;ii3++)
      for (int j=0;j<gsp;j++)
        val3a[ii3][j] = 1.f;

     #pragma acc parallel loop collapse(2) present(val4a[0:iN][0:gsp])
      for (int ii4=0;ii4<s4-s3;ii4++)
      for (int j=0;j<gsp;j++)
        val4a[ii4][j] = 1.f;

     //i3 on atom n
      for (int i3=s3;i3<s4;i3++)
      {
        int ii3 = i3-s3;
        vector<double> basis3 = basis[i3];
        int n3 = basis3[0]; int l3 = basis3[1]; int m3 = basis3[2]; double zeta3 = basis3[3];

        eval_sh(ii3,gsp,grid2c,val3a[ii3],n3,l3,m3,zeta3);
      }

     //i4 on atom n
      for (int i4=s3;i4<s4;i4++)
      {
        int ii4 = i4-s3;
        vector<double> basis4 = basis[i4];
        int n4 = basis4[0]; int l4 = basis4[1]; int m4 = basis4[2]; double zeta4 = basis4[3];

        eval_sh(ii4,gsp,grid2c,val4a[ii4],n4,l4,m4,zeta4);
      }

     //1122
      reduce_4c_1(s1,s2,s3,s4,gs,gsp,M,grid1a,grid2a,val1a,val2a,val3a,val4a,wtt1a,wtt2a,gt);
      reduce_4c_1(s1,s2,s3,s4,gs,gsp,M,grid1a,grid2b,val1a,val2a,val3b,val4b,wtt1a,wtt2b,gt);
      reduce_4c_1(s1,s2,s3,s4,gs,gsp,M,grid1b,grid2a,val1b,val2b,val3a,val4a,wtt1b,wtt2a,gt);
      reduce_4c_1(s1,s2,s3,s4,gs,gsp,M,grid1b,grid2b,val1b,val2b,val3b,val4b,wtt1b,wtt2b,gt);
      collect_4c_1(s1,s2,s3,s4,1,gs,gsp,M,N,gt,g);



    //1121
     #pragma acc parallel loop collapse(2) present(val4a[0:iN][0:gsp])
      for (int ii4=0;ii4<s2-s1;ii4++)
      for (int j=0;j<gsp;j++)
        val4a[ii4][j] = 1.f;

     //i4 on atom m
      for (int i4=s1;i4<s2;i4++)
      {
        int ii4 = i4-s1;
        vector<double> basis4 = basis[i4];
        int n4 = basis4[0]; int l4 = basis4[1]; int m4 = basis4[2]; double zeta4 = basis4[3];

        eval_sh(ii4,gsp,grid2a,val4a[ii4],n4,l4,m4,zeta4);
      }

     #pragma acc parallel loop collapse(2) present(val4b[0:iN][0:gsp])
      for (int ii4=0;ii4<s2-s1;ii4++)
      for (int j=0;j<gsp;j++)
        val4b[ii4][j] = 1.f;

     //i2 on atom m
      for (int i4=s1;i4<s2;i4++)
      {
        int ii4 = i4-s1;
        vector<double> basis4 = basis[i4];
        int n4 = basis4[0]; int l4 = basis4[1]; int m4 = basis4[2]; double zeta4 = basis4[3];

        eval_sh(ii4,gsp,grid2b,val4b[ii4],n4,l4,m4,zeta4);
      }

     //1121
      reduce_4c_1b(s1,s2,s3,s4,gs,gsp,M,grid1a,grid2a,val1a,val2a,val3a,val4a,wtt1a,wtt2a,gt);
      reduce_4c_1b(s1,s2,s3,s4,gs,gsp,M,grid1a,grid2b,val1a,val2a,val3b,val4b,wtt1a,wtt2b,gt);
      reduce_4c_1b(s1,s2,s3,s4,gs,gsp,M,grid1b,grid2a,val1b,val2b,val3a,val4a,wtt1b,wtt2a,gt);
      reduce_4c_1b(s1,s2,s3,s4,gs,gsp,M,grid1b,grid2b,val1b,val2b,val3b,val4b,wtt1b,wtt2b,gt);
      collect_4c_1b(s1,s2,s3,s4,gs,gsp,M,N,gt,g);



    //1221
     #pragma acc parallel loop collapse(2) present(val2a[0:iN][0:gs])
      for (int ii2=0;ii2<s4-s3;ii2++)
      for (int j=0;j<gs;j++)
        val2a[ii2][j] = 1.f;

     //i2 on atom n
      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1c,val2a[ii2],n2,l2,m2,zeta2);
      }

     #pragma acc parallel loop collapse(2) present(val2b[0:iN][0:gs])
      for (int ii2=0;ii2<s4-s3;ii2++)
      for (int j=0;j<gs;j++)
        val2b[ii2][j] = 1.f;

     //i2 on atom n
      for (int i2=s3;i2<s4;i2++)
      {
        int ii2 = i2-s3;
        vector<double> basis2 = basis[i2];
        int n2 = basis2[0]; int l2 = basis2[1]; int m2 = basis2[2]; double zeta2 = basis2[3];

        eval_sh(ii2,gs,grid1d,val2b[ii2],n2,l2,m2,zeta2);
      }

     //1221
      reduce_4c_1c(s1,s2,s3,s4,gs,gsp,M,grid1a,grid2a,val1a,val2a,val3a,val4a,wtt1a,wtt2a,gt);
      reduce_4c_1c(s1,s2,s3,s4,gs,gsp,M,grid1a,grid2b,val1a,val2a,val3b,val4b,wtt1a,wtt2b,gt);
      reduce_4c_1c(s1,s2,s3,s4,gs,gsp,M,grid1b,grid2a,val1b,val2b,val3a,val4a,wtt1b,wtt2a,gt);
      reduce_4c_1c(s1,s2,s3,s4,gs,gsp,M,grid1b,grid2b,val1b,val2b,val3b,val4b,wtt1b,wtt2b,gt);
      collect_4c_1c(s1,s2,s3,s4,gs,gsp,M,N,gt,g);


    //2221
     #pragma acc parallel loop collapse(2) present(val1a[0:iN][0:gs])
      for (int ii1=0;ii1<s4-s3;ii1++)
      for (int j=0;j<gs;j++)
        val1a[ii1][j] = 1.f;

     //i1 on atom n
      for (int i1=s3;i1<s4;i1++)
      {
        int ii1 = i1-s3;
        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_sh(ii1,gs,grid1c,val1a[ii1],n1,l1,m1,zeta1);
      }

     #pragma acc parallel loop collapse(2) present(val1b[0:iN][0:gs])
      for (int ii1=0;ii1<s4-s3;ii1++)
      for (int j=0;j<gs;j++)
        val1b[ii1][j] = 1.f;

     //i1 on atom n
      for (int i1=s3;i1<s4;i1++)
      {
        int ii1 = i1-s3;
        vector<double> basis1 = basis[i1];
        int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; double zeta1 = basis1[3];

        eval_sh(ii1,gs,grid1d,val1b[ii1],n1,l1,m1,zeta1);
      }

      reduce_4c_1d(s1,s2,s3,s4,gs,gsp,M,grid1a,grid2a,val1a,val2a,val3a,val4a,wtt1a,wtt2a,gt);
      reduce_4c_1d(s1,s2,s3,s4,gs,gsp,M,grid1a,grid2b,val1a,val2a,val3b,val4b,wtt1a,wtt2b,gt);
      reduce_4c_1d(s1,s2,s3,s4,gs,gsp,M,grid1b,grid2a,val1b,val2b,val3a,val4a,wtt1b,wtt2a,gt);
      reduce_4c_1d(s1,s2,s3,s4,gs,gsp,M,grid1b,grid2b,val1b,val2b,val3b,val4b,wtt1b,wtt2b,gt);
      collect_4c_1d(s1,s2,s3,s4,gs,gsp,M,N,gt,g);


    } //loop n over second atom

  } //loop m over natoms

 //normalization
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
  for (int l=0;l<N;l++)
    g[i*N3+j*N2+k*N+l] *= basis[i][4]*basis[j][4]*basis[k][4]*basis[l][4];

 //show the symmetry
  if (prl>1)
  for (int p=0;p<N;p++)
  for (int q=0;q<=p;q++)
  for (int r=0;r<=p;r++)
  for (int s=0;s<=r;s++)
  {
    double v1 = g[p*N3+q*N2+r*N+s];
    double v2 = g[p*N3+q*N2+s*N+r];
    double v3 = g[q*N3+p*N2+r*N+s];
    double v4 = g[q*N3+p*N2+s*N+r];

    double v5 = g[r*N3+s*N2+p*N+q];
    double v6 = g[s*N3+r*N2+p*N+q];
    double v7 = g[r*N3+s*N2+q*N+p];
    double v8 = g[s*N3+r*N2+q*N+p];

    printf("   pqrs: %i %i %i %i: %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f \n",p,q,r,s,v1,v2,v3,v4,v5,v6,v7,v8);
  }

  if ((prl>0 && N<10) || prl>1)
  {
    printf("\n g: \n");
    for (int i=0;i<N;i++)
    for (int j=0;j<N;j++)
    {
      printf(" ij: %2i %2i \n",i,j);
      print_square(N,&g[i*N3+j*N2]);
    }
  }


#if USE_ACC
  #pragma acc exit data delete(gt[0:M2*M2])

  #pragma acc exit data delete(grid1[0:4][0:gs6],wt1[0:4][0:gs])
  #pragma acc exit data delete(grid2[0:4][0:gsp6],wt2[0:4][0:gsp])
  //#pragma acc exit data delete(grid3[0:4][0:gsp6],wt3[0:4][0:gsp])
  //#pragma acc exit data delete(grid4[0:4][0:gsp6],wt4[0:4][0:gsp])
  for (int n=0;n<4;n++)
  {
    float** val1n = val1[n]; float** val2n = val2[n]; float** val3n = val3[n]; float** val4n = val4[n];
    #pragma acc exit data delete(val1n[0:iN][0:gs],val2n[0:iN][0:gs],val3n[0:iN][0:gsp],val4n[0:iN][0:gsp])
  }

  #pragma acc exit data delete(wtt[0:16][0:gsp])
  #pragma acc exit data delete(ang_g[0:3*nang],ang_w[0:nang])
  #pragma acc exit data delete(n2i[0:natoms])
#endif


  delete [] gt;
  delete [] n2i;

  delete [] ang_g;
  delete [] ang_w;

  for (int i=0;i<4;i++) 
  {
    for (int j=0;j<iN;j++)
      delete [] val1[i][j];
    delete [] val1[i];
    for (int j=0;j<iN;j++)
      delete [] val2[i][j];
    delete [] val2[i];
    for (int j=0;j<iN;j++)
      delete [] val3[i][j];
    delete [] val3[i];
    for (int j=0;j<iN;j++)
      delete [] val4[i][j];
    delete [] val4[i];
  }
 
  for (int i=0;i<4;i++)
  {
    delete [] grid1[i];
    delete [] grid2[i];
    delete [] wt1[i];
    delete [] wt2[i];
  }
  delete [] grid1;
  delete [] grid2;
  delete [] wt1;
  delete [] wt2;

  return;
}
