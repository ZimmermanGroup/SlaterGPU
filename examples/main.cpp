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

#include "fp_def.h"

using namespace std;

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
  int * atno = new int[100](); // 
  double * coords;
  vector< vector< double > > basis;
  vector< vector< double > > basis_aux;
  double Enn = 0;

  // seting up
  int natoms = initialize(false,basis,basis_aux,atno,coords,charge,nup,Enn,1);
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

    int prl = 2;
    #pragma acc enter data create(A[0:Naux2],C[0:N2a],S[0:N2],En[0:N2],T[0:N2],pVp[0:N2])
    #pragma acc enter data create(V[0:nc],dV[0:3*nc],Pao[0:N2],coordsc[0:3*nc])
    #pragma acc enter data create(g[0:N2*N2])
    printf("1e ints: %d\n2c2e ints: %d\n3c3e ints: %d\n",N2, Naux2, N2a);

    //core integrals
    auto t1 = chrono::high_resolution_clock::now();
    printf("Computing ST:\n");
    compute_ST(natoms, atno, coordsf, basis, nrad, size_ang, ang_g, ang_w, S, T, prl);

    auto t2 = chrono::high_resolution_clock::now();
    printf("Computing 2c_v2:\n");
    compute_all_2c_v2(0,natoms,atno,coordsf,basis_aux,nrad,size_ang,ang_g,ang_w,A,prl);

    auto t3 = chrono::high_resolution_clock::now();
    printf("Computing Enp:\n");
    compute_Enp(natoms,atno,coordsf,basis,nrad,size_ang,ang_g,ang_w,En,pVp,prl);

    auto t4 = chrono::high_resolution_clock::now();
    printf("Computing 3c_para:\n");
    compute_all_3c_para(ngpu,0,natoms,atno,coordsf,basis,basis_aux,nrad,size_ang,ang_g,ang_w,C,prl);

    auto t5 = chrono::high_resolution_clock::now();

    //Nate's Edits
    printf("Computing 3c_v2:\n");
    compute_all_3c_v2(0,natoms,atno,coordsf,basis,basis_aux,nrad,size_ang,ang_g,ang_w,C,prl);
    
    auto t6 = chrono::high_resolution_clock::now();
    printf("Computing VdV:\n");
    compute_VdV(natoms,atno,coordsf,basis,nrad,size_ang,ang_g,ang_w,nc,coordsc,Pao,V,dV,prl);
    
    auto t7 = chrono::high_resolution_clock::now();
    printf("Computing Enp_para:\n");
    compute_Enp_para(ngpu,natoms,atno,coordsf,basis,nrad,size_ang,ang_g,ang_w,En,pVp,prl);

    auto t8 = chrono::high_resolution_clock::now();
    printf("Computing 4c_v2:\n");
    compute_all_4c_v2(natoms,atno,coordsf,basis,nrad,size_ang,ang_g,ang_w,g,prl);
    
    auto t9 = chrono::high_resolution_clock::now();

    auto elapsed12 = chrono::duration_cast<chrono::nanoseconds>(t2-t1).count();
    auto elapsed23 = chrono::duration_cast<chrono::nanoseconds>(t3-t2).count();
    auto elapsed34 = chrono::duration_cast<chrono::nanoseconds>(t4-t3).count();
    auto elapsed45 = chrono::duration_cast<chrono::nanoseconds>(t5-t4).count();
    auto elapsed15 = chrono::duration_cast<chrono::nanoseconds>(t5-t1).count();

    auto elapsed56 = chrono::duration_cast<chrono::nanoseconds>(t6-t5).count();
    auto elapsed67 = chrono::duration_cast<chrono::nanoseconds>(t7-t6).count();    
    auto elapsed78 = chrono::duration_cast<chrono::nanoseconds>(t8-t7).count();
    auto elapsed89 = chrono::duration_cast<chrono::nanoseconds>(t9-t8).count();

    printf("Integral ST   time: %5.3e s\n",(double)elapsed12/1.e9);
    printf("Integral 2c2e time: %5.3e s\n",(double)elapsed23/1.e9);
    printf("Integral Vne  time: %5.3e s\n",(double)elapsed34/1.e9);
    printf("Integral 3c2e time: %5.3e s\n",(double)elapsed45/1.e9);
    printf("Integral tot  time: %5.3e s\n",(double)elapsed15/1.e9);
    
    printf("Integral 3c2e (v2)  time: %5.3e s\n",(double)elapsed56/1.e9);
    printf("Integral VdV        time: %5.3e s\n",(double)elapsed67/1.e9);
    printf("Integral Vne (para) time: %5.3e s\n",(double)elapsed78/1.e9);
    printf("Integral 4c (v2)    time: %5.3e s\n",(double)elapsed89/1.e9);

    #pragma acc exit data copyout(A[0:Naux2],C[0:N2a],S[0:N2],En[0:N2],T[0:N2],pVp[0:N2])
    #pragma acc exit data copyout(V[0:nc],dV[0:3*nc],Pao[0:N2],coordsc[0:3*nc],g[0:N2*N2])

    write_S_En_T(N,S,En,T);
    write_square(N,pVp,"pVp",2);

    // int mprint = min(N,10);
    // for (int i = 0; i < mprint; i++) {
    //   for (int j = 0; j < mprint; j++) {
    //     printf("%6.4f ",En[i*N+j]);
    //   }
    //   printf("\n");
    // }

    // for (int i=0;i<nc;i++)
    //   printf(" V[%i]: %8.5f \n",i,V[i]);
    // printf("\n");
    // printf(" dV: \n");
    // for (int i=0;i<nc;i++)
    //   printf(" %8.5f %8.5f %8.5f \n",dV[3*i+0],dV[3*i+1],dV[3*i+2]);
    // printf("\n");
   
    delete [] A, Anorm, C, S, En, T, pVp;
    delete [] V, dV, Pao, coordsc;

    delete [] g;
  }

  delete [] ang_g, ang_w;
  return 0;

}

