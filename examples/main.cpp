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
    int prl = 0;
    #pragma acc enter data \
      create(A[0:Naux2],C[0:N2a],S[0:N2],En[0:N2],T[0:N2],pVp[0:N2])
    printf("1e ints: %d\n2c2e ints: %d\n3c3e ints: %d\n",N2, Naux2, N2a);
    auto t1 = chrono::high_resolution_clock::now();
    compute_ST(natoms, atno, coordsf, basis, nrad, size_ang, ang_g, ang_w, S, T, prl);
    auto t2 = chrono::high_resolution_clock::now();
    compute_all_2c_v2(0,natoms,atno,coordsf,basis_aux,nrad,size_ang,ang_g,ang_w,A,prl);
    auto t3 = chrono::high_resolution_clock::now();
    compute_Enp(natoms,atno,coordsf,basis,nrad,size_ang,ang_g,ang_w,En,pVp,prl);
    auto t4 = chrono::high_resolution_clock::now();
    compute_all_3c_para(ngpu,0,natoms,atno,coordsf,basis,basis_aux,nrad,size_ang,ang_g,ang_w,C,prl);
    auto t5 = chrono::high_resolution_clock::now();
    auto elapsed12 = chrono::duration_cast<chrono::nanoseconds>(t2-t1).count();
    auto elapsed23 = chrono::duration_cast<chrono::nanoseconds>(t3-t2).count();
    auto elapsed34 = chrono::duration_cast<chrono::nanoseconds>(t4-t3).count();
    auto elapsed45 = chrono::duration_cast<chrono::nanoseconds>(t5-t4).count();
    auto elapsed15 = chrono::duration_cast<chrono::nanoseconds>(t5-t1).count();
    printf("Integral ST   time: %5.3e s\n",(double)elapsed12/1.e9);
    printf("Integral 2c2e time: %5.3e s\n",(double)elapsed23/1.e9);
    printf("Integral Vne  time: %5.3e s\n",(double)elapsed34/1.e9);
    printf("Integral 3c2e time: %5.3e s\n",(double)elapsed45/1.e9);
    printf("Integral tot  time: %5.3e s\n",(double)elapsed15/1.e9);
    #pragma acc exit data \
      copyout(A[0:Naux2],C[0:N2a],S[0:N2],En[0:N2],T[0:N2],pVp[0:N2])
    
    if (prl > 0) printf("Printing Standard Integral Files:\n");
    write_S_En_T(N,S,En,T);
    write_square(Naux,A,"A",2);
    write_square(N,pVp,"pVp",2);
    write_C(Naux, N2, C);

    int mprint = min(N,10);
    for (int i = 0; i < mprint; i++) {
      for (int j = 0; j < mprint; j++) {
        printf("%6.4f ",En[i*N+j]);
      }
      printf("\n");
    }
    delete [] A, Anorm, C, S, En, T, pVp;
  }

  delete [] ang_g, ang_w;
  return 0;
}
