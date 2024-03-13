#if USE_ACC
#include "accel.h"
#endif
#include <omp.h>
#include <cstdio>
#include <vector>
#include <string>
#include <chrono>

#include "read.h"
#include "lebedev2.h"
#include "integrals.h"
#include "grid_util.h"
#include "rel.h"

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
  FP2 * coords;
  vector< vector< FP2 > > basis;
  vector< vector< FP2 > > basis_aux;
  FP2 Enn = 0;

  // seting up
  int natoms = initialize(false,basis,basis_aux,atno,coords,charge,nup,Enn,1);
  FP1 * coordsf = new FP1[3*natoms];
  for (int i = 0; i < 3*natoms; i++) {
    coordsf[i] = coords[i];
  }
  
  FP2 SOL = 137.03599967994;
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

  FP2 * ang_g = new FP2[3*size_ang];
  FP2 * ang_w = new FP2[size_ang];
  get_angular_grid(size_ang,ang_g,ang_w);

  // Doing work
  if (true) {
    FP2 * A = new FP2[Naux2]();
    FP2 * Anorm = new FP2[Naux];
    FP2 * C = new FP2[N2a]();
    FP2 * S = new FP2[N2]();
    FP2 * En = new FP2[N2];
    FP2 * T = new FP2[N2];
    FP2 * pVp = new FP2[N2];
    FP2 * hcore = new FP2[N2];
    int prl = 0;
    #pragma acc enter data \
      create(A[0:Naux2],C[0:N2a],S[0:N2],En[0:N2],T[0:N2],pVp[0:N2],hcore[0:N2])
    printf("1e ints: %d\n2c2e ints: %d\n3c3e ints: %d\n",N2, Naux2, N2a);
    compute_ST(natoms, atno, coordsf, basis, nrad, size_ang, ang_g, ang_w, S, T, prl);
    compute_all_2c_v2(natoms,atno,coordsf,basis_aux,nrad,size_ang,ang_g,ang_w,A,prl);
    compute_Enp(natoms,atno,coordsf,basis,nrad,size_ang,ang_g,ang_w,En,pVp,prl);
    compute_all_3c_para2(ngpu,natoms,atno,coordsf,basis,basis_aux,nrad,size_ang,ang_g,ang_w,C,prl);
    compute_hcore(N, T, En, pVp, S, SOL, hcore, 0);
    #pragma acc exit data \
      copyout(A[0:Naux2],C[0:N2a],S[0:N2],En[0:N2],T[0:N2],pVp[0:N2],hcore[0:N2])

    printf("h: \n");
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        printf("%10.8f ",En[i*N+j]+T[i*N+j]);
      }
      printf("\n");
    }

    delete [] A, Anorm, C, S, En, T, pVp, hcore;
  }

  delete [] ang_g, ang_w;
  return 0;
}
