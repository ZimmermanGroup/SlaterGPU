#include "cintwrapper.h"
#include "elements.h"

#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <math.h>
#include "omp.h"

#if USEMKL
#include <mkl_lapacke.h>
#include <mkl_cblas.h>
#else
#include <lapacke.h>
#include <cblas.h>
#endif

using namespace std;

//extern "C" void dgetri(int *N, double **A, int *lda, int *IPIV,
//                       double *WORK, int *lwork, int *INFO);
//extern "C" void dgetrf(int *M, int *N, double **A, int *lda, int *IPIV,
//                       int *INFO);
//extern "C" void dsyev(char *jobz, char *uplo, int *n, double *a, int *lda,
//                      double *w, double * work, int *lwork, int* info);

//bool BT::DO_CART = true;
bool BT::DO_CART = false;

int binom(int n, int k) {
  return (tgamma(n+1))/(tgamma(n+1-k)*tgamma(k+1));
}

int calc_di(int i, int *bas) {
  int idx_i = 0;
  for (int j = 0; j < i; j++) {
    idx_i += CINTcgto_cart(j, bas);
  }
  return idx_i;
}

void get_overlap(double * overlap, int N,
                 int natm, int nbas, int nenv,
                 int *atm, int* bas, double *env) {
  int idx_i = 0;
  int di, dj;
  int shls[2];

  if (BT::DO_CART) {
    for (int i = 0; i < nbas; i++) {
      int idx_j = 0;
      shls[0] = i;
      di = CINTcgto_cart(i, bas);
      for (int j = 0; j <= i; j++) {
        dj = CINTcgto_cart(j, bas);
        shls[1] = j;
        double *buf = new double[di*dj];
        cint1e_ovlp_cart(buf,shls,atm,natm,bas,nbas,env);

        for (int j1 = 0; j1 < dj; j1++) {
          for (int i1 = 0; i1 < di; i1++) {
            int oi = i1 + idx_i;
            int oj = j1 + idx_j;
            overlap[oi * N + oj] = overlap[oj * N + oi] = buf[j1*di + i1];
          } // for i1
        } // for j1
        delete [] buf;
        idx_j += dj;
      } // for j
      idx_i += di;
    } // for i
  } // if BT::DO_CART
  else {
    for (int i = 0; i < nbas; i++) {
      int idx_j = 0;
      shls[0] = i;
      di = CINTcgto_spheric(i, bas);
      for (int j = 0; j <= i; j++) {
        dj = CINTcgto_spheric(j, bas);
        shls[1] = j;
        double *buf = new double[di*dj];
        cint1e_ovlp_sph(buf,shls,atm,natm,bas,nbas,env);

        for (int j1 = 0; j1 < dj; j1++) {
          for (int i1 = 0; i1 < di; i1++) {
            int oi = i1 + idx_i;
            int oj = j1 + idx_j;
            overlap[oi * N + oj] = overlap[oj * N + oi] = buf[j1*di + i1];
          } // for i1
        } // for j1
        delete [] buf;
        idx_j += dj;
      } // for j
      idx_i += di;
    }// for i
  } // else
}

void get_overlap_ri(double * overlap, int Naux,
                    int natm, int nbas, int nbas_ri, int nenv,
                    int *atm, int* bas, double *env) {
  int idx_i = 0;
  int di, dj;
  int shls[2];

  if (BT::DO_CART) {
    for (int i = 0; i < nbas_ri; i++) {
      int idx_j = 0;
      shls[0] = i + nbas;
      di = CINTcgto_cart(i + nbas, bas);
      for (int j = 0; j <= i; j++) {
        dj = CINTcgto_cart(j + nbas, bas);
        shls[1] = j + nbas;
        double *buf = new double[di*dj];
        cint1e_ovlp_cart(buf,shls,atm,natm,bas,nbas,env);

        for (int j1 = 0; j1 < dj; j1++) {
          for (int i1 = 0; i1 < di; i1++) {
            int oi = i1 + idx_i;
            int oj = j1 + idx_j;
            overlap[oi * Naux + oj] = overlap[oj * Naux + oi] = buf[j1*di + i1];
          } // for i1
        } // for j1
        delete [] buf;
        idx_j += dj;
      } // for j
      idx_i += di;
    }// for i
  } // if BT::DO_CART
  else {
    for (int i = 0; i < nbas_ri; i++) {
      int idx_j = 0;
      shls[0] = i + nbas;
      di = CINTcgto_spheric(i + nbas, bas);
      for (int j = 0; j <= i; j++) {
        dj = CINTcgto_spheric(j + nbas, bas);
        shls[1] = j + nbas;
        double *buf = new double[di*dj];
        cint1e_ovlp_sph(buf,shls,atm,natm,bas,nbas,env);

        for (int j1 = 0; j1 < dj; j1++) {
          for (int i1 = 0; i1 < di; i1++) {
            int oi = i1 + idx_i;
            int oj = j1 + idx_j;
            overlap[oi * Naux + oj] = overlap[oj * Naux + oi] = buf[j1*di + i1];
          } // for i1
        } // for j1
        delete [] buf;
        idx_j += dj;
      } // for j
      idx_i += di;
    }// for i
  } // else
}

void get_hcore(double *hcore, int N,
               int natm, int nbas, int nenv,
               int *atm, int *bas, double *env) {
  int idx_i = 0;
  int di, dj;
  int shls[2];

  double *tmp = new double[N * N];
  
  if (BT::DO_CART) {
    for (int i = 0; i < nbas; i++) {
      int idx_j = 0;
      shls[0] = i;
      di = CINTcgto_cart(i, bas);
      for (int j = 0; j <= i; j++) {
        dj = CINTcgto_cart(j, bas);
        shls[1] = j;
        double *buf = new double[di*dj];
        cint1e_nuc_cart(buf,shls,atm,natm,bas,nbas,env);

        for (int j1 = 0; j1 < dj; j1++) {
          for (int i1 = 0; i1 < di; i1++) {
            int oi = i1 + idx_i;
            int oj = j1 + idx_j;
            tmp[oi * N + oj] = tmp[oj * N + oi] = buf[j1*di + i1];
          } // for i1
        } // for j1
        delete [] buf;
        idx_j += dj;
      } // for j
      idx_i += di;
    }// for i

    for (int i = 0; i < N * N; i++) {
      hcore[i] = tmp[i];
    }

    idx_i = 0;
    for (int i = 0; i < nbas; i++) {
      int idx_j = 0;
      shls[0] = i;
      di = CINTcgto_cart(i, bas);
      for (int j = 0; j <= i; j++) {
        dj = CINTcgto_cart(j, bas);
        shls[1] = j;
        double *buf = new double[di*dj];
        cint1e_kin_cart(buf,shls,atm,natm,bas,nbas,env);

        for (int j1 = 0; j1 < dj; j1++) {
          for (int i1 = 0; i1 < di; i1++) {
            int oi = i1 + idx_i;
            int oj = j1 + idx_j;
            tmp[oi * N + oj] = tmp[oj * N + oi] = buf[j1*di + i1];
          } // for i1
        } // for j1
        delete [] buf;
        idx_j += dj;
      } // for j
      idx_i += di;
    }// for i
  } // if BT::DO_CART
  else {
    for (int i = 0; i < nbas; i++) {
      int idx_j = 0;
      shls[0] = i;
      di = CINTcgto_spheric(i, bas);
      for (int j = 0; j <= i; j++) {
        dj = CINTcgto_spheric(j, bas);
        shls[1] = j;
        double *buf = new double[di*dj];
        cint1e_nuc_sph(buf,shls,atm,natm,bas,nbas,env);

        for (int j1 = 0; j1 < dj; j1++) {
          for (int i1 = 0; i1 < di; i1++) {
            int oi = i1 + idx_i;
            int oj = j1 + idx_j;
            tmp[oi * N + oj] = tmp[oj * N + oi] = buf[j1*di + i1];
          } // for i1
        } // for j1
        delete [] buf;
        idx_j += dj;
      } // for j
      idx_i += di;
    }// for i

    for (int i = 0; i < N * N; i++) {
      hcore[i] = tmp[i];
    }

    idx_i = 0;
    for (int i = 0; i < nbas; i++) {
      int idx_j = 0;
      shls[0] = i;
      di = CINTcgto_spheric(i, bas);
      for (int j = 0; j <= i; j++) {
        dj = CINTcgto_spheric(j, bas);
        shls[1] = j;
        double *buf = new double[di*dj];
        cint1e_kin_sph(buf,shls,atm,natm,bas,nbas,env);

        for (int j1 = 0; j1 < dj; j1++) {
          for (int i1 = 0; i1 < di; i1++) {
            int oi = i1 + idx_i;
            int oj = j1 + idx_j;
            tmp[oi * N + oj] = tmp[oj * N + oi] = buf[j1*di + i1];
          } // for i1
        } // for j1
        delete [] buf;
        idx_j += dj;
      } // for j
      idx_i += di;
    }// for i
  }

  for (int i = 0; i < N * N; i++) {
    hcore[i] += tmp[i];
  }

  delete [] tmp;
}

void gen_eri(double **eri, int N, 
             int natm, int nbas, int nenv,
             int *atm, int *bas, double *env) {
  int N2 = N*N;

  int idxi = 0;
  int idxj = 0;
  int idxk = 0;
  int idxl = 0;
  int di = 0;
  int dj = 0;
  int dk = 0;
  int dl = 0;
  CINTOpt *no_opt = NULL;
  if (BT::DO_CART) {
    for (int i = 0; i < nbas; i++) {
      di = CINTcgto_cart(i,bas);
      for (int j = 0; j < nbas; j++) {
        dj = CINTcgto_cart(j,bas);
        for (int k = 0; k < nbas; k++) {
          dk = CINTcgto_cart(k,bas);
          for (int l = 0; l < nbas; l++) {
            dl = CINTcgto_cart(l,bas);
            double * buf = new double[di*dj*dk*dl];
            int shls[4];
            shls[0] = i;
            shls[1] = j;
            shls[2] = k;
            shls[3] = l;
            cint2e_cart(buf, shls, atm, natm, bas, nbas, env, no_opt);

            for (int l1 = 0; l1 < dl; l1++) {
              for (int k1 = 0; k1 < dk; k1++) {
                for (int j1 = 0; j1 < dj; j1++) {
                  for (int i1 = 0; i1 < di; i1++) {
                    int oi = i1 + idxi;
                    int oj = j1 + idxj;
                    int ok = k1 + idxk;
                    int ol = l1 + idxl;
                    int ij = oi * N + oj;
                    int kl = ok * N + ol;
                    int bdx = l1*(dk*dj*di) + k1*(dj*di) + j1*(di) + i1;
                    // printf("oi %2d oj %2d ok %2d ol %2d ij %6d kl %6d\n",oi, oj, ok, ol,ij,kl);
                    eri[ij][kl] = eri[kl][ij] = buf[bdx];
                  } // for i1
                } // for j1
              } // for k1
            } // for l1
            delete [] buf;
            idxl += dl;
          } // for l
          idxl = 0;
          idxk += dk;
        } // for k
        idxk = 0;
        idxj += dj;
      } // for j
      idxj = 0;
      idxi += di;
    } // for i
  } // if BT::DO_CART
  else {
    for (int i = 0; i < nbas; i++) {
      di = CINTcgto_spheric(i,bas);
      for (int j = 0; j < nbas; j++) {
        dj = CINTcgto_spheric(j,bas);
        for (int k = 0; k < nbas; k++) {
          dk = CINTcgto_spheric(k,bas);
          for (int l = 0; l < nbas; l++) {
            dl = CINTcgto_spheric(l,bas);
            double * buf = new double[di*dj*dk*dl];
            int shls[4];
            shls[0] = i;
            shls[1] = j;
            shls[2] = k;
            shls[3] = l;
            cint2e_sph(buf, shls, atm, natm, bas, nbas, env, no_opt);

            for (int l1 = 0; l1 < dl; l1++) {
              for (int k1 = 0; k1 < dk; k1++) {
                for (int j1 = 0; j1 < dj; j1++) {
                  for (int i1 = 0; i1 < di; i1++) {
                    int oi = i1 + idxi;
                    int oj = j1 + idxj;
                    int ok = k1 + idxk;
                    int ol = l1 + idxl;
                    int ij = oi * N + oj;
                    int kl = ok * N + ol;
                    int bdx = l1*(dk*dj*di) + k1*(dj*di) + j1*(di) + i1;
                    // printf("oi %2d oj %2d ok %2d ol %2d ij %6d kl %6d\n",oi, oj, ok, ol,ij,kl);
                    eri[ij][kl] = eri[kl][ij] = buf[bdx];
                  } // for i1
                } // for j1
              } // for k1
            } // for l1
            delete [] buf;
            idxl += dl;
          } // for l
          idxl = 0;
          idxk += dk;
        } // for k
        idxk = 0;
        idxj += dj;
      } // for j
      idxj = 0;
      idxi += di;
    } // for i
  } // else
}

void gen_ri_int(double *jSaux, double *jSinv, double **jCaux,
                double *jERI3c, double *jERI2c, int Naux, int N,
                int natm, int nbas, int nenv, int nbas_ri,
                int *atm, int *bas, double *env, double **eri) {
  int N2 = N*N;
  int Naux2 = Naux * Naux;

  gen_eri_3c(jERI3c, N, Naux, natm, nbas, nenv, nbas_ri, atm, bas, env);
  gen_eri_2c(jERI2c, Naux, natm, nbas, nenv, nbas_ri, atm, bas, env);

  double * jV = new double[Naux * Naux];
  for (int i = 0; i < Naux; i++) {
    for (int j = 0; j < Naux; j++) {
      jV[i*Naux+j] = jERI2c[i*Naux + j];
    }
  }
  
  double * eval = new double[Naux];
  char v = 'V';
  char u = 'U';
  double * work = new double[Naux];
  int info;
  // dsyev(&v,&u,&Naux,jV, &Naux, eval, work, &Naux, &info);
  // for (int i = 0; i < Naux; i++) {
  //   eval[i] = 1/sqrt(eval[i]);
  // }
//  LAPACKE_dsyev(LAPACK_ROW_MAJOR, v, u, Naux, jV, Naux, work);

  
  int ipiv[Naux];

  double *diag = new double[Naux2]();
  for (int i = 0; i < Naux; i++) {
    diag[i*Naux+i] = 1/sqrt(work[i]);
  }

  double *tmp = new double[Naux2]();
//  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
//              N,N,N,1.0,jV,N,diag,N,0.,tmp,N);
//  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
//              N,N,N,1.0,tmp,N,jV,N,0.,diag,N);
  
  for (int i = 0; i < N2; i++) {
    //can replace with cblas_dgemm later
    for (int mu = 0; mu < Naux; mu++) {
      jCaux[i][mu] = 0.;
      for (int nu = 0; nu < Naux; nu++) {
        jCaux[i][mu] += diag[nu*Naux+mu]*jERI3c[i*Naux+nu];
      }
    }
  }

  delete [] diag;
  delete [] tmp;
  delete [] work;
  delete [] jV;
  delete [] eval;
}

void gen_eri_ri(double *jSaux, double *jSinv, double **jCaux,
                double *jERI3c, double *jERI2c, int Naux, int N,
                int natm, int nbas, int nenv, int nbas_ri,
                int *atm, int *bas, double *env, double **eri) {
  int N2 = N*N;
  int Naux2 = Naux * Naux;
  gen_eri_3c(jERI3c, N, Naux, natm, nbas, nenv, nbas_ri, atm, bas, env);
  gen_eri_2c(jERI2c, Naux, natm, nbas, nenv, nbas_ri, atm, bas, env);
  double * jV = new double[Naux * Naux];
  for (int i = 0; i < Naux; i++) {
    for (int j = 0; j < Naux; j++) {
      jV[i*Naux+j] = jERI2c[i*Naux+j];
    }
  }
  int ipiv[Naux];
//  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, Naux, Naux, jV, Naux, ipiv);
//  LAPACKE_dgetri(LAPACK_ROW_MAJOR, Naux, jV, Naux, ipiv);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      for (int k = 0; k < N; k++) {
        for (int l = 0; l < N; l++) {
          int ij = i*N+j;
          int kl = k*N+l;
          eri[ij][kl] = 0.;
          for (int mu = 0; mu < Naux; mu++) {
            for (int nu = 0; nu < Naux; nu++) {
              eri[ij][kl] += jERI3c[ij*Naux + mu] * jERI3c[kl*Naux + nu] 
                          * jV[mu*Naux+nu];
            }
          }
        }
      }      
    }
  }
  delete [] jV;
}

void gen_eri_3c(double *eri, int N, int Naux,
                int natm, int nbas, int nenv, int nbas_ri,
                int *atm, int *bas, double *env) {
  int N2 = N*N;
  int idxi = 0;
  int idxj = 0;
  int idxmu = 0;
  int di = 0;
  int dj = 0;
  int dmu = 0;
  CINTOpt *no_opt = NULL;
  // printf("bas: %d %d\n",
  //        bas[(nbas+nbas_ri)*BAS_SLOTS+PTR_EXP],
  //        bas[(nbas+nbas_ri)*BAS_SLOTS+PTR_COEFF]);
  // printf("env: %7.4f %7.4f\n",
  //        env[bas[(nbas+nbas_ri)*BAS_SLOTS+PTR_EXP]],
  //        env[bas[(nbas+nbas_ri)*BAS_SLOTS+PTR_COEFF]]);
  // printf("dim: %d\n",CINTcgto_cart(nbas+nbas_ri,bas));
  if (BT::DO_CART) {
    for (int mu = 0; mu < nbas_ri; mu++) {
      dmu = CINTcgto_cart(mu+nbas,bas);
      for (int i = 0; i < nbas; i++) {
        di = CINTcgto_cart(i, bas);
        for (int j = 0; j < nbas; j++) {
          dj = CINTcgto_cart(j, bas);
          double * buf = new double[dmu * di * dj];
          int shls[4];
          shls[0] = mu + nbas;
          shls[1] = nbas + nbas_ri; // last bas has exp 0
          shls[2] = i;
          shls[3] = j;
          cint2e_cart(buf, shls, atm, natm, bas, nbas + nbas_ri + 1, env, no_opt);
          for (int mu1 = 0; mu1 < dmu; mu1++) {
            for (int i1 = 0; i1 < di; i1++) {
              for (int j1 = 0; j1 < dj; j1++) {
                int omu = mu1 + idxmu;
                int oi = i1 + idxi;
                int oj = j1 + idxj;
                int ij = oi * N + oj;
                int ji = oj * N + oi;
                int bdx = j1 * (di * dmu) + i1 * dmu + mu1;
                eri[ij*Naux + omu] = eri[ji*Naux + omu]= buf[bdx];
              } // for j1
            } // for i1
          } // for mu1
          delete [] buf;
          idxj += dj;
        } // for j
        idxj = 0;
        idxi += di;
      } // for i
      idxi = 0;
      idxmu += dmu;
    } // for mu
  } // if BT::DO_CART
  else {
    for (int mu = 0; mu < nbas_ri; mu++) {
      dmu = CINTcgto_spheric(mu+nbas,bas);
      for (int i = 0; i < nbas; i++) {
        di = CINTcgto_spheric(i, bas);
        for (int j = 0; j < nbas; j++) {
          dj = CINTcgto_spheric(j, bas);
          double * buf = new double[dmu * di * dj];
          int shls[4];
          shls[0] = mu + nbas;
          shls[1] = nbas + nbas_ri; // last bas has exp 0
          shls[2] = i;
          shls[3] = j;
          cint2e_sph(buf, shls, atm, natm, bas, nbas + nbas_ri + 1, env, no_opt);
          for (int mu1 = 0; mu1 < dmu; mu1++) {
            for (int i1 = 0; i1 < di; i1++) {
              for (int j1 = 0; j1 < dj; j1++) {
                int omu = mu1 + idxmu;
                int oi = i1 + idxi;
                int oj = j1 + idxj;
                int ij = oi * N + oj;
                int ji = oj * N + oi;
                int bdx = j1 * (di * dmu) + i1 * dmu + mu1;
                eri[ij*Naux + omu] = eri[ji*Naux + omu]= buf[bdx];
              } // for j1
            } // for i1
          } // for mu1
          delete [] buf;
          idxj += dj;
        } // for j
        idxj = 0;
        idxi += di;
      } // for i
      idxi = 0;
      idxmu += dmu;
    } // for mu
  } // else
} // gen_eri_3c

void gen_eri_2c(double *eri, int Naux,
                int natm, int nbas, int nenv, int nbas_ri,
                int *atm, int *bas, double *env) {
  int Naux2 = Naux * Naux;
  int idxmu = 0;
  int idxnu = 0;
  int dmu = 0;
  int dnu = 0;
  CINTOpt *no_opt = NULL;
  if (BT::DO_CART) {
    for (int mu = 0; mu < nbas_ri; mu++) {
      dmu = CINTcgto_cart(mu+nbas,bas);
      for (int nu = 0; nu <= mu; nu++) {
        dnu = CINTcgto_cart(nu+nbas,bas);
        double *buf = new double[dmu * dnu];
        int shls[4];
        shls[0] = mu + nbas;
        shls[1] = nbas + nbas_ri;
        shls[2] = nu + nbas;
        shls[3] = nbas + nbas_ri;
        cint2e_cart(buf, shls, atm, natm, bas, nbas + nbas_ri + 1, env, no_opt);
        for (int mu1 = 0; mu1 < dmu; mu1++) {
          for (int nu1 = 0; nu1 < dnu; nu1++) {
            int omu = mu1 + idxmu;
            int onu = nu1 + idxnu;
            int bdx = nu1 * dmu + mu1;
            eri[omu * Naux + onu] = eri[onu * Naux + omu] = buf[bdx];
          } // for nu1
        } // for mu1
        delete [] buf;
        idxnu += dnu;
      } // for nu
      idxnu = 0;
      idxmu += dmu;
    } // for mu
  } // if BT::DO_CART
  else {
    for (int mu = 0; mu < nbas_ri; mu++) {
      dmu = CINTcgto_spheric(mu+nbas,bas);
      for (int nu = 0; nu <= mu; nu++) {
        dnu = CINTcgto_spheric(nu+nbas,bas);
        double *buf = new double[dmu * dnu];
        int shls[4];
        shls[0] = mu + nbas;
        shls[1] = nbas + nbas_ri;
        shls[2] = nu + nbas;
        shls[3] = nbas + nbas_ri;
        cint2e_sph(buf, shls, atm, natm, bas, nbas + nbas_ri + 1, env, no_opt);
        for (int mu1 = 0; mu1 < dmu; mu1++) {
          for (int nu1 = 0; nu1 < dnu; nu1++) {
            int omu = mu1 + idxmu;
            int onu = nu1 + idxnu;
            int bdx = nu1 * dmu + mu1;
            eri[omu * Naux + onu] = eri[onu * Naux + omu] = buf[bdx];
          } // for nu1
        } // for mu1
        delete [] buf;
        idxnu += dnu;
      } // for nu
      idxnu = 0;
      idxmu += dmu;
    } // for mu
  } // else
}

void translate_basis(int nbas, int N, int natm, int* atm, int* bas, double* env,
                     vector< vector< double> > &basis) {
  vector< vector< double > >().swap(basis);

  for (int i = 0; i < nbas; i++) {
    int ang_a = bas[ANG_OF   + BAS_SLOTS*i];
    int idx_a = bas[ATOM_OF  + BAS_SLOTS*i];
    int atm_n = atm[idx_a * ATM_SLOTS + 0];
    int idx_c = atm[idx_a * ATM_SLOTS + 1];
    int dim = binom(3+ang_a-1,ang_a);
    double x = env[idx_c + 0];
    double y = env[idx_c + 1];
    double z = env[idx_c + 2];
    int start;
    int stop;
    if (BT::DO_CART) {
      start = 0;
      stop = dim;
    }
    else {
      start = -ang_a;
      stop = ang_a + 1;
    }
    for (int j = start; j < stop; j++) {
      vector< double > basvec;
      basvec.reserve(10);
      basvec.push_back(ang_a+1);
      basvec.push_back(ang_a);
      basvec.push_back(j);
      basvec.push_back(0);
      basvec.push_back(0);
      basvec.push_back(x);
      basvec.push_back(y);
      basvec.push_back(z);
      basvec.push_back(atm_n);
      basvec.push_back(atm_n);
      basis.push_back(basvec);
    }
  }
}

void compute_g_grad_ri(int natm, int N, int Naux, int nbas, int nbas_ri,
                       double *grad, double *GF, double *Pao, double *gQmunu,
                       double *gRS, int *atm, int *bas, double *env) {
  double * grad_term = new double[3*natm](); // zeroed in each contraction
  for (int i = 0; i < 3 * natm; i++) {
    grad[i] = 0.;
  }

  printf("Printing g(S)\n");
  contract_dS(natm, N, nbas, grad_term, GF, atm, bas, env);
  for (int i = 0; i < 3 * natm; i++) {
    grad_term[i] *= -1.;
    grad[i] += grad_term[i];
  }
  for (int i = 0; i < natm; i++) {
    for (int j = 0; j < 3; j++) {
      printf("%12.8f ",grad_term[i*3+j]);
    }
    printf("\n");
  }
  printf("\n");
  printf("Printing g(T)\n");
  contract_dK(natm, N, nbas, grad_term, Pao, atm, bas, env);
  for (int i = 0; i < 3 * natm; i++) {
    grad_term[i] *= -1;
    grad[i] += grad_term[i];
  }
  for (int i = 0; i < natm; i++) {
    for (int j = 0; j < 3; j++) {
      printf("%12.8f ",grad_term[i*3+j]);
    }
    printf("\n");
  }
  printf("\n");
  printf("Printing g(Vne)\n");
  contract_dVne(natm, N, nbas, grad_term, Pao, atm, bas, env);
  for (int i = 0; i < 3 * natm; i++) {
    grad_term[i] *= 1.;
    grad[i] += grad_term[i];
  }
  for (int i = 0; i < natm; i++) {
    for (int j = 0; j < 3; j++) {
      printf("%12.8f ",grad_term[i*3+j]);
    }
    printf("\n");
  }
  printf("\n");

  printf("Printing g(PQ)\n");
  contract_d2c2e(natm, N, nbas, Naux, nbas_ri, grad_term, gRS, atm, bas, env);
  for (int i = 0; i < 3 * natm; i++) {
    grad_term[i] *= 1.;
    grad[i] += grad_term[i];
  }
  for (int i = 0; i < natm; i++) {
    for (int j = 0; j < 3; j++) {
      printf("%12.8f ",grad_term[i*3+j]);
    }
    printf("\n");
  }
  printf("\n");
  printf("Printing g(gamma3c)\n");
  contract_d3c2e(natm, N, nbas, Naux, nbas_ri, grad_term, gQmunu, atm, bas, env);
  for (int i = 0; i < 3 * natm; i++) {
    grad_term[i] *= -2.;
    grad[i] += grad_term[i];
  }
  for (int i = 0; i < natm; i++) {
    for (int j = 0; j < 3; j++) {
      printf("%12.8f ",grad_term[i*3+j]);
    }
    printf("\n");
  }
  printf("\n");

  printf("Printing grad\n");
  for (int i = 0; i < natm; i++) {
    for (int j = 0; j < 3; j++) {
      printf("%12.8f ",grad[i*3+j]);
    }
    printf("\n");
  }
  printf("\n");

  delete [] grad_term;
}

void contract_dS(int natm, int N, int nbas, double *grad_term, double *GF,
                 int *atm, int *bas, double *env) {
  int N2 = N*N;
  double *jS = new double[3*N2];
  int idx_i = 0;
  int *atm_idx = new int[N];
  int di, dj;
  int shls[2];

  for (int i = 0; i < 3*natm; i++) {
    grad_term[i] = 0.;
  }

  for (int i = 0; i < nbas; i++) {
    int idx_j = 0;
    shls[0] = i;
    di = CINTcgto_spheric(i, bas);
    for (int j = 0; j < di; j++) {
      atm_idx[idx_i + j] = bas[ATOM_OF + BAS_SLOTS*i];
    }
    for (int j = 0; j < nbas; j++) {
      dj = CINTcgto_spheric(j, bas);
      shls[1] = j;
      double *buf = new double[di*dj*3];
      cint1e_ipovlp_sph(buf,shls,atm,natm,bas,nbas,env);

      for (int x = 0; x < 3; x++) {
        for (int j1 = 0; j1 < dj; j1++) {
          for (int i1 = 0; i1 < di; i1++) {
            int oi = i1 + idx_i;
            int oj = j1 + idx_j;
            jS[x*N2 + oi * N + oj] = buf[x * dj * di + j1*di + i1];
          } // for i1
        } // for j1
      } // for x
      delete [] buf;
      idx_j += dj;
    } // for j
    idx_i += di;
  }// for i

  for (int x = 0; x < 3; x++) {
    for (int i = 0; i < N; i++) {
      int atm_i = atm_idx[i];
      for (int j = 0; j < N; j++) {
        int atm_j = atm_idx[j];
        grad_term[atm_i*3 + x] += GF[i*N+j] * jS[x*N2 + i*N+j];
        grad_term[atm_j*3 + x] += GF[i*N+j] * jS[x*N2 + j*N+i];
      } // for j
    } // for i
  } // for x

  delete [] atm_idx;
  delete [] jS;
}


void contract_dK(int natm, int N, int nbas, double *grad_term, double *Pao,
                 int *atm, int *bas, double *env) {
  int N2 = N*N;
  double *jK1 = new double[3*N2];
  double *jK2 = new double[3*N2];
  int idx_i = 0;
  int *atm_idx = new int[N];
  int di, dj;
  int shls[2];

  for (int i = 0; i < 3*natm; i++) {
    grad_term[i] = 0.;
  }

  for (int i = 0; i < nbas; i++) {
    int idx_j = 0;
    shls[0] = i;
    di = CINTcgto_spheric(i, bas);
    for (int j = 0; j < di; j++) {
      atm_idx[idx_i + j] = bas[ATOM_OF + BAS_SLOTS*i];
    }
    for (int j = 0; j < nbas; j++) {
      dj = CINTcgto_spheric(j, bas);
      shls[1] = j;
      double *buf1 = new double[di*dj*3];
      double *buf2 = new double[di*dj*3];
      cint1e_ipkin_sph(buf1,shls,atm,natm,bas,nbas,env);
      cint1e_kinip_sph(buf2,shls,atm,natm,bas,nbas,env);

      for (int x = 0; x < 3; x++) {
        for (int j1 = 0; j1 < dj; j1++) {
          for (int i1 = 0; i1 < di; i1++) {
            int oi = i1 + idx_i;
            int oj = j1 + idx_j;
            jK1[x*N2 + oi * N + oj] = buf1[x * dj * di + j1*di + i1];
            jK2[x*N2 + oi * N + oj] = buf2[x * dj * di + j1*di + i1];
          } // for i1
        } // for j1
      } // for x
      delete [] buf1;
      delete [] buf2;
      idx_j += dj;
    } // for j
    idx_i += di;
  }// for i

  for (int x = 0; x < 3; x++) {
    for (int i = 0; i < N; i++) {
      int atm_i = atm_idx[i];
      for (int j = 0; j < N; j++) {
        int atm_j = atm_idx[j];
        grad_term[atm_i*3 + x] += Pao[i*N+j] * jK1[x*N2 + i*N+j];
        grad_term[atm_j*3 + x] += Pao[i*N+j] * jK2[x*N2 + i*N+j];
      } // for j
    } // for i
  } // for x

  delete [] atm_idx;
  delete [] jK1;
  delete [] jK2;
}

void contract_dVne(int natm, int N, int nbas, double *grad_term, double *Pao,
                   int *atm, int *bas, double *env) {
  int N2 = N*N;
  double *jVdi = new double[3*N2]();
  double *jVdR = new double[3*N2]();
  int idx_i = 0;
  int *atm_idx = new int[N];
  int di, dj;
  int shls[2];

  for (int i = 0; i < 3*natm; i++) {
    grad_term[i] = 0.;
  }

  for (int i = 0; i < nbas; i++) {
    di = CINTcgto_spheric(i, bas);
    for (int j = 0; j < di; j++) {
      atm_idx[idx_i + j] = bas[ATOM_OF + BAS_SLOTS*i];
    }
    idx_i += di;
  }// for i

  for (int a = 0; a < natm; a++) {
    int atomic_num = atm[CHARGE_OF + ATM_SLOTS * a];
    // PTR_RINV_ORIG = PTR_ENV_START + a * 3;
    for (int ra = 0; ra < 3; ra++) {
      env[PTR_RINV_ORIG + ra] = env[atm[a*ATM_SLOTS + PTR_COORD] + ra];
    }
    idx_i = 0;
    for (int i = 0; i < nbas; i++) {
      int idx_j = 0;
      shls[0] = i;
      di = CINTcgto_spheric(i, bas);
      for (int j = 0; j < nbas; j++) {
        dj = CINTcgto_spheric(j, bas);
        shls[1] = j;
        double *bufr = new double[di*dj*3];
        double *bufi = new double[di*dj*3];
        cint1e_drinv_sph(bufr,shls,atm,natm,bas,nbas,env);
        cint1e_iprinv_sph(bufi,shls,atm,natm,bas,nbas,env);

        for (int x = 0; x < 3; x++) {
          for (int j1 = 0; j1 < dj; j1++) {
            for (int i1 = 0; i1 < di; i1++) {
              int oi = i1 + idx_i;
              int oj = j1 + idx_j;
              jVdR[x*N2 + oi*N + oj] = bufr[x*dj*di + j1*di +i1];
              jVdi[x*N2 + oi*N + oj] = bufi[x*dj*di + j1*di +i1];
            } // for i1
          } // for j1
        } // for x
        delete [] bufr;
        delete [] bufi;
        idx_j += dj;
      } // for j
      idx_i += di;
    } // for i

    for (int x = 0; x < 3; x++) {
      for (int i = 0; i < N; i++) {
        int atm_i = atm_idx[i];
        for (int j = 0; j < N; j++) {
          int atm_j = atm_idx[j];
          grad_term[a*3 + x] -= atomic_num * Pao[i*N+j] * jVdR[x*N2 + j*N+i];
          grad_term[atm_i*3 + x] += atomic_num * Pao[i*N+j] * jVdi[x*N2 + i*N+j];
          grad_term[atm_j*3 + x] += atomic_num * Pao[i*N+j] * jVdi[x*N2 + j*N+i];
        } // for j
      } // for i
    } // for x
  } // for a

  // printf("Printing dVdR\n");
  // for (int i = 0; i < N; i++) {
  //   for (int j = 0; j < N; j++) {
  //     printf("%8.4f ",jVdR[i*N+j]);
  //   }
  //   printf("\n");
  // }
  // printf("\n");

  delete [] atm_idx;
  delete [] jVdR;
  delete [] jVdi;
}

void contract_dVne2(int natm, int N, int nbas, double *grad_term, double *Pao,
                   int *atm, int *bas, double *env) {
  int N2 = N*N;
  double *jVne = new double[3*N2];
  int idx_i = 0;
  int *atm_idx = new int[N];
  int di, dj;
  int shls[2];

  for (int i = 0; i < 3*natm; i++) {
    grad_term[i] = 0.;
  }

  for (int i = 0; i < nbas; i++) {
    int idx_j = 0;
    shls[0] = i;
    di = CINTcgto_spheric(i, bas);
    for (int j = 0; j < di; j++) {
      atm_idx[idx_i + j] = bas[ATOM_OF + BAS_SLOTS*i];
    }
    for (int j = 0; j < nbas; j++) {
      dj = CINTcgto_spheric(j, bas);
      shls[1] = j;
      double *buf = new double[di*dj*3];
      cint1e_ipnuc_sph(buf,shls,atm,natm,bas,nbas,env);

      for (int x = 0; x < 3; x++) {
        for (int j1 = 0; j1 < dj; j1++) {
          for (int i1 = 0; i1 < di; i1++) {
            int oi = i1 + idx_i;
            int oj = j1 + idx_j;
            jVne[x*N2 + oi*N + oj] = buf[x*dj*di + j1*di + i1];
          } // for i1
        } // for j1
      } // for x
      delete [] buf;
      idx_j += dj;
    } // for j
    idx_i += di;
  }// for i

  for (int x = 0; x < 3; x++) {
    for (int i = 0; i < N; i++) {
      int atm_i = atm_idx[i];
      for (int j = 0; j < N; j++) {
        int atm_j = atm_idx[j];
        grad_term[atm_i*3 + x] += Pao[i*N+j] * jVne[x*N2 + i*N+j];
        grad_term[atm_j*3 + x] += Pao[i*N+j] * jVne[x*N2 + j*N+i];
      } // for j
    } // for i
  } // for x

  // printf("Printing dVne\n");
  // for (int i = 0; i < N; i++) {
  //   for (int j = 0; j < N; j++) {
  //     printf("%8.4f ",jVne[i*N+j]);
  //   }
  //   printf("\n");
  // }
  // printf("\n");

  delete [] atm_idx;
  delete [] jVne;
}


void contract_d2c2e(int natm, int N, int nbas, int Naux, int nbas_ri,
                    double *grad_term, double *gRS,
                    int *atm, int *bas, double *env) {
  int Naux2 = Naux * Naux;
  double *jERI2c = new double[3*Naux2];
  int *atm_idx = new int[Naux];
  int ds = 0;
  int dr = 0;
  int idx_s = 0;
  int idx_r = 0;
  int shls[4];
  for (int i = 0; i < 3*natm; i++) {
    grad_term[i] = 0.;
  }

  for (int r = 0; r < nbas_ri; r++) {
    dr = CINTcgto_spheric(r+nbas, bas);
    shls[0] = r + nbas;
    shls[1] = nbas + nbas_ri;
    // shls[1] = nbas + nbas_ri;
    for (int s = 0; s < dr; s++) {
      atm_idx[idx_r + s] = bas[ATOM_OF + BAS_SLOTS*(r+nbas)];
    } // for s
    for (int s = 0; s < nbas_ri; s++) {
      ds = CINTcgto_spheric(s+nbas, bas);
      shls[2] = s + nbas;
      shls[3] = nbas + nbas_ri;
      // shls[3] = nbas + nbas_ri;
      double *buf = new double[dr * ds * 3];
      cint2e_ip1_sph(buf,shls,atm,natm,bas,nbas+nbas_ri+1,env,NULL);

      for (int x = 0; x < 3; x++) {
        for (int r1 = 0; r1 < dr; r1++) {
          for (int s1 = 0; s1 < ds; s1++) {
            int rr = r1 + idx_r;
            int os = s1 + idx_s;
            jERI2c[x*Naux2 + rr*Naux + os] = buf[x*dr*ds + s1*dr + r1];
          } // for s1
        } // for r1
      } // for x
      delete [] buf;
      idx_s += ds;
    } // for s
    idx_s = 0;
    idx_r += dr;
  } // for r

  for (int x = 0; x < 3; x++) {
    for (int r = 0; r < Naux; r++) {
      int atm_r = atm_idx[r];
      for (int s = 0; s < Naux; s++) {
        int atm_s = atm_idx[s];
        grad_term[atm_r*3 + x] += gRS[r*Naux+s] * jERI2c[x*Naux2 + r*Naux+s];
        grad_term[atm_s*3 + x] += gRS[r*Naux+s] * jERI2c[x*Naux2 + s*Naux+r];
      } // for s
    } // for r
  } // for x

  delete [] jERI2c;
  delete [] atm_idx;
}

void contract_d3c2e(int natm, int N, int nbas, int Naux, int nbas_ri, 
                    double *grad_term, double *gQmunu,
                    int *atm, int *bas, double *env) {
  int N2 = N*N;
  int N2a = N2*Naux;
  int *atm_idx = new int[N];
  int *atm_idx_ri = new int[Naux]();
  double *jERI3c_ip1 = new double[3*N2a]; // ( di j |  P )
  double *jERI3c_ip2 = new double[3*N2a]; // (  i j | dP )
  int dmu = 0;
  int dnu = 0;
  int dq = 0;
  int idx_mu = 0;
  int idx_nu = 0;
  int idx_q = 0;

  int shls[4];

  for (int i = 0; i < 3 * natm; i++) {
    grad_term[i] = 0.;
  }

  for (int mu = 0; mu < nbas; mu++) {
    dmu = CINTcgto_spheric(mu, bas);
    shls[0] = mu;
    for (int i = 0; i < dmu; i++) {
      atm_idx[idx_mu + i] = bas[ATOM_OF + BAS_SLOTS*mu];
    } // for i
    for (int nu = 0; nu < nbas; nu++) {
      dnu = CINTcgto_spheric(nu, bas);
      shls[1] = nu;
      for (int q = 0; q < nbas_ri; q++) {
        dq = CINTcgto_spheric(q + nbas, bas);
        shls[2] = nbas + q;
        shls[3] = nbas + nbas_ri;
        // shls[3] = nbas + nbas_ri;
        int dim = dmu * dnu * dq;
        double *buf_ip1 = new double[dim*3];
        double *buf_ip2 = new double[dim*3];
        cint2e_ip1_sph(buf_ip1, shls, atm, natm, bas, nbas + nbas_ri + 1, env, NULL);
        cint2e_ip2_sph(buf_ip2, shls, atm, natm, bas, nbas + nbas_ri + 1, env, NULL);
        for (int x = 0; x < 3; x++) {
          for (int mu1 = 0; mu1 < dmu; mu1++) {
            for (int nu1 = 0; nu1 < dnu; nu1++) {
              for (int q1 = 0; q1 < dq; q1++) {
                int omu = mu1 + idx_mu;
                int onu = nu1 + idx_nu;
                int oq = q1 + idx_q;
                int mn = omu * N + onu;
                int bdx = x*dim + q1 * (dmu * dnu) + nu1 * dmu + mu1;
                jERI3c_ip1[x * N2a + mn*Naux + oq] = buf_ip1[bdx];
                jERI3c_ip2[x * N2a + mn*Naux + oq] = buf_ip2[bdx];
              } // for q
            } // for nu1
          } // for mu1
        } // for x
        delete [] buf_ip1;
        delete [] buf_ip2;
        idx_q += dq;
      } // for q
      idx_q = 0;
      idx_nu += dnu;
    } // for nu
    idx_nu = 0;
    idx_mu += dmu;
  } // for mu

  idx_q = 0;
  for (int q = 0; q < nbas_ri; q++) {
    dq = CINTcgto_spheric(q+nbas, bas);
    shls[0] = q + nbas;
    for (int s = 0; s < dq; s++) {
      atm_idx_ri[idx_q + s] = bas[ATOM_OF + BAS_SLOTS*(q+nbas)];
    } // for s
    idx_q += dq;
  }

  for (int x = 0; x < 3; x++) {
    for (int mu = 0; mu < N; mu++) {
      int atm_mu = atm_idx[mu];
      for (int nu = 0; nu < N; nu++) {
        int atm_nu = atm_idx[nu];
        for (int q = 0; q < Naux; q++) {
          int atm_q = atm_idx_ri[q];
          int mn = mu * N + nu;
          int nm = nu * N + mu;
          int muatm = atm_mu * 3 + x;
          int nuatm = atm_nu * 3 + x;
          grad_term[muatm] += gQmunu[q*N*N + mn] * jERI3c_ip1[x*N2a + mn*Naux + q];
          grad_term[nuatm] += gQmunu[q*N*N + mn] * jERI3c_ip1[x*N2a + nm*Naux + q];
        } // for q
        for (int q = 0; q < Naux; q++) {
          int atm_q = atm_idx_ri[q];
          int mn = mu * N + nu;
          int qatm = atm_q * 3 + x;
          grad_term[qatm] += gQmunu[q*N*N + mn] * jERI3c_ip2[x*N2a + mn*Naux + q];
        } // for q
      } // for nu
    } // for mu
  } // for x

  delete [] atm_idx;
  delete [] atm_idx_ri;
  delete [] jERI3c_ip1;
  delete [] jERI3c_ip2;
}
