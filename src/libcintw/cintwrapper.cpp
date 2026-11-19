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
// commented out because it is not currently used; could remove if not needed
// #else
// #include <lapacke.h>
// #include <cblas.h>
#endif

using namespace std;

bool BT::DO_CART = true;

int binom(int n, int k)
{
  return (tgamma(n+1))/(tgamma(n+1-k)*tgamma(k+1));
}

int calc_di(int i, int *bas)
{
  int idx_i = 0;
  for (int j = 0; j < i; j++) {
    idx_i += CINTcgto_cart(j, bas);
  }
  return idx_i;
}

void get_vri(double** val, int gs, int N,
                 int natm, int nbas, int nbas_ri, int nenv,
                 int* atm, int* bas, double* env)
{
  int idx_i = 0;
  int di, dj;
  int shls[4] = {0, 0, 0, gs};

  int nbasx = nbas;
  int b0 = 0;
  int Nx = N;
  if (nbas_ri>0)
  {
    b0 = nbas;
    nbasx = nbas_ri;
    Nx = nbas_ri;

    printf("  using ri basis in get_vri \n");
  }

 /*
 //CPMZ something is off about the cache size estimate
  int cache_size = int1e_grids_sph(NULL, NULL, shls, atm, natm, bas, nbas+nbas_ri, env, NULL, NULL);
  //printf(" cache_size: %4i \n",cache_size);
  cache_size += gs;
  cache_size *= 10;
  //double* cache = new double[cache_size]();
 */

 //CPMZ just let libcint allocate for itself
  double* cache = NULL;

  if (BT::DO_CART)
  {
    for (int i=0;i<nbasx;i++)
    {
      int i0 = b0 + i;
      int idx_j = 0;
      shls[0] = i0;
      di = CINTcgto_cart(i0, bas);

      for (int j=0;j<=i;j++)
      {
        int j0 = b0 + j;
        dj = CINTcgto_cart(j0, bas);
        shls[1] = j0;

        double* buf = new double[di*dj*gs]();
        int1e_grids_cart(buf, NULL, shls, atm, natm, bas, nbas, env, NULL, cache);

        for (int j1=0;j1<dj;j1++)
        for (int i1=0;i1<di;i1++)
        {
          int oi = i1 + idx_i;
          int oj = j1 + idx_j;

          int ind = oi*Nx+oj;
          //double d12 = 0.5*(Pao[ind]+Pao[oj*Nx+oi]);

          for (int k=0;k<gs;k++)
            val[ind][k] = buf[(j1*di + i1)*gs + k];
        }

        delete [] buf;
        idx_j += dj;
      } // for j
      idx_i += di;
    } // for i

  } // if BT::DO_CART
  else
  {
    for (int i=0;i<nbasx;i++)
    {
      int i0 = b0+i;
      int idx_j = 0;
      shls[0] = i0;
      di = CINTcgto_spheric(i0, bas);

      for (int j=0;j<=i;j++)
      {
        int j0 = b0+j;
        dj = CINTcgto_spheric(j0, bas);
        shls[1] = j0;

        //printf("  dij: %2i %2i  gs: %4i \n",di,dj,gs);
        double* buf = new double[di*dj*gs]();
        int1e_grids_sph(buf, NULL, shls, atm, natm, bas, nbas, env, NULL, cache);

        for (int j1=0;j1<dj;j1++)
        for (int i1=0;i1<di;i1++)
        {
          int oi = i1 + idx_i;
          int oj = j1 + idx_j;

          int ind = oi*Nx+oj;
          //double d12 = 0.5*(Pao[ind]+Pao[oj*Nx+oi]);

          for (int k=0;k<gs;k++)
            val[ind][k] = buf[(j1*di + i1)*gs + k];
        }

        delete [] buf;
        idx_j += dj;
      } // for j

      idx_i += di;
    }// for i

  } // else

  //delete [] cache;

  return;
}

void get_position_squared_integrals(double* x2, double* y2, double* z2,
                                    int N, int natm, int nbas, int nenv, int* atm, int* bas, double* env)
{
  (void)nenv;

  int N2 = N*N;
  fill(x2,x2+N2,0.);
  fill(y2,y2+N2,0.);
  fill(z2,z2+N2,0.);

  int shls[2];
  int idx_i = 0;

  for (int i=0;i<nbas;i++)
  {
    shls[0] = i;
    const int di = BT::DO_CART ? CINTcgto_cart(i, bas) : CINTcgto_spheric(i, bas);

    int idx_j = 0;
    for (int j=0;j<=i;j++)
    {
      shls[1] = j;
      const int dj = BT::DO_CART ? CINTcgto_cart(j, bas) : CINTcgto_spheric(j, bas);
      const int block_size = di * dj;

      // xx, xy, xz, yx, yy, yz, zx, zy, zz
      vector<double> buf(9*block_size,0.);
      if (BT::DO_CART)
        cint1e_rr_cart(buf.data(),shls,atm,natm,bas,nbas,env);
      else
        cint1e_rr_sph(buf.data(),shls,atm,natm,bas,nbas,env);

      for (int j1=0;j1<dj;j1++)
      {
        for (int i1=0;i1<di;i1++)
        {
          const int oi = idx_i + i1;
          const int oj = idx_j + j1;
          const int ij = j1 * di + i1;

          // Diagonal components of r_alpha r_beta:
          const double x2ij = buf[ij]; // xx
          const double y2ij = buf[4*block_size+ij]; // yy
          const double z2ij = buf[8*block_size+ij]; // zz

          x2[oi*N+oj] = x2ij;
          y2[oi*N+oj] = y2ij;
          z2[oi*N+oj] = z2ij;
          x2[oj*N+oi] = x2ij;
          y2[oj*N+oi] = y2ij;
          z2[oj*N+oi] = z2ij;
        }
      }
      idx_j += dj;
    }
    idx_i += di;
  }

}

void get_position_integrals(double* rx, double* ry, double* rz,
                            int N, int natm, int nbas, int nenv, int* atm, int* bas, double* env)
{
  int N2 = N*N;
  fill(rx,rx+N2,0.);
  fill(ry,ry+N2,0.);
  fill(rz,rz+N2,0.);

  int shls[2];
  int idx_i = 0;

  for (int i=0;i<nbas;i++)
  {
    shls[0] = i;
    const int di = BT::DO_CART ? CINTcgto_cart(i, bas) : CINTcgto_spheric(i, bas);
    int idx_j = 0;

    for (int j=0;j<=i;j++)
    {
      shls[1] = j;
      const int dj = BT::DO_CART ? CINTcgto_cart(j, bas) : CINTcgto_spheric(j, bas);

      // Three components: x, y, and z.
      vector<double> buf(3*di*dj,0.);

      if (BT::DO_CART)
        cint1e_r_cart(buf.data(),shls,atm,natm,bas,nbas,env);
      else
        cint1e_r_sph(buf.data(),shls,atm,natm,bas,nbas,env);

      const int block_size = di*dj;
      for (int j1=0;j1<dj;j1++)
      {
        for (int i1=0;i1<di;i1++)
        {
          const int oi = idx_i + i1;
          const int oj = idx_j + j1;
          const int ij = j1 * di + i1;

          const double xij = buf[ij];
          const double yij = buf[block_size+ij];
          const double zij = buf[2*block_size+ij];

          rx[oi*N+oj] = xij;
          ry[oi*N+oj] = yij;
          rz[oi*N+oj] = zij;

          rx[oj*N+oi] = xij;
          ry[oj*N+oi] = yij;
          rz[oj*N+oi] = zij;
        }
      }
      idx_j += dj;
    }
    idx_i += di;
  }
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

void get_hcore(double *hcore, double *En, double *T, int N,
               int natm, int nbas, int nenv,
               int *atm, int *bas, double *env)
{
  int idx_i = 0;
  int di, dj;
  int shls[2];

  int N2 = N*N;
  double *tmp = new double[N2];

  if (BT::DO_CART)
  {
    for (int i = 0; i < nbas; i++)
    {
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

    for (int i = 0; i < N2; i++)
      hcore[i] = tmp[i];

    idx_i = 0;
    for (int i = 0; i < nbas; i++)
    {
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
  else
  {
    for (int i = 0; i < nbas; i++)
    {
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

    if (En!=NULL)
    for (int i = 0; i < N2; i++)
      En[i] = tmp[i];
    for (int i = 0; i < N2; i++)
      hcore[i] = tmp[i];

    idx_i = 0;
    for (int i = 0; i < nbas; i++)
    {
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

  if (T!=NULL)
  for (int i = 0; i < N2; i++)
    T[i] = tmp[i];
  for (int i = 0; i < N2; i++)
    hcore[i] += tmp[i];

  delete [] tmp;
}

void get_hcore(double *hcore, int N,
               int natm, int nbas, int nenv,
               int *atm, int *bas, double *env)
{
  return get_hcore(hcore,NULL,NULL,N,natm,nbas,nenv,atm,bas,env);
}

void get_tcore(double *tcore, int N,
               int natm, int nbas, int nenv,
               int *atm, int *bas, double *env) {
  int idx_i = 0;
  int di, dj;
  int shls[2];

  int N2 = N*N;
  double *tmp = new double[N2];

  if (BT::DO_CART)
  {
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
  else
  {
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

  for (int i = 0; i < N2; i++) {
    tcore[i] = tmp[i];
  }

  delete [] tmp;
}

//ADS
void gen_pvp(double *pvp, int N,
             int natm, int nbas, int nenv,
             int *atm, int *bas, double *env) {
  int N2 = N*N;
  double *tmp = new double[N2];
  int idx_i = 0;
  int idx_j = 0;
  int di = 0;
  int dj = 0;
  int shls[2];
  CINTOpt *no_opt = NULL;
  if (BT::DO_CART) {
    idx_i = 0;
    for (int i = 0; i < nbas; i++) {
      idx_j = 0;
      di = CINTcgto_cart(i,bas);
      shls[0] = i;
      for (int j = 0; j < nbas; j++) {
        dj = CINTcgto_cart(j,bas);
        shls[1] = j;
        double *buf = new double[di * dj * 9];
        cint1e_ipnucip_cart(buf, shls, atm, natm, bas, nbas, env);
        for (int i1 = 0; i1 < di; i++) {
          for (int j1 = 0; j1 < dj; j++) {
            int oi = i1 + idx_i;
            int oj = j1 + idx_j;
            tmp[oi * N + oj] = tmp[oj * N + oi] = -(buf[0*di*dj + j1*di + i1] + buf[4*di*dj + j1*di + i1] + buf[8*di*dj + j1*di + i1]);
          } //end i1
        } //end j1
        delete [] buf;
        idx_j += dj;
      } //end j
      idx_i += di;
    } //end i
  } //end BT::DO_CART
  else {
    idx_i = 0;
    for (int i = 0; i < nbas; i++) {
      idx_j = 0;
      di = CINTcgto_spheric(i,bas);
      shls[0] = i;
      for (int j = 0; j < nbas; j++) {
        dj = CINTcgto_spheric(j,bas);
        shls[1] = j;
        #if 0
        printf("di: %i dj: %i \n",di,dj);
        #endif
        double *buf = new double[di * dj * 9];
        cint1e_ipnucip_sph(buf, shls, atm, natm, bas, nbas, env);
        #if 0
        printf("buf:\n");
        for (int g1 = 0; g1 < di*dj*9; g1++)
        {
          printf("%6.8f ",buf[g1]);
        }
        printf("\n");
        #endif
        for (int i1 = 0; i1 < di; i1++) {
          for (int j1 = 0; j1 < dj; j1++) {
            int oi = i1 + idx_i;
            int oj = j1 + idx_j;
            tmp[oi * N + oj] = tmp[oj * N + oi] = -(buf[0*di*dj + j1*di + i1]+buf[4*di*dj + j1*di + i1]+buf[8*di*dj + j1*di + i1]);
          } //end i1
        } //end j1
        delete [] buf;
        idx_j += dj;
      } //end j
      idx_i += di;
    } //end i
  } //end else

  for (int i = 0; i < N2; i++) {
    pvp[i] = tmp[i];
  }

  delete [] tmp;
}

void gen_4c_overlap_m4(double* ovlp4, size_t N,
                    int natm, int nbas, int nenv,
                    int* atm, int* bas, double* env)
{
  #ifdef COMPILE_CINTW_4C
  size_t N2 = N*N;
  size_t Npairs = N*(N+1)/2;
  CINTOpt *no_opt = NULL;

 #define USE_PARA_4C 1

 #if USE_PARA_4C
 //when j==i, parallel execution could write to same place at the same time
  int nomp = 1;
 #pragma omp parallel
  nomp = omp_get_num_threads();
  if (nbas<nomp) nomp = nbas;

  printf("  gen_4c_overlap on %i threads \n",nomp);
  if (nbas<1) { printf("\n ERROR: no basis functions \n"); exit(-1); }
 #endif


  int dmax = 1;
  for (int i=0;i<nbas;i++)
  {
    int di = 0;
    if (BT::DO_CART)
      di = CINTcgto_cart(i, bas);
    else
      di = CINTcgto_spheric(i, bas);
    if (dmax<di)
      dmax = di;
  }

 #if USE_PARA_4C
 //will need to use these offsets in the parallel loops
  size_t* shell_offset = new size_t[nbas+1];
  shell_offset[0] = 0;
  if (BT::DO_CART)
  {
    for (int i=0;i<nbas;i++)
      shell_offset[i+1] = shell_offset[i] + CINTcgto_cart(i,bas);
  }
  else
  {
    for (int i=0;i<nbas;i++)
      shell_offset[i+1] = shell_offset[i] + CINTcgto_spheric(i,bas);
  }
 #endif

  size_t dm4 = 2*dmax*dmax*dmax*dmax;
  double* bufall = new double[dm4*nomp]();

  if (BT::DO_CART)
  {
    //size_t idxi = 0, idxj = 0, idxk = 0, idxl = 0;

   #pragma omp parallel for schedule(static,1) num_threads(nomp)
    for (size_t i=0;i<nbas;i++)
    {
      size_t di, dj, dk, dl;
      int tid = omp_get_thread_num();
      double* buf = &bufall[tid*dm4];
      //double* buf = bufall;

      di = CINTcgto_cart(i, bas);
      size_t idxi = shell_offset[i];
      for (size_t j=0;j<=i;j++) //j<=i
      {
        dj = CINTcgto_cart(j, bas);
        size_t idxj = shell_offset[j];
        for (size_t k=0;k<nbas;k++)
        {
          dk = CINTcgto_cart(k, bas);
          size_t idxk = shell_offset[k];
          for (size_t l=0;l<=k;l++) //l<=k
          {
            dl = CINTcgto_cart(l, bas);
            size_t idxl = shell_offset[l];

            int shls[4] = {i, j, k, l};
            cint4c1e_cart(buf, shls, atm, natm, bas, nbas, env, no_opt);

            for (size_t l1=0;l1<dl;l1++)
            for (size_t k1=0;k1<dk;k1++)
            for (size_t j1=0;j1<dj;j1++)
            for (size_t i1=0;i1<di;i1++)
            {
              size_t oia = i1 + idxi;
              size_t oja = j1 + idxj;
              size_t oi = max(oia,oja);
              size_t oj = min(oia,oja);

              size_t oka = k1 + idxk;
              size_t ola = l1 + idxl;
              size_t ok = max(oka,ola);
              size_t ol = min(oka,ola);

              size_t ij = oi*(oi+1)/2 + oj;
              size_t kl = ok*(ok+1)/2 + ol;
              size_t bdx = l1*(dk*dj*di) + k1*(dj*di) + j1*di + i1;

              ovlp4[ij*Npairs+kl] = buf[bdx];
            }

            //idxl += dl;
          }
          //idxl = 0;
          //idxk += dk;
        }
        //idxk = 0;
        //idxj += dj;
      }
      //idxj = 0;
      //idxi += di;
    }
  }
  else
  {
    //size_t idxi = 0, idxj = 0, idxk = 0, idxl = 0;

    // Same structure but use CINTcgto_spheric + cint4c1e_sph
   #pragma omp parallel for schedule(static,1) num_threads(nomp)
    for (size_t i=0;i<nbas;i++)
    {
      size_t di, dj, dk, dl;
      int tid = omp_get_thread_num();
      double* buf = &bufall[tid*dm4];
      //double* buf = bufall;

      di = CINTcgto_spheric(i, bas);
      size_t idxi = shell_offset[i];
      for (size_t j=0;j<=i;j++) //j<=i
      {
        dj = CINTcgto_spheric(j, bas);
        size_t idxj = shell_offset[j];
        for (size_t k=0;k<nbas;k++)
        {
          dk = CINTcgto_spheric(k, bas);
          size_t idxk = shell_offset[k];
          for (size_t l=0;l<=k;l++) //l<=k
          {
            dl = CINTcgto_spheric(l, bas);
            size_t idxl = shell_offset[l];

            int shls[4] = {i, j, k, l};
            cint4c1e_sph(buf, shls, atm, natm, bas, nbas, env, no_opt);

            for (size_t l1=0;l1<dl;l1++)
            for (size_t k1=0;k1<dk;k1++)
            for (size_t j1=0;j1<dj;j1++)
            for (size_t i1=0;i1<di;i1++)
            {
              size_t oia = i1 + idxi;
              size_t oja = j1 + idxj;
              size_t oi = max(oia,oja);
              size_t oj = min(oia,oja);

              size_t oka = k1 + idxk;
              size_t ola = l1 + idxl;
              size_t ok = max(oka,ola);
              size_t ol = min(oka,ola);

              size_t ij = oi*(oi+1)/2 + oj;
              size_t kl = ok*(ok+1)/2 + ol;
              size_t bdx = l1*(dk*dj*di) + k1*(dj*di) + j1*di + i1;

              ovlp4[ij*Npairs+kl] = buf[bdx];
            }
            //idxl += dl;
          }
          //idxl = 0;
          //idxk += dk;
        }
        //idxk = 0;
        //idxj += dj;
      }
      //idxj = 0;
      //idxi += di;
    }
  }

  delete [] bufall;
 #if USE_PARA_4C
  delete [] shell_offset;
 #endif

  #else
    printf("\n\n\n ERROR: this SlaterGPU is not set up for 4c center integrals gbasis \n");
    exit(-1);
  #endif

  return;
}


void gen_4c_overlap(double* ovlp4, size_t N,
                    int natm, int nbas, int nenv,
                    int* atm, int* bas, double* env)
{
  #ifdef COMPILE_CINTW_4C
  size_t N2 = N*N;
  size_t idxi = 0, idxj = 0, idxk = 0, idxl = 0;
  size_t di, dj, dk, dl;
  CINTOpt *no_opt = NULL;

  int dmax = 1;
  for (int i=0;i<nbas;i++)
  {
    di = CINTcgto_cart(i, bas);
    if (dmax<di)
      dmax = di;
    dj = CINTcgto_spheric(i, bas);
    if (dmax<dj)
      dmax = dj;
  }
  double* buf = new double[dmax*dmax*dmax*dmax]();

  if (BT::DO_CART)
  {
    for (size_t i=0;i<nbas;i++)
    {
      di = CINTcgto_cart(i, bas);
      for (size_t j=0;j<nbas;j++)
      {
        dj = CINTcgto_cart(j, bas);
        for (size_t k=0;k<nbas;k++)
        {
          dk = CINTcgto_cart(k, bas);
          for (size_t l=0;l<nbas;l++)
          {
            dl = CINTcgto_cart(l, bas);
            int shls[4] = {i, j, k, l};

            cint4c1e_cart(buf, shls, atm, natm, bas, nbas, env, no_opt);

            for (size_t l1=0;l1<dl;l1++)
            for (size_t k1=0;k1<dk;k1++)
            for (size_t j1=0;j1<dj;j1++)
            for (size_t i1=0;i1<di;i1++)
            {
              size_t oi = i1 + idxi;
              size_t oj = j1 + idxj;
              size_t ok = k1 + idxk;
              size_t ol = l1 + idxl;
              size_t ij = oi * N + oj;
              size_t kl = ok * N + ol;
              size_t bdx = l1*(dk*dj*di) + k1*(dj*di) + j1*di + i1;
              ovlp4[ij*N2+kl] = buf[bdx];
            }

            idxl += dl;
          }
          idxl = 0;
          idxk += dk;
        }
        idxk = 0;
        idxj += dj;
      }
      idxj = 0;
      idxi += di;
    }
  }
  else
  {
    // Same structure but use CINTcgto_spheric + cint4c1e_sph
    for (size_t i=0;i<nbas;i++)
    {
      di = CINTcgto_spheric(i, bas);
      for (size_t j=0;j<nbas;j++)
      {
        dj = CINTcgto_spheric(j, bas);
        for (size_t k=0;k<nbas;k++)
        {
          dk = CINTcgto_spheric(k, bas);
          for (size_t l=0;l<nbas;l++)
          {
            dl = CINTcgto_spheric(l, bas);
            int shls[4] = {i, j, k, l};

            cint4c1e_sph(buf, shls, atm, natm, bas, nbas, env, no_opt);

            for (size_t l1=0;l1<dl;l1++)
            for (size_t k1=0;k1<dk;k1++)
            for (size_t j1=0;j1<dj;j1++)
            for (size_t i1=0;i1<di;i1++)
            {
              size_t oi = i1 + idxi;
              size_t oj = j1 + idxj;
              size_t ok = k1 + idxk;
              size_t ol = l1 + idxl;
              size_t ij = oi * N + oj;
              size_t kl = ok * N + ol;
              size_t bdx = l1*(dk*dj*di) + k1*(dj*di) + j1*di + i1;
              ovlp4[ij*N2+kl] = buf[bdx];
            }
            idxl += dl;
          }
          idxl = 0;
          idxk += dk;
        }
        idxk = 0;
        idxj += dj;
      }
      idxj = 0;
      idxi += di;
    }
  }

  delete [] buf;

  #else
    printf("\n\n\n ERROR: this SlaterGPU is not set up for 4c center integrals gbasis \n");
    exit(-1);
  #endif
  return;
}

void gen_eri_4b(double** eri, int N1, int N2, int N3, int N,
             int natm, int nbas1, int nbas2, int nbas3, int nbas, int nenv,
             int* atm, int* bas, double* env, int nomp)
{
  int N4 = N-N1-N2-N3;
  int nbas12 = nbas1+nbas2;
  int nbas13 = nbas12+nbas3;
  if (nomp<1) nomp = 1;
  if (nomp>nbas3) nomp = nbas3;

  int idxi = 0;
  int idxj = 0;
  //int idxk = 0;
  CINTOpt* no_opt = NULL;

 //buffer size
  int mdi = 0; int mdj = 0; int mdk = 0; int mdl = 0;
  for (int i=0;i<nbas1;i++)
  {
    int di = CINTcgto_cart(i,bas);
    if (di>mdi) mdi = di;
  }
  for (int j=nbas1;j<nbas12;j++)
  {
    int dj = CINTcgto_cart(j,bas);
    if (dj>mdj) mdj = dj;
  }
  for (int k=nbas12;k<nbas13;k++)
  {
    int dk = CINTcgto_cart(k,bas);
    if (dk>mdk) mdk = dk;
  }
  for (int l=nbas13;l<nbas;l++)
  {
    int dl = CINTcgto_cart(l,bas);
    if (dl>mdl) mdl = dl;
  }

  int sizeb = mdi*mdj*mdk*mdl;
  double* bufall = new double[sizeb*nomp];
  if (BT::DO_CART)
  {
    for (int i=0;i<nbas1;i++)
    {
      int di = CINTcgto_cart(i,bas);
      for (int j=nbas1;j<nbas12;j++)
      {
        int dj = CINTcgto_cart(j,bas);
       #pragma omp parallel for num_threads(nomp)
        for (int k=nbas12;k<nbas13;k++)
        {
          //int tid = 0;
          int tid = omp_get_thread_num();
          double* buf = &bufall[sizeb*tid];

          int idxk = 0;
          for (int k1=nbas12;k1<k;k1++) idxk += CINTcgto_cart(k1,bas);
          int idxl = 0;

          int dk = CINTcgto_cart(k,bas);
          for (int l=nbas13;l<nbas;l++)
          {
            int dl = CINTcgto_cart(l,bas);
            int shls[4];
            shls[0] = i;
            shls[1] = j;
            shls[2] = k;
            shls[3] = l;
            cint2e_cart(buf, shls, atm, natm, bas, nbas, env, no_opt);

            for (int l1=0;l1<dl;l1++)
            {
              int ol = l1 + idxl;
              for (int k1=0;k1<dk;k1++)
              {
                int ok = k1 + idxk;
                for (int j1=0;j1<dj;j1++)
                {
                  int oj = j1 + idxj;
                  for (int i1=0;i1<di;i1++)
                  {
                    int oi = i1 + idxi;
                    int ij = oi * N2 + oj;
                    int kl = ok * N4 + ol;
                    int bdx = l1*(dk*dj*di) + k1*(dj*di) + j1*(di) + i1;
                    // printf("oi %2d oj %2d ok %2d ol %2d ij %6d kl %6d\n",oi, oj, ok, ol,ij,kl);
                    eri[ij][kl] = buf[bdx];
                  } // for i1
                } // for j1
              } // for k1
            } // for l1
            idxl += dl;
          } // for l
          idxl = 0;
          //idxk += dk;
        } // for k
        //idxk = 0;
        idxj += dj;
      } // for j
      idxj = 0;
      idxi += di;
    } // for i
  } // if BT::DO_CART
  else
  {
    for (int i=0;i<nbas1;i++)
    {
      int di = CINTcgto_spheric(i,bas);
      for (int j=nbas1;j<nbas12;j++)
      {
        int dj = CINTcgto_spheric(j,bas);
       #pragma omp parallel for num_threads(nomp)
        for (int k=nbas12;k<nbas13;k++)
        {
          int tid = omp_get_thread_num();
          double* buf = &bufall[sizeb*tid];

          int idxk = 0;
          for (int k1=nbas12;k1<k;k1++) idxk += CINTcgto_spheric(k1,bas);

          int idxl = 0;
          int dk = CINTcgto_spheric(k,bas);
          for (int l=nbas13;l<nbas;l++)
          {
            int dl = CINTcgto_spheric(l,bas);
            int shls[4];
            shls[0] = i;
            shls[1] = j;
            shls[2] = k;
            shls[3] = l;
            cint2e_sph(buf, shls, atm, natm, bas, nbas, env, no_opt);

            for (int l1=0;l1<dl;l1++)
            {
              int ol = l1 + idxl;
              for (int k1=0;k1<dk;k1++)
              {
                int ok = k1 + idxk;
                for (int j1=0;j1<dj;j1++)
                {
                  int oj = j1 + idxj;
                  for (int i1=0;i1<di;i1++)
                  {
                    int oi = i1 + idxi;
                    int ij = oi * N2 + oj;
                    int kl = ok * N4 + ol;
                    int bdx = l1*(dk*dj*di) + k1*(dj*di) + j1*(di) + i1;
                    // printf("oi %2d oj %2d ok %2d ol %2d ij %6d kl %6d\n",oi, oj, ok, ol,ij,kl);
                    eri[ij][kl] = buf[bdx];
                  } // for i1
                } // for j1
              } // for k1
            } // for l1
            idxl += dl;
          } // for l
          idxl = 0;
          //idxk += dk;
        } // for k
        //idxk = 0;
        idxj += dj;
      } // for j
      idxj = 0;
      idxi += di;
    } // for i
  } // else

  delete [] bufall;

  return;
}

double con_eri_4b(double* Pao1, double* Pao2, int N1, int N2, int N3, int N,
             int natm, int nbas1, int nbas2, int nbas3, int nbas, int nenv,
             int* atm, int* bas, double* env, int nomp)
{
  int N4 = N-N1-N2-N3;
  int nbas12 = nbas1+nbas2;
  int nbas13 = nbas12+nbas3;
  int nbas4  = nbas-nbas13;
  if (nomp<1) nomp = 1;
  CINTOpt* no_opt = NULL;

  if (BT::DO_CART)
  {
    printf("\n ERROR: con_eri_4b does not support CART functions \n");
    exit(-1);
  }

  // dims/offsets per block, computed once
  auto build = [&](int base, int n, vector<int>& nf, vector<int>& off)
  {
    nf.resize(n); off.resize(n);
    int cum=0;
    for (int a=0;a<n;a++){ nf[a] = CINTcgto_spheric(base+a,bas); off[a] = cum; cum += nf[a]; }
  };
  vector<int> nfi,offi,nfj,offj,nfk,offk,nfl,offl;
  build(0,       nbas1, nfi, offi);
  build(nbas1,   nbas2, nfj, offj);
  build(nbas12,  nbas3, nfk, offk);
  build(nbas13,  nbas4, nfl, offl);

  int mdi = *max_element(nfi.begin(),nfi.end());
  int mdj = *max_element(nfj.begin(),nfj.end());
  int mdk = *max_element(nfk.begin(),nfk.end());
  int mdl = *max_element(nfl.begin(),nfl.end());
  int sizeb = mdi*mdj*mdk*mdl;
  double* bufall = new double[sizeb*nomp];

  double vtot = 0.;

 #pragma omp parallel for collapse(4) num_threads(nomp) schedule(dynamic) reduction(+:vtot)
  for (int ii=0;ii<nbas1;ii++)
  for (int jj=0;jj<nbas2;jj++)
  for (int kk=0;kk<nbas3;kk++)
  for (int ll=0;ll<nbas4;ll++)
  {
    int tid = omp_get_thread_num();
    double* buf = &bufall[sizeb*tid];

    int i=ii, j=nbas1+jj, k=nbas12+kk, l=nbas13+ll;
    int di=nfi[ii], dj=nfj[jj], dk=nfk[kk], dl=nfl[ll];
    int idxi=offi[ii], idxj=offj[jj], idxk=offk[kk], idxl=offl[ll];

    int shls[4] = {i,j,k,l};
    cint2e_sph(buf, shls, atm, natm, bas, nbas, env, no_opt);

    double v1 = 0.;
    for (int l1=0;l1<dl;l1++)
    {
      int ol = l1+idxl;
      for (int k1=0;k1<dk;k1++)
      {
        int ok = k1+idxk;
        double pkl = Pao2[ok*N4+ol];
        for (int j1=0;j1<dj;j1++)
        {
          int oj = j1+idxj;
          for (int i1=0;i1<di;i1++)
          {
            int oi = i1+idxi;
            int bdx = l1*(dk*dj*di) + k1*(dj*di) + j1*di + i1;
            v1 += buf[bdx]*Pao1[oi*N2+oj]*pkl;
          }
        }
      }
    }
    vtot += v1;
  }

  delete [] bufall;
  return vtot;
}

void gen_eri_2b(double** eri, int N1, int N,
             int natm, int nbas1, int nbas, int nenv,
             int* atm, int* bas, double* env, int nomp)
{
  int N2 = N-N1;
  if (nomp<1) nomp = 1;

  int idxi = 0;
  int idxj = 0;
  int idxk = 0;
  CINTOpt* no_opt = NULL;

 //buffer size
  int mdi = 0; int mdk = 0;
  for (int i=0;i<nbas1;i++)
  {
    int di = CINTcgto_cart(i,bas);
    if (di>mdi) mdi = di;
  }
  for (int k=nbas1;k<nbas;k++)
  {
    int dk = CINTcgto_cart(k,bas);
    if (dk>mdk) mdk = dk;
  }

  int sizeb = mdi*mdi*mdk*mdk;
  double* bufall = new double[sizeb*nomp];
  if (BT::DO_CART)
  {
    for (int i=0;i<nbas1;i++)
    {
      int di = CINTcgto_cart(i,bas);
      for (int j=0;j<nbas1;j++)
      {
        int dj = CINTcgto_cart(j,bas);
       #pragma omp parallel for num_threads(nomp)
        for (int k=nbas1;k<nbas;k++)
        {
          //int tid = 0;
          int tid = omp_get_thread_num();
          double* buf = &bufall[sizeb*tid];

          int idxk = 0;
          for (int k1=nbas1;k1<k;k1++) idxk += CINTcgto_cart(k1,bas);

          int idxl = 0;
          int dk = CINTcgto_cart(k,bas);
          for (int l=nbas1;l<nbas;l++)
          {
            int dl = CINTcgto_cart(l,bas);
            int shls[4];
            shls[0] = i;
            shls[1] = j;
            shls[2] = k;
            shls[3] = l;
            cint2e_cart(buf, shls, atm, natm, bas, nbas, env, no_opt);

            for (int l1=0;l1<dl;l1++)
            {
              int ol = l1 + idxl;
              for (int k1=0;k1<dk;k1++)
              {
                int ok = k1 + idxk;
                for (int j1=0;j1<dj;j1++)
                {
                  int oj = j1 + idxj;
                  for (int i1=0;i1<di;i1++)
                  {
                    int oi = i1 + idxi;
                    int ij = oi * N1 + oj;
                    int kl = ok * N2 + ol;
                    int bdx = l1*(dk*dj*di) + k1*(dj*di) + j1*(di) + i1;
                    // printf("oi %2d oj %2d ok %2d ol %2d ij %6d kl %6d\n",oi, oj, ok, ol,ij,kl);
                    eri[ij][kl] = buf[bdx];
                  } // for i1
                } // for j1
              } // for k1
            } // for l1
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
  else
  {
    for (int i=0;i<nbas1;i++)
    {
      int di = CINTcgto_spheric(i,bas);
      for (int j=0;j<nbas1;j++)
      {
        int dj = CINTcgto_spheric(j,bas);
       #pragma omp parallel for num_threads(nomp)
        for (int k=nbas1;k<nbas;k++)
        {
          //int tid = 0;
          int tid = omp_get_thread_num();
          double* buf = &bufall[sizeb*tid];

          int idxk = 0;
          for (int k1=nbas1;k1<k;k1++) idxk += CINTcgto_spheric(k1,bas);

          int idxl = 0;
          int dk = CINTcgto_spheric(k,bas);
          for (int l=nbas1;l<nbas;l++)
          {
            int dl = CINTcgto_spheric(l,bas);
            int shls[4];
            shls[0] = i;
            shls[1] = j;
            shls[2] = k;
            shls[3] = l;
            cint2e_sph(buf, shls, atm, natm, bas, nbas, env, no_opt);

            for (int l1=0;l1<dl;l1++)
            {
              int ol = l1 + idxl;
              for (int k1=0;k1<dk;k1++)
              {
                int ok = k1 + idxk;
                for (int j1=0;j1<dj;j1++)
                {
                  int oj = j1 + idxj;
                  for (int i1=0;i1<di;i1++)
                  {
                    int oi = i1 + idxi;
                    int ij = oi * N1 + oj;
                    int kl = ok * N2 + ol;
                    int bdx = l1*(dk*dj*di) + k1*(dj*di) + j1*(di) + i1;
                    // printf("oi %2d oj %2d ok %2d ol %2d ij %6d kl %6d\n",oi, oj, ok, ol,ij,kl);
                    eri[ij][kl] = buf[bdx];
                  } // for i1
                } // for j1
              } // for k1
            } // for l1
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

  delete [] bufall;

  return;
}

double con_eri_2b(double* Pao1, double* Pao2, int N1, int N,
             int natm, int nbas1, int nbas, int nenv,
             int* atm, int* bas, double* env, int nomp)
{
  int N2 = N-N1;
  if (nomp<1) nomp = 1;
  CINTOpt* no_opt = NULL;

  if (BT::DO_CART)
  {
    printf("\n ERROR: con_eri_2b does not support CART functions \n");
    exit(-1);
  }

  // shell dims / offsets, precomputed once (was recomputed O(nbas) per k before)
  vector<int> nf1(nbas1), offs1(nbas1);
  { int cum=0; for (int i=0;i<nbas1;i++) { nf1[i] = CINTcgto_spheric(i,bas); offs1[i] = cum; cum += nf1[i]; } }

  int nbas2 = nbas-nbas1;
  vector<int> nf2(nbas2), offs2(nbas2);
  { int cum=0; for (int kk=0;kk<nbas2;kk++) { nf2[kk] = CINTcgto_spheric(nbas1+kk,bas); offs2[kk] = cum; cum += nf2[kk]; } }

  int mdi=0; for (int i=0;i<nbas1;i++) if (nf1[i]>mdi) mdi = nf1[i];
  int mdk=0; for (int kk=0;kk<nbas2;kk++) if (nf2[kk]>mdk) mdk = nf2[kk];
  int sizeb = mdi*mdi*mdk*mdk;
  double* bufall = new double[sizeb*nomp];

  // unique shell pairs per block, with degeneracy factor (1 diagonal, 2 off-diagonal)
  struct SP { int a,b; double fac; };
  vector<SP> pairs1; pairs1.reserve(nbas1*(nbas1+1)/2);
  for (int i=0;i<nbas1;i++) for (int j=0;j<=i;j++) pairs1.push_back({i,j,(i==j)?1.0:2.0});

  vector<SP> pairs2; pairs2.reserve(nbas2*(nbas2+1)/2);
  for (int k=0;k<nbas2;k++) for (int l=0;l<=k;l++) pairs2.push_back({k,l,(k==l)?1.0:2.0});

  int np1 = (int)pairs1.size(), np2 = (int)pairs2.size();
  double vtot = 0.;

 #pragma omp parallel for num_threads(nomp) collapse(2) schedule(dynamic) reduction(+:vtot)
  for (int p=0;p<np1;p++)
  for (int q=0;q<np2;q++)
  {
    int tid = omp_get_thread_num();
    double* buf = &bufall[sizeb*tid];

    int i = pairs1[p].a, j = pairs1[p].b; double fij = pairs1[p].fac;
    int di = nf1[i], dj = nf1[j], idxi = offs1[i], idxj = offs1[j];

    int kk = pairs2[q].a, ll = pairs2[q].b; double fkl = pairs2[q].fac;
    int dk = nf2[kk], dl = nf2[ll], idxk = offs2[kk], idxl = offs2[ll];
    int k = nbas1+kk, l = nbas1+ll;

    int shls[4] = {i,j,k,l};

    //could screen wrt Pao1,Pao2 here (Schwarz-type bound using fij*fkl)
    cint2e_sph(buf, shls, atm, natm, bas, nbas, env, no_opt);

    double v1 = 0.;
    for (int l1=0;l1<dl;l1++)
    {
      int ol = l1+idxl;
      for (int k1=0;k1<dk;k1++)
      {
        int ok = k1+idxk;
        double pkl = Pao2[ok*N2+ol];
        for (int j1=0;j1<dj;j1++)
        {
          int oj = j1+idxj;
          for (int i1=0;i1<di;i1++)
          {
            int oi = i1+idxi;
            int bdx = l1*(dk*dj*di) + k1*(dj*di) + j1*di + i1;
            v1 += buf[bdx]*Pao1[oi*N1+oj]*pkl;
          }
        }
      }
    }
    vtot += fij*fkl*v1;
  }

  delete [] bufall;
  return vtot;
}

void gen_ven_matrix_direct(double* V, bool clearV, double Z1, double* Pao,
                            int Nn, int Ne, int Nt,
                            int natoms, int nbas1, int nbas, int nenv,
                            int* atm, int* bas, double* env, int nomp)
{
  int Nn2 = Nn*Nn;
  if (clearV)
  for (int p=0;p<Nn2;p++)
    V[p] = 0.;

  CINTOpt* no_opt = NULL;
  if (nomp<1) nomp = 1;

  if (BT::DO_CART)
  { printf("\n ERROR: gen_ven_matrix_direct does not support Cartesian GTOs \n"); exit(-1); }

  int mdi = 0;
  for (int i=0;i<nbas;i++)
  {
    int di = CINTcgto_cart(i,bas);
    if (di>mdi) mdi = di;
  }
  int sbuf = mdi*mdi*mdi*mdi;

  vector<int> nf(nbas);
  vector<int> offs(nbas);
  int cum = 0;
  for (int i=0;i<nbas;i++)
  {
    nf[i] = CINTcgto_spheric(i, bas);
    offs[i] = cum; cum += nf[i];
  }
  if (cum!=Nt || Nt!=Nn+Ne)
  {
    printf("\n ERROR: gen_ven_matrix_direct AO count mismatch: shells sum to %d, Nt=%d, Nn+Ne=%d \n",cum,Nt,Nn+Ne);
    exit(-1);
  }
  bool aligned = (nbas1<nbas ? offs[nbas1]==Nn : cum==Nn);
  if (!aligned)
  {
    printf("\n ERROR: gen_ven_matrix_direct nuclear/electron block boundary doesn't align with shells \n");
    exit(-1);
  }

  //unique nuclear-nuclear bra shell pairs (i>=j), nuclear-block shells only
  struct SP { int i, j; };
  vector<SP> bra_pairs;
  bra_pairs.reserve(nbas1*(nbas1+1)/2);
  for (int i=0;i<nbas1;i++)
  for (int j=0;j<=i;j++)
    bra_pairs.push_back({i,j});
  int nbra = (int)bra_pairs.size();

  if (nomp>nbra && nbra>0) nomp = nbra;
  if (nomp<1) nomp = 1;

  double* bufall = new double[(size_t)sbuf*nomp];

  int shls[4];
  #pragma omp parallel for num_threads(nomp) schedule(dynamic) private(shls)
  for (int bp=0;bp<nbra;bp++)
  {
    int tid = omp_get_thread_num();
    double* buf = &bufall[(size_t)tid*sbuf];

    int i = bra_pairs[bp].i, j = bra_pairs[bp].j;
    int di = nf[i], dj = nf[j], idxi = offs[i], idxj = offs[j];
    bool fold = (i!=j);   //see gen_ven_direct: only stands in for the
                           //omitted (j,i) shell pair; a self-paired shell
                           //already visits both orderings as separate
                           //(i1,j1) iterations below

    for (int k=nbas1;k<nbas;k++)   //electron shells
    for (int l=nbas1;l<nbas;l++)   //electron shells, FULL sweep (Pao not
                                    //assumed symmetric, and no ket symmetry
                                    //is available here regardless)
    {
      int dk = nf[k], dl = nf[l], idxk = offs[k], idxl = offs[l];

      shls[0] = i; shls[1] = j; shls[2] = k; shls[3] = l;
      cint2e_sph(buf,shls,atm,natoms,bas,nbas,env,no_opt);

      for (int l1=0;l1<dl;l1++)
      for (int k1=0;k1<dk;k1++)
      {
        int r = (k1+idxk)-Nn;   //electron function index
        int s = (l1+idxl)-Nn;   //electron function index
        double pd = Pao[r*Ne+s];

        for (int j1=0;j1<dj;j1++)
        for (int i1=0;i1<di;i1++)
        {
          int oi = i1+idxi, oj = j1+idxj;
          int bdx = l1*(dk*dj*di) + k1*(dj*di) + j1*di + i1;
          double contrib = -Z1*pd*buf[bdx];

          V[oi*Nn+oj] += contrib;
          if (fold)
            V[oj*Nn+oi] += contrib;   //stands in for the (j,i) shell pair
        }
      }
    }
  }

  delete [] bufall;
  return;
}

void gen_eri_2(double** eri, int N,
             int natm, int nbas, int nenv,
             int* atm, int* bas, double* env, int nomp)
{
  CINTOpt* no_opt = NULL;
  if (nomp<1) nomp = 1;

  if (BT::DO_CART)
  { printf("\n ERROR: gen_eri_2 does not support Cartesian GTOs \n"); exit(-1); }

  int mdi = 0;
  for (int i=0;i<nbas;i++)
  {
    int di = CINTcgto_cart(i,bas);
    if (di>mdi) mdi = di;
  }

  int sbuf = mdi*mdi*mdi*mdi;
  double* bufall = new double[sbuf*nomp];

  vector<int> nf(nbas);
  vector<int> offs(nbas);
  int cum = 0;
  for (int i=0;i<nbas;i++)
  {
    nf[i] = CINTcgto_spheric(i, bas);
    offs[i] = cum; cum += nf[i];
  }

  //unique shell pairs
  struct SP { int i, j; };
  vector<SP> pairs;
  pairs.reserve(nbas*(nbas+1)/2);
  for (int i=0;i<nbas;i++)
  for (int j=0;j<=i;j++)
    pairs.push_back({i, j});
  int npair = (int)pairs.size();

  if (nomp>npair) nomp = npair;

  int shls[4];
  #pragma omp parallel for num_threads(nomp) schedule(dynamic) private(shls)
  for (int p=0;p<npair;p++)
  {
    int tid = omp_get_thread_num();
    double* buf = &bufall[tid*sbuf];

    int i = pairs[p].i, j = pairs[p].j;
    int di = nf[i], dj = nf[j], idxi = offs[i], idxj = offs[j];

    for (int q=0;q<=p;q++) // ket pair index <= bra pair index -> bra-ket swap folded in
    {
      int k = pairs[q].i; int l = pairs[q].j;
      int dk = nf[k]; int dl = nf[l]; int idxk = offs[k]; int idxl = offs[l];

      shls[0] = i; shls[1] = j; shls[2] = k; shls[3] = l;
      cint2e_sph(buf,shls,atm,natm,bas,nbas,env,no_opt);

      for (int l1=0;l1<dl;l1++)
      for (int k1=0;k1<dk;k1++)
      for (int j1=0;j1<dj;j1++)
      for (int i1=0;i1<di;i1++)
      {
        int oi = i1+idxi; int oj = j1+idxj; int ok = k1+idxk; int ol = l1+idxl;
        int bdx = l1*(dk*dj*di) + k1*(dj*di) + j1*di + i1;
        double val = buf[bdx];

        int ij = oi*N+oj; int ji = oj*N+oi;
        int kl = ok*N+ol; int lk = ol*N+ok;
        eri[ij][kl]=val; eri[ji][kl]=val; eri[ij][lk]=val; eri[ji][lk]=val;
        eri[kl][ij]=val; eri[kl][ji]=val; eri[lk][ij]=val; eri[lk][ji]=val;
      }
    }
  }

  delete [] bufall;
  return;
}

void gen_eri(double** eri, int N,
             int natm, int nbas, int nenv,
             int* atm, int* bas, double* env, int nomp)
{
  CINTOpt* no_opt = NULL;
  if (nomp<1) nomp = 1;

  int mdi = 0;
  for (int i=0;i<nbas;i++)
  {
    int di = CINTcgto_cart(i,bas);
    if (di>mdi) mdi = di;
  }

  int sbuf = mdi*mdi*mdi*mdi;
  double* bufall = new double[sbuf*nomp];

  if (BT::DO_CART)
  {
    int shls[4];
    int idxi = 0;
    for (int i=0;i<nbas;i++)
    {
      int di = CINTcgto_cart(i,bas);

     #pragma omp parallel for num_threads(nomp) private(shls)
      for (int j=0;j<nbas;j++)
      {
        int tid = omp_get_thread_num();
        double* buf = &bufall[tid*sbuf];

        int idxj = 0; for (int j1=0;j1<j;j1++) idxj += CINTcgto_cart(j1,bas);

        int dj = CINTcgto_cart(j,bas);
        int idxk = 0;
        for (int k=0;k<nbas;k++)
        {
          int dk = CINTcgto_cart(k,bas);
          int idxl = 0;
          for (int l=0;l<nbas;l++)
          {
            int dl = CINTcgto_cart(l,bas);
            shls[0] = i;
            shls[1] = j;
            shls[2] = k;
            shls[3] = l;
            cint2e_cart(buf, shls, atm, natm, bas, nbas, env, no_opt);

            for (int l1=0;l1<dl;l1++)
            for (int k1=0;k1<dk;k1++)
            for (int j1=0;j1<dj;j1++)
            for (int i1=0;i1<di;i1++)
            {
              int oi = i1 + idxi;
              int oj = j1 + idxj;
              int ok = k1 + idxk;
              int ol = l1 + idxl;
              int ij = oi * N + oj;
              int kl = ok * N + ol;
              int bdx = l1*(dk*dj*di) + k1*(dj*di) + j1*(di) + i1;
              // printf("oi %2d oj %2d ok %2d ol %2d ij %6d kl %6d\n",oi, oj, ok, ol,ij,kl);
              eri[ij][kl] = buf[bdx];
              //eri[ij][kl] = eri[kl][ij] = buf[bdx];

            } // for l1k1j1l1
            idxl += dl;
          } // for l
          idxl = 0;
          idxk += dk;
        } // for k
        idxk = 0;
        //idxj += dj;
      } // for j
      //idxj = 0;
      idxi += di;
    } // for i
  } // if BT::DO_CART
  else
  {
    int shls[4];
    int idxi = 0;
    for (int i=0;i<nbas;i++)
    {
      int di = CINTcgto_spheric(i,bas);
     #pragma omp parallel for num_threads(nomp) private(shls)
      for (int j=0;j<nbas;j++)
      {
        int tid = omp_get_thread_num();
        double* buf = &bufall[tid*sbuf];

        int idxj = 0; for (int j1=0;j1<j;j1++) idxj += CINTcgto_spheric(j1,bas);

        int dj = CINTcgto_spheric(j,bas);

        int idxk = 0;
        for (int k=0;k<nbas;k++)
        {
          int dk = CINTcgto_spheric(k,bas);

          int idxl = 0;
          for (int l=0;l<nbas;l++)
          {
            int dl = CINTcgto_spheric(l,bas);
            shls[0] = i;
            shls[1] = j;
            shls[2] = k;
            shls[3] = l;
            cint2e_sph(buf, shls, atm, natm, bas, nbas, env, no_opt);

            for (int l1=0;l1<dl;l1++)
            for (int k1=0;k1<dk;k1++)
            for (int j1=0;j1<dj;j1++)
            for (int i1=0;i1<di;i1++)
            {
              int oi = i1 + idxi;
              int oj = j1 + idxj;
              int ok = k1 + idxk;
              int ol = l1 + idxl;
              int ij = oi * N + oj;
              int kl = ok * N + ol;
              int bdx = l1*(dk*dj*di) + k1*(dj*di) + j1*(di) + i1;
              // printf("oi %2d oj %2d ok %2d ol %2d ij %6d kl %6d\n",oi, oj, ok, ol,ij,kl);
              eri[ij][kl] = buf[bdx];
              //eri[ij][kl] = eri[kl][ij] = buf[bdx];

            } // for l1k1j1i1
            idxl += dl;
          } // for l
          idxl = 0;
          idxk += dk;
        } // for k
        idxk = 0;
        //idxj += dj;
      } // for j
      //idxj = 0;
      idxi += di;
    } // for i
  } // else

  delete [] bufall;

  return;
}

void gen_eri(double** eri, int N,
             int natm, int nbas, int nenv,
             int* atm, int* bas, double* env)
{
  return gen_eri(eri,N,natm,nbas,nenv,atm,bas,env,1);
}

void gen_jMOI_gto(double** eri, int N,
                  int natm, int nbas, int nenv,
                  int *atm, int *bas, double *env)
{
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
  if (BT::DO_CART)
  {
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
                    int i = oi;
                    int jkl = oj * N * N + ok * N + ol;
                    int k = ok;
                    int lij = ol * N * N + oi * N + oj;
                    int bdx = l1*(dk*dj*di) + k1*(dj*di) + j1*(di) + i1;
                    // printf("oi %2d oj %2d ok %2d ol %2d ij %6d kl %6d\n",oi, oj, ok, ol,ij,kl);
                    eri[i][jkl] = eri[k][lij] = buf[bdx];
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
  else
  {
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
                    int i = oi;
                    int jkl = oj * N * N + ok * N + ol;
                    int k = ok;
                    int lij = ol * N * N + oi * N + oj;
                    int bdx = l1*(dk*dj*di) + k1*(dj*di) + j1*(di) + i1;
                    // printf("oi %2d oj %2d ok %2d ol %2d ij %6d kl %6d\n",oi, oj, ok, ol,ij,kl);
                    eri[i][jkl] = eri[k][lij] = buf[bdx];
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
  //int N2 = N*N;
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
  //int N2 = N*N;
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
                eri[(size_t)ij*Naux + omu] = eri[(size_t)ji*Naux + omu]= buf[bdx];
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
                eri[(size_t)ij*Naux + omu] = eri[(size_t)ji*Naux + omu]= buf[bdx];
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
                int *atm, int *bas, double *env)
{
  int Naux2 = Naux * Naux;
  int idxmu = 0;
  int idxnu = 0;
  int dmu = 0;
  int dnu = 0;
  CINTOpt *no_opt = NULL;
  if (BT::DO_CART)
  {
    for (int mu = 0; mu < nbas_ri; mu++)
    {
      dmu = CINTcgto_cart(mu+nbas,bas);
      for (int nu = 0; nu <= mu; nu++)
      {
        dnu = CINTcgto_cart(nu+nbas,bas);
        double *buf = new double[dmu * dnu];
        int shls[4];
        shls[0] = mu + nbas;
        shls[1] = nbas + nbas_ri;
        shls[2] = nu + nbas;
        shls[3] = nbas + nbas_ri;
        cint2e_cart(buf, shls, atm, natm, bas, nbas + nbas_ri + 1, env, no_opt);
        for (int mu1 = 0; mu1 < dmu; mu1++)
        {
          for (int nu1 = 0; nu1 < dnu; nu1++)
          {
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
                       double* grad, double* GF, double* Pao, double* gQmunu,
                       double* gRS, int* atm, int* bas, double* env)
{
  int prl = 2;

  double* grad_term = new double[3*natm](); // zeroed in each contraction
  for (int i=0;i<3*natm;i++)
    grad[i] = 0.;

  if (prl>0)
    printf("Printing g(S)\n");
  contract_dS(natm, N, nbas, grad_term, GF, atm, bas, env);
  for (int i=0;i<3*natm;i++)
  {
    grad_term[i] *= -1.;
    grad[i] += grad_term[i];
  }
  if (prl>0)
  for (int i=0;i<natm;i++)
  {
    for (int j=0;j<3;j++)
      printf("%12.8f ",grad_term[i*3+j]);
    printf("\n");
  }
  if (prl>0)
    printf("\n");
  if (prl>0)
    printf("Printing g(T)\n");
  contract_dK(natm, N, nbas, grad_term, Pao, atm, bas, env);
  for (int i = 0; i < 3 * natm; i++)
  {
    grad_term[i] *= -1;
    grad[i] += grad_term[i];
  }
  if (prl>0)
  for (int i=0;i<natm;i++)
  {
    for (int j=0;j<3;j++)
      printf("%12.8f ",grad_term[i*3+j]);
    printf("\n");
  }
  if (prl>0)
    printf("\n");
  if (prl>0)
    printf("Printing g(Vne)\n");
  contract_dVne(natm, N, nbas, grad_term, Pao, atm, bas, env);
  for (int i = 0; i < 3 * natm; i++)
  {
    grad_term[i] *= 1.;
    grad[i] += grad_term[i];
  }
  if (prl>0)
  for (int i = 0; i < natm; i++)
  {
    for (int j = 0; j < 3; j++)
      printf("%12.8f ",grad_term[i*3+j]);
    printf("\n");
  }
  printf("\n");

  if (prl>0)
    printf("Printing g(PQ)\n");
  contract_d2c2e(natm, N, nbas, Naux, nbas_ri, grad_term, gRS, atm, bas, env);
  for (int i = 0; i < 3 * natm; i++)
  {
    grad_term[i] *= 1.;
    grad[i] += grad_term[i];
  }
  if (prl>0)
  for (int i = 0; i < natm; i++) {
    for (int j = 0; j < 3; j++) {
      printf("%12.8f ",grad_term[i*3+j]);
    }
    printf("\n");
  }
  if (prl>0)
    printf("\n");
  if (prl>0)
    printf("Printing g(gamma3c)\n");
  contract_d3c2e(natm, N, nbas, Naux, nbas_ri, grad_term, gQmunu, atm, bas, env);
  for (int i=0;i<3*natm;i++)
  {
    grad_term[i] *= -2.;
    grad[i] += grad_term[i];
  }
  if (prl>0)
  for (int i=0;i<natm;i++)
  {
    for (int j=0;j<3;j++)
      printf("%12.8f ",grad_term[i*3+j]);
    printf("\n");
  }
  if (prl>0)
    printf("\n");

  if (prl>0)
    printf("Printing grad \n");
  if (prl>0)
  for (int i=0;i<natm;i++)
  {
    for (int j=0;j<3;j++)
      printf("%12.8f ",grad[i*3+j]);
    printf("\n");
  }
  if (prl>0)
    printf("\n");

  delete [] grad_term;

  return;
}

void contract_dS(int natm, int N, int nbas, double* grad_term, double* GF,
                 int* atm, int* bas, double* env)
{
  int N2 = N*N;
  double* jS = new double[3*N2];
  int idx_i = 0;
  int* atm_idx = new int[N];
  int di, dj;
  int shls[2];

  for (int i = 0; i < 3*natm; i++)
    grad_term[i] = 0.;

  for (int i=0;i<nbas;i++)
  {
    int idx_j = 0;
    shls[0] = i;
    di = CINTcgto_spheric(i, bas);
    for (int j = 0; j < di; j++)
      atm_idx[idx_i + j] = bas[ATOM_OF + BAS_SLOTS*i];
    for (int j = 0; j < nbas; j++)
    {
      dj = CINTcgto_spheric(j, bas);
      shls[1] = j;
      double *buf = new double[di*dj*3];
      cint1e_ipovlp_sph(buf,shls,atm,natm,bas,nbas,env);

      for (int x=0;x<3;x++)
      {
        for (int j1=0;j1<dj;j1++)
        {
          for (int i1=0;i1<di;i1++)
          {
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

  for (int x=0;x<3;x++)
  {
    for (int i=0;i<N;i++)
    {
      int atm_i = atm_idx[i];
      for (int j=0;j<N;j++)
      {
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
                 int *atm, int *bas, double *env)
{
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

void get_3center_overlap_sph(double* I, int Nao, int Ng, int natm, int nbas_ao, int nbas_g, int* atm, int* bas, double* env)
{
  if (I==nullptr || atm==nullptr || bas==nullptr || env==nullptr)
    throw invalid_argument("Null pointer in get_3center_overlap_sph");

  const int nbas_total = nbas_ao + nbas_g;

  // Construct AO offsets for both basis spaces.
  vector<int> ao_offset(nbas_ao + 1, 0);
  for (int ish=0;ish<nbas_ao;ish++)
  {
    ao_offset[ish+1] = ao_offset[ish] + CINTcgto_spheric(ish, bas);
  }

  vector<int> g_offset(nbas_g + 1, 0);
  for (int krel=0;krel<nbas_g;krel++)
  {
    const int ksh = nbas_ao + krel;
    g_offset[krel+1] = g_offset[krel] + CINTcgto_spheric(ksh, bas);
  }

  if (ao_offset[nbas_ao]!=Nao)
    throw invalid_argument("Nao is inconsistent with the orbital basis");

  if (g_offset[nbas_g]!=Ng)
    throw invalid_argument("Ng is inconsistent with the g basis");

  fill(I,I+Nao*Nao*Ng,0.);

  // Optional libcint optimizer.
  //CINTOpt* opt = nullptr;
  //int3c1e_optimizer(&opt,atm,natm,bas,nbas_total,env);

  int shls[3];

  for (int ish=0;ish<nbas_ao;ish++)
  {
    shls[0] = ish;

    const int di = CINTcgto_spheric(ish, bas);
    const int mu0 = ao_offset[ish];

    // Only calculate the lower shell triangle because
    // I[mu,nu,n] = I[nu,mu,n].
    for (int jsh=0;jsh<=ish;jsh++)
    {
      shls[1] = jsh;

      const int dj = CINTcgto_spheric(jsh, bas);
      const int nu0 = ao_offset[jsh];

      for (int krel=0;krel<nbas_g;krel++)
      {
        const int ksh = nbas_ao + krel;
        shls[2] = ksh;

        const int dk = CINTcgto_spheric(ksh, bas);
        const int n0 = g_offset[krel];

        /*
         * libcint ordering for a shell block is:
         *   buf[i + di*(j + dj*k)]
         * i.e. the first shell index varies fastest.
         */
        vector<double> buf(di*dj*dk);

        //const CACHE_SIZE_T rc = int3c1e_sph(buf.data(),nullptr,shls,atm,natm,bas,nbas_total,env,opt,nullptr);
        const int nonzero = cint3c1e_sph(buf.data(),shls,atm,natm,bas,nbas_total,env);

        /*
         * A zero return indicates a screened/zero shell block.
         * The output was initialized to zero, so no copy is needed.
         */
        if (nonzero==0)
          continue;

        for (int k=0;k<dk;k++)
        {
          const int n = n0 + k;
          for (int j=0;j<dj;j++)
          {
            const int nu = nu0 + j;
            for (int i=0;i<di;i++)
            {
              const int mu = mu0 + i;
              const size_t ibuf = i + di*(j+dj*k);
              const double value = buf[ibuf];

              const size_t out_mun = (mu*Nao+nu)*Ng + n;
              const size_t out_num = (nu*Nao+mu)*Ng + n;

              I[out_mun] = value;
              I[out_num] = value;
            }
          }
        }
      }
    }
  }

  //CINTdel_optimizer(&opt);
  return;
}

static inline size_t weight_index(int mu, int nu, int n, int Nao, int Naux)
{
  return (mu*Nao + nu) * Naux + n;
}

static inline size_t gradient_index(int atom, int xyz)
{
  return 3*atom + xyz;
}

/*
 * Contract
 *   E = sum_{mu,nu,n} W[mu,nu,n] I[mu,nu,n]
 * and return the explicit center gradient
 *   gradient[3*A+x] = dE/dR_{A,x}.
 */

#if 0
void contract_3center_overlap_gradient(double* gradient, const double* weight, int natoms, int nbas_ao, int nbas_ri, int nao, int naux, int* atm, int* bas, double* env)
{
  const int nbas = nbas_ao+nbas_ri;
  //const int nbas_ri = naux;
  //const int nbas_ao = nao;
  const int first_aux_shell = nbas_ao;

  vector<int> ao_offset(nbas_ao+1,0);
  for (int ish=0;ish<nbas_ao;ish++)
    ao_offset[ish+1] = ao_offset[ish] + CINTcgto_spheric(ish,bas);

  vector<int> aux_offset(nbas_ri+1,0);
  for (int krel=0;krel<nbas_ri;krel++)
  {
    const int ksh = first_aux_shell + krel;
    aux_offset[krel+1] = aux_offset[krel] + CINTcgto_spheric(ksh,bas);
  }

  if (ao_offset[nbas_ao]!=nao)
    throw invalid_argument("nao is inconsistent with bas");

  fill(gradient,gradient+3*natoms,0.);

  /*
   * Loop over the full mu,nu shell square.
   *
   * This supports a general, nonsymmetric W[mu,nu,n]. If W is known to
   * be symmetric in mu and nu, this can be reduced to a shell triangle.
   */
  for (int ish=0;ish<nbas_ao;ish++)
  {
    const int di = CINTcgto_spheric(ish,bas);
    const int mu0 = ao_offset[ish];
    const int atom_i = bas[ish * BAS_SLOTS + ATOM_OF];

    for (int jsh=0;jsh<nbas_ao;jsh++)
    {
      const int dj = CINTcgto_spheric(jsh,bas);
      const int nu0 = ao_offset[jsh];
      const int atom_j = bas[jsh * BAS_SLOTS + ATOM_OF];

      for (int krel=0;krel<nbas_ri;krel++)
      {
        const int ksh = first_aux_shell + krel;
        const int dk = CINTcgto_spheric(ksh,bas);
        const int n0 = aux_offset[krel];
        const int atom_k = bas[ksh * BAS_SLOTS + ATOM_OF];

        const size_t block_size = di * dj * dk;

        /*
         * Libcint derivative components are stored component-major:
         *   buf[xyz*block_size + i + di*(j + dj*k)]
         */
        vector<double> d_i(3*block_size);
        vector<double> d_j(3*block_size);
        vector<double> d_k(3*block_size);

        /*
         * d_i: electronic-coordinate gradient on chi_i.
         */
        {
          int shls[3] = {ish, jsh, ksh};

          const int nonzero = cint3c1e_ip1_sph(d_i.data(),shls,atm,natoms,bas,nbas,env);

          if (!nonzero)
            fill(d_i.begin(),d_i.end(),0.);
        }

        /*
         * d_j: place j first, giving gradient on chi_j.
         *
         * Returned dimensions are dj x di x dk, so index it with
         * j as the first local index and i as the second.
         */
        {
          int shls[3] = { jsh, ish, ksh };

          const int nonzero = cint3c1e_ip1_sph(d_j.data(),shls,atm,natoms,bas,nbas,env);

          if (!nonzero)
            fill(d_j.begin(), d_j.end(), 0.0);
        }

        /*
         * d_k: place k first, giving gradient on g_k.
         *
         * Returned dimensions are dk x di x dj.
         */
        {
          int shls[3] = {ksh, ish, jsh};

          const int nonzero = cint3c1e_ip1_sph(d_k.data(),shls,atm,natoms,bas,nbas,env);

          if (!nonzero)
            fill(d_k.begin(), d_k.end(), 0.0);
        }

        for (int k=0;k<dk;k++)
        {
          const int n = n0 + k;
          for (int j=0;j<dj;j++)
          {
            const int nu = nu0 + j;
            for (int i=0;i<di;i++)
            {
              const int mu = mu0 + i;
              const double w = weight[weight_index(mu,nu,n,nao,naux)];

              if (w==0.)
                continue;

              const size_t idx_i = i + di * (j + dj * k);
              const size_t idx_j = j + dj * (i + di * k);
              const size_t idx_k = k + dk * (i + di * j);

              for (int xyz=0;xyz<3;xyz++)
              {
                /*
                 * Center derivative:
                 *
                 *   d/dR chi(r-R) = -nabla_r chi(r-R).
                 */
                gradient[gradient_index(atom_i, xyz)] -= w * d_i[ xyz * block_size + idx_i];
                gradient[gradient_index(atom_j, xyz)] -= w * d_j[ xyz * block_size + idx_j];
                gradient[gradient_index(atom_k, xyz)] -= w * d_k[ xyz * block_size + idx_k];
              }
            }
          }
        }
      }
    }
  }
}
#endif

void contract_3center_overlap_gradient(double* gradient, const double* weight,
                                       int natoms, int nbas_ao, int nbas_ri,
                                       int nao, int naux,
                                       int* atm, int* bas, double* env)
{
  const int nbas = nbas_ao + nbas_ri;
  const int first_aux_shell = nbas_ao;

  vector<int> ao_offset(nbas_ao+1,0);
  for (int ish=0;ish<nbas_ao;ish++)
    ao_offset[ish+1] = ao_offset[ish] + CINTcgto_spheric(ish,bas);

  vector<int> aux_offset(nbas_ri+1,0);
  for (int krel=0;krel<nbas_ri;krel++)
  {
    const int ksh = first_aux_shell + krel;
    aux_offset[krel+1] = aux_offset[krel] + CINTcgto_spheric(ksh,bas);
  }

  if (ao_offset[nbas_ao]!=nao)
    throw invalid_argument("nao is inconsistent with bas");

  if (aux_offset[nbas_ri]!=naux)
    throw invalid_argument("naux is inconsistent with bas");

  fill(gradient,gradient+3*natoms,0.);

  /*
   * Loop over the full mu,nu shell square.
   * This supports a general, nonsymmetric W[mu,nu,n].
   */
  for (int ish=0;ish<nbas_ao;ish++)
  {
    const int di = CINTcgto_spheric(ish,bas);
    const int mu0 = ao_offset[ish];
    const int atom_i = bas[ish*BAS_SLOTS+ATOM_OF];

    for (int jsh=0;jsh<nbas_ao;jsh++)
    {
      const int dj = CINTcgto_spheric(jsh,bas);
      const int nu0 = ao_offset[jsh];
      const int atom_j = bas[jsh*BAS_SLOTS+ATOM_OF];

      for (int krel=0;krel<nbas_ri;krel++)
      {
        const int ksh = first_aux_shell + krel;
        const int dk = CINTcgto_spheric(ksh,bas);
        const int n0 = aux_offset[krel];
        const int atom_k = bas[ksh*BAS_SLOTS+ATOM_OF];

        const size_t block_size = (size_t)di*dj*dk;

        /*
         * d_i contains
         *   <(nabla chi_i) chi_j g_k>
         * d_j contains
         *   <(nabla chi_j) chi_i g_k>
         * The derivative with respect to the auxiliary center is obtained
         * from translational invariance:
         *   dI/dR_k = d_i + d_j.
         */
        vector<double> d_i(3*block_size,0.);
        vector<double> d_j(3*block_size,0.);

        {
          int shls[3] = {ish,jsh,ksh};

          const int nonzero = cint3c1e_ip1_sph(d_i.data(),shls,atm,natoms,bas,nbas,env);

          if (!nonzero)
            fill(d_i.begin(),d_i.end(),0.);
        }

        /*
         * Put j first so that ip1 differentiates chi_j.
         *
         * Returned dimensions are dj x di x dk.
         */
        {
          int shls[3] = {jsh,ish,ksh};

          const int nonzero =
            cint3c1e_ip1_sph(d_j.data(),shls,atm,natoms,bas,nbas,env);

          if (!nonzero)
            fill(d_j.begin(),d_j.end(),0.);
        }

        for (int k=0;k<dk;k++)
        {
          const int n = n0+k;
          for (int j=0;j<dj;j++)
          {
            const int nu = nu0+j;
            for (int i=0;i<di;i++)
            {
              const int mu = mu0+i;
              const double w = weight[weight_index(mu,nu,n,nao,naux)];

              if (w==0.)
                continue;

              const size_t idx_i = i + (size_t)di*(j + (size_t)dj*k);
              const size_t idx_j = j + (size_t)dj*(i + (size_t)di*k);

              for (int xyz=0;xyz<3;xyz++)
              {
                const double gi = d_i[(size_t)xyz*block_size+idx_i];
                const double gj = d_j[(size_t)xyz*block_size+idx_j];

                /*
                 * Libcint ip1 gives an electronic-coordinate derivative.
                 * For a Gaussian centered at R:
                 *   d/dR chi(r-R) = -nabla_r chi(r-R).
                 */
                gradient[gradient_index(atom_i,xyz)] -= w*gi;
                gradient[gradient_index(atom_j,xyz)] -= w*gj;

                /*
                 * Translational invariance:
                 *   dI/dR_i + dI/dR_j + dI/dR_k = 0
                 * and therefore
                 *   dI/dR_k = gi + gj.
                 */
                gradient[gradient_index(atom_k,xyz)] += w*(gi+gj);
              }
            }
          }
        }
      }
    }
  }
}
