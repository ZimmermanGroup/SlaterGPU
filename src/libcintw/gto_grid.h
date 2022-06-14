#ifndef _GTO_GRID_H_
#define _GTO_GRID_H_

#include <vector>
#include <string>
#include <map>
#include "fp_def.h"

using namespace std;

extern "C" {
  #include "cint_funcs.h"
  int cint1e_ovlp_cart(
    double *buf, int *shls,int *atm, int natm, 
    int *bas, int nbas, double *env
  );
  int cint1e_ovlp_sph(
    double *buf, int *shls,int *atm, int natm, 
    int *bas, int nbas, double *env
  ); 
  int cint1e_ipovlp_cart(
    double *buf, int *shls, int *atm, int natm, 
    int *bas, int nbas, double *env
  );
  int cint1e_ipovlp_sph(
    double *buf, int *shls, int *atm, int natm, 
    int *bas, int nbas, double *env
  );
  double CINTgto_norm(FINT n, double a);
  int CINTtot_cgto_spheric(const int *bas, const int nbas);
  FINT CINTcgto_spheric(const FINT n, const FINT *bas);
  int CINTtot_cgto_cart(const int *bas, const int nbas);
  FINT CINTcgto_cart(const FINT n, const FINT *bas);
}

void compute_ovlp_grid(
  int bas_idx, int * bas, double * env,
  FP1 * grid, double * gto_vals, int grid_size, int atom, double * at_coord
);

#endif // _GRO_GRID_H_
