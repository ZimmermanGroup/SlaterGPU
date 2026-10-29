#ifndef GAUSS_H
#define GAUSS_H

#include "integrals.h"
//#include "vxc.h"
#include "gamma.h"
#include "read.h"
#include "gto_grid.h"
#include "print.h"
#include "cuda_util.h"

void eval_gh(int gs, float* grid, float* val1, int l1, int m1, const float norm1, const float zeta1);
void eval_ghd(int gs, double* grid, double* val1, int l1, int m1, const double norm1, const double zeta1);
void eval_gh_spd(bool include_r, int gs, double* grid, double* val1, int l1, int m1);
void eval_gh_spd(int gs, double* grid, double* val1, int l1, int m1);

void eval_gh_ke(int gs, float* grid, float* val1, int n1, int l1, const float norm1, const float zeta1);
void eval_ghd_ke(int gs, double* grid, double* val1, int n1, int l1, const double norm1, const double zeta1);
//void eval_pdke_gh(int gs, double* grid, double* val1, int n1, int l1, int m1, double norm1, double zeta1);
void eval_p_gh(int gs, float* grid, float* val, int n1, int l1, int m1, float norm1, float zeta1);
void eval_pd_gh(int gs, double* grid, double* val, int n1, int l1, int m1, double norm1, double zeta1);

//Gaussian single-center potentials
void eval_vghd(int gsa, double* grid, double* val, int n, int l, double norm, double zt);

void eval_hess_ghd(int gs, double* grid, double* val, int n1, int l1, int m1, double norm1, double zeta1);
int eval_gh_full(int gs, float* grid, float** val1, int i1, int natoms, int nbas, int nenv, int N, int* atm, int* bas, double* env);

void wf_to_grid_gh(bool divide_by_rho, bool calc_rho, int natoms, int* atno, double* coords, vector<vector<double> > basis,
                   double* Pao, double* gfao, int gs, double* grid, double* wt,
                   double* rho, double* gf, double* Td, double* TL, int prl);
//void wf_to_grid_gh_ke(int natoms, int* atno, double* coords, vector<vector<double> > basis, double* jCA, int gs, float* grid, float* wt, double* TL, int prl);
//void wf_to_grid_gh_ke_2(int natoms, int* atno, double* coords, vector<vector<double> > basis, int nbas, int nenv, int N, int* atm, int* bas, double* env,
//                        double* jCA, int gs, float* grid, float* wt, double* Td, int prl);
//void integrate_hole_para_gh(double* rdm, bool full_rdm, bool hfx_on, int Nc, int No, int M, int natoms, int* atno, double* coords, int gs, int gsb, vector<vector<double> >& basis,
//                    double* Pao, double* Pmo, double* jCA, float* grid, float* gridb, float* wt, float* wtb, float* rho, float* vxch, int prl);

void gen_4cj_overlap(double* ol4j, int nrad, int nang, double* ang_g, double* ang_w,
                     int natoms, int* atno, double* coords, vector<vector<double> > basis, int prl);

#endif
