#ifndef JELLIUM_INTSH
#define JELLIUM_INTSH

//#include "hartree.h"
#include "read.h"
#include "write.h"
#include "cuda_util.h"
#include "cpu_util.h"
#include "integrals.h"
#include "becke.h"
#include "gauss.h"
//#include "vcf.h"
#include "scf_util.h"

#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <cstdio>
#include <algorithm>
#include <math.h>
#include <vector>

#include "accel.h"
#include <omp.h>
//#include "mpih.h"

using namespace std;


void print_murak_rmax(double Rc, int nrad, vector<vector<double> > basis);

int get_gs1(int nrad, int nang, float* gridf, double Rc);
int get_gs1(int nrad, int nang, double* grid, double Rc);

void compute_STEn_jellium(bool use_slater, int order, double Zg, double ztg, double rs, double Rc, int Ne, bool update_norm, vector<vector<double> >& basis,
         int nrad, int nang, double* ang_g, double* ang_w, double* S, double* T, double* En, int prl);

void compute_Vr_jellium(int order, double Zg, double ztg, double Rc, int Ne, int gs1, int gs2, double* grid, double* Vr);
void compute_Vr_jellium(double Zg, double ztg, double Rc, int Ne, int gs1, int gs2, double* grid, double* Vr);
void compute_Vr_jellium(double Zg, double ztg, double Rc, int Ne, int gs1, int gs2, float* grid, double* Vr);

void compute_C_jellium(bool slater, double Rc, vector<vector<double> > basis, vector<vector<double> > basis_aux, int nrad, int nang,
                       double* ang_g, double* ang_w, double* C);

void compute_4c_ol_jellium(bool slater, double Rc, vector<vector<double> > basis,
         int nrad, int nang, double* ang_g, double* ang_w, double* ol, int prl);

//overlap.cpp:
void orthonormalize_mos(int Nn, int N, double* S, double* jCA);
//soicas.cpp:
void fix_ortho_cpu(int N, double* jCA, double* S);

void symmetrize_MOs(int natoms, int N1, vector<vector<double> > basis, double* jCA, double* S);

#endif
