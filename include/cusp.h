#ifndef CUSPS
#define CUSPS

//#include "hartree.h"
#include "read.h"
#include "write.h"
#include "cuda_util.h"
#include "cpu_util.h"
#include "becke.h"
#include "integrals.h"
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


void compute_diatomic_symm(int natoms, int* atno, vector<vector<double> > basis, vector<double*>& pB_all, int prl);

void compute_cusp(int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* pB1, double* pB2, int prl);
void compute_cusp(int natoms, int* atno, double* coords, vector<vector<double> > &basis, vector<double*> pB, int prl);
void compute_Pc(int natoms, int N, double* pB, double* Pc);
void compute_Pc(int natoms, int N, vector<double*> pB, vector<double*> Pc);
int prepare_PSP(double thresh, int natoms, int N, double* S, double* Pc, double* Xp, cusolverDnHandle_t cu_hdl, int prl);
int prepare_PSP(double thresh, int natoms, int N, double* S, vector<double*> Pc, double* Xp, cusolverDnHandle_t cu_hdl, int prl);
int prepare_PSP(int natoms, int N, double* S, double* Pc, double* Xp, cusolverDnHandle_t cu_hdl, int prl);
int prepare_PSP(int natoms, int N, double* S, vector<double*> Pc, double* Xp, cusolverDnHandle_t cu_hdl, int prl);

//applies Pc to S matrix
void project_S(int N, double* S, vector<double*> Pc_all, double* X);

void check_cusp(int No, int natoms, int* atno, vector<vector<double> >& basis, double* pB, double* jCA);
void check_cusp(int No, int natoms, int* atno, vector<vector<double> >& basis, vector<double*> pB, double* jCA);

void evaluate_cusps(int nrad, int nang, int natoms, int* atno, double* coords, vector<vector<double> > basis, double* Pao, float* grid);

//vxc.cpp:
void rgrid_one_atom(int nrad, int Z1, float* r1);

#endif
