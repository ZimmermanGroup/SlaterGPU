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

void compute_cusp(int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* pb, int prl);
void compute_Pc(int natoms, int N, double* pB, double* Pc);
int prepare_PSP(int natoms, int N, double* S, double* Pc, double* Xp, cusolverDnHandle_t cu_hdl, int prl);

void check_cusp(int No, int natoms, int* atno, vector<vector<double> >& basis, double* pB, double* jCA);
void evaluate_cusps(int nrad, int nang, int natoms, int* atno, double* coords, vector<vector<double> > basis, double* Pao, float* grid);

//vxc.cpp:
void rgrid_one_atom(int nrad, int Z1, float* r1);

#endif
