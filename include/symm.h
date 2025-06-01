#ifndef SYMM
#define SYMM

#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <cstdio>
#include <algorithm>
#include <math.h>
#include <vector>

//#include "accel.h"
//#include <omp.h>
//#include "mpih.h"

#include "integrals.h"
#include "read.h"

using namespace std;


bool determine_mo_symmetry(int point_group, int natoms, int* atno, float* coords, vector<vector<double> > basis, double* jCA, double* jS, int* symm, int prl);
bool determine_mo_symmetry(int point_group, int natoms, int* atno, double* coords, vector<vector<double> > basis, double* jCA, double* jS, int* symm, int prl);

void determine_symmetry_atomic(int N, int* symmblocks, vector<vector<double> >& basis, double* jCA, int prl);

void determine_symmetry_aos(int N, int* symb, double* jCA, int prl);



#endif
