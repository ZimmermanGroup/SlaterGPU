#ifndef WRITEH
#define WRITEH

#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>
#include <vector>

#include "read.h"

using namespace std;

void write_iarray(short type1, short type2, short i1, int s1, int s2, FP1* A);
void write_iarray(short type1, short type2, short i1, int s1, int s2, FP2* A);

void save_geoms(int natoms, int* atno, vector<FP1> E, vector<FP1*> geoms, string fname);

void write_gridpts(int s1, int s2, FP1* A, string filename);

void write_square(int N, FP1* A, string fname, int prl);
void write_square(int N, FP2* A, string fname, int prl);

void write_molden(bool gbasis, int natoms, int* atno, FP2* coords, vector<vector<FP2> > &basis, FP2* jCA, int No, string fname);
void write_molden_g(int natoms, int* atno, FP2* coords, vector<vector<FP2> > &basis, FP2* jCA, int No, string fname);
string get_aname(int Z);

void write_C(int Naux, int N2, FP1* C);
void write_C(int Naux, int N2, FP2* C);

void write_S_En_T(int N, FP1* S, FP1* En, FP1* T);
void write_S_En_T(int N, FP2* S, FP2* En, FP2* T);

#endif
