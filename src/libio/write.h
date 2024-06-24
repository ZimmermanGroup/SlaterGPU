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
//#include "vcf.h"

using namespace std;

void write_iarray(short type1, short type2, short i1, int s1, int s2, float* A);
void write_iarray(short type1, short type2, short i1, int s1, int s2, double* A);

void save_geoms(int natoms, int* atno, vector<float> E, vector<float*> geoms, string fname);

void write_gridpts(int s1, int s2, float* A, string filename);

void write_float(double val, string filename);
void write_vector(int size, double* vals, string filename);
void write_vector(int size, float* vals, string filename);

void write_square(int N, float* A, string fname, int prl);
void write_square(int N, double* A, string fname, int prl);

void write_molden(bool gbasis, int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* jCA, int No, string fname);
void write_molden_g(int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* jCA, int No, string fname);
//void write_molden_vcf(int natoms, int* atno, float* coordsf, vector<vector<vcf> > &vcfs, string fname);
//void write_molden_vcf(int natoms, int* atno, double* coords, vector<vector<vcf> > &vcfs, string fname);

string get_aname(int Z);

void write_C(int Naux, int N2, float* C);
void write_C(int Naux, int N2, double* C);
void write_Cy(int Naux, int N2, double* C);

void write_S_En_T(int N, float* S, float* En, float* T);
void write_S_En_T(int N, double* S, double* En, double* T);

#endif
