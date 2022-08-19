#ifndef READH
#define READH

#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>
#include <algorithm>
#include <math.h>
#include <vector>
#include "fp_def.h"

using namespace std;

//#define A2B 1.8897261

string read_basis_text(string aname);
void read_rotate(int N, FP2* jCA);
FP2 read_FP1(string filename);
void read_thresh(FP1& no_thresh, FP1& occ_thresh);
int read_opt();
bool read_dft(int& dt1, int& dt2);
int read_ri();
int read_esci();
int read_mbe();
int read_lag();
int read_basis();
int read_cas();
bool read_cas_act(int& N, int& M);
bool read_nalpha_nbeta(int& Na, int& Nb);
int read_spinref();
int read_group();
int read_restart();
int read_nsteps();
void read_eps(FP2& eps1, FP2& eps2);
void read_eps(FP2& eps1, FP2& eps2, FP2& eps1s);
int read_tuples(vector<vector<int> >& tuples);
void read_SENT(string dirname, int N, FP2* S, FP2* T, FP2* En, int prl);
void read_integrals(string dirname, int N, int Naux, FP2* S, FP2* T, FP2* En, FP2* A, FP2* Ciap, int prl);
int read_iarray(short type1, short type2, short i1, int s1, int s2, FP2* A);
int read_gridpts(int s1, int s2, FP1* A, string filename);
int read_square(int N, FP2* Pao, string filename);
int read_square(vector<vector<FP2> > basis, FP2* Pao, string filename);
void print_rectangle_sm(int N1, int N2, FP2* S);
void print_square(int M, int N, FP2* S);
void print_square(int M, int N, FP1* S);
void print_square(int N, FP2* S);
void print_square(int N, FP1* S);

string get_aname(int Z);
int initialize(bool gbasis, vector<vector<FP2> >& basis, vector<vector<FP2> >& basis_aux, int* atno, FP2* &coords, int& charge, int& unpaired, FP2& Enn, int prl);

void read_nrad_nang(int& nrad, int& nang, int type);
FP2 nuclear_repulsion(int natoms, int* atno, FP2* coords);
FP2 nuclear_repulsion(int natoms, int* atno, FP1* coordsf);

#endif
