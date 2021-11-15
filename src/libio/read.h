#ifndef READH
#define READH

//#define A2B 1.8897261
#include <string>
#include <vector>
#include "fp_def.h"

using namespace std;

void read_thresh(FP1& no_thresh, FP1& occ_thresh);
int read_opt();
int read_esci();
int read_mbe();
int read_lag();
int read_basis();
int read_cas();
int read_restart();
int read_nsteps();
void read_eps(FP2& eps1, FP2& eps2);
void read_eps(FP2& eps1, FP2& eps2, FP2& eps1s);
int read_tuples(vector<vector<int> >& tuples);
void read_SENT(string dirname, int N, FP2* S, FP2* T, FP2* En, int prl);
void read_integrals(string dirname, int N, int Naux, FP2* S, FP2* T, FP2* En, FP2* A, FP2* Ciap, int prl);
// int read_iarray(short type1, short type2, short i1, int s1, int s2, FP2* A);
int read_square(int N, FP2* Pao, string filename);
int read_square(vector<vector<FP2> > basis, FP2* Pao, string filename);
void print_rectangle_sm(int N1, int N2, FP2* S);
int initialize(bool gbasis, vector<vector<FP2> >& basis, vector<vector<FP2> >& basis_aux, int* atno, FP2* &coords, int& charge, int& unpaired, FP2& Enn, int prl);
void read_nrad_nang(int& nrad, int& nang, int type);

string get_aname(int Z);

#endif
