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

using namespace std;

//#define A2B 1.8897261

// #define SSTRF( x ) ( ( static_cast<std::ostringstream&>(std::ostringstream() << std::fixed << std::setprecision(8) << std::scientific << (x)) ).str() )
// #define SSTRF2( x ) ( ( static_cast<std::ostringstream&>(std::ostringstream() << std::fixed << std::setprecision(14) << std::scientific << (x)) ).str() )
string SSTRF(float x);
string SSTRF2(double x);

vector<string> split1(const string &s, char delim);

int check_file(string filename);

int initialize(bool gbasis, vector<vector<double> >& basis, vector<vector<double> >& basis_aux, int* atno, double* &coords, int& charge, int& unpaired, double& Enn, int prl);
bool get_secondary_basis(string name, int natoms, int* atno, double* coords, vector<vector<double> >& basis, vector<vector<double> >& basis_aux, int prl);
void create_basis_aux(int natoms, vector<vector<double> >& basis_std, vector<vector<double> >& basis_aux);
string read_basis_text(string aname);
void print_basis(int natoms, vector<vector<double> >& basis, vector<vector<double> >& basis_aux, int prl);
int read_rotate(int N, double* jCA);
double read_float(string filename);
double read_float_2(string filename);
int read_int(string filename);
int read_int_2(string filename);
bool read_array(int size, double* A, string filename);
vector<double> read_vector(string filename);
void read_thresh(float& no_thresh, float& occ_thresh);
int read_opt();
bool read_dft(int& dt1, int& dt2, int type);
int read_ri();
int read_symm();
int read_mo_symmetry(int N, int* symm);
int read_pp();
int read_vnuc();
int read_hfx();
int read_cusp();
int read_esci();
int read_mbe();
int read_hf();
int read_lag();
int read_basis();
void read_T(int& use_td, int& use_tl, int& use_erf, int& use_th, int& use_mixed_t);
int read_cas();
void read_cas_size(int& Nc, int& Na, int& Nb);
bool read_cas_act(int& N, int& M);

void read_nrad_nang(int& nrad, int& nang, int type);
void read_gridps(int& nmu, int& nnu, int& nphi, int type);
void read_quad(int& qo1, int& qo2);

int read_spinref();
int read_group();
int read_restart();
int read_nsteps();
void read_eps(double& eps1, double& eps2);
void read_eps(double& eps1, double& eps2, double& eps1s);
int read_tuples(vector<vector<int> >& tuples);
void read_SENT(string dirname, int N, double* S, double* T, double* En, int prl);
void read_Ciap(int N, int Naux, double* Ciap, string filename);
bool read_integrals(string dirname, int N, int Naux, double* S, double* T, double* En, double* A, double* Ciap, int prl);
bool read_yukawa_potentials(int N, int Naux, double*& Ayd, double*& Cyd);
int read_iarray(short type1, short type2, short i1, int s1, int s2, double* A);
int read_gridpts(int s1, int s2, float* A, string filename);
int read_gridpts(int s1, int s2, double* A, string filename);
int read_square_check(int N, double* Pao, string filename);
int read_square(int N, double* Pao, string filename);
int read_square(vector<vector<double> > basis, double* Pao, string filename);

bool check_PS();

string get_aname(int Z);

#endif
