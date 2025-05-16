#ifndef SCF_UTILH
#define SCF_UTILH

#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>

#include <algorithm>
#include <math.h>

#include "cpu_util.h"
#include "cuda_util.h"

int count_zero_mo(int N, double* jCA);
int core_count(int natoms, int* atno);
int electron_count(int charge, int natoms, int* atno, int& Nc, int& Na, int& Nb);

void mo_coeff_to_pao_occ_complex(double* occ, int N, double* jCA, double* Pao, double* Pai);
void mo_coeff_to_pao_complex(int No, int N, double* jCA, double* Pao, double* Pai);
void ao_to_mo_complex(int No, int N, double* jAO, double* jMO, double* jCA, double* tmp);
void ao_to_mo_occ(int No, int N, double* jAO, double* jMO, double* jCA, double* tmp);
void ao_to_mo(int N, double* jAO, double* jMO, double* jCA, double* tmp);
void ao_to_mo_cpu(int N, double* jAO, double* jMO, double* jCA, double* tmp);

void compute_CSC(int N, double* jCA1, double* S, double* jCA2, double* O);

#endif
