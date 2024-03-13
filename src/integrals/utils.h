#include <chrono>
#include <string>
using namespace std;

int get_NOs(double thresh, int N, double* rho1, double* rho1e, double* jCA, double* jCA1, int prl);
void ao_to_mo_cpu(int N, double* jAO, double* jMO, double* jCA, double* tmp);
void mo_to_ao_cpu(int N, double* Pmo, double* Pao, double* jCA, double* tmp);

void mat_times_mat_cpu(double* C, double* A, double* B, int N);
void mat_times_mat_cpu(double* C, double* A, double* B, int M, int N, int K);
void mat_times_mat_at_cpu(double* C, double* A, double* B, int N);
void mat_times_mat_bt_cpu(double* C, double* A, double* B, int M, int N, int K, int LDAB);
void mat_times_mat_bt_cpu(double* C, double* A, double* B, int M, int N, int K);
void mat_times_mat_bt_cpu(double* C, double* A, double* B, int N);

void print_duration(chrono::high_resolution_clock::time_point t1, chrono::high_resolution_clock::time_point t2, string name);
void DiagonalizeP(double* A, double* Ae, int N);

void ao_to_mo_4c(double ** eri, double ** jMOI, double * jCA, int N);
