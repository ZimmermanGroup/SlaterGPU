#ifndef _PREP_CINT_H_
#define _PREP_CINT_H_

#include <string>
#include <vector>
#include <map>

extern "C" {
  #include "cint_funcs.h"
  int cint1e_nuc_cart(double *buf, int *shls,
                     int *atm, int natm, int *bas, int nbas, double *env);
  int cint1e_kin_cart(double *buf, int *shls,
                     int *atm, int natm, int *bas, int nbas, double *env);
  int cint1e_ovlp_cart(double *buf, int *shls,
                     int *atm, int natm, int *bas, int nbas, double *env);
  int cint2e_cart(double *buf, int *shls,
                  int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *no_opt);
  double CINTgto_norm(FINT n, double a);

  int CINTtot_cgto_cart(const int *bas, const int nbas);
  FINT CINTcgto_cart(const FINT n, const FINT *bas);

  int cint1e_nuc_sph(double *buf, int *shls,
                         int *atm, int natm, int *bas, int nbas, double *env);
  int cint1e_kin_sph(double *buf, int *shls,
                         int *atm, int natm, int *bas, int nbas, double *env);
  int cint1e_ovlp_sph(double *buf, int *shls,
                          int *atm, int natm, int *bas, int nbas, double *env);
  int cint1e_ovlpip_sph(double *buf, int *shls,
                          int *atm, int natm, int *bas, int nbas, double *env);
  int cint2e_sph(double *buf, int *shls,
                     int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *no_opt);
  
  // gradient integrals
  // < -ih∇ | Vnuc | -ih∇ >
  int cint1e_ipnucip_cart(double *buf, int *shls,
                     int *atm, int natm, int *bas, int nbas, double *env);
  int cint1e_ipnucip_sph(double *buf, int *shls,
                         int *atm, int natm, int *bas, int nbas, double *env);
  //
  int cint1e_ipnuc_sph(double *buf, int *shls,
                       int *atm, int natm, int *bas, int nbas, double *env);
  // < i | d 1 / R_1A | j >
  int cint1e_drinv_sph(double *buf, int *shls,
                       int *atm, int natm, int *bas, int nbas, double *env);
  // < di | 1 / R_1A | j >
  int cint1e_iprinv_sph(double *buf, int *shls,
                        int *atm, int natm, int *bas, int nbas, double *env);
  int cint1e_ipkin_sph(double *buf, int *shls,
                       int *atm, int natm, int *bas, int nbas, double *env);
  int cint1e_kinip_sph(double *buf, int *shls,
                       int *atm, int natm, int *bas, int nbas, double *env);
  // < di | j >
  int cint1e_ipovlp_sph(double *buf, int *shls,
                        int *atm, int natm, int *bas, int nbas, double *env);
  // ( di j | P )
  int cint3c2e_ip1_sph(double *buf, int *shls,
                       int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *no_opt);
  // ( i j | dP )
  int cint3c2e_ip2_sph(double *buf, int *shls,
                       int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *no_opt);
  // ( dP | Q )
  int cint2c2e_ip1_sph(double *buf, int *shls,
                       int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *no_opt);
  // ( di j | k l )
  int cint2e_ip1_sph(double *buf, int *shls,
                     int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *no_opt);
  // ( i j | dk l )
  int cint2e_ip2_sph(double *buf, int *shls,
                     int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *no_opt);

  int CINTtot_cgto_spheric(const int *bas, const int nbas);
  FINT CINTcgto_spheric(const FINT n, const FINT *bas);
}

using namespace std;

struct basis_t {
  int nuc; //nuclear charge
  vector< int > shells; // vector containing angular momentum numbers
  vector< vector< double > > exps; // vector containing exponents
  vector< vector< double > > coef; // vector containing coefficients
};

const double ANG2BOHR = 1.8897259886;
const double BOHR2ANG = 0.529177249;

class BT {
  public:
    static bool DO_CART;
};

int calc_di(int i, int *bas);

void get_overlap(double * overlap, int N,
                 int natm, int nbas, int nenv,
                 int *atm, int* bas, double *env);
void get_overlap_ri(double * overlap, int Naux,
                    int natm, int nbas, int nbas_ri, int nenv,
                    int *atm, int* bas, double *env);

void get_hcore(double *hcore, double *En, double *T, int N,
               int natm, int nbas, int nenv,
               int *atm, int *bas, double *env);

void get_hcore(double *hcore, int N,
               int natm, int nbas, int nenv,
               int *atm, int *bas, double *env);

void get_tcore(double *tcore, int N,
               int natm, int nbas, int nenv,
               int *atm, int *bas, double *env);

void gen_pvp(double *pvp, int N,
               int natm, int nbas, int nenv,
               int *atm, int *bas, double *env);

void gen_eri(double **eri, int N,
             int natm, int nbas, int nenv,
             int *atm, int *bas, double *env);

void gen_jMOI_gto(double **eri, int N,
                  int natm, int nbas, int nenv,
                  int *atm, int *bas, double *env);

void gen_ri_int(double *jSaux, double *jSinv, double **jCaux,
                double * jERI3c, double *jERI2c, int Naux, int N,
                int natm, int nbas, int nenv, int nbas_ri,
                int *atm, int *bas, double *env, double ** eri);

void gen_eri_ri(double *jSaux, double *jSinv, double **jCaux,
                double *jERI3c, double *jERI2c, int Naux, int N,
                int natm, int nbas, int nenv, int nbas_ri,
                int *atm, int *bas, double *env, double **eri);

void gen_eri_3c(double *eri, int N, int Naux,
                int natm, int nbas, int nenv, int nbas_ri,
                int *atm, int *bas, double *env);

void gen_eri_2c(double *eri, int Naux,
                int natm, int nbas, int nenv, int nbas_ri,
                int *atm, int *bas, double *env);

void translate_basis(int nbas, int N, int natm, int* atm, int* bas, double* env,
                     vector< vector< double> > &basis);

void compute_g_grad_ri(int natm, int N, int Naux, int nbas, int nbas_ri,
                       double *grad, double *GF, double *Pao, double *gQmunu,
                       double *gRS, int *atm, int *bas, double *env);
//contract functions also compute S,K,... integrals
void contract_dS(int natm, int N, int nbas, double *grad, double *GF,
                 int *atm, int *bas, double *env);
void contract_dK(int natm, int N, int nbas, double *grad_term, double *Pao,
                 int *atm, int *bas, double *env);
void contract_dVne(int natm, int N, int nbas, double *grad_term, double *Pao,
                   int *atm, int *bas, double *env);
void contract_d2c2e(int natm, int N, int nbas, int Naux, int nbas_ri,
                    double *grad_term, double *gRS,
                    int *atm, int *bas, double *env);
void contract_d3c2e(int natm, int N, int nbas, int Naux, int nbas_ri, 
                    double *grad_term, double *gQmunu,
                    int *atm, int *bas, double *env);
#endif

// Start over the gen_ri_int function
// use M_ik^mu formulation for jCaux
