#include "gto_grid.h"
#include "cintwrapper.h"
#include "cintprep.h"
#include "omp.h"
#include <math.h>

using namespace std;

void setup_env_bas_grid(
  int * bas_g, int * atm_g, double ** env_g, int & nenv_g,
  const int * bas, const int * atm, const double * env,
  int nbas, int natm, int nenv, int bas_idx
) {
  nenv_g = PTR_ENV_START;
  int offset = PTR_ENV_START;
  atm_g[CHARGE_OF ] = 1;
  atm_g[PTR_COORD ] = offset; // coordinates of grid point
  atm_g[NUC_MOD_OF] = 0;
  atm_g[PTR_ZETA  ] = 0;
  offset += 3;
  int bas_atm_of = bas[BAS_SLOTS * bas_idx + ATOM_OF];
  atm_g[CHARGE_OF  + ATM_SLOTS] = atm[ATM_SLOTS * bas_atm_of + CHARGE_OF];
  atm_g[PTR_COORD  + ATM_SLOTS] = offset;  // coordinates of atom
  atm_g[NUC_MOD_OF + ATM_SLOTS] = 0;
  atm_g[PTR_ZETA   + ATM_SLOTS] = 0;
  offset += 3;

  bas_g[ATOM_OF  ] = 0;
  bas_g[ANG_OF   ] = 0;
  bas_g[NPRIM_OF ] = 1;
  bas_g[NCTR_OF  ] = 1;
  bas_g[PTR_EXP  ] = offset;
  bas_g[PTR_COEFF] = offset + 1;
  offset += 2;
  bas_g[BAS_SLOTS + ATOM_OF  ] = 1;
  bas_g[BAS_SLOTS + ANG_OF   ] = bas[bas_idx * BAS_SLOTS + ANG_OF  ];
  bas_g[BAS_SLOTS + NPRIM_OF ] = bas[bas_idx * BAS_SLOTS + NPRIM_OF];
  bas_g[BAS_SLOTS + NCTR_OF  ] = bas[bas_idx * BAS_SLOTS + NCTR_OF ];
  bas_g[BAS_SLOTS + PTR_EXP  ] = offset;
  bas_g[BAS_SLOTS + PTR_COEFF] = offset + bas_g[BAS_SLOTS + NPRIM_OF];
  offset += bas_g[BAS_SLOTS + NPRIM_OF] + bas_g[BAS_SLOTS + NPRIM_OF] * bas_g[BAS_SLOTS + NCTR_OF];

  nenv_g += offset;
  *env_g = new double[nenv_g]();
  offset = PTR_ENV_START;
  offset += 3; // First three coordinates are set to grid point
  int atm_idx = PTR_COORD + ATM_SLOTS * bas[BAS_SLOTS * bas_idx + ATOM_OF];
  for (int i = 0; i < 3; i++) {
    (*env_g)[offset++] = env[atm[atm_idx]+i];
  }

  // Getting basis set specification
  // basis set for large decay
  const double DIRAC_EXP = 1.e15;
  (*env_g)[offset++] = DIRAC_EXP;
  // const double DIRAC_NORM = CINTgto_norm(0,DIRAC_EXP);
  // const double DIRAC_NORM = pow(DIRAC_EXP/M_PI,3./2.);
  const double DIRAC_NORM = 2*pow(DIRAC_EXP,1.5)/(M_PI);    // M_PI is the constant pi
  //const double DIRAC_NORM = pow((2*DIRAC_EXP/M_PI),3/4);
  (*env_g)[offset++] = DIRAC_NORM; // * DIRAC_NORM;

  int nprim = bas_g[BAS_SLOTS + NPRIM_OF];
  for (int i = 0; i < nprim; i++) {
    (*env_g)[offset + i] = env[bas[bas_idx * BAS_SLOTS + PTR_EXP] + i];
  }
  for (int i = 0; i < nprim; i++) {
    (*env_g)[offset + nprim + i] = env[bas[bas_idx * BAS_SLOTS + PTR_COEFF] + i];
  }
  offset += 2 * nprim;

  return;
}

template <class T>
void gen_gto_on_grid(
  const int nbas, const int natm, const int nenv, 
  const int * bas, const int * atm, const double * env,
  T * grid, int grid_size, int N, int bas_idx, T ** gto_on_grid,
  int & shl_size
) {
  int di = CINTcgto_spheric(bas_idx,bas);
  int dj = 1;
  shl_size = di;
  double * buf = new double[di * dj]();

  int shls[2] = {0,1};
  int bas_g[2 * BAS_SLOTS] = {0};
  int atm_g[2 * ATM_SLOTS] = {0};
  double * env_g = NULL;
  int nenv_g;

  setup_env_bas_grid(
    bas_g, atm_g, &env_g, nenv_g, bas, atm, env,
    nbas, natm, nenv, bas_idx
  );

  int offset = PTR_ENV_START;
  for (int i = 0; i < grid_size; i++) {
    env_g[offset + 0] = (double)grid[i * 6 + 0];
    env_g[offset + 1] = (double)grid[i * 6 + 1];
    env_g[offset + 2] = (double)grid[i * 6 + 2];

    cint1e_ovlp_sph(buf,shls,atm_g,2,bas_g,2,env_g);
    for (int j = 0; j < di * dj; j++) {
      gto_on_grid[j][i] = (T) buf[j];
    }
  }

  delete [] env_g;
  delete [] buf;
}

template <class T>
void gen_gto_grad_on_grid(
  const int nbas, const int natm, const int nenv, 
  const int * bas, const int * atm, const double * env,
  T * grid, int grid_size, int N, int bas_idx, T ** grad_on_grid,
  int & shl_size
) {
  int di = CINTcgto_spheric(bas_idx,bas);
  printf("di: %i\n",di);
  int dj = 1;
  shl_size = di;
  double * buf = new double[di * dj * 3]();
  float x = di * dj * 3;
  printf("size of buf : %8.5f\n",x);

  int shls[2] = {0,1};
  int bas_g[2 * BAS_SLOTS] = {0};
  int atm_g[2 * ATM_SLOTS] = {0};
  double * env_g = NULL;
  int nenv_g;

  setup_env_bas_grid(
    bas_g, atm_g, &env_g, nenv_g, bas, atm, env,
    nbas, natm, nenv, bas_idx
  );

  int offset = PTR_ENV_START;
  printf("gridsize: %i\n",grid_size);
  for (int i = 0; i < grid_size; i++) {
    env_g[offset + 0] = (double)grid[i * 6 + 0];
    env_g[offset + 1] = (double)grid[i * 6 + 1];
    env_g[offset + 2] = (double)grid[i * 6 + 2];
  
    cint1e_ovlpip_sph(buf,shls,atm_g,2,bas_g,2,env_g);
    printf("cint1e thing ran with i = %i\n",i);
    for (int j = 0; j < di * dj; j++) {
      for (int k = 0; k < 3; k++) {
        // printf("%3d %3d %3d\n",i,(j * grid_size * 3) + i * 3 + k,k*di+j);
        grad_on_grid[j][i * 3 + k] = (T) buf[k*di + j];
      }
    }
  }

  delete [] env_g;
  delete [] buf;
}


template void gen_gto_grad_on_grid<double>(
  const int nbas, const int natm, const int nenv, 
  const int * bas, const int * atm, const double * env,
  double * grid, int grid_size, int N, int bas_idx, double ** grad_on_grid,
  int & shl_size
);
template void gen_gto_grad_on_grid<float>(
  const int nbas, const int natm, const int nenv, 
  const int * bas, const int * atm, const double * env,
  float * grid, int grid_size, int N, int bas_idx, float ** grad_on_grid,
  int & shl_size
);
template void gen_gto_on_grid<double>(
  const int nbas, const int natm, const int nenv, 
  const int * bas, const int * atm, const double * env,
  double * grid, int grid_size, int N, int bas_idx, double ** gto_on_grid,
  int & shl_size
);
template void gen_gto_on_grid<float>(
  const int nbas, const int natm, const int nenv, 
  const int * bas, const int * atm, const double * env,
  float * grid, int grid_size, int N, int bas_idx, float ** gto_on_grid,
  int & shl_size
);
