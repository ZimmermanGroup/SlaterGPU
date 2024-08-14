#include "gto_grid.h"
#include "cintwrapper.h"

void compute_ovlp_grid(
  int bas_idx, int * bas0, double * env0,
  float * grid, double * gto_vals, int grid_size, int atom, double * at_coord
) {
  int shls[2] = {0,1};
  int di = BT::DO_CART ? CINTcgto_cart(bas_idx,bas0) : CINTcgto_spheric(bas_idx,bas0);
  int dj = 1;
  double * buf = new double[di*dj];

  int nenv = PTR_ENV_START;
  int offset = PTR_ENV_START;
  int bas[2 * BAS_SLOTS] = {0};
  int atm[2 * ATM_SLOTS] = {0};
  double * env;

  atm[0] = 1;
  atm[1] = offset; // coordinates of atom 1
  atm[2] = 0;
  atm[3] = 0;
  offset += 3;
  atm[0 + ATM_SLOTS] = atom;
  atm[1 + ATM_SLOTS] = offset;
  atm[2 + ATM_SLOTS] = 0;
  atm[3 + ATM_SLOTS] = 0;
  offset += 3;

  bas[ATOM_OF  ] = 0;
  bas[ANG_OF   ] = 0;
  bas[NPRIM_OF ] = 1;
  bas[NCTR_OF  ] = 1;
  bas[PTR_EXP  ] = offset;
  bas[PTR_COEFF] = offset + 1;
  offset += 2;
  bas[BAS_SLOTS + ATOM_OF  ] = 1;
  bas[BAS_SLOTS + ANG_OF   ] = bas0[bas_idx * BAS_SLOTS + ANG_OF  ];
  bas[BAS_SLOTS + NPRIM_OF ] = bas0[bas_idx * BAS_SLOTS + NPRIM_OF];
  bas[BAS_SLOTS + NCTR_OF  ] = bas0[bas_idx * BAS_SLOTS + NCTR_OF ];
  bas[BAS_SLOTS + PTR_EXP  ] = offset;
  bas[BAS_SLOTS + PTR_COEFF] = offset + bas[BAS_SLOTS + NPRIM_OF ];
  offset += bas[BAS_SLOTS + NPRIM_OF] + bas[BAS_SLOTS + NPRIM_OF] * bas[BAS_SLOTS + NCTR_OF];

  nenv += offset;
  env = new double[nenv]();
  offset = PTR_ENV_START;
  offset += 3; // need to insert grid coordinates to first 3 here
  for (int i = 0; i < 3; i++) {
    env[offset++] = at_coord[i];
  }
  
  // Getting basis set specification
  // basis set for large decay
  env[offset ++] = 1e9;
  env[offset ++] = 1.;

  for (int i = 0; i < bas[BAS_SLOTS + NPRIM_OF ]; i++) {
    env[offset++] = env0[bas0[bas_idx * BAS_SLOTS + PTR_EXP ]+i];
  }
  for (int i = 0; i < bas[BAS_SLOTS + NPRIM_OF] * bas[BAS_SLOTS + NCTR_OF]; i++) {
    env[offset++] = env0[bas0[bas_idx * BAS_SLOTS + PTR_COEFF] + i];
  }
  
  // assumes gto_grid has appropriate size for number of subshells
  for (int i = 0; i < grid_size; i ++) {
    env[PTR_ENV_START + 0] = grid[i*6 + 0];
    env[PTR_ENV_START + 1] = grid[i*6 + 1];
    env[PTR_ENV_START + 2] = grid[i*6 + 2];
    if (BT::DO_CART) {
      cint1e_ovlp_cart(buf,shls,atm,2,bas,2,env);
    }
    else {
      cint1e_ovlp_sph(buf,shls,atm,2,bas,2,env);
    }
    for (int j = 0; j < di; j++) {
      gto_vals[j * grid_size + i] = buf[j];
    }
  }

  delete [] env;
  delete [] bas;
  delete [] buf;
  delete [] atm;
  return;
}
