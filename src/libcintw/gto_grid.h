#ifndef __GTO_GRID_H__
#define __GTO_GRID_H__

void setup_env_bas_grid(
  int * bas_g, int * atm_g, double ** env_g, int & nenv_g,
  const int * bas, const int * atm, const double * env,
  int nbas, int natm, int nenv, int bas_idx
);
template <class T>
void gen_gto_on_grid(
  const int nbas, const int natm, const int nenv, 
  const int * bas, const int * atm, const double * env,
  T * grid, int grid_size, int N, int bas_idx, T ** gto_on_grid,
  int & shl_size
);
template <class T>
void gen_gto_grad_on_grid(
  const int nbas, const int natm, const int nenv, 
  const int * bas, const int * atm, const double * env,
  T * grid, int grid_size, int N, int bas_idx, T ** grad_on_grid,
  int & shl_size
);


#endif // __GTO_GRID_H__
