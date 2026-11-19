#ifndef SLATERGPU_CINT_SOC_H
#define SLATERGPU_CINT_SOC_H

#include <cstddef>
#include <vector>

namespace cintw {

struct SOCOptions {
  int threads = 0; // 0 uses the OpenMP runtime setting.
  std::size_t tile_cache_bytes = 2 * 1024 * 1024; // Per thread; excludes libcint scratch.
};

struct SOCIntegrals {
  int nao = 0;
  // Three real, row-major AO matrices, component order x,y,z. No i or alpha^2
  // prefactor. H_SO = i alpha^2/4 sum_k A[k] sigma[k].
  std::vector<double> nuclear;
  std::vector<double> somf;
  std::size_t shell_quartets = 0;
  double nuclear_antisymmetry_error = 0.;
  double somf_antisymmetry_error = 0.;
};

// All-electron Breit-Pauli nuclear SOC and spin-averaged SOMF. atm/bas/env use
// libcint's layout. nbas counts primary shells only (no auxiliary shells).
// cartesian is explicit and does not change BT::DO_CART. If reference_density
// is non-null it is a symmetric, spin-summed AO density P, including frozen
// cores, in the same AO ordering; the returned SOMF is J - 3/2 (K1 + K2).
// A null density requests nuclear SOC only. No four-index AO tensor is stored.
// Fractional/ECP nuclear charges and range-separated Coulomb operators are not
// supported. Calls must not mutate atm/bas/env/density during execution.
SOCIntegrals breit_pauli_soc(int natm, int nbas, int* atm, int* bas, double* env,
                            bool cartesian, const double* reference_density = nullptr,
                            const SOCOptions& options = SOCOptions());

} // namespace cintw
#endif
