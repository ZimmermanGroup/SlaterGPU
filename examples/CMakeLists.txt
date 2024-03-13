add_subdirectory(geom_timings)
add_subdirectory(rel_hamiltonian)

if(DO_GTO)
  add_subdirectory(test_gto_grid)
endif()
