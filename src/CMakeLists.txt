add_subdirectory(integrals)
add_subdirectory(libio)
add_subdirectory(librel)
add_subdirectory(libcint)

if(DO_GTO)
  add_subdirectory(libcintw)
endif()

