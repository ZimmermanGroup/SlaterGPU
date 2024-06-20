set(INTEGRALS_SOURCES
    integrals.cpp 
    integrals_d.cpp 
    gamma.cpp
    lebedev2.cpp     
    murak.cpp
    pVp.cpp
    reduce.cpp
    spherical.cpp
    Vinr.cpp
    grid_util.cpp
    )

add_library(SlaterGPU STATIC "${INTEGRALS_SOURCES}")
target_include_directories(SlaterGPU PUBLIC
    "${CMAKE_SOURCE_DIR}/include/SlaterGPU"
    "${CMAKE_SOURCE_DIR}/include"
    )

# set(INTEGRALS_SOURCES
#     integrals.cpp 
#     integrals_d.cpp 
#     gamma.cpp
#     lebedev2.cpp     
#     murak.cpp
#     pVp.cpp
#     reduce.cpp
#     spherical.cpp
#     Vinr.cpp
#     grid_util.cpp
#     vcf.cpp
#     opt.cpp
#     becke.cpp
#     soicas.cpp
#     nhbci.cpp
#     cusp.cpp
#     vfit.cpp
#     symm.cpp
#     cuda_util.cpp
#     cas.cpp
#     cpu_util.cpp
#     det.inc
#     mpih.cpp
#     dftj.cpp
#     gauss.cpp
#     orbsteps.cpp
#     overlap.cpp
#     x2c.cpp
#     vxc.cpp
#     hartree.cpp
#     )

# add_library(SlaterGPU STATIC "${INTEGRALS_SOURCES}")
# target_include_directories(SlaterGPU PUBLIC
#     "${CMAKE_SOURCE_DIR}/include/SlaterGPU"
#     "${CMAKE_SOURCE_DIR}/include"
#     "${CMAKE_SOURCE_DIR}/src/libxc/include"
#     )

# find_library(LIBXC_LIB NAMES libxc.a PATHS ${CMAKE_SOURCE_DIR}/src/libxc/lib NO_DEFAULT_PATH)
# if(NOT LIBXC_LIB)
#     message(FATAL_ERROR "libxc.a not found")
# endif()
# target_link_libraries(SlaterGPU PRIVATE ${LIBXC_LIB})


add_definitions("-w")

set_property(TARGET SlaterGPU PROPERTY CXX_STANDARD 11)
set_property(TARGET SlaterGPU PROPERTY POSITION_INDEPENDENT_CODE ON)

set(SLATER_COMP_FLAGS ${CMAKE_CXX_FLAGS} -DRED_DOUBLE=1 -O2
    -g -D__GCC_ATOMIC_TEST_AND_SET_TRUEVAL=1 -DGAMMA=1
    )
if(USE_ACC)
    target_compile_options(SlaterGPU PRIVATE
        ${SLATER_COMP_FLAGS} -DUSE_ACC=1 -acc=gpu
        )
else()
    target_compile_options(SlaterGPU PRIVATE
        ${SLATER_COMP_FLAGS} -DUSE_ACC=0
        )
endif()
