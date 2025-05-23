set(INTEGRALS_SOURCES
    integrals.cpp
    integrals_d.cpp
    integrals_aux.cpp
    integrals_ps.cpp
    spherical.cpp
    sphericald.cpp
    pVp.cpp
    pVpd.cpp
    becke.cpp
    gamma.cpp
    lebedev2.cpp
    murak.cpp
    reduce.cpp
    Vinr.cpp
    quad.cpp
    prosph.cpp
    ps_grid.cpp
    ps_first_second_order.cpp
    cusp.cpp
    cpu_util.cpp
    gpu_util.cpp
    cuda_util.cpp
    scf_util.cpp
    symm.cpp
    gauss.cpp
    hess.cpp
    opt.cpp
    jellium_ints.cpp
    )

add_library(SlaterGPU STATIC "${INTEGRALS_SOURCES}")
target_include_directories(SlaterGPU PUBLIC
    "${CMAKE_SOURCE_DIR}/include"
    "${CMAKE_SOURCE_DIR}/src/libio/"
    )

add_definitions("-w")

set_property(TARGET SlaterGPU PROPERTY CXX_STANDARD 14)
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
