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
    becke2.cpp
    sphericald.cpp
    pVpd.cpp
    integrals_aux.cpp
    )

add_library(SlaterGPU STATIC "${INTEGRALS_SOURCES}")
target_include_directories(SlaterGPU PUBLIC
    "${CMAKE_SOURCE_DIR}/include/SlaterGPU"
    "${CMAKE_SOURCE_DIR}/include"
    "${CMAKE_SOURCE_DIR}/src/libcintw/"
    "${CMAKE_SOURCE_DIR}/src/libio/"
    )

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
