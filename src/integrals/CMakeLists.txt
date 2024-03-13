set(INTEGRALS_SOURCES
    integrals.cpp integrals_d.cpp integrals.h
    braggslater.h 
    gamma.h gamma.cpp
    lebedev2.cpp lebedev2.h
    murak.cpp murak.h
    pVp.cpp pVp.h
    reduce.cpp reduce.h
    spherical.cpp spherical.h
    Vinr.cpp Vinr.h
    grid_util.cpp grid_util.h
    opt.h opt.cpp
    utils.cpp utils.h
    )

add_library(SlaterGPU STATIC "${INTEGRALS_SOURCES}")

find_package(BLAS REQUIRED)
message("BLAS: ${BLAS_LIBRARIES}")

target_include_directories(SlaterGPU PUBLIC
    "${CMAKE_SOURCE_DIR}/include/SlaterGPU"
    "${BLAS_HOME}/include/openblas"
    )

set_property(TARGET SlaterGPU PROPERTY CXX_STANDARD 11)
set_property(TARGET SlaterGPU PROPERTY POSITION_INDEPENDENT_CODE ON)

target_link_options(SlaterGPU PUBLIC -lgfortran) 

target_link_libraries(SlaterGPU PRIVATE
    "${BLAS_LIBRARIES}"
)

set(SLATER_COMP_FLAGS ${CMAKE_CXX_FLAGS} -DRED_DOUBLE=1 -O2
    -g -D__GCC_ATOMIC_TEST_AND_SET_TRUEVAL=1 -DGAMMA=1
    )
if(USE_ACC)
    target_compile_options(SlaterGPU PRIVATE
        ${SLATER_COMP_FLAGS} -DUSE_ACC=1 -acc=gpu -ta=tesla
        )
else()
    target_compile_options(SlaterGPU PRIVATE
        ${SLATER_COMP_FLAGS} -DUSE_ACC=0
        )
endif()
