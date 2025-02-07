set(CINT_SOURCES
    cint_funcs.h cint.h
    cintprep.cpp cintprep.h
    cintwrapper.cpp cintwrapper.h
    gto_grid.cpp gto_grid.h
    elements.h
    )

#find_package(BLAS REQUIRED)
#message("${BLAS_LIBRARIES}")

set(LIBCINT_PATH $ENV{LIBCINT_PATH})

add_library(cintw SHARED ${CINT_SOURCES})
set_property(TARGET cintw PROPERTY CXX_STANDARD 14)
set(CINTW_COMP_FLAGS ${CMAKE_CXX_FLAGS} 
    -O3 -g
    )

add_definitions("-w")

set(BLAS_INCLUDE_DIR $ENV{BLAS_INCLUDE_DIR})

target_compile_options(cintw PRIVATE ${CINTW_COMP_FLAGS})
target_include_directories(cintw PRIVATE
    "${LIBCINT_PATH}/include"
    "${CMAKE_SOURCE_DIR}/include/SlaterGPU"
    "{BLAS_INCLUDE_DIR}"
    "/usr/include/lapacke"
    "/usr/include/cblas"
    )
target_link_libraries(cintw PRIVATE
    "${LIBCINT_PATH}/lib64/libcint.so" 
    )
