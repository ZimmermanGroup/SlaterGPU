set(CINT_SOURCES
    cintprep.cpp
    cintwrapper.cpp
    gto_grid.cpp
    )

find_package(BLAS REQUIRED)
message(STATUS "${BLAS_LIBRARIES}")

set(LIBCINT_PATH $ENV{LIBCINT_PATH})

add_library(cintw SHARED ${CINT_SOURCES})
set_property(TARGET cintw PROPERTY CXX_STANDARD 14)
set(CINTW_COMP_FLAGS ${CMAKE_CXX_FLAGS}
    -O2 -g
    )

add_definitions("-w")

target_compile_options(cintw PRIVATE ${CINTW_COMP_FLAGS})
target_include_directories(cintw PRIVATE
    "${CMAKE_SOURCE_DIR}/include"
    "${LIBCINT_PATH}/include"
    "$ENV{BLAS_INCLUDE_DIR}"
    "$ENV{PREFIX}/include" # include pixi dependencies
    )

# Uncomment this block if you want to link against a manually built libcint
# outside of pixi
target_link_libraries(cintw PRIVATE
    "${LIBCINT_PATH}/lib64/libcint.so.5.3.0" 
    )
