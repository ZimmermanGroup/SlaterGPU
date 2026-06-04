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
if(COMPILE_CINTW_4C)
    target_compile_definitions(cintw PRIVATE COMPILE_CINTW_4C)
endif()
# Put the linked libcint's headers ahead of SlaterGPU/include so that the
# vendored cint.h does not shadow the cint.h that ships with the libcint we
# actually link against. A version/ABI skew between the two (e.g. I8 vs int
# FINT, or BAS_SLOTS drift) would otherwise corrupt integrals silently.
target_include_directories(cintw PRIVATE
    "${LIBCINT_PATH}/include"
    "$ENV{PREFIX}/include" # include pixi dependencies
    "${CMAKE_SOURCE_DIR}/include"
    "$ENV{BLAS_INCLUDE_DIR}"
    )

# Link against custom libcint if LIBCINT_PATH is set (manual build)
# Otherwise pixi-provided libcint is used
if(LIBCINT_PATH)
    message(STATUS "Using manual libcint installation from: ${LIBCINT_PATH}")
    target_link_libraries(cintw PRIVATE
        "${LIBCINT_PATH}/lib64/libcint.so.5.3.0"
        )
else()
    message(STATUS "Using pixi-provided libcint")
    find_library(LIBCINT_LIBRARY cint PATHS ${CMAKE_PREFIX_PATH}/lib)
    target_link_libraries(cintw ${LIBCINT_LIBRARY})

endif()