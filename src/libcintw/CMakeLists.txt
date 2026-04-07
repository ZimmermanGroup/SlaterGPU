set(CINT_SOURCES
    cintprep.cpp
    cintwrapper.cpp
    gto_grid.cpp
    )

find_package(BLAS REQUIRED)
message(STATUS "${BLAS_LIBRARIES}")

set(LIBCINT_PATH $ENV{LIBCINT_PATH})

# ---------------------------------------------------------------------------
# Resolve libcint dependency
#
# libcint (https://github.com/sunqm/libcint) provides the low-level integral
# engine for GTO-based calculations. We need it built with -DWITH_4C1E=ON to
# enable 4-center 1-electron overlap integrals (∫φ_iφ_jφ_kφ_l dr) used for
# DFT inversion. The conda-forge package does not include these, so we build
# from source by default.
#
# Two modes are supported:
#   1. FetchContent (default) — automatically downloads and builds libcint
#      from source as part of the SlaterGPU build. This is the portable path:
#      it works identically on any machine via `pixi install`.
#   2. LIBCINT_PATH (manual) — set the LIBCINT_PATH environment variable to
#      point to a pre-built libcint installation. This skips FetchContent and
#      links directly against the specified library. Useful for development
#      or when testing a patched libcint build.
#
# See: https://github.com/ZimmermanGroup/SlaterGPU/issues/82
# ---------------------------------------------------------------------------
if(LIBCINT_PATH)
    message(STATUS "Using manual libcint installation from: ${LIBCINT_PATH}")
    set(LIBCINT_INCLUDE_DIR "${LIBCINT_PATH}/include")
    set(LIBCINT_LINK_LIB "${LIBCINT_PATH}/lib64/libcint.so.5.3.0")
else()
    message(STATUS "Building libcint from source with WITH_4C1E=ON")
    include(FetchContent)

    # WITH_4C1E enables cint4c1e_cart/cint4c1e_sph functions in libcint.
    # This must be set before FetchContent_MakeAvailable so it's visible
    # to libcint's CMakeLists.txt during configure.
    set(WITH_4C1E ON CACHE BOOL "Enable 4-center 1-electron integrals" FORCE)

    FetchContent_Declare(
        libcint
        GIT_REPOSITORY https://github.com/sunqm/libcint.git
        GIT_TAG v5.3.0
    )
    FetchContent_MakeAvailable(libcint)

    # FetchContent populates libcint_SOURCE_DIR with the cloned source tree.
    # The headers we need (cint.h, etc.) live in its include/ directory.
    set(LIBCINT_INCLUDE_DIR "${libcint_SOURCE_DIR}/include")
    # Link against the CMake target so we get proper build-order dependencies.
    # Note: this target must also be added to the install export set in
    # src/CMakeLists.txt — see the TARGET cint guard there.
    set(LIBCINT_LINK_LIB cint)
endif()

# ---------------------------------------------------------------------------
# Build the cintw (libcint wrapper) shared library
# ---------------------------------------------------------------------------
add_library(cintw SHARED ${CINT_SOURCES})
set_property(TARGET cintw PROPERTY CXX_STANDARD 14)
set(CINTW_COMP_FLAGS ${CMAKE_CXX_FLAGS}
    -O2 -g
    )

add_definitions("-w")

target_compile_options(cintw PRIVATE ${CINTW_COMP_FLAGS})
target_include_directories(cintw PRIVATE
    "${CMAKE_SOURCE_DIR}/include"
    "${LIBCINT_INCLUDE_DIR}"
    "$ENV{BLAS_INCLUDE_DIR}"
    "$ENV{PREFIX}/include"
    )
target_link_libraries(cintw PRIVATE ${LIBCINT_LINK_LIB})
