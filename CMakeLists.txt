cmake_minimum_required(VERSION 3.15)

# Warning for potentially problematic CMake versions
if(CMAKE_VERSION VERSION_GREATER_EQUAL "4.1.0")
    message(WARNING 
        "CMake version 4.1.0 seemed to break CUDA library detection. "
        "If you encounter CUDART_LIBRARY errors, please downgrade to CMake <4.1.0. "
        "See: https://github.com/ZimmermanGroup/SlaterGPU/pull/50")
endif()

# Set the C and C++ compilers to NVIDIA HPC SDK compilers
# Ensure that the nvc and nvc++ compilers are available in the PATH
set(CMAKE_C_COMPILER nvc)
set(CMAKE_CXX_COMPILER nvc++)

project(SlaterGPU VERSION 0.1 
                  DESCRIPTION "Integrals over Slater type orbitals with GPU Acceleration"
                  LANGUAGES C CXX)

set(USE_ACC ON)
if("$ENV{USE_ACC}" STREQUAL "OFF")
    set(USE_ACC OFF)
endif()
message(STATUS "USE_ACC: ${USE_ACC}")
# NOTE: SlaterGPU contains no MPI code (no MPI_* calls, no <mpi.h>, no target
# links MPI::*). MPI is owned entirely by the ZEST/XCtera layers above, which
# find and link it themselves. So there is intentionally no find_package(MPI)
# or USE_MPI flag here — keeping the layer MPI-agnostic is what lets a
# USE_MPI=OFF stack configure without MPI present. See ZEST#92.
set(USE_OMP True)
set(RED_DOUBLE True)

if(NOT DO_GTO)
    set(DO_GTO True)
endif()

set(DO_GTO ${DO_GTO})

if(DO_GTO)
    # libcint's own headers (cint.h, cint_funcs.h) and library. LIBCINT_PATH
    # points at a manual install; otherwise the pixi environment is searched.
    set(LIBCINT_PATH $ENV{LIBCINT_PATH})
    find_path(LIBCINT_INCLUDE_DIR cint.h
        HINTS ${LIBCINT_PATH}/include
        REQUIRED)
    find_library(LIBCINT_LIBRARY cint
        HINTS ${LIBCINT_PATH}/lib ${LIBCINT_PATH}/lib64
        REQUIRED)
    message(STATUS "libcint headers: ${LIBCINT_INCLUDE_DIR}")
    message(STATUS "libcint library: ${LIBCINT_LIBRARY}")
endif()

# gen_4c_overlap / gen_4c_overlap_m4 need libcint built with WITH_4C1E=ON,
# which recipes/libcint provides. Set COMPILE_CINTW_4C=OFF in the environment
# when linking a libcint without 4c1e (e.g. a manual LIBCINT_PATH build);
# the functions then return zeros.
set(COMPILE_CINTW_4C ON)
if("$ENV{COMPILE_CINTW_4C}" STREQUAL "OFF")
    set(COMPILE_CINTW_4C OFF)
endif()
message(STATUS "COMPILE_CINTW_4C: ${COMPILE_CINTW_4C}")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# use this flag for verbose OpenACC compilation output to troubleshoot memory issues
# set(CMAKE_CXX_FLAGS " -Minfo=accel ")

set(CMAKE_CXX_FLAGS_RELEASE "-O2")

if (USE_OMP)
    find_package(OpenMP  REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
if (USE_ACC)
    find_package(OpenACC REQUIRED)
    # Provides CUDA::cublas / CUDA::cusolver imported targets (headers + libs),
    # so the GPU build does not depend on CPATH from the developer's shell
    # reaching the compiler (see issue #98)
    find_package(CUDAToolkit REQUIRED)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Install directories
include(GNUInstallDirs)
# Force lib directory instead of lib64 to match conda expectations
set(CMAKE_INSTALL_LIBDIR lib)

add_subdirectory(src)
add_subdirectory(examples)

# Install header files
install(DIRECTORY include/ 
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.h")

add_definitions("-w")
message(STATUS "${CMAKE_CXX_FLAGS}")
message(STATUS "DO_GTO: ${DO_GTO}")

# Generate package configuration files
include(CMakePackageConfigHelpers)

# Configure the config file
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/SlaterGPUConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/SlaterGPUConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SlaterGPU
    PATH_VARS CMAKE_INSTALL_INCLUDEDIR
)

# Configure the version file
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/SlaterGPUConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# Install the config files
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/SlaterGPUConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/SlaterGPUConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SlaterGPU
)
