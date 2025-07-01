cmake_minimum_required(VERSION 3.15)

# Set the C and C++ compilers to NVIDIA HPC SDK compilers
# Ensure that the nvc and nvc++ compilers are available in the PATH
set(CMAKE_C_COMPILER nvc)
set(CMAKE_CXX_COMPILER nvc++)

project(SlaterGPU VERSION 0.1 
                  DESCRIPTION "Integrals over Slater type orbitals with GPU Acceleration"
                  LANGUAGES C CXX)

set(USE_ACC True)
set(USE_MPI True)
set(USE_OMP True)
set(RED_DOUBLE True)

if(NOT DO_GTO)
    set(DO_GTO True)
endif()

set(DO_GTO ${DO_GTO})

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# use this flag for verbose OpenACC compilation output to troubleshoot memory issues
# set(CMAKE_CXX_FLAGS " -Minfo=accel ")

set(CMAKE_CXX_FLAGS_RELEASE "-O2")

if(USE_MPI)
    find_package(MPI     REQUIRED)
endif()
if (USE_OMP)
    find_package(OpenMP  REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
if (USE_ACC)
    find_package(OpenACC REQUIRED)
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
