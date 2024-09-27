cmake_minimum_required(VERSION 3.15)
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

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

add_subdirectory(src)
add_subdirectory(examples)

add_definitions("-w")
message(STATUS "${CMAKE_CXX_FLAGS}")
message(STATUS "DO_GTO: ${DO_GTO}")
