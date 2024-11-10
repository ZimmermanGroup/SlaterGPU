add_executable(sgpu.exe
    main.cpp
    )
target_include_directories(sgpu.exe PRIVATE 
    "${CMAKE_SOURCE_DIR}/src/integrals"
    "${CMAKE_SOURCE_DIR}/src/libio"
    "${CMAKE_SOURCE_DIR}/include"
    )
set(TEST_COMP_FLAGS ${CMAKE_CXX_FLAGS})
if(USE_ACC)
    target_link_options(sgpu.exe PRIVATE 
        -acc -ta=tesla
        )
    set(TEST_COMP_FLAGS ${TEST_COMP_FLAGS} 
        -DUSE_ACC=1 -acc -ta=tesla
        )
endif()

#LAPACK cmake not right.
#need to use this as a variable below
find_package(LAPACK)

find_package(BLAS)
enable_language(CUDA)
find_package(CUDA REQUIRED)
find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})

#enable_language(FORTRAN)
#find_library(FORTRANLIBS fortranlibs ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})

set_property(TARGET sgpu.exe PROPERTY CXX_STANDARD 14)
target_compile_options(sgpu.exe PRIVATE ${TEST_COMP_FLAGS} -O2)
#target_link_libraries(sgpu.exe PRIVATE SlaterGPU io "${CUDART_LIBRARY}" "${CUDA_LIBRARIES}" )
target_link_libraries(sgpu.exe PRIVATE SlaterGPU io "${FORTRANLIBS}" "${CUDART_LIBRARY}" "${CUDA_LIBRARIES}" "${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/lib64/libcublas.so" "${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/lib64/libcusolver.so" "-L/home/paulzim/integrate/gpu/lib/lapack-3.9.0/lib/ -llapack -lblas" "-L/export/apps/RockyOS8/nvidia_hpc/2024_249/Linux_x86_64/24.9/ -fortranlibs" )
file(COPY geom_1 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY geom_2 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY geom_3 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY geom_4 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY geom_5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY geom_ps DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY run.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
