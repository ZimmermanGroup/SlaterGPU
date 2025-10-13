add_executable(sgpu.exe
    main.cpp
    )
target_include_directories(sgpu.exe PRIVATE 
    "${CMAKE_SOURCE_DIR}/src/integrals"
    "${CMAKE_SOURCE_DIR}/src/libio"
    "${CMAKE_SOURCE_DIR}/src/libcintw"
    "${CMAKE_SOURCE_DIR}/include"
    )

set(TEST_COMP_FLAGS ${CMAKE_CXX_FLAGS})

# CUDA/GPU dependencies - only when USE_ACC is enabled
if(USE_ACC)
    enable_language(CUDA)
    find_package(CUDA REQUIRED)
    find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
    
    target_link_options(sgpu.exe PRIVATE 
        -acc=gpu
        )
    set(TEST_COMP_FLAGS ${TEST_COMP_FLAGS} 
        -DUSE_ACC=1 -acc=gpu
        )
    
    # GPU-specific link libraries
    set(GPU_LIBRARIES 
        "${CUDART_LIBRARY}" 
        "${CUDA_LIBRARIES}" 
        "${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/lib64/libcublas.so" 
        "${CUDA_TOOLKIT_ROOT_DIR}/../../math_libs/lib64/libcusolver.so"
        )
else()
    # CPU-only mode
    set(GPU_LIBRARIES "")
    message(STATUS "Building CPU-only version (USE_ACC=OFF)")
endif()

#LAPACK cmake not right.
#need to use this as a variable below
find_package(LAPACK)
find_package(BLAS)

set_property(TARGET sgpu.exe PROPERTY CXX_STANDARD 14)
target_compile_options(sgpu.exe PRIVATE ${TEST_COMP_FLAGS} -O2)

# Link libraries - GPU_LIBRARIES will be empty for CPU-only builds
target_link_libraries(sgpu.exe PRIVATE 
    SlaterGPU 
    io 
    cintw 
    "${FORTRANLIBS}" 
    ${GPU_LIBRARIES}
    "-L/home/paulzim/integrate/gpu/lib/lapack-3.9.0/lib/ -llapack -lblas" 
    "-L/export/apps/RockyOS8/nvidia_hpc/2024_249/Linux_x86_64/24.9/ -fortranlibs"
    )

file(COPY lih_VK1 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY lih_ri5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY h2_ri5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ch2_ri5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY beh2_ri5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY noplus_ri5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ch4_ri5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY run.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

install(TARGETS sgpu.exe
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})