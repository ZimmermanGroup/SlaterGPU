set(IO_SOURCES
    read.cpp read.h
    write.cpp write.h
    )

add_library(io SHARED "${IO_SOURCES}")
set_property(TARGET io PROPERTY CXX_STANDARD 11)

target_include_directories(io PRIVATE 
    "${CMAKE_SOURCE_DIR}/src/integrals"
)

target_link_libraries(io PRIVATE SlaterGPU)

set(IO_COMP_FLAGS ${CMAKE_CXX_FLAGS} 
    -O2 -g -D__GCC_ATOMIC_TEST_AND_SET_TRUEVAL=1
    )

if(USE_ACC)
    target_compile_options(io PRIVATE
        ${IO_COMP_FLAGS} -DUSE_ACC=1 -acc=gpu
        )
else()
    target_compile_options(io PRIVATE
        ${IO_COMP_FLAGS} -DUSE_ACC=0
        )
endif()
