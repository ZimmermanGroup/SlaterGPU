set(IO_SOURCES
    read.cpp
    write.cpp
    print.cpp
    )

add_library(io SHARED "${IO_SOURCES}")
set_property(TARGET io PROPERTY CXX_STANDARD 14)

target_include_directories(io PRIVATE 
    "${CMAKE_SOURCE_DIR}/src/integrals"
)

target_link_libraries(io PRIVATE SlaterGPU)

set(IO_COMP_FLAGS ${CMAKE_CXX_FLAGS} 
    -O2 -Minfo=accel -g -D__GCC_ATOMIC_TEST_AND_SET_TRUEVAL=1
    )

if(USE_ACC)
    target_compile_options(io PRIVATE
        ${IO_COMP_FLAGS} -DUSE_ACC=1 -acc=gpu -Minfo=accel
        )
    target_link_options(io PRIVATE
        -acc=gpu
        )
else()
    target_compile_options(io PRIVATE
        ${IO_COMP_FLAGS} -DUSE_ACC=0 -Minfo=accel
        )
endif()
