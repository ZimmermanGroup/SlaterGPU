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
        -acc=gpu
        )
    set(TEST_COMP_FLAGS ${TEST_COMP_FLAGS} 
        -DUSE_ACC=1 -acc=gpu
        )
endif()

set_property(TARGET sgpu.exe PROPERTY CXX_STANDARD 14)
target_compile_options(sgpu.exe PRIVATE ${TEST_COMP_FLAGS} -O2)
target_link_libraries(sgpu.exe PRIVATE SlaterGPU io)
file(COPY geom_1 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY geom_2 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY geom_3 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY geom_4 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY geom_5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY geom_ps DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY run.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
