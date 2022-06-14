add_executable(test_gto_grid
    main.cpp
    )
target_include_directories(test_gto_grid PRIVATE 
    "${CMAKE_SOURCE_DIR}/src/integrals"
    "${CMAKE_SOURCE_DIR}/src/libio"
    "${CMAKE_SOURCE_DIR}/src/libcintw"
    )
set(TEST_COMP_FLAGS ${CMAKE_CXX_FLAGS})
if(USE_ACC)
    target_link_options(test_gto_grid PRIVATE 
        -acc -ta=tesla
        )
    set(TEST_COMP_FLAGS ${TEST_COMP_FLAGS} 
        -DUSE_ACC=1 -acc -ta=tesla
        )
endif()

set_property(TARGET test_gto_grid PROPERTY CXX_STANDARD 11)
target_compile_options(test_gto_grid PRIVATE ${TEST_COMP_FLAGS} -O2)
target_link_libraries(test_gto_grid PRIVATE SlaterGPU io cintw)
file(COPY tests DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
