add_executable(test_rel_ham
    main.cpp
    )

target_include_directories(test_rel_ham PRIVATE 
    "${CMAKE_SOURCE_DIR}/src/librel"
    "${CMAKE_SOURCE_DIR}/src/integrals"
    "${CMAKE_SOURCE_DIR}/src/libio"
    )

set(TEST_COMP_FLAGS ${CMAKE_CXX_FLAGS})

if(USE_ACC)
    target_link_options(test_rel_ham PRIVATE 
        -acc -ta=tesla
        )
    set(TEST_COMP_FLAGS ${TEST_COMP_FLAGS} 
        -DUSE_ACC=1 -acc -ta=tesla
        )
endif()

set_property(TARGET test_rel_ham PROPERTY CXX_STANDARD 11)
target_compile_options(test_rel_ham PRIVATE ${TEST_COMP_FLAGS} -O2)
target_link_libraries(test_rel_ham PRIVATE SlaterGPU rel io)
file(COPY rel_ham.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
