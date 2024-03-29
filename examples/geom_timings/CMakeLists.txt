add_executable(test_geom_timings
    main.cpp
    )
target_include_directories(test_geom_timings PRIVATE 
    "${CMAKE_SOURCE_DIR}/src/integrals"
    "${CMAKE_SOURCE_DIR}/src/libio"
    )
set(TEST_COMP_FLAGS ${CMAKE_CXX_FLAGS})
if(USE_ACC)
    target_link_options(test_geom_timings PRIVATE 
        -acc -ta=tesla
        )
    set(TEST_COMP_FLAGS ${TEST_COMP_FLAGS} 
        -DUSE_ACC=1 -acc -ta=tesla
        )
endif()

set_property(TARGET test_geom_timings PROPERTY CXX_STANDARD 11)
target_compile_options(test_geom_timings PRIVATE ${TEST_COMP_FLAGS} -O2)
target_link_libraries(test_geom_timings PRIVATE SlaterGPU io)
file(COPY geom_1 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY geom_2 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY geom_3 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY geom_4 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY geom_5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY run_geom_timings.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
