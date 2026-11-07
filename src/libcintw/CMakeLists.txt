set(CINT_SOURCES
    cintprep.cpp
    cintwrapper.cpp
    gto_grid.cpp
    )

find_package(BLAS REQUIRED)
message(STATUS "${BLAS_LIBRARIES}")

add_library(cintw SHARED ${CINT_SOURCES})
set_property(TARGET cintw PROPERTY CXX_STANDARD 14)
set(CINTW_COMP_FLAGS ${CMAKE_CXX_FLAGS}
    -O2 -g
    )

add_definitions("-w")

target_compile_options(cintw PRIVATE ${CINTW_COMP_FLAGS})
if(COMPILE_CINTW_4C)
    target_compile_definitions(cintw PRIVATE COMPILE_CINTW_4C)
endif()
target_include_directories(cintw PRIVATE
    "${CMAKE_SOURCE_DIR}/include"
    "${LIBCINT_INCLUDE_DIR}"
    "$ENV{BLAS_INCLUDE_DIR}"
    "$ENV{PREFIX}/include" # include pixi dependencies
    )

# LIBCINT_LIBRARY is found in the top-level CMakeLists.txt
# (LIBCINT_PATH for a manual build, otherwise the pixi environment)
target_link_libraries(cintw ${LIBCINT_LIBRARY})
