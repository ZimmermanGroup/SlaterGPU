set(REL_SOURCES
    rel.cpp rel.h
    )

add_library(rel SHARED "${REL_SOURCES}")
set_property(TARGET rel PROPERTY CXX_STANDARD 11)

target_include_directories(rel PRIVATE 
    "${CMAKE_SOURCE_DIR}/src/integrals"
    "${CMAKE_SOURCE_DIR}/src/libio"
)

target_link_libraries(rel PRIVATE SlaterGPU)

