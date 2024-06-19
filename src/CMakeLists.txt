add_subdirectory(integrals)
add_subdirectory(libio)

if(DO_GTO)
    add_subdirectory(libcintw)
    message(STATUS "libcintw is enabled")
else()
    message(STATUS "libcintw is disabled")
endif()
