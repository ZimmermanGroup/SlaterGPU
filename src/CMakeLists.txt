add_subdirectory(integrals)
add_subdirectory(libio)

if(${DO_GTO})
    add_subdirectory(libcintw)
    message(STATUS "libcintw is enabled")
else()
    message(STATUS "libcintw is disabled")
endif()

# Install all library targets
install(TARGETS SlaterGPU io
        EXPORT SlaterGPUTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(${DO_GTO})
    set(_CINTW_INSTALL_TARGETS cintw)
    # When libcint is built from source via FetchContent (the default path in
    # libcintw/CMakeLists.txt), the "cint" CMake target exists and must be
    # included in the export set. CMake requires this because cintw links
    # against cint — exported targets must carry their dependencies.
    # When LIBCINT_PATH is set instead, cintw links to a file path (not a
    # CMake target), so there is no "cint" target and nothing extra to export.
    if(TARGET cint)
        list(APPEND _CINTW_INSTALL_TARGETS cint)
    endif()
    install(TARGETS ${_CINTW_INSTALL_TARGETS}
            EXPORT SlaterGPUTargets
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

# Export targets for use by other projects
install(EXPORT SlaterGPUTargets
        FILE SlaterGPUTargets.cmake
        NAMESPACE SlaterGPU::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SlaterGPU)
