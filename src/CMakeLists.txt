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
    install(TARGETS cintw
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
