option(DOXYGEN_ENABLED "Enable Doxygen documentation generation" ON)

find_package(Doxygen)

if(DOXYGEN_FOUND AND DOXYGEN_ENABLED)
    add_custom_target(doc_doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )
endif()