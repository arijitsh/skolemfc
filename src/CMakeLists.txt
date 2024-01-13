# Copyright (c) 2024, Arijit Shaw, Brendan Juba and Kuldeep S. Meel.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${GMP_INCLUDE_DIRS})
include_directories(${CRYPTOMINISAT5_INCLUDE_DIRS})


configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp" @ONLY)

add_library (skolemfc
    skolemfc.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp
)

add_executable (skolemfc-bin
    main.cpp
)

add_executable (skolemfc-test-bin
    test.cpp
)

target_link_libraries (skolemfc
    LINK_PUBLIC ${CRYPTOMINISAT5_LIBRARIES}
    ${GMPXX_LIB}
    ${GMP_LIB}
    ${ZLIB_LIBRARY}
)

target_link_libraries(skolemfc-bin
    ${Boost_LIBRARIES}
    skolemfc
)

target_link_libraries(skolemfc-test-bin
    ${Boost_LIBRARIES}
    skolemfc
)

set_target_properties(skolemfc PROPERTIES
    PUBLIC_HEADER "${skolemfc_public_headers}"
    VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
)

set_target_properties(skolemfc-bin PROPERTIES
    OUTPUT_NAME skolemfc
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

set_target_properties(skolemfc-test-bin PROPERTIES
    OUTPUT_NAME test
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
)

skolemfc_add_public_header(
    skolemfc ${CMAKE_CURRENT_SOURCE_DIR}/skolemfc.h )

# -----------------------------------------------------------------------------
# Copy public headers into build directory include directory.
# The skolemfcConfig.cmake we generate in the build directory depends on
# this.
# -----------------------------------------------------------------------------
set(HEADER_DEST "${PROJECT_BINARY_DIR}/include/skolemfc")
add_custom_target(CopyPublicHeaders ALL)
get_target_property(skolemfc_public_headers skolemfc PUBLIC_HEADER)
foreach(public_header ${skolemfc_public_headers})
    get_filename_component(HEADER_NAME ${public_header} NAME)
    add_custom_command(TARGET CopyPublicHeaders PRE_BUILD
                       COMMAND ${CMAKE_COMMAND} -E make_directory
                               "${HEADER_DEST}"
                       COMMAND ${CMAKE_COMMAND} -E echo
                       "Copying ${HEADER_NAME} to ${HEADER_DEST}"
                       COMMAND ${CMAKE_COMMAND} -E
                           copy_if_different
                           ${public_header}
                           "${HEADER_DEST}"
                      )
endforeach()


install(TARGETS skolemfc
    EXPORT ${SKOLEMFC_EXPORT_NAME}
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/skolemfc"
)

install(TARGETS skolemfc-bin
    EXPORT ${SKOLEMFC_EXPORT_NAME}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

