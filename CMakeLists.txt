#=============================================================================
#
#     Program: EPTlib
#     Author: Alessandro Arduino <a.arduino@inrim.it>
#
#  MIT License
#
#  Copyright (c) 2020  Alessandro Arduino
#  Istituto Nazionale di Ricerca Metrologica (INRiM)
#  Strada delle cacce 91, 10135 Torino
#  ITALY
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  SOFTWARE.
#
#=============================================================================

cmake_minimum_required(VERSION 3.12)

project(EPTlib VERSION 0.1.0 LANGUAGES C CXX)

# Disallow in-source build
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(FATAL_ERROR
        "FATAL: Out-of-source build required. Create a separate directory for build files.")
endif()

# Set some CMAKE_ variables
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# Configure files
string(TIMESTAMP BUILD_TIME "%Y-%m-%dT%H:%M:%SZ" UTC)
configure_file(${PROJECT_SOURCE_DIR}/include/eptlib/version.h.in
    ${PROJECT_BINARY_DIR}/include/eptlib/version.h
    @ONLY)

# Find packages
include(VerboseFindPackage)
# hdf5
set(HDF5_USE_STATIC_LIBRARIES ON)
VerboseFindPackage(HDF5 COMPONENTS C CXX REQUIRED)
include_directories("${HDF5_INCLUDE_DIRS}")
# boost
VerboseFindPackage(Boost REQUIRED)
include_directories("${Boost_INCLUDE_DIRS}")
# eigen3
VerboseFindPackage(Eigen3 REQUIRED)

# Enable testing with GoogleTest
include(GoogleTest)
VerboseFindPackage(GTest)
enable_testing()

# Set the library
add_subdirectory(src)

# Set the applications
set(EPTLIB_APP_NAME "EPTlib_app")
add_subdirectory(app)

# Test everything
if(GTest_FOUND)
    add_subdirectory(test)
endif()

# ----- Installing -----

install(TARGETS EPTlib
    EXPORT EPTlibTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include)

install(EXPORT EPTlibTargets
    FILE EPTlibTargets.cmake
    NAMESPACE EPTlib::
    DESTINATION lib/cmake/eptlib)

install(DIRECTORY
        ${PROJECT_BINARY_DIR}/include/
        ${PROJECT_SOURCE_DIR}/include/
    DESTINATION
        include
    PATTERN "*.in" EXCLUDE)
    
install(DIRECTORY
        ${PROJECT_SOURCE_DIR}/extern/
    DESTINATION
        share/eptlib/extern)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    EPTlibConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion)

configure_file(cmake/EPTlibConfig.cmake.in EPTlibConfig.cmake @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/EPTlibConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/EPTlibConfigVersion.cmake"
    DESTINATION lib/cmake/eptlib)

install(TARGETS ${EPTLIB_APP_NAME} RUNTIME DESTINATION bin)

# ----- Packaging -----

set(CPACK_PACKAGE_NAME "EPTlib")
set(CPACK_PACKAGE_CONTACT "Alessandro Arduino <a.arduino@inrim.it>")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "An open source, extensible collection of EPT methods. Visit https://eptlib.github.io for more details.")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.md")

set(CPACK_GENERATOR "WIX")
set(CPACK_WIX_UPGRADE_GUID "02ED789A-5006-44F2-AB33-C07C00DBEF6B")
configure_file(${PROJECT_SOURCE_DIR}/LICENSE ${PROJECT_BINARY_DIR}/LICENSE.txt COPYONLY)
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_BINARY_DIR}/LICENSE.txt")

set(CPACK_DEBIAN_FILE_NAME "DEB-DEFAULT")
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libhdf5-cpp-100")

set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
set(CPACK_SOURCE_IGNORE_FILES
    /.git
    /.vscode
    /.vs
    /*out*
    /*build*
    /*install*
)

include(CPack)
