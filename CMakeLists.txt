#=============================================================================
#
#     Program: EPTlib
#     Author: Alessandro Arduino <a.arduino@inrim.it>
#
#  MIT License
#
#  Copyright (c) 2020-2024  Alessandro Arduino
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

cmake_minimum_required(VERSION 3.16)

project(EPTlib VERSION 0.4.2 LANGUAGES C CXX)

# Define building options
option(build_docs "Build documentation (requires doxygen)" ON)
option(build_docs_latex "Build documentation with LaTeX output (requires latex)" ON)

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

# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Configure files
string(TIMESTAMP BUILD_TIME "%Y-%m-%dT%H:%M:%SZ" UTC)
configure_file(${PROJECT_SOURCE_DIR}/include/eptlib/version.h.in
    ${PROJECT_BINARY_DIR}/include/eptlib/version.h
    @ONLY)

# Find packages
include(VerboseFindPackage)
# hdf5
if(WIN32)
    set(HDF5_USE_STATIC_LIBRARIES ON)
endif()
VerboseFindPackage(HDF5 COMPONENTS C CXX REQUIRED)
include_directories("${HDF5_INCLUDE_DIRS}")
# boost
VerboseFindPackage(Boost REQUIRED)
include_directories("${Boost_INCLUDE_DIRS}")
# eigen3
VerboseFindPackage(Eigen3 REQUIRED)

# openmp
VerboseFindPackage(OpenMP)

# disable fast floating-point operations in Intel compiler
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel" OR CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fp:precise")
endif()

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

# Build the documentation
if(build_docs)
    VerboseFindPackage(Doxygen OPTIONAL_COMPONENTS dot)
    if(DOXYGEN_FOUND)
        # set the source files
        set(DOCS_SRC
            include/eptlib/differential_operator.h
            include/eptlib/ept_convreact.h
            include/eptlib/ept_gradient.h
            include/eptlib/ept_helmholtz_chi2.h
            include/eptlib/ept_helmholtz.h
            include/eptlib/ept_interface.h
            include/eptlib/ept_methods.h
            include/eptlib/image.h
            include/eptlib/shape.h
            include/eptlib/util.h
            include/eptlib/filter/anatomical_savitzky_golay.h
            include/eptlib/filter/moving_window.h
            include/eptlib/filter/postprocessing.h
            include/eptlib/filter/savitzky_golay.h
            include/eptlib/filter/weight_functions.h
            include/eptlib/io/io_hdf5.h
            include/eptlib/io/io_toml.h
            include/eptlib/io/io_util.h
            include/eptlib/linalg/matrix.h
            include/eptlib/linalg/qr.h
            include/eptlib/linalg/regression.h
            include/eptlib/linalg/vector.h
            include/eptlib/polynomial/fitting.h)
        # set custom doxygen options
        set(DOXYGEN_OUTPUT_DIRECTORY docs)
        set(DOXYGEN_PROJECT_BRIEF "An open-source, extensible C++ library of electric properties tomography methods.")
        set(DOXYGEN_PROJECT_LOGO "img/EPTlib_lib512.png")
        set(DOXYGEN_HTML_HEADER "docs/header.html")
        set(DOXYGEN_HTML_EXTRA_STYLESHEET "docs/style.css")
        set(DOXYGEN_EXTRACT_PRIVATE YES)
        set(DOXYGEN_EXTRACT_LOCAL_CLASSES NO)
        set(DOXYGEN_HIDE_UNDOC_CLASSES YES)
        if(build_docs_latex)
            set(DOXYGEN_GENERATE_LATEX YES)
        endif(build_docs_latex)
        # create the doxygen command
        doxygen_add_docs(docs ${DOCS_SRC} ALL #USE_STAMP_FILE
            COMMENT "Generating documentation with Doxygen")
    endif(DOXYGEN_FOUND)
endif(build_docs)

# ----- Installing -----
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    include(InstallRequiredSystemLibraries)
endif()

set(CMAKE_INSTALL_MFC_LIBRARIES ON)

install(TARGETS EPTlib
    EXPORT EPTlibTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include)

install(FILES ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}
    DESTINATION bin)

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
set(CPACK_PACKAGE_VENDOR "Alessandro Arduino")
set(CPACK_PACKAGE_CONTACT "Alessandro Arduino")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_PACKAGE_DESCRIPTION "An open source, extensible collection of EPT methods. Visit https://eptlib.github.io for more details.")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://eptlib.github.io")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.md")

set(CPACK_WIX_UPGRADE_GUID "02ED789A-5006-44F2-AB33-C07C00DBEF6B")
configure_file(${PROJECT_SOURCE_DIR}/LICENSE ${PROJECT_BINARY_DIR}/LICENSE.txt COPYONLY)
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_BINARY_DIR}/LICENSE.txt")
set(CPACK_WIX_PRODUCT_ICON "${PROJECT_SOURCE_DIR}/img/EPTlib_ico.ico")
set(CPACK_WIX_UI_BANNER "${PROJECT_SOURCE_DIR}/img/wix-banner.png")
set(CPACK_WIX_UI_DIALOG "${PROJECT_SOURCE_DIR}/img/wix-dialog.png")
set(CPACK_WIX_CMAKE_PACKAGE_REGISTRY "EPTlib")

set(CPACK_DEBIAN_FILE_NAME "DEB-DEFAULT")
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
set(CPACK_SOURCE_IGNORE_FILES
    /.git
    /.vscode
    /.vs
    /*out*
    /*build*
    /*install*
)

if(WIN32)
    set(CPACK_GENERATOR "WIX")
endif()

include(CPack)
