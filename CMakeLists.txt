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

cmake_minimum_required(VERSION 3.8)

project(EPTlib VERSION 1.0.0 LANGUAGES CXX)

# Disallow in-source build
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(FATAL_ERROR
        "FATAL: Out-of-source build required. Create a separate directory for build files.")
endif()

# Configure files
string(TIMESTAMP BUILD_TIME "%Y-%m-%dT%H:%M:%SZ" UTC)

configure_file(${PROJECT_SOURCE_DIR}/include/eptlib/version.h.in
    ${PROJECT_BINARY_DIR}/include/eptlib/version.h
    @ONLY)

# Find packages
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)

# Enable testing with GoogleTest
include(GoogleTest)
find_package(GTest)
enable_testing()

# Set other properties common to all building
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# Set the library
set(EPTLIB_SRC
    src/ept_helmholtz.cc
    src/ept_interface.cc
    src/finite_difference.cc
    src/shape.cc
    src/version.cc)

add_library(EPTlib SHARED ${EPTLIB_SRC})

target_include_directories(EPTlib
    PUBLIC
        ${Boost_INCLUDE_DIRS}
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(EPTlib PRIVATE Eigen3::Eigen)

target_compile_features(EPTlib PUBLIC cxx_std_11)

# Test the library
if(GTest_FOUND)
    add_subdirectory(test)
endif()

# Set the Helmholtz-based EPT executable
add_executable(EPTlib_HELM src/exec/ept_helmholtz_main.cc)

target_include_directories(EPTlib_HELM
    PRIVATE
        extern/tinytoml/include)

target_link_libraries(EPTlib_HELM PRIVATE EPTlib)

# ----- Packaging -----

