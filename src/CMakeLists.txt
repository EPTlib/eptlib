#=============================================================================
#
#     Program: EPTlib
#     Author: Alessandro Arduino <a.arduino@inrim.it>
#
#  MIT License
#
#  Copyright (c) 2020-2023  Alessandro Arduino
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

set(EPTLIB_SRC
    ept_convreact.cc
    ept_gradient.cc
    ept_helmholtz_chi2.cc
    ept_helmholtz.cc
    ept_interface.cc
    shape.cc
    filter/savitzky_golay.cc
    io/io_hdf5.cc
    io/io_toml.cc
    io/io_util.cc
    linalg/linalg_householder.cc
    linalg/linalg_qr.cc
    linalg/linalg_util.cc
)

add_library(EPTlib SHARED ${EPTLIB_SRC})

target_include_directories(EPTlib
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/tinytoml/include>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:include>
        $<INSTALL_INTERFACE:share/eptlib/extern/tinytoml/include>
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(EPTlib PRIVATE Eigen3::Eigen)

if(HDF5_LIBRARIES)
    target_link_libraries(EPTlib PUBLIC ${HDF5_LIBRARIES})
else()
    target_link_libraries(EPTlib PUBLIC hdf5 hdf5_cpp)
endif()

if(OpenMP_FOUND)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp:llvm")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -openmp:llvm")
    else()
        target_link_libraries(EPTlib PUBLIC OpenMP::OpenMP_CXX)
    endif()
endif()

target_compile_features(EPTlib PUBLIC cxx_std_17)

set_property(TARGET EPTlib
    PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
