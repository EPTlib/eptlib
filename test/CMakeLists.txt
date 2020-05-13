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

set(EPTlib_TEST_SRC
    main.cc
    test_ept_convreact.cc
    test_ept_helmholtz.cc
    test_ept_interface.cc
    test_finite_difference.cc
    test_io_hdf5.cc
    test_io_toml.cc
    test_median_filter.cc
    test_shape.cc)

add_executable(EPTlib_test ${EPTlib_TEST_SRC})

target_link_libraries(EPTlib_test PRIVATE GTest::GTest GTest::Main EPTlib)

gtest_add_tests(TARGET EPTlib_test
    WORKING_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input
    DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test)
