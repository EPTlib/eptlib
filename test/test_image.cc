/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2023  Alessandro Arduino
*  Istituto Nazionale di Ricerca Metrologica (INRiM)
*  Strada delle cacce 91, 10135 Torino
*  ITALY
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in all
*  copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*  SOFTWARE.
*
*****************************************************************************/

#include "gtest/gtest.h"

#include "eptlib/image.h"

TEST(ImageGTest,GetSize) {
    eptlib::Image<double> img(15,10,5);
    ASSERT_EQ(img.GetSize()[0], 15);
    ASSERT_EQ(img.GetSize()[1], 10);
    ASSERT_EQ(img.GetSize()[2], 5);
    ASSERT_EQ(img.GetSize(0), 15);
    ASSERT_EQ(img.GetSize(1), 10);
    ASSERT_EQ(img.GetSize(2), 5);
}

TEST(ImageGTest,GetNVox) {
    eptlib::Image<double> img(15,10,5);
    ASSERT_EQ(img.GetNVox(), 750);
}

TEST(ImageGTest,GetData) {
    eptlib::Image<double> img(15,10,5);
    std::iota(img.GetData().begin(),img.GetData().end(),1);
    for (int i = 0; i<img.GetNVox(); ++i) {
        ASSERT_EQ(img(i), i+1);
        ASSERT_EQ(img.At(i), i+1);
    }
    for (int k = 0; k<img.GetSize(2); ++k) {
        for (int j = 0; j<img.GetSize(1); ++j) {
            for (int i = 0; i<img.GetSize(0); ++i) {
                ASSERT_EQ(img(i,j,k), eptlib::IJKToIdx(i,j,k, 15,10)+1);
                ASSERT_EQ(img.At(i,j,k), eptlib::IJKToIdx(i,j,k, 15,10)+1);
            }
        }
    }
}

TEST(ImageGTest,AccessData) {
    eptlib::Image<double> img(15,10,5);
    for (int i = 0; i<img.GetNVox(); ++i) {
        img(i) = i+1;
    }
    for (int k = 0; k<img.GetSize(2); ++k) {
        for (int j = 0; j<img.GetSize(1); ++j) {
            for (int i = 0; i<img.GetSize(0); ++i) {
                ASSERT_EQ(img(i,j,k), eptlib::IJKToIdx(i,j,k, 15,10)+1);
                ASSERT_EQ(img.At(i,j,k), eptlib::IJKToIdx(i,j,k, 15,10)+1);
            }
        }
    }
}
