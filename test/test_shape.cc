/*****************************************************************************
*
*     Program: EPTlib
*     Author: Alessandro Arduino <a.arduino@inrim.it>
*
*  MIT License
*
*  Copyright (c) 2020-2023  Alessandro Arduino
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

#include "eptlib/shape.h"

TEST(ShapeGTest,GetSize) {
    eptlib::Shape shape(15,10,5);
    ASSERT_EQ(shape.GetSize()[0], 15);
    ASSERT_EQ(shape.GetSize()[1], 10);
    ASSERT_EQ(shape.GetSize()[2], 5);
    ASSERT_EQ(shape.GetSize(0), 15);
    ASSERT_EQ(shape.GetSize(1), 10);
    ASSERT_EQ(shape.GetSize(2), 5);
}

TEST(ShapeGTest,GetNVox) {
    eptlib::Shape shape(15,10,5);
    ASSERT_EQ(shape.GetNVox(), 750);
}

TEST(ShapeGTest,GetVolume) {
    eptlib::Shape shape(15,10,5);
    for (int idx = 0; idx<13; ++idx) {
        shape.GetData().set(idx);
    }
    ASSERT_EQ(shape.GetVolume(), 13);
}

TEST(ShapeGTest,GetData) {
    /*
     * shape = [
     *     1 1
     *     1 0
     *     0 1
     *     0 0
     * ];
     */
    eptlib::Shape shape(2,4,1);
    auto& data = shape.GetData();
    data.set(0);
    data.set(1);
    data.set(2);
    data.set(5);
    for (int idx = 0; idx<shape.GetNVox(); ++idx) {
        ASSERT_EQ(shape(idx), idx==0||idx==1||idx==2||idx==5);
    }
    for (int k = 0; k<shape.GetSize(2); ++k) {
        for (int j = 0; j<shape.GetSize(1); ++j) {
            for (int i = 0; i<shape.GetSize(0); ++i) {
                int idx = eptlib::IJKToIdx(i,j,k, 2,4);
                ASSERT_EQ(shape(i,j,k), idx==0||idx==1||idx==2||idx==5);
            }
        }
    }
}

TEST(ShapeGTest,AccessData) {
    eptlib::Shape shape(2,4,1);
    for (int idx = 0; idx<shape.GetNVox(); ++idx) {
        shape(idx) = idx==0||idx==1||idx==2||idx==5;
    }
    for (int k = 0; k<shape.GetSize(2); ++k) {
        for (int j = 0; j<shape.GetSize(1); ++j) {
            for (int i = 0; i<shape.GetSize(0); ++i) {
                int idx = eptlib::IJKToIdx(i,j,k, 2,4);
                ASSERT_EQ(shape(i,j,k), idx==0||idx==1||idx==2||idx==5);
            }
        }
    }
}

TEST(ShapeGTest,Operators) {
    eptlib::Shape cuboid1 = eptlib::shapes::Cuboid(3,5,1);
    cuboid1.Pad(1,0,2, 1,0,2);
    eptlib::Shape cuboid2 = eptlib::shapes::Cuboid(5,3,1);
    cuboid2.Pad(0,1,2, 0,1,2);
    eptlib::Shape cuboid3 = eptlib::shapes::Cuboid(5,1,5);
    cuboid3.Pad(0,2,0, 0,2,0);

    eptlib::Shape shape_union = cuboid1 + cuboid2 + cuboid3;
    boost::dynamic_bitset<> expected_union(std::string(
        "00000" "00000" "11111" "00000" "00000"
        "00000" "00000" "11111" "00000" "00000"
        "01110" "11111" "11111" "11111" "01110"
        "00000" "00000" "11111" "00000" "00000"
        "00000" "00000" "11111" "00000" "00000"
    ));
    ASSERT_TRUE(shape_union.GetData()==expected_union);

    eptlib::Shape shape_intersection = cuboid1 & cuboid2 & cuboid3;
    boost::dynamic_bitset<> expected_intersection(std::string(
        "00000" "00000" "00000" "00000" "00000"
        "00000" "00000" "00000" "00000" "00000"
        "00000" "00000" "01110" "00000" "00000"
        "00000" "00000" "00000" "00000" "00000"
        "00000" "00000" "00000" "00000" "00000"
    ));
    ASSERT_TRUE(shape_intersection.GetData()==expected_intersection);

    eptlib::Shape shape_difference = shape_union-shape_intersection;
    boost::dynamic_bitset<> expected_difference = expected_union - expected_intersection;
    ASSERT_TRUE(shape_difference.GetData()==expected_difference);
}

TEST(ShapeGTest,Pad) {
    eptlib::Shape cuboid = eptlib::shapes::Cuboid(3,2,1);
    cuboid.Pad(1,1,1, 2,2,2);
    boost::dynamic_bitset<> expected(std::string(
        "000000" "000000" "000000" "000000" "000000"
        "000000" "000000" "000000" "000000" "000000"
        "000000" "000000" "001110" "001110" "000000"
        "000000" "000000" "000000" "000000" "000000"
    ));
    ASSERT_TRUE(cuboid.GetData()==expected);
    return;
}

TEST(ShapeGTest,Shrink) {
    eptlib::Shape cuboid = eptlib::shapes::Cuboid(3,2,1);
    cuboid.Pad(1,1,1, 2,2,2);
    cuboid.Shrink(2,2,1, 1,1,1);
    boost::dynamic_bitset<> expected(std::string(
        "000" "000"
        "000" "011"
    ));
    ASSERT_TRUE(cuboid.GetData()==expected);
    return;
}

TEST(ShapeGTest,ShapesCuboid) {
    eptlib::Shape cuboid = eptlib::shapes::Cuboid(7,3,3);
    ASSERT_EQ(cuboid.GetSize(0), 7);
    ASSERT_EQ(cuboid.GetSize(1), 3);
    ASSERT_EQ(cuboid.GetSize(2), 3);
    ASSERT_TRUE(cuboid.GetData().all());
}

TEST(ShapeGTest,ShapesCuboidR) {
    eptlib::Shape cuboid = eptlib::shapes::CuboidR(3,1,1);
    ASSERT_EQ(cuboid.GetSize(0), 7);
    ASSERT_EQ(cuboid.GetSize(1), 3);
    ASSERT_EQ(cuboid.GetSize(2), 3);
    ASSERT_TRUE(cuboid.GetData().all());
}

TEST(ShapeGTest,ShapesEllipsoid) {
    eptlib::Shape ellipsoid = eptlib::shapes::Ellipsoid(3,3,1);
    ASSERT_EQ(ellipsoid.GetSize(0), 7);
    ASSERT_EQ(ellipsoid.GetSize(1), 7);
    ASSERT_EQ(ellipsoid.GetSize(2), 3);
    ASSERT_EQ(ellipsoid.GetVolume(), 31);
}

TEST(ShapeGTest,ShapesCross) {
    eptlib::Shape cross = eptlib::shapes::Cross(3,2,4);
    ASSERT_EQ(cross.GetSize(0), 7);
    ASSERT_EQ(cross.GetSize(1), 5);
    ASSERT_EQ(cross.GetSize(2), 9);
    ASSERT_EQ(cross.GetVolume(), 19);
}
