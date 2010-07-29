/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de
 ===========================================================================
  Copyright (C) 2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
 ===========================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
 ===========================================================================
  Tests for the banded DP programming algorithms in the seeds module.
 ===========================================================================
*/

#ifndef TEST_SEEDS_TEST_ALIGN_DYNPROG_BANDED_LINEAR_H_
#define TEST_SEEDS_TEST_ALIGN_DYNPROG_BANDED_LINEAR_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.

// Test matrix resizing.
SEQAN_DEFINE_TEST(test_align_dynprog_banded_linear_resize_matrix)
{
    using namespace seqan;

    Matrix<int, 2> matrix;

    _alignBanded_resizeMatrix(matrix, CharString("length 9"), CharString("length    11"), 2, NeedlemanWunsch());

    SEQAN_ASSERT_EQ(9u, length(matrix, 0));
    SEQAN_ASSERT_EQ(4u, length(matrix, 1));
}


// Test gutter initialization if gap costs are free.
SEQAN_DEFINE_TEST(test_align_dynprog_banded_linear_init_gutter_free)
{
    using namespace seqan;

    Matrix<int, 2> matrix;
    setLength(matrix, 0, 4);
    setLength(matrix, 1, 4);
    resize(matrix);

    _alignBanded_initGutter(matrix, Score<int, Simple>(1, -1, -2), AlignConfig<true, true, true, true>(), 2, NeedlemanWunsch());

    int inf = InfimumValue<int>::VALUE / 2;

    // "left" gutter
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 2, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 3, 0));
    // top gutter
    SEQAN_ASSERT_EQ(0, value(matrix, 0, 0));
    SEQAN_ASSERT_EQ(0, value(matrix, 0, 1));
    SEQAN_ASSERT_EQ(0, value(matrix, 0, 2));
    SEQAN_ASSERT_EQ(0, value(matrix, 0, 3));
    // right gutter must be very small so we never go there.
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 3));
    SEQAN_ASSERT_EQ(inf, value(matrix, 2, 3));
    SEQAN_ASSERT_EQ(inf, value(matrix, 3, 3));
}


// Test gutter initialization if gap costs are not free.
SEQAN_DEFINE_TEST(test_align_dynprog_banded_linear_init_gutter_not_free)
{
    using namespace seqan;

    Matrix<int, 2> matrix;
    setLength(matrix, 0, 4);
    setLength(matrix, 1, 4);
    resize(matrix);

    _alignBanded_initGutter(matrix, Score<int, Simple>(1, -1, -2), AlignConfig<false, false, true, true>(), 2, NeedlemanWunsch());

    int inf = InfimumValue<int>::VALUE / 2;

    // "left" gutter
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 2, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 3, 0));
    // top gutter
    SEQAN_ASSERT_EQ(0, value(matrix, 0, 0));
    SEQAN_ASSERT_EQ(-2, value(matrix, 0, 1));
    SEQAN_ASSERT_EQ(-4, value(matrix, 0, 2));
    SEQAN_ASSERT_EQ(-6, value(matrix, 0, 3));
    // right gutter must be very small so we never go there.
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 3));
    SEQAN_ASSERT_EQ(inf, value(matrix, 2, 3));
    SEQAN_ASSERT_EQ(inf, value(matrix, 3, 3));
}


// Test DP matrix filling
SEQAN_DEFINE_TEST(test_align_dynprog_banded_linear_fill_matrix)
{
    using namespace seqan;

    SEQAN_ASSERT_FAIL("Write me for banded!");

    Matrix<int, 2> matrix;
    DnaString const sequence0 = "CCA";
    DnaString const sequence1 = "CAA";
    Score<int, Simple> const scoringScheme(1, -1, -1);

    _align_resizeMatrix(matrix, sequence0, sequence1, NeedlemanWunsch());
    _align_initGutter(matrix, scoringScheme, AlignConfig<false, false, false, false>(), NeedlemanWunsch());
    _align_fillMatrix(matrix, sequence0, sequence1, scoringScheme, NeedlemanWunsch());

    int const expected[16] = {
         0, -1, -2, -3,
        -1,  1,  0, -1,
        -2,  0,  0, -1,
        -3, -1,  1,  1
    };

    for (unsigned j = 0; j < 4; ++j) {
        for (unsigned i = 0; i < 4; ++i) {
            SEQAN_ASSERT_EQ_MSG(expected[i * 4 + j], value(matrix, i, j), "i = %d, j = %d", i, j);
        }
    }
}

#endif  // TEST_SEEDS_TEST_ALIGN_DYNPROG_BANDED_LINEAR_H_
