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
  Tests for the classic DP programming algorithms in the seeds module.
 ===========================================================================
*/

#ifndef TEST_SEEDS_TEST_ALIGN_DYNPROG_AFFINE_H_
#define TEST_SEEDS_TEST_ALIGN_DYNPROG_AFFINE_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.

// Test matrix resizing.
SEQAN_DEFINE_TEST(test_align_dynprog_affine_resize_matrix)
{
    using namespace seqan;

    Matrix<int, 3> matrix;

    _align_resizeMatrix(matrix, CharString("length 9"), CharString("length    11"), Gotoh());

    SEQAN_ASSERT_EQ(9u, length(matrix, 0));
    SEQAN_ASSERT_EQ(13u, length(matrix, 1));
    SEQAN_ASSERT_EQ(3u, length(matrix, 2));
}

// Test gutter initialization if gap costs are free.
SEQAN_DEFINE_TEST(test_align_dynprog_affine_init_gutter_free)
{
    using namespace seqan;

    Matrix<int, 3> matrix;
    setLength(matrix, 0, 2);
    setLength(matrix, 1, 3);
    setLength(matrix, 2, 3);
    resize(matrix);

    _align_initGutter(matrix, Score<int, Simple>(1, -1, -2), AlignConfig<true, true, true, true>(), Gotoh());

    SEQAN_ASSERT_EQ(0, value(matrix, 0, 0));
    SEQAN_ASSERT_EQ(-2, value(matrix, 0, 1));
    SEQAN_ASSERT_EQ(-4, value(matrix, 0, 2));
    SEQAN_ASSERT_EQ(-2, value(matrix, 1, 0));
}

// Test gutter initialization if gap costs are not free.
SEQAN_DEFINE_TEST(test_align_dynprog_affine_init_gutter_not_free)
{
    using namespace seqan;

    Matrix<int, 3> matrix;
    setLength(matrix, 0, 2);
    setLength(matrix, 1, 3);
    setLength(matrix, 2, 3);
    resize(matrix);

    _align_initGutter(matrix, Score<int, Simple>(1, -1, -1, -2), AlignConfig<false, false, false, false>(), Gotoh());

    int inf = InfimumValue<int>::VALUE / 2;

    // Check M.
    SEQAN_ASSERT_EQ(0, value(matrix, 0, 0, 0));
    SEQAN_ASSERT_EQ(-1, value(matrix, 1, 0, 0));
    SEQAN_ASSERT_EQ(-1, value(matrix, 0, 1, 0));
    SEQAN_ASSERT_EQ(-2, value(matrix, 0, 2, 0));
    // Check I^a.
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 1, 1));
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 2, 1));
    // Check I^b.
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 0, 2));
}


// Test DP matrix filling
SEQAN_DEFINE_TEST(test_align_dynprog_affine_fill_matrix)
{
    using namespace seqan;

    Matrix<int, 3> matrix;
    DnaString const sequence0 = "CCA";
    DnaString const sequence1 = "CAA";
    Score<int, Simple> const scoringScheme(1, -1, -1, -2);

    _align_resizeMatrix(matrix, sequence0, sequence1, Gotoh());
    _align_initGutter(matrix, scoringScheme, AlignConfig<false, false, false, false>(), Gotoh());
    _align_fillMatrix(matrix, sequence0, sequence1, scoringScheme, Gotoh());

    int inf = InfimumValue<int>::VALUE / 2;

    // Check M.
    SEQAN_ASSERT_EQ(0, value(matrix, 0, 0, 0));
    SEQAN_ASSERT_EQ(-1, value(matrix, 0, 1, 0));
    SEQAN_ASSERT_EQ(-2, value(matrix, 0, 2, 0));
    SEQAN_ASSERT_EQ(-3, value(matrix, 0, 3, 0));
    SEQAN_ASSERT_EQ(-1, value(matrix, 1, 0, 0));
    SEQAN_ASSERT_EQ(1, value(matrix, 1, 1, 0));
    SEQAN_ASSERT_EQ(-1, value(matrix, 1, 2, 0));
    SEQAN_ASSERT_EQ(-2, value(matrix, 1, 3, 0));
    SEQAN_ASSERT_EQ(-2, value(matrix, 2, 0, 0));
    SEQAN_ASSERT_EQ(0, value(matrix, 2, 1, 0));
    SEQAN_ASSERT_EQ(0, value(matrix, 2, 2, 0));
    SEQAN_ASSERT_EQ(-2, value(matrix, 2, 3, 0));
    SEQAN_ASSERT_EQ(-3, value(matrix, 3, 0, 0));
    SEQAN_ASSERT_EQ(-2, value(matrix, 3, 1, 0));
    SEQAN_ASSERT_EQ(1, value(matrix, 3, 2, 0));
    SEQAN_ASSERT_EQ(1, value(matrix, 3, 3, 0));
    // Check I^a.
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 1, 1));
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 2, 1));
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 3, 1));
    SEQAN_ASSERT_EQ(-3, value(matrix, 1, 1, 1));
    SEQAN_ASSERT_EQ(-4, value(matrix, 1, 2, 1));
    SEQAN_ASSERT_EQ(-5, value(matrix, 1, 3, 1));
    SEQAN_ASSERT_EQ(-1, value(matrix, 2, 1, 1));
    SEQAN_ASSERT_EQ(-3, value(matrix, 2, 2, 1));
    SEQAN_ASSERT_EQ(-4, value(matrix, 2, 3, 1));
    SEQAN_ASSERT_EQ(-2, value(matrix, 3, 1, 1));
    SEQAN_ASSERT_EQ(-2, value(matrix, 3, 2, 1));
    SEQAN_ASSERT_EQ(-4, value(matrix, 3, 3, 1));
    // Check I^b.
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 0, 2));
    SEQAN_ASSERT_EQ(-3, value(matrix, 1, 1, 2));
    SEQAN_ASSERT_EQ(-1, value(matrix, 1, 2, 2));
    SEQAN_ASSERT_EQ(-2, value(matrix, 1, 3, 2));
    SEQAN_ASSERT_EQ(inf, value(matrix, 2, 0, 2));
    SEQAN_ASSERT_EQ(-4, value(matrix, 2, 1, 2));
    SEQAN_ASSERT_EQ(-2, value(matrix, 2, 2, 2));
    SEQAN_ASSERT_EQ(-4, value(matrix, 2, 3, 2));
    SEQAN_ASSERT_EQ(inf, value(matrix, 3, 0, 2));
    SEQAN_ASSERT_EQ(-5, value(matrix, 3, 1, 2));
    SEQAN_ASSERT_EQ(-4, value(matrix, 3, 2, 2));
    SEQAN_ASSERT_EQ(-3, value(matrix, 3, 3, 2));
}


SEQAN_DEFINE_TEST(test_align_dynprog_affine_traceback)
{
    using namespace seqan;

    typedef CharString TString;
    typedef Position<CharString>::Type TPosition;
    typedef Align<TString> TAlign;
    typedef Row<TAlign>::Type TAlignRow;
    typedef Iterator<TAlignRow, Standard>::Type TAlignRowIterator;
    typedef Iterator<TString, Standard>::Type TStringIterator;

    // Case: No free begin/end gaps.
    {
        // Fill the matrix with DP (tested by other tests).
        Matrix<int, 3> matrix;
        TString const sequence0 = "CCAAA";
        TString const sequence1 = "CAA";
        Score<int, Simple> const scoringScheme(1, -1, -2, -1);
        
        _align_resizeMatrix(matrix, sequence0, sequence1, Gotoh());
        _align_initGutter(matrix, scoringScheme, AlignConfig<false, false, false, false>(), Gotoh());
        _align_fillMatrix(matrix, sequence0, sequence1, scoringScheme, Gotoh());

        // Perform the traceback.
        Align<TString> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequence0);
        assignSource(row(alignment, 1), sequence1);

        size_t finalPos0 = 0;
        size_t finalPos1 = 0;
        TStringIterator seq0It = end(sequence0) - 1;
        TStringIterator seq1It = end(sequence1) - 1;
        TAlignRowIterator align0It = end(row(alignment, 0));
        TAlignRowIterator align1It = end(row(alignment, 1));
        int score = _align_traceBack(align0It, align1It, seq0It, seq1It, finalPos0, finalPos1, matrix, scoringScheme, 1, 1, true, AlignConfig<false, false, false, false>(), Gotoh());

        SEQAN_ASSERT_EQ(score, 1);
        SEQAN_ASSERT_TRUE(seq0It == begin(sequence0));
        SEQAN_ASSERT_TRUE(seq1It == begin(sequence1));
        // TODO(holtgrew): Why does this not work?
        // SEQAN_ASSERT_TRUE(align0It == begin(row(alignment, 0)));
        // SEQAN_ASSERT_TRUE(align1It == begin(row(alignment, 1)));
        SEQAN_ASSERT_EQ(finalPos0, static_cast<TPosition>(-1));
        SEQAN_ASSERT_EQ(finalPos1, 0u);
        // Expected alignment:
        //
        //  CCAAA
        //  --CAA
        SEQAN_ASSERT_TRUE(row(alignment, 0) == "CCAAA");
        // Leading gaps are not shown, so there should be two gaps in
        // front of the following.
        SEQAN_ASSERT_TRUE(row(alignment, 1) == "CAA");
    }
    // TODO(holtgrew): Case with free begin and end gaps.
}

#endif  // TEST_SEEDS_TEST_ALIGN_DYNPROG_AFFINE_H_
