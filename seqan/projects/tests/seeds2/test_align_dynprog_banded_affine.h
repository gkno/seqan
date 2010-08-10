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

#ifndef TEST_SEEDS_TEST_ALIGN_DYNPROG_BANDED_AFFINE_H_
#define TEST_SEEDS_TEST_ALIGN_DYNPROG_BANDED_AFFINE_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.

// Test matrix resizing.
SEQAN_DEFINE_TEST(test_align_dynprog_banded_affine_resize_matrix)
{
    using namespace seqan;

    Matrix<int, 3> matrix;

    _alignBanded_resizeMatrix(matrix, CharString("length   11"), CharString("length 8"), -3, 1, Gotoh());

    SEQAN_ASSERT_EQ(12u, length(matrix, 0));
    SEQAN_ASSERT_EQ(7u, length(matrix, 1));
    SEQAN_ASSERT_EQ(3u, length(matrix, 2));
}


// Test gutter initialization if gap costs are free.
SEQAN_DEFINE_TEST(test_align_dynprog_banded_affine_init_gutter_free)
{
    using namespace seqan;

    Matrix<int, 3> matrix;
    setLength(matrix, 0, 4);
    setLength(matrix, 1, 5);
    setLength(matrix, 2, 3);
    resize(matrix);
    // fill(matrix, 42);

    _alignBanded_initGutter(matrix, Score<int, Simple>(1, -1, -1, -2), -1, 1, AlignConfig<true, true, true, true>(), Gotoh());

    int inf = InfimumValue<int>::VALUE / 2;

    // Test matrix M
    //
    // Diagonal below the lower one.
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 0, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 0, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 2, 0, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 3, 0, 0));
    // Diagonal above the upper one.
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 4, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 4, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 2, 4, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 3, 4, 0));
    // Top gutter
    SEQAN_ASSERT_EQ(0, value(matrix, 0, 2, 0));
    SEQAN_ASSERT_EQ(0, value(matrix, 0, 3, 0));
    // Left gutter.
    SEQAN_ASSERT_EQ(0, value(matrix, 1, 1, 0));
    //
    // Test Matrix I^a
    //
    // Top gutter.
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 2, 1));
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 3, 1));
    // Diagonal above the upper one.
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 4, 1));
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 4, 1));
    SEQAN_ASSERT_EQ(inf, value(matrix, 2, 4, 1));
    SEQAN_ASSERT_EQ(inf, value(matrix, 3, 4, 1));
    //
    // Test Matrix I^b
    //
    // Left gutter.
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 1, 2));
    // Diagonal below the lower one.
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 0, 2));
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 0, 2));
    SEQAN_ASSERT_EQ(inf, value(matrix, 2, 0, 2));
    SEQAN_ASSERT_EQ(inf, value(matrix, 3, 0, 2));
}


// Test gutter initialization if gap costs are not free.
SEQAN_DEFINE_TEST(test_align_dynprog_banded_affine_init_gutter_not_free)
{
    using namespace seqan;

    Matrix<int, 3> matrix;
    setLength(matrix, 0, 4);
    setLength(matrix, 1, 5);
    setLength(matrix, 2, 3);
    resize(matrix);

    _alignBanded_initGutter(matrix, Score<int, Simple>(1, -1, -1, -2), -1, 1, AlignConfig<false, false, true, true>(), Gotoh());

    int inf = InfimumValue<int>::VALUE / 2;

    // Test matrix M
    //
    // Diagonal below the lower one.
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 0, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 0, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 2, 0, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 3, 0, 0));
    // Diagonal above the upper one.
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 4, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 4, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 2, 4, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 3, 4, 0));
    // Top gutter
    SEQAN_ASSERT_EQ(0, value(matrix, 0, 2, 0));
    SEQAN_ASSERT_EQ(-1, value(matrix, 0, 3, 0));
    // Left gutter.
    SEQAN_ASSERT_EQ(-1, value(matrix, 1, 1, 0));
    //
    // Test Matrix I^a
    //
    // Top gutter.
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 2, 1));
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 3, 1));
    // Diagonal above the upper one.
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 4, 1));
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 4, 1));
    SEQAN_ASSERT_EQ(inf, value(matrix, 2, 4, 1));
    SEQAN_ASSERT_EQ(inf, value(matrix, 3, 4, 1));
    //
    // Test Matrix I^b
    //
    // Left gutter.
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 1, 2));
    // Diagonal below the lower one.
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 0, 2));
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 0, 2));
    SEQAN_ASSERT_EQ(inf, value(matrix, 2, 0, 2));
    SEQAN_ASSERT_EQ(inf, value(matrix, 3, 0, 2));
}


// Test DP matrix filling
SEQAN_DEFINE_TEST(test_align_dynprog_banded_affine_fill_matrix)
{
    using namespace seqan;

    Matrix<int, 3> matrix;
    DnaString const sequence0 = "CACCCC";
    DnaString const sequence1 = "CCACCC";
    Score<int, Simple> const scoringScheme(1, -1, -1);

    _alignBanded_resizeMatrix(matrix, sequence0, sequence1, -1, 1, Gotoh());
    _alignBanded_initGutter(matrix, scoringScheme, -1, 1, AlignConfig<false, false, false, false>(), Gotoh());
    _alignBanded_fillMatrix(matrix, sequence0, sequence1, scoringScheme, -1, 1, Gotoh());

    int inf = InfimumValue<int>::VALUE / 2;

    // First, the gutters and border diagonals should not have been touched.
    //
    // Test matrix M
    //
    // Diagonal below the lower one.
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 0, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 0, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 2, 0, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 3, 0, 0));
    // Diagonal above the upper one.
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 4, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 4, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 2, 4, 0));
    SEQAN_ASSERT_EQ(inf, value(matrix, 3, 4, 0));
    // Top gutter
    SEQAN_ASSERT_EQ(0, value(matrix, 0, 2, 0));
    SEQAN_ASSERT_EQ(-1, value(matrix, 0, 3, 0));
    // Left gutter.
    SEQAN_ASSERT_EQ(-1, value(matrix, 1, 1, 0));
    //
    // Test Matrix I^a
    //
    // Top gutter.
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 2, 1));
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 3, 1));
    // Diagonal above the upper one.
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 4, 1));
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 4, 1));
    SEQAN_ASSERT_EQ(inf, value(matrix, 2, 4, 1));
    SEQAN_ASSERT_EQ(inf, value(matrix, 3, 4, 1));
    //
    // Test Matrix I^b
    //
    // Left gutter.
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 1, 2));
    // Diagonal below the lower one.
    SEQAN_ASSERT_EQ(inf, value(matrix, 0, 0, 2));
    SEQAN_ASSERT_EQ(inf, value(matrix, 1, 0, 2));
    SEQAN_ASSERT_EQ(inf, value(matrix, 2, 0, 2));
    SEQAN_ASSERT_EQ(inf, value(matrix, 3, 0, 2));

    // Second, check the fields in between.
    SEQAN_ASSERT_EQ(-1, value(matrix, 1, 1));
    SEQAN_ASSERT_EQ(1, value(matrix, 1, 2));
    SEQAN_ASSERT_EQ(0, value(matrix, 1, 3));
    // row 2
    SEQAN_ASSERT_EQ(0, value(matrix, 2, 1));
    SEQAN_ASSERT_EQ(0, value(matrix, 2, 2));
    SEQAN_ASSERT_EQ(1, value(matrix, 2, 3));
    // row 3
    SEQAN_ASSERT_EQ(1, value(matrix, 3, 1));
    SEQAN_ASSERT_EQ(0, value(matrix, 3, 2));
    SEQAN_ASSERT_EQ(2, value(matrix, 3, 3));
    // row 4
    SEQAN_ASSERT_EQ(0, value(matrix, 4, 1));
    SEQAN_ASSERT_EQ(1, value(matrix, 4, 2));
    SEQAN_ASSERT_EQ(3, value(matrix, 4, 3));
    // row 5
    SEQAN_ASSERT_EQ(1, value(matrix, 5, 1));
    SEQAN_ASSERT_EQ(2, value(matrix, 5, 2));
    SEQAN_ASSERT_EQ(4, value(matrix, 5, 3));
    // row 6
    SEQAN_ASSERT_EQ(2, value(matrix, 6, 1));
    SEQAN_ASSERT_EQ(3, value(matrix, 6, 2));
}


SEQAN_DEFINE_TEST(test_align_dynprog_banded_affine_traceback)
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
        Score<int, Simple> const scoringScheme(1, -1, -1);
        
        _alignBanded_resizeMatrix(matrix, sequence0, sequence1, -2, 1, Gotoh());
        _alignBanded_initGutter(matrix, scoringScheme, -2, 1, AlignConfig<false, false, false, false>(), Gotoh());
        _alignBanded_fillMatrix(matrix, sequence0, sequence1, scoringScheme, -2, 1, Gotoh());

        // TODO(holtgrew): Debug output, remove when not needed any more.
		{
			for (int k = 0; k < 3; ++k) {
				std::cout << ",-- *** filled banded alignment matrix " << k << std::endl;
				for (unsigned i = 0; i < length(matrix, 0); ++i) {
					std::cout << "| ";
					for (unsigned j = 0; j < i; ++j)
						std::cout << "\t";
					for (unsigned j = 0; j < length(matrix, 1); ++j) {
						if (value(matrix, i, j, k) < InfimumValue<int>::VALUE / 4)
							std::cout << "\tinf";
						else
							std::cout << "\t" << value(matrix, i, j, k);
					}
					std::cout << std::endl;
				}
				std::cout << "`--" << std::endl;
			}
		}

        // Perform the traceback.
        Align<TString> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequence0);
        assignSource(row(alignment, 1), sequence1);

        size_t finalPos0 = 4;
        size_t finalPos1 = 2;
        TStringIterator seq0It = end(sequence0) - 1;
        TStringIterator seq1It = end(sequence1) - 1;
        TAlignRowIterator align0It = end(row(alignment, 0));
        TAlignRowIterator align1It = end(row(alignment, 1));
        int score = _alignBanded_traceBack(align0It, align1It, seq0It, seq1It, finalPos0, finalPos1, matrix, scoringScheme, 5, 3, 2, 3, false, AlignConfig<false, false, false, false>(), Gotoh());

        // std::cout << alignment;
        SEQAN_ASSERT_EQ(score, 1);
        SEQAN_ASSERT_TRUE(seq0It == begin(sequence0));
        // std::cout << (begin(sequence1, Standard()) - seq1It) << std::endl;
        SEQAN_ASSERT_TRUE(seq1It + 1 == begin(sequence1));
        // TODO(holtgrew): Why does this not work?
        // SEQAN_ASSERT_TRUE(align0It == begin(row(alignment, 0)));
        // SEQAN_ASSERT_TRUE(align1It == begin(row(alignment, 1)));
        SEQAN_ASSERT_EQ(finalPos0, 0u);
        SEQAN_ASSERT_EQ(finalPos1, static_cast<TPosition>(-1));
        SEQAN_ASSERT_TRUE(row(alignment, 0) == "CCAAA");
        // Note that leading and trailing gaps are not saved.  The
        // following is postfixed with a gap on visual inspection.
        SEQAN_ASSERT_TRUE(row(alignment, 1) == "C-AA");
    }
    // TODO(holtgrew): Case with free begin and end gaps.
}

#endif  // TEST_SEEDS_TEST_ALIGN_DYNPROG_BANDED_AFFINE_H_
