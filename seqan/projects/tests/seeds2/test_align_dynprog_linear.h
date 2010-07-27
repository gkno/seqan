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

#ifndef TEST_SEEDS_TEST_ALIGN_DYNPROG_LINEAR_H_
#define TEST_SEEDS_TEST_ALIGN_DYNPROG_LINEAR_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.

// Test DP algorithm.
SEQAN_DEFINE_TEST(test_seeds_align_dynprog_linear_basic)
{
    typedef Matrix<int, 2> TMatrix;
    // TODO(holtgrew): Make matrices const-iterable, SeqAn is broken here.
    typedef Iterator<TMatrix>::Type TMatrixIterator;
    
    TMatrix dpMatrix;
    CharString sequence0 = "CAC";
    CharString sequence1 = "CCAA";
    Score<int, simple> scoringScheme(1, -1, -1);

    _alignDynprogNeedlemanWunsch(dpMatrix, sequence0, sequence1, AlignConfig<true, true, false, false>());

    SEQAN_ASSERT_EQ(4u, length(dpMatrix, 0));
    SEQAN_ASSERT_EQ(5u, length(dpMatrix, 1));
    int const expected[] = {
         0,  0,  0,  0,  0,
         0,  1,  1,  0, -1,
         0,  0,  0,  2,  1,
         0, -1, -1,  1,  1
    };
    int i = 0;
    for (TMatrixIterator it = begin(dpMatrix); it != end(dpMatrix); ++it, ++i)
        SEQAN_ASSERT_EQ_MSG(*it, expected[i], "i = %d", i);
}

#endif  // TEST_SEEDS_TEST_ALIGN_DYNPROG_LINEAR_H_
