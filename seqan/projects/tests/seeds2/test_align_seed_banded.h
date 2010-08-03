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
  Author: Carsten Kemena <carsten.kemena@crg.es>
 ===========================================================================
  Tests for the banded seed alignment.  Contains code that is based  on the
  test code by Carsten Kemena.
 ===========================================================================
*/

#ifndef TEST_SEEDS_TEST_ALIGN_SEED_BANDED_H_
#define TEST_SEEDS_TEST_ALIGN_SEED_BANDED_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.

// Test Needleman-Wunsch banded alignment around seeds.
SEQAN_DEFINE_TEST(test_align_seed_banded_needleman_wunsch)
{
    using namespace seqan;

    // In this test, we align CharStrings since this is easiest for
    // later comparison: The gap character is already embedded in the
    // alphabet.
    //
    // TODO(holtgrew): Maybe change to DnaString anyway? conversion to CharString could work well...
    {
        CharString query = "cgtacgtgagtga";
        CharString database = "cgatta";

        Seed<Simple> seed(0, 0, 9, 6);
        Score<int, Simple> scoringScheme(3, -3, -2);

        Align<CharString, ArrayGaps> alignment;
        Segment<CharString, InfixSegment> seg1(query, getBeginDim0(seed), getEndDim0(seed));
        Segment<CharString, InfixSegment> seg2(database, getBeginDim1(seed), getEndDim1(seed));
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), seg1);
        assignSource(row(alignment, 1), seg2);

        int result = bandedAlignment(alignment, seed, 1, scoringScheme);
        SEQAN_ASSERT_EQ(6, result);
        SEQAN_ASSERT_TRUE(row(alignment, 0) == "cgtacgtga");
        SEQAN_ASSERT_TRUE(row(alignment, 1) == "cg-at-t-a");
    }
    {
        CharString query = "cgtacgtgagtga";
        CharString database = "cgatta";

        Seed<Simple> seed(0, 0, 9, 6);
        Score<int, Simple> scoringScheme(3, -3, -2);
        setLowerDiagonal(seed, 1);

        Align<CharString, ArrayGaps> alignment;
        Segment<CharString, InfixSegment> seg1(query, getBeginDim0(seed), getEndDim0(seed));
        Segment<CharString, InfixSegment> seg2(database, getBeginDim1(seed), getEndDim1(seed));
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), seg1);
        assignSource(row(alignment, 1), seg2);

        int result = bandedAlignment(alignment, seed, 1, scoringScheme);
        SEQAN_ASSERT_EQ(result, 6);
        SEQAN_ASSERT_TRUE(row(alignment, 0) == "cgtacgtga");
        SEQAN_ASSERT_TRUE(row(alignment, 1) == "cg-at-t-a");
    }
    {
        CharString query = "ACTTTCATTTT";
        CharString database = "ACTGTTCAGGG";

        Seed<Simple> seed(0, 0, 4, 5);
        Score<int, Simple> scoringScheme(2, -1, -1);

        Align<CharString, ArrayGaps> alignment;
        Segment<CharString, InfixSegment> seg1(query, getBeginDim0(seed), getEndDim0(seed));
        Segment<CharString, InfixSegment> seg2(database, getBeginDim1(seed), getEndDim1(seed));
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), seg1);
        assignSource(row(alignment, 1), seg2);

        int result = bandedAlignment(alignment, seed, 2, scoringScheme);
        SEQAN_ASSERT_EQ(result, 7);
        SEQAN_ASSERT_TRUE(row(alignment, 0) == "ACT-T");
        SEQAN_ASSERT_TRUE(row(alignment, 1) == "ACTGT");
    }
}


// Test Gotoh banded alignment around seeds.
SEQAN_DEFINE_TEST(test_align_seed_banded_gotoh)
{
    using namespace seqan;

    // In this test, we align CharStrings since this is easiest for
    // later comparison: The gap character is already embedded in the
    // alphabet.
    //
    // TODO(holtgrew): Maybe change to DnaString anyway? conversion to CharString could work well...
    {
        CharString query = "cgtacgtgagtga";
        CharString database = "cgatta";

        Seed<Simple> seed(0, 0, 9, 6);
        setUpperDiagonal(seed, -3);
        Score<int, Simple> scoringScheme(3, -2, -1, -3);

        Align<CharString, ArrayGaps> alignment;
        Segment<CharString, InfixSegment> seg1(query, getBeginDim0(seed), getEndDim0(seed));
        Segment<CharString, InfixSegment> seg2(database, getBeginDim1(seed), getEndDim1(seed));
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), seg1);
        assignSource(row(alignment, 1), seg2);

        int result = bandedAlignment(alignment, seed, 1, scoringScheme);
        std::cout << alignment;
        SEQAN_ASSERT_EQ(result, 6);
        SEQAN_ASSERT_TRUE(row(alignment, 0) == "cgtacgtga");
        SEQAN_ASSERT_TRUE(row(alignment, 1) == "cg-a--tta");
    }
    {
        CharString query = "cgtacgtgagtga";
        CharString database = "cgatta";

        Seed<Simple> seed(0, 0, 9, 6);
        setUpperDiagonal(seed, -4);
        Score<int, Simple> scoringScheme(3, -2, -1, -3);

        Align<CharString, ArrayGaps> alignment;
        Segment<CharString, InfixSegment> seg1(query, getBeginDim0(seed), getEndDim0(seed));
        Segment<CharString, InfixSegment> seg2(database, getBeginDim1(seed), getEndDim1(seed));
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), seg1);
        assignSource(row(alignment, 1), seg2);

        SEQAN_ASSERT_EQ(bandedAlignment(alignment, seed, 1, scoringScheme), 6);
        SEQAN_ASSERT_TRUE(row(alignment, 0) == "cgtacgtga");
        SEQAN_ASSERT_TRUE(row(alignment, 1) == "cg-a--tta");
    }
    {
        CharString query = "ACTTTCATTTT";
        CharString database = "ACTGTTCAGGG";

        Seed<Simple> seed(0, 0, 4, 5);
        Score<int, Simple> scoringScheme(3, -2, -1, -3);

        Align<CharString, ArrayGaps> alignment;
        Segment<CharString, InfixSegment> seg1(query, getBeginDim0(seed), getEndDim0(seed));
        Segment<CharString, InfixSegment> seg2(database, getBeginDim1(seed), getEndDim1(seed));
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), seg1);
        assignSource(row(alignment, 1), seg2);

        SEQAN_ASSERT_EQ(bandedAlignment(alignment, seed, 2, scoringScheme), 9);
        SEQAN_ASSERT_TRUE(row(alignment, 0) == "ACT-T");
        SEQAN_ASSERT_TRUE(row(alignment, 1) == "ACTGT");
    }
}

#endif  // TEST_SEEDS_TEST_ALIGN_SEED_BANDED_H_
