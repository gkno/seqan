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
  Tests for the banded chain alignment.  Contains code that is based on
  the test code by Carsten Kemena.
 ===========================================================================
*/

#ifndef TEST_SEEDS_TEST_ALIGN_CHAIN_BANDED_H_
#define TEST_SEEDS_TEST_ALIGN_CHAIN_BANDED_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.

// Test overlap computation functions.
SEQAN_DEFINE_TEST(test_align_chain_banded_compute_upper_left_overlap)
{
    using namespace seqan;

    typedef Seed<Simple> TSimpleSeed;
    typedef Score<int, Simple> TScoringScheme;

    CharString sequence0 = "unimportant";
    CharString sequence1 = "also unimportant";
    _AlignmentChain<CharString, TScoringScheme, NeedlemanWunsch> alignmentChain(2, TScoringScheme(), sequence0, sequence1);

    // Test with a seed where start/end = lower/upper diagonal
    {
        TSimpleSeed seed(2, 1, 4, 5);
        // Diagonals: start/lower = -1, end/upper = 1;

        // Compute overlaps.
        unsigned overlap0, overlap1;
        _computeUpperLeftOverlap(overlap0, overlap1, seed, alignmentChain);
        SEQAN_ASSERT_EQ(2u, overlap0);
        SEQAN_ASSERT_EQ(4u, overlap1);
    }
    // Test with a seed where start/end = upper/lower diagonal
    {
        TSimpleSeed seed(1, 2, 5, 4);
        // Diagonals: start/upper = 1, end/lower = -1;

        // Compute overlaps.
        unsigned overlap0, overlap1;
        _computeUpperLeftOverlap(overlap0, overlap1, seed, alignmentChain);
        SEQAN_ASSERT_EQ(4u, overlap0);
        SEQAN_ASSERT_EQ(2u, overlap1);
    }
    // Test with a seed where start = end = lower = upper diagonal
    {
        TSimpleSeed seed(1, 1, 5, 5);
        // Diagonals: all are = 0

        // Compute overlaps.
        unsigned overlap0, overlap1;
        _computeUpperLeftOverlap(overlap0, overlap1, seed, alignmentChain);
        SEQAN_ASSERT_EQ(3u, overlap0);
        SEQAN_ASSERT_EQ(3u, overlap1);
    }
    // Test with a seed where {start, end} != {lower, upper diagonal}.
    {
        TSimpleSeed seed(1, 2, 5, 4);
        // Diagonals: start = 1, end = -1;
        setUpperDiagonal(seed, 3);
        setLowerDiagonal(seed, -2);

        // Compute overlaps.
        unsigned overlap0, overlap1;
        _computeUpperLeftOverlap(overlap0, overlap1, seed, alignmentChain);
        SEQAN_ASSERT_EQ(4u, overlap0);
        SEQAN_ASSERT_EQ(2u, overlap1);
    }
}


// Test overlap computation functions.
SEQAN_DEFINE_TEST(test_align_chain_banded_compute_lower_right_overlap)
{
    using namespace seqan;

    typedef Seed<Simple> TSimpleSeed;
    typedef Score<int, Simple> TScoringScheme;

    CharString sequence0 = "unimportant";
    CharString sequence1 = "also unimportant";
    _AlignmentChain<CharString, TScoringScheme, NeedlemanWunsch> alignmentChain(2, TScoringScheme(), sequence0, sequence1);

    // Test with a seed where start/end = lower/upper diagonal
    {
        TSimpleSeed seed(2, 1, 4, 5);
        // Diagonals: start/lower = -1, end/upper = 1;

        // Compute overlaps.
        unsigned overlap0, overlap1;
        _computeLowerRightOverlap(overlap0, overlap1, seed, alignmentChain);
        SEQAN_ASSERT_EQ(2u, overlap0);
        SEQAN_ASSERT_EQ(4u, overlap1);
    }
    // Test with a seed where start/end = upper/lower diagonal
    {
        TSimpleSeed seed(1, 2, 5, 4);
        // Diagonals: start/upper = 1, end/lower = -1;

        // Compute overlaps.
        unsigned overlap0, overlap1;
        _computeLowerRightOverlap(overlap0, overlap1, seed, alignmentChain);
        SEQAN_ASSERT_EQ(4u, overlap0);
        SEQAN_ASSERT_EQ(2u, overlap1);
    }
    // Test with a seed where start = end = lower = upper diagonal
    {
        TSimpleSeed seed(1, 1, 5, 5);
        // Diagonals: all are = 0

        // Compute overlaps.
        unsigned overlap0, overlap1;
        _computeLowerRightOverlap(overlap0, overlap1, seed, alignmentChain);
        SEQAN_ASSERT_EQ(3u, overlap0);
        SEQAN_ASSERT_EQ(3u, overlap1);
    }
    // Test with a seed where {start, end} != {lower, upper diagonal}.
    {
        TSimpleSeed seed(1, 2, 5, 4);
        // Diagonals: start = 1, end = -1;
        setUpperDiagonal(seed, 3);
        setLowerDiagonal(seed, -2);

        // Compute overlaps.
        unsigned overlap0, overlap1;
        _computeLowerRightOverlap(overlap0, overlap1, seed, alignmentChain);
        SEQAN_ASSERT_EQ(4u, overlap0);
        SEQAN_ASSERT_EQ(2u, overlap1);
    }
}

// Test banded alignment algorithm around a chain.  Linear gap cost
// case.
SEQAN_DEFINE_TEST(test_align_chain_banded_align_linear)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    
    // Test on whole strings.
    {
        // Resulting alignment should be something like this (seeds
        // are marked by < and >).
        //
        //     > <    > <    >  <
        //   GGCGATNNNCAT--GGCACA
        //   --CGA-ATCCATCCCACACA
        CharString sequence0 = "GGCGATNNNCATGGCACA";
        CharString sequence1 = "CGAATCCATCCCACACA";  // TODO(holtgrew): Switch back again.
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(2, 0, 6, 5));
        appendValue(seedChain, TSeed(9, 6, 12, 9));
        appendValue(seedChain, TSeed(14, 11, 16, 17));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequence0);
        assignSource(row(alignment, 1), sequence1);

        int result = bandedChainAlignment(alignment, seedChain, 1, scoringScheme, AlignConfig<false, false, false, false>());
        SEQAN_ASSERT_EQ(result, 5);

        // Compare alignment rows.
        SEQAN_ASSERT_TRUE(row(alignment, 0) == "GGCGATNNNCAT--GGCACA");
        // Note that leading gaps are not stored in the row itself,
        // when printed, there will be a "--" prepended.
        SEQAN_ASSERT_TRUE(row(alignment, 1) == "CGA-ATCCATCCCACACA");
    }
    // Test on infixes.
    {
        // The test data is the same as above but with a T and a TT prepended.
        CharString sequence0 = "TGGCGATNNNCATGGCACA";
        CharString sequence1 = "TTCGAATCCATCCCACACA";
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed > seedChain;
        appendValue(seedChain, TSeed(2, 0, 6, 5));
        appendValue(seedChain, TSeed(9, 6, 12, 9));
        appendValue(seedChain, TSeed(14, 11, 16, 17));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        // TODO(holtgrew): This infix assignment for alignments is a bit creepy, maybe one of the too many shortcuts?
        assignSource(row(alignment, 0), sequence0, 1, length(sequence0));
        assignSource(row(alignment, 1), sequence1, 2, length(sequence1));

        //cout << "Score: " << bandedChainAlignment(seedChain1, 2, alignment2, scoreMatrix) << endl;
        int result = bandedChainAlignment(alignment, seedChain, 1, scoringScheme, AlignConfig<false, false, false, false>());
        SEQAN_ASSERT_EQ(result, 5);

        // Compare alignment rows.
        SEQAN_ASSERT_TRUE(row(alignment, 0) == "GGCGATNNNCAT--GGCACA");
        // Note that leading gaps are not stored in the row itself,
        // when printed, there will be a "--" prepended.
        SEQAN_ASSERT_TRUE(row(alignment, 1) == "CGA-ATCCATCCCACACA");
    }
}


// Test banded alignment algorithm around a chain.  Affine gap cost
// case.
SEQAN_DEFINE_TEST(test_align_chain_banded_align_affine)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;

    /*
    // Test on whole strings.
    {
        CharString query = "ACGTCCTCGTACACCGTCTTAA";
        CharString database = "TACGATCCACACCGCGTCT";
        Score<int, Simple> scoringScheme(3, -2, -1, -3);
        
        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(1, 2, 6, 8));
        appendValue(seedChain, TSeed(12, 10, 19, 19));
        
        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), query);
        assignSource(row(alignment, 1), database);
        
        //cout << "Score: " << bandedChainAlignment(seedChain2, 2, alignment3, scoreMatrix2) << endl;
        int result = bandedChainAlignment(alignment, seedChain, 2, scoringScheme, AlignConfig<false, false, false, false>());
        SEQAN_ASSERT_EQ(result, 24);
        
        //cout << alignment3 << endl;
        SEQAN_ASSERT_TRUE(row(alignment, 0) == "ACG-TCCTCGTACAC--CGTCTTAA");
        SEQAN_ASSERT_TRUE(row(alignment, 1) == "TACGATCC----ACACCGCGTCT");
    }
    // Test on infixes.
    {
        CharString query = "ACGTCCTCGTACACCGTCTTAA";
        CharString database = "TACGATCCACACCGCGTCT";
        Score<int, Simple> scoringScheme(3, -2, -1, -3);
        
        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(1, 2, 6, 8));
        appendValue(seedChain, TSeed(12, 10, 19, 19));
        
        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), query, 1, length(query));
        assignSource(row(alignment, 1), database, 2, length(database));
        
        //cout << "Score: " << bandedChainAlignment(seedChain2, 2, alignment4, scoreMatrix2) << endl;
        int result = bandedChainAlignment(alignment, seedChain, 2, scoringScheme, AlignConfig<false, false, false, false>());
        SEQAN_ASSERT_EQ(result, 21);
        
        //cout << alignment4 << endl;
        SEQAN_ASSERT_TRUE(row(alignment, 0) == "CG-TCCTCGTACAC--CGTCTTAA");
        SEQAN_ASSERT_TRUE(row(alignment, 1) == "CGATCC----ACACCGCGTCT");
    }
    */
    // Test on whole strings -- linear gap costs.
    {
        // Resulting alignment should be something like this (seeds
        // are marked by < and >).
        //
        //     > <    > <    >  <
        //   GGCGATNNNCAT--GGCACA
        //   --CGA-ATCCATCCCACACA
        CharString sequence0 = "GGCGATNNNCATGGCACA";
        CharString sequence1 = "CGAATCCATCCCACACA";  // TODO(holtgrew): Switch back again.
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(2, 0, 6, 5));
        appendValue(seedChain, TSeed(9, 6, 12, 9));
        appendValue(seedChain, TSeed(14, 11, 16, 17));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequence0);
        assignSource(row(alignment, 1), sequence1);

        int result = bandedChainAlignment(alignment, seedChain, 1, scoringScheme, AlignConfig<false, false, false, false>());
        SEQAN_ASSERT_EQ(result, 5);

        // Compare alignment rows.
        SEQAN_ASSERT_TRUE(row(alignment, 0) == "GGCGATNNNCAT--GGCACA");
        // Note that leading gaps are not stored in the row itself,
        // when printed, there will be a "--" prepended.
        SEQAN_ASSERT_TRUE(row(alignment, 1) == "CGA-ATCCATCCCACACA");
    }
    // Test on infixes -- linear gap costs.
    {
        // The test data is the same as above but with a T and a TT prepended.
        CharString sequence0 = "TGGCGATNNNCATGGCACA";
        CharString sequence1 = "TTCGAATCCATCCCACACA";
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed > seedChain;
        appendValue(seedChain, TSeed(2, 0, 6, 5));
        appendValue(seedChain, TSeed(9, 6, 12, 9));
        appendValue(seedChain, TSeed(14, 11, 16, 17));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        // TODO(holtgrew): This infix assignment for alignments is a bit creepy, maybe one of the too many shortcuts?
        assignSource(row(alignment, 0), sequence0, 1, length(sequence0));
        assignSource(row(alignment, 1), sequence1, 2, length(sequence1));

        //cout << "Score: " << bandedChainAlignment(seedChain1, 2, alignment2, scoreMatrix) << endl;
        int result = bandedChainAlignment(alignment, seedChain, 1, scoringScheme, AlignConfig<false, false, false, false>());
        SEQAN_ASSERT_EQ(result, 5);

        // Compare alignment rows.
        SEQAN_ASSERT_TRUE(row(alignment, 0) == "GGCGATNNNCAT--GGCACA");
        // Note that leading gaps are not stored in the row itself,
        // when printed, there will be a "--" prepended.
        SEQAN_ASSERT_TRUE(row(alignment, 1) == "CGA-ATCCATCCCACACA");
    }
}

#endif  // TEST_SEEDS_TEST_ALIGN_CHAIN_BANDED_H_
