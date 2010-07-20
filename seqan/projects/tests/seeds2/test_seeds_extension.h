/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
  ============================================================================
  Copyright (C) 2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  ============================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ============================================================================
  Tests for the header seeds_extension.h.
  ==========================================================================*/

#ifndef TEST_SEEDS_TEST_SEEDS_EXTENSION_H_
#define TEST_SEEDS_TEST_SEEDS_EXTENSION_H_

template <typename TSeedSpec>
void
testSeedsExtensionMatchExtension(TSeedSpec const &)
{
    using namespace seqan;

    typedef Seed<TSeedSpec> TSeed;
    // Test extension to the left only.  Extension possible by 2.
    {
        DnaString query = "AAACCCCCCAA";
        DnaString database = "CCCCCCCCC";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_LEFT, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(3, 2, 4), seed);
    }
    // Test extension to the left only.  No extension is possible.
    {
        DnaString query = "AAAAACCCCAA";
        DnaString database = "CCCCCCCCC";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_LEFT, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 2), seed);
    }
    // Test extension to the right only.  Extension possible by 2.
    {
        DnaString query = "AAAAAAAAACC";
        DnaString database = "AAAAAAAAA";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_RIGHT, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 4), seed);
    }
    // Test extension to the right only.  No extension is possible.
    {
        DnaString query = "AAAAAAAAAA";
        DnaString database = "AAAAAACC";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_RIGHT, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 2), seed);
    }
    // Test extension in both directions.  Extension possible by 2 in both directions.
    {
        DnaString query = "AAACCCCCCC";
        DnaString database = "GGCCCCCCAA";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_BOTH, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(3, 2, 6), seed);
    }
    // Test extension in both directions.  Extension possible by 2 to the left.
    {
        DnaString query = "AAACCCCCCC";
        DnaString database = "GGCCCCAAAA";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_BOTH, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(3, 2, 4), seed);
    }
    // Test extension in both directions.  Extension possible by 2 to the right.
    {
        DnaString query = "AAAAACCCCC";
        DnaString database = "GGCCCCCCAA";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_BOTH, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 4), seed);
    }
    // Test extension in both directions.  No Extension possible in either direction.
    {
        DnaString query = "AAAAAAAAAAA";
        DnaString database = "CCCCAACCC";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_BOTH, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 2), seed);
    }
}

// Test the seed extension algorithm with match extension for simple seeds.
SEQAN_DEFINE_TEST(test_seeds_extension_match_extension_simple)
{
    using namespace seqan;
    testSeedsExtensionMatchExtension(Simple());
}

// Test the seed extension algorithm with match extension for chained seeds.
SEQAN_DEFINE_TEST(test_seeds_extension_match_extension_chained)
{
    using namespace seqan;
    testSeedsExtensionMatchExtension(ChainedSeed());
}

#endif  // TEST_SEEDS_TEST_SEEDS_EXTENSION_H_
