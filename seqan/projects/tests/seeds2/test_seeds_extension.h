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


template <typename TSeedSpec>
void
testSeedsExtensionUnGappedXDropExtension(TSeedSpec const &)
{
    using namespace seqan;

    typedef Seed<TSeedSpec> TSeed;

    // Test extension to the left only.  Extension possible by 2 over a gap of length 1.
    {
        DnaString query = "AAACACCCCAA";
        DnaString database = "CCCCCCCCC";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_LEFT, Score<int, EditDistance>(), 2, UnGappedXDrop());
        SEQAN_ASSERT_EQ(TSeed(3, 2, 4), seed);
    }
    // Test extension to the left only.  No extension is possible.
    {
        DnaString query = "AAAAACCACAA";
        DnaString database = "CCCCCCCCC";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_LEFT, Score<int, EditDistance>(), 2, UnGappedXDrop());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 2), seed);
    }
    // Test extension to the right only.  Extension possible by 2 over a gap of length 1.
    {
        DnaString query = "AAAAAAACACC";
        DnaString database = "AAAAAAAAA";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_RIGHT, Score<int, EditDistance>(), 2, UnGappedXDrop());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 4), seed);
    }
    // Test extension to the right only.  No extension is possible.
    {
        DnaString query = "AAAAAAAAAA";
        DnaString database = "AAAAAACC";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_RIGHT, Score<int, EditDistance>(), 2, UnGappedXDrop());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 2), seed);
    }
    // Test extension in both directions.  Extension possible by 2 in both directions over a gap of length 1 on each side.
    {
        DnaString query = "AAACCCCCCC";
        DnaString database = "GGCCCCCCAA";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_BOTH, Score<int, EditDistance>(), 2, UnGappedXDrop());
        SEQAN_ASSERT_EQ(TSeed(3, 2, 6), seed);
    }
    // Test extension in both directions.  Extension possible by 2 to the left over a gap of length 1.
    {
        DnaString query = "AAACCCCCCC";
        DnaString database = "GGCCCCAAAA";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_BOTH, Score<int, EditDistance>(), 2, UnGappedXDrop());
        SEQAN_ASSERT_EQ(TSeed(3, 2, 4), seed);
    }
    // Test extension in both directions.  Extension possible by 2 to the right over a gap of length 1.
    {
        DnaString query = "AAAAACCCCC";
        DnaString database = "GGCCCCCCAA";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_BOTH, Score<int, EditDistance>(), 2, UnGappedXDrop());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 4), seed);
    }
    // Test extension in both directions.  No Extension possible in either direction.
    {
        DnaString query = "AAAAAAAAAAA";
        DnaString database = "CCCCAACCC";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_BOTH, Score<int, EditDistance>(), 2, UnGappedXDrop());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 2), seed);
    }
}


template <typename TSeedSpec>
void
testSeedsExtensionGappedXDropExtension(TSeedSpec const &)
{
    using namespace seqan;

    typedef Seed<TSeedSpec> TSeed;
    SEQAN_ASSERT_FAIL("Write me!");
}


// Test the seed extension algorithm with match extension for simple seeds.
SEQAN_DEFINE_TEST(test_seeds_extension_match_extension_simple)
{
    using namespace seqan;
    testSeedsExtensionMatchExtension(Simple());
}

// Test the seed extension algorithm with ungapped x-drop extension for simple seeds.
SEQAN_DEFINE_TEST(test_seeds_extension_ungapped_xdrop_extension_simple)
{
    using namespace seqan;
    testSeedsExtensionUnGappedXDropExtension(Simple());
}

// Test the seed extension algorithm with gapped x-drop extension for simple seeds.
SEQAN_DEFINE_TEST(test_seeds_extension_gapped_xdrop_extension_simple)
{
    using namespace seqan;
    testSeedsExtensionGappedXDropExtension(Simple());
}

// Test the seed extension algorithm with ungapped x-drp extension for chained seeds.
SEQAN_DEFINE_TEST(test_seeds_extension_match_extension_chained)
{
    using namespace seqan;
    testSeedsExtensionMatchExtension(ChainedSeed());
}

// Test the seed extension algorithm with ungapped x-drop extension for simple seeds.
SEQAN_DEFINE_TEST(test_seeds_extension_ungapped_xdrop_extension_chained)
{
    using namespace seqan;
    testSeedsExtensionUnGappedXDropExtension(ChainedSeed());
}

// Test the seed extension algorithm with gapped x-drop extension for simple seeds.
SEQAN_DEFINE_TEST(test_seeds_extension_gapped_xdrop_extension_chained)
{
    using namespace seqan;
    testSeedsExtensionGappedXDropExtension(ChainedSeed());
}

#endif  // TEST_SEEDS_TEST_SEEDS_EXTENSION_H_
