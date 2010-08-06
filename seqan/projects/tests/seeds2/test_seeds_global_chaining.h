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
  Tests for the header seeds_global_chaining.h.
  ==========================================================================*/

#ifndef TEST_SEEDS_TEST_SEEDS_GLOBAL_CHAINING_H_
#define TEST_SEEDS_TEST_SEEDS_GLOBAL_CHAINING_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.


// Test global chaining weighting the seeds by their length only.
SEQAN_DEFINE_TEST(test_seeds_global_chaining_sparse_length)
{
    using namespace seqan;

    typedef SeedSet<Simple, Unordered, DefaultSeedSetConfig> TSeedSet;
    typedef Value<TSeedSet>::Type TSeed;
    typedef String<TSeed> TSeedChain;

    // Test with one seed.
    {
        TSeedSet seedSet;
        addSeed(seedSet, TSeed(1, 2, 3), Single());

        TSeedChain result;
        chainSeedsGlobally(result, seedSet, SparseChaining());

        SEQAN_ASSERT_EQ(1u, length(result));
        SEQAN_ASSERT_EQ(TSeed(1, 2, 3), front(result));
    }
    // Test with two seeds, both are part of the chain.
    {
        TSeedSet seedSet;
        addSeed(seedSet, TSeed(1, 2, 3), Single());
        addSeed(seedSet, TSeed(4, 5, 6), Single());

        TSeedChain result;
        chainSeedsGlobally(result, seedSet, SparseChaining());

        SEQAN_ASSERT_EQ(2u, length(result));
        SEQAN_ASSERT_EQ(TSeed(1, 2, 3), result[0]);
        SEQAN_ASSERT_EQ(TSeed(4, 5, 6), result[1]);
    }
    // Test with two seeds, only first one is part of the chain.
    {
        TSeedSet seedSet;
        addSeed(seedSet, TSeed(1, 2, 3), Single());
        addSeed(seedSet, TSeed(2, 1, 2), Single());

        TSeedChain result;
        chainSeedsGlobally(result, seedSet, SparseChaining());

        SEQAN_ASSERT_EQ(1u, length(result));
        SEQAN_ASSERT_EQ(TSeed(1, 2, 3), front(result));
    }
    // A bit larger example.
    {
        TSeedSet seedSet;
        addSeed(seedSet, TSeed(0, 0, 2), Single());
        addSeed(seedSet, TSeed(3, 5, 2), Single());
        addSeed(seedSet, TSeed(4, 1, 3), Single());
        addSeed(seedSet, TSeed(9, 9, 2), Single());

        TSeedChain result;
        chainSeedsGlobally(result, seedSet, SparseChaining());

        SEQAN_ASSERT_EQ(3u, length(result));
        SEQAN_ASSERT_EQ(TSeed(0, 0, 2), result[0]);
        SEQAN_ASSERT_EQ(TSeed(4, 1, 3), result[1]);
        SEQAN_ASSERT_EQ(TSeed(9, 9, 2), result[2]);
    }
    // Another non-trivial example that caused a problem earlier.  The
    // seeds are all overlapping here.
    {
        TSeedSet seedSet;
        addSeed(seedSet, TSeed(0, 93, 281, 342), Single());
        addSeed(seedSet, TSeed(3, 237, 127, 364), Single());
        addSeed(seedSet, TSeed(3, 284, 86, 368), Single());
        addSeed(seedSet, TSeed(5, 146, 239, 374), Single());

        TSeedChain result;
        chainSeedsGlobally(result, seedSet, SparseChaining());

        SEQAN_ASSERT_EQ(1u, length(result));
        SEQAN_ASSERT_EQ(TSeed(0, 93, 281, 342), result[0]);
    }
}

#endif  // TEST_SEEDS_TEST_SEEDS_GLOBAL_CHAINING_H_

