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
  Test the specialization Unordered SeedSet.
  ===========================================================================
*/

#ifndef TEST_SEEDS_TEST_SEEDS_SEED_SET_UNORDERED_H_
#define TEST_SEEDS_TEST_SEEDS_SEED_SET_UNORDERED_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.

template <typename TSeedSpec>
void testSeedsSeedSetConstructors(TSeedSpec const &)

// Test container functions for specialization unordered.
SEQAN_DEFINE_TEST(test_seeds_seed_set_container_functions_unordered)
{
    using namespace seqan;

    { // Construct with begin/end in both dimensions.
        // Define Seed type and declare a variable.
        typedef Seed<Simple> TSeed;
        TSeed s(1, 2, 3, 5);

        // Check values from construction.
        SEQAN_ASSERT_EQ(1u, getBeginDim0(s));
        SEQAN_ASSERT_EQ(2u, getBeginDim1(s));
        SEQAN_ASSERT_EQ(3u, getEndDim0(s));
        SEQAN_ASSERT_EQ(5u, getEndDim1(s));
        SEQAN_ASSERT_EQ(1, getLowerDiagonal(s));
        SEQAN_ASSERT_EQ(2, getUpperDiagonal(s));
        SEQAN_ASSERT_EQ(1, getStartDiagonal(s));
        SEQAN_ASSERT_EQ(2, getEndDiagonal(s));
    }
    { // Construct from ChainedSeed object.
        typedef Seed<ChainedSeed> TSeed2;
        TSeed2 s2(1, 2, 3);
        typedef Seed<Simple> TSeed;
        TSeed s(s2);

        // Check values from construction.
        SEQAN_ASSERT_EQ(1u, getBeginDim0(s));
        SEQAN_ASSERT_EQ(4u, getEndDim0(s));
        SEQAN_ASSERT_EQ(2u, getBeginDim1(s));
        SEQAN_ASSERT_EQ(5u, getEndDim1(s));
        SEQAN_ASSERT_EQ(1, getLowerDiagonal(s));
        SEQAN_ASSERT_EQ(1, getUpperDiagonal(s));
        SEQAN_ASSERT_EQ(1, getStartDiagonal(s));
        SEQAN_ASSERT_EQ(1, getEndDiagonal(s));
    }
}

#endif  // TEST_SEEDS_TEST_SEEDS_SEED_SET_UNORDERED_H_
