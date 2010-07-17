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
  Test the class SeedSet.
  ===========================================================================
*/

#ifndef TEST_SEEDS_TEST_SEEDS_SEED_SET_BASE_H_
#define TEST_SEEDS_TEST_SEEDS_SEED_SET_BASE_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.

// Test the container functions for the given SeedSet specialization.
template <typename TSeedSetSpec>
void testSeedsSeedSetContainerFunctions(TSeedSetSpec const &)
{
    using namespace seqan;

    // Define SeedSet type and declare a variable.
    typedef SeedSet<Simple, TSeedSetSpec> TSeedSet;
    TSeedSet s;

    // Test length/begin/end with empty set.
    SEQAN_ASSERT_TRUE(begin(s) == end(s));
    SEQAN_ASSERT_EQ(0u, length(s));
    {  // Test with const empty set.
        TSeedSet const & cs = s;
        SEQAN_ASSERT_TRUE(begin(cs) == end(cs));
        SEQAN_ASSERT_EQ(0u, length(cs));
    }

    // Insert one element, test basic accessor functions.
    typedef typename Value<TSeedSet>::Type TSeed;
    addSeed(s, TSeed(1, 2, 3), Single());

    // Test length/begin/end/front/back.
    SEQAN_ASSERT_EQ(1u, length(s));
    SEQAN_ASSERT_TRUE(begin(s) + 1 == end(s));
    SEQAN_ASSERT_TRUE(TSeed(1, 2, 3) == value(begin(s)));
    SEQAN_ASSERT_TRUE(TSeed(1, 2, 3) == front(s));
    SEQAN_ASSERT_TRUE(TSeed(1, 2, 3) == back(s));
    {  // Same tests with const seed set.
        TSeedSet const & cs = s;
        SEQAN_ASSERT_TRUE(begin(cs) + 1 == end(cs));
        SEQAN_ASSERT_TRUE(TSeed(1, 2, 3) == value(begin(cs)));
        SEQAN_ASSERT_TRUE(TSeed(1, 2, 3) == front(cs));
        SEQAN_ASSERT_TRUE(TSeed(1, 2, 3) == back(cs));
    }
}

// Test addSeed(..., Single) for the given SeedSet specialization.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSingle(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec> TSeedSet;
    TSeedSet s;

    // Test length/begin/end with empty set.
    SEQAN_ASSERT_TRUE(begin(s) == end(s));
    SEQAN_ASSERT_EQ(0u, length(s));
    {  // Test with const empty set.
        TSeedSet const & cs = s;
        SEQAN_ASSERT_TRUE(begin(cs) == end(cs));
        SEQAN_ASSERT_EQ(0u, length(cs));
    }

    // Insert one element, test basic accessor functions.
    typedef typename Value<TSeedSet>::Type TSeed;
    addSeed(s, TSeed(1, 2, 3), Single());

    // Test length/begin/end/front/back.
    SEQAN_ASSERT_EQ(1u, length(s));
    SEQAN_ASSERT_TRUE(begin(s) + 1 == end(s));
    SEQAN_ASSERT_TRUE(TSeed(1, 2, 3) == value(begin(s)));
    SEQAN_ASSERT_TRUE(TSeed(1, 2, 3) == front(s));
    SEQAN_ASSERT_TRUE(TSeed(1, 2, 3) == back(s));
    {  // Same tests with const seed set.
        TSeedSet const & cs = s;
        SEQAN_ASSERT_TRUE(begin(cs) + 1 == end(cs));
        SEQAN_ASSERT_TRUE(TSeed(1, 2, 3) == value(begin(cs)));
        SEQAN_ASSERT_TRUE(TSeed(1, 2, 3) == front(cs));
        SEQAN_ASSERT_TRUE(TSeed(1, 2, 3) == back(cs));
    }
}

// Test container functions for specialization Unordered SeedSet.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_functions_unordered)
{
    using namespace seqan;
    testSeedsSeedSetContainerFunctions(Unordered());
}

// Test addSeed(..., Single) for specialization Unordered SeedSet.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_single_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSingle(Unordered());
}

#endif  // TEST_SEEDS_TEST_SEEDS_SEED_SET_BASE_H_
