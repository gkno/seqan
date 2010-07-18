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


// Test addSeed(..., Merge) with the given SeedSet Specialization.
// Case: Seed left of added; Merging possible; Quality threshold
// reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdReached(TSeedSetSpec const &)
{
    using namespace seqan;
    SEQAN_ASSERT_FAIL("Write me!");
}


// Test addSeed(..., Merge) with the given SeedSet Specialization.
// Case: Seed right of added; Merging possible; Quality threshold
// reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeRightMergingPossibleThresholdReached(TSeedSetSpec const &)
{
    using namespace seqan;
    SEQAN_ASSERT_FAIL("Write me!");
}


// Test addSeed(..., Merge) with the given SeedSet Specialization.
// Case: Seed left of added; Merging impossible; Quality threshold
// reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeLeftMergingImpossibleThresholdReached(TSeedSetSpec const &)
{
    using namespace seqan;
    SEQAN_ASSERT_FAIL("Write me!");
}


// Test addSeed(..., Merge) with the given SeedSet Specialization.
// Case: Seed left of added; Merging possible; Quality threshold
// not reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdNotReached(TSeedSetSpec const &)
{
    using namespace seqan;
    SEQAN_ASSERT_FAIL("Write me!");
}


// Test addSeed(..., Merge) with the given SeedSet Specialization.
// Seeds have scores.  Case: Seed left of added; Merging possible;
// Quality threshold not reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdNotReachedScored(TSeedSetSpec const & )
{
    using namespace seqan;
    SEQAN_ASSERT_FAIL("Write me!");
}


// Test addSeed(..., Merge) with the given SeedSet Specialization.
// Seeds have scores.  Case: Seed left of added; Merging possible;
// Quality threshold reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdReachedScored(TSeedSetSpec const & )
{
    using namespace seqan;
    SEQAN_ASSERT_FAIL("Write me!");
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


// Test addSeed(..., Merge) for specialization Unordered SeedSet.
//
// Case: Seed in set is left of added;  Merging is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_merge_left_merging_possible_threshold_reached_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdReached(Unordered());
}


// Test addSeed(..., Merge) for specialization Unordered SeedSet.
//
// Case: Seed in set is right of added;  Merging is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_merge_right_merging_possible_threshold_reached_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeRightMergingPossibleThresholdReached(Unordered());
}


// Test addSeed(..., Merge) for specialization Unordered SeedSet.
//
// Case: Seed in set is left of added;  Merging is not possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_merge_left_merging_impossible_threshold_reached_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingImpossibleThresholdReached(Unordered());
}


// Test addSeed(..., Merge) for specialization Unordered SeedSet.
//
// Case: Seed in set is left of added;  Merging not possible;  Quality threshold not reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_merge_left_merging_possible_threshold_not_reached_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdNotReached(Unordered());
}


// Test addSeed(..., Merge) for specialization Unordered SeedSet.
// Seeds have scores.
//
// Case: Seed in set is left of added;  Merging is possible;  Quality threshold not reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_merge_left_merging_possible_threshold_not_reached_scored_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdNotReachedScored(Unordered());
}


// Test addSeed(..., Merge) for specialization Unordered SeedSet.
// Seeds have scores.
//
// Case: Seed in set is left of added;  Merging is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_merge_left_merging_possible_threshold_reached_scored_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdReachedScored(Unordered());
}

#endif  // TEST_SEEDS_TEST_SEEDS_SEED_SET_BASE_H_
