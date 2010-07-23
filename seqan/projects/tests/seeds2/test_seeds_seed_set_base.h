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
    {
        typedef typename Iterator<TSeedSet, Rooted>::Type TIterator; // TODO(holtgrew): Why explicit rooted necessary?
        TIterator it(begin(s));
        ++it;
        SEQAN_ASSERT_TRUE(it == end(s));
    }
    SEQAN_ASSERT_TRUE(TSeed(1, 2, 3) == value(begin(s)));
    SEQAN_ASSERT_TRUE(TSeed(1, 2, 3) == front(s));
    SEQAN_ASSERT_TRUE(TSeed(1, 2, 3) == back(s));
    {  // Same tests with const seed set.
        TSeedSet const & cs = s;
        {
            typedef typename Iterator<TSeedSet const, Rooted>::Type TIterator; // TODO(holtgrew): Why explicit rooted necessary?
            TIterator it(begin(cs));
            ++it;
            SEQAN_ASSERT_TRUE(it == end(cs));
        }
        SEQAN_ASSERT_TRUE(TSeed(1, 2, 3) == value(begin(cs)));
        SEQAN_ASSERT_TRUE(TSeed(1, 2, 3) == front(cs));
        SEQAN_ASSERT_TRUE(TSeed(1, 2, 3) == back(cs));
    }
}


// Test addSeed(..., Single) for the given SeedSet specialization.
//
// Case: No quality threshold.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSingleNoThreshold(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec, DefaultSeedSetConfig> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;

    TSeed seed(3, 3, 3);

    addSeed(set, seed, Single());
    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(seed, front(set));
}


// Test addSeed(..., Single) for the given SeedSet specialization.
//
// Case: Seed size threshold reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSingleThresholdReachedLength(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec, DefaultSeedSetConfigLength> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinSeedSizeThreshold(set, 1);

    TSeed seed(3, 3, 3);

    addSeed(set, seed, Single());
    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(seed, front(set));
}


// Test addSeed(..., Single) for the given SeedSet specialization.
//
// Case: Seed size threshold not reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSingleThresholdNotReachedLength(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec, DefaultSeedSetConfigLength> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinSeedSizeThreshold(set, 4);

    TSeed seed(3, 3, 3);

    addSeed(set, seed, Single());
    SEQAN_ASSERT_EQ(0u, length(set));
}


// Test addSeed(..., Single) for the given SeedSet specialization.
//
// Case: Seed score threshold reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSingleThresholdReachedScore(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec, DefaultSeedSetConfigScore> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinScoreThreshold(set, -1);

    TSeed seed(3, 3, 3);
    setScore(seed, 1);

    addSeed(set, seed, Single());
    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(seed, front(set));
}


// Test addSeed(..., Single) for the given SeedSet specialization.
//
// Case: Seed score threshold not reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSingleThresholdNotReachedScore(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec, DefaultSeedSetConfigScore> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinScoreThreshold(set, -1);

    TSeed seed(3, 3, 3);
    setScore(seed, -2);

    addSeed(set, seed, Single());
    SEQAN_ASSERT_EQ(0u, length(set));
}


// Test addSeed(..., Merge) with the given SeedSet Specialization.
// Case: Seed left of added; Merging possible; No quality threshold
// required.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeLeftMergingPossibleNoThreshold(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    addSeed(set, TSeed(3, 3, 3), Single());
    // Add seed with maximal diagonal distance 1, the two seeds are within this distance.
    bool ret = addSeed(set, TSeed(2, 2, 3), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), Merge());
    SEQAN_ASSERT_TRUE(ret);

    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(TSeed(2, 2, 4), front(set));
}


// Test addSeed(..., Merge) with the given SeedSet Specialization.
// Case: Seed right of added; Merging possible; No quality threshold
// required.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeRightMergingPossibleNoThreshold(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    addSeed(set, TSeed(2, 2, 3), Single());
    // Add seed with maximal diagonal distance 1, the two seeds are within this distance.
    bool ret = addSeed(set, TSeed(3, 3, 3), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), Merge());
    SEQAN_ASSERT_TRUE(ret);

    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(TSeed(2, 2, 4), front(set));
}


// Test addSeed(..., Merge) with the given SeedSet Specialization.
// Case: Seed left of added; Merging impossible; No quality threshold
// required.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeLeftMergingImpossibleNoThreshold(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    {  // Merging not possible because of diagonal distance
        TSeedSet set;
        addSeed(set, TSeed(0, 2, 3), Single());
        // Add seed with maximal diagonal distance 1, the two seeds
        // overlap but are are not within this distance.
        bool ret = addSeed(set, TSeed(3, 2, 3), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), Merge());
        SEQAN_ASSERT_NOT(ret);

        SEQAN_ASSERT_EQ(1u, length(set));
        SEQAN_ASSERT_EQ(TSeed(0, 2, 3), front(set));
    }
    {  // Merging not possible because non-overlapping
        TSeedSet set;
        addSeed(set, TSeed(0, 0, 2), Single());
        // Add seed with maximal diagonal distance 1, the two seeds
        // are within this diagonal distance but do not overlap.
        bool ret = addSeed(set, TSeed(3, 2, 3), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), Merge());
        SEQAN_ASSERT_NOT(ret);

        SEQAN_ASSERT_EQ(1u, length(set));
        SEQAN_ASSERT_EQ(TSeed(0, 0, 2), front(set));
    }
}


// Test addSeed(..., Merge) with the given SeedSet Specialization.
// Case: Seed left of added; Merging possible; Length quality
// threshold not reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdNotReachedLength(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec, DefaultSeedSetConfigLength> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinSeedSizeThreshold(set, 5);

    // Add a low-scoring seed.
    addSeed(set, TSeed(1, 1, 3), Single());
    // The seed is added to the set but not included in the set of
    // seeds above the quality threshold.
    SEQAN_ASSERT_EQ(0u, length(set));

    // Add a seed to merge into a low-quality seed.
    bool ret = addSeed(set, TSeed(0, 0, 2), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), Merge());
    SEQAN_ASSERT_TRUE(ret);
    // The seed is merged into the first one but the first one does
    // not exceed the quality threshold yet.
    SEQAN_ASSERT_EQ(0u, length(set));
}


// Test addSeed(..., Merge) with the given SeedSet Specialization.
// Case: Seed left of added; Merging possible; Length quality
// threshold reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdReachedLength(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec, DefaultSeedSetConfigLength> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinSeedSizeThreshold(set, 4);

    // Add a low-quality seed.
    addSeed(set, TSeed(1, 1, 3), Single());
    // The seed is added to the set but not included in the set of
    // seeds above the quality threshold.
    SEQAN_ASSERT_EQ(0u, length(set));

    // Add a seed to merge into a high-quality seed.
    bool ret = addSeed(set, TSeed(0, 0, 2), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), Merge());
    SEQAN_ASSERT_TRUE(ret);
    // The seed is merged into the first one and exceeds the quality
    // threshold.
    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(TSeed(0, 0, 4), front(set));
}


// Test addSeed(..., Merge) with the given SeedSet Specialization.
// Seeds have scores.  Case: Seed left of added; Merging possible;
// Score quality threshold not reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdNotReachedScored(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec, DefaultSeedSetConfigScore> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinScoreThreshold(set, -1);

    // Add a low-quality seed.
    TSeed s1(1, 1, 3);
    setScore(s1, -2);
    addSeed(set, s1, Single());
    SEQAN_ASSERT_EQ(0u, length(set));

    // Add a seed to merge into another low-quality seed.
    TSeed s2(0, 0, 2);
    setScore(s2, -1);
    bool ret = addSeed(set, s2, 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), Merge());
    SEQAN_ASSERT_TRUE(ret);
    // The seed is merged with the existing low-quality seed into one
    // of low quality.
    SEQAN_ASSERT_EQ(0u, length(set));
}


// Test addSeed(..., Merge) with the given SeedSet Specialization.
// Seeds have scores.  Case: Seed left of added; Merging possible;
// Score quality threshold reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdReachedScored(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec, DefaultSeedSetConfigScore> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinScoreThreshold(set, -2);

    // Add a low-quality seed.
    TSeed s1(1, 1, 3);
    setScore(s1, -3);
    addSeed(set, s1, Single());
    SEQAN_ASSERT_EQ(0u, length(set));

    // Add a seed to merge into a high-quality seed.
    TSeed s2(0, 0, 1);
    setScore(s2, 1);
    bool ret = addSeed(set, s2, 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), Merge());
    SEQAN_ASSERT_TRUE(ret);
    // The seed is merged with the existing low-quality seed into one
    // of sufficiently high quality.
    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(TSeed(0, 0, 4), front(set));
    SEQAN_ASSERT_EQ(-2, getScore(front(set)));
}


// Test addSeed(..., SimpleChain) with the given SeedSet Specialization.
// Case: Seed left of added; Chaining possible; No quality threshold
// required.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleNoThreshold(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    addSeed(set, TSeed(4, 5, 3), Single());
    // Add seed with maximal distance 2, the two seeds are within this distance.
    bool ret = addSeed(set, TSeed(1, 1, 3), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), SimpleChain());
    SEQAN_ASSERT_TRUE(ret);

    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(1u, getBeginDim0(front(set)));
    SEQAN_ASSERT_EQ(7u, getEndDim0(front(set)));
    SEQAN_ASSERT_EQ(1u, getBeginDim1(front(set)));
    SEQAN_ASSERT_EQ(8u, getEndDim1(front(set)));
}


// Test addSeed(..., SimpleChain) with the given SeedSet Specialization.
// Case: Seed right of added; Chaining possible; No quality threshold
// required.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSimpleChainRightChainingPossibleNoThreshold(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    addSeed(set, TSeed(1, 1, 3), Single());
    // Add seed with maximal distance 2, the two seeds are within this distance.
    bool ret = addSeed(set, TSeed(4, 5, 3), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), SimpleChain());
    SEQAN_ASSERT_TRUE(ret);

    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(1u, getBeginDim0(front(set)));
    SEQAN_ASSERT_EQ(7u, getEndDim0(front(set)));
    SEQAN_ASSERT_EQ(1u, getBeginDim1(front(set)));
    SEQAN_ASSERT_EQ(8u, getEndDim1(front(set)));
}


// Test addSeed(..., SimpleChain) with the given SeedSet Specialization.
// Case: Seed left of added; Chaining impossible; No quality threshold
// required.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSimpleChainLeftChainingImpossibleNoThreshold(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    {  // Chaining not possible because of distance
        TSeedSet set;
        addSeed(set, TSeed(0, 0, 3), Single());
        // Add seed with maximal distance 1, the two seeds do not
        // overlap but are are not within this distance.
        bool ret = addSeed(set, TSeed(5, 5, 3), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary*/, Nothing(), Nothing(), SimpleChain());
        SEQAN_ASSERT_NOT(ret);

        SEQAN_ASSERT_EQ(1u, length(set));
        SEQAN_ASSERT_EQ(TSeed(0, 0, 3), front(set));
    }
    {  // Chaining not possible because overlapping
        TSeedSet set;
        addSeed(set, TSeed(1, 2, 3), Single());
        // Add seed with maximal diagonal distance 1, the two seeds
        // are within this distance but do not overlap.
        bool ret = addSeed(set, TSeed(0, 0, 3), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary*/, Nothing(), Nothing(), SimpleChain());
        SEQAN_ASSERT_NOT(ret);

        SEQAN_ASSERT_EQ(1u, length(set));
        SEQAN_ASSERT_EQ(TSeed(1, 2, 3), front(set));
    }
}


// Test addSeed(..., SimpleChain) with the given SeedSet Specialization.
// Case: Seed left of added; Chaining possible; Length quality
// threshold not reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdNotReachedLength(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec, DefaultSeedSetConfigLength> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinSeedSizeThreshold(set, 6);

    // Add a low-quality seed.
    addSeed(set, TSeed(2, 2, 3), Single());
    // The seed is added to the set but not included in the set of
    // seeds above the quality threshold.
    SEQAN_ASSERT_EQ(0u, length(set));

    // Add a seed to chain into a low-quality seed.
    bool ret = addSeed(set, TSeed(0, 0, 2), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary*/, Nothing(), Nothing(), SimpleChain());
    SEQAN_ASSERT_TRUE(ret);
    // The seed is chained into the first one but the result does not
    // exceed the quality threshold yet.
    SEQAN_ASSERT_EQ(0u, length(set));
}


// Test addSeed(..., SimpleChain) with the given SeedSet Specialization.
// Case: Seed left of added; Chaining possible; Length quality
// threshold reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdReachedLength(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec, DefaultSeedSetConfigLength> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinSeedSizeThreshold(set, 5);

    // Add a low-quality seed.
    addSeed(set, TSeed(2, 2, 3), Single());
    // The seed is added to the set but not included in the set of
    // seeds above the quality threshold.
    SEQAN_ASSERT_EQ(0u, length(set));

    // Add a seed to chain into a low-quality seed.
    bool ret = addSeed(set, TSeed(0, 0, 2), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary */, Nothing(), Nothing(), SimpleChain());
    SEQAN_ASSERT_TRUE(ret);
    // The seed is chained into the first one and exceeds the quality
    // threshold.
    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(TSeed(0, 0, 5), front(set));
}


// Test addSeed(..., SimpleChain) with the given SeedSet Specialization.
// Seeds have scores.  Case: Seed left of added; Chaining possible;
// Score quality threshold not reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdNotReachedScored(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec, DefaultSeedSetConfigScore> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinScoreThreshold(set, 2);
    Score<int, Simple> scoringScheme(1, -1, -1);

    // Add a low-quality seed.
    TSeed s1(0, 0, 3);
    setScore(s1, 1);
    addSeed(set, s1, Single());
    SEQAN_ASSERT_EQ(0u, length(set));

    // Add a seed to chain into another low-quality seed.
    TSeed s2(4, 4, 2);
    setScore(s2, 1);
    bool ret = addSeed(set, s2, 1, Nothing(), scoringScheme, Nothing(), Nothing(), SimpleChain());
    SEQAN_ASSERT_TRUE(ret);
    // The seed is chained with the existing low-quality seed into one
    // of low quality.
    SEQAN_ASSERT_EQ(0u, length(set));
}


// Test addSeed(..., SimpleChain) with the given SeedSet Specialization.
// Seeds have scores.  Case: Seed left of added; Chaining possible;
// Score quality threshold reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdReachedScored(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec, DefaultSeedSetConfigScore> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinScoreThreshold(set, 3);
    Score<int, Simple> scoringScheme(1, -1, -1);

    // Add a low-quality seed.
    TSeed s1(0, 0, 3);
    setScore(s1, 2);
    addSeed(set, s1, Single());
    SEQAN_ASSERT_EQ(0u, length(set));

    // Add a seed to chain into a seed of sufficiently high quality.
    TSeed s2(4, 4, 2);
    setScore(s2, 2);
    bool ret = addSeed(set, s2, 1, Nothing(), scoringScheme, Nothing(), Nothing(), SimpleChain());
    SEQAN_ASSERT_TRUE(ret);
    // The seed is chained with the existing low-quality seed into one
    // of sufficiently high quality.
    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(TSeed(0, 0, 6), front(set));
    SEQAN_ASSERT_EQ(3, getScore(front(set)));
}


// Test addSeed(..., Chaos) with the given SeedSet Specialization.
// Case: Seed left of added; Chaining possible; No quality threshold
// required.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedChaosLeftChainingPossibleNoThreshold(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    DnaString sequence0 = "CCCCCCCCCC";
    DnaString sequence1 = "CCCCCCCCCC";

    TSeedSet set;
    addSeed(set, TSeed(4, 5, 3), Single());
    // Add seed with maximal distance 1, bandwidth 2, the two seeds are within this distance.
    bool ret = addSeed(set, TSeed(1, 1, 3), 1, 2, Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, sequence0, sequence1, Chaos());
    SEQAN_ASSERT_TRUE(ret);

    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(1u, getBeginDim0(front(set)));
    SEQAN_ASSERT_EQ(7u, getEndDim0(front(set)));
    SEQAN_ASSERT_EQ(1u, getBeginDim1(front(set)));
    SEQAN_ASSERT_EQ(8u, getEndDim1(front(set)));
}


// Test addSeed(..., Chaos) with the given SeedSet Specialization.
// Case: Seed right of added; Chaining possible; No quality threshold
// required.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedChaosRightChainingPossibleNoThreshold(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    DnaString sequence0 = "CCCCCCCCCC";
    DnaString sequence1 = "CCCCCCCCCC";

    TSeedSet set;
    addSeed(set, TSeed(1, 1, 3), Single());
    // Add seed with maximal distance 1, bandwidth 2, the two seeds are within this distance.
    bool ret = addSeed(set, TSeed(4, 5, 3), 1, 2, Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, sequence0, sequence1, Chaos());
    SEQAN_ASSERT_TRUE(ret);

    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(1u, getBeginDim0(front(set)));
    SEQAN_ASSERT_EQ(7u, getEndDim0(front(set)));
    SEQAN_ASSERT_EQ(1u, getBeginDim1(front(set)));
    SEQAN_ASSERT_EQ(8u, getEndDim1(front(set)));
}


// Test addSeed(..., Chaos) with the given SeedSet Specialization.
// Case: Seed left of added; Chaining impossible; No quality threshold
// required.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedChaosLeftChainingImpossibleNoThreshold(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    {  // Chaining not possible because of distance
        DnaString sequence0 = "CCCCCCCCCC";
        DnaString sequence1 = "CCCCCCCCCC";

        TSeedSet set;
        addSeed(set, TSeed(0, 0, 3), Single());
        // Add seed with maximal distance 1, bandwidth 2, the two
        // seeds do not overlap but are are not within this distance.
        bool ret = addSeed(set, TSeed(5, 5, 3), 1, 2, Score<int, Simple>()/*TODO(holtgrew): unnecessary*/, sequence0, sequence1, Chaos());
        SEQAN_ASSERT_NOT(ret);

        SEQAN_ASSERT_EQ(1u, length(set));
        SEQAN_ASSERT_EQ(TSeed(0, 0, 3), front(set));
    }
    {  // Chaining not possible because overlapping
        DnaString sequence0 = "CCCCCCCCCC";
        DnaString sequence1 = "CCCCCCCCCC";

        TSeedSet set;
        addSeed(set, TSeed(1, 2, 3), Single());
        // Add seed with maximal diagonal distance 1, bandwidth 2, the
        // two seeds are within this distance but do not overlap.
        bool ret = addSeed(set, TSeed(0, 0, 3), 1, 2, Score<int, Simple>()/*TODO(holtgrew): unnecessary*/, sequence0, sequence1, Chaos());
        SEQAN_ASSERT_NOT(ret);

        SEQAN_ASSERT_EQ(1u, length(set));
        SEQAN_ASSERT_EQ(TSeed(1, 2, 3), front(set));
    }
}


// Test addSeed(..., Chaos) with the given SeedSet Specialization.
// Case: Seed left of added; Chaining possible; Length quality
// threshold not reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdNotReachedLength(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec, DefaultSeedSetConfigLength> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinSeedSizeThreshold(set, 6);

    DnaString sequence0 = "CCCCCCCCCC";
    DnaString sequence1 = "CCCCCCCCCC";

    // Add a low-quality seed.
    addSeed(set, TSeed(2, 2, 3), Single());
    // The seed is added to the set but not included in the set of
    // seeds above the quality threshold.
    SEQAN_ASSERT_EQ(0u, length(set));

    // Add a seed to chain into a low-quality seed.
    bool ret = addSeed(set, TSeed(0, 0, 2), 1, 2, Score<int, Simple>()/*TODO(holtgrew): unnecessary*/, sequence0, sequence1, Chaos());
    SEQAN_ASSERT_TRUE(ret);
    // The seed is chained into the first one but the result does not
    // exceed the quality threshold yet.
    SEQAN_ASSERT_EQ(0u, length(set));
}


// Test addSeed(..., Chaos) with the given SeedSet Specialization.
// Case: Seed left of added; Chaining possible; Length quality
// threshold reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdReachedLength(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec, DefaultSeedSetConfigLength> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinSeedSizeThreshold(set, 5);

    DnaString sequence0 = "CCCCCCCCCC";
    DnaString sequence1 = "CCCCCCCCCC";

    // Add a low-quality seed.
    addSeed(set, TSeed(2, 2, 3), Single());
    // The seed is added to the set but not included in the set of
    // seeds above the quality threshold.
    SEQAN_ASSERT_EQ(0u, length(set));

    // Add a seed to chain into a low-quality seed.
    bool ret = addSeed(set, TSeed(0, 0, 2), 1, 2, Score<int, Simple>()/*TODO(holtgrew): unnecessary */, sequence0, sequence1, Chaos());
    SEQAN_ASSERT_TRUE(ret);
    // The seed is chained into the first one and exceeds the quality
    // threshold.
    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(TSeed(0, 0, 5), front(set));
}


// Test addSeed(..., Chaos) with the given SeedSet Specialization.
// Seeds have scores.  Case: Seed left of added; Chaining possible;
// Score quality threshold not reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdNotReachedScored(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec, DefaultSeedSetConfigScore> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinScoreThreshold(set, 2);
    Score<int, Simple> scoringScheme(1, -1, -1);

    DnaString sequence0 = "CCCCCCCCCC";
    DnaString sequence1 = "CCCCCCCCCC";

    // Add a low-quality seed.
    TSeed s1(0, 0, 3);
    setScore(s1, 1);
    addSeed(set, s1, Single());
    SEQAN_ASSERT_EQ(0u, length(set));

    // Add a seed to chain into another low-quality seed.
    TSeed s2(4, 4, 2);
    setScore(s2, 1);
    bool ret = addSeed(set, s2, 1, 2, scoringScheme, sequence0, sequence1, Chaos());
    SEQAN_ASSERT_TRUE(ret);
    // The seed is chained with the existing low-quality seed into one
    // of low quality.
    SEQAN_ASSERT_EQ(0u, length(set));
}


// Test addSeed(..., Chaos) with the given SeedSet Specialization.
// Seeds have scores.  Case: Seed left of added; Chaining possible;
// Score quality threshold reached.
template <typename TSeedSetSpec>
void testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdReachedScored(TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Simple, TSeedSetSpec, DefaultSeedSetConfigScore> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinScoreThreshold(set, 3);
    Score<int, Simple> scoringScheme(1, -1, -1);

    DnaString sequence0 = "CCCCCCCCCC";
    DnaString sequence1 = "CCCCCCCCCC";

    // Add a low-quality seed.
    TSeed s1(0, 0, 3);
    setScore(s1, 2);
    addSeed(set, s1, Single());
    SEQAN_ASSERT_EQ(0u, length(set));

    // Add a seed to chain into a seed of sufficiently high quality.
    TSeed s2(4, 4, 2);
    setScore(s2, 2);
    bool ret = addSeed(set, s2, 1, 2, scoringScheme, sequence0, sequence1, Chaos());
    SEQAN_ASSERT_TRUE(ret);
    // The seed is chained with the existing low-quality seed into one
    // of sufficiently high quality.
    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(TSeed(0, 0, 6), front(set));
    SEQAN_ASSERT_EQ(3, getScore(front(set)));
}

// Test container functions for specialization Unordered SeedSet.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_functions_unordered)
{
    using namespace seqan;
    testSeedsSeedSetContainerFunctions(Unordered());
}


// Test addSeed(..., Single) for specialization Unordered SeedSet.
//
// Case: No threshold.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_single_no_threshold_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSingleNoThreshold(Unordered());
}


// Test addSeed(..., Single) for specialization Unordered SeedSet.
//
// Case: Size threshold, threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_single_threshold_reached_length_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSingleThresholdReachedLength(Unordered());
}


// Test addSeed(..., Single) for specialization Unordered SeedSet.
//
// Case: Size threshold, threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_single_threshold_not_reached_length_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSingleThresholdNotReachedLength(Unordered());
}


// Test addSeed(..., Single) for specialization Unordered SeedSet.
//
// Case: Size threshold, threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_single_threshold_reached_score_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSingleThresholdReachedScore(Unordered());
}


// Test addSeed(..., Single) for specialization Unordered SeedSet.
//
// Case: Size threshold, threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_single_threshold_not_reached_score_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSingleThresholdNotReachedScore(Unordered());
}


// Test addSeed(..., Merge) for specialization Unordered SeedSet.
//
// Case: Seed in set is left of added;  Merging is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_merge_left_merging_possible_no_threshold_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingPossibleNoThreshold(Unordered());
}


// Test addSeed(..., Merge) for specialization Unordered SeedSet.
//
// Case: Seed in set is right of added;  Merging is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_merge_right_merging_possible_no_threshold_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeRightMergingPossibleNoThreshold(Unordered());
}


// Test addSeed(..., Merge) for specialization Unordered SeedSet.
//
// Case: Seed in set is left of added;  Merging is not possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_merge_left_merging_impossible_no_threshold_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingImpossibleNoThreshold(Unordered());
}


// Test addSeed(..., Merge) for specialization Unordered SeedSet.
//
// Case: Seed in set is left of added;  Merging not possible;  Length quality threshold not reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_merge_left_merging_possible_threshold_not_reached_length_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdNotReachedLength(Unordered());
}


// Test addSeed(..., Merge) for specialization Unordered SeedSet.
//
// Case: Seed in set is left of added;  Merging not possible;  Length quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_merge_left_merging_possible_threshold_reached_length_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdReachedLength(Unordered());
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

// Test addSeed(..., SimpleChain) for specialization Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_simple_chain_left_chaining_possible_no_threshold_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleNoThreshold(Unordered());
}


// Test addSeed(..., SimpleChain) for specialization Unordered SeedSet.
//
// Case: Seed in set is right of added;  Chaining is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_simple_chain_right_chaining_possible_no_threshold_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainRightChainingPossibleNoThreshold(Unordered());
}


// Test addSeed(..., SimpleChain) for specialization Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining is not possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_simple_chain_left_chaining_impossible_no_threshold_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainLeftChainingImpossibleNoThreshold(Unordered());
}


// Test addSeed(..., SimpleChain) for specialization Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining not possible;  Length quality threshold not reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_simple_chain_left_chaining_possible_threshold_not_reached_length_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdNotReachedLength(Unordered());
}


// Test addSeed(..., SimpleChain) for specialization Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining not possible;  Length quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_simple_chain_left_chaining_possible_threshold_reached_length_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdReachedLength(Unordered());
}


// Test addSeed(..., SimpleChain) for specialization Unordered SeedSet.
// Seeds have scores.
//
// Case: Seed in set is left of added;  Chaining is possible;  Quality threshold not reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_simple_chain_left_chaining_possible_threshold_not_reached_scored_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdNotReachedScored(Unordered());
}


// Test addSeed(..., SimpleChain) for specialization Unordered SeedSet.
// Seeds have scores.
//
// Case: Seed in set is left of added;  Chaining is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_simple_chain_left_chaining_possible_threshold_reached_scored_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdReachedScored(Unordered());
}


// Test addSeed(..., Chaos) for specialization Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_chaos_left_chaining_possible_no_threshold_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosLeftChainingPossibleNoThreshold(Unordered());
}


// Test addSeed(..., Chaos) for specialization Unordered SeedSet.
//
// Case: Seed in set is right of added;  Chaining is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_chaos_right_chaining_possible_no_threshold_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosRightChainingPossibleNoThreshold(Unordered());
}


// Test addSeed(..., Chaos) for specialization Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining is not possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_chaos_left_chaining_impossible_no_threshold_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosLeftChainingImpossibleNoThreshold(Unordered());
}


// Test addSeed(..., Chaos) for specialization Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining not possible;  Length quality threshold not reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_chaos_left_chaining_possible_threshold_not_reached_length_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdNotReachedLength(Unordered());
}


// Test addSeed(..., Chaos) for specialization Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining not possible;  Length quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_chaos_left_chaining_possible_threshold_reached_length_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdReachedLength(Unordered());
}


// Test addSeed(..., Chaos) for specialization Unordered SeedSet.
// Seeds have scores.
//
// Case: Seed in set is left of added;  Chaining is possible;  Quality threshold not reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_chaos_left_chaining_possible_threshold_not_reached_scored_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdNotReachedScored(Unordered());
}


// Test addSeed(..., Chaos) for specialization Unordered SeedSet.
// Seeds have scores.
//
// Case: Seed in set is left of added;  Chaining is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_add_seed_chaos_left_chaining_possible_threshold_reached_scored_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdReachedScored(Unordered());
}

#endif  // TEST_SEEDS_TEST_SEEDS_SEED_SET_BASE_H_
