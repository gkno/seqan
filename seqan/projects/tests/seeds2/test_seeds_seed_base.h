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
  Test the interface of the Seed class for specializations Simple Seed and
  Chained Seed.
  ===========================================================================
*/

#ifndef TEST_SEEDS_TEST_SEEDS_SEED_BASE_H_
#define TEST_SEEDS_TEST_SEEDS_SEED_BASE_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.

template <typename TSeedSpec>
void testSeedsSeedBaseConstructors(TSeedSpec const &)
{
    using namespace seqan;

    // Execute default constructor.
    {
        Seed<int, TSeedSpec> s;
    }
    // Execute with start position and length.
    {
        Seed<int, TSeedSpec> s(1, 2, 3);
    }
    // Execute with start and end position.
    {
        Seed<int, TSeedSpec> s(1, 2, 3, 4);
    }
}

template <typename TSeedSpec>
void testSeedsSeedBaseGettersSetters(TSeedSpec const &)
{
    using namespace seqan;

    // Define Seed type and declare a variable.
    typedef Seed<int, TSeedSpec> TSeed;
    TSeed s(1, 2, 3);

    // Check values from default construction.
    SEQAN_ASSERT_EQ(1, startDiagonal(s));
    SEQAN_ASSERT_EQ(2, endDiagonal(s));
    SEQAN_ASSERT_EQ(2, leftPosition(s, 0));
    SEQAN_ASSERT_EQ(2, rightPosition(s, 0));
    SEQAN_ASSERT_EQ(2, leftPosition(s, 1));
    SEQAN_ASSERT_EQ(2, rightPosition(s, 1));
    SEQAN_ASSERT_EQ(2u, dimension(s));
    SEQAN_ASSERT_EQ(2, leftDim0(s));
    SEQAN_ASSERT_EQ(2, rightDim0(s));
    SEQAN_ASSERT_EQ(2, leftDim1(s));
    SEQAN_ASSERT_EQ(2, rightDim1(s));
    SEQAN_ASSERT_EQ(2, leftDiagonal(s));
    SEQAN_ASSERT_EQ(2, rightDiagonal(s));
}

// Test constructors of the SimpleSeed specialization, as specified
// for the base class Seed.
SEQAN_DEFINE_TEST(test_seeds_seed_base_constructors_simple)
{
    using namespace seqan;
    testSeedsSeedBaseConstructors(Simple());
}

// Test the metafunctions of the SimpleSeed specialization, as
// specified for the base class Seed (none).
SEQAN_DEFINE_TEST(test_seeds_seed_base_metafunctions_simple)
{
    using namespace seqan;
    // No metafunctions in base, intentionally left blank.
}

// Test the getters and seeters of the SimpleSeed specialization as
// specified for the base class Seed.
SEQAN_DEFINE_TEST(test_seeds_seed_base_getters_setters_simple)
{
    using namespace seqan;
    testSeedsSeedBaseGettersSetters(Simple());
}

// Test constructors of the ChainedSeed specialization, as specified
// for the base class Seed.
SEQAN_DEFINE_TEST(test_seeds_seed_base_constructors_chained)
{
    using namespace seqan;
    testSeedsSeedBaseConstructors(ChainedSeed());
}

// Test the metafunctions of the ChainedSeed specialization, as
// specified for the base class Seed (none).
SEQAN_DEFINE_TEST(test_seeds_seed_base_metafunctions_chained)
{
    using namespace seqan;
    // No metafunctions in base, intentionally left blank.
}

// Test the getters and seeters of the ChainedSeed specialization as
// specified for the base class Seed.
SEQAN_DEFINE_TEST(test_seeds_seed_base_getters_setters_chained)
{
    using namespace seqan;
    testSeedsSeedBaseGettersSetters(ChainedSeed());
}

#endif  // TEST_SEEDS_TEST_SEEDS_SEED_BASE_H_
