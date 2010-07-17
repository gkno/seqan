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
        Seed<TSeedSpec> s;
    }
    // Execute with start position and length.
    {
        Seed<TSeedSpec> s(1, 2, 3);
    }
}

template <typename TSeedSpec>
void testSeedsSeedBaseGettersSetters(TSeedSpec const &)
{
    using namespace seqan;

    // Define Seed type and declare a variable.
    typedef Seed<TSeedSpec> TSeed;
    TSeed s(1, 2, 3);

    // Check values from construction.
    SEQAN_ASSERT_EQ(1u, getBeginDim0(s));
    SEQAN_ASSERT_EQ(4u, getEndDim0(s));
    SEQAN_ASSERT_EQ(2u, getBeginDim1(s));
    SEQAN_ASSERT_EQ(5u, getEndDim1(s));
    SEQAN_ASSERT_EQ(1, getLowerDiagonal(s));
    SEQAN_ASSERT_EQ(1, getUpperDiagonal(s));
    SEQAN_ASSERT_EQ(1, getStartDiagonal(s));
    SEQAN_ASSERT_EQ(1, getEndDiagonal(s));

    // Use setters from base class.
    setLowerDiagonal(s, 42);
    SEQAN_ASSERT_EQ(42, getLowerDiagonal(s));
    setUpperDiagonal(s, 5);
    SEQAN_ASSERT_EQ(5, getUpperDiagonal(s));
}

template <typename TSeedSpec>
void testSeedsSeedBaseBasicFunctions(TSeedSpec const &)
{
    using namespace seqan;

    // Define Seed type and declare a variable.
    typedef Seed<TSeedSpec> TSeed;
    TSeed s(1, 2, 3);

    {  // test assign
        TSeed x;
        assign(x, s);
        SEQAN_ASSERT_TRUE(x == s);
    }
    {  // test move
        TSeed s2(s);
        TSeed x;
        move(x, s);
        SEQAN_ASSERT_TRUE(x == s2);
    }
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

// Test the basic functions of the SimpleSeed specialization as
// specified for the base class Seed.
SEQAN_DEFINE_TEST(test_seeds_seed_base_basic_functions_simple)
{
    using namespace seqan;
    testSeedsSeedBaseBasicFunctions(Simple());
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

// Test the basic functions of the ChainedSeed specialization as
// specified for the base class Seed.
SEQAN_DEFINE_TEST(test_seeds_seed_base_basic_functions_chained)
{
    using namespace seqan;
    testSeedsSeedBaseBasicFunctions(ChainedSeed());
}

#endif  // TEST_SEEDS_TEST_SEEDS_SEED_BASE_H_
