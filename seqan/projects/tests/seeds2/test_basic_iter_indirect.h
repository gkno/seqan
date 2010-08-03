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

#ifndef TEST_SEEDS_TEST_BASIC_ITER_INDIRECT_H_
#define TEST_SEEDS_TEST_BASIC_ITER_INDIRECT_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.

#include <set>

// Test constructors of the indirect iterator.
SEQAN_DEFINE_TEST(test_seeds_basic_iter_indirect_constructors)
{
    using namespace seqan;

    SEQAN_ASSERT_FAIL("Write me!");

    // Default constructor.
    {
    }
    // Construct from wrapped iterator.
    {
    }
    // Copy constructor.
    {
    }
}

// Test the metafunctions.
SEQAN_DEFINE_TEST(test_seeds_basic_iter_indirect_metafunctions)
{
    using namespace seqan;

    SEQAN_ASSERT_FAIL("Write me!");

    // Iterator
    {
    }

    // Const Iterator
    {
    }
}

// Test the common iterator functions.
SEQAN_DEFINE_TEST(test_seeds_basic_iter_indirect_basic_functions)
{
    SEQAN_ASSERT_FAIL("Write me!");
}

#endif  // #ifndef TEST_SEEDS_TEST_BASIC_ITER_INDIRECT_H_
