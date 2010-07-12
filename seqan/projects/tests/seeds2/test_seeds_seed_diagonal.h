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
  Test the header seeds_seed_diagonal.h containing the SeedDiagonal class.
  ===========================================================================
*/

#ifndef TEST_SEEDS_TEST_SEEDS_SEED_DIAGONAL_H_
#define TEST_SEEDS_TEST_SEEDS_SEED_DIAGONAL_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.

// Test constructors of the SeedDiagonal class.  This also tests that
// all member variables are there.
SEQAN_DEFINE_TEST(test_seeds_seed_diagonal_constructors)
{
    using namespace seqan;

    typedef SeedDiagonal<int, int> TSeedDiagonal;

    // Default constructor.
    {
        TSeedDiagonal sd;
        SEQAN_ASSERT_EQ(0, sd.leftDim0);
        SEQAN_ASSERT_EQ(0, sd.leftDim1);
        SEQAN_ASSERT_EQ(0, sd.length);
    }
    // Only other constructor has all properties.
    {
        TSeedDiagonal sd(1, 2, 3);
        SEQAN_ASSERT_EQ(1, sd.leftDim0);
        SEQAN_ASSERT_EQ(2, sd.leftDim1);
        SEQAN_ASSERT_EQ(3, sd.length);
    }
}


// Test metafunctions of the SeedDiagonal class.
SEQAN_DEFINE_TEST(test_seeds_seed_diagonal_metafunctions)
{
    using namespace seqan;

    // Test the parametrization expected to be used most.
    {
        typedef SeedDiagonal<size_t, size_t> TSeedDiagonal;
        typedef Position<TSeedDiagonal>::Type TPosition;
        bool b1 = TYPECMP<size_t, TPosition>::VALUE;
        SEQAN_ASSERT_TRUE(b1);
        typedef Size<TSeedDiagonal>::Type TSize;
        bool b2 = TYPECMP<size_t, TSize>::VALUE;
        SEQAN_ASSERT_TRUE(b2);
    }
    // Test another parametrization.
    {
        typedef SeedDiagonal<double, int> TSeedDiagonal;
        typedef Position<TSeedDiagonal>::Type TPosition;
        bool b1 = TYPECMP<double, TPosition>::VALUE;
        SEQAN_ASSERT_TRUE(b1);
        typedef Size<TSeedDiagonal>::Type TSize;
        bool b2 = TYPECMP<int, TSize>::VALUE;
        SEQAN_ASSERT_TRUE(b2);
    }
}

#endif  // TEST_SEEDS_TEST_SEEDS_SEED_DIAGONAL_H_
