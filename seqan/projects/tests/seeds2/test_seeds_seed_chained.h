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
  Test the specialization Chained Seed.
  ===========================================================================
*/

#ifndef TEST_SEEDS_TEST_SEEDS_SEED_CHAINED_H_
#define TEST_SEEDS_TEST_SEEDS_SEED_CHAINED_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.

// Test the metafunctions of the ChainedSeed specialization.
SEQAN_DEFINE_TEST(test_seeds_seed_chained_metafunctions)
{
    using namespace seqan;

    // Test with the specialization expected in SeqAn.
    {
        typedef Seed<ChainedSeed> TSeed;
        typedef Value<TSeed>::Type TSeedDiagonal;
        bool b1 = TYPECMP<size_t, Position<TSeedDiagonal>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b1);
        bool b2 = TYPECMP<size_t, Size<TSeedDiagonal>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b2);
    }
    // Test with other specialization.
    {
        typedef Seed<ChainedSeed> TSeed;
        typedef Value<TSeed>::Type TSeedDiagonal;
        bool b1 = TYPECMP<int, Position<TSeedDiagonal>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b1);
        bool b2 = TYPECMP<int, Size<TSeedDiagonal>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b2);
    }
}

// Test the appendDiagonal() function for Chained Seed.
SEQAN_DEFINE_TEST(test_seeds_seed_chained_add_diagonal)
{
    using namespace seqan;

    typedef Seed<ChainedSeed> TSeed;
    typedef Value<TSeed>::Type TSeedDiagonal;
    TSeed s(1, 3, 4);

    SEQAN_ASSERT_EQ(1u, getLeftDim0(s));
    SEQAN_ASSERT_EQ(3u, getLeftDim1(s));
    SEQAN_ASSERT_EQ(4u, getRightDim0(s));
    SEQAN_ASSERT_EQ(6u, getRightDim1(s));
    SEQAN_ASSERT_EQ(2, getStartDiagonal(s));
    SEQAN_ASSERT_EQ(2, getEndDiagonal(s));
    SEQAN_ASSERT_EQ(1u, length(s._seedDiagonals));

    appendDiagonal(s, TSeedDiagonal(5, 5, 3));

    SEQAN_ASSERT_EQ(1u, getLeftDim0(s));
    SEQAN_ASSERT_EQ(3u, getLeftDim1(s));
    SEQAN_ASSERT_EQ(7u, getRightDim0(s));
    SEQAN_ASSERT_EQ(7u, getRightDim1(s));
    SEQAN_ASSERT_EQ(2, getStartDiagonal(s));
    SEQAN_ASSERT_EQ(0, getEndDiagonal(s));
    SEQAN_ASSERT_EQ(2u, length(s._seedDiagonals));
}

// Test the begin/end functions for chained seeds.
SEQAN_DEFINE_TEST(test_seeds_seed_chained_iterators)
{
    SEQAN_ASSERT_FAIL("Write me!");
}

#endif  // TEST_SEEDS_TEST_SEEDS_SEED_CHAINED_H_
