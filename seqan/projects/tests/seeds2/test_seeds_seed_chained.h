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

struct TestSmallSeedConfig
{
    typedef unsigned TPosition;
    typedef unsigned TSize;
    typedef int TDiagonal;
    typedef seqan::True THasScore;
    typedef int TScoreValue;
    typedef seqan::_ScoreMixin<int> TScoreMixin;
};

// Test the metafunctions of the ChainedSeed specialization.
SEQAN_DEFINE_TEST(test_seeds_seed_chained_metafunctions)
{
    using namespace seqan;

    // Test with the default configuration.
    {
        typedef Seed<ChainedSeed> TSeed;
        typedef Value<TSeed>::Type TSeedDiagonal;
        bool b;
        b = TYPECMP<size_t, Position<TSeedDiagonal>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = TYPECMP<size_t, Size<TSeedDiagonal>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);

        b = TYPECMP<size_t, Position<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = TYPECMP<size_t, Size<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = TYPECMP<_MakeSigned<size_t>::Type, Diagonal<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = TYPECMP<False, HasScore<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = TYPECMP<Nothing, SeedScore<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
    }
    // Test with other specialization.
    {
        typedef Seed<ChainedSeed, TestSmallSeedConfig> TSeed;
        typedef Value<TSeed>::Type TSeedDiagonal;
        bool b;
        b = TYPECMP<unsigned, Position<TSeedDiagonal>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = TYPECMP<unsigned, Size<TSeedDiagonal>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);

        b = TYPECMP<unsigned, Position<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = TYPECMP<unsigned, Size<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = TYPECMP<int, Diagonal<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = TYPECMP<True, HasScore<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = TYPECMP<int, SeedScore<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
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
    using namespace seqan;

    typedef Seed<ChainedSeed> TSeed;
    typedef Value<TSeed>::Type TSeedDiagonal;

    TSeed s(1, 2, 3);
    appendDiagonal(s, TSeedDiagonal(4, 5, 3));

    {  // non-const seed
        typedef Iterator<TSeed, Standard>::Type TIterator;
        TIterator it = begin(s);
        SEQAN_ASSERT_EQ(1u, it->leftDim0);
        SEQAN_ASSERT_EQ(2u, it->leftDim1);
        SEQAN_ASSERT_EQ(3u, it->length);
        ++it;
        SEQAN_ASSERT_EQ(4u, it->leftDim0);
        SEQAN_ASSERT_EQ(5u, it->leftDim1);
        SEQAN_ASSERT_EQ(3u, it->length);
        ++it;
        SEQAN_ASSERT_TRUE(it == end(s));
    }
    {  // const seed
        TSeed const & cs = s;
        typedef Iterator<TSeed const, Standard>::Type TIterator;
        TIterator it = begin(cs);
        SEQAN_ASSERT_EQ(1u, it->leftDim0);
        SEQAN_ASSERT_EQ(2u, it->leftDim1);
        SEQAN_ASSERT_EQ(3u, it->length);
        ++it;
        SEQAN_ASSERT_EQ(4u, it->leftDim0);
        SEQAN_ASSERT_EQ(5u, it->leftDim1);
        SEQAN_ASSERT_EQ(3u, it->length);
        ++it;
        SEQAN_ASSERT_TRUE(it == end(cs));
    }
}

#endif  // TEST_SEEDS_TEST_SEEDS_SEED_CHAINED_H_
