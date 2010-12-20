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
    typedef seqan::ScoreMixin_<int> TScoreMixin;
};

// Test assignment of chained seeds.
SEQAN_DEFINE_TEST(test_seeds_seed_chained_assign)
{
    using namespace seqan;

    typedef Seed<ChainedSeed> TSeed;

    TSeed s(0, 0, 3);

    TSeed s2 = s;
    s2 = s;
    assign(s2, s);
}

// Test the metafunctions of the ChainedSeed specialization.
SEQAN_DEFINE_TEST(test_seeds_seed_chained_metafunctions)
{
    using namespace seqan;

    // Test with the default configuration.
    {
        typedef Seed<ChainedSeed> TSeed;
        typedef Value<TSeed>::Type TSeedDiagonal;
        bool b;
        b = IsSameType<size_t, Position<TSeedDiagonal>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = IsSameType<size_t, Size<TSeedDiagonal>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);

        b = IsSameType<size_t, Position<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = IsSameType<size_t, Size<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = IsSameType<MakeSigned_<size_t>::Type, Diagonal<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = IsSameType<False, HasScore<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = IsSameType<Nothing, SeedScore<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
    }
    // Test with other specialization.
    {
        typedef Seed<ChainedSeed, TestSmallSeedConfig> TSeed;
        typedef Value<TSeed>::Type TSeedDiagonal;
        bool b;
        b = IsSameType<unsigned, Position<TSeedDiagonal>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = IsSameType<unsigned, Size<TSeedDiagonal>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);

        b = IsSameType<unsigned, Position<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = IsSameType<unsigned, Size<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = IsSameType<int, Diagonal<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = IsSameType<True, HasScore<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        b = IsSameType<int, SeedScore<TSeed>::Type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
    }
}

// Test the front() and back() functions for Chained Seed.
SEQAN_DEFINE_TEST(test_seeds_seed_chained_front_back)
{
    using namespace seqan;
    
    typedef Seed<ChainedSeed> TSeed;
    typedef Value<TSeed>::Type TSeedDiagonal;

    {
        TSeed s(1, 3, 4);
        SEQAN_ASSERT_EQ(TSeedDiagonal(1, 3, 4), front(s));
        SEQAN_ASSERT_EQ(TSeedDiagonal(1, 3, 4), back(s));
    }
    {
        TSeed const cs(1, 3, 4);
        SEQAN_ASSERT_EQ(TSeedDiagonal(1, 3, 4), front(cs));
        SEQAN_ASSERT_EQ(TSeedDiagonal(1, 3, 4), back(cs));
    }
}

// Test the appendDiagonal() function for Chained Seed.
SEQAN_DEFINE_TEST(test_seeds_seed_chained_append_diagonal)
{
    using namespace seqan;

    typedef Seed<ChainedSeed> TSeed;
    typedef Value<TSeed>::Type TSeedDiagonal;
    TSeed s(1, 3, 4);

    SEQAN_ASSERT_EQ(1u, getBeginDim0(s));
    SEQAN_ASSERT_EQ(3u, getBeginDim1(s));
    SEQAN_ASSERT_EQ(5u, getEndDim0(s));
    SEQAN_ASSERT_EQ(7u, getEndDim1(s));
    SEQAN_ASSERT_EQ(2, getStartDiagonal(s));
    SEQAN_ASSERT_EQ(2, getEndDiagonal(s));
    SEQAN_ASSERT_EQ(1u, length(s));

    appendDiagonal(s, TSeedDiagonal(5, 7, 3));

    SEQAN_ASSERT_EQ(1u, getBeginDim0(s));
    SEQAN_ASSERT_EQ(3u, getBeginDim1(s));
    SEQAN_ASSERT_EQ(8u, getEndDim0(s));
    SEQAN_ASSERT_EQ(10u, getEndDim1(s));
    SEQAN_ASSERT_EQ(2, getStartDiagonal(s));
    SEQAN_ASSERT_EQ(2, getEndDiagonal(s));
    SEQAN_ASSERT_EQ(2u, length(s));
}

// Test the truncateDiagonal() function for Chained Seed.
SEQAN_DEFINE_TEST(test_seeds_seed_chained_truncate_diagonals)
{
    using namespace seqan;

    typedef Seed<ChainedSeed> TSeed;
    typedef Value<TSeed>::Type TSeedDiagonal;
    TSeed s(1, 3, 4);

    appendDiagonal(s, TSeedDiagonal(5, 7, 3));
    SEQAN_ASSERT_EQ(2u, length(s));
    appendDiagonal(s, TSeedDiagonal(10, 10, 3));
    SEQAN_ASSERT_EQ(3u, length(s));
    typedef Iterator<TSeed, Standard>::Type TIterator;
    TIterator it = begin(s);
    ++it;
    truncateDiagonals(s, it);
    SEQAN_ASSERT_EQ(1u, length(s));
    TSeedDiagonal const diag = *begin(s, Standard());
    SEQAN_ASSERT_TRUE(diag == TSeedDiagonal(1, 3, 4));
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
        SEQAN_ASSERT_EQ(1u, it->beginDim0);
        SEQAN_ASSERT_EQ(2u, it->beginDim1);
        SEQAN_ASSERT_EQ(3u, it->length);
        ++it;
        SEQAN_ASSERT_EQ(4u, it->beginDim0);
        SEQAN_ASSERT_EQ(5u, it->beginDim1);
        SEQAN_ASSERT_EQ(3u, it->length);
        ++it;
        SEQAN_ASSERT_TRUE(it == end(s));
    }
    {  // const seed
        TSeed const & cs = s;
        typedef Iterator<TSeed const, Standard>::Type TIterator;
        TIterator it = begin(cs);
        SEQAN_ASSERT_EQ(1u, it->beginDim0);
        SEQAN_ASSERT_EQ(2u, it->beginDim1);
        SEQAN_ASSERT_EQ(3u, it->length);
        ++it;
        SEQAN_ASSERT_EQ(4u, it->beginDim0);
        SEQAN_ASSERT_EQ(5u, it->beginDim1);
        SEQAN_ASSERT_EQ(3u, it->length);
        ++it;
        SEQAN_ASSERT_TRUE(it == end(cs));
    }
}

#endif  // TEST_SEEDS_TEST_SEEDS_SEED_CHAINED_H_
