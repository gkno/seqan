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
   Specific tests for the combination of seeds.  This header should
   contain explicit tests for each combination type (i.e. merging,
   different chaining) for both Simple and Chained Seeds.  The code is
   already rudimentarily tested in test_seeds_seed_set_base.h.
   However, we only test the later begin/end dimensions there and not
   the diagonals for the seed sets.
  ===========================================================================
*/

#ifndef TEST_SEEDS_TEST_SEEDS_COMBINATION_H_
#define TEST_SEEDS_TEST_SEEDS_COMBINATION_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.

SEQAN_DEFINE_TEST(test_seeds_combination_seeds_combineable_merge_chained)
{
    using namespace seqan;

    typedef Seed<ChainedSeed> TSeed;
    typedef Value<Seed<ChainedSeed> >::Type TDiagonal;

    // Case: Is Mergeable, partial overlap.
    {
        TSeed left(0, 0, 3);
        TSeed right(2, 2, 3);
        SEQAN_ASSERT_TRUE(_seedsCombineable(left, right, 1, Nothing(), Merge()));
    }
    // Case: Is Mergeable, complete overlap.
    {
        TSeed left(0, 0, 3);
        TSeed right(0, 0, 3);
        SEQAN_ASSERT_TRUE(_seedsCombineable(left, right, 1, Nothing(), Merge()));
    }
    // Case: Is Mergeable, prefix overlap.
    {
        TSeed left(0, 0, 3);
        TSeed right(0, 0, 4);
        SEQAN_ASSERT_TRUE(_seedsCombineable(left, right, 1, Nothing(), Merge()));
    }
    // Case: Is Mergeable, suffix overlap.
    {
        TSeed left(0, 0, 4);
        TSeed right(1, 1, 3);
        SEQAN_ASSERT_TRUE(_seedsCombineable(left, right, 1, Nothing(), Merge()));
    }
    // Case: Not Mergeable, too large diagonal distance.
    {
        TSeed left(0, 0, 4);
        TSeed right(2, 0, 3);
        SEQAN_ASSERT_NOT(_seedsCombineable(left, right, 1, Nothing(), Merge()));
    }
    // Case: Not Mergeable, no overlap.
    {
        TSeed left(0, 0, 4);
        TSeed right(5, 5, 3);
        SEQAN_ASSERT_NOT(_seedsCombineable(left, right, 1, Nothing(), Merge()));
    }
    // Case: Not Mergeable, second starts right of first one.
    {
        TSeed left(2, 2, 3);
        TSeed right(0, 0, 3);
        SEQAN_ASSERT_NOT(_seedsCombineable(left, right, 1, Nothing(), Merge()));
    }
}


SEQAN_DEFINE_TEST(test_seeds_combination_seeds_combineable_simple_chaining_chained)
{
    using namespace seqan;

    typedef Seed<ChainedSeed> TSeed;
    typedef Value<Seed<ChainedSeed> >::Type TDiagonal;

    // Case: Chaineable.
    {
        TSeed left(0, 0, 3);
        TSeed right(4, 5, 3);
        SEQAN_ASSERT_TRUE(_seedsCombineable(left, right, 3, Nothing(), SimpleChain()));
    }
    // Case: Not Chaineable, overlap.
    {
        TSeed left(0, 0, 3);
        TSeed right(1, 0, 3);
        SEQAN_ASSERT_NOT(_seedsCombineable(left, right, 3, Nothing(), SimpleChain()));
    }
    // Case: Not Chainable, too far away.
    {
        TSeed left(0, 0, 3);
        TSeed right(6, 7, 3);
        SEQAN_ASSERT_NOT(_seedsCombineable(left, right, 3, Nothing(), SimpleChain()));
    }
}


SEQAN_DEFINE_TEST(test_seeds_combination_seeds_combineable_simple_chaos_chaining_chained)
{
    using namespace seqan;
    typedef Seed<ChainedSeed> TSeed;
    typedef Value<Seed<ChainedSeed> >::Type TDiagonal;

    // Case: Chaineable.
    {
        TSeed left(0, 0, 3);
        TSeed right(4, 5, 3);
        SEQAN_ASSERT_TRUE(_seedsCombineable(left, right, 3, Nothing(), SimpleChain()));
    }
    // Case: Not Chaineable, overlap.
    {
        TSeed left(0, 0, 3);
        TSeed right(1, 0, 3);
        SEQAN_ASSERT_NOT(_seedsCombineable(left, right, 3, Nothing(), SimpleChain()));
    }
    // Case: Not Chainable, too far away.
    {
        TSeed left(0, 0, 3);
        TSeed right(6, 7, 3);
        SEQAN_ASSERT_NOT(_seedsCombineable(left, right, 3, Nothing(), SimpleChain()));
    }
}

SEQAN_DEFINE_TEST(test_seeds_combination_combine_seeds_merge_chained)
{
    // TODO(holtgrew): Add scoring to the tests.
    using namespace seqan;

    typedef Seed<ChainedSeed> TSeed;
    typedef Value<Seed<ChainedSeed> >::Type TDiagonal;

    // Case: Each has one diagonal, first one is shortened.
    {
        TSeed left;
        appendDiagonal(left, TDiagonal(0, 0, 3));
        TSeed right;
        appendDiagonal(right, TDiagonal(2, 1, 3));
        _combineSeeds(left, right, Score<int, Simple>(0, -1, -1), Nothing(), Nothing(), Merge());
        TSeed expected;
        appendDiagonal(expected, TDiagonal(0, 0, 1));
        appendDiagonal(expected, TDiagonal(2, 1, 3));
        SEQAN_ASSERT_EQ(expected, left);
        // TODO(holtgrew): Check lower and upper diagonal.
    }
    // Case: First has two diagonals, one must be deleted, one shortened.
    {
        TSeed left;
        appendDiagonal(left, TDiagonal(0, 0, 3));
        appendDiagonal(left, TDiagonal(4, 5, 3));
        TSeed right;
        appendDiagonal(right, TDiagonal(6, 6, 3));
        _combineSeeds(left, right, Score<int, Simple>(0, -1, -1), Nothing(), Nothing(), Merge());
        TSeed expected;
        appendDiagonal(expected, TDiagonal(0, 0, 3));
        appendDiagonal(expected, TDiagonal(4, 5, 1));
        appendDiagonal(expected, TDiagonal(6, 6, 3));
        SEQAN_ASSERT_EQ(expected, left);
        // TODO(holtgrew): Check lower and upper diagonal.
    }
    // Case: First has two diagonals, one must be deleted, first one is kept.
    {
        TSeed left;
        appendDiagonal(left, TDiagonal(0, 0, 3));
        appendDiagonal(left, TDiagonal(4, 5, 3));
        TSeed right;
        appendDiagonal(right, TDiagonal(4, 5, 3));
        _combineSeeds(left, right, Score<int, Simple>(0, -1, -1), Nothing(), Nothing(), Merge());
        TSeed expected;
        appendDiagonal(expected, TDiagonal(0, 0, 3));
        appendDiagonal(expected, TDiagonal(4, 5, 3));
        SEQAN_ASSERT_EQ(expected, left);
        // TODO(holtgrew): Check lower and upper diagonal.
    }
    // Case: First has two diagonals, all are kept.
    {
        TSeed left;
        appendDiagonal(left, TDiagonal(0, 0, 3));
        appendDiagonal(left, TDiagonal(4, 5, 3));
        TSeed right;
        appendDiagonal(right, TDiagonal(7, 7, 3));
        _combineSeeds(left, right, Score<int, Simple>(0, -1, -1), Nothing(), Nothing(), Merge());
        TSeed expected;
        appendDiagonal(expected, TDiagonal(0, 0, 3));
        appendDiagonal(expected, TDiagonal(4, 5, 2));
        appendDiagonal(expected, TDiagonal(7, 7, 3));
        SEQAN_ASSERT_EQ(expected, left);
        // TODO(holtgrew): Check lower and upper diagonal.
    }
}


SEQAN_DEFINE_TEST(test_seeds_combination_combine_seeds_simple_chaining_chained)
{
    // TODO(holtgrew): Also test score updates...
    using namespace seqan;

    typedef Seed<ChainedSeed> TSeed;
    typedef Value<Seed<ChainedSeed> >::Type TDiagonal;
    
    // Simple chaining is easy, the diagonals are simply copied over.
    TSeed left;
    appendDiagonal(left, TDiagonal(1, 2, 3));
    appendDiagonal(left, TDiagonal(5, 6, 3));
    TSeed right;
    appendDiagonal(right, TDiagonal(10, 10, 3));
    appendDiagonal(right, TDiagonal(15, 14, 3));
    _combineSeeds(left, right, Score<int, Simple>(0, -1, -1), Nothing(), Nothing(), SimpleChain());
    TSeed expected;
    appendDiagonal(expected, TDiagonal(1, 2, 3));
    appendDiagonal(expected, TDiagonal(5, 6, 3));
    appendDiagonal(expected, TDiagonal(10, 10, 3));
    appendDiagonal(expected, TDiagonal(15, 14, 3));
    SEQAN_ASSERT_EQ(expected, left);
}


SEQAN_DEFINE_TEST(test_seeds_combination_combine_seeds_simple_chaos_chaining_chained)
{
    // TODO(holtgrew): Also test score updates...
    using namespace seqan;

    typedef Seed<ChainedSeed> TSeed;
    typedef Value<Seed<ChainedSeed> >::Type TDiagonal;

    // Simple case: Both seeds have one diagonal.    
    {
        DnaString sequence0 = "ACAAAC";
        DnaString sequence1 = "ACACCAAC";
        TSeed left(0, 0, 2);
        TSeed right(6, 4, 2);
        _combineSeeds(left, right, Score<int, Simple>(0, -1, -1), sequence0, sequence1, Chaos());
        TSeed expected;
        appendDiagonal(expected, TDiagonal(0, 0, 3));
        appendDiagonal(expected, TDiagonal(5, 3, 3));
        SEQAN_ASSERT_EQ(expected, left);
    }
    // Each seed has two diagonals.  The same example as above but CC
    // prepended and appended to both sequences.
    {
        DnaString sequence0 = "CCACAAACCC";
        DnaString sequence1 = "CCACACCAACCC";
        TSeed left;
        appendDiagonal(left, TDiagonal(0, 0, 2));
        appendDiagonal(left, TDiagonal(2, 2, 2));
        TSeed right;
        appendDiagonal(right, TDiagonal(8, 6, 2));
        appendDiagonal(right, TDiagonal(10, 8, 2));
        _combineSeeds(left, right, Score<int, Simple>(0, -1, -1), sequence0, sequence1, Chaos());
        TSeed expected;
        appendDiagonal(expected, TDiagonal(0, 0, 2));
        appendDiagonal(expected, TDiagonal(2, 2, 3));
        appendDiagonal(expected, TDiagonal(7, 5, 3));
        appendDiagonal(expected, TDiagonal(10, 8, 2));
        SEQAN_ASSERT_EQ(expected, left);
    }
}

#endif  // TEST_SEEDS_TEST_SEEDS_COMBINATION_H_
