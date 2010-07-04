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
  Tests for the random number generation code in seqan/random.
  ===========================================================================
*/

#ifndef TEST_RANDOM_TEST_RANDOM_RNG_H_
#define TEST_RANDOM_TEST_RANDOM_RNG_H_

// Construct MersenneTwister in all possible ways.
SEQAN_DEFINE_TEST(test_random_mt19937_constructors)
{
    using namespace seqan;

    {
        RNG<MersenneTwister> mt;
    }
    {
        RNG<MersenneTwister> mt(10);
    }
}

// Pick random numbers from the MT and make sure the same number is
// not returned twice in the first two picks.
SEQAN_DEFINE_TEST(test_random_mt19937_pick)
{
    using namespace seqan;

    RNG<MersenneTwister> mt(10);
    SEQAN_ASSERT_NEQ(pickRandomNumber(mt), pickRandomNumber(mt));
}

#endif  // TEST_RANDOM_TEST_RANDOM_RNG_H_
