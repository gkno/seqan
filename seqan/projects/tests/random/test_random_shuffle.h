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
  Tests for the shuffle() function.
  ===========================================================================
*/

#ifndef TEST_RANDOM_TEST_RANDOM_SHUFFLE_H_
#define TEST_RANDOM_TEST_RANDOM_SHUFFLE_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/random.h>

SEQAN_DEFINE_TEST(test_random_shuffle)
{
    using namespace seqan;

    RNG<MersenneTwister> mt(0);
    CharString container = "Hello!";
    CharString const before = container;
    
    shuffle(container, mt);
    SEQAN_ASSERT_NEQ(before, container);
}

#endif  // TEST_RANDOM_TEST_RANDOM_SHUFFLE_H_
