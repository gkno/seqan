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
  Tests for the probability distribution code in seqan/random.
  ===========================================================================
*/

#ifndef TEST_RANDOM_TEST_RANDOM_DISTS_H_
#define TEST_RANDOM_TEST_RANDOM_DISTS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

SEQAN_DEFINE_TEST(test_random_normal_constructors)
{
    using namespace seqan;
    PDF<Normal> pdf(1, 1);
}

SEQAN_DEFINE_TEST(test_random_normal_pick)
{
    using namespace seqan;

    
    // mean = 1.0, stddev = 1.0
    {
        RNG<MersenneTwister> mt(42);
        PDF<Normal> pdf(1, 1);
        
        double sum = 0;
        for (unsigned i = 0; i < 10000; ++i)
            sum += pickRandomNumber(mt, pdf);
        SEQAN_ASSERT_LT(fabs(sum / 10000 - 1), 0.02);
    }
    // mean = -3.4, stddev = 0.3
    {
        RNG<MersenneTwister> mt(42);
        PDF<Normal> pdf(-3.4, 0.3);
        
        double sum = 0;
        for (unsigned i = 0; i < 10000; ++i)
            sum += pickRandomNumber(mt, pdf);
        SEQAN_ASSERT_LT(fabs(sum / 10000 + 3.4), 0.02);
    }
    // mean = 4, stddev = 0.1
    {
        RNG<MersenneTwister> mt(42);
        PDF<Normal> pdf(4, 0.1);
        
        double sum = 0;
        for (unsigned i = 0; i < 10000; ++i)
            sum += pickRandomNumber(mt, pdf);
        SEQAN_ASSERT_LT(fabs(sum / 10000 - 4), 0.02);
    }
}

SEQAN_DEFINE_TEST(test_random_lognormal_constructors)
{
    using namespace seqan;

    {
        PDF<LogNormal> pdf(1, 1);
        SEQAN_ASSERT_EQ(pdf._normalDist._mu, 1);
        SEQAN_ASSERT_EQ(pdf._normalDist._sigma, 1);
    }
    {
        PDF<LogNormal> pdf(1, 1, MuSigma());
        SEQAN_ASSERT_EQ(pdf._normalDist._mu, 1);
        SEQAN_ASSERT_EQ(pdf._normalDist._sigma, 1);
    }
    {
        PDF<LogNormal> pdf(1, 1, MeanStdDev());
        SEQAN_ASSERT_NEQ(pdf._normalDist._mu, 1);
        SEQAN_ASSERT_NEQ(pdf._normalDist._sigma, 1);
    }
}

SEQAN_DEFINE_TEST(test_random_lognormal_pick)
{
    using namespace seqan;

    // mean = 1.0, stddev = 0.2
    {
        RNG<MersenneTwister> mt(42);
        PDF<LogNormal> pdf(1, 0.2, MeanStdDev());
        
        double sum = 0;
        for (unsigned i = 0; i < 10000; ++i)
            sum += pickRandomNumber(mt, pdf);
        SEQAN_ASSERT_LT(fabs(sum / 10000 - 1), 0.02);
    }
    // mean = 0.4, stddev = 0.1
    {
        RNG<MersenneTwister> mt(42);
        PDF<LogNormal> pdf(0.4, 0.1, MeanStdDev());
        
        double sum = 0;
        for (unsigned i = 0; i < 10000; ++i)
            sum += pickRandomNumber(mt, pdf);
        SEQAN_ASSERT_LT(fabs(sum / 10000 - 0.4), 0.02);
    }
    // mean = 4, stddev = 0.1
    {
        RNG<MersenneTwister> mt(42);
        PDF<LogNormal> pdf(4, 0.1, MeanStdDev());

        double sum = 0;
        for (unsigned i = 0; i < 10000; ++i)
            sum += pickRandomNumber(mt, pdf);
        SEQAN_ASSERT_LT(fabs(sum / 10000 - 4), 0.02);
    }
}

#endif  // TEST_RANDOM_TEST_RANDOM_DISTS_H_
