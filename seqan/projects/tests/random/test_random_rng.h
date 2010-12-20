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
        Rng<> mt;
    }
    {
        Rng<MersenneTwister> mt;
    }
    {
        Rng<MersenneTwister> mt(10);
    }
}

// Pick random numbers from the MT and make sure the same number is
// not returned twice in the first two picks.
SEQAN_DEFINE_TEST(test_random_mt19937_pick)
{
    using namespace seqan;

    Rng<MersenneTwister> mt(10);
    // Test pickRandomNumber().
    SEQAN_ASSERT_NEQ(pickRandomNumber(mt), pickRandomNumber(mt));
    // Test operator().
    SEQAN_ASSERT_NEQ(mt(), mt());
}

// Construct RngFunctor specialization in all possible ways.
SEQAN_DEFINE_TEST(test_random_rng_functor_constructors)
{
    using namespace seqan;
    
	typedef Rng<MersenneTwister> TMersenneTwister;
	typedef Pdf<Uniform<int> > TUniformPdf;
    
    {
        TMersenneTwister mt;
        TUniformPdf uniformPdf(0, 10);
        
        Rng<RngFunctor<TMersenneTwister, TUniformPdf> > rng(mt, uniformPdf);
    }
}

// Pick random number from RngFunctor and make sure it is the same as when
// directly using a MT and a Pdf.
SEQAN_DEFINE_TEST(test_random_rng_functor_pick)
{
    using namespace seqan;
    
    const int SEED = 10;

    typedef Rng<MersenneTwister> TMersenneTwister;
	typedef Pdf<Uniform<int> > TUniformPdf;
    typedef Rng<RngFunctor<TMersenneTwister, TUniformPdf> > TRngFunctor;

    // Compute by using the raw MT and uniform Pdf.
    String<int> rawInts;
    {
	    TMersenneTwister mt(SEED);
    	TUniformPdf uniform(10, 100);
        
        for (int i = 0; i < 100; ++i)
            appendValue(rawInts, pickRandomNumber(mt, uniform));
    }
    
    // Use RngFunctor with pickRandomNumber() and check equality.
    {
	    TMersenneTwister mt(SEED);
    	TUniformPdf uniform(10, 100);
	    TRngFunctor rngFunctor(mt, uniform);
        
        for (int i = 0; i < 100; ++i) {
            SEQAN_ASSERT_EQ_MSG(unsigned(rawInts[i]), pickRandomNumber(rngFunctor), "i = %d", i);
        }
    }
    
    // Use RngFunctor with operator() and check equality.
    {
	    TMersenneTwister mt(SEED);
    	TUniformPdf uniform(10, 100);
	    TRngFunctor rngFunctor(mt, uniform);
        
        for (int i = 0; i < 100; ++i) {
            SEQAN_ASSERT_EQ_MSG(rawInts[i], rngFunctor(), "i = %d", i);
        }
    }
}

#endif  // TEST_RANDOM_TEST_RANDOM_RNG_H_
