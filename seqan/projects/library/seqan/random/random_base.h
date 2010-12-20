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
  Basic definitions for the module random.
  ===========================================================================
*/

#ifndef SEQAN_RANDOM_RANDOM_BASE_H_
#define SEQAN_RANDOM_RANDOM_BASE_H_

namespace seqan {

// ===========================================================================
// Forwards, Tags.
// ===========================================================================

// Forward to MersenneTwister, really defined in random_mt19937.h.
struct MersenneTwister;

/**
.Class.RNG:
..summary:Random Number Generator
..signature:RNG<>
..signature:RNG<TSpec>
..cat:Random
..param.TSpec:Random Number Generator specialization.
...default:@Spec.Mersenne Twister RNG@
..include:seqan/random.h
*/
template <typename TSpec = MersenneTwister>
class RNG;

/**
.Class.PDF:
..summary:ProbabilityDensityFunction
..signature:PDF<TSpec>
..cat:Random
..param.TSpec:Specialization.
..include:seqan/random.h
*/
template <typename TSpec>
class PDF;

// ===========================================================================
// Classes
// ===========================================================================

/**
.Memfunc.RNG#operator()
..class:Class.RNG
..summary:Function call operator.
..signature:operator()
*/
        
// ===========================================================================
// Metafunctions
// ===========================================================================

///.Metafunction.Value.param.T.type:Class.PDF
// specification only

///.Metafunction.Value.param.T.type:Class.RNG
///.Metafunction.InfimumValue.param.T.type:Class.RNG
///.Metafunction.SupremumValue.param.T.type:Class.RNG
// specification only

// ===========================================================================
// Functions
// ===========================================================================

/**
.Function.pickRandomNumber
..summary:Pick a random number using a random number generator object, possibly following the given distribution.
..cat:Random
..include:seqan/random.h
..signature:pickRandomNumber(rng[, pdf])
..param.rng:Random number generator to use.
...type:Class.RNG
..param.pdf:Probability density function to use, if any.
...type:Class.PDF
..returns:Random number as specified in pdf, if any, or rng.
*/
// specification only

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_BASE_H_
