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
.Class.Rng:
..summary:Random Number Generator
..signature:Rng<>
..signature:Rng<TSpec>
..cat:Random
..param.TSpec:Random Number Generator specialization.
...default:@Spec.Mersenne Twister Rng@
..include:seqan/random.h
*/
template <typename TSpec = MersenneTwister>
class Rng;

/**
.Class.Pdf:
..summary:ProbabilityDensityFunction
..signature:Pdf<TSpec>
..cat:Random
..param.TSpec:Specialization.
..include:seqan/random.h
*/
template <typename TSpec>
class Pdf;

// ===========================================================================
// Classes
// ===========================================================================

/**
.Memfunc.Rng#operator()
..class:Class.Rng
..summary:Function call operator.
..signature:operator()
*/
        
// ===========================================================================
// Metafunctions
// ===========================================================================

///.Metafunction.Value.param.T.type:Class.Pdf
// specification only

///.Metafunction.Value.param.T.type:Class.Rng
///.Metafunction.MinValue.param.T.type:Class.Rng
///.Metafunction.MaxValue.param.T.type:Class.Rng
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
...type:Class.Rng
..param.pdf:Probability density function to use, if any.
...type:Class.Pdf
..returns:Random number as specified in pdf, if any, or rng.
*/
// specification only

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_BASE_H_
