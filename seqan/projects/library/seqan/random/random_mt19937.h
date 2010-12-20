 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
 ============================================================================
  An implementation of the Mersenne Twister 19937 random number generator.
 ==========================================================================*/

#ifndef SEQAN_RANDOM_RANDOM_MT19937_H_
#define SEQAN_RANDOM_RANDOM_MT19937_H_

#include <seqan/random/ext_MersenneTwister.h>

namespace seqan {

// ===========================================================================
// Forwards, Tags.
// ===========================================================================

// Tag for selecting a mersenne twister.
struct MersenneTwister {};

// ===========================================================================
// Classes
// ===========================================================================

/**
.Spec.Mersenne Twister Rng
..general:Class.Rng
..signature:Rng<MersenneTwister>
..summary:Mersenne Twister 19937 Random Number Generator
..cat:Random
..include:seqan/random.h
*/
template <>
class Rng<MersenneTwister>
{
public:
    ext::MTRand _mtRand;


/**
.Memfunc.Mersenne Twister Rng#Rng
..class:Spec.Mersenne Twister Rng
..summary:Constructor Mersenne Twister Rng.
..signature:Rng<MersenneTwister>([seed])
..param.seed:Seed for the initialization of the Mersenne Twister, defaults to 0.
...type:nolink:double.
*/
    Rng() : _mtRand(0lu)
    { SEQAN_CHECKPOINT; }

    Rng(unsigned seed) : _mtRand(seed)
    { SEQAN_CHECKPOINT; }
    
    inline
    unsigned
    operator()()
    {
        SEQAN_CHECKPOINT;
        return _mtRand.randInt();
    }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <>
struct Value<Rng<MersenneTwister> >
{
    typedef unsigned Type;
};

template <>
struct Value<const Rng<MersenneTwister> > : Value<Rng<MersenneTwister> > {};

template <>
struct MaxValue<Rng<MersenneTwister> >
{
    typedef Rng<MersenneTwister> _TMT;
    typedef Value<_TMT>::Type _TValue;
    static const _TValue VALUE;
};

const Value<Rng<MersenneTwister> >::Type MaxValue<Rng<MersenneTwister> >::VALUE = MaxValue<Value<Rng<MersenneTwister> >::Type>::VALUE;

template <>
struct MaxValue<const Rng<MersenneTwister> > : MaxValue<Rng<MersenneTwister> > {};

template <>
struct MinValue<Rng<MersenneTwister> >
{
    typedef Rng<MersenneTwister> _TMT;
    typedef Value<_TMT>::Type _TValue;
    static const _TValue VALUE;
};

const Value<Rng<MersenneTwister> >::Type MinValue<Rng<MersenneTwister> >::VALUE = MinValue<Value<Rng<MersenneTwister> >::Type>::VALUE;

template <>
struct MinValue<const Rng<MersenneTwister> > : MinValue<Rng<MersenneTwister> > {};

// ===========================================================================
// Functions
// ===========================================================================

inline unsigned
pickRandomNumber(Rng<MersenneTwister> & mt)
{
    SEQAN_CHECKPOINT;

    return mt._mtRand.randInt();
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_MT19937_H_
