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
.Spec.Mersenne Twister RNG
..general:Class.RNG
..summary:Mersenne Twister 19937 Random Number Generator
..cat:Random
..include:seqan/random.h
*/
template <>
class RNG<MersenneTwister>
{
public:
    ext::MTRand _mtRand;


/**
.Memfunc.Mersenne Twister RNG#RNG
..class:Spec.Mersenne Twister RNG
..summary:Constructor Mersenne Twister RNG.
..signature:RNG<MersenneTwister>([seed])
..param.seed:Seed for the initialization of the Mersenne Twister, defaults to 0.
...type:nolink:double.
*/
    RNG() : _mtRand(0lu)
    { SEQAN_CHECKPOINT; }

    RNG(unsigned seed) : _mtRand(seed)
    { SEQAN_CHECKPOINT; }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <>
struct Value<RNG<MersenneTwister> >
{
    typedef unsigned Type;
};

template <>
struct Value<const RNG<MersenneTwister> > : Value<RNG<MersenneTwister> > {};

template <>
struct SupremumValue<RNG<MersenneTwister> >
{
    typedef RNG<MersenneTwister> _TMT;
    typedef Value<_TMT>::Type _TValue;
    static const _TValue VALUE;
};

const Value<RNG<MersenneTwister> >::Type SupremumValue<RNG<MersenneTwister> >::VALUE = SupremumValue<Value<RNG<MersenneTwister> >::Type>::VALUE;

template <>
struct SupremumValue<const RNG<MersenneTwister> > : SupremumValue<RNG<MersenneTwister> > {};

template <>
struct InfimumValue<RNG<MersenneTwister> >
{
    typedef RNG<MersenneTwister> _TMT;
    typedef Value<_TMT>::Type _TValue;
    static const _TValue VALUE;
};

const Value<RNG<MersenneTwister> >::Type InfimumValue<RNG<MersenneTwister> >::VALUE = InfimumValue<Value<RNG<MersenneTwister> >::Type>::VALUE;

template <>
struct InfimumValue<const RNG<MersenneTwister> > : InfimumValue<RNG<MersenneTwister> > {};

// ===========================================================================
// Functions
// ===========================================================================

inline unsigned
pickRandomNumber(RNG<MersenneTwister> & mt)
{
    SEQAN_CHECKPOINT;

    return mt._mtRand.randInt();
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_MT19937_H_
