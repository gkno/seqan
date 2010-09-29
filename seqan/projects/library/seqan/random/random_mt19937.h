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

namespace seqan {

// ===========================================================================
// Forwards, Tags.
// ===========================================================================

// Tag for selecting a mersenne twister.
struct MersenneTwister {};

// Forward declaration for mersenne twister, needed for _initialize.
template <>
class RNG<MersenneTwister>;

// Forward declaration for _initialize, needed for RNG<MersenneTwister>'s
// constructor.
inline void _initialize(RNG<MersenneTwister> & mt, unsigned seed);

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
    unsigned long _buffer[SEQAN_MERSENNE_MT_LEN];
	  int _index;


/**
.Memfunc.Mersenne Twister RNG#RNG
..class:Spec.Mersenne Twister RNG
..summary:Constructor Mersenne Twister RNG.
..signature:RNG<MersenneTwister>([seed])
..param.seed:Seed for the initialization of the Mersenne Twister, defaults to 0.
...type:nolink:double.
*/
    RNG() : _index(0)
    {
        SEQAN_CHECKPOINT;
        _initialize(*this, 0);
    }

    RNG(unsigned seed) : _index(0)
    {
        SEQAN_CHECKPOINT;
        _initialize(*this, seed);
    }
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

const Value<RNG<MersenneTwister> >::Type SupremumValue<RNG<MersenneTwister> >::VALUE = SupremumValue<Value<RNG<MersenneTwister> >::Type>::VALUE / 2;

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

	unsigned long * b = mt._buffer;
	int idx = mt._index;
	unsigned long s;
	int i;
	
	if (idx == SEQAN_MERSENNE_MT_LEN*sizeof(unsigned long))
	{
		idx = 0;
		i = 0;
		for (; i < SEQAN_MERSENNE_MT_IB; i++) {
			s = SEQAN_MERSENNE_TWIST(b, i, i+1);
			b[i] = b[i + SEQAN_MERSENNE_MT_IA] ^ (s >> 1) ^ SEQAN_MERSENNE_MAGIC(s);
		}
		for (; i < SEQAN_MERSENNE_MT_LEN-1; i++) {
			s = SEQAN_MERSENNE_TWIST(b, i, i+1);
			b[i] = b[i - SEQAN_MERSENNE_MT_IB] ^ (s >> 1) ^ SEQAN_MERSENNE_MAGIC(s);
		}
        
		s = SEQAN_MERSENNE_TWIST(b, SEQAN_MERSENNE_MT_LEN-1, 0);
		b[SEQAN_MERSENNE_MT_LEN-1] = b[SEQAN_MERSENNE_MT_IA-1] ^ (s >> 1) ^ SEQAN_MERSENNE_MAGIC(s);
	}
	mt._index = idx + sizeof(unsigned long);
	return *(unsigned long *)((unsigned char *)b + idx);
}


/**
.Internal._initalize
..summary:Initialize Mersenne Twister's state with the given seed.
..cat:Random
..include:seqan/random.h
..signature:_initialize(mt, seed)
..param.mt:Mersenne Twister to initialize.
...type:Spec.Mersenne Twister RNG
..param.seed:Seed for the initialization.
...type:nolink:unsigned
*/
inline void
_initialize(RNG<MersenneTwister> & mt, unsigned seed)
{
    SEQAN_CHECKPOINT;

    // Initialize C stdlib standard number generator with the given seed.
    ::std::srand(seed);

	for (int i = 0; i < SEQAN_MERSENNE_MT_LEN; i++)
		mt._buffer[i] = ::std::rand();
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_MT19937_H_
