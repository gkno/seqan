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
  Code for uniformly distributed random number generation.
 ===========================================================================
*/

#ifndef SEQAN_RANDOM_RANDOM_UNIFORM_H_
#define SEQAN_RANDOM_RANDOM_UNIFORM_H_

namespace seqan {

// ===========================================================================
// Forwards, Tags.
// ===========================================================================

// Specialization tag for uniform distribution.
template <typename T>
struct Uniform;

// ===========================================================================
// Classes
// ===========================================================================

/**
.Spec.Uniform PDF
..signature:PDF<Uniform<T> >
..general:Class.PDF
..summary:Uniform distribution probability density function over a closed interval [min, max].
..param.T:Type of the values the PDF is defined on.
..cat:Random
..include:seqan/random.h
*/
template <typename T>
class PDF<Uniform<T> >
{
public:
    T _min;
    T _max;

// TODO(holtgrew): Switch to [begin, end) instead of [min, max] style?
/**
.MemfuncUniform PDF #PDF
..class:Spec.Uniform PDF
..summary:Constructor for uniform PDF.
..signature:PDF<Uniform<T> >(min, max)
..param.min:Smallest value of interval.
...type:nolink:T
..param.max:Largest value of interval.
...type:nolink:T
*/
    PDF(T min, T max)
            : _min(min), _max(max)
    {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_LEQ(_min, _max);
    }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <typename T>
struct Value<PDF<Uniform<T> > >
{
    typedef T Type;
};

template <typename T>
struct Value<const PDF<Uniform<T> > > : Value<PDF<Uniform<T> > > {};

// ===========================================================================
// Functions
// ===========================================================================

// Pick an integral random number uniformly distributed.
template <typename TRNG, typename T>
inline
typename Value<PDF<Uniform<T> > >::Type
_pickRandomNumber(TRNG & rng, PDF<Uniform<T> > const & pdf, True const &)
{
    SEQAN_CHECKPOINT;
    typename Value<TRNG>::Type limit = (SupremumValue<TRNG>::VALUE / (pdf._max - pdf._min)) * (pdf._max - pdf._min);
    typename Value<TRNG>::Type x;
    do {
        x = pickRandomNumber(rng);
    } while (x > limit);
    T y = x % (pdf._max - pdf._min + 1);
    return y + pdf._min;
}

// Pick a continuous random number uniformly distributed.
template <typename TRNG, typename T>
inline
typename Value<PDF<Uniform<T> > >::Type
_pickRandomNumber(TRNG & rng, PDF<Uniform<T> > const & pdf, False const &)
{
    SEQAN_CHECKPOINT;
    T x = static_cast<T>(pickRandomNumber(rng) - InfimumValue<TRNG>::VALUE);
    x /= static_cast<T>(SupremumValue<TRNG>::VALUE) - static_cast<T>(InfimumValue<TRNG>::VALUE);
    return pdf._min + x * (pdf._max - pdf._min);
}

template <typename TRNG, typename T>
inline
typename Value<PDF<Uniform<T> > >::Type
pickRandomNumber(TRNG & rng, PDF<Uniform<T> > const & pdf)
{
    SEQAN_CHECKPOINT;
    if (pdf._min == pdf._max)
        return pdf._min;
    return _pickRandomNumber(rng, pdf, typename IsIntegral<T>::Type());
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_UNIFORM_H_
