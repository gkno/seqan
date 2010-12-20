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
  Code for geometrically distributed random number generation, p=0.5.
 ===========================================================================
*/

#ifndef SEQAN_RANDOM_RANDOM_GEOMETRIC_H_
#define SEQAN_RANDOM_RANDOM_GEOMETRIC_H_

namespace seqan {

// ===========================================================================
// Forwards, Tags.
// ===========================================================================

// Specialization Tag for geometric distribution.
struct GeometricFairCoin {};

// ===========================================================================
// Classes
// ===========================================================================

/**
.Spec.Geometric Pdf
..signature:Pdf<GeometricFairCoin>
..general:Class.Pdf
..summary:Geometric probability density function with $p=0.5$.

This can be implemented efficiently not using any floating point arithmetics.
Just bit operations are needed.
..cat:Random
..include:seqan/random.h
*/
template <>
class Pdf<GeometricFairCoin>
{
public:

/**
.Memfunc.Geometric Pdf#Pdf
..class:Spec.Geometric Pdf
..summary:Constructor for geometric Pdf.
..signature:Pdf<Geometric>()
*/
    Pdf() { SEQAN_CHECKPOINT; }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <>
struct Value<Pdf<GeometricFairCoin> >
{
    typedef unsigned Type;
};

template <>
struct Value<const Pdf<GeometricFairCoin> > : Value<Pdf<GeometricFairCoin> > {};

// ===========================================================================
// Functions
// ===========================================================================

/*
..summary:Pick a geometricly distributed random number.
*/
template <typename TRNG>
inline
typename Value<Pdf<GeometricFairCoin> >::Type
pickRandomNumber(TRNG & rng, Pdf<GeometricFairCoin> const & /*pdf*/)
{
    SEQAN_CHECKPOINT;

    const int RG_IB1 = 1;
    const int RG_IB2 = 2;
    const int RG_IB5 = 16;
    const int RG_IB18 = 131072;
    const int RG_MASK = RG_IB1 + RG_IB2 + RG_IB5;
    
    typename Value<TRNG>::Type seed = pickRandomNumber(rng);
    typename Value<Pdf<GeometricFairCoin> >::Type value = 0;

    while (true) {
        if ((seed & RG_IB18)) {
            seed = ((seed ^ RG_MASK) << 1) | RG_IB1;
            ++value;
        } else {
            seed <<= 1;
            break;
        }
    }

    return value;
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_GEOMETRIC_H_
