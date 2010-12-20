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
  Functor for random number generation.
  ===========================================================================
*/

#ifndef SEQAN_RANDOM_RANDOM_RNG_FUNCTOR_H_
#define SEQAN_RANDOM_RANDOM_RNG_FUNCTOR_H_

namespace seqan {

// ===========================================================================
// Forwards, Tags.
// ===========================================================================

// Tag for selecting the RNG functor specialization.
template <typename TRng, typename TPdf>
struct RngFunctor {};
    
// ===========================================================================
// Classes
// ===========================================================================

/**
.Spec.RNG Functor
..general:Class.RNG
..signature:RNG<RngFunctor<TRng, TPdf> >
..summary:Functor wrapper for random number generation.
..cat:Random
..include:seqan/random.h
*/
template <typename TRng, typename TPdf>
class RNG<RngFunctor<TRng, TPdf> >
{
public:
    TRng & _rng;
    TPdf & _pdf;
    
/**
.Memfunc.RNG Functor#RNG
..class:Spec.RNG Functor
..summary:Constructor Functor RNG.
..signature:RNG<RngFunctor<TRng, TPdfd> >(rng, pdf)
..param.rng:@Class.RNG@ object to use.
..param.pdf:@Class.PDF@ object to use.
*/
    RNG(TRng & rng, TPdf & pdf)
	    : _rng(rng), _pdf(pdf)
    {}
    
    inline
    typename Value<TPdf>::Type
    operator()()
    {
        SEQAN_CHECKPOINT;
        return pickRandomNumber(_rng, _pdf);
    }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <typename TRng, typename TPdf>
struct Value<RngFunctor<TRng, TPdf> > : Value<TPdf> {};

template <typename TRng, typename TPdf>
struct Value<RngFunctor<TRng, TPdf> const> : Value<TPdf> {};

template <typename TRng, typename TPdf>
struct SupremumValue<RngFunctor<TRng, TPdf> > : SupremumValue<TPdf> {};

template <typename TRng, typename TPdf>
struct SupremumValue<RngFunctor<TRng, TPdf> const> : SupremumValue<TPdf> {};

template <typename TRng, typename TPdf>
struct InfimumValue<RngFunctor<TRng, TPdf> > : InfimumValue<TPdf> {};

template <typename TRng, typename TPdf>
struct InfimumValue<RngFunctor<TRng, TPdf> const> : InfimumValue<TPdf> {};
    
// ===========================================================================
// Functions
// ===========================================================================

template <typename TRng, typename TPdf>
inline unsigned
pickRandomNumber(RNG<RngFunctor<TRng, TPdf> > & rng)
{
    SEQAN_CHECKPOINT;
    
    return pickRandomNumber(rng._rng, rng._pdf);
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_RNG_FUNCTOR_H_
