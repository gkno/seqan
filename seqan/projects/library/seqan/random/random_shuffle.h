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
  Shuffling.
  ===========================================================================
*/

#ifndef SEQAN_RANDOM_RANDOM_SHUFFLE_H_
#define SEQAN_RANDOM_RANDOM_SHUFFLE_H_

namespace seqan {

// ===========================================================================
// Forwards, Tags.
// ===========================================================================

// ===========================================================================
// Classes
// ===========================================================================

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

/**
.Function.shuffle
..summary:Shuffle the given container.
..cat:Random
..include:seqan/random.h
..signature:shuffle(container, rng)
..param.container:Container to shuffle elements of.
..param.rng:Random number generator to use.
...type:Class.RNG
*/
template <typename TContainer, typename TRNG>
void shuffle(TContainer & container, TRNG & rng)
{
    SEQAN_CHECKPOINT;
    typedef typename Position<TContainer>::Type TPosition;
    typedef typename Value<TContainer>::Type TValue;

    TValue tmp;
    for (TPosition i = 0, iend = length(container); i < iend; ++i) {
        PDF<Uniform<TPosition> > uniformDist(i, iend - 1);
        TPosition j = pickRandomNumber(rng, uniformDist);
        // swap
        move(tmp, container[i]);
        move(container[i], container[j]);
        move(container[j], tmp);
    }
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_SHUFFLE_H_
