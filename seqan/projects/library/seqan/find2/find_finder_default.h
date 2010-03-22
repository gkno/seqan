/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007-2010

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
  Implementation of the default finder.
 ==========================================================================*/

#ifndef SEQAN_FIND2_FINDER_DEFAULT_H_
#define SEQAN_FIND2_FINDER_DEFAULT_H_

namespace seqan {

template <typename THaystack>
struct Finder2<THaystack, void> : _FindState {
    typedef typename Position<THaystack>::Type TPosition;
    typedef typename Iterator<THaystack>::Type TIterator;

    TState _state;

    Holder<THaystack> _holder;

    // TODO(holtgrew): Maybe switch to an iterator based implementation?
    TPosition _beginPosition;
    TPosition _endPosition;

    Finder2(THaystack & hstck)
        : _state(STATE_INITIAL),
          _holder(hstck) {
        SEQAN_CHECKPOINT;
    }
};


template <typename THaystack>
typename Position<THaystack>::Type & beginPosition(const Finder2<THaystack, void> & finder) {
    SEQAN_CHECKPOINT;
    return finder._beginPosition;
}


template <typename THaystack>
typename Position<THaystack>::Type & endPosition(const Finder2<THaystack, void> & finder) {
    SEQAN_CHECKPOINT;
    return finder._endPosition;
}


template <typename THaystack>
THaystack & haystack2(Finder2<THaystack, void> & finder) {
    SEQAN_CHECKPOINT;
    return value(finder._holder);
}

}  // namespace seqan
          
#endif  // SEQAN_FIND2_FINDER_DEFAULT_H_
