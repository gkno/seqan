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

template <typename _THaystack>
struct Finder<_THaystack, Default> : _FindState {
    typedef _THaystack THaystack;
    typedef typename Position<THaystack>::Type TPosition;
    typedef typename Iterator<THaystack>::Type TIterator;

    TState _state;

    Holder<THaystack> _holder;

    // TODO(holtgrew): Maybe switch to an iterator based implementation?
    TPosition _beginPosition;
    TPosition _endPosition;

    Finder(THaystack & haystack)
        : _state(STATE_INITIAL),
          _holder(haystack) {
        SEQAN_CHECKPOINT;
    }
};


template <typename THaystack>
THaystack const & haystack(Finder<THaystack, Default> const & finder) {
    SEQAN_CHECKPOINT;
    return value(finder._holder);
}


template <typename THaystack, typename TTag>
typename Iterator<THaystack const>::Type begin(Finder<THaystack, Default> const & finder,
                                               Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    SEQAN_ASSERT_TRUE(finder._state == TFinder::STATE_BEGIN_FOUND);
    return begin(haystack(finder), spec) + finder._beginPosition;
}


// TODO(holtgrew): Workaround for default begin implementation for non-const type.
template <typename THaystack, typename TTag>
typename Iterator<THaystack const>::Type begin(Finder<THaystack, Default> & finder,
                                               Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    return begin(const_cast<TFinder const &>(finder), spec);
}


template <typename THaystack, typename TTag>
typename Iterator<THaystack const>::Type end(Finder<THaystack, Default> const & finder,
                                             Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    SEQAN_ASSERT_TRUE(finder._state == TFinder::STATE_BEGIN_FOUND || finder._state == TFinder::STATE_FOUND);
    return begin(haystack(finder), spec) + finder._endPosition;
}


// TODO(holtgrew): Workaround for default end implementation for non-const type.
template <typename THaystack, typename TTag>
typename Iterator<THaystack const>::Type end(Finder<THaystack, Default> & finder,
                                             Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    return end(const_cast<TFinder const &>(finder), spec);
}


template <typename THaystack>
typename Position<THaystack const>::Type beginPosition(Finder<THaystack, Default> const & finder) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    SEQAN_ASSERT_TRUE(finder._state == TFinder::STATE_BEGIN_FOUND);
    return finder._beginPosition;
}


// TODO(holtgrew): Workaround for default beginPosition implementation for non-const type.
template <typename THaystack>
typename Position<THaystack>::Type beginPosition(Finder<THaystack, Default> & finder) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    return beginPosition(const_cast<TFinder const &>(finder));
}


template <typename THaystack>
typename Position<THaystack>::Type endPosition(Finder<THaystack, Default> const & finder) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    SEQAN_ASSERT_TRUE(finder._state == TFinder::STATE_BEGIN_FOUND || finder._state == TFinder::STATE_FOUND);
    return finder._endPosition;
}


// TODO(holtgrew): Workaround for default endPosition implementation for non-const type.
template <typename THaystack>
typename Position<THaystack>::Type endPosition(Finder<THaystack, Default> & finder) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    return endPosition(const_cast<TFinder const &>(finder));
}


template <typename THaystack>
Segment<THaystack const, InfixSegment> infix(Finder<THaystack, Default> & finder) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    SEQAN_ASSERT_TRUE(finder._state == TFinder::STATE_BEGIN_FOUND);
    return infix(haystack(finder), beginPosition(finder), endPosition(finder));
}

}  // namespace seqan
          
#endif  // SEQAN_FIND2_FINDER_DEFAULT_H_
