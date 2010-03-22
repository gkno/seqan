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
  Exact pattern matching using a naive implementation.
 ==========================================================================*/

#ifndef SEQAN_FIND2_FIND_SIMPLE_H_
#define SEQAN_FIND2_FIND_SIMPLE_H_

namespace seqan {

template <typename TNeedle>
struct Pattern2<TNeedle, Simple> : _FindState {
    // The pattern's state.
    TState _state;

    // The needle we store.
    Holder<TNeedle> _host;

    Pattern2() { SEQAN_CHECKPOINT; }

    Pattern2(TNeedle & ndl)
        : _state(STATE_INITIAL),
          _host(ndl)
    {}
};


template <typename TNeedle>
TNeedle & needle2(Pattern2<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}


template <typename TNeedle>
typename Position<TNeedle>::Type length(Pattern2<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    return length(needle2(pattern));
}


template <typename TNeedle>
Segment<TNeedle, InfixSegment> infix(Pattern2<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    return infix(pattern, 0u, length(needle2(pattern)));
}


template <typename TNeedle>
typename Iterator<TNeedle>::Type begin(Pattern2<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    return begin(needle2(pattern));
}


template <typename TNeedle>
typename Iterator<TNeedle>::Type end(Pattern2<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    return end(needle2(pattern));
}


template <typename TNeedle>
typename Position<TNeedle>::Type beginPosition(Pattern2<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    return 0u;
}


template <typename TNeedle>
typename Position<TNeedle>::Type endPosition(Pattern2<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    return length(needle2(pattern));
}


template <typename THaystack, typename TNeedle>
bool find(Finder2<THaystack, void> & finder,  // TODO(holtgrew): "Default" better than void?
          Pattern2<TNeedle, Simple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern2<THaystack, Simple> TFinder;
    typedef Pattern2<TNeedle, Simple> TPattern;
    typedef typename Position<TNeedle>::Type TPosition;

    // State of finder and pattern should be in sync.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);

    // Do not continue if the state is "not found".
    if (finder._state == TPattern::STATE_NOTFOUND)
        return false;
    // Initialize finder if state is "initial".
    if (finder._state == TPattern::STATE_INITIAL) {
        finder._beginPosition = 0u;
        finder._endPosition = length(needle(pattern));
    }

    // Search the needle in the haystack naively.
    while (finder._beginPosition < length(haystack2(finder)) - length(needle2(pattern))) {
        for (TPosition i = 0u; i < length(needle2(pattern)); ++i) {
            if (needle2(pattern)[i] != haystack2(finder)[finder._beginPosition + i]) {
                finder._beginPosition += 1;
                i = 0u;
                continue;
            }
        }
        finder._endPosition = finder._beginPosition + length(needle2(pattern));
        finder._state = TFinder::STATE_BEGIN_FOUND;
        pattern._state = TPattern::STATE_BEGIN_FOUND;
        return true;
    }

    finder._state = TFinder::STATE_NOTFOUND;
    pattern._state = TPattern::STATE_NOTFOUND;
    return false;
}


template <typename THaystack, typename TNeedle>
bool findBegin(Finder2<THaystack, void> & finder,  // TODO(holtgrew): "Default" better than void?
               Pattern2<TNeedle, Simple> & pattern) {
    SEQAN_CHECKPOINT;
    // State of finder and pattern should be in sync.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);
    typedef Pattern2<TNeedle, Simple> TPattern;
    return finder._state == TPattern::STATE_BEGIN_FOUND;
}


template <typename THaystack, typename TNeedle, typename TAlignSeq, typename TAlignSpec>
bool getAlignment(Finder2<THaystack, void> &finder,  // TODO(holtgrew): "Default" better than void?
                  Pattern2<TNeedle, Simple> &pattern,
                  Align<TAlignSeq, TAlignSpec> &outAlignment) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return false;
}

}  // namespace seqan
          
#endif  // SEQAN_FIND2_FIND_SIMPLE_H_
