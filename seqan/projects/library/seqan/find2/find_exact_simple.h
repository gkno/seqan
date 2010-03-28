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

#ifndef SEQAN_FIND2_FIND_EXACT_SIMPLE_H_
#define SEQAN_FIND2_FIND_EXACT_SIMPLE_H_

namespace seqan {

template <typename _TNeedle>
struct Pattern<_TNeedle, Simple> : _FindState {
    typedef _TNeedle TNeedle;

    // The pattern's state.
    TState _state;

    // The needle we store.
    Holder<TNeedle> _host;

    Pattern() : _state(STATE_EMPTY) {
        SEQAN_CHECKPOINT;
    }

    explicit
    Pattern(TNeedle & ndl)
        : _state(STATE_INITIAL),
          _host(ndl) {
        SEQAN_CHECKPOINT;
    }
};


template <typename TNeedle>
TNeedle & needle(Pattern<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}


template <typename TNeedle>
typename Position<TNeedle>::Type length(Pattern<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    return length(needle(pattern));
}


template <typename TNeedle>
Segment<TNeedle, InfixSegment> infix(Pattern<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    return infix(needle(pattern), 0u, length(needle(pattern)));
}


template <typename TNeedle, typename TTag>
    typename Iterator<TNeedle const, Tag<TTag> const>::Type begin(Pattern<TNeedle, Simple> const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    return begin(needle(pattern), spec);
}


template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type end(Pattern<TNeedle, Simple> const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    return end(needle(pattern), spec);
}


template <typename TNeedle>
typename Position<TNeedle>::Type beginPosition(Pattern<TNeedle, Simple> const &) {
    SEQAN_CHECKPOINT;
    return 0u;
}


template <typename TNeedle>
typename Position<TNeedle>::Type endPosition(Pattern<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    return length(needle(pattern));
}


template <typename THaystack, typename TNeedle>
bool find(Finder<THaystack, Default> & finder,  // TODO(holtgrew): "Default" better than void?
          Pattern<TNeedle, Simple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<THaystack, Simple> TFinder;
    typedef Pattern<TNeedle, Simple> TPattern;
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
    while (finder._beginPosition < length(haystack(finder)) - length(needle(pattern))) {
        for (TPosition i = 0u; i < length(needle(pattern)); ++i) {
            if (needle(pattern)[i] != haystack(finder)[finder._beginPosition + i]) {
                finder._beginPosition += 1;
                i = 0u;
                continue;
            }
        }
        finder._endPosition = finder._beginPosition + length(needle(pattern));
        finder._state = TFinder::STATE_BEGIN_FOUND;
        pattern._state = TPattern::STATE_BEGIN_FOUND;
        return true;
    }

    finder._state = TFinder::STATE_NOTFOUND;
    pattern._state = TPattern::STATE_NOTFOUND;
    return false;
}


template <typename THaystack, typename TNeedle>
bool findBegin(Finder<THaystack, Default> & finder,  // TODO(holtgrew): "Default" better than void?
               Pattern<TNeedle, Simple> & pattern) {
    SEQAN_CHECKPOINT;
    // State of finder and pattern should be in sync.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);
    typedef Pattern<TNeedle, Simple> TPattern;
    return finder._state == TPattern::STATE_BEGIN_FOUND;
}


template <typename THaystack, typename TNeedle, typename TAlignSeq, typename TAlignSpec>
bool getAlignment(Finder<THaystack, void> &finder,  // TODO(holtgrew): "Default" better than void?
                  Pattern<TNeedle, Simple> &pattern,
                  Align<TAlignSeq, TAlignSpec> &outAlignment) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return false;
}

}  // namespace seqan
          
#endif  // SEQAN_FIND2_FIND_EXACT_SIMPLE_H_
