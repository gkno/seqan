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
struct Needle<Pattern<TNeedle, Simple> > {
    typedef typename Value<TNeedle>::Type Value;
};


template <typename TNeedle>
TNeedle const & host(Pattern<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}


// TODO(holtgrew): Workaround for default host implementation for non-const type.
template <typename TNeedle>
TNeedle const & host(Pattern<TNeedle, Simple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    return host(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle>
TNeedle const & needle(Pattern<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}


template <typename TNeedle>
typename Position<TNeedle>::Type length(Pattern<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return length(needle(pattern));
}


template <typename TNeedle>
Segment<TNeedle const, InfixSegment> infix(Pattern<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return infix(needle(pattern), 0u, length(needle(pattern)));
}


// TODO(holtgrew): Workaround for default infix implementation for non-const type.
template <typename TNeedle>
Segment<TNeedle const, InfixSegment> infix(Pattern<TNeedle, Simple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    return infix(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type begin(Pattern<TNeedle, Simple> const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return begin(needle(pattern), spec);
}


// TODO(holtgrew): Workaround for default begin implementation for non-const type.
template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type begin(Pattern<TNeedle, Simple> & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    return begin(const_cast<TPattern const &>(pattern), spec);
}


template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type end(Pattern<TNeedle, Simple> const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return end(needle(pattern), spec);
}


// TODO(holtgrew): Workaround for default end implementation for non-const type.
template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type end(Pattern<TNeedle, Simple> & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    return end(const_cast<TPattern const &>(pattern), spec);
}


template <typename TNeedle>
typename Position<TNeedle>::Type beginPosition(Pattern<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return 0u;
}


// TODO(holtgrew): Workaround for default beginPosition implementation for non-const type.
template <typename TNeedle>
typename Position<TNeedle>::Type beginPosition(Pattern<TNeedle, Simple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    return beginPosition(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle>
typename Position<TNeedle>::Type endPosition(Pattern<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return length(needle(pattern));
}


// TODO(holtgrew): Workaround for default endPosition implementation for non-const type.
template <typename TNeedle>
typename Position<TNeedle>::Type endPosition(Pattern<TNeedle, Simple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    return endPosition(const_cast<TPattern const &>(pattern));
}


template <typename THaystack, typename TNeedle>
bool find(Finder<THaystack, Default> & finder,
          Pattern<TNeedle, Simple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, Simple> TPattern;
    typedef typename Position<TNeedle>::Type TPosition;

    // State of finder and pattern should be in sync.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);

    // Do not continue if the state is "not found".
    if (finder._state == TPattern::STATE_NOTFOUND)
        return false;
    // Initialize finder if state is "initial".  Otherwise advance at
    // least by one (if not set to of haystack with setEndPosition()).
    if (finder._state == TPattern::STATE_INITIAL) {
        finder._beginPosition = 0u;
        finder._endPosition = length(needle(pattern));
    } else if (finder._state == TPattern::STATE_NO_HIT) {
        // Only advance if not at end if set manually to a "no hit" position.
        if (finder._endPosition == length(haystack(finder)))
            return false;
        finder._beginPosition += 1;
    } else {
        finder._beginPosition += 1;
    }

    // Search the needle in the haystack naively.
    for (TPosition i = 0u; i < length(needle(pattern));) {
        // Break out of loop if no more match is possible.
        if (finder._beginPosition >= length(haystack(finder)) - length(needle(pattern))) {
            finder._state = TFinder::STATE_NOTFOUND;
            pattern._state = TPattern::STATE_NOTFOUND;
            return false;
        }
        // Otherwise, go on searching.
        if (needle(pattern)[i] != haystack(finder)[finder._beginPosition + i]) {
            finder._beginPosition += 1;
            i = 0u;
            continue;
        }
        i += 1;
    }
    finder._endPosition = finder._beginPosition + length(needle(pattern));
    finder._state = TFinder::STATE_BEGIN_FOUND;
    pattern._state = TPattern::STATE_BEGIN_FOUND;
    return true;
}


template <typename THaystack, typename TNeedle>
bool findBegin(Finder<THaystack, Default> & finder,
               Pattern<TNeedle, Simple> & pattern) {
    SEQAN_CHECKPOINT;
    // State of finder and pattern should be in sync.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);
    typedef Pattern<TNeedle, Simple> TPattern;
    return finder._state == TPattern::STATE_BEGIN_FOUND;
}


template <typename THaystack, typename TNeedle, typename TPosition>
bool setEndPosition(Finder<THaystack, Default> & finder,
                    Pattern<TNeedle, Simple> & pattern,
                    TPosition const & pos) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, Simple> TPattern;
    // State of finder and pattern should be in sync.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);
    // End position must not be right of the end of the haystack.
    SEQAN_ASSERT_LEQ(static_cast<typename _MakeUnsigned<TPosition>::Type>(pos), length(haystack(finder)));
    // Begin position must not be left of the beginning of the haystack.
    SEQAN_ASSERT_GEQ(static_cast<typename _MakeUnsigned<TPosition>::Type>(pos), length(needle(pattern)));

    // Set the end position.
    finder._endPosition = pos;
    finder._beginPosition = pos - length(needle(pattern));

    // Check whether there is a hit at this position and update the
    // state accordingly.
    typedef typename Position<THaystack>::Type THaystackPos;
    for (THaystackPos i = 0u; i < length(needle(pattern)); ++i) {
        if (needle(pattern)[i] != haystack(finder)[finder._beginPosition + i]) {
            finder._state = TPattern::STATE_NO_HIT;
            pattern._state = TFinder::STATE_NO_HIT;
            return false;
        }
    }
    finder._state = TPattern::STATE_BEGIN_FOUND;
    pattern._state = TPattern::STATE_BEGIN_FOUND;
    return true;
}


/*
  Build the alignment resulting from the search result as specified by the
  finder and the pattern.  If the state is not "begin found" then no alignment
  is built and false is returned.
*/
template <typename THaystack, typename TNeedle, typename TAlignSeq, typename TAlignSpec>
bool buildAlignment(Finder<THaystack, Default> &finder,
                    Pattern<TNeedle, Simple> &pattern,
                    Align<TAlignSeq, TAlignSpec> &outAlignment) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, Simple> TPattern;
    typedef Align<TAlignSeq, TAlignSpec> TAlign;
    typedef typename Row<TAlign>::Type TRow;

    // Both finder and pattern must be in the "found begin position" state.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);
    // Can only build alignment if the state is "begin found".
    if (finder._state != TFinder::STATE_BEGIN_FOUND)
        return false;

    // Initialize alignment with the two sequences.
    resize(rows(outAlignment), 2);
    assignSource(row(outAlignment, 0), haystack(finder));
    assignSource(row(outAlignment, 1), needle(pattern));
    // Insert gap into the needle.
    insertGaps(row(outAlignment, 1), 0, beginPosition(finder));

    return true;
}

}  // namespace seqan
          
#endif  // SEQAN_FIND2_FIND_EXACT_SIMPLE_H_
