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
  Approximate string matching with Hamming distance (also known as "k-Mismatch
  Problem") using a naive algorithm.
 ==========================================================================*/

#ifndef SEQAN_FIND2_FIND_HAMMING_SIMPLE_H_
#define SEQAN_FIND2_FIND_HAMMING_SIMPLE_H_

namespace seqan {

struct HammingSimple_;
typedef Tag<HammingSimple_> HammingSimple;
    
template <typename _TNeedle>
struct Pattern<_TNeedle, HammingSimple> : _FindState {
    typedef _TNeedle TNeedle;

    // The pattern's state.
    TState _state;

    // The needle we store.
    Holder<TNeedle> _host;

    // The lowest score to return hits for.
    int _scoreLimit;
    
    // The score of the needle at the current position.
    int _score;
    
    Pattern() : _state(STATE_EMPTY) {
        SEQAN_CHECKPOINT;
    }

    Pattern(TNeedle & ndl, int scoreLimit)
        : _state(STATE_INITIAL),
          _host(ndl),
          _scoreLimit(scoreLimit) {
        SEQAN_CHECKPOINT;
        }
};


template <typename TNeedle>
struct Needle<Pattern<TNeedle, HammingSimple> > {
    typedef typename Value<TNeedle>::Type Value;
};


template <typename TNeedle>
int getScoreLimit(Pattern<TNeedle, HammingSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    return pattern._scoreLimit;
}


template <typename TNeedle>
int score(Pattern<TNeedle, HammingSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, HammingSimple> TPattern;
    // State of pattern should be in "found", "found begin" or "found
    // no begin" state.
    SEQAN_ASSERT_TRUE(pattern._state == TPattern::STATE_FOUND ||
                      pattern._state == TPattern::STATE_BEGIN_FOUND ||
                      pattern._state == TPattern::STATE_BEGIN_NOTFOUND);
    return pattern._score;
}


template <typename TNeedle>
TNeedle const & host(Pattern<TNeedle, HammingSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}


// TODO(holtgrew): Workaround for default host implementation for non-const type.
template <typename TNeedle>
TNeedle const & host(Pattern<TNeedle, HammingSimple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, HammingSimple> TPattern;
    return host(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle>
TNeedle const & needle(Pattern<TNeedle, HammingSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}


template <typename TNeedle>
typename Position<TNeedle>::Type length(Pattern<TNeedle, HammingSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, HammingSimple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return length(needle(pattern));
}


template <typename TNeedle>
Segment<TNeedle const, InfixSegment> infix(Pattern<TNeedle, HammingSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, HammingSimple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return infix(needle(pattern), 0u, length(needle(pattern)));
}


// TODO(holtgrew): Workaround for default infix implementation for non-const type.
template <typename TNeedle>
Segment<TNeedle const, InfixSegment> infix(Pattern<TNeedle, HammingSimple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, HammingSimple> TPattern;
    return infix(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type begin(Pattern<TNeedle, HammingSimple> const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, HammingSimple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return begin(needle(pattern), spec);
}


// TODO(holtgrew): Workaround for default begin implementation for non-const type.
template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type begin(Pattern<TNeedle, HammingSimple> & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, HammingSimple> TPattern;
    return begin(const_cast<TPattern const &>(pattern), spec);
}


template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type end(Pattern<TNeedle, HammingSimple> const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, HammingSimple> TPattern;
    // State of pattern should be in "found", "found begin" or "found
    // no begin" state.
    SEQAN_ASSERT_TRUE(pattern._state == TPattern::STATE_FOUND ||
                      pattern._state == TPattern::STATE_BEGIN_FOUND ||
                      pattern._state == TPattern::STATE_BEGIN_NOTFOUND);
    return end(needle(pattern), spec);
}


// TODO(holtgrew): Workaround for default end implementation for non-const type.
template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type end(Pattern<TNeedle, HammingSimple> & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, HammingSimple> TPattern;
    return end(const_cast<TPattern const &>(pattern), spec);
}


template <typename TNeedle>
typename Position<TNeedle>::Type beginPosition(Pattern<TNeedle, HammingSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, HammingSimple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return 0u;
}


// TODO(holtgrew): Workaround for default beginPosition implementation for non-const type.
template <typename TNeedle>
typename Position<TNeedle>::Type beginPosition(Pattern<TNeedle, HammingSimple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, HammingSimple> TPattern;
    return beginPosition(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle>
typename Position<TNeedle>::Type endPosition(Pattern<TNeedle, HammingSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, HammingSimple> TPattern;
    // State of pattern should be in "found", "found begin" or "found
    // no begin" state.
    SEQAN_ASSERT_TRUE(pattern._state == TPattern::STATE_FOUND ||
                      pattern._state == TPattern::STATE_BEGIN_FOUND ||
                      pattern._state == TPattern::STATE_BEGIN_NOTFOUND);
    return length(needle(pattern));
}


// TODO(holtgrew): Workaround for default endPosition implementation for non-const type.
template <typename TNeedle>
typename Position<TNeedle>::Type endPosition(Pattern<TNeedle, HammingSimple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, HammingSimple> TPattern;
    return endPosition(const_cast<TPattern const &>(pattern));
}


template <typename THaystack, typename TNeedle>
bool find(Finder<THaystack, Default> & finder,
          Pattern<TNeedle, HammingSimple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, HammingSimple> TPattern;
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
    int mismatchCount = 0;
    for (TPosition i = 0u; i < length(needle(pattern));) {
        // Break out of loop if no more match is possible.
        if (finder._beginPosition >= length(haystack(finder)) - length(needle(pattern))) {
            finder._state = TFinder::STATE_NOTFOUND;
            pattern._state = TPattern::STATE_NOTFOUND;
            return false;
        }
        mismatchCount += needle(pattern)[i] != haystack(finder)[finder._beginPosition + i];
        // Otherwise, go on searching.  We count mismatches until we find too many.
        if (mismatchCount > -getScoreLimit(pattern)) {
            finder._beginPosition += 1;
            mismatchCount = 0;
            i = 0u;
            continue;
        }
        i += 1;
    }
    finder._endPosition = finder._beginPosition + length(needle(pattern));
    finder._state = TFinder::STATE_FOUND;
    pattern._state = TPattern::STATE_FOUND;
    pattern._score = -mismatchCount;
    return true;
}


template <typename THaystack, typename TNeedle>
bool findBegin(Finder<THaystack, Default> & finder,
               Pattern<TNeedle, HammingSimple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, HammingSimple> TPattern;
    // State of finder and pattern should be in sync and in "found" or
    // "found begin" state.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);
    SEQAN_ASSERT_TRUE(pattern._state == TPattern::STATE_FOUND ||
                      pattern._state == TPattern::STATE_BEGIN_FOUND);
    if (pattern._state == TPattern::STATE_FOUND) {
        finder._state = TFinder::STATE_BEGIN_FOUND;
        pattern._state = TPattern::STATE_BEGIN_FOUND;
        return true;
    }
    finder._state = TFinder::STATE_BEGIN_NOTFOUND;
    pattern._state = TPattern::STATE_BEGIN_NOTFOUND;
    return false;
}


template <typename THaystack, typename TNeedle, typename TPosition>
bool setEndPosition(Finder<THaystack, Default> & finder,
                    Pattern<TNeedle, HammingSimple> & pattern,
                    TPosition const & pos) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, HammingSimple> TPattern;
    // State of finder and pattern should be in sync.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);
    // End position must not be right of the end of the haystack.
    SEQAN_ASSERT_LEQ(static_cast<typename MakeUnsigned_<TPosition>::Type>(pos), length(haystack(finder)));
    // Begin position must not be left of the beginning of the haystack.
    SEQAN_ASSERT_GEQ(static_cast<typename MakeUnsigned_<TPosition>::Type>(pos), length(needle(pattern)));

    // Set the end position.
    finder._endPosition = pos;
    finder._beginPosition = pos - length(needle(pattern));

    // Check whether there is a hit at this position and update the
    // state accordingly.
    int mismatchCount = 0;
    typedef typename Position<THaystack>::Type THaystackPos;
    for (THaystackPos i = 0u; i < length(needle(pattern)); ++i) {
        mismatchCount += (needle(pattern)[i] != haystack(finder)[finder._beginPosition + i]);
        if (mismatchCount > -getScoreLimit(pattern)) {
            finder._state = TPattern::STATE_NO_HIT;
            pattern._state = TFinder::STATE_NO_HIT;
            return false;
        }
    }
    finder._state = TPattern::STATE_FOUND;
    pattern._state = TPattern::STATE_FOUND;
    pattern._score = -mismatchCount;
    return true;
}


template <typename THaystack, typename TNeedle, typename TPosition>
bool setBeginPosition(Finder<THaystack, Default> & finder,
                      Pattern<TNeedle, HammingSimple> & pattern,
                      TPosition const & pos) {
    SEQAN_CHECKPOINT;
    return setEndPosition(finder, pattern, pos + length(needle(pattern)));
}


// Build the alignment resulting from the search result as specified by the
// finder and the pattern.  If the state is not "begin found" then no alignment
// is built and false is returned.
template <typename THaystack, typename TNeedle, typename TAlignSeq, typename TAlignSpec>
bool buildAlignment(Finder<THaystack, Default> &finder,
                    Pattern<TNeedle, HammingSimple> &pattern,
                    Align<TAlignSeq, TAlignSpec> &outAlignment) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, HammingSimple> TPattern;
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
    insertGaps(row(outAlignment, 1), 0, beginPosition(finder));

    return true;
}

}  // namespace seqan
          
#endif  // SEQAN_FIND2_FIND_HAMMING_SIMPLE_H_
