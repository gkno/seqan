/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 207-010

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
  Exact pattern matching for multiple needles at once with a naive
  implementation.
 ==========================================================================*/

#ifndef SEQAN_FIND2_FIND_HAMMING_SIMPLE_H_
#define SEQAN_FIND2_FIND_HAMMING_SIMPLE_H_

namespace seqan {

struct _MultipleSimple;
typedef Tag<_MultipleSimple> MultipleSimple;


template <typename _TNeedle, typename TStringSetSpec>
struct Pattern<StringSet<_TNeedle, TStringSetSpec>, MultipleSimple> : _FindState {
    typedef StringSet<_TNeedle, TStringSetSpec> TNeedleContainer;
    typedef typename Position<TNeedleContainer>::Type TNeedleIndex;
    typedef _TNeedle TNeedle;
    typedef typename Position<TNeedle>::Type TPosition;

    // The needle set we work on.
    Holder<TNeedleContainer> _host;

    // The index of the needle that was searche for last.
    TNeedleIndex _needleIndex;

    // The pattern's state.
    TState _state;

    Pattern() : _state(STATE_EMPTY) { SEQAN_CHECKPOINT; }

    Pattern(TNeedleContainer & needles)
        : _host(needles),
          _state(STATE_INITIAL)
    { SEQAN_CHECKPOINT; }
};


template <typename TNeedle, typename TStringSetSpec>
struct Needle<Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> > {
    typedef TNeedle Value;
};


template <typename TNeedle, typename TStringSetSpec>
StringSet<TNeedle, TStringSetSpec> const & host(Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}


// TODO(holtgrew): Workaround for default host implementation for non-const type.
template <typename TNeedle, typename TStringSetSpec>
StringSet<TNeedle, TStringSetSpec> const & host(Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> TPattern;
    return host(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle, typename TStringSetSpec>
StringSet<TNeedle, TStringSetSpec> const & needles(Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    return host(pattern);
}


template <typename TNeedle, typename TStringSetSpec>
typename Position<StringSet<TNeedle, TStringSetSpec> >::Type needleIndex(Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return pattern._needleIndex;
}


template <typename TNeedle, typename TStringSetSpec>
typename Position<StringSet<TNeedle, TStringSetSpec> >::Type length(Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return length(needles(pattern)[pattern._needleIndex]);
}


template <typename TNeedle, typename TStringSetSpec>
Segment<TNeedle const, InfixSegment>
infix(Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return infix(needles(pattern)[needleIndex(pattern)], 0, length(pattern));
}


// TODO(holtgrew): Workaround for default infix implementation for non-const type.
template <typename TNeedle, typename TStringSetSpec>
Segment<TNeedle const, InfixSegment>
infix(Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> TPattern;
    return infix(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle, typename TStringSetSpec, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type begin(Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<StringSet<TNeedle, TStringSetSpec>, Simple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return begin(needles(pattern)[needleIndex(pattern)], spec);
}


// TODO(holtgrew): Workaround for default begin implementation for non-const type.
template <typename TNeedle, typename TStringSetSpec, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type begin(Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> TPattern;
    return begin(const_cast<TPattern const &>(pattern), spec);
}


template <typename TNeedle, typename TStringSetSpec, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type end(Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<StringSet<TNeedle, TStringSetSpec>, Simple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return end(needles(pattern)[needleIndex(pattern)], spec);
}


// TODO(holtgrew): Workaround for default end implementation for non-const type.
template <typename TNeedle, typename TStringSetSpec, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type end(Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> TPattern;
    return end(const_cast<TPattern const &>(pattern), spec);
}


template <typename TNeedle, typename TStringSetSpec>
typename Position<TNeedle>::Type beginPosition(Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<StringSet<TNeedle, TStringSetSpec>, Simple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return 0u;
}


// TODO(holtgrew): Workaround for default beginPosition implementation for non-const type.
template <typename TNeedle, typename TStringSetSpec>
typename Position<TNeedle>::Type beginPosition(Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> TPattern;
    return beginPosition(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle, typename TStringSetSpec>
typename Position<TNeedle>::Type endPosition(Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> TPattern;
    return length(needles(pattern)[needleIndex(pattern)]);
}


// TODO(holtgrew): Workaround for default endPosition implementation for non-const type.
template <typename TNeedle, typename TStringSetSpec>
typename Position<TNeedle>::Type endPosition(Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> TPattern;
    return endPosition(const_cast<TPattern const &>(pattern));
}


template <typename THaystack, typename TNeedle, typename TStringSetSpec>
bool find(Finder<THaystack, Default> & finder,
          Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef StringSet<TNeedle, TStringSetSpec> TStringSet;
    typedef Pattern<TStringSet, MultipleSimple> TPattern;
    typedef typename TPattern::TNeedleIndex TNeedleIndex;
    typedef typename Position<TNeedle>::Type TPosition;

    // State of finder and pattern should be in sync.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);

    // Get some shortcuts.
    TNeedleIndex const kNeedleCount = length(needles(pattern));
    THaystack const & hstck = haystack(finder);
    TStringSet const & ndls = needles(pattern);
    
    // Do not continue if the state is "not found".
    if (finder._state == TFinder::STATE_NOTFOUND)
        return false;
    // Initialize finder if state is "initial".  Otherwise advance at
    // least by one (if not set to of haystack with setEndPosition()).
    if (finder._state == TFinder::STATE_INITIAL) {
        pattern._needleIndex = 0;
        finder._beginPosition = 0u;
        // The end position is set to the length of the shortest needle.
        TPosition minNeedleLength = length(ndls[0]);
        for (TNeedleIndex nidx = 0; nidx < length(ndls); ++nidx)
            minNeedleLength = _min(minNeedleLength, length(ndls[nidx]));
        finder._endPosition = minNeedleLength;
    } else if (finder._state == TFinder::STATE_NO_HIT) {
        // Only advance if not at end if set manually to a "no hit" position.
        if (finder._endPosition == length(haystack(finder)))
            return false;
        // Go to next needle if possible and handle "overflow" in the
        // needle index by incrementing the end position.
        pattern._needleIndex += 1;
        if (pattern._needleIndex >= kNeedleCount) {
            pattern._needleIndex = 0;
            finder._endPosition += 1;
        }
    } else {
        // Go to next needle if possible and handle "overflow" in the
        // needle index by incrementing the end position.
        pattern._needleIndex += 1;
        SEQAN_ASSERT_LEQ(pattern._needleIndex, kNeedleCount);
        if (pattern._needleIndex == kNeedleCount) {
            pattern._needleIndex = 0;
            finder._endPosition += 1;
        }
    }

    // Search the needle in the haystack naively.
    for (; finder._endPosition <= length(hstck); ++finder._endPosition) {
        for (; pattern._needleIndex < length(ndls); ++pattern._needleIndex) {
            if (length(ndls[pattern._needleIndex]) > finder._endPosition)
                continue;  // Needle does not fit left of current end position.
            // Compare from needle and haystack from the left to the right.
            TNeedle const & ndl = ndls[pattern._needleIndex];
            TNeedleIndex needleLength = length(ndl);
            bool match = true;
            for (TPosition i = 0u; i < length(ndl); ++i) {
                if (ndl[i] != hstck[finder._endPosition - needleLength + i]) {
                    match = false;
                    break;
                }
            }
            // On success, update the finder and return true.
            if (match) {
                finder._beginPosition = finder._endPosition - needleLength;
                finder._state = TFinder::STATE_BEGIN_FOUND;
                pattern._state = TPattern::STATE_BEGIN_FOUND;
                return true;
            }
            // If no matching location was found, this loop is left and the next
            // needle is tried.
        }
        // Start with first needle for the next end position.
        pattern._needleIndex = 0;
    }
    return false;
}


template <typename THaystack, typename TNeedle, typename TStringSetSpec>
bool findBegin(Finder<THaystack, Default> & finder,
               Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> & pattern) {
    SEQAN_CHECKPOINT;
    // State of finder and pattern should be in sync.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);
    typedef Pattern<TNeedle, Simple> TPattern;
    return finder._state == TPattern::STATE_BEGIN_FOUND;
}


template <typename THaystack, typename TNeedle, typename TStringSetSpec, typename TPos>
bool setEndPosition(Finder<THaystack, Default> & finder,
                    Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> & pattern,
                    TPos const & pos) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef StringSet<TNeedle, TStringSetSpec> TStringSet;
    typedef Pattern<TStringSet, MultipleSimple> TPattern;
    typedef typename TPattern::TNeedleIndex TNeedleIndex;
    typedef typename Position<TNeedle>::Type TPosition;

    // State of finder and pattern should be in sync.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);

    // Get some shortcuts.
    TNeedleIndex const kNeedleCount = length(needles(pattern));
    THaystack const & hstck = haystack(finder);
    TStringSet const & ndls = needles(pattern);

    // Set the end position.
    finder._endPosition = pos;

    // Find the first needle that matches here, if any.
    bool found = false;
    for (pattern._needleIndex = 0; pattern._needleIndex < kNeedleCount; ++pattern._needleIndex) {
        TNeedle const & ndl = ndls[pattern._needleIndex];
        TPosition kNeedleLength = length(ndl);
        // Skip needle if it does not fit left of current end position.
        if (kNeedleLength > finder._endPosition)
            continue;
        bool mismatch = false;
        for (TPosition i = 0; i < kNeedleLength; ++i) {
            if (ndl[i] != hstck[finder._endPosition - kNeedleLength + i]) {
                mismatch = true;
                break;
            }
        }
        // If we found a mismatch, try the next needle.
        if (mismatch)
            continue;
        // Otherwise, we can leave this loop and have found a match
        // for the needle pattern._needleIndex.
        found = true;
        break;
    }
    
    if (found) {
        finder._beginPosition = finder._endPosition - length(ndls[pattern._needleIndex]);
        finder._state = TFinder::STATE_BEGIN_FOUND;
        pattern._state = TPattern::STATE_BEGIN_FOUND;
        return true;
    } else {
        finder._state = TFinder::STATE_NO_HIT;
        pattern._state = TPattern::STATE_NO_HIT;
        return false;
    }
}


/*
  Build the alignment resulting from the search result as specified by the
  finder and the pattern.  If the state is not "begin found" then no alignment
  is built and false is returned.
*/
template <typename THaystack, typename TNeedle, typename TStringSetSpec, typename TAlignSeq, typename TAlignSpec>
bool buildAlignment(Finder<THaystack, Default> &finder,
                    Pattern<StringSet<TNeedle, TStringSetSpec>, MultipleSimple> &pattern,
                    Align<TAlignSeq, TAlignSpec> &outAlignment) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef StringSet<TNeedle, TStringSetSpec> TStringSet;
    typedef Pattern<TStringSet, Simple> TPattern;
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
    assignSource(row(outAlignment, 1), needles(pattern)[needleIndex(pattern)]);
    insertGaps(row(outAlignment, 1), 0, beginPosition(finder));

    return true;
}

}  // namespace seqan
          
#endif  // SEQAN_FIND2_FIND_HAMMING_SIMPLE_H_
