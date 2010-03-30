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


template <typename _TNeedleContainer>
struct Pattern<_TNeedleContainer, MultipleSimple> : _FindState {
    typedef _TNeedleContainer TNeedleContainer;
    typedef typename Position<TNeedleContainer>::Type TNeedleIndex;
    typedef typename Value<TNeedleContainer>::Type TNeedle;

    // The needle set we work on.
    Holder<TNeedleContainer> _host;

    // The index of the current needle to search for.
    TNeedleIndex _needleIndex;

    // The pattern's state.
    TState _state;

    Pattern() : _state(STATE_EMPTY) { SEQAN_CHECKPOINT; }

    Pattern(TNeedleContainer & needles)
        : _host(needles),
          _state(STATE_INITIAL)
    { SEQAN_CHECKPOINT; }
};


template <typename TNeedleContainer>
TNeedleContainer const & host(Pattern<TNeedleContainer, MultipleSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}


// TODO(holtgrew): Workaround for default host implementation for non-const type.
template <typename TNeedleContainer>
TNeedleContainer const & host(Pattern<TNeedleContainer, MultipleSimple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedleContainer, MultipleSimple> TPattern;
    return host(const_cast<TPattern const &>(pattern));
}


template <typename TNeedleContainer>
TNeedleContainer const & needles(Pattern<TNeedleContainer, MultipleSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}


template <typename TNeedleContainer>
typename Position<TNeedleContainer>::Type needleIndex(Pattern<TNeedleContainer, MultipleSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedleContainer, MultipleSimple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return pattern._needleIndex;
}


template <typename TNeedleContainer>
typename Position<TNeedleContainer>::Type length(Pattern<TNeedleContainer, MultipleSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedleContainer, MultipleSimple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return length(needles(pattern)[pattern._needleIndex]);
}


template <typename TNeedleContainer>
Segment<typename TNeedleContainer::TNeedle const, InfixSegment> infix(Pattern<TNeedleContainer, MultipleSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedleContainer, MultipleSimple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return infix(needles(pattern)[needleIndex(pattern)], 0, length(pattern));
}


// TODO(holtgrew): Workaround for default infix implementation for non-const type.
template <typename TNeedleContainer>
Segment<typename TNeedleContainer::TNeedle const, InfixSegment> infix(Pattern<TNeedleContainer, MultipleSimple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedleContainer, MultipleSimple> TPattern;
    return infix(const_cast<TPattern const &>(pattern));
}


template <typename TNeedleContainer, typename TTag>
typename Iterator<typename TNeedleContainer::TNeedle const, Tag<TTag> const>::Type begin(Pattern<TNeedleContainer, MultipleSimple> const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


// TODO(holtgrew): Workaround for default begin implementation for non-const type.
template <typename TNeedleContainer, typename TTag>
typename Iterator<typename TNeedleContainer::TNeedle const, Tag<TTag> const>::Type begin(Pattern<TNeedleContainer, MultipleSimple> & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedleContainer, MultipleSimple> TPattern;
    return begin(const_cast<TPattern const &>(pattern), spec);
}


template <typename TNeedleContainer, typename TTag>
typename Iterator<typename TNeedleContainer::TNeedle const, Tag<TTag> const>::Type end(Pattern<TNeedleContainer, MultipleSimple> const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


// TODO(holtgrew): Workaround for default end implementation for non-const type.
template <typename TNeedleContainer, typename TTag>
typename Iterator<typename TNeedleContainer::TNeedle const, Tag<TTag> const>::Type end(Pattern<TNeedleContainer, MultipleSimple> & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedleContainer, MultipleSimple> TPattern;
    return end(const_cast<TPattern const &>(pattern), spec);
}


template <typename TNeedleContainer>
typename Position<typename TNeedleContainer::TNeedle>::Type beginPosition(Pattern<TNeedleContainer, MultipleSimple> const &) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


// TODO(holtgrew): Workaround for default beginPosition implementation for non-const type.
template <typename TNeedleContainer>
typename Position<typename TNeedleContainer::TNeedle>::Type beginPosition(Pattern<TNeedleContainer, MultipleSimple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedleContainer, MultipleSimple> TPattern;
    return beginPosition(const_cast<TPattern const &>(pattern));
}


template <typename TNeedleContainer>
typename Position<typename TNeedleContainer::TNeedle>::Type endPosition(Pattern<TNeedleContainer, MultipleSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


// TODO(holtgrew): Workaround for default endPosition implementation for non-const type.
template <typename TNeedleContainer>
typename Position<typename TNeedleContainer::TNeedle>::Type endPosition(Pattern<TNeedleContainer, MultipleSimple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedleContainer, MultipleSimple> TPattern;
    return beginPosition(const_cast<TPattern const &>(pattern));
}


template <typename THaystack, typename TNeedleContainer>
bool find(Finder<THaystack, Default> & finder,
          Pattern<TNeedleContainer, MultipleSimple> & pattern) {
    (void) finder;
    (void) pattern;
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return false;
}


template <typename THaystack, typename TNeedleContainer>
bool findBegin(Finder<THaystack, Default> & finder,
               Pattern<TNeedleContainer, MultipleSimple> & pattern) {
    (void) finder;
    (void) pattern;
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return false;
}


template <typename THaystack, typename TNeedleContainer, typename TPosition>
bool setEndPosition(Finder<THaystack, Default> & finder,
                    Pattern<TNeedleContainer, MultipleSimple> & pattern,
                    TPosition const & pos) {
    (void) finder;
    (void) pattern;
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return false;
}


/*
  Build the alignment resulting from the search result as specified by the
  finder and the pattern.  If the state is not "begin found" then no alignment
  is built and false is returned.
*/
template <typename THaystack, typename TNeedleContainer, typename TAlignSeq, typename TAlignSpec>
bool buildAlignment(Finder<THaystack, Default> &finder,
                    Pattern<TNeedleContainer, MultipleSimple> &pattern,
                    Align<TAlignSeq, TAlignSpec> &outAlignment) {
    (void) finder;
    (void) pattern;
    (void) outAlignment;
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return false;
}

}  // namespace seqan
          
#endif  // SEQAN_FIND2_FIND_HAMMING_SIMPLE_H_
