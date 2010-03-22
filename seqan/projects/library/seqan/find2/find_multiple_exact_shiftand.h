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
  Exact pattern matching for multiple needles at once with a variant of the
  Shift-And Algorithm.
 ==========================================================================*/

#ifndef SEQAN_FIND2_FIND_HAMMING_SIMPLE_H_
#define SEQAN_FIND2_FIND_HAMMING_SIMPLE_H_

namespace seqan {

template <typename TNeedleContainer>
struct Pattern2<TNeedleContainer, MultipleShiftAnd> {
    typedef int TScoreValue;

    // The needle set we work on.
    Holder<TNeedleContainer> _host;

    Pattern2() {}
};


template <typename TNeedleSet>
TNeedleSet const & host(Pattern2<TNeedleSet, MultipleShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedleSet>
TNeedleSet const & needles(Pattern2<TNeedleSet, MultipleShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedleSet>
typename Position<TNeedleSet>::Type needleIdentifier(Pattern2<TNeedleSet, MultipleShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedleSet>
typename Position<typename Value<TNeedleSet>::Type>::Type length(Pattern2<TNeedleSet, MultipleShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedleSet>
Segment<typename Value<TNeedleSet>::Type, InfixSegment> infix(Pattern2<TNeedleSet, MultipleShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedleSet>
typename Iterator<typename Value<TNeedleSet>::Type>::Type begin(Pattern2<TNeedleSet, MultipleShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedleSet>
typename Iterator<typename Value<TNeedleSet>::Type>::Type end(Pattern2<TNeedleSet, MultipleShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedleSet>
typename Position<typename Value<TNeedleSet>::Type>::Type beginPosition(Pattern2<TNeedleSet, MultipleShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedleSet>
typename Position<typename Value<TNeedleSet>::Type>::Type endPosition(Pattern2<TNeedleSet, MultipleShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename THaystack, typename TNeedleSet>
bool find(Finder2<THaystack, void> & finder,  // TODO(holtgrew): "Default" better than void?
          Pattern2<TNeedleSet, MultipleShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return false;
}


template <typename THaystack, typename TNeedleSet>
bool findBegin(Finder2<THaystack, void> & finder,  // TODO(holtgrew): "Default" better than void?
               Pattern2<TNeedleSet, MultipleShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return false;
}


template <typename THaystack, typename TNeedleSet, typename TAlignSeq, typename TAlignSpec>
bool getAlignment(Finder2<THaystack, void> &finder,  // TODO(holtgrew): "Default" better than void?
                  Pattern2<TNeedleSet, MultipleShiftAnd> &pattern,
                  Align<TAlignSeq, TAlignSpec> &outAlignment) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return false;
}

}  // namespace seqan
          
#endif  // SEQAN_FIND2_FIND_HAMMING_SIMPLE_H_
