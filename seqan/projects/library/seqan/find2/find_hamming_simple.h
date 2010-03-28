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
  Approximate string matching with Hamming distance (also known as "k-Mismatch
  Problem") using a naive algorithm.
 ==========================================================================*/

#ifndef SEQAN_FIND2_FIND_HAMMING_SIMPLE_H_
#define SEQAN_FIND2_FIND_HAMMING_SIMPLE_H_

namespace seqan {

struct _HammingSimple;
typedef Tag<_HammingSimple> HammingSimple;
    

template <typename TNeedle>
struct Pattern<TNeedle, HammingSimple> {
    typedef int TScoreValue;

    // The needle we work on.
    Holder<TNeedle> _host;

    // The minimal score of a match.
    TScoreValue _scoreLimit;

    // The current score of a match.
    TScoreValue _currentScore;

    Pattern() {}
};


template <typename TNeedle>
int score(Pattern<TNeedle, HammingSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
int getScoreLimit(Pattern<TNeedle, HammingSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}

template <typename TNeedle>
TNeedle const & host(Pattern<TNeedle, HammingSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
TNeedle const & needle(Pattern<TNeedle, HammingSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
typename Position<TNeedle>::Type length(Pattern<TNeedle, HammingSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
Segment<TNeedle, InfixSegment> infix(Pattern<TNeedle, HammingSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
typename Iterator<TNeedle>::Type begin(Pattern<TNeedle, HammingSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
typename Iterator<TNeedle>::Type end(Pattern<TNeedle, HammingSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
typename Position<TNeedle>::Type beginPosition(Pattern<TNeedle, HammingSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
typename Position<TNeedle>::Type endPosition(Pattern<TNeedle, HammingSimple> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename THaystack, typename TNeedle>
bool find(Finder<THaystack, void> & finder,  // TODO(holtgrew): "Default" better than void?
          Pattern<TNeedle, HammingSimple> & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return false;
}


template <typename THaystack, typename TNeedle>
bool findBegin(Finder<THaystack, void> & finder,  // TODO(holtgrew): "Default" better than void?
               Pattern<TNeedle, HammingSimple> & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return false;
}


template <typename THaystack, typename TNeedle, typename TAlignSeq, typename TAlignSpec>
bool getAlignment(Finder<THaystack, void> &finder,  // TODO(holtgrew): "Default" better than void?
                  Pattern<TNeedle, HammingSimple> &pattern,
                  Align<TAlignSeq, TAlignSpec> &outAlignment) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return false;
}

}  // namespace seqan
          
#endif  // SEQAN_FIND2_FIND_HAMMING_SIMPLE_H_
