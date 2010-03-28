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
  Wildcard pattern matching using a modification of the Shift-And Algorithm.
 ==========================================================================*/

#ifndef SEQAN_FIND2_FIND_PATTERN_WILD_SHIFTAND_H_
#define SEQAN_FIND2_FIND_PATTERN_WILD_SHIFTAND_H_

namespace seqan {

struct _WildShiftAnd;
typedef Tag<_WildShiftAnd> WildShiftAnd;

template <typename TNeedle>
struct Pattern<TNeedle, WildShiftAnd> {
    Holder<TNeedle> _host;

    Pattern() {}
};


template <typename TNeedle>
TNeedle const & host(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
TNeedle const & needle(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
typename Position<TNeedle>::Type length(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
Segment<TNeedle, InfixSegment> infix(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
typename Iterator<TNeedle>::Type begin(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
typename Iterator<TNeedle>::Type end(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
typename Position<TNeedle>::Type beginPosition(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
typename Position<TNeedle>::Type endPosition(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename THaystack, typename TNeedle>
bool find(Finder<THaystack, void> & finder,  // TODO(holtgrew): "Default" better than void?
          Pattern<TNeedle, WildShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return false;
}


template <typename THaystack, typename TNeedle>
bool findBegin(Finder<THaystack, void> & finder,  // TODO(holtgrew): "Default" better than void?
               Pattern<TNeedle, WildShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return false;
}


template <typename THaystack, typename TNeedle, typename TAlignSeq, typename TAlignSpec>
bool getAlignment(Finder<THaystack, void> &finder,  // TODO(holtgrew): "Default" better than void?
                  Pattern<TNeedle, WildShiftAnd> &pattern,
                  Align<TAlignSeq, TAlignSpec> &outAlignment) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return false;
}

}  // namespace seqan
          
#endif  // SEQAN_FIND2_FIND_PATTERN_WILD_SHIFTAND_H_
