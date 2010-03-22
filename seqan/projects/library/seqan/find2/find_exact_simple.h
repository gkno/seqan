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
struct Pattern2<TNeedle, Simple> {
    Holder<TNeedle> _host;

    Pattern2() {}
};


template <typename TNeedle>
TNeedle const & host(Pattern2<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
TNeedle const & needle(Pattern2<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
typename Position<TNeedle>::Type length(Pattern2<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
Segment<TNeedle, InfixSegment> infix(Pattern2<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
typename Iterator<TNeedle>::Type begin(Pattern2<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
typename Iterator<TNeedle>::Type end(Pattern2<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
typename Position<TNeedle>::Type beginPosition(Pattern2<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename TNeedle>
typename Position<TNeedle>::Type endPosition(Pattern2<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename THaystack, typename TNeedle>
bool find(Finder2<THaystack, void> & finder,  // TODO(holtgrew): "Default" better than void?
          Pattern2<TNeedle, Simple> & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename THaystack, typename TNeedle>
bool findBegin(Finder2<THaystack, void> & finder,  // TODO(holtgrew): "Default" better than void?
               Pattern2<TNeedle, Simple> & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename THaystack, typename TNeedle, typename TAlignSeq, typename TAlignSpec>
bool alignment(Finder2<THaystack, void> &finder,  // TODO(holtgrew): "Default" better than void?
               Pattern2<TNeedle, Simple> &pattern,
               Align<TAlignSeq, TAlignSpec> &outAlignment) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}

}  // namespace seqan
          
#endif  // SEQAN_FIND2_FIND_SIMPLE_H_
