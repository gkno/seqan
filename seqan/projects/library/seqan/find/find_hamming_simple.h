/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  A finder that implements a simple/naive/brute force approach to approximate
  string matching with hamming distance.
 ============================================================================
  Status: Testing.  Should work, some performance improvements possible,
          see TODOs.
 ============================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
 ==========================================================================*/

#ifndef SEQAN_FIND_FIND_SIMPLE_H_
#define SEQAN_FIND_FIND_SIMPLE_H_

#include <algorithm>

namespace SEQAN_NAMESPACE_MAIN {

/**
.Spec.HammingSimpleFinder:
..summary:A brute force online searching algorithm for approximate string matching with hamming distance.
..general:Class.Pattern
..cat:Searching
..signature:Pattern<TNeedle, HammingSimple>
..param.TNeedle:The needle type.
...type:Class.String
..remarks:This specialization should only be used if no other is applicable or for verification purposes.
..include:seqan/find.h
*/

///.Class.Pattern.param.TSpec.type:Spec.HammingSimpleFinder

struct _HammingSimple;
typedef Tag<_HammingSimple> HammingSimple;


template <typename TNeedle>
class Pattern<TNeedle, HammingSimple> {
public:
    // The holder for the needle.
    Holder<TNeedle> data_host;

    // The maximal distance.  Must be >= 0, i.e. -score.
    int maxDistance;

    // The current distance, >= 0, i.e. -current score.
    int distance;

    Pattern() : maxDistance(-1), distance(0) {}

    template <typename TNeedle2>
    Pattern(const TNeedle2 &ndl, int k = -1) {
        SEQAN_CHECKPOINT;
        setHost(*this, ndl, k);
    }
};


template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, HammingSimple> & me, 
              const TNeedle2 & needle, int k) {
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_NOT(empty(needle));
    SEQAN_ASSERT_LEQ_MSG(k, 0, "Are you confusing distances and scores?");

    setValue(me.data_host, needle);
    me.maxDistance = -k;
}


template <typename TNeedle, typename TNeedle2>
void
setHost(Pattern<TNeedle, HammingSimple> &horsp, TNeedle2 &ndl, int k) {
    SEQAN_CHECKPOINT;
    setHost(horsp, reinterpret_cast<const TNeedle2&>(ndl), k);
}


template <typename TNeedle>
inline void _finderInit(Pattern<TNeedle, HammingSimple> & me) {
    SEQAN_CHECKPOINT;
    (void) me;  // Suppress unused variable warning.
}


template <typename TNeedle>
inline int score(const Pattern<TNeedle, HammingSimple> &me) {
    SEQAN_CHECKPOINT;
    return -me.distance;
}


template <typename TNeedle>
inline int getScore(const Pattern<TNeedle, HammingSimple> &me) {
    SEQAN_CHECKPOINT;
    return -me.distance;
}


template <typename TNeedle>
inline void setScoreLimit(Pattern<TNeedle, HammingSimple> & me, int _limit) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_LEQ(_limit, 0);
    me.maxDistance = -_limit;
}


template <typename TFinder, typename TNeedle>
inline bool find(TFinder &finder, 
                 Pattern<TNeedle, HammingSimple> &me) {
    SEQAN_CHECKPOINT;

    typedef typename Haystack<TFinder>::Type THaystack;
    typedef typename Size<THaystack>::Type TSize;

    // Shortcuts to haystack and needle.
    const THaystack &hstk = haystack(finder);
    const TNeedle &ndl = needle(me);

    // If the needle is longer than the haystack then we cannot find anything.
    if (length(hstk) < length(ndl))
        return false;

    // Initialize or advance finder, depending whether it has been
    // initialized before.
    if (empty(finder)) {
        _finderInit(me);
        _setFinderLength(finder, length(needle(me)));
        _finderSetNonEmpty(finder);
    } else {
        finder += 1;
    }

    // Check whether we are beyond the last possible match position.
    if (position(finder) > length(hstk) - length(ndl))
        return false;

    // TODO(holtgrew): Switch from indices to iterators to improve performance.

    // Perform a naive search for the needle in the haystack such that
    // the difference is <= me.maxDistance.
    TSize i;
    for (i = position(finder); i <= length(hstk) - length(ndl); ++i) {
        me.distance = 0;  // Reset mismatch count.
        for (TSize j = 0; j < length(ndl); ++j) {
            me.distance += (ndl[j] != hstk[i + j]);
            if (me.distance > me.maxDistance)
                break;
        }
        if (me.distance <= me.maxDistance)
            break;
    }

    // Return false if we did not break out of the for-loop but it
    // stopped normally.
    if (i > length(hstk) - length(ndl))
        return false;

    _setFinderEnd(finder, i + length(ndl));
    setPosition(finder, beginPosition(finder));
    return true; 
}

}  // namespace SEQAN_NAMESPACE_MAIN

#endif  // SEQAN_FIND_FIND_SIMPLE_H_
