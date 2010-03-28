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
  Approximate string matching with linear and affine gap costs based on
  Seller's algorithm, Gotoh's.
 ==========================================================================*/

#ifndef SEQAN_FIND2_FIND_APPROX_DPSEARCH_H_
#define SEQAN_FIND2_FIND_APPROX_DPSEARCH_H_

namespace seqan {

struct _FindInfix;
typedef Tag<_FindInfix> FindInfix;

struct _FindPrefix;
typedef Tag<_FindPrefix> FindPrefix;

template <typename TScore, typename TSpec = FindInfix, typename TFindBeginPatternSpec = True>
struct DPSearch;


// TODO(holtgrew): Rename TScore to TScoringScheme?  Also: Rename Score class to ScoringScheme?
template <typename TNeedle, typename TScore, typename TSpec, typename TSupportFindBegin>
struct Pattern<TNeedle, DPSearch<TScore, TSpec, TSupportFindBegin> > {
    typedef typename Value<TScore>::Value TScoreValue;

    // The needle we work on.
    Holder<TNeedle> _host;

    // The scoring scheme to use.
    TScore _scoringScheme;

    // The minimal score of a match.
    TScoreValue _scoreLimit;

    // The current score of a match.
    TScoreValue _currentScore;

    Pattern() {}
};


template <typename THaystack, typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
int score(Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return pattern._score;
}


template <typename THaystack, typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
int getScoreLimit(Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return pattern._scoreLimit;
}


template <typename THaystack, typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
TScore getScoringScheme(Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> >  const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename THaystack, typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
TScore const & scoringScheme(Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> >  const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename THaystack, typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
TNeedle const & host(Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename THaystack, typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
TNeedle const & needle(Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename THaystack, typename TNeedle, typename TScore, typename TSpec>
typename Position<TNeedle>::Type length(Pattern<TNeedle, DPSearch<TScore, TSpec, True> > const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename THaystack, typename TNeedle, typename TScore, typename TSpec>
Segment<TNeedle, InfixSegment> infix(Pattern<TNeedle, DPSearch<TScore, TSpec, True> > const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename THaystack, typename TNeedle, typename TScore, typename TSpec>
typename Iterator<TNeedle>::Type begin(Pattern<TNeedle, DPSearch<TScore, TSpec, True> > const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename THaystack, typename TNeedle, typename TScore, typename TSpec>
typename Iterator<TNeedle>::Type end(Pattern<TNeedle, DPSearch<TScore, TSpec, True> > const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename THaystack, typename TNeedle, typename TScore, typename TSpec>
typename Position<TNeedle>::Type beginPosition(Pattern<TNeedle, DPSearch<TScore, TSpec, True> > const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
}


template <typename THaystack, typename TNeedle, typename TScore, typename TSpec>
typename Position<TNeedle>::Type endPosition(Pattern<TNeedle, DPSearch<TScore, TSpec, True> > const & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return false;
}


template <typename THaystack, typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
bool find(Finder<THaystack, void> & finder,  // TODO(holtgrew): "Default" better than void?
          Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return false;
}


// The function "findBegin" is only supported if the TFindBeginPatternSpec is True.
template <typename THaystack, typename TNeedle, typename TScore, typename TSpec>
bool findBegin(Finder<THaystack, void> & finder,  // TODO(holtgrew): "Default" better than void?
               Pattern<TNeedle, DPSearch<TScore, TSpec, True> > & pattern) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return false;
}


// The function "alignment" is only supported if the TFindBeginPatternSpec is True.
template <typename THaystack, typename TNeedle, typename TAlignSeq, typename TAlignSpec, typename TScore, typename TSpec>
bool getAlignment(Finder<THaystack, void> &finder,  // TODO(holtgrew): "Default" better than void?
                  Pattern<TNeedle, DPSearch<TScore, TSpec, True> > & pattern,
                  Align<TAlignSeq, TAlignSpec> & outAlignment) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return false;
}

}  // namespace seqan
          
#endif  // SEQAN_FIND2_FIND_APPROX_DPSEARCH_H_
