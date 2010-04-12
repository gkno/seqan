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
  Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
 ==========================================================================*/

// TODO(holtgrew): Should the public interface for the class Score not be defined here?

#ifndef SEQAN_SCORE_SCORE_BASE_H_
#define SEQAN_SCORE_SCORE_BASE_H_

namespace SEQAN_NAMESPACE_MAIN {

/**
.Class.Score:
..cat:Miscellaneous
..summary:A scoring scheme.
..signature:Score<TValue, TSpec>
..param.TValue:The value type.
...default:int
..param.TSpec:The specializing type.
...default:@Tag.Simple@
..include:seqan/score.h
*/
template <typename TValue = int, typename TSpec = Simple>
struct Score;


/**
.Metafunction.Value.param.T.type:Class.Score
 */
template <typename TValue, typename TSpec>
struct Value<Score<TValue, TSpec> > {
    typedef TValue Type;
};


/**
.Function.scoreGapOpenHorizontal
..cat:Scoring
..summary:TODO(holtgrew): David wrote this, he should comment it.
..signature:scoreGapOpenHorizontal(score, pos1, pos2, seq1, seq2)
..include:seqan/score.h
 */
template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenHorizontal(
    Score<TValue, TSpec> const & me,
    TPos1,
    TPos2,
    TSeq1 const &,
    TSeq2 const &) {
    SEQAN_CHECKPOINT;
    return scoreGapOpen(me);
}


/**
.Function.scoreGapOpenVertical
..cat:Scoring
..summary:TODO(holtgrew): David wrote this, he should comment it.
..signature:scoreGapOpenVertical(score, pos1, pos2, seq1, seq2)
..include:seqan/score.h
 */
template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenVertical(
    Score<TValue, TSpec> const & me,
    TPos1,
    TPos2,
    TSeq1 const &,
    TSeq2 const &) {
    SEQAN_CHECKPOINT;
    return scoreGapOpen(me);
}


/**
.Function.scoreGapExtendHorizontal
..cat:Scoring
..summary:TODO(holtgrew): David wrote this, he should comment it.
..signature:scoreGapExtendHorizontal(score, pos1, pos2, seq1, seq2)
..include:seqan/score.h
 */
template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
    Score<TValue, TSpec> const & me,
    TPos1,
    TPos2,
    TSeq1 const &,
    TSeq2 const &) {
    SEQAN_CHECKPOINT;
    return scoreGapExtend(me);
}


/**
.Function.scoreGapExtendVertical
..cat:Scoring
..summary:TODO(holtgrew): David wrote this, he should comment it.
..signature:scoreGapExtendVertical(score, pos1, pos2, seq1, seq2)
..include:seqan/score.h
 */
template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
    Score<TValue, TSpec> const & me,
    TPos1,
    TPos2,
    TSeq1 const &,
    TSeq2 const &) {
    SEQAN_CHECKPOINT;
    return scoreGapExtend(me);
}


/**
.Function.scoreGapHorizontal
..cat:Scoring
..summary:TODO(holtgrew): David wrote this, he should comment it.
..signature:scoreGapHorizontal(score, pos1, pos2, seq1, seq2)
..include:seqan/score.h
 */
template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapHorizontal(
    Score<TValue, TSpec> const & me,
    TPos1,
    TPos2,
    TSeq1 const &,
    TSeq2 const &) {
    SEQAN_CHECKPOINT;
    return scoreGap(me);
}


/**
.Function.scoreGapVertical
..cat:Scoring
..summary:TODO(holtgrew): David wrote this, he should comment it.
..signature:scoreGapVertical(score, pos1, pos2, seq1, seq2)
..include:seqan/score.h
 */
template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapVertical(
    Score<TValue, TSpec> const & me,
    TPos1,
    TPos2,
    TSeq1 const &,
    TSeq2 const &) {
    SEQAN_CHECKPOINT;
    return scoreGap(me);
}


/**
.Function.score
..cat:Scoring
..signature:score(score, pos1, pos2, seq1, seq2)
..remark:ATTENTION -- score(TScore, TVal1, TVal2) is deprecated.  Better use this function.
..remark:TODO(holtgrew): David wrote this, he should comment it.
..include:seqan/score.h
 */
template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, TSpec> const & me,
      TPos1 pos1,
      TPos2 pos2,
      TSeq1 const &seq1,
      TSeq2 const &seq2) {
    SEQAN_CHECKPOINT;
    return score(me, seq1[pos1], seq2[pos2]);
}

}  // namespace SEQAN_NAMESPACE_MAIN

#endif  // SEQAN_SCORE_SCORE_BASE_H_
