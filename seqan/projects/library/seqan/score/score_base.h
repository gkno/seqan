// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================

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
class Score;


///.Metafunction.Value.param.T.type:Class.Score
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
