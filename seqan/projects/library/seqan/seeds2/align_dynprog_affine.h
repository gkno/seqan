/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2010

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
 ==========================================================================*/

#ifndef SEQAN_SEEDS_ALIGN_DYNPROG_AFFINE_H_
#define SEQAN_SEEDS_ALIGN_DYNPROG_AFFINE_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

template <typename TScoreValue, typename TSequence>
inline void
_align_resizeMatrix(Matrix<TScoreValue, 2> & matrix, TSequence const & sequence0, TSequence const & sequence1, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_FAIL("Not implemented!");
}


template <typename TScoreValue, bool BEGIN1_FREE, bool BEGIN0_FREE, bool END1_FREE, bool END0_FREE>
inline void
_align_initGutter(Matrix<TScoreValue, 2> & matrix, Score<TScoreValue, Simple> const scoringScheme, AlignConfig<BEGIN1_FREE, BEGIN0_FREE, END1_FREE, END0_FREE> const &, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Iterator<TMatrix>::Type TIterator;
    typedef typename Position<TMatrix>::Type TPosition;

    SEQAN_ASSERT_EQ_MSG(scoreGapOpen(scoringScheme), scoreGapExtend(scoringScheme),
                        "Only linear gap costs allowed for Needleman-Wunsch.");

    // Init left gutter with zeroes if begin gaps are in dimension 0
    // free or with gap scores otherwise.
    if (BEGIN0_FREE) {
        TIterator it = begin(matrix);
        for (TPosition i = 0, iend = length(matrix, 0); i < iend; ++i) {
            *it = 0;
            goNext(it, 0);
        }
    } else {
        TScoreValue gapScore = scoreGap(scoringScheme);
        TScoreValue x = 0;
        TIterator it = begin(matrix);
        for (TPosition i = 0, iend = length(matrix, 0); i < iend; ++i) {
            *it = x;
            x += gapScore;
            goNext(it, 0);
        }
    }
    SEQAN_ASSERT_FAIL("Not implemented!");
}


template <typename TScoreValue, typename TSequence>
inline void
_align_fillMatrix(Matrix<TScoreValue, 2> & matrix, TSequence const & sequence0, TSequence const & sequence1, Score<TScoreValue, Simple> const & scoringScheme, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_FAIL("Not implemented!");
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_DYNPROG_AFFINE_H_

