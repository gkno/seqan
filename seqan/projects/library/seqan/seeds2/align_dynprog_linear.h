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
  Classic linear programming for sequence alignment.  Needleman Wunsch
  implementation, i.e. for linear gap costs only.

  The alignment code could use some optimization.  However, we cannot
  use the same optimization as in the graph alignment since we want to
  compute globally optimal trace through multiple connected alignment
  matrices.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_ALIGN_DYNPROG_LINEAR_H_
#define SEQAN_SEEDS_ALIGN_DYNPROG_LINEAR_H_

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
_align_resizeMatrix(Matrix<TScoreValue, 2> & matrix, TSequence const & sequence0, TSequence const & sequence1, NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;

    setLength(matrix, 0, length(sequence0) + 1);
    setLength(matrix, 1, length(sequence1) + 1);
    resize(matrix);
}


template <typename TScoreValue, bool BEGIN1_FREE, bool BEGIN0_FREE, bool END1_FREE, bool END0_FREE>
inline void
_align_initGutter(Matrix<TScoreValue, 2> & matrix, Score<TScoreValue, Simple> const scoringScheme, AlignConfig<BEGIN1_FREE, BEGIN0_FREE, END1_FREE, END0_FREE> const &, NeedlemanWunsch const &)
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

    // Init top gutter with zeroes if begin gaps are in dimension 1
    // free or with gap scores otherwise.
    if (BEGIN1_FREE) {
        TIterator it = begin(matrix);
        for (TPosition i = 0, iend = length(matrix, 1); i < iend; ++i) {
            *it = 0;
            goNext(it, 1);
        }
    } else {
        TScoreValue gapScore = scoreGap(scoringScheme);
        TScoreValue x = 0;
        TIterator it = begin(matrix);
        for (TPosition i = 0, iend = length(matrix, 1); i < iend; ++i) {
            *it = x;
            x += gapScore;
            goNext(it, 1);
        }
    }
}


template <typename TScoreValue, typename TSequence>
inline void
_align_fillMatrix(Matrix<TScoreValue, 2> & matrix, TSequence const & sequence0, TSequence const & sequence1, Score<TScoreValue, Simple> const & scoringScheme, NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_EQ_MSG(scoreGapOpen(scoringScheme), scoreGapExtend(scoringScheme),
                        "Only linear gap costs allowed for Needleman-Wunsch.");

    // Fill matrix with standard NW DP programming.

    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Iterator<TMatrix>::Type TMatrixIterator;
    typedef typename Iterator<TSequence, Standard>::Type TSequenceIterator;

    // We need three iterators in the alignment matrix to fill it.
    // itTop points to the cell in the top row of the current column.
    // itLeft points to the column to the top left of the current
    // cell.  itAbove points to the cell above the current cell.  We
    // can use itAbove for value assignment.
    TMatrixIterator itTop = begin(matrix);
    TMatrixIterator itAbove;
    TMatrixIterator itLeft;

    // Compute score values so they are available locally and it maybe
    // is easier to the compiler to save some indirect adressing.
    TScoreValue matchScore = scoreMatch(scoringScheme);
    TScoreValue mismatchScore = scoreMismatch(scoringScheme);
    TScoreValue gapScore = scoreGap(scoringScheme);
    
    // Perform the Needleman-Wunsch dynamic programming.
    for (TSequenceIterator it1 = begin(sequence1), it1end = end(sequence1); it1 != it1end; ++it1) {
        itLeft = itTop;
        goNext(itTop, 1);
        itAbove = itTop;
        for (TSequenceIterator it0 = begin(sequence0), it0end = end(sequence0); it0 != it0end; ++it0) {
            TScoreValue scoreMoveDiagonal = *itLeft + ((*it1 == *it0) ? matchScore : mismatchScore);
            goNext(itLeft, 0);
            TScoreValue scoreMoveRight = *itLeft + gapScore;
            TScoreValue scoreMoveDown = *itAbove + gapScore;
            goNext(itAbove, 0);
            *itAbove = _max(scoreMoveDiagonal, _max(scoreMoveRight, scoreMoveDown));
        }
    }
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_DYNPROG_LINEAR_H_

