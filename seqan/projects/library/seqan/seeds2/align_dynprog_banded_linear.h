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
 ============================================================================
  Banded linear programming for sequence alignment.  Needleman Wunsch
  implementation, i.e. for linear gap costs only.

  The alignment code could use some optimization.  However, we cannot
  use the same optimization as in the graph alignment since we want to
  compute globally optimal trace through multiple connected alignment
  matrices.

  The banding is done in dimension 1, dimension 0 is there completely.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_ALIGN_DYNPROG_BANDED_LINEAR_H_
#define SEQAN_SEEDS_ALIGN_DYNPROG_BANDED_LINEAR_H_

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

template <typename TScoreValue, typename TSequence, typename TBandwidth>
inline void
_alignBanded_resizeMatrix(Matrix<TScoreValue, 2> & matrix, TSequence const & sequence0, TSequence const & sequence1, TBandwidth bandwidth, NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;

    // We need space length(sequence0) x bandwidth for the alignment,
    // the top gutter and the left and right gutters.
    setLength(matrix, 0, length(sequence0) + 1);
    setLength(matrix, 1, bandwidth + 2);
    resize(matrix);
}


template <typename TScoreValue, bool BEGIN1_FREE, bool BEGIN0_FREE, bool END1_FREE, bool END0_FREE, typename TBandwidth>
inline void
_alignBanded_initGutter(Matrix<TScoreValue, 2> & matrix, Score<TScoreValue, Simple> const scoringScheme, AlignConfig<BEGIN1_FREE, BEGIN0_FREE, END1_FREE, END0_FREE> const &, TBandwidth bandwidth, NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_EQ_MSG(scoreGapOpen(scoringScheme), scoreGapExtend(scoringScheme),
                        "Only linear gap costs allowed for Needleman-Wunsch.");

    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Iterator<TMatrix>::Type TIterator;
    typedef typename Position<TMatrix>::Type TPosition;
    
    // Initialize the left gutter with infima.
    TIterator it = begin(matrix);
    goNext(it, 0);
    for (TPosition i = 0, iend = length(matrix, 0); i < iend; ++i, goNext(it, 0))
        *it = InfimumValue<TScoreValue>::VALUE / 2;
    // Initialize the top gutter according to the AlignConfig.
    it = begin(matrix);
    if (BEGIN0_FREE) {
        for (TPosition i = 0, iend = length(matrix, 1); i < iend; ++i, goNext(it, 1))
            *it = 0;
    } else {
        TScoreValue gapScore = scoreGap(scoringScheme);
        TScoreValue x = 0;
        for (TPosition i = 0, iend = length(matrix, 1); i < iend; ++i, goNext(it, 1)) {
            *it = x;
            x += gapScore;
        }
    }
    goPrevious(it, 1);
    goNext(it, 0);
    // Initialize the right gutter with infima.
    for (TPosition i = 0, iend = length(matrix, 0) - 1; i < iend; ++i, goNext(it, 0))
        *it = InfimumValue<TScoreValue>::VALUE / 2;
}


template <typename TScoreValue, typename TSequence, typename TBandwidth>
inline void
_alignBanded_fillMatrix(Matrix<TScoreValue, 2> & matrix, TSequence const & sequence0, TSequence const & sequence1, Score<TScoreValue, Simple> const & scoringScheme, TBandwidth bandwidth, NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_EQ_MSG(scoreGapOpen(scoringScheme), scoreGapExtend(scoringScheme),
                        "Only linear gap costs allowed for Needleman-Wunsch.");

    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Iterator<TMatrix>::Type TMatrixIterator;
    typedef typename Iterator<TSequence, Standard>::Type TSequenceIterator;

    // We use three iterators for filling the matrix.  itLeftmost
    // always points to the leftmost cell of the current row, itLeft
    // points to the field left of the current cell, itTop points to
    // the cell to the upper left of the current cell.
    TMatrixIterator itLeftmost = begin(matrix);
    TMatrixIterator itLeft;
    TMatrixIterator itTop;

    // Compute score values so they are available locally and it maybe
    // is easier to the compiler to save some indirect adressing.
    TScoreValue matchScore = scoreMatch(scoringScheme);
    TScoreValue mismatchScore = scoreMismatch(scoringScheme);
    TScoreValue gapScore = scoreGap(scoringScheme);
    
    // TODO(holtgrew): Because we band in dimension 0, we have to fill the alignment matrix row-wise.  This, however is not cache-efficient.  Maybe tranpose matrix here?

    for (TSequenceIterator it0 = begin(sequence0), it0End = end(sequence0); it0 != it0End; ++it0) {
        itTop = itLeftmost;
        goNext(itLeftmost, 0);
        itLeft = itLeftmost;
        for (TSequenceIterator it1 = begin(sequence0), it1End = end(sequence0); it1 != it1End; ++it1) {
            TScoreValue scoreDiagonal = *itTop + ((*it0 == *it1) ? matchScore : mismatchScore);
            goNext(itTop, 1);
            TScoreValue scoreVertical = *itTop + gapScore;
            TScoreValue scoreHorizontal = *itLeft + gapScore;
            goNext(itLeft, 0);
            *itLeft = _max(scoreDiagonal, _max(scoreVertical, scoreHorizontal));
        }
    }
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_DYNPROG_BANDED_LINEAR_H_

