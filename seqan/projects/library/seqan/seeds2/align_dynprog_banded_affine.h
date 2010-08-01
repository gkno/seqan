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
 ==========================================================================*/

#ifndef SEQAN_SEEDS_ALIGN_DYNPROG_BANDED_AFFINE_H_
#define SEQAN_SEEDS_ALIGN_DYNPROG_BANDED_AFFINE_H_

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

template <typename TScoreValue, typename TSequence, typename TDiagonal>
inline void
_alignBanded_resizeMatrix(Matrix<TScoreValue, 2> & matrix, TSequence const & sequence0, TSequence const & sequence1, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    (void) matrix;
    (void) sequence0;
    (void) sequence1;
    (void) lowerDiagonal;
    (void) upperDiagonal;

    SEQAN_ASSERT_FAIL("Write me!");
}


template <typename TScoreValue, typename TDiagonal, bool BEGIN1_FREE, bool BEGIN0_FREE, bool END1_FREE, bool END0_FREE>
inline void
_alignBanded_initGutter(Matrix<TScoreValue, 2> & matrix, Score<TScoreValue, Simple> const & scoringScheme, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, AlignConfig<BEGIN1_FREE, BEGIN0_FREE, END1_FREE, END0_FREE> const &, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    (void) matrix;
    (void) scoringScheme;
    (void) lowerDiagonal;
    (void) upperDiagonal;

    SEQAN_ASSERT_FAIL("Write me!");
}


template <typename TScoreValue, typename TDiagonal, typename TOverlap>
inline void
_alignBanded_initGutterFromUnbanded(Matrix<TScoreValue, 2> & matrix, Score<TScoreValue, Simple> const & scoringScheme, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, Matrix<TScoreValue, 2> /*const*/ & otherMatrix, TOverlap overlap0, TOverlap overlap1, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    (void) matrix;
    (void) scoringScheme;
    (void) lowerDiagonal;
    (void) upperDiagonal;
    (void) otherMatrix;
    (void) overlap0;
    (void) overlap1;

    SEQAN_ASSERT_FAIL("Write me!");
}


template <typename TScoreValue, typename TSequence, typename TDiagonal>
inline void
_alignBanded_fillMatrix(Matrix<TScoreValue, 2> & matrix, TSequence const & sequence0, TSequence const & sequence1, Score<TScoreValue, Simple> const & scoringScheme, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    (void) matrix;
    (void) sequence0;
    (void) sequence1;
    (void) scoringScheme;
    (void) lowerDiagonal;
    (void) upperDiagonal;

    SEQAN_ASSERT_FAIL("Write me!");
}


// Compute traceback in the given banded alignment matrix, starting at
// the lower right, shifted by the given overlap to the upper left.
// The matrix was filled with the Needleman-Wunschla lgorithm, end
// gaps are free as configured by the AlignConfig object.  Returns the
// best score.
template <typename TAlignmentIterator, typename TSequenceIterator, typename TPosition, typename TScoreValue, typename TScoringScheme, typename TOverlap, bool START0_FREE, bool START1_FREE, bool END0_FREE, bool END1_FREE>
TScoreValue
_alignBanded_traceBack(TAlignmentIterator & alignmentIt0, TAlignmentIterator & alignmentIt1, TSequenceIterator & sourceIt0, TSequenceIterator & sourceIt1, TPosition & finalPos0, TPosition & finalPos1, Matrix<TScoreValue, 2> /*const*/ & matrix, TScoringScheme const & scoringScheme, TOverlap overlap0, TOverlap overlap1, bool goToTopLeft, AlignConfig<START1_FREE, START0_FREE, END1_FREE, END0_FREE> const &, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    (void) alignmentIt0;
    (void) alignmentIt1;
    (void) sourceIt0;
    (void) sourceIt1;
    (void) matrix;
    (void) scoringScheme;
    (void) overlap0;
    (void) overlap1;
    (void) goToTopLeft;
    (void) finalPos0;
    (void) finalPos1;

    SEQAN_ASSERT_FAIL("Write me!");
    return 0;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_DYNPROG_BANDED_AFFINE_H_

