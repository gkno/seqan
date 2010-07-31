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

    (void)matrix;
    (void)sequence0;
    (void)sequence1;

    SEQAN_ASSERT_FAIL("Not implemented!");
}


template <typename TScoreValue, bool BEGIN1_FREE, bool BEGIN0_FREE, bool END1_FREE, bool END0_FREE>
inline void
_align_initGutter(Matrix<TScoreValue, 2> & matrix, Score<TScoreValue, Simple> const scoringScheme, AlignConfig<BEGIN1_FREE, BEGIN0_FREE, END1_FREE, END0_FREE> const &, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    (void)matrix;
    (void)scoringScheme;

    SEQAN_ASSERT_FAIL("Not implemented!");
}


template <typename TScoreValue, typename TDiagonal, typename TOverlap>
inline void
_align_initGutterFromBanded(Matrix<TScoreValue, 2> & matrix, Score<TScoreValue, Simple> const & scoringScheme, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, Matrix<TScoreValue, 2> /*const*/ & otherMatrix, TOverlap overlap0, TOverlap overlap1, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    (void)matrix;
    (void)scoringScheme;
    (void)lowerDiagonal;
    (void)upperDiagonal;
    (void)otherMatrix;
    (void)overlap0;
    (void)overlap1;

    SEQAN_ASSERT_FAIL("Not implemented!");
}


template <typename TScoreValue, typename TSequence>
inline void
_align_fillMatrix(Matrix<TScoreValue, 2> & matrix, TSequence const & sequence0, TSequence const & sequence1, Score<TScoreValue, Simple> const & scoringScheme, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    (void)matrix;
    (void)sequence0;
    (void)sequence1;
    (void)scoringScheme;

    SEQAN_ASSERT_FAIL("Not implemented!");
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_DYNPROG_AFFINE_H_

