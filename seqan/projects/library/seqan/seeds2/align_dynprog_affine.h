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

// Fill the Gotoh DP matrices for the alignment of the two given
// sequences.  The alignment configuration gives whether the first
// column and row are initialized with value 0 or with the gap scores.
//
// This is the most basic case.
template <typename TScoreValue, typename TSequence, typename TScoringScheme, typename TAlignmentConfig>
void
_align_dynProg(
        Matrix<TScoreValue, 2> & matrix,
        TSequence const & seq0,
        TSequence const & seq1,
        TScoringScheme const & scoringScheme,
        TAlignmentConfig const &,
        Gotoh const &)
{
    SEQAN_CHECKPOINT;
}


// Fill the NW DP matrix for the alignment of the two given sequences.
//
// This is the case where data is copied over from a banded seed
// alignment DP matrix.
//
// The data of the lower right rectangle of the DP seed alignment
// matrix is copied over into the upper right corner of the alignment
// matrix and the alignment border is filled with infimum values.
template <typename TScoreValue, typename TSequence, typename TScoringScheme>
void
_align_dynProg(
        Matrix<TScoreValue, 2> & matrix,
        TSequence const & seq0,
        TSequence const & seq1,
        TScoringScheme const & scoringScheme,
        _DPMatrixRectangle<TScoreValue> /*const*/ & lowerRightRectangle,
        Gotoh const &)
{
    SEQAN_CHECKPOINT;
}


// Perform a traceback in the Gotoh alignment matrices beginning with
// the characters at pos0/pos1 of the segments of this matrix.
template <typename TTarget, typename TScoreValue, typename TPosition>
void
_align_traceBack(
        TTarget gaps0It,
        TTarget gaps1It,
        Matrix<TScoreValue, 2> const & matrix,
        TPosition pos0,
        TPosition pos1,
        Gotoh const &)
{
    SEQAN_CHECKPOINT;
}   

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_DYNPROG_AFFINE_H_

