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

#ifndef SEQAN_SEEDS_ALIGN_DYNPROG_LINEAR_H_
#define SEQAN_SEEDS_ALIGN_DYNPROG_LINEAR_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// Describes a rectangle from a DP seed alignment matrix.  This is
// used in the chain alignment to describe DP matrix parts to copy
// into other alignment matrices.
template <typename TScoreValue>
struct _DPMatrixRectangle
{
    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Size<TMatrix>::Type TSize;
    typedef typename Position<TMatrix>::Type TPosition;
    
    TMatrix & matrix;
    TSize length0;
    TSize length1;

    _DPMatrixRectangle(TMatrix & matrix_, TSize length0_, TSize length1_)
            : matrix(matrix_), length0(length0_), length1(length1_)
    { SEQAN_CHECKPOINT; }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

// Fill the Needleman-Wunsch DP matrix for the alignment of the two
// given sequences.  The alignment configuration gives whether the
// first column and row are initialized with value 0 or with the gap
// scores.
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
        NeedlemanWunsch const &)
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
        NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;
}


// Perform a traceback in the Needleman-Wunsch alignment matrix
// beginning with the characters at pos0/pos1 of the segments of this
// matrix.
template <typename TTarget, typename TScoreValue, typename TPosition>
void
_align_traceBack(
        TTarget gaps0It,
        TTarget gaps1It,
        Matrix<TScoreValue, 2> const & matrix,
        TPosition pos0,
        TPosition pos1,
        NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;
}   

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_DYNPROG_LINEAR_H_

