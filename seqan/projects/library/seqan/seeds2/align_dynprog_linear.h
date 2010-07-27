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

template <typename TScoreValue, typename TScoringScheme>
inline void
_initDPMatrixEdgeDim(
        Matrix<TScoreValue, 2> & matrix,
        TScoringScheme const & scoringScheme,
        unsigned dim,
        bool fillWithZeros,
        NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;
    
    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Iterator<TMatrix>::Type TMatrixIterator;
    typedef typename Position<TMatrix>::Type TPosition;

    TMatrixIterator it = begin(matrix);
    if (fillWithZeros) {
        for (TPosition i = 0; i < length(matrix, dim); ++i) {
            *it = 0;
            goNext(it, dim);
        }
    } else {
        TScoreValue x = 0;
        for (TPosition i = 0; i < length(matrix, dim); ++i) {
            *it = x;
            x += scoreGap(scoringScheme);
            goNext(it, dim);
        }
    }
}

// Fill the Needleman-Wunsch DP matrix for the alignment of the two
// given sequences.  The alignment configuration gives whether the
// first column and row are initialized with value 0 or with the gap
// scores.
//
// This is the most basic case.
template <typename TScoreValue, typename TSequence, typename TScoringScheme, bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TAlignConfigSpec>
// TODO(holtgrew): Which bool for AlignConfig is for which sequence? We assume LEFT/RIGHT is dim0, TOP/BOTTOM is dim1.
void
_align_dynProg(
        Matrix<TScoreValue, 2> & matrix,
        TSequence const & sequence0,
        TSequence const & sequence1,
        TScoringScheme const & scoringScheme,  // TODO(holtgrew): Enforce Simple Score?
        AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TAlignConfigSpec> const &,
        NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;

    // In this function, we assume sequence 0 is written along
    // dimension 0 of the matrix and this is the vertical direction.
    // sequence 1 corresponds to dimension 0 and the horizontal
    // direction.

    SEQAN_ASSERT_EQ_MSG(scoreGapOpen(scoringScheme), scoreGapExtend(scoringScheme),
                        "Needleman-Wunsch DP only support linear gap costs.");

    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Iterator<TMatrix>::Type TMatrixIterator;
    typedef typename Position<TSequence>::Type TPosition;
    typedef typename Iterator<TSequence, Standard>::Type TSequenceIterator;

    // Resize the matrix, first row and column are gap scores.
    setLength(matrix, 0, length(sequence0) + 1);
    setLength(matrix, 1, length(sequence0) + 1);
    resize(matrix);

    // Initialize the left and top edges of the matrix.
    _initDPMatrixEdgeDim(matrix, scoringScheme, 0, LEFT, NeedlemanWunsch());
    _initDPMatrixEdgeDim(matrix, scoringScheme, 1, TOP, NeedlemanWunsch());

    // We need thre iterators in the alignment matrix to fill it.
    // itTop points to the cell in the top row of the current column.
    // itLeft points to the column to the top left of the current
    // cell.  itAbove points to the cell above the current cell.  We
    // can use itAbove for value assignment.
    TMatrixIterator itTop = begin(matrix);
    TMatrixIterator itLeft;
    TMatrixIterator itAbove;

    // Perform the Needleman-Wunsch dynamic programming.
    for (TSequenceIterator it1 = begin(sequence1); it1 != end(sequence1); ++it1) {
        itLeft = itTop;
        goNext(itTop, 1);
        itAbove = itTop;
        for (TSequenceIterator it0 = begin(sequence0); it0 != end(sequence0); ++it0) {
            TScoreValue scoreMoveDiagonal = *itLeft + ((*it1 == *it0) ? scoreMatch(scoringScheme) : scoreMismatch(scoringScheme));
            goNext(itLeft, 0);
            TScoreValue scoreMoveRight = *itLeft + scoreGap(scoringScheme);
            TScoreValue scoreMoveDown = *itAbove + scoreGap(scoringScheme);
            goNext(itAbove, 0);
            *itAbove = _max(scoreMoveDiagonal, _max(scoreMoveRight, scoreMoveDown));
        }
    }
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

    SEQAN_ASSERT_EQ_MSG(scoreGapOpen(scoringScheme), scoreGapExtend(scoringScheme),
                        "Needleman-Wunsch DP only support linear gap costs.");

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

