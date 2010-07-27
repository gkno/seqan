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

// Fill the NW DP matrix for the banded alignment of the two given
// sequences.
//
// This is the case where data is copied over from an alignment DP
// matrix.
//
// The data of the lower right rectangle of the DP alignment matrix is
// copied over into the upper right corner of the alignment matrix and
// the alignment border is filled with infimum values.
template <typename TScoreValue, typename TSequence, typename TDiagonal, typename TScoringScheme>
void
_align_banded_dynProg(
        Matrix<TScoreValue, 2> & matrix,
        TSequence const & seq0,
        TSequence const & seq1,
        TDiagonal lowerDiagonal,
        TDiagonal upperDiagonal,
        TScoringScheme const & scoringScheme,
        _DPMatrixRectangle<TScoreValue> /*const*/ & lowerRightRectangle,
        NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_GT(lowerRightRectangle.length0, 0u);
    SEQAN_ASSERT_GT(lowerRightRectangle.length1, 0u);
    
    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Iterator<TMatrix>::Type TMatrixIterator;
    typedef typename Size<TSequence>::Type TSize;
    typedef typename Position<TSequence>::Type TPosition;

    // Compute bandwidth in direction of dimension 1.  Band is
    // "jolted" to the right to give a rectangle for the matrix.  This
    // matrix has size length(seq0) x bandwidth.
    TSize bandwidth = upperDiagonal - lowerDiagonal;
    // The width of the overlapping rectangle is in relation to the bandwidth.
    SEQAN_ASSERT_EQ(lowerRightRectangle.length1 + 1, bandwidth);
    // Length of the empty triangles.
    TSize triangleSideLength = lowerRightRectangle.length0;

    // Allocate memory for the alignment matrix.
    setLength(matrix, 0, length(seq0) + 1);
    setLength(matrix, 1, bandwidth + 2);
    resize(matrix);

    // Initialize matrix.
    //
    // First, fill the diagonal left of the lower one with infimas (we
    // actually need to pick a little higher value so no overflow
    // happens).  The diagonal corresponds to the first column.
    {
        TMatrixIterator it = begin(matrix);
        for (TPosition pos = lowerRightRectangle.length0; pos < length(seq0) + 1; ++pos) {
            *it = InfimumValue<TScoreValue>::VALUE / 2;
            goNext(it, 0);
        }
    }
    // Second, copy in the data from the rectangle.
    //
    // TODO(holtgrew): Actually, we could copy in less, only L-shaped border is required. Or this part could be ignored in the alignment in the upper left rectangle...
    TMatrixIterator destLeft = begin(matrix);
    TMatrixIterator srcLeft = begin(lowerRightRectangle.matrix);
    setPosition(srcLeft, length(lowerRightRectangle.matrix, 0)  - lowerRightRectangle.length0 + (length(lowerRightRectangle.matrix, 1) - lowerRightRectangle.length1) * _dataFactors(lowerRightRectangle.matrix)[1]);  // TODO(holtgrew): Extend matrix interface for something like this...
    unsigned skipped = 1;  // Number of cells skipped because of band/rectangle difference.
    for (TPosition offset0 = 0; offset0 < lowerRightRectangle.length0; ++offset0) {
        TMatrixIterator destIt = destLeft;
        TMatrixIterator srcIt = srcLeft;
        for (TPosition offset1 = skipped; offset1 < lowerRightRectangle.length1; ++offset1) {
            *destIt = *srcIt;
            goNext(destIt, 1);
            goNext(srcIt, 1);
        }
        goNext(destLeft, 0);
        goNext(srcLeft, 0);
        skipped += 1;
    }
}


// Perform a traceback in the Needleman-Wunsch banded alignment matrix
// beginning with the positions pos0/pos1 in the sequences.
template <typename TTarget, typename TScoreValue, typename TPosition>
void
_align_banded_traceBack(
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

#endif  // #ifndef SEQAN_SEEDS_ALIGN_DYNPROG_BANDED_LINEAR_H_

