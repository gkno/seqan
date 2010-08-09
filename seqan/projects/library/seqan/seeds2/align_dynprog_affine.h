/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  ============================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ============================================================================
  Classic Gotoh DP algorithm for affine gap costs.  The algorithm is
  split into matrix resizing, filling the gutter, filling the rest of
  the matrix and doing the traceback.  Otherwise, it is a text book
  implementation without any tricks.

  We implement the three matrices as one 3-dimensional matrix.  The
  first two dimensions are for positions in sequences 0 and 1.  The
  third dimension differentiates between the three matrices.  They are
  -- in order -- M, I^a, I^b.
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
_align_resizeMatrix(Matrix<TScoreValue, 3> & matrix, TSequence const & sequence0, TSequence const & sequence1, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    setLength(matrix, 0, length(sequence0) + 1);
    setLength(matrix, 1, length(sequence1) + 1);
    setLength(matrix, 2, 3);
    resize(matrix);
}


template <typename TScoreValue, bool BEGIN1_FREE, bool BEGIN0_FREE, bool END1_FREE, bool END0_FREE>
inline void
_align_initGutter(Matrix<TScoreValue, 3> & matrix, Score<TScoreValue, Simple> const scoringScheme, AlignConfig<BEGIN1_FREE, BEGIN0_FREE, END1_FREE, END0_FREE> const &, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_EQ(length(matrix, 2), 3u);

    typedef Matrix<TScoreValue, 3> TMatrix;
    typedef typename Iterator<TMatrix>::Type TIterator;
    typedef typename Position<TMatrix>::Type TPosition;

    // TODO(holtgrew): Support free begin gaps.
    SEQAN_ASSERT_NOT_MSG(BEGIN0_FREE, "Free begin gaps are not supported in seeds align module yet.");
    SEQAN_ASSERT_NOT_MSG(BEGIN1_FREE, "Free begin gaps are not supported in seeds align module yet.");

    // We do not take the real minimum here because of overflows.
    TScoreValue inf = InfimumValue<TScoreValue>::VALUE / 2;
    // Get shortcuts to gap extension score.
    TScoreValue gapExtendScore = scoreGapExtend(scoringScheme);

    // Get iterators to the top left corners of the matrices.
    TIterator itMBegin = begin(matrix);
    TIterator itIABegin = itMBegin;
    goNext(itIABegin, 2);
    TIterator itIBBegin = itIABegin;
    goNext(itIBBegin, 2);

    // Init top and left gutters.
    {
        // Left gutters...
        TIterator itM = itMBegin;
        *itM = 0;
        goNext(itM, 0);
        TIterator itIB = itIBBegin;
        goNext(itIB, 0);
        for (TPosition pos = 1, posEnd = length(matrix, 0); pos < posEnd; ++pos) {
            *itM = pos * gapExtendScore;
            *itIB = inf;
            goNext(itM, 0);
            goNext(itIB, 0);
        }
        // Top gutters...
        itM = itMBegin;
        goNext(itM, 1);
        TIterator itIA = itIABegin;
        goNext(itIA, 1);
        for (TPosition pos = 1, posEnd = length(matrix, 1); pos < posEnd; ++pos) {
            *itM = pos * gapExtendScore;
            *itIA = inf;
            goNext(itM, 1);
            goNext(itIA, 1);
        }
    }
}


template <typename TScoreValue, typename TDiagonal, typename TOverlap>
inline void
_align_initGutterFromBanded(Matrix<TScoreValue, 3> & matrix, Score<TScoreValue, Simple> const & scoringScheme, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, Matrix<TScoreValue, 3> /*const*/ & otherMatrix, TOverlap overlap0, TOverlap overlap1, Gotoh const &)
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
_align_fillMatrix(Matrix<TScoreValue, 3> & matrix, TSequence const & sequence0, TSequence const & sequence1, Score<TScoreValue, Simple> const & scoringScheme, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_EQ(length(matrix, 2), 3u);

    typedef Matrix<TScoreValue, 3> TMatrix;
    typedef typename Position<TMatrix>::Type TPosition;
    typedef typename Iterator<TMatrix>::Type TMatrixIterator;
    typedef typename Iterator<TSequence, Standard>::Type TSequenceIterator;

    // We need three iterators in each alignment matrix to fill it.
    // it*Top point to the cell in the top row of the current column.
    // it*Left points to the column to the top left of the current
    // cell.  it*Above points to the cell above the current cell.  We
    // can use it*Above for value assignment.
    TMatrixIterator itMTop = begin(matrix);
    TMatrixIterator itMAbove;
    TMatrixIterator itMLeft;
    TMatrixIterator itIATop = itMTop;
    goNext(itIATop, 2);
    TMatrixIterator itIAAbove;
    TMatrixIterator itIALeft;
    TMatrixIterator itIBTop = itIATop;
    goNext(itIBTop, 2);
    TMatrixIterator itIBAbove;
    TMatrixIterator itIBLeft;

    // Compute score values so they are available locally and it maybe
    // is easier to the compiler to save some indirect adressing.
    TScoreValue matchScore = scoreMatch(scoringScheme);
    TScoreValue mismatchScore = scoreMismatch(scoringScheme);
    TScoreValue gapOpenScore = scoreGapOpen(scoringScheme);
    TScoreValue gapExtendScore = scoreGapExtend(scoringScheme);

    // Perform the matrix filling, column-wise for cache efficiency.
    for (TSequenceIterator it1 = begin(sequence1), it1End = end(sequence1); it1 != it1End; ++it1) {
        // std::cout << "iteration it1" << std::endl;
        itMLeft = itMTop;
        itIALeft = itIATop;
        itIBLeft = itIBTop;
        goNext(itMTop, 1);
        goNext(itIATop, 1);
        goNext(itIBTop, 1);
        itMAbove = itMTop;
        itIAAbove = itIATop;
        itIBAbove = itIBTop;
        for (TSequenceIterator it0 = begin(sequence0), it0End = end(sequence0); it0 != it0End; ++it0) {
            // std::cout << "  iteration it0" << std::endl;
            // Compute M_{i-1,j-1} + match/mismatch
            TScoreValue scoreMoveDiagonal = *itMLeft + ((*it0 == *it1) ? matchScore : mismatchScore);
            goNext(itMLeft, 0);
            // Compute I^a_{i,j}
            TScoreValue scoreIA = _max(*itMAbove + gapOpenScore, *itIAAbove + gapExtendScore);
            goNext(itIALeft, 0);  // TODO(holtgrew): Remove IALeft? Not necessary!
            goNext(itIAAbove, 0);
            *itIAAbove = scoreIA;
            // Compute I^b_{i,j}
            goNext(itIBLeft, 0);
            TScoreValue scoreIB = _max(*itMLeft + gapOpenScore, *itIBLeft + gapExtendScore);
            goNext(itIBAbove, 0);
            *itIBAbove = scoreIB;
            // Assign M_{i,j}
            goNext(itMAbove, 0);
            *itMAbove = _max(scoreMoveDiagonal, _max(*itIAAbove, *itIBAbove));
        }
    }
}

// Compute traceback in the given normal DP alignment matrix, starting
// at the lower right, shifted by the given overlap to the upper left.
// The matrix was filled with the Needleman-Wunschla lgorithm, end
// gaps are free as configured by the AlignConfig object.  Returns the
// best score.
template <typename TAlignmentIterator, typename TSequenceIterator, typename TPosition, typename TScoreValue, typename TScoringScheme, typename TOverlap, bool START0_FREE, bool START1_FREE, bool END0_FREE, bool END1_FREE>
TScoreValue
_align_traceBack(TAlignmentIterator & alignmentIt0, TAlignmentIterator & alignmentIt1, TSequenceIterator & sourceIt0, TSequenceIterator & sourceIt1, TPosition & finalPos0, TPosition & finalPos1, Matrix<TScoreValue, 3> /*const*/ & matrix, TScoringScheme const & scoringScheme, TOverlap overlap0, TOverlap overlap1, bool goToTopLeft, AlignConfig<START1_FREE, START0_FREE, END1_FREE, END0_FREE> const &, Gotoh const &)
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
    
    SEQAN_ASSERT_FAIL("Not implemented!");

    return 0;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_DYNPROG_AFFINE_H_

