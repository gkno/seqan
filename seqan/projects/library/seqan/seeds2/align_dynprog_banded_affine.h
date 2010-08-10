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

// TODO(holtgrew): Maybe adjust rausch's code so things can be copied in and use it instead since it is heavily tuned?

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
_alignBanded_resizeMatrix(Matrix<TScoreValue, 3> & matrix, TSequence const & sequence0, TSequence const & sequence1, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    (void) sequence1;  // In case we do not run in debug mode.

    SEQAN_ASSERT_GEQ(upperDiagonal, 0);
    SEQAN_ASSERT_LEQ(lowerDiagonal, 0);
    SEQAN_ASSERT_GEQ_MSG(length(sequence1) - lowerDiagonal, length(sequence0), "Lower diagonal is not low enough.");
    SEQAN_ASSERT_GEQ_MSG(length(sequence0) + upperDiagonal, length(sequence1), "Upper diagonal is not high enough.");

    // We need space length(sequence0) x bandwidth for the alignment,
    // the top gutter and the left and right gutters.
    setLength(matrix, 0, length(sequence0) + 1);
    setLength(matrix, 1, upperDiagonal - lowerDiagonal + 3);
    setLength(matrix, 2, 3);
    resize(matrix);
    //fill(matrix, -42);
}


template <typename TScoreValue, typename TDiagonal, bool BEGIN1_FREE, bool BEGIN0_FREE, bool END1_FREE, bool END0_FREE>
inline void
_alignBanded_initGutter(Matrix<TScoreValue, 3> & matrix, Score<TScoreValue, Simple> const & scoringScheme, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, AlignConfig<BEGIN1_FREE, BEGIN0_FREE, END1_FREE, END0_FREE> const &, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    typedef Matrix<TScoreValue, 3> TMatrix;
    typedef typename Iterator<TMatrix>::Type TIterator;
    typedef typename Position<TMatrix>::Type TPosition;
    
    // Initialize the diagonal below the lower one with infimas.
    TIterator diagonalIt = begin(matrix);
    TIterator verticalIt = diagonalIt;
    goNext(verticalIt, 2);
    TIterator horizontalIt = verticalIt;
    goNext(horizontalIt, 2);
    for (TPosition i = 0, iend = length(matrix, 0); i < iend; ++i, goNext(diagonalIt, 0), goNext(horizontalIt, 0)) {
        *diagonalIt = InfimumValue<TScoreValue>::VALUE / 2;
        *horizontalIt = InfimumValue<TScoreValue>::VALUE / 2;
    }
    // Initialize the left gutter according to the AlignConfig.
    setPosition(diagonalIt, (1 - lowerDiagonal) * _dataFactors(matrix)[1]); // TODO(holtgrew): There should be a function that accepts two coordinates for the Matrix class.
    verticalIt = diagonalIt;
    goNext(verticalIt, 2);
    horizontalIt = verticalIt;
    goNext(horizontalIt, 2);
    if (BEGIN0_FREE) {
        for (TPosition i = 0, iend = -lowerDiagonal; i <= iend; ++i, goNext(diagonalIt, 0), goPrevious(diagonalIt, 1)) {
            *diagonalIt = 0;
        }
    } else {
        TScoreValue gapScore = scoreGap(scoringScheme);
        TScoreValue x = 0;
        for (TPosition i = 0, iend = -lowerDiagonal; i <= iend; ++i, goNext(diagonalIt, 0), goPrevious(diagonalIt, 1)) {
            *diagonalIt = x;
            x += gapScore;
        }
    }
    for (TPosition i = 0, iend = -lowerDiagonal; i <= iend; ++i, goNext(verticalIt, 0), goPrevious(verticalIt, 1), goNext(horizontalIt, 0), goPrevious(horizontalIt, 1)) {
        *verticalIt = InfimumValue<TScoreValue>::VALUE / 2;
        *horizontalIt = InfimumValue<TScoreValue>::VALUE / 2;
    }
    // Initialize the top gutter according to the AlignConfig.
    setPosition(diagonalIt, (1 - lowerDiagonal) * _dataFactors(matrix)[1]); // TODO(holtgrew): There should be a function that accepts two coordinates for the Matrix class.
    verticalIt = diagonalIt;
    goNext(verticalIt, 2);
    horizontalIt = verticalIt;
    goNext(horizontalIt, 2);
    if (BEGIN1_FREE) {
        for (TPosition i = 0, iend = upperDiagonal; i <= iend; ++i, goNext(diagonalIt, 1))
            *diagonalIt = 0;
    } else {
        TScoreValue gapScore = scoreGap(scoringScheme);
        TScoreValue x = 0;
        for (TPosition i = 0, iend = upperDiagonal; i <= iend; ++i, goNext(diagonalIt, 1)) {
            *diagonalIt = x;
            x += gapScore;
        }
    }
    for (TPosition i = 0, iend = upperDiagonal; i <= iend; ++i, goNext(verticalIt, 1), goNext(horizontalIt, 1)) {
        *verticalIt = InfimumValue<TScoreValue>::VALUE / 2;
        *horizontalIt = InfimumValue<TScoreValue>::VALUE / 2;
    }
    // Initialize the diagonal above the upper one with infimas.
    for (TPosition i = 0, iend = length(matrix, 0); i < iend; ++i, goNext(diagonalIt, 0), goNext(verticalIt, 0)) {
        *diagonalIt = InfimumValue<TScoreValue>::VALUE / 2;
        *verticalIt = InfimumValue<TScoreValue>::VALUE / 2;
    }
}


template <typename TScoreValue, typename TDiagonal, typename TOverlap>
inline void
_alignBanded_initGutterFromUnbanded(Matrix<TScoreValue, 3> & matrix, Score<TScoreValue, Simple> const & scoringScheme, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, Matrix<TScoreValue, 3> /*const*/ & otherMatrix, TOverlap overlap0, TOverlap overlap1, Gotoh const &)
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
_alignBanded_fillMatrix(Matrix<TScoreValue, 3> & matrix, TSequence const & sequence0, TSequence const & sequence1, Score<TScoreValue, Simple> const & scoringScheme, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_GEQ_MSG(upperDiagonal, 0, "Upper diagonal must not lie below main diagonal.");
    SEQAN_ASSERT_LEQ_MSG(lowerDiagonal, 0, "Upper diagonal must not lie above main diagonal.");

    typedef Matrix<TScoreValue, 3> TMatrix;
    typedef typename Iterator<TMatrix>::Type TMatrixIterator;
    typedef typename Iterator<TSequence, Standard>::Type TSequenceIterator;
    typedef typename Size<TSequence>::Type TSize;
    typedef typename Position<TSequence>::Type TPosition;

    // We need three iterators in each alignment matrix to fill it.
    // it*Top points to the cell in the top row of the current column.
    // it*Left points to the column to the top left of the current
    // cell.  it*Above points to the cell above the current cell.  We
    // can use itAbove for value assignment.
    TMatrixIterator itMTop = begin(matrix);
    TMatrixIterator itMAbove;
    TMatrixIterator itMLeft;
    TMatrixIterator itIATop;
    TMatrixIterator itIAAbove;
    TMatrixIterator itIALeft;
    TMatrixIterator itIBTop;
    TMatrixIterator itIBAbove;
    TMatrixIterator itIBLeft;

    // Compute score values so they are available locally and it maybe
    // is easier to the compiler to save some indirect adressing.
    TScoreValue matchScore = scoreMatch(scoringScheme);
    TScoreValue mismatchScore = scoreMismatch(scoringScheme);
    TScoreValue gapOpenScore = scoreGapOpen(scoringScheme);
    TScoreValue gapExtendScore = scoreGapExtend(scoringScheme);

    // Fill matrix column wise for cache efficiency.  For this, we
    // need the length of the current column in the sheared alignment
    // matrix which complicates things.
    setPosition(itMTop, (1 - lowerDiagonal) * _dataFactors(matrix)[1]); // TODO(holtgrew): There should be a function that accepts two coordinates for the Matrix class.
    itIATop = itMTop;
    goNext(itIATop, 2);
    itIBTop = itIATop;
    goNext(itIBTop, 2);
    TPosition seq1Pos = 0;
    TSequenceIterator it0Begin = begin(sequence0);
    for (TSequenceIterator it1 = begin(sequence1), it1end = end(sequence1); it1 != it1end; ++it1, ++seq1Pos) {
        if (seq1Pos <= static_cast<TPosition>(upperDiagonal)) {
            itMLeft = itMTop;
            goNext(itMTop, 1);
            itIALeft = itIATop;
            goNext(itIATop, 1);
            itIBLeft = itIBTop;
            goNext(itIBTop, 1);
        } else {
            goNext(itMTop, 0);
            itMLeft = itMTop;
            goPrevious(itMLeft, 1);
            goNext(itIATop, 0);
            itIALeft = itIATop;
            goPrevious(itIALeft, 1);
            goNext(itIBTop, 0);
            itIBLeft = itIBTop;
            goPrevious(itIBLeft, 1);
        }
        itMAbove = itMTop;
        itIAAbove = itIATop;
        itIBAbove = itIBTop;
        TDiagonal from = _max(static_cast<TDiagonal>(seq1Pos) - upperDiagonal, TDiagonal(0));
        TDiagonal to = _min(seq1Pos - lowerDiagonal, length(sequence0) - 1);
        TDiagonal tmp = 1 + to - from;
        TSize runLength = _max(TDiagonal(0), tmp);
        // std::cout << "RUN LENGTH == " << runLength << std::endl;
        // std::cout << "from == " << from << std::endl << "to == " << to << std::endl;
        // std::cout << "seq1Pos = " << seq1Pos << std::endl;
        // std::cout << "lowerDiagonal = " << lowerDiagonal << std::endl;
        // std::cout << "length(sequence0) = " << length(sequence0) << std::endl;
        // std::cout << "seq1Pos = " << seq1Pos << std::endl;
        // std::cout << "upperDiagonal = " << upperDiagonal << std::endl;
        TSize i = 0;
        for (TSequenceIterator it0 = it0Begin; i < runLength; ++it0, ++i) {
            // Compute M_{i-1,j-1} + match/mismatch
            TScoreValue scoreMoveDiagonal = *itMLeft + ((*it0 == *it1) ? matchScore : mismatchScore);
            goNext(itMLeft, 0);
            goPrevious(itMLeft, 1);
            // Compute I^a_{i,j}
            TScoreValue scoreIA = _max(*itMAbove + gapOpenScore, *itIAAbove + gapExtendScore);
            goNext(itIALeft, 0);  // TODO(holtgrew): Remove IALeft? Not necessary!
            goPrevious(itIALeft, 1);
            goNext(itIAAbove, 0);
            goPrevious(itIAAbove, 1);
            *itIAAbove = scoreIA;
            // Compute I^b_{i,j}
            goNext(itIBLeft, 0);
            goPrevious(itIBLeft, 1);
            TScoreValue scoreIB = _max(*itMLeft + gapOpenScore, *itIBLeft + gapExtendScore);
            goNext(itIBAbove, 0);
            goPrevious(itIBAbove, 1);
            *itIBAbove = scoreIB;
            // Assign M_{i,j}
            goNext(itMAbove, 0);
            goPrevious(itMAbove, 1);
            *itMAbove = _max(scoreMoveDiagonal, _max(*itIAAbove, *itIBAbove));
        }
        // std::cout << "++it1" << std::endl;
        // We only need an offset for it0 when all diagonals up to the
        // upper one are aligned.
        if (seq1Pos >= static_cast<TPosition>(upperDiagonal))
            it0Begin += 1;
    }
}


// Compute traceback in the given banded alignment matrix, starting at
// the lower right, shifted by the given overlap to the upper left.
// The matrix was filled with the Needleman-Wunschla lgorithm, end
// gaps are free as configured by the AlignConfig object.  Returns the
// best score.
template <typename TAlignmentIterator, typename TSequenceIterator, typename TPosition, typename TScoreValue, typename TScoringScheme, typename TOverlap, bool START0_FREE, bool START1_FREE, bool END0_FREE, bool END1_FREE>
TScoreValue
_alignBanded_traceBack(TAlignmentIterator & alignmentIt0, TAlignmentIterator & alignmentIt1, TSequenceIterator & sourceIt0, TSequenceIterator & sourceIt1, TPosition & finalPos0, TPosition & finalPos1, Matrix<TScoreValue, 3> /*const*/ & matrix, TScoringScheme const & scoringScheme, TOverlap overlap0, TOverlap overlap1, TOverlap upperTriangleEdgeWidth, TOverlap lowerTriangleEdgeWidth, bool goToTopLeft, AlignConfig<START1_FREE, START0_FREE, END1_FREE, END0_FREE> const &, Gotoh const &)
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
    (void) upperTriangleEdgeWidth;
    (void) lowerTriangleEdgeWidth;

    SEQAN_ASSERT_FAIL("Write me!");
    return 0;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_DYNPROG_BANDED_AFFINE_H_

