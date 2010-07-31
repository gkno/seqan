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

template <typename TScoreValue, typename TSequence, typename TDiagonal>
inline void
_alignBanded_resizeMatrix(Matrix<TScoreValue, 2> & matrix, TSequence const & sequence0, TSequence const & /*sequence1*/, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;

    // We need space length(sequence0) x bandwidth for the alignment,
    // the top gutter and the left and right gutters.
    setLength(matrix, 0, length(sequence0) + 1);
    setLength(matrix, 1, upperDiagonal - lowerDiagonal + 3);
    fill(matrix, -42);
    // resize(matrix);
}


template <typename TScoreValue, typename TDiagonal, bool BEGIN1_FREE, bool BEGIN0_FREE, bool END1_FREE, bool END0_FREE>
inline void
_alignBanded_initGutter(Matrix<TScoreValue, 2> & matrix, Score<TScoreValue, Simple> const & scoringScheme, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, AlignConfig<BEGIN1_FREE, BEGIN0_FREE, END1_FREE, END0_FREE> const &, NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_EQ_MSG(scoreGapOpen(scoringScheme), scoreGapExtend(scoringScheme),
                        "Only linear gap costs allowed for Needleman-Wunsch.");

    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Iterator<TMatrix>::Type TIterator;
    typedef typename Position<TMatrix>::Type TPosition;
    
    // Initialize the diagonal below the lower one with infimas.
    TIterator it = begin(matrix);
    for (TPosition i = 0, iend = length(matrix, 0); i < iend; ++i, goNext(it, 0))
        *it = InfimumValue<TScoreValue>::VALUE / 2;
    // Initialize the left gutter according to the AlignConfig.
    setPosition(it, (1 - lowerDiagonal) * _dataFactors(matrix)[1]); // TODO(holtgrew): There should be a function that accepts two coordinates for the Matrix class.
    if (BEGIN0_FREE) {
        for (TPosition i = 0, iend = -lowerDiagonal; i <= iend; ++i, goNext(it, 0), goPrevious(it, 1))
            *it = 0;
    } else {
        TScoreValue gapScore = scoreGap(scoringScheme);
        TScoreValue x = 0;
        for (TPosition i = 0, iend = -lowerDiagonal; i <= iend; ++i, goNext(it, 0), goPrevious(it, 1)) {
            *it = x;
            x += gapScore;
        }
    }
    // Initialize the top gutter according to the AlignConfig.
    setPosition(it, (1 - lowerDiagonal) * _dataFactors(matrix)[1]); // TODO(holtgrew): There should be a function that accepts two coordinates for the Matrix class.
    if (BEGIN1_FREE) {
        for (TPosition i = 0, iend = upperDiagonal; i <= iend; ++i, goNext(it, 1))
            *it = 0;
    } else {
        TScoreValue gapScore = scoreGap(scoringScheme);
        TScoreValue x = 0;
        for (TPosition i = 0, iend = upperDiagonal; i <= iend; ++i, goNext(it, 1)) {
            *it = x;
            x += gapScore;
        }
    }
    // Initialize the diagonal above the upper one with infimas.
    for (TPosition i = 0, iend = length(matrix, 0); i < iend; ++i, goNext(it, 0))
        *it = InfimumValue<TScoreValue>::VALUE / 2;
}


template <typename TScoreValue, typename TDiagonal, typename TOverlap>
inline void
_alignBanded_initGutterFromUnbanded(Matrix<TScoreValue, 2> & matrix, Score<TScoreValue, Simple> const & scoringScheme, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, Matrix<TScoreValue, 2> /*const*/ & otherMatrix, TOverlap overlap0, TOverlap overlap1, NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;

    // TODO(holtgrew): Really unnecessary? Remove along with all other unused parameters in all align_*.h files.
    (void) upperDiagonal;

    SEQAN_ASSERT_EQ_MSG(scoreGapOpen(scoringScheme), scoreGapExtend(scoringScheme),
                        "Only linear gap costs allowed for Needleman-Wunsch.");

    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Iterator<TMatrix>::Type TIterator;
    typedef typename Position<TMatrix>::Type TPosition;
    
    // Initialize the diagonal below the lower one with infimas.
    TIterator it = begin(matrix);
    for (TPosition i = 0, iend = length(matrix, 0); i < iend; ++i, goNext(it, 0))
        *it = InfimumValue<TScoreValue>::VALUE / 2;
    // // TODO(holtgrew): Debug code.
    // std::cout << ",-- NW Matrix to copy in... " << length(otherMatrix, 0) << " x " << length(otherMatrix, 1) << std::endl;
    // for (unsigned i = 0; i < length(otherMatrix, 0); ++i) {
    //     std::cout << "|\t";
    //     for (unsigned j = 0; j < length(otherMatrix, 1); ++j) {
    //         std::cout << value(otherMatrix, i, j) << "\t";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "`--" << std::endl;
    // Copy over a column from the other matrix into the left gutter.
    // std::cout << "overlap0 = " << overlap0 << " overlap1 = " << overlap1 << std::endl;
    // std::cout << "pos = " << (length(otherMatrix, 0) - overlap0 - 1) + (length(otherMatrix, 1) - overlap1 - 1) * _dataFactors(otherMatrix)[1] << std::endl;
    TIterator otherIt = begin(otherMatrix);
    setPosition(otherIt, (length(otherMatrix, 0) - overlap0 - 1) + (length(otherMatrix, 1) - overlap1 - 1) * _dataFactors(otherMatrix)[1]);
    TIterator srcIt = otherIt;
    setPosition(it, (1 - lowerDiagonal) * _dataFactors(matrix)[1]); // TODO(holtgrew): There should be a function that accepts two coordinates for the Matrix class.
    for (TOverlap i = 0; i < overlap0; ++i) {
        // std::cout << "*srcIt = " << *srcIt << std::endl;
        *it = *srcIt;
        goNext(srcIt, 0);
        goNext(it, 0);
        goPrevious(it, 1);
    }
    // Copy over a row from the other matrix into the top gutter.
    srcIt = otherIt;
    setPosition(it, (1 - lowerDiagonal) * _dataFactors(matrix)[1]); // TODO(holtgrew): There should be a function that accepts two coordinates for the Matrix class.
    for (TOverlap i = 0; i < overlap1; ++i) {
        // std::cout << "*srcIt = " << *srcIt << std::endl;
        *it = *srcIt;
        goNext(srcIt, 1);
        goNext(it, 1);
    }
    // Initialize the diagonal above the upper one with infimas.
    for (TPosition i = 0, iend = length(matrix, 0); i < iend; ++i, goNext(it, 0))
        *it = InfimumValue<TScoreValue>::VALUE / 2;
}


template <typename TScoreValue, typename TSequence, typename TDiagonal>
inline void
_alignBanded_fillMatrix(Matrix<TScoreValue, 2> & matrix, TSequence const & sequence0, TSequence const & sequence1, Score<TScoreValue, Simple> const & scoringScheme, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_EQ_MSG(scoreGapOpen(scoringScheme), scoreGapExtend(scoringScheme),
                        "Only linear gap costs allowed for Needleman-Wunsch.");
    SEQAN_ASSERT_GEQ_MSG(upperDiagonal, 0, "Upper diagonal must not lie below main diagonal.");
    SEQAN_ASSERT_LEQ_MSG(lowerDiagonal, 0, "Upper diagonal must not lie above main diagonal.");

    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Iterator<TMatrix>::Type TMatrixIterator;
    typedef typename Iterator<TSequence, Standard>::Type TSequenceIterator;
    typedef typename Size<TSequence>::Type TSize;
    typedef typename Position<TSequence>::Type TPosition;

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

    // Fill matrix column wise for cache efficiency.  For this, we
    // need the length of the current column in the sheared alignment
    // matrix which complicates things.
    setPosition(itTop, (1 - lowerDiagonal) * _dataFactors(matrix)[1]); // TODO(holtgrew): There should be a function that accepts two coordinates for the Matrix class.
    TPosition seq1Pos = 0;
    TSequenceIterator it0Begin = begin(sequence0);
    for (TSequenceIterator it1 = begin(sequence1), it1end = end(sequence1); it1 != it1end; ++it1, ++seq1Pos) {
        if (seq1Pos <= static_cast<TPosition>(upperDiagonal)) {
            itLeft = itTop;
            goNext(itTop, 1);
        } else {
            goNext(itTop, 0);
            itLeft = itTop;
            goPrevious(itLeft, 1);
        }
        itAbove = itTop;
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
            // std::cout << "Compare " << *it0 << " and " << *it1 << std::endl;
            TScoreValue scoreMoveDiagonal = *itLeft + ((*it0 == *it1) ? matchScore : mismatchScore);
            // std::cout << "*itLeft == " << *itLeft << std::endl;
            goNext(itLeft, 0);
            goPrevious(itLeft, 1);
            TScoreValue scoreMoveRight = *itLeft + gapScore;
            TScoreValue scoreMoveDown = *itAbove + gapScore;
            // std::cout << "scoreMoveDiagonal == " << scoreMoveDiagonal << ", scoreMoveRight == " << scoreMoveRight << ", scoreMoveDown == " << scoreMoveDown << std::endl;
            // std::cout << "*itLeft == " << *itLeft << ", *itAbove == " << *itAbove << std::endl;
            goNext(itAbove, 0);
            goPrevious(itAbove, 1);
            *itAbove = _max(scoreMoveDiagonal, _max(scoreMoveRight, scoreMoveDown));
            // // TODO(holtgrew): Debug output, remove when not needed any more.
            // {
            //     std::cout << ",-- matrix while filling" << std::endl;
            //     for (unsigned i = 0; i < length(matrix, 0); ++i) {
            //         std::cout << "| ";
            //         for (unsigned j = 0; j < i; ++j)
            //             std::cout << "\t";
            //         for (unsigned j = 0; j < length(matrix, 1); ++j) {
            //             if (value(matrix, i, j) == InfimumValue<int>::VALUE / 2)
            //                 std::cout << "\tinf";
            //             else
            //                 std::cout << "\t" << value(matrix, i, j);
            //         }
            //         std::cout << std::endl;
            //     }
            //     std::cout << "`--" << std::endl;
            // }
        }
        // std::cout << "++it1" << std::endl;
        // We only need an offset for it0 when all diagonals up to the
        // upper one are aligned.
        if (seq1Pos >= static_cast<TPosition>(upperDiagonal))
            it0Begin += 1;
    }
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_DYNPROG_BANDED_LINEAR_H_

