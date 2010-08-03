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
_alignBanded_resizeMatrix(Matrix<TScoreValue, 2> & matrix, TSequence const & sequence0, TSequence const & sequence1, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, NeedlemanWunsch const &)
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
    resize(matrix);
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


// Compute traceback in the given banded alignment matrix, starting at
// the lower right, shifted by the given overlap to the upper left.
// The matrix was filled with the Needleman-Wunschl algorithm, end
// gaps are free as configured by the AlignConfig object.  Returns the
// best score.
template <typename TAlignmentIterator, typename TSequenceIterator, typename TPosition, typename TScoreValue, typename TScoringScheme, typename TOverlap, bool START0_FREE, bool START1_FREE, bool END0_FREE, bool END1_FREE>
TScoreValue
_alignBanded_traceBack(TAlignmentIterator & alignmentIt0, TAlignmentIterator & alignmentIt1, TSequenceIterator & sourceIt0, TSequenceIterator & sourceIt1, TPosition & finalPos0, TPosition & finalPos1, Matrix<TScoreValue, 2> /*const*/ & matrix, TScoringScheme const & scoringScheme, TOverlap overlap0, TOverlap overlap1, TOverlap upperTriangleEdgeLength, TOverlap lowerTriangleEdgeLength, bool goToTopLeft, AlignConfig<START1_FREE, START0_FREE, END1_FREE, END0_FREE> const &, NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;

    // std::cout << "trace back banded" << std::endl;
    // // TODO(holtgrew): Debug output, remove when not needed any more.
    // {
    //     std::cout << ",-- filled banded alignment matrix " << length(matrix, 0) << " x " << length(matrix, 1) << std::endl;
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
    // std::cout << "length(matrix, 0) = " << length(matrix, 0) << std::endl;
    // std::cout << "overlap0 = " << overlap0 << std::endl;
    // std::cout << "finalPos0 = " << finalPos0 << std::endl;
    // std::cout << "length(matrix, 1) = " << length(matrix, 1) << std::endl;
    // std::cout << "overlap1 = " << overlap1 << std::endl;
    // std::cout << "finalPos1 = " << finalPos1 << std::endl;
    // std::cout << "upperTriangleEdgeLength = " << upperTriangleEdgeLength << std::endl;
    // std::cout << "lowerTriangleEdgeLength = " << lowerTriangleEdgeLength << std::endl;
    
    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Iterator<TMatrix>::Type TMatrixIterator;

    // Initialization
    //
    // Precomputation of the score difference between mismatch and gap.
	TScoreValue scoreDifference = scoreMismatch(scoringScheme) - scoreGap(scoringScheme);
    // Current position in the matrix.  Note that this is the position
    // in the matrix including the gutter.  When writing this out, we
    // want to have the position in the matrix, excluding the gutter,
    // i.e. everything is shifted to the upper left.
    TPosition pos0 = length(matrix, 0) - overlap0 + finalPos0;
    TPosition pos1 = length(matrix, 1) - overlap1 - lowerTriangleEdgeLength - 1 + (length(matrix, 0) - pos0 - 1) + finalPos1;
    // std::cout << "length(matrix, 0) - pos0 - 1 == " << (length(matrix, 0) - pos0 - 1) << std::endl;
    // std::cout << "begin pos0 = " << pos0 << std::endl;
    // std::cout << "begin pos1 = " << pos1 << std::endl;
    TPosition origPos0 = pos0;
    (void) origPos0;  // In case we run without assertions.
    // TPosition origPos1 = pos1;
    // Iterator to current entry in the matrix.
    TMatrixIterator matrixIt = begin(matrix);
    setPosition(matrixIt, pos0 + pos1 * _dataFactors(matrix)[1]);  // TODO(holtgrew): Matrix class should have setPosition with coordinates.
    TMatrixIterator origMatrixIt = matrixIt;

    // std::cout << "STARTING AT " << pos0 << ", " << pos1 << std::endl;

    // Search for starting point of the trace.
    //
    // Insert free end gaps if any.
    //
    // TODO(holtgrew): Currently not supported.

    // Now, perform the traceback.
    while (pos0 > static_cast<TPosition>(1)) {
        SEQAN_ASSERT_GT(pos1, static_cast<TPosition>(0));
        // Flags for "go horizontal" and "go vertical".
        bool gh = false;
        bool gv = false;

        // Determine whether to go vertical/horizontal.
        // std::cout << "Compare: " << *sourceIt0 << ", " << *sourceIt1 << std::endl;
        // std::cout << "  Score: " << *matrixIt << std::endl;
        if (*sourceIt0 == *sourceIt1) {
            gh = true;
            gv = true;
        } else {
			TMatrixIterator it = matrixIt;

			goPrevious(it, 1);
			TScoreValue h = *it;

			it = matrixIt;
			goPrevious(it, 1);
			TScoreValue d = *it;

			goNext(it, 1);
			TScoreValue v = *it;

			gv = (v >= h) || (d + scoreDifference >= h);
			gh = (h > v) || (d + scoreDifference >= v);
        }

        // if (gv && gh) {
        //     std::cout << "GO DIAGONAL" << std::endl;
        // } else if (gv) {
        //     std::cout << "GO VERTICAL" << std::endl;
        // } else if (gh) {
        //     std::cout << "GO HORIZONTAL" << std::endl;
        // }

        // Move iterators in source sequence, alignment rows, matrix
        // and possibly insert gaps.
        if (gv && gh) {
            pos0 -= 1;
            goPrevious(matrixIt, 0);
        } else if (gv) {
            pos0 -= 1;
            pos1 += 1;
            goPrevious(matrixIt, 0);
            goNext(matrixIt, 1);
        } else if (gh) {
            pos1 -= 1;
            goPrevious(matrixIt, 1);
        }
        if (gv) {
            // std::cout << "  moving in dimension 0" << std::endl;
            goPrevious(sourceIt0);
            goPrevious(alignmentIt0);
            // std::cout << "  *alignmentIt0 == " << convert<char>(*alignmentIt0) << std::endl;
        } else {
            // std::cout << __FILE__ << ":" << __LINE__ << "--  gap in dimension 0" << std::endl;
            insertGap(alignmentIt0);
        }
        if (gh) {
            // std::cout << "  moving in dimension 1" << std::endl;
            goPrevious(sourceIt1);
            goPrevious(alignmentIt1);
            // std::cout << "  *alignmentIt1 == " << convert<char>(*alignmentIt1) << std::endl;
        } else {
            // std::cout << __FILE__ << ":" << __LINE__ << "--  gap in dimension 0" << std::endl;
            insertGap(alignmentIt1);
        }
    }

    // Go to the top left of the matrix if configured to do so.
    //
    // TODO(holtgrew): Currently not supported.
    SEQAN_ASSERT_NOT_MSG(goToTopLeft, "goToTopLeft is currently not supported in banded alignment.");

    // Write out the final positions in the alignment matrix.  Convert
    // from position in the current alignment matrix with gutter to
    // positioin without gutter by shifting the position to the upper
    // left.  Also adjust for the shear of the matrix.
    // std::cout << "@and: pos0 = " << pos0 << ", pos1 = " << pos1 << std::endl;
    finalPos0 = pos0 - 1;
    finalPos1 = pos1 - 2 - (upperTriangleEdgeLength - 1);
    SEQAN_ASSERT_EQ_MSG(finalPos0, 0u, "Must reach top of matrix in banded alignment traceback.");
    // std::cout << "finalPos0 = " << finalPos0 << ", finalPos1 = " << finalPos1 << std::endl;

    return *origMatrixIt;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_DYNPROG_BANDED_LINEAR_H_

