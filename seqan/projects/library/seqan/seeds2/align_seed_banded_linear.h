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
  Author: Carsten Kemena <carsten.kemena@crg.es>
 ============================================================================
  Banded alignment around seeds with linear gap costs.

  This header defines the function bandedAlignment() with tag
  NeedlemanWunsch.  This function first calls
  _bandedAlignment_NW_align() which performs the banded alignment in
  an alignment matrix.  Then, the traceback is computed by
  _bandedAlignment_NW_traceback()
 ==========================================================================*/

#ifndef SEQAN_SEEDS_ALIGN_SEED_BANDED_LINEAR_H_
#define SEQAN_SEEDS_ALIGN_SEED_BANDED_LINEAR_H_

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

// The function bandedAlignment(..., NeedlemanWunsch()) is the public
// interface to the code in this header.
template <typename TSource, typename TAlignSpec, typename TSpecSeed, typename TSeedConfig, typename TBandwidth, typename TScoreValue>
TScoreValue
bandedAlignment(
        Align<TSource, TAlignSpec> & alignment,
        Seed<TSpecSeed, TSeedConfig> const & seed,
        TBandwidth k,
        Score<TScoreValue, Simple> const & scoringScheme,
        NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;

    typedef Matrix<TScoreValue> TMatrix;

    clearGaps(row(alignment, 0));
    clearGaps(row(alignment, 1));

    // First, compute the banded alignment matrix.
    TMatrix matr;
    TScoreValue ret = _bandedAlignment_NW_align(
            matr, seed, k, sourceSegment(row(alignment, 0)),
            sourceSegment(row(alignment, 1)), scoringScheme);

    // Then, compute the traceback through the matrix.
    //
    // We give an iterator into the banded alignment matrix so the
    // traceback function does not have to know the bandwidth and
    // seed.
    Iter<TMatrix, PositionIterator> iter = begin(matr);
    setPosition(iter, 1 + k + getLowerDiagonal(seed) - getStartDiagonal(seed));
    _bandedAlignment_NW_traceback(alignment, matr, iter, scoringScheme);

    return ret;
}


// Compute traceback through the alignment matrix.  Step 1 for
// bandedAlignment(..., NeedlemanWunsch()).  Step 2 is
// _bandedAlignment_NW_traceback().
//
// If this function is used for the banded alignment around one seed,
// initialValues is empty, if it is used for banded chain alignment,
// it is not.

template <typename TScoreValue, typename TSeedConfig, unsigned DIMENSION, typename TSequence, typename TValue, typename TSpecSeed>
TScoreValue
_bandedAlignment_NW_align(
        Matrix<TScoreValue, DIMENSION> & matrix,
        Seed<TSpecSeed, TSeedConfig> const & seed,
        TValue k,
        TSequence const & sequence0,
        TSequence const & sequence1,
        Score<TScoreValue, Simple> const & scoringScheme)
{
    // Case: Not banded chain alignment, plain banded alignment.
    SEQAN_CHECKPOINT;
    String<TScoreValue> emptyInitialValues;
    return _bandedAlignment_NW_align(matrix, seed, k, sequence0, sequence1, scoringScheme, emptyInitialValues);
}


template <typename TScoreValue, typename TSeedConfig, unsigned DIMENSION, typename TSequence, typename TValue, typename TSpecSeed, typename TInitialValuesSpec>
TScoreValue
_bandedAlignment_NW_align(
        Matrix<TScoreValue, DIMENSION> & matrix,
        Seed<TSpecSeed, TSeedConfig> const & seed,
        TValue k,
        TSequence const & sequence0,
        TSequence const & sequence1,
        Score<TScoreValue, Simple> const & scoringScheme,
        String<TScoreValue, TInitialValuesSpec> const & initialValues)
{
    // Case: Banded chain alignment if length(initialValues) > 0.
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_EQ_MSG(scoreGapOpen(scoringScheme), scoreGapExtend(scoringScheme),
                        "Banded Needleman-Wunsch only works with linear gap costs.");

    typedef Matrix<TScoreValue, DIMENSION> TMatrix;

    typedef typename Size<TMatrix>::Type TSize;
    typedef typename Position<TMatrix>::Type TPosition;
    typedef typename Iterator<TMatrix, Standard>::Type TMatrixIterator;

    typedef typename Iterator<TSequence const, Rooted>::Type TSequenceIterator;
    typedef typename Value<TSequence const>::Type TAlphabet;

    // -----------------------------------------------------------------------
    // Declare variables, compute properties of the matrix band
    // -----------------------------------------------------------------------

    // The height of the area between the upper and lower diagonal of the seed.
    TSize height = getLowerDiagonal(seed) - getUpperDiagonal(seed) + 1;

    // The dimensions of the alignment matrix we use for the computation.
    //
    // TODO(holtgrew): The name of the variables *Length are highly misleading.
    TSize sequence0Length = height + 2 * k;
    TSize sequence1Length = length(sequence1);

    // Iterators to the first and last characters of the sequences.
    //
    // TODO(holtgrew): The _end variables should really be called _last!
    // TODO(holtgrew): Also, they should probably be called seq0{Begin,Last}, seq1{Begin,Last}.
    TSequenceIterator x_begin = begin(sequence0);
    TSequenceIterator x_end = end(sequence0) - 1;
    TSequenceIterator y_begin = begin(sequence1);
    --y_begin;
    TSequenceIterator y_end = end(sequence1) - 1;

    // The height of the upper empty triangle.
    TValue upperEmptyTriangleHeight = getLowerDiagonal(seed) - getStartDiagonal(seed) + k;
    // The height and width of the lower empty triangle.
    TValue lowerEmpyTriangleHeight = getEndDiagonal(seed) - getUpperDiagonal(seed) + k;
    TValue lowerEmptyTriangleWidth = sequence0Length - lowerEmpyTriangleHeight;
    // Length of the upper and lower diagonals.
    TSize upperDiagonalLength = sequence1Length - lowerEmpyTriangleHeight;
    TSize lowerDiagonalLength = sequence1Length - upperEmptyTriangleHeight;

    // Get shortcuts to match/mismatch/gap scores.
    TScoreValue matchScore = scoreMatch(scoringScheme);
    TScoreValue mismatchScore = scoreMismatch(scoringScheme);
    TScoreValue gapScore = scoreGapExtend(scoringScheme);

    // TODO(holtgrew): This looks like values for handling the borders of the band.
    TScoreValue horizontalValue = 0;
    TScoreValue borderScore = gapScore;
    TScoreValue verticalValue = borderScore;

    // Set the size of the alignment matrix.
    //
    // TODO(holtgrew): The interface looks not too good, should it not be resizeMatrix(matrix, width, height)?
    setDimension(matrix, 2);
    setLength(matrix, 0, sequence0Length + 2);
    setLength(matrix, 1, sequence1Length + 1);
    fill(matrix, 0);

    // Get iterators into the matrix.
    //
    // TODO(holtgrew): Which of the iterators is for the match/mismatch, horizontal gap, vertical gap case?
    TMatrixIterator finger1;
    TMatrixIterator finger2 = begin(matrix);
    TMatrixIterator finger3;

    // TODO(holtgrew): seq1It and seq0It would be better names.
    TSequenceIterator x = x_begin;
    TSequenceIterator y = y_begin;

    // -----------------------------------------------------------------------
    // Matrix Initialization
    // -----------------------------------------------------------------------

    // TODO(holtgrew): Using a real infimum here can use problems below if anything is subtracted.  This should not be the case...
    TScoreValue inf = InfimumValue<TScoreValue>::VALUE / 2;

    // Fill parts outside the upper-right diagonal.
    setPosition(finger2, length(matrix, 0) - 1);
    for (unsigned i = 1; i != upperDiagonalLength; ++i) {
        *finger2 = inf;
        goNext(finger2, 1);
    }

    TPosition pos = 0;
    borderScore = lowerEmpyTriangleHeight * gapScore;

    // Fill along the cut-off towards the lower right empty triangle.
    *finger2 = inf;
    for (int i = -1; i != lowerEmpyTriangleHeight; ++i){
        goPrevious(finger2, 0);
        goNext(finger2,1);
        if (length(initialValues) > 0) {
            *finger2 = initialValues[pos++];
        } else {
            *finger2 = borderScore;
            borderScore -= gapScore;
        }
    }

    borderScore = gapScore;
    finger3 = finger2;

    // Fill along bottom.
    for (int i = 0; i != lowerEmptyTriangleWidth; ++i){
        goPrevious(finger2, 0);
        if (length(initialValues) > 0) {
            *finger2 = initialValues[pos++];
        } else {
            *finger2 = borderScore;
            borderScore += gapScore;
        }
    }

    *finger2 = inf;
    TScoreValue tmp;

    // -----------------------------------------------------------------------
    // Perform Banded Alignment
    // -----------------------------------------------------------------------
    //
    // TODO(holtgrew): It appears that the alignment is computed from the lower right to the upper left for some reason.  Find out why.
    //
    // Compute part next to the lower empty triangle.
    //
    y = y_end;
    TSize runLength = lowerEmptyTriangleWidth;
    TSize measure = 0;

    // This first loop appears to compute the lower part of the band
    // until the lower empty triangle height is reached.  The lower
    // right corner is empty so handling this separately make sure we
    // do not access this uninitialized memory.
    for (int i = -1; i != lowerEmpyTriangleHeight; ++i) {
        verticalValue = *finger3;
        finger1 = finger3;
        goPrevious(finger1, 0);
        goPrevious(finger3, 1);
        finger2 = finger3;
        goNext(finger3, 0);
        horizontalValue = *finger3;
        x = x_end;
        for (unsigned j = 0; j != runLength; ++j) {
            if (*x==*y) {
                *finger2 = verticalValue + matchScore;
            } else {
                tmp = *finger1;
                TScoreValue s1 = verticalValue + mismatchScore;
                TScoreValue s2 = gapScore + _max(horizontalValue, tmp);
                *finger2 = _max(s1, s2);
            }
            horizontalValue = *finger2;
            goPrevious(finger2);
            verticalValue = *finger1;
            goPrevious(finger1);
            --x;
        }
        *finger2 = inf;
        ++measure;
        if (measure < lowerDiagonalLength)
            ++runLength;
        --y;
    }
    --runLength;

    //
    // Compute part above the lower empty triangle.
    //
    goPrevious(finger3);

    // The second loop now computes the band up from the lower right
    // triangle.  The upper left triangle is not computed.
    while (y != y_begin) {
		verticalValue = *finger3;
		finger1 = finger3;
		goPrevious(finger1, 0);
		goPrevious(finger3, 1);
		finger2 = finger3;
		x = --x_end;
		horizontalValue = inf;
		for (unsigned j = 0; j != runLength; ++j) {
			if (*x == *y) {
				*finger2 = verticalValue + matchScore;
			} else {
				tmp = *finger1;
				TScoreValue s1 = verticalValue + mismatchScore;
				TScoreValue s2 = gapScore + _max(horizontalValue, tmp);
				*finger2 = _max(s1, s2);
			}
			horizontalValue = *finger2;
			goPrevious(finger2);
			verticalValue = *finger1;
			goPrevious(finger1);
			--x;
		}
		*finger2 = inf;
		++measure;
		if (measure >= lowerDiagonalLength)
			--runLength;
		--y;
	}

    // finger2 now points left of the cell with the result.
    finger2 += 1;
	return *finger2;
}


// Compute traceback through the alignment matrix.  Step 2 for
// bandedAlignment(..., NeedlemanWunsch()).  Step 1 is
// _bandedAlignment_NW_align().

template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, unsigned DIMENSION>
typename Size<Matrix<TScoreValue, DIMENSION> >::Type
_bandedAlignment_NW_traceback(
        Align<TTargetSource, TTargetSpec> & alignment,
        Matrix<TScoreValue, DIMENSION> const & matrix,
        Iter<Matrix<TScoreValue, DIMENSION>, PositionIterator> & alignmentMatrixIt,
        Score<TScoreValue, Simple> const & scoringScheme)
{
    // Case: Not banded chain alignment, plain banded alignment.
    SEQAN_CHECKPOINT;
    return _bandedAlignment_NW_traceback(alignment, matrix, alignmentMatrixIt, scoringScheme, false);
}

template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, unsigned DIMENSION>
typename Size<Matrix<TScoreValue, DIMENSION> >::Type
_bandedAlignment_NW_traceback(
        Align<TTargetSource, TTargetSpec> & alignment,
        Matrix<TScoreValue, DIMENSION> const & matrix,
        Iter<Matrix<TScoreValue, DIMENSION>, PositionIterator> & alignmentMatrixIt,
        Score<TScoreValue, Simple> const & scoringScheme,
        bool bandedChainAlignment)
{
    // Case: Banded chain alignment if bandedChainAlignment == true.
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_EQ_MSG(scoreGapOpen(scoringScheme), scoreGapExtend(scoringScheme),
                        "Banded Needleman-Wunsch only works with linear gap costs.");

	typedef Iter<Matrix<TScoreValue, DIMENSION>, PositionIterator > TMatrixIterator;

    typedef typename Infix<TTargetSource>::Type TTargetSourceSegment;
	typedef typename Iterator<TTargetSourceSegment, Rooted>::Type TSequenceIterator;

	typedef Align<TTargetSource, TTargetSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TTargetIterator;

    // -----------------------------------------------------------------------
    // Compute Values Used Below
    // -----------------------------------------------------------------------

	TTargetSourceSegment sequence0 = sourceSegment(row(alignment, 0));
	TTargetSourceSegment sequence1 = sourceSegment(row(alignment, 1));

	TTargetIterator targetSequence0It = iter(row(alignment, 0), 0);
	TTargetIterator targetSequence1It = iter(row(alignment, 1), 0);

	TSequenceIterator itSeq0 = begin(sequence0);
	TSequenceIterator itSeq0End = end(sequence0);

	TSequenceIterator itSeq1 = begin(sequence1);
	TSequenceIterator itSeq1End = end(sequence1);

	TScoreValue scoreDiff = scoreMismatch(scoringScheme) - scoreGap(scoringScheme);

    // -----------------------------------------------------------------------
	// Follow Trace
    // -----------------------------------------------------------------------
    //
    // The trace is followed until the border is reached.
	while (itSeq0 != itSeq0End && itSeq1 != itSeq1End) {
        // The variables gv/gh are set to true if there is a gap to be
        // introduced in sequence1/sequence0.
        //
        // First, compute gv and gh.
		bool gv;
		bool gh;
		if (*itSeq0 == *itSeq1) {
			gv = gh = true;
		} else {
			TMatrixIterator it = alignmentMatrixIt;

			goNext(it, 0);
			TScoreValue h = *it;

			it = alignmentMatrixIt;
			goNext(it, 1);
			TScoreValue d = *it;

			goPrevious(it, 0);
			TScoreValue v = *it;

			gv = (v >= h) || (d + scoreDiff >= h);
			gh = (h >  v) || (d + scoreDiff >= v);
		}
        SEQAN_ASSERT_TRUE(gv || gh);

        // Then, introduce gaps in one of the sequences if necessary
        // and increment the iterators.
		if (gv) {
			++itSeq1;
			goNext(alignmentMatrixIt, 1);
			goPrevious(alignmentMatrixIt, 0);
		} else {
			insertGap(targetSequence1It);
		}
		if (gh) {
			++itSeq0;
			goNext(alignmentMatrixIt, 0);
		} else {
			insertGap(targetSequence0It);
		}

		++targetSequence0It;
		++targetSequence1It;
	}

    if (bandedChainAlignment) {
        setSourceEndPosition(row(alignment, 0), position(itSeq0));
        setSourceEndPosition(row(alignment, 1), position(itSeq1));
        return length(matrix, 0) - coordinate(alignmentMatrixIt, 0) - 2;
    } else {
        return 0;
    }
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_SEED_BANDED_LINEAR_H_
