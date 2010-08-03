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
  Banded alignment around seeds with affine gap costs.

  This header defines the function bandedAlignment() with tag Gotoh.
  This function first calls _bandedAlignment_Gotoh_align() which
  performs the banded alignment in an alignment matrix.  Then, the
  traceback is computed by _bandedAlignment_Gotoh_traceback().
 ==========================================================================*/

#ifndef SEQAN_SEEDS_ALIGN_SEED_BANDED_AFFINE_H_
#define SEQAN_SEEDS_ALIGN_SEED_BANDED_AFFINE_H_

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

// The function bandedAlignment(..., Gotoh()) is the public interface
// to the code in this header.
template <typename TSource, typename TAlignSpec, typename TSpecSeed, typename TSeedConfig, typename TBandwidth, typename TScoreValue>
TScoreValue
bandedAlignment(
        Align<TSource, TAlignSpec> & alignment,
        Seed<TSpecSeed, TSeedConfig> const & seed,
        TBandwidth k,
        Score<TScoreValue, Simple> const & scoringScheme,
        Gotoh const &)
{
    SEQAN_CHECKPOINT;

    typedef Matrix<TScoreValue> TMatrix;

	clearGaps(row(alignment, 0));
	clearGaps(row(alignment, 1));

	TMatrix d_matr;
	TMatrix v_matr;
	TMatrix h_matr;
	TScoreValue ret = _bandedAlignment_Gotoh_align(d_matr, v_matr, h_matr, seed, k, sourceSegment(row(alignment, 0)), sourceSegment(row(alignment, 1)), scoringScheme);
	_bandedAlignment_Gotoh_traceback(alignment, d_matr, v_matr, h_matr, 1 + k + getLowerDiagonal(seed) - getStartDiagonal(seed));
	return ret;
}


// Compute traceback through the alignment matrix.  Step 1 for
// bandedAlignment(..., Gotoh()).  Step 2 is
// _bandedAlignment_Gotoh_traceback().
template <typename TScoreValue, unsigned DIMENSION, typename TString, typename TValue, typename TSpecSeed, typename TSeedConfig>
TScoreValue
_bandedAlignment_Gotoh_align(
        Matrix<TScoreValue, DIMENSION> & diag_matrix_,
        Matrix<TScoreValue, DIMENSION> & vert_matrix_,
        Matrix<TScoreValue, DIMENSION> & hori_matrix_,
        Seed<TSpecSeed, TSeedConfig> const & seed,
        TValue k,
        TString const & str1_,
        TString const & str2_,
        Score<TScoreValue, Simple> const & scoringScheme)
{
    SEQAN_CHECKPOINT;

    // TODO(holtgrew): This function could use some polishing.
    // TODO(holtgrew): More comments! Seed NW alignment for an example.
    // TODO(holtgrew): Variable names should follow the SeqAn conventions.

	typedef Matrix<TScoreValue, DIMENSION> TMatrix;

	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Iterator<TMatrix, Rooted>::Type TMatrixIterator;

	typedef typename Iterator<TString const, Rooted>::Type TStringIterator;
	typedef typename Value<TString const>::Type TAlphabet;

    // -----------------------------------------------------------------------
    // Declare variables, compute properties of the matrix band
    // -----------------------------------------------------------------------

	TSize height = getLowerDiagonal(seed) - getUpperDiagonal(seed) + 1;
	TSize str1_length = height+2*k;
	TSize str2_length = length(str2_);

    // The height of the upper empty triangle.
	TValue up_height = getLowerDiagonal(seed) - getStartDiagonal(seed) + k;
    // The height and width of the lower empty triangle.
	TValue down_height = getEndDiagonal(seed) - getUpperDiagonal(seed) + k;
	TValue down_width = str1_length - down_height;
    // Length of the upper and lower diagonals.
	TValue length_right_diag = str2_length - down_height;
	TValue length_left_diag = str2_length - up_height;

    // TODO(holtgrew): Again, _end should probably be called _last
	TStringIterator x_begin = begin(str1_) - 1;
	TStringIterator x_end = end(str1_) - 1;
	TStringIterator y_begin = begin(str2_) - 1;
	TStringIterator y_end = end(str2_) - 1;

    // Get shortcuts to the match, mismatch, gap opening and extension
    // scores.
	TScoreValue matchScore = scoreMatch(scoringScheme);
	TScoreValue mismatchScore = scoreMismatch(scoringScheme);
	TScoreValue gapOpenScore = scoreGapOpen(scoringScheme);
	TScoreValue gapExtendScore = scoreGapExtend(scoringScheme);

    // Set the size of the alignment matrix.
    //
    // TODO(holtgrew): The interface looks not too good, should it not be resizeMatrix(matrix, width, height)?
	setDimension(diag_matrix_, 2);
	setLength(diag_matrix_, 0, str1_length + 2);
	setLength(diag_matrix_, 1, str2_length + 1);
	resize(diag_matrix_);
	setDimension(vert_matrix_, 2);
	setLength(vert_matrix_, 0, str1_length + 2);
	setLength(vert_matrix_, 1, str2_length + 1);
	resize(vert_matrix_);
	setDimension(hori_matrix_, 2);
	setLength(hori_matrix_, 0, str1_length + 2);
	setLength(hori_matrix_, 1, str2_length + 1);
	resize(hori_matrix_);

    // Get iterators into the matrices.
	TSize pos = length(diag_matrix_, 0) - 1;
	TMatrixIterator diag_col_ = begin(diag_matrix_);
	TMatrixIterator diag_finger1(diag_matrix_,pos);
	TMatrixIterator diag_finger2;
	TMatrixIterator diag_finger3;
	TMatrixIterator vert_col_ = begin(vert_matrix_);
	TMatrixIterator vert_finger1(vert_matrix_,pos);
	TMatrixIterator vert_finger2;
	TMatrixIterator vert_finger3;
	TMatrixIterator hori_col_ = begin(hori_matrix_);
	TMatrixIterator hori_finger1(hori_matrix_,pos);
	TMatrixIterator hori_finger2;
	TMatrixIterator hori_finger3;

    // TODO(holtgrew): seq1It and seq0It would be better names.
	TStringIterator x = x_end;
	TStringIterator y;

    // -----------------------------------------------------------------------
    // Matrix Initialization
    // -----------------------------------------------------------------------

    // TODO(holtgrew): Using a real infimum here can use problems below if anything is subtracted.
	TScoreValue inf = InfimumValue<TScoreValue>::VALUE / 2;
	for (int i = 1; i != length_right_diag; ++i){
		*diag_finger1 = inf;
		*hori_finger1 = inf;
		*vert_finger1 = inf;
		goNext(diag_finger1, 1);
		goNext(hori_finger1, 1);
		goNext(vert_finger1, 1);
	}

	*diag_finger1 = inf;
	*hori_finger1 = inf;
	*vert_finger1 = inf;

	TScoreValue border_ = gapOpenScore + (down_height-1) * gapExtendScore;

	for (int i = -1; i != down_height; ++i){
		goPrevious(diag_finger1, 0);
		goNext(diag_finger1,1);
		goPrevious(vert_finger1, 0);
		goNext(vert_finger1,1);
		goPrevious(hori_finger1, 0);
		goNext(hori_finger1,1);
		*diag_finger1 = border_;
		*vert_finger1 = border_;
		*hori_finger1 = inf;
		border_ -= gapExtendScore;
	}

	*diag_finger1 = 0;
	*vert_finger1 = inf;
	*hori_finger1 = inf;

	diag_finger3 = diag_finger1;
	hori_finger3 = hori_finger1;
	vert_finger3 = vert_finger1;

	for (int i = 0; i < down_width; ++i){
		goPrevious(diag_finger1, 0);
		goPrevious(hori_finger1, 0);
		goPrevious(vert_finger1, 0);
		border_ += gapExtendScore;
		*diag_finger1 = border_;
		*hori_finger1 = border_;
		*vert_finger1 = inf;
	}

	*diag_finger1 = inf;
	*hori_finger1 = inf;

    // -----------------------------------------------------------------------
    // Perform Banded Alignment
    // -----------------------------------------------------------------------
    //
    // TODO(holtgrew): It appears that the alignment is computed from the lower right to the upper left for some reason.  Find out why.
    //
    // Compute part next to the lower empty triangle.
    //
	TScoreValue hori_value, vert_value, diag_value, tmp_diag, tmp_vert, diag_front, diag_match, diag_back;
	y=y_end;
	TSize run_length = down_width;
	TSize measure = 0;

    // TODO(holtgrew): This first loop appears to compute the lower part of the band until the lower empty triangle height is reached.
	for (int i = -1; i != down_height; ++i) {
		goPrevious(hori_finger3,1);
		hori_finger1 = hori_finger3;
		goNext(hori_finger3);
		hori_value = *hori_finger3;

		vert_finger2 = vert_finger3;
		goPrevious(vert_finger2);
		goPrevious(vert_finger3,1);
		vert_finger1 = vert_finger3;
		goNext(vert_finger3);

		diag_finger2 = diag_finger3;
		diag_match=*diag_finger2;
		goPrevious(diag_finger2);
		goPrevious(diag_finger3,1);
		diag_finger1 = diag_finger3;
		goNext(diag_finger3);
		diag_back=*diag_finger3;

		x = x_end;
		for (unsigned int j = 0; j != run_length; ++j){
			diag_front= *diag_finger2;
			hori_value = (hori_value+gapExtendScore > diag_back+gapOpenScore) ? hori_value+gapExtendScore : diag_back+gapOpenScore;
			*hori_finger1 = hori_value;
			goPrevious(hori_finger1);

			tmp_vert = *vert_finger2;
			vert_value = (tmp_vert+gapExtendScore > diag_front+gapOpenScore) ? tmp_vert+gapExtendScore : diag_front+gapOpenScore;
			*vert_finger1 = vert_value;
			goPrevious(vert_finger1);
			goPrevious(vert_finger2);

			tmp_diag = (vert_value > hori_value) ? vert_value  : hori_value;
			diag_value = diag_match + ((*x == *y) ? matchScore : mismatchScore);
			if (diag_value < tmp_diag)
				diag_value = tmp_diag;
			diag_back = diag_value;
			*diag_finger1 = diag_value;
			goPrevious(diag_finger1);
			diag_match = *diag_finger2;
			goPrevious(diag_finger2);
			--x;
		}
		*diag_finger1 = inf;
		*vert_finger1 = inf;
		++measure;
		if (static_cast<TValue>(measure) < length_left_diag)
			++run_length;
		--y;
	}
	--run_length;

    //
    // Compute part above the lower empty triangle.
    //
    // TODO(holtgrew): The second loop now computes the band up from the lower right triangle. It is not clear to me whether the upper left empty triangle is computed or not.
	while (y != y_begin) {
		goPrevious(diag_finger3);
		goPrevious(hori_finger3);
		goPrevious(vert_finger3);

		goPrevious(hori_finger3,1);
		hori_finger1 = hori_finger3;
		goNext(hori_finger3);
		hori_value = inf;

		vert_finger2 = vert_finger3;
		goPrevious(vert_finger2);
		goPrevious(vert_finger3,1);
		vert_finger1 = vert_finger3;
		goNext(vert_finger3);

		diag_finger2 = diag_finger3;
		diag_match=*diag_finger2;
		goPrevious(diag_finger2);
		goPrevious(diag_finger3,1);
		diag_finger1 = diag_finger3;
		goNext(diag_finger3);
		diag_back= inf;
		x = --x_end;
		for (unsigned int j = 0; j != run_length; ++j){
			diag_front= *diag_finger2;
			hori_value = (hori_value+gapExtendScore > diag_back+gapOpenScore) ? hori_value+gapExtendScore : diag_back+gapOpenScore;
			*hori_finger1 = hori_value;
			goPrevious(hori_finger1);

			tmp_vert = *vert_finger2;
			vert_value = (tmp_vert+gapExtendScore > diag_front+gapOpenScore) ? tmp_vert+gapExtendScore : diag_front+gapOpenScore;
			*vert_finger1 = vert_value;
			goPrevious(vert_finger1);
			goPrevious(vert_finger2);

			tmp_diag = (vert_value > hori_value) ? vert_value  : hori_value;
			diag_value = diag_match + ((*x == *y) ? matchScore : mismatchScore);
			if (diag_value < tmp_diag)
				diag_value = tmp_diag;
			diag_back = diag_value;
			*diag_finger1 = diag_value;
			goPrevious(diag_finger1);
			diag_match = *diag_finger2;
			goPrevious(diag_finger2);
			--x;
		}
		*vert_finger1 = inf;
		*diag_finger1 = inf;
		++measure;
		if (static_cast<TValue>(measure) >= length_left_diag)
			--run_length;
		--y;
	}
	--run_length;
	++diag_finger1;

	return *diag_finger1;
}


// Compute traceback through the alignment matrix.  Step 2 for
// bandedAlignment(..., Gotoh()).  Step 1 is
// _bandedAlignment_Gotoh_align().
template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, unsigned DIMENSION, typename TPos>
void
_bandedAlignment_Gotoh_traceback(
        Align<TTargetSource, TTargetSpec> & alignment,
        // TODO(holtgrew): Why can't the matrices be const? Something within getValue() does not compile if made const.
        Matrix<TScoreValue, DIMENSION> /*const*/ & diagMatrix,
        Matrix<TScoreValue, DIMENSION> /*const*/ & vertMatrix,
        Matrix<TScoreValue, DIMENSION> /*const*/ & horiMatrix,
        TPos position)
{
    SEQAN_CHECKPOINT;

	typedef typename Position<Matrix<TScoreValue, DIMENSION> >::Type TPosition;

	typedef Align<TTargetSource, TTargetSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TTargetIterator;

	typedef typename Iterator<TTargetSource, Standard>::Type TStringIterator;
    typedef typename Infix<TTargetSource>::Type TTargetSourceSegment;
    typedef typename Size<TTargetSourceSegment>::Type TSize;

    // -----------------------------------------------------------------------
    // Compute Values Used Below
    // -----------------------------------------------------------------------

	TTargetSourceSegment sequence0 = sourceSegment(row(alignment, 0));
	TTargetSourceSegment sequence1 = sourceSegment(row(alignment, 1));
	TSize lengthDim0 = length(horiMatrix, 0);

	TPosition pos = position;

	TTargetIterator targetSequence0It = iter(row(alignment, 0), 0);
	TTargetIterator targetSequence1It = iter(row(alignment, 1), 0);

	TStringIterator itSeq0 = iter(sequence0, 0);
	TStringIterator itSeq0End = end(sequence0);

	TStringIterator itSeq1 = iter(sequence1, 0);
	TStringIterator it1End = end(sequence1);

    // -----------------------------------------------------------------------
	// Follow Trace
    // -----------------------------------------------------------------------
    //
    // The trace is followed until the border is reached.
	while (itSeq0 != itSeq0End && itSeq1 != it1End) {
		if (getValue(diagMatrix, pos) > getValue(horiMatrix, pos)) {
			if (getValue(diagMatrix, pos) > getValue(vertMatrix, pos)) {
				++itSeq0;
				++itSeq1;
				pos += lengthDim0;
			} else {
				++itSeq1;
				insertGap(targetSequence0It);
				pos += lengthDim0 - 1;
			}
		} else {
			if (getValue(horiMatrix, pos) > getValue(vertMatrix, pos)) {
				++itSeq0;
				insertGap(targetSequence1It);
				++pos;
			} else {
				++itSeq1;
				insertGap(targetSequence0It);
				pos += lengthDim0 - 1;
			}
		}
		++targetSequence0It;
		++targetSequence1It;
	}
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_SEED_BANDED_AFFINE_H_
