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
	TScoreValue ret = _banded_gotoh(d_matr, v_matr, h_matr, seed, k, sourceSegment(row(alignment, 0)), sourceSegment(row(alignment, 1)), scoringScheme);
	_banded_gotoh_trace(alignment, d_matr, v_matr, h_matr, 1 + k + getLowerDiagonal(seed) - getStartDiagonal(seed));	
	return ret;
}
 
   
// Compute traceback through the alignment matrix.  Step 1 for
// bandedAlignment(..., NeedlemanWunsch()).  Step 2 is
// _bandedAlignment_NW_traceback().
template <typename TScoreValue, unsigned DIMENSION, typename TString, typename TValue, typename TSpecSeed, typename TSeedConfig>
TScoreValue
_banded_gotoh(
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
	typedef Matrix<TScoreValue, DIMENSION> TMatrix;

	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Iterator<TMatrix, Rooted>::Type TMatrixIterator;

	typedef typename Iterator<TString const, Rooted>::Type TStringIterator;
	typedef typename Value<TString const>::Type TAlphabet;

	//-------------------------------------------------------------------------
	//define some variables

	TSize height = getLowerDiagonal(seed) - getUpperDiagonal(seed) + 1;
	TSize str1_length = height+2*k;
	TSize str2_length = length(str2_);

	TValue up_height = getLowerDiagonal(seed) - getStartDiagonal(seed) + k; //equals length of empty lower triangle
	//TValue up_width = str1_length - up_height;
	TValue down_height = getEndDiagonal(seed) - getUpperDiagonal(seed) + k; //equals length of empty lower triangle
	TValue down_width = str1_length - down_height;
	TValue length_right_diag = str2_length - down_height;
	TValue length_left_diag = str2_length - up_height;


	TStringIterator x_begin = begin(str1_) - 1;
	TStringIterator x_end = end(str1_) - 1;
	TStringIterator y_begin = begin(str2_) - 1;
	TStringIterator y_end = end(str2_) - 1;

	TStringIterator x = x_end;
	TStringIterator y;

	TScoreValue scoringSchemematch = scoreMatch(scoringScheme);
	TScoreValue scoringSchememismatch = scoreMismatch(scoringScheme);
	TScoreValue scoringSchemegap_open = scoreGapOpen(scoringScheme);
	TScoreValue scoringSchemegap_extend = scoreGapExtend(scoringScheme);

	//TScoreValue v;

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

	TSize pos = length(diag_matrix_, 0)-1;
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

	
	//-------------------------------------------------------------------------
	// init
	
    // TODO(holtgrew): Using a real infimum here can use problems below if anything is subtracted.  This should not be the case...
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

	TScoreValue border_ = scoringSchemegap_open + (down_height-1) * scoringSchemegap_extend;

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
		border_ -= scoringSchemegap_extend;
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
		border_ += scoringSchemegap_extend;
		*diag_finger1 = border_;
		*hori_finger1 = border_;
		*vert_finger1 = inf;
	}

	*diag_finger1 = inf;
	*hori_finger1 = inf;
	
	//-------------------------------------------------------------------------
	//first section
	TScoreValue hori_value, vert_value, diag_value, tmp_diag, tmp_vert, diag_front, diag_match, diag_back;
	y=y_end;
	TSize run_length = down_width;
	TSize measure = 0;

	for (int i = -1; i != down_height; ++i){	
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
			hori_value = (hori_value+scoringSchemegap_extend > diag_back+scoringSchemegap_open) ? hori_value+scoringSchemegap_extend : diag_back+scoringSchemegap_open;
			*hori_finger1 = hori_value;
			goPrevious(hori_finger1);

			tmp_vert = *vert_finger2;
			vert_value = (tmp_vert+scoringSchemegap_extend > diag_front+scoringSchemegap_open) ? tmp_vert+scoringSchemegap_extend : diag_front+scoringSchemegap_open;
			*vert_finger1 = vert_value;
			goPrevious(vert_finger1);
			goPrevious(vert_finger2);

			
			tmp_diag = (vert_value > hori_value) ? vert_value  : hori_value;
			diag_value = diag_match + ((*x == *y) ? scoringSchemematch : scoringSchememismatch);
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

	while(y!= y_begin)//for (int i = 0; i != main; ++i){
	{
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
			hori_value = (hori_value+scoringSchemegap_extend > diag_back+scoringSchemegap_open) ? hori_value+scoringSchemegap_extend : diag_back+scoringSchemegap_open;
			*hori_finger1 = hori_value;
			goPrevious(hori_finger1);

			tmp_vert = *vert_finger2;
			vert_value = (tmp_vert+scoringSchemegap_extend > diag_front+scoringSchemegap_open) ? tmp_vert+scoringSchemegap_extend : diag_front+scoringSchemegap_open;
			*vert_finger1 = vert_value;
			goPrevious(vert_finger1);
			goPrevious(vert_finger2);

			tmp_diag = (vert_value > hori_value) ? vert_value  : hori_value;
			diag_value = diag_match + ((*x == *y) ? scoringSchemematch : scoringSchememismatch);
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
// bandedAlignment(..., NeedlemanWunsch()).  Step 1 is
// _bandedAlignment_NW_align().
template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, unsigned DIMENSION, typename TValue>
void
_banded_gotoh_trace(
        Align<TTargetSource, TTargetSpec> & target_,
        Matrix<TScoreValue, DIMENSION> & diag_matrix_,
        Matrix<TScoreValue, DIMENSION> & vert_matrix_,
        Matrix<TScoreValue, DIMENSION> & hori_matrix_,
        TValue position_)
{
    SEQAN_CHECKPOINT;
	typedef typename Position<Matrix<TScoreValue, DIMENSION> >::Type TPosition;

	typedef Align<TTargetSource, TTargetSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TTargetIterator;
    typedef typename Infix<TTargetSource>::Type TTargetSourceSegment;
	typedef typename Iterator<TTargetSource, Standard>::Type TStringIterator;

	TTargetSourceSegment str_0 = sourceSegment(row(target_, 0));
	TTargetSourceSegment str_1 = sourceSegment(row(target_, 1));
	typename Size<TTargetSourceSegment>::Type dim_0_len = length(hori_matrix_,0);

	TPosition pos = position_;

	TTargetIterator target_0 = iter(row(target_, 0), 0, Standard());
	TTargetIterator target_1 = iter(row(target_, 1), 0, Standard());

	TStringIterator it_0 = iter(str_0, 0, Standard());
	TStringIterator it_0_end = end(str_0);

	TStringIterator it_1 = iter(str_1, 0, Standard());
	TStringIterator it_1_end = end(str_1);
	
	//-------------------------------------------------------------------------
	//follow the trace until the border is reached
	while ((it_0 != it_0_end) && (it_1 != it_1_end))
	{
		if (getValue(diag_matrix_,pos) > getValue(hori_matrix_,pos))
		{
			if (getValue(diag_matrix_,pos) > getValue(vert_matrix_,pos))
			{
				++it_0;
				++it_1;
				pos += dim_0_len;
			} else{
				++it_1;
				insertGap(target_0);
				pos += dim_0_len-1;
			}
		}
		else
		{
			if (getValue(hori_matrix_,pos) > getValue(vert_matrix_,pos)) 
			{					
				++it_0;
				insertGap(target_1);
				++pos;
			} 
			else
			{
				++it_1;
				insertGap(target_0);
				pos += dim_0_len-1;
			}
		}
		++target_0;
		++target_1;
	}
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_SEED_BANDED_AFFINE_H_
