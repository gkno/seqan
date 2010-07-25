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
  Banded chain alignment around a chain of seeds; Case with affine gap costs.
  Based on the original code by Carsten Kemena, adapted to the new seeds
  interface.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_ALIGN_CHAIN_BANDED_AFFINE_H_
#define SEQAN_SEEDS_ALIGN_CHAIN_BANDED_AFFINE_H_

namespace seqan {

// TODO(bkehr): large parts of this file are identical to banded_chain_align.h

// ===========================================================================
// Forwards
// ===========================================================================

template <typename TScoreValue, unsigned DIMENSION, typename TSeedSpec, typename TSeedConfig, typename TString, typename TValue2>
void
_gotoh_rectangle(Matrix<TScoreValue, DIMENSION> & matrix_diag,
				 Matrix<TScoreValue, DIMENSION> & matrix_vert,
				 Matrix<TScoreValue, DIMENSION> & matrix_hori,	//edit matrix
				 Seed<TSeedSpec, TSeedConfig> const &seed1,				//Seed nearer to the end
				 Seed<TSeedSpec, TSeedConfig> const &seed2,				//Seed nearer to the start
				 TValue2 k_begin,									//upper diagonal extension
				 TValue2 k_end,										//lower diagonal extension
				 TString const & str1_,								//first sequence
				 TString const & str2_,								//secondSequence
				 Score<TScoreValue, Simple> const & score_,			//score matrix
				 String<TScoreValue> & init_diag,					//Values for initialisation
				 String<TScoreValue> & init_vert,					//Values for initialisation
				 String<TScoreValue> & init_hori);					//Values for initialisation

template <typename TContainer, typename TValue, typename TScore, typename TAlign>
TScore
chain_to_alignment_gotoh(TContainer const &seedChain, 
						TValue k,
						TAlign & whole_alignment, 
						Score< TScore, Simple> const &scoreMatrix);

template <typename TScoreValue, unsigned DIMENSION, typename TSeedSpec, typename TSeedConfig, typename TString, typename TValue2>
TScoreValue
_banded_gotoh(Matrix<TScoreValue, DIMENSION> & matrix_diag,
			  Matrix<TScoreValue, DIMENSION> & matrix_vert,
			  Matrix<TScoreValue, DIMENSION> & matrix_hori,
			  Seed<TSeedSpec, TSeedConfig> const &seed,
			  TValue2 k,
			  TString const & str1_,
			  TString const & str2_,
			  Score<TScoreValue, Simple> const & score_,
			  String<TScoreValue> & init_diag,
			  String<TScoreValue> & init_vert,
			  String<TScoreValue> & init_hori);

template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, unsigned DIMENSION>
typename Size<Matrix<TScoreValue, DIMENSION> >::Type
_banded_gotoh_trace2(Align<TTargetSource, TTargetSpec> & target_,
					 Matrix<TScoreValue, DIMENSION> & diag_matrix_,
					 Matrix<TScoreValue, DIMENSION> & vert_matrix_,
					 Matrix<TScoreValue, DIMENSION> & hori_matrix_,
					 Iter< Matrix<TScoreValue, DIMENSION>, PositionIterator > source_
					);

template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, unsigned DIMENSION>
TScoreValue
_gotoh_trace_lastRectangle(Align<TTargetSource, TTargetSpec> & target_,
						Matrix<TScoreValue, DIMENSION> & diag_matrix_,
						Matrix<TScoreValue, DIMENSION> & vert_matrix_,
						Matrix<TScoreValue, DIMENSION> & hori_matrix_,
						Iter< Matrix<TScoreValue, DIMENSION>, PositionIterator > source_);

template<typename TSeed, typename TString, typename TDiff, typename TMatrix, typename TScoreString, typename TValue, typename TAlign, typename TScoreMatrix>
void
_calculateBandedSeedGotoh(TSeed const &seed,
					 TDiff k,
					 TMatrix &matrix_diag,
					 TMatrix &matrix_vert,
					 TMatrix &matrix_hori,
					 TString *p_seq1,
					 TString *p_seq2,
					 TScoreString &score_str_diag,
					 TScoreString &score_str_vert,
					 TScoreString &score_str_hori,
					 TValue &score_length,
					 ::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > &alignmentVector,
					 TScoreMatrix const &scoreMatrix);

template<typename TSeed, typename TString, typename TDiff, typename TMatrix, typename TScoreString, typename TValue, typename TAlign, typename TScoreMatrix>
void
_calculateFirstRectangleGotoh(TSeed const &seed,
						 TDiff k,
						 TMatrix &matrix_diag,
						 TMatrix &matrix_vert,
						 TMatrix &matrix_hori,
						 TString *p_seq1,
						 TString *p_seq2,
						 TScoreString &score_str_diag,
						 TScoreString &score_str_vert,
						 TScoreString &score_str_hori,
						 TValue &score_length,
						 ::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > &alignmentVector,
						 TScoreMatrix const &scoreMatrix);
template<typename TSeed, typename TString, typename TDiff, typename TMatrix, typename TScoreString, typename TValue, typename TAlign, typename TScoreMatrix>
void
_calculateLastRectangleGotoh(TSeed const &seed,
						TDiff k,
						TMatrix & matrix_diag,
						TMatrix & matrix_vert,
						TMatrix & matrix_hori,
						TString *p_seq1,
						TString *p_seq2,
						TScoreString & score_str_diag,
						TScoreString & score_str_vert,
						TScoreString & score_str_hori,
						TValue &score_length,
						::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > &alignmentVector,
						TScoreMatrix const &scoreMatrix);

template<typename TSeed, typename TString, typename TDiff, typename TMatrix, typename TScoreString, typename TValue, typename TAlign, typename TScoreMatrix>
void
_calculateRectangleGotoh(TSeed const &seed,
					TSeed const &seed2,
					TDiff k_begin,
					TDiff k_end,
					TMatrix & matrix_diag,
					TMatrix & matrix_vert,
					TMatrix & matrix_hori,
					TString *p_seq1,
					TString *p_seq2,
					TScoreString &score_str_diag,
					TScoreString &score_str_vert,
					TScoreString &score_str_hori,
					TValue &score_length,
					::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > & alignmentVector,
					TScoreMatrix const &scoreMatrix);

template <typename TScoreValue, unsigned DIMENSION, typename TSeedSpec, typename TSeedConfig, typename TString, typename TValue2>
void
_banded_gotoh_rectangle_first(Matrix<TScoreValue, DIMENSION> & matrix_diag,	//edit matrix
							  Matrix<TScoreValue, DIMENSION> & matrix_vert,
							  Matrix<TScoreValue, DIMENSION> & matrix_hori,
							  Seed<TSeedSpec, TSeedConfig> const &seed,				//Seed
							  TValue2 k,											//diagonal extension
							  TString const & str1_,							//first sequence
							  TString const & str2_,							//secondSequence
							  Score<TScoreValue, Simple> const & score_,		//score matrix
							  String<TScoreValue> & init_diag,					//Values for initialisation
							  String<TScoreValue> & init_vert,
							  String<TScoreValue> & init_hori);

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

template <typename TContainer, typename TValue, typename TScore, typename TAlign>
TScore
chain_to_alignment_gotoh(TContainer const &seedChain, 
						TValue k,
						TAlign & whole_alignment, 
						Score< TScore, Simple> const &scoreMatrix)
{
    typedef typename Infix<typename Source<TAlign>::Type>::Type TString;    //Sequence in an align object
	typedef typename Size<TString>::Type TSize;								//Size of the string
	typedef typename Iterator<const TContainer, Standard>::Type TIterator;	//Iterator for the seed chain
	typedef typename Value<TContainer>::Type TSeed;							//Type of Seed
	typedef typename Value<TSeed>::Type TPosition;
	typedef String<TScore> TScoreString;
	typedef Matrix<TScore> TMatrix;
	typedef Iter<TMatrix, PositionIterator > TMatrixIterator;
	typedef ::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > TAlignVector;
	typedef typename ::std::map<TValue, Pair<TValue, TAlign> >::iterator TMapIterator;
	
	TScoreString score_str_diag;
	TScoreString score_str_vert;
	TScoreString score_str_hori;
	TAlignVector alignmentVector;
	TString seq1 = sourceSegment(row(whole_alignment, 0));
	TString seq2 = sourceSegment(row(whole_alignment, 1));
	TString *p_seq1 = &seq1;
	TString *p_seq2 = &seq2;
	
	TValue score_length = 0;
	TMatrix matrix_diag;
	TMatrix matrix_vert;
	TMatrix matrix_hori;
	setDimension(matrix_vert, 2);
	setDimension(matrix_hori, 2);
	setDimension(matrix_diag, 2);
	TIterator it = begin(seedChain);

	//calculation begins at the end
	//rectangle between last seed and end
	TValue k_begin = k;
	TValue k_end = k;
	if ((length(*it) <= k) ||(((getEndDim1(*it) - 1)-getBeginDim1(*it)) <=k))
 		k_end = length(*it) -1;

	
	_calculateLastRectangleGotoh(*it, k_end, matrix_diag, matrix_vert, matrix_hori, p_seq1, p_seq2, score_str_diag, score_str_vert, score_str_hori, score_length, alignmentVector, scoreMatrix);
	/*
	cout <<"INIT: " << endl;
	cout << length(score_str_diag) << endl;
	for ( int i = 0; i < length(score_str_diag); ++i)
		cout << score_str_diag[i] << " ";
	cout << endl;*/


	_calculateBandedSeedGotoh(*it, k_end, matrix_diag, matrix_vert, matrix_hori, p_seq1, p_seq2, score_str_diag, score_str_vert, score_str_hori, score_length, alignmentVector, scoreMatrix);
	
	TIterator it_begin = end(seedChain);
	TIterator it2 = it;
	++it2;
	while (it2 != it_begin)
	{
		k_begin = k_end;
		if ((length(*it2) <= k) || (((getEndDim1(*it2) - 1)-getBeginDim1(*it2)) <=k))
			k_end = length(*it2) -2;
		else 
			k_end = k;
		_calculateRectangleGotoh(*it, *it2, k_begin, k_end, matrix_diag, matrix_vert, matrix_hori, p_seq1, p_seq2, score_str_diag, score_str_vert, score_str_hori, score_length, alignmentVector, scoreMatrix);	
		
		_calculateBandedSeedGotoh(*it2, k_end, matrix_diag, matrix_vert, matrix_hori, p_seq1, p_seq2, score_str_diag, score_str_vert, score_str_hori, score_length, alignmentVector, scoreMatrix);
		
		++it;
		++it2;
	}
	//cout << "last" << endl;
	_calculateFirstRectangleGotoh(*it, k_end, matrix_diag, matrix_vert, matrix_hori, p_seq1, p_seq2, score_str_diag, score_str_vert, score_str_hori, score_length, alignmentVector, scoreMatrix);
	//cout << "last" << endl;
	_constructAlignment(alignmentVector, whole_alignment);
	//cout << "last" << endl;
	return score_str_diag[0];
}

template <typename TScoreValue, unsigned DIMENSION, typename TSeedSpec, typename TSeedConfig, typename TString, typename TValue2>
TScoreValue
_banded_gotoh(Matrix<TScoreValue, DIMENSION> & matrix_diag,
			  Matrix<TScoreValue, DIMENSION> & matrix_vert,
			  Matrix<TScoreValue, DIMENSION> & matrix_hori,
			  Seed<TSeedSpec, TSeedConfig> const &seed,
			  TValue2 k,
			  TString const & str1_,
			  TString const & str2_,
			  Score<TScoreValue, Simple> const & score_,
			  String<TScoreValue> & init_diag,
			  String<TScoreValue> & init_vert,
			  String<TScoreValue> & init_hori)
{
	typedef Matrix<TScoreValue, DIMENSION> TMatrix;

	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Iterator<TMatrix, Rooted>::Type TMatrixIterator;

	typedef typename Iterator<TString const, Rooted>::Type TStringIterator;
	typedef typename Value<TString const>::Type TAlphabet;

	//-------------------------------------------------------------------------
	//define some variables

	TSize height = getLowerDiagonal(seed) - getUpperDiagonal(seed)+1;
	TSize str1_length = height+2*k;
	TSize str2_length = length(str2_);

	TSize up_height = getLowerDiagonal(seed) - getStartDiagonal(seed) + k; //equals length of empty lower triangle
	//TValue up_width = str1_length - up_height;
	TSize down_height = getEndDiagonal(seed) - getUpperDiagonal(seed) + k; //equals length of empty lower triangle
	TSize down_width = str1_length - down_height;
	TSize length_right_diag = str2_length - down_height;
	TSize length_left_diag = str2_length - up_height;


	TStringIterator x_begin = begin(str1_) - 1;
	TStringIterator x_end = end(str1_) - 1;
	TStringIterator y_begin = begin(str2_) - 1;
	TStringIterator y_end = end(str2_) - 1;

	TStringIterator x = x_end;
	TStringIterator y;

	TScoreValue score_match = scoreMatch(score_);
	TScoreValue score_mismatch = scoreMismatch(score_);
	TScoreValue score_gap_open = scoreGapOpen(score_);
	TScoreValue score_gap_extend = scoreGapExtend(score_);

	setLength(matrix_diag, 0, str1_length + 2);
	setLength(matrix_diag, 1, str2_length + 1);
	resize(matrix_diag);
	setLength(matrix_vert, 0, str1_length + 2);
	setLength(matrix_vert, 1, str2_length + 1);
	resize(matrix_vert);
	setLength(matrix_hori, 0, str1_length + 2);
	setLength(matrix_hori, 1, str2_length + 1);
	resize(matrix_hori);

	TSize pos = length(matrix_diag, 0)-1;
	TMatrixIterator diag_col_ = begin(matrix_diag);
	TMatrixIterator diag_finger1(matrix_diag,pos);
	TMatrixIterator diag_finger2;
	TMatrixIterator diag_finger3;
	TMatrixIterator vert_col_ = begin(matrix_vert);
	TMatrixIterator vert_finger1(matrix_vert,pos);
	TMatrixIterator vert_finger2;
	TMatrixIterator vert_finger3;
	TMatrixIterator hori_col_ = begin(matrix_hori);
	TMatrixIterator hori_finger1(matrix_hori,pos);
	TMatrixIterator hori_finger2;
	TMatrixIterator hori_finger3;
	

	//-------------------------------------------------------------------------
	// init
	

	TScoreValue inf = -1000000;
	for (TSize i = 1; i != length_right_diag; ++i){
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

	int z = -1;

	for (int i = -1; i != down_height; ++i){
		goPrevious(diag_finger1, 0);
		goNext(diag_finger1,1);
		goPrevious(vert_finger1, 0);
		goNext(vert_finger1,1);
		goPrevious(hori_finger1, 0);
		goNext(hori_finger1,1);
		*diag_finger1 = init_diag[++z];
		*vert_finger1 = init_vert[z];
		*hori_finger1 = init_hori[z];
	}

	//*diag_finger1 = init_diag[++z];
	//*vert_finger1 = init_vert[z];
	//*hori_finger1 = init_hori[z];
	
	diag_finger3 = diag_finger1;
	hori_finger3 = hori_finger1;
	vert_finger3 = vert_finger1;

	while (z < static_cast<int>(length(init_diag))-1)
	{
		goPrevious(diag_finger1, 0);
		goPrevious(hori_finger1, 0);
		goPrevious(vert_finger1, 0);
		*diag_finger1 = init_diag[++z];
		*hori_finger1 = init_hori[z];
		*vert_finger1 = init_vert[z];
	}

	TSize arg = length(matrix_diag,0) - 2;
	while ( z < arg)
	{
	goPrevious(diag_finger1, 0);
	goPrevious(hori_finger1, 0);
	goPrevious(vert_finger1, 0);
	*diag_finger1 = inf;
	*vert_finger1 = inf;
	*hori_finger1 = inf;
	++z;
	}


	//printMatrix(matrix_diag);
	//printMatrix(matrix_vert);
	//printMatrix(matrix_hori);

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
			hori_value = (hori_value+score_gap_extend > diag_back+score_gap_open) ? hori_value+score_gap_extend : diag_back+score_gap_open;
			*hori_finger1 = hori_value;
			goPrevious(hori_finger1);

			tmp_vert = *vert_finger2;
			vert_value = (tmp_vert+score_gap_extend > diag_front+score_gap_open) ? tmp_vert+score_gap_extend : diag_front+score_gap_open;
			*vert_finger1 = vert_value;
			goPrevious(vert_finger1);
			goPrevious(vert_finger2);

			
			tmp_diag = (vert_value > hori_value) ? vert_value  : hori_value;
			diag_value = diag_match + ((*x == *y) ? score_match : score_mismatch);
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
		if (measure < length_left_diag)
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
			hori_value = (hori_value+score_gap_extend > diag_back+score_gap_open) ? hori_value+score_gap_extend : diag_back+score_gap_open;
			*hori_finger1 = hori_value;
			goPrevious(hori_finger1);

			tmp_vert = *vert_finger2;
			vert_value = (tmp_vert+score_gap_extend > diag_front+score_gap_open) ? tmp_vert+score_gap_extend : diag_front+score_gap_open;
			*vert_finger1 = vert_value;
			goPrevious(vert_finger1);
			goPrevious(vert_finger2);

			tmp_diag = (vert_value > hori_value) ? vert_value  : hori_value;
			diag_value = diag_match + ((*x == *y) ? score_match : score_mismatch);
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
		if (measure >= length_left_diag)
			--run_length;
		--y;
	}
	--run_length;
	++diag_finger1;
	//cout << "ENDE" << endl << endl << endl;
	return *diag_finger1;
}


//////////////////////////////////////////////////////////////////////////////
//traceback through needleman wunsch matrix


//Position berechnen!
template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, unsigned DIMENSION>
typename Size<Matrix<TScoreValue, DIMENSION> >::Type
_banded_gotoh_trace2(Align<TTargetSource, TTargetSpec> & target_,
					 Matrix<TScoreValue, DIMENSION> & diag_matrix_,
					 Matrix<TScoreValue, DIMENSION> & vert_matrix_,
					 Matrix<TScoreValue, DIMENSION> & hori_matrix_,
					 Iter< Matrix<TScoreValue, DIMENSION>, PositionIterator > source_
					)
{
	SEQAN_CHECKPOINT
	
	typedef typename Position<Matrix<TScoreValue, DIMENSION> >::Type TPosition;

	typedef Align<TTargetSource, TTargetSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TTargetIterator;
    typedef typename Infix<TTargetSource>::Type TTargetSourceSegment;
	typedef typename Iterator<TTargetSourceSegment, Rooted>::Type TStringIterator;

	TTargetSourceSegment str_0 = sourceSegment(row(target_, 0));
	TTargetSourceSegment str_1 = sourceSegment(row(target_, 1));
	typename Size<TTargetSourceSegment>::Type dim_0_len = length(diag_matrix_,0);

	TPosition pos = position(source_);

	TTargetIterator target_0 = iter(row(target_, 0), 0, Standard());
	TTargetIterator target_1 = iter(row(target_, 1), 0, Standard());

	TStringIterator it_0 = iter(str_0, 0);
	TStringIterator it_0_end = end(str_0);

	TStringIterator it_1 = iter(str_1, 0);
	TStringIterator it_1_end = end(str_1);
	
	//-------------------------------------------------------------------------
	//follow the trace until the border is reached

	while ((it_0 != it_0_end) && (it_1 != it_1_end))
	{
		if (getValue(diag_matrix_,pos) == getValue(vert_matrix_,pos))
		{
			++it_1;
			insertGap(target_0);
			pos += dim_0_len-1;
		}
		else
		{
			if (getValue(hori_matrix_,pos) == getValue(diag_matrix_,pos)) 
			{					
				++it_0;
				insertGap(target_1);
				++pos;
			} 
			else
			{
				++it_0;
				++it_1;
				pos += dim_0_len;
			}
		}
		++target_0;
		++target_1;
	}

	setSourceEndPosition(row(target_,1),position(it_1));
	setSourceEndPosition(row(target_,0),position(it_0));
	
	return length(diag_matrix_,0) - (pos % dim_0_len)-2;
}


//changed version of usual trace-back
template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, unsigned DIMENSION>
TScoreValue
_gotoh_trace_lastRectangle(Align<TTargetSource, TTargetSpec> & target_,
						Matrix<TScoreValue, DIMENSION> & diag_matrix_,
						Matrix<TScoreValue, DIMENSION> & vert_matrix_,
						Matrix<TScoreValue, DIMENSION> & hori_matrix_,
						Iter< Matrix<TScoreValue, DIMENSION>, PositionIterator > source_)
{
	SEQAN_CHECKPOINT
	typedef typename Position<Matrix<TScoreValue, DIMENSION> >::Type TPosition;

	typedef Align<TTargetSource, TTargetSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TTargetIterator;
    typedef typename Infix<TTargetSource>::Type TTargetSourceSegment;
	typedef typename Iterator<TTargetSourceSegment, Rooted>::Type TStringIterator;

	TTargetSourceSegment str_0 = sourceSegment(row(target_, 0));
	TTargetSourceSegment str_1 = sourceSegment(row(target_, 1));
	typename Size<TTargetSourceSegment>::Type dim_0_len = length(hori_matrix_,0);

	TPosition pos = position(source_);

	TTargetIterator target_0 = iter(row(target_, 0), 0, Standard());
	TTargetIterator target_1 = iter(row(target_, 1), 0, Standard());

	TStringIterator it_0 = iter(str_0, 0);
	TStringIterator it_0_end = end(str_0);

	TStringIterator it_1 = iter(str_1, 0);
	TStringIterator it_1_end = end(str_1);
	
	//-------------------------------------------------------------------------
	//follow the trace until the border is reached
	while ((it_0 != it_0_end) && (it_1 != it_1_end))
	{
		if (getValue(diag_matrix_,pos) == getValue(hori_matrix_,pos))
		{
			++it_0;
			insertGap(target_1);
			++pos;
		} 
		else
		{
			if (getValue(diag_matrix_,pos) == getValue(vert_matrix_,pos))
			{
				++it_1;
				insertGap(target_0);
				pos += dim_0_len;
			} else{
				++it_0;
				++it_1;
				pos += dim_0_len+1;
			}
		}
		++target_0;
		++target_1;
	}
	return 0;
}

//calculation and backtracking of the banded alignment of a seed
template<typename TSeed, typename TString, typename TDiff, typename TMatrix, typename TScoreString, typename TValue, typename TAlign, typename TScoreMatrix>
void
_calculateBandedSeedGotoh(TSeed const &seed,
					 TDiff k,
					 TMatrix &matrix_diag,
					 TMatrix &matrix_vert,
					 TMatrix &matrix_hori,
					 TString *p_seq1,
					 TString *p_seq2,
					 TScoreString &score_str_diag,
					 TScoreString &score_str_vert,
					 TScoreString &score_str_hori,
					 TValue &score_length,
					 ::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > &alignmentVector,
					 TScoreMatrix const &scoreMatrix)
{
	typedef typename ::std::map<TValue,Pair<TValue, TAlign> >::iterator TMapIterator;
	typedef Iter<TMatrix, PositionIterator> TMatrixIterator;
    typedef typename Infix<TString>::Type TSegment;
	TSegment seg1_align = infix(host(*p_seq1), getBeginDim0(seed), (getEndDim0(seed) - 1)+1);
	TSegment seg2_align = infix(host(*p_seq2), getBeginDim1(seed), (getEndDim1(seed) - 1)+1);
	_banded_gotoh(matrix_diag, matrix_vert, matrix_hori, seed, k, seg1_align, seg2_align, scoreMatrix, score_str_diag, score_str_vert, score_str_hori);

	//printMatrix(matrix_diag);
	TValue height_diag = getLowerDiagonal(seed)-getStartDiagonal(seed)+k;
	TValue width_diag = getStartDiagonal(seed)-getUpperDiagonal(seed)+k;
	TValue overall = height_diag + width_diag + 1;
	
	resize(score_str_diag,overall);
	resize(score_str_vert,overall);
	resize(score_str_hori,overall);

	TMatrixIterator matr_it = begin(matrix_diag);
	setPosition(matr_it, length(matrix_diag,0)-2);

	alignmentVector.push_back(::std::map<TValue,Pair<TValue, TAlign> >());
	TValue width_align = getBeginDim0(seed) + width_diag;
	TValue height_align = getBeginDim1(seed);

	//diagonal back_track
	TValue new_connect;
	TValue old_connect = -1;

	for(TValue j = 0; j<=width_diag; ++j)
	{
		alignmentVector.back().insert(std::make_pair(j,Pair<TValue, TAlign> ()));
		TMapIterator mapIt = --(alignmentVector.back().end());
		TSegment seg1_align = infix(host(*p_seq1), width_align, (getEndDim0(seed) - 1)+1);
		TSegment seg2_align = infix(host(*p_seq2), height_align,(getEndDim1(seed) - 1)+1);
		resize(rows(mapIt->second.i2),2);
		assignSource(row(mapIt->second.i2,0), seg1_align);
		assignSource(row(mapIt->second.i2,1), seg2_align);
		new_connect = _banded_gotoh_trace2(mapIt->second.i2, matrix_diag, matrix_vert, matrix_hori, matr_it);
		score_str_diag[j] = *matr_it;
		score_str_vert[j] = getValue(matrix_vert, position(matr_it));
		score_str_hori[j] = getValue(matrix_hori, position(matr_it));
		mapIt->second.i1 = new_connect;
		goPrevious(matr_it,0);
		--width_align;
		if (old_connect != new_connect)
			_deleteAlignment(alignmentVector, old_connect, new_connect);
		old_connect = new_connect;
	}
			
	++width_align;
	++height_align;

	for(TValue j = width_diag+1; j < overall; ++j)
	{
		goNext(matr_it,1);
		alignmentVector.back().insert(std::make_pair(j,Pair<TValue, TAlign> ()));
		TMapIterator mapIt = --(alignmentVector.back().end());
		TSegment seg1_align = infix(host(*p_seq1), width_align, (getEndDim0(seed) - 1)+1);
		TSegment seg2_align = infix(host(*p_seq2), height_align,(getEndDim1(seed) - 1)+1);
		resize(rows(mapIt->second.i2),2);
		assignSource(row(mapIt->second.i2,0), seg1_align);
		assignSource(row(mapIt->second.i2,1), seg2_align);
		
		new_connect = _banded_gotoh_trace2(mapIt->second.i2, matrix_diag, matrix_vert, matrix_hori,  matr_it);
		
		score_str_diag[j] = *matr_it;
		score_str_vert[j] = getValue(matrix_vert, position(matr_it));
		score_str_hori[j] = getValue(matrix_hori, position(matr_it));


		mapIt->second.i1 = new_connect;
		goPrevious(matr_it,0);
		++height_align;
		if (old_connect != new_connect)
			_deleteAlignment(alignmentVector, old_connect, new_connect);
		old_connect = new_connect;
	}
	_deleteAlignment(alignmentVector, old_connect, score_length);
	score_length = length(score_str_diag);
}

template<typename TSeed, typename TString, typename TDiff, typename TMatrix, typename TScoreString, typename TValue, typename TAlign, typename TScoreMatrix>
void
_calculateFirstRectangleGotoh(TSeed const &seed,
						 TDiff k,
						 TMatrix &matrix_diag,
						 TMatrix &matrix_vert,
						 TMatrix &matrix_hori,
						 TString *p_seq1,
						 TString *p_seq2,
						 TScoreString &score_str_diag,
						 TScoreString &score_str_vert,
						 TScoreString &score_str_hori,
						 TValue &score_length,
						 ::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > &alignmentVector,
						 TScoreMatrix const &scoreMatrix)
{
	typedef typename ::std::map<TValue,Pair<TValue, TAlign> >::iterator TMapIterator;
	typedef Iter<TMatrix, PositionIterator> TMatrixIterator;
	TValue new_connect;
    typedef typename Infix<TString>::Type TSegment;
	TSegment seg1b_align = infix(host(*p_seq1), beginPosition(*p_seq1), getBeginDim0(seed) + getStartDiagonal(seed) - getUpperDiagonal(seed) + k);
	TSegment seg2b_align = infix(host(*p_seq2), beginPosition(*p_seq2), getBeginDim1(seed) + getLowerDiagonal(seed)  - getStartDiagonal(seed) + k);

	_banded_gotoh_rectangle_first(matrix_diag, matrix_vert, matrix_hori, seed, k, seg1b_align, seg2b_align, scoreMatrix, score_str_diag, score_str_vert, score_str_hori);

	TValue w_d2 = getStartDiagonal(seed) - getUpperDiagonal(seed) + k;
	TValue h_d2 = getLowerDiagonal(seed) - getStartDiagonal(seed) + k;

	resize(score_str_diag,1);
	resize(score_str_vert,1);
	resize(score_str_hori,1);

	TMatrixIterator matr_it = begin(matrix_diag);

	alignmentVector.push_back(::std::map<TValue,Pair<TValue, TAlign> >());

	TValue width_stop = getBeginDim0(seed) - beginPosition(*p_seq1);
	TValue height_stop = getBeginDim1(seed) - beginPosition(*p_seq2);

	alignmentVector.back().insert(std::make_pair(0,Pair<TValue, TAlign> ()));
	TMapIterator mapIt = --(alignmentVector.back().end());
	TSegment seg1_align = infix(host(*p_seq1), beginPosition(*p_seq1), getBeginDim0(seed) + w_d2);
	TSegment seg2_align = infix(host(*p_seq2), beginPosition(*p_seq2), getBeginDim1(seed)+ h_d2);
	
	resize(rows(mapIt->second.i2),2);
	assignSource(row(mapIt->second.i2,0), seg1_align);
	assignSource(row(mapIt->second.i2,1), seg2_align);

	new_connect = _gotoh_trace_rectangle(mapIt->second.i2,  matr_it, matrix_diag, matrix_vert, matrix_hori, width_stop, height_stop);

	score_str_diag[0] = *matr_it;
	mapIt->second.i1 = new_connect;
	_deleteAlignment(alignmentVector, -1, new_connect);
	_deleteAlignment(alignmentVector, new_connect, score_length);
}

template<typename TSeed, typename TString, typename TDiff, typename TMatrix, typename TScoreString, typename TValue, typename TAlign, typename TScoreMatrix>
void
_calculateLastRectangleGotoh(TSeed const &seed,
						TDiff k,
						TMatrix & matrix_diag,
						TMatrix & matrix_vert,
						TMatrix & matrix_hori,
						TString *p_seq1,
						TString *p_seq2,
						TScoreString & score_str_diag,
						TScoreString & score_str_vert,
						TScoreString & score_str_hori,
						TValue &score_length,
						::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > &alignmentVector,
						TScoreMatrix const &scoreMatrix)
{
	typedef typename ::std::map<TValue,Pair<TValue, TAlign> >::iterator TMapIterator;
	typedef Iter<TMatrix, PositionIterator> TMatrixIterator;
    typedef typename Infix<TString>::Type TSegment;

	TValue seq1_end = endPosition(*p_seq1);
	TValue seq2_end = endPosition(*p_seq2);
	TValue width_diag=getLowerDiagonal(seed) - getEndDiagonal(seed)+k;
	TValue width =  width_diag + seq1_end - (getEndDim0(seed) - 1);
	TValue height_diag = getEndDiagonal(seed) - getUpperDiagonal(seed)+k;
	TValue height = height_diag + seq2_end - (getEndDim1(seed) - 1);
	
	resize(score_str_diag, width_diag+height_diag+1);
	resize(score_str_vert, width_diag+height_diag+1);
	resize(score_str_hori, width_diag+height_diag+1);
	
	TSegment seg1 = infix(host(*p_seq1),seq1_end-width+1,seq1_end);
	TSegment seg2 = infix(host(*p_seq2),seq2_end-height+1,seq2_end);

	//_gotoh(matrix_diag, matrix_vert, matrix_hori, seg1, seg2, scoreMatrix);
	_gotoh2(matrix_diag, matrix_vert, matrix_hori, seg1, seg2, scoreMatrix);



	TValue width_align = seq1_end - width  + width_diag +1;
	TValue height_align = seq2_end - height+1;

	TMatrixIterator iter_ = begin(matrix_diag);

	TValue x = width_diag;
	alignmentVector.push_back(::std::map<TValue,Pair<TValue, TAlign> >());
	
	//last rectangle
	for(TValue i = 0; i<height_diag; ++i)
	{
		alignmentVector[0].insert(std::make_pair(i,Pair<TValue, TAlign> ()));
		TMapIterator mapIt = --alignmentVector[0].end();
		TSegment seg1_align = infix(host(*p_seq1), width_align, seq1_end);
		TSegment seg2_align = infix(host(*p_seq2), height_align, seq2_end);
		resize(rows(mapIt->second.i2),2);
		assignSource(row(mapIt->second.i2,0),seg1_align);
		assignSource(row(mapIt->second.i2,1), seg2_align);
		setPosition(iter_, x);
		score_str_diag[i] = *iter_;
		score_str_vert[i] = getValue(matrix_vert,position(iter_));
		score_str_hori[i] = getValue(matrix_hori,position(iter_));

		mapIt->second.i1 = -1;
		_gotoh_trace_lastRectangle(mapIt->second.i2, matrix_diag, matrix_vert, matrix_hori, iter_);
		x+=width;
		++height_align;
	}
	TValue a = height_diag + width_diag+1;
	
	for(TValue i = height_diag; i<a; ++i)
	{
		alignmentVector[0].insert(std::make_pair(i,Pair<TValue, TAlign> ()));
		TMapIterator mapIt = --alignmentVector[0].end();
		TSegment seg1_align = infix(host(*p_seq1), width_align, seq1_end);
	    TSegment seg2_align = infix(host(*p_seq2), height_align, seq2_end);
		resize(rows(mapIt->second.i2),2);
		assignSource(row(mapIt->second.i2,0),seg1_align);
		assignSource(row(mapIt->second.i2,1), seg2_align);
		setPosition(iter_, x);
		score_str_diag[i] = *iter_;
		score_str_vert[i] = getValue(matrix_vert,position(iter_));
		score_str_hori[i] = getValue(matrix_hori,position(iter_));
		mapIt->second.i1 = -1;
		_gotoh_trace_lastRectangle(mapIt->second.i2, matrix_diag, matrix_vert, matrix_hori, iter_);
		--x;
		--width_align;
	}
	score_length = length(score_str_diag);
}

//calculation and backtracking of the edit matrix between two seeds
template<typename TSeed, typename TString, typename TDiff, typename TMatrix, typename TScoreString, typename TValue, typename TAlign, typename TScoreMatrix>
void
_calculateRectangleGotoh(TSeed const &seed,
					TSeed const &seed2,
					TDiff k_begin,
					TDiff k_end,
					TMatrix & matrix_diag,
					TMatrix & matrix_vert,
					TMatrix & matrix_hori,
					TString *p_seq1,
					TString *p_seq2,
					TScoreString &score_str_diag,
					TScoreString &score_str_vert,
					TScoreString &score_str_hori,
					TValue &score_length,
					::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > & alignmentVector,
					TScoreMatrix const &scoreMatrix)
{
	typedef typename ::std::map<TValue,Pair<TValue, TAlign> >::iterator TMapIterator;
	typedef Iter<TMatrix, PositionIterator> TMatrixIterator;
    typedef typename Infix<TString>::Type TSegment;
	TSegment seg1b_align = infix(host(*p_seq1), (getEndDim0(seed2) - 1)-(getLowerDiagonal(seed2) - getEndDiagonal(seed2)   + k_end) + 1, getBeginDim0(seed) + getStartDiagonal(seed) - getUpperDiagonal(seed) + k_begin);
	TSegment seg2b_align = infix(host(*p_seq2), (getEndDim1(seed2) - 1)-(getEndDiagonal(seed2) -  getUpperDiagonal(seed2) + k_end) + 1, getBeginDim1(seed) + getLowerDiagonal(seed)  - getStartDiagonal(seed) + k_begin);

	_gotoh_rectangle(matrix_diag, matrix_vert, matrix_hori, seed, seed2, k_begin, k_end, seg1b_align, seg2b_align, scoreMatrix, score_str_diag, score_str_vert, score_str_hori);


//	printMatrix(matrix_diag);
//	cout << endl;
	TValue width_diag = getLowerDiagonal(seed2) - getEndDiagonal(seed2)+k_end;
	TValue height_diag = getEndDiagonal(seed2) - getUpperDiagonal(seed2)+k_end;
	
	TValue w_d2 = getStartDiagonal(seed) - getUpperDiagonal(seed) + k_begin;
	TValue h_d2 = getLowerDiagonal(seed) - getStartDiagonal(seed) + k_begin;

	TValue overall = height_diag + width_diag + 1;

	resize(score_str_diag, overall);
	resize(score_str_vert, overall);
	resize(score_str_hori, overall);

	TMatrixIterator matr_it = begin(matrix_diag);
	setPosition(matr_it, width_diag);

	alignmentVector.push_back(::std::map<TValue,Pair<TValue, TAlign> >());
	TValue width_align = (getEndDim0(seed2) - 1)+1;
	TValue height_align = (getEndDim1(seed2) - 1) - height_diag + 1;

	TValue width_stop = getBeginDim0(seed) - (getEndDim0(seed2) - 1)-1;
	TValue height_stop = getBeginDim1(seed) - ((getEndDim1(seed2) - 1) - height_diag)-1;

	TValue old_connect = -1;
	TValue new_connect;
	for(TValue j = 0; j<height_diag; ++j)
	{
		alignmentVector.back().insert(std::make_pair(j,Pair<TValue, TAlign> ()));
		TMapIterator mapIt = --(alignmentVector.back().end());
		TSegment seg1_align = infix(host(*p_seq1), width_align, getBeginDim0(seed) + w_d2);
		TSegment seg2_align = infix(host(*p_seq2), height_align, getBeginDim1(seed)+ h_d2);

		resize(rows(mapIt->second.i2),2);
		assignSource(row(mapIt->second.i2,0), seg1_align);
		assignSource(row(mapIt->second.i2,1), seg2_align);

		new_connect = _gotoh_trace_rectangle(mapIt->second.i2,  matr_it, matrix_diag, matrix_vert, matrix_hori, width_stop, height_stop);
		score_str_diag[j] = *matr_it;
		score_str_vert[j] = getValue(matrix_vert,position(matr_it));
		score_str_hori[j] = getValue(matrix_hori,position(matr_it));

		mapIt->second.i1 = new_connect;
		goNext(matr_it, 1);
		++height_align;
		--height_stop;
		if (old_connect != new_connect){
			_deleteAlignment(alignmentVector, old_connect, new_connect);
		}
		old_connect = new_connect;
	}
	for(TValue j = height_diag; j < overall; ++j)
	{
		alignmentVector.back().insert(std::make_pair(j,Pair<TValue, TAlign> ()));
		TMapIterator mapIt = --(alignmentVector.back().end());
		TSegment seg1_align = infix(host(*p_seq1), width_align, getBeginDim0(seed) + w_d2);
		TSegment seg2_align = infix(host(*p_seq2), height_align, getBeginDim1(seed)+ h_d2);
		resize(rows(mapIt->second.i2),2);
		assignSource(row(mapIt->second.i2,0), seg1_align);
		assignSource(row(mapIt->second.i2,1), seg2_align);
		new_connect = _gotoh_trace_rectangle(mapIt->second.i2,  matr_it, matrix_diag, matrix_vert, matrix_hori, width_stop, height_stop);
		score_str_diag[j] = *matr_it;
		score_str_vert[j] = getValue(matrix_vert,position(matr_it));
		score_str_hori[j] = getValue(matrix_hori,position(matr_it));
		mapIt->second.i1 = new_connect;
		goPrevious(matr_it,0);
		--width_align;
		++width_stop;
		if (old_connect != new_connect)
			_deleteAlignment(alignmentVector, old_connect, new_connect);
		old_connect = new_connect;
	}
	_deleteAlignment(alignmentVector, old_connect, score_length);
	score_length = length(score_str_diag);
}

// Rectangle calculation between two seeds
template <typename TScoreValue, unsigned DIMENSION, typename TSeedSpec, typename TSeedConfig, typename TString, typename TValue2>
void
_gotoh_rectangle(Matrix<TScoreValue, DIMENSION> & matrix_diag,
				 Matrix<TScoreValue, DIMENSION> & matrix_vert,
				 Matrix<TScoreValue, DIMENSION> & matrix_hori,	//edit matrix
				 Seed<TSeedSpec, TSeedConfig> const &seed1,				//Seed nearer to the end
				 Seed<TSeedSpec, TSeedConfig> const &seed2,				//Seed nearer to the start
				 TValue2 k_begin,									//upper diagonal extension
				 TValue2 k_end,										//lower diagonal extension
				 TString const & str1_,								//first sequence
				 TString const & str2_,								//secondSequence
				 Score<TScoreValue, Simple> const & score_,			//score matrix
				 String<TScoreValue> & init_diag,					//Values for initialisation
				 String<TScoreValue> & init_vert,					//Values for initialisation
				 String<TScoreValue> & init_hori)					//Values for initialisation
{
SEQAN_CHECKPOINT


	typedef Matrix<TScoreValue, DIMENSION> TMatrix;

	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Iterator<TMatrix, PositionIterator>::Type TMatrixIterator;

	typedef typename Iterator<TString const, Rooted>::Type TStringIterator;
	typedef typename Value<TString const>::Type TAlphabet;

	//-------------------------------------------------------------------------
	//define some variables

	TScoreValue inf = -1000000;

	TSize diag_height1 = getLowerDiagonal(seed1) - getStartDiagonal(seed1) + k_begin;
	TSize diag_width1 = getStartDiagonal(seed1) - getUpperDiagonal(seed1) + k_begin;

	TSize diag_height2 = getEndDiagonal(seed2) - getUpperDiagonal(seed2) + k_end;
	TSize diag_width2 = getLowerDiagonal(seed2) - getEndDiagonal(seed2) + k_end;

	TSize rectangle_height = getBeginDim1(seed1) - (getEndDim1(seed2) - 1) -1;
	TSize rectangle_width = length(str1_);

	TScoreValue score_match = scoreMatch(score_);
	TScoreValue score_mismatch = scoreMismatch(score_);
	TScoreValue score_gap_extend = scoreGapExtend(score_);
	TScoreValue score_gap_open = scoreGapOpen(score_);

	TSize str1_length = length(str1_);
	TSize str2_length = length(str2_);
	setLength(matrix_diag, 0, str1_length + 1);
	setLength(matrix_diag, 1, str2_length + 1);
	setLength(matrix_vert, 0, str1_length + 1);
	setLength(matrix_vert, 1, str2_length + 1);
	setLength(matrix_hori, 0, str1_length + 1);
	setLength(matrix_hori, 1, str2_length + 1);
	
	//cout << str1_length << " "<<str2_length << endl;
	resize(matrix_diag);
	//cout << "ende" << endl;
	resize(matrix_vert);
	resize(matrix_hori);

	TMatrixIterator col_diag = begin(matrix_diag);
	TMatrixIterator col_hori = begin(matrix_hori);
	TMatrixIterator col_vert = begin(matrix_vert);
	//Initialisierung

	setPosition(col_diag, str1_length);
	setPosition(col_hori, str1_length);
	setPosition(col_vert, str1_length);
	for (int i = 0; i!= diag_height2+rectangle_height;++i)
	{
		*col_diag = inf;
		goNext(col_diag,1);
		*col_vert = inf;
		goNext(col_vert,1);
		*col_hori = inf;
		goNext(col_hori,1);
	}


	setPosition(col_diag, (diag_height2 + rectangle_height+1)*(str1_length+1)-1);
	setPosition(col_hori, (diag_height2 + rectangle_height+1)*(str1_length+1)-1);
	setPosition(col_vert, (diag_height2 + rectangle_height+1)*(str1_length+1)-1);
	
	for (int i = 0; i != diag_width1; ++i)
	{		
		*col_diag = init_diag[i];
		--col_diag;
		*col_vert = init_vert[i];
		--col_vert;
		*col_hori = init_hori[i];
		--col_hori;
	}

	TSize len = length(init_diag);
	for (int i = diag_width1; i != len; ++i)
	{
		*col_diag=init_diag[i];
		goNext(col_diag,1);
		*col_vert=init_vert[i];
		goNext(col_vert,1);
		*col_hori=init_hori[i];
		goNext(col_hori,1);
	}

	goPrevious(col_diag,1);
	TMatrixIterator finger2_diag =col_diag;
	--col_diag;

	goPrevious(col_vert,1);
	TMatrixIterator finger2_vert =col_vert;
	--col_vert;
	
	goPrevious(col_hori,1);
	TMatrixIterator finger2_hori =col_hori;
	--col_hori;


	TSize width_align = str1_length - diag_width1;
	for (int i = 0; i!= width_align; ++i)
	{
		*col_diag = inf;
		--col_diag;
		*col_vert = inf;
		--col_vert;
		*col_hori = inf;
		--col_hori;
	}

	//MAIN
	TStringIterator x_begin = begin(str1_);
	--x_begin;
	TStringIterator x_end = end(str1_) -1;
	setPosition(x_end,width_align-1);
	TStringIterator y_begin = begin(str2_);
	TStringIterator y_end = end(str2_)-1;
	TStringIterator x,y;

	col_diag = finger2_diag;	
	col_vert = finger2_vert;
	col_hori = finger2_hori;

	TMatrixIterator finger1;
	
	
	TScoreValue hori_value, vert_value, diag_value_down, diag_value_right, diag_value_tmp;
	TMatrixIterator finger_diag_down, finger_vert_down;

	for (int i = 0; i < diag_height1; ++i)
	{
		finger2_hori = col_hori;
		goPrevious(finger2_hori,1);
		hori_value = *finger2_hori;
		--finger2_hori;
		goPrevious(col_hori,1);

		finger2_vert = col_vert;
		--finger2_vert;
		finger_vert_down = finger2_vert;
		goPrevious(finger2_vert,1);
		goPrevious(col_vert,1);

		finger2_diag = col_diag;
		diag_value_down = *finger2_diag;
		--finger2_diag;
		finger_diag_down = finger2_diag;
		goPrevious(finger2_diag,1);
		goPrevious(col_diag,1);
		diag_value_right = *col_diag;

		for (x = x_end; x != x_begin; --x)
		{
			diag_value_tmp = *finger_diag_down;

			hori_value = (hori_value+score_gap_extend > diag_value_right+score_gap_open) ? hori_value+score_gap_extend : diag_value_right+score_gap_open;
			*finger2_hori = hori_value;
			--finger2_hori;
			
			vert_value = *finger_vert_down;
			vert_value = (vert_value+score_gap_extend > diag_value_tmp+score_gap_open) ? vert_value+score_gap_extend : diag_value_tmp+score_gap_open;
			*finger2_vert = vert_value;
			--finger2_vert;
			--finger_vert_down;

			TScoreValue tmp_diag = hori_value > vert_value ? hori_value : vert_value;
			diag_value_down = diag_value_down + ((*x == *y_end) ? score_match : score_mismatch);
			if (diag_value_down < tmp_diag)
				diag_value_down = tmp_diag;
			diag_value_right= diag_value_down;
			*finger2_diag = diag_value_down;
			
			
			diag_value_down = diag_value_tmp;
			--finger2_diag;
			--finger_diag_down;

		}
		--y_end;
	}

	setPosition(col_diag,(diag_height2 + rectangle_height+1)*(rectangle_width+1)-1);
	setPosition(col_vert,(diag_height2 + rectangle_height+1)*(rectangle_width+1)-1);
	setPosition(col_hori,(diag_height2 + rectangle_height+1)*(rectangle_width+1)-1);

	x_end = end(str1_) -1;
	for (int i = 0; i < rectangle_height; ++i)
	{
		finger2_hori = col_hori;
		goPrevious(finger2_hori,1);
		hori_value = *finger2_hori;
		--finger2_hori;
		goPrevious(col_hori,1);

		finger2_vert = col_vert;
		--finger2_vert;
		finger_vert_down = finger2_vert;
		goPrevious(finger2_vert,1);
		goPrevious(col_vert,1);

		finger2_diag = col_diag;
		diag_value_down = *finger2_diag;
		--finger2_diag;
		finger_diag_down = finger2_diag;
		goPrevious(finger2_diag,1);
		goPrevious(col_diag,1);
		diag_value_right = *col_diag;

		for (x = x_end; x != x_begin; --x)
		{
			diag_value_tmp = *finger_diag_down;

			hori_value = (hori_value+score_gap_extend > diag_value_right+score_gap_open) ? hori_value+score_gap_extend : diag_value_right+score_gap_open;
			*finger2_hori = hori_value;
			--finger2_hori;
			
			vert_value = *finger_vert_down;
			vert_value = (vert_value+score_gap_extend > diag_value_tmp+score_gap_open) ? vert_value+score_gap_extend : diag_value_tmp+score_gap_open;
			*finger2_vert = vert_value;
			--finger2_vert;
			--finger_vert_down;

			TScoreValue tmp_diag = hori_value > vert_value ? hori_value : vert_value;
			diag_value_down = diag_value_down + ((*x == *y_end) ? score_match : score_mismatch);
			if (diag_value_down < tmp_diag)
				diag_value_down = tmp_diag;
			diag_value_right= diag_value_down;
			*finger2_diag = diag_value_down;
			
			
			diag_value_down = diag_value_tmp;
			--finger2_diag;
			--finger_diag_down;
		}
		--y_end;
	}

	setPosition(x_begin,diag_width2-1);

	for (int i = 0; i < diag_height2; ++i)
	{
		finger2_hori = col_hori;
		goPrevious(finger2_hori,1);
		hori_value = *finger2_hori;
		--finger2_hori;
		goPrevious(col_hori,1);

		finger2_vert = col_vert;
		--finger2_vert;
		finger_vert_down = finger2_vert;
		goPrevious(finger2_vert,1);
		goPrevious(col_vert,1);

		finger2_diag = col_diag;
		diag_value_down = *finger2_diag;
		--finger2_diag;
		finger_diag_down = finger2_diag;
		goPrevious(finger2_diag,1);
		goPrevious(col_diag,1);
		diag_value_right = *col_diag;

		for (x = x_end; x != x_begin; --x)
		{
			diag_value_tmp = *finger_diag_down;

			hori_value = (hori_value+score_gap_extend > diag_value_right+score_gap_open) ? hori_value+score_gap_extend : diag_value_right+score_gap_open;
			*finger2_hori = hori_value;
			--finger2_hori;
			
			vert_value = *finger_vert_down;
			vert_value = (vert_value+score_gap_extend > diag_value_tmp+score_gap_open) ? vert_value+score_gap_extend : diag_value_tmp+score_gap_open;
			*finger2_vert = vert_value;
			--finger2_vert;
			--finger_vert_down;

			TScoreValue tmp_diag = hori_value > vert_value ? hori_value : vert_value;
			diag_value_down = diag_value_down + ((*x == *y_end) ? score_match : score_mismatch);
			if (diag_value_down < tmp_diag)
				diag_value_down = tmp_diag;
			diag_value_right= diag_value_down;
			*finger2_diag = diag_value_down;
			
			
			diag_value_down = diag_value_tmp;
			--finger2_diag;
			--finger_diag_down;
		}
		--y_end;
	}
}


template <typename TScoreValue, unsigned DIMENSION, typename TString>
void
_gotoh2(Matrix<TScoreValue, DIMENSION> & matrix_diag,		//edit matrix
		Matrix<TScoreValue, DIMENSION> & matrix_vert,		//edit matrix
		Matrix<TScoreValue, DIMENSION> & matrix_hori,		//edit matrix
		TString const & str1_,								//first sequence
		TString const & str2_,								//secondSequence
		Score<TScoreValue, Simple> const & score_)			//score matrix
{
	SEQAN_CHECKPOINT

	typedef Matrix<TScoreValue, DIMENSION> TMatrix;
	typedef typename Size<TMatrix> ::Type TValue;

	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Iterator<TMatrix, PositionIterator>::Type TMatrixIterator;

	typedef typename Iterator<TString const, Rooted>::Type TStringIterator;
	typedef typename Value<TString const>::Type TAlphabet;

	//-------------------------------------------------------------------------
	//define some variables

	TScoreValue inf = -1000000;

	TScoreValue score_match = scoreMatch(score_);
	TScoreValue score_mismatch = scoreMismatch(score_);
	TScoreValue score_gap_extend = scoreGapExtend(score_);
	TScoreValue score_gap_open = scoreGapOpen(score_);

	

	TValue str1_length = length(str1_);
	TValue str2_length = length(str2_);
	setLength(matrix_diag, 0, str1_length + 1);
	setLength(matrix_diag, 1, str2_length + 1);
	setLength(matrix_vert, 0, str1_length + 1);
	setLength(matrix_vert, 1, str2_length + 1);
	setLength(matrix_hori, 0, str1_length + 1);
	setLength(matrix_hori, 1, str2_length + 1);
	
		//cout << str1_length << " "<<str2_length << endl;
	resize(matrix_diag);
	//cout << "ende" << endl;
	resize(matrix_vert);
	resize(matrix_hori);

	TMatrixIterator col_diag = begin(matrix_diag);
	TMatrixIterator col_hori = begin(matrix_hori);
	TMatrixIterator col_vert = begin(matrix_vert);
	TMatrixIterator finger2_diag, finger2_vert, finger2_hori;
	//Initialisierung

	setPosition(col_diag, str1_length);
	setPosition(col_hori, str1_length);
	setPosition(col_vert, str1_length);
	
	TScoreValue border = (str2_length-1) * score_gap_extend + score_gap_open; 
	
	for (int i = str2_length; i!= 0 ;--i)
	{
		*col_diag = border;
		goNext(col_diag,1);
		*col_vert = border;
		goNext(col_vert,1);
		*col_hori = inf;
		goNext(col_hori,1);
		border -= score_gap_extend;
	}


	*col_diag = 0;
	*col_hori = inf;
	*col_hori = inf;
	border = score_gap_open;

	for (unsigned int i = 0; i!= str1_length; ++i)
	{
		--col_diag;
		*col_diag = border;
		--col_vert;
		*col_vert = inf;
		--*col_hori;
		*col_hori = border;
		border += score_gap_extend;
	}

	col_diag = end(matrix_diag);
	--col_diag;
	col_vert = end(matrix_vert);
	--col_vert;
	col_hori = end(matrix_hori);
	--col_hori;

	//MAIN
	TStringIterator x_begin = begin(str1_);
	--x_begin;
	TStringIterator x_end = end(str1_) -1;
	
	TStringIterator y_begin = begin(str2_);
	TStringIterator y_end = end(str2_)-1;
	TStringIterator x,y;

	finger2_diag = col_diag;	
	finger2_vert = col_vert;
	finger2_hori = col_hori;

	TMatrixIterator finger1;
	
	
	TScoreValue hori_value, vert_value, diag_value_down, diag_value_right, diag_value_tmp;
	TMatrixIterator finger_diag_down, finger_vert_down;

	x_end = end(str1_) -1;
	for (unsigned int i = 0; i < str2_length; ++i)
	{
		finger2_hori = col_hori;
		goPrevious(finger2_hori,1);
		hori_value = *finger2_hori;
		--finger2_hori;
		goPrevious(col_hori,1);

		finger2_vert = col_vert;
		--finger2_vert;
		finger_vert_down = finger2_vert;
		goPrevious(finger2_vert,1);
		goPrevious(col_vert,1);

		finger2_diag = col_diag;
		diag_value_down = *finger2_diag;
		--finger2_diag;
		finger_diag_down = finger2_diag;
		goPrevious(finger2_diag,1);
		goPrevious(col_diag,1);
		diag_value_right = *col_diag;

		for (x = x_end; x != x_begin; --x)
		{
			diag_value_tmp = *finger_diag_down;

			hori_value = (hori_value+score_gap_extend > diag_value_right+score_gap_open) ? hori_value+score_gap_extend : diag_value_right+score_gap_open;
			*finger2_hori = hori_value;
			--finger2_hori;
			
			vert_value = *finger_vert_down;
			vert_value = (vert_value+score_gap_extend > diag_value_tmp+score_gap_open) ? vert_value+score_gap_extend : diag_value_tmp+score_gap_open;
			*finger2_vert = vert_value;
			--finger2_vert;
			--finger_vert_down;

			TScoreValue tmp_diag = hori_value > vert_value ? hori_value : vert_value;
			diag_value_down = diag_value_down + ((*x == *y_end) ? score_match : score_mismatch);
			if (diag_value_down < tmp_diag)
				diag_value_down = tmp_diag;
			diag_value_right= diag_value_down;
			*finger2_diag = diag_value_down;
			
			
			diag_value_down = diag_value_tmp;
			--finger2_diag;
			--finger_diag_down;
		}
		--y_end;
	}
}

//changed version of usual trace-back
template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, unsigned DIMENSION, typename TValue, typename TMatrix>
TValue
_gotoh_trace_rectangle(Align<TTargetSource, TTargetSpec> & target_,
						Iter< Matrix<TScoreValue, DIMENSION>, PositionIterator > & diag_source,
						TMatrix & diag_matrix_, 
						TMatrix & vert_matrix_, 
						TMatrix & hori_matrix_, 
						TValue width_stop, 
						TValue height_stop)
{
	typedef typename Position<Matrix<TScoreValue, DIMENSION> >::Type TPosition;

	typedef Align<TTargetSource, TTargetSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TTargetIterator;
    typedef typename Infix<TTargetSource>::Type TTargetSourceSegment;
	typedef typename Iterator<TTargetSourceSegment, Rooted>::Type TStringIterator;

	TTargetSourceSegment str_0 = sourceSegment(row(target_, 0));
	TTargetSourceSegment str_1 = sourceSegment(row(target_, 1));
	typename Size<TTargetSourceSegment>::Type dim_0_len = length(hori_matrix_,0);

	TPosition pos = position(diag_source);

	TTargetIterator target_0 = iter(row(target_, 0), 0, Standard());
	TTargetIterator target_1 = iter(row(target_, 1), 0, Standard());

	TStringIterator it_0 = iter(str_0, 0);
	TStringIterator it_1 = iter(str_1, 0);

	while ((static_cast<TValue>(position(it_0)) < width_stop) || (static_cast<TValue>(position(it_1)) < height_stop))
	{
		//cout << "arg" << endl;
		if (getValue(diag_matrix_,pos) == getValue(hori_matrix_,pos))
		{
			++it_0;
			insertGap(target_1);
			++pos;
		}
		else
		{
			if (getValue(diag_matrix_,pos) != getValue(vert_matrix_,pos)) 
			{					
				++it_0;
				++it_1;
				pos += dim_0_len+1;
			} 
			else
			{
				++it_1;
				insertGap(target_0);
				pos += dim_0_len;
			}
		}
		++target_0;
		++target_1;
	}

	//-------------------------------------------------------------------------
	//follow the trace until the border is reached




	setSourceEndPosition(row(target_,1),position(it_1));
	setSourceEndPosition(row(target_,0),position(it_0));

	return length(diag_matrix_,0) - (pos % dim_0_len) + position(it_1) - height_stop -1;
}


//Rectangle calculation between two seeds
template <typename TScoreValue, unsigned DIMENSION, typename TSeedSpec, typename TSeedConfig, typename TString, typename TValue2>
void
_banded_gotoh_rectangle_first(Matrix<TScoreValue, DIMENSION> & matrix_diag,	//edit matrix
							  Matrix<TScoreValue, DIMENSION> & matrix_vert,
							  Matrix<TScoreValue, DIMENSION> & matrix_hori,
							  Seed<TSeedSpec, TSeedConfig> const &seed,				//Seed
							  TValue2 k,											//diagonal extension
							  TString const & str1_,							//first sequence
							  TString const & str2_,							//secondSequence
							  Score<TScoreValue, Simple> const & score_,		//score matrix
							  String<TScoreValue> & init_diag,					//Values for initialisation
							  String<TScoreValue> & init_vert,
							  String<TScoreValue> & init_hori)					
{
SEQAN_CHECKPOINT


	typedef Matrix<TScoreValue, DIMENSION> TMatrix;

	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Iterator<TMatrix, PositionIterator>::Type TMatrixIterator;

	typedef typename Iterator<TString const, Rooted>::Type TStringIterator;
	typedef typename Value<TString const>::Type TAlphabet;

	//-------------------------------------------------------------------------
	//define some variables
	TScoreValue inf = -1000000;

	TSize diag_height1 = getLowerDiagonal(seed) - getStartDiagonal(seed) + k;
	TSize diag_width1 = getStartDiagonal(seed) - getUpperDiagonal(seed) + k;

	TSize rectangle_height = getBeginDim1(seed) - beginPosition(str1_);
	TSize rectangle_width = length(str1_);

	TScoreValue score_match = scoreMatch(score_);
	TScoreValue score_mismatch = scoreMismatch(score_);
	TScoreValue score_gap_extend = scoreGapExtend(score_);
	TScoreValue score_gap_open = scoreGapOpen(score_);

	TSize str1_length = length(str1_);
	TSize str2_length = length(str2_);
	setLength(matrix_diag, 0, str1_length + 1);
	setLength(matrix_diag, 1, str2_length + 1);
	setLength(matrix_vert, 0, str1_length + 1);
	setLength(matrix_vert, 1, str2_length + 1);
	setLength(matrix_hori, 0, str1_length + 1);
	setLength(matrix_hori, 1, str2_length + 1);

	resize(matrix_diag);
	resize(matrix_vert);
	resize(matrix_hori);

	TMatrixIterator col_diag = begin(matrix_diag);
	TMatrixIterator col_hori = begin(matrix_hori);
	TMatrixIterator col_vert = begin(matrix_vert);

	//Initialisierung

	setPosition(col_diag, str1_length);
	setPosition(col_hori, str1_length);
	setPosition(col_vert, str1_length);
	for (int i = 0; i!= rectangle_height+1;++i)
	{
		*col_diag = inf;
		goNext(col_diag,1);
		*col_vert = inf;
		goNext(col_vert,1);
		*col_hori = inf;
		goNext(col_hori,1);
	}

	setPosition(col_diag, (rectangle_height+1)*(str1_length+1)-1);
	setPosition(col_hori, (rectangle_height+1)*(str1_length+1)-1);
	setPosition(col_vert, (rectangle_height+1)*(str1_length+1)-1);
	
	for (int i = 0; i != diag_width1; ++i)
	{		
		*col_diag = init_diag[i];
		--col_diag;
		*col_vert = init_vert[i];
		--col_vert;
		*col_hori = init_hori[i];
		--col_hori;
	}

	TSize len = length(init_diag);
	for (int i = diag_width1; i != len; ++i)
	{
		*col_diag=init_diag[i];
		goNext(col_diag,1);
		*col_vert=init_vert[i];
		goNext(col_vert,1);
		*col_hori=init_hori[i];
		goNext(col_hori,1);
	}

	goPrevious(col_diag,1);
	TMatrixIterator finger2_diag =col_diag;
	--col_diag;

	goPrevious(col_vert,1);
	TMatrixIterator finger2_vert =col_vert;
	--col_vert;
	
	goPrevious(col_hori,1);
	TMatrixIterator finger2_hori =col_hori;
	--col_hori;


	TSize width_align = str1_length - diag_width1;
	for (int i = 0; i!= width_align; ++i)
	{
		*col_diag = inf;
		--col_diag;
		*col_vert = inf;
		--col_vert;
		*col_hori = inf;
		--col_hori;
	}


	//MAIN
	TStringIterator x_begin = begin(str1_);
	--x_begin;
	TStringIterator x_end = end(str1_) -1;
	setPosition(x_end,width_align-1);
	TStringIterator y_begin = begin(str2_);
	TStringIterator y_end = end(str2_)-1;
	TStringIterator x,y;

	col_diag = finger2_diag;	
	col_vert = finger2_vert;
	col_hori = finger2_hori;

	TMatrixIterator finger1;
	
	
	TScoreValue hori_value, vert_value, diag_value_down, diag_value_right, diag_value_tmp;
	TMatrixIterator finger_diag_down, finger_vert_down;

	for (int i = 0; i < diag_height1; ++i)
	{
		finger2_hori = col_hori;
		goPrevious(finger2_hori,1);
		hori_value = *finger2_hori;
		--finger2_hori;
		goPrevious(col_hori,1);

		finger2_vert = col_vert;
		--finger2_vert;
		finger_vert_down = finger2_vert;
		goPrevious(finger2_vert,1);
		goPrevious(col_vert,1);

		finger2_diag = col_diag;
		diag_value_down = *finger2_diag;
		--finger2_diag;
		finger_diag_down = finger2_diag;
		goPrevious(finger2_diag,1);
		goPrevious(col_diag,1);
		diag_value_right = *col_diag;

		for (x = x_end; x != x_begin; --x)
		{
			diag_value_tmp = *finger_diag_down;

			hori_value = (hori_value+score_gap_extend > diag_value_right+score_gap_open) ? hori_value+score_gap_extend : diag_value_right+score_gap_open;
			*finger2_hori = hori_value;
			--finger2_hori;
			
			vert_value = *finger_vert_down;
			vert_value = (vert_value+score_gap_extend > diag_value_tmp+score_gap_open) ? vert_value+score_gap_extend : diag_value_tmp+score_gap_open;
			*finger2_vert = vert_value;
			--finger2_vert;
			--finger_vert_down;

			TScoreValue tmp_diag = hori_value > vert_value ? hori_value : vert_value;
			diag_value_down = diag_value_down + ((*x == *y_end) ? score_match : score_mismatch);
			if (diag_value_down < tmp_diag)
				diag_value_down = tmp_diag;
			diag_value_right= diag_value_down;
			*finger2_diag = diag_value_down;
			
			
			diag_value_down = diag_value_tmp;
			--finger2_diag;
			--finger_diag_down;

		}
		--y_end;
	}

	setPosition(col_diag,(rectangle_height+1)*(rectangle_width+1)-1);
	setPosition(col_vert,(rectangle_height+1)*(rectangle_width+1)-1);
	setPosition(col_hori,(rectangle_height+1)*(rectangle_width+1)-1);

	x_end = end(str1_) -1;
	for (int i = 0; i < rectangle_height; ++i)
	{
		finger2_hori = col_hori;
		goPrevious(finger2_hori,1);
		hori_value = *finger2_hori;
		--finger2_hori;
		goPrevious(col_hori,1);

		finger2_vert = col_vert;
		--finger2_vert;
		finger_vert_down = finger2_vert;
		goPrevious(finger2_vert,1);
		goPrevious(col_vert,1);

		finger2_diag = col_diag;
		diag_value_down = *finger2_diag;
		--finger2_diag;
		finger_diag_down = finger2_diag;
		goPrevious(finger2_diag,1);
		goPrevious(col_diag,1);
		diag_value_right = *col_diag;

		for (x = x_end; x != x_begin; --x)
		{
			diag_value_tmp = *finger_diag_down;

			hori_value = (hori_value+score_gap_extend > diag_value_right+score_gap_open) ? hori_value+score_gap_extend : diag_value_right+score_gap_open;
			*finger2_hori = hori_value;
			--finger2_hori;
			
			vert_value = *finger_vert_down;
			vert_value = (vert_value+score_gap_extend > diag_value_tmp+score_gap_open) ? vert_value+score_gap_extend : diag_value_tmp+score_gap_open;
			*finger2_vert = vert_value;
			--finger2_vert;
			--finger_vert_down;

			TScoreValue tmp_diag = hori_value > vert_value ? hori_value : vert_value;
			diag_value_down = diag_value_down + ((*x == *y_end) ? score_match : score_mismatch);
			if (diag_value_down < tmp_diag)
				diag_value_down = tmp_diag;
			diag_value_right= diag_value_down;
			*finger2_diag = diag_value_down;
			
			
			diag_value_down = diag_value_tmp;
			--finger2_diag;
			--finger_diag_down;
		}
		--y_end;
	}
	//printMatrix(matrix_diag);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_CHAIN_BANDED_AFFINE_H_
