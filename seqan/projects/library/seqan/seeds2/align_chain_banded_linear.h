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
  Author: Birte Kehr <bkehr@fu-berlin.de>
  Author: Carsten Kemena <carsten.kemena@crg.es>
 ============================================================================
  Banded chain alignment around a chain of seeds; Case with linear gap costs.
  Based on the original code by Carsten Kemena, adapted to the new seeds
  interface.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_ALIGN_CHAIN_BANDED_LINEAR_H_
#define SEQAN_SEEDS_ALIGN_CHAIN_BANDED_LINEAR_H_

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

template <typename TContainer, typename TBandwidth, typename TScoreValue, typename TAlign>
TScoreValue
_bandedChainAlignment_NW(
        TAlign & alignment,
        TContainer const & seedChain, 
        TBandwidth k,
        Score<TScoreValue, Simple> const & scoringScheme)
{
	SEQAN_CHECKPOINT;
    // TODO(holtgrew): The alignment matrix is computed from the lower right to the upper left. Why?
    
    // -----------------------------------------------------------------------
    // Function-Global Typedefs
    // -----------------------------------------------------------------------

    typedef typename Source<TAlign>::Type TSequence;
    typedef typename Infix<TSequence>::Type TInfix;
	typedef typename Size<TSequence>::Type TSize;
	typedef typename Position<TSequence>::Type TPosition;

	typedef typename Iterator<TContainer const, Standard>::Type TIterator;

	typedef String<TScoreValue> TScoreString;

	typedef Matrix<TScoreValue> TMatrix;

    // TODO(holtgrew): Using maps appears to be pretty heavyweight...
    // TODO(holtgrew): Is TBandwidth correct here?
    // TODO(holtgrew): What is the map used for here?
	typedef ::std::vector< ::std::map<TBandwidth, Pair<TBandwidth, TAlign> > > TAlignVector;
	
    // -----------------------------------------------------------------------
    // Declare variables, Simple initialization
    // -----------------------------------------------------------------------
	TScoreString score_str;
	TAlignVector alignmentVector;

	TInfix seq1 = sourceSegment(row(alignment, 0));
	TInfix seq2 = sourceSegment(row(alignment, 1));
    TInfix *p_seq1 = &seq1;  // TODO(holtgrew): Why not copy infixes? Cheap...
	TInfix *p_seq2 = &seq2;
	
	TBandwidth score_length = 0;
	TMatrix matrix_;
	TIterator it = begin(seedChain);

    // -----------------------------------------------------------------------
    // Banded Alignment
    // -----------------------------------------------------------------------
    //
    // The calculation begins at the end.  We start with computing the
    // last rectangle and the last seed.  Then, we perform the banded
    // alignment in the gaps and along the seeds.
	TBandwidth k_begin = k;
	TBandwidth k_end = k;
	if (length(*it) <= k || (getEndDim1(*it) - 1 - getBeginDim1(*it)) <= k)
 		k_end = length(*it) - 1;

	_calculateLastRectangle(*it, k_end, matrix_, p_seq1, p_seq2, score_str, score_length, alignmentVector, scoringScheme);
	_calculateBandedSeed(*it, k_end, matrix_, p_seq1, p_seq2, score_str, score_length, alignmentVector, scoringScheme);

    // Compute alignment around all but the first and last rectangle
    // and last seed.
	TIterator it_begin = end(seedChain);
	TIterator it2 = it;
	++it2;
	while (it2 != it_begin) {
		k_begin = k_end;
		if (length(*it2) <= k || (getEndDim1(*it2) - 1 - getBeginDim1(*it2)) <= k)
			k_end = length(*it2) - 2;
		else 
			k_end = k;
		
		_calculateRectangle(*it, *it2, k_begin, k_end, matrix_, p_seq1, p_seq2, score_str, score_length, alignmentVector, scoringScheme);	
		_calculateBandedSeed(*it2, k_end, matrix_, p_seq1, p_seq2, score_str, score_length, alignmentVector, scoringScheme);
		++it;
		++it2;
	}

    // Compute alignment of first rectangle.
	_calculateFirstRectangle(*it, k_end, matrix_, p_seq1, p_seq2, score_str, score_length, alignmentVector, scoringScheme);

    // -----------------------------------------------------------------------
    // Construct Alignment, Return Score
    // -----------------------------------------------------------------------
	_constructAlignment(alignmentVector, alignment);
	return score_str[0];
}

// TODO(holtgrew): _constructAlignment is shared with affine gap costs.
//
//"Glues" single alignments together
template<typename TValue, typename TAlign, typename TAlign2>
void
// TODO(holtgrew): wholeAlignment should be named alignment and be the first alignment.
_constructAlignment(::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > const & me,
					TAlign2 & wholeAlignment)
{
	SEQAN_CHECKPOINT;

	typedef typename ::std::map<TValue,Pair<TValue, TAlign> >::const_iterator TIterator;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TTargetIterator;
	typedef typename Row<TAlign2>::Type TRow2;
	typedef typename Iterator<TRow2, Standard>::Type TTargetIterator2;

	TTargetIterator2 align_it0 = iter(row(wholeAlignment, 0), 0);
	TTargetIterator2 align_it1 = iter(row(wholeAlignment, 1), 0);
	//TValue position = 0;
	int length_ = me.size()-1;
	TIterator it;

	typedef typename Position<TAlign>::Type TPosition;
	for (int i = length_; i != -1; --i)
	{
		//cout << "LENGHT: " << me[i].size() << endl;
		it = me[i].begin();//find(position);
		//position = it->second.i1;
		TRow tmp_row0 = row(it->second.i2, 0);
		TRow tmp_row1 = row(it->second.i2, 1);
		TPosition end_ = endPosition(cols(it->second.i2));
		unsigned int j = 0;
		for (TPosition it_ = beginPosition(cols(it->second.i2)); it_ != end_; ++it_)
		{
			if(isGap(tmp_row0, j))
				insertGap(align_it0);
			if(isGap(tmp_row1, j))
				insertGap(align_it1);
			++align_it0;
			++align_it1;
			++j;
		}
	}
}


// TODO(holtgrew): This is more or less the same as _bandedAlignment_NW_align() from banded seed alignment! Make a diff of the original code and create a more general function.
template <typename TScoreValue, unsigned DIMENSION, typename TSeedSpec, typename TSeedConfig, typename TString, typename TValue2>
TScoreValue
_banded_needleman_wunsch(Matrix<TScoreValue, DIMENSION> & matrix_,
						 Seed<TSeedSpec, TSeedConfig> const & seed,
						 TValue2 k,
						 TString const & str1_,
						 TString const & str2_,
						 Score<TScoreValue, Simple> const & score_,
						 String<TScoreValue> init)
{
	SEQAN_CHECKPOINT;

	typedef Matrix<TScoreValue, DIMENSION> TMatrix;

	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Iterator<TMatrix, Standard>::Type TMatrixIterator;

	typedef typename Iterator<TString const, Rooted>::Type TStringIterator;
	typedef typename Value<TString const>::Type TAlphabet;

	//-------------------------------------------------------------------------
	//define some variables
	TSize height = getLowerDiagonal(seed) - getUpperDiagonal(seed)+1;
	TSize str1_length = height+2*k;
	TSize str2_length = length(str2_);

	TStringIterator x_begin = begin(str1_);
	TStringIterator x_end = end(str1_) -1;
	TStringIterator y_begin = begin(str2_);
	--y_begin;
	TStringIterator y_end = end(str2_)-1;
	
	TSize up_height = getLowerDiagonal(seed) - getStartDiagonal(seed) + k; //equals length of empty lower triangle
	//TValue up_width = str1_length - up_height;
	TSize down_height = getEndDiagonal(seed) - getUpperDiagonal(seed) + k; //equals length of empty lower triangle
	TSize down_width = str1_length - down_height;
	TSize length_right_diag = str2_length - down_height;
	TSize length_left_diag = str2_length - up_height;

	TScoreValue score_match = scoreMatch(score_);
	TScoreValue score_mismatch = scoreMismatch(score_);
	TScoreValue score_gap = scoreGapExtend(score_);

	TScoreValue horizontalValue = 0;
	TScoreValue border_ = score_gap;
	TScoreValue verticalValue = border_;

	setDimension(matrix_, 2);
	setLength(matrix_, 0, str1_length + 2);
	setLength(matrix_, 1, str2_length + 1);
	resize(matrix_);
	
	TMatrixIterator col_ = begin(matrix_);
	TMatrixIterator finger1;
	TMatrixIterator finger2;
	TMatrixIterator finger3;

	TStringIterator x = x_begin;
	TStringIterator y = y_begin;

	//-------------------------------------------------------------------------
	// init

	TScoreValue inf = -1000000;
	finger2 = col_;

	setPosition(finger2, length(matrix_, 0)-1);
	for (unsigned int i = 1; i != length_right_diag; ++i){
		*finger2 = inf;
		goNext(finger2, 1);
	}

	TScoreValue pos = -1;  // TODO(holtgrew): Shouldn't this be TPosition? Why negative?

	*finger2 = inf;
	for (int i = -1; i != down_height; ++i){
		goPrevious(finger2, 0);
		goNext(finger2,1);
		*finger2 = init[++pos];
	}
		
	finger3 = finger2;

	for (int i = 0; i != down_width; ++i){
		goPrevious(finger2, 0);
		*finger2 = init[++pos];
	}

	*finger2 = inf;
	TScoreValue tmp;
	//-------------------------------------------------------------------------
	//first section
	
	y=y_end;
	TSize run_length = down_width;
	TSize measure = 0;

	for (int i = -1; i != down_height; ++i){	
		verticalValue = *finger3;
		finger1 = finger3;
		goPrevious(finger1,0);
		goPrevious(finger3,1);
		finger2 = finger3;
		goNext(finger3,0);
		horizontalValue = *finger3;
		x = x_end;
		for (unsigned int j = 0; j != run_length; ++j){
			if (*x==*y){
				*finger2 = verticalValue + score_match;
			} else {
				tmp = *finger1;
				TScoreValue s1 = verticalValue + score_mismatch;
				TScoreValue s2 = score_gap + ((horizontalValue > tmp) ? horizontalValue : tmp);
				*finger2 = (s1 > s2) ? s1 : s2;
			}
			horizontalValue = *finger2;
			goPrevious(finger2);
			verticalValue = *finger1;
			goPrevious(finger1);
			--x;
		}
		*finger2 = inf;
		++measure;
		if (measure < length_left_diag)
			++run_length;
		--y;
	}
	--run_length;

	goPrevious(finger3);

	while (y!= y_begin)
	{	
		verticalValue = *finger3;
		finger1 = finger3;
		goPrevious(finger1,0);
		goPrevious(finger3,1);
		finger2 = finger3;
		x = --x_end;
		horizontalValue = inf;
		for (unsigned int j = 0; j != run_length; ++j){
			if (*x==*y){
				*finger2 = verticalValue + score_match;
			} else {
				tmp = *finger1;
				TScoreValue s1 = verticalValue + score_mismatch;
				TScoreValue s2 = score_gap + ((horizontalValue > tmp) ? horizontalValue : tmp);
				*finger2 = (s1 > s2) ? s1 : s2;
			}
			horizontalValue = *finger2;
			goPrevious(finger2);
			verticalValue = *finger1;
			goPrevious(finger1);
			--x;
		}
		*finger2 = inf;
		++measure;
		if (measure >= length_left_diag)
			--run_length;
		--y;
	}
	return *(++finger2);

}


//////////////////////////////////////////////////////////////////////////////
//traceback through needleman wunsch matrix


//Position berechnen!
template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, unsigned DIMENSION>
typename Size<Matrix<TScoreValue, DIMENSION> >::Type
_banded_needleman_wunsch_trace2(Align<TTargetSource, TTargetSpec> & target_,
								Matrix<TScoreValue, DIMENSION> & matrix,
								Iter< Matrix<TScoreValue, DIMENSION>, PositionIterator > source_,
								Score<TScoreValue, Simple> const & score_)
{
	SEQAN_CHECKPOINT;

	typedef Iter<Matrix<TScoreValue, DIMENSION>, PositionIterator > TMatrixIterator;
    typedef typename Infix<TTargetSource>::Type TTargetSourceSegment;

	TTargetSourceSegment str_0 = sourceSegment(row(target_, 0));
	TTargetSourceSegment str_1 = sourceSegment(row(target_, 1));


	typedef Align<TTargetSource, TTargetSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TTargetIterator;
	TTargetIterator target_0 = iter(row(target_, 0), 0);
	TTargetIterator target_1 = iter(row(target_, 1), 0);



	typedef typename Iterator<TTargetSourceSegment, Rooted>::Type TStringIterator;
	TStringIterator it_0 = begin(str_0);
	TStringIterator it_0_end = end(str_0);

	TStringIterator it_1 = begin(str_1);
	TStringIterator it_1_end = end(str_1);

	TScoreValue score_diff = scoreMismatch(score_) - scoreGapExtend(score_);


	//-------------------------------------------------------------------------
	//follow the trace until the border is reached
	while ((it_0 != it_0_end) && (it_1 != it_1_end))
	{
		bool gv;
		bool gh;

		if (*it_0 == *it_1)
		{
			gv = gh = true;
		}
		else
		{
			TMatrixIterator it_ = source_;

			goNext(it_, 0);
			TScoreValue h = *it_;

			it_ = source_;
			goNext(it_, 1);
			TScoreValue d = *it_;
		
			goPrevious(it_, 0);
			TScoreValue v = *it_;

			gv = (v >= h) | (d + score_diff >= h);
			gh = (h >  v) | (d + score_diff >= v);

		}

	
		if (gv){
			++it_1;
			goNext(source_, 1);
			goPrevious(source_, 0);
		}else{
			insertGap(target_1);
		}

		if (gh) {
			++it_0;
			goNext(source_, 0);
	
		}else{
			insertGap(target_0);
		}
		++target_0;
		++target_1;
	}

	setSourceEndPosition(row(target_,1),position(it_1));
	setSourceEndPosition(row(target_,0),position(it_0));

	return length(matrix,0) - coordinate(source_,0) -2;
}


//changed version of usual trace-back
template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, unsigned DIMENSION>
int
_needleman_wunsch_trace_lastRectangle(Align<TTargetSource, TTargetSpec> & target_,
						Iter< Matrix<TScoreValue, DIMENSION>, PositionIterator > source_,
						Score<TScoreValue, Simple> const & score_)
{
	SEQAN_CHECKPOINT;

	typedef Iter<Matrix<TScoreValue, DIMENSION>, PositionIterator > TMatrixIterator;
    typedef typename Infix<TTargetSource>::Type TTargetSourceSegment;

	TTargetSourceSegment str_0 = sourceSegment(row(target_, 0));
	TTargetSourceSegment str_1 = sourceSegment(row(target_, 1));

	typedef Align<TTargetSource, TTargetSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TTargetIterator;
	TTargetIterator target_0 = iter(row(target_, 0), 0);
	TTargetIterator target_1 = iter(row(target_, 1), 0);

	typedef typename Iterator<TTargetSourceSegment, Standard>::Type TStringIterator;
	TStringIterator it_0 = iter(str_0, 0);
	TStringIterator it_0_end = end(str_0);

	TStringIterator it_1 = iter(str_1, 0);
	TStringIterator it_1_end = end(str_1);

	TScoreValue score_diff = scoreMismatch(score_) - scoreGapExtend(score_);

	//-------------------------------------------------------------------------
	//follow the trace until the border is reached

	while ((it_0 != it_0_end) && (it_1 != it_1_end))
	{
		bool gv;
		bool gh;

		if (*it_0 == *it_1)
		{
			gv = gh = true;
		}
		else
		{

			TMatrixIterator it_ = source_;

			goNext(it_, 0);
			TScoreValue v = *it_;

			goNext(it_, 1);
			TScoreValue d = *it_;

			it_ = source_;
			goNext(it_, 1);
			TScoreValue h = *it_;

			gv = (v >= h) | (d + score_diff >= h);
			gh = (h >  v) | (d + score_diff >= v);

		}

		if (gv)
		{
			++it_0;
			goNext(source_, 0);
		}
		else
		{
			insertGap(target_0);
		}

		if (gh) 
		{
			++it_1;
			goNext(source_, 1);
		}
		else
		{
			insertGap(target_1);
		}

		++target_0;
		++target_1;
	}
	return 0;
}

template<typename TAlign, typename TValue>
void
_rec_delete(::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > &vec,	//alignment vector
		   TValue index,										//position im vector
		   TValue position)										//alignment to delete
{
	SEQAN_CHECKPOINT;

	if (position != -1)
	{
		typename ::std::map<TValue,Pair<TValue, TAlign> >::iterator it, it1, it2;
		it = vec[index].find(position);
		bool x = true; //true = successor of it can be deleted
		if (it != vec[index].begin())
		{
			it1 = it;
			--it1;
			x = (it1->second.i1 != it->second.i1);
		}
		it2 = it;
		++it2;
		if (x && (it2 != vec[index].end()))
		{
			x = (it2->second.i1 != it->second.i1);
		}
		if (x)
			_rec_delete(vec, index-1, it->second.i1);
		vec[index].erase(position);
	}
}

template<typename TValue, typename TAlign, typename TSize>
void
_deleteAlignment(::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > &me,
				TSize old_end,
				TSize new_end)
{
	SEQAN_CHECKPOINT;

	//cout << "me: " << me[0].size() << endl;
	int length = me.size()-2;
	typename ::std::map<TValue,Pair<TValue, TAlign> >::iterator it, it1, it2;
	for(TSize i = old_end+1; i < new_end; ++i)
	{
	
		it = me[length].find(i);
		bool x = true; //true = successor of it can be deleted
		if (it != me[length].begin())
		{
			it1 = it;
			--it1;
			x = (it1->second.i1 != it->second.i1);
		}
		it2 = it;
		++it2;
		if (x && (it2 != me[length].end()))
		{
			x = (it2->second.i1 != it->second.i1);
		}
		if (x)
			_rec_delete(me, length-1, it->second.i1);
		me[length].erase(i);
		//cout << "me2: " << me[0].size() << endl;
	}
	
}

//calculation and backtracking of the banded alignment of a seed
template<typename TSeed, typename TString, typename TDiff, typename TMatrix, typename TScoreString, typename TValue, typename TAlign, typename TScoreMatrix>
void
_calculateBandedSeed(TSeed const &seed,
					 TDiff k,
					 TMatrix &matrix_,
					 TString *p_seq1,
					 TString *p_seq2,
					 TScoreString &score_str,
					 TValue &score_length,
					 ::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > &alignmentVector,
					 TScoreMatrix const &scoringScheme)
{
	SEQAN_CHECKPOINT;

	typedef typename ::std::map<TValue,Pair<TValue, TAlign> >::iterator TMapIterator;
	typedef Iter<TMatrix, PositionIterator> TMatrixIterator;
    typedef typename Infix<TString>::Type TSegment;
	TSegment seg1_align = infix(host(*p_seq1), getBeginDim0(seed), (getEndDim0(seed) - 1)+1);
	TSegment seg2_align = infix(host(*p_seq2), getBeginDim1(seed), (getEndDim1(seed) - 1)+1);
	_banded_needleman_wunsch(matrix_, seed, k, seg1_align, seg2_align, scoringScheme, score_str);

	TValue height_diag = getLowerDiagonal(seed)-getStartDiagonal(seed)+k;
	TValue width_diag = getStartDiagonal(seed)-getUpperDiagonal(seed)+k;
	TValue overall = height_diag + width_diag + 1;
	
	resize(score_str,overall);
	TMatrixIterator matr_it = begin(matrix_);
	setPosition(matr_it, length(matrix_,0)-2);

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
		new_connect = _banded_needleman_wunsch_trace2(mapIt->second.i2, matrix_,  matr_it, scoringScheme);
		score_str[j] = *matr_it;
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
		
		new_connect = _banded_needleman_wunsch_trace2(mapIt->second.i2, matrix_,  matr_it, scoringScheme);
		score_str[j] = *matr_it;
		mapIt->second.i1 = new_connect;
		goPrevious(matr_it,0);
		++height_align;
		if (old_connect != new_connect)
			_deleteAlignment(alignmentVector, old_connect, new_connect);
		old_connect = new_connect;
	}
	_deleteAlignment(alignmentVector, old_connect, score_length);
	score_length = length(score_str);
}

template<typename TSeed, typename TString, typename TDiff, typename TMatrix, typename TScoreString, typename TValue, typename TAlign, typename TScoreMatrix>
void
_calculateFirstRectangle(TSeed const &seed,
						 TDiff k,
						 TMatrix &matrix_,
						 TString *p_seq1,
						 TString *p_seq2,
						 TScoreString &score_str,
						 TValue &score_length,
						 ::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > &alignmentVector,
						 TScoreMatrix const &scoringScheme)
{
	SEQAN_CHECKPOINT;

	typedef typename ::std::map<TValue,Pair<TValue, TAlign> >::iterator TMapIterator;
	typedef Iter<TMatrix, PositionIterator> TMatrixIterator;
	TValue new_connect;
    typedef typename Infix<TString>::Type TSegment;
	TSegment seg1b_align = infix(host(*p_seq1), beginPosition(*p_seq1), getBeginDim0(seed) + getStartDiagonal(seed) - getUpperDiagonal(seed) + k);
	TSegment seg2b_align = infix(host(*p_seq2), beginPosition(*p_seq2), getBeginDim1(seed) + getLowerDiagonal(seed) - getStartDiagonal(seed) + k);

	_banded_needleman_wunsch_rectangle_first(matrix_, seed, k, seg1b_align, seg2b_align, scoringScheme,score_str);

	TValue w_d2 = getStartDiagonal(seed) - getUpperDiagonal(seed) + k;
	TValue h_d2 = getLowerDiagonal(seed) - getStartDiagonal(seed) + k;

	resize(score_str,1);

	TMatrixIterator matr_it = begin(matrix_);

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

	new_connect = _needleman_wunsch_trace_rectangle(mapIt->second.i2,  matr_it, scoringScheme, matrix_, width_stop, height_stop);
	score_str[0] = *matr_it;
	mapIt->second.i1 = new_connect;
	_deleteAlignment(alignmentVector, -1, new_connect);
	_deleteAlignment(alignmentVector, new_connect, score_length);
}

// TODO(holtgrew): The result is alignmentVector so this should be the first argument.
template<typename TSeed, typename TString, typename TDiff, typename TMatrix, typename TScoreString, typename TValue, typename TAlign, typename TScoreMatrix>
void
_calculateLastRectangle(TSeed const &seed,
						TDiff k,
						TMatrix &matrix_,
						TString *p_seq1,
						TString *p_seq2,
						TScoreString &score_str,
						TValue &score_length,
						::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > &alignmentVector,
						TScoreMatrix const & scoringScheme)
{
	SEQAN_CHECKPOINT;

	typedef typename ::std::map<TValue,Pair<TValue, TAlign> >::iterator TMapIterator;
	typedef Iter<TMatrix, PositionIterator> TMatrixIterator;

    TValue seq1_end = endPosition(*p_seq1);
	TValue seq2_end = endPosition(*p_seq2);
	TValue width_diag=getLowerDiagonal(seed) - getEndDiagonal(seed)+k;
	TValue width =  width_diag + seq1_end - (getEndDim0(seed) - 1);
	TValue height_diag = getEndDiagonal(seed) - getUpperDiagonal(seed)+k;
	TValue height = height_diag + seq2_end - (getEndDim1(seed) - 1);
	
	resize(score_str, width_diag+height_diag+1);
	
    typedef typename Infix<TString>::Type TSegment;
	TSegment seg1 = infix(host(*p_seq1),seq1_end-width+1,seq1_end);
	TSegment seg2 = infix(host(*p_seq2),seq2_end-height+1,seq2_end);

	_needleman_wunsch(matrix_, seg1, seg2, scoringScheme);

	TValue width_align = seq1_end - width  + width_diag +1;
	TValue height_align = seq2_end - height+1;

	TMatrixIterator iter_ = begin(matrix_);

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
		score_str[i] = *iter_;
		mapIt->second.i1 = -1;
		_needleman_wunsch_trace_lastRectangle(mapIt->second.i2, iter_, scoringScheme);
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
		score_str[i] = *iter_;
		mapIt->second.i1 = -1;
		_needleman_wunsch_trace_lastRectangle(mapIt->second.i2, iter_, scoringScheme);
		--x;
		--width_align;
	}
	score_length = length(score_str);
}

//calculation and backtracking of the edit matrix between two seeds
template<typename TSeed, typename TString, typename TDiff, typename TMatrix, typename TScoreString, typename TValue, typename TAlign, typename TScoreMatrix>
void
_calculateRectangle(TSeed const &seed,
					TSeed const &seed2,
					TDiff k_begin,
					TDiff k_end,
					TMatrix &matrix_,
					TString *p_seq1,
					TString *p_seq2,
					TScoreString &score_str,
					TValue &score_length,
					::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > &alignmentVector,
					TScoreMatrix const &scoringScheme)
{
	SEQAN_CHECKPOINT;

	typedef typename ::std::map<TValue,Pair<TValue, TAlign> >::iterator TMapIterator;
	typedef Iter<TMatrix, PositionIterator> TMatrixIterator;
    typedef typename Infix<TString>::Type TSegment;
	TSegment seg1b_align = infix(host(*p_seq1), (getEndDim0(seed2) - 1)-(getLowerDiagonal(seed2) - getEndDiagonal(seed2)   + k_end) + 1, getBeginDim0(seed) + getStartDiagonal(seed) - getUpperDiagonal(seed) + k_begin);
	TSegment seg2b_align = infix(host(*p_seq2), (getEndDim1(seed2) - 1)-(getEndDiagonal(seed2) -  getUpperDiagonal(seed2) + k_end) + 1, getBeginDim1(seed) + getLowerDiagonal(seed)  - getStartDiagonal(seed) + k_begin);

	_needleman_wunsch_rectangle(matrix_, seed, seed2, k_begin, k_end, seg1b_align, seg2b_align, scoringScheme,score_str);

	TValue width_diag = getLowerDiagonal(seed2) - getEndDiagonal(seed2)+k_end;
	TValue height_diag = getEndDiagonal(seed2) - getUpperDiagonal(seed2)+k_end;
	
	TValue w_d2 = getStartDiagonal(seed) - getUpperDiagonal(seed) + k_begin;
	TValue h_d2 = getLowerDiagonal(seed) - getStartDiagonal(seed) + k_begin;

	TValue overall = height_diag + width_diag + 1;

	resize(score_str, overall);
	TMatrixIterator matr_it = begin(matrix_);
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

		//cout << "ALL: " << length(seg1_align) << " " << length(seg2_align) << " " << width_stop << " " <<height_stop << endl;
		new_connect = _needleman_wunsch_trace_rectangle(mapIt->second.i2,  matr_it, scoringScheme, matrix_, width_stop, height_stop);
		score_str[j] = *matr_it;
		mapIt->second.i1 = new_connect;
		goNext(matr_it, 1);
		++height_align;
		--height_stop;
		if (old_connect != new_connect){
			//cout << "ups "<< old_connect << " " << new_connect << endl;
			_deleteAlignment(alignmentVector, old_connect, new_connect);
			//cout << "ups_ENDE" << endl;
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
		new_connect = _needleman_wunsch_trace_rectangle(mapIt->second.i2,  matr_it, scoringScheme, matrix_, width_stop, height_stop);
		score_str[j] = *matr_it;
		mapIt->second.i1 = new_connect;
		goPrevious(matr_it,0);
		--width_align;
		++width_stop;
		if (old_connect != new_connect)
			_deleteAlignment(alignmentVector, old_connect, new_connect);
		old_connect = new_connect;
	}
	_deleteAlignment(alignmentVector, old_connect, score_length);
	score_length = length(score_str);
}



//Rectangle calculation between two seeds
template <typename TScoreValue, unsigned DIMENSION, typename TSeedSpec, typename TSeedConfig, typename TString, typename TValue2>
void
_needleman_wunsch_rectangle(Matrix<TScoreValue, DIMENSION> & matrix_,			//edit matrix
							Seed<TSeedSpec, TSeedConfig> const &seed1,		//Seed nearer to the end
							Seed<TSeedSpec, TSeedConfig> const &seed2,		//Seed nearer to the start
							TValue2 k_begin,							//upper diagonal extension
							TValue2 k_end,								//lower diagonal extension
							TString const & str1_,						//first sequence
							TString const & str2_,						//secondSequence
							Score<TScoreValue, Simple> const & score_,	//score matrix
							String<TScoreValue> init)//Values for initialisation
{
	SEQAN_CHECKPOINT;

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
	TScoreValue score_gap = scoreGapExtend(score_);

	TSize str1_length = length(str1_);
	TSize str2_length = length(str2_);
	setLength(matrix_, 0, str1_length + 1);
	setLength(matrix_, 1, str2_length + 1);
	
	resize(matrix_);

	TMatrixIterator col_ = begin(matrix_);
	//Initialisierung

	setPosition(col_, (diag_height2 + rectangle_height+1)*(str1_length+1)-1);

	for (int i = 0; i != diag_width1; ++i)
	{
		
		*col_ = init[i];
		--col_;
	}

	TSize len = length(init);
	for (int i = diag_width1; i != len; ++i)
	{
		*col_=init[i];
		goNext(col_,1);
	}

	goPrevious(col_,1);
	TMatrixIterator finger2 =col_;
	
	--col_;



	TSize width_align = str1_length - diag_width1;
	for (int i = 0; i!= width_align; ++i)
	{
		*col_ = inf;
		--col_;
	}


	TStringIterator x_begin = begin(str1_);
	--x_begin;
	TStringIterator x_end = end(str1_) -1;
	setPosition(x_end,width_align-1);
	TStringIterator y_begin = begin(str2_);
	TStringIterator y_end = end(str2_)-1;
	TStringIterator x,y;

	col_ = finger2;
	TMatrixIterator finger1;
	
	TScoreValue h, v;
	TScoreValue s1;
	for (int i = 0; i < diag_height1; ++i)
	{
		finger2 = col_;		//points to last column
		goPrevious(col_, 1);	//points to this column
		finger1 = col_;
		v = *finger1;
		h = *finger2;
		for (x = x_end; x != x_begin; --x)
		{
			goPrevious(finger1, 0);
			goPrevious(finger2, 0);
			s1 = h + ((*x == *y_end) ? score_match : score_mismatch);
			h = *finger2;
			TScoreValue s2 = score_gap + ((h > v) ? h : v);
			v = (s1 > s2) ? s1 : s2;
			*finger1 = v;
		}
		--y_end;
	}

	setPosition(col_,(diag_height2 + rectangle_height+1)*(rectangle_width+1)-1);

	h =*col_;
	x_end = end(str1_) -1;
	for (int i = 0; i < rectangle_height; ++i)
	{
		v = inf;
		finger2 = col_;		//points to last column
		goPrevious(col_, 1);	//points to this column
		finger1 = col_;

		*finger1 = v;

		for (x = x_end; x != x_begin; --x)
		{
			goPrevious(finger1, 0);
			goPrevious(finger2, 0);
			s1 = h + ((*x == *y_end) ? score_match : score_mismatch);
			h = *finger2;
			TScoreValue s2 = score_gap + ((h > v) ? h : v);
			v = (s1 > s2) ? s1 : s2;
			//cout << position(finger1,0) << endl;
			*finger1 = v;
		}
		h = inf;
		--y_end;
	}

	setPosition(x_begin,diag_width2-1);

	for (int i = 0; i < diag_height2; ++i)
	{
		v = inf;
		finger2 = col_;		//points to last column
		goPrevious(col_, 1);	//points to this column
		finger1 = col_;

		*finger1 = v;

		for (x = x_end; x != x_begin; --x)
		{
			goPrevious(finger1, 0);
			goPrevious(finger2, 0);
			s1 = h + ((*x == *y_end) ? score_match : score_mismatch);
			h = *finger2;
			TScoreValue s2 = score_gap + ((h > v) ? h : v);
			v = (s1 > s2) ? s1 : s2;
			*finger1 = v;
			
		}
		h = inf;
		--y_end;
	}
}

//changed version of usual trace-back
template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, unsigned DIMENSION, typename TValue, typename TMatrix>
TScoreValue
_needleman_wunsch_trace_rectangle(Align<TTargetSource, TTargetSpec> & target_,
						Iter< Matrix<TScoreValue, DIMENSION>, PositionIterator > source_,
						Score<TScoreValue, Simple> const & score_,
						TMatrix matrix, 
						TValue width_stop, 
						TValue height_stop)
{
	SEQAN_CHECKPOINT;

	typedef Iter<Matrix<TScoreValue, DIMENSION>, PositionIterator > TMatrixIterator;
    typedef typename Infix<TTargetSource>::Type TTargetSourceSegment;

	TTargetSourceSegment str_0 = sourceSegment(row(target_, 0));
	TTargetSourceSegment str_1 = sourceSegment(row(target_, 1));
	typedef Align<TTargetSource, TTargetSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TTargetIterator;
	TTargetIterator target_0 = iter(row(target_, 0), 0);
	TTargetIterator target_1 = iter(row(target_, 1), 0);

	typedef typename Iterator<TTargetSourceSegment, Rooted>::Type TStringIterator;
	TStringIterator it_0 = iter(str_0, 0);
	TStringIterator it_0_end = end(str_0);

	TStringIterator it_1 = iter(str_1, 0);
	TStringIterator it_1_end = end(str_1);

	TScoreValue score_diff = scoreMismatch(score_) - scoreGapExtend(score_);
	//TScoreValue score_match = scoreMatch(score_);

	//-------------------------------------------------------------------------
	//follow the trace until the border is reached

	while ((static_cast<TValue>(position(it_0)) < width_stop) || (static_cast<TValue>(position(it_1)) < height_stop))
	{
		
		bool gv;
		bool gh;
		TMatrixIterator it3_ = source_;
		goNext(it3_,1);
		goNext(it3_,0);
		if ((*it_0 == *it_1)&&(*it3_ != -1000000))
		{
			gv = gh = true;
		}
		else
		{

			TMatrixIterator it_ = source_;

			goNext(it_, 0);
			TScoreValue v = *it_;

			goNext(it_, 1);
			TScoreValue d = *it_;

			it_ = source_;
			goNext(it_, 1);
			TScoreValue h = *it_;

			gv = (v >= h) | (d + score_diff >= h);
			gh = (h >  v) | (d + score_diff >= v);

		}
	
		if (gv)
		{
			++it_0;
			goNext(source_, 0);
		}
		else
		{
			insertGap(target_0);
		}

		if (gh) 
		{
			++it_1;
			goNext(source_, 1);
		}
		else
		{
			insertGap(target_1);
		}
		++target_0;
		++target_1;
	}

	setSourceEndPosition(row(target_,1),position(it_1));
	setSourceEndPosition(row(target_,0),position(it_0));

	return length(matrix,0) - coordinate(source_,0) + position(it_1) - height_stop -1;
}


//Rectangle calculation between two seeds
template <typename TScoreValue, unsigned DIMENSION, typename TSeedSpec, typename TSeedConfig, typename TString, typename TValue2>
void
_banded_needleman_wunsch_rectangle_first(Matrix<TScoreValue, DIMENSION> & matrix_,	//edit matrix
								   Seed<TSeedSpec, TSeedConfig> const &seed,		//Seed
								   TValue2 k,									//diagonal extension
								   TString const & str1_,						//first sequence
								   TString const & str2_,						//secondSequence
								   Score<TScoreValue, Simple> const & score_,	//score matrix
								   String<TScoreValue> init)					//Values for initialisation
{
	SEQAN_CHECKPOINT;

	typedef Matrix<TScoreValue, DIMENSION> TMatrix;

	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Iterator<TMatrix, Standard>::Type TMatrixIterator;

	typedef typename Iterator<TString const, Rooted>::Type TStringIterator;
	typedef typename Value<TString const>::Type TAlphabet;

	//-------------------------------------------------------------------------
	//define some variables

	TScoreValue inf = -1000000;

	TSize diag_height1 = getLowerDiagonal(seed) - getStartDiagonal(seed) + k;
	TSize diag_width1 = getStartDiagonal(seed) - getUpperDiagonal(seed) + k;

	TSize rectangle_height = getBeginDim1(seed) - beginPosition(str2_);
	TSize rectangle_width = length(str1_);
	


	TScoreValue score_match = scoreMatch(score_);
	TScoreValue score_mismatch = scoreMismatch(score_);
	TScoreValue score_gap = scoreGapExtend(score_);

	TSize str1_length = length(str1_);
	TSize str2_length = length(str2_);
	setLength(matrix_, 0, str1_length + 1);
	setLength(matrix_, 1, str2_length + 1);
	
	resize(matrix_);
	TMatrixIterator col_ = begin(matrix_);

	//Initialisierung
	setPosition(col_, (rectangle_height+1)*(rectangle_width+1)-1);

	for (int i = 0; i != diag_width1; ++i)
	{
		*col_ = init[i];
		--col_;
	}

	TSize len = length(init);
	for (int i = diag_width1; i != len; ++i)
	{
		*col_=init[i];
		goNext(col_,1);
	}

	goPrevious(col_,1);
	TMatrixIterator finger2 =col_;
	--col_;

	TSize width_align = str1_length - diag_width1;
	for (int i = 0; i!= width_align; ++i)
	{
		*col_ = inf;
		--col_;
	}

	//calculation of matrix

	TStringIterator x_begin = begin(str1_);
	--x_begin;
	TStringIterator x_end = end(str1_) -1;
	
	setPosition(x_end, width_align-1);
	TStringIterator y_begin = begin(str2_);
	TStringIterator y_end = end(str2_)-1;
	TStringIterator x,y;
	col_ = finger2;
	TMatrixIterator finger1;
	
	TScoreValue h, v;
	TScoreValue s1;
	for (int i = 0; i < diag_height1; ++i)
	{
		finger2 = col_;		//points to last column
		goPrevious(col_, 1);	//points to this column
		finger1 = col_;
		v = *finger1;
		h = *finger2;
		for (x = x_end; x != x_begin; --x)
		{
			goPrevious(finger1, 0);
			goPrevious(finger2, 0);
			s1 = h + ((*x == *y_end) ? score_match : score_mismatch);
			h = *finger2;
			TScoreValue s2 = score_gap + ((h > v) ? h : v);
			v = (s1 > s2) ? s1 : s2;
			*finger1 = v;
		}
		--y_end;
	}

	setPosition(col_, (rectangle_height+1)*(rectangle_width+1)-1);


	h =*col_;
	x_end = end(str1_) -1;
	for (int i = 0; i < rectangle_height; ++i)
	{
		v = inf;
		finger2 = col_;		//points to last column
		goPrevious(col_, 1);	//points to this column
		finger1 = col_;

		*finger1 = v;

		for (x = x_end; x != x_begin; --x)
		{
			goPrevious(finger1, 0);
			goPrevious(finger2, 0);
			s1 = h + ((*x == *y_end) ? score_match : score_mismatch);
			h = *finger2;
			TScoreValue s2 = score_gap + ((h > v) ? h : v);
			v = (s1 > s2) ? s1 : s2;
			*finger1 = v;
			
		}
		h = inf;
		--y_end;
	}	
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_CHAIN_BANDED_LINEAR_H_
