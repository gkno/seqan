 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: graph_align_banded_needleman_wunsch.h 1812 2008-03-31 15:54:55Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_BANDED_NEEDLEMAN_WUNSCH_H
#define SEQAN_HEADER_GRAPH_ALIGN_BANDED_NEEDLEMAN_WUNSCH_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TStringSet, typename TScore, typename TColumn, typename TValPair, typename TIndexPair, typename TDiagonal>
inline void
_align_banded_nw_trace(TAlign& align,
					   TStringSet const& str,
					   TScore const& sc,
					   TColumn const& mat,
					   TValPair const& overallMaxValue,
					   TIndexPair const& overallMaxIndex,
					   TDiagonal const diagL,
					   TDiagonal const diagU)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Id<TStringSet>::Type TId;
	typedef typename Size<TColumn>::Type TSize;
	typedef unsigned char TTraceValue;

	// Gotoh back-trace values
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;

	// Initialization
	TScoreValue infValue = InfimumValue<TScoreValue>::VALUE;
	TString const& str1 = str[0];
	TString const& str2 = str[1];	
	TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
	TId id2 = positionToId(const_cast<TStringSet&>(str), 1);
	TSize len1 = length(str1) + 1;
	TSize len2 = length(str2) + 1;
	TSize lo_row = 0;
	if (diagU <= 0) lo_row = -1 * diagU;
	TSize diagonalWidth = (TSize) (diagU - diagL + 1);
	
	//// Debug stuff
	//TColumn originalMat;
	//resize(originalMat, len1 * len2);
	//TSize count = 0;
	//for(TSize i=0; i<len2; ++i) {
	//	for(TSize j=0; j<len1; ++j) {
	//		value(originalMat, i * len1 + j) = count;
	//		std::cout << count << ',';
	//		++count;
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;

	// Start the trace from the cell with the max value
	TSize row = 0;
	TSize col = 0;
	if (overallMaxValue.i1 > overallMaxValue.i2) {
		row = overallMaxIndex.i1.i1;
		col = overallMaxIndex.i1.i2;
	} else {
		row = overallMaxIndex.i2.i1;
		col = overallMaxIndex.i2.i2;
	}

	// Handle the skipped sequence parts
	TSize actualRow = row + lo_row;
	TSize actualCol = col + diagL + actualRow;
	if (actualCol + 1 < len1) _align_trace_print(align, str, id1, actualCol, id2, actualRow, (len1 - (actualCol + 1)),  Horizontal);
	if (actualRow + 1 < len2) _align_trace_print(align, str, id1, actualCol, id2, actualRow, (len2 - (actualRow + 1)),  Vertical);

	// Find initial direction
	TTraceValue tv = Diagonal;
	TScoreValue diag_val = mat[row * diagonalWidth + col];
	TScoreValue verti_val = (col < diagonalWidth - 1) ?	mat[(row - 1) * diagonalWidth + (col + 1)] + scoreGapExtendVertical(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2) : infValue;
	TScoreValue hori_val = (col > 0) ? mat[row * diagonalWidth + (col - 1)] + scoreGapExtendHorizontal(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2) : infValue;
	if (hori_val == diag_val) {
		tv = Horizontal;
		--col;
	} else if (verti_val == diag_val) {
		tv = Vertical;
		--row; ++col;
	} else {
		--row;
	}
	
	// Walk until we hit a border
	TTraceValue oldTraceValue = tv;
	TSize seqLen = 1;
	while(true) {
		actualRow = row + lo_row;
		actualCol = col + diagL + actualRow;

		diag_val = mat[row * diagonalWidth + col];

		/*
		// Direction changed, so make aligned segments
		if (oldTraceValue != tv) {
			_align_trace_print(align, str, id1, actualCol, id2, actualRow, seqLen, oldTraceValue);
			seqLen = 0;
		}
		*/

		// Check if we hit a border
		if ((actualRow == 0) || (actualCol == 0)) break;
		else {
			//std::cout << row << ',' << col << ':' << value(originalMat, actualRow * len1 + actualCol) << std::endl; 
			
			// Last value was diagonal
			oldTraceValue = tv;
			if (tv == Diagonal) {
				verti_val = (col < diagonalWidth - 1) ?	mat[(row - 1) * diagonalWidth + (col + 1)] + scoreGapExtendVertical(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2) : infValue;
				hori_val = (col > 0) ? mat[row * diagonalWidth + (col - 1)] + scoreGapExtendHorizontal(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2) : infValue;
				if (hori_val == mat[row * diagonalWidth + col]) {
					tv = Horizontal;
					_align_trace_print(align, str, id1, actualCol, id2, actualRow, seqLen, oldTraceValue);
					--col; seqLen = 1;
				} else if (verti_val == mat[row * diagonalWidth + col]) {
					tv = Vertical;
					_align_trace_print(align, str, id1, actualCol, id2, actualRow, seqLen, oldTraceValue);
					--row; ++col; seqLen = 1;
				} else {
					--row; ++seqLen;
				}
			} else {
				diag_val = mat[(row - 1) * diagonalWidth + col] + score(const_cast<TScore&>(sc), ((int) actualCol - 1), ((int) actualRow - 1), str1, str2);
				if (tv == Horizontal) { // Last value was horizontal
					if (diag_val == mat[row * diagonalWidth + col]) {
						tv = Diagonal;
						_align_trace_print(align, str, id1, actualCol, id2, actualRow, seqLen, oldTraceValue);
						--row; seqLen = 1;
					} else {
						--col; ++seqLen;
					}
				} else { // Vertical
					if (diag_val == mat[row * diagonalWidth + col]) {
						tv = Diagonal;
						_align_trace_print(align, str, id1, actualCol, id2, actualRow, seqLen, oldTraceValue);
						--row; seqLen = 1;
					} else {
						--row; ++col; ++seqLen;
					}
				}
			}
		}
	}
	
	// Align left overs
	if (seqLen) _align_trace_print(align, str, id1, actualCol, id2, actualRow, seqLen, tv);

	// Handle the remaining sequence
	if (actualCol != 0) _align_trace_print(align, str, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) actualCol,  Horizontal);
	else if (actualRow != 0) _align_trace_print(align, str, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) actualRow,  Vertical);

}

////////////////////////////////////////////////////////////////////////////


template <typename TColumn, typename TStringSet, typename TScore, typename TValPair, typename TIndexPair, typename TDiagonal, typename TAlignConfig>
inline typename Value<TScore>::Type
_align_banded_nw(TColumn& mat,
				 TStringSet const& str,
				 TScore const & sc,
				 TValPair& overallMaxValue,
				 TIndexPair& overallMaxIndex,
				 TDiagonal diagL,
				 TDiagonal diagU,
				 TAlignConfig const)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Size<TColumn>::Type TSize;

	// Initialization
	TString const& str1 = str[0];
	TString const& str2 = str[1];		
	TSize len1 = length(str1) + 1;
	TSize len2 = length(str2) + 1;
	TSize diagonalWidth = (TSize) (diagU - diagL + 1);
	TSize hi_diag = diagonalWidth;
	TSize lo_diag = 0;
	if (diagL > 0) lo_diag = 0;
	else if (diagU < 0) lo_diag = hi_diag;
	else lo_diag = (TSize) (1-1 * diagL); 
	TSize lo_row = 0;
	if (diagU <= 0) lo_row = -1 * diagU;
	TSize hi_row = len2;
	if (len1 - diagL < hi_row) hi_row = len1 - diagL;
	TSize height = hi_row - lo_row;
	TScoreValue infValue = InfimumValue<TScoreValue>::VALUE;
	fill(mat, height * diagonalWidth, infValue);
	overallMaxValue = Pair<TScoreValue, TScoreValue>(infValue, infValue);
	overallMaxIndex = Pair<Pair<TSize,TSize>, Pair<TSize,TSize> >(Pair<TSize,TSize>(diagonalWidth, height), Pair<TSize,TSize>(diagonalWidth, height));
	
	//// Debug stuff
	//TColumn originalMat;
	//resize(originalMat, len1 * len2);
	//TSize count = 0;
	//for(TSize i=0; i<len2; ++i) {
	//	for(TSize j=0; j<len1; ++j) {
	//		value(originalMat, i * len1 + j) = count;
	//		std::cout << count << ',';
	//		++count;
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;

	// Classical DP with affine gap costs
	TSize actualCol = 0;
	TScoreValue verti_val = 0;
	TScoreValue hori_val = 0;
	for(TSize row = 0; row < height; ++row) {
		TSize actualRow = row + lo_row;
		if (lo_diag > 0) --lo_diag;
		if (row + lo_row >= len1 - diagU) --hi_diag;
		for(TSize col = lo_diag; col<hi_diag; ++col) {
			actualCol = col + diagL + actualRow;
			//std::cout << row << ',' << col << ':' << value(originalMat, actualRow * len1 + actualCol) << std::endl;

			// Usual initialization for first row and column
			if (actualRow == 0) {
				_initFirstRow(TAlignConfig(), value(mat, col), (TScoreValue) actualCol * scoreGapExtendHorizontal(sc, ((int) actualCol - 1), -1, str1, str2));
				//std::cout << row << ',' << col << ':' << value(mat, row * diagonalWidth + col) << std::endl;
				continue;
			} else if (actualCol == 0) {
				_initFirstColumn(TAlignConfig(), value(mat, row * diagonalWidth + col), (TScoreValue) actualRow * scoreGapExtendVertical(sc, -1, ((int) actualRow - 1), str1, str2));
				//std::cout << row << ',' << col << ':' << value(mat, row * diagonalWidth + col) << std::endl;
				continue;
			}

			// Get the new maximum for vertical
			verti_val = (col < diagonalWidth - 1) ?	mat[(row - 1) * diagonalWidth + (col + 1)] + scoreGapExtendVertical(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2) : infValue;

			// Get the new maximum for horizontal
			hori_val = (col > 0) ? mat[row * diagonalWidth + (col - 1)] + scoreGapExtendHorizontal(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2) : infValue;

			// Get the new maximum for mat
			TScoreValue& maxVal = mat[row * diagonalWidth + col];
			maxVal = mat[(row - 1) * diagonalWidth + col] + score(const_cast<TScore&>(sc), ((int) actualCol - 1), ((int) actualRow - 1), str1, str2);
			if (verti_val > maxVal) 
				maxVal = verti_val;
			if (hori_val > maxVal)
				maxVal = hori_val;
			
			// Store the maximum
			if (actualCol == len1 - 1) __processLastColumn(TAlignConfig(), overallMaxValue, overallMaxIndex, value(mat, row * diagonalWidth + col), row, col);
			if (actualRow == len2 - 1) __processLastRow(TAlignConfig(), overallMaxValue, overallMaxIndex, value(mat, row * diagonalWidth + col), row, col);
			//std::cout << row << ',' << col << ':' << value(mat, row * diagonalWidth + col) << std::endl;
		}
	}
	if (overallMaxValue.i1 > overallMaxValue.i2) return overallMaxValue.i1;
	else return overallMaxValue.i2;
}

////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScore, typename TAlignConfig, typename TDiagonal>
inline typename Value<TScore>::Type
_globalAlignment(TAlign& align,
				 TStringSet const& str,
				 TScore const& sc,
				 TAlignConfig const,
				 TDiagonal diag1,
				 TDiagonal diag2,
				 BandedNeedlemanWunsch)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	  
	// Maximum value
	Pair<TScoreValue, TScoreValue> overallMaxValue;
	Pair<Pair<TSize,TSize>, Pair<TSize,TSize> > overallMaxIndex;

	// Create the trace
	typedef typename Value<TScore>::Type TScoreValue;
	typedef String<TScoreValue> TColumn;
	TColumn mat;
	TScoreValue maxScore = _align_banded_nw(mat, str, sc, overallMaxValue, overallMaxIndex, (int) diag1, (int) diag2, TAlignConfig());
	
	// Follow the trace and create the graph
	_align_banded_nw_trace(align, str, sc, mat, overallMaxValue, overallMaxIndex, (int) diag1, (int) diag2);

	return maxScore;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScore, typename TAlignConfig, typename TDiagonal>
inline typename Value<TScore>::Type
_globalAlignment(TStringSet const& str,
				 TScore const& sc,
				 TAlignConfig const,
				 TDiagonal diag1,
				 TDiagonal diag2,
				 BandedNeedlemanWunsch)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	  
	// Maximum value
	Pair<TScoreValue, TScoreValue> overallMaxValue;
	Pair<Pair<TSize,TSize>, Pair<TSize,TSize> > overallMaxIndex;

	// Calculate the score
	typedef typename Value<TScore>::Type TScoreValue;
	typedef String<TScoreValue> TColumn;
	TColumn mat;
	return _align_banded_nw(mat, str, sc, overallMaxValue, overallMaxIndex, (int) diag1, (int) diag2, TAlignConfig());
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
