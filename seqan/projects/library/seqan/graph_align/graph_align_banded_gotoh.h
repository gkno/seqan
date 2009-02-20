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
  $Id: graph_align_banded_gotoh.h 1812 2008-03-31 15:54:55Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_BANDED_GOTOH_H
#define SEQAN_HEADER_GRAPH_ALIGN_BANDED_GOTOH_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Banded Gotoh Alignment, affine gap cost
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TRight, typename TSpec, typename TValue1, typename TIndex1, typename TValue2, typename TIndex2>
inline void
__processLastRow(AlignConfig<TTop, TLeft, TRight, false, TSpec> const,
				 TValue1& maxValue,
				 TIndex1& maxIndex,
				 TValue2 const val,
				 TIndex2 const row,
				 TIndex2 const col)
{
	SEQAN_CHECKPOINT
	maxValue.i1 = val;
	maxIndex.i1 = Pair<TIndex2,TIndex2>(row, col);
}

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TRight, typename TSpec, typename TValue1, typename TIndex1, typename TValue2, typename TIndex2>
inline void
__processLastRow(AlignConfig<TTop, TLeft, TRight, true, TSpec> const,
				 TValue1& maxValue,
				 TIndex1& maxIndex,
				 TValue2 const val,
				 TIndex2 const row,
				 TIndex2 const col)
{
	SEQAN_CHECKPOINT
	if (val > maxValue.i1) {
		maxValue.i1 = val;
		maxIndex.i1 = Pair<TIndex2,TIndex2>(row, col);
	}
}



//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TBottom, typename TSpec, typename TValue1, typename TIndex1, typename TValue2, typename TIndex2>
inline void
__processLastColumn(AlignConfig<TTop, TLeft, false, TBottom, TSpec> const,
					TValue1& maxValue,
					TIndex1& maxIndex,
					TValue2 const val,	
					TIndex2 const row,
					TIndex2 const col)
{
	SEQAN_CHECKPOINT
	maxValue.i2 = val;
	maxIndex.i2 = Pair<TIndex2,TIndex2>(row, col);
}

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TBottom, typename TSpec, typename TValue1, typename TIndex1, typename TValue2, typename TIndex2>
inline void
__processLastColumn(AlignConfig<TTop, TLeft, true, TBottom, TSpec> const,
					TValue1& maxValue,
					TIndex1& maxIndex,
					TValue2 const val,
					TIndex2 const row,
					TIndex2 const col)
{
	SEQAN_CHECKPOINT
	if (val > maxValue.i2) {
		maxValue.i2 = val;
		maxIndex.i2 = Pair<TIndex2,TIndex2>(row, col);
	}
}



//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TStringSet, typename TScore, typename TColumn, typename TValPair, typename TIndexPair, typename TDiagonal>
inline void
_align_banded_gotoh_trace(TAlign& align,
						  TStringSet const& str,
						  TScore const& sc,
						  TColumn const& mat,
						  TColumn const& horizontal,
						  TColumn const& vertical,
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
	TString const& str1 = str[0];
	TString const& str2 = str[1];	
	TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
	TId id2 = positionToId(const_cast<TStringSet&>(str), 1);
	TSize len1 = length(str1) + 1;
	TSize len2 = length(str2) + 1;
	TScoreValue gapOpen = scoreGapOpen(sc);
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
	TScoreValue horiOrigMax = value(horizontal, row * diagonalWidth + col);
	TScoreValue vertMax = value(vertical, row * diagonalWidth + col);
	if (horiOrigMax == value(mat, row * diagonalWidth + col)) tv = Horizontal;
	else if (vertMax == value(mat, row * diagonalWidth + col)) tv = Vertical;
	else tv = Diagonal;
	
	// Walk until we hit a border
	TTraceValue oldTraceValue = tv;
	TSize seqLen = 0;
	bool hitBorder = false;
	do {
		actualRow = row + lo_row;
		actualCol = col + diagL + actualRow;

		// Direction changed, so make aligned segments
		if (oldTraceValue != tv) {
			_align_trace_print(align, str, id1, actualCol, id2, actualRow, seqLen, oldTraceValue);
			seqLen = 0;
		}

		// Check if we hit a border
		if ((actualRow == 0) || (actualCol == 0)) hitBorder = true;
		else {
			//std::cout << row << ',' << col << ':' << value(originalMat, actualRow * len1 + actualCol) << std::endl; 
			
			// Last value was diagonal
			if (tv == Diagonal) {
				oldTraceValue = Diagonal;
				if ((value(horizontal, row * diagonalWidth + col)) == value(mat, row * diagonalWidth + col)) {
					tv = Horizontal;
				} else if ((value(vertical, row * diagonalWidth + col)) == value(mat, row * diagonalWidth + col)) {
					tv = Vertical;
				} else {
					--row; ++seqLen;
				}
			} else if (tv == Horizontal) { // Last value was horizontal
				oldTraceValue = Horizontal;
				if ((value(mat, row * diagonalWidth + (col - 1)) + gapOpen) == value(horizontal, row * diagonalWidth + col)) {
					tv = Diagonal;
				}
				--col; ++seqLen;
			} else { // Vertical
				oldTraceValue = Vertical;
				if ((value(mat, (row - 1) * diagonalWidth + (col + 1)) + gapOpen) == value(vertical, row * diagonalWidth + col)) {
					tv = Diagonal;
				}
				--row; ++col; ++seqLen;
			}	
		}
	} while(!hitBorder);
	
	// Align left overs
	if (seqLen) _align_trace_print(align, str, id1, actualCol, id2, actualRow, seqLen, tv);

	// Handle the remaining sequence
	if (actualCol != 0) _align_trace_print(align, str, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) actualCol,  Horizontal);
	else if (actualRow != 0) _align_trace_print(align, str, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) actualRow,  Vertical);

}


//////////////////////////////////////////////////////////////////////////////

template <typename TColumn, typename TStringSet, typename TScore, typename TValPair, typename TIndexPair, typename TDiagonal, typename TAlignConfig>
inline typename Value<TScore>::Type
_align_banded_gotoh(TColumn& mat,
					TColumn& horizontal,
					TColumn& vertical,
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
	TScoreValue gap = scoreGapExtend(sc);
	TScoreValue gapOpen = scoreGapOpen(sc);
	TScoreValue tmp = 0;
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
	TScoreValue infValue = -1 * _getInfinityDistance<TScoreValue>();
	fill(mat, height * diagonalWidth, infValue);
	fill(horizontal, height * diagonalWidth, infValue);
	fill(vertical, height * diagonalWidth, infValue);
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
	for(TSize row = 0; row < height; ++row) {
		TSize actualRow = row + lo_row;
		if (lo_diag > 0) --lo_diag;
		if (row + lo_row >= len1 - diagU) --hi_diag;
		for(TSize col = lo_diag; col<hi_diag; ++col) {
			TSize actualCol = col + diagL + actualRow;
			//std::cout << row << ',' << col << ':' << value(originalMat, actualRow * len1 + actualCol) << std::endl;

			// Usual initialization for first row and column
			if (actualRow == 0) {
				if (actualCol != 0) {
					_initFirstRow(TAlignConfig(), value(mat, col), gapOpen + (actualCol - 1) * gap);
					value(vertical, col) = value(mat, col) + gapOpen - gap;
				} else value(mat, col) = 0;
				//std::cout << row << ',' << col << ':' << value(mat, row * diagonalWidth + col) << std::endl;
				//std::cout << row << ',' << col << ':' << value(horizontal, row * diagonalWidth + col) << std::endl;
				//std::cout << row << ',' << col << ':' << value(vertical, row * diagonalWidth + col) << std::endl;
				continue;
			} else if (actualCol == 0) {
				_initFirstColumn(TAlignConfig(), value(mat, row * diagonalWidth + col), gapOpen + (actualRow - 1) * gap);
				value(horizontal, row * diagonalWidth + col) = value(mat, row * diagonalWidth + col) + gapOpen - gap;
				//std::cout << row << ',' << col << ':' << value(mat, row * diagonalWidth + col) << std::endl;
				//std::cout << row << ',' << col << ':' << value(horizontal, row * diagonalWidth + col) << std::endl;
				//std::cout << row << ',' << col << ':' << value(vertical, row * diagonalWidth + col) << std::endl;
				continue;
			}

			// Get the new maximum for vertical
			if (col < diagonalWidth - 1) {
				if ((tmp = value(mat, (row - 1) * diagonalWidth + (col + 1)) + gapOpen) >  (value(vertical, (row - 1) * diagonalWidth + (col + 1)) + gap)) {
					value(vertical, row * diagonalWidth + col) = tmp;
				} else {
					value(vertical, row * diagonalWidth + col) = (value(vertical, (row - 1) * diagonalWidth + (col + 1)) + gap);
				}
			}

			// Get the new maximum for horizontal
			if (col > 0) {
				if ((tmp = value(mat, row * diagonalWidth + (col - 1)) + gapOpen) > value(horizontal, row * diagonalWidth + (col - 1)) + gap) {
					value(horizontal, row * diagonalWidth + col) = tmp;
				} else {
					value(horizontal, row * diagonalWidth + col) = value(horizontal, row * diagonalWidth + (col - 1)) + gap;
				}
			}

			// Get the new maximum for mat
			TScoreValue sc_ = score(const_cast<TScore&>(sc), actualCol-1, actualRow-1, str1, str2);
			tmp = value(mat, (row - 1) * diagonalWidth + col) + sc_;
			if (value(vertical, row * diagonalWidth + col) > tmp) {
				tmp = value(vertical, row * diagonalWidth + col);
			}
			if (value(horizontal, row * diagonalWidth + col) > tmp) {
				tmp = value(horizontal, row * diagonalWidth + col);
			}
			value(mat, row * diagonalWidth + col) = tmp;

			// Store the maximum
			if (actualCol == len1 - 1) __processLastColumn(TAlignConfig(), overallMaxValue, overallMaxIndex, value(mat, row * diagonalWidth + col), row, col);
			if (actualRow == len2 - 1) __processLastRow(TAlignConfig(), overallMaxValue, overallMaxIndex, value(mat, row * diagonalWidth + col), row, col);
			//std::cout << row << ',' << col << ':' << value(mat, row * diagonalWidth + col) << std::endl;
			//std::cout << row << ',' << col << ':' << value(horizontal, row * diagonalWidth + col) << std::endl;
			//std::cout << row << ',' << col << ':' << value(vertical, row * diagonalWidth + col) << std::endl;
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
				 BandedGotoh)
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
	TColumn mat; TColumn horizontal; TColumn vertical;
	TScoreValue maxScore = _align_banded_gotoh(mat, horizontal, vertical, str, sc, overallMaxValue, overallMaxIndex, (int) diag1, (int) diag2, TAlignConfig());
	
	// Follow the trace and create the graph
	_align_banded_gotoh_trace(align, str, sc, mat, horizontal, vertical, overallMaxValue, overallMaxIndex, (int) diag1, (int) diag2);

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
				 BandedGotoh)
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
	TColumn mat; TColumn horizontal; TColumn vertical;
	return _align_banded_gotoh(mat, horizontal, vertical, str, sc, overallMaxValue, overallMaxIndex, (int) diag1, (int) diag2, TAlignConfig());
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
