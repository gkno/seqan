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

template<bool TTop, bool TLeft, bool TBottom, typename TSpec, typename TValue1, typename TIndex1, typename TValue2, typename TIndex2>
inline void
_lastColumn(AlignConfig<TTop, TLeft, false, TBottom, TSpec> const,
			TValue1& maxValue,
			TIndex1& maxIndex,
			TValue2 const val,	
			TIndex2 const row,
			TIndex2 const col)
{
	SEQAN_CHECKPOINT
	maxValue[1] = val; maxIndex[2] = row; maxIndex[3] = col;
}

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TBottom, typename TSpec, typename TValue1, typename TIndex1, typename TValue2, typename TIndex2>
inline void
_lastColumn(AlignConfig<TTop, TLeft, true, TBottom, TSpec> const,
			TValue1& maxValue,
			TIndex1& maxIndex,
			TValue2 const val,
			TIndex2 const row,
			TIndex2 const col)
{
	SEQAN_CHECKPOINT
	if (val > maxValue[1]) {maxValue[1] = val; maxIndex[2] = row; maxIndex[3] = col; }
}

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TRight, typename TSpec, typename TValue1, typename TIndex1, typename TValue2, typename TIndex2>
inline void
_lastRow(AlignConfig<TTop, TLeft, TRight, false, TSpec> const,		
		 TValue1& maxValue,
		 TIndex1& maxIndex,
		 TValue2 const val,
		 TIndex2 const row,
		 TIndex2 const col)
{
	SEQAN_CHECKPOINT
	maxValue[0] = val; maxIndex[0] = row; maxIndex[1] = col;
}

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TRight, typename TSpec, typename TValue1, typename TIndex1, typename TValue2, typename TIndex2>
inline void
_lastRow(AlignConfig<TTop, TLeft, TRight, true, TSpec> const,
		 TValue1& maxValue,
		 TIndex1& maxIndex,
		 TValue2 const val,
		 TIndex2 const row,
		 TIndex2 const col)
{
	SEQAN_CHECKPOINT
	if (val > maxValue[0]) {maxValue[0] = val; maxIndex[0] = row; maxIndex[1] = col; }
}



//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TStringSet, typename TTrace, typename TValPair, typename TIndexPair, typename TDiagonal>
inline void
_align_banded_nw_trace(TAlign& align,
					   TStringSet const& str,
					   TTrace const& trace,
					   TValPair const& overallMaxValue,
					   TIndexPair const& overallMaxIndex,
					   TDiagonal const diagL,
					   TDiagonal const diagU)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Id<TStringSet>::Type TId;
	typedef typename Size<TTrace>::Type TSize;
	typedef typename Value<TTrace>::Type TTraceValue;

	// Traceback values
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;

	// Initialization
	TString const& str1 = str[0];
	TString const& str2 = str[1];	
	TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
	TId id2 = positionToId(const_cast<TStringSet&>(str), 1);
	TSize len1 = length(str1) + 1;
	TSize len2 = length(str2) + 1;
	TSize lo_row = (diagU <= 0) ? -1 * diagU : 0;
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
	TSize row = (overallMaxValue[0] > overallMaxValue[1]) ? overallMaxIndex[0] : overallMaxIndex[2];
	TSize col = (overallMaxValue[0] > overallMaxValue[1]) ? overallMaxIndex[1] : overallMaxIndex[3];

	// Handle the skipped sequence parts
	TSize actualRow = row + lo_row;
	TSize actualCol = col + diagL + actualRow;
	if (actualCol + 1 < len1) _align_trace_print(align, str, id1, actualCol, id2, actualRow, (len1 - (actualCol + 1)),  Horizontal);
	if (actualRow + 1 < len2) _align_trace_print(align, str, id1, actualCol, id2, actualRow, (len2 - (actualRow + 1)),  Vertical);

	// Find initial direction
	TTraceValue tv = trace[row * diagonalWidth + col];
	if (tv == Horizontal) --col;
	else if (tv == Vertical) {--row; ++col;} 
	else --row;
	
	// Walk until we hit a border
	TSize seqLen = 1;
	TTraceValue newTv = tv;
	while(true) {
		actualRow = row + lo_row;
		actualCol = col + diagL + actualRow;
		newTv = trace[row * diagonalWidth + col];

		// Check if we hit a border
		if ((actualRow == 0) || (actualCol == 0)) break;
		else {
			//std::cout << row << ',' << col << ':' << value(originalMat, actualRow * len1 + actualCol) << std::endl; 
			if (tv == Diagonal) {
				if (newTv == Horizontal) {
					_align_trace_print(align, str, id1, actualCol, id2, actualRow, seqLen, tv);
					--col; seqLen = 1;
				} else if (newTv == Vertical) {
					_align_trace_print(align, str, id1, actualCol, id2, actualRow, seqLen, tv);
					--row; ++col; seqLen = 1;
				} else {
					--row; ++seqLen;
				}
			} else {
				if (tv == Horizontal) { 
					if (newTv == Diagonal) {
						_align_trace_print(align, str, id1, actualCol, id2, actualRow, seqLen, tv);
						--row; seqLen = 1;
					} else {
						--col; ++seqLen;
					}
				} else { 
					if (newTv == Diagonal) {
						_align_trace_print(align, str, id1, actualCol, id2, actualRow, seqLen, tv);
						--row; seqLen = 1;
					} else {
						--row; ++col; ++seqLen;
					}
				}
			}
			tv = newTv;
		}
	}
	
	// Align left overs
	if (seqLen) _align_trace_print(align, str, id1, actualCol, id2, actualRow, seqLen, tv);

	// Handle the remaining sequence
	if (actualCol != 0) _align_trace_print(align, str, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) actualCol,  Horizontal);
	else if (actualRow != 0) _align_trace_print(align, str, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) actualRow,  Vertical);

}

////////////////////////////////////////////////////////////////////////////


template <typename TColumn, typename TTrace, typename TStringSet, typename TScore, typename TValPair, typename TIndexPair, typename TDiagonal, typename TAlignConfig>
inline typename Value<TScore>::Type
_align_banded_nw(TColumn& mat,
				 TTrace& trace,
				 TStringSet const& str,
				 TScore const & sc,
				 TValPair& overallMaxValue,
				 TIndexPair& overallMaxIndex,
				 TDiagonal diagL,
				 TDiagonal diagU,
				 TAlignConfig const)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Size<TColumn>::Type TSize;

	// Initialization
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;
	TString const& str1 = str[0];
	TString const& str2 = str[1];		
	TSize len1 = length(str1) + 1;
	TSize len2 = length(str2) + 1;
	TSize diagonalWidth = (TSize) (diagU - diagL + 1);
	TSize hi_diag = diagonalWidth;
	TSize lo_diag = 0;
	if (diagL > 0) lo_diag = 0;
	else lo_diag = (diagU < 0) ? hi_diag : (TSize) (1-1 * diagL); 
	TSize lo_row = (diagU <= 0) ? -1 * diagU : 0;
	TSize hi_row = len2;
	if (len1 - diagL < hi_row) hi_row = len1 - diagL;
	TSize height = hi_row - lo_row;
	resize(mat, height * diagonalWidth);
	resize(trace, height * diagonalWidth);
	overallMaxValue[0] = InfimumValue<TScoreValue>::VALUE;
	overallMaxValue[1] = InfimumValue<TScoreValue>::VALUE;
	overallMaxIndex[0] = diagonalWidth;	overallMaxIndex[1] = height;
	overallMaxIndex[2] = diagonalWidth;	overallMaxIndex[3] = height;
	
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
	typedef typename Iterator<TColumn, Standard>::Type TColIter;
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	TSize actualCol = 0;
	TSize actualRow = 0;
	TScoreValue verti_val = 0;
	TScoreValue hori_val = 0;
	for(TSize row = 0; row < height; ++row) {
		actualRow = row + lo_row;
		if (lo_diag > 0) --lo_diag;
		if (row + lo_row >= len1 - diagU) --hi_diag;
		TTraceIter traceIt = begin(trace, Standard()) + row * diagonalWidth + lo_diag;
		TColIter matIt = begin(mat, Standard()) + row * diagonalWidth + lo_diag;
		TColIter matItVerti = begin(mat, Standard()) + (row - 1) * diagonalWidth + lo_diag + 1;
		for(TSize col = lo_diag; col<hi_diag; ++col, ++matIt, ++matItVerti, ++traceIt) {
			actualCol = col + diagL + actualRow;
			//std::cout << row << ',' << col << ':' << value(originalMat, actualRow * len1 + actualCol) << std::endl;

			if ((actualRow != 0) && (actualCol != 0)) {
				// Get the new maximum for mat
				*matIt = *(matItVerti - 1) + score(const_cast<TScore&>(sc), ((int) actualCol - 1), ((int) actualRow - 1), str1, str2);
				*traceIt = Diagonal;
				if ((verti_val = (col < diagonalWidth - 1) ?	*matItVerti + scoreGapExtendVertical(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2) : InfimumValue<TScoreValue>::VALUE) > *matIt) {
					*matIt = verti_val;
					*traceIt = Vertical;
				}						
				if ((hori_val = (col > 0) ? *(matIt - 1) + scoreGapExtendHorizontal(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2) : InfimumValue<TScoreValue>::VALUE) > *matIt) {
					*matIt = hori_val;
					*traceIt = Horizontal;
				}
			
				// Store the maximum
				if (actualCol == len1 - 1) _lastColumn(TAlignConfig(), overallMaxValue, overallMaxIndex, *matIt, row, col);
				if (actualRow == len2 - 1) _lastRow(TAlignConfig(), overallMaxValue, overallMaxIndex, *matIt, row, col);
				//std::cout << row << ',' << col << ':' << value(mat, row * diagonalWidth + col) << std::endl;
			
			} else {
				
				// Usual initialization for first row and column
				if (actualRow == 0) _initFirstRow(TAlignConfig(), *matIt, (TScoreValue) actualCol * scoreGapExtendHorizontal(sc, ((int) actualCol - 1), -1, str1, str2));
				else _initFirstColumn(TAlignConfig(), *matIt, (TScoreValue) actualRow * scoreGapExtendVertical(sc, -1, ((int) actualRow - 1), str1, str2));			
			}
		}
	}
	return (overallMaxValue[0] > overallMaxValue[1]) ? overallMaxValue[0] : overallMaxValue[1];
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
	TScoreValue overallMaxValue[2];
	TSize overallMaxIndex[4];
	
	// Create the trace
	typedef typename Value<TScore>::Type TScoreValue;
	typedef String<TScoreValue> TColumn;
	TColumn mat;
	String<TraceBack> trace;
	TScoreValue maxScore = _align_banded_nw(mat, trace, str, sc, overallMaxValue, overallMaxIndex, (int) diag1, (int) diag2, TAlignConfig());
	
	// Follow the trace and create the graph
	_align_banded_nw_trace(align, str, trace, overallMaxValue, overallMaxIndex, (int) diag1, (int) diag2);

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
	TScoreValue overallMaxValue[2];
	TSize overallMaxIndex[4];
	
	// Calculate the score
	typedef typename Value<TScore>::Type TScoreValue;
	typedef String<TScoreValue> TColumn;
	TColumn mat;
	return _align_banded_nw(mat, str, sc, overallMaxValue, overallMaxIndex, (int) diag1, (int) diag2, TAlignConfig());
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
