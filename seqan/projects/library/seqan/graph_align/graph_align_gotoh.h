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
  $Id: graph_align_gotoh.h 1812 2008-03-31 15:54:55Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_GOTOH_H
#define SEQAN_HEADER_GRAPH_ALIGN_GOTOH_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Gotoh Alignment, affine gap cost
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TStringSet, typename TTrace, typename TIndexPair, typename TVal>
inline void
_align_gotoh_trace(TAlign& align,		 
				   TStringSet const& str,
				   TTrace const& trace,
				   TIndexPair const& overallMaxIndex,
				   TVal const initialDir)
{
	SEQAN_CHECKPOINT
	
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Id<TStringSet>::Type TId;

	// TraceBack values for Gotoh
	TTraceValue Diagonal = 0;
	TTraceValue Horizontal = 1;
	TTraceValue Vertical = 2;

	TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
	TId id2 = positionToId(const_cast<TStringSet&>(str), 1);
	TSize len1 = overallMaxIndex[0];
	TSize len2 = overallMaxIndex[1];
	if (len1 < length(str[0])) {
		_align_trace_print(align, str, id1, len1, id2, len2, length(str[0]) - len1,  Horizontal);
	} else if (len2 < length(str[1])) {
		_align_trace_print(align, str, id1, len1, id2, len2, length(str[1]) - len2,  Vertical);
	}
	TSize numRows = length(str[1]);

	// Initialize everything	
	TTraceValue nextTraceValue = value(trace, (len1 - 1)*numRows + (len2 - 1));
	TTraceValue tv = 0;
	if (initialDir == Diagonal) {
		if ( ((unsigned char) nextTraceValue >> 2) == 0) tv = Diagonal;
		else if ( ((unsigned char) nextTraceValue >> 2) == 1) tv =  Horizontal;
		else tv =  Vertical;
	} else if (initialDir == Horizontal) {
		if (((unsigned char) nextTraceValue >> 1) % (unsigned char) 2 == 0) {
		  _align_trace_print(align, str, id1, --len1, id2, len2, (TSize) 1,  Horizontal);
		  tv =  Diagonal;
		}
		else tv =  Horizontal;
	} else if (initialDir == Vertical) {
		if ((unsigned char) nextTraceValue % (unsigned char) 2 == 0) {
		  _align_trace_print(align, str, id1, len1, id2, --len2, (TSize) 1,  Vertical);
		  tv =  Diagonal;
		}
		else tv =  Vertical;
	}
	TSize segLen = 0;
	TTraceValue tvOld = tv;

	// Now follow the trace
	do {
		nextTraceValue = value(trace, (len1 - 1)*numRows + (len2 - 1));
		if (tv == Diagonal) {
			if (((unsigned char) nextTraceValue >> 2) == 0) tv = Diagonal;
			else if (((unsigned char) nextTraceValue >> 2) == 1) tv = Horizontal;
			else tv =  Vertical;
		} else if (tv == Horizontal) {
			if ((((unsigned char) nextTraceValue >> 1) % (unsigned char) 2) == 0) tv =  Diagonal;
		    else tv =  Horizontal;
		} else if (tv == Vertical) {
			if (((unsigned char) nextTraceValue %  (unsigned char) 2) == 0) tv =  Diagonal;
			else tv =  Vertical;
		}
		if (tv == Diagonal) {
			if (tv != tvOld) {
				if (tvOld ==  Vertical) --len2;
				else --len1;
				_align_trace_print(align, str, id1, len1, id2, len2, ++segLen, tvOld);
				tvOld = tv; segLen = 0;
			} else {
				++segLen;
				--len1; --len2;
			}
		} else if(tv == Horizontal) {
			if (tv != tvOld) {
				_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
				if ((((unsigned char) nextTraceValue >> 1) % (unsigned char) 2) == 0) {
					_align_trace_print(align, str, id1, --len1, id2, len2, (TSize) 1,  Horizontal);
					tv =  Diagonal; segLen = 0;
				} else {
					tvOld = tv; segLen = 1;
					--len1;
				}
			} else {
				++segLen;
				--len1;
			}
		} else if (tv == Vertical) {
			if (tv != tvOld) {
				_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
				if (((unsigned char) nextTraceValue %  (unsigned char) 2) == 0) {
					_align_trace_print(align, str, id1, len1, id2, --len2, (TSize) 1,  Vertical);
					tv =  Diagonal; segLen = 0;
				} else {
					tvOld = tv; segLen = 1;
					--len2;
				}
			} else {
				++segLen;
				--len2;
			}
		}
	} while ((len1 != 0) && (len2 !=0));

	// Process left-overs
	if (segLen) _align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);

	// Handle the remaining sequence
	if (len1 != 0) _align_trace_print(align, str, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) len1,  Horizontal);
	else if (len2 != 0) _align_trace_print(align, str, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) len2,  Vertical);
}



//////////////////////////////////////////////////////////////////////////////

template <typename TTrace, typename TStringSet, typename TScore, typename TValPair, typename TIndexPair, typename TAlignConfig>
inline typename Value<TScore>::Type
_align_gotoh(TTrace& trace,	     
			 TStringSet const& str,
			 TScore const & sc,
			 TValPair& overallMaxValue,
			 TIndexPair& overallMaxIndex,
			 typename Value<TTrace>::Type& initialDir,
			 TAlignConfig const)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Value<TTrace>::Type TTraceValue;

	// TraceBack values for Gotoh
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;

	// The DP Matrix for diagonal walks
	typedef typename Value<TScore>::Type TScoreValue;
	typedef String<TScoreValue> TColumn;
	typedef typename Size<TColumn>::Type TSize;
	TColumn mat;
	// The DP Matrix for gaps from the left
	TColumn horizontal;
	// The DP Matrix for gaps from the top
	TScoreValue vert = 0;

	// Initialization
	TString const& str1 = str[0];
	TString const& str2 = str[1];		
	TSize len1 = length(str1);
	TSize len2 = length(str2);
	TScoreValue gap = scoreGapExtend(sc);
	TScoreValue gapOpen = scoreGapOpen(sc);
	resize(mat, (len2+1));   // One column for the diagonal matrix
	resize(horizontal, (len2+1));   // One column for the horizontal matrix
	resize(trace, len1*len2);
	TTraceValue tvMat, tvHorizontal, tvVertical;
	
	// Classical DP
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	TTraceIter it = begin(trace, Standard() );
	overallMaxValue[0] = InfimumValue<TScoreValue>::VALUE;
	overallMaxValue[1] = InfimumValue<TScoreValue>::VALUE;
	overallMaxIndex[0] = len1;
	overallMaxIndex[1] = len2;
	mat[0] = 0;

	TScoreValue a = 0;
	TScoreValue b = 0;
	TScoreValue max_val = 0;
	for(TSize row = 1; row <= len2; ++row) {
		_initFirstColumn(TAlignConfig(), mat[row], gapOpen + (row - 1) * gap);
		horizontal[row] = mat[row] + gapOpen - gap;
		//std::cout << getValue(mat, row) << std::endl;
		//std::cout << getValue(horizontal, row) << std::endl;
		//std::cout << "====" << std::endl;
	}
	for(TSize col = 1; col <= len1; ++col) {
		TScoreValue diagValMat = mat[0];
		_initFirstRow(TAlignConfig(), mat[0], gapOpen + (col - 1) * gap);
		vert = mat[0] + gapOpen - gap;
		for(TSize row = 1; row <= len2; ++row, ++it) {
			// Get the new maximum for vertical
			a = mat[row - 1] + gapOpen;
			b = vert + gap;
			if (a > b) {vert = a;	tvVertical = Diagonal;}
			else { vert = b; tvVertical = Vertical;}

			// Get the new maximum for left
			a = mat[row] + gapOpen;
			b = horizontal[row] + gap;
			if (a > b) {horizontal[row] = a; tvHorizontal = Diagonal;}
			else {horizontal[row] = b; tvHorizontal = Horizontal;}

			// Get the new maximum for mat
			max_val = diagValMat + score(const_cast<TScore&>(sc), col-1, row-1, str1, str2);
			tvMat = Diagonal;
			if (vert > max_val) {
				max_val = vert;
				tvMat = Vertical;
			}
			if (horizontal[row] > max_val) {
				max_val = horizontal[row];
				tvMat = Horizontal;
			}

			// Assign the new diagonal values
			diagValMat = mat[row];
			mat[row] = max_val;

			// Assign the right trace value
			if (tvMat == Diagonal) {
				if (tvHorizontal == Diagonal) *it = (tvVertical == Diagonal) ? 0 : 1;
				else if (tvHorizontal == Horizontal) *it = (tvVertical == Diagonal) ? 2 : 3;
			} else if (tvMat ==  Horizontal) {
				if (tvHorizontal ==  Diagonal) *it = (tvVertical ==  Diagonal) ? 4 : 5;
				else if (tvHorizontal ==  Horizontal) *it = (tvVertical ==  Diagonal) ? 6 : 7;
			} else if (tvMat ==  Vertical) {
				if (tvHorizontal ==  Diagonal) *it = (tvVertical ==  Diagonal) ? 8 : 9;
				else if (tvHorizontal ==  Horizontal) *it =  (tvVertical ==  Diagonal) ? 10 : 11;
			}
		}
		_lastRow(TAlignConfig(), overallMaxValue, overallMaxIndex, mat[len2], col);
		// If we got a new index, store direction
		if (overallMaxIndex[0] == col) initialDir = tvMat;
	}
	_lastColumn(TAlignConfig(), overallMaxValue, overallMaxIndex, mat);
	
	// If we got a new index, store direction
	if ((overallMaxIndex[1] != len2)  && (overallMaxValue[1] > overallMaxValue[0])) 
		initialDir = (horizontal[overallMaxIndex[1]] == mat[overallMaxIndex[1]]) ? Horizontal : Diagonal;

	// If we end up in the bottom right corner, get direction
	if ((overallMaxIndex[0] == len1) && (overallMaxIndex[1] == len2)) {
		initialDir =  Diagonal;
		if (horizontal[len2] == mat[len2]) initialDir =  Horizontal;
		else if (vert == mat[len2]) initialDir =  Vertical;
	}

	//// Debug code
	//for(TSize i= 0; i<len2;++i) {
	//	for(TSize j= 0; j<len1;++j) {
	//		std::cout << (TSize) getValue(trace, j*len2 + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << (TSize) initialDir << std::endl;

	return _maxOfAlignment<TScoreValue>(TAlignConfig(), overallMaxValue, overallMaxIndex, len1, len2);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScore, typename TAlignConfig>
inline typename Value<TScore>::Type
_globalAlignment(TAlign& align,
				 TStringSet const& str,
				 TScore const& sc,
				 TAlignConfig const,
				 Gotoh)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	  
	// Trace
	String<TraceBackGotoh> trace;
	TraceBackGotoh initialDir;
	TScoreValue overallMaxValue[2];
	TSize overallMaxIndex[2];

	// Create the trace
	TScoreValue maxScore = _align_gotoh(trace, str, sc, overallMaxValue, overallMaxIndex, initialDir, TAlignConfig());
	// Follow the trace and create the graph
	_align_gotoh_trace(align, str, trace, overallMaxIndex, initialDir);

	return maxScore;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScore, typename TAlignConfig>
inline typename Value<TScore>::Type
_globalAlignment(TStringSet const& str,
				 TScore const& sc,
				 TAlignConfig const,
				 Gotoh)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	TraceBackGotoh initialDir;
	String<TraceBackGotoh> trace;
	TScoreValue overallMaxValue[2];
	TSize overallMaxIndex[2];
	return _align_gotoh(trace, str, sc, overallMaxValue, overallMaxIndex, initialDir, TAlignConfig());	
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
