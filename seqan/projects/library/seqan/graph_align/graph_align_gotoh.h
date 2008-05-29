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
	TSize len1 = overallMaxIndex.first;
	TSize len2 = overallMaxIndex.second;
	if (len1 < length(str[0])) {
		_align_trace_print(align, str, id1, len1, id2, len2, length(str[0]) - len1,  Horizontal);
	} else if (len2 < length(str[1])) {
		_align_trace_print(align, str, id1, len1, id2, len2, length(str[1]) - len2,  Vertical);
	}
	TSize numRows = length(str[1]);

	// Initialize everything	
	TTraceValue nextTraceValue = getValue(trace, (len1 - 1)*numRows + (len2 - 1));
	TTraceValue tv = 0;
	if (initialDir == Diagonal) {
		if ( (Byte) nextTraceValue /  (Byte) 4 == 0) tv =  Diagonal;
		else if ( (Byte) nextTraceValue /  (Byte) 4 == 1) tv =  Horizontal;
		else tv =  Vertical;
	} else if (initialDir == Horizontal) {
		if (( (Byte) nextTraceValue /  (Byte) 2) % (Byte) 2 == 0) {
		  _align_trace_print(align, str, id1, --len1, id2, len2, (TSize) 1,  Horizontal);
		  tv =  Diagonal;
		}
		else tv =  Horizontal;
	} else if (initialDir == Vertical) {
		if (  (Byte) nextTraceValue % (Byte) 2 == 0) {
		  _align_trace_print(align, str, id1, len1, id2, --len2, (TSize) 1,  Vertical);
		  tv =  Diagonal;
		}
		else tv =  Vertical;
	}
	TSize segLen = 0;
	TTraceValue tvOld = tv;

	// Now follow the trace
	do {
		nextTraceValue = getValue(trace, (len1 - 1)*numRows + (len2 - 1));
		if (tv == Diagonal) {
			if (  (Byte) nextTraceValue /  (Byte) 4 == 0) tv =  Diagonal;
			else if ( (Byte) nextTraceValue /  (Byte) 4 == 1) tv =  Horizontal;
			else tv =  Vertical;
		} else if (tv == Horizontal) {
			if (( (Byte) nextTraceValue /  (Byte) 2) % (Byte) 2 == 0) tv =  Diagonal;
		    else tv =  Horizontal;
		} else if (tv == Vertical) {
			if (  (Byte) nextTraceValue %  (Byte) 2 == 0) tv =  Diagonal;
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
				if (( (Byte) nextTraceValue /  (Byte) 2) % (Byte) 2 == 0) {
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
				if (  (Byte) nextTraceValue %  (Byte) 2 == 0) {
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
	TTraceValue Diagonal = 0;
	TTraceValue Horizontal = 1;
	TTraceValue Vertical = 2;

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
	TScoreValue tmp = 0;
	resize(mat, (len2+1));   // One column for the diagonal matrix
	resize(horizontal, (len2+1));   // One column for the horizontal matrix
	resize(trace, len1*len2);
	TTraceValue tvMat, tvHorizontal, tvVertical;
	
	// Classical DP
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	TTraceIter it = begin(trace, Standard() );
	// Max values: overall.first = last column, overall.second = last row
	TScoreValue infValue = -1 * (_getInfinity<TScoreValue>() / 2);
	overallMaxValue = std::make_pair(infValue, infValue);
	//overallMaxValue = std::make_pair(InfimumValue<TScoreValue>::VALUE, InfimumValue<TScoreValue>::VALUE);
	overallMaxIndex = std::make_pair(len1, len2);
	assignValue(mat, 0, 0);
	for(TSize row = 1; row <= len2; ++row) {
		_initFirstColumn(TAlignConfig(), value(mat, row), gapOpen + (row - 1) * gap);
		assignValue(horizontal, row, getValue(mat, row) + gapOpen - gap);
		//std::cout << getValue(mat, row) << std::endl;
		//std::cout << getValue(horizontal, row) << std::endl;
		//std::cout << "====" << std::endl;
	}
	for(TSize col = 1; col <= len1; ++col) {
		TScoreValue diagValMat = getValue(mat, 0);
		_initFirstRow(TAlignConfig(), value(mat, 0), gapOpen + (col - 1) * gap);
		vert = getValue(mat, 0) + gapOpen - gap;
		for(TSize row = 1; row <= len2; ++row) {
			// Get the new maximum for vertical
			if ((tmp = getValue(mat, row - 1) + gapOpen) > vert + gap) {
				vert = tmp;
				tvVertical = Diagonal;
			} else {
				vert = vert + gap;
				tvVertical = Vertical;
			}

			// Get the new maximum for left
			if ((tmp = getValue(mat, row) + gapOpen) > getValue(horizontal, row) + gap) {
				assignValue(horizontal, row, tmp);
				tvHorizontal = Diagonal;
			} else {
				assignValue(horizontal, row, getValue(horizontal, row) + gap);
				tvHorizontal = Horizontal;
			}

			// Get the new maximum for mat
			TScoreValue sc_ = score(const_cast<TScore&>(sc), str1[col-1], str2[row-1]);
			tmp = diagValMat + sc_;
			tvMat = Diagonal;
			if (vert > tmp) {
				tmp = vert;
				tvMat = Vertical;
			}
			if (getValue(horizontal, row) > tmp) {
				tmp = getValue(horizontal,row);
				tvMat = Horizontal;
			}

			// Assign the new diagonal values
			diagValMat = getValue(mat, row);
			assignValue(mat, row, tmp);

			// Assign the right trace value
			if (tvMat == Diagonal) {
				if (tvHorizontal == Diagonal) {
					if (tvVertical == Diagonal) assignValue(it, 0);
					else assignValue(it, 1);
				} else if (tvHorizontal == Horizontal) {
					if (tvVertical == Diagonal) assignValue(it, 2);
					else assignValue(it, 3);
				}
			} else if (tvMat ==  Horizontal) {
				if (tvHorizontal ==  Diagonal) {
					if (tvVertical ==  Diagonal) assignValue(it, 4);
					else assignValue(it, 5);
				} else if (tvHorizontal ==  Horizontal) {
					if (tvVertical ==  Diagonal) assignValue(it, 6);
					else assignValue(it, 7);
				}
			} else if (tvMat ==  Vertical) {
				if (tvHorizontal ==  Diagonal) {
					if (tvVertical ==  Diagonal) assignValue(it, 8);
					else assignValue(it, 9);
				} else if (tvHorizontal ==  Horizontal) {
					if (tvVertical ==  Diagonal) assignValue(it, 10);
					else assignValue(it, 11);
				}
			}
			goNext(it);
		}
		_processLastRow(TAlignConfig(), overallMaxValue, overallMaxIndex, getValue(mat, len2), col);
		// If we got a new index, store direction
		if (overallMaxIndex.first == col) {
			tmp = getValue(mat, len2);
			initialDir =  Diagonal;
			if (vert == tmp) {
				initialDir =  Vertical;
			}
		}
	}
	_processLastColumn(TAlignConfig(), overallMaxValue, overallMaxIndex, mat);
	// If we got a new index, store direction
	if (overallMaxIndex.second != len2) {
		if (overallMaxValue.second > overallMaxValue.first) {
			overallMaxIndex.first = len1;
			tmp = getValue(mat, overallMaxIndex.second);
			initialDir =  Diagonal;
			if (getValue(horizontal, overallMaxIndex.second) ==  tmp) {
				initialDir =  Horizontal;
			}
		} else if (overallMaxIndex.first != len1) {
			overallMaxIndex.second = len2;
		}
	}

	// If we end up in the bottom right corner, get direction
	if ((overallMaxIndex.first == len1) &&
		(overallMaxIndex.second == len2)) {
		tmp = getValue(mat, len2);
		initialDir =  Diagonal;
		if (getValue(horizontal, len2) ==  tmp) {
			initialDir =  Horizontal;
		}
		else if (vert == tmp) {
			initialDir =  Vertical;
		}
	}

	//// Debug code
	//for(TSize i= 0; i<len2;++i) {
	//	for(TSize j= 0; j<len1;++j) {
	//		std::cout << (TSize) getValue(trace, j*len2 + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << (TSize) initialDir << std::endl;

	return _retrieveMaxOfAlignment(TAlignConfig(), overallMaxValue);
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
	std::pair<TScoreValue, TScoreValue> overallMaxValue;
	std::pair<TSize, TSize> overallMaxIndex;

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
	std::pair<TScoreValue, TScoreValue> overallMaxValue;
	std::pair<TSize, TSize> overallMaxIndex;
	return _align_gotoh(trace, str, sc, overallMaxValue, overallMaxIndex, initialDir, TAlignConfig());	
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
