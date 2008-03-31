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
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_NEEDLEMAN_WUNSCH_H
#define SEQAN_HEADER_GRAPH_NEEDLEMAN_WUNSCH_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Needleman-Wunsch Alignment, constant gap cost
// Gap extension score is taken as the constant gap score!!!
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


template <typename TAlign, typename TStringSet, typename TTrace, typename TIndexPair>
void
_align_needleman_wunsch_trace(TAlign& align,
							  TStringSet const& str,
							  TTrace const& trace,
							  TIndexPair const& overallMaxIndex)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Id<TStringSet>::Type TId;
	typedef typename Value<TTrace>::Type TTraceValue;


	// TraceBack values
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};
	
	// Initialization
	TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
	TId id2 = positionToId(const_cast<TStringSet&>(str), 1);
	TSize len1 = overallMaxIndex.first;
	TSize len2 = overallMaxIndex.second;
	if (len1 < length(str[0])) {
		_align_trace_print(align, str, id1, len1, id2, len2, length(str[0]) - len1, (Byte) Horizontal);
	} else if (len2 < length(str[1])) {
		_align_trace_print(align, str, id1, len1, id2, len2, length(str[1]) - len2, (Byte) Vertical);
	}
	TSize numRows = length(str[1]);
		
	// Initialize everything
	TTraceValue tv = getValue(trace, (len1-1)*numRows + (len2-1));
	TTraceValue tvOld = tv;  // We need to know when the direction of the trace changes

	TSize segLen = 1;
	if (tv == (Byte) Diagonal) {
		--len1; --len2;
	}
	else if (tv == (Byte) Horizontal) --len1;
	else if (tv == (Byte) Vertical) --len2;

	// Now follow the trace
	do {
		tv = getValue(trace, (len1-1)*numRows + (len2-1));
		if (tv == (Byte) Diagonal) {
			if (tv != tvOld) {
				_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
				tvOld = tv; segLen = 1;
			} else ++segLen;
			--len1; --len2;
		} else if (tv == (Byte) Horizontal) {
			//std::cout << '(' << ((*str)[0])[len1 - 1] << ',' << '-' << ')' << std::endl;
			if (tv != tvOld) {
				_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
				tvOld = tv; segLen = 1;
			} else ++segLen;
			--len1;
		} else if (tv == (Byte) Vertical) {
			//std::cout << '(' << '-' << ',' << ((*str)[1])[len2-1] << ')' << std::endl;
			if (tv != tvOld) {
				_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
				tvOld = tv; segLen = 1;
			} else ++segLen;
			--len2;
		}
	} while ((len1 != 0) && (len2 !=0));
	// Process left-overs
	_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);

	// Handle the remaining sequence
	if (len1 != 0) _align_trace_print(align, str, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) len1, (Byte) Horizontal);
	else if (len2 != 0) _align_trace_print(align, str, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) len2, (Byte) Vertical);
}

	

//////////////////////////////////////////////////////////////////////////////

template <typename TTrace, typename TStringSet, typename TScore, typename TValPair, typename TIndexPair, typename TAlignConfig>
inline typename Value<TScore>::Type
_align_needleman_wunsch(TTrace & trace,
						TStringSet const & str,
						TScore const & sc,
						TValPair & overallMaxValue,
						TIndexPair & overallMaxIndex,
						TAlignConfig const) 
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;

	// TraceBack values
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};

	TScore const & _sc = sc;

	// One DP Matrix column
	typedef String<TScoreValue> TColumn;
	typedef typename Size<TColumn>::Type TSize;
	TColumn column;

	// Initialization
	typedef typename Value<TStringSet>::Type TString;
	TString const& str1 = str[0];
	TString const& str2 = str[1];
	TSize len1 = length(str1);
	TSize len2 = length(str2);
	TScoreValue gap = scoreGapExtend(sc);
	resize(column, len2 + 1);  
	resize(trace, len1*len2);
	for(TSize row = 0; row <= len2; ++row) _initFirstColumn(TAlignConfig(), value(column, row), row*gap);

	typedef typename Iterator<TColumn, Standard>::Type TColIterator;
	TColIterator col_end = end(column, Standard());

	// Classical DP
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	typedef typename Value<TTrace>::Type TTraceValue;
	TTraceIter it = begin(trace, Standard() );
	// Max values: overall.first = last column, overall.second = last row
	overallMaxValue = std::make_pair(InfimumValue<TScoreValue>::VALUE, InfimumValue<TScoreValue>::VALUE);
	overallMaxIndex = std::make_pair(len1, len2);

	for(TSize col = 1; col <= len1; ++col) 
	{
		TScoreValue diagVal = getValue(column, 0);
		_initFirstRow(TAlignConfig(), value(column, 0), col * gap);
		TScoreValue vertiVal = getValue(column, 0);

		typedef typename Value<TString>::Type TStringValue;
		TStringValue char1 = str1[col-1];

		TColIterator coit = begin(column, Standard());

		typedef typename Iterator<TString, Standard>::Type TStringIterator;
		TStringIterator s2_it = begin(str2, Standard());

		while (true)
		{
			++coit;
			if (coit == col_end) break;

			TScoreValue max_hori_verti;
			TTraceValue tv_hori_verti;
			if (*coit >= vertiVal)
			{
				max_hori_verti = *coit;
				tv_hori_verti = (Byte) Horizontal;
			}
			else
			{
				max_hori_verti = vertiVal;
				tv_hori_verti = (Byte) Vertical;
			}
			max_hori_verti += gap;

			vertiVal = diagVal + score(_sc, char1, *s2_it); //compute the maximum in vertiVal 
			++s2_it;

			if (vertiVal >= max_hori_verti)
			{
				*it = (Byte) Diagonal;
				++it;
			}
			else
			{
				vertiVal = max_hori_verti;

				*it = tv_hori_verti;
				++it;
			}

			diagVal = *coit;
			*coit = vertiVal;
		}
		_processLastRow(TAlignConfig(), overallMaxValue, overallMaxIndex, vertiVal, col);
	}
	_processLastColumn(TAlignConfig(), overallMaxValue, overallMaxIndex, column);

	//for(TSize i= 0; i<len2;++i) {
	//	for(TSize j= 0; j<len1;++j) {
	//		std::cout << (TSize) getValue(trace, j*len2 + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}

	return _retrieveMaxOfAlignment(TAlignConfig(), overallMaxValue);
}

//*/
/*
template <typename TTrace, typename TStringSet, typename TScore, typename TValPair, typename TIndexPair, typename TAlignConfig>
inline typename Value<TScore>::Type
_align_needleman_wunsch(TTrace& trace,
						TStringSet const& str,
						TScore const& sc,
						TValPair& overallMaxValue,
						TIndexPair& overallMaxIndex,
						TAlignConfig const) 
{
	typedef typename Value<TScore>::Type TScoreValue;

	// TraceBack values
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};

	// One DP Matrix column
	typedef String<TScoreValue> TColumn;
	typedef typename Size<TColumn>::Type TSize;
	TColumn column;
		
	// Initialization
	typedef typename Value<TStringSet>::Type TString;
	TString const& str1 = str[0];
	TString const& str2 = str[1];
	TSize len1 = length(str1);
	TSize len2 = length(str2);
	TScoreValue gap = scoreGapExtend(sc);
	resize(column, len2 + 1);  
	resize(trace, len1*len2);
	for(TSize row = 0; row <= len2; ++row) _initFirstColumn(TAlignConfig(), value(column, row), row*gap);

	// Classical DP
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	typedef typename Value<TTrace>::Type TTraceValue;
	TTraceIter it = begin(trace, Standard() );
	// Max values: overall.first = last column, overall.second = last row
	overallMaxValue = std::make_pair(InfimumValue<TScoreValue>::VALUE, InfimumValue<TScoreValue>::VALUE);
	overallMaxIndex = std::make_pair(len1, len2);
	for(TSize col = 1; col <= len1; ++col) {
		TScoreValue diagVal = getValue(column, 0);
		_initFirstRow(TAlignConfig(), value(column, 0), col * gap);
		TScoreValue maxVal = 0;
		for(TSize row = 1; row <= len2; ++row) {	
			// Get the new maximum	
			TScoreValue tmp = 0;
			maxVal = diagVal + score(const_cast<TScore&>(sc), str1[col-1], str2[row-1]);
			TTraceValue tv = (Byte) Diagonal;
			if ((tmp = getValue(column, row) + gap) > maxVal) {
				maxVal = tmp;
				tv = (Byte) Horizontal;
			}
			if ((tmp = getValue(column, (row - 1)) + gap) > maxVal) {
				maxVal = tmp;
				tv = (Byte) Vertical;
			}
			diagVal = getValue(column, row);
			// Assign the new value
			assignValue(column, row, maxVal);
			
			// Remember the trace
			assignValue(it, tv);
			goNext(it);
		}
		_processLastRow(TAlignConfig(), overallMaxValue, overallMaxIndex, maxVal, col);
	}
	_processLastColumn(TAlignConfig(), overallMaxValue, overallMaxIndex, column);

	//for(TSize i= 0; i<len2;++i) {
	//	for(TSize j= 0; j<len1;++j) {
	//		std::cout << (TSize) getValue(trace, j*len2 + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}

	return _retrieveMaxOfAlignment(TAlignConfig(), overallMaxValue);
}
//*/
//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScore, typename TAlignConfig>
inline typename Value<TScore>::Type
_globalAlignment(TAlign& align,
				 TStringSet const& str,
				 TScore const& sc,
				 TAlignConfig const,
				 NeedlemanWunsch)
				 
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;

	// Gap extension score is taken as the constant gap score!!!
	SEQAN_TASSERT(scoreGapOpen(sc) == 0)

	// Create the trace
	String<TraceBack> trace;
	std::pair<TScoreValue, TScoreValue> overallMaxValue;
	std::pair<TSize, TSize> overallMaxIndex;
	TScoreValue	maxScore = _align_needleman_wunsch(trace, str, sc, overallMaxValue, overallMaxIndex, TAlignConfig());	

	// Follow the trace and create the graph
	_align_needleman_wunsch_trace(align, str, trace, overallMaxIndex);	

	return maxScore;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScore, typename TAlignConfig>
inline typename Value<TScore>::Type
_globalAlignment(TStringSet const& str,
				 TScore const& sc,
				 TAlignConfig const,
				 NeedlemanWunsch)
{
	SEQAN_CHECKPOINT

	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	
	// Gap extension score is taken as the constant gap score!!!
	SEQAN_TASSERT(scoreGapOpen(sc) == 0)
	String<TraceBack> trace;
	std::pair<TScoreValue, TScoreValue> overallMaxValue;
	std::pair<TSize, TSize> overallMaxIndex;
	return _align_needleman_wunsch(trace, str, sc, overallMaxValue, overallMaxIndex, TAlignConfig());	
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
