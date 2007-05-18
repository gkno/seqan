#ifndef SEQAN_HEADER_GRAPH_NEEDLEMAN_WUNSCH_H
#define SEQAN_HEADER_GRAPH_NEEDLEMAN_WUNSCH_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Needleman-Wunsch Alignment, constant gap cost
// Gap extension score is taken as the constant gap score!!!
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


template <typename TAlign, typename TStringSet, typename TTrace>
void
_align_needleman_wunsch_trace(TAlign& align,
							  TStringSet const& str,
							  TTrace const& trace)
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
	TSize len1 = length(str[0]);
	TSize len2 = length(str[1]);
	TSize numRows = len2;
		
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

template <typename TTrace, typename TStringSet, typename TScoreValue>
TScoreValue
_align_needleman_wunsch(TTrace& trace,
			TStringSet const& str,
			Score<TScoreValue, Simple> const& sc,
			bool createTrace) 
{
	SEQAN_CHECKPOINT
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	typedef String<TScoreValue> TColumn;
	typedef typename Size<TColumn>::Type TSize;
	typedef typename Value<TStringSet>::Type TString;

	// TraceBack values
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};

	// One DP Matrix column
	TColumn column;
		
	// Initialization
	TString const& str1 = str[0];
	TString const& str2 = str[1];
	TSize len1 = length(str1);
	TSize len2 = length(str2);
	TScoreValue gap = scoreGapExtend(sc);
	resize(column, len2 + 1);  
	if (createTrace) resize(trace, len1*len2);
	for(TSize row = 0; row <= len2; ++row) assignValue(column, row, row * gap);

	// Classical DP
	TTraceIter it = begin(trace, Standard() );
	TTraceValue tv;
	TScoreValue maxVal = 0;
	TScoreValue tmp = 0;
	for(TSize col = 1; col <= len1; ++col) {
		TScoreValue diagVal = getValue(column, 0);
		assignValue(column, 0, col * gap);
		for(TSize row = 1; row <= len2; ++row) {
			// Get the new maximum	
			maxVal = diagVal + score(const_cast<Score<TScoreValue, Simple>&>(sc), str1[col-1], str2[row-1]);
			tv = (Byte) Diagonal;
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
			if (createTrace) {
			  assignValue(it, tv);
			  goNext(it);
			}
		}
	}

	//for(unsigned int i= 0; i<len2;++i) {
	//	for(unsigned int j= 0; j<len1;++j) {
	//		std::cout << (unsigned int) getValue(trace, j*len2 + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}

	return getValue(column, len2);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TTrace, typename TStringSet, typename TScoreValue>
TScoreValue
_align_needleman_wunsch(TTrace& trace,
			TStringSet const& str,
			Score<TScoreValue, Simple> const& sc) 
{
	SEQAN_CHECKPOINT
	return _align_needleman_wunsch(trace,str,sc,true);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScoreValue>
TScoreValue
_globalAlignment(TAlign& align,
				 TStringSet const& str,
				 Score<TScoreValue, Simple> const& sc,
				 NeedlemanWunsch)
				 
{
	SEQAN_CHECKPOINT
	// Gap extension score is taken as the constant gap score!!!
	SEQAN_TASSERT(scoreGapOpen(sc) == 0)

	typedef typename Size<TStringSet>::Type TSize;
	  
	TScoreValue maxScore;
	TSize maxLen = length(str[0]);
	TSize tmp;
	if ((tmp = length(str[1])) > maxLen) maxLen = tmp;

	if (maxLen > 100000) {
		String<TraceBack, External<> > trace;
		//open(trace, "D:\\seqan.dat");

		// Create the trace
		maxScore = _align_needleman_wunsch(trace, str, sc);	

		// Follow the trace and create the graph
		_align_needleman_wunsch_trace(align, str, trace);	
	} else {
		String<TraceBack> trace;

		// Create the trace
		maxScore = _align_needleman_wunsch(trace, str, sc);	

		// Follow the trace and create the graph
		_align_needleman_wunsch_trace(align, str, trace);	
	}
	return maxScore;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScoreValue>
TScoreValue
_globalAlignment(TStringSet const& str,
				 Score<TScoreValue, Simple> const& sc,
				 NeedlemanWunsch)
{
	SEQAN_CHECKPOINT
	// Gap extension score is taken as the constant gap score!!!
	SEQAN_TASSERT(scoreGapOpen(sc) == 0)
	String<TraceBack> trace;
	return _align_needleman_wunsch(trace, str, sc, false);	
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
