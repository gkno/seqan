#ifndef SEQAN_HEADER_GRAPH_ALIGN_GOTOH_H
#define SEQAN_HEADER_GRAPH_ALIGN_GOTOH_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Gotoh Alignment, affine gap cost
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TStringSet, typename TTrace, typename TVal>
void
_align_gotoh_trace(TAlign& align,		 
				   TStringSet const& str,
				   TTrace const& trace,
				   TVal const initialDir)
{
	SEQAN_CHECKPOINT
	
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Id<TStringSet>::Type TId;

	// TraceBack values for Gotoh
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};

	TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
	TId id2 = positionToId(const_cast<TStringSet&>(str), 1);	 
	TSize len1 = length(str[0]);
	TSize len2 = length(str[1]);
	TSize numRows = len2;

	// Initialize everything	
	TTraceValue tv = getValue(trace, (len1 - 1)*numRows + (len2 - 1));
	TTraceValue tvOld = initialDir;
	switch( (Byte) initialDir) {
		case Diagonal:
			if ((Byte) tv / (Byte) 4 == 0) tv = (Byte) Diagonal;
			else if ((Byte) tv / (Byte) 4 == 1) tv = (Byte) Horizontal;
			else tv = (Byte) Vertical;
			break;
		case Horizontal:
			if (((Byte) tv / (Byte) 2) % 2 == 0) tv = (Byte) Diagonal;
			else tv = (Byte) Horizontal;
			break;
		case Vertical:
			if ( (Byte) tv % 2 == 0) tv = (Byte) Diagonal;
			else tv = (Byte) Vertical;
			break;
	}
	TSize segLen = 0;

	// Now follow the trace
	do {
		switch( (Byte) tvOld) {
			case Diagonal: 
				//std::cout << '(' << ((*str)[0])[len1 - 1] << ',' << ((*str)[1])[len2-1] << ')' << std::endl;
				--len1; --len2;
				++segLen;
				if (tv != tvOld) {
					_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
					tvOld = tv; 
					segLen = 0;
				}
				break;
			case Horizontal:
				//std::cout << '(' << ((*str)[0])[len1 - 1] << ',' << '-' << ')' << std::endl;
				--len1;
				++segLen;
				if (tv != tvOld) {
					_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
					tvOld = tv;
					segLen = 0;
				}
				break;
			case Vertical:
				//std::cout << '(' << '-' << ',' << ((*str)[1])[len2-1] << ')' << std::endl;
				--len2;
				++segLen;
				if (tv != tvOld) {
					_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
					tvOld = tv; 
					segLen = 0;
				}
				break;
		}
		if ((len1 != 0) && (len2 !=0)) {
			TTraceValue nextTraceValue = getValue(trace, (len1 - 1)*numRows + (len2 - 1));
			switch( (Byte) tv) {
				case Diagonal:
					if ( (Byte) nextTraceValue / (Byte) 4 == 0) tv = (Byte) Diagonal;
					else if ((Byte) nextTraceValue / (Byte) 4 == 1) tv = (Byte) Horizontal;
					else tv = (Byte) Vertical;
					break;
				case Horizontal:
					if (((Byte) nextTraceValue / (Byte) 2) % 2 == 0) tv = (Byte) Diagonal;
					else tv = (Byte) Horizontal;
					break;
				case Vertical:
					if ( (Byte) nextTraceValue % (Byte) 2 == 0) tv = (Byte) Diagonal;
					else tv = (Byte) Vertical;
					break;
			}
		}
	} while ((len1 != 0) && (len2 !=0));
	// Process left-overs
	if (segLen > 0) _align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);

	// Handle the remaining sequence
	if (len1 != 0) _align_trace_print(align, str, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) len1, (Byte) Horizontal);
	else if (len2 != 0) _align_trace_print(align, str, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) len2, (Byte) Vertical);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TTrace, typename TStringSet, typename TScoreValue>
TScoreValue
_align_gotoh(TTrace& trace,
	     TStringSet const& str,
	     Score<TScoreValue, Simple> const & sc,
	     typename Value<TTrace>::Type& initialDir,
	     bool createTrace)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	typedef String<TScoreValue> TColumn;
	typedef typename Size<TColumn>::Type TSize;
	typedef typename Value<TStringSet>::Type TString;

	// TraceBack values for Gotoh
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};

	// The DP Matrix for diagonal walks
	TColumn mat;
	// The DP Matrix for gaps from the left
	TColumn horizontal;
	// The DP Matrix for gaps from the top
	TColumn vertical;

	// Initialization
	TString const& str1 = str[0];
	TString const& str2 = str[1];		
	TSize len1 = length(str1);
	TSize len2 = length(str2);
	TScoreValue gap = scoreGapExtend(sc);
	TScoreValue gapOpen = scoreGapOpen(sc);
	TScoreValue maxVal = 0;
	TScoreValue tmp = 0;
	resize(mat, (len2+1));   // One column for the diagonal matrix
	resize(horizontal, (len2+1));   // One column for the horizontal matrix
	resize(vertical, (len2+1));   // One column for the vertical matrix
	if (createTrace) resize(trace, len1*len2);

	TTraceValue tvMat, tvHorizontal, tvVertical;
	
	TScoreValue inf = infimumValue<TScoreValue>() / 2;

	// Classical DP
	TTraceIter it = begin(trace, Standard() );
	assignValue(mat, 0, 0);
	assignValue(horizontal, 0, inf);
	assignValue(vertical, 0, inf);
	for(TSize row = 1; row <= len2; ++row) {
		assignValue(mat, row, gapOpen + (row - 1) * gap);
		assignValue(horizontal, row, inf);
		assignValue(vertical, row, gapOpen + (row - 1) * gap);
	}
	for(TSize col = 1; col <= len1; ++col) {
		TScoreValue diagValMat = getValue(mat, 0);
		TScoreValue diagValHori = getValue(horizontal, 0);
		TScoreValue diagValVert = getValue(vertical, 0);
		TScoreValue diagValVertTmp;
		TScoreValue diagValHoriTmp;
		assignValue(mat, 0, gapOpen + (col - 1) * gap);
		assignValue(horizontal, 0, gapOpen + (col - 1) * gap);
		assignValue(vertical, 0, inf);
		for(TSize row = 1; row <= len2; ++row) {
			// Get the new maximum for vertical
			maxVal = getValue(mat, row - 1) + gapOpen;
			tvVertical = (Byte) Diagonal;
			if ((tmp = getValue(vertical, row - 1) + gap) > maxVal) {
				maxVal = tmp;
				tvVertical = (Byte) Vertical;
			}
			diagValVertTmp = getValue(vertical, row);
			assignValue(vertical, row, maxVal);

			// Get the new maximum for left
			maxVal = getValue(mat, row) + gapOpen;
			tvHorizontal = (Byte) Diagonal;
			if ((tmp = getValue(horizontal, row) + gap) > maxVal) {
				maxVal = tmp;
				tvHorizontal = (Byte) Horizontal;
			}
			diagValHoriTmp = getValue(horizontal, row);
			assignValue(horizontal, row, maxVal);

			// Get the new maximum for mat
			TScoreValue sc_ = score(const_cast<Score<TScoreValue, Simple>&>(sc), str1[col-1], str2[row-1]);
			maxVal = diagValMat + sc_;
			tvMat = (Byte) Diagonal;
			if ((tmp = diagValVert + sc_) > maxVal) {
				maxVal = tmp;
				tvMat = (Byte) Vertical;
			}
			if ((tmp = diagValHori + sc_) > maxVal) {
				maxVal = tmp;
				tvMat = (Byte) Horizontal;
			}

			// Assign the new diagonal values
			diagValMat = getValue(mat, row);
			diagValHori = diagValHoriTmp;
			diagValVert = diagValVertTmp;
			assignValue(mat, row, maxVal);

			// Assign the right trace value
			if (createTrace) {
			  if (tvMat == (Byte) Diagonal) {
			    if (tvHorizontal == (Byte) Diagonal) {
			      if (tvVertical == (Byte) Diagonal) assignValue(it, 0);
			      else assignValue(it, 1);
			    } else if (tvHorizontal == (Byte) Horizontal) {
			      if (tvVertical == (Byte) Diagonal) assignValue(it, 2);
			      else assignValue(it, 3);
			    }
			  } else if (tvMat == (Byte) Horizontal) {
			    if (tvHorizontal == (Byte) Diagonal) {
			      if (tvVertical == (Byte) Diagonal) assignValue(it, 4);
			      else assignValue(it, 5);
			    } else if (tvHorizontal == (Byte) Horizontal) {
			      if (tvVertical == (Byte) Diagonal) assignValue(it, 6);
			      else assignValue(it, 7);
			    }
			  } else if (tvMat == (Byte) Vertical) {
			    if (tvHorizontal == (Byte) Diagonal) {
			      if (tvVertical == (Byte) Diagonal) assignValue(it, 8);
			      else assignValue(it, 9);
			    } else if (tvHorizontal == (Byte) Horizontal) {
			      if (tvVertical == (Byte) Diagonal) assignValue(it, 10);
			      else assignValue(it, 11);
			    }
			  }
			  goNext(it);
			}
		}
	}

	maxVal = getValue(mat, len2);
	initialDir = (Byte) Diagonal;
	if ((tmp = getValue(horizontal, len2)) > maxVal) {
		maxVal = tmp;
		initialDir = (Byte) Horizontal;
	}
	if ((tmp = getValue(vertical, len2)) > maxVal) {
		maxVal = tmp;
		initialDir = (Byte) Vertical;
	}

	/*
	for(unsigned int i= 0; i<len2;++i) {
		for(unsigned int j= 0; j<len1;++j) {
			std::cout << (unsigned int) getValue(trace, j*len2 + i) << ',';
		}
		std::cout << std::endl;
	}
	*/

	return maxVal;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TTrace, typename TStringSet, typename TScoreValue>
TScoreValue
_align_gotoh(TTrace& trace,
	     TStringSet const& str,
	     Score<TScoreValue, Simple> const & sc,
	     typename Value<TTrace>::Type& initialDir)
{
	SEQAN_CHECKPOINT
	return _align_gotoh(trace, str, sc, initialDir, true);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScoreValue>
TScoreValue
_globalAlignment(TAlign& align,
				 TStringSet const& str,
				 Score<TScoreValue, Simple> const& sc,
				 Gotoh)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TStringSet>::Type TSize;
	  
	TScoreValue maxScore;
	TSize maxLen = length(str[0]);
	TSize tmp;
	if ((tmp = length(str[1])) > maxLen) maxLen = tmp;

	if (maxLen > 100000) {
		// Trace
		String<TraceBackGotoh, External<> > trace;
		//open(trace, "D:\\seqan.dat");
		TraceBackGotoh initialDir;

		// Create the trace
		maxScore = _align_gotoh(trace, str, sc, initialDir);	
		// Follow the trace and create the graph
		_align_gotoh_trace(align, str, trace, initialDir);
	} else {
		// Trace
		String<TraceBackGotoh> trace;
		TraceBackGotoh initialDir;

		// Create the trace
		maxScore = _align_gotoh(trace, str, sc, initialDir);	
		// Follow the trace and create the graph
		_align_gotoh_trace(align, str, trace, initialDir);
	}
	return maxScore;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScoreValue>
TScoreValue
_globalAlignment(TStringSet const& str,
		 Score<TScoreValue, Simple> const& sc,
		 Gotoh)
{
	SEQAN_CHECKPOINT
	TraceBackGotoh initialDir;
	String<TraceBackGotoh> trace;
	return _align_gotoh(trace, str, sc, initialDir, false);	
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
