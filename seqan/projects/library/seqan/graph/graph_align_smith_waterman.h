#ifndef SEQAN_HEADER_GRAPH_SMITH_WATERMAN_H
#define SEQAN_HEADER_GRAPH_SMITH_WATERMAN_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Smith Waterman Alignment, affine gap cost
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TStringSet, typename TTrace, typename TVal, typename TSize, typename TForbidden>
inline void
_align_smith_waterman_trace(TAlign& align,
							TStringSet const& str,
							TTrace const& trace,
							TVal const initialDir,
							TSize best_row,
							TSize best_col,
							TForbidden& forbidden)
{
	SEQAN_CHECKPOINT
	
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Id<TStringSet>::Type TId;

	// TraceBack values for Gotoh
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2, Stop = 12};

	TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
	TId id2 = positionToId(const_cast<TStringSet&>(str), 1);	 
	TSize len1 = best_col;
	TSize len2 = best_row;
	if ((best_col == 0) || (best_row == 0)) return;
	TSize numRows = length(str[1]);

	// Initialize everything	
	TTraceValue nextTraceValue = getValue(trace, (len1 - 1)*numRows + (len2 - 1));
	TTraceValue tv = 0;
	switch( (Byte) initialDir) {
		case Diagonal:
			if ((Byte) nextTraceValue / (Byte) 4 == 0) tv = (Byte) Diagonal;
			else if ((Byte) nextTraceValue / (Byte) 4 == 1) tv = (Byte) Horizontal;
			else tv = (Byte) Vertical;
			break;
		case Horizontal:
			if (((Byte) nextTraceValue / (Byte) 2) % 2 == 0) {
			  _align_trace_print(align, str, id1, --len1, id2, len2, (TSize) 1, (Byte) Horizontal);
			  tv = (Byte) Diagonal;
			}
			else tv = (Byte) Horizontal;
			break;
		case Vertical:
			if ( (Byte) nextTraceValue % 2 == 0) {
			  _align_trace_print(align, str, id1, len1, id2, --len2, (TSize) 1, (Byte) Vertical);
			  tv = (Byte) Diagonal;
			}
			else tv = (Byte) Vertical;
			break;
	}
	TSize segLen = 0;
	TTraceValue tvOld = tv;

	// Now follow the trace
	do {
		nextTraceValue = getValue(trace, (len1 - 1)*numRows + (len2 - 1));
		if ((Byte) nextTraceValue == Stop) break;
		if (!empty(forbidden)) assignValue(forbidden, (len1 - 1)*numRows + (len2 - 1), true);
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
	  switch( (Byte) tv) {
		case Diagonal: 
			if (tv != tvOld) {
				if (tvOld == (Byte) Vertical) --len2;
				else --len1;
				_align_trace_print(align, str, id1, len1, id2, len2, ++segLen, tvOld);
				tvOld = tv; segLen = 0;
			} else {
				++segLen;
				--len1; --len2;
			}
			break;
		case Horizontal:
			if (tv != tvOld) {
				_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
				if (((Byte) nextTraceValue / (Byte) 2) % 2 == 0) {
					_align_trace_print(align, str, id1, --len1, id2, len2, (TSize) 1, (Byte) Horizontal);
					tv = (Byte) Diagonal; segLen = 0;
				} else {
					tvOld = tv; segLen = 1;
					--len1;
				}
			} else {
				++segLen;
				--len1;
			}
			break;
		case Vertical:
			if (tv != tvOld) {
				_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
				if ( (Byte) nextTraceValue % (Byte) 2 == 0) {
					_align_trace_print(align, str, id1, len1, id2, --len2, (TSize) 1, (Byte) Vertical);
					tv = (Byte) Diagonal; segLen = 0;
				} else {
					tvOld = tv; segLen = 1;
					--len2;
				}
			} else {
				++segLen;
				--len2;
			}
			break;
	  }
	} while ((len1 != 0) && (len2 !=0));
	// Process left-overs
	if (segLen) _align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TTrace, typename TStringSet, typename TScore, typename TSize, typename TForbidden>
inline typename Value<TScore>::Type
_align_smith_waterman(TTrace& trace,
					  TStringSet const& str,
					  TScore const & sc,
					  typename Value<TTrace>::Type& initialDir,
					  TSize& best_row,
					  TSize& best_col,
					  TForbidden& forbidden)
{
	SEQAN_CHECKPOINT
	// TraceBack values for Smith Waterman
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2, Stop = 12};

	// The DP Matrix for diagonal walks
	typedef typename Value<TScore>::Type TScoreValue;
	typedef String<TScoreValue> TColumn;
	TColumn mat;
	// The DP Matrix for gaps from the left
	TColumn horizontal;
	// The DP Matrix for gaps from the top
	TScoreValue vert = 0;

	// Initialization
	typedef typename Value<TStringSet>::Type TString;
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
	typedef typename Value<TTrace>::Type TTraceValue;
	TTraceValue tvMat=0, tvHorizontal=0, tvVertical=0;

	// Record the max score
	TScoreValue score_max = 0;
	best_row = 0;
	best_col = 0;
	initialDir = (Byte) Stop;
	bool emptyForbidden = empty(forbidden);

	
	// Classical DP
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	TTraceIter it = begin(trace, Standard() );
	assignValue(mat, 0, 0);
	for(TSize row = 1; row <= len2; ++row) {
		assignValue(mat, row, 0);
		assignValue(horizontal, row, gapOpen);
	}
	for(TSize col = 1; col <= len1; ++col) {
		TScoreValue diagValMat = getValue(mat, 0);
		assignValue(mat, 0, 0);
		vert = gapOpen;
		for(TSize row = 1; row <= len2; ++row) {
			if ((!emptyForbidden) && (getValue(forbidden, (col-1) * len2 + (row-1)) == true)) {
				tmp = 0;
				vert = 0;
				tvMat = (Byte) Stop;
				assignValue(horizontal, row, 0);
			} else {
				// Get the new maximum for vertical
				if ((tmp = getValue(mat, row - 1) + gapOpen) > vert + gap) {
					vert = tmp;
					tvVertical = (Byte) Diagonal;
				} else {
					vert = vert + gap;
					tvVertical = (Byte) Vertical;
				}
	
				// Get the new maximum for horizontal
				if ((tmp = getValue(mat, row) + gapOpen) > getValue(horizontal, row) + gap) {
					assignValue(horizontal, row, tmp);
					tvHorizontal = (Byte) Diagonal;
				} else {
					assignValue(horizontal, row, getValue(horizontal, row) + gap);
					tvHorizontal = (Byte) Horizontal;
				}
	
				// Get the new maximum for mat
				TScoreValue sc_ = score(const_cast<TScore&>(sc), str1[col-1], str2[row-1]);
				tmp = diagValMat + sc_;
				tvMat = (Byte) Diagonal;
				if (vert > tmp) {
					tmp = vert;
					tvMat = (Byte) Vertical;
				}
				if (getValue(horizontal, row) > tmp) {
					tmp = getValue(horizontal,row);
					tvMat = (Byte) Horizontal;
				}
				if (0 >= tmp) {
					tmp = 0;
					tvMat = (Byte) Stop;
				}
			}

			// Assign the new diagonal values
			diagValMat = getValue(mat, row);
			assignValue(mat, row, tmp);

			// Record the new best score
			if (tmp > score_max) {
				best_row = row;
				best_col = col;
				score_max = tmp;
				initialDir = (Byte) Diagonal;
				if (getValue(horizontal, row) ==  tmp) {
					initialDir = (Byte) Horizontal;
				}
				else if (vert == tmp) {
					initialDir = (Byte) Vertical;
				}
			}

			// Assign the right trace value
			if (tvMat == (Byte) Stop) {
				assignValue(it, 12);
			} else if (tvMat == (Byte) Diagonal) {
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

	//// Debug code
	//for(unsigned int i= 0; i<len2;++i) {
	//	for(unsigned int j= 0; j<len1;++j) {
	//		std::cout << (unsigned int) getValue(trace, j*len2 + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << "Max score: " << best_row << ',' << best_col << ':' << score_max << " (" << (unsigned int) initialDir << ")" << std::endl;

	return score_max;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScore>
inline typename Value<TScore>::Type
_localAlignment(TAlign& align,
				TStringSet const& str,
				TScore const& sc,
				SmithWaterman)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	  
	TScoreValue maxScore;
	TSize best_row = 0;
	TSize best_col = 0;
	String<bool> forbidden;

	// Trace
	String<TraceBackGotoh> trace;
	TraceBackGotoh initialDir;

	// Create the trace
	maxScore = _align_smith_waterman(trace, str, sc, initialDir, best_row, best_col, forbidden);	
	// Follow the trace and create the graph
	_align_smith_waterman_trace(align, str, trace, initialDir, best_row, best_col, forbidden);
	
	return maxScore;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScore>
inline typename Value<TScore>::Type
_localAlignment(TStringSet const& str,
				TScore const& sc,
				SmithWaterman)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TStringSet>::Type TSize;
	TraceBackGotoh initialDir;
	TSize best_row = 0;
	TSize best_col = 0;
	String<bool> forbidden;
	String<TraceBackGotoh> trace;
	return _align_smith_waterman(trace, str, sc, initialDir, best_row, best_col, forbidden);	
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
