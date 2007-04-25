#ifndef SEQAN_HEADER_GRAPH_ALGORITHM_ALIGNMENT_H
#define SEQAN_HEADER_GRAPH_ALGORITHM_ALIGNMENT_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Needleman-Wunsch Alignment, constant gap cost
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////


template <typename TStringSet, typename TCargo, typename TSpec, typename TTrace>
void
needleman_wunsch_trace(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TTrace const& trace)
{
	SEQAN_CHECKPOINT
	enum {Diagonal = 0, Left = 1, Top=2};
	
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Id<TStringSet>::Type TId;
	 
	TStringSet* str = &stringSet(g);
	TId id1 = positionToId(*str, 0);
	TId id2 = positionToId(*str, 1);
	TSize len1 = length(((*str)[0]));
	TSize cols = len1 + 1;
	TSize len2 = length(((*str)[1]));

	// Initialize everything
	TTraceValue tv = getValue(trace, len2*cols + len1);
	TTraceValue tvOld = tv;  // We need to know when the direction of the trace changes

	TSize segLen = 1;
	switch(tv) {
		case Diagonal:
			//std::cout << '(' << ((*str)[0])[len1 - 1] << ',' << ((*str)[1])[len2-1] << ')' << std::endl;
			--len1; --len2;
			break;
		case Left:
			//std::cout << '(' << ((*str)[0])[len1 - 1] << ',' << '-' << ')' << std::endl;
			--len1;
			break;
		case Top:
			//std::cout << '(' << '-' << ',' << ((*str)[1])[len2-1] << ')' << std::endl;
			--len2;
			break;
	}
	// Now follow the trace
	do {
		tv = getValue(trace, len2*cols + len1);
		switch(tv) {
			case Diagonal: 
				//std::cout << '(' << ((*str)[0])[len1 - 1] << ',' << ((*str)[1])[len2-1] << ')' << std::endl;
				if (tv != tvOld) {
					switch(tvOld) {
						case Left:	
							addVertex(g, id1, len1, segLen);
							break;
						case Top:
							addVertex(g, id2, len2, segLen);
							break;
					}
					tvOld = tv; segLen = 1;
				} else ++segLen;
				--len1; --len2;
				break;
			case Left:
				//std::cout << '(' << ((*str)[0])[len1 - 1] << ',' << '-' << ')' << std::endl;
				if (tv != tvOld) {
					switch(tvOld) {
						case Diagonal: 
							addEdge(g, addVertex(g, id1, len1, segLen), addVertex(g, id2, len2, segLen));	
							break;
						case Top:	
							addVertex(g, id2, len2, segLen);
							break;
					}
					tvOld = tv; segLen = 1;
				} else ++segLen;
				--len1;
				break;
			case Top:
				//std::cout << '(' << '-' << ',' << ((*str)[1])[len2-1] << ')' << std::endl;
				if (tv != tvOld) {
					switch(tvOld) {
						case Diagonal: 
							addEdge(g, addVertex(g, id1, len1, segLen), addVertex(g, id2, len2, segLen));
							break;
						case Left:
							addVertex(g, id1, len1, segLen);
							break;
					}
					tvOld = tv; segLen = 1;
				} else ++segLen;
				--len2;
				break;
		}
	} while ((len1 != 0) && (len2 !=0));
	// Process left-overs
	switch(tvOld) {
		case Diagonal: 
			addEdge(g, addVertex(g, id1, len1, segLen), addVertex(g, id2, len2, segLen));
			break;
		case Left:
			addVertex(g, id1, len1, segLen);
			break;
		case Top:
			addVertex(g, id2, len2, segLen);
			break;
	}
	// Handle the remaining sequence
	if (len1 != 0) {
		addVertex(g, id1, 0, len1);
	} else if (len2 != 0) {
		addVertex(g, id2, 0, len2);
	} 
}

//////////////////////////////////////////////////////////////////////////////

template <typename TScoreValue, typename TValue, typename TSpec, typename TString>
TScoreValue
needleman_wunsch(String<TValue, TSpec>& trace,
				 TString const & str1,
				 TString const & str2,
				 Score<TScoreValue, Simple> const& sc) 
{
	SEQAN_CHECKPOINT

	// For the DP Matrix: Where did we come from?
	enum {Diagonal = 0, Left = 1, Top=2};
	typedef String<TValue, TSpec> TTrace;
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Size<String<TScoreValue> >::Type TSize;

	// The DP Matrix: A simple string
	String<TScoreValue> mat;
		
	TSize len1 = length(str1);
	TSize cols = len1 + 1;
	TSize len2 = length(str2);
	TScoreValue gap = scoreGapExtend(sc);
	TScoreValue maxVal = 0;
	TScoreValue tmp = 0;
	resize(mat, (len1+1)*(len2+1));   // One additional row and column
	resize(trace, (len1+1)*(len2+1));

	// Classical DP
	TTraceValue tv;
	for(TSize col = 0; col <= len1; ++col) assignValue(mat, col, col * gap);
	for(TSize row = 0; row < len2; ++row) {
		assignValue(mat, (row + 1) * cols, (row+1) * gap);
		for(TSize col = 0; col < len1; ++col) {
			// Get the new maximum	
			maxVal = getValue(mat, row*cols+col) + score(const_cast<Score<TScoreValue, Simple>&>(sc), str1[col], str2[row]);
			tv = Diagonal;
			if ((tmp = getValue(mat, row * cols + (col + 1)) + gap) > maxVal) {
				maxVal = tmp;
				tv = Top;
			}
			if ((tmp = getValue(mat, (row + 1) * cols + col) + gap) > maxVal) {
				maxVal = tmp;
				tv = Left;
			}
			// Assign the new value
			assignValue(mat, (row+1) * cols + (col+1), maxVal);
			assignValue(trace, (row+1) * cols + (col+1), tv);
		}
	}
	return getValue(mat, len2*cols + len1);
}

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Alignment: Gotoh Alignment, affine gap cost
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet, typename TCargo, typename TSpec, typename TTrace, typename TVal>
void
gotoh_trace(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
			TTrace const& trace,
			TVal const initialDir)
{
	SEQAN_CHECKPOINT
	enum {Diagonal = 0, Left = 1, Top=2};
	
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Id<TStringSet>::Type TId;
	 
	TStringSet* str = &stringSet(g);
	TId id1 = positionToId(*str, 0);
	TId id2 = positionToId(*str, 1);
	TSize len1 = length(((*str)[0]));
	TSize cols = len1 + 1;
	TSize len2 = length(((*str)[1]));

	// Initialize everything
	
	TTraceValue tv = 0;
	TTraceValue tvOld = initialDir;
	switch(initialDir) {
		case Diagonal:
			tv = ((getValue(trace, len2*cols + len1) / 9) % 3);
			break;
		case Left:
			tv = ((getValue(trace, len2*cols + len1) / 3) % 3);
			break;
		case Top:
			tv = (getValue(trace, len2*cols + len1) % 3);
			break;
	}
	TSize segLen = 0;

	// Now follow the trace
	do {
		switch(tvOld) {
			case Diagonal: 
				//std::cout << '(' << ((*str)[0])[len1 - 1] << ',' << ((*str)[1])[len2-1] << ')' << std::endl;
				--len1; --len2;
				++segLen;
				if (tv != tvOld) {
					addEdge(g, addVertex(g, id1, len1, segLen), addVertex(g, id2, len2, segLen));
					tvOld = tv; 
					segLen = 0;
				}
				break;
			case Left:
				//std::cout << '(' << ((*str)[0])[len1 - 1] << ',' << '-' << ')' << std::endl;
				--len1;
				++segLen;
				if (tv != tvOld) {
					addVertex(g, id1, len1, segLen);
					tvOld = tv;
					segLen = 0;
				}
				break;
			case Top:
				//std::cout << '(' << '-' << ',' << ((*str)[1])[len2-1] << ')' << std::endl;
				--len2;
				++segLen;
				if (tv != tvOld) {
					addVertex(g, id2, len2, segLen);
					tvOld = tv; 
					segLen = 0;
				}
				break;
		}
		switch(tv) {
			case Diagonal:
				tv = ((getValue(trace, len2*cols + len1) / 9) % 3);
				break;
			case Left:
				tv = ((getValue(trace, len2*cols + len1) / 3) % 3);
				break;
			case Top:
				tv = (getValue(trace, len2*cols + len1) % 3);
				break;
		}
	} while ((len1 != 0) && (len2 !=0));
	// Process left-overs
	if (segLen > 0) {
		switch(tvOld) {
			case Diagonal: 
				addEdge(g, addVertex(g, id1, len1, segLen), addVertex(g, id2, len2, segLen));
				break;
			case Left:
				addVertex(g, id1, len1, segLen);
				break;
			case Top:
				addVertex(g, id2, len2, segLen);
				break;
		}
	}
	// Handle the remaining sequence
	if (len1 != 0) {
		addVertex(g, id1, 0, len1);
	} else if (len2 != 0) {
		addVertex(g, id2, 0, len2);
	} 
}

//////////////////////////////////////////////////////////////////////////////

template <typename TScoreValue, typename TValue, typename TSpec, typename TString>
TScoreValue
gotoh(String<TValue, TSpec>& trace,
	  TString const & str1,
	  TString const & str2,
	  Score<TScoreValue, Simple> const & sc,
	  TValue& initialDir)
{
SEQAN_CHECKPOINT
	// DiagLeftTop means diagonal in mat matrix, left in left matrix and left in top matrix
	enum {Diagonal = 0, Left = 1, Top=2};
	typedef String<TValue, TSpec> TTrace;
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Size<String<TScoreValue> >::Type TSize;

	// The DP Matrix for diagonal walks
	String<TScoreValue> mat;
	// The DP Matrix for gaps from the left
	String<TScoreValue> left;
	// The DP Matrix for gaps from the top
	String<TScoreValue> top;
		
	TSize len1 = length(str1);
	TSize cols = len1 + 1;
	TSize len2 = length(str2);
	TScoreValue gap = scoreGapExtend(sc);
	TScoreValue gapOpen = scoreGapOpen(sc);
	TScoreValue maxVal = 0;
	TScoreValue tmp = 0;
	resize(mat, (len1+1)*(len2+1));   // One additional row and column
	resize(left, (len1+1)*(len2+1));   // One additional row and column
	resize(top, (len1+1)*(len2+1));   // One additional row and column
	resize(trace, (len1+1)*(len2+1));

	TTraceValue tvMat, tvLeft, tvTop;
	
	TScoreValue inf = infimumValue<TScoreValue>();

	// Classical DP
	assignValue(mat, 0, 0);
	assignValue(left, 0, inf);
	assignValue(top, 0, inf);
	for(TSize col = 0; col < len1; ++col) {
		assignValue(mat, (col + 1), gapOpen + col * gap);
		assignValue(top, (col + 1), gapOpen + col * gap);
		assignValue(left, (col + 1), inf);
	}
	for(TSize row = 0; row < len2; ++row) {
		assignValue(mat, (row + 1) * cols, gapOpen + row * gap);
		assignValue(left, (row + 1) * cols, gapOpen + row * gap);
		assignValue(top, (row + 1) * cols, inf);
		for(TSize col = 0; col < len1; ++col) {
			// Get the new maximum for top
			maxVal = getValue(mat, row * cols + (col + 1)) + gapOpen;
			tvTop = Diagonal;
			if ((tmp = getValue(top, row * cols + (col + 1)) + gap) > maxVal) {
				maxVal = tmp;
				tvTop = Top;
			}
			assignValue(top, (row+1) * cols + (col+1), maxVal);

			// Get the new maximum for left
			maxVal = getValue(mat, (row + 1) * cols + col) + gapOpen;
			tvLeft = Diagonal;
			if ((tmp = getValue(left, (row + 1) * cols + col) + gap) > maxVal) {
				maxVal = tmp;
				tvLeft = Left;
			}
			assignValue(left, (row+1) * cols + (col+1), maxVal);

			// Get the new maximum for mat
			TScoreValue sc_ = score(const_cast<Score<TScoreValue, Simple>&>(sc), str1[col], str2[row]);
			maxVal = getValue(mat, row*cols+col) + sc_;
			tvMat = Diagonal;
			if ((tmp = getValue(top, row * cols + col) + sc_) > maxVal) {
				maxVal = tmp;
				tvMat = Top;
			}
			if ((tmp = getValue(left, row * cols + col) + sc_) > maxVal) {
				maxVal = tmp;
				tvMat = Left;
			}
			// Assign the new value
			assignValue(mat, (row+1) * cols + (col+1), maxVal);
			assignValue(trace, (row+1) * cols + (col+1), tvMat * 3 * 3 + tvLeft * 3 + tvTop);
		}

	}

	maxVal = getValue(mat, len2*cols + len1);
	initialDir = Diagonal;
	if ((tmp = getValue(left, len2 * cols + len1)) > maxVal) {
		maxVal = tmp;
		initialDir = Left;
	}
	if ((tmp = getValue(top, len2 * cols + len1)) > maxVal) {
		maxVal = tmp;
		initialDir = Top;
	}


	//for(unsigned int i= 0; i<=len2;++i) {
	//	for(unsigned int j= 0; j<=len1;++j) {
	//		std::cout << getValue(mat, i*(len1 + 1) + j) << ',';
	//	}
	//	std::cout << std::endl;
	//}

	//for(unsigned int i= 0; i<=len2;++i) {
	//	for(unsigned int j= 0; j<=len1;++j) {
	//		std::cout << getValue(left, i*(len1 + 1) + j) << ',';
	//	}
	//	std::cout << std::endl;
	//}

	//for(unsigned int i= 0; i<=len2;++i) {
	//	for(unsigned int j= 0; j<=len1;++j) {
	//		std::cout << getValue(top, i*(len1 + 1) + j) << ',';
	//	}
	//	std::cout << std::endl;
	//}
	//
	//for(unsigned int i= 0; i<=len2;++i) {
	//	for(unsigned int j= 0; j<=len1;++j) {
	//		std::cout << (getValue(trace, i*(len1 + 1) + j) % 3) << ',';
	//	}
	//	std::cout << std::endl;
	//}

	//for(unsigned int i= 0; i<=len2;++i) {
	//	for(unsigned int j= 0; j<=len1;++j) {
	//		std::cout << ((getValue(trace, i*(len1 + 1) + j) / 3) % 3) << ',';
	//	}
	//	std::cout << std::endl;
	//}

	//for(unsigned int i= 0; i<=len2;++i) {
	//	for(unsigned int j= 0; j<=len1;++j) {
	//		std::cout << ((getValue(trace, i*(len1 + 1) + j) / 9) % 3) << ',';
	//	}
	//	std::cout << std::endl;
	//}

	return maxVal;
}
 
/**
.Function.globalAlignment:
..summary:Computes the best global alignment of the two sequences given in the StringSet of graph g.
..cat:Alignments
..signature:globalAlignment(g, score)
..param.g:The alignment graph having 2 sequences.
...type:Class.Graph Alignment
..param.score:The score values to be used for computing the alignment.
...type:Class.Score
..returns:The score value of the best scoring global alignment.
*/
template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue>
TScoreValue
globalAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				Score<TScoreValue, Simple> const& sc)
{
	SEQAN_CHECKPOINT
	// So far only simple gap score
	TScoreValue maxScore;
	if(scoreGapOpen(sc)==0) {	
		// Needleman wunsch	
		String<short> trace;
		maxScore = needleman_wunsch(trace, stringSet(g)[0], stringSet(g)[1], sc);	
	
		// Follow the trace and create the graph
		needleman_wunsch_trace(g, trace);	
	} else {
		// Gotoh
		String<short> trace;
		short initialDir;
		maxScore = gotoh(trace, stringSet(g)[0], stringSet(g)[1], sc, initialDir);	

		// Follow the trace and create the graph
		gotoh_trace(g, trace, initialDir);
	}

	return maxScore;
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
