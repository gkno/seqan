#ifndef SEQAN_HEADER_GRAPH_ALGORITHM_ALIGNMENT_H
#define SEQAN_HEADER_GRAPH_ALGORITHM_ALIGNMENT_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
//needleman wunsch alignment

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


template <typename TScoreValue, typename TTrace, typename TSpec, typename TString>
TScoreValue
needleman_wunsch(String<TScoreValue, TSpec> & mat,
				 TTrace& trace,
				 TString const & str1,
				 TString const & str2,
				 Score<TScoreValue, Simple> const& sc) 
{
	SEQAN_CHECKPOINT

	// For the DP Matrix: Where did we come from?
	enum {Diagonal = 0, Left = 1, Top=2};
	typedef String<TScoreValue, TSpec> TMatrix;
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Size<TMatrix>::Type TSize;
		
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
 
/**
.Function.needlemanWunsch:
..summary:Computes the best global alignment of the two sequences given in the StringSet of graph g.
..cat:Alignments
..signature:needlemanWunsch(g, score)
..param.g:The alignment graph having 2 sequences.
...type:Class.Graph Alignment
..param.score:The score values to be used for computing the alignment.
...type:Class.Score
..returns:The score value of the best scoring global alignment.
*/
template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue>
TScoreValue
needlemanWunsch(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				Score<TScoreValue, Simple> const& sc)
{
	SEQAN_CHECKPOINT
	// So far only simple gap score
	SEQAN_TASSERT(scoreGapOpen(sc)==0)

	// For the DP Matrix: Where did we come from?
	enum {Diagonal = 0, Left = 1, Top = 2};
	String<short> trace;

	// The DP Matrix: A simple string
	String<TScoreValue> mat;
	// Create the DP Matrix and the trace
	TScoreValue maxScore = needleman_wunsch(mat, trace, stringSet(g)[0], stringSet(g)[1], sc);	
	
	// Follow the trace and create the graph
	needleman_wunsch_trace(g, trace);	

	return maxScore;
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
