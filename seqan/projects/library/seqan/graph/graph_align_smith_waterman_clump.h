#ifndef SEQAN_HEADER_GRAPH_SMITH_WATERMAN_CLUMP_H
#define SEQAN_HEADER_GRAPH_SMITH_WATERMAN_CLUMP_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Smith Waterman Alignment with Clumping, affine gap cost
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScoreValue, typename TSpec>
TScoreValue
_localAlignment(TAlign& align,
				TStringSet& str,
				Score<TScoreValue, TSpec> const& sc,
				SmithWatermanClump)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename VertexDescriptor<TAlign>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TAlign>::Type TEdgeDescriptor;
	typedef typename Id<TAlign>::Type TId;
	typedef typename Cargo<TAlign>::Type TCargo;
	typedef typename Iterator<TAlign, EdgeIterator>::Type TEdgeIterator;
  
	// Score of the best alignment
	TScoreValue maxScore = 0;

	// For clumpping remember the used positions
	TSize len0 = length(str[0]);
	TSize len1 = length(str[1]);
	String<bool> forbidden;
	fill(forbidden, len0 * len1, false);

	// Get the length of the shortest string
	TSize minLen = len0;
	if (len1  < minLen) minLen = len1;
	minLen = (TSize) 0.1 * minLen;
	// Local matches must be at least 20 characters long
	if (minLen < 20) minLen = 20;

	// String of fragments to combine all found local alignments into one alignment graph
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;
	String<TScoreValue, Block<> > score_values;

	// Stop looking for local alignments, if there are to short
	TSize local_len;
	TSize count = 0;
	do {
		local_len = 0;
		typedef Graph<Alignment<TStringSet, void> > TPairGraph;
		typedef typename VertexDescriptor<TPairGraph>::Type TVD;
		typedef typename Iterator<TPairGraph, EdgeIterator>::Type TEI;
		
		// Trace
		String<TraceBackGotoh> trace;	// Trace-back path
		TraceBackGotoh initialDir;		// First direction on the path
		TSize best_row = 0;				// Where to start
		TSize best_col = 0;

		// The local alignment graph
		TPairGraph pGraph(str);

		// Create the local alignment
		TScoreValue tmpScore = _align_smith_waterman(trace, str, sc, initialDir, best_row, best_col, forbidden);
		if (tmpScore > maxScore) maxScore = tmpScore;
		
		// Follow the trace and create the graph
		_align_smith_waterman_trace(pGraph, str, trace, initialDir, best_row, best_col, forbidden);
		//std::cout << pGraph << std::endl;

		// Extract the matches
		TEI it(pGraph);
		for(;!atEnd(it);++it) {
			TVD sV = sourceVertex(it);
			TVD tV = targetVertex(it);
			push_back(matches, TFragment( (unsigned int) sequenceId(pGraph, sV), (unsigned int) fragmentBegin(pGraph,sV), (unsigned int) sequenceId(pGraph, tV),  (unsigned int)  fragmentBegin(pGraph,tV),  (unsigned int)  fragmentLength(pGraph,tV)));
			push_back(score_values, tmpScore);
			local_len += fragmentLength(pGraph,tV);
		}

		//// Debug code
		//for(unsigned int i= 0; i<length(str[1]);++i) {
		//	for(unsigned int j= 0; j<length(str[0]);++j) {
		//		std::cout << (unsigned int) getValue(forbidden, j*length(str[1]) + i) << ',';
		//	}
		//	std::cout << std::endl;
		//}
		//std::cout << std::endl;
		++count;
	} while ((local_len > minLen) && (count < 10));

	//// Debug Code
	//// Print all the matches
	//std::cout << "The sequences:" << std::endl;
	//for(TSize i = 0;i<length(str);++i) {
	//	std::cout << positionToId(str,i) << ':' << str[i] << std::endl;
	//}
	//std::cout << "The matches:" << std::endl;
	//for(TSize i = 0;i<length(matches);++i) {
	//	TId tmp_id1 = sequenceId(matches[i],0);
	//	std::cout << tmp_id1 << ',' << fragmentBegin(matches[i],tmp_id1) << ',';
	//	for(TSize j = fragmentBegin(matches[i],tmp_id1); j < fragmentBegin(matches[i],tmp_id1) + fragmentLength(matches[i],tmp_id1); ++j) {
	//		std::cout << str[idToPosition(str, tmp_id1)][j];
	//	}
	//	TId tmp_id2 = sequenceId(matches[i],1);
	//	std::cout << ',' <<	tmp_id2 << ',' << fragmentBegin(matches[i],tmp_id2) << ',';
	//	for(TSize j = fragmentBegin(matches[i],tmp_id2); j < fragmentBegin(matches[i],tmp_id2) + fragmentLength(matches[i],tmp_id2); ++j) {
	//		std::cout << str[idToPosition(str, tmp_id2)][j];
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << "=====" << std::endl;
	
	// Refine all matches and create multiple alignment
	matchRefinement(matches,str,align);

	// Adapt edge weights
	TFragmentStringIter endIt = end(matches);
	for(TFragmentStringIter it = begin(matches); it != endIt; ++it) {
		TId id1 = sequenceId(*it,0);
		TId id2 = sequenceId(*it,1);
		TSize pos1 = fragmentBegin(*it, id1);
		TSize pos2 = fragmentBegin(*it, id2);
		TSize end1 = pos1 + fragmentLength(*it, id1);
		while(pos1 < end1) {
			TVertexDescriptor p1 = findVertex(align, id1, pos1);
			TVertexDescriptor p2 = findVertex(align, id2, pos2);
			TEdgeDescriptor e = findEdge(align, p1, p2);
			cargo(e) = (TCargo) getValue(score_values, position(it));
			SEQAN_TASSERT(fragmentLength(align, p1) == fragmentLength(align, p2))
			pos1 += fragmentLength(align, p1);
			pos2 += fragmentLength(align, p2);
		}
	}

	return maxScore;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
