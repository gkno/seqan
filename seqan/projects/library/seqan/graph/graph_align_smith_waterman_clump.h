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
	typedef typename Id<TAlign>::Type TId;
	typedef typename Cargo<TAlign>::Type TCargo;
	typedef typename Iterator<TAlign, EdgeIterator>::Type TEdgeIterator;
  
	TScoreValue maxScore = 0;
	TSize minLen = length(str[0]);
	TSize best_row = 0;
	TSize best_col = 0;
	String<bool> forbidden;
	fill(forbidden, length(str[1])* length(str[0]), false);
	TSize tmp;
	if ((tmp = length(str[1])) < minLen) minLen = tmp;

	// String of fragments to combine all pairwise alignments into a multiple alignment
	typedef Fragment<> TFragment;
	typedef String<TFragment, External<> > TFragmentString;
	typedef Iterator<TFragmentString> TFragmentStringIter;
	TFragmentString matches;
	String<TScoreValue> scores;


	minLen /= 20;
	if (minLen < 2) minLen = 2;
	for(TSize i = 0; i<minLen; ++i) {
		typedef Graph<Alignment<TStringSet, void> > TPairGraph;
		typedef typename VertexDescriptor<TPairGraph>::Type TVD;
		typedef typename Iterator<TPairGraph, EdgeIterator>::Type TEI;
		
		// Trace
		String<TraceBackGotoh> trace;
		TraceBackGotoh initialDir;

		TPairGraph pGraph(str);

		// Create the trace
		TScoreValue tmpScore = _align_smith_waterman(trace, str, sc, initialDir, best_row, best_col, forbidden);
		if (i == 0) maxScore = tmpScore;
		
		// Follow the trace and create the graph
		_align_smith_waterman_trace(pGraph, str, trace, initialDir, best_row, best_col, forbidden);


		TEI it(pGraph);
		for(;!atEnd(it);++it) {
			TVD sV = sourceVertex(it);
			TVD tV = targetVertex(it);
			push_back(matches, TFragment( (unsigned int) sequenceId(pGraph, sV), (unsigned int) fragmentBegin(pGraph,sV), (unsigned int) sequenceId(pGraph, tV),  (unsigned int)  fragmentBegin(pGraph,tV),  (unsigned int)  fragmentLength(pGraph,tV)));
			appendValue(scores, tmpScore);
		}

		//// Debug code
		//for(unsigned int i= 0; i<length(str[1]);++i) {
		//	for(unsigned int j= 0; j<length(str[0]);++j) {
		//		std::cout << (unsigned int) getValue(forbidden, j*length(str[1]) + i) << ',';
		//	}
		//	std::cout << std::endl;
		//}
		//std::cout << std::endl;
	}

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
	matchRefinement(matches,str,const_cast<Score<TScoreValue, TSpec>&>(sc),align);
	return maxScore;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
