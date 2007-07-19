#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_PROGRESSIVE_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_PROGRESSIVE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// T-Coffee - Progressive Alignment
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TSequenceIn, typename TSequence, typename TTag>
inline void 
_alignStringsAccordingToGraph(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
							  TSequenceIn const& seq1,
							  TSequenceIn const& seq2,
							  TSequence& alignSeq,
							  TTag)
{
	SEQAN_CHECKPOINT
	// Clear output parameter
	clear(alignSeq);

	// Just heaviest common subsequence
	heaviestCommonSubsequence(g, seq1, seq2, alignSeq);

	//// Debug Code
	//std::cout << "The sequences:" << std::endl;
	//for(unsigned int i = 0;i<length(stringSet(g));++i) {
	//	std::cout << positionToId(stringSet(g),i) << ':' << (stringSet(g))[i] << std::endl;
	//}
	//for(unsigned int i = 0; i<length(alignSeq);++i) {
	//	std::cout << '(';
	//	for(unsigned int j=0; j<length(alignSeq[i]);++j) {
	//		TVertexDescriptor v = getValue(alignSeq[i], j);
	//		std::cout << v << ',';
	//		if (v != getNil<TVertexDescriptor>()) {
	//			std::cout << label(g, getValue(alignSeq[i], j)) << ',';
	//		}
	//	}
	//	std::cout << ')' << ',';
	//}
	//std::cout << std::endl;

}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TVertexDescriptor, typename TSequence, typename TTag>
inline void 
_recursiveProgressiveAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
							   TGuideTree& tree,
							   TVertexDescriptor const root,
							   TSequence& alignSeq,
							   TTag)
{
	SEQAN_CHECKPOINT

	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjacencyIterator;
	typedef typename Iterator<TGuideTree, DfsPreorder>::Type TDfsPreorderIterator;

	TStringSet& str = stringSet(g);
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

	if(isLeaf(tree, root)) {
		TId seqId = positionToId(str, root);
		TSize i = 0;
		while(i<length(str[root])) {
			TVertexDescriptor nextVertex = findVertex(g, seqId, i);
			if (nextVertex == nilVertex) {
				//std::cout << "Nil vertex" << std::endl;
				TSize j = i + 1;
				while ((j < length(str[root])) && (findVertex(g, seqId, j) == nilVertex)) ++j;
				nextVertex = addVertex(g, seqId, i, j-i);
			}
			appendValue(alignSeq, String<TVertexDescriptor>(nextVertex));
			i += fragmentLength(g, nextVertex);
		}
	} else {
		// Align the two children
		typedef String<String<TVertexDescriptor> > TSegmentString;
		TSegmentString seq1;
		TSegmentString seq2;
		TAdjacencyIterator adjIt(tree, root);
		_recursiveProgressiveAlignment(g,tree, *adjIt, seq1, TTag());
		goNext(adjIt);
		_recursiveProgressiveAlignment(g,tree, *adjIt, seq2, TTag());
		_alignStringsAccordingToGraph(g,seq1,seq2,alignSeq, TTag());

		//// Debug Code
		//for(unsigned int i = 0; i<length(alignSeq);++i) {
		//	std::cout << '(';
		//	for(unsigned int j=0; j<length(alignSeq[i]);++j) {
		//		std::cout << getValue(alignSeq[i], j) << ',';
		//	}
		//	std::cout << ')' << ',';
		//}
		//std::cout << std::endl;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TOutGraph, typename TTag>
inline void 
progressiveAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					 TGuideTree& tree,
					 TOutGraph& gOut,
					 TTag)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjacencyIterator;
	typedef typename Iterator<TGuideTree, DfsPreorder>::Type TDfsPreorderIterator;

	// Perform progressive alignment
	typedef String<String<TVertexDescriptor> > TSegmentString;
	TSegmentString alignSeq;
	_recursiveProgressiveAlignment(g,tree,getRoot(tree),alignSeq, TTag());

	// Create the alignment graph
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	for(TSize i = 0; i<length(alignSeq);++i) {
		for(TSize j=0; j<length(alignSeq[i]);++j) {
			TVertexDescriptor v = getValue(alignSeq[i], j);
			if (v == nilVertex) continue;
			SEQAN_TASSERT(fragmentBegin(g,v) < length(getValueById(stringSet(g), sequenceId(g,v))))
			SEQAN_TASSERT(fragmentLength(g,v) > 0)
			SEQAN_TASSERT(fragmentBegin(g,v) + fragmentLength(g,v) <= length(getValueById(stringSet(g), sequenceId(g,v))))
			TVertexDescriptor l = addVertex(gOut, sequenceId(g, v), fragmentBegin(g,v), fragmentLength(g,v));
			//std::cout << l << label(gOut, l) << ',';
			TSize count = 1;
			for(TSize k = j; k>0; --k) {
				if (getValue(alignSeq[i], k - 1) != nilVertex) {
					SEQAN_TASSERT(fragmentLength(gOut,l) == fragmentLength(gOut,l - count))
					addEdge(gOut, l - count, l);
					++count;
				}
			}
		}
		//std::cout << std::endl;
	}

}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
