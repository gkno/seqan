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

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_PROGRESSIVE_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_PROGRESSIVE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Progressive Alignment
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TPosition, typename TSequence>
inline void 
_buildLeafString(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
				 TPosition const pos,
				 TSequence& alignSeq)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

	//TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	TStringSet& str = stringSet(g);
	TSize lenRoot = length(str[pos]);
	TId seqId = positionToId(str, pos);
	TSize i = 0;
	while(i<lenRoot) {
		TVertexDescriptor nextVertex = findVertex(const_cast<TGraph&>(g), seqId, i);
		SEQAN_TASSERT(nextVertex != getNil<TVertexDescriptor>())
		if (nextVertex == nilVertex) {
			std::cout << "Warning: Nil Vertex" << std::endl;
			TSize j = i + 1;
			while ((j < lenRoot) && (findVertex(const_cast<TGraph&>(g), seqId, j) == nilVertex)) ++j;
			nextVertex = addVertex(const_cast<TGraph&>(g), seqId, i, j-i);
		}
		appendValue(alignSeq, String<TVertexDescriptor>(nextVertex));
		i += fragmentLength(g, nextVertex);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TSegmentString, typename TOutGraph>
inline void 
_createAlignmentGraph(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
					  TSegmentString& alignSeq,
					  TOutGraph& gOut)			 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef String<TVertexDescriptor> TVertexString;

	// Initialization
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

	// Create the alignment graph
	TSize alignSeqLen = length(alignSeq);
	for(TSize i = 0; i<alignSeqLen;++i) {
		TVertexString& alignSeq_i = alignSeq[i];
		TSize len_i = length(alignSeq_i);
		for(TSize j=0; j<len_i; ++j) {
			TVertexDescriptor v = getValue(alignSeq_i, j);
			if (v == nilVertex) continue;
			SEQAN_TASSERT(fragmentBegin(g,v) < length(getValueById(stringSet(g), sequenceId(g,v))))
			SEQAN_TASSERT(fragmentLength(g,v) > 0)
			SEQAN_TASSERT(fragmentBegin(g,v) + fragmentLength(g,v) <= length(getValueById(stringSet(g), sequenceId(g,v))))
			TVertexDescriptor l = addVertex(gOut, sequenceId(g, v), fragmentBegin(g,v), fragmentLength(g,v));
			//std::cout << l << label(gOut, l) << ',';
			TSize count = 1;
			for(TSize k = j; k>0; --k) {
				if (getValue(alignSeq_i, k - 1) != nilVertex) {
					SEQAN_TASSERT(fragmentLength(gOut,l) == fragmentLength(gOut,l - count))
					addEdge(gOut, (TVertexDescriptor) (l - count), (TVertexDescriptor) l);
					++count;
				}
			}
		}
		//std::cout << std::endl;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TVertexDescriptor, typename TSequence>
inline void 
_recursiveProgressiveAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
							   TGuideTree& tree,
							   TVertexDescriptor const root,
							   TSequence& alignSeq)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjacencyIterator;
	typedef typename Iterator<TGuideTree, DfsPreorder>::Type TDfsPreorderIterator;

	if(isLeaf(tree, root)) {
		_buildLeafString(g, root, alignSeq);
	} else {
		// Align the two children (Binary tree)
		typedef String<String<TVertexDescriptor> > TSegmentString;
		TSegmentString seq1;
		TSegmentString seq2;
		TAdjacencyIterator adjIt(tree, root);
		_recursiveProgressiveAlignment(g,tree, *adjIt, seq1);
		goNext(adjIt);
		_recursiveProgressiveAlignment(g,tree, *adjIt, seq2);
		heaviestCommonSubsequence(g,seq1,seq2,alignSeq);

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

/**
.Function.progressiveAlignment:
..summary:Performs a progressive alignment.
..cat:Graph
..signature:
combineGraphs(inputGraph, guideTree, outputGraph [,seqsPerGroup])
..param.graph:The alignment graph with multiple sequence information.
...type:Spec.Alignment Graph
..param.guideTree:A guide tree.
...type:Spec.Tree
..param.outputGraph:An alignment graph for final MSA.
...type:Spec.Alignment Graph
..param.seqsPerGroup:Integer indicating the amount of sequences per group.
...remarks:The sequence set is split into groups of seqsPerGroup members.
Then a triplet extension is performed solely on the group members.
Thus, you do not need to do a separate triplet extension if you use this option!
..returns:void
*/
template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TOutGraph>
inline void 
progressiveAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					 TGuideTree& tree,
					 TOutGraph& gOut)			 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef String<TVertexDescriptor> TVertexString;
	typedef String<TVertexString> TSegmentString;

	// Perform progressive alignment
	TSegmentString alignSeq;
	_recursiveProgressiveAlignment(g,tree,getRoot(tree),alignSeq);

	// Create the alignment graph
	_createAlignmentGraph(g, alignSeq, gOut);
}


				

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TOutGraph, typename TSize>
inline void
progressiveAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					 TGuideTree& tree,
					 TOutGraph& gOut,
					 TSize seqPerGroup)			 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Id<TStringSet>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGuideTree, DfsPreorder>::Type TDfsPreorderIterator;
	typedef typename Iterator<TGuideTree, BfsIterator>::Type TBfsIterator;
	typedef String<TVertexDescriptor> TVertexString;
	typedef String<TVertexString> TSegmentString;
	typedef std::set<TSize> TChildrenSet;

	// Build groups of sequences
	TGuideTree copy_tree(tree);
	String<std::set<TSize> > sequenceGroups;
	String<TId> sequenceGroupRoots;
	while (numChildren(copy_tree, getRoot(copy_tree)) > 0) {
		TBfsIterator bfsIt(copy_tree, getRoot(copy_tree));
		for(;!atEnd(bfsIt);goNext(bfsIt)) {
			TChildrenSet sequenceSet;
			collectLeaves(copy_tree, *bfsIt, sequenceSet);
			if ((TSize) sequenceSet.size() <= (TSize) seqPerGroup) {
				TChildrenSet finalSequenceSet;
				for(typename TChildrenSet::iterator pos = sequenceSet.begin(); pos != sequenceSet.end(); ++pos) {
					finalSequenceSet.insert(positionToId(stringSet(g), *pos));
				}
				appendValue(sequenceGroups, finalSequenceSet);
				// Find the proper root
				TDfsPreorderIterator dfsItRootFinder(copy_tree, *bfsIt);
				for(;!atEnd(dfsItRootFinder);goNext(dfsItRootFinder)) {
					if (numChildren(copy_tree, *dfsItRootFinder) == 1) continue;
					appendValue(sequenceGroupRoots, *dfsItRootFinder);
					break;
				}
				// Cut the tree
				if (getRoot(copy_tree) != *bfsIt) {
					removeChild(copy_tree, parentVertex(tree, *bfsIt), *bfsIt);
				} else {
					removeAllChildren(copy_tree, *bfsIt);
				}
				break;
			}
		}
	}
	clear(copy_tree);

	//// Debug code
	//for(TSize i = 0;i<length(sequenceGroups);++i) {
	//	std::cout << "Sequences: ";
	//	copy((sequenceGroups[i]).begin(), (sequenceGroups[i]).end(), std::ostream_iterator<int>(std::cout, " "));
	//	std::cout << std::endl;
	//}

	// Align the sequence groups and flip the order for profile alignment
	String<TSegmentString> profileStrings;
	resize(profileStrings, length(sequenceGroups));
	TSize numGroups = length(sequenceGroups);
	for(TSize i = 0;i<numGroups;++i) {
		// Perform progressive alignment
		TGraph copy_graph(g);
		if ((sequenceGroups[i]).size() > 1) tripletLibraryExtension(copy_graph, sequenceGroups[i]);
		TSegmentString alignSeq;
		_recursiveProgressiveAlignment(copy_graph,tree,sequenceGroupRoots[i],alignSeq);
		profileStrings[i] = alignSeq;
		//std::cout << "One mini tree finished " << (sequenceGroups[i]).size() << std::endl;
		clear(copy_graph);
	}

	// Align the profiles
	TSize lastIndex = numGroups - 1;
	for(TSize y = 1; y < numGroups;++y) {
		//TGraph copy_graph(g);
		//restrictedTripletLibraryExtension(copy_graph, sequenceGroups[lastIndex], sequenceGroups[lastIndex-y]);
		TSegmentString alignSeq;
		//heaviestCommonSubsequence(copy_graph,profileStrings[lastIndex],profileStrings[lastIndex-y],alignSeq);
		heaviestCommonSubsequence(g,profileStrings[lastIndex],profileStrings[lastIndex-y],alignSeq);
		profileStrings[lastIndex] = alignSeq;
		(sequenceGroups[lastIndex]).insert((sequenceGroups[lastIndex-y]).begin(),(sequenceGroups[lastIndex-y]).end());
		//std::cout << "Two mini trees joined " << (sequenceGroups[lastIndex]).size() << std::endl;
		//clear(copy_graph);
	}

	// Create the alignment graph
	_createAlignmentGraph(g, profileStrings[lastIndex], gOut);
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
