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
  $Id: graph_align_tcoffee_progressive.h 1919 2008-05-02 15:54:46Z rausch@PCPOOL.MI.FU-BERLIN.DE $
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
	typedef typename Value<TSequence>::Type TVertexString;

	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	TStringSet& str = stringSet(g);
	TSize lenRoot = length(value(str, pos));
	TId seqId = positionToId(str, pos);
	TSize i = 0;
	while(i<lenRoot) {
		TVertexDescriptor nextVertex = findVertex(const_cast<TGraph&>(g), seqId, i);
		SEQAN_TASSERT(nextVertex != nilVertex)
		if (nextVertex == nilVertex) {
			std::cout << "Warning: Nil Vertex" << std::endl;
			TSize j = i + 1;
			while ((j < lenRoot) && (findVertex(const_cast<TGraph&>(g), seqId, j) == nilVertex)) ++j;
			nextVertex = addVertex(const_cast<TGraph&>(g), seqId, i, j-i);
		}
		appendValue(alignSeq, TVertexString(nextVertex));
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
		TVertexString& alignSeq_i = value(alignSeq, i);
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

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TVertexDescriptor, typename TSegmentString>
inline void 
_recursiveProgressiveAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
							   TGuideTree& tree,
							   TVertexDescriptor const root,
							   TSegmentString& alignSeq)
{
	SEQAN_CHECKPOINT
	typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjacencyIterator;

	if(isLeaf(tree, root)) {
		_buildLeafString(g, root, alignSeq);
	} else {
		// Align the two children (Binary tree)
		TSegmentString seq1;
		TSegmentString seq2;
		TAdjacencyIterator adjIt(tree, root);
		_recursiveProgressiveAlignment(g,tree, *adjIt, seq1);
		goNext(adjIt);
		_recursiveProgressiveAlignment(g,tree, *adjIt, seq2);
		heaviestCommonSubsequence(g,seq1,seq2,alignSeq);
		clear(seq1);
		clear(seq2);

		//// Debug Code
		//for(TSize i = 0; i<length(alignSeq);++i) {
		//	std::cout << '(';
		//	for(TSize j=0; j<length(alignSeq[i]);++j) {
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
progressiveAlignment(inputGraph, guideTree, outputGraph [,seqsPerGroup])
..param.inputGraph:The alignment graph with multiple sequence information.
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
		value(profileStrings, i) = alignSeq;
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




















// Testing stuff




//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TScore, typename TOutGraph>
inline TCargo 
iterativeProgressiveAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
							  TGuideTree& tree,
							  TScore const& score_type,
							  TOutGraph& gOut)			 
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjacencyIterator;
	typedef String<TVertexDescriptor> TVertexString;
	typedef String<TVertexString> TSegmentString;

	// Perform initial progressive alignment
	TSegmentString alignSeq;
	_recursiveProgressiveAlignment(g,tree,getRoot(tree),alignSeq);
	clearVertices(gOut);
	_createAlignmentGraph(g, alignSeq, gOut);
	TScoreValue maxSumScore = sumOfPairsScore(gOut, score_type);
	std::cout << maxSumScore << std::endl;

	// Build cutting order of edges
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	String<TVertexDescriptor, Block<> > finalOrder;
	std::deque<TVertexDescriptor> toDoList;
	toDoList.push_back(getRoot(tree));
	while(!toDoList.empty()) {
		TVertexDescriptor v = toDoList[0];
		toDoList.pop_front();
		push_back(finalOrder, v);
		if (!isLeaf(tree, v)) {
			TAdjacencyIterator adjIt(tree, v);
			toDoList.push_back(*adjIt);
			goNext(adjIt);
			toDoList.push_back(*adjIt);
		}
	}

	// Iterative alignment of profiles
	TSize nseq = length(stringSet(g));
	TSize len = length(finalOrder) - 1;  // Don't touch the root again
	TSize iterationsWithoutImprovement = 0;
	for(TSize edge_count=0; ((edge_count<len) && (iterationsWithoutImprovement<5));++edge_count) {
		TVertexDescriptor subtree_root = top(finalOrder);
		pop(finalOrder);

		// Collect all vertex descriptors that belong to the subtree
		std::set<TVertexDescriptor> subtree;
		toDoList.clear();
		toDoList.push_back(subtree_root);
		TSize numSeqs2 = 0;
		while(!toDoList.empty()) {
			TVertexDescriptor v = toDoList[0];
			toDoList.pop_front();
			if (!isLeaf(tree, v)) {
				TAdjacencyIterator adjIt(tree, v);
				toDoList.push_back(*adjIt);
				goNext(adjIt);
				toDoList.push_back(*adjIt);
			} else {
				// Insert the sequence id into the subtree set
				subtree.insert(positionToId(stringSet(g), v));
				++numSeqs2;
			}
		}

		// Build the 2 profile strings
		TSegmentString seq1;
		TSegmentString seq2;
		TSize numSeqs1 = nseq - numSeqs2;
		TSize alignSeqLen = length(alignSeq);
		for(TSize i = 0; i<alignSeqLen;++i) {
			TVertexString set1; TVertexString set2;
			TSize count1 = 0; TSize count2 = 0;
			TVertexString& alignSeq_i = alignSeq[i];
			TSize vertexSetLen = length(alignSeq_i);
			for(TSize j=0; j<vertexSetLen; ++j) {
				fill(set1, numSeqs1, nilVertex);
				fill(set2, numSeqs2, nilVertex);
				TVertexDescriptor v = getValue(alignSeq_i, j);
				if (v ==nilVertex) continue;
				else if (subtree.find(sequenceId(g,v)) == subtree.end()) set1[count1++] = v;
				else set2[count2++] = v;
			}
			if (count1 != 0) appendValue(seq1, set1);
			if (count2 != 0) appendValue(seq2, set2);
		}

		// Align profile strings
		TSegmentString localAlignSeq;
		heaviestCommonSubsequence(g,seq1,seq2,localAlignSeq);
		clearVertices(gOut);
		_createAlignmentGraph(g, localAlignSeq, gOut);
		TScoreValue localSumScore = sumOfPairsScore(gOut, score_type);
		std::cout << localSumScore << std::endl;

		// New maximum?
		if (localSumScore > maxSumScore) {
			iterationsWithoutImprovement = 0;
			maxSumScore = localSumScore;
			alignSeq = localAlignSeq;
		} else {
			++iterationsWithoutImprovement;
		}
	}

	// Create the alignment graph
	clearVertices(gOut);
	_createAlignmentGraph(g, alignSeq, gOut);
	maxSumScore = sumOfPairsScore(gOut, score_type);
	std::cout << maxSumScore << std::endl;

	// Return alignment score
	return (TCargo) maxSumScore;
}







//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TSegmentString, typename TEdgeMap, typename TOutGraph>
inline void 
_createMatchingGraph(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
					 TSegmentString& alignSeq,
					 TEdgeMap& edgeMap,
					 TOutGraph& gOut,
					 TEdgeMap& edgeMapOut)			 
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef String<TVertexDescriptor> TVertexString;

	// Initialization
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

	// Create the matching graph
	clear(edgeMapOut);
	TSize alignSeqLen = length(alignSeq);
	for(TSize i = 0; i<alignSeqLen;++i) {
		TVertexString& alignSeq_i = alignSeq[i];
		TSize len_i = length(alignSeq_i);
		for(TSize j=0; j<len_i - 1; ++j) {
			for(TSize k=j+1; k<len_i; ++k) {
				TVertexDescriptor v1 = getValue(alignSeq_i, j);
				TVertexDescriptor v2 = getValue(alignSeq_i, k);
				if ((v1 == nilVertex) || (v2 == nilVertex)) continue;
				
				TVertexDescriptor v1New = findVertex(gOut, sequenceId(g, v1), fragmentBegin(g,v1));
				if (v1New == nilVertex) v1New = addVertex(gOut, sequenceId(g, v1), fragmentBegin(g,v1), fragmentLength(g,v1));
				TVertexDescriptor v2New = findVertex(gOut, sequenceId(g, v2), fragmentBegin(g,v2));
				if (v2New == nilVertex) v2New = addVertex(gOut, sequenceId(g, v2), fragmentBegin(g,v2), fragmentLength(g,v2));
			
				TEdgeDescriptor e = findEdge(g, v1, v2);
				addEdge(gOut, v1New, v2New, cargo(e));
				appendValue(edgeMapOut, property(edgeMap, e));
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TString, typename TOutString>
inline TCargo
heaviestMatching(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
				 TString const& str1, 
				 TString const& str2,
				 TOutString& align) 
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef __int64 TLargeSize;
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename Value<TString>::Type TVertexSet;

	TSize m = length(str1);  // How many sets of vertex descriptors in seq1
	TSize n = length(str2);  // How many sets of vertex descriptors in seq1

	// Size of the sequences
	// Note for profile alignments every member of the sequence is a String!!! of vertex descriptors
	TCargo divider = (TCargo) length(str1[0]) * (TCargo) length(str2[0]);
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	
	// Fill the vertex to position map for str1
	// Remember for each vertex descriptor the position in the sequence
	typedef std::map<TVertexDescriptor, TSize> TVertexToPosMap;
	typedef typename TVertexToPosMap::const_iterator TVertexToPosMapIter;
	TVertexToPosMap map;
	typedef typename Iterator<TString const>::Type TStringIterConst;
	typedef typename Iterator<TVertexSet const>::Type TVertexSetIterConst;
	TStringIterConst itStrEnd1 = end(str1);
	TSize pos = 0;
	for(TStringIterConst itStr1 = begin(str1);itStr1 != itStrEnd1;++itStr1, ++pos) {
		TVertexSetIterConst itVEnd = end(getValue(itStr1));	
		for(TVertexSetIterConst itV = begin(getValue(itStr1));itV != itVEnd;++itV) {
			if (*itV != nilVertex) map.insert(std::make_pair(*itV, pos));
		}
	}

	// Given a full graph, what positions are occupied?
	typedef std::set<TLargeSize> TOccupiedPositions;
	TOccupiedPositions occupiedPositions;
	TStringIterConst itStrEnd2 = end(str2);
	TSize posItStr2 = 0;
	for(TStringIterConst itStr2 = begin(str2);itStr2 != itStrEnd2;++itStr2, ++posItStr2) {
		TVertexSetIterConst itVEnd = end(getValue(itStr2));
		for(TVertexSetIterConst itV = begin(getValue(itStr2));itV != itVEnd;++itV) {
			if (*itV != nilVertex) {
				TOutEdgeIterator itOut(g, *itV);
				for(;!atEnd(itOut); ++itOut) {
					// Target vertex must be in the map
					TVertexToPosMapIter pPos = map.find(targetVertex(itOut));
					if (pPos != map.end()) occupiedPositions.insert( (TLargeSize) (pPos->second) * (TLargeSize) n + (TLargeSize) (n - posItStr2 - 1) );
				}
			}
		}
	}
	// Map the occupied positions to slots
	typedef std::map<TLargeSize, TSize> TPositionToSlotMap;
	TPositionToSlotMap posToSlotMap;
	TSize counter = 0;
	for(typename TOccupiedPositions::const_iterator setIt = occupiedPositions.begin();setIt != occupiedPositions.end(); ++setIt, ++counter) {
		posToSlotMap.insert(std::make_pair(*setIt, counter));
	}
	occupiedPositions.clear();

	// Walk through str2 and fill in the weights of the actual edges
	typedef String<TCargo> TWeights;
	typedef typename Iterator<TWeights>::Type TWeightsIter;
	TWeights weights;
	fill(weights, posToSlotMap.size(),0);
	posItStr2 = 0;
	for(TStringIterConst itStr2 = begin(str2);itStr2 != itStrEnd2;++itStr2, ++posItStr2) {
		TVertexSetIterConst itVEnd = end(getValue(itStr2));
		for(TVertexSetIterConst itV = begin(getValue(itStr2));itV != itVEnd;++itV) {
			if (*itV != nilVertex) {
				TOutEdgeIterator itOut(g, *itV);
				for(;!atEnd(itOut); ++itOut) {
					// Target vertex must be in the map
					TVertexToPosMapIter pPos = map.find(targetVertex(itOut));
					if (pPos != map.end()) weights[posToSlotMap[ (TLargeSize) (pPos->second) * (TLargeSize) n + (TLargeSize) (n - posItStr2 - 1) ]] += (TCargo) cargo(*itOut);
				}
			}
		}
	}
	map.clear();
	// Average weights
	TWeightsIter itWeights = begin(weights);
	TWeightsIter itWeightsEnd = begin(weights);
	for(;itWeights != itWeightsEnd; ++itWeights) *itWeights /= divider;

	// Create the match graph
	Graph<Undirected<void> > matchGraph;
	typedef typename VertexDescriptor<Graph<Undirected<void> > >::Type TVD;
	String<bool> vertexMap;
	fill(vertexMap, m + n, false);
	for(TSize i = 0; i < m; ++i) addVertex(matchGraph);
	for(TSize i = 0; i < n; ++i) value(vertexMap, addVertex(matchGraph)) = true;
	for(typename TPositionToSlotMap::const_iterator mapIt = posToSlotMap.begin();mapIt != posToSlotMap.end(); ++mapIt) {
		TLargeSize pos = mapIt->first;
		TSize i = (TSize) (pos / (TLargeSize) n);   // Get the index in str1
		TSize j = n - 1 - (TSize) (pos % (TLargeSize) n); // Get the index in str2
		addEdge(matchGraph, i, m+j);
	}

	// Compute the best matching
	typedef String<Pair<TVD, TVD> > TEdges;
	TEdges edges;

	//std::fstream strm;
	//strm.open("Z:\\my_graph.dot", std::ios_base::out | std::ios_base::trunc);
	//write(strm,matchGraph,DotDrawing());
	//strm.close();
	
	TCargo val = weighted_bipartite_matching(matchGraph, vertexMap, weights, edges);
	
	// Retrieve the aligned segments
	TSize seqsInStr1 = length(str1[0]);	 
	TSize seqsInStr2 = length(str2[0]);	
	typedef typename Value<TString>::Type TVertexSet;
	typedef typename Iterator<TString const, Rooted>::Type TStringIter;
	typedef typename Iterator<TString, Rooted>::Type TSIter;
	typedef typename Iterator<TVertexSet const, Rooted>::Type TVertexSetIter;
	typedef typename Iterator<TVertexSet, Rooted>::Type TIter;	
	clear(align);
	// Retrieve all matches
	String<bool> matchedVertices;
	fill(matchedVertices, n + m, false);
	typedef typename Iterator<TEdges>::Type TEdgesIter;
	TEdgesIter itEdges = begin(edges);
	TEdgesIter itEdgesEnd = end(edges);
	for(;itEdges != itEdgesEnd; ++itEdges) {
		TSize i = (value(itEdges)).i1;
		TSize j = (value(itEdges)).i2 - m; 
		value(matchedVertices, i) = true;
		value(matchedVertices, m+j) = true;
		TVertexSet tmp;
		fill(tmp, (seqsInStr1 + seqsInStr2), nilVertex);
		TVertexSetIter itVEnd = end(value(str1,i));
		TSize count = 0;
		for(TVertexSetIter itV = begin(value(str1,i));itV != itVEnd;++itV) {
			tmp[count] = *itV;
			++count;
		}
		TVertexSetIter itVEnd2 = end(value(str2,j));
		for(TVertexSetIter itV2 = begin(value(str2,j));itV2 != itVEnd2;++itV2) {
			tmp[count] = *itV2;
			++count;
		}
		appendValue(align, tmp);
	}
	typedef typename Iterator<Graph<Undirected<void> >, VertexIterator>::Type TVertexIterator;
	TVertexIterator itVertex(matchGraph);
	for(;!atEnd(itVertex);++itVertex) {
		if (getProperty(matchedVertices, *itVertex) == true) continue;
		if (*itVertex < m) {
			TSize i = *itVertex;
			TVertexSet tmp;
			fill(tmp, (seqsInStr1 + seqsInStr2), nilVertex);
			TVertexSetIter itVEnd = end(value(str1, i));
			TSize count = 0;
			for(TVertexSetIter itV = begin(value(str1, i));itV != itVEnd;++itV) {
				tmp[count] = *itV;
				++count;
			}
			appendValue(align, tmp);
		} else {
			TSize j = *itVertex - m;
			TVertexSet tmp;
			fill(tmp, (seqsInStr1 + seqsInStr2), nilVertex);
			TVertexSetIter itVEnd = end(value(str2, j));
			TSize count = 0;
			for(TVertexSetIter itV = begin(value(str2, j));itV != itVEnd;++itV) {
				tmp[count] = *itV;
				++count;
			}
			appendValue(align, tmp);
		}
	}

	return val;
}


template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TVertexDescriptor, typename TSequence>
inline void 
_recursiveProgressiveMatching(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
							  TGuideTree& tree,
							  TVertexDescriptor const root,
							  TSequence& alignSeq)
{
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
		_recursiveProgressiveMatching(g,tree, *adjIt, seq1);
		goNext(adjIt);
		_recursiveProgressiveMatching(g,tree, *adjIt, seq2);

		heaviestMatching(g, seq1, seq2, alignSeq);

		//// Debug Code
		//for(TSize i = 0; i<length(alignSeq);++i) {
		//	std::cout << '(';
		//	for(TSize j=0; j<length(alignSeq[i]);++j) {
		//		std::cout << getValue(alignSeq[i], j) << ',';
		//	}
		//	std::cout << ')' << ',';
		//}
		//std::cout << std::endl;
	}
}


template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TEdgeMap, typename TOutGraph>
inline void 
progressiveMatching(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					TGuideTree& tree,
					TEdgeMap& edgeMap,
					TOutGraph& gOut,
					TEdgeMap& edgeMapOut)			 
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef String<TVertexDescriptor> TVertexString;
	typedef String<TVertexString> TSegmentString;

	// Perform progressive alignment
	TSegmentString alignSeq;
	_recursiveProgressiveMatching(g,tree,getRoot(tree),alignSeq);

	// Create the alignment graph
	_createMatchingGraph(g, alignSeq, edgeMap, gOut, edgeMapOut);
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
