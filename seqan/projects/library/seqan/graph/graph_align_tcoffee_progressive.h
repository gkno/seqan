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
// T-Coffee - Progressive Alignment
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
					addEdge(gOut, l - count, l);
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


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TVertexDescriptor, typename TSequence, typename TTag>
inline void 
_recursiveProgressiveAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
							   TGuideTree& tree,
							   TVertexDescriptor const root,
							   TSequence& alignSeq,
							   Tag<TTag>)
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
		typedef StringSet<TSegmentString, Owner<> > TSegmentStringSet;
		TSegmentStringSet strSet;
		resize(strSet,2);
		TAdjacencyIterator adjIt(tree, root);
		_recursiveProgressiveAlignment(g,tree, *adjIt, value(strSet, 0), Tag<TTag>());
		goNext(adjIt);
		_recursiveProgressiveAlignment(g,tree, *adjIt, value(strSet, 1), Tag<TTag>());
		Score<int, ScoreAlignmentGraph<TGraph> > score_type = Score<int, ScoreAlignmentGraph<TGraph> >(g);
		TSequence tmp;
		globalAlignment(tmp, strSet, score_type, Tag<TTag>());
		TSize len = length(tmp);
		resize(alignSeq, len);
		for(TSize i = len; i>0;--i) alignSeq[len-i] = tmp[i-1];

		//// Debug Code
		//for(unsigned int i = 0; i<length(alignSeq);++i) {
		//	int fragLen = 0;
		//	std::cout << '(';
		//	for(unsigned int j=0; j<length(alignSeq[i]);++j) {
		//		std::cout << getValue(alignSeq[i], j) << ',';
		//		if ((fragLen>0) &&
		//			(getNil<typename VertexDescriptor<TGraph>::Type>() != getValue(alignSeq[i], j))) {
		//			SEQAN_TASSERT(fragmentLength(g, getValue(alignSeq[i], j)) == fragLen)
		//		} else {
		//			if (getNil<typename VertexDescriptor<TGraph>::Type>() != getValue(alignSeq[i], j)) {
		//				fragLen = fragmentLength(g, getValue(alignSeq[i], j));
		//			}
		//		}
		//	}
		//	std::cout << ')' << ',';
		//}
		//std::cout << std::endl;
	}
}

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TOutGraph, typename TTag>
inline void 
progressiveAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					 TGuideTree& tree,
					 TOutGraph& gOut,
					 Tag<TTag>)			 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef String<TVertexDescriptor> TVertexString;
	typedef String<TVertexString> TSegmentString;

	// Perform initial progressive alignment
	TSegmentString alignSeq;
	_recursiveProgressiveAlignment(g,tree,getRoot(tree),alignSeq, Tag<TTag>());

	// Create the alignment graph
	_createAlignmentGraph(g, alignSeq, gOut);
}



/*
//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TScore, typename TOutGraph>
inline TCargo 
iterativeProgressiveAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
							  TGuideTree& tree,
							  TScore const& score_type,
							  TOutGraph& gOut)			 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjacencyIterator;
	typedef String<TVertexDescriptor> TVertexString;
	typedef String<TVertexString> TSegmentString;

	// Perform initial progressive alignment
	TSegmentString alignSeq;
	//TCargo maxSumScore = _recursiveProgressiveAlignment(g,tree,getRoot(tree),alignSeq);
	_recursiveProgressiveAlignment(g,tree,getRoot(tree),alignSeq);
	TScoreValue maxSumScore = sumOfPairsScore(g, alignSeq, score_type);

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
		//TCargo localSumScore = _alignStringsAccordingToGraph(g,seq1,seq2,localAlignSeq);
		heaviestCommonSubsequence(g,seq1,seq2,localAlignSeq);
		TScoreValue localSumScore = sumOfPairsScore(g, localAlignSeq, score_type);

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
	_createAlignmentGraph(g, alignSeq, gOut);

	// Return alignment score
	return (TCargo) maxSumScore;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TOutGraph>
inline void 
highestScoreFirstAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
						   TOutGraph& gOut)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;

	// Initialization
	String<double> distanceMatrix;
	TStringSet& str = stringSet(g);	
	TSize nseq = length(str);
	resize(distanceMatrix, nseq * nseq);

	// Build initial distance matrix
	typedef String<String<TVertexDescriptor> > TSegmentString;
	String<TSegmentString> segmentStr;
	String<bool> active;
	resize(segmentStr, nseq);
	fill(active, nseq, true);
	for(TSize i=0; i<nseq; ++i) _buildLeafString(g, i, value(segmentStr,i));
	for(TSize i=0; i<nseq; ++i) {
		TSize len1 = length(value(segmentStr,i));
		for(TSize j=i+1; j<nseq; ++j) {
			// Align the 2 strings
			TSize len2 = length(value(segmentStr,j));
			double score = hcsPairwiseScore(g,value(segmentStr,i),value(segmentStr,j));
			
			// Normalize by distance
			score /= ((len1 + len2) / 2);
						
			// Remember the value
			assignValue(distanceMatrix, i*nseq+j, score);
		}
	}

	bool oneLeft = true;
	do {
		oneLeft = true;

		// Find highest scoring pair
		TSize index_i = 0;
		TSize index_j = 0;
		double maxScore = 0;
		for(TSize i=0; i<nseq; ++i) {
			if (!active[i]) continue;
			for(TSize j=i+1; j<nseq; ++j) {
				if (!active[j]) continue;
				oneLeft = false;
				//std::cout << getValue(distanceMatrix, i*nseq+j) << ',';
				double tmp;
				if ((tmp = getValue(distanceMatrix, i*nseq+j)) > maxScore) {
					maxScore = tmp;
					index_i = i;
					index_j = j;
				}
			}
			//std::cout << std::endl;
		}
		if (oneLeft) break;
		//std::cout << index_i << ',' << index_j << ':' << maxScore << std::endl;

		// Align this sequence pair
		TSegmentString alignSeq;
		heaviestCommonSubsequence(g,segmentStr[index_i],segmentStr[index_j],alignSeq);
		active[index_j] = false;
		clear(value(segmentStr, index_j));
		segmentStr[index_i] = alignSeq;

		// Recalculate distances
		for(TSize i=0; i<nseq; ++i) {
			if ((!active[i]) || (index_i == i)) continue;
			TSize len1 = length(segmentStr[index_i]);
			TSize len2 = length(segmentStr[i]);
			double score = hcsPairwiseScore(g,segmentStr[index_i],segmentStr[i]);
			
			// Normalize by distance
			score /= ((len1 + len2) / 2);
					
			// Remember the value
			if (index_i < i) assignValue(distanceMatrix, index_i*nseq+i, score);
			else assignValue(distanceMatrix, i*nseq+index_i, score);
		}
	} while (!oneLeft);

	// Create alignment graph
	for(TSize i=0; i<nseq; ++i) {
		if (!active[i]) continue;
		_createAlignmentGraph(g, segmentStr[i], gOut);
		break;
	}
}
*/



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
