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

	//TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	TStringSet& str = stringSet(g);
	TSize lenRoot = length(str[pos]);
	TId seqId = positionToId(str, pos);
	TSize i = 0;
	while(i<lenRoot) {
		TVertexDescriptor nextVertex = findVertex(const_cast<TGraph&>(g), seqId, i);
		SEQAN_TASSERT(nextVertex != getNil<TVertexDescriptor>())
		//if (nextVertex == nilVertex) {
		//	std::cout << "Nil Vertex" << std::endl;
		//	TSize j = i + 1;
		//	while ((j < lenRoot) && (findVertex(const_cast<TGraph&>(g), seqId, j) == nilVertex)) ++j;
		//	nextVertex = addVertex(const_cast<TGraph&>(g), seqId, i, j-i);
		//}
		appendValue(alignSeq, String<TVertexDescriptor>(nextVertex));
		i += fragmentLength(g, nextVertex);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TVertexDescriptor, typename TSequence>
inline TCargo 
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

	TCargo score = 0;
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
		score = heaviestCommonSubsequence(g,seq1,seq2,alignSeq);

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
	return score;
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

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TOutGraph>
inline TCargo 
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

	// Perform initial progressive alignment
	TSegmentString alignSeq;
	TCargo maxSumScore = _recursiveProgressiveAlignment(g,tree,getRoot(tree),alignSeq);

	// Create the alignment graph
	_createAlignmentGraph(g, alignSeq, gOut);

	// Return alignment score
	return (TCargo) maxSumScore;
}


///////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TVertexDescriptor, typename TSet>
inline void
_getAllChildren(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				TGuideTree& tree,
				TVertexDescriptor const root,
				TSet& set1)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TGuideTree>::Type TSize;
	typedef typename Id<TGuideTree>::Type TId;
	typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjacencyIterator;

	if (isLeaf(tree, root)) {
		set1.insert(positionToId(stringSet(g), root));
	} else {
		TAdjacencyIterator adjIt(tree, root);
		for(;!atEnd(adjIt);goNext(adjIt)) {
			_getAllChildren(g, tree, *adjIt, set1);
		}
	}
}
				

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TVertexDescriptor, typename TSet, typename TSequence>
inline TCargo 
_recursiveProgressiveAlignmentWithTriplet(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
							   TGuideTree& tree,
							   TVertexDescriptor const root,
							   TSet& id_set,
							   bool tripletDone,
							   TSequence& alignSeq)
{
	SEQAN_CHECKPOINT

	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjacencyIterator;
	typedef typename Iterator<TGuideTree, DfsPreorder>::Type TDfsPreorderIterator;

	TCargo score = 0;
	if(isLeaf(tree, root)) {
		_buildLeafString(g, root, alignSeq);
		id_set.insert(positionToId(stringSet(g), root));
	} else {
		if (tripletDone) {
			// Align the two children (Binary tree)
			typedef String<String<TVertexDescriptor> > TSegmentString;
			TSegmentString seq1;
			TSegmentString seq2;
			TAdjacencyIterator adjIt(tree, root);
			std::set<TId> set1;
			_recursiveProgressiveAlignmentWithTriplet(g,tree, *adjIt, set1, tripletDone, seq1);
			goNext(adjIt);
			std::set<TId> set2;
			_recursiveProgressiveAlignmentWithTriplet(g,tree, *adjIt, set2, tripletDone, seq2);
			id_set.insert(set1.begin(),set1.end());
			id_set.insert(set2.begin(),set2.end());
			score = heaviestCommonSubsequence(g,seq1,seq2,alignSeq);
		} else {
			std::set<TId> children;
			_getAllChildren(g, tree, root, children);
			if (children.size() < 5) {
				TGraph copy_graph(g);
				tripletLibraryExtensionOnSet(copy_graph, children);
				tripletDone = true;
				_recursiveProgressiveAlignmentWithTriplet(copy_graph,tree, root, id_set, tripletDone, alignSeq);
				std::cout << "One mini tree finished " << id_set.size() << std::endl;
				clear(copy_graph);
			} else {
				// Align the two children (Binary tree)
				typedef String<String<TVertexDescriptor> > TSegmentString;
				TSegmentString seq1;
				TSegmentString seq2;
				TAdjacencyIterator adjIt(tree, root);
				std::set<TId> set1;
				_recursiveProgressiveAlignmentWithTriplet(g,tree, *adjIt, set1, tripletDone, seq1);
				goNext(adjIt);
				std::set<TId> set2;
				_recursiveProgressiveAlignmentWithTriplet(g,tree, *adjIt, set2, tripletDone, seq2);
				
				TGraph copy_graph(g);
				restrictedTripletLibraryExtension(copy_graph, set1, set2);
				score = heaviestCommonSubsequence(g,seq1,seq2,alignSeq);
				id_set.insert(set1.begin(),set1.end());
				id_set.insert(set2.begin(),set2.end());
				std::cout << "Two mini trees joined " << id_set.size() << std::endl;
			}
		} 
		
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
	return score;
}



//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TOutGraph, typename TSize>
inline void
progressiveAlignmentWithTriplet(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
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


	//tripletLibraryExtension(g);
	//progressiveAlignment(g, tree, gOut);
	//return;

	// Build groups of sequences
	TGuideTree copy_tree(tree);
	String<std::set<TSize> > sequenceGroups;
	String<TId> sequenceGroupRoots;
	while (numChildren(copy_tree, getRoot(copy_tree)) > 0) {
		TBfsIterator bfsIt(copy_tree, getRoot(copy_tree));
		for(;!atEnd(bfsIt);goNext(bfsIt)) {
			std::set<TSize> sequenceSet;
			_getAllChildren(g, copy_tree, *bfsIt, sequenceSet);
			if ((TSize) sequenceSet.size() <= (TSize) seqPerGroup) {
				appendValue(sequenceGroups, sequenceSet);
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
		if ((sequenceGroups[i]).size() > 1) tripletLibraryExtensionOnSet(copy_graph, sequenceGroups[i]);
		TSegmentString alignSeq;
		_recursiveProgressiveAlignment(copy_graph,tree,sequenceGroupRoots[i],alignSeq);
		profileStrings[i] = alignSeq;
		std::cout << "One mini tree finished " << (sequenceGroups[i]).size() << std::endl;
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
		std::cout << "Two mini trees joined " << (sequenceGroups[lastIndex]).size() << std::endl;
		//clear(copy_graph);
	}

	// Create the alignment graph
	_createAlignmentGraph(g, profileStrings[lastIndex], gOut);
}

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
/*
template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TScore, typename TOutGraph>
inline void 
iterativeProgressiveAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
							  TGuideTree& tree,
							  TScore const& score_type,
							  TOutGraph& gOut)			 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Infix<TString>::Type TInfix;
	typedef typename Iterator<TGuideTree, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjacencyIterator;
	typedef typename Iterator<TGuideTree, DfsPreorder>::Type TDfsPreorderIterator;

	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	TSize nseq = length(stringSet(g));


	// Perform initial progressive alignment
	typedef String<TVertexDescriptor> TVertexString;
	typedef String<TVertexString> TSegmentString;
	TSegmentString alignSeq;
	_recursiveProgressiveAlignment(g,tree,getRoot(tree),alignSeq);
	TScoreValue maxSumScore = sumOfPairsScore(g, alignSeq, score_type);

	// Build cutting order of edges
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
	TSize len = length(finalOrder) - 1;  // Don't touch the root again
	for(TSize edge_count=0; edge_count<len;++edge_count) {
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
		for(TSize i = 0; i<length(alignSeq);++i) {
			TVertexString set1;
			TVertexString set2;
			TSize count1 = 0;
			TSize count2 = 0;
			for(TSize j=0; j<length(alignSeq[i]);++j) {
				fill(set1, numSeqs1, nilVertex);
				fill(set2, numSeqs2, nilVertex);
				TVertexDescriptor v = (alignSeq[i])[j];
				if (v ==nilVertex) continue;
				else if (subtree.find(sequenceId(g,v)) == subtree.end()) {
					set1[count1] = v;
					++count1;
				} else {
					set2[count2] = v;
					++count2;
				}
			}
			if (count1 != 0) {
				appendValue(seq1, set1);
			}
			if (count2 != 0) {
				appendValue(seq2, set2);
			}
		}

		// Align profile strings
		TSegmentString localAlignSeq;
		_alignStringsAccordingToGraph(g,seq1,seq2,localAlignSeq);
		TScoreValue localSumScore = sumOfPairsScore(g, localAlignSeq, score_type);

		if (localSumScore > maxSumScore) {
			maxSumScore = localSumScore;
			alignSeq = localAlignSeq;
		}
	}

	if (nseq > 10) {
		// Now try a random 3-cut
		mtRandInit();
		TSize repeatsWithoutImprove = 0;
		while (repeatsWithoutImprove < 5) {
			// Build the vertexPool
			std::set<TVertexDescriptor> vertexPool;
			TVertexIterator itV(tree);
			for(;!atEnd(itV);++itV) vertexPool.insert(*itV);
			vertexPool.erase(getRoot(tree));

			// Pick a random vertex
			typename std::set<TVertexDescriptor>::const_iterator pos = vertexPool.begin();
			TSize limit = ((Byte) mtRand() % (vertexPool.size()));
			for(TSize i=0; i < limit; ++i) ++pos;
			TSize subtree_root1 = *pos;
				
			// Collect all vertex descriptors that belong to the subtree
			std::set<TVertexDescriptor> subtree1;
			std::deque<TVertexDescriptor> toDoList;
			toDoList.push_back(subtree_root1);
			TSize numSeqs1 = 0;
			while(!toDoList.empty()) {
				TVertexDescriptor v = toDoList[0];
				toDoList.pop_front();
				vertexPool.erase(v);
				if (!isLeaf(tree, v)) {
					TAdjacencyIterator adjIt(tree, v);
					toDoList.push_back(*adjIt);
					goNext(adjIt);
					toDoList.push_back(*adjIt);
				} else {
					// Insert the sequence id into the subtree set
					subtree1.insert(positionToId(stringSet(g), v));
					++numSeqs1;
				}
			}
	
			// Pick a second random vertex
			pos = vertexPool.begin();
			limit = ((Byte) mtRand() % (vertexPool.size()));
			for(TSize i=0; i < limit; ++i) ++pos;
			TSize subtree_root2 = *pos;

			// Collect all vertex descriptors that belong to the subtree
			std::set<TVertexDescriptor> subtree2;
			toDoList.clear();
			toDoList.push_back(subtree_root2);
			TSize numSeqs2 = 0;
			while(!toDoList.empty()) {
				TVertexDescriptor v = toDoList[0];
				toDoList.pop_front();
				if (v == subtree_root1) {
					// A child is the root of the other tree -> Break!
					vertexPool.clear();
					break;
				}
				vertexPool.erase(v);
				if (!isLeaf(tree, v)) {
					TAdjacencyIterator adjIt(tree, v);
					toDoList.push_back(*adjIt);
					goNext(adjIt);
					toDoList.push_back(*adjIt);
				} else {
					// Insert the sequence id into the subtree set
					subtree2.insert(positionToId(stringSet(g), v));
					++numSeqs2;
				}
			}

			if (vertexPool.empty()) continue;

			// Build the 3 profile strings
			TSegmentString seq0;
			TSegmentString seq1;
			TSegmentString seq2;
			TSize numSeqs0 = nseq - numSeqs1 - numSeqs2;
			if (numSeqs0 == 0) continue;
	
			for(TSize i = 0; i<length(alignSeq);++i) {
				TVertexString set0;
				TVertexString set1;
				TVertexString set2;
				TSize count0 = 0;
				TSize count1 = 0;
				TSize count2 = 0;
				for(TSize j=0; j<length(alignSeq[i]);++j) {
					fill(set0, numSeqs0, nilVertex);
					fill(set1, numSeqs1, nilVertex);
					fill(set2, numSeqs2, nilVertex);
					TVertexDescriptor v = (alignSeq[i])[j];
					if (v ==nilVertex) continue;
					else if (subtree1.find(sequenceId(g,v)) != subtree1.end()) {
						set1[count1] = v;
						++count1;
					} else if (subtree2.find(sequenceId(g,v)) != subtree2.end()) {
						set2[count2] = v;
						++count2;
					} else {
						set0[count0] = v;
						++count0;
					}
				}
				if (count0 != 0) {
					appendValue(seq0, set0);
				}
				if (count1 != 0) {
					appendValue(seq1, set1);
				}
				if (count2 != 0) {
					appendValue(seq2, set2);
				}
			}

			// Align all possible combinations
			TSegmentString tmp;
			TSegmentString localAlignSeq0;
			_alignStringsAccordingToGraph(g,seq0,seq1,tmp);
			_alignStringsAccordingToGraph(g,tmp,seq2,localAlignSeq0);
			TScoreValue score0 = sumOfPairsScore(g, localAlignSeq0, score_type);
			clear(tmp);
			TSegmentString localAlignSeq1;
			_alignStringsAccordingToGraph(g,seq0,seq2,tmp);
			_alignStringsAccordingToGraph(g,tmp,seq1,localAlignSeq1);
			TScoreValue score1 = sumOfPairsScore(g, localAlignSeq0, score_type);
			clear(tmp);
			TSegmentString localAlignSeq2;
			_alignStringsAccordingToGraph(g,seq1,seq2,tmp);
			_alignStringsAccordingToGraph(g,tmp,seq0,localAlignSeq2);
			TScoreValue score2 = sumOfPairsScore(g, localAlignSeq0, score_type);

		
			if ((score0 > maxSumScore) &&
				(score0 > score1) &&
				(score0 > score2)) {
				maxSumScore = score0;
				alignSeq = localAlignSeq0;
				repeatsWithoutImprove = 0;
			} else if ((score1 > maxSumScore) &&
					(score1 > score2)) {
				maxSumScore = score1;
				alignSeq = localAlignSeq1;
				repeatsWithoutImprove = 0;
			} else if (score2 > maxSumScore) {
				maxSumScore = score2;
				alignSeq = localAlignSeq2;
				repeatsWithoutImprove = 0;
			} else {
				++repeatsWithoutImprove;
			}
		}
	}

	// Create the alignment graph
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
*/


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
