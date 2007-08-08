#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_PROGRESSIVE_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_PROGRESSIVE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// T-Coffee - Progressive Alignment
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TSequenceIn, typename TSequence>
inline TCargo 
_alignStringsAccordingToGraph(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
							  TSequenceIn const& seq1,
							  TSequenceIn const& seq2,
							  TSequence& alignSeq)
{
	SEQAN_CHECKPOINT
	// Clear output parameter
	clear(alignSeq);

	// Just heaviest common subsequence
	return heaviestCommonSubsequence(g, seq1, seq2, alignSeq);
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

	TStringSet& str = stringSet(g);
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	TCargo score = 0;

	if(isLeaf(tree, root)) {
		TId seqId = positionToId(str, root);
		TSize i = 0;
		TSize lenRoot = length(str[root]);
		while(i<lenRoot) {
			TVertexDescriptor nextVertex = findVertex(g, seqId, i);
			SEQAN_TASSERT(nextVertex != getNil<TVertexDescriptor>())
			if (nextVertex == nilVertex) {
				std::cout << "Nil Vertex" << std::endl;
				TSize j = i + 1;
				while ((j < lenRoot) && (findVertex(g, seqId, j) == nilVertex)) ++j;
				nextVertex = addVertex(g, seqId, i, j-i);
			}
			appendValue(alignSeq, String<TVertexDescriptor>(nextVertex));
			i += fragmentLength(g, nextVertex);
		}
	} else {
		// Align the two children (Binary tree)
		typedef String<String<TVertexDescriptor> > TSegmentString;
		TSegmentString seq1;
		TSegmentString seq2;
		TAdjacencyIterator adjIt(tree, root);
		_recursiveProgressiveAlignment(g,tree, *adjIt, seq1);
		goNext(adjIt);
		_recursiveProgressiveAlignment(g,tree, *adjIt, seq2);
		score = _alignStringsAccordingToGraph(g,seq1,seq2,alignSeq);

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
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Infix<TString>::Type TInfix;
	typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjacencyIterator;
	typedef typename Iterator<TGuideTree, DfsPreorder>::Type TDfsPreorderIterator;

	// Initialization
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	TSize nseq = length(stringSet(g));

	// Perform initial progressive alignment
	typedef String<TVertexDescriptor> TVertexString;
	typedef String<TVertexString> TSegmentString;
	TSegmentString alignSeq;
	//TCargo maxSumScore = _recursiveProgressiveAlignment(g,tree,getRoot(tree),alignSeq);
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
		_alignStringsAccordingToGraph(g,seq1,seq2,localAlignSeq);
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
