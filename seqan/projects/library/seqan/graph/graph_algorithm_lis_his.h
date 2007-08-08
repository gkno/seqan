#ifndef SEQAN_HEADER_GRAPH_ALGORITHM_LIS_HIS_H
#define SEQAN_HEADER_GRAPH_ALGORITHM_LIS_HIS_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// LIS: Longest Increasing Subsequence
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSortedSequence, typename TKey>
inline typename TSortedSequence::const_iterator
_previousInSortedSequence(TSortedSequence const& list, TKey const key) {
	SEQAN_CHECKPOINT
	typedef typename TSortedSequence::const_iterator TSortedSequenceIter;

	TSortedSequenceIter a_k_it = list.lower_bound(key);
	// Now we need to move one to the front

	if (a_k_it != list.end()) {
		// If we are at the beginning, no predecessor
		if (a_k_it == list.begin()) a_k_it = list.end();
		else --a_k_it;
	} else {
		// If we are at the end, the predecessor is the last element of the list
		TSortedSequenceIter tmp = list.begin();
		if (tmp != list.end()) {
			do {
				a_k_it = tmp;
			} while(++tmp != list.end());
		}
	}
	return a_k_it;
}


//////////////////////////////////////////////////////////////////////////////


template<typename TSortedSequence, typename TIterator>
inline typename TSortedSequence::const_iterator
_nextInSortedSequence(TSortedSequence const& list, TIterator const& prev) {
	SEQAN_CHECKPOINT
	typedef typename TSortedSequence::const_iterator TSortedSequenceIter;
		
	TSortedSequenceIter b_l_it;
	if (prev == list.end()) b_l_it = list.begin();
	else b_l_it = list.upper_bound(*prev);

	return b_l_it;
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.longestIncreasingSubsequence:
..summary:Computes the longest increasing subsequence.
..cat:Alignments
..signature:longestIncreasingSubsequence(str, pos)
..param.str:In-parameter: An arbitrary string.
...type:Class.String
..param.pos:Out-parameter: A String with the positions that belong to the longest increasing subsequence.
...remarks:
The last position in pos indicates the first element in the longest increasing subsequence.
That's why pos should be a Block-String (Stack).
*/
template<typename TString, typename TPositions>
inline void
longestIncreasingSubsequence(TString const& str, TPositions& pos) {
	SEQAN_CHECKPOINT

	// The list of decreasing covers, only the smallest number must be remembered
	// See Gusfield
	typedef std::pair<typename Value<TString>::Type, typename Value<TPositions>::Type> TKey;
	typedef std::set<TKey, std::less<TKey> > TSortedSequence;
	typedef typename TSortedSequence::const_iterator TSortedSequenceIter;
	TSortedSequence list;

	// The trace-back graph
	typedef Graph<Directed<void, WithoutEdgeId> > TGraph;
	typedef Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TGraph g;

	// Walk through the sequence and build the decreasing covers
	typedef typename Iterator<TString const>::Type TStringIter;
	TStringIter endIt = end(str);
	for(TStringIter it = begin(str); it != endIt; ++it) {
		// Get previous element
		TSortedSequenceIter a_k_it = _previousInSortedSequence(list, std::make_pair(*it, 0)); 
		
		// Get next element
		TSortedSequenceIter b_l_it = _nextInSortedSequence(list, a_k_it);
		
		// Delete from list
		if (b_l_it != list.end()) list.erase(*b_l_it);

		// Insert new list element
		list.insert(std::make_pair(*it, position(it)));

		// Create the corresponding node
		// Note: The VertexDescriptor == position(it)
		addVertex(g);

		// Connect to predecessor
		if (a_k_it != list.end()) addEdge(g, position(it), a_k_it->second);
	}

	// Trace-back
	if (list.rbegin() == list.rend()) return;
	else {
		bool finished = false;
		// Start with the maximal position in the list == Vertex Descriptor
		TVertexDescriptor v = list.rbegin()->second;
		while (!finished) {
			push_back(pos, v);
			TOutEdgeIterator it(g, v);
			if (atEnd(it)) finished = true;
			else v = targetVertex(it);
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
// LCS: Longest Common Subsequence
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Function.longestCommonSubsequence:
..summary:Computes the longest common subsequence.
..cat:Alignments
..signature:longestCommonSubsequence(str1, str2, pos)
..param.str1:In-parameter: An arbitrary string.
...type:Class.String
..param.str2:In-parameter: An arbitrary string.
...type:Class.String
..param.pos:Out-parameter: A String with pairs of positions that indicate the longest common subsequence.
...remarks:
The last pair of positions in pos indicates the first pair in the longest common subsequence.
That's why pos should be a Block-String (Stack).
*/
template<typename TString1, typename TString2, typename TFinalPos>
inline void
longestCommonSubsequence(TString1 const& str1,
						 TString2 const& str2,
						 TFinalPos& pos) 
{
	SEQAN_CHECKPOINT
	typedef typename Value<TString1>::Type TValue;
	typedef typename Size<TString1>::Type TSize;
	typedef typename Position<TString1>::Type TPos;
	TSize alphabet_size = ValueSize<TValue>::VALUE;

	// The occurences of each letter in the second string
	typedef String<TPos, Block<> > TPositions;
	String<TPositions> occ;
	fill(occ, alphabet_size, TPositions());
	typedef typename Iterator<TString2 const>::Type TStringIter;
	TStringIter endIt = end(str2);
	for(TStringIter it = begin(str2); it != endIt; ++it) {
		push_back(value(occ, (unsigned int) *it), position(it));
	}

	// Build the combined string
	String<TPos, Block<> > finalSeq;
	String<TPos, Block<> > mapping;
	TStringIter endIt1 = end(str1);
	for(TStringIter it = begin(str1); it != endIt1; ++it) {
		for(int i = length(occ[(unsigned int) *it])-1; i>=0; --i) {
			TPos source = position(it);
			push_back(finalSeq, (occ[(unsigned int) *it])[i]);
			push_back(mapping, source);
		}
	}

	// Call longest increasing subsequence
	typedef String<unsigned int, Block<> > TResult;
	TResult result;
	longestIncreasingSubsequence(finalSeq, result);

	// Insert the common pairs
	typedef typename Iterator<TResult>::Type TResultIter;
	TResultIter endResult = end(result);
	for(TResultIter it = begin(result); it != endResult; ++it) {
		push_back(pos, std::make_pair(mapping[*it], finalSeq[*it]));
	}

	//// Debug code
	//for(int i=0; i<length(pos);++i) {
	//	std::cout << pos[i].first << ',' << pos[i].second << ';';
	//}
	//std::cout << std::endl;
}


//////////////////////////////////////////////////////////////////////////////
// HIS: Heaviest Increasing Subsequence
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Function.heaviestIncreasingSubsequence:
..summary:Computes the heaviest increasing subsequence.
..cat:Alignments
..signature:heaviestIncreasingSubsequence(str, weights, pos)
..param.str:In-parameter: An arbitrary string.
...type:Class.String
..param.weights:In-parameter: A weight for each position in the string.
..param.pos:Out-parameter: A String of positions that indicate the members of the heaviest increasing subsequence.
...remarks:
The last position in pos indicates the first member of the heaviest increasing subsequence.
That's why pos should be a Block-String (Stack).
Note that only members that contribute a weight are selected, that is, positions with associated weight=0 are ignored.
*/
template<typename TString, typename TWeightMap, typename TPositions>
inline typename Value<TWeightMap>::Type
heaviestIncreasingSubsequence(TString const& str, 
							  TWeightMap const& weights, 
							  TPositions& pos) 
{
	SEQAN_CHECKPOINT
	typedef typename Value<TString>::Type TValue;
	typedef typename Value<TPositions>::Type TPos;
	typedef typename Value<TWeightMap>::Type TWeight;

	// The list of decreasing covers, only the smallest element of each member must be remembered
	typedef std::pair<TValue, std::pair<TWeight, TPos> > TKey;
	typedef std::set<TKey, std::less<TKey> > TSortedSequence;
	typedef typename TSortedSequence::const_iterator TSortedSequenceIter;
	TSortedSequence list;
	
	// The trace-back graph
	typedef Graph<Directed<void, WithoutEdgeId> > TGraph;
	typedef Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TGraph g;

	// Walk through the sequence and build the decreasing covers
	typedef typename Iterator<TString const>::Type TStringIter;
	TStringIter endIt = end(str);
	for(TStringIter it = begin(str); it != endIt; ++it) {
		TWeight w = getProperty(weights, position(it));
		// Letters that do not contribute a weight (e.g., w = 0) are excluded!
		// Weights must increase!
		if (w <= 0) {
			addVertex(g);
		}


		// Get previous element
		TSortedSequenceIter a_k_it = _previousInSortedSequence(list, std::make_pair(*it, std::make_pair(0, 0))); 
		
		// Get next element
		TSortedSequenceIter b_l_it = _nextInSortedSequence(list, a_k_it);

		// Determine new weight
		if (a_k_it != list.end()) w += a_k_it->second.first;

		// Delete from list
		while ((b_l_it != list.end()) && 
				(w >= b_l_it->second.first)) {
					TSortedSequenceIter tmp = b_l_it;
					b_l_it = _nextInSortedSequence(list, b_l_it);
					list.erase(*tmp);
		}

		// Insert new list element
		if ((b_l_it == list.end()) ||
			(*it < b_l_it->first)) {
				list.insert(std::make_pair(*it, std::make_pair(w, position(it))));
		}

		// Create the corresponding node, position(it) == Vertex Descriptor
		addVertex(g);

		// Connect to predecessor
		if (a_k_it != list.end()) addEdge(g, position(it), a_k_it->second.second);
	}

	// Trace-back
	TWeight w = 0;
	if (list.rbegin() == list.rend()) return 0;
	else {
		bool finished = false;
		// Last vertex is end of heaviest increasing subsequence
		TVertexDescriptor v = list.rbegin()->second.second;
		while (!finished) {
			// Exclude edges with weight 0 !!!
			// Note: Very important, do not delete this check!!!
			if (getProperty(weights, v) > 0) {
				appendValue(pos, v);
				w+=getProperty(weights, v);
			}
			TOutEdgeIterator it(g, v);
			if (atEnd(it)) finished = true;
			else v = targetVertex(it);
		}
	}
	return w;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TString>
inline TCargo
heaviestCommonSubsequence(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
						  TString const& str1, 
						  TString const& str2,
						  TString& align) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TString const>::Type TStringIter;
	typedef typename Iterator<TString>::Type TSIter;
	typedef typename Value<TString>::Type TVertexSet;
	typedef typename Iterator<TVertexSet const>::Type TVertexSetIter;
	typedef typename Iterator<TVertexSet>::Type TIter;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	
	// Size of the sequences
	// Note for profile alignments every member of the sequence is a String!!! of vertex descriptors
	TSize m = length(str1);  // How many sets of vertex descriptors in seq1
	TSize n = length(str2);  // How many sets of vertex descriptors in seq1
	TSize seqsInStr1 = length(str1[0]);	 // #Vertex descriptors per node
	TSize seqsInStr2 = length(str2[0]);	
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	
	// Fill the vertex to position map for str1
	// Remember for each vertex descriptor the position in the sequence
	typedef std::map<TVertexDescriptor, TSize> TVertexToPosMap;
	typedef typename TVertexToPosMap::const_iterator TVertexToPosMapIter;
	TVertexToPosMap map;
	TStringIter itStrEnd1 = end(str1);
	for(TStringIter itStr1 = begin(str1);itStr1 != itStrEnd1;++itStr1) {
		TVertexSetIter itVEnd = end(getValue(itStr1));
		for(TVertexSetIter itV = begin(getValue(itStr1));itV != itVEnd;++itV) {
			if (*itV != nilVertex) map.insert(std::make_pair(*itV, position(itStr1)));
		}
	}

	// We create the full graph
	// For a given node number the edges in decreasing order, so
	// during increasing subsequence computation no two edges are selected for a given node
	// We do this for every node 0 ... m-1
	String<TSize> seq;
	resize(seq, n*m);
	typedef typename Iterator<String<TSize> >::Type TSeqIter;
	TSeqIter itSeq = begin(seq);
	for(TSize i=0; i<m;++i) {
		for(int j=n-1;j>=0;--j) {
			*itSeq = (TSize) j;
			++itSeq;
		}
	}
	// Initially every edge receives weight=0
	String<double> weights;
	// For profile alignments, take the average weight
	double divider = (double) seqsInStr1 * (double) seqsInStr2;
	fill(weights, n*m, 0);

	// Walk through str2 and fill in the weights of the actual edges
	TStringIter itStrEnd2 = end(str2);
	for(TStringIter itStr2 = begin(str2);itStr2 != itStrEnd2;++itStr2) {
		TVertexSetIter itVEnd = end(getValue(itStr2));
		for(TVertexSetIter itV = begin(getValue(itStr2));itV != itVEnd;++itV) {
			if (*itV != nilVertex) {
				TOutEdgeIterator itOut(g, *itV);
				for(;!atEnd(itOut); ++itOut) {
					// Target vertex must be in the map
					TVertexToPosMapIter pPos = map.find(targetVertex(itOut));
					if (pPos != map.end()) {
						// Calculate the edge index
						TSize index = pPos->second * n + (n - position(itStr2) - 1);
						weights[index] += (double) cargo(*itOut) / divider;	
					}
				}
			}
		}
	}
	map.clear();

	// Calculate the heaviest increasing subsequence
	// Note edges with weight=0 are ignored!
	String<unsigned int> pos;
	TCargo score = (TCargo) heaviestIncreasingSubsequence(seq, weights, pos);
	
	// Create the alignment sequence
	TSize numMatches = length(pos);
	TSize alignLength = numMatches + (n - numMatches) + (m - numMatches);
	clear(align);
	resize(align, alignLength);
	TSIter pointerAlign = begin(align);
	TSIter pointerAlignEnd = end(align);
	TStringIter pointerStr1 = begin(str1);
	TStringIter pointerStr2 = begin(str2);
	int p = length(pos)-1;
	while(pointerAlign != pointerAlignEnd) {
		TSize i = m;
		TSize j = n;
		if (p>=0) {
			i = pos[p] / n;   // Get the index in str1
			j = n - 1 - (pos[p] % n); // Get the index in str2
		};
		// Gaps in seq 2
		while (i != position(pointerStr1)) {
			TVertexSet tmp;
			fill(tmp, (seqsInStr1 + seqsInStr2), nilVertex);
			TVertexSetIter itVEnd = end(value(pointerStr1));
			TSize count = 0;
			for(TVertexSetIter itV = begin(value(pointerStr1));itV != itVEnd;++itV) {
				tmp[count] = *itV;
				++count;
			}
			value(pointerAlign) = tmp;
			++pointerAlign;
			++pointerStr1;
		}
		// Gaps in seq 1
		while (j != position(pointerStr2)) {
			TVertexSet tmp;
			fill(tmp, (seqsInStr1 + seqsInStr2), nilVertex);
			TVertexSetIter itVEnd = end(value(pointerStr2));
			TSize count = 0;
			for(TVertexSetIter itV = begin(value(pointerStr2));itV != itVEnd;++itV) {
				tmp[count] = *itV;
				++count;
			}
			value(pointerAlign) = tmp;
			++pointerAlign;
			++pointerStr2;
		}
		if (p>=0) {
			// Matches
			TVertexSet tmp;
			fill(tmp, (seqsInStr1 + seqsInStr2), nilVertex);
			TVertexSetIter itVEnd = end(value(pointerStr1));
			TSize count = 0;
			for(TVertexSetIter itV = begin(value(pointerStr1));itV != itVEnd;++itV) {
				tmp[count] = *itV;
				++count;
			}
			TVertexSetIter itVEnd2 = end(value(pointerStr2));
			for(TVertexSetIter itV2 = begin(value(pointerStr2));itV2 != itVEnd2;++itV2) {
				tmp[count] = *itV2;
				++count;
			}
			value(pointerAlign) = tmp;
			++pointerAlign;
			++pointerStr1;
			++pointerStr2;
			--p;
		}
	}
	SEQAN_TASSERT(position(pointerAlign) == length(align))

	return score;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
