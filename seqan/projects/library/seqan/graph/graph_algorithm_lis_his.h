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
	typedef typename Iterator<TString const, Rooted>::Type TStringIter;
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
	typedef typename Iterator<TString2 const, Rooted>::Type TStringIter;
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
	typedef typename Iterator<TResult, Rooted>::Type TResultIter;
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
	typedef typename Size<TString>::Type TSize;
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
	TSize pos_of_iterator = 0;
	for(TStringIter it = begin(str); it != endIt; ++it, ++pos_of_iterator) {
		TWeight w = getValue(weights, pos_of_iterator);
		// Letters that do not contribute a weight (e.g., w = 0) are excluded!
		// Weights must increase!
		if (w <= 0) {
			addVertex(g);  // Note: The vertex id corresponds to the position
			continue;
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
				list.insert(std::make_pair(*it, std::make_pair(w, pos_of_iterator)));
		}

		// Create the corresponding node, pos_of_iterator == Vertex Descriptor
		addVertex(g);

		// Connect to predecessor
		if (a_k_it != list.end()) addEdge(g, pos_of_iterator, a_k_it->second.second);
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
			//if (getProperty(weights, v) > 0) {
				appendValue(pos, v);
				w+=getValue(weights, v);
			//}
			TOutEdgeIterator it(g, v);
			if (atEnd(it)) finished = true;
			else v = targetVertex(it);
		}
	}
	return w;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue, typename TPositions, typename TSize, typename TVertexDescriptor, typename TString>
inline void
__heaviestCommonSubsequence(std::map<TKey, TValue>&,
							TPositions const&,
							TSize const,
							TSize const,
							TVertexDescriptor const,
							TString const&, 
							TString const&,
							Nothing&) 
{
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue, typename TPositions, typename TSize, typename TVertexDescriptor, typename TString, typename TOutString>
inline void
__heaviestCommonSubsequence(std::map<TKey, TValue>& posToSlotMap,
							TPositions const& positions,
							TSize const m,
							TSize const n,
							TVertexDescriptor const nilVertex,
							TString const& str1, 
							TString const& str2,
							TOutString& align) 
{
	SEQAN_CHECKPOINT
	typedef __int64 TLargeSize;
	typedef std::map<TKey, TValue> TPositionToSlotMap;
	typedef typename Value<TString>::Type TVertexSet;
	typedef typename Iterator<TString const, Rooted>::Type TStringIter;
	typedef typename Iterator<TString, Rooted>::Type TSIter;
	typedef typename Iterator<TVertexSet const, Rooted>::Type TVertexSetIter;
	typedef typename Iterator<TVertexSet, Rooted>::Type TIter;	

	// #Vertex descriptors per node
	TSize seqsInStr1 = length(str1[0]);	 
	TSize seqsInStr2 = length(str2[0]);	

	// Reverse the map
	String<TLargeSize> slotToPos;
	resize(slotToPos, posToSlotMap.size());
	unsigned int counter = 0;
	for(typename TPositionToSlotMap::const_iterator mapIt = posToSlotMap.begin();mapIt != posToSlotMap.end(); ++mapIt, ++counter) {
		value(slotToPos, counter) = mapIt->first;
	}
	posToSlotMap.clear();

	// Create the alignment sequence
	TSize numMatches = length(positions);
	TSize alignLength = numMatches + (n - numMatches) + (m - numMatches);
	clear(align);
	resize(align, alignLength);
	TSIter pointerAlign = begin(align);
	TSIter pointerAlignEnd = end(align);
	TStringIter pointerStr1 = begin(str1);
	TStringIter pointerStr2 = begin(str2);
	int p = length(positions) - 1;
	while(pointerAlign != pointerAlignEnd) {
		TSize i = m;
		TSize j = n;
		if (p>=0) {
			i = (TSize) (slotToPos[positions[p]] / (TLargeSize) n);   // Get the index in str1
			j = n - 1 - (TSize) (slotToPos[positions[p]] % (TLargeSize) n); // Get the index in str2
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
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TString, typename TOutString>
inline TCargo
heaviestCommonSubsequence(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
						  TString const& str1, 
						  TString const& str2,
						  TOutString& align) 
{
	SEQAN_CHECKPOINT
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

	// We could create the full graph -> too expensive
	// Remember which edges are actually present
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


	// Now the tough part: Find the right number for a given position
	typedef String<TSize> TSequenceString;
	typedef typename Iterator<TSequenceString>::Type TSeqIter;
	TSequenceString seq;
	resize(seq, posToSlotMap.size());
	TSeqIter itSeq = begin(seq);
	for(typename TPositionToSlotMap::const_iterator mapIt = posToSlotMap.begin();mapIt != posToSlotMap.end(); ++mapIt, ++counter) {
		*itSeq = n - 1 - (TSize) (mapIt->first % (TLargeSize) n); ++itSeq;
	}

	// Calculate the heaviest increasing subsequence
	String<TSize> positions;
	TCargo score = (TCargo) heaviestIncreasingSubsequence(seq, weights, positions);

	// Retrieve the alignment sequence
	__heaviestCommonSubsequence(posToSlotMap, positions, m, n, nilVertex, str1, str2, align);

	return score;

}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TString>
inline TCargo
heaviestCommonSubsequence(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
						  TString const& str1, 
						  TString const& str2) 
{
	SEQAN_CHECKPOINT
	Nothing noth;
	return heaviestCommonSubsequence(g, str1, str2, noth);
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
