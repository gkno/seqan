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
  $Id: graph_algorithm_lis_his.h 1811 2008-03-31 15:38:54Z rausch@PCPOOL.MI.FU-BERLIN.DE $
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
		if (a_k_it != list.end()) addEdge(g, (TVertexDescriptor) position(it), (TVertexDescriptor) a_k_it->second);
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
template<typename TString1, typename TString2, typename TNeighborhoodSize, typename TFinalPos>
inline void
longestCommonSubsequence(TString1 const& str1,
						 TString2 const& str2,
						 TNeighborhoodSize nSize,
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
	TPos current_pos = 0;
	for(TStringIter it = begin(str2); it != endIt; ++it, ++current_pos) {
		push_back(value(occ, ordValue(value(it))), current_pos);
	}

	// Build the combined string
	String<TPos, Block<> > finalSeq;
	String<TPos, Block<> > mapping;
	TStringIter endIt1 = end(str1);
	current_pos = 0;
	for(TStringIter it = begin(str1); it != endIt1; ++it, ++current_pos) {
		TPositions& current_occ = occ[ordValue(value(it))];
		for(int i = length(current_occ)-1; i>=0; --i) {
			// Do we have a neighborhood
			if (nSize != 0) {
				TPos diff = current_pos - current_occ[i];
				if (current_pos < current_occ[i]) diff = current_occ[i] - current_pos;
				if ((TPos) diff > (TPos) nSize) continue;
			}
			push_back(finalSeq, current_occ[i]);
			push_back(mapping, current_pos);
		}
	}

	// Call longest increasing subsequence
	typedef String<TSize, Block<> > TResult;
	TResult result;
	longestIncreasingSubsequence(finalSeq, result);
	
	// Insert the common pairs
	typedef typename Iterator<TResult>::Type TResultIter;
	TResultIter endResult = end(result);
	for(TResultIter it = begin(result); it != endResult; ++it) {
		push_back(pos, std::make_pair(mapping[*it], finalSeq[*it]));
	}

	//// Debug code
	//for(TSize i=0; i<length(pos);++i) {
	//	std::cout << pos[i].first << ',' << pos[i].second << ';';
	//}
	//std::cout << std::endl;
}

template<typename TString1, typename TString2, typename TFinalPos>
inline void
longestCommonSubsequence(TString1 const& str1,
						 TString2 const& str2,
						 TFinalPos& pos) 
{
	SEQAN_CHECKPOINT
	longestCommonSubsequence(str1, str2, 0, pos);
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
	typedef typename Iterator<TString const, Standard>::Type TStringIter;
	TStringIter it = begin(str, Standard());
	TStringIter endIt = end(str, Standard());
	TSize pos_of_iterator = 0;
	for(; it != endIt; ++it, ++pos_of_iterator) {
		TWeight w = weights[pos_of_iterator];
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
		if (a_k_it != list.end()) addEdge(g, (TVertexDescriptor) pos_of_iterator, (TVertexDescriptor) a_k_it->second.second);
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



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
