#ifndef SEQAN_HEADER_GRAPH_ALGORITHM_LIS_HIS_H
#define SEQAN_HEADER_GRAPH_ALGORITHM_LIS_HIS_H

namespace SEQAN_NAMESPACE_MAIN
{

template<typename TSortedSequence, typename TKey>
inline typename TSortedSequence::const_iterator
_previousInSortedSequence(TSortedSequence const& list, TKey const key) {
	SEQAN_CHECKPOINT
	typedef typename TSortedSequence::const_iterator TSortedSequenceIter;

	TSortedSequenceIter a_k_it = list.lower_bound(key);
	if (a_k_it != list.end()) {
		if (a_k_it == list.begin()) a_k_it = list.end();
		else --a_k_it;
	} else {
		TSortedSequenceIter tmp = list.begin();
		if (tmp != list.end()) {
			do {
				a_k_it = tmp;
			} while(++tmp != list.end());
		}
	}
	return a_k_it;
}

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

template<typename TString, typename TPositions>
inline void
longestIncreasingSubsequence(TString const& str, TPositions& pos) {
	SEQAN_CHECKPOINT
	typedef typename Size<TString>::Type TSize;
	typedef typename Value<TString>::Type TValue;
	typedef typename Value<TPositions>::Type TPos;
	typedef typename Iterator<TString const>::Type TStringIter;
	typedef std::pair<TValue, TPos> TKey;
	typedef std::set<TKey, std::less<TKey> > TSortedSequence;
	typedef typename TSortedSequence::const_iterator TSortedSequenceIter;
	typedef Graph<Directed<void, WithoutEdgeId> > TGraph;
	typedef Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	
	TSortedSequence list;
	TGraph g;
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
		addVertex(g);

		// Connect to predecessor
		if (a_k_it != list.end()) addEdge(g, position(it), a_k_it->second);
	}

	// Trace-back
	if (list.rbegin() == list.rend()) return;
	else {
		bool finished = false;
		TVertexDescriptor v = list.rbegin()->second;
		while (!finished) {
			TOutEdgeIterator it(g, v);
			appendValue(pos, v);
			if (atEnd(it)) finished = true;
			else v = targetVertex(it);
		}
	}
}

template<typename TString, typename TWeightMap, typename TPositions>
inline void
heaviestIncreasingSubsequence(TString const& str, TWeightMap const& weights, TPositions& pos) {
	SEQAN_CHECKPOINT
	typedef typename Size<TString>::Type TSize;
	typedef typename Value<TString>::Type TValue;
	typedef typename Value<TPositions>::Type TPos;
	typedef typename Value<TWeightMap>::Type TWeight;
	typedef typename Iterator<TString const>::Type TStringIter;
	typedef std::pair<TValue, std::pair<TWeight, TPos> > TKey;
	typedef std::set<TKey, std::less<TKey> > TSortedSequence;
	typedef typename TSortedSequence::const_iterator TSortedSequenceIter;
	typedef Graph<Directed<void, WithoutEdgeId> > TGraph;
	typedef Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;

	TSortedSequence list;
	TGraph g;
	TStringIter endIt = end(str);
	for(TStringIter it = begin(str); it != endIt; ++it) {
		// Get previous element
		TSortedSequenceIter a_k_it = _previousInSortedSequence(list, std::make_pair(*it, std::make_pair(0, 0))); 
		// Get next element
		TSortedSequenceIter b_l_it = _nextInSortedSequence(list, a_k_it);

		// Determine new weight
		TWeight w = getProperty(weights, position(it));
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

		// Create the corresponding node
		addVertex(g);

		// Connect to predecessor
		if (a_k_it != list.end()) addEdge(g, position(it), a_k_it->second.second);
	}

	// Trace-back
	if (list.rbegin() == list.rend()) return;
	else {
		bool finished = false;
		TVertexDescriptor v = list.rbegin()->second.second;
		while (!finished) {
			TOutEdgeIterator it(g, v);
			appendValue(pos, v);
			if (atEnd(it)) finished = true;
			else v = targetVertex(it);
		}
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
