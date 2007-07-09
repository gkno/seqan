#ifndef SEQAN_HEADER_GRAPH_ALGORITHM_LIS_HIS_H
#define SEQAN_HEADER_GRAPH_ALGORITHM_LIS_HIS_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TWeightMap, typename TPositions>
inline typename Value<TWeightMap>::Type
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
	TWeight w = 0;
	if (list.rbegin() == list.rend()) return 0;
	else {
		bool finished = false;
		TVertexDescriptor v = list.rbegin()->second.second;
		while (!finished) {
			TOutEdgeIterator it(g, v);
			appendValue(pos, v);
			w+=getProperty(weights, v);
			if (atEnd(it)) finished = true;
			else v = targetVertex(it);
		}
	}
	return w;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TString>
inline void
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

	typedef std::map<TVertexDescriptor, TSize> TVertexToPosMap;
	typedef typename TVertexToPosMap::const_iterator TVertexToPosMapIter;
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	TSize m = length(str1);
	TSize n = length(str2);
	TSize seqsInStr1 = length(str1[0]);	
	TSize seqsInStr2 = length(str2[0]);	

	// Fill the vertex to position map for str1
	TVertexToPosMap map;
	TStringIter itStrEnd1 = end(str1);
	for(TStringIter itStr1 = begin(str1);itStr1 != itStrEnd1;++itStr1) {
		TVertexSetIter itVEnd = end(getValue(itStr1));
		for(TVertexSetIter itV = begin(getValue(itStr1));itV != itVEnd;++itV) {
			if (*itV != nilVertex) {
				map.insert(std::make_pair(*itV, position(itStr1)));
			}
		}
	}

	// Create the sequence and the corresponding weights
	String<TSize> seq;
	resize(seq, n*m);
	typedef typename Iterator<String<TSize> >::Type TSeqIter;
	TSeqIter itSeq = begin(seq);
	for(TSize i=0; i<m;++i) {
		for(TSize j=n;j>0;--j) {
			*itSeq = j - 1;
			++itSeq;
		}
	}
	String<double> weights;
	double divider = (double) seqsInStr1 * (double) seqsInStr2;
	fill(weights, n*m, 0);

	// Walk through str2 and fill in the weights
	TStringIter itStrEnd2 = end(str2);
	for(TStringIter itStr2 = begin(str2);itStr2 != itStrEnd2;++itStr2) {
		TVertexSetIter itVEnd = end(getValue(itStr2));
		for(TVertexSetIter itV = begin(getValue(itStr2));itV != itVEnd;++itV) {
			if (*itV != nilVertex) {
				TOutEdgeIterator itOut(g, *itV);
				for(;!atEnd(itOut); ++itOut) {
					TVertexToPosMapIter pPos = map.find(targetVertex(itOut));
					if (pPos != map.end()) {
						TSize index = pPos->second * n + (n - position(itStr2) - 1);
						weights[index] += (double) cargo(*itOut) / divider;
					}
				}
			}
		}
	}
	map.clear();

	String<unsigned int> pos;
	heaviestIncreasingSubsequence(seq, weights, pos);
	
	//// Debug code
	//for(TSize i = 0; i<n*m; ++i) {
	//	std::cout << seq[i] << ':' << weights[i] << ',';
	//}
	//std::cout << std::endl;
	//for(int i = length(pos)-1; i>=0; --i) {
	//	std::cout << pos[i] <<  ',';
	//}
	//std::cout << std::endl;
	//std::cout << w << std::endl;

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
			i = pos[p] / n;
			j = n - 1 - (pos[p] % n); 
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

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
