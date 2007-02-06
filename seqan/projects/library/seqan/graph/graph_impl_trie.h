#ifndef SEQAN_HEADER_GRAPH_IMPL_TRIE_H
#define SEQAN_HEADER_GRAPH_IMPL_TRIE_H

namespace SEQAN_NAMESPACE_MAIN
{

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TPos, typename TChar>
inline void
createTrie(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g,
		   String<String<TPos> >& terminalStateMap,
		   String<String<TChar> > const& keywords)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<String<TChar> >::Type TSize;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	TVertexDescriptor root = addVertex(g);
	assignRoot(g,root);
	typename Iterator<String<String<TChar> > const>::Type it = begin(keywords);
	for(;!atEnd(it);goNext(it)) {
		TVertexDescriptor current = root;
		typename Iterator<String<TChar> const>::Type sIt = begin(*it);
		for(;!atEnd(sIt);goNext(sIt)) {
			if (getSuccessor(g, current, *sIt) == nilVal) break;
			current = getSuccessor(g, current, *sIt);
		}
		for(;!atEnd(sIt);goNext(sIt)) {
			TVertexDescriptor newState = addVertex(g);
			addEdge(g,current,newState,*sIt);
			current = newState;
		}
		if (length(terminalStateMap) <= current) {
			resize(terminalStateMap,current+1);
		}
		if (empty(getProperty(terminalStateMap,current))) {
			assignProperty(terminalStateMap,current,String<TPos>(position(it)));
		} else {
			String<TPos> tmp = getProperty(terminalStateMap,current);
			appendValue(tmp, position(it));
			assignProperty(terminalStateMap,current,tmp);
		}
	}
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TNodeMap>
inline void
_createTrieNodeNames(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
					 String<String<unsigned int> > pos,
					 TNodeMap& nodeMap)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	initVertexMap(g, nodeMap);
	typedef typename Iterator<TGraph, VertexIterator<> >::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		std::stringstream s;
		s << "\"" << *it << " {";
		String<unsigned int> endPositions = getProperty(pos,*it);
		if (!empty(endPositions)) {
			typename Iterator<String<unsigned int> >::Type itP = begin(endPositions);
			typename Iterator<String<unsigned int> >::Type beginP = itP;
			for(;!atEnd(itP);goNext(itP)) {
				if (beginP != itP) s << ", ";
				s << *itP;
			}
		}
		s << "}" << "\"";
		assignProperty(nodeMap, *it, String<char>(s.str().c_str()));
	}

}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
