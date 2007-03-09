#ifndef SEQAN_HEADER_GRAPH_IMPL_TRIE_H
#define SEQAN_HEADER_GRAPH_IMPL_TRIE_H

namespace SEQAN_NAMESPACE_MAIN
{

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TTerminalStateMap, typename TKeyword, typename TPos>
inline void
_addStringToTrie(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g,
				 TTerminalStateMap& terminalStateMap,
				 TKeyword const& str,
				 TPos const& keywordIndex)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TKeyword>::Type TSize;

	TVertexDescriptor current = getRoot(g);
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	typename Iterator<TKeyword const>::Type sIt = begin(str);
	for(;!atEnd(sIt);goNext(sIt)) {
		if (getSuccessor(g, current, *sIt) == nilVal) break;
		current = getSuccessor(g, current, *sIt);
	}
	for(;!atEnd(sIt);goNext(sIt)) {
		TVertexDescriptor newState = addVertex(g);
		resize(terminalStateMap, numVertices(g), Generous());
		assignProperty(terminalStateMap,newState,String<TPos>());
		addEdge(g,current,newState,*sIt);
		current = newState;
	}
	String<TPos> tmp = getProperty(terminalStateMap,current);
	appendValue(tmp, keywordIndex);
	assignProperty(terminalStateMap,current,tmp);
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TTerminalStateMap, typename TKeywords>
inline void
createTrie(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g,
		   TTerminalStateMap& terminalStateMap,
		   TKeywords const& keywords)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Position<TKeywords>::Type TPos;
	TVertexDescriptor root = addVertex(g);
	assignRoot(g,root);
	resize(terminalStateMap, numVertices(g), Generous());
	assignProperty(terminalStateMap,root,String<TPos>());
	typename Iterator<TKeywords const>::Type it = begin(keywords);
	for(;!atEnd(it);goNext(it)) {
		_addStringToTrie(g,terminalStateMap,*it,position(it));
	}
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TTerminalStateMap, typename TKeywords>
inline void
createTrieOnReverse(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g,
					TTerminalStateMap& terminalStateMap,
					TKeywords const& keywords)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Position<TKeywords>::Type TPos;
	TVertexDescriptor root = addVertex(g);
	assignRoot(g,root);
	resize(terminalStateMap, numVertices(g), Generous());
	assignProperty(terminalStateMap,root,String<TPos>());
	typename Iterator<TKeywords const>::Type it = begin(keywords);
	for(;!atEnd(it);goNext(it)) {
		typedef typename Value<TKeywords>::Type TKeyword;
		TKeyword tmp;
		typename Iterator<TKeyword const>::Type sIt = end(*it);
		while(!atBegin(sIt)) {
			goPrevious(sIt);
			appendValue(tmp,*sIt);
		}
		_addStringToTrie(g,terminalStateMap,tmp,position(it));
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
