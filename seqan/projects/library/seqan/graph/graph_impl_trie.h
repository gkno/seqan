#ifndef SEQAN_HEADER_GRAPH_IMPL_TRIE_H
#define SEQAN_HEADER_GRAPH_IMPL_TRIE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph - Trie
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Trie:
..cat:Graph
..general:Class.Graph
..summary:A keyword trie.
..remarks:A keyword trie is a special automaton and thus, it is not implemented in its own class.
It solely provides create functions where based upon a set of strings a keyword trie is created.
..signature:Graph<Automaton<TAlphabet, TCargo, TSpec> > 
..param.TAlphabet:The alphabet type that is used for the transition labels.
...metafunction:Metafunction.Alphabet
...remarks:Use @Metafunction.Alphabet@ to get the type of the labels in an automaton.
...default:$char$
..param.TCargo:The cargo type that can be attached to the edges.
...metafunction:Metafunction.Cargo
...remarks:Use @Metafunction.Cargo@ to get the cargo type of an undirected graph.
...default:$void$
..param.TSpec:The specializing type for the graph.
...metafunction:Metafunction.Spec
...remarks:Use WithoutEdgeId here to omit edge ids.
Note: If edges do not store ids external property maps do not work.
...default:$Default$, see @Tag.Default@.
..include:graph.h
*/


//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TTerminalStateMap, typename TKeyword, typename TPos>
inline void
_addStringToTrie(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
				 TTerminalStateMap& terminalStateMap,
				 TKeyword const& str,
				 TPos const& keywordIndex)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
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

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TNodeMap>
inline void
_createTrieNodeNames(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
					 String<String<unsigned int> > pos,
					 TNodeMap& nodeMap)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	initVertexMap(g, nodeMap);
	typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
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

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.createTrie:
..cat:Spec.Trie
..summary:Creates a trie.
..signature:createTrie(g, terminalStateMap, keywords)
..param.g:Out-parameter: An automaton.
...type:Spec.Trie
..param.terminalStateMap:Out-parameter: An external property map.
...type:Class.External Property Map
...remarks:The external property map must be a String<String<unsigned int> > because
in every vertex of the trie a number of keywords can end. 
Hence, we must associate with every vertex a set of indices that correspond to keywords.
..param.text:In-parameter: A set of strings.
...type:Class.String
..returns:void
..see:Function.createTrieOnReverse
*/
template<typename TAlphabet, typename TCargo, typename TSpec, typename TTerminalStateMap, typename TKeywords>
inline void
createTrie(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
		   TTerminalStateMap& terminalStateMap,
		   TKeywords const& keywords)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
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

//////////////////////////////////////////////////////////////////////////////

/**
.Function.createTrieOnReverse:
..cat:Spec.Trie
..summary:Creates a trie for all reversed keywords.
..signature:createTrieOnReverse(g, terminalStateMap, keywords)
..param.g:Out-parameter: An automaton.
...type:Spec.Trie
..param.terminalStateMap:Out-parameter: An external property map.
...type:Class.External Property Map
...remarks:The external property map must be a String<String<unsigned int> > because
in every vertex of the trie a number of keywords can end. 
Hence, we must associate with every vertex a set of indices that correspond to keywords.
..param.text:In-parameter: A set of strings.
...type:Class.String
..returns:void
..see:Function.createTrie
*/
template<typename TAlphabet, typename TCargo, typename TSpec, typename TTerminalStateMap, typename TKeywords>
inline void
createTrieOnReverse(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
					TTerminalStateMap& terminalStateMap,
					TKeywords const& keywords)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
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

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
