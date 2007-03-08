#ifndef SEQAN_HEADER_GRAPH_IMPL_ORACLE_H
#define SEQAN_HEADER_GRAPH_IMPL_ORACLE_H

namespace SEQAN_NAMESPACE_MAIN
{

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TPropertyMap, typename TChar>
inline void
addLetterToOracle(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g,
				  TPropertyMap& supplyState,
				  TChar const c)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	TVertexDescriptor newState = addVertex(g);
	TVertexDescriptor pred = newState - 1;
	addEdge(g, pred, newState, c);
	TVertexDescriptor k = getProperty(supplyState, pred);
	while ((k!=nilVal) &&
			(getTarget(&g.data_vertex[k].data_edge[(TSize) TAlphabet(c)])==nilVal))
	{
		addEdge(g,k,newState,c);
		k = getProperty(supplyState, k);
	}
	TVertexDescriptor s;
	if (k==nilVal) s=0;
	else s = getTarget(&g.data_vertex[k].data_edge[(TSize) TAlphabet(c)]);
	assignProperty(supplyState, newState, s);
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TText>
inline void
createOracle(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g,
			 TText const text)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Size<TText>::Type TSize;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	TSize len = length(text);
	String<TVertexDescriptor> supplyState;
	resize(supplyState, len+1);
	TVertexDescriptor v1 = addVertex(g);
	assignProperty(supplyState, v1, nilVal);
	for(TSize i = 0; i<len; ++i) {
		addLetterToOracle(g, supplyState, getValue(text,i));
	}
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TText>
inline void
createOracleOnReverse(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g,
			 TText const text)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Size<TText>::Type TSize;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	TSize len = length(text);
	String<TVertexDescriptor> supplyState;
	resize(supplyState, len+1);
	TVertexDescriptor v1 = addVertex(g);
	assignProperty(supplyState, v1, nilVal);
	for(TSize i = len-1; i>0; --i) {
		addLetterToOracle(g, supplyState, getValue(text,i));
	}
	addLetterToOracle(g, supplyState, getValue(text,0));
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
