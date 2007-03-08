#ifndef SEQAN_HEADER_GRAPH_IMPL_WORDGRAPH_H
#define SEQAN_HEADER_GRAPH_IMPL_WORDGRAPH_H

namespace SEQAN_NAMESPACE_MAIN
{

template <typename TSpec = Default>
struct WordGraph;

////////////////
// WordGraph
////////////////

template<typename TAlphabet, typename TSpec, typename TGraphSpec>
class Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> >, TGraphSpec> 
{
	public:
		typedef typename Id<Graph>::Type TIdType;
		typedef typename EdgeType<Graph>::Type TEdge;
		typedef typename VertexDescriptor<Graph>::Type TVertexDescriptor;
		typedef typename IdHandler<TEdge, TIdType>::Type TEdgeIdManager;

		String<AutomatonEdgeArray<TEdge, TAlphabet> > data_vertex;		// List of tables
		IdManager<TIdType> data_id_managerV;
		TEdgeIdManager data_id_managerE;
		TVertexDescriptor root;


//____________________________________________________________________________


		Graph() : root(0) {
			SEQAN_CHECKPOINT
		}


		~Graph() {
			SEQAN_CHECKPOINT
			clear(*this);
		}

		Graph(Graph const & _other) 
		{
			SEQAN_CHECKPOINT
			_copyGraph(_other, *this);
		}
	
		Graph const& operator = (Graph const & _other) {
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			_copyGraph(_other, *this);
			return *this;
		}
};

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> >::Type 
addEdge(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec>& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		String<TAlphabet> const& label) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, target) == true)
	
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Size<TAlphabet>::Type TSize;
	TAlphabet firstChar = getValue(label, 0);
	TEdgeDescriptor e = &g.data_vertex[source].data_edge[(TSize) firstChar];
	TId id = obtainId(g.data_id_managerE);
	_assignId(e, id);
	assignTarget(e, target);
	String<TAlphabet> suf(suffix(label,1));
	assignCargo(e, suf);
	return e;
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec, typename TVertexDescriptor, typename TChars>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> >::Type 
addEdge(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec>& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		TChars const* chars) 
{
	SEQAN_CHECKPOINT
	return addEdge(g,source,target,String<TAlphabet>(chars));
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec, typename TVertexDescriptor, typename TLabel, typename TEdgeCargo>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> >::Type 
addEdge(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec>& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		TLabel const label,
		TEdgeCargo const cargo)
{
	SEQAN_CHECKPOINT
	// No additional cargo allowed. Cargo is used for string
	// Use external property map
	SEQAN_ASSERT(false)
}


template<typename TFile, typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec, typename TIDString>
inline void
write(TFile & target,
	  Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> const& g,
	  TIDString const &,
	  Raw)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TAlphabet>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();

	_streamWrite(target,"WordGraph - EdgeList:\n");
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const>::Type TIterConst;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		if (!idInUse(g.data_id_managerV, position(it))) continue;
		TVertexDescriptor sourceVertex = position(it);
		for(TSize i=0;i<table_length;++i) {
			TEdgeDescriptor ed = &g.data_vertex[sourceVertex].data_edge[i];
			if (getTarget(ed) ==  nilVal) continue;
			_streamPutInt(target, sourceVertex);
			_streamWrite(target,"->");
			_streamPutInt(target, getTarget(ed));
			_streamPut(target, ' ');
			_streamPut(target, ' ');
			_streamWrite(target, "Label: ");
			_streamPut(target, TAlphabet(i));
			_streamWrite(target, getCargo(ed));
			_streamPut(target, '\n');
		}
	}
}


template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec, typename TVertexDescriptor, typename TCharacters>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> >::Type 
getSuccessor(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> const& g,
			 TVertexDescriptor vertex,
			 TCharacters const& chars)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> const TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TAlphabet>::Type TSize;
	TAlphabet letter(getValue(chars, 0));
	TEdgeDescriptor ed = &g.data_vertex[vertex].data_edge[(TSize) letter];
	if (getCargo(ed) == suffix(chars, 1)) {
		return getTarget(ed);
	} else {
		return _get_nil<TVertexDescriptor>();
	}
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec, typename TVertexDescriptor, typename TCharacters>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> >::Type 
getSuccessor(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> const& g,
			 TVertexDescriptor vertex,
			 TCharacters const* chars)
{
	SEQAN_CHECKPOINT
	return getSuccessor(g,vertex,String<TAlphabet>(chars));
}


template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec, typename TVertexDescriptor, typename TLabel>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> >::Type 
getPredecessor(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> const& g,
			   TVertexDescriptor vertex,
			   TLabel const label) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TGraph>::Type TSize;
	TAlphabet letter(getValue(label, 0));
	String<TAlphabet> suf(suffix(label, 1));
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const>::Type TIterConst;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TEdgeDescriptor ed = (TEdgeDescriptor) &(getValue(it)).data_edge[(TSize) letter];
		if ( (TVertexDescriptor) getTarget(ed)==nilVal) continue;
		String<TAlphabet> tmp(getCargo(ed));
		if (( (TVertexDescriptor) getTarget(ed)==vertex) &&
			(tmp == suf))
		{
			return position(it);
	    }
	}
	return _get_nil<TVertexDescriptor>();
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec, typename TVertexDescriptor, typename TCharacters>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> >::Type 
getPredecessor(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> const& g,
			   TVertexDescriptor vertex,
			   TCharacters const* chars)
{
	SEQAN_CHECKPOINT
	return getPredecessor(g,vertex,String<TAlphabet>(chars));
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec, typename TVertexDescriptor, typename TIterator>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> >::Type 
parseString(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> const& g,
			TVertexDescriptor const vertex,
			TIterator beginIt,
			TIterator endIt)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> TGraph;
	typedef typename Size<TGraph>::Type TSize;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	TVertexDescriptor succ = vertex;
	while (beginIt!=endIt) {
		String<TAlphabet> label(*beginIt);
		TSize range = 1;
		TVertexDescriptor tmp = getSuccessor(g,succ,label);
		while ((tmp == nilVal) &&
				(beginIt+range != endIt))
		{
			appendValue(label, *(beginIt + range));
			tmp = getSuccessor(g,succ,label);
			++range;
		}
		if (tmp == nilVal) break;
		succ = tmp;
		beginIt = beginIt+range;
	}
	return succ;
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
