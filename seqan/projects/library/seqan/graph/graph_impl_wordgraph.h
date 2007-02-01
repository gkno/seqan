#ifndef SEQAN_HEADER_GRAPH_IMPL_WORDGRAPH_H
#define SEQAN_HEADER_GRAPH_IMPL_WORDGRAPH_H

namespace SEQAN_NAMESPACE_MAIN
{

template <typename TSpec = Default>
struct WordGraph;

////////////////
// WordGraph
////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec>
class Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> 
{
	public:
		typedef typename Id<Graph>::Type TIdType;
		typedef typename EdgeType<Graph>::Type TEdge;

		String<AutomatonEdgeArray<TEdge, TAlphabet> > data_vertex;		// List of tables
		IdManager<TIdType> data_id_managerV;
		String<String<TAlphabet> > data_edge_label;
	

//____________________________________________________________________________


		Graph() {
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

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec>
inline void
_copyGraph(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> const& source,
		   Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec>& dest,
		   bool transpose = false)
{
	_copyAutomatonGraph(source, dest, transpose);
	resize(dest.data_edge_label,length(source.data_edge_label));
	if (transpose == false) {
		dest.data_edge_label=source.data_edge_label;
	} else {
		typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> TGraph;
		typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
		typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
		typedef typename EdgeType<TGraph>::Type TEdge;
		typedef typename Size<TAlphabet>::Type TSize;
		TSize table_length = ValueSize<TAlphabet>::VALUE;
		TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
		typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const>::Type TIterConst;
		for(TIterConst it = begin(source.data_vertex);!atEnd(it);goNext(it)) {
			if (!idInUse(source.data_id_managerV, position(it))) continue;
			TVertexDescriptor sourceVertex = position(it);
			for(TSize i=0;i<table_length;++i) {
				TVertexDescriptor targetVertex = (getValue(it)).data_edge[i].data_target;
				if (targetVertex==nilVal) continue;
				assignProperty(dest.data_edge_label, targetVertex * table_length + i, getProperty(source.data_edge_label, sourceVertex * table_length +i));
			}
		}
	}
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec>
inline void
clearEdges(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec>& g)
{
	SEQAN_CHECKPOINT
	_clearAutomatonEdges(g);
	clear(g.data_edge_label);
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec>
inline void
clearVertices(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec>& g)
{
	SEQAN_CHECKPOINT
	_clearAutomatonVertices(g);
	clear(g.data_edge_label);
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> >::Type 
addVertex(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec>& g) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> TGraph;
	typedef typename Size<TGraph>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	resize(g.data_edge_label, (numVertices(g)+1) * table_length);
	return _addAutomatonVertex(g);
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> >::Type 
addEdge(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec>& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		String<TAlphabet> const& label) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TAlphabet>::Type TSize;
	TAlphabet firstChar = getValue(label, 0);
	TEdgeDescriptor e = _addAutomatonEdge(g,source,target, firstChar);
	assignProperty(g.data_edge_label, e, suffix(label,1));
	return e;
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec, typename TVertexDescriptor, typename TLabel>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> >::Type 
addEdge(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec>& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		TLabel const label) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TAlphabet>::Type TSize;
	String<TAlphabet> tmp(label);
	return addEdge(g, source, target, tmp);
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec, typename TVertexDescriptor, typename TEdgeCargo>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> >::Type 
addEdge(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec>& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		String<TAlphabet> const& label,
		TEdgeCargo const cargo)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TAlphabet>::Type TSize;
	TAlphabet firstChar = getValue(label, 0);
	TEdgeDescriptor e = _addAutomatonEdge(g,source,target, firstChar, (TCargo) cargo);
	assignProperty(g.data_edge_label, e, suffix(label,1));
	return e;
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
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TAlphabet>::Type TSize;
	String<TAlphabet> tmp(label);
	return addEdge(g,source,target,tmp,cargo);
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec, typename TVertexDescriptor, typename TLabel>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> >::Type 
getSuccessor(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> const& g,
			 TVertexDescriptor vertex,
			 TLabel const label) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TAlphabet>::Type TSize;
	TAlphabet letter(getValue(label, 0));
	if (getProperty(g.data_edge_label, TEdgeDescriptor(vertex, letter)) == suffix(label, 1)) {
		return g.data_vertex[vertex].data_edge[(TSize) letter].data_target;
	} else {
		return _get_nil<TVertexDescriptor>();
	}
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
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const>::Type TIterConst;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		if ((getValue(it)).data_edge[(TSize) letter].data_target==nilVal) continue;
		String<TAlphabet> edgeLabel(getProperty(g.data_edge_label, TEdgeDescriptor(position(it), letter)));
		if (((getValue(it)).data_edge[(TSize) letter].data_target==vertex) &&
			(edgeLabel == suf))
		{
			return position(it);
	    }
	}
	return _get_nil<TVertexDescriptor>();
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
			if (g.data_vertex[sourceVertex].data_edge[i].data_target ==  nilVal) continue;
			_streamPutInt(target, sourceVertex);
			_streamWrite(target,"->");
			_streamPutInt(target, g.data_vertex[sourceVertex].data_edge[i].data_target);
			_streamPut(target, ' ');
			_streamPut(target, ' ');
			_streamWrite(target, "Label: ");
			_streamPut(target, TAlphabet(i));
			_streamWrite(target, getProperty(g.data_edge_label, TEdgeDescriptor(sourceVertex, TAlphabet(i))));
			_streamPut(target, '\n');
		}
	}
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
