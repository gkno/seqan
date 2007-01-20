#ifndef SEQAN_HEADER_GRAPH_IMPL_AUTOMATON_H
#define SEQAN_HEADER_GRAPH_IMPL_AUTOMATON_H

namespace SEQAN_NAMESPACE_MAIN
{

template<typename TEdge, typename TAlphabet>
class AutomatonEdgeArray {
public:
	TEdge data_edge[ValueSize<TAlphabet>::VALUE];

	AutomatonEdgeArray() 
	{
		typedef typename VertexDescriptor<TEdge>::Type TVertexDescriptor;
		typedef typename Size<TAlphabet>::Type TSize;
		TSize table_length = ValueSize<TAlphabet>::VALUE;
		TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
		for(TSize i=0;i<table_length;++i) {
			data_edge[i].data_target=nilVal;
		}
	}
};

////////////////
// Automaton
////////////////
template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
class Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> 
{
	public:
		typedef typename Id<Graph>::Type TIdType;
		typedef typename EdgeType<Graph>::Type TEdge;

		String<AutomatonEdgeArray<TEdge, TAlphabet> > data_vertex;		// List of tables
		IdManager<TIdType> data_id_managerV;
	

//____________________________________________________________________________


		Graph() {
			SEQAN_CHECKPOINT
		}


		~Graph() {
			SEQAN_CHECKPOINT
		}

		Graph(Graph const & _other) 
		{
			SEQAN_CHECKPOINT
	
		}
	
		Graph const& operator = (Graph const & _other) {
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			return *this;
		}
};

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
addVertex(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g) 
{
	SEQAN_CHECKPOINT
	typedef typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type TVertexDescriptor;
	typedef typename EdgeType<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type TEdge;

	TVertexDescriptor vd = obtainId(g.data_id_managerV);
	if (vd == length(g.data_vertex)) {
		appendValue(g.data_vertex, AutomatonEdgeArray<TEdge, TAlphabet>()); 
	} else {
		value(g.data_vertex, vd) =  AutomatonEdgeArray<TEdge, TAlphabet>();
	}
	return vd;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor, typename TLabel>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
addEdge(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		TLabel const label) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, target) == true)
	typedef typename Size<TAlphabet>::Type TSize;
	typedef typename EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type TEdgeDescriptor;
	TAlphabet letter = label;
	g.data_vertex[source].data_edge[(TSize) letter].data_target = target;
	return TEdgeDescriptor(source,label);
}

template<typename TFile, typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TIDString>
inline void
write(TFile & target,
	  Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
	  TIDString const &,
	  Raw)
{
	SEQAN_CHECKPOINT
	typedef typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type TVertexDescriptor;
	typedef typename EdgeType<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type TEdge;
	typedef typename Size<TAlphabet>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;

	_streamWrite(target,"Automaton - State: (Input / NextState)\n");
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const>::Type TIterConst;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TVertexDescriptor sourceVertex = position(it);
		_streamPutInt(target, sourceVertex);
		_streamWrite(target,": ");
		for(TSize i=0;i<table_length;++i) {
			_streamPut(target, ' ');
			_streamPut(target, '(');
			_streamPut(target, TAlphabet(i));
			_streamPut(target, ' ');
			_streamPut(target, '/');
			_streamPut(target, ' ');
			if (g.data_vertex[sourceVertex].data_edge[i].data_target ==  _get_nil<TVertexDescriptor>()) _streamWrite(target,"nil");
			else _streamPutInt(target, g.data_vertex[sourceVertex].data_edge[i].data_target);
			_streamPut(target, ')');
			_streamPut(target, ' ');
		}
		_streamPut(target, '\n');
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
