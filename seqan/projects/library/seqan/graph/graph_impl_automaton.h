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
			clear(*this);
		}

		Graph(Graph const & _other) 
		{
			SEQAN_CHECKPOINT
			_copy(_other, *this);
		}
	
		Graph const& operator = (Graph const & _other) {
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			_copy(_other, *this);
			return *this;
		}
};

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
_copy(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& source,
	  Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& dest,
	  bool transpose = false)
{
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TAlphabet>::Type TSize;
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const>::Type TIterConst;
	
	clear(dest);
	dest.data_id_managerV = source.data_id_managerV;
	resize(dest.data_vertex, length(source.data_vertex));
	if (!transpose) {
		dest.data_vertex = source.data_vertex;
	} else {
		TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
		for(TIterConst it = begin(source.data_vertex);!atEnd(it);goNext(it)) {
			TSize table_length = ValueSize<TAlphabet>::VALUE;
			TVertexDescriptor sourceVertex = position(it);
			for(TSize i=0;i<table_length;++i) {
				TVertexDescriptor targetVertex = source.data_vertex[sourceVertex].data_edge[i].data_target;
				if (targetVertex == nilVal) continue;
				dest.data_vertex[targetVertex].data_edge[i].data_target = sourceVertex;
				assignCargo(dest.data_vertex[targetVertex].data_edge[i], getCargo(dest.data_vertex[sourceVertex].data_edge[i]));
			}
		}
	}
} 

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
clearEdges(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g)
{
	SEQAN_CHECKPOINT
	clear(g.data_vertex);
	resize(g.data_vertex, getIdUpperBound(g.data_id_managerV));
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
clearVertices(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g)
{
	SEQAN_CHECKPOINT
	releaseAll(g.data_id_managerV);
	clear(g.data_vertex);
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
outDegree(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
		  TVertexDescriptor const vertex)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename Size<TGraph>::Type TSize;
	TSize count=0;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	for(TSize i=0;i<table_length;++i) {
		if (g.data_vertex[vertex].data_edge[i].data_target!=nilVal) ++count;
	}
	return count;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
inDegree(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
		 TVertexDescriptor const vertex)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TGraph>::Type TSize;
	TSize count=0;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const>::Type TIterConst;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		if (!idInUse(g.data_id_managerV, position(it))) continue;
		for(TSize i=0;i<table_length;++i) {
			if (((getValue(it)).data_edge[i].data_target!=nilVal) &&
				((getValue(it)).data_edge[i].data_target==vertex))
			{
				++count;
			}
		}
	}
	return count;
}


template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
inline typename Size<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
numEdges(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TGraph>::Type TSize;
	TSize count=0;
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const>::Type TIterConst;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		if (!idInUse(g.data_id_managerV, position(it))) continue;
		count += outDegree(g, position(it));
	}
	return count;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
addVertex(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph >::Type TEdge;

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
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TAlphabet>::Type TSize;
	TAlphabet letter(label);
	g.data_vertex[source].data_edge[(TSize) letter].data_target = target;
	return TEdgeDescriptor(source,label);
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
addEdge(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		TAlphabet const label,
		TCargo const cargo)
{
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TEdgeDescriptor e = addEdge(g,source,target,label);
	assignCargo(e,cargo);
	return e;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
sourceVertex(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, edge.i1) == true)
	return edge.i1;
}

//Slow
template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
inline void
removeEdge(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g,
	       TVertexDescriptor const source,
	       TVertexDescriptor const target)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, target) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename Size<TGraph>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	for(TSize i=0;i<table_length;++i) {
		if (g.data_vertex[source].data_edge[i].data_target==target) {
			g.data_vertex[source].data_edge[i].data_target=nilVal;
			// Don't return here: Maybe multiple edges with different labels but the same target
		}
	}
}

//Fast
template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeDescriptor>
inline void
removeEdge(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g,
	   TEdgeDescriptor const edge)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, edge.i1) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	TAlphabet letter(edge.i2);
	g.data_vertex[edge.i1].data_edge[(TSize) letter].data_target=nilVal;
}


template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
inline void
removeOutEdges(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g,
			   TVertexDescriptor const vertex)
{
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename Size<TGraph>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	for(TSize i=0;i<table_length;++i) {
		g.data_vertex[vertex].data_edge[i].data_target=nilVal;
	}
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
inline void
removeInEdges(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g,
			  TVertexDescriptor const vertex)
{
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TGraph>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > >::Type TIter;
	for(TIter it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		if (!idInUse(g.data_id_managerV, position(it))) continue;
		for(TSize i=0;i<table_length;++i) {
			if (((getValue(it)).data_edge[i].data_target!=nilVal) &&
				((getValue(it)).data_edge[i].data_target==vertex))
			{
				(getValue(it)).data_edge[i].data_target=nilVal;
			}
		}
	}
}


template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
targetVertex(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, edge.i1) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename Size<TAlphabet>::Type TSize;
	TAlphabet letter(edge.i2);
	return g.data_vertex[edge.i1].data_edge[(TSize) letter].data_target;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TMatrix>
inline void
getAdjacencyMatrix(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
				   TMatrix& mat) 
{
	SEQAN_CHECKPOINT

	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Size<TMatrix>::Type TMatrixSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	TMatrixSize len = getIdUpperBound(g.data_id_managerV);
	setDimension(mat, 2);
	setLength(mat, 0, len);
	setLength(mat, 1, len);
	resize(mat);
	for (TMatrixSize i=0;i<len*len;++i) value(mat,i) = 0;
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const>::Type TIterConst;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		if (!idInUse(g.data_id_managerV, position(it))) continue;
		for(TSize i=0;i<table_length;++i) {
			if (((getValue(it)).data_edge[i].data_target!=nilVal))
			{
				TVertexDescriptor const source = position(it);
				TVertexDescriptor const target = (getValue(it)).data_edge[i].data_target;
				assignValue(mat, source*len+target, getValue(mat, source*len+target)+1);
			}
		}
	}
}


template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor, typename TChar>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
getSuccessorVertex(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
				   TVertexDescriptor vertex,
				   TChar const c) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	typedef Graph<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> > TGraph;
	typedef typename Size<TAlphabet>::Type TSize;
	TAlphabet letter(c);
	return g.data_vertex[vertex].data_edge[(TSize) letter].data_target;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor, typename TChar>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
getPredecessorVertex(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
					 TVertexDescriptor vertex,
					 TChar const c) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TGraph>::Type TSize;
	TAlphabet letter(c);
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const>::Type TIterConst;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		if (((getValue(it)).data_edge[(TSize) letter].data_target!=nilVal) &&
			((getValue(it)).data_edge[(TSize) letter].data_target==vertex))
		{
			return position(it);
	    }
	}
	// We should never reach this point
	SEQAN_ASSERT(false)
	return 0;
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
		if (!idInUse(g.data_id_managerV, position(it))) continue;
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
