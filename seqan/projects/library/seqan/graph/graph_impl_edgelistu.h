#ifndef SEQAN_HEADER_GRAPH_IMPL_EDGELISTU_H
#define SEQAN_HEADER_GRAPH_IMPL_EDGELISTU_H

namespace SEQAN_NAMESPACE_MAIN
{
////////////////
// Undirected graph which stores the edges in a list
////////////////
template<typename TCargo, typename TEdgeSpec, typename TSpec>
class Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> 
{
	public:
		typedef typename Id<Graph>::Type TIdType;
		typedef typename EdgeType<Graph>::Type TEdgeStump;	
		typedef typename IdHandler<TEdgeStump, TIdType>::Type TEdgeIdManager;
		typedef Allocator<SinglePool<sizeof(TEdgeStump)> > TAllocator;
		
		String<TEdgeStump*> data_vertex;			// Pointers to EdgeStump lists
		IdManager<TIdType> data_id_managerV;
		TEdgeIdManager data_id_managerE;		
		TAllocator data_allocator;

//____________________________________________________________________________


		Graph() {
			SEQAN_CHECKPOINT
		}


		template<typename TEdgeArray, typename TSize>
		Graph(TEdgeArray edges, TSize size) {
			SEQAN_CHECKPOINT
			_copyGraph(edges,size,*this);
		}


		~Graph() {
			SEQAN_CHECKPOINT
			clear(*this);
		}

		Graph(Graph const & _other) :
			data_allocator(_other.data_allocator)
		{
			SEQAN_CHECKPOINT
			_copyGraph(_other, *this);		
		}
	
		Graph const& operator = (Graph const & _other) {
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			clear(*this);
			data_allocator = _other.data_allocator;
			_copyGraph(_other, *this);
			return *this;
		}
};


//////////////////////////////////////////////////////////////////////////////
// EdgeListU specific graph functions
//////////////////////////////////////////////////////////////////////////////
template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
_copyGraph(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& source,
		   Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& dest,
		   bool transpose) 
{
	SEQAN_CHECKPOINT
	// Never transpose because both graphs are equal (for undirected graphs)
	_copyEdgeListGraph(source,dest,false);
}


template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
_copyGraph(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& source,
		   Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& dest) 
{
	_copyGraph(source, dest, false);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
transpose(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& g)
{
	SEQAN_CHECKPOINT
	// Nothing to do
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
transpose(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& g)
{
	// Nothing to do
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline typename Size<Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> >::Type 
numEdges(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& g) 
{
	SEQAN_CHECKPOINT
	return idCount(g.data_id_managerE);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
clearEdges(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& g) 
{
	SEQAN_CHECKPOINT
	typedef Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> >::Type TIter;
	for(TIter it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TEdgeStump* current = getValue(it);
		if(current == (TEdgeStump*) 0) continue;
		TVertexDescriptor sourceV = position(it);
		while(getValue(g.data_vertex, sourceV) != (TEdgeStump*) 0) {
			TVertexDescriptor targetV = targetVertex(g,(getValue(g.data_vertex, sourceV)));
			removeEdge(g,sourceV,targetV);
		}
	}
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
clearVertices(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& g) 
{
	SEQAN_CHECKPOINT
	_clearEdgeListVertices(g);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> >::Type 
outDegree(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& g, 
		  TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	return (_inEdgeListDegree(g,vertex) + _outEdgeListDegree(g,vertex));
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> >::Type 
inDegree(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& g, 
		 TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	return (_inEdgeListDegree(g,vertex) + _outEdgeListDegree(g,vertex));
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> >::Type 
degree(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& g,
	   TVertexDescriptor const vertex) 
{
	return outDegree(g,vertex);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> >::Type 
degree(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& g,
	   TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	return outDegree(g,vertex);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec> 
inline typename VertexDescriptor<Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> >::Type 
addVertex(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& g) 
{
	SEQAN_CHECKPOINT
	return _addEdgeListVertex(g);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename EdgeDescriptor<Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> >::Type 
addEdge(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(source != target)
	if (source > target) return _addDirectedEdge(g,target,source);
	else return _addDirectedEdge(g,source,target);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename EdgeDescriptor<Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> >::Type 
addEdge(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		TCargo const cargo) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(source != target)
	typedef Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TEdgeDescriptor e;
	if (source > target) {
		e = _addDirectedEdge(g,target,source);
		assignCargo(e,cargo);
	} else {
		e = _addDirectedEdge(g,source,target);
		assignCargo(e,cargo);
	}
	return e;
}


template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
void 
removeEdge(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& g, 
		   TVertexDescriptor const source, 
		   TVertexDescriptor const target) 
{
	SEQAN_CHECKPOINT
	if (source > target) _removeDirectedEdge(g,target,source);
	else _removeDirectedEdge(g,source,target);
	
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeDescriptor>
void 
removeEdge(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& g,
		   TEdgeDescriptor const edge)
{
	SEQAN_CHECKPOINT
	removeEdge(g, sourceVertex(g,edge), targetVertex(g,edge));
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline void 
removeOutEdges(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& g, 
			   TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	_removeEdgeListInEdges(g,v);
	_removeEdgeListOutEdges(g,v);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline void 
removeInEdges(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& g, 
			   TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	_removeEdgeListInEdges(g,v);
	_removeEdgeListOutEdges(g,v);
}





template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> >::Type 
targetVertex(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	return getTarget(edge);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> >::Type 
sourceVertex(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	return _sourceEdgeListVertex(g,edge);
}

template <typename TFile, typename TCargo, typename TEdgeSpec, typename TSpec, typename TIDString>
inline void
write(TFile & target,
	  Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& g,
	  TIDString const &id,
	  Raw)
{
	SEQAN_CHECKPOINT
	_writeEdgeList(target, g, id, Raw());
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TMatrix>
inline void
getAdjacencyMatrix(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& g, 
				   TMatrix& mat) 
{
	SEQAN_CHECKPOINT
	typedef typename Size<TMatrix>::Type TSize;
	_getEdgeListAdjacencyMatrix(g,mat);
	TSize len = getIdUpperBound(g.data_id_managerV);
	for (TSize row=0;row<len;++row) {
		for (TSize col=0;col<row+1;++col) {
			value(mat,row*len+col) = getValue(mat,col*len+row);
		}
	}
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
