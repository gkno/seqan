#ifndef SEQAN_HEADER_GRAPH_IMPL_EDGELIST_H
#define SEQAN_HEADER_GRAPH_IMPL_EDGELIST_H

namespace SEQAN_NAMESPACE_MAIN
{
////////////////
// Graph which stores the edges in a list
////////////////
template<typename TCargo, typename TEdgeSpec, typename TSpec>
class Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> 
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
// EdgeList specific graph functions
//////////////////////////////////////////////////////////////////////////////
template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
_copyGraph(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> const& source,
		   Graph<EdgeList<TCargo, TEdgeSpec>, TSpec>& dest,
		   bool transpose) 
{
	SEQAN_CHECKPOINT
	typedef Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> const>::Type TIterConst;
	typedef typename Iterator<String<TEdgeStump*> >::Type TIter;
	clear(dest);
	resize(dest.data_vertex, length(source.data_vertex));
	TIter itInit = begin(dest.data_vertex);
	while(!atEnd(itInit)) {
		*itInit = (TEdgeStump*) 0;
		goNext(itInit);
	}
	TIterConst it = begin(source.data_vertex);
	while(!atEnd(it)) {
		TEdgeStump const* current = getValue(it);
		TVertexDescriptor sourceVertex = position(it);
		while(current != (TEdgeStump*) 0) {
			TVertexDescriptor targetVertex = current->data_target;
			// Create missing vertices
			if (sourceVertex>targetVertex) _createVertices(dest,sourceVertex);
			else _createVertices(dest,targetVertex);
			// Add edge
			TEdgeDescriptor e;
			if (!transpose) e = addEdge(dest, sourceVertex, targetVertex);
			else e = addEdge(dest, targetVertex, sourceVertex);
			_assignId(e, _getId(current));
			assignCargo(e, getCargo(current));
			current = current->data_next;
		}
		goNext(it);
	}
	dest.data_id_managerV = source.data_id_managerV;
	dest.data_id_managerE = source.data_id_managerE;
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
_copyGraph(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> const& source,
				   Graph<EdgeList<TCargo, TEdgeSpec>, TSpec>& dest) 
{
	_copyGraph(source, dest, false);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline typename Size<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> >::Type 
numEdges(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> const& g) 
{
	SEQAN_CHECKPOINT
	return idCount(g.data_id_managerE);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
clearEdges(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec>& g) 
{
	SEQAN_CHECKPOINT
	typedef Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> >::Type TIter;
	for(TIter it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TEdgeStump* current = getValue(it);
		if(current == (TEdgeStump*) 0) continue;
		TVertexDescriptor sourceVertex = position(it);
		removeOutEdges(g, sourceVertex);
	}
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
clearVertices(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec>& g) 
{
	SEQAN_CHECKPOINT
	clearEdges(g);
	releaseAll(g.data_id_managerV);
	clear(g.data_vertex);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> >::Type 
outDegree(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> const& g, 
		  TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> const TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Size<TGraph>::Type TSize;
	TSize count=0;
	TEdgeStump const* current = getValue(g.data_vertex, vertex);
	while(current!=0) {
		current = current->data_next;
		++count;
	}
	return count;
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> >::Type 
inDegree(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> const& g, 
		 TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Iterator<String<TEdgeStump*> const>::Type TIterConst;
	TSize count=0;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TEdgeStump const* current = getValue(it);
		while(current!=0) {
			if ( (TVertexDescriptor) current->data_target==vertex) ++count;
			current = current->data_next;
		}
	}
	return count;
}

template<typename TCargo, typename TEdgeSpec, typename TSpec> 
inline typename VertexDescriptor<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> >::Type 
addVertex(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec>& g) 
{
	SEQAN_CHECKPOINT	
	typedef Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	TVertexDescriptor vd = obtainId(g.data_id_managerV);
	if (vd == length(g.data_vertex)) {
		appendValue(g.data_vertex, (TEdgeStump*) 0); 
	} else {
		value(g.data_vertex, vd) = (TEdgeStump*) 0;
	}
	return vd;
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename EdgeDescriptor<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> >::Type 
addEdge(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec>& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, target) == true)

	typedef Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Id<TGraph>::Type TId;

	TEdgeStump* edge_ptr;
	allocate(g.data_allocator, edge_ptr, 1);
	valueConstruct(edge_ptr);
	edge_ptr->data_target = target;
	edge_ptr->data_next = (TEdgeStump*) 0;
	TId id = obtainId(g.data_id_managerE);
	_assignId(edge_ptr, id);
	if (getValue(g.data_vertex, source)!=0) {
		edge_ptr->data_next = getValue(g.data_vertex, source);
	}
	value(g.data_vertex, source)=edge_ptr;
	return edge_ptr;
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename EdgeDescriptor<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> >::Type 
addEdge(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec>& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		TCargo const cargo) 
{
	SEQAN_CHECKPOINT
	typedef Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TEdgeDescriptor e = addEdge(g,source,target);
	assignCargo(e,cargo);
	return e;
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
void 
removeEdge(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec>& g, 
		   TVertexDescriptor const source, 
		   TVertexDescriptor const target) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, target) == true)
	
	typedef Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;

	// Find edge and predecessor
	TEdgeStump* pred = 0;
	TEdgeStump* current = getValue(g.data_vertex, source);
	while(current != (TEdgeStump*) 0) {
		if ( (TVertexDescriptor) current->data_target == target) break;
		pred = current;
		current = current->data_next;
	}
	
	// Not found?
	if (current == (TEdgeStump*) 0) return;
	
	// Relink the next pointer of predecessor
	if (pred != (TEdgeStump*) 0) pred->data_next = current->data_next;
	else value(g.data_vertex, source) = current->data_next;
	
	// Deallocate
	releaseId(g.data_id_managerE, _getId(current));
	valueDestruct(current);
	deallocate(g.data_allocator, current, 1);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeDescriptor>
void 
removeEdge(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec>& g,
		   TEdgeDescriptor const edge)
{
	SEQAN_CHECKPOINT
	removeEdge(g, sourceVertex(g,edge), targetVertex(g,edge));
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline void 
removeOutEdges(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec>& g, 
			   TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)

	typedef Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	while(getValue(g.data_vertex, v) != (TEdgeStump*) 0) {
		TVertexDescriptor target = targetVertex(g,(getValue(g.data_vertex, v)));
		removeEdge(g,v,target);
	}
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline void 
removeInEdges(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec>& g, 
			   TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)

	typedef Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> >::Type TIter;
	for(TIter it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TEdgeStump* current = getValue(it);
		TVertexDescriptor const sourceVertex = position(it);
		while(current!=0) {
			if ( (TVertexDescriptor) current->data_target==v) {
				removeEdge(g, sourceVertex, v);
				current = getValue(g.data_vertex, sourceVertex);
			} else {
				current = current->data_next;
			}
		}
	}
}


template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> >::Type 
targetVertex(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> const& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	return getTarget(edge);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> >::Type 
sourceVertex(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> const& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	typedef Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> const>::Type TIterConst;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TEdgeDescriptor current = getValue(it);
		while(current!=(TEdgeDescriptor) 0) {
			if (current == edge) return position(it);
			current=current->data_next;
		}
	}
	// We should never reach this point
	SEQAN_ASSERT(false)
	return 0;
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor, typename TChar>
inline typename VertexDescriptor<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> >::Type 
getSuccessor(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> const& g,
				   TVertexDescriptor vertex,
				   TChar const c) 
{
	SEQAN_CHECKPOINT
	typedef Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Cargo<TEdgeStump>::Type TInternalCargo;
	InternalMap<TInternalCargo> eMap;
	initEdgeMap(g,eMap);
	return getSuccessor(g,vertex,eMap,c);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeLabelMap, typename TVertexDescriptor, typename TChar>
inline typename VertexDescriptor<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> >::Type 
getSuccessor(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> const& g,
				   TVertexDescriptor vertex,
				   TEdgeLabelMap const& eMap,
				   TChar const c) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	TEdgeStump const* current = getValue(g.data_vertex, vertex);
	while(current!=(TEdgeStump*) 0) {
		if (getProperty(eMap,current) == c) return targetVertex(g,current);
		current=current->data_next;
	}
	return _get_nil<TVertexDescriptor>();
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor, typename TChar>
inline typename VertexDescriptor<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> >::Type 
getPredecessor(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> const& g,
				   TVertexDescriptor vertex,
				   TChar const c) 
{
	SEQAN_CHECKPOINT
	typedef Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
   	typedef typename Cargo<TEdgeStump>::Type TInternalCargo;
	InternalMap<TInternalCargo> eMap;
	initEdgeMap(g,eMap);
	return getPredecessor(g,vertex,eMap,c);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeLabelMap, typename TVertexDescriptor, typename TChar>
inline typename VertexDescriptor<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> >::Type 
getPredecessor(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> const& g,
				   TVertexDescriptor vertex,
				   TEdgeLabelMap const& eMap,
				   TChar const c) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	
	typedef Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> const>::Type TIterConst;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TEdgeStump const* current = getValue(it);
		while(current!=(TEdgeStump*) 0) {
			if ((vertex = targetVertex(g,current)) &&
				(getProperty(eMap,current) == c)) {
					return position(it);
			}
			current=current->data_next;
		}
	}
	// We should never reach this point
	SEQAN_ASSERT(false)
	return 0;
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TMatrix>
inline void
getAdjacencyMatrix(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> const& g, 
				   TMatrix& mat) 
{
	SEQAN_CHECKPOINT

	typedef Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Size<TMatrix>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<String<TEdgeStump*> const>::Type TIterConst;
	TSize len = getIdUpperBound(g.data_id_managerV);
	setDimension(mat, 2);
	setLength(mat, 0, len);
	setLength(mat, 1, len);
	resize(mat);
	for (TSize i=0;i<len*len;++i) value(mat,i) = 0;
	
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TEdgeStump const* current = getValue(it);
		TVertexDescriptor const source = position(it);
		while(current != (TEdgeStump*) 0) {
			TVertexDescriptor target = targetVertex(g,current);
			assignValue(mat, source*len+target, getValue(mat, source*len+target)+1);
			current = current->data_next;
		}
	}
}

template <typename TFile, typename TCargo, typename TEdgeSpec, typename TSpec, typename TIDString>
inline void
write(TFile & target,
	  Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> const& g,
	  TIDString const &,
	  Raw)
{
	SEQAN_CHECKPOINT
	typedef Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> const>::Type TIterConst;
	_streamWrite(target,"Adjacency list:\n");
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		_streamPutInt(target, position(it));
		_streamWrite(target," -> ");
		TEdgeStump const* current = getValue(it);
		while(current!=0) {
			_streamPutInt(target, current->data_target);
			_streamPut(target, ',');
			current=current->data_next;
		}
		_streamPut(target, '\n');
	}
	_streamWrite(target,"Edge list:\n");
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TEdgeStump const* current = getValue(it);
		while(current!=0) {
			_streamWrite(target,"Source: ");
			_streamPutInt(target, position(it));		
			_streamPut(target, ',');
			_streamWrite(target,"Target: ");
			_streamPutInt(target, current->data_target);
			_streamPut(target, ' ');
			_streamWrite(target,"(Id: ");
			_streamPutInt(target, _getId(current));
			_streamPut(target, ',');
			_streamWrite(target," Cargo-Type: ");
			_streamWrite(target, typeid(getCargo(current)).name());
			_streamPut(target, ')');
			_streamPut(target, '\n');
			current=current->data_next;
		}
	}
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
