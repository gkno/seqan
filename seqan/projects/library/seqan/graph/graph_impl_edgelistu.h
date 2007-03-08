#ifndef SEQAN_HEADER_GRAPH_IMPL_EDGELISTU_H
#define SEQAN_HEADER_GRAPH_IMPL_EDGELISTU_H

namespace SEQAN_NAMESPACE_MAIN
{
////////////////
// Undirected graph which stores the edges in a list
////////////////

/**
.Spec.EdgeListU:
..cat:Graph
..general:Class.Graph
..summary:Undirected graph that stores the edges in an adjacency list.
..signature:Graph<EdgeListU<TCargo, TSpec>, TGraphSpec>
..param.TCargo:The cargo type that can be attached to the edges.
...metafunction:Metafunction.Cargo
...remarks:Use Cargo<EdgeListU<TCargo, TSpec> >@ to get the cargo type of an undirected graph.
...default:$void$
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Default$, see @Tag.Default@.
..include:graph.h
*/
template<typename TCargo, typename TEdgeSpec, typename TSpec>
class Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> 
{
	public:
		typedef typename Id<Graph>::Type TIdType;
		typedef typename EdgeType<Graph>::Type TEdgeStumpU;	
		typedef typename IdHandler<TEdgeStumpU, TIdType>::Type TEdgeIdManager;
		typedef Allocator<SinglePool<sizeof(TEdgeStumpU)> > TAllocator;
		
		String<TEdgeStumpU*> data_vertex;			// Pointers to EdgeStumpU lists
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
// EdgeListU specific graph functions
//////////////////////////////////////////////////////////////////////////////
template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
_copyGraph(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& source,
		   Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& dest,
		   bool transpose) 
{
	SEQAN_CHECKPOINT
	typedef Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStumpU;
	typedef typename Iterator<String<TEdgeStumpU*> const>::Type TIterConst;
	typedef typename Iterator<String<TEdgeStumpU*> >::Type TIter;
	clear(dest);
	resize(dest.data_vertex, length(source.data_vertex));
	TIter itInit = begin(dest.data_vertex);
	while(!atEnd(itInit)) {
		*itInit = (TEdgeStumpU*) 0;
		goNext(itInit);
	}
	TIterConst it = begin(source.data_vertex);
	while(!atEnd(it)) {
		TEdgeStumpU const* current = getValue(it);
		TVertexDescriptor sourceVertex = position(it);
		while(current != (TEdgeStumpU*) 0) {
			TVertexDescriptor targetVertex = getTarget(current);

			if (targetVertex != sourceVertex) {
				// Create missing vertices, targetVertex is always bigger than sourceVertex
				_createVertices(dest,targetVertex);
			
				// Add edge
				TEdgeDescriptor e = addEdge(dest, sourceVertex, targetVertex);
				_assignId(e, _getId(current));
				assignCargo(e, getCargo(current));
				current = current->data_next_source;  // Follow the link belonging to the source id
			} else {
				// Do nothing here because we don't want to create edges twice!!!
				current = current->data_next_target;
			}
		}
		goNext(it);
	}
	dest.data_id_managerV = source.data_id_managerV;
	dest.data_id_managerE = source.data_id_managerE;
}


template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
_copyGraph(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& source,
		   Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& dest) 
{
	_copyGraph(source, dest, false); // Never transpose because undirected and transposed undirected graph are equal
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
transpose(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& source,
		  Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& dest)
{
	SEQAN_CHECKPOINT
	// Undirected graph, no transpose just copy
	_copyGraph(source, dest, false);
	
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
inline typename Size<Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> >::Type 
numVertices(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& g) 
{
	SEQAN_CHECKPOINT
	return idCount(g.data_id_managerV);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline bool 
empty(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& g) 
{
	SEQAN_CHECKPOINT
	return (!(idCount(g.data_id_managerV)));
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
clearEdges(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& g) 
{
	SEQAN_CHECKPOINT
	typedef Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStumpU;
	typedef typename Iterator<String<TEdgeStumpU*> >::Type TIter;
	for(TIter it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		removeOutEdges(g,position(it)); // Remove all outgoing edges
	}
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
clearVertices(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& g) 
{
	SEQAN_CHECKPOINT
	clearEdges(g);
	releaseAll(g.data_id_managerV);
	clear(g.data_vertex);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void 
clear(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& g) 
{
	SEQAN_CHECKPOINT
	clearVertices(g);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> >::Type 
outDegree(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& g, 
		  TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStumpU;
	typedef typename Size<TGraph>::Type TSize;
	TSize count=0;
	TEdgeStumpU const* current = getValue(g.data_vertex, vertex);
	while(current!=0) {
		if ( (TVertexDescriptor) getTarget(current)==vertex) current = current->data_next_target;
		else current = current->data_next_source;
		++count;
	}
	return count;
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> >::Type 
inDegree(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& g, 
		 TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	// In-degree and out-degree are equal for undirected graphs
	return outDegree(g,vertex);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> >::Type 
degree(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& g,
	   TVertexDescriptor const vertex) 
{
	return outDegree(g,vertex);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec> 
inline typename VertexDescriptor<Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> >::Type 
addVertex(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& g) 
{
	SEQAN_CHECKPOINT
	typedef Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStumpU;
	TVertexDescriptor vd = obtainId(g.data_id_managerV);
	if (vd == length(g.data_vertex)) {
		appendValue(g.data_vertex, (TEdgeStumpU*) 0); 
	} else {
		value(g.data_vertex, vd) = (TEdgeStumpU*) 0;
	}
	return vd;
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
inline void 
removeVertex(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& g, 
			 TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)

	removeOutEdges(g,v); // Remove all outgoing edges
	releaseId(g.data_id_managerV, v); // Release id
}



template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename EdgeDescriptor<Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> >::Type 
addEdge(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& g, 
		TVertexDescriptor source, 
		TVertexDescriptor target) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(source != target)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, target) == true)

	typedef Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStumpU;
	typedef typename Id<TGraph>::Type TId;

	// Source must be the smaller vertex id
	if (source > target) {
		TVertexDescriptor tmp = target;
		target = source;
		source = tmp;
	}

	TEdgeStumpU* edge_ptr;
	allocate(g.data_allocator, edge_ptr, 1);
	valueConstruct(edge_ptr);
	assignSource(edge_ptr, source);
	assignTarget(edge_ptr, target);
	edge_ptr->data_next_source = (TEdgeStumpU*) 0;
	edge_ptr->data_next_target = (TEdgeStumpU*) 0;
	TId id = obtainId(g.data_id_managerE);
	_assignId(edge_ptr, id);
	if (getValue(g.data_vertex, source)!=0) {
		edge_ptr->data_next_source = getValue(g.data_vertex, source);
	}
	if (getValue(g.data_vertex, target)!=0) {
		edge_ptr->data_next_target = getValue(g.data_vertex, target);
	}
	value(g.data_vertex, source)=edge_ptr;
	value(g.data_vertex, target)=edge_ptr;
	return edge_ptr;
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
	TEdgeDescriptor e = addEdge(g,source,target);
	assignCargo(e,cargo);
	return e;
}


template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
void 
removeEdge(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& g, 
		   TVertexDescriptor const source, 
		   TVertexDescriptor const target) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, target) == true)
	
	typedef Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStumpU;

	// Find edge and source predecessor
	TEdgeStumpU* predSource = 0;
	TEdgeStumpU* current = getValue(g.data_vertex, source);
	while(current != (TEdgeStumpU*) 0) {
		TVertexDescriptor adjV = getTarget(current);
		if (adjV != source) {
			if ( adjV == target) break;
			predSource = current;
			current = current->data_next_source;
		} else {
			adjV = getSource(current);
			if ( adjV == target) break;
			predSource = current;
			current = current->data_next_target;
		}
	}
	
	// Not found?
	if (current == (TEdgeStumpU*) 0) return;

	// Find edge and target predecessor
	TEdgeStumpU* predTarget = 0;
	current = getValue(g.data_vertex, target);
	while(current != (TEdgeStumpU*) 0) {
		TVertexDescriptor adjV = getTarget(current);
		if (adjV != target) {
			if ( adjV == source) break;
			predTarget = current;
			current = current->data_next_source;
		} else {
			adjV = getSource(current);
			if ( adjV == source) break;
			predTarget = current;
			current = current->data_next_target;
		}
	}

	
	// Relink the next pointer of source predecessor
	if (predSource != (TEdgeStumpU*) 0) {
		if (source != (TVertexDescriptor) getTarget(current)) {
			if (source != (TVertexDescriptor) getTarget(predSource)) predSource->data_next_source = current->data_next_source;
			else predSource->data_next_target = current->data_next_source;
		} else {
			if (source != (TVertexDescriptor) getTarget(predSource)) predSource->data_next_source = current->data_next_target;
			else predSource->data_next_target = current->data_next_target;
		}
	}
	else {
		if (source != (TVertexDescriptor) getTarget(current)) value(g.data_vertex, source) = current->data_next_source;
		else value(g.data_vertex, source) = current->data_next_target;
	}

	// Relink the next pointer of target predecessor
	if (predTarget != (TEdgeStumpU*) 0) {
		if (target != (TVertexDescriptor) getTarget(current)) {
			if (target != (TVertexDescriptor) getTarget(predTarget)) predTarget->data_next_source = current->data_next_source;
			else predTarget->data_next_target = current->data_next_source;
		} else {
			if (target != (TVertexDescriptor) getTarget(predTarget)) predTarget->data_next_source = current->data_next_target;
			else predTarget->data_next_target = current->data_next_target;
		}
	}
	else {
		if (target != (TVertexDescriptor) getTarget(current)) value(g.data_vertex, target) = current->data_next_source;
		else value(g.data_vertex, target) = current->data_next_target;
	}

	// Deallocate
	releaseId(g.data_id_managerE, _getId(current));
	valueDestruct(current);
	deallocate(g.data_allocator, current, 1);	
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
	typedef Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStumpU;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TEdgeDescriptor eD = getValue(g.data_vertex, v);
	while(eD != (TEdgeStumpU*) 0) {
		TVertexDescriptor target = getTarget(eD);
		if (v == target) target = getSource(eD);
		removeEdge(g,v,target);
		eD = getValue(g.data_vertex, v);
	}
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline void 
removeInEdges(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec>& g, 
			   TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	removeOutEdges(g,v);
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
	return getSource(edge);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TMatrix>
inline void
getAdjacencyMatrix(Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& g, 
				   TMatrix& mat) 
{
	SEQAN_CHECKPOINT
	typedef Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStumpU;
	typedef typename Size<TMatrix>::Type TSize;
	TSize len = getIdUpperBound(g.data_id_managerV);
	setDimension(mat, 2);
	setLength(mat, 0, len);
	setLength(mat, 1, len);
	resize(mat);
	for (TSize i=0;i<len*len;++i) value(mat,i) = 0;

	typedef typename Iterator<String<TEdgeStumpU*> const>::Type TIterConst;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TVertexDescriptor sourceV = position(it);
		TEdgeStumpU const* current = getValue(it);
		while(current!=0) {
			TVertexDescriptor adjV = getTarget(current);
			if (adjV != sourceV) {
				assignValue(mat,sourceV*len+adjV, getValue(mat,sourceV*len+adjV)+1);
				current=current->data_next_source;
			} else {
				adjV = getSource(current);
				assignValue(mat,sourceV*len+adjV, getValue(mat,sourceV*len+adjV)+1);
				current=current->data_next_target;
			}
		}
	}
}


template <typename TFile, typename TCargo, typename TEdgeSpec, typename TSpec, typename TIDString>
inline void
write(TFile & target,
	  Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const& g,
	  TIDString const &,
	  Raw)
{
	SEQAN_CHECKPOINT
	typedef Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStumpU;
	typedef typename Iterator<String<TEdgeStumpU*> const>::Type TIterConst;
	_streamWrite(target,"Adjacency list:\n");
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TVertexDescriptor sourceV = position(it);
		_streamPutInt(target, sourceV);
		_streamWrite(target," -> ");
		TEdgeStumpU const* current = getValue(it);
		while(current!=0) {
			TVertexDescriptor adjV = getTarget(current);
			if (adjV != sourceV) {
				_streamPutInt(target, adjV);
				_streamPut(target, ',');
				current=current->data_next_source;
			} else {
				adjV = getSource(current);
				_streamPutInt(target, adjV);
				_streamPut(target, ',');
				current=current->data_next_target;
			}
		}
		_streamPut(target, '\n');
	}
	_streamWrite(target,"Edge list:\n");
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TVertexDescriptor sourceV = position(it);
		TEdgeStumpU const* current = getValue(it);
		while(current!=0) {
			TVertexDescriptor targetV = getTarget(current);
			if (sourceV != targetV) {
				_streamWrite(target,"Source: ");
				_streamPutInt(target, sourceV);		
				_streamPut(target, ',');
				_streamWrite(target,"Target: ");
				_streamPutInt(target, targetV);
				_streamPut(target, ' ');
				_streamWrite(target,"(Id: ");
				_streamPutInt(target, _getId(current));
				_streamPut(target, ',');
				_streamWrite(target," Cargo-Type: ");
				_streamWrite(target, typeid(getCargo(current)).name());
				_streamPut(target, ')');
				_streamPut(target, '\n');
				current=current->data_next_source;
			} else {
				current=current->data_next_target;
			}
		}
	}
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
