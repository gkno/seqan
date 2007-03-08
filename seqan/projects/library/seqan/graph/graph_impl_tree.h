#ifndef SEQAN_HEADER_GRAPH_IMPL_TREE_H
#define SEQAN_HEADER_GRAPH_IMPL_TREE_H

namespace SEQAN_NAMESPACE_MAIN
{
////////////////
// Tree that stores the tree edges in a list
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
class Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> 
{
	public:
		typedef typename Id<Graph>::Type TIdType;
		typedef typename EdgeType<Graph>::Type TEdgeStumpT;	
		typedef typename VertexDescriptor<Graph>::Type TVertexDescriptor;
		typedef Allocator<SinglePool<sizeof(TEdgeStumpT)> > TAllocator;
		
		String<TEdgeStumpT*> data_vertex;			// Pointers to EdgeStumpT lists
		IdManager<TIdType> data_id_managerV;
		TAllocator data_allocator;
		TVertexDescriptor root;

//____________________________________________________________________________


		Graph() : root(0) {
			SEQAN_CHECKPOINT
			addVertex(*this); // Add the root node
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
// EdgeListT specific graph functions
//////////////////////////////////////////////////////////////////////////////
template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
_copyGraph(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const& source,
		   Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec>& dest,
		   bool transpose) 
{
	SEQAN_CHECKPOINT
	typedef Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStumpT;
	typedef typename Iterator<String<TEdgeStumpT*> const>::Type TIterConst;
	typedef typename Iterator<String<TEdgeStumpT*> >::Type TIter;
	clear(dest);
	resize(dest.data_vertex, length(source.data_vertex));
	TIter itInit = begin(dest.data_vertex);
	while(!atEnd(itInit)) {
		*itInit = (TEdgeStumpT*) 0;
		goNext(itInit);
	}
	TIterConst it = begin(source.data_vertex);
	while(!atEnd(it)) {
		TEdgeStumpT const* current = getValue(it);
		TVertexDescriptor parentVertex = position(it);
		while(current != (TEdgeStumpT*) 0) {
			TVertexDescriptor childVertex = current->data_child;
			// Create missing vertices
			if (parentVertex > childVertex) _createVertices(dest,parentVertex);
			else _createVertices(dest,childVertex);
			// Add edge
			TEdgeDescriptor e = addEdge(dest, parentVertex, childVertex);
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
_copyGraph(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const& source,
		   Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec>& dest) 
{
	// To transpose a tree doesn't make sense
	_copyGraph(source, dest, false); 
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
transpose(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const& source,
		  Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec>& dest)
{
	SEQAN_CHECKPOINT
	// On a tree transpose makes no sense, just copy
	_copyGraph(source, dest, false);
	
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
transpose(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const& g)
{
	// Nothing to do
}


template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline typename Size<Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> >::Type 
numEdges(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const& g) 
{
	SEQAN_CHECKPOINT
	return (idCount(g.data_id_managerV) - 1);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline typename Size<Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> >::Type 
numVertices(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const& g) 
{
	SEQAN_CHECKPOINT
	return idCount(g.data_id_managerV);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline bool 
empty(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const& g) 
{
	SEQAN_CHECKPOINT
	// Just the root
	return (idCount(g.data_id_managerV) == 1);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
clearEdges(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec>& g) 
{
	SEQAN_CHECKPOINT
	removeAllChildren(g, getRoot(g));
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
clearVertices(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec>& g) 
{
	SEQAN_CHECKPOINT
	// No need to release ids because clearEdges removes the children
	clearEdges(g);
	clear(g.data_vertex);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline void 
clear(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec>& g) 
{
	SEQAN_CHECKPOINT
	clearVertices(g);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> >::Type 
outDegree(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const& g, 
		  TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	return (numChildren(g,vertex) + 1);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> >::Type 
inDegree(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const& g, 
		 TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	
	// Trees are considered as undirected graphs
	return (outDegree(g,vertex));
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> >::Type 
degree(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const& g,
	   TVertexDescriptor const vertex) 
{
	return outDegree(g,vertex);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec> 
inline typename VertexDescriptor<Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> >::Type 
addVertex(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec>& g) 
{
	SEQAN_CHECKPOINT	
	typedef Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStumpT;
	TVertexDescriptor vd = obtainId(g.data_id_managerV);
	if (vd == length(g.data_vertex)) {
		appendValue(g.data_vertex, (TEdgeStumpT*) 0); 
	} else {
		value(g.data_vertex, vd) = (TEdgeStumpT*) 0;
	}
	return vd;
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
inline void 
removeVertex(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec>& g, 
			 TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)
	removeChild(g,getParent(v),v);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename EdgeDescriptor<Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> >::Type 
addEdge(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec>& g, 
		TVertexDescriptor const parent, 
		TVertexDescriptor const child) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, parent) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, child) == true)

	typedef Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStumpT;
	typedef typename Id<TGraph>::Type TId;

	TEdgeStumpT* edge_ptr;
	allocate(g.data_allocator, edge_ptr, 1);
	valueConstruct(edge_ptr);
	assignParent(edge_ptr, parent);
	assignChild(edge_ptr, child);
	edge_ptr->data_next = (TEdgeStumpT*) 0;
	if (getValue(g.data_vertex, parent)!=0) {
		edge_ptr->data_next = getValue(g.data_vertex, parent);
	}
	value(g.data_vertex, parent)=edge_ptr;
	return edge_ptr;
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename EdgeDescriptor<Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> >::Type 
addEdge(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec>& g, 
		TVertexDescriptor const parent, 
		TVertexDescriptor const child,
		TCargo const cargo) 
{
	SEQAN_CHECKPOINT
	typedef Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TEdgeDescriptor e = addEdge(g,parent,child);
	assignCargo(e,cargo);
	return e;
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
void 
removeEdge(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec>& g,
		   TVertexDescriptor const parent,
		   TVertexDescriptor const child) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, parent) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, child) == true)
	
	typedef Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStumpT;

	// Find edge and predecessor
	TEdgeStumpT* pred = 0;
	TEdgeStumpT* current = getValue(g.data_vertex, parent);
	while(current != (TEdgeStumpT*) 0) {
		if ( (TVertexDescriptor) current->data_child == child) break;
		pred = current;
		current = current->data_next;
	}
	
	// Not found?
	if (current == (TEdgeStumpT*) 0) return;
	
	// Relink the next pointer of predecessor
	if (pred != (TEdgeStumpT*) 0) pred->data_next = current->data_next;
	else value(g.data_vertex, parent) = current->data_next;
	
	// Deallocate
	valueDestruct(current);
	deallocate(g.data_allocator, current, 1);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeDescriptor>
void 
removeEdge(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec>& g,
		   TEdgeDescriptor const edge)
{
	SEQAN_CHECKPOINT
	removeEdge(g, parentVertex(g,edge), childVertex(g,edge));
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline void 
removeOutEdges(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec>& g, 
			   TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)
	removeAllChildren(g, getParent(v));
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline void 
removeInEdges(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec>& g,
			  TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	removeOutEdges(g,v);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> >::Type 
targetVertex(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	return childVertex(g,edge);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> >::Type 
sourceVertex(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	return parentVertex(g,edge);
}

template <typename TFile, typename TCargo, typename TEdgeSpec, typename TSpec, typename TIDString>
inline void
write(TFile & target,
	  Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const& g,
	  TIDString const &,
	  Raw)
{
	SEQAN_CHECKPOINT
	typedef Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> const>::Type TIterConst;
	_streamWrite(target,"Adjacency list:\n");
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		_streamPutInt(target, position(it));
		_streamWrite(target," -> ");
		TEdgeStump const* current = getValue(it);
		while(current!=0) {
			_streamPutInt(target, getChild(current));
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
			_streamPutInt(target, getChild(current));
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

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline typename VertexDescriptor<Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> >::Type 
getRoot(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const& g)
{
	SEQAN_CHECKPOINT
	return g.root;
}

template<typename TCargo, typename TEdgeSpec, typename TSpec>
inline typename VertexDescriptor<Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> >::Type&
root(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec>& g)
{
	SEQAN_CHECKPOINT
	return g.root;
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
inline bool
isRoot(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const& g,
	   TVertexDescriptor v)
{
	SEQAN_CHECKPOINT
	return ( (TVertexDescriptor) g.root == v);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
inline bool
isLeaf(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const& g,
	   TVertexDescriptor v)
{
	SEQAN_CHECKPOINT
	typedef Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStumpT;
	return (value(g.data_vertex, v) ==  (TEdgeStumpT*) 0);
}


template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename VertexDescriptor<Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> >::Type 
addChild(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec>& g,
		 TVertexDescriptor parent) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, parent) == true)
	TVertexDescriptor child = addVertex(g);
	addEdge(g,parent,child);
	return child;
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename VertexDescriptor<Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> >::Type 
addChild(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec>& g,
		 TVertexDescriptor const parent,
		 TCargo const cargo) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, parent) == true)
	TVertexDescriptor child = addVertex(g);
	addEdge(g,parent,child,cargo);
	return child;
}



template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline void 
removeAllChildren(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec>& g, 
				  TVertexDescriptor const parent) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, parent) == true)

	typedef Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStumpT;
	while(getValue(g.data_vertex, parent) != (TEdgeStumpT*) 0) {
		TVertexDescriptor child = childVertex(g,(getValue(g.data_vertex, parent)));
		if (!isLeaf(g,child)) removeAllChildren(g,child);
		removeEdge(g,parent, child);
		releaseId(g.data_id_managerV, child); // Release id
	}
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline void 
removeChild(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec>& g,
			TVertexDescriptor const parent,
			TVertexDescriptor const child)
{
	SEQAN_CHECKPOINT
	if (!isLeaf(g,child)) removeAllChildren(g,child);
	removeEdge(g,parent, child);
	releaseId(g.data_id_managerV, child); // Release id
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> >::Type 
childVertex(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const& g,
			TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	return getChild(edge);
}

template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> >::Type 
parentVertex(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	return getParent(edge);
}



template<typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> >::Type 
numChildren(Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const& g,
			TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStumpT;
	typedef typename Size<TGraph>::Type TSize;
	TSize count=0;
	TEdgeStumpT const* current = getValue(g.data_vertex, vertex);
	while(current!=0) {
		current = current->data_next;
		++count;
	}
	return count;
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
