#ifndef SEQAN_HEADER_GRAPH_IMPL_TREE_H
#define SEQAN_HEADER_GRAPH_IMPL_TREE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph - Tree
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Tree:
..cat:Graph
..general:Class.Graph
..summary:A tree that stores the edges in an adjacency list.
..signature:Graph<Tree<TCargo, TSpec> >
..param.TCargo:The cargo type that can be attached to the edges.
...metafunction:Metafunction.Cargo
...remarks:Use @Metafunction.Cargo@ to get the cargo type of the tree.
...default:$void$
..param.TSpec:The specializing type for the graph.
...metafunction:Metafunction.Spec
...default:$Default$, see @Tag.Default@.
..include:graph.h
*/
template<typename TCargo, typename TSpec>
class Graph<Tree<TCargo, TSpec> > 
{
	public:
		typedef typename VertexIdHandler<Graph>::Type TVertexIdManager;
		typedef typename VertexDescriptor<Graph>::Type TVertexDescriptor;
		typedef typename EdgeType<Graph>::Type TEdgeStump;	
		typedef Allocator<SinglePool<sizeof(TEdgeStump)> > TAllocator;
		
		TVertexDescriptor data_root;
		String<TEdgeStump*> data_vertex;			// Pointers to EdgeStumpT lists
		String<TVertexDescriptor> data_parent;		// Map to the parents of each node
		TVertexIdManager data_id_managerV;
		TAllocator data_allocator;
		

//____________________________________________________________________________


		Graph() : data_root(getNil<TVertexDescriptor>()) {
			SEQAN_CHECKPOINT
		}


		~Graph() {
			SEQAN_CHECKPOINT
			clear(*this);
		}

		Graph(Graph const & _other) :
			data_root(getNil<TVertexDescriptor>()),
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
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline String<typename EdgeType<Graph<Tree<TCargo, TSpec> > >::Type*>&
_getVertexString(Graph<Tree<TCargo, TSpec> > const& g) {
	SEQAN_CHECKPOINT
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	return const_cast<String<TEdgeStump*>&>(g.data_vertex);
}

/////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline IdManager<typename Id<Graph<Tree<TCargo, TSpec> > >::Type, Default> const &
_getVertexIdManager(Graph<Tree<TCargo, TSpec> > const& g) {
	return g.data_id_managerV;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline IdManager<typename Id<Graph<Tree<TCargo, TSpec> > >::Type, Default> &
_getVertexIdManager(Graph<Tree<TCargo, TSpec> >& g) {
	return g.data_id_managerV;
}

/////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline IdManager<typename Id<Graph<Tree<TCargo, TSpec> > >::Type, Default> const &
_getEdgeIdManager(Graph<Tree<TCargo, TSpec> > const& g) {
	return g.data_id_managerV;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline IdManager<typename Id<Graph<Tree<TCargo, TSpec> > >::Type, Default> &
_getEdgeIdManager(Graph<Tree<TCargo, TSpec> >& g) {
	return g.data_id_managerV;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec> 
inline void
_rebuildParentMap(Graph<Tree<TCargo, TSpec> > const& g) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> const>::Type TIterConst;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TVertexDescriptor parent = position(it);
		TEdgeStump* current = getValue(it);
		while(current!=0) {
			assignValue(g.data_parent, getTarget(current), parent);
			current = getNextT(current);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Tree<TCargo, TSpec> > const& source,
		   Graph<Tree<TCargo, TSpec> >& dest,
		   bool transpose) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> const>::Type TIterConst;
	typedef typename Iterator<String<TEdgeStump*> >::Type TIter;
	clear(dest);
	resize(dest.data_vertex, length(source.data_vertex));
	resize(dest.data_parent, length(source.data_parent));
	TIter itInit = begin(dest.data_vertex);
	while(!atEnd(itInit)) {
		*itInit = (TEdgeStump*) 0;
		goNext(itInit);
	}
	TIterConst it = begin(source.data_vertex);
	while(!atEnd(it)) {
		TEdgeStump* current = getValue(it);
		TVertexDescriptor parentVertex = position(it);
		while(current != (TEdgeStump*) 0) {
			TVertexDescriptor childVertex = getTarget(current);
			// Create missing vertices
			if (parentVertex > childVertex) _createVertices(dest,parentVertex);
			else _createVertices(dest,childVertex);
			// Add edge
			TEdgeDescriptor e;
			if (!transpose) e = addEdge(dest, parentVertex, childVertex);
			else e = addEdge(dest, childVertex, parentVertex);
			assignCargo(e, getCargo(current));
			current = getNextT(current);
		}
		goNext(it);
	}
	dest.data_id_managerV = source.data_id_managerV;
	dest.data_root = source.data_root;
}

template<typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Tree<TCargo, TSpec> > const& source,
		   Graph<Tree<TCargo, TSpec> >& dest) 
{
	_copyGraph(source, dest, false); 
}

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
transpose(Graph<Tree<TCargo, TSpec> > const& source,
		  Graph<Tree<TCargo, TSpec> >& dest)
{
	SEQAN_CHECKPOINT
	_copyGraph(source, dest, true);
	
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
transpose(Graph<Tree<TCargo, TSpec> >& g)
{
	Graph<Tree<TCargo, TSpec> > dest;
	_copyGraph(g, dest, true);
	g = dest;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Size<Graph<Tree<TCargo, TSpec> > >::Type 
numEdges(Graph<Tree<TCargo, TSpec> > const& g) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Iterator<String<TEdgeStump*> const>::Type TIterConst;
	TSize count=0;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TEdgeStump* current = getValue(it);
		while(current!=0) {
			++count;
			current = getNextT(current);
		}
	}
	return count;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Size<Graph<Tree<TCargo, TSpec> > >::Type 
numVertices(Graph<Tree<TCargo, TSpec> > const& g) 
{
	SEQAN_CHECKPOINT
	return idCount(g.data_id_managerV);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline bool 
empty(Graph<Tree<TCargo, TSpec> > const& g) 
{
	SEQAN_CHECKPOINT
	return (!idCount(g.data_id_managerV));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
clearEdges(Graph<Tree<TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	//removeAllChildren(g, getRoot(g));
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> >::Type TIter;
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	for(TIter it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TEdgeStump* current = getValue(it);
		if(current == (TEdgeStump*) 0) continue;
		TVertexDescriptor sourceVertex = position(it);
		assignValue(g.data_parent, sourceVertex, nilVertex);
		removeOutEdges(g, sourceVertex);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
clearVertices(Graph<Tree<TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	clearEdges(g);
	releaseAll(g.data_id_managerV);
	clear(g.data_vertex);
	clear(g.data_parent);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void 
clear(Graph<Tree<TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	clearVertices(g);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<Tree<TCargo, TSpec> > >::Type 
outDegree(Graph<Tree<TCargo, TSpec> > const& g, 
		  TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Size<TGraph>::Type TSize;
	TSize count=0;
	TEdgeStump* current = getValue(g.data_vertex, vertex);
	while(current!=0) {
		current = getNextT(current);
		++count;
	}
	return count;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<Tree<TCargo, TSpec> > >::Type 
inDegree(Graph<Tree<TCargo, TSpec> > const& g, 
		 TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Iterator<String<TEdgeStump*> const>::Type TIterConst;
	TSize count=0;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TEdgeStump* current = getValue(it);
		while(current!=0) {
			if ( (TVertexDescriptor) getTarget(current)==vertex) ++count;
			current = getNextT(current);
		}
	}
	return count;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<Tree<TCargo, TSpec> > >::Type 
degree(Graph<Tree<TCargo, TSpec> > const& g,
	   TVertexDescriptor const vertex) 
{
	return (inDegree(g,vertex)+outDegree(g,vertex));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec> 
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type 
addVertex(Graph<Tree<TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT	
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	TVertexDescriptor vd;
	if (empty(g)) {
		vd = obtainId(g.data_id_managerV);
		g.data_root = vd;
	} else {
		vd = obtainId(g.data_id_managerV);
	}
	if (vd == length(g.data_vertex)) {
		appendValue(g.data_vertex, (TEdgeStump*) 0);
		fill(g.data_parent, vd + 1, nilVertex, Generous());
	} else {
		value(g.data_vertex, vd) = (TEdgeStump*) 0;
		assignValue(g.data_parent, vd, nilVertex);
	}
	return vd;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
removeVertex(Graph<Tree<TCargo, TSpec> >& g, 
			 TVertexDescriptor const v) 
{
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)

	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	TEdgeStump* current = getValue(g.data_vertex, v);
	while(current!=0) {
		assignValue(g.data_parent, childVertex(g, current), nilVertex);
		current = getNextT(current);
	}
	assignValue(g.data_parent, v, nilVertex);
	removeOutEdges(g,v); // Remove all outgoing edges
	removeInEdges(g,v); // Remove all incoming edges
	releaseId(g.data_id_managerV, v); // Release id
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename EdgeDescriptor<Graph<Tree<TCargo, TSpec> > >::Type 
addEdge(Graph<Tree<TCargo, TSpec> >& g, 
		TVertexDescriptor const parent, 
		TVertexDescriptor const child) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, parent) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, child) == true)

	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Id<TGraph>::Type TId;

	TEdgeStump* edge_ptr;
	allocate(g.data_allocator, edge_ptr, 1);
	valueConstruct(edge_ptr);
	assignTarget(edge_ptr, child);
	assignValue(g.data_parent, child, parent);
	assignNextT(edge_ptr, (TEdgeStump*) 0);
	if (getValue(g.data_vertex, parent)!=0) {
		assignNextT(edge_ptr, getValue(g.data_vertex, parent));
	}
	value(g.data_vertex, parent)=edge_ptr;
	return edge_ptr;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename EdgeDescriptor<Graph<Tree<TCargo, TSpec> > >::Type 
addEdge(Graph<Tree<TCargo, TSpec> >& g, 
		TVertexDescriptor const parent, 
		TVertexDescriptor const child,
		TCargo const cargo) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TEdgeDescriptor e = addEdge(g,parent,child);
	assignCargo(e,cargo);
	return e;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
void 
removeEdge(Graph<Tree<TCargo, TSpec> >& g,
		   TVertexDescriptor const parent,
		   TVertexDescriptor const child) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, parent) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, child) == true)
	
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;

	// Find edge and predecessor
	TEdgeStump* pred = 0;
	TEdgeStump* current = getValue(g.data_vertex, parent);
	while(current != (TEdgeStump*) 0) {
		if ( (TVertexDescriptor) getTarget(current) == child) break;
		pred = current;
		current = getNextT(current);
	}
	
	// Not found?
	if (current == (TEdgeStump*) 0) return;
	assignValue(g.data_parent, child, getNil<TVertexDescriptor>());
	
	// Relink the next pointer of predecessor
	if (pred != (TEdgeStump*) 0) assignNextT(pred, getNextT(current));
	else value(g.data_vertex, parent) = getNextT(current);
	
	// Deallocate
	valueDestruct(current);
	deallocate(g.data_allocator, current, 1);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TEdgeDescriptor>
void 
removeEdge(Graph<Tree<TCargo, TSpec> >& g,
		   TEdgeDescriptor const edge)
{
	SEQAN_CHECKPOINT
	removeEdge(g, parentVertex(g,edge), childVertex(g,edge));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline void 
removeOutEdges(Graph<Tree<TCargo, TSpec> >& g, 
			   TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)

	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	while(getValue(g.data_vertex, v) != (TEdgeStump*) 0) {
		TVertexDescriptor target = targetVertex(g,(getValue(g.data_vertex, v)));
		removeEdge(g,v,target);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline void 
removeInEdges(Graph<Tree<TCargo, TSpec> >& g,
			  TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)

	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> >::Type TIter;
	for(TIter it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		if (!idInUse(g.data_id_managerV, position(it))) continue;
		TEdgeStump* current = getValue(it);
		TVertexDescriptor const sourceVertex = position(it);
		while(current!=0) {
			if ( (TVertexDescriptor) current->data_target==v) {
				removeEdge(g, sourceVertex, v);
				current = getValue(g.data_vertex, sourceVertex);
			} else {
				current = getNextT(current);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type 
targetVertex(Graph<Tree<TCargo, TSpec> > const& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	return getTarget(edge);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type 
sourceVertex(Graph<Tree<TCargo, TSpec> > const& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> const>::Type TIterConst;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TEdgeDescriptor current = getValue(it);
		while(current!=(TEdgeDescriptor) 0) {
			if (current == edge) return position(it);
			current=getNextT(current);
		}
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TMatrix>
inline void
getAdjacencyMatrix(Graph<Tree<TCargo, TSpec> > const& g, 
				   TMatrix& mat) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Size<TMatrix>::Type TSize;
	TSize len = getIdUpperBound(g.data_id_managerV);
	setDimension(mat, 2);
	setLength(mat, 0, len);
	setLength(mat, 1, len);
	resize(mat);
	for (TSize i=0;i<len*len;++i) value(mat,i) = 0;

	typedef typename Iterator<String<TEdgeStump*> const>::Type TIterConst;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TVertexDescriptor parentV = position(it);
		TEdgeStump* current = getValue(it);
		while(current!=0) {
			TVertexDescriptor childV = getTarget(current);
			assignValue(mat,parentV*len+childV, getValue(mat,parentV*len+childV)+1);
			current=getNextT(current);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Tree<TCargo, TSpec> > >::Type 
findEdge(Graph<Tree<TCargo, TSpec> >& g,
		 TVertexDescriptor const v,
		 TVertexDescriptor const w)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, w) == true)
	
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;

	if (getValue(g.data_parent, w) == v) {
		TEdgeStump* current = getValue(g.data_vertex, v);
		while((TEdgeStump*) current != 0) {
			if (getTarget(current) == w) return current;
			current = getNextT(current);
		}
	} else if (getValue(g.data_parent, v) == w) {
		TEdgeStump* current = getValue(g.data_vertex, w);
		while((TEdgeStump*) current != 0) {
			if (getTarget(current) == v) return current;
			current = getNextT(current);
		}
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec, typename TIDString>
inline void
write(TFile & target,
	  Graph<Tree<TCargo, TSpec> > const& g,
	  TIDString const &,
	  Raw)
{
	SEQAN_CHECKPOINT
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> const>::Type TIterConst;
	_streamWrite(target,"Adjacency list:\n");
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		_streamPutInt(target, position(it));
		_streamWrite(target," -> ");
		TEdgeStump* current = getValue(it);
		while(current!=0) {
			_streamPutInt(target, getTarget(current));
			_streamPut(target, ',');
			current=getNextT(current);
		}
		_streamPut(target, '\n');
	}
	_streamWrite(target,"Edge list:\n");
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		TEdgeStump* current = getValue(it);
		while(current!=0) {
			_streamWrite(target,"Source: ");
			_streamPutInt(target, position(it));		
			_streamPut(target, ',');
			_streamWrite(target,"Target: ");
			_streamPutInt(target, getTarget(current));
			_streamPut(target, ' ');
			_streamWrite(target,"(Id: ");
			_streamPutInt(target, _getId(current));
			_streamPut(target, ',');
			_streamWrite(target," Cargo-Type: ");
			_streamWrite(target, typeid(getCargo(current)).name());
			_streamPut(target, ')');
			_streamPut(target, '\n');
			current=getNextT(current);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
assignRoot(Graph<Tree<TCargo, TSpec> >& g,
		   TVertexDescriptor const vertex)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	
	g.data_root = vertex;
	assignValue(g.data_parent, vertex, getNil<TVertexDescriptor>());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
createRoot(Graph<Tree<TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TVertexDescriptor r = addVertex(g);
	g.data_root = r;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type&
root(Graph<Tree<TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	return g.data_root;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type 
getRoot(Graph<Tree<TCargo, TSpec> > const& g)
{
	SEQAN_CHECKPOINT
	return g.data_root;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline bool
isRoot(Graph<Tree<TCargo, TSpec> > const& g,
	   TVertexDescriptor v)
{
	SEQAN_CHECKPOINT
	return ( (TVertexDescriptor) g.data_root == v);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.isLeaf:
..cat:Spec.Tree
..summary:Tests whether a given vertex is a leaf or not.
..signature:isLeaf(g, v)
..param.g:A tree.
...type:Spec.Tree
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:True if vertex is a leaf.
*/

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline bool
isLeaf(Graph<Tree<TCargo, TSpec> > const& g,
	   TVertexDescriptor v)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)
	
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	return (value(g.data_vertex, v) ==  (TEdgeStump*) 0);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.numTreeEdges:
..cat:Spec.Tree
..summary:Number of tree edges.
..signature:numTreeEdges(g)
..param.g:A tree.
...type:Spec.Tree
..returns:Number of tree edges. Faster than numEdges for trees.
*/

template<typename TCargo, typename TSpec>
inline typename Size<Graph<Tree<TCargo, TSpec> > >::Type 
numTreeEdges(Graph<Tree<TCargo, TSpec> > const& g) 
{
	SEQAN_CHECKPOINT
	if (empty(g)) return 0;
	else return numVertices(g) - 1;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.numChildren:
..cat:Spec.Tree
..summary:Number of children of a given tree vertex.
..signature:numChildren(g, v)
..param.g:A tree.
...type:Spec.Tree
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:Number of children
*/

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<Tree<TCargo, TSpec> > >::Type 
numChildren(Graph<Tree<TCargo, TSpec> > const& g,
			TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	
	return outDegree(g, vertex);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.addChild:
..cat:Spec.Tree
..summary:Adds a new child vertex to a parent vertex.
Optionally a cargo can be attached to the parent-child edge.
..signature:addChild(g, parent [, cargo])
..param.g:A tree.
...type:Spec.Tree
..param.parent:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..param.cargo:A cargo object.
...type:Metafunction.Cargo
..returns:A new vertex descriptor.
...type:Metafunction.VertexDescriptor
..see:Function.removeChild
..see:Function.removeAllChildren
*/

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type 
addChild(Graph<Tree<TCargo, TSpec> >& g,
		 TVertexDescriptor parent) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, parent) == true)
	TVertexDescriptor child = addVertex(g);
	addEdge(g,parent,child);
	return child;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type 
addChild(Graph<Tree<TCargo, TSpec> >& g,
		 TVertexDescriptor const parent,
		 TCargo const cargo) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, parent) == true)
	TVertexDescriptor child = addVertex(g);
	addEdge(g,parent,child,cargo);
	return child;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.removeChild:
..cat:Spec.Tree
..summary:Removes a child from the tree given a parent.
..signature:removeChild(g, parent, child)
..param.g:A tree.
...type:Spec.Tree
..param.parent:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..param.child:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:void
..see:Function.addChild
..see:Function.removeAllChildren
*/

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline void 
removeChild(Graph<Tree<TCargo, TSpec> >& g,
			TVertexDescriptor const parent,
			TVertexDescriptor const child)
{
	SEQAN_CHECKPOINT
	if (!isLeaf(g,child)) removeAllChildren(g,child);
	removeEdge(g,parent, child);  // Parent map is cleared in removeEdge
	releaseId(g.data_id_managerV, child); // Release id
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.removeAllChildren:
..cat:Spec.Tree
..summary:Removes all children from the tree given a parent.
..signature:removeChild(g, parent)
..param.g:A tree.
...type:Spec.Tree
..param.parent:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:void
..see:Function.addChild
..see:Function.removeChild
*/

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline void 
removeAllChildren(Graph<Tree<TCargo, TSpec> >& g, 
				  TVertexDescriptor const parent) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, parent) == true)
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;

	while(getValue(g.data_vertex, parent) != (TEdgeStump*) 0) {
		TVertexDescriptor child = childVertex(g,(getValue(g.data_vertex, parent)));
		if (!isLeaf(g,child)) removeAllChildren(g,child);
		removeEdge(g,parent, child);  // Parent map is cleared in removeEdge
		releaseId(g.data_id_managerV, child); // Release id
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.childVertex:
..cat:Spec.Tree
..summary:Returns the child vertex of an edge.
..signature:childVertex(g, e)
..param.g:A tree.
...type:Spec.Tree
..param.e:An edge descriptor.
...type:Metafunction.EdgeDescriptor
..returns:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..see:Function.parentVertex
*/
template<typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type 
childVertex(Graph<Tree<TCargo, TSpec> > const& g,
			TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	return getTarget(edge);
}

/**
.Function.parentVertex:
..cat:Spec.Tree
..summary:Returns the parent vertex of an edge.
..signature:parentVertex(g, e)
..param.g:A tree.
...type:Spec.Tree
..param.e:An edge descriptor.
...type:Metafunction.EdgeDescriptor
..returns:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..see:Function.parentVertex
*/
template<typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type 
parentVertex(Graph<Tree<TCargo, TSpec> > const& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	return getValue(g.data_parent, getTarget(edge));
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
