#ifndef SEQAN_HEADER_GRAPH_INTERFACE_H
#define SEQAN_HEADER_GRAPH_INTERFACE_H

namespace SEQAN_NAMESPACE_MAIN
{

// Default directed graph
template<typename TCargo = void, typename TSpec = Default>
struct Directed;

// Default undirected graph
template<typename TCargo = void, typename TSpec = Default>
struct Undirected;

// Default Tree
template<typename TCargo = void, typename TSpec = Default>
struct Tree;

// Default Automaton
template<typename TAlphabet = char, typename TCargo = void, typename TSpec = Default>
struct Automaton;


//////////////////////////////////////////////////////////////////////////////
// Graph
//////////////////////////////////////////////////////////////////////////////

/**
.Class.Graph:
..cat:Graph
..summary:Generic graph.
..signature:Graph<TEdges, TSpec>
..param.TEdges:The storage type for edges, e.g., edge list (adjacency list) or automaton (transition table).
...remarks:The default corresponds to a directed graph.
The @Metafunction.EdgeType@ can be used to retrieve the type of the edge stump.
...default:Directed<>
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Default$, see @Tag.Default@.
..include:graph.h
*/
template<typename TEdges = Directed<>, typename TSpec = Default>
class Graph;


//////////////////////////////////////////////////////////////////////////////
// General Graph Metafunction
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Spec.param.T.type:Class.Graph

template<typename TEdges, typename TSpec>
struct Spec<Graph<TEdges, TSpec> > 
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TEdges, typename TSpec>
struct Spec<Graph<TEdges, TSpec> const>
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Size.param.T.type:Class.Graph

template<typename TEdges, typename TSpec>
struct Size<Graph<TEdges, TSpec> > 
{
	typedef size_t Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TEdges, typename TSpec>
struct Size<Graph<TEdges, TSpec> const>
{
	typedef size_t Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.EdgeDescriptor.param.T.type:Class.Graph

template<typename TEdges, typename TSpec>
struct EdgeDescriptor<Graph<TEdges, TSpec> > 
{
	typedef typename EdgeType<Graph<TEdges, TSpec> >::Type* Type;
};

template<typename TEdges, typename TSpec>
struct EdgeDescriptor<Graph<TEdges, TSpec> const>
{
	typedef typename EdgeType<Graph<TEdges, TSpec> const>::Type* Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.VertexDescriptor.param.T.type:Class.Graph

template<typename TEdges, typename TSpec>
struct VertexDescriptor<Graph<TEdges, TSpec> > 
{
	typedef typename Id<Graph<TEdges, TSpec> >::Type Type;
};

template<typename TEdges, typename TSpec>
struct VertexDescriptor<Graph<TEdges, TSpec> const>
{
	typedef typename Id<Graph<TEdges, TSpec> >::Type Type;
};


//////////////////////////////////////////////////////////////////////////////

///.Metafunction.EdgeType.param.T.type:Class.Graph

template<typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<Directed<TCargo, TEdgeSpec>, TSpec> > {
	typedef EdgeStump<TCargo, true, false, true, TEdgeSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<Directed<TCargo, TEdgeSpec>, TSpec> const> {
	typedef EdgeStump<TCargo, true, false, true, TEdgeSpec> const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
struct EdgeType<Graph<Directed<TCargo, WithoutEdgeId>, TSpec> > {
	typedef EdgeStump<TCargo, true, false, false, TSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
struct EdgeType<Graph<Directed<TCargo, WithoutEdgeId>, TSpec> const> {
	typedef EdgeStump<TCargo, true, false, false, TSpec> const Type;
};


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<Tree<TCargo, TEdgeSpec>, TSpec> > {
	typedef EdgeStump<TCargo, true, true, false, TEdgeSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<Tree<TCargo, TEdgeSpec>, TSpec> const> {
	typedef EdgeStump<TCargo, true, true, false, TEdgeSpec> const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<Undirected<TCargo, TEdgeSpec>, TSpec> > {
	typedef EdgeStump<TCargo, true, true, true, TEdgeSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<Undirected<TCargo, TEdgeSpec>, TSpec> const> {
	typedef EdgeStump<TCargo, true, true, true, TEdgeSpec> const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
struct EdgeType<Graph<Undirected<TCargo, WithoutEdgeId>, TSpec> > {
	typedef EdgeStump<TCargo, true, true, false, TSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
struct EdgeType<Graph<Undirected<TCargo, WithoutEdgeId>, TSpec> const> {
	typedef EdgeStump<TCargo, true, true, false, TSpec> const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> > {
	typedef EdgeStump<TCargo, false, false, true, TEdgeSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const> {
	typedef EdgeStump<TCargo, false, false, true, TSpec> const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
struct EdgeType<Graph<Automaton<TAlphabet, TCargo, WithoutEdgeId>, TSpec> > {
	typedef EdgeStump<TCargo, false, false, false, TSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
struct EdgeType<Graph<Automaton<TAlphabet, TCargo, WithoutEdgeId>, TSpec> const> {
	typedef EdgeStump<TCargo, false, false, false, TSpec> const Type;
};


//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Cargo.param.T.type:Class.Graph

template<typename TEdges, typename TSpec>
struct Cargo<Graph<TEdges, TSpec> > {
	typedef typename Cargo<typename EdgeType<Graph<TEdges, TSpec> >::Type>::Type TType;
};


template<typename TEdges, typename TSpec>
struct Cargo<Graph<TEdges, TSpec> const> {
	typedef typename Cargo<typename EdgeType<Graph<TEdges, TSpec> const>::Type>::Type TType;
};



//////////////////////////////////////////////////////////////////////////////

///.Metafunction.EdgeType.param.T.type:Class.Graph

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
struct Alphabet<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> > {
	typedef TAlphabet Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
struct Alphabet<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const> {
	typedef TAlphabet Type;
};


//////////////////////////////////////////////////////////////////////////////
// Generic Graph Functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getInfinity:
..cat:Graph
..summary:Utility function returning a value that represents infinity.
Useful for various graph algorithms, e.g., Dijkstra.
..signature:getInfinity<T>()
..returns:Pseudo infinity value for type T.
..see:Function.getNil
*/
template <typename T>
inline T const
getInfinity()
{
SEQAN_CHECKPOINT
	T * _tag = 0;
	return supremumValueImpl(_tag);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getNil:
..cat:Graph
..summary:Utility function returning a value that represents nil.
Useful for various graph algorithms, e.g., missing predecessors, vertices that have not been visited, etc.
..signature:getNil<T>()
..returns:Pseudo nil value for type T.
..see:Function.getInfinity
*/
template <typename T>
inline T const
getNil(T *)
{
	return ~0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline T const
getNil()
{
SEQAN_CHECKPOINT
	T * _tag = 0;
	return getNil(_tag);
}

//////////////////////////////////////////////////////////////////////////////

// Simple _getId function to get the id for a vertex descriptor which is the id!
template<typename TId>
inline TId
_getId(TId const id)
{
	SEQAN_CHECKPOINT
	return id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TEdges, typename TSpec, typename TVertexDescriptor>
inline void
_createVertices(Graph<TEdges, TSpec>& g,
				TVertexDescriptor const maxId) 
{
		// Create missing vertices
		while (maxId >= getIdUpperBound(g.data_id_managerV))
		{
			addVertex(g);
		}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.addEdges:
..cat:Graph
..summary:Shortcut to add multiple edges at once.
Creates vertices implicitly.
..signature:addEdge(g, edges, size)
..param.g:A graph.
...type:Class.Graph
..param.edges:An array of vertex descriptors. It is assumed that the
edges are stored in the following way: Source1, Target1, Source2, Target2, Source3, ...
...type:Metafunction.VertexDescriptor
..param.size:Size of the array. Must be a multiple of 2.
...type:Metafunction.Size
..returns:void
..see:Function.addEdge
*/
template<typename TEdges, typename TSpec, typename TEdgeArray, typename TSize>
inline void
addEdges(Graph<TEdges, TSpec>& dest,
		 TEdgeArray const edges,
		 TSize const size) 
{
	typedef typename VertexDescriptor<Graph<TEdges, TSpec> >::Type TVertexDescriptor;
	for(TSize i=0;i<size;++i) {
		TVertexDescriptor source = edges[2*i];
		TVertexDescriptor target = edges[2*i+1];
		// Create missing vertices
		if (source>target) _createVertices(dest,source);
		else _createVertices(dest,target);
		// Add edge
		SEQAN_ASSERT(idInUse(dest.data_id_managerV, source) == true)
		SEQAN_ASSERT(idInUse(dest.data_id_managerV, target) == true)
		addEdge(dest, source, target);
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename TEdges, typename TSpec>
inline TStream &
operator << (TStream & target, 
			 Graph<TEdges, TSpec> const& source)
{
SEQAN_CHECKPOINT
	write(target, source);
	return target;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
