#ifndef SEQAN_HEADER_GRAPH_INTERFACE_H
#define SEQAN_HEADER_GRAPH_INTERFACE_H

namespace SEQAN_NAMESPACE_MAIN
{

// Default EdgeList
template<typename TCargo = void, typename TSpec = Default>
struct EdgeList;

// Default EdgeListU
template<typename TCargo = void, typename TSpec = Default>
struct EdgeListU;

// Default Tree
template<typename TCargo = void, typename TSpec = Default>
struct EdgeListT;

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
...default:EdgeList<>
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Default$, see @Tag.Default@.
..include:graph.h
*/
template<typename TEdges = EdgeList<>, typename TSpec = Default>
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

///.Metafunction.EdgeType.param.T.type:Class.Graph

template<typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> > {
	typedef EdgeStump<TCargo, TEdgeSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> const> {
	typedef EdgeStump<TCargo, TEdgeSpec> const Type;
};


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> > {
	typedef EdgeStumpT<TCargo, TEdgeSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> const> {
	typedef EdgeStumpT<TCargo, TEdgeSpec> const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> > {
	typedef EdgeStumpU<TCargo, TEdgeSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> const> {
	typedef EdgeStumpU<TCargo, TEdgeSpec> const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> > {
	typedef EdgeStumpA<TCargo, TEdgeSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const> {
	typedef EdgeStumpA<TCargo, TEdgeSpec> const Type;
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
// Utility functions
//////////////////////////////////////////////////////////////////////////////
template <typename T>
inline T const
_get_nil(T *)
{
	return ~0;
}

template <typename T>
inline T const
_get_nil()
{
SEQAN_CHECKPOINT
	T * _tag = 0;
	return _get_nil(_tag);
}

template <typename T>
inline T const
_get_infinity()
{
SEQAN_CHECKPOINT
	T * _tag = 0;
	return supremumValueImpl(_tag);
}

// Simple _getId function to get the id for a vertex descriptor which is the id!
template<typename TId>
inline TId
_getId(TId const id)
{
	SEQAN_CHECKPOINT
	return id;
}



//////////////////////////////////////////////////////////////////////////////
// Generic Graph Functions
//////////////////////////////////////////////////////////////////////////////
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

template<typename TEdges, typename TSpec, typename TEdgeArray, typename TSize, typename TPropertyMap, typename TProperties>
inline void
addEdges(Graph<TEdges, TSpec>& dest,
		 TEdgeArray const edges,
		 TSize const size,
		 TPropertyMap& pm,
		 TProperties const& prop)
{
	SEQAN_CHECKPOINT
	typedef Graph<TEdges, TSpec> TGraph;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	for(TSize i=0;i<size;++i) {
		TVertexDescriptor source = edges[2*i];
		TVertexDescriptor target = edges[2*i+1];
		// Create missing vertices
		if (source>target) _createVertices(dest,source);
		else _createVertices(dest,target);
		// Add edge
		SEQAN_ASSERT(idInUse(dest.data_id_managerV, source) == true)
		SEQAN_ASSERT(idInUse(dest.data_id_managerV, target) == true)
		TEdgeDescriptor e = addEdge(dest, source, target);
		initEdgeMap(dest, pm);
		assignProperty(pm,e,prop[i]);
	}
}

template<typename TEdges, typename TSpec, typename TEdgeArray, typename TSize, typename TPropertyMap1, typename TProperties1, typename TPropertyMap2, typename TProperties2>
inline void
addEdges(Graph<TEdges, TSpec>& dest,
		 TEdgeArray const edges,
		 TSize const size,
		 TPropertyMap1& pm1,
		 TProperties1 const& prop1,
		 TPropertyMap2& pm2,
		 TProperties2 const& prop2)
{
	SEQAN_CHECKPOINT
	typedef Graph<TEdges, TSpec> TGraph;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	for(TSize i=0;i<size;++i) {
		TVertexDescriptor source = edges[2*i];
		TVertexDescriptor target = edges[2*i+1];
		// Create missing vertices
		if (source>target) _createVertices(dest,source);
		else _createVertices(dest,target);
		// Add edge
		SEQAN_ASSERT(idInUse(dest.data_id_managerV, source) == true)
		SEQAN_ASSERT(idInUse(dest.data_id_managerV, target) == true)
		TEdgeDescriptor e = addEdge(dest, source, target);
		initEdgeMap(dest, pm1);
		assignProperty(pm1,e,prop1[i]);
		initEdgeMap(dest, pm2);
		assignProperty(pm2,e,prop2[i]);
	}
}

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
