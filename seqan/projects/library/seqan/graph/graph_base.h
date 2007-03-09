#ifndef SEQAN_HEADER_GRAPH_BASE_H
#define SEQAN_HEADER_GRAPH_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// General Graph Metafunction
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

	
/**
.Metafunction.Id:
..summary:Type of an object that represents an id.
..signature:Id<T>::Type
..param.T:Type for which a suitable id type is determined.
..returns.param.Type:Id type.
..remarks.text:The id type of a container is the type that is used to uniquely identify its elements.
In most cases this type is unsigned int.
..example.code:Id<Graph<> >::Type id; //id has type unsigned int
*/
template<typename T>
struct Id {
	typedef unsigned int Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename T>
struct Id<T const> {
	typedef unsigned int Type;
};

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.VertexDescriptor:
..summary:Type of an object that represents a vertex descriptor.
..signature:VertexDescriptor<T>::Type
..param.T:Type T must be a graph. All graphs currently use ids as vertex descriptors.
..returns.param.Type:VertexDescriptor type.
..remarks.text:The vertex descriptor is a unique handle to a vertex in a graph.
It is used in various graph functions, e.g., to add edges, to create OutEdge Iterators or to remove a vertex.
It is also used to attach properties to vertices.
..example.code:VertexDescriptor<Graph<> >::Type vD; //vD is a vertex descriptor
*/
template<typename T>
struct VertexDescriptor;



//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.EdgeDescriptor:
..summary:Type of an object that represents an edge descriptor.
..signature:EdgeDescriptor<T>::Type
..param.T:Type T must be a graph. All graphs currently use pointer to edge stumps as edge descriptors.
..returns.param.Type:EdgeDescriptor type.
..remarks.text:The edge descriptor is a unique handle to a given edge in a graph.
It is used in various graph functions, e.g., to remove edges, to assign a cargo to an edge or to get the endpoints of an edge.
It is also used to attach properties to edges.
..example.code:EdgeDescriptor<Graph<> >::Type eD; //eD is an edge descriptor
*/
template<typename T>
struct EdgeDescriptor;

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Cargo:
..summary:Type of the cargo of an edge.
..signature:Cargo<T>::Type
..param.T:Edge type for which the cargo type is determined.
..returns.param.Type:Cargo type.
..remarks.text:The cargo type of an edge indicates the kind of information that is stored with the edge.
..example.code:Cargo<EdgeStump<int, TSpec> >::Type c; //c has type int
*/
template<typename T>
struct Cargo;


//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.EdgeType:
..summary:Edge type of a graph object.
..signature:EdgeType<T>::Type
..param.T:Type T must be a graph.
..returns.param.Type:Edge type.
..remarks.text:The specific edge stump type that is used in a graph.
..example.code:EdgeType<TGraph>::Type e; //e is an edge in TGraph
*/
template<typename T>
struct EdgeType;

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Alphabet:
..summary:Access to the Alphabet type.
..signature:Alphabet<T>::Type
..param.T:Type T must be a type that uses some kind of alphabet internally.
..returns.param.Type:Alphabet type.
..remarks.text:Type T can be for example an automaton where the alphabet type describes the domain of the transition labels.
..example.code:Alphabet<Graph<Automaton<Dna> > >::Type alph; //alph is of type Dna
*/
template<typename T>
struct Alphabet;


//////////////////////////////////////////////////////////////////////////////
// General Graph Tags
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.WithoutEdgeId
..summary:Indicates whether an edge stump has an edge id or not.
..value.WithoutEdgeId:No edge id is stored in the edge stump. 
*/
struct WithoutEdgeId_;
typedef Tag<WithoutEdgeId_> const WithoutEdgeId;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.DotDrawing
..summary:Switch to trigger drawing in dot format.
..value.DotDrawing:Graphs in dot format.
*/

struct TagDotDrawing_;
typedef Tag<TagDotDrawing_> const DotDrawing;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.EmptyMap
..summary:Switch to indicate whether a property map is present or not.
..value.EmptyMap:The property map is not present. 
*/
struct EmptyMap_;
typedef Tag<EmptyMap_> const EmptyMap;


//////////////////////////////////////////////////////////////////////////////
// Graph - Default edge stumps
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

// Default Edge Stump: No cargo but edge id
template<typename TCargo = void, typename TSpec = Default>
class EdgeStump;

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.VertexDescriptor.param.T.type:Class.EdgeStump

template<typename TCargo, typename TSpec>
struct VertexDescriptor<EdgeStump<TCargo, TSpec> > 
{
	typedef typename Id<EdgeStump<TCargo, TSpec> >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
struct VertexDescriptor<EdgeStump<TCargo, TSpec> const> 
{
	typedef typename Id<EdgeStump<TCargo, TSpec> >::Type Type;
};


//////////////////////////////////////////////////////////////////////////////

// Default Edge Stump Undirected: No cargo but edge id
template<typename TCargo = void, typename TSpec = Default>
class EdgeStumpU;


//////////////////////////////////////////////////////////////////////////////

///.Metafunction.VertexDescriptor.param.T.type:Class.EdgeStumpU

template<typename TCargo, typename TSpec>
struct VertexDescriptor<EdgeStumpU<TCargo, TSpec> > 
{
	typedef typename Id<EdgeStumpU<TCargo, TSpec> >::Type Type;
};

template<typename TCargo, typename TSpec>
struct VertexDescriptor<EdgeStumpU<TCargo, TSpec> const> 
{
	typedef typename Id<EdgeStumpU<TCargo, TSpec> >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

// Default Edge Stump Tree: No cargo but edge id
template<typename TCargo = void, typename TSpec = Default>
class EdgeStumpT;

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.VertexDescriptor.param.T.type:Class.EdgeStumpT

template<typename TCargo, typename TSpec>
struct VertexDescriptor<EdgeStumpT<TCargo, TSpec> > 
{
	typedef typename Id<EdgeStumpT<TCargo, TSpec> >::Type Type;
};

template<typename TCargo, typename TSpec>
struct VertexDescriptor<EdgeStumpT<TCargo, TSpec> const> 
{
	typedef typename Id<EdgeStumpT<TCargo, TSpec> >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

// Default Automaton Edge: No cargo but edge id
template<typename TCargo = void, typename TSpec = Default>
class EdgeStumpA;

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.VertexDescriptor.param.T.type:Class.EdgeStumpA

template<typename TCargo, typename TSpec>
struct VertexDescriptor<EdgeStumpA<TCargo, TSpec> > 
{
	typedef typename Id<EdgeStumpA<TCargo, TSpec> >::Type Type;
};

template<typename TCargo, typename TSpec>
struct VertexDescriptor<EdgeStumpA<TCargo, TSpec> const> 
{
	typedef typename Id<EdgeStumpA<TCargo, TSpec> >::Type Type;
};


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
