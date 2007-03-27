#ifndef SEQAN_HEADER_GRAPH_BASE_H
#define SEQAN_HEADER_GRAPH_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// General Graph Metafunction
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.EdgeDescriptor:
..summary:Type of an object that represents an edge descriptor.
..signature:EdgeDescriptor<T>::Type
..param.T:Type T must be a graph. All graphs use a pointer to an edge stump as an edge descriptor.
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
..example.code:Cargo<Graph<Directed<int> > >::Type c; //c has type int
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
..example.code:Graph<Directed<void, WithoutEdgeId> > g; //g stores no edge ids
*/
struct WithoutEdgeId_;
typedef Tag<WithoutEdgeId_> const WithoutEdgeId;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.TreeTag
..summary:Indicates whether an edge stump belongs to a tree or not.
..remarks:Only for internal use.
..value.TreeTag:EdgeStump belongs to a tree.
*/
struct TreeTag_;
typedef Tag<TreeTag_> const TreeTag;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.DotDrawing
..summary:Switch to trigger drawing in dot format.
..value.DotDrawing:Graphs in dot format.
*/

struct TagDotDrawing_;
typedef Tag<TagDotDrawing_> const DotDrawing;



//////////////////////////////////////////////////////////////////////////////
// Graph Iterator Tags
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.VertexIterator
..summary:Tag for a VertexIterator
..value.VertexIterator:Specifies the kind of graph iterator, i.e., VertexIterator.
*/
struct VertexIterator_;
typedef Tag<VertexIterator_> const VertexIterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.EdgeIterator
..summary:Tag for an EdgeIterator
..value.EdgeIterator:Specifies the kind of graph iterator, i.e., EdgeIterator.
*/
struct EdgeIterator_;
typedef Tag<EdgeIterator_> const EdgeIterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.OutEdgeIterator
..summary:Tag for an OutEdgeIterator
..value.OutEdgeIterator:Specifies the kind of graph iterator, i.e., OutEdgeIterator.
*/
struct OutEdgeIterator_;
typedef Tag<OutEdgeIterator_> const OutEdgeIterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.AdjacencyIterator
..summary:Tag for an AdjacencyIterator
..value.AdjacencyIterator:Specifies the kind of graph iterator, i.e., AdjacencyIterator.
*/
struct AdjacencyIterator_;
typedef Tag<AdjacencyIterator_> const AdjacencyIterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.BfsIterator
..summary:Tag for a breath-first search iterator.
..value.BfsIterator:Specifies the kind of graph iterator, i.e., BfsIterator.
*/
struct BfsIterator_;
typedef Tag<BfsIterator_> const BfsIterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.DfsPreorder
..summary:Tag for a depth-first postorder search iterator.
..value.DfsPreorder:Specifies the kind of graph iterator, i.e., DfsPreorder.
*/
struct DfsPreorder_;
typedef Tag<DfsPreorder_> const DfsPreorder;

//////////////////////////////////////////////////////////////////////////////
// Graph - Default edge stump
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo = void, bool TList = true, bool TSource = false, bool TId = true, typename TSpec = Default>
class EdgeStump;

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.VertexDescriptor.param.T.type:Class.EdgeStump

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
struct VertexDescriptor<EdgeStump<TCargo, TList, TSource, TId, TSpec> > 
{
	typedef typename Id<EdgeStump<TCargo, TList, TSource, TId, TSpec> >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
struct VertexDescriptor<EdgeStump<TCargo, TList, TSource, TId, TSpec> const> 
{
	typedef typename Id<EdgeStump<TCargo, TList, TSource, TId, TSpec> >::Type Type;
};

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
