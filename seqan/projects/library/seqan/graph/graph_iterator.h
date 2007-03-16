#ifndef SEQAN_HEADER_GRAPH_ITERATOR_H
#define SEQAN_HEADER_GRAPH_ITERATOR_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph Iterators
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec>
struct GraphIterator;

	// Vertex iterator
	template <typename TSpec = Default>
	struct InternalVertexIterator;

	// Edge iterator
	template <typename TSpec = Default>
	struct InternalEdgeIterator;

	// OutEdge iterator
	template <typename TSpec = Default>
	struct InternalOutEdgeIterator;

	// Adjacency iterator
	template <typename TSpec = Default>
	struct InternalAdjacencyIterator;

	// Bfs iterator
	template <typename TSpec = Default>
	struct InternalBfsIterator;

	// Dfs iterator
	template <typename TSpec = Default>
	struct InternalDfsIterator;

//////////////////////////////////////////////////////////////////////////////
// Graph Iterators - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Host.param.T.type:Class.Graph

template<typename TGraph, typename TIteratorSpec>
struct Host<Iter<TGraph, GraphIterator<TIteratorSpec> > >
{	
	typedef TGraph Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct Host<Iter<TGraph const, GraphIterator<TIteratorSpec> > >
{	
	typedef TGraph const Type;
};


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
