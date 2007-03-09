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
	struct VertexIterator;

	// Edge iterator
	template <typename TSpec = Default>
	struct EdgeIterator;

	// OutEdge iterator
	template <typename TSpec = Default>
	struct OutEdgeIterator;

	// Adjacency iterator
	template <typename TSpec = Default>
	struct AdjacencyIterator;

	// Bfs iterator
	template <typename TSpec = Default>
	struct BfsIterator;

	// Dfs iterator
	template <typename TSpec = Default>
	struct DfsIterator;

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
