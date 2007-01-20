#ifndef SEQAN_HEADER_GRAPH_INTERFACE_H
#define SEQAN_HEADER_GRAPH_INTERFACE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Generic Graph Functions
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
	
template<typename TEdgeArray, typename TSize, typename TEdges, typename TSpec>
inline void
_copy(TEdgeArray const edges, 
	  TSize const size, 
	  Graph<TEdges, TSpec>& dest) 
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

template<typename TEdges, typename TSpec>
inline typename Size<Graph<TEdges, TSpec> >::Type 
numVertices(Graph<TEdges, TSpec> const& g) 
{
	SEQAN_CHECKPOINT
	return idCount(g.data_id_managerV);
}


template<typename TEdges, typename TSpec>
inline typename Size<Graph<TEdges, TSpec> >::Type 
numEdges(Graph<TEdges, TSpec> const& g) 
{
	SEQAN_CHECKPOINT
	return idCount(g.data_id_managerE);
}

template<typename TEdges, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<TEdges, TSpec> >::Type 
degree(Graph<TEdges, TSpec> const& g, 
	   TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	return (inDegree(g,vertex)+outDegree(g,vertex));
}

template<typename TEdges, typename TSpec>
inline bool 
empty(Graph<TEdges, TSpec> const& g) 
{
	SEQAN_CHECKPOINT
	return (!(idCount(g.data_id_managerV)));
}


template<typename TEdges, typename TSpec>
inline void 
clear(Graph<TEdges, TSpec>& g) 
{
	SEQAN_CHECKPOINT
	clearVertices(g);
}

template<typename TEdges, typename TSpec>
inline void 
transpose(Graph<TEdges, TSpec> const& source,
		  Graph<TEdges, TSpec>& dest)
{
	SEQAN_CHECKPOINT
	_copy(source, dest, true);
}

template<typename TEdges, typename TSpec, typename TVertexDescriptor>
inline void 
removeVertex(Graph<TEdges, TSpec>& g, 
			 TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)

	removeOutEdges(g,v); // Remove all outgoing edges
	removeInEdges(g,v); // Remove all incoming edges
	releaseId(g.data_id_managerV, v); // Release id
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

//////////////////////////////////////////////////////////////////////////////
// Graph Functions relying on graph data structure
//////////////////////////////////////////////////////////////////////////////
// _copy(source,dest,transpose) 
// removeOutEdges(graph, vertex);
// removeInEdges(graph, vertex);
// removeVertex(graph, vertex);
// clearVertices(graph)
// clearEdges(graph)
// outDegree(graph, vertex)
// inDegree(graph, vertex)
// addEdge(graph,source,target) 
// addEdge(graph,source,target,cargo)
// removeEdge(graph,source,target) 
// targetVertex(graph, edge);
// sourceVertex(graph, edge);
// getAdjacencyMatrix(graph, mat);
// getPredecessorVertex(graph, vertex, edegMap, char)
// getPredecessorVertex(graph, vertex, char)
// getSuccessorVertex(graph, vertex, edegMap, char)
// getSuccessorVertex(graph, vertex, char)



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
