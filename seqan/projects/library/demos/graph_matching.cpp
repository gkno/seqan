/// This code example illustrates how to perform a weighted bipartite matching
#include <seqan/graph.h>
#include <iostream>

using namespace seqan;

int main() {
	typedef Graph<Undirected<void> > TGraph;
	typedef Size<TGraph>::Type TSize;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;

/// Create a simple bipartite graph
	TGraph matchGraph;
	TVertexDescriptor v0 = addVertex(matchGraph);
	TVertexDescriptor v1 = addVertex(matchGraph);
	TVertexDescriptor v2 = addVertex(matchGraph);
	TVertexDescriptor v3 = addVertex(matchGraph);
	TVertexDescriptor v4 = addVertex(matchGraph);
	TVertexDescriptor v5 = addVertex(matchGraph);
	TEdgeDescriptor e0 = addEdge(matchGraph, v0, v3);
	TEdgeDescriptor e1 = addEdge(matchGraph, v0, v5);
	TEdgeDescriptor e2 = addEdge(matchGraph, v1, v3);
	TEdgeDescriptor e3 = addEdge(matchGraph, v1, v4);
	TEdgeDescriptor e4 = addEdge(matchGraph, v2, v4);
	TEdgeDescriptor e5 = addEdge(matchGraph, v2, v5);

/// Differentiate the two vertex sets with true and false
	String<bool> vertexMap;
	fill(vertexMap, numVertices(matchGraph), false);
	value(vertexMap, v0) = true;
	value(vertexMap, v1) = true;
	value(vertexMap, v2) = true;

/// Add the weights for the edges
	String<TSize> weights;
	resize(weights, numEdges(matchGraph));
	property(weights, e0) = 1;
	property(weights, e1) = 4;
	property(weights, e2) = 6;
	property(weights, e3) = 8;
	property(weights, e4) = 6;
	property(weights, e5) = 1;

/// Compute a weighted bipartite matching
	typedef String<Pair<TVertexDescriptor, TVertexDescriptor> > TEdges;
	TEdges edges;
	TSize val =  weighted_bipartite_matching(matchGraph, vertexMap, weights, edges);

/// Output the weight and the selected edges
	std::cout << "Weight: " << val << std::endl;
	typedef Iterator<TEdges>::Type TEdgeIter;
	for(TEdgeIter itE = begin(edges); itE != end(edges); ++itE) 
		std::cout << (value(itE)).i1 << ',' << (value(itE)).i2 << ": Weight = " << property(weights, findEdge(matchGraph, (value(itE)).i1, (value(itE)).i2)) << std::endl;

	return 0;
}

