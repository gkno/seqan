/// This code example illustrates the simple path growing heuristic for matching
#include <seqan/graph.h>
#include <iostream>

using namespace seqan;

int main() {
	typedef Graph<Undirected<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef Size<TGraph>::Type TSize;
/// Graph creation: 14 undirected edges {0,7}, {0,5}, ...
	TSize numEdges = 14;
	TVertexDescriptor edges[] = {0,7, 0,5, 0,8, 1,5, 1,6, 2,6, 2,5, 2,8, 3,8, 3,4, 4,5, 4,6, 4,7, 4,8};
	TGraph g;
	addEdges(g,edges, numEdges);
	std::cout << g << std::endl;
/// One external property map: Edge weights
	unsigned int weights[] =    {20,  5,   19,  6,   7,   10,  3,   9,   12,  11,  13,  12,  9,   12};
	String<unsigned int> weightMap;
	resizeEdgeMap(g, weightMap, weights);
/// Out-parameters: Selected edges
	String<bool> edgeMap;	
/// Path growing algorithm
	unsigned int weight = path_growing_algorithm(g, weightMap, edgeMap);
/// Console output
	std::cout << "Found matching of weight: " << weight << std::endl;
	std::cout << "Selected edges are: " << std::endl;
	TEdgeIterator it(g);
	for(;!atEnd(it);++it) {
		if (getProperty(edgeMap, *it) == true) {
			std::cout << '{' << sourceVertex(it) << ',' << targetVertex(it) << '}' << " - Weight: " << getProperty(weightMap, *it) << std::endl;
		}
	}
	return 0;
}
