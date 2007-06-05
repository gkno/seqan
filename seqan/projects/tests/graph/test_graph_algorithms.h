#ifndef SEQAN_HEADER_TEST_GRAPH_ALGORITHMS_H
#define SEQAN_HEADER_TEST_GRAPH_ALGORITHMS_H

using namespace std;
using namespace seqan;

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
void Test_BreadthFirstSearch() {
//____________________________________________________________________________
// Breadth-First Search
	typedef Graph<Undirected<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,4, 1,5, 2,5, 2,6, 2,3, 3,6, 3,7, 5,6, 6,7};
	
	//Create the graph
	TGraph g;
	addEdges(g,edges, numEdges);

	// Predecessor and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Bfs
	breadth_first_search(g, 1, predMap, distMap);
	
	SEQAN_TASSERT(getProperty(distMap, 0) == 1)
	SEQAN_TASSERT(getProperty(distMap, 1) == 0)
	SEQAN_TASSERT(getProperty(distMap, 2) == 2)
	SEQAN_TASSERT(getProperty(distMap, 3) == 3)
	SEQAN_TASSERT(getProperty(distMap, 4) == 2)
	SEQAN_TASSERT(getProperty(distMap, 5) == 1)
	SEQAN_TASSERT(getProperty(distMap, 6) == 2)
	SEQAN_TASSERT(getProperty(distMap, 7) == 3)
	SEQAN_TASSERT(getProperty(predMap, 0) == 1)
	SEQAN_TASSERT(getProperty(predMap, 1) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 2) == 5)
	SEQAN_TASSERT(getProperty(predMap, 3) == 6)
	SEQAN_TASSERT(getProperty(predMap, 4) == 0)
	SEQAN_TASSERT(getProperty(predMap, 5) == 1)
	SEQAN_TASSERT(getProperty(predMap, 6) == 5)
	SEQAN_TASSERT(getProperty(predMap, 7) == 6)
}

//////////////////////////////////////////////////////////////////////////////

void Test_DepthFirstSearch() {
//____________________________________________________________________________
// Depth-First Search
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 8;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,3, 0,1, 1,4, 2,4, 2,5, 3,1, 4,3, 5,5};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);

	// Predecessor and distance map
	String<unsigned int> predMap;
	String<unsigned int> discoveryTimeMap;
	String<unsigned int> finishingTimeMap;

	// Dfs
	depth_first_search(g, predMap, discoveryTimeMap, finishingTimeMap);

	SEQAN_TASSERT(getProperty(discoveryTimeMap, 0) == 1)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 1) == 2)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 2) == 9)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 3) == 4)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 4) == 3)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 5) == 10)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 0) == 8)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 1) == 7)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 2) == 12)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 3) == 5)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 4) == 6)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 5) == 11)
	SEQAN_TASSERT(getProperty(predMap, 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 1) == 0)
	SEQAN_TASSERT(getProperty(predMap, 2) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 3) == 4)
	SEQAN_TASSERT(getProperty(predMap, 4) == 1)
	SEQAN_TASSERT(getProperty(predMap, 5) == 2)
}

//////////////////////////////////////////////////////////////////////////////

void Test_TopologicalSort() {
//____________________________________________________________________________
// Topological Sort
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 9;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,3, 0,1, 1,2, 3,2, 5,7, 5,6, 6,7, 6,3, 8,7};
	
	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);

	// Predecessor and distance map
	String<TVertexDescriptor> order;
	
	// Topological sort
	topological_sort(g, order);

	SEQAN_TASSERT(getValue(order, 0) == 8)
	SEQAN_TASSERT(getValue(order, 1) == 5)
	SEQAN_TASSERT(getValue(order, 2) == 6)
	SEQAN_TASSERT(getValue(order, 3) == 7)
	SEQAN_TASSERT(getValue(order, 4) == 4)
	SEQAN_TASSERT(getValue(order, 5) == 0)
	SEQAN_TASSERT(getValue(order, 6) == 3)
	SEQAN_TASSERT(getValue(order, 7) == 1)
	SEQAN_TASSERT(getValue(order, 8) == 2)
}

//////////////////////////////////////////////////////////////////////////////

void Test_StronglyConnectedComponents() {
//____________________________________________________________________________
// Strongly-Connected-Components

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 14;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {1,0, 0,4, 2,1, 4,1, 5,1, 6,2, 3,2, 2,3, 7,3, 5,4, 6,5, 5,6, 7,6, 7,7};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);

	// Predecessor and distance map
	String<unsigned int> component;
	
	// Strongly Connected Components
	strongly_connected_components(g, component);

	SEQAN_TASSERT(getValue(component, 0) == 3)
	SEQAN_TASSERT(getValue(component, 1) == 3)
	SEQAN_TASSERT(getValue(component, 2) == 2)
	SEQAN_TASSERT(getValue(component, 3) == 2)
	SEQAN_TASSERT(getValue(component, 4) == 3)
	SEQAN_TASSERT(getValue(component, 5) == 1)
	SEQAN_TASSERT(getValue(component, 6) == 1)
	SEQAN_TASSERT(getValue(component, 7) == 0)
}

//////////////////////////////////////////////////////////////////////////////

void Test_PrimsAlgorithm() {
//____________________________________________________________________________
// Prim's algorithm
	typedef Graph<Undirected<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 14;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,6, 1,2, 1,6, 2,3, 2,4, 2,8, 3,5, 3,8, 4,6, 4,7, 5,8, 6,7, 7,8};
	unsigned int weights[] =    {4,   8,   8,   11,  7,   2,   4,   9,   14,  7,   6,   10,  1,   2  };

	//Create the graph
	TGraph g;
	addEdges(g,edges, numEdges);
	String<unsigned int> weightMap;
	resizeEdgeMap(g, weightMap, weights);

	// Tree and predecessor map 
	String<TVertexDescriptor> predMap;

	prims_algorithm(g, 0, weightMap, predMap);

	SEQAN_TASSERT(getProperty(predMap, 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 1) == 0)
	SEQAN_TASSERT(getProperty(predMap, 2) == 1)
	SEQAN_TASSERT(getProperty(predMap, 3) == 2)
	SEQAN_TASSERT(getProperty(predMap, 4) == 2)
	SEQAN_TASSERT(getProperty(predMap, 5) == 3)
	SEQAN_TASSERT(getProperty(predMap, 6) == 7)
	SEQAN_TASSERT(getProperty(predMap, 7) == 8)
	SEQAN_TASSERT(getProperty(predMap, 8) == 2)
}


//////////////////////////////////////////////////////////////////////////////

void Test_KruskalsAlgorithm() {
//____________________________________________________________________________
// Kruskal's algorithm
	typedef Graph<Undirected<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 14;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,6, 1,2, 1,6, 2,3, 2,4, 2,8, 3,5, 3,8, 4,6, 4,7, 5,8, 6,7, 7,8};
	unsigned int weights[] =    {4,   8,   8,   11,  7,   2,   4,   9,   14,  7,   6,   10,  1,   2  };

	//Create the graph
	TGraph g;
	addEdges(g,edges, numEdges);
	String<unsigned int> weightMap;
	resizeEdgeMap(g, weightMap, weights);

	// Tree edges
	String<TVertexDescriptor> treeEdges;
	kruskals_algorithm(g, 0, weightMap, treeEdges);

	SEQAN_TASSERT(getValue(treeEdges, 0) == 6)
	SEQAN_TASSERT(getValue(treeEdges, 1) == 7)
	SEQAN_TASSERT(getValue(treeEdges, 2) == 2)
	SEQAN_TASSERT(getValue(treeEdges, 3) == 4)
	SEQAN_TASSERT(getValue(treeEdges, 4) == 7)
	SEQAN_TASSERT(getValue(treeEdges, 5) == 8)
	SEQAN_TASSERT(getValue(treeEdges, 6) == 0)
	SEQAN_TASSERT(getValue(treeEdges, 7) == 1)
	SEQAN_TASSERT(getValue(treeEdges, 8) == 2)
	SEQAN_TASSERT(getValue(treeEdges, 9) == 8)
	SEQAN_TASSERT(getValue(treeEdges, 10) == 2)
	SEQAN_TASSERT(getValue(treeEdges, 11) == 3)
	SEQAN_TASSERT(getValue(treeEdges, 12) == 0)
	SEQAN_TASSERT(getValue(treeEdges, 13) == 6)
	SEQAN_TASSERT(getValue(treeEdges, 14) == 3)
	SEQAN_TASSERT(getValue(treeEdges, 15) == 5)
}

//////////////////////////////////////////////////////////////////////////////

void Test_DagShortestPath() {
//____________________________________________________________________________
// DAG-Shortest Paths
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,2, 0,1, 1,3, 1,2, 2,5, 2,4, 2,3, 3,5, 3,4, 4,5};
	int weights[] =             {3,   5,   6,   2,   2,   4,   7,   1,   -1,  -2};
	
	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);

	String<int> weightMap;
	resizeEdgeMap(g, weightMap, weights);

	// Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// DAG-Shortest path(Graph, sourceVertex_vertex, weightMap, predMap, distMap)
	dag_shortest_path(g,1,weightMap,predMap,distMap);
	
	SEQAN_TASSERT(getProperty(distMap, 0) == (unsigned int) getInfinityDistance(weightMap))
	SEQAN_TASSERT(getProperty(distMap, 1) == 0)
	SEQAN_TASSERT(getProperty(distMap, 2) == 2)
	SEQAN_TASSERT(getProperty(distMap, 3) == 6)
	SEQAN_TASSERT(getProperty(distMap, 4) == 5)
	SEQAN_TASSERT(getProperty(distMap, 5) == 3)
	SEQAN_TASSERT(getProperty(predMap, 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 1) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 2) == 1)
	SEQAN_TASSERT(getProperty(predMap, 3) == 1)
	SEQAN_TASSERT(getProperty(predMap, 4) == 3)
	SEQAN_TASSERT(getProperty(predMap, 5) == 4)
}

//////////////////////////////////////////////////////////////////////////////

void Test_BellmanFord() {
//____________________________________________________________________________
// Bellman-Ford

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,3, 1,2, 1,3, 2,4, 3,1, 3,2, 3,4, 4,0, 4,2};
	unsigned int weights[] =    {10,  5,   1,   2,   4,   3,   9,   2,   7,   6};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);
	
	String<unsigned int> weightMap;
	resizeEdgeMap(g,weightMap, weights);

	// Out parameters of Bellman-Ford: Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Bellman-Ford
	bool noNegativeCycle = bellman_ford_algorithm(g,0,weightMap,predMap,distMap);

	SEQAN_TASSERT(getProperty(distMap, 0) == 0)
	SEQAN_TASSERT(getProperty(distMap, 1) == 8)
	SEQAN_TASSERT(getProperty(distMap, 2) == 9)
	SEQAN_TASSERT(getProperty(distMap, 3) == 5)
	SEQAN_TASSERT(getProperty(distMap, 4) == 7)
	SEQAN_TASSERT(getProperty(predMap, 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 1) == 3)
	SEQAN_TASSERT(getProperty(predMap, 2) == 1)
	SEQAN_TASSERT(getProperty(predMap, 3) == 0)
	SEQAN_TASSERT(getProperty(predMap, 4) == 3)
	SEQAN_TASSERT(noNegativeCycle == true)
}

//////////////////////////////////////////////////////////////////////////////

void Test_Dijkstra() {
//____________________________________________________________________________
// Dijkstra

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,3, 1,2, 1,3, 2,4, 3,1, 3,2, 3,4, 4,0, 4,2};
	unsigned int weights[] =    {10,  5,   1,   2,   4,   3,   9,   2,   7,   6};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);

	String<unsigned int> weightMap;
	resizeEdgeMap(g , weightMap, weights);

	// Out parameters of Dijkstra: Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Dijkstra
	dijkstra(g,0,weightMap,predMap,distMap);

	SEQAN_TASSERT(getProperty(distMap, 0) == 0)
	SEQAN_TASSERT(getProperty(distMap, 1) == 8)
	SEQAN_TASSERT(getProperty(distMap, 2) == 9)
	SEQAN_TASSERT(getProperty(distMap, 3) == 5)
	SEQAN_TASSERT(getProperty(distMap, 4) == 7)
	SEQAN_TASSERT(getProperty(predMap, 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 1) == 3)
	SEQAN_TASSERT(getProperty(predMap, 2) == 1)
	SEQAN_TASSERT(getProperty(predMap, 3) == 0)
	SEQAN_TASSERT(getProperty(predMap, 4) == 3)
}

//////////////////////////////////////////////////////////////////////////////

void Test_AllPairsShortestPath() {
//____________________________________________________________________________
// All-Pairs Shortest Path

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 9;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,2, 0,4, 1,3, 1,4, 2,1, 3,0, 3,2, 4,3};
	int weights[] =    {3,   8,   -4,  1,   7,   4,   2,   -5,  6};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);
	
	String<int> weightMap;	
	resizeEdgeMap(g,weightMap, weights);

	// Out parameter
	Matrix<int> distMat;
	Matrix<TVertexDescriptor> predMat;

	// All-Pairs shortest path
	all_pairs_shortest_path(g,weightMap, distMat, predMat);

	unsigned int len = length(distMat, 0);
	SEQAN_TASSERT(getValue(distMat, 0*len + 0) == 0)
	SEQAN_TASSERT(getValue(distMat, 0*len + 1) == 1)
	SEQAN_TASSERT(getValue(distMat, 0*len + 2) == -3)
	SEQAN_TASSERT(getValue(distMat, 0*len + 3) == 2)
	SEQAN_TASSERT(getValue(distMat, 0*len + 4) == -4)
	SEQAN_TASSERT(getValue(distMat, 1*len + 0) == 3)
	SEQAN_TASSERT(getValue(distMat, 1*len + 1) == 0)
	SEQAN_TASSERT(getValue(distMat, 1*len + 2) == -4)
	SEQAN_TASSERT(getValue(distMat, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(distMat, 1*len + 4) == -1)
	SEQAN_TASSERT(getValue(distMat, 2*len + 0) == 7)
	SEQAN_TASSERT(getValue(distMat, 2*len + 1) == 4)
	SEQAN_TASSERT(getValue(distMat, 2*len + 2) == 0)
	SEQAN_TASSERT(getValue(distMat, 2*len + 3) == 5)
	SEQAN_TASSERT(getValue(distMat, 2*len + 4) == 3)
	SEQAN_TASSERT(getValue(distMat, 3*len + 0) == 2)
	SEQAN_TASSERT(getValue(distMat, 3*len + 1) == -1)
	SEQAN_TASSERT(getValue(distMat, 3*len + 2) == -5)
	SEQAN_TASSERT(getValue(distMat, 3*len + 3) == 0)
	SEQAN_TASSERT(getValue(distMat, 3*len + 4) == -2)
	SEQAN_TASSERT(getValue(distMat, 4*len + 0) == 8)
	SEQAN_TASSERT(getValue(distMat, 4*len + 1) == 5)
	SEQAN_TASSERT(getValue(distMat, 4*len + 2) == 1)
	SEQAN_TASSERT(getValue(distMat, 4*len + 3) == 6)
	SEQAN_TASSERT(getValue(distMat, 4*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 0*len + 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 0*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 0*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 0*len + 3) == 4)
	SEQAN_TASSERT(getValue(predMat, 0*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 1*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 1*len + 1) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 1*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(predMat, 1*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 2*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 2*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 2*len + 2) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 2*len + 3) == 1)
	SEQAN_TASSERT(getValue(predMat, 2*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 3*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 3*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 3*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 3*len + 3) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 3*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 4*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 4*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 4*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 4*len + 3) == 4)
	SEQAN_TASSERT(getValue(predMat, 4*len + 4) == getNilPredecessor(g))
}

//////////////////////////////////////////////////////////////////////////////

void Test_FloydWarshall() {
//____________________________________________________________________________
// Floyd Warshall

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 9;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,2, 0,4, 1,3, 1,4, 2,1, 3,0, 3,2, 4,3};
	int weights[] =    {3,   8,   -4,  1,   7,   4,   2,   -5,  6};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);
	
	String<int> weightMap;	
	resizeEdgeMap(g,weightMap, weights);

	// Out parameter
	Matrix<int> distMat;
	Matrix<TVertexDescriptor> predMat;

	// Floyd-Warshall
	floyd_warshall(g,weightMap, distMat, predMat);

	unsigned int len = length(distMat, 0);
	SEQAN_TASSERT(getValue(distMat, 0*len + 0) == 0)
	SEQAN_TASSERT(getValue(distMat, 0*len + 1) == 1)
	SEQAN_TASSERT(getValue(distMat, 0*len + 2) == -3)
	SEQAN_TASSERT(getValue(distMat, 0*len + 3) == 2)
	SEQAN_TASSERT(getValue(distMat, 0*len + 4) == -4)
	SEQAN_TASSERT(getValue(distMat, 1*len + 0) == 3)
	SEQAN_TASSERT(getValue(distMat, 1*len + 1) == 0)
	SEQAN_TASSERT(getValue(distMat, 1*len + 2) == -4)
	SEQAN_TASSERT(getValue(distMat, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(distMat, 1*len + 4) == -1)
	SEQAN_TASSERT(getValue(distMat, 2*len + 0) == 7)
	SEQAN_TASSERT(getValue(distMat, 2*len + 1) == 4)
	SEQAN_TASSERT(getValue(distMat, 2*len + 2) == 0)
	SEQAN_TASSERT(getValue(distMat, 2*len + 3) == 5)
	SEQAN_TASSERT(getValue(distMat, 2*len + 4) == 3)
	SEQAN_TASSERT(getValue(distMat, 3*len + 0) == 2)
	SEQAN_TASSERT(getValue(distMat, 3*len + 1) == -1)
	SEQAN_TASSERT(getValue(distMat, 3*len + 2) == -5)
	SEQAN_TASSERT(getValue(distMat, 3*len + 3) == 0)
	SEQAN_TASSERT(getValue(distMat, 3*len + 4) == -2)
	SEQAN_TASSERT(getValue(distMat, 4*len + 0) == 8)
	SEQAN_TASSERT(getValue(distMat, 4*len + 1) == 5)
	SEQAN_TASSERT(getValue(distMat, 4*len + 2) == 1)
	SEQAN_TASSERT(getValue(distMat, 4*len + 3) == 6)
	SEQAN_TASSERT(getValue(distMat, 4*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 0*len + 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 0*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 0*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 0*len + 3) == 4)
	SEQAN_TASSERT(getValue(predMat, 0*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 1*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 1*len + 1) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 1*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(predMat, 1*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 2*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 2*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 2*len + 2) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 2*len + 3) == 1)
	SEQAN_TASSERT(getValue(predMat, 2*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 3*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 3*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 3*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 3*len + 3) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 3*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 4*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 4*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 4*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 4*len + 3) == 4)
	SEQAN_TASSERT(getValue(predMat, 4*len + 4) == getNilPredecessor(g))
}


//////////////////////////////////////////////////////////////////////////////

void Test_TransitiveClosure() {
//____________________________________________________________________________
// Transitive Closure

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 5;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {3,0, 1,2, 2,1, 1,3, 3,2};
		
	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);

	// Transitive-Closure
	Matrix<bool> closure;
	transitive_closure(g,closure);
	
	unsigned int len = length(closure, 0);
	SEQAN_TASSERT(getValue(closure, 0*len + 0) == 1)
	SEQAN_TASSERT(getValue(closure, 0*len + 1) == 0)
	SEQAN_TASSERT(getValue(closure, 0*len + 2) == 0)
	SEQAN_TASSERT(getValue(closure, 0*len + 3) == 0)
	SEQAN_TASSERT(getValue(closure, 1*len + 0) == 1)
	SEQAN_TASSERT(getValue(closure, 1*len + 1) == 1)
	SEQAN_TASSERT(getValue(closure, 1*len + 2) == 1)
	SEQAN_TASSERT(getValue(closure, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(closure, 2*len + 0) == 1)
	SEQAN_TASSERT(getValue(closure, 2*len + 1) == 1)
	SEQAN_TASSERT(getValue(closure, 2*len + 2) == 1)
	SEQAN_TASSERT(getValue(closure, 2*len + 3) == 1)
	SEQAN_TASSERT(getValue(closure, 3*len + 0) == 1)
	SEQAN_TASSERT(getValue(closure, 3*len + 1) == 1)
	SEQAN_TASSERT(getValue(closure, 3*len + 2) == 1)
	SEQAN_TASSERT(getValue(closure, 3*len + 3) == 1)
}

//////////////////////////////////////////////////////////////////////////////

void Test_FordFulkerson() {
//____________________________________________________________________________
// Ford-Fulkerson
	typedef Graph<Directed<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,4, 1,2, 1,4, 2,3, 2,4, 4,1, 4,5, 5,2, 5,3};
	unsigned int capacity[] =    {16,  13,  12,  10,  20,  9,   4,   14,  7,   4};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);
	String<unsigned int> capMap;	
	resizeEdgeMap(g,capMap, capacity);

	// Out-parameter
	String<unsigned int> flow;	
	unsigned int valF = ford_fulkerson(g, 0, 3, capMap, flow);
	
	SEQAN_TASSERT(valF == 23)
	TEdgeIterator itEdge(g);
	for(;!atEnd(itEdge);goNext(itEdge)) {
		SEQAN_TASSERT(getProperty(flow, getValue(itEdge)) <= getProperty(capMap, getValue(itEdge)))
	}
}

//////////////////////////////////////////////////////////////////////////////

void Test_PathGrowingAlgorithm() {
//____________________________________________________________________________
// Path growing algorithm
	typedef Graph<Undirected<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 14;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,7, 0,5, 0,8, 1,5, 1,6, 2,6, 2,5, 2,8, 3,8, 3,4, 4,5, 4,6, 4,7, 4,8};
	unsigned int weights[] =    {20,  5,   19,  6,   7,   10,  3,   9,   12,  11,  13,  12,  9,   12};

	//Create the graph
	TGraph g;
	addEdges(g,edges, numEdges);
	String<unsigned int> weightMap;
	resizeEdgeMap(g, weightMap, weights);

	// Path growing algorithm
	String<bool> edgeMap;	
	resizeEdgeMap(g, edgeMap);

	// EdgeMap indicates whether an edge is selected or not
	unsigned int weight = path_growing_algorithm(g, weightMap, edgeMap);

	// Iterate over all the edges and output the ones that are selected
	std::cout << "Found matching of weight: " << weight << std::endl;
	std::cout << "Selected edges are: " << std::endl;
	TEdgeIterator it(g);
	for(;!atEnd(it);++it) {
		if (getProperty(edgeMap, *it) == true) {
			std::cout << '{' << sourceVertex(it) << ',' << targetVertex(it) << '}' << " - Weight: " << getProperty(weightMap, *it) << std::endl;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

void Test_GraphAlgorithms() {
//____________________________________________________________________________
// Graph Algorithms
	// Elementary graph algorithms
	Test_BreadthFirstSearch();
	Test_DepthFirstSearch();
	Test_TopologicalSort();
	Test_StronglyConnectedComponents();

	// Minimum Spanning Trees
	Test_PrimsAlgorithm();
	Test_KruskalsAlgorithm();

	// Single-Source shortest paths
	Test_DagShortestPath();
	Test_BellmanFord();
	Test_Dijkstra();

	// All-Pairs Shortest paths
	Test_AllPairsShortestPath();
	Test_FloydWarshall();
	Test_TransitiveClosure();

	//Maximum Flow
	Test_FordFulkerson();

	//Matching
	Test_PathGrowingAlgorithm();
}


}

#endif

