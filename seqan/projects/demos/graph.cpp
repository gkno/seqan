#include <seqan/graph.h>
#include <iostream>
#include <fstream>

#define TEST_PATH "projects/tests/graph/"
#define LIB_PATH "projects/library/seqan/graph/"

using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

void BreadthFirstSearch() {
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
	// Vertex names
	char names[] = {'r', 's', 't', 'u', 'v', 'w', 'x', 'y'};
	
	//Create the graph
	TGraph g;
	addEdges(g, edges, numEdges);
	std::cout << g << ::std::endl;
	
	String<char> nameMap;
	resizeVertexMap(g,nameMap, names);

	// Predecessor and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Bfs
	breadth_first_search(g, 1, predMap, distMap);
	
	// Output
	std::cout << "Breadth-First search: " << ::std::endl;
	typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Vertex " << getProperty(nameMap, getValue(it)) << ": ";
		if (getProperty(distMap, getValue(it))== getInfinityDistance(distMap)) {
			std::cout << "Not reachable!";
		} else {
			std::cout << "Level = " << getProperty(distMap, getValue(it));
		}
		typedef Value<String<unsigned int> >::Type TPredVal;
		TPredVal pre = getProperty(predMap, getValue(it));
		if (pre != getNilPredecessor(g)) {
			std::cout << ", Predecessor = " << getProperty(nameMap, pre) << ::std::endl;
		} else {
			std::cout << ", Predecessor = nil" << ::std::endl;
		}
		goNext(it);
	}
}

//////////////////////////////////////////////////////////////////////////////

void DepthFirstSearch() {
//____________________________________________________________________________
// Depth-First Search
	typedef Graph<> TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 8;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,3, 0,1, 1,4, 2,4, 2,5, 3,1, 4,3, 5,5};
	// Vertex names
	char names[] = {'u', 'v', 'w', 'x', 'y', 'z'};

	//Create the graph
	TGraph g;
	addEdges(g, edges, numEdges);
	std::cout << g << ::std::endl;
	String<char> nameMap;
	resizeVertexMap(g,nameMap, names);

	// Predecessor and distance map
	String<unsigned int> predMap;
	String<unsigned int> discoveryTimeMap;
	String<unsigned int> finishingTimeMap;

	// Dfs
	depth_first_search(g, predMap, discoveryTimeMap, finishingTimeMap);

	// Output
	std::cout << "Depth-First search: " << ::std::endl;
	typedef Iterator<Graph<>, VertexIterator>::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Vertex " << getProperty(nameMap, getValue(it)) << ": ";
		std::cout << "Discovery time = " << getProperty(discoveryTimeMap, getValue(it)) << ",";
		std::cout << "Finishing time = " << getProperty(finishingTimeMap, getValue(it)) << ",";
		typedef Value<String<unsigned int> >::Type TPredVal;
		TPredVal pre = getProperty(predMap, getValue(it));
		if (pre != getNilPredecessor(g)) {
			std::cout << "Predecessor = " << getProperty(nameMap, pre) << ::std::endl;
		} else {
			std::cout << "Predecessor = nil" << ::std::endl;
		}
		goNext(it);
	}
}

//////////////////////////////////////////////////////////////////////////////

void TopologicalSort() {
//____________________________________________________________________________
// Topological Sort
	typedef Graph<> TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 9;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,3, 0,1, 1,2, 3,2, 5,7, 5,6, 6,7, 6,3, 8,7};
	std::string names[] = {"shirt", "tie", "jacket", "belt", "watch", "undershorts", "pants", "shoes", "socks"};
	
	//Create the graph
	TGraph g;
	addEdges(g, edges, numEdges);
	std::cout << g << ::std::endl;
	String<std::string> nameMap;
	resizeVertexMap(g,nameMap, names);

	// Predecessor and distance map
	String<TVertexDescriptor> order;
	
	// Topological sort
	topological_sort(g, order);

	// Output
	std::cout << "Topological sort: " << ::std::endl;
	typedef Iterator<String<TVertexDescriptor> >::Type TStringIterator;
	TStringIterator it = begin(order);
	while(!atEnd(it)) {
		std::cout << getProperty(nameMap, getValue(it)) << ",";
		goNext(it);
	}
	std::cout << ::std::endl;
}

//////////////////////////////////////////////////////////////////////////////

void StronglyConnectedComponents() {
//____________________________________________________________________________
// Strongly-Connected-Components

	typedef Graph<> TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 14;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {1,0, 0,4, 2,1, 4,1, 5,1, 6,2, 3,2, 2,3, 7,3, 5,4, 6,5, 5,6, 7,6, 7,7};
	// Vertex names
	char names[] = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'};

	//Create the graph
	TGraph g;
	addEdges(g, edges, numEdges);
	std::cout << g << ::std::endl;
	String<char> nameMap;
	resizeVertexMap(g,nameMap, names);

	// Predecessor and distance map
	String<unsigned int> component;
	
	// Strongly Connected Components
	strongly_connected_components(g, component);

	// Output
	std::cout << "Strongly Connected Components: " << ::std::endl;
	typedef Iterator<Graph<>, VertexIterator>::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Vertex " << getProperty(nameMap, getValue(it)) << ": ";
		std::cout << "Component = " << getProperty(component, getValue(it)) << ::std::endl;
		goNext(it);
	}
}

//////////////////////////////////////////////////////////////////////////////

void PrimsAlgorithm() {
//____________________________________________________________________________
// Prim's algorithm
	typedef Graph<Undirected<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	// Number of edges
	TSize numEdges = 14;
	// Edges
	TVertexDescriptor edges[] = {0,1, 0,6, 1,2, 1,6, 2,3, 2,4, 2,8, 3,5, 3,8, 4,6, 4,7, 5,8, 6,7, 7,8};
	unsigned int weights[] =    {4,   8,   8,   11,  7,   2,   4,   9,   14,  7,   6,   10,  1,   2  };
	// Vertex names
	char names[] = {'a', 'b', 'c', 'd', 'i', 'e', 'h', 'g', 'f'};

	//Create the graph
	TGraph g;
	addEdges(g,edges, numEdges);
	String<int> weightMap;
	resizeEdgeMap(g, weightMap, weights);
	String<char> nameMap;
	resizeVertexMap(g,nameMap, names);
	std::cout << g << std::endl;

	// Tree and predecessor map 
	String<TVertexDescriptor> predMap;

	prims_algorithm(g, 0, weightMap, predMap);

	// Output
	std::cout << "Minimum Spanning Tree (Prim's algorithm): " << ::std::endl;
	typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Path from " << getProperty(nameMap, 0) << " to " << getProperty(nameMap, getValue(it)) << ": ";
		_print_path(g,predMap,(TVertexDescriptor) 0, getValue(it), nameMap);
		std::cout << ::std::endl;
		goNext(it);
	}
}

//////////////////////////////////////////////////////////////////////////////

void KruskalsAlgorithm() {
//____________________________________________________________________________
// Kruskal's algorithm
	typedef Graph<Undirected<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	// Number of edges
	TSize numEdges = 14;
	// Edges
	TVertexDescriptor edges[] = {0,1, 0,6, 1,2, 1,6, 2,3, 2,4, 2,8, 3,5, 3,8, 4,6, 4,7, 5,8, 6,7, 7,8};
	unsigned int weights[] =    {4,   8,   8,   11,  7,   2,   4,   9,   14,  7,   6,   10,  1,   2  };
	// Vertex names
	char names[] = {'a', 'b', 'c', 'd', 'i', 'e', 'h', 'g', 'f'};

	//Create the graph
	TGraph g;
	addEdges(g,edges, numEdges);
	String<int> weightMap;
	resizeEdgeMap(g, weightMap, weights);
	String<char> nameMap;
	resizeVertexMap(g,nameMap, names);
	std::cout << g << std::endl;

	// Tree edges
	String<TVertexDescriptor> treeEdges;
	kruskals_algorithm(g, 0, weightMap, treeEdges);

	
	// Output
	std::cout << "Minimum Spanning Tree (Kruskal's algorithm): " << ::std::endl;
	std::cout << "Tree Edges: ";
	typedef Iterator<String<TVertexDescriptor> >::Type TStrIterator;
	TStrIterator it = begin(treeEdges);
	while(!atEnd(it)) {
		std::cout << "(" << getProperty(nameMap,getValue(it)) << ",";
		goNext(it);
		std::cout << getProperty(nameMap,getValue(it)) << "), ";
		goNext(it);
	}
	std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

void DagShortestPath() {
//____________________________________________________________________________
// DAG-Shortest Paths
	typedef Graph<> TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,2, 0,1, 1,3, 1,2, 2,5, 2,4, 2,3, 3,5, 3,4, 4,5};
	int weights[] =             {3,   5,   6,   2,   2,   4,   7,   1,   -1,  -2};
	
	//Create the graph
	TGraph g;
	addEdges(g, edges, numEdges);
	String<int> weightMap;
	resizeEdgeMap(g,weightMap, weights);
	std::cout << g << ::std::endl;


	// Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// DAG-Shortest path(Graph, sourceVertex_vertex, weightMap, predMap, distMap)
	dag_shortest_path(g,1,weightMap,predMap,distMap);

	// Output
	std::cout << "Single-Source Shortest Paths in DAG: " << ::std::endl;
	typedef Iterator<Graph<>, VertexIterator>::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Path from 1 to " << getValue(it) << ": ";
		_print_path(g,predMap,(TVertexDescriptor) 1, getValue(it));
		std::cout << " (Distance: " << getProperty(distMap, getValue(it)) << ")" << ::std::endl;
		goNext(it);
	}
}

//////////////////////////////////////////////////////////////////////////////

void BellmanFord() {
//____________________________________________________________________________
// Bellman-Ford

	typedef Graph<> TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,3, 1,2, 1,3, 2,4, 3,1, 3,2, 3,4, 4,0, 4,2};
	unsigned int weights[] =    {10,  5,   1,   2,   4,   3,   9,   2,   7,   6};

	//Create the graph
	TGraph g;
	addEdges(g, edges, numEdges);
	String<unsigned int> weightMap;
	resizeEdgeMap(g,weightMap, weights);
	std::cout << g << ::std::endl;

	// Out parameters of Bellman-Ford: Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Bellman-Ford
	bool noNegativeCycle = bellman_ford_algorithm(g,0,weightMap,predMap,distMap);

	// Output
	std::cout << "Single-Source Shortest Paths: " << ::std::endl;
	std::cout << "Graph without negative cycles? " << noNegativeCycle << ::std::endl;
	typedef Iterator<Graph<>, VertexIterator>::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Path from 0 to " << getValue(it) << ": ";
		_print_path(g,predMap,(TVertexDescriptor) 0, getValue(it));
		std::cout << " (Distance: " << getProperty(distMap, getValue(it)) << ")" << ::std::endl;
		goNext(it);
	}
}

//////////////////////////////////////////////////////////////////////////////

void Dijkstra() {
//____________________________________________________________________________
// Dijkstra with external edge map

	typedef Graph<> TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,3, 1,2, 1,3, 2,4, 3,1, 3,2, 3,4, 4,0, 4,2};
	unsigned int weights[] =    {10,  5,   1,   2,   4,   3,   9,   2,   7,   6};

	//Create the graph
	TGraph g;
	addEdges(g, edges, numEdges);
	String<unsigned int> weightMap;
	resizeEdgeMap(g,weightMap, weights);
	std::cout << g << ::std::endl;

	// Out parameters of Dijkstra: Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Dijkstra
	dijkstra(g,0,weightMap,predMap,distMap);

	// Output
	std::cout << "Single-Source Shortest Paths: " << ::std::endl;
	typedef Iterator<Graph<>, VertexIterator>::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Path from 0 to " << getValue(it) << ": ";
		_print_path(g,predMap,(TVertexDescriptor) 0, getValue(it));
		std::cout << " (Distance: " << getProperty(distMap, getValue(it)) << ")" << ::std::endl;
		goNext(it);
	}
}

//////////////////////////////////////////////////////////////////////////////

void DijkstraInternalMap() {
//____________________________________________________________________________
// Dijkstra with internal property map
	//typedef SimpleTypeToMember<unsigned int> TEdgeCargo;
	typedef unsigned int TEdgeCargo;
	typedef Directed<TEdgeCargo> TEdges;
	//No edge ids require manual creation of the graph!!!
	//typedef Directed<TEdgeCargo, WithoutEdgeId> TEdges;

	typedef Graph<TEdges> TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,3, 1,2, 1,3, 2,4, 3,1, 3,2, 3,4, 4,0, 4,2};
	unsigned int weights[] =    {10,  5,   1,   2,   4,   3,   9,   2,   7,   6};

	//Create the graph
	// (1) Internal Map using member ids
	InternalMap<unsigned int> intMap;
	// (2) Internal Map using class
	//InternalPointerMap<unsigned int TEdgeCargo:: *, &TEdgeCargo::member> intMap;
	// (3) Internal Map using raw pointers to members
	//unsigned int TEdgeCargo:: * intMap = &TEdgeCargo::member;
	TGraph g;
	addEdges(g, edges, numEdges);
	resizeEdgeMap(g,intMap, weights);
	std::cout << g << ::std::endl;

	// Out parameters of Dijkstra: Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Dijkstra
	dijkstra(g,0,intMap,predMap,distMap);

	// Output
	std::cout << "Single-Source Shortest Paths: " << ::std::endl;
	typedef Iterator<Graph<TEdges> , VertexIterator>::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Path from 0 to " << getValue(it) << ": ";
		_print_path(g,predMap,(TVertexDescriptor) 0, getValue(it));
		std::cout << " (Distance: " << getProperty(distMap, getValue(it)) << ")" << ::std::endl;
		goNext(it);
	}
}

//////////////////////////////////////////////////////////////////////////////

void AllPairsShortestPath() {
//____________________________________________________________________________
// All-Pairs Shortest Path
	typedef Graph<> TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 9;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,2, 0,4, 1,3, 1,4, 2,1, 3,0, 3,2, 4,3};
	int weights[] =    {3,   8,   -4,  1,   7,   4,   2,   -5,  6};

	//Create the graph
	TGraph g;
	addEdges(g, edges, numEdges);
	String<int> weightMap;
	resizeEdgeMap(g,weightMap, weights);
	std::cout << g << ::std::endl;

	// Out parameter
	Matrix<int> distMat;
	Matrix<TVertexDescriptor> predMat;

	// All-Pairs shortest path
	all_pairs_shortest_path(g,weightMap, distMat, predMat);

	// Output
	unsigned int len = length(distMat, 0);
	for (TSize row=0;row < len;++row) {
		for (TSize col=0;col < len;++col) {
			std::cout << row << "," << col << " (Distance=" << getValue(distMat, row*len + col) << "): "; 
			_print_all_pairs_shortest_path(g,predMat,row,col);
			std::cout << ::std::endl;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

void FloydWarshall() {
//____________________________________________________________________________
// Floyd Warshall

	typedef Graph<> TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 9;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,2, 0,4, 1,3, 1,4, 2,1, 3,0, 3,2, 4,3};
	int weights[] =    {3,   8,   -4,  1,   7,   4,   2,   -5,  6};

	//Create the graph
	TGraph g;
	addEdges(g,edges, numEdges);
	String<int> weightMap;
	resizeEdgeMap(g,weightMap, weights);
	std::cout << g << ::std::endl;

	// Out parameter
	Matrix<int> distMat;
	Matrix<TVertexDescriptor> predMat;

	// Floyd-Warshall
	floyd_warshall(g,weightMap, distMat, predMat);

	// Output
	unsigned int len = length(distMat, 0);
	for (TSize row=0;row < len;++row) {
		for (TSize col=0;col < len;++col) {
			std::cout << row << "," << col << " (Distance=" << getValue(distMat, row*len + col) << "): "; 
			_print_all_pairs_shortest_path(g,predMat,row,col);
			std::cout << ::std::endl;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

void TransitiveClosure() {
//____________________________________________________________________________
// Transitive Closure
	typedef Graph<> TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 5;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {3,0, 1,2, 2,1, 1,3, 3,2};
	
	//Create the graph
	TGraph g;
	addEdges(g, edges, numEdges);
	std::cout << g << ::std::endl;

	// Transitive-Closure
	Matrix<bool> closure;
	transitive_closure(g,closure);

	TSize len = length(closure, 0);
	for (TSize row=0;row < len;++row) {
		for (TSize col=0;col < len;++col) {
			std::cout << getValue(closure, row*len+col) << ",";
		}
		std::cout << ::std::endl;
	}
}

//////////////////////////////////////////////////////////////////////////////

void FordFulkerson() {
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
	std::cout << g << ::std::endl;

	// Out-parameter
	String<unsigned int> flow;	
	unsigned int valF = ford_fulkerson(g, 0, 3, capMap, flow);
	
	// Output
	std::cout << "Ford-Fulkerson (Value of the flow = " << valF << ")" << ::std::endl;
	TEdgeIterator itEdge(g);
	for(;!atEnd(itEdge);goNext(itEdge)) {
		std::cout << "(" << sourceVertex(itEdge) << "," << targetVertex(itEdge) << "): ";
		std::cout << "Flow: " << getProperty(flow, getValue(itEdge)) << ", Capacity: " << getProperty(capMap, getValue(itEdge)) << ::std::endl;
	}
}

//////////////////////////////////////////////////////////////////////////////

void PathGrowingAlgorithm() {
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

void LongestIncreasingSubsequence() {
//____________________________________________________________________________
// Longest Increasing Subsequence
	
	// The sequence of numbers (Gusfield example)
	String<unsigned int> seq;
	appendValue(seq, 5); appendValue(seq, 3); appendValue(seq, 4);
	appendValue(seq, 9); appendValue(seq, 6); appendValue(seq, 2);
	appendValue(seq, 1); appendValue(seq, 8); appendValue(seq, 7);
	appendValue(seq, 10);

	// The string of positions belonging to the lis
	String<unsigned int, Block<> > pos;
	longestIncreasingSubsequence(seq,pos);

	// Output the lis
	for(int i = 0; i<(int) length(seq); ++i) {
		std::cout << seq[i] << ',';
	}
	std::cout << std::endl;
	std::cout << "Lis: " << std::endl;
	for(int i = length(pos)-1; i>=0; --i) {
		std::cout << seq[pos[i]] <<  ',';
	}
	std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

void LongestCommonSubsequence() {
//____________________________________________________________________________
// Longest Common Subsequence

	String<char> seq1("abacx");
	String<char> seq2("baabca");
	String<std::pair<unsigned int, unsigned int>, Block<> > pos;
	longestCommonSubsequence(seq1, seq2, pos);

	std::cout << seq1 << std::endl;
	std::cout << seq2 << std::endl;
	std::cout << "Lcs:" << std::endl;

	// Output the lcs
	for(int i = length(pos)-1; i>=0; --i) {
		std::cout << seq1[pos[i].first] <<  ',';
	}
	std::cout << std::endl;
	for(int i = length(pos)-1; i>=0; --i) {
		std::cout << seq2[pos[i].second] <<  ',';
	}
	std::cout << std::endl;
}


//////////////////////////////////////////////////////////////////////////////

void HeaviestIncreasingSubsequence() {
//____________________________________________________________________________
// Heaviest Increasing Subsequence
	String<char> seq("zeitgeist");
	String<unsigned int> weights;
	String<unsigned int> pos;
	fill(weights, length(seq), 1);
	assignProperty(weights, 2, 10);
	unsigned int w = heaviestIncreasingSubsequence(seq, weights, pos);

	// Output
	for(int i = 0; i< (int) length(seq); ++i) {
		std::cout << seq[i] << "(Weight=" << getProperty(weights, i) << "),";
	}
	std::cout << std::endl;
	std::cout << "His: " << std::endl;
	for(int i = length(pos)-1; i>=0; --i) {
		std::cout << seq[pos[i]] <<  ',';
	}
	std::cout << "(Weight=" << w << ')' << std::endl;
}


//////////////////////////////////////////////////////////////////////////////

void NeedlemanWunschDemo() {
	typedef String<char> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;

	// Ordinary DP alignment
	TStringSet str;
	TString str0("Myannealing");appendValue(str, str0);
	TString str1("annual"); appendValue(str, str1);
	TGraph g(str);
	Score<int> score_type = Score<int>(0,-1,-1,0);
	int score = globalAlignment(g, score_type, NeedlemanWunsch() );
	std::cout << "Scoring schema: Match=0, Mismatch=-1, Gap=-1" << std::endl;
	std::cout << g << std::endl;
	std::cout << "Score: " << score << std::endl;

	// Ends free-space alignment
	// AlignConfig<TTop, TLeft, TRight, TBottom>
	AlignConfig<true,false,false,true> ac;	// Gaps are free at the top and maximum is also searched for in the last row
	int score2 = globalAlignment(stringSet(g), score_type, ac, NeedlemanWunsch() );
	std::cout << "Score with ends free-space alignment: " << score2 << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

void GotohDemo() {
	typedef String<char> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	
	// Ordinary DP alignment
	TStringSet str;
	TString str0("TarfieldandGarfieldarestupid.");appendValue(str, str0);
	TString str1("Garfield");appendValue(str, str1);
	TGraph g(str);
	Score<int> score_type = Score<int>(2,-1,-1,-4);
	int score = globalAlignment(g, score_type, Gotoh());
	std::cout << "Scoring schema: Match=2, Mismatch=-1, Gap-extension=-1, Gap-opening=-4" << std::endl;
	std::cout << g << std::endl;
	std::cout << "Score: " << score << std::endl;

	// Ends free-space alignment
	// AlignConfig<TTop, TLeft, TRight, TBottom>
	AlignConfig<true,false,false,true> ac;	// Gaps are free at the top and maximum is also searched for in the last row
	int score2 = globalAlignment(g, score_type, ac, Gotoh() );
	std::cout << g << std::endl;
	std::cout << "Score with ends free-space alignment: " << score2 << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

void HirschbergDemo() {
	typedef String<char> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	
	// Ordinary DP alignment
	TStringSet str;
	TString str0("TarfieldandGarfieldarestupid.");appendValue(str, str0);
	TString str1("Garfield");appendValue(str, str1);
	TGraph g(str);
	Score<int> score_type = Score<int>(2,-1,-1,-4);
	int score = globalAlignment(g, score_type, Hirschberg());
	std::cout << "Scoring schema: Match=2, Mismatch=-1, Gap-extension=-1, Gap-opening=-4" << std::endl;
	std::cout << g << std::endl;
	std::cout << "Score: " << score << std::endl;

	// Hirschberg: Currently no ends free-space alignments possible
	// Always AlignConfig<false, false, false, false>
}

//////////////////////////////////////////////////////////////////////////////

void SmithWatermanDemo() {
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, Dna> > TGraph;

	
	TStringSet str;
	TString str0("TTGACACCCTCCCAATTGTA"); appendValue(str, str0);
	TString str1("ACCCCAGGCTTTACACAT"); appendValue(str, str1);

	// Well, as usual, we can use a graph, a file, or like now a fragment string
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	TFragmentString matches;
	Score<int> score_type = Score<int>(2,-1,-1,-2);
	int score = localAlignment(matches, str, score_type, SmithWaterman() );

	std::cout << "Scoring schema: Match=2, Mismatch=-1, Gap-extension=-1, Gap-opening=-2" << std::endl;
	std::cout << str0 << std::endl;
	std::cout << str1 << std::endl;
	std::cout << "Local match: " << std::endl;
	std::cout << label(matches[0], str, 0) << std::endl;
	std::cout << label(matches[0], str, 1) << std::endl;
	std::cout << "Score: " << score << std::endl;

}


//////////////////////////////////////////////////////////////////////////////

void AutomatonTest() {
//____________________________________________________________________________
// Automaton
	typedef Graph<Automaton<Dna> > DnaAutomaton;
	typedef VertexDescriptor<DnaAutomaton>::Type VertexDescriptorType;
	typedef EdgeDescriptor<DnaAutomaton>::Type EdgeDescriptorType;
		
	//Create the automaton
	DnaAutomaton automaton;

	// Property maps used during graph creation
	// Do not use this way unless you really need to
	String<char> propMap;
	reserve(propMap, 6);		// Maximum of 6 vertices
	String<char> propMapEdges; 
	reserve(propMapEdges, 100); // Maximum of edges does not work here 
								// because ids depend on size of alphabet and number of vertices
								// Size of alphabet * number of vertices is safe
	
	// Add vertices and edges
	VertexDescriptorType rootVertex = addVertex(automaton); // 0
	assignProperty(propMap, rootVertex, 'a');
	VertexDescriptorType v = addVertex(automaton); // 1
	assignProperty(propMap, v, 'b');
	v = addVertex(automaton); // 2
	assignProperty(propMap, v, 'c');
	v = addVertex(automaton); // 3
	assignProperty(propMap, v, 'd');
	v = addVertex(automaton); // 4
	assignProperty(propMap, v, 'e');
	v = addVertex(automaton); // 5
	assignProperty(propMap, v, 'f');
	EdgeDescriptorType edge1 = addEdge(automaton,0,1,'C');
	assignProperty(propMapEdges, edge1, 'a');
	EdgeDescriptorType edge3 = addEdge(automaton,1,0,'A');
	assignProperty(propMapEdges, edge3, 'c');
	EdgeDescriptorType e = addEdge(automaton,4,0,'G');
	assignProperty(propMapEdges, e, 'h');
	e = addEdge(automaton,0,3,'T');
	assignProperty(propMapEdges, e, 'b');
	e = addEdge(automaton,1,1,'G');
	assignProperty(propMapEdges, e, 'd');
	e = addEdge(automaton,1,2,'T');
	assignProperty(propMapEdges, e, 'e');
	e = addEdge(automaton,5,1,'T');
	assignProperty(propMapEdges, e, 'i');
	e = addEdge(automaton,2,5,'C');
	assignProperty(propMapEdges, e, 'f');
	e = addEdge(automaton,3,4,'G');
	assignProperty(propMapEdges, e, 'g');
	e = addEdge(automaton,5,3,'A');
	assignProperty(propMapEdges, e, 'j');

	// Property maps used after graph creation -> Much better :-)
	String<std::string> nameMap;
	resizeVertexMap(automaton, nameMap);
	typedef Iterator<DnaAutomaton, VertexIterator>::Type TVertexIterator;
	TVertexIterator it(automaton);
	std::string str = "Hallo";
	for(;!atEnd(it);goNext(it)) {
		assignProperty(nameMap, getValue(it), ::std::string(str));
		str = str.append("o");
	}

	// Print automaton
	std::cout << automaton << ::std::endl;

	// Perform walks on it
	std::cout << "A simple forward walk through the automaton:" << ::std::endl;
	std::cout << rootVertex << " -T-> ";
	VertexDescriptorType succ = getSuccessor(automaton,rootVertex,'T');
	std::cout << succ << " -G-> ";
	succ = getSuccessor(automaton,succ,'G');
	std::cout << succ << " -G-> ";
	succ = getSuccessor(automaton,succ,'G');
	std::cout << succ << " -C-> ";
	succ = getSuccessor(automaton,succ,'C');
	std::cout << succ << ::std::endl;
	std::cout << "------------------------------" << ::std::endl;
	std::cout << "Multiple forward transitions at once:" << ::std::endl;
	std::cout << rootVertex << " -TGGC-> ";
	succ = parseString(automaton,rootVertex,"TGGC");
	std::cout << succ << ::std::endl;

	goBegin(it);
	std::cout << "All vertices: ";
	for(;!atEnd(it);goNext(it)) {
		std::cout << getValue(it) << "(" << getProperty(nameMap, getValue(it)) << "," << getProperty(propMap, getValue(it)) << "), ";
	}
	std::cout << ::std::endl;
	std::cout << "Some edge properties:" << ::std::endl;
	std::cout << getProperty(propMapEdges, edge1) << ::std::endl;
	std::cout << getProperty(propMapEdges, edge3) << ::std::endl;
}

//////////////////////////////////////////////////////////////////////////////

void NeighborJoining() {
//____________________________________________________________________________
// Neighbor Joining

	// Create a distance matrix
	String<double> mat;
	fill(mat, 8*8, 0);
	assignValue(mat, 0*8+0, 0);assignValue(mat, 0*8+1, 7);assignValue(mat, 0*8+2, 8);assignValue(mat, 0*8+3, 11);
	assignValue(mat, 0*8+4, 13);assignValue(mat, 0*8+5, 16);assignValue(mat, 0*8+6, 13);assignValue(mat, 0*8+7, 17);
	assignValue(mat, 1*8+0, 7);assignValue(mat, 1*8+1, 0);assignValue(mat, 1*8+2, 5);assignValue(mat, 1*8+3, 8);
	assignValue(mat, 1*8+4, 10);assignValue(mat, 1*8+5, 13);assignValue(mat, 1*8+6, 10);assignValue(mat, 1*8+7, 14);
	assignValue(mat, 2*8+0, 8);assignValue(mat, 2*8+1, 5);assignValue(mat, 2*8+2, 0);assignValue(mat, 2*8+3, 5);
	assignValue(mat, 2*8+4, 7);assignValue(mat, 2*8+5, 10);assignValue(mat, 2*8+6, 7);assignValue(mat, 2*8+7, 11);
	assignValue(mat, 3*8+0, 11);assignValue(mat, 3*8+1, 8);assignValue(mat, 3*8+2, 5);assignValue(mat, 3*8+3, 0);
	assignValue(mat, 3*8+4, 8);assignValue(mat, 3*8+5, 11);assignValue(mat, 3*8+6, 8);assignValue(mat, 3*8+7, 12);
	assignValue(mat, 4*8+0, 13);assignValue(mat, 4*8+1, 10);assignValue(mat, 4*8+2, 7);assignValue(mat, 4*8+3, 8);
	assignValue(mat, 4*8+4, 0);assignValue(mat, 4*8+5, 5);assignValue(mat, 4*8+6, 6);assignValue(mat, 4*8+7, 10);
	assignValue(mat, 5*8+0, 16);assignValue(mat, 5*8+1, 13);assignValue(mat, 5*8+2, 10);assignValue(mat, 5*8+3, 11);
	assignValue(mat, 5*8+4, 5);assignValue(mat, 5*8+5, 0);assignValue(mat, 5*8+6, 9);assignValue(mat, 5*8+7, 13);
	assignValue(mat, 6*8+0, 13);assignValue(mat, 6*8+1, 10);assignValue(mat, 6*8+2, 7);assignValue(mat, 6*8+3, 8);
	assignValue(mat, 6*8+4, 6);assignValue(mat, 6*8+5, 9);assignValue(mat, 6*8+6, 0);assignValue(mat, 6*8+7, 8);
	assignValue(mat, 7*8+0, 17);assignValue(mat, 7*8+1, 14);assignValue(mat, 7*8+2, 11);assignValue(mat, 7*8+3, 12);
	assignValue(mat, 7*8+4, 10);assignValue(mat, 7*8+5, 13);assignValue(mat, 7*8+6, 8);assignValue(mat, 7*8+7, 0);

	typedef Graph<Tree<double> > TGraph;
	TGraph njTreeOut;
	slowNjTree(mat, njTreeOut);

	// Print the tree
	// std::cout << njTreeOut << std::endl;
	// or:
	std::fstream strm; 
	strm.open(TEST_PATH "my_nj_tree.dot", std::ios_base::out | std::ios_base::trunc);
	write(strm,njTreeOut,DotDrawing());
	strm.close();


	// Dfs Iterator
	typedef Iterator<TGraph, DfsPreorder>::Type TDfsPreorder;
	TDfsPreorder dfsIt(njTreeOut,getRoot(njTreeOut));
	std::cout << "Dfs: ";
	for(;!atEnd(dfsIt);goNext(dfsIt)) {
		std::cout << *dfsIt << ", ";
	}
	std::cout << std::endl;

	typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIter;
	TEdgeIter edIt(njTreeOut);
	std::cout << "Edges: " << std::endl;
	for(;!atEnd(edIt);goNext(edIt)) {
		std::cout << sourceVertex(edIt) << " -- " << targetVertex(edIt) << ": Weight = " << getCargo(*edIt) << std::endl;
	}

}

//////////////////////////////////////////////////////////////////////////////

void Upgma() {
//____________________________________________________________________________
// Upgma

	// Create a distance matrix
	String<double> mat;
	fill(mat, 8*8, 0);
	assignValue(mat, 0*8+0, 0);assignValue(mat, 0*8+1, 32);assignValue(mat, 0*8+2, 48);assignValue(mat, 0*8+3, 51);
	assignValue(mat, 0*8+4, 50);assignValue(mat, 0*8+5, 48);assignValue(mat, 0*8+6, 98);assignValue(mat, 0*8+7, 148);
	assignValue(mat, 1*8+0, 32);assignValue(mat, 1*8+1, 0);assignValue(mat, 1*8+2, 26);assignValue(mat, 1*8+3, 34);
	assignValue(mat, 1*8+4, 29);assignValue(mat, 1*8+5, 33);assignValue(mat, 1*8+6, 84);assignValue(mat, 1*8+7, 136);
	assignValue(mat, 2*8+0, 48);assignValue(mat, 2*8+1, 26);assignValue(mat, 2*8+2, 0);assignValue(mat, 2*8+3, 42);
	assignValue(mat, 2*8+4, 44);assignValue(mat, 2*8+5, 44);assignValue(mat, 2*8+6, 92);assignValue(mat, 2*8+7, 152);
	assignValue(mat, 3*8+0, 51);assignValue(mat, 3*8+1, 34);assignValue(mat, 3*8+2, 42);assignValue(mat, 3*8+3, 0);
	assignValue(mat, 3*8+4, 44);assignValue(mat, 3*8+5, 38);assignValue(mat, 3*8+6, 86);assignValue(mat, 3*8+7, 142);
	assignValue(mat, 4*8+0, 50);assignValue(mat, 4*8+1, 29);assignValue(mat, 4*8+2, 44);assignValue(mat, 4*8+3, 44);
	assignValue(mat, 4*8+4, 0);assignValue(mat, 4*8+5, 24);assignValue(mat, 4*8+6, 89);assignValue(mat, 4*8+7, 142);
	assignValue(mat, 5*8+0, 48);assignValue(mat, 5*8+1, 33);assignValue(mat, 5*8+2, 44);assignValue(mat, 5*8+3, 38);
	assignValue(mat, 5*8+4, 24);assignValue(mat, 5*8+5, 0);assignValue(mat, 5*8+6, 90);assignValue(mat, 5*8+7, 142);
	assignValue(mat, 6*8+0, 98);assignValue(mat, 6*8+1, 84);assignValue(mat, 6*8+2, 92);assignValue(mat, 6*8+3, 86);
	assignValue(mat, 6*8+4, 89);assignValue(mat, 6*8+5, 90);assignValue(mat, 6*8+6, 0);assignValue(mat, 6*8+7, 148);
	assignValue(mat, 7*8+0, 148);assignValue(mat, 7*8+1, 136);assignValue(mat, 7*8+2, 152);assignValue(mat, 7*8+3, 142);
	assignValue(mat, 7*8+4, 142);assignValue(mat, 7*8+5, 142);assignValue(mat, 7*8+6, 148);assignValue(mat, 7*8+7, 0);

	typedef Graph<Tree<double> > TGraph;
	TGraph njTreeOut;
	upgmaTree(mat, njTreeOut);

	// Print the tree
	// std::cout << njTreeOut << std::endl;
	// or:
	std::fstream strm; 
	strm.open(TEST_PATH "my_nj_tree.dot", std::ios_base::out | std::ios_base::trunc);
	write(strm,njTreeOut,DotDrawing());
	strm.close();


	// Dfs Iterator
	typedef Iterator<TGraph, DfsPreorder>::Type TDfsPreorder;
	TDfsPreorder dfsIt(njTreeOut,getRoot(njTreeOut));
	std::cout << "Dfs: ";
	for(;!atEnd(dfsIt);goNext(dfsIt)) {
		std::cout << *dfsIt << ", ";
	}
	std::cout << std::endl;

	typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIter;
	TEdgeIter edIt(njTreeOut);
	std::cout << "Edges: " << std::endl;
	for(;!atEnd(edIt);goNext(edIt)) {
		std::cout << sourceVertex(edIt) << " -- " << targetVertex(edIt) << ": Weight = " << getCargo(*edIt) << std::endl;
	}

}

//////////////////////////////////////////////////////////////////////////////

void TCoffee() {
//____________________________________________________________________________
// T-Coffee
	typedef String<AminoAcid> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int, Default> > TGraph;
	TString str1 = "GARFIELDTHELASTFATCAT";
	TString str2 = "GARFIELDTHEFASTCAT";
	TString str3 = "GARFIELDTHEVERYFASTCAT";
	TString str4 = "THEFATCAT";
	TStringSet strSet;
	assignValueById(strSet, str1);
	assignValueById(strSet, str2);
	assignValueById(strSet, str3);
	assignValueById(strSet, str4);

	// Score object
	Score<int> score_type_global = Score<int>(2,-1,-1,-2);
	Score<int> score_type_local = Score<int>(2,-1,-1,-2);

	// Generate a primary library, i.e., all global pairwise alignments
	TGraph lib1(strSet);
	String<double> distanceMatrix; 
	generatePrimaryLibrary(lib1, distanceMatrix, score_type_global, GlobalPairwise_Library() );

	// Generate a primary library, i.e., all local pairwise alignments
	TGraph lib2(strSet);
	generatePrimaryLibrary(lib2, score_type_local, LocalPairwise_Library() );

	// Weighting of libraries (Signal addition)
	TGraph g(strSet);
	String<TGraph*> libs;
	appendValue(libs, &lib1);
	appendValue(libs, &lib2);
	combineGraphs(g, libs);

	// Triplet library extension
	tripletLibraryExtension(g);

	// Build the neighbor joining tree
	Graph<Tree<double> > njTreeOut;
	slowNjTree(distanceMatrix, njTreeOut);

	// Perform a progressive alignment
	Graph<Alignment<TStringSet, void> > gOut(strSet);
	iterativeProgressiveAlignment(g, njTreeOut, score_type_local, gOut);

	// Print the alignment
	std::cout << gOut << std::endl;
}


//////////////////////////////////////////////////////////////////////////////

int main ()
{
//____________________________________________________________________________
// Graph Algorithms
	// Elementary graph algorithms
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Breadth-Frist Search-----------" << ::std::endl;
	BreadthFirstSearch();
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Depth-Frist Search-------------" << ::std::endl;
	DepthFirstSearch();
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Topological Sort---------------" << ::std::endl;
	TopologicalSort();
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Strongly-Connected-Components--" << ::std::endl;
	StronglyConnectedComponents();

	// Minimum Spanning Trees
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Prim's algorithm---------------" << ::std::endl;
	PrimsAlgorithm();
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Kruskal's algorithm---------------" << ::std::endl;
	KruskalsAlgorithm();


	// Single-Source shortest paths
	std::cout << "===================================" << ::std::endl;
	std::cout << "----DAG-Shortest Path--------------" << ::std::endl;
	DagShortestPath();
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Bellman-Ford-------------------" << ::std::endl;
	BellmanFord();
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Dijkstra-----------------------" << ::std::endl;
	Dijkstra();
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Dijkstra (Internal Map)--------" << ::std::endl;
	DijkstraInternalMap();

	// All-Pairs Shortest paths
	std::cout << "===================================" << ::std::endl;
	std::cout << "----All-Pairs Shortest Path--------" << ::std::endl;
	AllPairsShortestPath();
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Floyd-Warshall-----------------" << ::std::endl;
	FloydWarshall();
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Transitive-Closure-------------" << ::std::endl;
	TransitiveClosure();

	// Maximum Flow
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Ford-Fulkerson-----------------" << ::std::endl;
	FordFulkerson();

	// Matching
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Path growing algorithm---------" << ::std::endl;
	PathGrowingAlgorithm();

	// Lis, lcs, his
	std::cout << "===================================" << ::std::endl;
	std::cout << "--Longest Increasing Subsequence---" << ::std::endl;
	LongestIncreasingSubsequence();
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Longest Common Subsequence-----" << ::std::endl;
	LongestCommonSubsequence();
	std::cout << "===================================" << ::std::endl;
	std::cout << "--Heaviest Increasing Subsequence--" << ::std::endl;
	HeaviestIncreasingSubsequence();

	// Global Alignments
	std::cout << "===================================" << ::std::endl;
	std::cout << "--------Needleman Wunsch-----------" << ::std::endl;
	NeedlemanWunschDemo();
	std::cout << "===================================" << ::std::endl;
	std::cout << "-------------Gotoh-----------------" << ::std::endl;
	GotohDemo();
	std::cout << "===================================" << ::std::endl;
	std::cout << "----------Hirschberg---------------" << ::std::endl;
	HirschbergDemo();

	// Local Alignments
	std::cout << "===================================" << ::std::endl;
	std::cout << "--------Smith Waterman-------------" << ::std::endl;
	SmithWatermanDemo();


	// Automaton
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Automaton----------------------" << ::std::endl;
	AutomatonTest();

	// Guide Tree Construction Methods
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Neighbor Joining---------------" << ::std::endl;
	NeighborJoining();
	std::cout << "===================================" << ::std::endl;
	std::cout << "---------UPGMA---------------------" << ::std::endl;
	Upgma();

	// T-Coffee
	std::cout << "===================================" << ::std::endl;
	std::cout << "----T-Coffee-----------------------" << ::std::endl;
	TCoffee();

	return 0;
}
