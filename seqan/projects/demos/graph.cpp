#include <seqan/graph.h>
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

void BreadthFirstSearch() {
//____________________________________________________________________________
// Breadth-First Search

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 20;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 1,0, 0,4, 4,0, 1,5, 5,1, 2,5, 5,2, 2,6, 6,2, 2,3, 3,2, 3,6, 6,3, 3,7, 7,3, 5,6, 6,5, 6,7, 7,6};
	// Vertex names
	char names[] = {'r', 's', 't', 'u', 'v', 'w', 'x', 'y'};
	
	//Create the graph
	Graph<> g(edges, numEdges);
	std::cout << g << std::endl;
	String<char> nameMap;
	initVertexMap(g,nameMap, names);

	// Predecessor and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Bfs
	breadth_first_search(g, 1, predMap, distMap);
	
	// Output
	std::cout << "Breadth-First search: " << std::endl;
	typedef Iterator<Graph<>, VertexIterator<> >::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Vertex " << getProperty(nameMap, getValue(it)) << ": ";
		if (getProperty(distMap, getValue(it))== getInfinityDistance(EmptyMap())) {
			std::cout << "Not reachable!";
		} else {
			std::cout << "Level = " << getProperty(distMap, getValue(it));
		}
		typedef Value<String<unsigned int> >::Type TPredVal;
		TPredVal pre = getProperty(predMap, getValue(it));
		if (pre != getNilPredecessor(g)) {
			std::cout << ", Predecessor = " << getProperty(nameMap, pre) << std::endl;
		} else {
			std::cout << ", Predecessor = nil" << std::endl;
		}
		goNext(it);
	}
}

//////////////////////////////////////////////////////////////////////////////

void DepthFirstSearch() {
//____________________________________________________________________________
// Depth-First Search

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 8;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,3, 0,1, 1,4, 2,4, 2,5, 3,1, 4,3, 5,5};
	// Vertex names
	char names[] = {'u', 'v', 'w', 'x', 'y', 'z'};

	//Create the graph
	Graph<> g(edges, numEdges);
	std::cout << g << std::endl;
	String<char> nameMap;
	initVertexMap(g,nameMap, names);

	// Predecessor and distance map
	String<unsigned int> predMap;
	String<unsigned int> discoveryTimeMap;
	String<unsigned int> finishingTimeMap;

	// Dfs
	depth_first_search(g, predMap, discoveryTimeMap, finishingTimeMap);

	// Output
	std::cout << "Depth-First search: " << std::endl;
	typedef Iterator<Graph<>, VertexIterator<> >::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Vertex " << getProperty(nameMap, getValue(it)) << ": ";
		std::cout << "Discovery time = " << getProperty(discoveryTimeMap, getValue(it)) << ",";
		std::cout << "Finishing time = " << getProperty(finishingTimeMap, getValue(it)) << ",";
		typedef Value<String<unsigned int> >::Type TPredVal;
		TPredVal pre = getProperty(predMap, getValue(it));
		if (pre != getNilPredecessor(g)) {
			std::cout << "Predecessor = " << getProperty(nameMap, pre) << std::endl;
		} else {
			std::cout << "Predecessor = nil" << std::endl;
		}
		goNext(it);
	}
}

//////////////////////////////////////////////////////////////////////////////

void TopologicalSort() {
//____________________________________________________________________________
// Topological Sort
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 9;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,3, 0,1, 1,2, 3,2, 5,7, 5,6, 6,7, 6,3, 8,7};
	std::string names[] = {"shirt", "tie", "jacket", "belt", "watch", "undershorts", "pants", "shoes", "socks"};
	
	//Create the graph
	Graph<> g(edges, numEdges);
	std::cout << g << std::endl;
	String<std::string> nameMap;
	initVertexMap(g,nameMap, names);

	// Predecessor and distance map
	String<TVertexDescriptor> order;
	
	// Topological sort
	topological_sort(g, order);

	// Output
	std::cout << "Topological sort: " << std::endl;
	typedef Iterator<String<TVertexDescriptor> >::Type TStringIterator;
	TStringIterator it = begin(order);
	while(!atEnd(it)) {
		std::cout << getProperty(nameMap, getValue(it)) << ",";
		goNext(it);
	}
	std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

void StronglyConnectedComponents() {
//____________________________________________________________________________
// Strongly-Connected-Components

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 14;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {1,0, 0,4, 2,1, 4,1, 5,1, 6,2, 3,2, 2,3, 7,3, 5,4, 6,5, 5,6, 7,6, 7,7};
	// Vertex names
	char names[] = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'};

	//Create the graph
	Graph<> g(edges, numEdges);
	std::cout << g << std::endl;
	String<char> nameMap;
	initVertexMap(g,nameMap, names);

	// Predecessor and distance map
	String<unsigned int> component;
	
	// Strongly Connected Components
	strongly_connected_components(g, component);

	// Output
	std::cout << "Strongly Connected Components: " << std::endl;
	typedef Iterator<Graph<>, VertexIterator<> >::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Vertex " << getProperty(nameMap, getValue(it)) << ": ";
		std::cout << "Component = " << getProperty(component, getValue(it)) << std::endl;
		goNext(it);
	}
}

//////////////////////////////////////////////////////////////////////////////

void DagShortestPath() {
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
	Graph<> g(edges, numEdges);
	std::cout << g << std::endl;
	String<int> weightMap;
	initEdgeMap(g,weightMap, weights);

	// Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// DAG-Shortest path(Graph, sourceVertex_vertex, weightMap, predMap, distMap)
	dag_shortest_path(g,1,weightMap,predMap,distMap);

	// Output
	std::cout << "Single-Source Shortest Paths in DAG: " << std::endl;
	typedef Iterator<Graph<>, VertexIterator<> >::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Path from 1 to " << getValue(it) << ": ";
		_print_path(g,predMap,(TVertexDescriptor) 1, getValue(it));
		std::cout << " (Distance: " << getProperty(distMap, getValue(it)) << ")" << std::endl;
		goNext(it);
	}
}

//////////////////////////////////////////////////////////////////////////////

void BellmanFord() {
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
	Graph<> g(edges, numEdges);
	std::cout << g << std::endl;

	// In parameters of Bellman-Ford
	// Define the distances or weights of the edges
	String<unsigned int> weightMap;
	initEdgeMap(g,weightMap, weights);

	// Out parameters of Bellman-Ford: Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Bellman-Ford
	bool noNegativeCycle = bellman_ford_algorithm(g,0,weightMap,predMap,distMap);

	// Output
	std::cout << "Single-Source Shortest Paths: " << std::endl;
	std::cout << "Graph without negative cycles? " << noNegativeCycle << std::endl;
	typedef Iterator<Graph<>, VertexIterator<> >::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Path from 0 to " << getValue(it) << ": ";
		_print_path(g,predMap,(TVertexDescriptor) 0, getValue(it));
		std::cout << " (Distance: " << getProperty(distMap, getValue(it)) << ")" << std::endl;
		goNext(it);
	}
}

//////////////////////////////////////////////////////////////////////////////

void Dijkstra() {
//____________________________________________________________________________
// Dijkstra with external edge map

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,3, 1,2, 1,3, 2,4, 3,1, 3,2, 3,4, 4,0, 4,2};
	unsigned int weights[] =    {10,  5,   1,   2,   4,   3,   9,   2,   7,   6};

	//Create the graph
	Graph<> g(edges, numEdges);
	std::cout << g << std::endl;

	// In parameters of Dijkstra
	// Define the distances or weights of the edges
	String<unsigned int> weightMap;
	initEdgeMap(g,weightMap, weights);

	// Out parameters of Dijkstra: Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Dijkstra
	dijkstra(g,0,weightMap,predMap,distMap);

	// Output
	std::cout << "Single-Source Shortest Paths: " << std::endl;
	typedef Iterator<Graph<>, VertexIterator<> >::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Path from 0 to " << getValue(it) << ": ";
		_print_path(g,predMap,(TVertexDescriptor) 0, getValue(it));
		std::cout << " (Distance: " << getProperty(distMap, getValue(it)) << ")" << std::endl;
		goNext(it);
	}
}

//////////////////////////////////////////////////////////////////////////////

void DijkstraInternalMap() {
//____________________________________________________________________________
// Dijkstra with internal property map
	//typedef SimpleTypeToMember<unsigned int> TEdgeCargo;
	typedef unsigned int TEdgeCargo;
	typedef EdgeList<TEdgeCargo> TEdges;
	//No edge ids require manual creation of the graph!!!
	//typedef EdgeList<TEdgeCargo, WithoutEdgeId> TEdges;
	typedef VertexDescriptor<Graph<TEdges> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<TEdges> >::Type TEdgeDescriptor;
	typedef Size<Graph<TEdges> >::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,3, 1,2, 1,3, 2,4, 3,1, 3,2, 3,4, 4,0, 4,2};
	unsigned int weights[] =    {10,  5,   1,   2,   4,   3,   9,   2,   7,   6};

	//Create the graph
	Graph<TEdges> g(edges, numEdges);
	std::cout << g << std::endl;

	// In parameters of Dijkstra
	// Define the distances or weights of the edges
	// (1) Internal Map using member ids
	InternalMap<unsigned int> intMap;
	// (2) Internal Map using class
	//InternalPointerMap<unsigned int TEdgeCargo:: *, &TEdgeCargo::member> intMap;
	// (3) Internal Map using raw pointers to members
	//unsigned int TEdgeCargo:: * intMap = &TEdgeCargo::member;
	initEdgeMap(g,intMap, weights);

	// Out parameters of Dijkstra: Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Dijkstra
	dijkstra(g,0,intMap,predMap,distMap);

	// Output
	std::cout << "Single-Source Shortest Paths: " << std::endl;
	typedef Iterator<Graph<TEdges> , VertexIterator<> >::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Path from 0 to " << getValue(it) << ": ";
		_print_path(g,predMap,(TVertexDescriptor) 0, getValue(it));
		std::cout << " (Distance: " << getProperty(distMap, getValue(it)) << ")" << std::endl;
		goNext(it);
	}
}

//////////////////////////////////////////////////////////////////////////////

void AllPairsShortestPath() {
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
	Graph<> g(edges, numEdges);
	std::cout << g << std::endl;

	// In parameters
	String<int> weightMap;
	initEdgeMap(g,weightMap, weights);
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
			std::cout << std::endl;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

void FloydWarshall() {
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
	Graph<> g(edges, numEdges);
	std::cout << g << std::endl;

	// In parameters
	String<int> weightMap;
	initEdgeMap(g,weightMap, weights);
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
			std::cout << std::endl;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

void TransitiveClosure() {
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
	Graph<> g(edges, numEdges);
	std::cout << g << std::endl;

	// Transitive-Closure
	Matrix<bool> closure;
	transitive_closure(g,closure);

	TSize len = length(closure, 0);
	for (TSize row=0;row < len;++row) {
		for (TSize col=0;col < len;++col) {
			std::cout << getValue(closure, row*len+col) << ",";
		}
		std::cout << std::endl;
	}
}

//////////////////////////////////////////////////////////////////////////////

void AutomatonTest() {
//____________________________________________________________________________
// Automaton
	typedef char CargoType;
	typedef EdgeList<CargoType> EdgeType;
	typedef VertexDescriptor<Graph<EdgeType> >::Type VertexDescriptorType;
	typedef EdgeDescriptor<Graph<EdgeType> >::Type EdgeDescriptorType;
	typedef Size<Graph<EdgeType> >::Type TSize;

	//Create the graph
	Graph<EdgeType> automaton;

	// Add vertices and edges
	VertexDescriptorType rootVertex = addVertex(automaton); // A = 0
	addVertex(automaton); // B = 1
	addVertex(automaton); // C = 2
	addVertex(automaton); // D = 3
	addVertex(automaton); // E = 4
	addVertex(automaton); // F = 5
	addEdge(automaton,0,1,'2');
	addEdge(automaton,1,0,'1');
	addEdge(automaton,4,0,'6');
	addEdge(automaton,0,3,'7');
	addEdge(automaton,1,1,'3');
	addEdge(automaton,1,2,'4');
	addEdge(automaton,5,1,'8');
	addEdge(automaton,2,5,'5');
	addEdge(automaton,3,4,'2');
	addEdge(automaton,5,3,'7');

	// Vertex names are handled in an external property map
	char names[] = {'A', 'B', 'C', 'D', 'E', 'F'};
	String<char> nameMap;
	initVertexMap(automaton,nameMap, names);

	std::cout << "A simple forward walk through the automaton:" << std::endl;
	std::cout << getProperty(nameMap, rootVertex) << " -7-> ";
	VertexDescriptorType succ = getSuccessorVertex(automaton,rootVertex,'7');
	std::cout << getProperty(nameMap, succ) << " -2-> ";
	succ = getSuccessorVertex(automaton,succ,'2');
	std::cout << getProperty(nameMap, succ) << " -6-> ";
	succ = getSuccessorVertex(automaton,succ,'6');
	std::cout << getProperty(nameMap, succ) << " -2-> ";
	succ = getSuccessorVertex(automaton,succ,'2');
	std::cout << getProperty(nameMap, succ) << std::endl;
	std::cout << "...and now let's walk backwards:" << std::endl;
	VertexDescriptorType pred = succ;
	std::cout << getProperty(nameMap, pred) << " <-8- ";
	pred = getPredecessorVertex(automaton,pred,'8');
	std::cout << getProperty(nameMap, pred) << " <-5- ";
	pred = getPredecessorVertex(automaton,pred,'5');
	std::cout << getProperty(nameMap, pred) << " <-4- ";
	pred = getPredecessorVertex(automaton,pred,'4');
	std::cout << getProperty(nameMap, pred) << " <-3- ";
	pred = getPredecessorVertex(automaton,pred,'3');
	std::cout << getProperty(nameMap, pred) << std::endl;

	// Print automaton
	std::cout << "Automaton - State: (Input, NextState)" << std::endl;
	typedef Iterator<Graph<EdgeType> , VertexIterator<> >::Type TVertexIterator;
	TVertexIterator it(automaton);
	while(!atEnd(it)) {
		VertexDescriptorType source = getValue(it);
		std::cout << getProperty(nameMap, source) << ": ";
		typedef Iterator<Graph<EdgeType>, OutEdgeIterator<> >::Type TOutEdgeIterator;
		TOutEdgeIterator itout(automaton,source);
		for(;!atEnd(itout);goNext(itout)) {
			VertexDescriptorType sink = targetVertex(automaton, getValue(itout));
			std::cout << "(" << getCargo(getValue(itout)) << "/" << getProperty(nameMap, sink) << "), ";
		}
		std::cout << std::endl;
		goNext(it);
	}
}

int main ()
{
//____________________________________________________________________________
// Graph Algorithms
	// Elementary graph algorithms
	std::cout << "===================================" << std::endl;
	std::cout << "----Breadth-Frist Search-----------" << std::endl;
	BreadthFirstSearch();
	std::cout << "===================================" << std::endl;
	std::cout << "----Depth-Frist Search-------------" << std::endl;
	DepthFirstSearch();
	std::cout << "===================================" << std::endl;
	std::cout << "----Topological Sort---------------" << std::endl;
	TopologicalSort();
	std::cout << "===================================" << std::endl;
	std::cout << "----Strongly-Connected-Components--" << std::endl;
	StronglyConnectedComponents();

	// Single-Source shortest paths
	std::cout << "===================================" << std::endl;
	std::cout << "----DAG-Shortest Path--------------" << std::endl;
	DagShortestPath();
	std::cout << "===================================" << std::endl;
	std::cout << "----Bellman-Ford-------------------" << std::endl;
	BellmanFord();
	std::cout << "===================================" << std::endl;
	std::cout << "----Dijkstra-----------------------" << std::endl;
	Dijkstra();
	std::cout << "===================================" << std::endl;
	std::cout << "----Dijkstra (Internal Map)--------" << std::endl;
	DijkstraInternalMap();

	// All-Pairs Shortest paths
	std::cout << "===================================" << std::endl;
	std::cout << "----All-Pairs Shortest Path--------" << std::endl;
	AllPairsShortestPath();
	std::cout << "===================================" << std::endl;
	std::cout << "----Floyd-Warshall-----------------" << std::endl;
	FloydWarshall();
	std::cout << "===================================" << std::endl;
	std::cout << "----Transitive-Closure-------------" << std::endl;
	TransitiveClosure();

	// Automaton
	std::cout << "===================================" << std::endl;
	std::cout << "----Automaton----------------------" << std::endl;
	AutomatonTest();

	return 0;
}
