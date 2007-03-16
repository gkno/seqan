#ifndef SEQAN_HEADER_GRAPH_ALGORITHM_H
#define SEQAN_HEADER_GRAPH_ALGORITHM_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph - Algorithms
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////







//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Elementary graph algorithms
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Breadth-first search
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Function.breadth_first_search:
..cat:Graph
..summary:Implements a breadth-first search on a graph.
..remarks:Breadth-first search computes the distance from source to all reachable
vertices. It also produces a breath-first tree where each node has a predecessor / parent.
..signature:breadth_first_search(g, source, predecessor, distance)
..param.g:In-parameter:A graph.
...type:Class.Graph
..param.source:In-parameter:A vertex descriptor.
...type:Metafunction.VertexDescriptor
...remarks:The breadth-first search is started from this vertex.
..param.predecessor:Out-parameter:A property map.
...remarks:The predecessor map stores implicitly the breadth-first tree.
..param.distance:Out-parameter:A property map.
...remarks:The distance map indicates at what depth a vertex was discovered.
..returns:void.
..see:Function.depth_first_search
*/
template<typename TSpec, typename TVertexDescriptor, typename TPredecessorMap, typename TDistanceMap>
void
breadth_first_search(Graph<TSpec> const& g,
					 TVertexDescriptor const source,
					 TPredecessorMap& predecessor, 
					 TDistanceMap& distance)
{
	typedef typename Iterator<Graph<TSpec>, EdgeIterator<> >::Type TEdgeIterator;
	typedef typename Iterator<Graph<TSpec>, VertexIterator<> >::Type TVertexIterator;
	typedef typename Value<TPredecessorMap>::Type TPredVal;
	typedef typename Value<TDistanceMap>::Type TDistVal;

	// Initialization
	initVertexMap(g,predecessor);
	initVertexMap(g,distance);
	TPredVal nilPred = getNilPredecessor(g);
	TDistVal infDist = getInfinityDistance(distance);
	
	String<bool> tokenMap;
	initVertexMap(g, tokenMap);
	TVertexIterator it(g);
	for(;!atEnd(it);goNext(it)) {
		assignProperty(tokenMap, getValue(it), false);
		assignProperty(distance, getValue(it), infDist);
		assignProperty(predecessor, getValue(it), nilPred);
	}
	assignProperty(tokenMap, source, true);
	assignProperty(distance, source, 0);
	assignProperty(predecessor, source, nilPred);
	std::deque<TVertexDescriptor> queue;
	queue.push_back(source);
	
	// Bfs
	while (!queue.empty()) {
		TVertexDescriptor u = queue.front();
		queue.pop_front();
		typedef typename Iterator<Graph<TSpec>, OutEdgeIterator<> >::Type TOutEdgeIterator;
		TOutEdgeIterator itout(g,u);
		for(;!atEnd(itout);goNext(itout)) {
			TVertexDescriptor v = targetVertex(itout);
			if (getProperty(tokenMap, v) == false) {
				assignProperty(tokenMap, v, true);
				assignProperty(distance, v, getProperty(distance,u) + 1);
				assignProperty(predecessor, v, u);
				queue.push_back(v);
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
// Depth-first search
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TVertexDescriptor, typename TTokenMap, typename TPredecessorMap, typename TDiscoveryTimeMap, typename TFinishingTimeMap, typename TVal>
void
_dfs_visit(Graph<TSpec> const& g,
		   TVertexDescriptor const u,
		   TTokenMap& tokenMap,
		   TPredecessorMap& predecessor,
		   TDiscoveryTimeMap& disc,
		   TFinishingTimeMap& finish,
		   TVal& time)
{
	typedef typename Iterator<Graph<TSpec>, AdjacencyIterator<> >::Type TAdjacencyIterator;

	assignProperty(tokenMap, u, true);
	++time;
	assignProperty(disc, u, time);
	TAdjacencyIterator itad(g,u);
	for(;!atEnd(itad);goNext(itad)) {
		TVertexDescriptor v = getValue(itad);
		if (getProperty(tokenMap, v) == false) {
			assignProperty(predecessor, v, u);
			_dfs_visit(g, v, tokenMap, predecessor, disc, finish, time);
		}
	}
	++time;
	assignProperty(finish, u, time);
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.depth_first_search:
..cat:Graph
..summary:Implements a depth-first search on a graph.
..remarks:In contrast to a breadth-first search the depth-first search is repeated from multiple sources if the graph is not connected.
Hence, depth-first search produces a depth-first forest. To ensure each vertex ends up in exactly one tree we need not just a distance but a
discovery and finishing time.
..signature:depth_first_search(g, predecessor, discovery, finish)
..param.g:In-parameter:A graph.
...type:Class.Graph
..param.predecessor:Out-parameter:A property map.
...remarks:Predecessor subgraph produced by the depth-first search.
..param.discovery:Out-parameter:A property map.
...remarks:The discovery time of a vertex v.
..param.finish:Out-parameter:A property map.
...remarks:The time when v's adjacency list has been fully explored.
..returns:void.
..see:Function.breadth_first_search
*/
template<typename TSpec, typename TPredecessorMap, typename TDiscoveryTimeMap, typename TFinishingTimeMap>
void
depth_first_search(Graph<TSpec> const& g,
				   TPredecessorMap& predecessor,
				   TDiscoveryTimeMap& disc,
				   TFinishingTimeMap& finish)
{
	typedef typename Iterator<Graph<TSpec>, EdgeIterator<> >::Type TEdgeIterator;
	typedef typename Iterator<Graph<TSpec>, VertexIterator<> >::Type TVertexIterator;
	typedef typename VertexDescriptor<Graph<TSpec> >::Type TVertexDescriptor;
	typedef typename Value<TPredecessorMap>::Type TPredVal;

	// Initialization
	initVertexMap(g,predecessor);
	initVertexMap(g,disc);
	initVertexMap(g,finish);
	TPredVal nilPred = getNilPredecessor(g);
		
	String<bool> tokenMap;
	initVertexMap(g, tokenMap);
	TVertexIterator it(g);
	for(;!atEnd(it);goNext(it)) {
		assignProperty(tokenMap, getValue(it), false);
		assignProperty(predecessor, getValue(it), nilPred);
	}

	unsigned int time = 0;

	goBegin(it);
	for(;!atEnd(it);goNext(it)) {
		TVertexDescriptor u = getValue(it);
		if (getProperty(tokenMap, u) == false) {
			_dfs_visit(g, u, tokenMap, predecessor, disc, finish, time);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// Topological sort
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.topological_sort:
..cat:Graph
..summary:Performs a topological sort on a directed acyclic graph (DAG).
..remarks:A topological sort is a linear ordering of all its vertices such that if the graph contains an edge (u,v) then u appears before v in the ordering.
..signature:topological_sort(g, topSort)
..param.g:In-parameter:A directed acyclic graph.
...type:Spec.Directed graph
..param.topSort:Out-parameter:A linear ordering of the vertices.
...type:Class.String
..returns:void.
*/
template<typename TSpec, typename TVertexDescriptor>
void
topological_sort(Graph<TSpec> const& g,
				 String<TVertexDescriptor>& topSort)
{
	// Initialization
	String<unsigned int> predMap;
	String<unsigned int> discoveryTimeMap;
	String<unsigned int> finishingTimeMap;
	
	// Dfs
	depth_first_search(g, predMap, discoveryTimeMap, finishingTimeMap);

	// Order vertices
	typedef ::std::pair<unsigned int, TVertexDescriptor> TTimeVertexPair;
	std::priority_queue<TTimeVertexPair> q;
	typedef typename Iterator<Graph<TSpec>, VertexIterator<> >::Type TVertexIterator;
	TVertexIterator it(g);
	for(;!atEnd(it);++it) {
		q.push(std::make_pair(getProperty(finishingTimeMap, getValue(it)), getValue(it)));
	}

	// Create topological order
	resize(topSort,numVertices(g));
	unsigned int count=0;
	while(!q.empty()) {
		assignValue(topSort, count, q.top().second);
		q.pop();
		++count;
	}
}

//////////////////////////////////////////////////////////////////////////////
// Strongly connected components
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Function.strongly_connected_components:
..cat:Graph
..summary:Decomposes a directed graph into its strongly connected components.
..signature:strongly_connected_components(g, components)
..param.g:In-parameter:A directed graph.
...type:Spec.Directed graph
..param.components:Out-parameter:A property map.
...remarks:Each vertex is mapped to a component id. If two vertices share the same id they are in the same component.
..returns:void.
*/

template<typename TSpec, typename TComponents>
void
strongly_connected_components(Graph<TSpec> const& g_source,
							  TComponents& components)
{
	// Initialization
	typedef typename VertexDescriptor<Graph<TSpec> >::Type TVertexDescriptor;
	typedef typename Iterator<Graph<TSpec>, EdgeIterator<> >::Type TEdgeIterator;
	typedef typename Iterator<Graph<TSpec>, VertexIterator<> >::Type TVertexIterator;
	typedef typename Value<TComponents>::Type TCompVal;
	initVertexMap(g_source,components);
	String<unsigned int> predMap;
	String<unsigned int> discoveryTimeMap;
	String<unsigned int> finishingTimeMap;
	
	// Dfs
	depth_first_search(g_source, predMap, discoveryTimeMap, finishingTimeMap);

	Graph<TSpec> g;
	transpose(g_source, g);

	// Second Dfs
	String<unsigned int> predecessor;
	String<unsigned int> disc;
	String<unsigned int> finish;
	initVertexMap(g,predecessor);
	initVertexMap(g,disc);
	initVertexMap(g,finish);
	TCompVal nilPred = getNilPredecessor(g);
	String<bool> tokenMap;
	initVertexMap(g, tokenMap);
	TVertexIterator it(g);
	for(;!atEnd(it);goNext(it)) {
		assignProperty(components, getValue(it), nilPred);
		assignProperty(tokenMap, getValue(it), false);
		assignProperty(predecessor, getValue(it), nilPred);
	}

	// Order vertices
	typedef ::std::pair<unsigned int, TVertexDescriptor> TTimeVertexPair;
	std::priority_queue<TTimeVertexPair> q;
	goBegin(it);
	for(;!atEnd(it);++it) {
		q.push(std::make_pair(getProperty(finishingTimeMap, getValue(it)), getValue(it)));
	}

	unsigned int time = 0;
	unsigned int label = 0;
	while(!q.empty()) {
		TVertexDescriptor u = q.top().second;
		q.pop();
		if (getProperty(tokenMap, u) == false) {
			_dfs_visit(g, u, tokenMap, predecessor, disc, finish, time);
			TVertexIterator it_label(g);
			for(;!atEnd(it_label);goNext(it_label)) {
				if ((getProperty(tokenMap, getValue(it_label)) == true) &&
					(getProperty(components, getValue(it_label)) == nilPred)) {
					assignProperty(components, getValue(it_label), label);
				}
			}
			++label;
		}
	}
}














//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Single-Source Shortest Paths
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TPredecessorMap, typename TVertexDescriptor>
inline void
_print_path(Graph<TSpec> const& g,
			TPredecessorMap const& predecessor,
			TVertexDescriptor const source,
			TVertexDescriptor const v)
{
	if (source == v) {
		std::cout << source;
	} else if (getProperty(predecessor, v) == getNilPredecessor(g)) {
		std::cout << "No path from " << source << " to " << v << " exists.";
	} else {
		_print_path(g,predecessor, source, getProperty(predecessor, v));
		std::cout << "," << v;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
inline void 
_initialize_single_source(Graph<TSpec> const& g,
						  TVertexDescriptor const source,
						  TWeightMap const& weight,
						  TPredecessorMap& predecessor, 
						  TDistanceMap& distance)
{
	typedef typename Iterator<Graph<TSpec>, VertexIterator<> >::Type TVertexIterator;
	typedef typename Value<TPredecessorMap>::Type TPredVal;
	typedef typename Value<TWeightMap>::Type TDistVal;
	TPredVal nilPred = getNilPredecessor(g);
	TDistVal infDist = getInfinityDistance(weight);
	
	TVertexIterator it(g);
	while(!atEnd(it)) {
		assignProperty(distance, getValue(it), infDist);
		assignProperty(predecessor, getValue(it), nilPred);
		goNext(it);
	}
	assignProperty(distance, source, 0);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap, typename TVertexDescriptor, typename TEdgeDescriptor>
inline void 
_relax(Graph<TSpec> const& g,
	    TWeightMap const& weight,
		TPredecessorMap& predecessor, 
		TDistanceMap& distance,
		TVertexDescriptor const u,
		TEdgeDescriptor const e)
{
	TVertexDescriptor v = targetVertex(g,e);
	if (getProperty(distance, v) > getProperty(distance,u) + getProperty(weight,e)) {
		assignProperty(distance, v, getProperty(distance,u) + getProperty(weight,e));
		assignProperty(predecessor, v, u);
	}
}



//////////////////////////////////////////////////////////////////////////////
// DAG Shortest Path
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Function.dag_shortest_path:
..cat:Graph
..summary:Computes shortest paths from a single source in a directed acyclic graph (DAG).
..signature:dag_shortest_path(g, source, weight, predecessor, distance)
..param.g:In-parameter:A directed acyclic graph.
...type:Spec.Directed graph
..param.source:In-parameter:A source vertex.
...type:Metafunction.VertexDescriptor
..param.weight:In-parameter:A weight map.
...remarks:In a directed acyclic graph edge weights can be negative because no cycles do exist.
..param.predecessor:Out-parameter:A property map.
...remarks:A property map that represents predecessor relationships among vertices. It determines a shortest-paths tree.
..param.distance:Out-parameter:A property map.
...remarks:Indicates for each vertex the distance from the source.
..returns:void.
..see:Function.bellman_ford_algorithm
..see:Function.dijkstra
*/
template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
void
dag_shortest_path(Graph<TSpec> const& g,
				  TVertexDescriptor const source,
				  TWeightMap const& weight,
				  TPredecessorMap& predecessor,
				  TDistanceMap& distance)
{
	typedef typename EdgeDescriptor<Graph<TSpec> >::Type TEdgeDescriptor;
	typedef typename Iterator<Graph<TSpec>, EdgeIterator<> >::Type TEdgeIterator;
	typedef typename Iterator<Graph<TSpec>, OutEdgeIterator<> >::Type TOutEdgeIterator;
	typedef typename Iterator<String<TVertexDescriptor> >::Type TStringIterator;
	
	// Initialization
	initVertexMap(g,predecessor);
	initVertexMap(g,distance);

	// Topological sort
	String<TVertexDescriptor> order;
	topological_sort(g, order);

	_initialize_single_source(g, source, weight, predecessor, distance);

	//DAG Shortest Paths
	TStringIterator it = begin(order);
	while(!atEnd(it)) {
		TOutEdgeIterator itout(g, getValue(it));
		for(;!atEnd(itout);++itout) {
			_relax(g,weight,predecessor, distance, getValue(it), getValue(itout));
		}
		goNext(it);
	}
}


//////////////////////////////////////////////////////////////////////////////
// Bellman-Ford
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Function.bellman_ford_algorithm:
..cat:Graph
..summary:Computes shortest paths from a single source in a graph.
..remarks:Edge weights may be negative in the Bellman-Ford algorithm.
The out parameters are only valid if the algorithm returns true.
..signature:bellman_ford_algorithm(g, source, weight, predecessor, distance)
..param.g:In-parameter:A graph.
...type:Class.Graph
..param.source:In-parameter:A source vertex.
...type:Metafunction.VertexDescriptor
..param.weight:In-parameter:A weight map.
...remarks:A property map with edge weights. Edge weights may be negative.
..param.predecessor:Out-parameter:A property map.
...remarks:A property map that represents predecessor relationships among vertices. It determines a shortest-paths tree.
..param.distance:Out-parameter:A property map.
...remarks:Indicates for each vertex the distance from the source.
..returns:True if the graph has no negative weight cycles, false otherwise.
..see:Function.dag_shortest_path
..see:Function.dijkstra
*/
template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
bool 
bellman_ford_algorithm(Graph<TSpec> const& g,
					   TVertexDescriptor const source,
					   TWeightMap const& weight,
					   TPredecessorMap& predecessor, 
					   TDistanceMap& distance)
{
	// Initialization
	typedef typename EdgeDescriptor<Graph<TSpec> >::Type TEdgeDescriptor;
	typedef typename Iterator<Graph<TSpec>, VertexIterator<> >::Type TVertexIterator;
	typedef typename Iterator<Graph<TSpec>, OutEdgeIterator<> >::Type TOutEdgeIterator;
	initVertexMap(g,predecessor);
	initVertexMap(g,distance);
	_initialize_single_source(g, source, weight, predecessor, distance);

	// Run Bellman-Ford
	for(unsigned int i=0; i<numVertices(g) - 1; ++i) {
		TVertexIterator it(g);
		for(;!atEnd(it);goNext(it)) {
			TVertexDescriptor u = getValue(it);
			TOutEdgeIterator itout(g, u);
			for(;!atEnd(itout);++itout) {
				_relax(g,weight,predecessor, distance, u, getValue(itout));
			}
		}
	}

	TVertexIterator it(g);
	for(;!atEnd(it);goNext(it)) {
		TVertexDescriptor u = getValue(it);
		TOutEdgeIterator itout(g, u);
		for(;!atEnd(itout);++itout) {
			TVertexDescriptor v = targetVertex(g, getValue(itout));
			if (getProperty(distance, v) > getProperty(distance,u) + getProperty(weight,getValue(itout))) {
				return false;
			}	
		}
	}
	return true;
}


//////////////////////////////////////////////////////////////////////////////
// Dijkstra
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Function.dijkstra:
..cat:Graph
..summary:Computes shortest paths from a single source in a graph.
..remarks:Edge weights have to be nonnegative.
..signature:dijkstra(g, source, weight, predecessor, distance)
..param.g:In-parameter:A graph.
...type:Class.Graph
..param.source:In-parameter:A source vertex.
...type:Metafunction.VertexDescriptor
..param.weight:In-parameter:A weight map.
...remarks:A property map with edge weights. Edge weights have to be nonnegative.
..param.predecessor:Out-parameter:A property map.
...remarks:A property map that represents predecessor relationships among vertices. It determines a shortest-paths tree.
..param.distance:Out-parameter:A property map.
...remarks:Indicates for each vertex the distance from the source.
..returns:void
..see:Function.dag_shortest_path
..see:Function.bellman_ford_algorithm
*/
template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
void 
dijkstra(Graph<TSpec> const& g,
		 TVertexDescriptor const source,
		 TWeightMap const& weight,
		 TPredecessorMap& predecessor, 
		 TDistanceMap& distance)
{
	typedef typename Value<TDistanceMap>::Type TDistVal;
	typedef typename Iterator<Graph<TSpec>, VertexIterator<> >::Type TVertexIterator;
	typedef typename Iterator<Graph<TSpec>, OutEdgeIterator<> >::Type TOutEdgeIterator;
	
	// Initialization
	initVertexMap(g,predecessor);
	initVertexMap(g,distance);

	_initialize_single_source(g, source, weight, predecessor, distance);
	
	String<bool> setS;
	initVertexMap(g, setS);
	TVertexIterator it(g);
	for(;!atEnd(it);++it) {
		assignProperty(setS, getValue(it), false);
	}
	TDistVal infDist = getInfinityDistance(weight);
	TVertexDescriptor nilVertex = getNilPredecessor(g);

	// Run Dijkstra
	unsigned int count = numVertices(g);
	while (count > 0) {
		// Extract min
		TDistVal min = infDist;
		TVertexDescriptor u = nilVertex;
		TVertexIterator it_find(g);
		for(;!atEnd(it_find);++it_find) {
			if(getProperty(setS,getValue(it_find))==true) continue;
			if ((u == nilVertex) ||
				(getProperty(distance,getValue(it_find))<getProperty(distance,u))) {
					u = getValue(it_find);
					min = getProperty(distance,getValue(it_find));
			}
		}
		assignProperty(setS, u, true);
		TOutEdgeIterator itout(g, u);
		for(;!atEnd(itout);++itout) {
			_relax(g,weight,predecessor, distance, u, getValue(itout));
		}
		--count;
	}
}








//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// All-Pairs shortest paths
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TPredecessor, typename TVertexDescriptor>
inline void
_print_all_pairs_shortest_path(Graph<TSpec> const& g,
							   TPredecessor& predecessor, 
							   TVertexDescriptor const i,
							   TVertexDescriptor const j)
{
	typedef typename Size<TPredecessor>::Type TSize;
	TSize len = getIdUpperBound(g.data_id_managerV);
	if (i==j) {
		std::cout << i;
	} else if (getValue(predecessor, i*len+j) == getNilPredecessor(g)) {
		std::cout << "No path from " << i << " to " << j << " exists.";
	} else {
		_print_all_pairs_shortest_path(g,predecessor, i, (TVertexDescriptor) getValue(predecessor, i*len+j));
		std::cout << "," << j;
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TWeightMap, typename TMatrix, typename TPredecessor>
void 
_initialize_all_pairs(Graph<TSpec> const& g,
						TWeightMap const& weight,
						TMatrix& matrix,
						TPredecessor& predecessor)
{
	typedef typename VertexDescriptor<Graph<TSpec> >::Type TVertexDescriptor;
	typedef typename Iterator<Graph<TSpec>, VertexIterator<> >::Type TVertexIterator;
	typedef typename Iterator<Graph<TSpec>, OutEdgeIterator<> >::Type TOutEdgeIterator;
	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Value<TWeightMap>::Type TWeightVal;
	typedef typename Value<TPredecessor>::Type TPredVal;
	
	// Create adjacency-like matrix
	TSize len = getIdUpperBound(g.data_id_managerV);
	setDimension(matrix, 2);
	setLength(matrix, 0, len);
	setLength(matrix, 1, len);
	resize(matrix);
	setDimension(predecessor, 2);
	setLength(predecessor, 0, len);
	setLength(predecessor, 1, len);
	resize(predecessor);
	TWeightVal infWeight = getInfinityDistance(weight);
	TPredVal nilPred = getNilPredecessor(g);
	for (TSize row=0;row < len;++row) {
		for (TSize col=0;col < len;++col) {
			if (row != col) assignValue(matrix, row*len + col, infWeight);
			else assignValue(matrix, row*len + col, 0);
			assignValue(predecessor, row*len + col, nilPred);
		}
	}

	// Include edge weights and initial predecessors
	TVertexIterator it(g);
	for(;!atEnd(it);goNext(it)) {
		TVertexDescriptor u = getValue(it);
		TOutEdgeIterator itout(g, u);
		for(;!atEnd(itout);++itout) {
			TVertexDescriptor v = targetVertex(g,getValue(itout));
			assignValue(matrix, u*len + v, getProperty(weight, getValue(itout)));
			assignValue(predecessor, u*len + v, u);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TMatrix, typename TPredecessor, typename TInfDist>
void 
_extend_shortest_paths(TMatrix& local,
					   TMatrix& w,
					   TPredecessor& predecessor,
					   TInfDist const infDist)
{
	typedef typename Value<TMatrix>::Type TMatrixVal;
	typedef typename Value<TPredecessor>::Type TPredVal;
	typedef typename Size<TMatrix>::Type TSize;
	TMatrix oldLocal = local;
	TPredecessor oldPredecessor = predecessor;
	TSize len = length(oldLocal, 0);
	for(TSize i = 0; i<len;++i) {
		for(TSize j = 0; j<len;++j) {
			if (i==j) continue;
			assignValue(local, i*len+j,infDist);
			TPredVal ind;
			for(TSize k = 0; k<len;++k) {
				TMatrixVal min1 = getValue(local, i*len+j);
				TMatrixVal min2 = getValue(oldLocal, i*len+k) + getValue(w, k*len + j);
				if (min2 < min1) {
					assignValue(local, i*len+j,min2);
					ind = k;
				}
			}
			if (getValue(oldLocal, i*len+j) > getValue(local, i*len+j)) {
				assignValue(predecessor, i*len+j,ind);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// All-Pairs shortest path
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
/**
.Function.all_pairs_shortest_path:
..cat:Graph
..summary:Finds shortest paths between all pairs of vertices in a graph.
..signature:all_pairs_shortest_path(g, weight, distance, predecessor)
..param.g:In-parameter:A graph.
...type:Class.Graph
..param.weight:In-parameter:A weight map.
...remarks:A property map with edge weights. Edge weights may be negative.
..param.distance:Out-parameter:A matrix with distances.
...type:Class.Matrix
...remarks:Entry (i,j) in this matrix indicates the distance from vertex i to vertex j.
..param.predecessor:Out-parameter:A matrix with predecessors.
...type:Class.Matrix
...remarks:Entry (i,j) in this matrix indicates the predecessor of j on a shortest path from vertex i to vertex j.
You can use _print_all_pairs_shortest_path(g, predecessor, i, j) to print the shortest path from i to j.
..returns:void
..see:Function.floyd_warshall
*/
template<typename TSpec, typename TWeightMap, typename TMatrix, typename TPredecessor>
void 
all_pairs_shortest_path(Graph<TSpec> const& g,
						TWeightMap const& weight,
						TMatrix& distMatrix,
						TPredecessor& predecessor)
{
	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Value<TWeightMap>::Type TWeightVal;
	TWeightVal infWeight = getInfinityDistance(weight);

	// Initialize first distance matrix
	_initialize_all_pairs(g,weight,distMatrix,predecessor);

	TSize len = length(distMatrix, 0);
	TMatrix local = distMatrix;
	for(TSize m=2;m<len;++m) {
		_extend_shortest_paths(local,distMatrix,predecessor, infWeight);
	}
	distMatrix = local;
}


//////////////////////////////////////////////////////////////////////////////
// Floyd-Warshall
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
/**
.Function.floyd_warshall:
..cat:Graph
..summary:Finds shortest paths between all pairs of vertices in a graph.
..signature:floyd_warshall(g, weight, distance, predecessor)
..remarks:The graph must be free of negative-weight cycles.
..param.g:In-parameter:A graph.
...type:Class.Graph
..param.weight:In-parameter:A weight map.
...remarks:A property map with edge weights. Edge weights may be negative.
..param.distance:Out-parameter:A matrix with distances.
...type:Class.Matrix
...remarks:Entry (i,j) in this matrix indicates the distance from vertex i to vertex j.
..param.predecessor:Out-parameter:A matrix with predecessors.
...type:Class.Matrix
...remarks:Entry (i,j) in this matrix indicates the predecessor of j on a shortest path from vertex i to vertex j.
You can use _print_all_pairs_shortest_path(g, predecessor, i, j) to print the shortest path from i to j.
..returns:void
..see:Function.all_pairs_shortest_path
*/
template<typename TSpec, typename TWeightMap, typename TMatrix, typename TPredecessor>
void 
floyd_warshall(Graph<TSpec> const& g,
			   TWeightMap const& weight,
			   TMatrix& distMatrix,
			   TPredecessor& predecessor)
{
	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Value<TMatrix>::Type TMatrixVal;

	// Initialize first distance matrix
	_initialize_all_pairs(g,weight,distMatrix,predecessor);

	// Floyd-Warshall
	TSize len = length(distMatrix, 0);
	TMatrix local = distMatrix;
	for(TSize k=0;k<len;++k) {
		for(TSize i=0;i<len;++i) {
			for(TSize j=0;j<len;++j) {
				TMatrixVal min1 = getValue(distMatrix, i*len+j);
				TMatrixVal min2 = getValue(distMatrix, i*len+k) + getValue(distMatrix, k*len + j);
				if (min2 < min1) {
					assignValue(local, i*len+j,min2);
					assignValue(predecessor, i*len+j,getValue(predecessor, k*len+j));
				} else {
					assignValue(local, i*len+j,min1);
					assignValue(predecessor, i*len+j, getValue(predecessor, i*len+j));
				}
			}
		}
		distMatrix=local;
	}
}

//////////////////////////////////////////////////////////////////////////////
// Transitive Closure
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
/**
.Function.transitive_closure:
..cat:Graph
..summary:Determines whether there is a path between any two given vertices or not.
..signature:transitive_closure(g, closure)
..param.g:In-parameter:A graph.
...type:Class.Graph
..param.closure:Out-parameter:A matrix which indicates the closure.
...type:Class.Matrix
...remarks:Entry (i,j) in this matrix indicates whether there is a path from i to j in the graph or not.
..returns:void
*/
template<typename TSpec, typename TMatrix>
void 
transitive_closure(Graph<TSpec> const& g,
				   TMatrix& closure)
{
	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Value<TMatrix>::Type TMatrixVal;

	// Initialize first closure matrix
	getAdjacencyMatrix(g,closure);
	TSize len = length(closure, 0);
	for (TSize diag=0;diag < len;++diag) assignValue(closure, diag*len+diag,1);

	// Transitive Closure
	TMatrix local = closure;
	for (TSize k=0;k<len;++k) {
		for(TSize i=0;i<len;++i) {
			for(TSize j=0;j<len;++j) {
				TMatrixVal t_ij = getValue(closure, i*len+j);
				TMatrixVal t_ik = getValue(closure, i*len+k);
				TMatrixVal t_kj = getValue(closure, k*len+j);
				assignValue(local, i*len+j, t_ij | (t_ik & t_kj));
			}
		}
		closure = local;
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
