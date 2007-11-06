 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id$
 ==========================================================================*/

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
	SEQAN_CHECKPOINT
	typedef typename Iterator<Graph<TSpec>, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<Graph<TSpec>, VertexIterator>::Type TVertexIterator;
	typedef typename Value<TPredecessorMap>::Type TPredVal;
	typedef typename Value<TDistanceMap>::Type TDistVal;

	// Initialization
	resizeVertexMap(g,predecessor);
	resizeVertexMap(g,distance);
	TPredVal nilPred = getNilPredecessor(g);
	TDistVal infDist = getInfinityDistance(distance);
	
	String<bool> tokenMap;
	resizeVertexMap(g, tokenMap);
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
		typedef typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type TOutEdgeIterator;
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
	SEQAN_CHECKPOINT

	typedef typename Iterator<Graph<TSpec>, AdjacencyIterator>::Type TAdjacencyIterator;

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
	SEQAN_CHECKPOINT

	typedef typename Iterator<Graph<TSpec>, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<Graph<TSpec>, VertexIterator>::Type TVertexIterator;
	typedef typename VertexDescriptor<Graph<TSpec> >::Type TVertexDescriptor;
	typedef typename Value<TPredecessorMap>::Type TPredVal;

	// Initialization
	resizeVertexMap(g,predecessor);
	resizeVertexMap(g,disc);
	resizeVertexMap(g,finish);
	TPredVal nilPred = getNilPredecessor(g);
		
	String<bool> tokenMap;
	resizeVertexMap(g, tokenMap);
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
	SEQAN_CHECKPOINT

	// Initialization
	String<unsigned int> predMap;
	String<unsigned int> discoveryTimeMap;
	String<unsigned int> finishingTimeMap;
	
	// Dfs
	depth_first_search(g, predMap, discoveryTimeMap, finishingTimeMap);

	// Order vertices
	typedef ::std::pair<unsigned int, TVertexDescriptor> TTimeVertexPair;
	std::priority_queue<TTimeVertexPair> q;
	typedef typename Iterator<Graph<TSpec>, VertexIterator>::Type TVertexIterator;
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
	SEQAN_CHECKPOINT

	// Initialization
	typedef typename VertexDescriptor<Graph<TSpec> >::Type TVertexDescriptor;
	typedef typename Iterator<Graph<TSpec>, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<Graph<TSpec>, VertexIterator>::Type TVertexIterator;
	typedef typename Value<TComponents>::Type TCompVal;
	resizeVertexMap(g_source,components);
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
	resizeVertexMap(g,predecessor);
	resizeVertexMap(g,disc);
	resizeVertexMap(g,finish);
	TCompVal nilPred = getNilPredecessor(g);
	String<bool> tokenMap;
	resizeVertexMap(g, tokenMap);
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
// Minimum Spanning Trees
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Prim's algorithm
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
.Function.prims_algorithm:
..cat:Graph
..summary:Computes a minimum spanning tree on a graph.
..signature:prims_algorithm(g, source, weight, predecessor)
..param.g:In-parameter:A graph.
...type:Class.Graph
..param.source:In-parameter:A source vertex.
...type:Metafunction.VertexDescriptor
..param.weight:In-parameter:Edge weights.
..param.predecessor:Out-parameter:A property map.
...remarks:A property map that represents predecessor relationships among vertices. It determines a minimum spanning tree.
..returns:void.
..see:Function.kruskals_algorithm
*/
template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap>
void
prims_algorithm(Graph<TSpec> const& g,
				TVertexDescriptor const source,
				TWeightMap const& weight,
				TPredecessorMap& predecessor)
{
	SEQAN_CHECKPOINT
	typedef Graph<TSpec> TGraph;
	typedef typename Iterator<Graph<TSpec>, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename Value<TPredecessorMap>::Type TPred;
	typedef typename Value<TWeightMap>::Type TWeight;

	typedef ::std::pair<TWeight, TVertexDescriptor> TWeightVertexPair;
	std::priority_queue<TWeightVertexPair, std::vector<TWeightVertexPair>, std::greater<TWeightVertexPair> > q;
	
	// Initialization
	String<bool> tokenMap;
	String<TWeight> key;
	TPred nilPred = getNilPredecessor(g);
	TWeight infWeight = getInfinityDistance(weight);
	resizeVertexMap(g,predecessor);
	resizeVertexMap(g,tokenMap);
	resizeVertexMap(g,key);

	TVertexIterator it(g);
	while(!atEnd(it)) {
		TVertexDescriptor u = getValue(it);
		if (u == source) q.push(std::make_pair(0, u));
		assignProperty(predecessor, u, nilPred);
		assignProperty(key, u, infWeight);
		assignProperty(tokenMap, u, false);
		goNext(it);
	}

	assignProperty(key, source, 0);
	while(!q.empty()) {
		TVertexDescriptor u = q.top().second;
		q.pop();
		if (getProperty(tokenMap, u)) continue;
		assignProperty(tokenMap, u, true);
		TOutEdgeIterator itOut(g,u);
		while(!atEnd(itOut)) {
			TVertexDescriptor v = targetVertex(itOut);
			TWeight w = getProperty(weight, getValue(itOut));
			if ((!getProperty(tokenMap, v)) &&
				(w < getProperty(key, v))) {
					assignProperty(predecessor, v, u);
					assignProperty(key, v, w);
					q.push(std::make_pair(w, v));
			}
			goNext(itOut);
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
// Kruskal's algorithm
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
.Function.kruskals_algorithm:
..cat:Graph
..summary:Computes a minimum spanning tree on a graph.
..signature:kruskals_algorithm(g, source, weight, edges)
..param.g:In-parameter:A graph.
...type:Class.Graph
..param.source:In-parameter:A source vertex.
...type:Metafunction.VertexDescriptor
..param.weight:In-parameter:Edge weights.
..param.edges:Out-parameter:Array of vertex descriptors.
...remarks:Array or string where two consecutive entries are an edge.
..returns:void.
..see:Function.prims_algorithm
*/
template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TEdges>
void
kruskals_algorithm(Graph<TSpec> const& g,
				   TVertexDescriptor const,
				   TWeightMap const& weight,
				   TEdges& edges)
{
	SEQAN_CHECKPOINT
	typedef Graph<TSpec> TGraph;
	typedef typename Iterator<Graph<TSpec>, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Value<TWeightMap>::Type TWeight;

	typedef ::std::pair<TWeight, std::pair<TVertexDescriptor, TVertexDescriptor> > TWeightEdgePair;
	std::priority_queue<TWeightEdgePair, std::vector<TWeightEdgePair>, std::greater<TWeightEdgePair> > q;


	resize(edges, 2 * (numVertices(g) - 1));
	String<String<TVertexDescriptor> > set;
	String<TVertexDescriptor> id;
	resizeVertexMap(g, set);
	resizeVertexMap(g, id);
	
	// Make the sets
	TVertexIterator it(g);
	while(!atEnd(it)) {
		TVertexDescriptor v = getValue(it);
		appendValue(property(set, v), v);
		assignProperty(id, v, v);
		goNext(it);
	}

	// Sort the edges
	TEdgeIterator itE(g);
	while(!atEnd(itE)) {
		TVertexDescriptor x = sourceVertex(itE);
		TVertexDescriptor y = targetVertex(itE);
		TWeight w = getProperty(weight, getValue(itE));
		q.push(std::make_pair(w, std::make_pair(x,y)));
		goNext(itE);
	}

	// Process each edge
	unsigned int index = 0;
	while(!q.empty()) {
		TVertexDescriptor x = q.top().second.first;
		TVertexDescriptor y = q.top().second.second;
		q.pop();
		if (getProperty(id, x) != getProperty(id,y)) {
			TVertexDescriptor owner = getProperty(id, x);
			assignValue(edges, index, x);
			assignValue(edges, index+1, y);
			index = index + 2;
			typedef typename Iterator<String<TVertexDescriptor> >::Type TStrIterator;
			TStrIterator strIt = begin(property(set,getProperty(id, y)));
			TStrIterator strItEnd = end(property(set,getProperty(id, y)));
			for(;strIt != strItEnd;goNext(strIt)) {				
				TVertexDescriptor setMember = getValue(strIt);
				appendValue(property(set, owner), setMember);
				assignProperty(id, setMember, owner);
			}
		
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

template<typename TSpec, typename TPredecessorMap, typename TVertexDescriptor, typename TNameMap>
inline void
_print_path(Graph<TSpec> const& g,
			TPredecessorMap const& predecessor,
			TVertexDescriptor const source,
			TVertexDescriptor const v,
			TNameMap const& nameMap)
{
	if (source == v) {
		std::cout << getProperty(nameMap, source);
	} else if (getProperty(predecessor, v) == getNilPredecessor(g)) {
		std::cout << "No path from " << getProperty(nameMap, source) << " to " << getProperty(nameMap, v) << " exists.";
	} else {
		_print_path(g,predecessor, source, getProperty(predecessor, v), nameMap);
		std::cout << "," << getProperty(nameMap, v);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TPredecessorMap, typename TVertexDescriptor1, typename TVertexDescriptor2>
inline void
_print_path(Graph<TSpec> const& g,
			TPredecessorMap const& predecessor,
			TVertexDescriptor1 const source,
			TVertexDescriptor2 const v)
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
	typedef typename Iterator<Graph<TSpec>, VertexIterator>::Type TVertexIterator;
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
	SEQAN_CHECKPOINT
	typedef typename EdgeDescriptor<Graph<TSpec> >::Type TEdgeDescriptor;
	typedef typename Iterator<Graph<TSpec>, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename Iterator<String<TVertexDescriptor>, Rooted>::Type TStringIterator;
	
	// Initialization
	resizeVertexMap(g,predecessor);
	resizeVertexMap(g,distance);

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
	SEQAN_CHECKPOINT
	// Initialization
	typedef typename EdgeDescriptor<Graph<TSpec> >::Type TEdgeDescriptor;
	typedef typename Iterator<Graph<TSpec>, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type TOutEdgeIterator;
	resizeVertexMap(g,predecessor);
	resizeVertexMap(g,distance);
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
	SEQAN_CHECKPOINT
	typedef typename Value<TDistanceMap>::Type TDistVal;
	typedef typename Iterator<Graph<TSpec>, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type TOutEdgeIterator;
	
	// Initialization
	resizeVertexMap(g,predecessor);
	resizeVertexMap(g,distance);

	_initialize_single_source(g, source, weight, predecessor, distance);
	
	String<bool> setS;
	resizeVertexMap(g, setS);
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
	typedef typename Iterator<Graph<TSpec>, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type TOutEdgeIterator;
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
	SEQAN_CHECKPOINT
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
	SEQAN_CHECKPOINT
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
	SEQAN_CHECKPOINT
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



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Maximum Flow
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
template<typename TSpec, typename TCapMap, typename TFlowMap, typename TResidualGraph>
void
_build_residual_graph(Graph<TSpec> const& g,
					  TCapMap const& capacity,
					  TFlowMap const& flow,
					  TResidualGraph& rG)
{
	SEQAN_CHECKPOINT
	typedef Graph<TSpec> TGraph;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Value<TFlowMap>::Type TFlow;
	typedef typename Value<TCapMap>::Type TCap;

	clear(rG);
	TVertexIterator itV(g);
	for(;!atEnd(itV);goNext(itV)) {
		_createVertices(rG, getValue(itV));
	}

	TEdgeIterator itE(g);
	for(;!atEnd(itE);goNext(itE)) {
		typedef typename EdgeDescriptor<TResidualGraph>::Type TEdgeDescriptor;
		TFlow f = getProperty(flow, getValue(itE));
		TCap cap = getProperty(capacity, getValue(itE));
		if (f > 0) {
			TEdgeDescriptor e_rG = findEdge(rG, targetVertex(itE), sourceVertex(itE));
			if (e_rG == 0) addEdge(rG, targetVertex(itE), sourceVertex(itE), f);
			else cargo(e_rG) += f;
		}
		if (f < cap) {
			TEdgeDescriptor e_rG = findEdge(rG, sourceVertex(itE), targetVertex(itE));
			if (e_rG == 0) addEdge(rG, sourceVertex(itE), targetVertex(itE), cap - f);
			else cargo(e_rG) += cap - f;			
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TPredecessorMap, typename TVertexDescriptor>
inline unsigned int
_get_minimum_aug(Graph<TSpec> const& rG,
				 TPredecessorMap& predecessor,
				 TVertexDescriptor const source,
				 TVertexDescriptor sink)
{
	typedef unsigned int TFlow;
	typedef typename Iterator<String<TVertexDescriptor>, Rooted>::Type TIterator;
	
	// Build secondary predecessor map just containing the path
	TVertexDescriptor nilPred = getNilPredecessor(rG);
	String<TVertexDescriptor> predMap;
	resizeVertexMap(rG, predMap);
	TIterator it = begin(predMap);
	for(;!atEnd(it);goNext(it)) {
		*it = nilPred;
	}

	// Find minimum flow
	TVertexDescriptor pred = getProperty(predecessor, sink);
	TFlow f = getCargo(findEdge(rG, pred,sink));
	assignProperty(predMap, sink, pred);
	while(pred != source) {
		sink = pred;
		pred = getProperty(predecessor, sink);
		TFlow f2 = getCargo(findEdge(rG, pred,sink));
		assignProperty(predMap, sink, pred);
		if (f2 < f) f = f2;
	}

	// Just return the augmenting path
	predecessor = predMap;
	return f;
}

//////////////////////////////////////////////////////////////////////////////
// Ford Fulkerson
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
.Function.ford_fulkerson:
..cat:Graph
..summary:Computes a maximum flow in a directed graph.
..signature:ford_fulkerson(g, source, sink, capacity, flow)
..param.g:In-parameter:A graph.
...type:Spec.Directed graph
..param.source:In-parameter:A source vertex.
...type:Metafunction.VertexDescriptor
..param.sink:In-parameter:A sink vertex.
...type:Metafunction.VertexDescriptor
..param.capacity:In-parameter:A property map of edge capacities.
..param.flow:Out-parameter:A property map with the flow of each edge.
..returns:The value of the flow.
*/
template<typename TSpec, typename TVertexDescriptor, typename TCapMap, typename TFlowMap>
typename Value<TFlowMap>::Type
ford_fulkerson(Graph<TSpec> const& g,
			   TVertexDescriptor const source,
			   TVertexDescriptor const sink,
			   TCapMap const& capacity,
			   TFlowMap& flow)
{
	SEQAN_CHECKPOINT
	typedef Graph<TSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename Value<TFlowMap>::Type TFlow;

	// Initialization
	TVertexDescriptor nilPred = getNilPredecessor(g);
	resizeEdgeMap(g,flow);
	TEdgeIterator itE(g);
	for(;!atEnd(itE);goNext(itE)) {
		assignProperty(flow, getValue(itE), 0);
	}

	// Build the residual graph
	Graph<Directed<TFlow> > rG;
	_build_residual_graph(g,capacity, flow, rG);

		
	// Determine whether the sink is reachable
	String<TVertexDescriptor> predMap;
	String<TVertexDescriptor> distMap;
	breadth_first_search(rG, source, predMap, distMap);
	
	while (getProperty(predMap, sink) != nilPred) {
		TFlow inc = _get_minimum_aug(rG, predMap, source, sink);
		TEdgeIterator itEdge(g);
		for(;!atEnd(itEdge);goNext(itEdge)) {
			TVertexDescriptor u = sourceVertex(itEdge);
			TVertexDescriptor v = targetVertex(itEdge);
			TEdgeDescriptor e = getValue(itEdge);
			if (getProperty(predMap, v) == u) assignProperty(flow, e, getProperty(flow, e) + inc);
			if (getProperty(predMap, u) == v) assignProperty(flow, e, getProperty(flow, e) - inc);
		}
		// Build the residual graph
		_build_residual_graph(g,capacity, flow, rG);
		// Determine whether the sink is reachable
		clear(predMap);
		clear(distMap);
		breadth_first_search(rG, source, predMap, distMap);
	}

	TFlow valF = 0;
	TOutEdgeIterator itOutEdge(g, source);
	for(;!atEnd(itOutEdge);goNext(itOutEdge)) {
		valF += getProperty(flow, getValue(itOutEdge));
	}
	return valF;
}


//////////////////////////////////////////////////////////////////////////////
// Matching
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Path Growing Algorithm
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Function.path_growing_algorithm:
*/

template<typename TSpec, typename TWeightMap, typename TEdgeMap>
typename Value<TWeightMap>::Type
path_growing_algorithm(Graph<TSpec>& g,
					   TWeightMap const& weightMap,
					   TEdgeMap& edgeMap1)
{
	SEQAN_CHECKPOINT
	typedef Graph<TSpec> TGraph;
	typedef typename Value<TWeightMap>::Type TValue;
	typedef typename Size<Graph<TSpec> >::Type TSize;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	// Make a copy of the graph
	TGraph mutant(g);

	// Initialy not a single edge is selected
	fill(edgeMap1, getIdUpperBound(_getEdgeIdManager(g)), false);
	TEdgeMap edgeMap2 = edgeMap1;
	TValue edgeMap1Sum = 0;
	TValue edgeMap2Sum = 0;
	
	// Run the algorithm
	TSize i = 1;
	while (numEdges(mutant) > 0) {
		TVertexIterator itVert(mutant);
		while (outDegree(mutant, *itVert) < 1) goNext(itVert);
		TVertexDescriptor x = *itVert;
		TVertexDescriptor y;
		while (outDegree(mutant, x) >= 1) {
			TOutEdgeIterator itOut(mutant, x);
			TEdgeDescriptor e = *itOut;
			TValue max = getProperty(weightMap, e);
			y = targetVertex(itOut);
			goNext(itOut);
			for(;!atEnd(itOut);++itOut) {
				if (getProperty(weightMap, *itOut) > max) {
					e = *itOut;
					max = getProperty(weightMap, e);
					y = targetVertex(itOut);
				}
			}
			if (i == 1) {
				// Mark the edge for m1
				assignProperty(edgeMap1, e, true);
				edgeMap1Sum += max;
			} else {
				// Mark the edge for m2
				assignProperty(edgeMap2, e, true);
				edgeMap2Sum += max;
			}
			i = 3 - i;
			removeVertex(mutant, x);
			x = y;
		}
	}
	

	// Check whether we have to swap bool arrays
	if (edgeMap2Sum > edgeMap1Sum) {
		edgeMap1Sum = edgeMap2Sum;
		edgeMap1 = edgeMap2;
	}

	return edgeMap1Sum;
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
