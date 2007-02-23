#ifndef SEQAN_HEADER_GRAPH_ALGORITHM_H
#define SEQAN_HEADER_GRAPH_ALGORITHM_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph Algorithms
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Utility functions
//////////////////////////////////////////////////////////////////////////////
template<typename TEdges, typename TSpec>
inline typename VertexDescriptor<Graph<TEdges, TSpec> >::Type
getNilPredecessor(Graph<TEdges, TSpec> const& g)
{
	return _get_nil<typename VertexDescriptor<Graph<TEdges, TSpec> >::Type>();
}

template<typename TWeightMap>
inline typename Value<TWeightMap>::Type
getInfinityDistance(TWeightMap const& weight)
{
	// ToDo
	// We need to divide by 2 because of addition: infinity + something
	return (_get_infinity<typename Value<TWeightMap>::Type>()/2);
}

inline unsigned int
getInfinityDistance(EmptyMap)
{
	// ToDo
	// We need to divide by 2 because of addition: infinity + something
	return (_get_infinity<unsigned int>() / 2);
}

template<typename TEdges, typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
inline void 
_initialize_single_source(Graph<TEdges, TSpec> const& g,
							TVertexDescriptor const source,
							TWeightMap const& weight,
							TPredecessorMap& predecessor, 
							TDistanceMap& distance)
{
	typedef typename Value<TPredecessorMap>::Type TPredVal;
	typedef typename Value<TWeightMap>::Type TDistVal;
	TPredVal nilPred = getNilPredecessor(g);
	TDistVal infDist = getInfinityDistance(weight);
	typedef typename Iterator<Graph<TEdges, TSpec>, VertexIterator<> >::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		assignProperty(distance, getValue(it), infDist);
		assignProperty(predecessor, getValue(it), nilPred);
		goNext(it);
	}
	assignProperty(distance, source, 0);
}

template<typename TEdges, typename TSpec, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap, typename TVertexDescriptor, typename TEdgeDescriptor>
inline void 
_relax(Graph<TEdges, TSpec> const& g,
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

template<typename TEdges, typename TSpec, typename TPredecessorMap, typename TVertexDescriptor>
inline void
_print_path(Graph<TEdges, TSpec> const& g,
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

template<typename TEdges, typename TSpec, typename TPredecessor, typename TVertexDescriptor>
inline void
_print_all_pairs_shortest_path(Graph<TEdges, TSpec> const& g,
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
// Breadth-first search
//////////////////////////////////////////////////////////////////////////////
template<typename TEdges, typename TSpec, typename TVertexDescriptor, typename TPredecessorMap, typename TDistanceMap>
void
breadth_first_search(Graph<TEdges, TSpec> const& g,
					 TVertexDescriptor const source,
					 TPredecessorMap& predecessor, 
					 TDistanceMap& distance)
{
	typedef typename Iterator<Graph<TEdges, TSpec>, EdgeIterator<> >::Type TEdgeIterator;
	typedef typename Iterator<Graph<TEdges, TSpec>, VertexIterator<> >::Type TVertexIterator;
	typedef typename Value<TPredecessorMap>::Type TPredVal;
	typedef typename Value<TDistanceMap>::Type TWeightVal;

	// Initialization
	initVertexMap(g,predecessor);
	initVertexMap(g,distance);
	TPredVal nilPred = getNilPredecessor(g);
	TWeightVal infDist = getInfinityDistance(EmptyMap());
	
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
		typedef typename Iterator<Graph<TEdges, TSpec>, OutEdgeIterator<> >::Type TOutEdgeIterator;
		TOutEdgeIterator itout(g,u);
		for(;!atEnd(itout);goNext(itout)) {
			TVertexDescriptor v = targetVertex(g, getValue(itout));
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
template<typename TEdges, typename TSpec, typename TVertexDescriptor, typename TTokenMap, typename TPredecessorMap, typename TDiscoveryTimeMap, typename TFinishingTimeMap, typename TVal>
void
_dfs_visit(Graph<TEdges, TSpec> const& g,
		   TVertexDescriptor const u,
		   TTokenMap& tokenMap,
		   TPredecessorMap& predecessor,
		   TDiscoveryTimeMap& disc,
		   TFinishingTimeMap& finish,
		   TVal& time)
{
	assignProperty(tokenMap, u, true);
	++time;
	assignProperty(disc, u, time);
	typedef typename Iterator<Graph<TEdges, TSpec>, AdjacencyIterator<> >::Type TAdjacencyIterator;
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


template<typename TEdges, typename TSpec, typename TPredecessorMap, typename TDiscoveryTimeMap, typename TFinishingTimeMap>
void
depth_first_search(Graph<TEdges, TSpec> const& g,
					 TPredecessorMap& predecessor, 
					 TDiscoveryTimeMap& disc,
					 TFinishingTimeMap& finish)
{
	typedef typename Iterator<Graph<TEdges, TSpec>, EdgeIterator<> >::Type TEdgeIterator;
	typedef typename Iterator<Graph<TEdges, TSpec>, VertexIterator<> >::Type TVertexIterator;
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
		typedef typename VertexDescriptor<Graph<TEdges, TSpec> >::Type TVertexDescriptor;
		TVertexDescriptor u = getValue(it);
		if (getProperty(tokenMap, u) == false) {
			_dfs_visit(g, u, tokenMap, predecessor, disc, finish, time);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// Topological sort
//////////////////////////////////////////////////////////////////////////////
template<typename TEdges, typename TSpec, typename TVertexDescriptor>
void
topological_sort(Graph<TEdges, TSpec> const& g,
				 String<TVertexDescriptor>& topSort)
{
	// Initialization
	String<unsigned int> predMap;
	String<unsigned int> discoveryTimeMap;
	String<unsigned int> finishingTimeMap;
	
	// Dfs
	depth_first_search(g, predMap, discoveryTimeMap, finishingTimeMap);

	// Order vertices
	typedef ::std::pair<unsigned int, unsigned int> TTimeVertexPair;
	std::priority_queue<TTimeVertexPair> q;
	typedef typename Iterator<Graph<TEdges, TSpec>, VertexIterator<> >::Type TVertexIterator;
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
template<typename TEdges, typename TSpec, typename TComponents>
void
strongly_connected_components(Graph<TEdges, TSpec> const& g_source,
				 				TComponents& components)
{
	// Initialization
	typedef typename Iterator<Graph<TEdges, TSpec>, EdgeIterator<> >::Type TEdgeIterator;
	typedef typename Iterator<Graph<TEdges, TSpec>, VertexIterator<> >::Type TVertexIterator;
	typedef typename Value<TComponents>::Type TCompVal;
	initVertexMap(g_source,components);
	String<unsigned int> predMap;
	String<unsigned int> discoveryTimeMap;
	String<unsigned int> finishingTimeMap;
	
	// Dfs
	depth_first_search(g_source, predMap, discoveryTimeMap, finishingTimeMap);

	Graph<TEdges, TSpec> g;
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
	typedef ::std::pair<unsigned int, unsigned int> TTimeVertexPair;
	std::priority_queue<TTimeVertexPair> q;
	goBegin(it);
	for(;!atEnd(it);++it) {
		q.push(std::make_pair(getProperty(finishingTimeMap, getValue(it)), getValue(it)));
	}



	unsigned int time = 0;
	unsigned int label = 0;
	while(!q.empty()) {
		typedef typename VertexDescriptor<Graph<TEdges, TSpec> >::Type TVertexDescriptor;
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
// DAG Shortest Path
//////////////////////////////////////////////////////////////////////////////
template<typename TEdges, typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
void
dag_shortest_path(Graph<TEdges, TSpec> const& g,
					TVertexDescriptor const source,
					TWeightMap const& weight,
					TPredecessorMap& predecessor, 
					TDistanceMap& distance)
{
	// Initialization
	typedef typename EdgeDescriptor<Graph<TEdges, TSpec> >::Type TEdgeDescriptor;
	typedef typename Iterator<Graph<TEdges, TSpec>, EdgeIterator<> >::Type TEdgeIterator;
	initVertexMap(g,predecessor);
	initVertexMap(g,distance);

	// Topological sort
	String<TVertexDescriptor> order;
	topological_sort(g, order);

	_initialize_single_source(g, source, weight, predecessor, distance);

	//DAG Shortest Paths
	typedef typename Iterator<String<TVertexDescriptor> >::Type TStringIterator;
	TStringIterator it = begin(order);
	while(!atEnd(it)) {
		typedef typename Iterator<Graph<TEdges, TSpec>, OutEdgeIterator<> >::Type TOutEdgeIterator;
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
template<typename TEdges, typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
bool 
bellman_ford_algorithm(Graph<TEdges, TSpec> const& g,
						TVertexDescriptor const source,
						TWeightMap const& weight,
						TPredecessorMap& predecessor, 
						TDistanceMap& distance)
{
	// Initialization
	typedef typename EdgeDescriptor<Graph<TEdges, TSpec> >::Type TEdgeDescriptor;
	typedef typename Iterator<Graph<TEdges, TSpec>, VertexIterator<> >::Type TVertexIterator;
	typedef typename Iterator<Graph<TEdges, TSpec>, OutEdgeIterator<> >::Type TOutEdgeIterator;
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
template<typename TEdges, typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
void 
dijkstra(Graph<TEdges, TSpec> const& g,
			TVertexDescriptor const source,
			TWeightMap const& weight,
			TPredecessorMap& predecessor, 
			TDistanceMap& distance)
{
	// Initialization
	typedef typename Value<TDistanceMap>::Type TDistVal;
	typedef typename Iterator<Graph<TEdges, TSpec>, VertexIterator<> >::Type TVertexIterator;
	typedef typename Iterator<Graph<TEdges, TSpec>, OutEdgeIterator<> >::Type TOutEdgeIterator;
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
// All-Pairs shortest path
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
template<typename TEdges, typename TSpec, typename TWeightMap, typename TMatrix, typename TPredecessor>
void 
_initialize_all_pairs(Graph<TEdges, TSpec> const& g,
						TWeightMap const& weight,
						TMatrix& matrix,
						TPredecessor& predecessor)
{
	typedef typename VertexDescriptor<Graph<TEdges, TSpec> >::Type TVertexDescriptor;
	typedef typename Iterator<Graph<TEdges, TSpec>, VertexIterator<> >::Type TVertexIterator;
	typedef typename Iterator<Graph<TEdges, TSpec>, OutEdgeIterator<> >::Type TOutEdgeIterator;
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

template<typename TEdges, typename TSpec, typename TWeightMap, typename TMatrix, typename TPredecessor>
void 
all_pairs_shortest_path(Graph<TEdges, TSpec> const& g,
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

template<typename TEdges, typename TSpec, typename TWeightMap, typename TMatrix, typename TPredecessor>
void 
floyd_warshall(Graph<TEdges, TSpec> const& g,
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

template<typename TEdges, typename TSpec, typename TMatrix>
void 
transitive_closure(Graph<TEdges, TSpec> const& g,
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
