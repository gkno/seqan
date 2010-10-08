#ifndef SEQAN_HEADER_TEST_GRAPH_ITERATORS_H
#define SEQAN_HEADER_TEST_GRAPH_ITERATORS_H

using namespace seqan;

//////////////////////////////////////////////////////////////////////////////
template <typename TGraphType>
void Test_VertexIterator() {
//____________________________________________________________________________
// Graph InternalVertexIterator
	typedef Graph<TGraphType> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	
	TGraph g;
	TVertexDescriptor v0 = addVertex(g);
	TVertexDescriptor v1 = addVertex(g);
	addEdge(g,v0,v1,'t');
	addVertex(g); //2
	addVertex(g); //3
	addVertex(g); //4

	//Tricky case -> id 0 is released
	removeVertex(g,v0);
	removeVertex(g,3);
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	TVertexIterator it(g);
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
	SEQAN_ASSERT_TRUE(getValue(it)==1);
	SEQAN_ASSERT_TRUE(value(it)==1);
	SEQAN_ASSERT_TRUE(getValue(it)==1);
	goNext(it);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	SEQAN_ASSERT_TRUE(getValue(it)==2);
	++it;
	SEQAN_ASSERT_TRUE(getValue(it)==4);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	it++; 
	SEQAN_ASSERT_TRUE(atEnd(it)==true);
	++it; 
	SEQAN_ASSERT_TRUE(atEnd(it)==true);
	goPrevious(it);
	// No assignment to vertex iterators
	// *it = 3;
	SEQAN_ASSERT_TRUE((*it)==4);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	it--;
	SEQAN_ASSERT_TRUE(getValue(it)==2);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	--it;
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
	SEQAN_ASSERT_TRUE(getValue(it)==1);
	--it;
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
	SEQAN_ASSERT_TRUE(getValue(it)==1);
	TVertexIterator it2(g);
	TVertexIterator it3;
	it3 = it;
	SEQAN_ASSERT_TRUE(it == it2);
	SEQAN_ASSERT_TRUE(it2 == it3);
	goEnd(it);
	SEQAN_ASSERT_TRUE(it2 != it);
	goEnd(it2);
	SEQAN_ASSERT_TRUE(it2 == it);
	goBegin(it2);
	SEQAN_ASSERT_TRUE(it2 != it);
	SEQAN_ASSERT_TRUE(&hostGraph(it) == &g);
}

//////////////////////////////////////////////////////////////////////////////
void Test_TreeInternalVertexIterator() {
//____________________________________________________________________________
// Tree InternalVertexIterator
	typedef Graph<Tree<void> > TTree;
	typedef VertexDescriptor<TTree>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TTree>::Type TEdgeDescriptor;
	typedef Size<TTree>::Type TSize;
	
	TTree gV;
	TSize numEdges = 8;
	//Parent, Child, Parent, Child, ...
	//The root must be the first vertex
	TVertexDescriptor edges[] = {0,8, 0,3, 0,2, 0,1, 2,4, 4,5, 5,7, 5,6};
	addEdges(gV,edges, numEdges);

	typedef Iterator<TTree, VertexIterator>::Type TVertexIterator;
	TVertexIterator itV(gV);
	SEQAN_ASSERT_TRUE(atBegin(itV)==true);
	SEQAN_ASSERT_TRUE(getValue(itV)==0);
	SEQAN_ASSERT_TRUE(value(itV)==0);
	SEQAN_ASSERT_TRUE(getValue(itV)==0);
	goNext(itV);
	SEQAN_ASSERT_TRUE(atBegin(itV)==false);
	SEQAN_ASSERT_TRUE(getValue(itV)==1);
	++itV;
	SEQAN_ASSERT_TRUE(getValue(itV)==2);
	SEQAN_ASSERT_TRUE(atEnd(itV)==false);
	goPrevious(itV);
	SEQAN_ASSERT_TRUE((*itV)==1);
	SEQAN_ASSERT_TRUE(atEnd(itV)==false);
	itV--;
	SEQAN_ASSERT_TRUE(getValue(itV)==0);
	SEQAN_ASSERT_TRUE(atBegin(itV)==true);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TGraphType>
void Test_OutEdgeIterator() {
//____________________________________________________________________________
// Graph InternalOutEdgeIterator
	typedef Graph<TGraphType> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	
	TGraph g;
	TVertexDescriptor v0 = addVertex(g);
	TVertexDescriptor v1 = addVertex(g);
	addVertex(g); //2
	addEdge(g,0,2,'t');
	addEdge(g,v0,v1,'a');
	addVertex(g); //3
	addVertex(g); //4

	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	TOutEdgeIterator it(g, v0);
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(it))==0);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(it))==1);
	SEQAN_ASSERT_TRUE(sourceVertex(g, value(it))==0);
	SEQAN_ASSERT_TRUE(targetVertex(g, *it)==1);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
	goNext(it);
	// Slow
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(it))==0);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(it))==2);
	// Fast
	SEQAN_ASSERT_TRUE(sourceVertex(it)==0);
	SEQAN_ASSERT_TRUE(targetVertex(it)==2);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	++it;
	SEQAN_ASSERT_TRUE(atEnd(it)==true);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	goPrevious(it);
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(it))==0);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(it))==2);
	--it;
	// Slow
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(it))==0);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(it))==1);
	// Fast
	SEQAN_ASSERT_TRUE(sourceVertex(it)==0);
	SEQAN_ASSERT_TRUE(targetVertex(it)==1);
	it++; 
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(it))==0);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(it))==2);
	it--;
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(it))==0);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(it))==1);
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
	TOutEdgeIterator it2(g, v0);
	TOutEdgeIterator it3;
	it3 = it;
	SEQAN_ASSERT_TRUE(it == it2);
	SEQAN_ASSERT_TRUE(it2 == it3);
	goEnd(it);
	SEQAN_ASSERT_TRUE(it2 != it);
	goEnd(it2);
	SEQAN_ASSERT_TRUE(it2 == it);
	goBegin(it2);
	SEQAN_ASSERT_TRUE(it2 != it);
	SEQAN_ASSERT_TRUE(&g == &hostGraph(it));
}


//////////////////////////////////////////////////////////////////////////////

	
template <typename TGraphType>
void Test_EdgeIterator() {
//____________________________________________________________________________
// Graph InternalEdgeIterator
	typedef Graph<TGraphType> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	
	TGraph g;
	TVertexDescriptor v0 = addVertex(g);
	TVertexDescriptor v1 = addVertex(g);
	addVertex(g); //2
	addEdge(g,0,2,'t');
	addEdge(g,v0,v1,'a');
	addVertex(g); //3
	addVertex(g); //4
	addVertex(g); //5
	addEdge(g,3,4,'t');
	addEdge(g,2,3,'a');
	addEdge(g,4,5,'a');

	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	TEdgeIterator it(g);
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(it))==0);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(it))==1);
	SEQAN_ASSERT_TRUE(sourceVertex(g, value(it))==0);
	SEQAN_ASSERT_TRUE(targetVertex(g, *it)==1);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
	goNext(it);
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(it))==0);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(it))==2);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	++it;
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	// Slow
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(it))==2);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(it))==3);
	// Fast
	SEQAN_ASSERT_TRUE(sourceVertex(it)==2);
	SEQAN_ASSERT_TRUE(targetVertex(it)==3);
	it++;
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(it))==3);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(it))==4);
	it++;
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(it))==4);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(it))==5);
	++it;
	SEQAN_ASSERT_TRUE(atEnd(it)==true);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	++it;
	++it;
	goPrevious(it);
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(it))==4);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(it))==5);
	--it;
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(it))==3);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(it))==4);
	it--;
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(it))==2);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(it))==3);
	TEdgeIterator it2(g);
	TEdgeIterator it3;
	goBegin(it);
	it3 = it;
	SEQAN_ASSERT_TRUE(it == it2);
	SEQAN_ASSERT_TRUE(it2 == it3);
	goEnd(it);
	SEQAN_ASSERT_TRUE(it2 != it);
	goEnd(it2);
	SEQAN_ASSERT_TRUE(it2 == it);
	goBegin(it2);
	SEQAN_ASSERT_TRUE(it2 != it);
	SEQAN_ASSERT_TRUE(&hostGraph(it) == &g);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TGraphType>
void Test_AdjacencyIterator() {
//____________________________________________________________________________
// Graph InternalAdjacencyIterator
	typedef Graph<TGraphType> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	
	TGraph g;
	TVertexDescriptor v0 = addVertex(g);
	TVertexDescriptor v1 = addVertex(g);
	addVertex(g); //2
	addEdge(g,0,2,'t');
	addEdge(g,v0,v1,'a');
	addVertex(g); //3
	addVertex(g); //4
	addVertex(g); //5
	addEdge(g,3,2,'t');
	addEdge(g,3,4,'a');
	addEdge(g,4,5,'t');

	typedef typename Iterator<TGraph, AdjacencyIterator>::Type TAdjacencyIterator;
	TAdjacencyIterator it(g,3);
	SEQAN_ASSERT_TRUE(getValue(it) == 4);
	SEQAN_ASSERT_TRUE(&hostGraph(it) == &g);
	SEQAN_ASSERT_TRUE(value(it) == 4);
	SEQAN_ASSERT_TRUE(*it == 4);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
	goNext(it);
	SEQAN_ASSERT_TRUE(getValue(it)==2);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	++it;
	SEQAN_ASSERT_TRUE(atEnd(it)==true);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	goBegin(it);
	it++;
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	it++;
	goPrevious(it);
	SEQAN_ASSERT_TRUE(getValue(it)==2);
	--it;
	SEQAN_ASSERT_TRUE(getValue(it)==4);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
	goEnd(it);
	it--;
	SEQAN_ASSERT_TRUE(getValue(it)==2);
	goBegin(it);
	TAdjacencyIterator it2(g, 3);
	TAdjacencyIterator it3;
	it3 = it;
	SEQAN_ASSERT_TRUE(it == it2);
	SEQAN_ASSERT_TRUE(it2 == it3);
	goEnd(it);
	SEQAN_ASSERT_TRUE(it2 != it);
	goEnd(it2);
	SEQAN_ASSERT_TRUE(it2 == it);
	goBegin(it2);
	SEQAN_ASSERT_TRUE(it2 != it);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TGraphType>
void Test_BfsIter() {
//____________________________________________________________________________
// Graph BfsIterator
	typedef Graph<TGraphType> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;

	//Create the graph
	TGraph g;
	addVertex(g);addVertex(g);addVertex(g);addVertex(g);
	addEdge(g,0,2,'t');
	addEdge(g,0,1,'a');
	addEdge(g,1,3,'a');
	
	typedef typename Iterator<TGraph, BfsIterator>::Type TBfsIterator;
	TBfsIterator it(g,0);
	SEQAN_ASSERT_TRUE(getValue(it) == 0);
	SEQAN_ASSERT_TRUE(&hostGraph(it) == &g);
	SEQAN_ASSERT_TRUE(value(it) == 0);
	SEQAN_ASSERT_TRUE(*it == 0);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
	goNext(it);
	SEQAN_ASSERT_TRUE(getValue(it)==1);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	goNext(it);
	SEQAN_ASSERT_TRUE(getValue(it)==2);
	goNext(it);
	SEQAN_ASSERT_TRUE(getValue(it)==3);
	goNext(it);
	SEQAN_ASSERT_TRUE(atEnd(it)==true);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	goBegin(it);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
}

//////////////////////////////////////////////////////////////////////////////

void Test_BfsIterator() {
//____________________________________________________________________________
// Graph InternalBfsIterator
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Create the graph
	Graph<> g;
	addVertex(g);addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);addVertex(g);addVertex(g);
	addEdge(g,1,0);
	addEdge(g,1,5);
	addEdge(g,0,4);
	addEdge(g,5,2);
	addEdge(g,5,6);
	addEdge(g,2,6);
	addEdge(g,2,3);
	addEdge(g,6,3);
	addEdge(g,6,7);
	addEdge(g,3,7);
	typedef Iterator<Graph<>, BfsIterator>::Type TBfsIterator;
	TBfsIterator it(g,1);
	SEQAN_ASSERT_TRUE(getValue(it) == 1);
	SEQAN_ASSERT_TRUE(&hostGraph(it) == &g);
	SEQAN_ASSERT_TRUE(value(it) == 1);
	SEQAN_ASSERT_TRUE(*it == 1);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
	goNext(it);
	SEQAN_ASSERT_TRUE(getValue(it)==5);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	goNext(it);
	SEQAN_ASSERT_TRUE(getValue(it)==0);
	goNext(it);
	SEQAN_ASSERT_TRUE(getValue(it)==6);
	goNext(it);
	SEQAN_ASSERT_TRUE(getValue(it)==2);
	++it;
	SEQAN_ASSERT_TRUE(getValue(it)==4);
	it++;
	SEQAN_ASSERT_TRUE(getValue(it)==7);
	goNext(it);
	SEQAN_ASSERT_TRUE(getValue(it)==3);
	goNext(it);
	SEQAN_ASSERT_TRUE(atEnd(it)==true);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	goBegin(it);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
	TBfsIterator it2(g, 1);
	TBfsIterator it3;
	it3 = it;
	SEQAN_ASSERT_TRUE(it == it2);
	SEQAN_ASSERT_TRUE(it2 == it3);
	goEnd(it);
	SEQAN_ASSERT_TRUE(it2 != it);
	goEnd(it2);
	SEQAN_ASSERT_TRUE(it2 == it);
	goBegin(it2);
	SEQAN_ASSERT_TRUE(it2 != it);
}

//////////////////////////////////////////////////////////////////////////////


template <typename TGraphType>
void Test_DfsPreorderIter() {
//____________________________________________________________________________
// Graph DfsIterator
	typedef Graph<TGraphType> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;

	//Create the graph
	TGraph g;
	addVertex(g);addVertex(g);addVertex(g);addVertex(g);
	addEdge(g,0,1,'t');
	addEdge(g,0,2,'a');
	addEdge(g,1,3,'a');
	
	typedef typename Iterator<TGraph, DfsPreorder>::Type TDfsIterator;
	TDfsIterator it(g,0);
	SEQAN_ASSERT_TRUE(getValue(it) == 0);
	SEQAN_ASSERT_TRUE(&hostGraph(it) == &g);
	SEQAN_ASSERT_TRUE(value(it) == 0);
	SEQAN_ASSERT_TRUE(*it == 0);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
	goNext(it);
	SEQAN_ASSERT_TRUE(getValue(it)==1);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	goNext(it);
	SEQAN_ASSERT_TRUE(getValue(it)==3);
	goNext(it);
	SEQAN_ASSERT_TRUE(getValue(it)==2);
	goNext(it);
	SEQAN_ASSERT_TRUE(atEnd(it)==true);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	goBegin(it);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
}


//////////////////////////////////////////////////////////////////////////////

void Test_DfsPreorderIterator() {
//____________________________________________________________________________
// Graph DfsIterator
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Create the graph
	Graph<> g;
	addVertex(g);addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);addVertex(g);addVertex(g);
	addEdge(g,1,0);
	addEdge(g,1,5);
	addEdge(g,0,4);
	addEdge(g,5,2);
	addEdge(g,5,6);
	addEdge(g,2,6);
	addEdge(g,2,3);
	addEdge(g,6,3);
	addEdge(g,6,7);
	addEdge(g,3,7);

	typedef Iterator<Graph<>, DfsPreorder>::Type TDfsPreorder;
	TDfsPreorder it(g,1);
	SEQAN_ASSERT_TRUE(getValue(it) == 1);
	SEQAN_ASSERT_TRUE(&hostGraph(it) == &g);
	SEQAN_ASSERT_TRUE(value(it) == 1);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
	goNext(it);
	SEQAN_ASSERT_TRUE(getValue(it)==0);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	goNext(it);
	SEQAN_ASSERT_TRUE(getValue(it)==4);
	goNext(it);
	SEQAN_ASSERT_TRUE(getValue(it)==5);
	SEQAN_ASSERT_TRUE(value(it)==5);
	SEQAN_ASSERT_TRUE(*it == 5);
	goNext(it);
	SEQAN_ASSERT_TRUE(getValue(it)==2);
	++it;
	SEQAN_ASSERT_TRUE(getValue(it)==3);
	it++;
	SEQAN_ASSERT_TRUE(getValue(it)==7);
	goNext(it);
	SEQAN_ASSERT_TRUE(getValue(it)==6);
	goNext(it);
	SEQAN_ASSERT_TRUE(atEnd(it)==true);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	goBegin(it);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
	TDfsPreorder it2(g, 1);
	TDfsPreorder it3;
	it3 = it;
	SEQAN_ASSERT_TRUE(it == it2);
	SEQAN_ASSERT_TRUE(it2 == it3);
	goEnd(it);
	SEQAN_ASSERT_TRUE(it2 != it);
	goEnd(it2);
	SEQAN_ASSERT_TRUE(it2 == it);
	goBegin(it2);
	SEQAN_ASSERT_TRUE(it2 != it);
}


//////////////////////////////////////////////////////////////////////////////

void Test_GraphIterators() {
	Test_VertexIterator<Directed<char> >();
	Test_VertexIterator<Undirected<char> >();
	Test_VertexIterator<Automaton<char> >();
	Test_VertexIterator<Hmm<Dna, char> >();
	Test_TreeInternalVertexIterator();
	
	Test_OutEdgeIterator<Directed<char> >();
	Test_OutEdgeIterator<Undirected<char> >();
	Test_OutEdgeIterator<Tree<char> >();
	Test_OutEdgeIterator<Automaton<char> >();
	Test_OutEdgeIterator<Hmm<Dna, char> >();
	
	Test_EdgeIterator<Directed<char> >();
	Test_EdgeIterator<Undirected<char> >();
	Test_EdgeIterator<Tree<char> >();
	Test_EdgeIterator<Automaton<char> >();
	Test_EdgeIterator<Hmm<Dna, char> >();
	
	Test_AdjacencyIterator<Directed<char> >();
	Test_AdjacencyIterator<Undirected<char> >();
	Test_AdjacencyIterator<Tree<char> >();
	Test_AdjacencyIterator<Automaton<char> >();
	Test_AdjacencyIterator<Hmm<Dna, char> >();
	
	// Test bfs and dfs iterator
	Test_BfsIter<Directed<char> >();
	Test_BfsIter<Undirected<char> >();
	Test_BfsIter<Tree<char> >();
	Test_BfsIter<Automaton<char> >();
	Test_BfsIter<Hmm<Dna, char> >();
	Test_BfsIterator();
	
	Test_DfsPreorderIter<Directed<char> >();
	Test_DfsPreorderIter<Undirected<char> >();
	Test_DfsPreorderIter<Tree<char> >();
	Test_DfsPreorderIter<Automaton<char> >();
	Test_DfsPreorderIter<Hmm<Dna, char> >();
	Test_DfsPreorderIterator();
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(test_graph_iterators)
{
	Test_GraphIterators();
}

#endif

