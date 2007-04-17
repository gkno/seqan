#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>
#include <time.h>
#include <string>

#define SEQAN_DEBUG
#define SEQAN_TEST

#define TEST_PATH "projects/tests/graph/"
#define LIB_PATH "projects/library/seqan/graph/"

#include <seqan/graph.h>
#include "test_graph_basic.h"
#include "test_graph_types.h"
#include "test_graph_iterators.h"
#include "test_graph_properties.h"
#include "test_graph_derived.h"
#include "test_graph_algorithms.h"

using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	// Test Graph Basic
	Test_IdManager();	// Test Id Manager
	Test_EdgeStump();	// Test EdgeStumps


	// Test Graph Types
	Test_Directed();	// Directed graphs
	Test_Undirected();	// Undirected graphs
	Test_Automaton();	// Automatons
	Test_WordGraph();	// Word Graph
	Test_Tree();		// Trees
	Test_Alignment();	// Alignment graph

	// Test Graph Iterators
	Test_VertexIterator<Directed<char> >();
	Test_VertexIterator<Undirected<char> >();
	Test_VertexIterator<Automaton<char> >();
	Test_TreeInternalVertexIterator();
	Test_OutEdgeIterator<Directed<char> >();
	Test_OutEdgeIterator<Undirected<char> >();
	Test_OutEdgeIterator<Tree<char> >();
	Test_OutEdgeIterator<Automaton<char> >();
	Test_EdgeIterator<Directed<char> >();
	Test_EdgeIterator<Undirected<char> >();
	Test_EdgeIterator<Tree<char> >();
	Test_EdgeIterator<Automaton<char> >();
	Test_AdjacencyIterator<Directed<char> >();
	Test_AdjacencyIterator<Undirected<char> >();
	Test_AdjacencyIterator<Tree<char> >();
	Test_AdjacencyIterator<Automaton<char> >();
	// Test bfs and dfs iterator
	Test_BfsIter<Directed<char> >();
	Test_BfsIter<Undirected<char> >();
	Test_BfsIter<Tree<char> >();
	Test_BfsIter<Automaton<char> >();
	Test_BfsIterator();
	Test_DfsPreorderIter<Directed<char> >();
	Test_DfsPreorderIter<Undirected<char> >();
	Test_DfsPreorderIter<Tree<char> >();
	Test_DfsPreorderIter<Automaton<char> >();
	Test_DfsPreorderIterator();

	// Test Graph Properties
	Test_ExternalProperty<Directed<char> >();
	Test_ExternalProperty<Undirected<char> >();
	Test_ExternalProperty<Tree<char> >();
	Test_ExternalProperty<Automaton<char> >();	
	Test_Property();

	// Test Graph Derived
	Test_Oracle();
	Test_Trie();

	// Test Graph Algorithms
	Test_Algorithms();


	// T-Coffee
	Test_TCoffee();


//____________________________________________________________________________
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_base.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_idmanager.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_interface.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_edgestump.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_directed.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_undirected.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_automaton.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_wordgraph.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_tree.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_align.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_iterator.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_iterator_vertex.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_iterator_outedge.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_iterator_adjacency.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_iterator_edge.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_property.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_oracle.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_trie.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_drawing.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_iterator_bfs.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_iterator_dfs.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_algorithm.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_algorithm_tcoffee.h");

	SEQAN_TREPORT("TEST END")

	return 0;
}
