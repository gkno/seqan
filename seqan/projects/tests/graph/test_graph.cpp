#define SEQAN_DEBUG
#define SEQAN_TEST
//#define SEQAN_VVERBOSE

#define TEST_PATH "projects/tests/graph/"
#define LIB_PATH "projects/library/seqan/graph/"

#include <seqan/graph.h>


#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>
#include <time.h>
#include <string>
#include <ctime>


#include "test_graph_basic.h"
#include "test_graph_types.h"
#include "test_graph_iterators.h"
#include "test_graph_properties.h"
#include "test_graph_derived.h"
#include "test_graph_algorithms.h"
#include "test_graph_alignment.h"
#include "test_graph_tcoffee.h"
#include "test_graph_match_refinement.h"

using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	Test_GraphBasics();			// Test Graph Basic
	Test_GraphTypes();			// Test Graph Types
	Test_GraphIterators();		// Test Graph Iterators
	Test_GraphProperties();		// Test internal and external property maps
	Test_GraphDerivedTypes();	// Test Additional graph types, e.g., oracle, trie,...
	Test_GraphAlgorithms();		// Test Graph Algorithms
	Test_GraphAlignment();		// Test Graph Alignment
	Test_GraphTCoffee();		// Test T-Coffee
	Test_GraphMatchRefinement();// Test Match Refinement


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
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_utility_alphabets.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_utility_parsing.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_algorithm.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_algorithm_lis_his.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_base.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_gotoh.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_hirschberg.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_interface.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_myers.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_needleman_wunsch.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_tcoffee.h");


	SEQAN_TREPORT("TEST END")

	return 0;
}
