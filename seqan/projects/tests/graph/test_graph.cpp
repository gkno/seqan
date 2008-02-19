#define SEQAN_DEBUG
#define SEQAN_TEST
//#define SEQAN_VVERBOSE

#define TEST_PATH "projects/tests/graph/"
#define LIB_PATH "projects/library/seqan/graph/"


//#include <seqan/index.h> // For MUMPairwise_Library

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
//#include "test_graph_consensus.h"
//#include "test_graph_match_refinement.h"
//#include "test_graph_interval_tree.h"
#include "test_graph_folding.h"

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
	//Test_GraphConsensus();		// Test Graph Consensus
	//Test_GraphMatchRefinement();// Test Match Refinement
	//Test_GraphIntervalTree();	// Test Interval Tree
	Test_GraphFolding();		// Test Folding

	debug::verifyCheckpoints("projects/library/seqan/graph/graph_interface.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_drawing.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_utility_parsing.h");


	SEQAN_TREPORT("TEST END")

	return 0;
}
