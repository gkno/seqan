#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_VVERBOSE

// Seqan
#include <seqan/graph_types.h>

// Test files
#include "test_graph_basic.h"
#include "test_graph_types.h"
#include "test_graph_iterators.h"
#include "test_graph_properties.h"
#include "test_graph_derived.h"

//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	Test_GraphBasics();			// Test Graph Basic
	Test_GraphTypes();			// Test Graph Types
	Test_GraphIterators();		// Test Graph Iterators
	Test_GraphProperties();		// Test internal and external property maps
	Test_GraphDerivedTypes();	// Test Additional graph types, e.g., oracle, trie,...

	debug::verifyCheckpoints("projects/library/seqan/graph_types/graph_interface.h");

	
	SEQAN_TREPORT("TEST END")

	return 0;
}
