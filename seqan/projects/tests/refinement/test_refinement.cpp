#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_VERBOSE

// SeqAn Includes
#include <seqan/refinement.h>
#include "test_graph_match_refinement.h"
#include "test_graph_interval_tree.h"

//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	//Test_GraphMatchRefinement();// Test Match Refinement
	//Test_GraphIntervalTree();	// Test Interval Tree

	debug::verifyCheckpoints("projects/library/seqan/refinement/graph_impl_interval_types.h");
	debug::verifyCheckpoints("projects/library/seqan/refinement/graph_impl_interval_tree.h");
	debug::verifyCheckpoints("projects/library/seqan/refinement/graph_algorithm_refine.h");
	debug::verifyCheckpoints("projects/library/seqan/refinement/graph_algorithm_refine_scoring.h");
	debug::verifyCheckpoints("projects/library/seqan/refinement/graph_algorithm_refine_fragment.h");
	debug::verifyCheckpoints("projects/library/seqan/refinement/graph_algorithm_refine_aligngraph.h");
	debug::verifyCheckpoints("projects/library/seqan/refinement/graph_algorithm_refine_align.h");
	debug::verifyCheckpoints("projects/library/seqan/refinement/graph_algorithm_refine_exact.h");
	debug::verifyCheckpoints("projects/library/seqan/refinement/graph_algorithm_refine_annotation.h");
	debug::verifyCheckpoints("projects/library/seqan/refinement/graph_algorithm_refine_inexact.h");

	SEQAN_TREPORT("TEST END")

	return 0;
}
