#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_VERBOSE

// SeqAn
#include <seqan/graph_algorithms.h>
#include "test_graph_algorithms.h"

// SeqAn Namespace
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

SEQAN_BEGIN_TESTSUITE(test_graph_algorithms) {
	mtRandInit();

	SEQAN_CALL_TEST(test_heap_tree);
	SEQAN_CALL_TEST(test_breadth_first_search);
	SEQAN_CALL_TEST(test_depth_first_search);
	SEQAN_CALL_TEST(test_topological_sort);
	SEQAN_CALL_TEST(test_strongly_connected_components);
	SEQAN_CALL_TEST(test_connected_components);
	SEQAN_CALL_TEST(test_prims_algorithm);
	SEQAN_CALL_TEST(test_kruskals_algorithm);
	SEQAN_CALL_TEST(test_mst_all);
	SEQAN_CALL_TEST(test_dag_shortest_path);
	SEQAN_CALL_TEST(test_bellmann_ford);
	SEQAN_CALL_TEST(test_dijkstra);
	SEQAN_CALL_TEST(test_all_pairs_shortest_path);
	SEQAN_CALL_TEST(test_floyd_warshall);
	SEQAN_CALL_TEST(test_transitive_closure);
	SEQAN_CALL_TEST(test_ford_fulkerson);
	SEQAN_CALL_TEST(test_path_growing_algorithm);
	SEQAN_CALL_TEST(test_longest_increasing_subsequence);
	SEQAN_CALL_TEST(test_longest_common_subsequence);
	SEQAN_CALL_TEST(test_heaviest_increasing_subsequence);
	SEQAN_CALL_TEST(test_hmm_algorithm);	

    // Verify Checkpoints.
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_algorithms/graph_algorithm.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_algorithms/graph_algorithm_hmm.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_algorithms/graph_algorithm_lis_his.h");
}
SEQAN_END_TESTSUITE

