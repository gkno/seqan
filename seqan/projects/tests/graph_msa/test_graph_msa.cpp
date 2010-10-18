#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_VERBOSE


// Test path
#define TEST_PATH "projects/tests/graph_msa/"

// External / STL
#include <iostream>
#include <fstream>
#include <string>

//SeqAn
#include <seqan/graph_msa.h>
#include <seqan/misc/misc_random.h>

// Test files
#include "test_graph_tcoffee.h"


using namespace seqan;

//////////////////////////////////////////////////////////////////////////////


SEQAN_BEGIN_TESTSUITE(test_graph_msa) {
    // Call Tests.
    SEQAN_CALL_TEST(test_guide_tree_neighbour_joining);
    SEQAN_CALL_TEST(test_guide_tree_upgma_weight_avg);
    SEQAN_CALL_TEST(test_guide_tree_upgma_avg);
    SEQAN_CALL_TEST(test_guide_tree_upgma_min);
    SEQAN_CALL_TEST(test_guide_tree_upgma_max);
    SEQAN_CALL_TEST(test_distances);
	SEQAN_CALL_TEST(test_libraries);
	SEQAN_CALL_TEST(test_external_libraries);
	SEQAN_CALL_TEST(test_triplet_extension);
	SEQAN_CALL_TEST(test_sop);
	SEQAN_CALL_TEST(test_progressive);
	SEQAN_CALL_TEST(test_reversable_fragments);	
	
    // Verify Checkpoints.
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_msa/graph_align_tcoffee_kmer.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_msa/graph_align_tcoffee_distance.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_msa/graph_align_tcoffee_guidetree.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_msa/graph_align_tcoffee_library.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_msa/graph_align_tcoffee_io.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_msa/graph_align_tcoffee_base.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_msa/graph_align_tcoffee_progressive.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_msa/graph_align_tcoffee_msa.h");
}
SEQAN_END_TESTSUITE

