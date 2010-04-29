// Test path
#define TEST_PATH "projects/tests/graph_align/"

#include <seqan/basic.h>
#include <seqan/map.h>
#include <seqan/graph_align.h>
#include "test_graph_align.h"

SEQAN_BEGIN_TESTSUITE(test_graph_align) {
    // global alignments
    SEQAN_CALL_TEST(test_graph_align_needleman_wunsch);
    SEQAN_CALL_TEST(test_graph_align_gotoh);
    SEQAN_CALL_TEST(test_graph_align_hirschberg);
    SEQAN_CALL_TEST(test_graph_align_allAgainstAll);
    SEQAN_CALL_TEST(test_graph_align_gotohVsBandedGotoh);
    
    // local alignments
    SEQAN_CALL_TEST(test_graph_align_smith_waterman);
    SEQAN_CALL_TEST(test_graph_align_smith_waterman_clump);
    SEQAN_CALL_TEST(test_graph_align_banded_smith_waterman_clump);
    
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_align/graph_align_base.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_align/graph_align_config.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_align/graph_align_interface.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_align/graph_align_needleman_wunsch.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_align/graph_align_gotoh.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_align/graph_align_hirschberg.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_align/graph_align_smith_waterman.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_align/graph_align_smith_waterman_clump.h");
}
SEQAN_END_TESTSUITE
