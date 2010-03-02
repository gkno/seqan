#include <seqan/basic.h>
#include <seqan/consensus.h>

SEQAN_BEGIN_TESTSUITE(test_consensus) {
    // Yep, no tests here :(

    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/consensus/consensus_base.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/consensus/consensus_generated_forwards.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/consensus/consensus_library.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/consensus/consensus_realign.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/consensus/consensus_score.h");
}
SEQAN_END_TESTSUITE
