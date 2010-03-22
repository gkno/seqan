#include <seqan/basic.h>
#include <seqan/find2.h>

SEQAN_DEFINE_TEST(test_find2_foo) {
    SEQAN_ASSERT_TRUE(false, "Write me!");
}

SEQAN_BEGIN_TESTSUITE(test_find2) {
    SEQAN_CALL_TEST(test_find2_foo);

    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_approx_dpsearch.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_exact_simple.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_hamming_simple.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_multiple_exact_shiftand.h");
}
SEQAN_END_TESTSUITE
