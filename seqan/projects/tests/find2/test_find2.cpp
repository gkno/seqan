#include <seqan/basic.h>
#include <seqan/find2.h>

#include "test_find_approx_dpsearch.h"
#include "test_find_exact_simple.h"
#include "test_find_finder_default.h"
#include "test_find_hamming_simple.h"
#include "test_find_multiple_exact_simple.h"
#include "test_find_pattern_wild_shiftand.h"

SEQAN_BEGIN_TESTSUITE(test_find2) {
    std::cout << "Running Tests" << std::endl;
    std::cout << "=============" << std::endl;
    std::cout << std::endl;

    SEQAN_CALL_TEST(test_find2_find_finder_default_interface);

    SEQAN_CALL_TEST(test_find2_find_exact_simple_pattern_find_easy);
    SEQAN_CALL_TEST(test_find2_find_exact_simple_pattern_find_harder);
    SEQAN_CALL_TEST(test_find2_find_exact_simple_pattern_find_nomatch);
    SEQAN_CALL_TEST(test_find2_find_exact_simple_pattern_interface);
    SEQAN_CALL_TEST(test_find2_find_exact_simple_pattern_set_end_position);

    SEQAN_CALL_TEST(test_find2_find_multiple_exact_simple_pattern_find_easy);
    SEQAN_CALL_TEST(test_find2_find_multiple_exact_simple_pattern_find_harder);
    SEQAN_CALL_TEST(test_find2_find_multiple_exact_simple_pattern_find_nomatch);
    SEQAN_CALL_TEST(test_find2_find_multiple_exact_simple_pattern_interface);
    SEQAN_CALL_TEST(test_find2_find_multiple_exact_simple_pattern_set_end_position);

    SEQAN_CALL_TEST(test_find2_find_approx_dpsearch_pattern_find_easy_score_limit_0);
    SEQAN_CALL_TEST(test_find2_find_approx_dpsearch_pattern_find_easy_score_limit_1);
    SEQAN_CALL_TEST(test_find2_find_approx_dpsearch_pattern_find_easy_score_limit_0_prefix);
    SEQAN_CALL_TEST(test_find2_find_approx_dpsearch_pattern_find_easy_score_limit_1_prefix);
    SEQAN_CALL_TEST(test_find2_find_approx_dpsearch_pattern_find_easy_score_limit_1_use_score_limit);
    SEQAN_CALL_TEST(test_find2_find_approx_dpsearch_pattern_find_harder_score_limit_0);
    SEQAN_CALL_TEST(test_find2_find_approx_dpsearch_pattern_find_harder_score_limit_1);
    SEQAN_CALL_TEST(test_find2_find_approx_dpsearch_pattern_interface);
    SEQAN_CALL_TEST(test_find2_find_approx_dpsearch_pattern_set_end_position_score_limit_0);
    SEQAN_CALL_TEST(test_find2_find_approx_dpsearch_pattern_set_end_position_score_limit_1);
    SEQAN_CALL_TEST(test_find2_find_approx_dpsearch_pattern_set_end_position_score_limit_0_prefix);
    SEQAN_CALL_TEST(test_find2_find_approx_dpsearch_pattern_set_end_position_score_limit_1_prefix);
    SEQAN_CALL_TEST(test_find2_find_approx_dpsearch_pattern_find_begin);

    SEQAN_CALL_TEST(test_find2_find_hamming_simple_pattern_find_easy_score_limit_0);
    SEQAN_CALL_TEST(test_find2_find_hamming_simple_pattern_find_easy_score_limit_1);
    SEQAN_CALL_TEST(test_find2_find_hamming_simple_pattern_find_harder_score_limit_0);
    SEQAN_CALL_TEST(test_find2_find_hamming_simple_pattern_find_harder_score_limit_1);
    SEQAN_CALL_TEST(test_find2_find_hamming_simple_pattern_interface);
    SEQAN_CALL_TEST(test_find2_find_hamming_simple_pattern_set_end_position_score_limit_0);
    SEQAN_CALL_TEST(test_find2_find_hamming_simple_pattern_set_end_position_score_limit_1);

    SEQAN_CALL_TEST(test_find2_find_pattern_wild_shiftand_is_unsigned);
    SEQAN_CALL_TEST(test_find2_find_pattern_wild_shiftand_is_valid);
    SEQAN_CALL_TEST(test_find2_find_pattern_wild_shiftand_length_without_wildcards);
    SEQAN_CALL_TEST(test_find2_find_pattern_wild_shiftand_get_character_class);
    SEQAN_CALL_TEST(test_find2_find_pattern_wild_shiftand_pattern_find_easy);
    SEQAN_CALL_TEST(test_find2_find_pattern_wild_shiftand_pattern_find_harder);
    SEQAN_CALL_TEST(test_find2_find_pattern_wild_shiftand_pattern_find_nomatch);
    SEQAN_CALL_TEST(test_find2_find_pattern_wild_shiftand_pattern_interface);
    SEQAN_CALL_TEST(test_find2_find_pattern_wild_shiftand_pattern_set_end_position);

    std::cout << std::endl;
    std::cout << "Verifying Check Points" << std::endl;
    std::cout << "======================" << std::endl;
    std::cout << std::endl;

    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_approx_dpsearch.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_approx_find_begin.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_exact_simple.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_finder_default.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_hamming_simple.h");
    //SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_multiple_exact_simple.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_pattern_wild_shiftand.h");

    std::cout << std::endl;
}
SEQAN_END_TESTSUITE

