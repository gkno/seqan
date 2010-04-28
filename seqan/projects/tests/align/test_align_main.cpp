#include <iostream>
#include <cstdio>
#include <vector>

#include <seqan/basic.h>

// Definitions to externally defined tests.
// TODO(holtgrew): Remove, each test .cpp should be a test suite.
void SEQAN_TEST_test_align_align_char_string_array_gaps();
void SEQAN_TEST_test_align_align_dna_array_gaps();
void SEQAN_TEST_test_align_align_char_string_sumlist_gaps();
void SEQAN_TEST_test_align_align_dna_sumlist_gaps();

void SEQAN_TEST_test_align_gaps_base_char_string_array_gaps();
void SEQAN_TEST_test_align_gaps_base_char_string_sumlist_gaps();
void SEQAN_TEST_test_align_gaps_test_gaps_iterator();
void SEQAN_TEST_test_align_gaps_test_gap_manipulation_char_string_array_gaps();
void SEQAN_TEST_test_align_gaps_test_gap_manipulation_char_string_sumlist_gaps();
void SEQAN_TEST_test_align_gaps_test_sequence_gaps_base();

void SEQAN_TEST_test_align_myers_test_short();
void SEQAN_TEST_test_align_myers_test_long();
void SEQAN_TEST_test_align_hirschberger();

void SEQAN_TEST_testLocalAlign();
void SEQAN_TEST_testLocalAlign2();
void SEQAN_TEST_testBandedLocalAlign();


SEQAN_BEGIN_TESTSUITE("test_align") {
    // Call tests from test_align_align.cpp.
    SEQAN_CALL_TEST(test_align_align_char_string_array_gaps);
    SEQAN_CALL_TEST(test_align_align_dna_array_gaps);
    SEQAN_CALL_TEST(test_align_align_char_string_sumlist_gaps);
    SEQAN_CALL_TEST(test_align_align_dna_sumlist_gaps);

    SEQAN_CALL_TEST(test_align_gaps_base_char_string_array_gaps);
    SEQAN_CALL_TEST(test_align_gaps_base_char_string_sumlist_gaps);
    SEQAN_CALL_TEST(test_align_gaps_test_gaps_iterator);
    SEQAN_CALL_TEST(test_align_gaps_test_gap_manipulation_char_string_array_gaps);
    SEQAN_CALL_TEST(test_align_gaps_test_gap_manipulation_char_string_sumlist_gaps);
    SEQAN_CALL_TEST(test_align_gaps_test_sequence_gaps_base);

    SEQAN_CALL_TEST(testLocalAlign);
    SEQAN_CALL_TEST(testLocalAlign2);
    SEQAN_CALL_TEST(testBandedLocalAlign);

    SEQAN_CALL_TEST(test_align_myers_test_short);
    SEQAN_CALL_TEST(test_align_myers_test_long);
    SEQAN_CALL_TEST(test_align_hirschberger);

    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/align/align_algorithms.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/align/align_base.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/align/align_cols_base.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/align/align_dynprog.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/align/align_hirschberg.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/align/align_iterator_base.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/align/align_local_dynprog.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/align/align_myers.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/align/align_trace.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/align/gaps_array.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/align/gaps_base.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/align/gaps_iterator_base.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/align/gaps_simple.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/align/gaps_sumlist.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/align/hirschberg_set.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/align/matrix_base.h");
}
SEQAN_END_TESTSUITE
