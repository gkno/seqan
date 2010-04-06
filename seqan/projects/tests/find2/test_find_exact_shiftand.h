/*
  Test the find2/find_exact_shiftand.h header.
*/
#ifndef TESTS_FIND2_TEST_FIND_EXACT_SHIFTAND_H_
#define TESTS_FIND2_TEST_FIND_EXACT_SHIFTAND_H_

#include "test_find_exact_tests.h"

using namespace seqan;

// Test the basic interface for the ShiftAnd pattern class.
SEQAN_DEFINE_TEST(test_find2_find_exact_shiftand_pattern_interface) {
    testFind2_ExactPattern_Interface<ShiftAnd>();
}


// "Easy" (= short, few hits) test of find() using the ShiftAnd pattern.
SEQAN_DEFINE_TEST(test_find2_find_exact_shiftand_pattern_find_easy) {
    testFind2_ExactPattern_FindEasy<ShiftAnd>();
}


// "Harder" (= longer, more hits hits) test of find() using the exact pattern.
SEQAN_DEFINE_TEST(test_find2_find_exact_shiftand_pattern_find_harder) {
    testFind2_ExactPattern_FindHarder<ShiftAnd>();
}


// "Easy" (= short, few hits) test of find() using the exact pattern.
SEQAN_DEFINE_TEST(test_find2_find_exact_shiftand_pattern_find_nomatch) {
    testFind2_ExactPattern_FindNoMatch<ShiftAnd>();
}


// "Harder" (= longer, more hits hits) test of find() using the exact pattern.
SEQAN_DEFINE_TEST(test_find2_find_exact_shiftand_pattern_set_end_position) {
    testFind2_ExactPattern_SetEndPosition<ShiftAnd>();
}


// Tests for setBeginPosition() with ShiftAnd pattern.
SEQAN_DEFINE_TEST(test_find2_find_exact_shiftand_pattern_set_begin_position) {
    testFind2_ExactPattern_SetBeginPosition<ShiftAnd>();
}

#endif  // TESTS_FIND2_TEST_FIND_EXACT_SHIFTAND_H_
