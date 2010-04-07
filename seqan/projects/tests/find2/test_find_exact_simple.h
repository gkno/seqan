/*
  Test the find2/find_exact_simple.h header.
*/
#ifndef TESTS_FIND2_TEST_FIND_EXACT_SIMPLE_H_
#define TESTS_FIND2_TEST_FIND_EXACT_SIMPLE_H_

#include "test_find_exact_tests.h"

using namespace seqan;

// Test the basic interface for the Simple pattern class.
SEQAN_DEFINE_TEST(test_find2_find_exact_simple_pattern_interface) {
    testFind2_ExactPattern_Interface<Simple>();
}


// "Easy" (= short, few hits) test of find() using the Simple pattern.
SEQAN_DEFINE_TEST(test_find2_find_exact_simple_pattern_find_easy) {
    testFind2_ExactPattern_FindEasy<Simple>();
}


// "Harder" (= longer, more hits hits) test of find() using the exact pattern.
SEQAN_DEFINE_TEST(test_find2_find_exact_simple_pattern_find_harder) {
    testFind2_ExactPattern_FindHarder<Simple>();
}


// "Easy" (= short, few hits) test of find() using the exact pattern.
SEQAN_DEFINE_TEST(test_find2_find_exact_simple_pattern_find_nomatch) {
    testFind2_ExactPattern_FindNoMatch<Simple>();
}


// Setting end position of exact Simple finder.
SEQAN_DEFINE_TEST(test_find2_find_exact_simple_pattern_set_end_position) {
    testFind2_ExactPattern_SetEndPosition<Simple>();
}


// Setting end position of exact Simple finder with long needle.
SEQAN_DEFINE_TEST(test_find2_find_exact_simple_pattern_set_end_position_long_needle) {
    testFind2_ExactPattern_SetEndPosition<Simple>();
}


// Test find() with long needle.
SEQAN_DEFINE_TEST(test_find2_find_exact_simple_pattern_find_long_needle) {
    testFind2_ExactPattern_FindLongNeedle<Simple>();
}


// Tests for setBeginPosition() with Simple pattern.
SEQAN_DEFINE_TEST(test_find2_find_exact_simple_pattern_set_begin_position) {
    testFind2_ExactPattern_SetBeginPosition<Simple>();
}


// Tests for setBeginPosition() with Simple pattern and a long
// needle.
SEQAN_DEFINE_TEST(test_find2_find_exact_simple_pattern_set_begin_position_long_needle) {
    testFind2_ExactPattern_SetBeginPositionLongNeedle<Simple>();
}

#endif  // TESTS_FIND2_TEST_FIND_EXACT_SIMPLE_H_
