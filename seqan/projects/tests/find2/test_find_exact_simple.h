// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Test the find2/find_exact_simple.h header.
// ==========================================================================

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
