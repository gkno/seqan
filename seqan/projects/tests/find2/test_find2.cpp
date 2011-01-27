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
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/find2.h>

#include "test_find_approx_dpsearch.h"
#include "test_find_exact_shiftand.h"
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
    SEQAN_CALL_TEST(test_find2_find_exact_simple_pattern_find_long_needle);
    SEQAN_CALL_TEST(test_find2_find_exact_simple_pattern_interface);
    SEQAN_CALL_TEST(test_find2_find_exact_simple_pattern_set_end_position);
    SEQAN_CALL_TEST(test_find2_find_exact_simple_pattern_set_end_position_long_needle);
    SEQAN_CALL_TEST(test_find2_find_exact_simple_pattern_set_begin_position);
    SEQAN_CALL_TEST(test_find2_find_exact_simple_pattern_set_begin_position_long_needle);

    SEQAN_CALL_TEST(test_find2_find_exact_shiftand_pattern_find_easy);
    SEQAN_CALL_TEST(test_find2_find_exact_shiftand_pattern_find_harder);
    SEQAN_CALL_TEST(test_find2_find_exact_shiftand_pattern_find_nomatch);
    SEQAN_CALL_TEST(test_find2_find_exact_shiftand_pattern_find_long_needle);
    SEQAN_CALL_TEST(test_find2_find_exact_shiftand_pattern_interface);
    SEQAN_CALL_TEST(test_find2_find_exact_shiftand_pattern_set_end_position);
    SEQAN_CALL_TEST(test_find2_find_exact_shiftand_pattern_set_end_position_long_needle);
    SEQAN_CALL_TEST(test_find2_find_exact_shiftand_pattern_set_begin_position);
    SEQAN_CALL_TEST(test_find2_find_exact_shiftand_pattern_set_begin_position_long_needle);

    SEQAN_CALL_TEST(test_find2_find_multiple_exact_simple_pattern_find_easy);
    SEQAN_CALL_TEST(test_find2_find_multiple_exact_simple_pattern_find_harder);
    SEQAN_CALL_TEST(test_find2_find_multiple_exact_simple_pattern_find_nomatch);
    SEQAN_CALL_TEST(test_find2_find_multiple_exact_simple_pattern_interface);
    SEQAN_CALL_TEST(test_find2_find_multiple_exact_simple_pattern_set_end_position);
    SEQAN_CALL_TEST(test_find2_find_multiple_exact_simple_pattern_set_begin_position);

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
    SEQAN_CALL_TEST(test_find2_find_hamming_simple_pattern_set_begin_position_score_limit_1);

    // This needs some work...
//     SEQAN_CALL_TEST(test_find2_find_pattern_wild_shiftand_is_unsigned);
//     SEQAN_CALL_TEST(test_find2_find_pattern_wild_shiftand_is_valid);
//     SEQAN_CALL_TEST(test_find2_find_pattern_wild_shiftand_length_without_wildcards);
//     SEQAN_CALL_TEST(test_find2_find_pattern_wild_shiftand_get_character_class);
//     SEQAN_CALL_TEST(test_find2_find_pattern_wild_shiftand_pattern_find_easy);
//     SEQAN_CALL_TEST(test_find2_find_pattern_wild_shiftand_pattern_find_harder);
//     SEQAN_CALL_TEST(test_find2_find_pattern_wild_shiftand_pattern_find_nomatch);
//     SEQAN_CALL_TEST(test_find2_find_pattern_wild_shiftand_pattern_interface);
//     SEQAN_CALL_TEST(test_find2_find_pattern_wild_shiftand_pattern_set_end_position);

    std::cout << std::endl;
    std::cout << "Verifying Check Points" << std::endl;
    std::cout << "======================" << std::endl;
    std::cout << std::endl;

    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_approx_dpsearch.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_approx_find_begin.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_exact_shiftand.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_exact_simple.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_finder_default.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_hamming_simple.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_multiple_exact_simple.h");
//     SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_pattern_wild_shiftand.h");

    std::cout << std::endl;
}
SEQAN_END_TESTSUITE

