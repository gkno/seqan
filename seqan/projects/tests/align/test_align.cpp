/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
  ===========================================================================
  Copyright (C) 2010
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.
  
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  
  ===========================================================================
  Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
  ===========================================================================
  Tests for the align-module.
  ===========================================================================
*/

#include <seqan/basic.h>  
#include <seqan/file.h>   

#include "test_align_align.h"
#include "test_align_gaps.h"
#include "test_align_myers.h"
#include "test_local_align.h"

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
    SEQAN_CALL_TEST(test_align_gaps_test_trailing_gaps_char_string_array_gaps);
    SEQAN_CALL_TEST(test_align_gaps_test_count_characters_char_string_array_gaps);

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
