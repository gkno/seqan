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
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ===========================================================================
  Tests for the seeds module.
  ===========================================================================
*/

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.

#include "test_seeds_seed_diagonal.h"
#include "test_seeds_seed_base.h"
#include "test_seeds_seed_chained.h"
#include "test_seeds_seed_simple.h"
#include "test_seeds_seed_set_base.h"

SEQAN_BEGIN_TESTSUITE(test_seeds) {
    // Call tests.
    SEQAN_CALL_TEST(test_seeds_seed_diagonal_constructors);
    SEQAN_CALL_TEST(test_seeds_seed_diagonal_metafunctions);

    SEQAN_CALL_TEST(test_seeds_seed_base_constructors_simple);
    SEQAN_CALL_TEST(test_seeds_seed_base_metafunctions_simple);
    SEQAN_CALL_TEST(test_seeds_seed_base_getters_setters_simple);
    SEQAN_CALL_TEST(test_seeds_seed_base_basic_functions_simple);
    SEQAN_CALL_TEST(test_seeds_seed_base_constructors_chained);
    SEQAN_CALL_TEST(test_seeds_seed_base_metafunctions_chained);
    SEQAN_CALL_TEST(test_seeds_seed_base_getters_setters_chained);
    SEQAN_CALL_TEST(test_seeds_seed_base_basic_functions_chained);

    SEQAN_CALL_TEST(test_seeds_seed_chained_metafunctions);
    SEQAN_CALL_TEST(test_seeds_seed_chained_append_diagonal);
    SEQAN_CALL_TEST(test_seeds_seed_chained_truncate_diagonals);
    SEQAN_CALL_TEST(test_seeds_seed_chained_iterators);

    SEQAN_CALL_TEST(test_seeds_seed_simple_constructors);
    SEQAN_CALL_TEST(test_seeds_seed_simple_setters);

    SEQAN_CALL_TEST(test_seeds_seed_set_base_container_functions_unordered);

    SEQAN_CALL_TEST(test_seeds_seed_set_base_container_add_seed_merge_left_merging_possible_no_threshold_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_container_add_seed_merge_right_merging_possible_no_threshold_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_container_add_seed_merge_left_merging_impossible_no_threshold_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_container_add_seed_merge_left_merging_possible_threshold_not_reached_length_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_container_add_seed_merge_left_merging_possible_threshold_reached_length_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_container_add_seed_merge_left_merging_possible_threshold_not_reached_scored_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_container_add_seed_merge_left_merging_possible_threshold_reached_scored_unordered);

    // Verify checkpoints.
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/seeds2/seeds_seed_diagonal.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/seeds2/seeds_seed_base.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/seeds2/seeds_seed_chained.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/seeds2/seeds_seed_simple.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/seeds2/seeds_seed_set_base.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/seeds2/seeds_seed_set_unordered.h");
}
SEQAN_END_TESTSUITE
