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
  Tests for the sequence_journal module.
  ===========================================================================
*/

#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_sequence_journal.h"


SEQAN_BEGIN_TESTSUITE(test_sequence_journal) {
    // Call tests of the sequence journal with unbalanced tree journal.
    SEQAN_CALL_TEST(test_sequence_journal_unbalanced_tree_host);
    SEQAN_CALL_TEST(test_sequence_journal_unbalanced_tree_clear);
    SEQAN_CALL_TEST(test_sequence_journal_unbalanced_tree_erase_position);
    SEQAN_CALL_TEST(test_sequence_journal_unbalanced_tree_erase_begin_end);
    SEQAN_CALL_TEST(test_sequence_journal_unbalanced_tree_insert);
    SEQAN_CALL_TEST(test_sequence_journal_unbalanced_tree_insert_value);
    SEQAN_CALL_TEST(test_sequence_journal_unbalanced_tree_assign_value);
    SEQAN_CALL_TEST(test_sequence_journal_unbalanced_tree_assign_infix);
    SEQAN_CALL_TEST(test_sequence_journal_unbalanced_tree_length);
    SEQAN_CALL_TEST(test_sequence_journal_unbalanced_tree_copy_constructor);
    SEQAN_CALL_TEST(test_sequence_journal_unbalanced_tree_begin_end_iterator);
    SEQAN_CALL_TEST(test_sequence_journal_unbalanced_tree_begin_end_const_iterator);
    SEQAN_CALL_TEST(test_sequence_journal_unbalanced_tree_subscript_operator);
    SEQAN_CALL_TEST(test_sequence_journal_unbalanced_tree_fuzzying);

    // Call tests of the sequence journal with sorted array_ journals.
    SEQAN_CALL_TEST(test_sequence_journal_sorted_array_host);
    SEQAN_CALL_TEST(test_sequence_journal_sorted_array_clear);
    SEQAN_CALL_TEST(test_sequence_journal_sorted_array_erase_position);
    SEQAN_CALL_TEST(test_sequence_journal_sorted_array_erase_begin_end);
    SEQAN_CALL_TEST(test_sequence_journal_sorted_array_insert);
    SEQAN_CALL_TEST(test_sequence_journal_sorted_array_insert_value);
    SEQAN_CALL_TEST(test_sequence_journal_sorted_array_assign_value);
    SEQAN_CALL_TEST(test_sequence_journal_sorted_array_assign_infix);
    SEQAN_CALL_TEST(test_sequence_journal_sorted_array_length);
    SEQAN_CALL_TEST(test_sequence_journal_sorted_array_copy_constructor);
    SEQAN_CALL_TEST(test_sequence_journal_sorted_array_begin_end_iterator);
    SEQAN_CALL_TEST(test_sequence_journal_sorted_array_begin_end_const_iterator);
    SEQAN_CALL_TEST(test_sequence_journal_sorted_array_subscript_operator);
    SEQAN_CALL_TEST(test_sequence_journal_sorted_array_fuzzying);

    // Verify checkpoints.
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence_journal/journal_entry.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence_journal/journal_entries_unbalanced_tree.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence_journal/journal_entries_unbalanced_tree_iterator.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence_journal/journal_entries_unbalanced_tree_node.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence_journal/journal_entries_sorted_array.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence_journal/sequence_journal.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence_journal/sequence_journal_iterator.h");
}
SEQAN_END_TESTSUITE
