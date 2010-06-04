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
    // Call tests of the sequence journal.
    SEQAN_CALL_TEST(test_sequence_journal_host);
    SEQAN_CALL_TEST(test_sequence_clear);
    SEQAN_CALL_TEST(test_sequence_erase_position);
    SEQAN_CALL_TEST(test_sequence_erase_begin_end);
    SEQAN_CALL_TEST(test_sequence_insert);
    SEQAN_CALL_TEST(test_sequence_insert_value);
    SEQAN_CALL_TEST(test_sequence_assign_value);
    SEQAN_CALL_TEST(test_sequence_assign_infix);
    SEQAN_CALL_TEST(test_sequence_length);
    SEQAN_CALL_TEST(test_sequence_copy_constructor);
    SEQAN_CALL_TEST(test_sequence_begin_end_iterator);
    // TODO(holtgrew): Broken.
//     SEQAN_CALL_TEST(test_sequence_begin_end_const_iterator);
    SEQAN_CALL_TEST(test_sequence_journal_fuzzying);

    // Verify checkpoints.
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence_journal/journal_entry.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence_journal/journal_tree_unbalanced.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence_journal/journal_tree_unbalanced_iterator.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence_journal/journal_tree_unbalanced_node.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence_journal/sequence_journal.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence_journal/sequence_journal_iterator.h");
}
SEQAN_END_TESTSUITE
