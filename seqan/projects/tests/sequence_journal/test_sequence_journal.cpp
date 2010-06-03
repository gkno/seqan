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
    // Call tests.
    SEQAN_CALL_TEST(test_sequence_journal_simple_demo);

    // Verify checkpoints.
//     SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence_journal/journal_node.h");
//     SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence_journal/journal_tree.h");
//     SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence_journal/journal_iterator.h");
//     SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence_journal/journal_st_node.h");
//     SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence_journal/journal_tree_iterator.h");
//     SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence_journal/string_journal.h");
}
SEQAN_END_TESTSUITE
