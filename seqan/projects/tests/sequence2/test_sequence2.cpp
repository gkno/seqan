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
  Tests for the sequence module.

  Despite the name "test_sequence2", the tests are for the original
  sequence module.  However, this test suite uses new style tests.
 ===========================================================================
*/

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/sequence.h>  // Include module under test.

#include "test_sequence_std_adaptions.h"

SEQAN_BEGIN_TESTSUITE(test_sequence) {
    // Test adaptions for std::string.
    SEQAN_CALL_TEST(test_sequence_adaptions_metafunctions_std_string);
    SEQAN_CALL_TEST(test_sequence_adaptions_iterators_std_string);
    SEQAN_CALL_TEST(test_sequence_adaptions_sequence_interface_std_string);
    SEQAN_CALL_TEST(test_sequence_adaptions_sequence_memory_std_string);

    // Test adaptions for std::list.
    SEQAN_CALL_TEST(test_sequence_adaptions_metafunctions_std_list);
    SEQAN_CALL_TEST(test_sequence_adaptions_iterators_std_list);
    SEQAN_CALL_TEST(test_sequence_adaptions_sequence_interface_std_list);

    // TODO(holtgrew): Rename to sequence_adapt_std_{string,list}.
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence/std_string.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence/adapt_std_list.h");
}
SEQAN_END_TESTSUITE

