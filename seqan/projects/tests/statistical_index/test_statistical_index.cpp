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
  Author: Jonathan Goeke <goeke@molgen.mpg.de>
  ===========================================================================
  Tests for SeqAn's module statistical_index (markov model)
  ===========================================================================
*/

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/statistical_index.h>  // The module under test.

#include "test_markov_model.h"
#include "test_statistics.h"


SEQAN_BEGIN_TESTSUITE(test_statistical_index) {
    // Call Tests.
    SEQAN_CALL_TEST(test_statistical_index_markov_model);
    SEQAN_CALL_TEST(test_statistical_index_statistics);

    // Verify Checkpoints.
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/statistical_index.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/statistical_index/statistical_index_markov_model.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/statistical_index/statistical_index_statistics.h");
}
SEQAN_END_TESTSUITE
