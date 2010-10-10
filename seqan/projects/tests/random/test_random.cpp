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
  Tests for SeqAn's module random.
  ===========================================================================
*/

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/random.h>  // The module under test.

#include "test_random_rng.h"
#include "test_random_dists.h"
#include "test_random_shuffle.h"


SEQAN_BEGIN_TESTSUITE(test_random) {
    // Call Tests.
    SEQAN_CALL_TEST(test_random_mt19937_constructors);
    SEQAN_CALL_TEST(test_random_mt19937_pick);

    SEQAN_CALL_TEST(test_random_normal_constructors);
    SEQAN_CALL_TEST(test_random_normal_pick);

    SEQAN_CALL_TEST(test_random_geometric_fair_coin_constructors);
    SEQAN_CALL_TEST(test_random_geometric_fair_coin_pick);

    SEQAN_CALL_TEST(test_random_lognormal_constructors);
    SEQAN_CALL_TEST(test_random_lognormal_pick);

    SEQAN_CALL_TEST(test_random_uniform_int_constructors);
    SEQAN_CALL_TEST(test_random_uniform_int_pick);

    SEQAN_CALL_TEST(test_random_uniform_double_constructors);
    SEQAN_CALL_TEST(test_random_uniform_double_pick);

    SEQAN_CALL_TEST(test_random_shuffle);

    // Verify Checkpoints.
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/random.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/random/random_lognormal.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/random/random_mt19937.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/random/random_normal.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/random/random_shuffle.h");
}
SEQAN_END_TESTSUITE
