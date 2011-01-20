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
// Author: Manuel Holtgrewe <manuel.holtgrew@fu-berlin.de>
// ==========================================================================
// Tests for the module parallel.
// ==========================================================================

#include <seqan/basic.h>     // Testing infrastructure.
#include <seqan/file.h>      // Required to print strings in tests.
#include <seqan/parallel.h>  // Header under test.

#if defined(_OPENMP)
#include <omp.h>
#endif  // #if defined(_OPENMP)

#include "test_parallel_atomic_primitives.h"
#include "test_parallel_atomic_misc.h"

SEQAN_BEGIN_TESTSUITE(test_parallel) {
#if defined(_OPENMP)
    // Set number of threads to >=2 so there actually is parallelism.
    if (omp_get_max_threads() < 2)
        omp_set_num_threads(2);
#endif  // #if defined(_OPENMP)

    // Tests for atomic primitives.
    SEQAN_CALL_TEST(test_parallel_atomic_inc);
    SEQAN_CALL_TEST(test_parallel_atomic_dec);
    SEQAN_CALL_TEST(test_parallel_atomic_add);
    SEQAN_CALL_TEST(test_parallel_atomic_or);
    SEQAN_CALL_TEST(test_parallel_atomic_xor);
    SEQAN_CALL_TEST(test_parallel_atomic_cas);

    // Tests for misc simpmle atomic operations.
    SEQAN_CALL_TEST(test_parallel_atomic_min);
    SEQAN_CALL_TEST(test_parallel_atomic_max);
}
SEQAN_END_TESTSUITE
