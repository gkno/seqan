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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Tests for atomic primitives.
// ==========================================================================

#ifndef TEST_PARALLEL_TEST_PARALLEL_ATOMIC_PRIMITIVES_H_
#define TEST_PARALLEL_TEST_PARALLEL_ATOMIC_PRIMITIVES_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/parallel.h>
#include <seqan/random.h>

template <typename T>
void atomicIncTestImpl(T const &)
{
    using namespace seqan;

    T const ITERATIONS = 4 * 1024;

    T x = 0;
    #pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < static_cast<int>(ITERATIONS); ++i)
        atomicInc(x);

    SEQAN_ASSERT_EQ(x, ITERATIONS);
}

template <typename T>
void atomicDecTestImpl(T const &)
{
    using namespace seqan;

    T const ITERATIONS = 4 * 1024;

    T x = 2 * ITERATIONS;
    #pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < static_cast<int>(ITERATIONS); ++i)
        atomicDec(x);

    SEQAN_ASSERT_EQ(x, ITERATIONS);
}

template <typename T>
void atomicAddTestImpl(T const &)
{
    using namespace seqan;

    T const ITERATIONS = 4 * 1024;

    T x = 0;
    #pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < static_cast<int>(ITERATIONS); ++i)
        atomicAdd(x, 1);

    SEQAN_ASSERT_EQ(x, ITERATIONS);
}

template <typename T>
void atomicOrTestImpl(T const &)
{
    using namespace seqan;

    T const ITERATIONS = 4 * 1024;

    T expected = 0;
    for (T i = 0; i < ITERATIONS; ++i)
        expected |= i;
    
    T x = 0;
    #pragma omp parallel for schedule(static, 1)
    for (T i = 0; i < ITERATIONS; ++i)
        atomicOr(x, i);

    SEQAN_ASSERT_EQ(expected, x);
}

template <typename T>
void atomicXorTestImpl(T const &)
{
    using namespace seqan;

    T const ITERATIONS = 4 * 1024;

    T expected = 0;
    for (T i = 0; i < ITERATIONS; ++i)
        expected ^= i;
    
    T x = 0;
    #pragma omp parallel for schedule(static, 1)
    for (T i = 0; i < ITERATIONS; ++i)
        atomicXor(x, i);

    SEQAN_ASSERT_EQ(expected, x);
}

template <typename T>
void atomicCasTestImpl(T const &)
{
    using namespace seqan;

    // Re-implementation of a parallel max with test.

    int const ARR_SIZE = 4 * 1024;

    // Generate pseudorandom numbers to compute minimum of.
    T expectedMin = MaxValue<T>::VALUE;
    T arr[ARR_SIZE];
    T x = 3;
    for (int i = 0; i < ARR_SIZE; ++i) {
        x = 32456 * x + 9874;
        arr[i] = x;
        expectedMin = _min(expectedMin, arr[i]);
    }

    // Compute minimum.
    x = MaxValue<T>::VALUE;
    #pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < ARR_SIZE; ++i) {
        T val = x;
        while (val > arr[i]) {
            T m = _min(val, arr[i]);
            val = atomicCas(x, val, m);
        }
    }

    SEQAN_ASSERT_EQ(expectedMin, x);
}

SEQAN_DEFINE_TEST(test_parallel_atomic_inc)
{
    using namespace seqan;

    atomicIncTestImpl(__int16());
    atomicIncTestImpl(__int32());
    atomicIncTestImpl(__int64());
}

SEQAN_DEFINE_TEST(test_parallel_atomic_dec)
{
    using namespace seqan;

    atomicDecTestImpl(__int16());
    atomicDecTestImpl(__int32());
    atomicDecTestImpl(__int64());
}

SEQAN_DEFINE_TEST(test_parallel_atomic_add)
{
    using namespace seqan;

    atomicAddTestImpl(__int16());
    atomicAddTestImpl(__int32());
    atomicAddTestImpl(__int64());
}

SEQAN_DEFINE_TEST(test_parallel_atomic_or)
{
    using namespace seqan;

    atomicOrTestImpl(__int16());
    atomicOrTestImpl(__int32());
    atomicOrTestImpl(__int64());
}

SEQAN_DEFINE_TEST(test_parallel_atomic_xor)
{
    using namespace seqan;

    atomicXorTestImpl(__int16());
    atomicXorTestImpl(__int32());
    atomicXorTestImpl(__int64());
}

SEQAN_DEFINE_TEST(test_template_others_cas)
{
    using namespace seqan;

    atomicCasTestImpl(__int16());
    atomicCasTestImpl(__int32());
    atomicCasTestImpl(__int64());
}

#endif  // TEST_PARALLEL_TEST_PARALLEL_ATOMIC_PRIMITIVES_H_
