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
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Tests for the SeqAn metaprogramming header.
// ==========================================================================

#ifndef TEST_BASIC_TEST_BASIC_METAPROGRAMMING_H_
#define TEST_BASIC_TEST_BASIC_METAPROGRAMMING_H_

#include <seqan/basic.h>

using namespace seqan;

SEQAN_DEFINE_TEST(test_basic_metaprogramming_true)
{
    bool b = True::VALUE;
    SEQAN_ASSERT_EQ(b, true);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_false)
{
    bool b = False::VALUE;
    SEQAN_ASSERT_EQ(b, false);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_eval)
{
    bool b = Eval<true>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    b = Eval<false>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, false);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_or)
{
    bool b = Or<True, True>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    b = Or<True, False>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    b = Or<False, True>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    b = Or<False, False>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, false);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_and)
{
    bool b = And<True, True>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    b = And<True, False>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, false);
    b = And<False, True>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, false);
    b = And<False, False>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, false);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_if)
{
    bool b = If<true, False, True>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, false);
    b = If<false, False, True>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, true);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_is_same_type)
{
    bool b = IsSameType<True, False>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, false);
    b = IsSameType<True, True>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, true);
}

// Helper functions and struct for the test for metaprogrammign Switch.

int switchTest(Nothing const &) { return -1; }
int switchTest(False const &) { return 0; }
int switchTest(True const &) { return 1; }
int switchTest(NilCase const &) { return 2; }

template <int X>
struct SwitchTest
{
    typedef typename Switch<
        X,
        Case<-1, Nothing,
        Case<0, False,
        Case<1, True
        > > > >::Type Type;
};

SEQAN_DEFINE_TEST(test_basic_metaprogramming_switch)
{
    using namespace seqan;

    typedef typename SwitchTest<-1>::Type T1;
    SEQAN_ASSERT_EQ(switchTest(T1()), -1);

    typedef typename SwitchTest<0>::Type T2;
    SEQAN_ASSERT_EQ(switchTest(T2()), 0);

    typedef typename SwitchTest<1>::Type T3;
    SEQAN_ASSERT_EQ(switchTest(T3()), 1);

    typedef typename SwitchTest<2>::Type T4;
    SEQAN_ASSERT_EQ(switchTest(T4()), 2);
}

struct LoopTestWorker_
{
    template <typename TArg>
    static inline void body(TArg & arg, int I)
    {
        arg[I - 1] += I;
    }
};

SEQAN_DEFINE_TEST(test_basic_metaprogramming_loop)
{
    int i[5] = { 1, 2, 3, 4, 5 };
    Loop<LoopTestWorker_, 5>::run(i);

    SEQAN_ASSERT_EQ(i[0], 2);
    SEQAN_ASSERT_EQ(i[1], 4);
    SEQAN_ASSERT_EQ(i[2], 6);
    SEQAN_ASSERT_EQ(i[3], 8);
    SEQAN_ASSERT_EQ(i[4], 10);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_loop_reverse)
{
    int i[5] = { 1, 2, 3, 4, 5 };
    LoopReverse<LoopTestWorker_, 5>::run(i);

    SEQAN_ASSERT_EQ(i[0], 2);
    SEQAN_ASSERT_EQ(i[1], 4);
    SEQAN_ASSERT_EQ(i[2], 6);
    SEQAN_ASSERT_EQ(i[3], 8);
    SEQAN_ASSERT_EQ(i[4], 10);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_log2)
{
    __uint64 x = Log2<1>::VALUE;
    SEQAN_ASSERT_EQ(x, 0u);
    x = Log2<2>::VALUE;
    SEQAN_ASSERT_EQ(x, 1u);
    x = Log2<3>::VALUE;
    SEQAN_ASSERT_EQ(x, 2u);
    x = Log2<4>::VALUE;
    SEQAN_ASSERT_EQ(x, 2u);
    x = Log2<5>::VALUE;
    SEQAN_ASSERT_EQ(x, 3u);
    x = Log2<15>::VALUE;
    SEQAN_ASSERT_EQ(x, 4u);
    x = Log2<16>::VALUE;
    SEQAN_ASSERT_EQ(x, 4u);
    x = Log2<17>::VALUE;
    SEQAN_ASSERT_EQ(x, 5u);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_log2_floor)
{
    __uint64 x = Log2Floor<1>::VALUE;
    SEQAN_ASSERT_EQ(x, 0u);
    x = Log2Floor<2>::VALUE;
    SEQAN_ASSERT_EQ(x, 1u);
    x = Log2Floor<3>::VALUE;
    SEQAN_ASSERT_EQ(x, 1u);
    x = Log2Floor<4>::VALUE;
    SEQAN_ASSERT_EQ(x, 2u);
    x = Log2Floor<5>::VALUE;
    SEQAN_ASSERT_EQ(x, 2u);
    x = Log2Floor<15>::VALUE;
    SEQAN_ASSERT_EQ(x, 3u);
    x = Log2Floor<16>::VALUE;
    SEQAN_ASSERT_EQ(x, 4u);
    x = Log2Floor<17>::VALUE;
    SEQAN_ASSERT_EQ(x, 4u);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_power)
{
    __uint64 x = Power<2, 2>::VALUE;
    SEQAN_ASSERT_EQ(x, 4u);
    x = Power<3, 2>::VALUE;
    SEQAN_ASSERT_EQ(x, 9u);
    x = Power<10, 0>::VALUE;
    SEQAN_ASSERT_EQ(x, 1u);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_make_unsigned)
{
    bool b = IsSameType<unsigned int, typename MakeUnsigned_<int>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<__uint64, typename MakeUnsigned_<__int64>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_make_signed)
{
    bool b = IsSameType<int, typename MakeSigned_<unsigned>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<__int64, typename MakeSigned_<__uint64>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_make_remove_const)
{
    bool b = IsSameType<int, const int>::Type::VALUE;
    SEQAN_ASSERT_NOT(b);
    b = IsSameType<int, RemoveConst_<const int>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<int, RemoveConst_<int>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_make_copy_const)
{
    bool b = IsSameType<int, CopyConst_<int, int>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<int const, CopyConst_<int, int const>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<int const, CopyConst_<int const, int const>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<int const, CopyConst_<int const, int>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_make_is_const)
{
    bool b = IsConst_<int>::Type::VALUE;
    SEQAN_ASSERT_NOT(b);
    b = IsConst_<int const>::Type::VALUE;
    SEQAN_ASSERT(b);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_class_identifier)
{
    void * id1 = ClassIdentifier_<int>::getID();
    void * id2 = ClassIdentifier_<signed int>::getID();
    void * id3 = ClassIdentifier_<double>::getID();
    SEQAN_ASSERT_EQ(id1, id2);
    SEQAN_ASSERT_NEQ(id1, id3);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_memset)
{
    int i[] = { 1, 2, 3 };

    memset<3 * sizeof(int), '\0'>(i);

    SEQAN_ASSERT_EQ(i[0], 0);
    SEQAN_ASSERT_EQ(i[1], 0);
    SEQAN_ASSERT_EQ(i[2], 0);
}

#endif  // TEST_BASIC_TEST_BASIC_METAPROGRAMMING_H_