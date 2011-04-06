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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================

#ifndef TESTS_BASIC_TEST_BASIC_PROXY_H_
#define TESTS_BASIC_TEST_BASIC_PROXY_H_

// TODO(holtgrew): Use the helper struct from test construct/destruct.

SEQAN_DEFINE_TEST(test_basic_proxy_iterator)
{
	int i1[] = {10, 20, 30};
	int * pi1 = i1;
	Proxy<IteratorProxy<int *> > px(pi1);
	SEQAN_ASSERT_EQ(px, 10);

//assign
	px = 11;
	SEQAN_ASSERT_EQ(i1[0], 11);

	int i2 = 12;
	px = i2;
	SEQAN_ASSERT_EQ(i1[0], 12);

	int * pi2 = i1 + 1;
	Proxy<IteratorProxy<int *> > px2(pi2);
	px = px2;
	SEQAN_ASSERT_EQ(i1[0], 20);
	SEQAN_ASSERT_EQ(px, 20);

//copy ctor
	Proxy<IteratorProxy<int *> > px3(px2);
	SEQAN_ASSERT_EQ(px3, 20);


//assign
	char s1[100] = "";
	char * it1 = s1;
	Proxy<IteratorProxy<char *> > px4(it1);

	assign(px4, 'X');
	SEQAN_ASSERT_EQ(px4, 'X');

	char c1 = 'a';
	assign(px4, c1);
	SEQAN_ASSERT_EQ(px4, 'a');

   SEQAN_ASSERT_FAIL("Move me to other tests in this header!");
}

SEQAN_DEFINE_TEST(test_basic_proxy_iterator_constructors)
{
    SEQAN_ASSERT_FAIL("Write me!");
}

SEQAN_DEFINE_TEST(test_basic_proxy_iterator_assign)
{
    SEQAN_ASSERT_FAIL("Write me!");
}

SEQAN_DEFINE_TEST(test_basic_proxy_iterator_set)
{
    SEQAN_ASSERT_FAIL("Write me!");
}

SEQAN_DEFINE_TEST(test_basic_proxy_iterator_move)
{
    SEQAN_ASSERT_FAIL("Write me!");
}

SEQAN_DEFINE_TEST(test_basic_proxy_iterator_getValue)
{
    SEQAN_ASSERT_FAIL("Write me!");
}

SEQAN_DEFINE_TEST(test_basic_proxy_iterator_comparators)
{
    SEQAN_ASSERT_FAIL("Write me!");
}

SEQAN_DEFINE_TEST(test_basic_proxy_iterator_stream_read)
{
    SEQAN_ASSERT_FAIL("Write me!");
}

SEQAN_DEFINE_TEST(test_basic_proxy_iterator_stream_write)
{
    SEQAN_ASSERT_FAIL("Write me!");
}

#endif  // #ifndef TESTS_BASIC_TEST_BASIC_PROXY_H_
