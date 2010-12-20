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

#include <iostream>
#include <cstdio>
#include <vector>

#define SEQAN_DEBUG
#define SEQAN_TEST

#include "seqan/basic.h"
#include "seqan/sequence.h"

using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////


SEQAN_DEFINE_TEST(Test_Definition) {
	SEQAN_ASSERT_EQ(ClassIdentifier_<int>::getID(), ClassIdentifier_<int>::getID());
	SEQAN_ASSERT(ClassIdentifier_<char>::getID() != ClassIdentifier_<int>::getID());
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(Test_Type) {
	int i;
	int const ci = 99;
	int a[10];

	_toParameter<int>(& i) = 10;
	SEQAN_ASSERT_EQ(i, 10);

	*_toParameter<int *>(& i) = 20;
	SEQAN_ASSERT_EQ(i, 20);

	Pointer_<int>::Type p1 = _toPointer(i);
	*p1 = 30;
	SEQAN_ASSERT_EQ(i, 30);

	Pointer_<int *>::Type p2 = _toPointer(p1);
	*p2 = 40;
	SEQAN_ASSERT_EQ(i, 40);

	Pointer_<int[10]>::Type p3 = _toPointer(a);
	p3[1] = 50;
	SEQAN_ASSERT_EQ(a[1], 50);

	Pointer_<int const *>::Type p4 = _toPointer(ci);
	SEQAN_ASSERT_EQ(*p4, 99);
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(Test_Iterator_Adapt_Std) {
//test SeqAn iterator to fulfill std iterator 

	typedef ::std::iterator_traits<Iterator<char *, Rooted>::Type>::value_type T1;
	bool b1 = _isSameType<T1, char>();
	SEQAN_ASSERT_TRUE(b1);
}
