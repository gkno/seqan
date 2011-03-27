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

#ifndef TESTS_BASIC_TEST_BASIC_ITERATOR_H_
#define TESTS_BASIC_TEST_BASIC_ITERATOR_H_

struct Test_Iterator_1
{
	int data_dat_1;
	mutable int data_dat_2;

	Test_Iterator_1(int data_):
		data_dat_1(data_ + 1),
		data_dat_2(data_ + 2)
	{
	}
	Test_Iterator_1(Test_Iterator_1 const & other_):
		data_dat_1(other_.data_dat_1 + 10),
		data_dat_2(other_.data_dat_2 + 10)
	{
	}
	~Test_Iterator_1() 
	{
	}

	Test_Iterator_1 & operator = (Test_Iterator_1 const & other_)
	{
		data_dat_1 = other_.data_dat_1 + 20;
		data_dat_2 = other_.data_dat_2 + 20;
		return *this;
	}

	int & operator * ()
	{
		return data_dat_1;
	}
	int & operator * () const
	{
		return data_dat_2;
	}
};

namespace seqan {

template <>
struct Value<Test_Iterator_1>
{
	typedef int Type;
};

}  // namespace seqan

SEQAN_DEFINE_TEST(test_basic_iterator_basic) {
//test default iterator functions

	//value, getValue, operator *
	int i1 = 10;
	int * it1 = & i1;
	SEQAN_ASSERT_EQ(value(it1), 10);
	SEQAN_ASSERT_EQ(getValue(it1), 10);

	Test_Iterator_1 it2(10);
	SEQAN_ASSERT_EQ(value(it2), 11);
	SEQAN_ASSERT_EQ(getValue(it2), 11);
/*
	Test_Iterator_1 const it3(10);
	SEQAN_ASSERT_EQ(value(it3), 12);
	SEQAN_ASSERT_EQ(getValue(it3), 12);

	//assign to value reference
	value(it2) = 15;
	SEQAN_ASSERT_EQ(getValue(it2), 15);

	//moveValue
	moveValue(it2, 50);
	SEQAN_ASSERT_EQ(value(it2), 50);

	//defaults of some advanced functions
	SEQAN_ASSERT_EQ(position(it1), 0);
	SEQAN_ASSERT_EQ(container(it1), *it1);
*/
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
void Test_Iter()
{
	typedef Iter<char *, TSpec> TIterator;

	char arr1[] = "XYZ";
	char arr2[] = "abcdefg";

	TIterator it1 = seqan::begin(arr1);
	SEQAN_ASSERT_EQ(container(it1), arr1);
	SEQAN_ASSERT_EQ(*it1, 'X');

	setContainer(it1, arr2);
	SEQAN_ASSERT_EQ(*it1, 'a');

	TIterator it2(it1);
	SEQAN_ASSERT_EQ(*it2, 'a');

	TIterator it3;
	it3 = it2 + 1;
	SEQAN_ASSERT_EQ(*it3, 'b');

	++it3;
	SEQAN_ASSERT_EQ(*it3, 'c');

	it3++;
	SEQAN_ASSERT_EQ(*it3, 'd');

	SEQAN_ASSERT_EQ((int)position(it3), 3);

	setPosition(it3, 2);
	SEQAN_ASSERT_EQ(*it3, 'c');

	--it3;
	SEQAN_ASSERT_EQ(*it3, 'b');

	it3--;
	SEQAN_ASSERT_EQ(*it3, 'a');
	SEQAN_ASSERT(it3 == it2);
	SEQAN_ASSERT(!(it3 != it2));

	++it3;
	assignValue(it3, 'z');
	SEQAN_ASSERT_EQ(*it3, 'z');

	moveValue(it3, 'y');
	SEQAN_ASSERT_EQ(*it3, 'y');

	TIterator const it4 = it3;
	assignValue(it4, 'x');
	SEQAN_ASSERT_EQ(*it4, 'x');
	SEQAN_ASSERT_EQ(*it3, 'x');

	moveValue(it4, 'w');
	SEQAN_ASSERT_EQ(*it4, 'w');
	SEQAN_ASSERT_EQ(*it3, 'w');


	it3 = 1 + it4;
	SEQAN_ASSERT_EQ(*it3, 'c');

	it3 = it4 - 1;
	SEQAN_ASSERT_EQ(*it3, 'a');

	it3 += 4;
	SEQAN_ASSERT_EQ(*it3, 'e');

	it3 -= 2;
	SEQAN_ASSERT_EQ(*it3, 'c');

	SEQAN_ASSERT_EQ(it3 - it4, 1);

}

SEQAN_DEFINE_TEST(test_basic_iterator_adaptor)
{
	typedef AdaptorIterator<char *> TSpec;
	typedef Iter<char *, TSpec> TIterator;

	Test_Iter<TSpec>();

	char arr1[] = "abc";
	TIterator it1 = seqan::begin(arr1);
	char * ptr1 = it1;
	SEQAN_ASSERT_EQ(*ptr1, 'a');
}

SEQAN_DEFINE_TEST(test_basic_iterator_position)
{
	Test_Iter<PositionIterator>();
}

#endif  // #ifndef TESTS_BASIC_TEST_BASIC_ITERATOR_H_
