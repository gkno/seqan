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

#include <iostream>
#include <cstdio>
#include <vector>

#define SEQAN_DEBUG
#define SEQAN_TEST

#include "seqan/basic.h"
#include "seqan/sequence.h"

#include "test_allocator.h"
#include "test_common.h"
#include "test_alphabet.h"


using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////
//test iterator proxy

SEQAN_DEFINE_TEST(Test_Proxy_Iterator) {
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

}


//////////////////////////////////////////////////////////////////////////////
// Helper class for constructor/destructor call counting.
// This is needed for test of Holder.

struct CtorDtorCounter
{
	static int static_ctors;
	static int static_dtors;

	int data_value;

	CtorDtorCounter():
		data_value(0)
	{
		++static_ctors;
	}
	CtorDtorCounter(CtorDtorCounter const & other_):
		data_value(other_.data_value)
	{
		++static_ctors;
	}
	~CtorDtorCounter()
	{
		++static_dtors;
	}
	CtorDtorCounter & operator = (CtorDtorCounter const & other_)
	{
		data_value = other_.data_value;
		return *this;
	}
};

int CtorDtorCounter::static_ctors = 0;
int CtorDtorCounter::static_dtors = 0;

//____________________________________________________________________________
// test for holder class

SEQAN_DEFINE_TEST(Test_Holder) {
	{
//ctors
		Holder<CtorDtorCounter> ho1;
		SEQAN_ASSERT(empty(ho1));
		SEQAN_ASSERT(!dependent(ho1));

		create(ho1);
		SEQAN_ASSERT(!empty(ho1));
		SEQAN_ASSERT(!dependent(ho1));

		Holder<CtorDtorCounter> ho2(ho1);

		Holder<CtorDtorCounter> ho3(value(ho1));
		SEQAN_ASSERT(!dependent(ho1));


//create
		CtorDtorCounter rco1;
		create(ho3, rco1);

//setValue
		setValue(ho3, rco1);
		SEQAN_ASSERT(dependent(ho3));

		rco1.data_value = 10;
		create(ho3);
		SEQAN_ASSERT_EQ(value(ho3).data_value, 10);

		CtorDtorCounter rco2;
		rco2.data_value = 20;

//operator = (value) => assignValue
		ho2 = rco2;
		SEQAN_ASSERT_EQ(value(ho2).data_value, 20);
		SEQAN_ASSERT(!dependent(ho2));

		rco2.data_value = 30;
		SEQAN_ASSERT_EQ(value(ho2).data_value, 20);

//operator = (holder) => assign
		setValue(ho1, rco1);
		ho1 = ho2;
		SEQAN_ASSERT_EQ(value(ho1).data_value, 20);

//clear
		clear(ho3);
		SEQAN_ASSERT(empty(ho3));

		assign(ho2, ho3);
		SEQAN_ASSERT(empty(ho2));

//conversion operator
		rco1 = ho1;
		SEQAN_ASSERT_EQ(rco1.data_value, 20);

//moveValue
		moveValue(ho1, rco2);
		SEQAN_ASSERT_EQ(rco1.data_value, 30);

	}

	SEQAN_ASSERT_EQ(CtorDtorCounter::static_ctors, CtorDtorCounter::static_dtors);


//test default implementations of addRef and releaseRef

//test const object holders
/*
	typedef char Bla[100];
	Holder<Bla const> cho1 = "test";*/
}

//////////////////////////////////////////////////////////////////////////////
// Test Iterator, basic functions

namespace SEQAN_NAMESPACE_MAIN
{

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

template <>
struct Value<Test_Iterator_1>
{
	typedef int Type;
};

}
//____________________________________________________________________________


SEQAN_DEFINE_TEST(Test_Iterator_Basic) {
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


SEQAN_DEFINE_TEST(Test_Iterator_Adaptor) {
	typedef AdaptorIterator<char *> TSpec;
	typedef Iter<char *, TSpec> TIterator;

	Test_Iter<TSpec>();

	char arr1[] = "abc";
	TIterator it1 = seqan::begin(arr1);
	char * ptr1 = it1;
	SEQAN_ASSERT_EQ(*ptr1, 'a');
}

SEQAN_DEFINE_TEST(Test_Iterator_Position) {
	Test_Iter<PositionIterator>();
}

//////////////////////////////////////////////////////////////////////////////
//class for testing move operations

struct MoveObj
{
	mutable int data_dat;

	MoveObj(int dat = 0): data_dat(dat) {}
	MoveObj(MoveObj const & other_): data_dat(other_.data_dat) 
	{ 
		other_.data_dat = 0; 
	}
	MoveObj const & operator = (MoveObj const & other_) const
	{ 
		data_dat = other_.data_dat; 
		other_.data_dat = 0; 
		return *this;
	}
	~MoveObj() {}
};

//____________________________________________________________________________

SEQAN_DEFINE_TEST(Test_Transport) {
	MoveObj o1(10);
	MoveObj o2;
	move(o2, o1);
	SEQAN_ASSERT_EQ(o1.data_dat, 0);
	SEQAN_ASSERT_EQ(o2.data_dat, 10);

	MoveObj const o3;
	move(o3, o2);
	SEQAN_ASSERT_EQ(o2.data_dat, 0);
	SEQAN_ASSERT_EQ(o3.data_dat, 10);

	MoveObj const o4;
	move(o4, o3);
	SEQAN_ASSERT_EQ(o3.data_dat, 0);
	SEQAN_ASSERT_EQ(o4.data_dat, 10);

	move(o1, o4);
	SEQAN_ASSERT_EQ(o4.data_dat, 0);
	SEQAN_ASSERT_EQ(o1.data_dat, 10);
}


SEQAN_DEFINE_TEST(test_basic_suprema_infima)
{
  using namespace seqan;

    // These tests are only here to instantiate the MaxValue and
    // MinValue Metafunctions for double and float.
    {
        double x = MaxValue<double>::VALUE;
        SEQAN_ASSERT_GT(x, 0);
    }
    {
        double x = MinValue<double>::VALUE;
        SEQAN_ASSERT_LT(x, 0);
    }
    {
        float x = MaxValue<float>::VALUE;
        SEQAN_ASSERT_GT(x, 0);
    }
    {
        float x = MinValue<float>::VALUE;
        SEQAN_ASSERT_LT(x, 0);
    }
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_BEGIN_TESTSUITE(test_basic) {
	SEQAN_CALL_TEST(test_basic_suprema_infima);
	SEQAN_CALL_TEST(Test_Proxy_Iterator);
	SEQAN_CALL_TEST(Test_Holder);
	SEQAN_CALL_TEST(Test_Iterator_Basic);
	SEQAN_CALL_TEST(Test_Iterator_Adaptor);
	SEQAN_CALL_TEST(Test_Iterator_Position);
	SEQAN_CALL_TEST(Test_Transport);

    // Tests from test_common.cpp.
    SEQAN_CALL_TEST(Test_Definition);
    SEQAN_CALL_TEST(Test_Type);
    SEQAN_CALL_TEST(Test_Iterator_Adapt_Std);
    // Tests from test_alphabet.cpp.
    SEQAN_CALL_TEST(TestAlphabetInterface);
    SEQAN_CALL_TEST(TestSimpleTypeConversions);
    SEQAN_CALL_TEST(TestExtremeValues);
    SEQAN_CALL_TEST(Test_Simple_Types);
    SEQAN_CALL_TEST(Test_Array_Functions);
    // Tests from test_allocator.cpp.
    SEQAN_CALL_TEST(testSimpleAllocator);
    SEQAN_CALL_TEST(testPoolAllocator);
    SEQAN_CALL_TEST(testMultiPoolAllocator);

    // Verify all check points.
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_aggregates.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_allocator_chunkpool.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_allocator_interface.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_allocator_multipool.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_allocator_simple.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_allocator_singlepool.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_allocator_to_std.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_alphabet_interface.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_alphabet_interface2.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_alphabet_simple.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_alphabet_simple_tabs.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_alphabet_trait_basic.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_compare.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_converter.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_debug.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_definition.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_forwards.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_holder.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_host.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_iterator.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_iterator_adapt_std.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_iterator_adaptor.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_iterator_base.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_iterator_position.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_logvalue.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_metaprogramming.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_pointer.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_profchar.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_profile.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_proxy.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_sse2.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_tag.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_testing.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_transport.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_type.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_volatile_ptr.h");
}
SEQAN_END_TESTSUITE

