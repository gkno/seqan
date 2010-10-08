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
//helper class for reference counting
//this is needed for test of Holder

struct RefCountObj
{
	static int static_addrefs;
	static int static_releaserefs;
	static int static_ctors;
	static int static_dtors;

	mutable int data_addrefs;
	mutable int data_releaserefs;

	int data_value;

	RefCountObj():
		data_addrefs(0),
		data_releaserefs(0),
		data_value(0)
	{
		++static_ctors;
	}
	RefCountObj(RefCountObj const & other_):
		data_addrefs(0),
		data_releaserefs(0),
		data_value(other_.data_value)
	{
		++static_ctors;
	}
	~RefCountObj()
	{
SEQAN_ASSERT_EQ(data_addrefs, data_releaserefs);

		++static_dtors;
	}
	RefCountObj & operator = (RefCountObj const & other_)
	{
		data_addrefs = data_releaserefs = 0;
		data_value = other_.data_value;
		return *this;
	}
};

int RefCountObj::static_addrefs = 0;
int RefCountObj::static_releaserefs = 0;
int RefCountObj::static_ctors = 0;
int RefCountObj::static_dtors = 0;

void addRef(RefCountObj & me)
{
	++me.data_addrefs;
	++RefCountObj::static_addrefs;
}
void addRef(RefCountObj const & me)
{
	++me.data_addrefs;
	++RefCountObj::static_addrefs;
}
void releaseRef(RefCountObj & me)
{
	++me.data_releaserefs;
	++RefCountObj::static_releaserefs;
}
void releaseRef(RefCountObj const & me)
{
	++me.data_releaserefs;
	++RefCountObj::static_releaserefs;
}

//____________________________________________________________________________
// test for holder class

SEQAN_DEFINE_TEST(Test_Holder) {
	{
//ctors
		Holder<RefCountObj> ho1;
		SEQAN_ASSERT_TRUE(empty(ho1));
		SEQAN_ASSERT_TRUE(!dependent(ho1));

		create(ho1);
		SEQAN_ASSERT_TRUE(!empty(ho1));
		SEQAN_ASSERT_TRUE(!dependent(ho1));

		Holder<RefCountObj> ho2(ho1);

		Holder<RefCountObj> ho3(value(ho1));
		SEQAN_ASSERT_EQ(value(ho1).data_addrefs, 1);
		SEQAN_ASSERT_EQ(value(ho1).data_releaserefs, 0);
		SEQAN_ASSERT_TRUE(!dependent(ho1));


//create
		RefCountObj rco1;
		create(ho3, rco1);
		SEQAN_ASSERT_EQ(value(ho3).data_addrefs, 0);
		SEQAN_ASSERT_EQ(value(ho1).data_addrefs, 1);
		SEQAN_ASSERT_EQ(value(ho1).data_releaserefs, 1);

//setValue
		setValue(ho3, rco1);
		SEQAN_ASSERT_TRUE(dependent(ho3));
		SEQAN_ASSERT_EQ(value(ho3).data_addrefs, 1);

		rco1.data_value = 10;
		create(ho3);
		SEQAN_ASSERT_EQ(value(ho3).data_value, 10);

		RefCountObj rco2;
		rco2.data_value = 20;

//operator = (value) => assignValue
		ho2 = rco2;
		SEQAN_ASSERT_EQ(value(ho2).data_value, 20);
		SEQAN_ASSERT_EQ(rco2.data_addrefs, 0);
		SEQAN_ASSERT_TRUE(!dependent(ho2));

		rco2.data_value = 30;
		SEQAN_ASSERT_EQ(value(ho2).data_value, 20);

//operator = (holder) => assign
		setValue(ho1, rco1);
		ho1 = ho2;
		SEQAN_ASSERT_EQ(value(ho1).data_value, 20);
		SEQAN_ASSERT_EQ(rco2.data_addrefs, 0);

//clear
		clear(ho3);
		SEQAN_ASSERT_TRUE(empty(ho3));

		assign(ho2, ho3);
		SEQAN_ASSERT_TRUE(empty(ho2));

//conversion operator
		rco1 = ho1;
		SEQAN_ASSERT_EQ(rco1.data_value, 20);

//moveValue
		moveValue(ho1, rco2);
		SEQAN_ASSERT_EQ(rco1.data_value, 30);

	}

	SEQAN_ASSERT_EQ(RefCountObj::static_addrefs, RefCountObj::static_releaserefs);
	SEQAN_ASSERT_EQ(RefCountObj::static_ctors, RefCountObj::static_dtors);


//test default implementations of addRef and releaseRef

	int i1;
	int const i2 = 0;

	addRef(i1);
	addRef(i2);
	releaseRef(i1);
	releaseRef(i2);

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

	TIterator it1 = begin(arr1);
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
	SEQAN_ASSERT_TRUE(it3 == it2);
	SEQAN_ASSERT_TRUE(!(it3 != it2));

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
	TIterator it1 = begin(arr1);
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

//////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew): This is highly temporary!
// Forward-declarations of tests defined elsewhere.
void SEQAN_TEST_Test_Definition();
void SEQAN_TEST_Test_Type();
void SEQAN_TEST_Test_Iterator_Adapt_Std();

void SEQAN_TEST_TestAlphabetInterface();
void SEQAN_TEST_TestSimpleTypeConversions();
void SEQAN_TEST_TestExtremeValues();
void SEQAN_TEST_Test_Simple_Types();
void SEQAN_TEST_Test_Array_Functions();

void SEQAN_TEST_testSimpleAllocator();
void SEQAN_TEST_testPoolAllocator();
void SEQAN_TEST_testMultiPoolAllocator();

SEQAN_DEFINE_TEST(test_basic_suprema_infima)
{
  using namespace seqan;

  {
    double x = SupremumValue<double>::VALUE;
    SEQAN_ASSERT_EQ(x, 1);
  }
  {
    double x = InfimumValue<double>::VALUE;
    SEQAN_ASSERT_EQ(x, 1);
  }
  {
    float x = SupremumValue<float>::VALUE;
    SEQAN_ASSERT_EQ(x, 1);
  }
  {
    float x = InfimumValue<float>::VALUE;
    SEQAN_ASSERT_EQ(x, 1);
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
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_counted_ptr.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_debug.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_definition.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_forwards.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_generated_forwards.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_holder.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_holder_dynamic.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_host.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_iterator.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_iterator_adapt_std.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_iterator_adaptor.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_iterator_base.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_iterator_position.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_iterator_simple.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_logvalue.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_metaprogramming.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/basic/basic_operator.h");
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

