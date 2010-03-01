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
	SEQAN_ASSERT_EQ(_ClassIdentifier<int>::getID(), _ClassIdentifier<int>::getID());;
	SEQAN_ASSERT(_ClassIdentifier<char>::getID() != _ClassIdentifier<int>::getID());
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(Test_Type) {
	int i;
	int const ci = 99;
	int a[10];

	_toParameter<int>(& i) = 10;
	SEQAN_ASSERT_EQ(i, 10);;

	*_toParameter<int *>(& i) = 20;
	SEQAN_ASSERT_EQ(i, 20);;

	_Pointer<int>::Type p1 = _toPointer(i);
	*p1 = 30;
	SEQAN_ASSERT_EQ(i, 30);;

	_Pointer<int *>::Type p2 = _toPointer(p1);
	*p2 = 40;
	SEQAN_ASSERT_EQ(i, 40);;

	_Pointer<int[10]>::Type p3 = _toPointer(a);
	p3[1] = 50;
	SEQAN_ASSERT_EQ(a[1], 50);;

	_Pointer<int const *>::Type p4 = _toPointer(ci);
	SEQAN_ASSERT_EQ(*p4, 99);;
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(Test_Iterator_Adapt_Std) {
//test SeqAn iterator to fulfill std iterator 

	typedef ::std::iterator_traits<Iterator<char *, Rooted>::Type>::value_type T1;
	bool b1 = _isSameType<T1, char>();
	SEQAN_ASSERT_TRUE(b1);
}
