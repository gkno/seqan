#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>

#define SEQAN_DEBUG
#define SEQAN_TEST

#include "seqan/sequence.h"

using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(Infix)
{
    SEQAN_CHECKPOINT;
//____________________________________________________________________________
// infix

	Infix<String<char> >::Type infix_1;
	String<char> str_1 = "this is a string";
	setHost(infix_1, str_1);
	SEQAN_ASSERT_TRUE(length(infix_1) == 0);
	SEQAN_ASSERT_TRUE(id(infix_1) == id(str_1));

	setEnd(infix_1, end(str_1));
	SEQAN_ASSERT_TRUE(infix_1 == str_1);
	SEQAN_ASSERT_TRUE(length(infix_1) == length(infix_1));

	setEndPosition(infix_1, 9);
	SEQAN_ASSERT_TRUE(infix_1 == infix(str_1, 0, 9));

	Infix<String<char> >::Type infix_2(infix_1);
	SEQAN_ASSERT_TRUE(infix_2 == infix(str_1, 0, 9));
	SEQAN_ASSERT_TRUE(infix_2 == infix_1);
	SEQAN_ASSERT_TRUE(id(infix_1) == id(infix_2));

	setBeginPosition(infix_2, 5);
	SEQAN_ASSERT_TRUE(infix_2 == "is a");

	setBegin(infix_2, begin(str_1));
	SEQAN_ASSERT_TRUE(infix_2 == "this is a");

	Infix<String<char> >::Type infix_3(str_1);
	SEQAN_ASSERT_TRUE(infix_3 == str_1);
	SEQAN_ASSERT_TRUE(id(infix_3) == id(str_1));

	Infix<String<char> >::Type infix_4(str_1, 5, 9);
	SEQAN_ASSERT_TRUE(infix_4 == "is a");

	SEQAN_ASSERT_TRUE(capacity(infix_4) == capacity(str_1) - length(str_1) + length(infix_4));

	Infix<String<char> >::Type infix_5(str_1, begin(str_1), end(str_1));
	SEQAN_ASSERT_TRUE(infix_5 == str_1);

	SEQAN_ASSERT_TRUE(begin(infix_5) == begin(str_1));
	SEQAN_ASSERT_TRUE(beginPosition(infix_5) == 0);
	SEQAN_ASSERT_TRUE(end(infix_5) == end(str_1));
	SEQAN_ASSERT_TRUE(endPosition(infix_5) == length(str_1));

	SEQAN_ASSERT_TRUE(begin(infix(str_1, 0, length(str_1))) == begin(str_1));
	SEQAN_ASSERT_TRUE(end(infix(str_1, 0, length(str_1))) == end(str_1));
	SEQAN_ASSERT_TRUE(length(infix(str_1, 0, length(str_1))) == length(str_1));

	str_1 = "begin middle end";
	assign(infix(str_1, 6, 12),  "to");
	SEQAN_ASSERT_TRUE(str_1 == "begin to end");

	assign(infix(str_1, 6, 8), "the test", 14);
	SEQAN_ASSERT_TRUE(str_1 == "begin the test");

//	setEnd(infix_1);
//	SEQAN_ASSERT_TRUE(infix_1 == "");

//____________________________________________________________________________
// test infix iteration

	str_1 = "begin middle end";
	goBegin(infix_1);
	SEQAN_ASSERT_TRUE(infix_1 == "b");

	goBegin(infix_1, str_1);
	SEQAN_ASSERT_TRUE(infix_1 == "b");

	goPrevious(infix_1);
	SEQAN_ASSERT_TRUE(atBegin(infix_1));

	goEnd(infix_1);
	SEQAN_ASSERT_TRUE(infix_1 == str_1);

	goEnd(infix_1, str_1);
	SEQAN_ASSERT_TRUE(infix_1 == str_1);

	goNext(infix_1);
	SEQAN_ASSERT_TRUE(atEnd(infix_1));
//____________________________________________________________________________
// compare operators 

	str_1 = "hello";
	Infix<String<char> >::Type infix_6(str_1, 0, 5);

	SEQAN_ASSERT_TRUE(infix_6 == str_1);

	infix_6 += str_1;
	SEQAN_ASSERT_TRUE(isEqual(infix_6, "hellohello"));

	SEQAN_ASSERT_TRUE(infix_6 != "bla");
	SEQAN_ASSERT_TRUE(!isNotEqual(infix_6, "hellohello"));

	SEQAN_ASSERT_TRUE(!(infix_6 < "hello"));
	SEQAN_ASSERT_TRUE(!isLess(infix_6, "hello"));

	SEQAN_ASSERT_TRUE(!(infix_6 <= "hello"));
	SEQAN_ASSERT_TRUE(!isLessOrEqual(infix_6, "hello"));

	SEQAN_ASSERT_TRUE(infix_6 > "hello");
	SEQAN_ASSERT_TRUE(isGreater(infix_6, "hello"));

	SEQAN_ASSERT_TRUE(infix_6 >= "hello");
	SEQAN_ASSERT_TRUE(isGreaterOrEqual(infix_6, "hello"));
//____________________________________________________________________________

	clear(infix_6);
	SEQAN_ASSERT_TRUE(infix_6 == "");
//____________________________________________________________________________
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(Suffix)
{
    SEQAN_CHECKPOINT;
//____________________________________________________________________________
// suffix

	String<char> str_1 = "this is a string";

	Suffix<String<char> >::Type suffix_1;
	setHost(suffix_1, str_1);
	SEQAN_ASSERT_TRUE(length(suffix_1) == length(str_1));
	SEQAN_ASSERT_TRUE(id(suffix_1) == id(str_1));

	Suffix<String<char> >::Type suffix_2(suffix_1);
	SEQAN_ASSERT_TRUE(suffix_2 == suffix_1);
	SEQAN_ASSERT_TRUE(id(suffix_1) == id(suffix_2));

	setBeginPosition(suffix_2, 5);
	SEQAN_ASSERT_TRUE(suffix_2 == "is a string");

	setBegin(suffix_2, begin(str_1));
	SEQAN_ASSERT_TRUE(suffix_2 == "this is a string");

	Suffix<String<char> >::Type suffix_3(str_1);
	SEQAN_ASSERT_TRUE(suffix_3 == str_1);
	SEQAN_ASSERT_TRUE(id(suffix_3) == id(str_1));

	Suffix<String<char> >::Type suffix_4(str_1, 5);
	SEQAN_ASSERT_TRUE(suffix_4 == "is a string");

	SEQAN_ASSERT_TRUE(capacity(suffix_4) == capacity(str_1) - length(str_1) + length(suffix_4));

	Suffix<String<char> >::Type suffix_5(str_1, begin(str_1));
	SEQAN_ASSERT_TRUE(suffix_5 == str_1);

	SEQAN_ASSERT_TRUE(begin(suffix_5) == begin(str_1));
	SEQAN_ASSERT_TRUE(beginPosition(suffix_5) == 0);
	SEQAN_ASSERT_TRUE(end(suffix_5) == end(str_1));
	SEQAN_ASSERT_TRUE(endPosition(suffix_5) == length(str_1));

	SEQAN_ASSERT_TRUE(begin(suffix(str_1, 0)) == begin(str_1));
	SEQAN_ASSERT_TRUE(end(suffix(str_1, 3)) == end(str_1));
	SEQAN_ASSERT_TRUE(length(suffix(str_1, 0)) == length(str_1));

	str_1 = "begin middle end";
	assign(suffix(str_1, 6), "to panic");
	SEQAN_ASSERT_TRUE(str_1 == "begin to panic");

	assign(suffix(str_1, 6), "the test", 9);
	SEQAN_ASSERT_TRUE(str_1 == "begin the");

	char str_2[200] = "begin middle end";
	assign(suffix(str_2, 6), "to panic");
	SEQAN_ASSERT_TRUE(isEqual(str_2, "begin to panic"));

	assign(suffix(str_2, 6), "the test", 9);
	SEQAN_ASSERT_TRUE(isEqual(str_2, "begin the"));

//____________________________________________________________________________
// test suffix iteration
/*
	str_1 = "begin middle end";
	goBegin(suffix_1);
	SEQAN_ASSERT_TRUE(suffix_1 == str_1);

	goBegin(suffix_1, str_1);
	SEQAN_ASSERT_TRUE(suffix_1 == str_1);
	SEQAN_ASSERT_TRUE(atBegin(suffix_1));

	goEnd(suffix_1);
	SEQAN_ASSERT_TRUE(suffix_1 == "d");

	goEnd(suffix_1, str_1);
	SEQAN_ASSERT_TRUE(suffix_1 == "d");

	goNext(suffix_1);
	SEQAN_ASSERT_TRUE(atEnd(suffix_1));

	goPrevious(suffix_1);
	SEQAN_ASSERT_TRUE(atEnd(suffix_1));
*/

}

SEQAN_DEFINE_TEST(Ticket317)
{
    SEQAN_CHECKPOINT;
    // http://trac.mi.fu-berlin.de/seqan/ticket/317
    
    CharString text = "thisisatext";
    Infix<CharString>::Type sub1 = infix(text, 2, 8);
    Infix<Infix<CharString>::Type>::Type sub2 = infix(sub1, 1, 3);

    // Check our assumptions.
    SEQAN_ASSERT_EQ(sub1, "isisat");
    SEQAN_ASSERT_EQ(sub2, "si");

    // operator - (triggered the ticket)
    SEQAN_ASSERT_EQ(begin(sub1) - begin(text), 2u);
    SEQAN_ASSERT_EQ(begin(sub2) - begin(text), 3u);
    SEQAN_ASSERT_EQ(begin(sub2) - begin(sub1), 1u);

    // operator == and !=
#   define TEST_IT(a, b) \
    SEQAN_ASSERT_EQ(a, b); \
    SEQAN_ASSERT_LEQ(a, b); \
    SEQAN_ASSERT_GEQ(a, b);

    TEST_IT(begin(sub1), begin(text) + 2);
    TEST_IT(begin(sub2), begin(text) + 3);
    TEST_IT(begin(sub2), begin(sub1) + 1);
#   undef TEST_IT

    SEQAN_ASSERT_NEQ(begin(sub1), begin(text));
    SEQAN_ASSERT_NEQ(begin(sub2), begin(text));
    SEQAN_ASSERT_NEQ(begin(sub2), begin(sub1));

    // operator <, >, <=, >=
#   define TEST_IT(a, b) \
    SEQAN_ASSERT_LT(a, b); \
    SEQAN_ASSERT_GT(b, a); \
    SEQAN_ASSERT_LEQ(a, b); \
    SEQAN_ASSERT_GEQ(b, a);

    TEST_IT(begin(sub1), end(text));
    TEST_IT(begin(text), begin(sub1));
    TEST_IT(begin(sub2), end(sub1));
    TEST_IT(begin(sub1), begin(sub2));
    TEST_IT(begin(sub2), end(text));
    TEST_IT(begin(text), begin(sub2));
#   undef TEST_IT
}
