#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>

#define SEQAN_DEBUG
#define SEQAN_TEST

#include "seqan/sequence.h"
#include "seqan/file.h"

using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////
//test sequence default interface:
//non-container objects are treated like containers of length 1


template <typename TSpec>
void Test_StringSet()
{	
	typedef StringSet<CharString, TSpec> TStringSet;
	TStringSet set;

	resize(set, 3);
	set[0] = "Hallo ";
	set[1] = "schlauer ";
	set[2] = "Hamster!";

	SEQAN_ASSERT_TRUE(length(set) == 3);

	SEQAN_ASSERT_TRUE(isEqual(set[0], "Hallo "));
	SEQAN_ASSERT_TRUE(isEqual(set[1], "schlauer "));
	SEQAN_ASSERT_TRUE(isEqual(set[2], "Hamster!"));
/*
	// currently, this won't work for Owner<ConcatDirect<..> > StringSets
	// to fix it, we need to introduce Modifiers for Segments
	// which propagate their resize events to their StringSets
	resize(set[0], 9);
	infix(set[0], 6, 9) = "du ";
	SEQAN_ASSERT_TRUE(isEqual(set[0], "Hallo du "));
*/

	//StringSet iterators
	typedef typename Iterator<TStringSet>::Type TIterator;
	int i = 0;
	for (TIterator it = begin(set); it != end(set); goNext(it))
	{
		SEQAN_ASSERT_TRUE(*it == set[i]);
		++i;
	}

	TIterator itBegin = begin(set);
	SEQAN_ASSERT_TRUE(atBegin(itBegin) == true);
	SEQAN_ASSERT_TRUE(atEnd(itBegin) == false);
	TIterator itEnd = end(set);
	SEQAN_ASSERT_TRUE(atBegin(itEnd) == false);
	SEQAN_ASSERT_TRUE(atEnd(itEnd) == true);
	SEQAN_ASSERT_TRUE(i == 3);
}

//____________________________________________________________________________


template <typename TSpec>
void Test_StringSet_Concat()
{
	StringSet<CharString, TSpec> set;

	CharString s1 = "Hallo ";

	appendValue(set, s1);
	appendValue(set, "schlauer ");
	appendValue(set, "Hamster!");

	SEQAN_ASSERT_TRUE(length(set) == 3);

	CharString all = concat(set);

	SEQAN_ASSERT_TRUE(concat(set)[10] == 'a');
	SEQAN_ASSERT_TRUE(isEqual(set[0], "Hallo "));
	SEQAN_ASSERT_TRUE(isEqual(set[1], "schlauer "));
	SEQAN_ASSERT_TRUE(isEqual(set[2], "Hamster!"));
	SEQAN_ASSERT_TRUE(isEqual(all, "Hallo schlauer Hamster!"));

	SEQAN_ASSERT_TRUE(stringSetLimits(set)[0] == 0);
	SEQAN_ASSERT_TRUE(stringSetLimits(set)[1] == 6);
	SEQAN_ASSERT_TRUE(stringSetLimits(set)[2] == 15);
	SEQAN_ASSERT_TRUE(stringSetLimits(set)[3] == 23);

	StringSet<CharString, TSpec> const &cset = set;
	
	all = concat(cset);
	SEQAN_ASSERT_TRUE(concat(cset)[10] == 'a');
	SEQAN_ASSERT_TRUE(isEqual(all, "Hallo schlauer Hamster!"));
}


//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet>
void Test_StringSetIdHolder() {
	typedef	typename Id<TStringSet>::Type TId;
	typedef StringSet<String<char>, Dependent<Tight> > TSetTight;
	typedef StringSet<String<char>, Dependent<Generous> > TSetGenerous;

	TStringSet str;
	String<char> bla("a");
	TId id0 = assignValueById(str, bla);
	SEQAN_ASSERT_TRUE(id0 == 0);
	SEQAN_ASSERT_TRUE(idToPosition(str, id0) == 0);
	SEQAN_ASSERT_TRUE(positionToId(str, 0) == id0);
	SEQAN_ASSERT_TRUE(length(str) == 1);
	SEQAN_ASSERT_TRUE(str[0] == "a");
	SEQAN_ASSERT_TRUE(getValueById(str, id0) == "a");
	String<char> bla1("b");
	TId id1 = assignValueById(str, bla1);
	SEQAN_ASSERT_TRUE(id1 == 1);
	SEQAN_ASSERT_TRUE(idToPosition(str, id1) == 1);
	SEQAN_ASSERT_TRUE(positionToId(str, 1) == id1);
	SEQAN_ASSERT_TRUE(str[1] == "b");
	SEQAN_ASSERT_TRUE(length(str) == 2);
	SEQAN_ASSERT_TRUE(getValueById(str, id1) == "b");
	String<char> bla2("c");
	TId id2 = assignValueById(str, bla2);
	SEQAN_ASSERT_TRUE(id2 == 2);
	SEQAN_ASSERT_TRUE(str[2] == "c");
	SEQAN_ASSERT_TRUE(length(str) == 3);
	SEQAN_ASSERT_TRUE(getValueById(str, id2) == "c");
	String<char> bla3("d");
	TId id3 = assignValueById(str, bla3);
	SEQAN_ASSERT_TRUE(id3 == 3);
	SEQAN_ASSERT_TRUE(str[3] == "d");
	SEQAN_ASSERT_TRUE(length(str) == 4);
	SEQAN_ASSERT_TRUE(getValueById(str, id3) == "d");
	removeValueById(str,id1);
	SEQAN_ASSERT_TRUE(getValueById(str, id0) == "a");
	SEQAN_ASSERT_TRUE(getValueById(str, id2) == "c");
	SEQAN_ASSERT_TRUE(getValueById(str, id3) == "d");
	if (TYPECMP<TStringSet, TSetTight>::VALUE)
	{
		SEQAN_ASSERT_TRUE(length(str) == 3);
	}
	else if (TYPECMP<TStringSet, TSetGenerous>::VALUE)
	{
		SEQAN_ASSERT_TRUE(length(str) == 4);
	}
	removeValueById(str,id2);
	SEQAN_ASSERT_TRUE(getValueById(str, id0) == "a");
	SEQAN_ASSERT_TRUE(getValueById(str, id3) == "d");
	if (TYPECMP<TStringSet, TSetTight>::VALUE)
	{
		SEQAN_ASSERT_TRUE(length(str) == 2);
	}
	else if (TYPECMP<TStringSet, TSetGenerous>::VALUE)
	{
		SEQAN_ASSERT_TRUE(length(str) == 4);
	}

	String<char> bla4("e");
	TId id4 = assignValueById(str, bla4, 100);
	SEQAN_ASSERT_TRUE(id4 == 100);
	SEQAN_ASSERT_TRUE(getValueById(str, id4) == "e");
	SEQAN_ASSERT_TRUE(getValueById(str, id0) == "a");
	SEQAN_ASSERT_TRUE(getValueById(str, id3) == "d");
	removeValueById(str,id3);
	SEQAN_ASSERT_TRUE(getValueById(str, id0) == "a");
	SEQAN_ASSERT_TRUE(getValueById(str, id4) == "e");
	if (TYPECMP<TStringSet, TSetTight>::VALUE)
	{
		SEQAN_ASSERT_TRUE(length(str) == 2);
	}
	else if (TYPECMP<TStringSet, TSetGenerous>::VALUE)
	{
		SEQAN_ASSERT_TRUE(length(str) == 101);
	}
	String<char> bla5("f");
	TId id5 = assignValueById(str, bla5); 
	SEQAN_ASSERT_TRUE(getValueById(str, id0) == "a");
	SEQAN_ASSERT_TRUE(getValueById(str, id4) == "e");
	SEQAN_ASSERT_TRUE(getValueById(str, id5) == "f");
	assignValueById(str, bla5, id4); 
	SEQAN_ASSERT_TRUE(getValueById(str, id0) == "a");
	SEQAN_ASSERT_TRUE(getValueById(str, id4) == "f");
	SEQAN_ASSERT_TRUE(getValueById(str, id5) == "f");
	removeValueById(str,id4);
	SEQAN_ASSERT_TRUE(getValueById(str, id0) == "a");
	SEQAN_ASSERT_TRUE(getValueById(str, id5) == "f");
	if (TYPECMP<TStringSet, TSetTight>::VALUE)
	{
		SEQAN_ASSERT_TRUE(length(str) == 2);
	}
	else if (TYPECMP<TStringSet, TSetGenerous>::VALUE)
	{
		SEQAN_ASSERT_TRUE(length(str) == 102);
	}
	clear(str);
	id1 = assignValueById(str, bla1);
	id2 = assignValueById(str, bla2);
	id3 = assignValueById(str, bla3);	
	SEQAN_ASSERT_TRUE(getValueById(str, id1) == "b");
	SEQAN_ASSERT_TRUE(getValueById(str, id2) == "c");
	SEQAN_ASSERT_TRUE(getValueById(str, id3) == "d");
	if (TYPECMP<TStringSet, TSetTight>::VALUE)
	{
		SEQAN_ASSERT_TRUE(length(str) == 3);
	}
	else if (TYPECMP<TStringSet, TSetGenerous>::VALUE)
	{
		SEQAN_ASSERT_TRUE(length(str) == 3);
	}
}

//____________________________________________________________________________


template <typename TSpec>
void Test_StringSet_Id()
{	
	StringSet<CharString, Owner<Default> > origin;
	StringSet<CharString, TSpec> set;

	resize(origin, 3);
	origin[0] = "Hallo ";
	origin[1] = "schlauer ";
	origin[2] = "Hamster!";

	appendValue(set, origin[0]);
	appendValue(set, origin[1]);
	appendValue(set, origin[2]);

	SEQAN_ASSERT_TRUE(length(set) == 3);

	CharString all = concat(set);

	SEQAN_ASSERT_TRUE(concat(set)[10] == 'a');
	SEQAN_ASSERT_TRUE(isEqual(set[0], "Hallo "));
	SEQAN_ASSERT_TRUE(isEqual(set[1], "schlauer "));
	SEQAN_ASSERT_TRUE(isEqual(set[2], "Hamster!"));
	SEQAN_ASSERT_TRUE(isEqual(all, "Hallo schlauer Hamster!"));

	SEQAN_ASSERT_TRUE(stringSetLimits(set)[0] == 0);
	SEQAN_ASSERT_TRUE(stringSetLimits(set)[1] == 6);
	SEQAN_ASSERT_TRUE(stringSetLimits(set)[2] == 15);
	SEQAN_ASSERT_TRUE(stringSetLimits(set)[3] == 23);

	StringSet<CharString, TSpec> const &cset = set;
	
	all = concat(cset);
	SEQAN_ASSERT_TRUE(concat(cset)[10] == 'a');
	SEQAN_ASSERT_TRUE(isEqual(all, "Hallo schlauer Hamster!"));
}

SEQAN_DEFINE_TEST(StringSet_Owner_Default) {
    SEQAN_CHECKPOINT;
	Test_StringSet< Owner<Default> >();
}

SEQAN_DEFINE_TEST(StringSet_Concat_Owner_Default) {
    SEQAN_CHECKPOINT;
	Test_StringSet_Concat< Owner<Default> >();
}

SEQAN_DEFINE_TEST(StringSet_Concat_Owner_ConcatDirect) {
    SEQAN_CHECKPOINT;
	Test_StringSet_Concat< Owner<ConcatDirect<> > >();
}

SEQAN_DEFINE_TEST(StringSet_Id_Dependent_Tight) {
    SEQAN_CHECKPOINT;
	Test_StringSet_Id< Dependent<Tight> >();
}

SEQAN_DEFINE_TEST(StringSet_Id_Dependent_Generous) {
    SEQAN_CHECKPOINT;
	Test_StringSet_Id< Dependent<Generous> >();
}

SEQAN_DEFINE_TEST(StringSetIdHolder_Char_Dependent_Tight) {
    SEQAN_CHECKPOINT;
	Test_StringSetIdHolder<StringSet<String<char>, Dependent<Tight> > >();
}

SEQAN_DEFINE_TEST(StringSetIdHolder_Char_Dependent_Generous) {
    SEQAN_CHECKPOINT;
	Test_StringSetIdHolder<StringSet<String<char>, Dependent<Generous> > >();
}
