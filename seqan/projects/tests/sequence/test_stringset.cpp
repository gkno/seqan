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
	StringSet<CharString, TSpec> set;

	resize(set, 3);
	set[0] = "Hallo ";
	set[1] = "schlauer ";
	set[2] = "Hamster!";

	SEQAN_TASSERT(length(set) == 3)

	SEQAN_TASSERT(isEqual(set[0], "Hallo "))
	SEQAN_TASSERT(isEqual(set[1], "schlauer "))
	SEQAN_TASSERT(isEqual(set[2], "Hamster!"))
/*
	// currently, this won't work for Owner<ConcatDirect<..> > StringSets
	// to fix it, we need to introduce Modifiers for Segments
	// which propagate their resize events to their StringSets
	resize(set[0], 9);
	infix(set[0], 6, 9) = "du ";
	SEQAN_TASSERT(isEqual(set[0], "Hallo du "))
*/
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

	SEQAN_TASSERT(length(set) == 3)

	CharString all = concat(set);

	SEQAN_TASSERT(concat(set)[10] == 'a');
	SEQAN_TASSERT(isEqual(set[0], "Hallo "))
	SEQAN_TASSERT(isEqual(set[1], "schlauer "))
	SEQAN_TASSERT(isEqual(set[2], "Hamster!"))
	SEQAN_TASSERT(isEqual(all, "Hallo schlauer Hamster!"))

	SEQAN_TASSERT(stringSetLimits(set)[0] == 0)
	SEQAN_TASSERT(stringSetLimits(set)[1] == 6)
	SEQAN_TASSERT(stringSetLimits(set)[2] == 15)
	SEQAN_TASSERT(stringSetLimits(set)[3] == 23)

	StringSet<CharString, TSpec> const &cset = set;
	
	all = concat(cset);
	SEQAN_TASSERT(concat(cset)[10] == 'a');
	SEQAN_TASSERT(isEqual(all, "Hallo schlauer Hamster!"))
}


//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet>
void Test_StringSetIdHolder() {
	typedef	typename Id<TStringSet>::Type TId;

	TStringSet str;
	String<char> bla("a");
	TId id0 = assignValueById(str, bla);
	SEQAN_TASSERT(id0 == 0)
	SEQAN_TASSERT(length(str) == 1)
	SEQAN_TASSERT(str[0] == "a")
	SEQAN_TASSERT(getValueById(str, id0) == "a")
	String<char> bla1("b");
	TId id1 = assignValueById(str, bla1);
	SEQAN_TASSERT(id1 == 1)
	SEQAN_TASSERT(str[1] == "b")
	SEQAN_TASSERT(length(str) == 2)
	SEQAN_TASSERT(getValueById(str, id1) == "b")
	String<char> bla2("c");
	TId id2 = assignValueById(str, bla2);
	SEQAN_TASSERT(id2 == 2)
	SEQAN_TASSERT(str[2] == "c")
	SEQAN_TASSERT(length(str) == 3)
	SEQAN_TASSERT(getValueById(str, id2) == "c")
	String<char> bla3("d");
	TId id3 = assignValueById(str, bla3);
	SEQAN_TASSERT(id3 == 3)
	SEQAN_TASSERT(str[3] == "d")
	SEQAN_TASSERT(length(str) == 4)
	SEQAN_TASSERT(getValueById(str, id3) == "d")
	removeValueById(str,id1);
	SEQAN_TASSERT(getValueById(str, id0) == "a")
	SEQAN_TASSERT(getValueById(str, id1) == "")
	SEQAN_TASSERT(getValueById(str, id2) == "c")
	SEQAN_TASSERT(getValueById(str, id3) == "d")
	SEQAN_TASSERT(length(str) == 3)
	removeValueById(str,id2);
	SEQAN_TASSERT(getValueById(str, id0) == "a")
	SEQAN_TASSERT(getValueById(str, id1) == "")
	SEQAN_TASSERT(getValueById(str, id2) == "")
	SEQAN_TASSERT(getValueById(str, id3) == "d")
	SEQAN_TASSERT(length(str) == 2)
	removeValueById(str,id3);
	SEQAN_TASSERT(getValueById(str, id0) == "a")
	SEQAN_TASSERT(getValueById(str, id1) == "")
	SEQAN_TASSERT(getValueById(str, id2) == "")
	SEQAN_TASSERT(getValueById(str, id3) == "")
	SEQAN_TASSERT(length(str) == 1)
	id1 = assignValueById(str, bla1);
	id2 = assignValueById(str, bla2);
	id3 = assignValueById(str, bla3);	
	SEQAN_TASSERT(getValueById(str, id0) == "a")
	SEQAN_TASSERT(getValueById(str, id1) == "b")
	SEQAN_TASSERT(getValueById(str, id2) == "c")
	SEQAN_TASSERT(getValueById(str, id3) == "d")
	SEQAN_TASSERT(length(str) == 4)
	TStringSet subSet;
	String<unsigned int> ids;
	appendValue(ids, id1);
	appendValue(ids, id3);
	subset(str, subSet, ids);
	SEQAN_TASSERT(getValueById(subSet, id0) == "")
	SEQAN_TASSERT(getValueById(subSet, id1) == "b")
	SEQAN_TASSERT(getValueById(subSet, id2) == "")
	SEQAN_TASSERT(getValueById(subSet, id3) == "d")
	SEQAN_TASSERT(subSet[0] == "b")
	SEQAN_TASSERT(subSet[1] == "d")
	SEQAN_TASSERT(length(subSet) == 2)
	TStringSet subSet2;
	String<unsigned int> ids2;
	appendValue(ids2, id3);
	subset(subSet, subSet2, ids2);
	SEQAN_TASSERT(getValueById(subSet2, id0) == "")
	SEQAN_TASSERT(getValueById(subSet2, id1) == "")
	SEQAN_TASSERT(getValueById(subSet2, id2) == "")
	SEQAN_TASSERT(getValueById(subSet2, id3) == "d")
	SEQAN_TASSERT(length(subSet2) == 1)
	clear(subSet);
	assignValueById(subSet,subSet2, id3);
	assignValueById(subSet,str, id0);
	SEQAN_TASSERT(valueById(subSet, id0) == "a")
	SEQAN_TASSERT(getValueById(subSet, id1) == "")
	SEQAN_TASSERT(getValueById(subSet, id2) == "")
	SEQAN_TASSERT(valueById(subSet, id3) == "d")
	SEQAN_TASSERT(length(subSet) == 2)
	String<char> bla5("f");
	assignValueById(subSet, bla5, id3);
	SEQAN_TASSERT(getValueById(subSet, id3) == "f")
	SEQAN_TASSERT(length(subSet) == 2)
	assignValueById(subSet, bla3, 20);
	SEQAN_TASSERT(getValueById(subSet, 20) == "d")
	SEQAN_TASSERT(length(subSet) == 3)
}

//____________________________________________________________________________


template <typename TSpec>
void Test_StringSet_IdHolder()
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

	SEQAN_TASSERT(length(set) == 3)

	CharString all = concat(set);

	SEQAN_TASSERT(concat(set)[10] == 'a');
	SEQAN_TASSERT(isEqual(set[0], "Hallo "))
	SEQAN_TASSERT(isEqual(set[1], "schlauer "))
	SEQAN_TASSERT(isEqual(set[2], "Hamster!"))
	SEQAN_TASSERT(isEqual(all, "Hallo schlauer Hamster!"))

	SEQAN_TASSERT(stringSetLimits(set)[0] == 0)
	SEQAN_TASSERT(stringSetLimits(set)[1] == 6)
	SEQAN_TASSERT(stringSetLimits(set)[2] == 15)
	SEQAN_TASSERT(stringSetLimits(set)[3] == 23)

	StringSet<CharString, TSpec> const &cset = set;
	
	all = concat(cset);
	SEQAN_TASSERT(concat(cset)[10] == 'a');
	SEQAN_TASSERT(isEqual(all, "Hallo schlauer Hamster!"))

	Test_StringSetIdHolder<StringSet<String<char>, TSpec> >();
}

//____________________________________________________________________________


int mainTestStringSet()
{
	SEQAN_TREPORT("TEST STRINGSET BEGIN")

	Test_StringSet< Owner<Default> >();
	Test_StringSet_Concat< Owner<Default> >();
	Test_StringSet_Concat< Owner<ConcatDirect<> > >();

	Test_StringSet_IdHolder< Dependent<Tight> >();
	Test_StringSet_IdHolder< Dependent<Generous> >();

	debug::verifyCheckpoints("projects/library/seqan/sequence/sequence_multiple.h");

	SEQAN_TREPORT("TEST STRINGSET END")

	return 0;
}
