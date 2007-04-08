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
	// currently, this won't work for ConcatDirect<..> StringSets
	// to fix it, we need to introduce Modifiers for Segments
	// which propagate their resize events to their StringSets
	resize(set[0], 9);
	infix(set[0], 6, 9) = "du ";
	SEQAN_TASSERT(isEqual(set[0], "Hallo du "))
*/

//____________________________________________________________________________

}

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

	StringSet<CharString, TSpec> const &cset = set;
	
	all = concat(cset);
	SEQAN_TASSERT(concat(cset)[10] == 'a');
	SEQAN_TASSERT(isEqual(all, "Hallo schlauer Hamster!"))

//____________________________________________________________________________

}

int mainTestStringSet()
{
	SEQAN_TREPORT("TEST STRINGSET BEGIN")

	Test_StringSet< ConcatVirtual<> >();
	Test_StringSet_Concat< ConcatVirtual<> >();
	Test_StringSet_Concat< ConcatDirect<> >();

	debug::verifyCheckpoints("projects/library/seqan/sequence/sequence_multiple.h");

	SEQAN_TREPORT("TEST STRINGSET END")

	return 0;
}
