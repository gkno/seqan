#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>
#include <time.h>

#define SEQAN_DEBUG
#define SEQAN_TEST

#define TEST_PATH "projects/tests/find/"
#define LIB_PATH "projects/library/seqan/find/"

#include "seqan/find.h"

using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

template <typename TAlgorithmSpec>
void Test_OnlineAlg()
{
	String<unsigned int> pos;

//____________________________________________________________________________
// Test1 - small needle

	String<char> haystack("Dies ist ein Haystack. Ja, das ist wirklich einer!");
	Finder<String<char> > finder(haystack);

	String<char> needle("ist");
	Pattern<String<char>, TAlgorithmSpec> pattern(needle);

	while (find(finder, pattern))
		append(pos,position(finder));

	SEQAN_TASSERT(pos[0] == 5);
	SEQAN_TASSERT(pos[1] == 31);

//____________________________________________________________________________
// Test2 - large needle

	haystack = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefgaabcdef";
	setHost(finder, haystack);
	clear(finder);

	needle = "abcdefghijklmnopqrstuvwxyzabcdefg";
	setHost(pattern, needle);

	clear(pos);
	while (find(finder, pattern))
		append(pos,position(finder));

	SEQAN_TASSERT(pos[0] == 0);
	SEQAN_TASSERT(pos[1] == 26);

//____________________________________________________________________________
// Test3 - different alphabet

	String<Dna> hstk = "aaaaaaacaa";
	Finder<String<Dna> > finderDna(hstk);

	String<Dna> ndl = "aa";
	setHost(pattern, ndl);

	clear(pos);
	while (find(finderDna, pattern))
		append(pos,position(finderDna));

	SEQAN_TASSERT(pos[0] == 0);
	SEQAN_TASSERT(pos[1] == 1);
	SEQAN_TASSERT(pos[2] == 2);
	SEQAN_TASSERT(pos[3] == 3);
	SEQAN_TASSERT(pos[4] == 4);
	SEQAN_TASSERT(pos[5] == 5);
	SEQAN_TASSERT(pos[6] == 8);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlgorithmSpec>
void Test_OnlineAlgMulti()
{
	String<unsigned int> pos;

//____________________________________________________________________________
// Test1 - Single keyword
	String<char> haystack("Dies ist ein Haystack. Ja, das ist wirklich einer!");
	Finder<String<char> > finder(haystack);

	typedef String<String<char> > TNeedle;
	TNeedle keywords;
	appendValue(keywords, String<char>("ist"));
	Pattern<TNeedle, TAlgorithmSpec> pattern(keywords);

	while (find(finder, pattern))
		append(pos,position(finder));

	SEQAN_TASSERT(pos[0] == 5);
	SEQAN_TASSERT(pos[1] == 31);

	haystack = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefgaabcdef";
	setHost(finder, haystack);
	clear(finder);

	clear(keywords);
	appendValue(keywords, String<char>("abcdefghijklmnopqrstuvwxyzabcdefg"));
	setHost(pattern, keywords);
	clear(pos);

	while (find(finder, pattern))
		append(pos,position(finder));

	SEQAN_TASSERT(pos[0] == 0);
	SEQAN_TASSERT(pos[1] == 26);


	String<Dna> hstk = "aaaaaaacaa";
	Finder<String<Dna> > finderDna(hstk);

	typedef String<String<Dna> > TDnaNeedle;
	TDnaNeedle dna_keywords;
	appendValue(dna_keywords, String<Dna>("aa"));
	setHost(pattern, dna_keywords);

	clear(pos);
	while (find(finderDna, pattern))
		append(pos,position(finderDna));

	SEQAN_TASSERT(pos[0] == 0);
	SEQAN_TASSERT(pos[1] == 1);
	SEQAN_TASSERT(pos[2] == 2);
	SEQAN_TASSERT(pos[3] == 3);
	SEQAN_TASSERT(pos[4] == 4);
	SEQAN_TASSERT(pos[5] == 5);
	SEQAN_TASSERT(pos[6] == 8);

//____________________________________________________________________________
// Test2 - Multiple keywords
	String<char> hst("annual_announce");
	Finder<String<char> > fd(hst);

	typedef String<String<char> > TN;
	TN kyw;
	appendValue(kyw, String<char>("announce"));
	appendValue(kyw, String<char>("annual"));
	appendValue(kyw, String<char>("annually"));
	Pattern<TN, TAlgorithmSpec> pt(kyw);

	String<unsigned int> finderPos;
	String<unsigned int> keywordIndex;
	while (find(fd, pt)) {
		append(finderPos,position(fd));
		append(keywordIndex,position(pt));
	}

	SEQAN_TASSERT(finderPos[0] == 0);
	SEQAN_TASSERT(keywordIndex[0] == 1);
	SEQAN_TASSERT(finderPos[1] == 7);
	SEQAN_TASSERT(keywordIndex[1] == 0);

	String<Dna> hstDna("AGATACGATATATAC");
	Finder<String<Dna> > fdDna(hstDna);

	typedef String<String<Dna> > TNDna;
	TNDna kywDna;
	appendValue(kywDna, String<Dna>("ATATATA"));
	appendValue(kywDna, String<Dna>("TATAT"));
	appendValue(kywDna, String<Dna>("ACGATAT"));
	Pattern<TNDna, TAlgorithmSpec> ptDna(kywDna);

	clear(finderPos);
	clear(keywordIndex);
	while (find(fdDna, ptDna)) {
		append(finderPos,position(fdDna));
		append(keywordIndex,position(ptDna));
	}

	SEQAN_TASSERT(finderPos[0] == 4);
	SEQAN_TASSERT(keywordIndex[0] == 2);
	SEQAN_TASSERT(finderPos[1] == 8);
	SEQAN_TASSERT(keywordIndex[1] == 1);
	SEQAN_TASSERT(finderPos[2] == 7);
	SEQAN_TASSERT(keywordIndex[2] == 0);
}



void Test_Various()
{
	//Horspool
	String<char> haystk("Dies ist ein Haystack. Ja, das ist wirklich einer!");
	String<char> ndl("ist");

//	Pattern<String<char>, Horspool> horsp(ndl);
	Iterator<String<char>, Rooted>::Type it = begin(haystk);
/*
	while (find(it, horsp))
		printf("%i\n", position(it));


	//Multi
*/
	String<String<char> > ndl_2;
	resize(ndl_2, 2);
	ndl_2[0] = "ist";
	ndl_2[1] = "ein";

/*
	Finder<String<char>, MultipatternFinder> finder;
	it = begin(haystk, Rooted());

	while (find(finder, it, ndl_2))
		printf("%i at position %i\n", needle(finder), position(it));
*/
	
	//Score

	String<char> ndl_3("des");

	Pattern<String<char>, SimpleScore> fnd(-1);
	setHost(fnd, ndl_3);
	it = begin(haystk, Rooted());

	while (find(it, fnd))
		printf("position %i\n", position(it));
}

//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	Test_OnlineAlg<Horspool>();	
	Test_OnlineAlg<ShiftAnd>();
	Test_OnlineAlg<ShiftOr>();
	Test_OnlineAlg<BndmAlgo>();
	Test_OnlineAlg<BomAlgo>();
	
	Test_OnlineAlgMulti<AhoCorasick>();

	Test_Various();

//	testMyersUkkonen("accagaatatggagatctagggatcca", "agata", -2);
//	testMyersUkkonen("actacctttatctatcatcggattcgcgatctctcgcgatcgatggcttcgagtacgtcacacagtgcatctagccggattcgcgatctctcgcgatcgatggcttcgtgtacgtcac", 
//		"cggattcgcgatctctcgcgatcgatggcttcgtgtacgtcacacagtgcatctagc", -1);

	
//	debug::verifyCheckpoints("projects/library/seqan/find/find_myers_ukkonen.h");

	debug::verifyCheckpoints("projects/library/seqan/find/find_horspool.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_base.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_shiftand.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_shiftor.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_bndm.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_bom.h");



	SEQAN_TREPORT("TEST END")
	return 0;
}
