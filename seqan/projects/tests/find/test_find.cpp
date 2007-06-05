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

	SEQAN_TASSERT(host(pattern) == needle);
	SEQAN_TASSERT(host(reinterpret_cast<Pattern<String<char>, TAlgorithmSpec> const &>(pattern)) == needle);
	SEQAN_TASSERT(pos[0] == 5);
	SEQAN_TASSERT(pos[1] == 31);
	SEQAN_TASSERT(length(pos) == 2);

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
	SEQAN_TASSERT(length(pos) == 2);

//____________________________________________________________________________
// Test3 - different alphabet, small needle

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
	SEQAN_TASSERT(length(pos) == 7);

//____________________________________________________________________________
// Test3b - different alphabet, small needle, jumping finder

	goBegin(finderDna); // That's a repositioning
	clear(finderDna);	// That's why, clear state
	clear(pos);

	bool firstHit = true;
	while (find(finderDna, pattern)) {
		if (firstHit) {
			firstHit = false;
			finderDna += 2;
			clear(finderDna);  // clear the state of the finder
		} else {
			//unsigned int p = position(finderDna);
			append(pos,position(finderDna));
		}
	}

	SEQAN_TASSERT(pos[0] == 2);
	SEQAN_TASSERT(pos[1] == 3);
	SEQAN_TASSERT(pos[2] == 4);
	SEQAN_TASSERT(pos[3] == 5);
	SEQAN_TASSERT(pos[4] == 8);
	SEQAN_TASSERT(length(pos) == 5);

//____________________________________________________________________________
// Test4 - different alphabet, large needle
	String<Dna> text = "taaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaat";
	Finder<String<Dna> > finderText(text);

	String<Dna> query = "taaaataaaataaaataaaataaaataaaataaaataaaataaaat";
	setHost(pattern, query);

	clear(pos);
	while (find(finderText, pattern)) 
		append(pos,position(finderText));

	SEQAN_TASSERT(pos[0] == 0);
	SEQAN_TASSERT(pos[1] == 5);
	SEQAN_TASSERT(pos[2] == 10);
	SEQAN_TASSERT(pos[3] == 15);
	SEQAN_TASSERT(pos[4] == 20);
	SEQAN_TASSERT(pos[5] == 25);
	SEQAN_TASSERT(length(pos) == 6);
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

	SEQAN_TASSERT(host(pattern) == keywords);
	SEQAN_TASSERT(host(reinterpret_cast<Pattern<TNeedle, TAlgorithmSpec> const &>(pattern)) == keywords);
	SEQAN_TASSERT(pos[0] == 5);
	SEQAN_TASSERT(pos[1] == 31);
	SEQAN_TASSERT(length(pos) == 2);

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
	SEQAN_TASSERT(length(pos) == 2);


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
	SEQAN_TASSERT(length(pos) == 7);

	String<Dna> text = "taaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaat";
	Finder<String<Dna> > finderText(text);

	clear(dna_keywords);
	appendValue(dna_keywords, String<Dna>("taaaataaaataaaataaaataaaataaaataaaataaaataaaat"));
	setHost(pattern, dna_keywords);

	clear(pos);
	while (find(finderText, pattern)) 
		append(pos,position(finderText));

	SEQAN_TASSERT(pos[0] == 0);
	SEQAN_TASSERT(pos[1] == 5);
	SEQAN_TASSERT(pos[2] == 10);
	SEQAN_TASSERT(pos[3] == 15);
	SEQAN_TASSERT(pos[4] == 20);
	SEQAN_TASSERT(pos[5] == 25);
	SEQAN_TASSERT(length(pos) == 6);

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
	SEQAN_TASSERT(length(finderPos) == 2);
	SEQAN_TASSERT(length(keywordIndex) == 2);

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
	SEQAN_TASSERT(length(finderPos) == 3);
	SEQAN_TASSERT(length(keywordIndex) == 3);

//____________________________________________________________________________
// Test2 - Multiple keywords that do not fit into a machine word
	String<Dna> my_haystack("AGATACGATATATACAGATACGATATATACAGATACGATATATACAGATACGATATATACAGATACGATATATAC");
	Finder<String<Dna> > my_finder(my_haystack);

	typedef String<String<Dna> > TNeedle_My;
	TNeedle_My my_keywords;
	appendValue(my_keywords, String<Dna>("ATATATA"));
	appendValue(my_keywords, String<Dna>("ACCGATCCAT"));
	appendValue(my_keywords, String<Dna>("TATAT"));
	appendValue(my_keywords, String<Dna>("ACCGAT"));
	appendValue(my_keywords, String<Dna>("ACGATAT"));
	appendValue(my_keywords, String<Dna>("CCAA"));
	Pattern<TNeedle_My, TAlgorithmSpec> my_pattern(my_keywords);

	clear(finderPos);
	clear(keywordIndex);
	while (find(my_finder, my_pattern)) {
		//std::cout << position(my_finder) << "-" << position(my_pattern) << ::std::endl;
		append(finderPos,position(my_finder));
		append(keywordIndex,position(my_pattern));
	}

	SEQAN_TASSERT(finderPos[0] == 4);
	SEQAN_TASSERT(keywordIndex[0] == 4);
	SEQAN_TASSERT(finderPos[1] == 8);
	SEQAN_TASSERT(keywordIndex[1] == 2);
	SEQAN_TASSERT(finderPos[2] == 7);
	SEQAN_TASSERT(keywordIndex[2] == 0);
	SEQAN_TASSERT(finderPos[3] == 19);
	SEQAN_TASSERT(keywordIndex[3] == 4);
	SEQAN_TASSERT(finderPos[4] == 23);
	SEQAN_TASSERT(keywordIndex[4] == 2);
	SEQAN_TASSERT(finderPos[5] == 22);
	SEQAN_TASSERT(keywordIndex[5] == 0);
	SEQAN_TASSERT(finderPos[6] == 34);
	SEQAN_TASSERT(keywordIndex[6] == 4);
	SEQAN_TASSERT(finderPos[7] == 38);
	SEQAN_TASSERT(keywordIndex[7] == 2);
	SEQAN_TASSERT(finderPos[8] == 37);
	SEQAN_TASSERT(keywordIndex[8] == 0);
	SEQAN_TASSERT(finderPos[9] == 49);
	SEQAN_TASSERT(keywordIndex[9] == 4);
	SEQAN_TASSERT(finderPos[10] == 53);
	SEQAN_TASSERT(keywordIndex[10] == 2);
	SEQAN_TASSERT(finderPos[11] == 52);
	SEQAN_TASSERT(keywordIndex[11] == 0);
	SEQAN_TASSERT(finderPos[12] == 64);
	SEQAN_TASSERT(keywordIndex[12] == 4);
	SEQAN_TASSERT(finderPos[13] == 68);
	SEQAN_TASSERT(keywordIndex[13] == 2);
	SEQAN_TASSERT(finderPos[14] == 67);
	SEQAN_TASSERT(keywordIndex[14] == 0);


//____________________________________________________________________________
// Multiple keywords with overlapping matches
	String<Dna> my2_haystack("aaaacaaa");
	Finder<String<Dna> > my2_finder(my2_haystack);

	typedef String<String<Dna> > TNeedle_My2;
	TNeedle_My2 my2_keywords;
	appendValue(my2_keywords, String<Dna>("aa"));
	appendValue(my2_keywords, String<Dna>("aaa"));
	appendValue(my2_keywords, String<Dna>("ac"));
	appendValue(my2_keywords, String<Dna>("aac"));
	appendValue(my2_keywords, String<Dna>("gctccacctgacctagcccatggggcccaaatttccggccttaattcccattt"));
	Pattern<TNeedle_My2, TAlgorithmSpec> my2_pattern(my2_keywords);

	clear(finderPos);
	clear(keywordIndex);
	while (find(my2_finder, my2_pattern)) {
		//std::cout << position(my2_finder) << ":" << position(my2_pattern) << ::std::endl;
		append(finderPos,position(my2_finder));
		append(keywordIndex,position(my2_pattern));
	}

	SEQAN_TASSERT(finderPos[0] == 0);
	SEQAN_TASSERT(keywordIndex[0] == 0);
	SEQAN_TASSERT(finderPos[1] == 1);
	SEQAN_TASSERT(keywordIndex[1] == 0);
	SEQAN_TASSERT(finderPos[2] == 0);
	SEQAN_TASSERT(keywordIndex[2] == 1);
	SEQAN_TASSERT(finderPos[3] == 2);
	SEQAN_TASSERT(keywordIndex[3] == 0);
	SEQAN_TASSERT(finderPos[4] == 1);
	SEQAN_TASSERT(keywordIndex[4] == 1);
	SEQAN_TASSERT(finderPos[5] == 3);
	SEQAN_TASSERT(keywordIndex[5] == 2);
	SEQAN_TASSERT(finderPos[6] == 2);
	SEQAN_TASSERT(keywordIndex[6] == 3);
	SEQAN_TASSERT(finderPos[7] == 5);
	SEQAN_TASSERT(keywordIndex[7] == 0);
	SEQAN_TASSERT(finderPos[8] == 6);
	SEQAN_TASSERT(keywordIndex[8] == 0);
	SEQAN_TASSERT(finderPos[9] == 5);
	SEQAN_TASSERT(keywordIndex[9] == 1);

//____________________________________________________________________________
// Multiple duplicated keywords with overlapping matches, jumping finder
	goBegin(my2_finder); // That's a repositioning
	clear(my2_finder);	// That's why, clear state
	clear(finderPos);
	clear(keywordIndex);

	unsigned int hits = 0;
	while (find(my2_finder, my2_pattern)) {
		if (hits == 2) break;
		//std::cout << position(my2_finder) << ":" << position(my2_pattern) << ::std::endl;
		append(finderPos,position(my2_finder));
		append(keywordIndex,position(my2_pattern));
		++hits;
	}
	goBegin(my2_finder);
	my2_finder+=5;
	clear(my2_finder);
	clear(finderPos);
	clear(keywordIndex);
	while (find(my2_finder, my2_pattern)) {
		//std::cout << position(my2_finder) << ":" << position(my2_pattern) << ::std::endl;
		append(finderPos,position(my2_finder));
		append(keywordIndex,position(my2_pattern));
	}

	SEQAN_TASSERT(finderPos[0] == 5);
	SEQAN_TASSERT(keywordIndex[0] == 0);
	SEQAN_TASSERT(finderPos[1] == 6);
	SEQAN_TASSERT(keywordIndex[1] == 0);
	SEQAN_TASSERT(finderPos[2] == 5);
	SEQAN_TASSERT(keywordIndex[2] == 1);
}



void Test_Various()
{
	String<char> haystk("Dies ist ein Haystack. Ja, das ist wirklich einer!");
	String<char> ndl("des");

	Finder<String<char> > fnd(haystk);
	Pattern<String<char>, SimpleScore> pat(ndl, -1);

	while (find(fnd, pat))
		printf("-- position %i\n", position(fnd));
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
	Test_OnlineAlgMulti<MultipleShiftAnd>();
	Test_OnlineAlgMulti<SetHorspool>();
//	Test_OnlineAlgMulti<WuManber>();

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
	//debug::verifyCheckpoints("projects/library/seqan/find/find_quasar.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_ahocorasick.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_multiple_shiftand.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_set_horspool.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_wumanber.h");


	SEQAN_TREPORT("TEST END")
	return 0;
}
