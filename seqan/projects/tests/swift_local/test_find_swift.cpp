#ifndef SEQAN_HEADER_TEST_FIND_SWIFT_H
#define SEQAN_HEADER_TEST_FIND_SWIFT_H

#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seeds.h>

#include "../../library/apps/swift_local/swift_local.h"

using namespace seqan;

// calls local swift and compares swift hits to expected swift hits
template<typename THaystack, typename TIndex>
void testLocalSwift(Finder<THaystack, Swift<SwiftLocal> > & finder,
					Pattern<TIndex, Swift<SwiftLocal> > & pattern,
					double epsilon,
					int minLength,
                    String<String<char> > & expectedPositions) {
	typedef typename Position<String<String<char> > >::Type TPosition;
	TPosition i = 0;
	while (find(finder, pattern, epsilon, minLength)) {
        // write swift hit into string and compare to expected output
        std::ostringstream pos;
		pos << "< " << positionRange(finder) << " , " << positionRange(pattern) << " >";
		SEQAN_ASSERT_EQ(pos.str().c_str(), value(expectedPositions,i));
		++i;
	}
}

void testOneLocalSwiftHit() {
	// a single pattern and a single hit 
    typedef Finder<DnaString, Swift<SwiftLocal> > TFinder;
	DnaString text = "aaaaaaacgatcgatgcaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
	TFinder finder_swift(text);

    typedef Index<DnaString, Index_QGram<UngappedShape<4> > > TQGramIndexSimple;
	TQGramIndexSimple index_4gram("tttcgatcgatgctttttttttttttttttttttttttttttt");
    Pattern<TQGramIndexSimple, Swift<SwiftLocal> > pattern_swift(index_4gram);

    String<String<char> > expectedPositions;
    appendValue(expectedPositions, "< < 7 , 17 > , < 0 , 17 > >");

	testLocalSwift(finder_swift, pattern_swift, 0.1, 6, expectedPositions);
}

void testOneLocalSwiftHit2() {
    // a single pattern and a single hit 2
    typedef Finder<DnaString, Swift<SwiftLocal> > TFinder;
	DnaString text = "aaaaagagaccccccagagaaaaaa";
	TFinder finder_swift(text);

    typedef Index<DnaString, Index_QGram<UngappedShape<4> > > TQGramIndexSimple;
	TQGramIndexSimple index_4gram("tttttggccccccggtttttt");
    Pattern<TQGramIndexSimple, Swift<SwiftLocal> > pattern_swift(index_4gram);

    String<String<char> > expectedPositions;
    appendValue(expectedPositions, "< < 9 , 15 > , < 0 , 15 > >");

	testLocalSwift(finder_swift, pattern_swift, 0.1, 6, expectedPositions);
}

void testOneLocalSwiftHitStringSet() {
	// a single hit, pattern in a StringSet
    typedef Finder<DnaString, Swift<SwiftLocal> > TFinder;
	DnaString text = "aaaaaaacgatcgatgcaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
	TFinder finder_swift(text);

	typedef Index<StringSet<DnaString>, Index_QGram< UngappedShape<4> > > TQGramIndex;
	TQGramIndex index_4gram;
	appendValue(indexText(index_4gram), "tttcgatcgatgctttttttttttttttttttttttttttttt");
    Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift(index_4gram);

    String<String<char> > expectedPositions;
    appendValue(expectedPositions, "< < 7 , 17 > , < < 0 , 0 > , < 0 , 17 > > >");

	testLocalSwift(finder_swift, pattern_swift, 0.1, 6, expectedPositions);
}

void testOneLocalSwiftHitBucketBorder() {
	// hit at bucket border
    typedef Finder<DnaString, Swift<SwiftLocal> > TFinder;
	DnaString text = "aaacgatcgatgcaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
	TFinder finder_swift(text);

	typedef Index<StringSet<DnaString>, Index_QGram< UngappedShape<4> > > TQGramIndex;
	TQGramIndex index_4gram;
	appendValue(indexText(index_4gram), "tttcgatcgatgctttttttttttttttttttttttttttttt");
    Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift(index_4gram);

    String<String<char> > expectedPositions;
    appendValue(expectedPositions, "< < 3 , 13 > , < < 0 , 0 > , < 0 , 13 > > >");
    appendValue(expectedPositions, "< < 3 , 13 > , < < 0 , 2 > , < 0 , 29 > > >");

	testLocalSwift(finder_swift, pattern_swift, 0.1, 6, expectedPositions);
}

void testOneLocalSwiftHitNegDiag() {

	// hit at lower haystack position -> diagonal start position negative
    typedef Finder<DnaString, Swift<SwiftLocal> > TFinder;
	DnaString text = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaacgatcgatgc";
	TFinder finder_swift(text);

	typedef Index<StringSet<DnaString>, Index_QGram< UngappedShape<4> > > TQGramIndex;
	TQGramIndex index_4gram;
	appendValue(indexText(index_4gram), "ttttcgatcgatgctttttttttttttttttttttttttttttt");
    Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift(index_4gram);

    String<String<char> > expectedPositions;
    appendValue(expectedPositions, "< < 36 , 46 > , < < 0 , 0 > , < 0 , 14 > > >");
    appendValue(expectedPositions, "< < 36 , 46 > , < < 0 , 3 > , < 0 , 30 > > >");

	testLocalSwift(finder_swift, pattern_swift, 0.1, 10, expectedPositions);
}

void testLocalSwiftTwoPatterns() {
	// two pattern sequences
    typedef Finder<DnaString, Swift<SwiftLocal> > TFinder;
    DnaString text = "aaaacgttccaaaaaaa";
	TFinder finder_swift(text);

	typedef Index<StringSet<DnaString>, Index_QGram< UngappedShape<4> > > TQGramIndex;
	TQGramIndex index_4gram;
	appendValue(indexText(index_4gram), "ttttcgttcctttttttttttttttttttttt"); // length = 32
    appendValue(indexText(index_4gram), "ttcgttcctt"); // length = 10
	Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift(index_4gram);

    String<String<char> > expectedPositions;
    appendValue(expectedPositions, "< < 4 , 10 > , < < 0 , 0 > , < 0 , 10 > > >");
    appendValue(expectedPositions, "< < 4 , 10 > , < < 0 , 3 > , < 0 , 26 > > >");
    appendValue(expectedPositions, "< < 4 , 10 > , < < 1 , 0 > , < 1 , 10 > > >");

	testLocalSwift(finder_swift, pattern_swift, 0.1, 6, expectedPositions);
}

void testLocalSwiftLongPatterns() {

	// two longer pattern sequences
    typedef Finder<DnaString, Swift<SwiftLocal> > TFinder;
	DnaString text = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaacgatcagtgacaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
	TFinder finder_swift(text);

	typedef Index<StringSet<DnaString>, Index_QGram< UngappedShape<4> > > TQGramIndex;
	TQGramIndex index_4gram;
	appendValue(indexText(index_4gram), "ttttttttttttttcgatcagtgacttttttttttttttttttttttttttttttttatcagt"); // length = 63
    appendValue(indexText(index_4gram), "ttttttttttttttcgatcagtgacttttttttttttttttttttttttttttttttatcagt"); // length = 63
	Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift(index_4gram);

    String<String<char> > expectedPositions;
    appendValue(expectedPositions, "< < 88 , 99 > , < < 0 , 7 > , < 0 , 35 > > >");
    appendValue(expectedPositions, "< < 88 , 99 > , < < 1 , 7 > , < 1 , 35 > > >");
    appendValue(expectedPositions, "< < 90 , 96 > , < < 0 , 41 > , < 0 , 63 > > >");
    appendValue(expectedPositions, "< < 90 , 96 > , < < 1 , 41 > , < 1 , 63 > > >");

	testLocalSwift(finder_swift, pattern_swift, 0.1, 6, expectedPositions);
}

SEQAN_DEFINE_TEST(test_find_swift) {
    testOneLocalSwiftHit();
    testOneLocalSwiftHit2();
    testOneLocalSwiftHitStringSet();
    testOneLocalSwiftHitBucketBorder();
    testOneLocalSwiftHitNegDiag();
    testLocalSwiftTwoPatterns();
    testLocalSwiftLongPatterns();
}

template<typename TString, typename TPos>
Align<TString>
testLongestEpsMatch(TString const & seq1, TString const & seq2, TPos seedBegin, TPos seedEnd) {
    Align<TString> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), seq1);
    assignSource(row(alignment, 1), seq2);

    Score<int> score(1, -1, -1);
    globalAlignment(alignment, score);

    longestEpsMatch(alignment, seedBegin, seedEnd, 6, 0.125);

    return alignment;
}

SEQAN_DEFINE_TEST(test_longest_epsMatch) {
    DnaString seq1 = "ACCTTTGCCCCCCCCCCTAAAAAAAATTAAAA";
    DnaString seq2 = "ACGTTTACCCCCCCCCCGAAAAAAAAGAAAA";
    Align<DnaString> alignment = testLongestEpsMatch(seq1, seq2, 7, 17);
    SEQAN_ASSERT_TRUE(row(alignment, 0) == "ACCTTTGCCCCCCCCCCTAAAAAAAA");
    SEQAN_ASSERT_TRUE(row(alignment, 1) == "ACGTTTACCCCCCCCCCGAAAAAAAA");

    seq1 = "ACCTTTGCCCCCCCCCCTAAAAAAAATTAAAA";
    seq2 = "ACGTTTCCCCCCCCCCGAAAAAAAAGAAAA";
    alignment = testLongestEpsMatch(seq1, seq2, 7, 17);
    SEQAN_ASSERT_TRUE(row(alignment, 0) == "ACCTTTGCCCCCCCCCCTAAAAAAAA");
    SEQAN_ASSERT_TRUE(row(alignment, 1) == "ACGTTT-CCCCCCCCCCGAAAAAAAA");
    
    seq1 = "AAAATTAAAAAAAATCCCCCCCCCCGTTTCCA";
    seq2 = "AAAAGAAAAAAAAGCCCCCCCCCCTTTGCA";
    alignment = testLongestEpsMatch(seq1, seq2, 15, 25);
    SEQAN_ASSERT_TRUE(row(alignment, 0) == "AAAAAAAATCCCCCCCCCCGTTTCCA");
    SEQAN_ASSERT_TRUE(row(alignment, 1) == "AAAAAAAAGCCCCCCCCCC-TTTGCA");
}

template<typename TString, typename TScore, typename TAliString, typename TScoreValue>
void testXDropAlign(TString const & seq1, TString const & seq2, TScore scoring, TScoreValue scoreDropOff, TScoreValue minScore, TAliString & aliString) {
    Align<TString> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq1);
    assignSource(row(align, 1), seq2);

    globalAlignment(align, scoring);
    //std::cout << align << std::endl;

    _splitAtXDrops(align, scoring, scoreDropOff, minScore, aliString);
}

SEQAN_DEFINE_TEST(test_split_xDrop_align) {
    Score<int> scoring(2,-1,-2);

    DnaString seq1 = "cgataagctcttggacta";
    DnaString seq2 = "cgataatatggactagg";
    String<Align<DnaString> > aliString;
    testXDropAlign(seq1, seq2, scoring, 3/*scoreDropOff*/, 10/*minScore*/, aliString);
    SEQAN_ASSERT_EQ(length(aliString), 2u);
    SEQAN_ASSERT_TRUE(row(value(aliString, 0), 0) == "cgataa");
    SEQAN_ASSERT_TRUE(row(value(aliString, 0), 1) == "cgataa");
    SEQAN_ASSERT_TRUE(row(value(aliString, 1), 0) == "tggacta");
    SEQAN_ASSERT_TRUE(row(value(aliString, 1), 1) == "tggacta");

    seq1 = "cgataagctcttggacta";
    seq2 = "cgataatatggactagggg";
    clear(aliString);
    testXDropAlign(seq1, seq2, scoring, 3/*scoreDropOff*/, 14/*minScore*/, aliString);
    SEQAN_ASSERT_EQ(length(aliString), 1u);
    SEQAN_ASSERT_TRUE(row(value(aliString, 0), 0) == "tggacta");
    SEQAN_ASSERT_TRUE(row(value(aliString, 0), 1) == "tggacta");

    seq1 = "cgataagctcagttggacta";
    seq2 = "cgataatcactggactagggg";
    clear(aliString);
    testXDropAlign(seq1, seq2, scoring, 2/*scoreDropOff*/, 6/*minScore*/, aliString);
    SEQAN_ASSERT_EQ(length(aliString), 3u);
    SEQAN_ASSERT_TRUE(row(value(aliString, 0), 0) == "cgataa");
    SEQAN_ASSERT_TRUE(row(value(aliString, 0), 1) == "cgataa");
    SEQAN_ASSERT_TRUE(row(value(aliString, 1), 0) == "tca");
    SEQAN_ASSERT_TRUE(row(value(aliString, 1), 1) == "tca");
    SEQAN_ASSERT_TRUE(row(value(aliString, 2), 0) == "tggacta");
    SEQAN_ASSERT_TRUE(row(value(aliString, 2), 1) == "tggacta");

    clear(aliString);
    testXDropAlign(seq1, seq2, scoring, 3/*scoreDropOff*/, 6/*minScore*/, aliString);
    SEQAN_ASSERT_EQ(length(aliString), 2u);
    SEQAN_ASSERT_TRUE(row(value(aliString, 0), 0) == "cgataa");
    SEQAN_ASSERT_TRUE(row(value(aliString, 0), 1) == "cgataa");
    SEQAN_ASSERT_TRUE(row(value(aliString, 1), 0) == "tcagttggacta");
    SEQAN_TASSERT(row(value(aliString, 2), 1) == "tca-ctggacta"); // TODO: Warum geht SEQAN_ASSERT_TRUE hier nicht?    
    
    seq1 = "aaaaaa";
    seq2 = "ccaaaaaa";
    clear(aliString);
    testXDropAlign(seq1, seq2, scoring, 3/*scoreDropOff*/, 12/*minScore*/, aliString);
    SEQAN_ASSERT_EQ(length(aliString), 1u);
    SEQAN_ASSERT_TRUE(row(value(aliString, 0), 0) == "aaaaaa");
    SEQAN_ASSERT_TRUE(row(value(aliString, 0), 1) == "aaaaaa"); 

    seq1 = "aaaaaa";
    seq2 = "aaaaaa";
    clear(aliString);
    testXDropAlign(seq1, seq2, scoring, 3/*scoreDropOff*/, 12/*minScore*/, aliString);
    SEQAN_ASSERT_EQ(length(aliString), 1u);
    SEQAN_ASSERT_TRUE(row(value(aliString, 0), 0) == "aaaaaa");
    SEQAN_ASSERT_TRUE(row(value(aliString, 0), 1) == "aaaaaa");
    
    seq1 = "CCCCAGGGGGGACAAAAAAGAAACCCAGGGGGGACCCAGGG";
    seq2 = "CCCCTGGGGGGTCTTTTTGTCCCTGGGGGGTCCCTGGG";
    clear(aliString);
    testXDropAlign(seq1, seq2, Score<int>(1, -1, -1), 3/*scoreDropOff*/, 7/*minScore*/, aliString);
    SEQAN_ASSERT_EQ(length(aliString), 2u);

    seq1 = "GCTGTAGCTGGAG";
    seq2 = "GCTGTAGCTGGAG";
    clear(aliString);
    testXDropAlign(seq1, seq2, Score<int>(1, -9, -9), 5/*scoreDropOff*/, 7/*minScore*/, aliString);
    SEQAN_ASSERT_EQ(length(aliString), 1u);
}

SEQAN_BEGIN_TESTSUITE(test_find_swift) {
    SEQAN_CALL_TEST(test_find_swift);
    SEQAN_CALL_TEST(test_longest_epsMatch);
    SEQAN_CALL_TEST(test_split_xDrop_align);

    //SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/index/find_swift.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/apps/swift_local/swift_local.h");
}
SEQAN_END_TESTSUITE

#endif
