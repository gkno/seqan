#include <iostream>

#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/basic/basic_testing.h>
#include <seqan/index.h>
#include <seqan/seeds.h>

#include <demos/swift_local.h>

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
    appendValue(expectedPositions, "< < 36 , 46 > , < < 0 , 3 > , < 0 , 30 > > >");
    appendValue(expectedPositions, "< < 36 , 46 > , < < 0 , 0 > , < 0 , 14 > > >");

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
    appendValue(expectedPositions, "< < 90 , 96 > , < < 0 , 41 > , < 0 , 63 > > >");
    appendValue(expectedPositions, "< < 88 , 99 > , < < 0 , 7 > , < 0 , 35 > > >");
    appendValue(expectedPositions, "< < 90 , 96 > , < < 1 , 41 > , < 1 , 63 > > >");
    appendValue(expectedPositions, "< < 88 , 99 > , < < 1 , 7 > , < 1 , 35 > > >");

	testLocalSwift(finder_swift, pattern_swift, 0.1, 6, expectedPositions);
}

void testShrink1() {
    DnaString seq1 = "ACCTTTGCCCCCCCCCCTAAAAAAAATTAAAA";
    DnaString seq2 = "ACGTTTACCCCCCCCCCGAAAAAAAAGAAAA";

    Align<DnaString> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), seq1);
    assignSource(row(alignment, 1), seq2);

    Score<int> score(1, -1, -1);
    globalAlignment(alignment, score);

    shrinkToMaxEpsMatch(alignment, 6, 0.125);

    SEQAN_TASSERT(row(alignment, 0) == "ACCTTTGCCCCCCCCCCTAAAAAAAA");
    SEQAN_TASSERT(row(alignment, 1) == "ACGTTTACCCCCCCCCCGAAAAAAAA");
}

void testShrink2() {
    DnaString seq1 = "ACCTTTGCCCCCCCCCCTAAAAAAAATTAAAA";
    DnaString seq2 = "ACGTTTCCCCCCCCCCGAAAAAAAAGAAAA";

    Align<DnaString> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), seq1);
    assignSource(row(alignment, 1), seq2);

    Score<int> score(1, -1, -1);
    globalAlignment(alignment, score);
    
    shrinkToMaxEpsMatch(alignment, 6, 0.125);
    
    SEQAN_TASSERT(row(alignment, 0) == "ACCTTTGCCCCCCCCCCTAAAAAAAA");
    SEQAN_TASSERT(row(alignment, 1) == "ACGTTT-CCCCCCCCCCGAAAAAAAA");
}

void testShrink3() {
    DnaString seq1 = "AAAATTAAAAAAAATCCCCCCCCCCGTTTCCA";
    DnaString seq2 = "AAAAGAAAAAAAAGCCCCCCCCCCTTTGCA";

    Align<DnaString> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), seq1);
    assignSource(row(alignment, 1), seq2);

    Score<int> score(1, -1, -1);
    globalAlignment(alignment, score);

    shrinkToMaxEpsMatch(alignment, 6, 0.125);

    SEQAN_TASSERT(row(alignment, 0) == "AAAAAAAATCCCCCCCCCCGTTTCCA");
    SEQAN_TASSERT(row(alignment, 1) == "AAAAAAAAGCCCCCCCCCC-TTTGCA");
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

SEQAN_DEFINE_TEST(test_find_longest_epsMatch) {
    testShrink1();
    testShrink2();
    testShrink3();
}

SEQAN_BEGIN_TESTSUITE(test_find_swift) {
    SEQAN_CALL_TEST(test_find_swift);
    SEQAN_CALL_TEST(test_find_longest_epsMatch);
}
SEQAN_END_TESTSUITE
