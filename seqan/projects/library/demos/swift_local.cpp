#define SEQAN_DEBUG
#define SEQAN_TEST

///Some tests for the swift algorithm used as filter for local alignments.
#include <iostream>
#include <seqan/index.h>
#include <seqan/find/find_swift.h>
#include <seqan/seeds.h>
#include <demos/swift_local.h>

using namespace seqan;


// prints the positions and ranges of swift hits
void printOutput(String<String<char> > & positions,
                 String<String<char> > & ranges) {
	::std::cout << "hits at " << ::std::endl;
	for (unsigned i = 0; i < length(positions); i++) {
		::std::cout << value(positions, i) << ::std::endl;
		//::std::cout << value(ranges, i) << ::std::endl;
	}
	::std::cout << length(positions) << " hits found." << ::std::endl;
}

// compares the positions to the expected positions
// returns true only if both strings contain the same positions in the same order!!!
bool verifyOutput(String<String<char> > & positions,
                  String<String<char> > & expectedPositions) {
    if (length(positions) != length(expectedPositions)) return false;

    for (unsigned i = 0; i < length(positions); i++) {
        if (value(positions, i) != value(expectedPositions, i)) return false;
    }
    return true;
}

// calls swift and verifies the output
// if output differs from expected output, output is printed
template<typename TIndex>
void testLocalSwift(Finder<DnaString, Swift<SwiftLocal> > & finder,
					Pattern<TIndex, Swift<SwiftLocal> > & pattern,
					double epsilon,
					int minLength,
                    String<String<char> > & expectedPositions) {
    typedef typename Infix<DnaString>::Type TFinderInfix;
    typedef typename Infix<typename GetSequenceByNo<TIndex const >::Type >::Type TPatternInfix;

    typedef Pair<TFinderInfix, TPatternInfix> TRangePair;

    StringSet<String<Align<TPatternInfix> > > matches;
    resize(matches, countSequences(needle(pattern)));

    String<String<char> > positions;
    String<String<char> > ranges;

	while (find(finder, pattern, epsilon, minLength)) {
        // append output to positions and ranges string
        std::ostringstream pos, ran;
		pos << "< " << positionRange(finder) << " , " << positionRange(pattern) << " >";
		ran << range(finder) << "\n" << range(pattern);
        appendValue(positions, pos.str().c_str());
        appendValue(ranges, ran.str().c_str());

        // verification
        verifySwiftHitByLocalAlign(range(finder), range(pattern),
            epsilon, minLength, value(matches, pattern.curSeqNo));

        //unsigned thresh = _swiftBucketParams(pattern, pattern.curSeqNo).threshold;
        //unsigned w = _swiftBucketParams(pattern, pattern.curSeqNo).distanceCut;
        //unsigned q = length(indexShape(needle(pattern)));
        //verifySwiftHitByChains/*MaxChain*/(range(finder), range(pattern),
        //    epsilon, minLength, q, thresh, w, 1/*bandwidth*/,
        //    value(matches, pattern.curSeqNo));
	}

    std::cout << std::endl;
    for (unsigned i = 0; i < length(matches); i++) {
        std::cout << "Pattern sequence " << i << ":" <<  std::endl;
        for (unsigned j = 0; j < length(value(matches, i)); j++) {
            Align<TPatternInfix> m = value(value(matches, i), j);
            std::cout << "< " << toSourcePosition(row(m, 0), beginPosition(row(m, 0))) + beginPosition(source(row(m, 0)));
            std::cout << " , " << toSourcePosition(row(m, 0), endPosition(row(m, 0))) + beginPosition(source(row(m, 0))); 
            std::cout << " >< " << toSourcePosition(row(m, 1), beginPosition(row(m, 1))) + beginPosition(source(row(m, 1)));
            std::cout << " , " << toSourcePosition(row(m, 1), endPosition(row(m, 1))) + beginPosition(source(row(m, 1))) << " >" << std::endl;

            std::cout << value(value(matches, i), j) << std::endl;
        }
        std::cout << length(value(matches, i)) << " eps-matches" << std::endl;
    }

    if (verifyOutput(positions, expectedPositions))
        ::std::cout << "OK." << ::std::endl;
    else {
        ::std::cout << length(positions) << " swift hits" << ::std::endl;
        //printOutput(positions, ranges);
    }
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

void testOneLocalSwiftHit1() {
    typedef Index<DnaString, Index_QGram<UngappedShape<4> > > TQGramIndexSimple;
    typedef Finder<DnaString, Swift<SwiftLocal> > TFinder;

	// a single pattern and a single hit 
	::std::cout << "A single pattern and a single hit: ";

	DnaString text = "aaaaaaacgatcgatgcaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
	TFinder finder_swift(text);

	TQGramIndexSimple index_4gram("tttcgatcgatgctttttttttttttttttttttttttttttt");
    Pattern<TQGramIndexSimple, Swift<SwiftLocal> > pattern_swift(index_4gram);

    String<String<char> > expectedPositions;
    appendValue(expectedPositions, "< < 7 , 17 > , < 0 , 17 > >");

	testLocalSwift(finder_swift, pattern_swift, 0.1, 6, expectedPositions);
	::std::cout << ::std::endl;
}

void testOneLocalSwiftHit2() {
    typedef Index<DnaString, Index_QGram<UngappedShape<4> > > TQGramIndexSimple;
    typedef Finder<DnaString, Swift<SwiftLocal> > TFinder;

    // a single pattern and a single hit 2
	::std::cout << "A single pattern and a single hit 2: ";

	DnaString text = "aaaaagagaccccccagagaaaaaa";
	TFinder finder_swift(text);

	TQGramIndexSimple index_4gram("tttttggccccccggtttttt");
    Pattern<TQGramIndexSimple, Swift<SwiftLocal> > pattern_swift(index_4gram);

    String<String<char> > expectedPositions;
    appendValue(expectedPositions, "< < 9 , 15 > , < 0 , 15 > >");

	testLocalSwift(finder_swift, pattern_swift, 0.1, 6, expectedPositions);
	::std::cout << ::std::endl;
}

void testOneLocalSwiftHitStringSet() {
	typedef Index<StringSet<DnaString>, Index_QGram< UngappedShape<4> > > TQGramIndex;
    typedef Finder<DnaString, Swift<SwiftLocal> > TFinder;

	// a single hit, pattern in a StringSet
	::std::cout << "A single hit, pattern in a StringSet: ";

	DnaString text = "aaaaaaacgatcgatgcaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
	TFinder finder_swift(text);

	TQGramIndex index_4gram;
	appendValue(indexText(index_4gram), "tttcgatcgatgctttttttttttttttttttttttttttttt");
    Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift(index_4gram);

    String<String<char> > expectedPositions;
    appendValue(expectedPositions, "< < 7 , 17 > , < < 0 , 0 > , < 0 , 17 > > >");

	testLocalSwift(finder_swift, pattern_swift, 0.1, 6, expectedPositions);
	::std::cout << ::std::endl;
}

void testOneLocalSwiftHitBucketBorder() {
	typedef Index<StringSet<DnaString>, Index_QGram< UngappedShape<4> > > TQGramIndex;
    typedef Finder<DnaString, Swift<SwiftLocal> > TFinder;

	// hit at bucket border
	::std::cout << "Hit at bucket border: ";

	DnaString text = "aaacgatcgatgcaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
	TFinder finder_swift(text);

	TQGramIndex index_4gram;
	appendValue(indexText(index_4gram), "tttcgatcgatgctttttttttttttttttttttttttttttt");
    Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift(index_4gram);

    String<String<char> > expectedPositions;
    appendValue(expectedPositions, "< < 3 , 13 > , < < 0 , 0 > , < 0 , 13 > > >");
    appendValue(expectedPositions, "< < 3 , 13 > , < < 0 , 2 > , < 0 , 29 > > >");

	testLocalSwift(finder_swift, pattern_swift, 0.1, 6, expectedPositions);
	::std::cout << ::std::endl;
}

void testOneLocalSwiftHitNegDiag() {
	typedef Index<StringSet<DnaString>, Index_QGram< UngappedShape<4> > > TQGramIndex;
    typedef Finder<DnaString, Swift<SwiftLocal> > TFinder;

	// hit at lower haystack position -> diagonal start position negative
	::std::cout << "Hit at lower haystack position -> diagonal start position negative: ";

	DnaString text = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaacgatcgatgc";
	TFinder finder_swift(text);

	TQGramIndex index_4gram;
	appendValue(indexText(index_4gram), "ttttcgatcgatgctttttttttttttttttttttttttttttt");
    Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift(index_4gram);

    String<String<char> > expectedPositions;
    appendValue(expectedPositions, "< < 36 , 46 > , < < 0 , 3 > , < 0 , 30 > > >");
    appendValue(expectedPositions, "< < 36 , 46 > , < < 0 , 0 > , < 0 , 14 > > >");

	testLocalSwift(finder_swift, pattern_swift, 0.1, 10, expectedPositions);
	::std::cout << ::std::endl;
}

void testLocalSwiftTwoPatterns() {
	typedef Index<StringSet<DnaString>, Index_QGram< UngappedShape<4> > > TQGramIndex;
    typedef Finder<DnaString, Swift<SwiftLocal> > TFinder;

	// two pattern sequences
	::std::cout << "Two pattern sequences: ";

	DnaString text = "aaaacgttccaaaaaaa";
	TFinder finder_swift(text);

	TQGramIndex index_4gram;
	appendValue(indexText(index_4gram), "ttttcgttcctttttttttttttttttttttt"); // length = 32
    appendValue(indexText(index_4gram), "ttcgttcctt"); // length = 10
	Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift(index_4gram);

    String<String<char> > expectedPositions;
    appendValue(expectedPositions, "< < 4 , 10 > , < < 0 , 0 > , < 0 , 10 > > >");
    appendValue(expectedPositions, "< < 4 , 10 > , < < 0 , 3 > , < 0 , 26 > > >");
    appendValue(expectedPositions, "< < 4 , 10 > , < < 1 , 0 > , < 1 , 10 > > >");

	testLocalSwift(finder_swift, pattern_swift, 0.1, 6, expectedPositions);
	::std::cout << ::std::endl;
}

void testLocalSwiftLongPatterns() {
	typedef Index<StringSet<DnaString>, Index_QGram< UngappedShape<4> > > TQGramIndex;
    typedef Finder<DnaString, Swift<SwiftLocal> > TFinder;

	// two longer pattern sequences
	::std::cout << "Two longer pattern sequences: ";

	DnaString text = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaacgatcagtgacaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
	TFinder finder_swift(text);

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
	::std::cout << ::std::endl;
}

void testLocalSwiftArgSeqs(const char *argv[]) {
	typedef Index<StringSet<DnaString>, Index_QGram< UngappedShape<5> > > TQGramIndex;
    typedef Finder<DnaString, Swift<SwiftLocal> > TFinder;

	// adenovirus sequences
	::std::cout << "Adenovirus sequences: ";

	DnaString text = String<Dna, FileReader<Fasta> >(argv[1]);
	TFinder finder_swift(text);

	TQGramIndex index_qgram;
	appendValue(indexText(index_qgram), String<Dna, FileReader<Fasta> >(argv[2]));
	Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift(index_qgram);

    String<String<char> > expectedPositions; 

	testLocalSwift(finder_swift, pattern_swift, 0.1, 30, expectedPositions);
}

int main(int argc, const char *argv[]) {
    // Get rid of warnings for unused variables.
    (void)argc;
    (void)argv;

    testOneLocalSwiftHit1();
    testOneLocalSwiftHit2();
    testOneLocalSwiftHitStringSet();
    testOneLocalSwiftHitBucketBorder();
    testOneLocalSwiftHitNegDiag();
    testLocalSwiftTwoPatterns();
    testLocalSwiftLongPatterns();
    
	time_t startTime = time(0);
	unsigned iterations = 1;
	for(unsigned t = 0; t < iterations; t++) {
        testLocalSwiftArgSeqs(argv);
	}
	double runningTime = (time(0)-startTime)/(double)iterations;

    std::cout << "Running time: " << runningTime << "s" << std::endl;

    testShrink1();
    testShrink2();
    testShrink3();

	return 0;
}