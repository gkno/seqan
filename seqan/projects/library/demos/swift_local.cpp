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
		::std::cout << value(ranges, i) << ::std::endl;
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
        epsilon = pattern._currentErrorRate;
        unsigned thresh = _swiftBucketParams(pattern, pattern.curSeqNo).threshold;
        unsigned w = _swiftBucketParams(pattern, pattern.curSeqNo).distanceCut;
        unsigned q = length(indexShape(needle(pattern)));
        verifySwiftHit(range(finder), range(pattern),
            epsilon, minLength, q, thresh, w, 1/*bandwidth*/,
            value(matches, pattern.curSeqNo));

        //verifySwiftHitByLocalAlign(range(finder), range(pattern),
        //    epsilon, minLength, value(matches, pattern.curSeqNo));
	}

    std::cout << std::endl;
    for (unsigned i = 0; i < length(matches); i++) {
        std::cout << "Pattern sequence " << i << ":" <<  std::endl;
        for (unsigned j = 0; j < length(value(matches, i)); j++) {
            std::cout << value(value(matches, i), j) << std::endl;
        }
    }

    if (verifyOutput(positions, expectedPositions))
        ::std::cout << "OK." << ::std::endl;
    else {
        ::std::cout << length(positions) << " hits found." << ::std::endl;
        printOutput(positions, ranges);
    }
}

int main(int argc, const char *argv[]) {
    // Get rid of warnings for unused variables.
    (void)argc;
    (void)argv;
    
    typedef Index<DnaString, Index_QGram<UngappedShape<4> > > TQGramIndexSimple;
	typedef Index<StringSet<DnaString>, Index_QGram< UngappedShape<4> > > TQGramIndex;
    typedef Finder<DnaString, Swift<SwiftLocal> > TFinder;
    
    typedef Pair<Position<TFinder>::Type> TFinderPosition;
    typedef Pair<SAValue<TQGramIndex>::Type> TPatternPosition;
    typedef Pair<TFinderPosition, TPatternPosition> TPosPair;

	// a single pattern and a single hit 
	::std::cout << "A single pattern and a single hit: ";

	DnaString text_0 = "aaaaaaacgatcgatgcaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
	TFinder finder_swift_0(text_0);

	TQGramIndexSimple index_4gram_0("tttcgatcgatgctttttttttttttttttttttttttttttt");
    Pattern<TQGramIndexSimple, Swift<SwiftLocal> > pattern_swift_0(index_4gram_0);

    String<String<char> > expectedPositions_0;
    appendValue(expectedPositions_0, "< < 7 , 17 > , < 0 , 17 > >");

	testLocalSwift(finder_swift_0, pattern_swift_0, 0.1, 6, expectedPositions_0);
	::std::cout << ::std::endl;

 //   // a single pattern and a single hit 2
	//::std::cout << "A single pattern and a single hit 2: ";

	//DnaString text_0a = "aaaaagagaccccccagagaaaaaa";
	//TFinder finder_swift_0a(text_0a);

	//TQGramIndexSimple index_4gram_0a("tttttggccccccggtttttt");
 //   Pattern<TQGramIndexSimple, Swift<SwiftLocal> > pattern_swift_0a(index_4gram_0a);

 //   String<String<char> > expectedPositions_0a;
 //   appendValue(expectedPositions_0a, "< < 9 , 15 > , < 0 , 15 > >");

	//testLocalSwift(finder_swift_0a, pattern_swift_0a, 0.1, 6, expectedPositions_0a);
	//::std::cout << ::std::endl;

	//// a single hit, pattern in a StringSet
	//::std::cout << "A single hit, pattern in a StringSet: ";

	//DnaString text_1 = "aaaaaaacgatcgatgcaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
	//TFinder finder_swift_1(text_1);

	//TQGramIndex index_4gram_1;
	//appendValue(indexText(index_4gram_1), "tttcgatcgatgctttttttttttttttttttttttttttttt");
 //   Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift_1(index_4gram_1);

 //   String<String<char> > expectedPositions_1;
 //   appendValue(expectedPositions_1, "< < 7 , 17 > , < < 0 , 0 > , < 0 , 17 > > >");

	//testLocalSwift(finder_swift_1, pattern_swift_1, 0.1, 6, expectedPositions_1);
	//::std::cout << ::std::endl;


	//// hit at bucket border
	//::std::cout << "Hit at bucket border: ";

	//DnaString text_2 = "aaacgatcgatgcaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
	//TFinder finder_swift_2(text_2);

	//TQGramIndex index_4gram_2;
	//appendValue(indexText(index_4gram_2), "tttcgatcgatgctttttttttttttttttttttttttttttt");
 //   Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift_2(index_4gram_2);

 //   String<String<char> > expectedPositions_2;
 //   appendValue(expectedPositions_2, "< < 3 , 13 > , < < 0 , 0 > , < 0 , 13 > > >");
 //   appendValue(expectedPositions_2, "< < 3 , 13 > , < < 0 , 2 > , < 0 , 29 > > >");

	//testLocalSwift(finder_swift_2, pattern_swift_2, 0.1, 6, expectedPositions_2);
	//::std::cout << ::std::endl;


	//// hit at lower haystack position -> diagonal start position negative
	//::std::cout << "Hit at lower haystack position -> diagonal start position negative: ";

	//DnaString text_3 = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaacgatcgatgc";
	//TFinder finder_swift_3(text_3);

	//TQGramIndex index_4gram_3;
	//appendValue(indexText(index_4gram_3), "ttttcgatcgatgctttttttttttttttttttttttttttttt");
 //   Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift_3(index_4gram_3);

 //   String<String<char> > expectedPositions_3;
 //   appendValue(expectedPositions_3, "< < 36 , 46 > , < < 0 , 3 > , < 0 , 30 > > >");
 //   appendValue(expectedPositions_3, "< < 36 , 46 > , < < 0 , 0 > , < 0 , 14 > > >");

	//testLocalSwift(finder_swift_3, pattern_swift_3, 0.2, 10, expectedPositions_3);
	//::std::cout << ::std::endl;


	//// two pattern sequences
	//::std::cout << "Two pattern sequences: ";

	//DnaString text_4 = "aaaacgttccaaaaaaa";
	//TFinder finder_swift_4(text_4);

	//TQGramIndex index_4gram_4;
	//appendValue(indexText(index_4gram_4), "ttttcgttcctttttttttttttttttttttt"); // length = 32
 //   appendValue(indexText(index_4gram_4), "ttcgttcctt"); // length = 10
	//Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift_4(index_4gram_4);

 //   String<String<char> > expectedPositions_4;
 //   appendValue(expectedPositions_4, "< < 4 , 10 > , < < 0 , 0 > , < 0 , 10 > > >");
 //   appendValue(expectedPositions_4, "< < 4 , 10 > , < < 0 , 3 > , < 0 , 26 > > >");
 //   appendValue(expectedPositions_4, "< < 4 , 10 > , < < 1 , 0 > , < 1 , 10 > > >");

	//testLocalSwift(finder_swift_4, pattern_swift_4, 0.1, 6, expectedPositions_4);
	//::std::cout << ::std::endl;


	//// two longer pattern sequences
	//::std::cout << "Two longer pattern sequences: ";

	//DnaString text_5 = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaacgatcagtgacaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
	//TFinder finder_swift_5(text_5);

	//TQGramIndex index_4gram_5;
	//appendValue(indexText(index_4gram_5), "ttttttttttttttcgatcagtgacttttttttttttttttttttttttttttttttatcagt"); // length = 63
 //   appendValue(indexText(index_4gram_5), "ttttttttttttttcgatcagtgacttttttttttttttttttttttttttttttttatcagt"); // length = 63
	//Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift_5(index_4gram_5);

 //   String<String<char> > expectedPositions_5;
 //   appendValue(expectedPositions_5, "< < 90 , 96 > , < < 0 , 41 > , < 0 , 63 > > >");
 //   appendValue(expectedPositions_5, "< < 88 , 99 > , < < 0 , 7 > , < 0 , 35 > > >");
 //   appendValue(expectedPositions_5, "< < 90 , 96 > , < < 1 , 41 > , < 1 , 63 > > >");
 //   appendValue(expectedPositions_5, "< < 88 , 99 > , < < 1 , 7 > , < 1 , 35 > > >");

	//testLocalSwift(finder_swift_5, pattern_swift_5, 0.1, 6, expectedPositions_5);
	//::std::cout << ::std::endl;

	//// adenovirus sequences
	//::std::cout << "Adenovirus sequences: ";

	//DnaString text_6 = String<Dna, FileReader<Fasta> >(argv[1]);
	//TFinder finder_swift_6(text_6);

	//TQGramIndex index_4gram_6;
	//appendValue(indexText(index_4gram_6), String<Dna, FileReader<Fasta> >(argv[2]));
	//Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift_6(index_4gram_6);

 //   String<String<char> > expectedPositions_6; 

	//testLocalSwift(finder_swift_6, pattern_swift_6, 0.1, 6, expectedPositions_6);


	return 0;
}
