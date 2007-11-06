/// This code example illustrates a graph-based Smith-Waterman alignment
#include <seqan/graph.h>
#include <iostream>

using namespace seqan;

int main() {
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, Dna> > TGraph;
/// Alignments are carried out on a StringSet that holds the sequences
	TStringSet str;
	TString str0("TTGACACCCTCCCAATTGTA"); appendValue(str, str0);
	TString str1("ACCCCAGGCTTTACACAT"); appendValue(str, str1);
/// Configuration of alignment algorithm: Scoring (Match = 2, Mismatch = -1, Gap-extension = -1, Gap-opening = -2)
	Score<int> score_type = Score<int>(2,-1,-1,-2);
/// We could again use a graph to retrieve the local alignment, but as with all other alignment algorithms we can also use a string of fragments or simply std::cout.
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	TFragmentString matches;
/// Local alignment with Smith-Waterman
	int score = localAlignment(matches, str, score_type, SmithWaterman() );
/// Console output
	std::cout << "Scoring schema: Match=2, Mismatch=-1, Gap-extension=-1, Gap-opening=-2" << std::endl;
	std::cout << str0 << std::endl;
	std::cout << str1 << std::endl;
	std::cout << "Local match: " << std::endl;
	std::cout << label(matches[0], str, 0) << std::endl;
	std::cout << label(matches[0], str, 1) << std::endl;
	std::cout << "Score: " << score << std::endl;
	return 0;
}

