/// This code example illustrates a graph-based Gotoh alignment in linear space (Hirschberg)
#include <seqan/graph_align.h>
#include <iostream>

using namespace seqan;

int main() {
	typedef String<char> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
/// Alignments are carried out on a StringSet that holds the sequences
	TStringSet str;
	TString str0("TarfieldandGarfieldarestupid.");appendValue(str, str0);
	TString str1("Garfield");appendValue(str, str1);
/// Configuration of alignment algorithm: Scoring (Match = 2, Mismatch = -1, Gap-extension = -1, Gap-opening = -4)	
	Score<int> score_type = Score<int>(2,-1,-1,-4);
/// Out-parameter: An alignment graph
	TGraph g(str);
/// Global alignment with Gotoh in linear space (Hirschberg)
	int score = globalAlignment(g, score_type, Hirschberg());
/// Console output
	std::cout << "Scoring schema: Match=2, Mismatch=-1, Gap-extension=-1, Gap-opening=-4" << std::endl;
	std::cout << g << std::endl;
	std::cout << "Score: " << score << std::endl;
	return 0;
}
