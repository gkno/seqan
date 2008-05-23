/// This code example illustrates a graph-based Needleman-Wunsch alignment
#include <seqan/graph_align.h>
#include <iostream>

using namespace seqan;

int main() {
	typedef String<char> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
/// Alignments are carried out on a StringSet that holds the sequences
	TStringSet str;
	TString str0("Myannealing");appendValue(str, str0);
	TString str1("annual"); appendValue(str, str1);
/// Configuration of alignment algorithm: Scoring (Match = 0, Mismatch = -1, Gap = -1)
	Score<int> score_type = Score<int>(0,-1,-1,0);
/// Out-parameter: An alignment graph	
	TGraph g(str);
/// Global alignment with Needleman-Wunsch
	int score = globalAlignment(g, score_type, NeedlemanWunsch() );
/// Console output
	std::cout << "Scoring schema: Match=0, Mismatch=-1, Gap=-1" << std::endl;
	std::cout << g << std::endl;
	std::cout << "Score: " << score << std::endl;
/// The initialization of the matrix and the traceback start can also be controlled.
/// Configuration:	AlignConfig<TTop, TLeft, TRight, TBottom>.
/// If TTop, TLeft, TRight, or TBottom is true the corresponding row / column is initialized with 0's or searched for the maximum, respectively.
	AlignConfig<true,false,false,true> ac;
/// Global alignment with Needleman-Wunsch and special initialization
	int score2 = globalAlignment(stringSet(g), score_type, ac, NeedlemanWunsch() );
/// Console output
	std::cout << "Score with ends free-space alignment: " << score2 << std::endl;
	return 0;
}
