#include <iostream>
#include <seqan/align.h>
#include <seqan/graph.h>

using namespace std;
using namespace seqan;

int main()
{
///Here are some sequences to align:
	DnaString seq1 = "atcgaatgcgga";
	DnaString seq2 = "actcgttgca";
///Now we choose a scoring scheme with affine gap costs ("gap open" == -2, "gap extend" == -1).
	Score<int> score(0, -1, -2, -1);
///Example 1: We use @Class.Align@ to align the two sequences. 
///Since we does not specify an @Tag.Global Alignment Algorithms|algorithm tag@ when we call @Function.globalAlignment@, 
///a suitable algorithm (@Tag.Global Alignment Algorithms|Gotoh@) is automatically choosen.
	Align<DnaString, ArrayGaps> align;
	resize(rows(align), 2);
	assignSource(row(align, 0), seq1);
	assignSource(row(align, 1), seq2);

	cout << "Score = " << globalAlignment(align, score) << endl;
	cout << align << endl;
///Example 2: We now choose explicitely the algorithm @Tag.Global Alignment Algorithms|MyersHirschberg@.
///Since this algorithm always works on Levenshtein distance, $score$ is ignored here.
///Therefore, this algorithm computes a different alignment and returns a different score.
	cout << "Score = " << globalAlignment(align, score, MyersHirschberg()) << endl;
	cout << align << endl;
///Example 3: We now do the same as in case 1, but now we use an @Spec.Alignment Graph@ for storing the alignment.
///Here we use @Tag.Global Alignment Algorithms|Hirschberg's algorithm@.
	typedef StringSet<DnaString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;

	TStringSet string_set;
	appendValue(string_set, seq1);
	appendValue(string_set, seq2);
	TAlignmentGraph alignment_graph(string_set);

	cout << "Score = " << globalAlignment(alignment_graph, score, Hirschberg()) << endl;
	cout << alignment_graph << endl;

	return 0;
}
