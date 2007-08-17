#include <fstream>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/align.h>
#include <seqan/graph.h>

using namespace std;
using namespace seqan;

const unsigned loops = 2;

int main(int argc, char* argv[])
{
	Align<DnaString, ArrayGaps> align;
	resize(rows(align), 2);

	// levenshtein-score (match, mismatch, gapOpen, gapExtend)
	Score<int> score = Score<int>(0, -1, -2, -1);

	DnaString seq1 = "atcgaatgcgga";
	DnaString seq2 = "actcgttgca";
/*
	ifstream file;
	file.open(argv[1], ios_base::in | ios_base::binary);
	read(file, seq1, Fasta());
	file.close();

	file.open(argv[2], ios_base::in | ios_base::binary);
	read(file, seq2, Fasta());
	file.close();
*/	

	assignSource(row(align, 0), seq1);
	assignSource(row(align, 1), seq2);
/*
	needlemanWunsch(align, score);

	cout << align << endl;
*/

	typedef StringSet<DnaString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;

	// Ordinary DP alignment
	TStringSet str;
	appendValue(str, seq1);
	appendValue(str, seq2);
	TGraph g(str);

	clock_t start, finish;
	start = clock();
	for(unsigned loop = 0; loop < loops; ++loop)
		globalAlignment(g, score, Hirschberg());
	finish = clock();

	double duration = (double)(finish - start) / (double)(loops * CLOCKS_PER_SEC);
	cout << "aligning took " << duration << " seconds" << endl;

//	cout << align << endl;

	return 0;
}
