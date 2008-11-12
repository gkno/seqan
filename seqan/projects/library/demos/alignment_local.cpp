#include <iostream>
#include <seqan/align.h>

using namespace std;
using namespace seqan;

int main()
{
///Example 1: This program applies the Smith-Waterman algorithm to compute the best local alignment between two given sequences:
	Align< String<char> > ali;
	appendValue(rows(ali), "aphilologicaltheorem");
	appendValue(rows(ali), "bizarreamphibology");
    cout << "Score = " << localAlignment(ali, Score<int>(3,-3,-2), SmithWaterman()) << endl;
	cout << ali;
	cout << "Aligns Seq1[" << sourceBeginPosition(row(ali, 0)) << ":" << (sourceEndPosition(row(ali, 0))-1) << "]";
	cout << " and Seq2[" << sourceBeginPosition(row(ali, 1)) << ":" <<  (sourceEndPosition(row(ali, 1))-1) << "]" << endl << endl;


///Example 2: This program applies the Waterman-Eggert algorithm to compute all non-overlapping local alignments with score better than 2:
	Align< String<Dna> > ali2;
	appendValue(rows(ali2), "ataagcgtctcg");
	appendValue(rows(ali2), "tcatagagttgc");

	LocalAlignmentFinder<> finder(ali2);
	Score<int> scoring(2, -1, -2, 0);
	while (localAlignment(ali2, finder, scoring, 2, WatermanEggert()))
	{
		cout << "Score = " << getScore(finder) << endl;
		cout << ali2;
		cout << "Aligns Seq1[" << sourceBeginPosition(row(ali2, 0)) << ":" << (sourceEndPosition(row(ali2, 0))-1) << "]";
		cout << " and Seq2[" << sourceBeginPosition(row(ali2, 1)) << ":" <<  (sourceEndPosition(row(ali2, 1))-1) << "]" << endl << endl;
	}

	return 0;
}
