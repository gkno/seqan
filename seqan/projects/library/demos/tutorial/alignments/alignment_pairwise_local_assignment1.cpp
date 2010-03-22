//FRAGMENT(main)
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
//FRAGMENT(init)
	Align< String<AminoAcid> > ali;
	appendValue(rows(ali), "PNCFDAKQRTASRPL");
	appendValue(rows(ali), "CFDKQKNNRTATRDTA");


//FRAGMENT(ali)
	LocalAlignmentFinder<> finder(ali);
	Score<int> sc(3,-2,-5,-1);
	unsigned count = 0;
	while (localAlignment(ali, finder, sc, 0, WatermanEggert()) && count < 3) {
		::std::cout << "Score = " << getScore(finder) << ::std::endl;
		::std::cout << ali;
		::std::cout << "Aligns Seq1[" << sourceBeginPosition(row(ali, 0)) << ":" << (sourceEndPosition(row(ali, 0))-1) << "]";
		::std::cout << " and Seq2[" << sourceBeginPosition(row(ali, 1)) << ":" <<  (sourceEndPosition(row(ali, 1))-1) << "]" << ::std::endl << ::std::endl;
		++count;
	}
	return 0;
}
