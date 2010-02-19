// FRAGMENT(includes)
#include <iostream>
#include "seqan/find_motif.h"

using namespace seqan;

// FRAGMENT(typedefs)
int main() 
{
    typedef MotifFinder<Dna, EPatternBranching> TMotifFinder;
    typedef String<DnaString> TString;
    typedef Size<TString>::Type TSize;

// FRAGMENT(sequences)
	TString dataset;
	appendValue(dataset, DnaString("ACAGCA"));
	appendValue(dataset, DnaString("AGGCAG"));
	appendValue(dataset, DnaString("TCAGTC"));
    TSize seqCount = length(dataset);

// FRAGMENT(initialization)
	::std::srand((unsigned) time(NULL));

	unsigned int motifLength = 2;		//length of motif
	unsigned int mm = 1;	        	//number of mismatches
	bool is_exact = false;	            //occurences of motif need to have exactly mm mismatches
	unsigned int h = 0;                 //size of the neighborhood considering at first 

	TMotifFinder finder_epb1(seqCount, motifLength, mm, is_exact, h);

// FRAGMENT(search)
	findMotif(finder_epb1,dataset,OMOPS());

	for (int i = 0; i < (int) motifCount(finder_epb1); ++i)
		::std::cout << i << ": " << getMotif(finder_epb1, i) << ::std::endl;

	return 0;
}

