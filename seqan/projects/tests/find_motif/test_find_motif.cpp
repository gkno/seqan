#include <iostream>

#define SEQAN_DEBUG
#define SEQAN_TEST

#define TEST_PATH "/projects/tests/find_motif/"

#include <seqan/find_motif.h>
#include <seqan/align.h>

using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

void Test_approximationAlgorithms()
{
	unsigned int t = 0;		//number of input sequences
	unsigned int n = 0;		//length of sequence
	unsigned int l = 0;		//length of motif
	unsigned int d = 0;		//number of substitutions
	bool is_exact = false;	//size of Hamming distance
	unsigned int m =0;		//total number of possible l-mers
	unsigned int h = 0;		//size of the neighborhood considering at first
	unsigned int i = 0;

	srand((unsigned) time(NULL));

//____________________________________________________________________________
// Test1 - Search for OOPS motifs on a small set of nucleotide sequences
//         given the exact Hamming distance (=d)
	t = 3;		
	n = 6;		
	l = 4;		
	d = 1;		
	is_exact = true;	
	m = t*(n-l+1);

	String<DnaString> dataset1;
	appendValue(dataset1,DnaString("ACAGCA"));
	appendValue(dataset1,DnaString("AGGCAG"));
	appendValue(dataset1,DnaString("TCAGTC"));

	//Application of PROJECTION-OOPS
    MotifFinder<Dna, Projection> motif_finder1(t,l,m,d,is_exact);
	findMotif(motif_finder1, dataset1, OOPS());
	SEQAN_TASSERT(motif_finder1.consensus_pattern=="AGCC");

	//Application of ePatternBranching-OOPS
    MotifFinder<Dna, EPatternBranching> motif_finder2(t,l,d,is_exact,h);
	findMotif(motif_finder2, dataset1, OOPS());
	SEQAN_TASSERT(motif_finder2.set_of_motifs[0]=="AGCC");

//____________________________________________________________________________
// Test2 - Search for OMOPS motifs on a small set of nucleotide sequences
//         given the exact Hamming distance (<=d)

	is_exact = false;

	//Application of PROJECTION-OMOPS
    MotifFinder<Dna, Projection> motif_finder3(t,l,m,d,is_exact);
	findMotif(motif_finder3, dataset1, OMOPS());
	displayResult(motif_finder3);
	//SEQAN_TASSERT(motif_finder3.consensus_pattern=="AGCC");

	//Application of ePatternBranching-OMOPS
	MotifFinder<Dna, EPatternBranching> motif_finder4(t,l,d,is_exact,h);
	findMotif(motif_finder4, dataset1, OMOPS());
	displayResult(motif_finder4);
	//SEQAN_TASSERT(motif_finder4.set_of_motifs[0]=="AGCC");

//____________________________________________________________________________
// Test3 - Search for ZOOPS motifs on a set of small nucleotide sequences
//         given the inexact Hamming distance (<=d)

	//Application of PROJECTION-ZOOPS
    MotifFinder<Dna, Projection> motif_finder5(t,l,m,d,is_exact);
	findMotif(motif_finder5, dataset1, ZOOPS());
	SEQAN_TASSERT(motif_finder5.consensus_pattern=="TCAG");

//____________________________________________________________________________
// Test4 - Search for TCM motifs on a set of small nucleotide sequences
//         given the exact Hamming distance (=d)

	is_exact = true;

	//Application of PROJECTION-TCM
    MotifFinder<Dna, Projection> motif_finder6(t,l,m,d,is_exact);
	findMotif(motif_finder6, dataset1, ZOOPS());
	SEQAN_TASSERT(motif_finder6.consensus_pattern=="TCAG");
}

void Test_exactAlgorithms()
{
	unsigned int l = 0;		//length of motif
	unsigned int d = 0;		//number of substitutions
	bool is_exact = false;	//size of Hamming distance
	unsigned int i = 0;

//____________________________________________________________________________
// Test1 - Search for OOPS motifs on a small set of nucleotide sequences
//         given the exact Hamming distance (=d)

	l = 8;		
	d = 2;		
	is_exact = true;	

	String<DnaString> dataset1;
	appendValue(dataset1,DnaString("CAGAAGGCTCTAAACAGGTA"));
	appendValue(dataset1,DnaString("CCACAAATCTTTCTCCGGCG"));
	appendValue(dataset1,DnaString("TCGACTGAAATGGAGAACAG"));
	appendValue(dataset1,DnaString("GCTTACGTACTGAACGTGCG"));
	appendValue(dataset1,DnaString("TAACTTACACTCAACAGGTG"));

	//Application of PMS1-OOPS
	MotifFinder<Dna, PMS1> motif_finder1(l,d,is_exact);
	findMotif(motif_finder1,dataset1,OOPS());

	//Application of PMSP-OOPS
	MotifFinder<Dna, PMSP> motif_finder2(l,d,is_exact);
	findMotif(motif_finder2,dataset1,OOPS());

	SEQAN_TASSERT(length(motif_finder1.set_of_motifs)==length(motif_finder2.set_of_motifs));

	for(i=0; i<length(motif_finder1.set_of_motifs); ++i)
	{
		SEQAN_TASSERT(motif_finder1.set_of_motifs[i]==motif_finder2.set_of_motifs[i]);
	}

//____________________________________________________________________________
// Test2 - Search for OMOPS motifs on a small set of nucleotide sequences
//         given the inexact Hamming distance (<=d)

	l = 6;		
	d = 2;		
	is_exact = false;	

	String<DnaString> dataset2;
	appendValue(dataset2,DnaString("TTTAGGCCACTCGATAAGAG"));
	appendValue(dataset2,DnaString("GAGCTACAACCGAGGAATTG"));
	appendValue(dataset2,DnaString("CGAATGTTGTGACTCGATAA"));
	appendValue(dataset2,DnaString("GTCCGAAAATACCTGGAGCC"));
	appendValue(dataset2,DnaString("CCGCTTAAAAATGTATTATC"));

	//Application of PMS1-OMOPS
	MotifFinder<Dna, PMS1> motif_finder3(l,d,is_exact);
	findMotif(motif_finder3,dataset2,OMOPS());

	//Application of PMSP-OMOPS
	MotifFinder<Dna, PMSP> motif_finder4(l,d,is_exact);
	findMotif(motif_finder4,dataset2,OMOPS());

	SEQAN_TASSERT(length(motif_finder3.set_of_motifs)==length(motif_finder4.set_of_motifs));

	for(i=0; i<length(motif_finder3.set_of_motifs); ++i)
	{
		SEQAN_TASSERT(motif_finder3.set_of_motifs[i]==motif_finder4.set_of_motifs[i]);
	}

//____________________________________________________________________________
// Test3 - Search for ZOOPS motifs on a small set of nucleotide sequences
//         given the exact Hamming distance (=d)

	l = 6;		
	d = 1;		
	is_exact = true;	

	String<DnaString> dataset3;
	appendValue(dataset3,DnaString("CTATTACCCAGAGTGCCAAG"));
	appendValue(dataset3,DnaString("AACGGACCCAACGGCAAAGA"));
	appendValue(dataset3,DnaString("TCTTGTCTTCACAGTCAGAT"));

	//Application of PMS1-ZOOPS
	MotifFinder<Dna, PMS1> motif_finder5(l,d,is_exact);
	findMotif(motif_finder5,dataset3,ZOOPS());

	//Application of PMSP-ZOOPS
	MotifFinder<Dna, PMSP> motif_finder6(l,d,is_exact);
	findMotif(motif_finder6,dataset3,ZOOPS());

	SEQAN_TASSERT(length(motif_finder5.set_of_motifs)==length(motif_finder6.set_of_motifs));

	for(i=0; i<length(motif_finder5.set_of_motifs); ++i)
	{
		SEQAN_TASSERT(motif_finder5.set_of_motifs[i]==motif_finder6.set_of_motifs[i]);
	}
}

int main() 
{
	SEQAN_TREPORT("TEST FIND MOTIF BEGIN")

	Test_approximationAlgorithms();
	Test_exactAlgorithms();

	debug::verifyCheckpoints("projects/library/seqan/find_motif/find_motif_base.h");
	debug::verifyCheckpoints("projects/library/seqan/find_motif/find_motif_projection.h");
	debug::verifyCheckpoints("projects/library/seqan/find_motif/find_motif_epatternbranching.h");
	debug::verifyCheckpoints("projects/library/seqan/find_motif/find_motif_pms1.h");
	debug::verifyCheckpoints("projects/library/seqan/find_motif/find_motif_pmsp.h");

	SEQAN_TREPORT("TEST FIND MOTIF END");

	return 0;
}


/*
#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/find_motif.h>

using namespace seqan;

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	srand((unsigned) time(NULL));

///Motif search on a small set of nucleotide sequences.
	unsigned int t = 3;		//number of input sequences
	unsigned int n = 6;		//length of sequence
	unsigned int l = 4;		//length of motif
	unsigned int d = 1;		//number of substitutions
	bool is_exact = true;	//size of Hamming distance
	unsigned int h = 0;		//size of the neighborhood considering at first

	String<DnaString> dataset;
	appendValue(dataset,DnaString("ACAGCA"));
	appendValue(dataset,DnaString("AGGCAG"));
	appendValue(dataset,DnaString("TCAGTC"));
	
///Application of ePatternBranching (h=0)
	MotifFinder<Dna, EPatternBranching> finder_epb1(t,l,d,is_exact,h);
	findMotif(finder_epb1,dataset,OMOPS());
	displayResult(finder_epb1); 

///Application of ePatternBranching (h=0)
	MotifFinder<Dna, EPatternBranching> finder_epb2(t,l,d,is_exact,h);
	findMotif(finder_epb2,dataset,OOPS());
	displayResult(finder_epb2); 

///Application of PMS1-ZOOPS 
	MotifFinder<Dna, PMS1> finder_pms1(l,d,is_exact);
	findMotif(finder_pms1,dataset,ZOOPS());
	displayResult(finder_pms1); 

///Application of PMSP-TCM
	MotifFinder<Dna, PMSP> finder_pmsp(l,d,is_exact);
	findMotif(finder_pmsp,dataset,TCM());
	displayResult(finder_pmsp); 
	

///Application of PROJECTION-OOPS
	unsigned int m = t*(n-l+1);
    MotifFinder<Dna, Projection> finder_proj(t,l,m,d,is_exact);
	findMotif(finder_proj, dataset, OOPS());
	displayResult(finder_proj);

///Application of PROJECTION-OMOPS
    MotifFinder<Dna, Projection> finder_proj_omops(t,l,m,d,is_exact);
	findMotif(finder_proj_omops, dataset, OMOPS());
	displayResult(finder_proj_omops);

///Application of PROJECTION-ZOOPS
	MotifFinder<Dna, Projection> finder_proj_zoops(t,l,m,d,is_exact);
	findMotif(finder_proj_zoops, dataset, ZOOPS());
	displayResult(finder_proj_zoops);
	
///Application of PROJECTION-TCM
    MotifFinder<Dna, Projection> finder_proj_tcm(t,l,m,d,is_exact);
	findMotif(finder_proj_tcm, dataset, TCM());
	displayResult(finder_proj_tcm);

	SEQAN_TREPORT("TEST END")
	return 0;
}
*/