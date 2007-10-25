#define SEQAN_DEBUG
#define SEQAN_TEST

#include "seqan/find_motif.h"

using namespace seqan;

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	/*typedef FrequencyDistribution<Dna> TFrequencyDistribution;
	TFrequencyDistribution fd('a');
	//Iterator<TFrequencyDistribution>::Type fd_begin = begin(fd);
	//Iterator<TFrequencyDistribution>::Type fd_end = end(fd);
	if(std::find(begin(fd),end(fd),static_cast<double>(0.5))!=end(fd))
	{
		std::cout << "TEST\n";
	}*/


	//Motif search on a small set of nucleotide sequences.
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

	//Application of ePatternBranching (h=0)
	MotifFinder<Dna, EPatternBranching> finder_epb1(t,l,d,is_exact,h);
	findMotif(finder_epb1,dataset,OMOPS());
	displayResult(finder_epb1); 
	//output: 
	//[0] AGCC
	//[1] CCAG

	//Application of ePatternBranching (h=0)
	MotifFinder<Dna, EPatternBranching> finder_epb2(t,l,d,is_exact,h);
	findMotif(finder_epb2,dataset,OOPS());
	displayResult(finder_epb2); 
	//output: 
	//[0] AGCC
	//[1] CCAG

	
	//Application of PMS1-ZOOPS 
	MotifFinder<Dna, PMS1> finder_pms1(l,d,is_exact);
	findMotif(finder_pms1,dataset,ZOOPS());
	displayResult(finder_pms1); 
	//output: 
	//[0] AAGC
	//...
	//[13] TGCA

	//Application of PMSP-TCM
	MotifFinder<Dna, PMSP> finder_pmsp(l,d,is_exact);
	findMotif(finder_pmsp,dataset,TCM());
	displayResult(finder_pmsp); 
	//output: 
	//[0] AAGC
	//...
	//[13] TGCA


	/*
	//Application of PROJECTION-OOPS
	unsigned int m_total = t*(n-l+1);
    MotifFinder<Dna, Projection> finder_proj(t,l,m_total,d);
	findMotif(finder_proj, dataset, OOPS());
	displayResult(finder_proj);

	//Application of PROJECTION-TCM
    MotifFinder<Dna, Projection> finder_proj_tcm(t,l,m_total,d);
	findMotif(finder_proj_tcm, dataset, TCM());
	displayResult(finder_proj);
	*/

	
	SEQAN_TREPORT("TEST END")
	return 0;
}