#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_VERBOSE

// Test path
#define TEST_PATH "projects/tests/graph_align/"

// External / STL
#include <iostream>
#include <fstream>
#include <string>

#include <seqan/map.h>


// SeqAn Includes
#include <seqan/graph_align.h>


// Test files
#include "test_graph_align.h"


using namespace seqan;


////////////////////////////////////////////////////////////////////////////////
//
//#include <seqan/graph_msa.h>
//void simulateSequences() 
//{
//	SEQAN_TREPORT("TEST BEGIN")
//
//	typedef String<Dna> TString;
//	typedef StringSet<TString, Dependent<> > TStringSet;
//	typedef Graph<Alignment<TStringSet, int> > TGraph;
//	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
//	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
//	typedef	Id<TStringSet>::Type TId;
//	Score<int> score_type = Score<int>(1,-1,-1,-2);
//	TStringSet str;
//	TString str0("TGCA");
//	appendValue(str, str0);
//	TString str1("TTCAGTTA");
//	appendValue(str, str1);
//
//	
//	mtRandInit();
//	while(true) {
//		clear(str);
//		clear(str0);
//		clear(str1);
//
//		unsigned int stop = (unsigned int) mtRand() % 40 + 1;
//		for (unsigned int i=0; i<stop ; ++i) {
//			append(str0, (Byte) mtRand() % 4 );
//		}
//		stop = (unsigned int) mtRand() % 40 + 1;
//		for (unsigned int i=0; i<stop; ++i) {
//			append(str1,  (Byte) mtRand() % 4);
//		}
//		std::cout << str0 << std::endl;
//		std::cout << str1 << std::endl;
//		appendValue(str, str0);
//		appendValue(str, str1);
//		std::cout << "Length Seq0: " << length(str0) << std::endl;
//		std::cout << "Length Seq1: " << length(str1) << std::endl;
//		
//		TGraph g(str);
//		int score = globalAlignment(g, score_type, Gotoh());
//		std::cout << g << std::endl;
//		std::cout << score << std::endl;
//		std::cout << sumOfPairsScore(g, score_type) << std::endl;
//		
//		int diag2 = length(str0);
//		int diag1 = -1 * length(str1);
//		TGraph g2(str);
//		score = globalAlignment(g2, score_type, diag1, diag2, BandedGotoh());
//		std::cout << g2 << std::endl;
//		std::cout << score << std::endl;
//		std::cout << sumOfPairsScore(g2, score_type) << std::endl;
//
//		if (mtRand() % 2 == 0) {
//			diag2 = mtRand() % length(str0);
//			if (diag2 - 5 > diag1) diag1 = diag2 - 5;
//		} else {
//			diag1 = -1 * (mtRand() % length(str1));
//			if (diag1 + 5 < diag2) diag2 = diag1 + 5;
//		}
//		std::cout << diag1 << ',' << diag2 << std::endl;
//		TGraph gOverlap(str);
//		typedef Fragment<> TFragment;
//		typedef String<TFragment> TFragmentString;
//		typedef Iterator<TFragmentString>::Type TFragmentStringIter;
//		TFragmentString matches;
//		globalAlignment(matches, stringSet(gOverlap), score_type, AlignConfig<true, true, true, true>(), diag1, diag2, BandedGotoh());
//		matchRefinement(matches,stringSet(gOverlap),score_type,gOverlap);
//		std::cout << gOverlap << std::endl;
//		String<char> mat;
//		if (convertAlignment(gOverlap, mat) == false) {
//			std::cout << "Cannot convert alignment" << std::endl;
//			exit(0);
//		}
//		if (sumOfPairsScore(g, score_type) !=  sumOfPairsScore(g2, score_type)) {
//			std::cout << "Different sum of pairs score" << std::endl;
//			exit(0);
//		}
//	}
//	
//	//TGraph g(str);
//	//int score = globalAlignment(g, score_type, AlignConfig<true, true, true, true>(), Gotoh());
//	//std::cout << g << std::endl;
//	//std::cout << score << std::endl;
//	//int diag2 = length(str0);
//	//int diag1 = -1 * length(str1);
//	//diag2 = -2;
//	//diag1 = -4;
//	//clearVertices(g);
//	//score = globalAlignment(g, score_type, AlignConfig<true, true, true, true>(), diag1, diag2, BandedGotoh());
//	//std::cout << g << std::endl;
//	//std::cout << score << std::endl;
//}



//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")
	//simulateSequences();

	// Redirect std::cout
    std::ofstream file(TEST_PATH "redirect.txt");
    std::streambuf* strm_puffer = std::cout.rdbuf();
    std::cout.rdbuf(file.rdbuf());
	
	Test_GraphAlignment();		// Test Graph Alignment

	// Restore std::cout
	std::cout.rdbuf(strm_puffer);

	SEQAN_TREPORT("TEST END")

	return 0;
}
