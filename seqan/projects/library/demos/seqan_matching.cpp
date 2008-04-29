#include <iostream>
#include <fstream>
#include <cstdio>
#include <seqan/graph.h>

using namespace seqan;

/*
template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TConfigOptions>
inline void
globalMatching(StringSet<TString, TSpec> const& seqSet,
			   StringSet<TName, TSpec2> const& nameSet,
			   TConfigOptions& cfgOpt)
{
	SEQAN_CHECKPOINT
	
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;

	// Score objects
	Score<int> score_type_global = Score<int>(5,-4,-4,-14);
	Score<int> score_type_local = Score<int>(5,-4,-4,-14);

	TGraph g(seqSet);
	String<Pair<TId, TId> > pList;
	selectPairsForLibraryGeneration(g, pList);

	std::fstream strm_lib;
	strm_lib.open(toCString(value(cfgOpt, "matches")), std::ios_base::in | std::ios_base::binary);
	read(strm_lib, g, nameSet, BlastLib());	// Read library
	strm_lib.close();
	//generatePrimaryLibrary(g, score_type_global, 5, Kmer_Library() );

	String<double> distanceMatrix;
	getDistanceMatrix(g,distanceMatrix,LibraryDistance() );

	tripletLibraryExtension(g);

	
	// ... and normal progressive alignment with guide tree
	Graph<Tree<double> > guideTree;
	//upgmaTree(distanceMatrix, guideTree);	
	slowNjTree(distanceMatrix, guideTree);

	TGraph gOut(seqSet);
	progressiveAlignment(g, guideTree, gOut);

	TGraph gOut2(seqSet);
	progressiveMatching(g, guideTree, gOut2);


	std::fstream strm_lib1;
	strm_lib1.open("/home/takifugu/rausch/matches/reads/library.txt", std::ios_base::out | std::ios_base::trunc);
	write(strm_lib1, g, nameSet, BlastLib());	// Write library
	strm_lib1.close();

	std::fstream strm_lib2;
	strm_lib2.open("/home/takifugu/rausch/matches/reads/alignment.txt", std::ios_base::out | std::ios_base::trunc);
	write(strm_lib2, gOut, nameSet, BlastLib());	// Write library
	strm_lib2.close();

	std::fstream strm_lib3;
	strm_lib3.open("/home/takifugu/rausch/matches/reads/matching.txt", std::ios_base::out | std::ios_base::trunc);
	write(strm_lib3, gOut2, nameSet, BlastLib());	// Write library
	strm_lib3.close();

	clear(guideTree);
	clear(distanceMatrix);
	clear(g);
}

*/


int main(/*int argc, const char *argv[]*/) {

	return 0;
}
