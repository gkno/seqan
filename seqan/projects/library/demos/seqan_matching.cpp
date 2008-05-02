#include <iostream>
#include <fstream>
#include <cstdio>
#include <seqan/graph.h>

using namespace seqan;


template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TConfigOptions, typename TOutGraph>
inline void
globalMatching(StringSet<TString, TSpec> const& seqSet,
			   StringSet<TName, TSpec2> const& nameSet,
			   TConfigOptions& cfgOpt,
			   TOutGraph& gOut)
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

	progressiveAlignment(g, guideTree, gOut);

	TGraph gOut2(seqSet);
	progressiveMatching(g, guideTree, gOut2);


	std::fstream strm_lib1;
	strm_lib1.open("/home/takifugu/rausch/matches/reads/library.txt", std::ios_base::out | std::ios_base::trunc);
	write(strm_lib1, g, nameSet, BlastLib());	// Write library
	strm_lib1.close();

	std::fstream strm_lib3;
	strm_lib3.open("/home/takifugu/rausch/matches/reads/matching.txt", std::ios_base::out | std::ios_base::trunc);
	write(strm_lib3, gOut2, nameSet, BlastLib());	// Write library
	strm_lib3.close();

	clear(guideTree);
	clear(distanceMatrix);
	clear(g);
}

int main(int argc, const char *argv[]) {
	////////////////////////////////////////////////////////////////////////////////
	//// Command line parsing
	////////////////////////////////////////////////////////////////////////////////
	//
	//// Set the keys
	//typedef String<char> TKey;
	//typedef String<char> TValue;
	//typedef Size<TKey>::Type TSize;
	//ConfigOptions<TKey, TValue> cfgOpt;
	//TKey keys[] = {"seq","outfile","output", "matches"};
	//assignKeys(cfgOpt, keys, 4);
	//// Set default options
	//assign(cfgOpt, "output", "fasta");
	//assign(cfgOpt, "outfile", "out.fasta");
	//// Help Message
	//String<char> helpMsg;
	//append(helpMsg, "Usage: seqan_matching -seq <FASTA Sequence File> [ARGUMENTS]\n");
	//assignHelp(cfgOpt, helpMsg);
	//
	//if (argc < 2) {	std::cerr << valueHelp(cfgOpt) << std::endl; return -1; }
	//if (!parseCmdLine(argc, argv, cfgOpt)) return -1;

	////////////////////////////////////////////////////////////////////////////////
	//// Read the sequences
	////////////////////////////////////////////////////////////////////////////////
	//
	//typedef String<Dna> TSequence;
	//typedef String<char> TName;
	//StringSet<TSequence, Owner<> > origStrSet;
	//StringSet<TName> names;
	//_alignImportSequences(value(cfgOpt, "seq"), origStrSet, names);
	//typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
	//typedef Size<TDepSequenceSet>::Type TSize;
	//TDepSequenceSet strSet(origStrSet);
	//
	////////////////////////////////////////////////////////////////////////////////
	//// Alignment
	////////////////////////////////////////////////////////////////////////////////

	//typedef Graph<Alignment<TDepSequenceSet, TSize> > TGraph;
	//Graph<Alignment<TDepSequenceSet, void, WithoutEdgeId> > gOut(strSet);
	//globalMatching(strSet, names, cfgOpt, gOut);

	//////////////////////////////////////////////////////////////////////////////
	// Alignment output
	//////////////////////////////////////////////////////////////////////////////

	//std::fstream strm;
	//strm.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::trunc);
	//write(strm,gOut,names,XMFA());
	//strm.close();
	
	return 0;
}
