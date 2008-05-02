#include <iostream>
#include <fstream>
#include <cstdio>
#include <seqan/graph.h>

using namespace seqan;


template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TConfigOptions, typename TOutGraph, typename TEdgeMap>
inline void
globalMatching(StringSet<TString, TSpec> const& seqSet,
			   StringSet<TName, TSpec2> const& nameSet,
			   TConfigOptions& cfgOpt,
			   TOutGraph& gOut,
			   TEdgeMap& edgeMapOut)
{
	SEQAN_CHECKPOINT
	
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef TOutGraph TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;

	// Read all matches and reversed matches
	TGraph g(seqSet);
	TEdgeMap edgeMap;
	std::fstream strm_lib;
	strm_lib.open(toCString(value(cfgOpt, "matches")), std::ios_base::in | std::ios_base::binary);
	read(strm_lib, g, nameSet, edgeMap, MummerLib());	// Read library
	strm_lib.close();

	// Write out that library
	std::fstream strm_lib1;
	strm_lib1.open(toCString(value(cfgOpt, "library")), std::ios_base::out | std::ios_base::trunc);
	write(strm_lib1, g, nameSet, edgeMap, BlastLib());	// Write library
	strm_lib1.close();

	// Compute a rough and roudy distance matrix
	String<double> distanceMatrix;
	getDistanceMatrix(g,distanceMatrix,LibraryDistance() );

	// Perform the triplet extension
	tripletLibraryExtension(g);
	
	// Get a guide tree
	Graph<Tree<double> > guideTree;
	slowNjTree(distanceMatrix, guideTree);

	// Progressive matching
	progressiveMatching(g, guideTree, edgeMap, gOut, edgeMapOut);

	// Clean-up
	clear(guideTree);
	clear(distanceMatrix);
	clear(g);
}

int main(int argc, const char *argv[]) {
	//////////////////////////////////////////////////////////////////////////////
	// Command line parsing
	//////////////////////////////////////////////////////////////////////////////
	
	// Set the keys
	typedef String<char> TKey;
	typedef String<char> TValue;
	typedef Size<TKey>::Type TSize;
	ConfigOptions<TKey, TValue> cfgOpt;
	TKey keys[] = {"seq","outfile","output", "matches", "library"};
	assignKeys(cfgOpt, keys, 5);
	// Set default options
	assign(cfgOpt, "output", "fasta");
	assign(cfgOpt, "library", "library.txt");
	assign(cfgOpt, "outfile", "matching.txt");
	// Help Message
	String<char> helpMsg;
	append(helpMsg, "Usage: seqan_matching -seq <FASTA Sequence File> [ARGUMENTS]\n");
	assignHelp(cfgOpt, helpMsg);
	
	if (argc < 2) {	std::cerr << valueHelp(cfgOpt) << std::endl; return -1; }
	if (!parseCmdLine(argc, argv, cfgOpt)) return -1;

	//////////////////////////////////////////////////////////////////////////////
	// Read the sequences
	//////////////////////////////////////////////////////////////////////////////
	
	typedef String<Dna> TSequence;
	typedef String<char> TName;
	StringSet<TSequence, Owner<> > origStrSet;
	StringSet<TName> names;
	_alignImportSequences(value(cfgOpt, "seq"), origStrSet, names);
	typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
	typedef Size<TDepSequenceSet>::Type TSize;
	TDepSequenceSet strSet(origStrSet);
	
	//////////////////////////////////////////////////////////////////////////////
	// Alignment
	//////////////////////////////////////////////////////////////////////////////

	typedef Graph<Alignment<TDepSequenceSet, TSize> > TGraph;
	Graph<Alignment<TDepSequenceSet, TSize> > gOut(strSet);
	String<bool> edgeMap;
	globalMatching(strSet, names, cfgOpt, gOut, edgeMap);

	////////////////////////////////////////////////////////////////////////////
	// Alignment output
	////////////////////////////////////////////////////////////////////////////

	std::fstream strm_lib3;
	strm_lib3.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::trunc);
	write(strm_lib3, gOut, names, edgeMap, BlastLib());	// Write library
	strm_lib3.close();

	//std::fstream strm;
	//strm.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::trunc);
	//write(strm,gOut,names,XMFA());
	//strm.close();
	
	return 0;
}
