#include <iostream>
#include <fstream>
#include <cstdio>
#include <seqan/graph_msa.h>

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


template<typename TGraph, typename TVertexSet, typename TCliques>
inline void
enumNextMaximalClique(TGraph& g, 
					  TVertexSet& setC, 
					  TVertexSet& setN, 
					  TVertexSet& setP,
					  TCliques& cliques) 
{
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, AdjacencyIterator>::Type TAdjacencyIterator;
	typedef typename TVertexSet::const_iterator TSetIter;
	// Check the bound
	bool oldClique = false;
	for(TSetIter itP = setP.begin(); itP != setP.end(); ++itP) {
		TVertexSet copySetN = setN;
		for(TAdjacencyIterator itA(g, *itP); !atEnd(itA); ++itA) copySetN.erase(*itA);
		if (copySetN.empty()) {
			oldClique = true;
			break;
		}
	}
	if (!oldClique) {
		while (!setN.empty()) {
			TVertexDescriptor v = *(setN.begin());
			setC.insert(v);
			setN.erase(v);
			TVertexSet setNN;
			TVertexSet setNP;
			for(TAdjacencyIterator itA(g,v);!atEnd(itA); ++itA) {
				if (setN.find(*itA) != setN.end()) setNN.insert(*itA);
				else if (setP.find(*itA) != setP.end()) setNP.insert(*itA);
			}
			if ((setNN.empty()) && (setNP.empty())) appendValue(cliques, setC);
			else enumNextMaximalClique(g, setC, setNN, setNP, cliques);
			setC.erase(v);
			setP.insert(v);
		}
	}
}

template<typename TGraph, typename TCliques>
inline void
enumMaximalCliques(TGraph& g,
				   TCliques& cliques) {
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef std::set<TVertexDescriptor> TVertexSet;
	TVertexSet setN;
	TVertexSet setP;
	TVertexSet setC;
	for(TVertexIterator itV(g);!atEnd(itV);++itV) setN.insert(value(itV));
	enumNextMaximalClique(g, setC, setN, setP, cliques);
}

int main(int argc, const char *argv[]) {
	typedef Graph<Undirected<> > TUndirectedGraph;
	typedef Size<TUndirectedGraph>::Type TSize;
	typedef VertexDescriptor<TUndirectedGraph>::Type TVertexDescriptor;
	typedef std::set<TVertexDescriptor> TVertexSet;

	TUndirectedGraph g;
	addVertex(g);addVertex(g);addVertex(g);addVertex(g);addVertex(g);addVertex(g);addVertex(g);
	addEdge(g, 0, 1);addEdge(g, 0, 2);addEdge(g, 0, 3);addEdge(g, 0, 4);
	addEdge(g, 1, 2);addEdge(g, 1, 3);addEdge(g, 1, 4);
	addEdge(g, 3, 4);
	addEdge(g, 3, 5);
	addEdge(g, 5, 6);
	addEdge(g, 4, 6);

	String<TVertexSet> cliques;
	enumMaximalCliques(g, cliques);
	
	std::cout << g << std::endl;
	std::cout << "Cliques" << std::endl;
	for(TSize i = 0; i < length(cliques); ++i) {
		for(TVertexSet::const_iterator vIt = (value(cliques, i)).begin(); vIt != (value(cliques, i)).end(); ++vIt) {
			std::cout << *vIt << ',';
		}
		std::cout << std::endl;
	}
	


	exit(0);






















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
