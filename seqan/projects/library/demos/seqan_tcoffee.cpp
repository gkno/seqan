#include <iostream>
#include <fstream>
#include <seqan/graph.h>

using namespace seqan;





template<typename TConfigOptions>
inline int
dnaAlignment(TConfigOptions& cfgOpt) {
	// Read the sequences
	typedef String<Dna> TString;
	StringSet<TString, Owner<> > origStrSet;
	String<String<char> > names;
	// Read the sequences from sequencefile
	_alignImportSequences(value(cfgOpt, "seq"), "", "", origStrSet, names);
	
	// Make a dependent string set
	typedef StringSet<TString, Dependent<> > TStringSet;
	TStringSet strSet;
	for(unsigned int i = 0; i<length(origStrSet); ++i) appendValue(strSet, origStrSet[i]);

	// Start the alignment
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	typedef Id<TGraph>::Type TId;
	Graph<Alignment<TStringSet, void, WithoutEdgeId> > gOut(strSet);
	if (value(cfgOpt, "method") == "dna") tCoffeeDnaAlignment(strSet, gOut);
	else {
		String<char> lib = value(cfgOpt, "matches");
		tCoffeeLongDnaAlignment(strSet, names, lib, gOut);
	}

	// Output alignment
	if (value(cfgOpt, "output") == "fasta") {
		std::fstream strm;
		strm.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::trunc);
		write(strm,gOut,names,FastaFormat());
		strm.close();
	} else if (value(cfgOpt, "output") == "msf") {
		std::fstream strm;
		strm.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::trunc);
		write(strm,gOut,names,MsfFormat());
		strm.close();
	}

	return 0;
}


int main(int argc, const char *argv[]) {
	// Parse the command line options
	typedef String<char> TKey;
	typedef String<char> TValue;
	ConfigOptions<TKey, TValue> cfgOpt;
	TKey keys[] = {"seq","lib","matches","usetree","outfile","method","output"};
	assignKeys(cfgOpt, keys, 7);
	assign(cfgOpt, "output", "fasta");
	assign(cfgOpt, "outfile", "out.fasta");
	assign(cfgOpt, "method", "protein");
	assignHelp(cfgOpt, "Usage: seqan_tcoffee -seq <sequence file> [ARGUMENTS]");
	if (argc < 2) {
		std::cerr << valueHelp(cfgOpt) << std::endl;
		return -1;
	}
	if (!parseCmdLine(argc, argv, cfgOpt)) return -1;
	if ((value(cfgOpt, "method") == "dna") || (value(cfgOpt, "method") == "genome")) return dnaAlignment(cfgOpt);

	// Read the sequences
	typedef String<AminoAcid> TString;
	StringSet<TString, Owner<> > origStrSet;
	String<String<char> > names;
	if ((!length(value(cfgOpt, "seq"))) && (length(value(cfgOpt, "lib")))) {
		// Read sequences from library
		std::fstream strm;
		strm.open(toCString(value(cfgOpt, "lib")), std::ios_base::in | std::ios_base::binary);
		read(strm,origStrSet,names,TCoffeeLib());	// Read identifiers and strings
		strm.close();			
	} else if (length(value(cfgOpt, "seq"))) {
		// Read the sequences from sequencefile
		_alignImportSequences(value(cfgOpt, "seq"), "", "", origStrSet, names);
	} else {
		std::cerr << valueHelp(cfgOpt) << std::endl;
		return -1;
	}

	// Make a dependent string set
	typedef StringSet<TString, Dependent<> > TStringSet;
	TStringSet strSet;
	for(unsigned int i = 0; i<length(origStrSet); ++i) appendValue(strSet, origStrSet[i]);

	// Start the alignment
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	typedef Id<TGraph>::Type TId;
	Graph<Alignment<TStringSet, void, WithoutEdgeId> > gOut(strSet);
	if (value(cfgOpt, "output") == "tc_lib") {
		Blosum62 score_type_global(-1,-11);
		Blosum62 score_type_local(-2,-8);
		TGraph lib(strSet);
		String<Pair<TId, TId> > pList;
		selectPairsForLibraryGeneration(lib, pList);
		if (value(cfgOpt, "method") == "global") {
			generatePrimaryLibrary(lib, pList, score_type_global, GlobalPairwise_Library() );
		} else if (value(cfgOpt, "method") == "local") {
			generatePrimaryLibrary(lib, pList, score_type_local, LocalPairwise_Library() );
		} else if (value(cfgOpt, "method") == "blast") {
			if (length(value(cfgOpt, "matches"))) {
				std::fstream strm_lib;
				strm_lib.open(toCString(value(cfgOpt, "matches")), std::ios_base::in | std::ios_base::binary);
				read(strm_lib, lib, names, BlastLib());	// Read library
				strm_lib.close();
			}
		} else {
			TGraph lib1(strSet);
			generatePrimaryLibrary(lib1, pList, score_type_global, GlobalPairwise_Library() );
			TGraph lib2(strSet);
			generatePrimaryLibrary(lib2, pList, score_type_local, LocalPairwise_Library() );
	
			// Weighting of libraries (Signal addition)
			TGraph g(strSet);
			String<TGraph*> libs;
			appendValue(libs, &lib1);
			appendValue(libs, &lib2);
			combineGraphs(lib, libs);
			//// Only topology
			//combineGraphs(g, true, libs);
		}

		// Write the library
		std::fstream strm;
		strm.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::trunc);
		write(strm,lib,names,TCoffeeLib());
		strm.close();
	} 
	else if (length(value(cfgOpt, "lib"))) {
		TGraph g(strSet);
		std::fstream strm_lib;
		strm_lib.open(toCString(value(cfgOpt, "lib")), std::ios_base::in | std::ios_base::binary);
		read(strm_lib,g,TCoffeeLib());	// Read library
		strm_lib.close();

		//unsigned int nucs = 0;
		//for(unsigned int i = 0; i<length(strSet); ++i) nucs += length(strSet[i]);
		//std::cerr << "Avg. segment length: " << (double) nucs / (double) numVertices(g) << std::endl;

		// Build the guide tree
		Graph<Tree<double> > guideTree;
		if (length(strSet) < 3) {
			typedef VertexDescriptor<Graph<Tree<double> > >::Type TVertexDescriptor;
			TVertexDescriptor v1 = addVertex(guideTree);
			TVertexDescriptor v2 = addVertex(guideTree);
			TVertexDescriptor internalVertex = addVertex(guideTree);
			addEdge(guideTree, internalVertex, v1, 1.0);
			addEdge(guideTree, internalVertex, v2, 1.0);
			guideTree.data_root = internalVertex;
		} else if (length(value(cfgOpt, "usetree"))) {
			std::fstream strm_tree;
			strm_tree.open(toCString(value(cfgOpt, "usetree")), std::ios_base::in | std::ios_base::binary);
			read(strm_tree,guideTree,names,NewickFormat());	// Read newick tree
			strm_tree.close();
		} else {
			// Pairwise Distances
			String<double> distanceMatrix;
			getDistanceMatrix(g, distanceMatrix, LibraryDistance());
			upgmaTree(distanceMatrix, guideTree);
		}

		// Fast progressive alignment
		unsigned int nSeq = length(strSet);
		unsigned int threshold = 30;
		if (nSeq < threshold) {
			// Full triplet...
			tripletLibraryExtension(g);

			progressiveAlignment(g, guideTree, gOut);
		} else {
			// Triplet only on groups of sequences
			progressiveAlignment(g, guideTree, gOut, threshold);
		}
	} 
	//else if (length(cfgOpt.libfile)) {
	//	TGraph lib1(strSet);
	//	std::fstream strm_lib;
	//	strm_lib.open(toCString(cfgOpt.libfile), std::ios_base::in | std::ios_base::binary);
	//	read(strm_lib,lib1,TCoffeeLib());	// Read library
	//	strm_lib.close();

	//	Blosum62 score_type_global(-1,-11);
	//	TGraph lib2(strSet);
	//	String<Pair<TId, TId> > pList;
	//	selectPairsForLibraryGeneration(lib2, pList);
	//	String<double> distanceMatrix;
	//	generatePrimaryLibrary(lib2, pList, distanceMatrix, score_type_global, GlobalPairwise_Library() );
	//
	//	TGraph g(strSet);
	//	String<TGraph*> libs;
	//	appendValue(libs, &lib1);
	//	appendValue(libs, &lib2);
	//	combineGraphs(g, true, libs);

	//	// Triplet library extension
	//	tripletLibraryExtension(g);
	//	// Build the guide tree
	//	Graph<Tree<double> > guideTree;
	//	upgmaTree(distanceMatrix, guideTree);

	//	// Progressive alignment
	//	progressiveAlignment(g, guideTree, gOut);
	//} 
	else {
		// Full-blown alignment
		tCoffeeProteinAlignment(strSet, names, value(cfgOpt, "matches"), gOut);
		//tCoffeeProteinAlignment(strSet, names, gOut);
	}

	if (value(cfgOpt, "output") == "fasta") {
		// Output alignment
		std::fstream strm;
		strm.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::trunc);
		write(strm,gOut,names,FastaFormat());
		strm.close();
	} else if (value(cfgOpt, "output") == "msf") {
		std::fstream strm;
		strm.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::trunc);
		write(strm,gOut,names,MsfFormat());
		strm.close();
	}

	return 0;
}
