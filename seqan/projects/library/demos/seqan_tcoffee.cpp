#include <iostream>
#include <fstream>
#include <seqan/graph.h>

using namespace seqan;


template<typename TName>
class ConfigOptions {
public:
	TName sequencefile;
	TName libfile;
	TName matchfile;
	TName treefile;
	TName outfile;
	TName method;
	TName output;

	ConfigOptions() : outfile("out.fasta"), output("fasta") {}
};

inline bool
printErrorForCmd()
{
	std::cerr << "Usage: seqan_tcoffee -seq <sequence file> [OPTIONS]" << std::endl;
	return false;
}

template<typename TConfigOptions>
inline bool
parseCommandLine(int argc, const char *argv[], TConfigOptions& cfgOpt) {
	if (argc <= 2) return printErrorForCmd();
	for(int i = 1; i < argc; ++i) {
		if (argv[i][0] == '-') {
			// parse option
			if (strcmp(argv[i], "-seq")==0) {
				if (i + 1 == argc) return printErrorForCmd();
				++i;
				cfgOpt.sequencefile = argv[i];
				continue;
			}
			if (strcmp(argv[i], "-lib")==0) {
				if (i + 1 == argc) return printErrorForCmd();
				++i;
				cfgOpt.libfile = argv[i];
				continue;
			}
			if (strcmp(argv[i], "-matches")==0) {
				if (i + 1 == argc) return printErrorForCmd();
				++i;
				cfgOpt.matchfile = argv[i];
				continue;
			}
			if (strcmp(argv[i], "-usetree")==0) {
				if (i + 1 == argc) return printErrorForCmd();
				++i;
				cfgOpt.treefile = argv[i];
				continue;
			}
			if (strcmp(argv[i], "-outfile")==0) {
				if (i + 1 == argc) return printErrorForCmd();
				++i;
				cfgOpt.outfile = argv[i];
				continue;
			}
			if (strcmp(argv[i], "-method")==0) {
				if (i+1 == argc) return printErrorForCmd();
				++i;
				cfgOpt.method = argv[i];
				continue;
			}
			if (strcmp(argv[i], "-output")==0) {
				if (i+1 == argc) return printErrorForCmd();
				++i;
				cfgOpt.output = argv[i];
				continue;
			}
		}
	}
	return true;
}

template<typename TConfigOptions>
inline int
dnaAlignment(TConfigOptions& cfgOpt) {
	// Read the sequences
	typedef String<Dna> TString;
	StringSet<TString, Owner<> > origStrSet;
	String<String<char> > names;
	if (length(cfgOpt.sequencefile)) {
		// Read the sequences from sequencefile
		_alignImportSequences(cfgOpt.sequencefile, "", "", origStrSet, names);
	} else {
		std::cerr << "No sequences provided." << std::endl;
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
	if (cfgOpt.method == "dna") tCoffeeDnaAlignment(strSet, gOut);
	else {
		String<char> lib = cfgOpt.matchfile;
		tCoffeeLongDnaAlignment(strSet, names, lib, gOut);
	}

	// Output alignment
	if (cfgOpt.output == "fasta") {
		std::fstream strm;
		strm.open(toCString(cfgOpt.outfile), std::ios_base::out | std::ios_base::trunc);
		write(strm,gOut,names,FastaFormat());
		strm.close();
	} else if (cfgOpt.output == "msf") {
		std::fstream strm;
		strm.open(toCString(cfgOpt.outfile), std::ios_base::out | std::ios_base::trunc);
		write(strm,gOut,names,MsfFormat());
		strm.close();
	}

	return 0;
}


int main(int argc, const char *argv[]) {
	// Parse the command line options
	ConfigOptions<String<char> > cfgOpt;
	if (!parseCommandLine(argc, argv, cfgOpt)) return -1;
	if ((cfgOpt.method == "dna") || (cfgOpt.method == "genome")) return dnaAlignment(cfgOpt);
	
	// Read the sequences
	typedef String<AminoAcid> TString;
	StringSet<TString, Owner<> > origStrSet;
	String<String<char> > names;
	if ((!length(cfgOpt.sequencefile)) && (length(cfgOpt.libfile))) {
		// Read sequences from library
		std::fstream strm;
		strm.open(toCString(cfgOpt.libfile), std::ios_base::in | std::ios_base::binary);
		read(strm,origStrSet,names,TCoffeeLib());	// Read identifiers and strings
		strm.close();			
	} else if (length(cfgOpt.sequencefile)) {
		// Read the sequences from sequencefile
		_alignImportSequences(cfgOpt.sequencefile, "", "", origStrSet, names);
	} else {
		printErrorForCmd();
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
	if (cfgOpt.output == "tc_lib") {
		Blosum62 score_type_global(-1,-11);
		Blosum62 score_type_local(-2,-8);
		TGraph lib(strSet);
		String<Pair<TId, TId> > pList;
		selectPairsForLibraryGeneration(lib, pList);
		if (cfgOpt.method == "global") {
			generatePrimaryLibrary(lib, pList, score_type_global, GlobalPairwise_Library() );
		} else if (cfgOpt.method == "local") {
			generatePrimaryLibrary(lib, pList, score_type_local, LocalPairwise_Library() );
		} else if (cfgOpt.method == "motif") {
			if (length(cfgOpt.matchfile)) {
				std::fstream strm_lib;
				strm_lib.open(toCString(cfgOpt.matchfile), std::ios_base::in | std::ios_base::binary);
				read(strm_lib, lib, names, MemeMotif());	// Read library
				strm_lib.close();
			}
		} else if (cfgOpt.method == "blast") {
			if (length(cfgOpt.matchfile)) {
				std::fstream strm_lib;
				strm_lib.open(toCString(cfgOpt.matchfile), std::ios_base::in | std::ios_base::binary);
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
		strm.open(toCString(cfgOpt.outfile), std::ios_base::out | std::ios_base::trunc);
		write(strm,lib,names,TCoffeeLib());
		strm.close();
	} 
	else if (length(cfgOpt.libfile)) {
		TGraph g(strSet);
		std::fstream strm_lib;
		strm_lib.open(toCString(cfgOpt.libfile), std::ios_base::in | std::ios_base::binary);
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
		} else if (length(cfgOpt.treefile)) {
			std::fstream strm_tree;
			strm_tree.open(toCString(cfgOpt.treefile), std::ios_base::in | std::ios_base::binary);
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
		tCoffeeProteinAlignment(strSet, names, cfgOpt.matchfile, gOut);
		//tCoffeeProteinAlignment(strSet, names, gOut);
	}

	if (cfgOpt.output == "fasta") {
		// Output alignment
		std::fstream strm;
		strm.open(toCString(cfgOpt.outfile), std::ios_base::out | std::ios_base::trunc);
		write(strm,gOut,names,FastaFormat());
		strm.close();
	} else if (cfgOpt.output == "msf") {
		std::fstream strm;
		strm.open(toCString(cfgOpt.outfile), std::ios_base::out | std::ios_base::trunc);
		write(strm,gOut,names,MsfFormat());
		strm.close();
	}

	return 0;
}
