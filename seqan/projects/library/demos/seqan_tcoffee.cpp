#include <iostream>
#include <fstream>
#include <seqan/graph.h>

using namespace seqan;

int main(int argc, const char *argv[]) {
	// The original string set
	typedef String<AminoAcid> TString;
	StringSet<TString, Owner<> > origStrSet;
	typedef String<char> TFileName;
	String<TFileName> names;
	TFileName outfile = "out.fasta";
	TFileName sequencefile = "";
	TFileName libfile = "";
	TFileName treefile = "";
	String<char> method = "";
	String<char> output = "";

	// Command-line parsing
	if (argc <= 1) return 0;
	for(int i = 1; i < argc; ++i) {
		//std::cout << argv[i] << std::endl;
		if (argv[i][0] == '-') {
			// parse option
			if (strcmp(argv[i], "-seq")==0) {
				if (i + 1 == argc) {
					//printHelp(argc, argv);
					return 0;
				}
				++i;
				sequencefile = argv[i];
				continue;
			}
			if (strcmp(argv[i], "-lib")==0) {
				if (i + 1 == argc) {
					//printHelp(argc, argv);
					return 0;
				}
				++i;
				libfile = argv[i];
				continue;
			}
			if (strcmp(argv[i], "-usetree")==0) {
				if (i + 1 == argc) {
					//printHelp(argc, argv);
					return 0;
				}
				++i;
				treefile = argv[i];
				continue;
			}
			if (strcmp(argv[i], "-outfile")==0) {
				if (i + 1 == argc) {
					//printHelp(argc, argv);
					return 0;
				}
				++i;
				outfile = argv[i];
				continue;
			}
			if (strcmp(argv[i], "-method")==0) {
				if (i+1 == argc) {
					//printHelp(argc, argv);
					return 0;
				}
				++i;
				method = argv[i];
				continue;
			}
			if (strcmp(argv[i], "-output")==0) {
				if (i+1 == argc) {
					//printHelp(argc, argv);
					return 0;
				}
				++i;
				output = argv[i];
				continue;
			}
		}
	}

	if (output == "tc_lib") {
		std::cout << "Sequence File: " << sequencefile << std::endl;
		std::cout << "Output File: " << outfile << std::endl;
		std::cout << "Method: " << method << std::endl;
		std::cout << "Output format: " << output << std::endl;
	} else {
		std::cout << "Sequence File: " << sequencefile << std::endl;
		std::cout << "Library file: " << libfile << std::endl;
		std::cout << "Tree file: " << treefile << std::endl;
		std::cout << "Output File: " << outfile << std::endl;
	}

	if ((length(libfile)) &&
		(!length(sequencefile))) {
			// Read sequences from library
			std::fstream strm;
			strm.open(toCString(libfile), std::ios_base::in | std::ios_base::binary);
			read(strm,origStrSet,names,TCoffeeLib());	// Read identifiers and strings
			strm.close();			
	} else {
		// Read the sequences from sequencefile
		_alignImportSequences(sequencefile, "", "", origStrSet, names);
	}

	// Make dependent string set
	typedef StringSet<TString, Dependent<> > TStringSet;
	TStringSet strSet;
	for(unsigned int i = 0; i<length(origStrSet); ++i) appendValue(strSet, origStrSet[i]);


	// Start the alignment
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	Graph<Alignment<TStringSet, void, WithoutEdgeId> > gOut(strSet);
	if (output == "tc_lib") {
		Blosum62 score_type_global(-1,-11);
		Blosum62 score_type_local(-2,-8);
		TGraph lib(strSet);
		if (method == "global") {
			generatePrimaryLibrary(lib, score_type_global, GlobalPairwise_Library() );
		} else if (method == "local") {
			generatePrimaryLibrary(lib, score_type_local, LocalPairwise_Library() );
		} else {
			TGraph lib1(strSet);
			generatePrimaryLibrary(lib1, score_type_global, GlobalPairwise_Library() );
			TGraph lib2(strSet);
			generatePrimaryLibrary(lib2, score_type_local, LocalPairwise_Library() );
	
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
		strm.open(toCString(outfile), std::ios_base::out | std::ios_base::trunc);
		write(strm,lib,names,TCoffeeLib());
		strm.close();
	} else if (length(libfile)) {
		TGraph g(strSet);
		std::fstream strm_lib;
		strm_lib.open(toCString(libfile), std::ios_base::in | std::ios_base::binary);
		read(strm_lib,g,TCoffeeLib());	// Read library
		strm_lib.close();

		//// Debug code
		//typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
		//typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
		//TEdgeIterator it(g);
		//for(;!atEnd(it);++it) {
		//	TVertexDescriptor sV = sourceVertex(it);
		//	TVertexDescriptor tV = targetVertex(it);
		//	std::cout << label(g, sV) << ',' << label(g,tV) << std::endl;	
		//}

		// Triplet library extension
		tripletLibraryExtension(g);

		// Build the guide tree
		Graph<Tree<double> > guideTree;
		if (length(treefile)) {
			std::fstream strm_tree;
			strm_tree.open(toCString(treefile), std::ios_base::in | std::ios_base::binary);
			read(strm_tree,guideTree,names,NewickFormat());	// Read newick tree
			strm_tree.close();
		} else {
			// Pairwise Distances
			String<double> distanceMatrix;
			getDistanceMatrix(g, distanceMatrix, LibraryDistance());
			upgmaTree(distanceMatrix, guideTree);
		}

		// Progressive alignment
		progressiveAlignment(g, guideTree, gOut);
	} else {
		if (length(treefile)) {
			// Score objects
			Blosum62 score_type_global(-1,-11);
			Blosum62 score_type_local(-2,-8);
	
			// Generate a primary library, i.e., all global pairwise alignments
			TGraph lib1(strSet);
			generatePrimaryLibrary(lib1, score_type_global, GlobalPairwise_Library() );
	
			TGraph lib2(strSet);
			generatePrimaryLibrary(lib2, score_type_local, LocalPairwise_Library() );
	
			// Weighting of libraries (Signal addition)
			TGraph g(strSet);
			String<TGraph*> libs;
			appendValue(libs, &lib1);
			appendValue(libs, &lib2);
			combineGraphs(g, libs);

			// Clear the old libraries
			clear(lib1);
			clear(lib2);

			//Retrieve the guide tree
			Graph<Tree<double> > guideTree;
			std::fstream strm_tree;
			strm_tree.open(toCString(treefile), std::ios_base::in | std::ios_base::binary);
			read(strm_tree,guideTree,names,NewickFormat());	// Read newick tree
			strm_tree.close();

			unsigned int nSeq = length(strSet);
			unsigned int threshold = 30;
			if (nSeq < threshold) {
				tripletLibraryExtension(g);
				progressiveAlignment(g, guideTree, gOut);		
			} else {
				progressiveAlignment(g, guideTree, gOut, threshold);
			}

			//write(std::cout,guideTree,names,NewickFormat());	// Write newick tree

			clear(guideTree);
		} else {
			// Full-blown alignment
			tCoffeeProteinAlignment(strSet, gOut);
		}
	}

	if (output != "tc_lib") {

		// Output alignment
		std::fstream strm; // Alignment graph as fasta
		strm.open(toCString(outfile), std::ios_base::out | std::ios_base::trunc);
		write(strm,gOut,names,FastaFormat());
		strm.close();
	}

	return 0;
}
