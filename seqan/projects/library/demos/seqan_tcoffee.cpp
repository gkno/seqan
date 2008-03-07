#include <iostream>
#include <fstream>
#include <cstdio>
#include <seqan/graph.h>

using namespace seqan;



template<typename TConfigOptions>
inline int
dnaAlignment(TConfigOptions& cfgOpt) {
	//////////////////////////////////////////////////////////////////////////////
	// Read the sequences
	//////////////////////////////////////////////////////////////////////////////
	
	typedef String<Dna> TSequence;
	StringSet<TSequence, Owner<> > origStrSet;
	typedef String<char> TName;
	StringSet<TName> names;
	_alignImportSequences(value(cfgOpt, "seq"), origStrSet, names);
	typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
	TDepSequenceSet strSet(origStrSet);
	
	//////////////////////////////////////////////////////////////////////////////
	// Alignment
	//////////////////////////////////////////////////////////////////////////////
	typedef Graph<Alignment<TDepSequenceSet, unsigned int> > TGraph;
	typedef Id<TGraph>::Type TId;
	Graph<Alignment<TDepSequenceSet, void, WithoutEdgeId> > gOut(strSet);
	if (value(cfgOpt, "method") == "dna") globalAlignment(strSet, names, value(cfgOpt, "matches"), gOut, MSA_Dna() );
	else globalAlignment(strSet, names, value(cfgOpt, "matches"), gOut, MSA_Genome() );

	//////////////////////////////////////////////////////////////////////////////
	// Output
	//////////////////////////////////////////////////////////////////////////////
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
	//////////////////////////////////////////////////////////////////////////////
	// Command line parsing
	//////////////////////////////////////////////////////////////////////////////
	
	// Set the keys
	typedef String<char> TKey;
	typedef String<char> TValue;
	ConfigOptions<TKey, TValue> cfgOpt;
	TKey keys[] = {"seq","lib","matches","usetree","outfile","method","output"};
	assignKeys(cfgOpt, keys, 7);
	// Set default options
	assign(cfgOpt, "output", "fasta");
	assign(cfgOpt, "outfile", "out.fasta");
	assign(cfgOpt, "method", "protein");
	// Help Message
	String<char> helpMsg;
	append(helpMsg, "Usage: seqan_tcoffee -seq <FASTA Sequence File> [ARGUMENTS]\n");
	append(helpMsg, "Arguments\n");
	append(helpMsg, "-matches <File with Segment Matches>\n");
	append(helpMsg, "\tSpecifies a file with gapless segment matches in BLAST Tabular Format (-m 8).\n");
	append(helpMsg, "-method [protein | dna | genome]\n");
	append(helpMsg, "\tSpecifies the type of sequence data, default is protein.\n");
	append(helpMsg, "-outfile <Alignment Filename>\n");
	append(helpMsg, "\tSpecifies the name of the output file, default is out.fasta.\n");
	append(helpMsg, "-output [fasta| msf]\n");
	append(helpMsg, "\tSpecifies the output format, default is fasta.\n");
	assignHelp(cfgOpt, helpMsg);
	if (argc < 2) {
		std::cerr << valueHelp(cfgOpt) << std::endl;
		return -1;
	}
	if (!parseCmdLine(argc, argv, cfgOpt)) return -1;
	if ((value(cfgOpt, "method") == "dna") || (value(cfgOpt, "method") == "genome")) return dnaAlignment(cfgOpt);



	//////////////////////////////////////////////////////////////////////////////
	// Read the sequences
	//////////////////////////////////////////////////////////////////////////////

	typedef String<AminoAcid> TSequence;
	typedef String<char> TName;
	StringSet<TSequence, Owner<> > origStrSet;
	StringSet<TName> names;

	// Read the sequences from a T-Coffee library file or a regular sequence file in FASTA format
	if ((!length(value(cfgOpt, "seq"))) && (length(value(cfgOpt, "lib")))) {
		std::fstream strm;
		strm.open(toCString(value(cfgOpt, "lib")), std::ios_base::in | std::ios_base::binary);
		read(strm,origStrSet,names,TCoffeeLib());	
		strm.close();			
	} else if (length(value(cfgOpt, "seq"))) {
		_alignImportSequences(value(cfgOpt, "seq"), origStrSet, names);
	} else {
		std::cerr << valueHelp(cfgOpt) << std::endl;
		return -1;
	}
	typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
	TDepSequenceSet strSet(origStrSet);


	//////////////////////////////////////////////////////////////////////////////
	// Alignment of the sequences
	//////////////////////////////////////////////////////////////////////////////
	
	typedef Graph<Alignment<TDepSequenceSet, unsigned int> > TGraph;
	typedef Id<TGraph>::Type TId;
	Graph<Alignment<TDepSequenceSet, void, WithoutEdgeId> > gOut(strSet);
	
	// Case 1: Create a T-Coffee library
	if (value(cfgOpt, "output") == "tc_lib") {
		Blosum62 score_type_global(-1,-11);
		Blosum62 score_type_local(-2,-8);
		TGraph lib(strSet);
		if (value(cfgOpt, "method") == "global") generatePrimaryLibrary(lib, score_type_global, GlobalPairwise_Library() );
		else if (value(cfgOpt, "method") == "local") generatePrimaryLibrary(lib, score_type_local, LocalPairwise_Library() );
		else if (value(cfgOpt, "method") == "blast") {
			if (length(value(cfgOpt, "matches"))) {
				std::fstream strm_lib;
				strm_lib.open(toCString(value(cfgOpt, "matches")), std::ios_base::in | std::ios_base::binary);
				read(strm_lib, lib, names, BlastLib());	// Read library
				strm_lib.close();
			}
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
		}

		// Write the library
		std::fstream strm;
		strm.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::trunc);
		write(strm,lib,names,TCoffeeLib());
		strm.close();
	} 
	// Case 2: A library is send to us and we just align the library (SeqAn::M-Coffee)
	else if (length(value(cfgOpt, "lib"))) {
		TGraph g(strSet);
		std::fstream strm_lib;
		strm_lib.open(toCString(value(cfgOpt, "lib")), std::ios_base::in | std::ios_base::binary);
		read(strm_lib,g,TCoffeeLib());	// Read library
		strm_lib.close();

		//unsigned int nucs = 0;
		//for(unsigned int i = 0; i<length(strSet); ++i) nucs += length(strSet[i]);
		//std::cerr << "Avg. segment length: " << (double) nucs / (double) numVertices(g) << std::endl;

		// Read or build the guide tree
		Graph<Tree<double> > guideTree;
		if (length(value(cfgOpt, "usetree"))) {
			std::fstream strm_tree;
			strm_tree.open(toCString(value(cfgOpt, "usetree")), std::ios_base::in | std::ios_base::binary);
			read(strm_tree,guideTree,names,NewickFormat());	// Read newick tree
			strm_tree.close();
		} else {
			String<double> distanceMatrix;
			getDistanceMatrix(g, distanceMatrix, LibraryDistance());
			upgmaTree(distanceMatrix, guideTree);
		}

		//Progressive alignment
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
	// Case 3: Full protein alignment
	else {
		globalAlignment(strSet, names, value(cfgOpt, "matches"), gOut, MSA_Protein() );
		//testProfiles(strSet, names, value(cfgOpt, "matches"), gOut);
	}



	//////////////////////////////////////////////////////////////////////////////
	// Alignment output
	//////////////////////////////////////////////////////////////////////////////

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
