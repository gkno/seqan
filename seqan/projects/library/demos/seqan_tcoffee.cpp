#include <iostream>
#include <fstream>
#include <cstdio>
#include <seqan/graph.h>

using namespace seqan;


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
	else if (value(cfgOpt, "method") == "matching") globalMatching(strSet, names, cfgOpt);
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

template<typename TConfigOptions>
inline int
evaluateAlignment(TConfigOptions& cfgOpt) {
	typedef unsigned int TSize;
	typedef String<AminoAcid> TSequence;
	typedef String<char> TName;
	StringSet<TSequence, Owner<> > origStrSet;
	StringSet<TName> names;

	if (!length(value(cfgOpt, "infile"))) {
		std::cerr << valueHelp(cfgOpt) << std::endl;
		return -1;
	}

	// Read the sequences
	std::fstream strm;
	strm.open(toCString(value(cfgOpt, "infile")), std::ios_base::in | std::ios_base::binary);
	read(strm,origStrSet,names,FastaAlign());	
	strm.close();
	typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
	TDepSequenceSet strSet(origStrSet);
	
	//// Debug code
	//std::cout << "Sequences:" << std::endl;
	//for(TSize i = 0; i<length(origStrSet); ++i) {
	//	std::cout << '>' << names[i] << std::endl;
	//	std::cout << strSet[i] << std::endl;
	//}

	// Read the scoring matrix
	double gopening = 0; 
	double gextend = 0; 
	std::stringstream ssStream1(toCString(value(cfgOpt, "gop")));
	ssStream1 >> gopening; 
	gopening *= -1.0;
	std::stringstream ssStream2(toCString(value(cfgOpt, "gex")));
	ssStream2 >> gextend;
	gextend *= -1.0;
	typedef Score<double, ScoreMatrix<> > TScore;
	TScore scType;
	if (!length(value(cfgOpt, "matrix"))) {
		Blosum62 sc(-1 , -11);
		convertScoringMatrix(sc, scType, gextend, gopening);
	} else {
		loadScoreMatrix(scType, value(cfgOpt, "matrix"));
		scType.data_gap_extend = (double) (gextend);
		scType.data_gap_open = (double) (gopening);
	}

	// Read the alignment
	typedef Graph<Alignment<TDepSequenceSet, unsigned int> > TGraph;
	TGraph g(strSet);
	std::fstream strm_lib;
	strm_lib.open(toCString(value(cfgOpt, "infile")), std::ios_base::in | std::ios_base::binary);
	read(strm_lib,g,names,scType,FastaAlign());	
	strm_lib.close();

	//// Debug code
	//std::cout << g << std::endl;
	//std::cout << "Scoring parameters:" << std::endl;
	//std::cout << "Gap opening: " << gopening << std::endl;
	//std::cout << "Gap extension: " << gextend << std::endl;
	//typedef typename Value<TSequence>::Type TAlphabet;
	//TSize alphSize = ValueSize<TAlphabet>::VALUE;
	//for(TSize col = 0; col<alphSize; ++col) std::cout << TAlphabet(col) << ',';
	//std::cout << std::endl;
	//for(TSize row = 0; row<alphSize; ++row) {
	//	for(TSize col = 0; col<alphSize; ++col) {
	//		std::cout << score(scType, TAlphabet(row), TAlphabet(col)) << ',';
	//	}
	//	std::cout << std::endl;
	//}

	unsigned int numGapEx = 0;
	unsigned int numGap = 0;
	unsigned int numPairs = 0;
	unsigned int alignLen = 0;
	String<unsigned int> pairCount;
	double alignScore = alignmentEvaluation(g, scType, numGapEx, numGap, numPairs, pairCount, alignLen);
	std::cout << "Alignment Score: " << alignScore << std::endl;
	std::cout << "Alignment Length: " << alignLen << std::endl;
	std::cout << "#Match/Mismatch pairs: " << numPairs << std::endl;
	std::cout << "Score contribution by match/mismatch pairs: " << (alignScore - ((numGap * gopening) + (numGapEx * gextend))) << std::endl;
	std::cout << "#Gap extensions: " << numGapEx << std::endl;
	std::cout << "Score contribution by gap extensions: " << (numGapEx * gextend) << std::endl;
	std::cout << "#Gap openings: " << numGap << std::endl;
	std::cout << "Score contribution by gap openings: " << (numGap * gopening) << std::endl;
	typedef typename Value<TSequence>::Type TAlphabet;
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	std::cout << "#Pairs: " << std::endl;
	for(TSize col = 0; col<alphSize; ++col) std::cout << TAlphabet(col) << ',';
	std::cout << std::endl;
	for(TSize row = 0; row<alphSize; ++row) {
		for(TSize col = 0; col<alphSize; ++col) {
			std::cout << value(pairCount, row * alphSize + col) << ',';
		}
		std::cout << std::endl;
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
	TKey keys[] = {"seq","infile", "aln", "lib","matches","usetree","outfile","method","output","gop", "gex", "matrix", "score"};
	assignKeys(cfgOpt, keys, 13);
	// Set default options
	assign(cfgOpt, "output", "fasta");
	assign(cfgOpt, "outfile", "out.fasta");
	assign(cfgOpt, "method", "protein");
	assign(cfgOpt, "score", "false");
	assign(cfgOpt, "gop", "11");
	assign(cfgOpt, "gex", "1");
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
	if ((value(cfgOpt, "method") == "dna") || 
		(value(cfgOpt, "method") == "genome") ||
		(value(cfgOpt, "method") == "matching")) {
			return dnaAlignment(cfgOpt);
	}

	//////////////////////////////////////////////////////////////////////////////
	// Evaluation mode?
	//////////////////////////////////////////////////////////////////////////////

	if (value(cfgOpt, "score") == "true") {
		return evaluateAlignment(cfgOpt);
	}

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
		//write(strm,lib,names,BlastLib());
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
		if (length(value(cfgOpt, "matrix"))) testChristian(strSet, names, value(cfgOpt, "matrix"), value(cfgOpt, "gex"), value(cfgOpt, "gop"), gOut);
		else if (!length(value(cfgOpt, "aln"))) globalAlignment(strSet, names, value(cfgOpt, "matches"), gOut, MSA_Protein() );
		else metaGlobalAlignment(strSet, names, value(cfgOpt, "aln"), gOut, MSA_Protein() );
		//testChristian(strSet, names, value(cfgOpt, "matrix"), value(cfgOpt, "gex"), value(cfgOpt, "gop"), gOut);
		//testProfiles(strSet, names, value(cfgOpt, "matches"), gOut);
		//testIslands(strSet, names, value(cfgOpt, "matches"), gOut);
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
