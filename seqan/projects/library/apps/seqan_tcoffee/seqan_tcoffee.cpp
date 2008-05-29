#include <iostream>
#include <fstream>
#include <cstdio>
#include <seqan/graph_msa.h>

using namespace seqan;

template<typename TConfigOptions>
inline int
dnaAlignment(TConfigOptions& cfgOpt) {
	//////////////////////////////////////////////////////////////////////////////
	// Read the sequences
	//////////////////////////////////////////////////////////////////////////////
	
	typedef String<Dna> TSequence;
	typedef String<char> TName;
	StringSet<TSequence, Owner<> > origStrSet;
	StringSet<TName> names;
	_loadSequences(value(cfgOpt, "seq"), origStrSet, names);
	typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
	typedef typename Size<TDepSequenceSet>::Type TSize;
	TDepSequenceSet strSet(origStrSet);
	
	//////////////////////////////////////////////////////////////////////////////
	// Alignment
	//////////////////////////////////////////////////////////////////////////////

	typedef Graph<Alignment<TDepSequenceSet, TSize> > TGraph;
	typedef Id<TGraph>::Type TId;
	Graph<Alignment<TDepSequenceSet, void, WithoutEdgeId> > gOut(strSet);
	if (value(cfgOpt, "method") == "dna") globalAlignment(strSet, names, cfgOpt, gOut, MSA_Dna() );
	else globalAlignment(strSet, names, cfgOpt, gOut, MSA_Genome() );

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
	typedef typename Size<TConfigOptions>::Type TSize;
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

	// Make a dependent StringSet
	typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
	TDepSequenceSet strSet(origStrSet);
	
	// Read the scoring matrix
	typedef Score<double, ScoreMatrix<> > TScore;
	typedef typename Value<TScore>::Type TScoreValue;
	TScoreValue gopening = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gop"));
	TScoreValue gextend = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gex"));	
	TScore scType;
	if (!length(value(cfgOpt, "matrix"))) {
		Blosum62 sc(-1 , -11);
		convertScoringMatrix(sc, scType, gextend, gopening);
	} else {
		loadScoreMatrix(scType, value(cfgOpt, "matrix"));
		scType.data_gap_extend = (TScoreValue) (gextend);
		scType.data_gap_open = (TScoreValue) (gopening);
	}

	// Read the alignment
	typedef Graph<Alignment<TDepSequenceSet, TSize> > TGraph;
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

	TSize numGapEx = 0;
	TSize numGap = 0;
	TSize numPairs = 0;
	TSize alignLen = 0;
	String<TSize> pairCount;
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
	typedef Size<TKey>::Type TSize;
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
	append(helpMsg, "\nArguments\n===================\n");
	append(helpMsg, "-seq <FASTA Sequence File>\n");
	append(helpMsg, "\tA file with multiple sequences in FASTA format.\n");
	append(helpMsg, "-method [protein | dna | genome]\n");
	append(helpMsg, "\tSpecifies the type of sequence data, default is protein.\n");
	append(helpMsg, "-matches <File with Segment Matches>\n");
	append(helpMsg, "\tA file with gapless segment matches in BLAST Tabular Format (-m 8).\n");
	append(helpMsg, "-aln <FASTA Alignment File>,<FASTA Alignment File>,<FASTA Alignment File>\n");
	append(helpMsg, "\tTriggers a meta-alignment of all given alignment files.\n");
	append(helpMsg, "-infile <FASTA Alignment File>\n");
	append(helpMsg, "\tSwitch to evaluate a given multiple alignment file in FASTA format.\n");
	append(helpMsg, "-outfile <Alignment Filename>\n");
	append(helpMsg, "\tSpecifies the name of the output file, default is out.fasta.\n");
	append(helpMsg, "-output [fasta | msf]\n");
	append(helpMsg, "\tSpecifies the output format, default is fasta.\n");
	append(helpMsg, "-gop <Number>\n");
	append(helpMsg, "\tSpecifies the gap open penalty, default is 11.\n");
	append(helpMsg, "-gex <Number>\n");
	append(helpMsg, "\tSpecifies the gap extension penalty, default is 1.\n");
	append(helpMsg, "-matrix <Matrix File>\n");
	append(helpMsg, "\tSpecifies a matrix file, default for protein is Blosum62.\n");
	append(helpMsg, "-usetree <Newick Guide Tree>\n");
	append(helpMsg, "\tA guide tree in newick format.\n");
	append(helpMsg, "-lib <T-Coffee library>\n");
	append(helpMsg, "\tA T-Coffee library that is used to align the sequences.\n");

	append(helpMsg, "\n\n\nCommon Tasks\n===================\n");
	append(helpMsg, "\nProtein Alignment:\n");
	append(helpMsg, "seqan_tcoffee -seq seq.fasta\n");
	append(helpMsg, "seqan_tcoffee -seq seq.fasta -usetree guideTree.dnd -output fasta -outfile align.fasta\n");
	append(helpMsg, "seqan_tcoffee -seq seq.fasta -usetree guideTree.dnd -matrix BLOSUM45 -gop 11 -gex 1 -output msf -outfile align.msf\n");
	append(helpMsg, "\nDna / Genome Alignment:\n");
	append(helpMsg, "seqan_tcoffee -seq seq.fasta -method dna\n");
	append(helpMsg, "seqan_tcoffee -seq seq.fasta -method genome -output fasta -outfile genomeAlign.fasta\n");
	append(helpMsg, "seqan_tcoffee -seq seq.fasta -matches blastHit.tab -method genome\n");
	append(helpMsg, "\nMeta-Alignment:\n");
	append(helpMsg, "seqan_tcoffee -seq seq.fasta -aln sub1.fasta,sub2.fasta\n");
	append(helpMsg, "seqan_tcoffee -seq seq.fasta -usetree guideTree.dnd -aln sub1.fasta,sub2.fasta -output msf -outfile final.msf\n");
	append(helpMsg, "\nAlignment Evaluation:\n");
	append(helpMsg, "seqan_tcoffee -infile align.fasta\n");
	append(helpMsg, "seqan_tcoffee -infile align.fasta -gop 2 -gex 1 -matrix BLOSUM45\n");
	append(helpMsg, "\nAlignment from a T-Coffee Library:\n");
	append(helpMsg, "seqan_tcoffee -lib tcoffee.tc_lib\n");
	append(helpMsg, "seqan_tcoffee -lib tcoffee.tc_lib -usetree tcTree.dnd -output msf -outfile align.msf\n");
	assignHelp(cfgOpt, helpMsg);
	
	if (argc < 2) {	std::cerr << valueHelp(cfgOpt) << std::endl; return -1; }
	if (!parseCmdLine(argc, argv, cfgOpt)) return -1;

	//////////////////////////////////////////////////////////////////////////////
	// Alignment of Dna Sequences
	//////////////////////////////////////////////////////////////////////////////

	if ((value(cfgOpt, "method") == "dna") || (value(cfgOpt, "method") == "genome")) {
		return dnaAlignment(cfgOpt);
	}


	//////////////////////////////////////////////////////////////////////////////
	// Evaluation mode?
	//////////////////////////////////////////////////////////////////////////////

	if (!empty(value(cfgOpt, "infile"))) assign(cfgOpt, "score", "true");
	if (value(cfgOpt, "score") == "true") return evaluateAlignment(cfgOpt);

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
		_loadSequences(value(cfgOpt, "seq"), origStrSet, names);
	} else {
		std::cerr << valueHelp(cfgOpt) << std::endl;
		return -1;
	}
	typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
	TDepSequenceSet strSet(origStrSet);

	//////////////////////////////////////////////////////////////////////////////
	// Alignment of the sequences
	//////////////////////////////////////////////////////////////////////////////
	
	typedef Graph<Alignment<TDepSequenceSet, TSize> > TGraph;
	typedef Id<TGraph>::Type TId;
	Graph<Alignment<TDepSequenceSet, void, WithoutEdgeId> > gOut(strSet);
	
	// Case 1: Interaction with T-Coffee
	// A T-Coffee library is send to us and we just align according to the library (SeqAn::M-Coffee)
	if (length(value(cfgOpt, "lib"))) {
		TGraph g(strSet);
		std::fstream strm_lib;
		strm_lib.open(toCString(value(cfgOpt, "lib")), std::ios_base::in | std::ios_base::binary);
		read(strm_lib,g,TCoffeeLib());	// Read library
		strm_lib.close();

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
		TSize nSeq = length(strSet);
		TSize threshold = 30;
		if (nSeq < threshold) {
			// Full triplet...
			tripletLibraryExtension(g);

			progressiveAlignment(g, guideTree, gOut);
		} else {
			// Triplet only on groups of sequences
			progressiveAlignment(g, guideTree, gOut, threshold);
		}
	} 
	// Case 2: Alignment via SeqAn
	else {
		// Normal Alignment
		if ((!length(value(cfgOpt, "aln"))) && (!length(value(cfgOpt, "matrix")))) globalAlignment(strSet, names, cfgOpt, gOut, MSA_Protein() );
		// Meta-Alignment
		else if (length(value(cfgOpt, "aln"))) metaGlobalAlignment(strSet, names, cfgOpt, gOut, MSA_Protein() );
		// Custom Alignment
		else globalAlignment(strSet, names, cfgOpt, gOut, MSA_ProteinCustom() );
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
