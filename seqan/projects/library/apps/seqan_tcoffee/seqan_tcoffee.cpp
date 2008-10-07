#define SEQAN_PROFILE

#include <seqan/graph_msa.h>

#include <iostream>
#include <fstream>
#include <cstdio>

using namespace seqan;


// Profiling
SEQAN_PROTIMESTART(__myProfileTime);


//////////////////////////////////////////////////////////////////////////////////

template<typename TConfigOptions, typename TScore, typename TAlphabet>
inline int
evaluateAlignment(TConfigOptions const& cfgOpt, TScore const& scType, TAlphabet) {
	typedef typename Value<TScore>::Type TScoreValue;
	typedef String<TAlphabet> TSequence;
	typedef typename Size<TSequence>::Type TSize;
	typedef String<char> TName;
	StringSet<TSequence, Owner<> > origStrSet;
	StringSet<TName> names;

	// Read the sequences
	std::fstream strm;
	strm.open(toCString(value(cfgOpt, "infile")), std::ios_base::in | std::ios_base::binary);
	read(strm,origStrSet,names,FastaAlign());	
	strm.close();

	// Make a dependent StringSet
	typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
	TDepSequenceSet strSet(origStrSet);
	
	// Read the alignment
	typedef String<Fragment<> > TFragmentString;
	String<TScoreValue> scores;
	TFragmentString matches;
	std::fstream strm_lib;
	strm_lib.open(toCString(value(cfgOpt, "infile")), std::ios_base::in | std::ios_base::binary);
	read(strm_lib,matches, scores, names, FastaAlign());	
	strm_lib.close();

	// Build the alignment graph
	typedef Graph<Alignment<TDepSequenceSet, TSize> > TGraph;
	TGraph g(strSet);
	buildAlignmentGraph(matches, g, FrequencyCounting() );

	// Print the scoring information
	TScoreValue gop = scoreGapOpen(scType);
	TScoreValue gex = scoreGap(scType);
	std::cout << "Scoring parameters:" << std::endl;
	std::cout << "*Gap opening: " << gop << std::endl;
	std::cout << "*Gap extension: " << gex << std::endl;
	std::cout << "*Scoring matrix: " << std::endl;
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	std::cout << "   ";
	for(TSize col = 0; col<alphSize; ++col) std::cout << TAlphabet(col) << ',';
	std::cout << std::endl;
	for(TSize row = 0; row<alphSize; ++row) {
		for(TSize col = 0; col<alphSize; ++col) {
			if (col == 0) std::cout << TAlphabet(row) << ": ";
			std::cout << score(scType, TAlphabet(row), TAlphabet(col));
			if (col < alphSize - 1) std::cout << ',';
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	// Print the alignment information
	TSize numGapEx = 0;
	TSize numGap = 0;
	TSize numPairs = 0;
	TSize alignLen = 0;
	String<TSize> pairCount;
	TScoreValue alignScore = alignmentEvaluation(g, scType, numGapEx, numGap, numPairs, pairCount, alignLen);
	std::cout << "Alignment Score: " << alignScore << std::endl;
	std::cout << "Alignment Length: " << alignLen << std::endl;
	std::cout << "#Match-Mismatch pairs: " << numPairs << std::endl;
	std::cout << "Score contribution by match-mismatch pairs: " << (alignScore - (((TScoreValue) numGap * gop) + ((TScoreValue) numGapEx * gex))) << std::endl;
	std::cout << "#Gap extensions: " << numGapEx << std::endl;
	std::cout << "Score contribution by gap extensions: " << ((TScoreValue) numGapEx * gex) << std::endl;
	std::cout << "#Gap openings: " << numGap << std::endl;
	std::cout << "Score contribution by gap openings: " << ((TScoreValue) numGap * gop) << std::endl;
	std::cout << std::endl;
	std::cout << "#Pairs: " << std::endl;
	std::cout << "   ";
	for(TSize col = 0; col<alphSize; ++col) std::cout << TAlphabet(col) << ',';
	std::cout << std::endl;
	for(TSize row = 0; row<alphSize; ++row) {
		for(TSize col = 0; col<alphSize; ++col) {
			if (col == 0) std::cout << TAlphabet(row) << ": ";
			std::cout << value(pairCount, row * alphSize + col);
			if (col < alphSize - 1) std::cout << ',';
		}
		std::cout << std::endl;
	}

	return 0;

}

//////////////////////////////////////////////////////////////////////////////////

template<typename TConfigOptions>
inline int
evaluateAlignment(TConfigOptions const& cfgOpt, AminoAcid) 
{
	if (!length(value(cfgOpt, "matrix"))) {
		typedef typename Value<Blosum62>::Type TScoreValue;
		Blosum62 sc(-1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gex")) , -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gop")));
		return evaluateAlignment(cfgOpt, sc, AminoAcid() );
	} else {
		typedef Score<int, ScoreMatrix<> > TScore;
		typedef typename Value<TScore>::Type TScoreValue;
		TScore sc;
		loadScoreMatrix(sc, value(cfgOpt, "matrix"));
		sc.data_gap_extend = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gex")); 
		sc.data_gap_open = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gop"));
		return evaluateAlignment(cfgOpt, sc, AminoAcid() );
	}
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TConfigOptions>
inline int
evaluateAlignment(TConfigOptions const& cfgOpt, Dna) 
{
	typedef Score<int> TScore;
	typedef typename Value<TScore>::Type TScoreValue;
	TScore scType(_stringToNumber<TScoreValue>(value(cfgOpt, "msc")),_stringToNumber<TScoreValue>(value(cfgOpt, "mmsc")),-1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gex")),-1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gop")));
	return evaluateAlignment(cfgOpt, scType, Dna() );
}



//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TConfigOptions, typename TAlignmentGraph>
inline void
predefinedGlobalAlignment(StringSet<TString, TSpec> const& seqSet,
						  StringSet<TName, TSpec2> const& nameSet,
						  TConfigOptions& cfgOpt,
						  TAlignmentGraph& gOut)
{
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef typename Size<TStringSet>::Type TSize;
	typedef Graph<Alignment<TStringSet, TSize> > TGraph;
	typedef typename Id<TGraph>::Type TId;
	TSize nSeq = length(seqSet);
	TSize threshold = 30;

	// Select all pairs
	String<Pair<TId, TId> > pList;
	selectPairs(seqSet, pList);

	// Set-up a distance matrix
	typedef String<double> TDistanceMatrix;
	typedef typename Value<TDistanceMatrix>::Type TDistanceValue;
	TDistanceMatrix distanceMatrix;

	// Set-up alignment scoring matrices
	typedef Blosum62 TScore;
	typedef typename Value<TScore>::Type TScoreValue;
	TScore score_type_global(-1,-11);
	TScore score_type_local(-2,-8);

	// Containers for segment matches and corresponding scores 
	typedef String<Fragment<> > TFragmentString;
	TFragmentString matches;
	typedef String<TScoreValue> TScoreValues;
	TScoreValues scores;

	// Compute segment matches from global pairwise alignments
	appendSegmentMatches(seqSet, pList, score_type_global, matches, scores, distanceMatrix, GlobalPairwise_Library() );
	
	// Append segment matches from local pairwise alignments
	appendSegmentMatches(seqSet, pList, score_type_local, matches, scores, LocalPairwise_Library() );

	// Select a subset of pairs for end-gaps free
	typedef std::multimap<TDistanceValue, Pair<TId, TId> > TBestPairs;
	TBestPairs bestPairs;
	TDistanceValue dist = 0;
	for(TSize i=0;i<nSeq-1;++i) {
		for(TSize j=i+1;j<nSeq;++j) {
			TDistanceValue d = value(distanceMatrix, i*nSeq+j);
			bestPairs.insert(std::make_pair(d, Pair<TId, TId>(positionToId(seqSet, i),positionToId(seqSet, j)))); 
			dist+=d;
		}
	}
	TSize numPairs = nSeq * (nSeq - 1) / 2;
	dist /= numPairs;
	String<Pair<TId, TId> > pListLocal;
	pListLocal = pList;
	if (nSeq > threshold) {
		clear(pListLocal);
		typename TBestPairs::reverse_iterator itBestR = bestPairs.rbegin();
		typename TBestPairs::iterator itBestF = bestPairs.begin();
		TSize limit = 4 * nSeq;
		for(TSize i=0;i<limit;++i, ++itBestR, ++itBestF) {
			appendValue(pListLocal, itBestF->second);
			appendValue(pListLocal, itBestR->second);
		}
	}

	// Append segment matches from end-gaps free alignments
	Nothing noth;
	if (dist > 0.75) {
		Blosum30 sT(-4,-20);
		AlignConfig<true,true,true,true> ac;
		appendSegmentMatches(seqSet, pListLocal, sT, matches, scores, noth, ac, GlobalPairwise_Library() );
	} else if (dist < 0.50) {
		Blosum80 sT(-2, -12);
		AlignConfig<true,true,true,true> ac;
		appendSegmentMatches(seqSet, pListLocal, sT, matches, scores, noth, ac, GlobalPairwise_Library() );
	} else {
		Blosum62 sT(-3,-14);
		AlignConfig<true,true,true,true> ac;
		appendSegmentMatches(seqSet, pListLocal, sT, matches, scores, noth, ac, GlobalPairwise_Library() );
	}

	// Include external segment matches
	if (length(value(cfgOpt, "blast"))) {
		std::fstream strm_lib;
		strm_lib.open(toCString(value(cfgOpt, "blast")), std::ios_base::in | std::ios_base::binary);
		read(strm_lib, matches, scores, nameSet, BlastLib());
		strm_lib.close();
	}

	// Score the matches
	scoreMatches(seqSet, score_type_global, matches, scores);
	
	// Use these segment matches for the initial alignment graph
	TGraph g(seqSet);
	buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
	
	// Guide tree
	Graph<Tree<double> > guideTree;
	if (length(value(cfgOpt, "usetree"))) {
		std::fstream strm_tree;
		strm_tree.open(toCString(value(cfgOpt, "usetree")), std::ios_base::in | std::ios_base::binary);
		read(strm_tree,guideTree,nameSet,NewickFormat());	// Read newick tree
		strm_tree.close();
	} else {
		//upgmaTree(distanceMatrix, guideTree);
		slowNjTree(distanceMatrix, guideTree);
	}

	// Triplet extension and progressive alignment
	if (nSeq < threshold) {
		// Full triplet...
		tripletLibraryExtension(g);

		// ... and normal progressive alignment with guide tree
		progressiveAlignment(g, guideTree, gOut);
	} else {
		// Triplet only on groups of sequences
		progressiveAlignment(g, guideTree, gOut, threshold);
	}

	clear(guideTree);
	clear(distanceMatrix);
	clear(g);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TScoreMatrix, typename TConfigOptions, typename TAlignmentGraph>
inline void
globalAlignment(StringSet<TString, TSpec> const& seqSet,
				StringSet<TName, TSpec2> const& nameSet,
				TConfigOptions const& cfgOpt,
				TScoreMatrix const& scType,
				TAlignmentGraph& gOut)
{
	typedef typename Value<TScoreMatrix>::Type TScoreValue;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef typename Size<TStringSet>::Type TSize;
	typedef Graph<Alignment<TStringSet, TSize> > TGraph;
	typedef typename Id<TGraph>::Type TId;
	TSize nSeq = length(seqSet);
	TSize threshold = 30;

	// Select all pairs
	String<Pair<TId, TId> > pList;
	selectPairs(seqSet, pList);

	// Set-up a distance matrix
	typedef String<double> TDistanceMatrix;
	typedef typename Value<TDistanceMatrix>::Type TDistanceValue;
	TDistanceMatrix distanceMatrix;

	// Containers for segment matches and corresponding scores 
	typedef String<Fragment<> > TFragmentString;
	TFragmentString matches;
	typedef String<TScoreValue> TScoreValues;
	TScoreValues scores;

	// Include segment matches from other alignments
	if (length(value(cfgOpt, "aln"))) {
		std::cout << "Parsing external alignment files:" << std::endl;
		String<char> alnFiles = value(cfgOpt, "aln");
		String<char> alignmentFile;
		for(TSize i = 0; i<=length(alnFiles); ++i) {
			if ((i == length(alnFiles) || (value(alnFiles, i) == ','))) {		
				std::cout << "*Alignment file: " << alignmentFile << std::endl;
				std::stringstream input;
				input << alignmentFile;
				std::ifstream strm_lib;
				strm_lib.open(input.str().c_str(), std::ios_base::in | std::ios_base::binary);
				read(strm_lib, matches, scores, nameSet, FastaAlign());
				strm_lib.close();
				clear(alignmentFile);
			} else {
				if ((value(alnFiles, i) != ' ') && (value(alnFiles, i) != '\t')) appendValue(alignmentFile, value(alnFiles, i));
			}
		}	
		std::cout << "External segment matches done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
	}

	// Include computed segment matches
	if (length(value(cfgOpt, "method"))) {
		std::cout << "Computing segment matches:" << std::endl;
		String<char> methodNames = value(cfgOpt, "method");
		String<char> currentMethod;
		for(TSize i = 0; i<=length(methodNames); ++i) {
			if ((i == length(methodNames) || (value(methodNames, i) == ','))) {		
				std::cout << "*Method: " << currentMethod << std::endl;
				if (currentMethod == "global") {
					// Compute segment matches from global pairwise alignments
					appendSegmentMatches(seqSet, pList, scType, matches, scores, distanceMatrix, GlobalPairwise_Library() );
				} else if (currentMethod == "local") {
					// Compute segment matches from local pairwise alignments
					appendSegmentMatches(seqSet, pList, scType, matches, scores, LocalPairwise_Library() );
				} else if (currentMethod == "overlap") {
					// Compute segment matches from overlap alignments
					Nothing noth;
					AlignConfig<true,true,false,false> ac;
					appendSegmentMatches(seqSet, pList, scType, matches, scores, noth, ac, GlobalPairwise_Library() );
				} else if (currentMethod == "lcs") {
					// Compute segment matches from the longest common subsequence
					appendSegmentMatches(seqSet, matches, scores, Lcs_Library() );
				} else {
					std::cout << "*Unknown method!!!" << std::endl;
				}
				std::cout << "*Done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
				clear(currentMethod);
			} else {
				if ((value(methodNames, i) != ' ') && (value(methodNames, i) != '\t')) appendValue(currentMethod, value(methodNames, i));
			}
		}	
	}

	// Include a T-Coffee library
	if (length(value(cfgOpt, "lib"))) {
		std::cout << "Parsing a T-Coffee Library:" << std::endl;
		String<char> libNames = value(cfgOpt, "lib");
		String<char> currentLib;
		for(TSize i = 0; i<=length(libNames); ++i) {
			if ((i == length(libNames) || (value(libNames, i) == ','))) {
				std::cout << "*T-Coffee library: " << currentLib << std::endl;
				std::stringstream input;
				input << currentLib;
				std::ifstream strm_lib;
				strm_lib.open(input.str().c_str(), std::ios_base::in | std::ios_base::binary);
				read(strm_lib, matches, scores, nameSet, TCoffeeLib());
				strm_lib.close();
				clear(currentLib);
			} else {
				if ((value(libNames, i) != ' ') && (value(libNames, i) != '\t')) appendValue(currentLib, value(libNames, i));
			}
		}
		std::cout << "Parsing done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
	}
	
	// Include MUMmer segment matches
	if (length(value(cfgOpt, "mummer"))) {
		std::cout << "Parsing MUMmer segment matches:" << std::endl;
		String<char> mummerFiles = value(cfgOpt, "mummer");
		String<char> currentMumFile;
		for(TSize i = 0; i<=length(mummerFiles); ++i) {
			if ((i == length(mummerFiles) || (value(mummerFiles, i) == ','))) {		
				std::cout << "*MUMmer file: " << currentMumFile << std::endl;
				std::stringstream input;
				input << currentMumFile;
				std::ifstream strm_lib;
				strm_lib.open(input.str().c_str(), std::ios_base::in | std::ios_base::binary);
				read(strm_lib, matches, scores, seqSet, nameSet, MummerLib());		
				strm_lib.close();
				clear(currentMumFile);
			} else {
				if ((value(mummerFiles, i) != ' ') && (value(mummerFiles, i) != '\t')) appendValue(currentMumFile, value(mummerFiles, i));
			}
		}	
		std::cout << "Parsing done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
	}

	// Include BLAST segment matches
	if (length(value(cfgOpt, "blast"))) {
		std::cout << "Parsing BLAST segment matches:" << std::endl;
		String<char> blastFiles = value(cfgOpt, "blast");
		String<char> currentBlastFile;
		for(TSize i = 0; i<=length(blastFiles); ++i) {
			if ((i == length(blastFiles) || (value(blastFiles, i) == ','))) {		
				std::cout << "*BLAST file: " << currentBlastFile << std::endl;
				std::stringstream input;
				input << currentBlastFile;
				std::ifstream strm_lib;
				strm_lib.open(input.str().c_str(), std::ios_base::in | std::ios_base::binary);
				read(strm_lib, matches, scores, nameSet, BlastLib());
				strm_lib.close();
				clear(currentBlastFile);
			} else {
				if ((value(blastFiles, i) != ' ') && (value(blastFiles, i) != '\t')) appendValue(currentBlastFile, value(blastFiles, i));
			}
		}
		std::cout << "Parsing done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
	}

	std::cout << "Total number of segment matches: " << length(matches) << std::endl;
	// Score the matches
	if ((value(cfgOpt, "rescore") == "true") && (length(value(cfgOpt, "seq")))) {
		std::cout << "Scoring method: Re-Score" << std::endl;
		scoreMatches(seqSet, scType, matches, scores);
		std::cout << "Scoring done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;\
	}
	
	// Use these segment matches for the initial alignment graph
	TGraph g(seqSet);
	std::cout << "Construction of alignment graph: FractionalScore" << std::endl;
	buildAlignmentGraph(matches, scores, g, FractionalScore() );
	clear(matches);
	clear(scores);
	std::cout << "Alignment graph construction done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;

	// Guide tree
	Graph<Tree<double> > guideTree;
	if (length(value(cfgOpt, "usetree"))) {
		std::cout << "Guide Tree: " << value(cfgOpt, "usetree") << std::endl;
		std::fstream strm_tree;
		strm_tree.open(toCString(value(cfgOpt, "usetree")), std::ios_base::in | std::ios_base::binary);
		read(strm_tree,guideTree,nameSet,NewickFormat());	// Read newick tree
		strm_tree.close();
	} else {
		std::cout << "Guide Tree: Neighbor Joining" << std::endl;
		// Check if we have a valid distance matrix
		if (empty(distanceMatrix)) getDistanceMatrix(g, distanceMatrix, LibraryDistance());
		slowNjTree(distanceMatrix, guideTree);
	}
	std::cout << "Guide tree done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;

	// Triplet extension and progressive alignment
	if (nSeq < threshold) {
		std::cout << "Progressive Alignment: Default" << std::endl;
		// Full triplet...
		tripletLibraryExtension(g);

		// ... and normal progressive alignment with guide tree
		progressiveAlignment(g, guideTree, gOut);
	} else {
		std::cout << "Progressive Alignment: Group-based" << std::endl;
		// Triplet only on groups of sequences
		progressiveAlignment(g, guideTree, gOut, threshold);
	}
	std::cout << "Progressive alignment done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
	clear(guideTree);
	clear(distanceMatrix);
	clear(g);
	std::cout << "Clean-up done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TConfigOptions, typename TAlignmentGraph>
inline void
customizedMsaAlignment(StringSet<TString, TSpec> const& seqSet,
					   StringSet<TName, TSpec2> const& nameSet,
					   TConfigOptions& cfgOpt,
					   TAlignmentGraph& gOut,
					   AminoAcid)
{
	std::cout << "Scoring: " << std::endl;
	if (!length(value(cfgOpt, "matrix"))) {
		typedef typename Value<Blosum62>::Type TScoreValue;
		TScoreValue gopening = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gop")); 
		TScoreValue gextend = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gex")); 
		std::cout << "*Gap open penalty: " << gopening << std::endl;
		std::cout << "*Gap extension penalty: " << gextend << std::endl;
		std::cout << "*Scoring Matrix: Blosum62" << std::endl;
		Blosum62 sc(gextend , gopening);
		globalAlignment(seqSet, nameSet, cfgOpt, sc, gOut);
	} else {
		typedef Score<int, ScoreMatrix<> > TScore;
		typedef typename Value<TScore>::Type TScoreValue;
		TScore sc;
		TScoreValue gopening = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gop")); 
		TScoreValue gextend = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gex")); 
		std::cout << "*Gap open penalty: " << gopening << std::endl;
		std::cout << "*Gap extension penalty: " << gextend << std::endl;
		std::cout << "*Scoring Matrix: " << value(cfgOpt, "matrix") << std::endl;
		loadScoreMatrix(sc, value(cfgOpt, "matrix"));
		sc.data_gap_extend = (TScoreValue) (gextend);
		sc.data_gap_open = (TScoreValue) (gopening);
		globalAlignment(seqSet, nameSet, cfgOpt, sc, gOut);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TConfigOptions, typename TAlignmentGraph>
inline void
customizedMsaAlignment(StringSet<TString, TSpec> const& seqSet,
					   StringSet<TName, TSpec2> const& nameSet,
					   TConfigOptions& cfgOpt,
					   TAlignmentGraph& gOut,
					   Dna)
{
	typedef Score<int> TScore;
	typedef typename Value<TScore>::Type TScoreValue;
	TScoreValue gop = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gop")); 
	TScoreValue gex = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gex")); 
	TScoreValue match = _stringToNumber<TScoreValue>(value(cfgOpt, "msc")); 
	TScoreValue mismatch = _stringToNumber<TScoreValue>(value(cfgOpt, "mmsc")); 
	std::cout << "Scoring: " << std::endl;
	std::cout << "*Gap open penalty: " << gop << std::endl;
	std::cout << "*Gap extension penalty: " << gex << std::endl;
	std::cout << "*Match score: " << match << std::endl;
	std::cout << "*Mismatch score: " << mismatch << std::endl;
	TScore scType(match,mismatch,gex,gop);
	globalAlignment(seqSet, nameSet, cfgOpt, scType, gOut);
}



//////////////////////////////////////////////////////////////////////////////////

template<typename TConfigOptions, typename TAlphabet>
inline int
customizedMsaAlignment(TConfigOptions const& cfgOpt, TAlphabet) {

	//////////////////////////////////////////////////////////////////////////////
	// Read the sequences
	//////////////////////////////////////////////////////////////////////////////
	typedef String<TAlphabet> TSequence;
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
	typedef typename Size<TDepSequenceSet>::Type TSize;
	TDepSequenceSet strSet(origStrSet);
	std::cout << "Number of sequences: " << length(strSet) << std::endl;
	std::cout << "Import of sequences done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
	
	//////////////////////////////////////////////////////////////////////////////
	// Alignment of the sequences
	//////////////////////////////////////////////////////////////////////////////
	typedef Graph<Alignment<TDepSequenceSet, TSize> > TGraph;
	Graph<Alignment<TDepSequenceSet, void, WithoutEdgeId> > gOut(strSet);
	
	// Predefined alignment pipeline
	//predefinedGlobalAlignment(strSet, names, cfgOpt, gOut);
	// Custom Alignment
	customizedMsaAlignment(strSet, names, cfgOpt, gOut, TAlphabet() );

	//////////////////////////////////////////////////////////////////////////////
	// Alignment output
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
	std::cout << "Output done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;

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
	TKey keys[] = {"seq", "alphabet", "aln", "lib","blast","mummer", "usetree","outfile","method","output","gop", "gex", "matrix", "infile", "msc", "mmsc", "rescore"};
	assignKeys(cfgOpt, keys, 17);

	// Set default options
	assign(cfgOpt, "rescore", "true");
	assign(cfgOpt, "output", "fasta");
	assign(cfgOpt, "outfile", "out.fasta");
	assign(cfgOpt, "alphabet", "protein");
	assign(cfgOpt, "gop", "11");
	assign(cfgOpt, "gex", "1");
	assign(cfgOpt, "msc", "5");
	assign(cfgOpt, "mmsc", "-4");

	// Help Message
	String<char> helpMsg;
	append(helpMsg, "Usage: seqan_tcoffee -seq <FASTA Sequence File> [Options]\n");
	append(helpMsg, "\nOptions\n");
	// Main options
	append(helpMsg, "\nMain Options\n------------\n");
	append(helpMsg, "-seq <FASTA Sequence File>\n");
	append(helpMsg, "\tFile with multiple sequences in FASTA format.\n\n");
	append(helpMsg, "-alphabet [protein | dna ]\n");
	append(helpMsg, "\tSequence alphabet, default is protein.\n\n");
	append(helpMsg, "-outfile <Alignment Filename>\n");
	append(helpMsg, "\tName of the output file, default is out.fasta.\n\n");
	append(helpMsg, "-output [fasta | msf]\n");
	append(helpMsg, "\tOutput format, default is fasta.\n\n");
	// Segment match generation options
	append(helpMsg, "\nSegment-Match Generation Options\n------------\n");
	append(helpMsg, "-method [global], [local], [overlap], [lcs]\n");
	append(helpMsg, "\tMethods to generate segment matches, default is global, local.\n\n");
	append(helpMsg, "-blast <BLAST matches>, <BLAST matches>, ...\n");
	append(helpMsg, "\tFiles with gapless segment matches in BLAST tabular format (-m 8 -g F).\n\n");
	append(helpMsg, "-mummer <MUMmer matches>, <MUMmer matches>, ...\n");
	append(helpMsg, "\tFiles with gapless segment matches in MUMmer format.\n\n");
	append(helpMsg, "-aln <FASTA Alignment File>, <FASTA Alignment File>, ...\n");
	append(helpMsg, "\tSegment-matches are extracted from given alignments (Meta-Alignment).\n\n");
	append(helpMsg, "-lib <T-Coffee library>, <T-Coffee library>, ...\n");
	append(helpMsg, "\tSegment-matches are extracted from T-Coffee libraries.\n\n");
	// Scoring
	append(helpMsg, "\nScoring Options\n------------\n");
	append(helpMsg, "-gop <Number>\n");
	append(helpMsg, "\tGap open penalty, default is 11.\n\n");
	append(helpMsg, "-gex <Number>\n");
	append(helpMsg, "\tGap extension penalty, default is 1.\n\n");
	append(helpMsg, "-matrix <Score-Matrix File>\n");
	append(helpMsg, "\tSpecifies a score-matrix file, default is Blosum62.\n\n");
	append(helpMsg, "-msc <Number>\n");
	append(helpMsg, "\tMatch score for Dna alphabet, default is 5.\n\n");
	append(helpMsg, "-mmsc <Number>\n");
	append(helpMsg, "\tMismatch penalty for Dna alphabet, default is 4.\n\n");
	append(helpMsg, "-rescore [ true | false ]\n");
	append(helpMsg, "\tRe-score all segment-matches, default is true.\n\n");
	// Guide Tree
	append(helpMsg, "Guide Tree Options\n------------\n");
	append(helpMsg, "-usetree <Newick Guide Tree>\n");
	append(helpMsg, "\tA guide tree in newick format.\n\n");
	// Alignment evaluation
	append(helpMsg, "Alignment Evaluation Options\n------------\n");
	append(helpMsg, "-infile <FASTA Alignment File>\n");
	append(helpMsg, "\tEvaluate the given multiple alignment.\n\n");


	append(helpMsg, "\n\n\nExamples\n");
	append(helpMsg, "\nProtein Alignment:\n");
	append(helpMsg, "\t./seqan_tcoffee -seq seq.fasta\n");
	append(helpMsg, "\t./seqan_tcoffee -seq seq.fasta -method global, local\n");
	append(helpMsg, "\nDna Alignment:\n");
	append(helpMsg, "\t./seqan_tcoffee -seq seq.fasta -alphabet dna\n");
	append(helpMsg, "\nGenome Alignment:\n");
	append(helpMsg, "\t./seqan_tcoffee -seq seq.fasta -method lcs -alphabet dna\n");
	append(helpMsg, "\t./seqan_tcoffee -seq seq.fasta -mummer my.mums -blast my1, my2 -alphabet dna\n");
	append(helpMsg, "\nMeta-Alignment:\n");
	append(helpMsg, "\t./seqan_tcoffee -seq seq.fasta -aln sub1.fasta, sub2.fasta\n");
	append(helpMsg, "\nAlignment Evaluation:\n");
	append(helpMsg, "\t./seqan_tcoffee -infile align.fasta\n");
	append(helpMsg, "\t./seqan_tcoffee -infile align.fasta -gop 5 -gex 4 -alphabet dna\n");
	append(helpMsg, "\nAlignment from a T-Coffee Library:\n");
	append(helpMsg, "\t./seqan_tcoffee -lib tcoffee.tc_lib\n");
	assignHelp(cfgOpt, helpMsg);

	std::cout << "*********************************************" << std::endl;
	std::cout << "* Segment-based multiple sequence alignment *" << std::endl;
	std::cout << "*                                           *" << std::endl;
	std::cout << "* SeqAn::T-Coffee                           *" << std::endl;
	std::cout << "* Version: 1.1 (10. September 2008)         *" << std::endl;
	std::cout << "*********************************************" << std::endl;
	std::cout << std::endl;
	
	if (argc < 2) {	std::cerr << valueHelp(cfgOpt) << std::endl; return -1; }
	if (!parseCmdLine(argc, argv, cfgOpt)) return -1;

	// Check if any segment-match generation procedure is selected, otherwise set the default
	if ((!length(value(cfgOpt, "blast"))) && (!length(value(cfgOpt, "mummer"))) && (!length(value(cfgOpt, "method"))) && (!length(value(cfgOpt, "aln"))) && (!length(value(cfgOpt, "lib")))) {
		assign(cfgOpt, "method", "global, local");
	}

	//////////////////////////////////////////////////////////////////////////////
	// Evaluation mode?
	//////////////////////////////////////////////////////////////////////////////

	if (length(value(cfgOpt, "infile"))) {
		if ((value(cfgOpt, "alphabet") == "dna")) {
			std::cout << "Alignment evaluation for Dna sequences" << std::endl;
			return evaluateAlignment(cfgOpt, Dna() );
		} else {
			std::cout << "Alignment evaluation for AminoAcid sequences" << std::endl;
			return evaluateAlignment(cfgOpt, AminoAcid() );
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	// Alignment of Dna or Amino Acid Sequences
	//////////////////////////////////////////////////////////////////////////////

	if ((value(cfgOpt, "alphabet") == "dna")) {
		std::cout << "Multiple Sequence Alignment of Dna sequences" << std::endl;
		return customizedMsaAlignment(cfgOpt, Dna() );
	} else if ((value(cfgOpt, "alphabet") == "protein")) {
		std::cout << "Multiple Sequence Alignment of AminoAcid sequences" << std::endl;
		return customizedMsaAlignment(cfgOpt, AminoAcid() );
	}
	else {	std::cerr << valueHelp(cfgOpt) << std::endl; return -1; }
}
