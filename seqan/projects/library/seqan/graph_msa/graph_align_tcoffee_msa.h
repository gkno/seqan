 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: graph_align_tcoffee_msa.h 1911 2008-05-02 09:28:04Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_MSA_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_MSA_H

namespace SEQAN_NAMESPACE_MAIN
{

template<typename TAlphabet, typename TScore>
struct MsaOptions {
public:
	// Rescore segment matches after refinement
	bool rescore;

	// Output format
	// 0: Fasta
	// 1: Msf
	unsigned int outputFormat;

	// Scoring object
	TScore sc;

	// All methods to compute segment matches
	// 0: global alignment
	// 1: local alignments
	// 2: overlap alignments
	// 3: longest common subsequence
	String<unsigned int> method;
	
	// Various input and output file names
	String<std::string> alnfiles;		// External alignment files
	String<std::string> libfiles;		// T-Coffee library files
	String<std::string> blastfiles;		// Blast match files
	String<std::string> mummerfiles;	// MUMmer files
	std::string outfile;				// Output file name
	std::string seqfile;				// Sequence file name
	std::string infile;					// Alignment file for alignment evaluation
	std::string treefile;				// Guide tree file
	
	// Initialization
	MsaOptions() : rescore(true), outputFormat(0) {}
};

//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TScore>
inline void
evaluateAlignment(MsaOptions<TAlphabet, TScore> const& msaOpt) {
	typedef typename Value<TScore>::Type TScoreValue;
	typedef String<TAlphabet> TSequence;
	typedef typename Size<TSequence>::Type TSize;
	StringSet<TSequence, Owner<> > origStrSet;
	StringSet<String<char> > names;

	// Read the sequences
	std::fstream strm;
	strm.open(msaOpt.infile.c_str(), ::std::ios_base::in | ::std::ios_base::binary);
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
	strm_lib.open(msaOpt.infile.c_str(), ::std::ios_base::in | ::std::ios_base::binary);
	read(strm_lib,matches, scores, names, FastaAlign());	
	strm_lib.close();

	// Build the alignment graph
	typedef Graph<Alignment<TDepSequenceSet, TSize> > TGraph;
	TGraph g(strSet);
	buildAlignmentGraph(matches, g, FrequencyCounting() );

	// Print the scoring information
	TScoreValue gop = msaOpt.sc.data_gap_open;
	TScoreValue gex = msaOpt.sc.data_gap_extend;
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
			std::cout << score(msaOpt.sc, TAlphabet(row), TAlphabet(col));
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
	String<char> mat;
	if (convertAlignment(g, mat)) {
		TScoreValue alignScore = alignmentEvaluation(g, msaOpt.sc, numGapEx, numGap, numPairs, pairCount, alignLen);
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
	} else {
		std::cout << "No valid alignment!" << std::endl;
	}
}


//////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TStringSet1, typename TNames, typename TAlphabet, typename TScore>
inline void
globalMsaAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& gAlign, 
				   TStringSet1& sequenceSet,
				   TNames& sequenceNames,
				   MsaOptions<TAlphabet, TScore> const& msaOpt)
{
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	typedef Graph<Alignment<TStringSet, TSize> > TGraph;
	typedef typename Id<TGraph>::Type TId;
	
	// Initialize alignment object
	clear(gAlign);
	assignStringSet(gAlign, sequenceSet);

#ifdef SEQAN_PROFILE
	std::cout << "Scoring parameters:" << std::endl;
	std::cout << "*Gap opening: " << msaOpt.sc.data_gap_open << std::endl;
	std::cout << "*Gap extension: " << msaOpt.sc.data_gap_extend << std::endl;
	//std::cout << "*Scoring matrix: " << std::endl;
	//TSize alphSize = ValueSize<TAlphabet>::VALUE;
	//std::cout << "   ";
	//for(TSize col = 0; col<alphSize; ++col) std::cout << TAlphabet(col) << ',';
	//std::cout << std::endl;
	//for(TSize row = 0; row<alphSize; ++row) {
	//	for(TSize col = 0; col<alphSize; ++col) {
	//		if (col == 0) std::cout << TAlphabet(row) << ": ";
	//		std::cout << score(msaOpt.sc, TAlphabet(row), TAlphabet(col));
	//		if (col < alphSize - 1) std::cout << ',';
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;
#endif

	// Some alignment constants
	TStringSet& seqSet = stringSet(gAlign);
	TSize nSeq = length(seqSet);
	TSize threshold = 30;

	// Select all possible pairs for global and local alignments
	String<TSize> pList;
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

	// Include segment matches from subalignments
	if (!empty(msaOpt.alnfiles)) {
#ifdef SEQAN_PROFILE
		std::cout << "Parsing external alignment files:" << std::endl;
#endif
		typedef typename Iterator<String<std::string>, Standard>::Type TIter;
		TIter begIt = begin(msaOpt.alnfiles, Standard() );
		TIter begItEnd = end(msaOpt.alnfiles, Standard() );
		for(;begIt != begItEnd; goNext(begIt)) {
#ifdef SEQAN_PROFILE
			std::cout << "*External file " << (*begIt) << std::endl;
#endif
			std::ifstream strm_lib;
			strm_lib.open((*begIt).c_str(), ::std::ios_base::in | ::std::ios_base::binary);
			read(strm_lib, matches, scores, sequenceNames, FastaAlign());
			strm_lib.close();
		}
#ifdef SEQAN_PROFILE
		std::cout << "External segment matches done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
	}

	// Include computed segment matches
	if (!empty(msaOpt.method)) {
#ifdef SEQAN_PROFILE
		std::cout << "Computing segment matches:" << std::endl;
#endif
		typedef typename Iterator<String<unsigned int>, Standard>::Type TIter;
		TIter begIt = begin(msaOpt.method, Standard() );
		TIter begItEnd = end(msaOpt.method, Standard() );
		for(;begIt != begItEnd; goNext(begIt)) {
#ifdef SEQAN_PROFILE
			if (*begIt == 0) std::cout << "*Method: global" << std::endl;
			else if (*begIt == 1) std::cout << "*Method: local" << std::endl;
			else if (*begIt == 2) std::cout << "*Method: overlap" << std::endl;
			else if (*begIt == 3) std::cout << "*Method: lcs" << std::endl;
#endif
			if (*begIt == 0) appendSegmentMatches(seqSet, pList, msaOpt.sc, matches, scores, distanceMatrix, GlobalPairwise_Library() );
			else if (*begIt == 1) appendSegmentMatches(seqSet, pList, msaOpt.sc, matches, scores, LocalPairwise_Library() );
			else if (*begIt == 2) {
				Nothing noth;
				AlignConfig<true,true,true, true> ac;
				appendSegmentMatches(seqSet, pList, msaOpt.sc, matches, scores, noth, ac, GlobalPairwise_Library() );
			} 
			else if (*begIt == 3) appendSegmentMatches(seqSet, pList, matches, scores, Lcs_Library() );
			else {
#ifdef SEQAN_PROFILE
				std::cout << "*Unknown method!!!" << std::endl;
#endif
			}
#ifdef SEQAN_PROFILE
			std::cout << "*Done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
		}	
	}

	// Include a T-Coffee library
	if (!empty(msaOpt.libfiles)) {
#ifdef SEQAN_PROFILE
		std::cout << "Parsing a T-Coffee Library:" << std::endl;
#endif
		typedef typename Iterator<String<std::string>, Standard>::Type TIter;
		TIter begIt = begin(msaOpt.libfiles, Standard() );
		TIter begItEnd = end(msaOpt.libfiles, Standard() );
		for(;begIt != begItEnd; goNext(begIt)) {
#ifdef SEQAN_PROFILE
			std::cout << "*T-Coffee library: " << (*begIt) << std::endl;
#endif
			std::ifstream strm_lib;
			strm_lib.open((*begIt).c_str(), std::ios_base::in | std::ios_base::binary);
			read(strm_lib, matches, scores, sequenceNames, TCoffeeLib());
			strm_lib.close();
		}
#ifdef SEQAN_PROFILE
		std::cout << "Parsing done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
	}

	// Include MUMmer segment matches
	if (!empty(msaOpt.mummerfiles)) {
#ifdef SEQAN_PROFILE
		std::cout << "Parsing MUMmer segment matches:" << std::endl;
#endif
		typedef typename Iterator<String<std::string>, Standard>::Type TIter;
		TIter begIt = begin(msaOpt.mummerfiles, Standard() );
		TIter begItEnd = end(msaOpt.mummerfiles, Standard() );
		for(;begIt != begItEnd; goNext(begIt)) {
#ifdef SEQAN_PROFILE
			std::cout << "*MUMmer file: " << (*begIt) << std::endl;
#endif
			std::ifstream strm_lib;
			strm_lib.open((*begIt).c_str(), std::ios_base::in | std::ios_base::binary);
			read(strm_lib, matches, scores, seqSet, sequenceNames, MummerLib());		
			strm_lib.close();
		}
#ifdef SEQAN_PROFILE
		std::cout << "Parsing done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
	}

	// Include BLAST segment matches
	if (!empty(msaOpt.blastfiles)) {
#ifdef SEQAN_PROFILE
		std::cout << "Parsing BLAST segment matches:" << std::endl;
#endif
		typedef typename Iterator<String<std::string>, Standard>::Type TIter;
		TIter begIt = begin(msaOpt.blastfiles, Standard() );
		TIter begItEnd = end(msaOpt.blastfiles, Standard() );
		for(;begIt != begItEnd; goNext(begIt)) {
#ifdef SEQAN_PROFILE
			std::cout << "*BLAST file: " << (*begIt) << std::endl;
#endif
			std::ifstream strm_lib;
			strm_lib.open((*begIt).c_str(), std::ios_base::in | std::ios_base::binary);
			read(strm_lib, matches, scores, sequenceNames, BlastLib());
			strm_lib.close();
		}
#ifdef SEQAN_PROFILE
		std::cout << "Parsing done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
	}

#ifdef SEQAN_PROFILE
	std::cout << "Total number of segment matches: " << length(matches) << std::endl;
#endif
	// Score the matches
	if (msaOpt.rescore) {
#ifdef SEQAN_PROFILE	
		std::cout << "Scoring method: Re-Score" << std::endl;
#endif
		scoreMatches(seqSet, msaOpt.sc, matches, scores);
#ifdef SEQAN_PROFILE
		std::cout << "Scoring done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
	}


	// Use these segment matches for the initial alignment graph
	TGraph g(seqSet);
#ifdef SEQAN_PROFILE
	std::cout << "Construction of alignment graph: FractionalScore" << std::endl;
#endif
	buildAlignmentGraph(matches, scores, g, FractionalScore() );
	clear(matches);
	clear(scores);
#ifdef SEQAN_PROFILE
	std::cout << "Alignment graph construction done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	// Guide tree
	Graph<Tree<double> > guideTree;
	if (!msaOpt.treefile.empty()) {
#ifdef SEQAN_PROFILE
		std::cout << "Guide Tree: " << msaOpt.treefile << std::endl;
#endif
		std::fstream strm_tree;
		strm_tree.open(msaOpt.treefile.c_str(), ::std::ios_base::in | ::std::ios_base::binary);
		read(strm_tree,guideTree,sequenceNames,NewickFormat());	// Read newick tree
		strm_tree.close();
	} else {
#ifdef SEQAN_PROFILE
		std::cout << "Guide Tree: Neighbor Joining" << std::endl;
#endif
		// Check if we have a valid distance matrix
		if (empty(distanceMatrix)) getDistanceMatrix(g, distanceMatrix, KmerDistance());
		slowNjTree(distanceMatrix, guideTree);
	}
#ifdef SEQAN_PROFILE
	std::cout << "Guide tree done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
		
	// Triplet extension
	if (nSeq < threshold) tripletLibraryExtension(g);
	else tripletLibraryExtension(g, guideTree, threshold / 2);
#ifdef SEQAN_PROFILE
	std::cout << "Triplet extension done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	// Progressive Alignment
	progressiveAlignment(g, guideTree, gAlign);
#ifdef SEQAN_PROFILE
	std::cout << "Progressive alignment done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	clear(guideTree);
	clear(distanceMatrix);
	clear(g);
#ifdef SEQAN_PROFILE
	std::cout << "Clean-up done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

}


//////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TScore>
inline void
globalMsaAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& gAlign,
				   TScore const& scoreObject)
{
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Value<TString>::Type TAlphabet;
	TStringSet sequenceSet = stringSet(gAlign);
	String<String<char> > sequenceNames;
	fill(sequenceNames, length(sequenceSet), String<char>("tmpName"));
	MsaOptions<AminoAcid, TScore> msaOpt;
	msaOpt.sc = scoreObject;
	appendValue(msaOpt.method, 0);  // Global pairwise
	appendValue(msaOpt.method, 1);	// Local pairwise
	globalMsaAlignment(gAlign, sequenceSet, sequenceNames, msaOpt);
}

//////////////////////////////////////////////////////////////////////////////


template<typename TSource, typename TSpec, typename TScore>
inline void
globalMsaAlignment(Align<TSource, TSpec>& align,
				   TScore const& scoreObject)
{
	typedef StringSet<TSource, Dependent<> > TStringSet;
	TStringSet sequenceSet = stringSet(align);
	Graph<Alignment<TStringSet, void, WithoutEdgeId> > gAlign(sequenceSet);
	globalMsaAlignment(gAlign, scoreObject);
	
	// Pipe into Align data structure
	String<char> mat;
	convertAlignment(gAlign, mat);
	typedef Align<TSource, TSpec> TAlign;
	typedef typename Size<TAlign>::Type TSize;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow>::Type TRowIterator;
	clearGaps(align);
	TSize nseq = length(sequenceSet);
	String<TRowIterator> rowIter;
	resize(rowIter, nseq);
	for(TSize i = 0; i<nseq; ++i) value(rowIter, i) = begin(row(align, i));
	TSize lenMat = length(mat);
	TSize colLen = lenMat / nseq;
	TSize gapCount = 0;
	char gapChar = gapValue<char>();
	for(TSize alignRow = 0; alignRow < nseq; ++alignRow) {
		for(TSize pos = alignRow * colLen; pos < (alignRow + 1) * colLen; ++pos) {
			if (value(mat, pos) != gapChar) {
				if (gapCount) {
					insertGaps(value(rowIter, alignRow), gapCount);
					goFurther(value(rowIter, alignRow), gapCount);
					gapCount = 0;
				}
				goNext(value(rowIter,alignRow));
			} else ++gapCount;
		}
		if (gapCount) {
			insertGaps(value(rowIter, alignRow), gapCount);
			goFurther(value(rowIter, alignRow), gapCount);
			gapCount = 0;
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
// The 2 most useful functions in SeqAn!!!
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TMatches>
void
_debugMatches(TStringSet& str, 
			  TMatches& matches)
{
	typedef typename Id<TStringSet>::Type TId;
	typedef typename Size<TStringSet>::Type TSize;

	// Print all the matches
	std::cout << "The sequences:" << std::endl;
	for(TSize i = 0;i<length(str);++i) {
		std::cout << positionToId(str,i) << ':' << str[i] << std::endl;
	}
	std::cout << "The matches:" << std::endl;
	for(TSize i = 0;i<length(matches);++i) {
		TId tmp_id1 = sequenceId(matches[i],0);
		std::cout << tmp_id1 << ',' << fragmentBegin(matches[i],tmp_id1) << ',';
		for(TSize j = fragmentBegin(matches[i],tmp_id1); j < fragmentBegin(matches[i],tmp_id1) + fragmentLength(matches[i],tmp_id1); ++j) {
			std::cout << str[idToPosition(str, tmp_id1)][j];
		}
		TId tmp_id2 = sequenceId(matches[i],1);
		std::cout << ',' <<	tmp_id2 << ',' << fragmentBegin(matches[i],tmp_id2) << ',';
		for(TSize j = fragmentBegin(matches[i],tmp_id2); j < fragmentBegin(matches[i],tmp_id2) + fragmentLength(matches[i],tmp_id2); ++j) {
			std::cout << str[idToPosition(str, tmp_id2)][j];
		}
		std::cout << std::endl;

		SEQAN_TASSERT(sequenceId(matches[i],0) != sequenceId(matches[i],1))
		SEQAN_TASSERT(fragmentBegin(matches[i],tmp_id1) < length(str[idToPosition(str, tmp_id1)]))
		SEQAN_TASSERT(fragmentBegin(matches[i],tmp_id1) + fragmentLength(matches[i],tmp_id1) <= length(str[idToPosition(str, tmp_id1)]))
		SEQAN_TASSERT(fragmentBegin(matches[i],tmp_id2) < length(str[idToPosition(str, tmp_id2)]))
		SEQAN_TASSERT(fragmentBegin(matches[i],tmp_id2) + fragmentLength(matches[i],tmp_id2) <= length(str[idToPosition(str, tmp_id2)]))
		SEQAN_TASSERT(fragmentLength(matches[i],tmp_id2) == fragmentLength(matches[i],tmp_id1))
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph>
void
_debugRefinedMatches(TGraph& g)
{
	typedef typename Id<TGraph>::Type TId;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	
	std::cout << "Refined matches" << std::endl;
	TEdgeIterator it_tmp(g);
	for(;!atEnd(it_tmp);++it_tmp) {
		TId id1 = sequenceId(g,sourceVertex(it_tmp));
		TId id2 = sequenceId(g,targetVertex(it_tmp));
		std::cout << id1 << ',' << fragmentBegin(g,sourceVertex(it_tmp)) << ',';
		std::cout << label(g,sourceVertex(it_tmp));
		std::cout << ',' <<	id2 << ',' << fragmentBegin(g,targetVertex(it_tmp)) << ',';
		std::cout << label(g,targetVertex(it_tmp));
		std::cout << " (" << cargo(*it_tmp) << ")";
		std::cout << std::endl;	

		SEQAN_TASSERT(sequenceId(g,sourceVertex(it_tmp)) != sequenceId(g,targetVertex(it_tmp)))
		SEQAN_TASSERT(fragmentBegin(g,sourceVertex(it_tmp)) < length((stringSet(g))[idToPosition((stringSet(g)), id1)]))
		SEQAN_TASSERT(fragmentBegin(g,sourceVertex(it_tmp)) + fragmentLength(g,sourceVertex(it_tmp)) <= length((stringSet(g))[idToPosition((stringSet(g)), id1)]))
		SEQAN_TASSERT(fragmentBegin(g,targetVertex(it_tmp)) < length((stringSet(g))[idToPosition((stringSet(g)), id2)]))
		SEQAN_TASSERT(fragmentBegin(g,targetVertex(it_tmp)) + fragmentLength(g,targetVertex(it_tmp)) <= length((stringSet(g))[idToPosition((stringSet(g)), id2)]))
		SEQAN_TASSERT(fragmentLength(g,sourceVertex(it_tmp)) == fragmentLength(g,targetVertex(it_tmp)))

	}
}













////// Testing Stuff
//
////////////////////////////////////////////////////////////////////////////////
//
//template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TScoreMatrixFile, typename TGapOpen,  typename TGapEx, typename TAlignmentGraph>
//inline void
//testChristian(StringSet<TString, TSpec> const& seqSet,
//			  StringSet<TName, TSpec2> const&,
//			  TScoreMatrixFile& scoreMatrix,
//			  TGapEx& gex,
//			  TGapOpen& gop,
//			  TAlignmentGraph& gOut)
//{
//	SEQAN_CHECKPOINT
//
//	typedef StringSet<TString, Dependent<> > TStringSet;
//	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
//	typedef typename Size<TGraph>::Type TSize;
//	typedef typename Id<TGraph>::Type TId;
//	TSize nSeq = length(seqSet);
//	TSize threshold = 30;
//	
//	// Select informative pairs
//	TGraph g(seqSet);
//	String<Pair<TId, TId> > pList;
//	selectPairs(g, pList);
//
//	// Generate a primary library, i.e., all global pairwise alignments
//	String<double> distanceMatrix;
//	double gopening = 0; 
//	double gextend = 0; 
//	std::stringstream ssStream1(toCString(gop));
//	ssStream1 >> gopening; 
//	std::stringstream ssStream2(toCString(gex));
//	ssStream2 >> gextend;
//	typedef Score<double, ScoreMatrix<> > TScore;
//	TScore score_type_global;
//	if (!length(scoreMatrix)) {
//		Blosum62 sc(-1 , -11);
//		convertScoringMatrix(sc, score_type_global, -1.0 * gextend, -1.0 * gopening);
//	} else {
//		loadScoreMatrix(score_type_global, scoreMatrix);
//		score_type_global.data_gap_extend = (double) (-1.0 * gextend);
//		score_type_global.data_gap_open = (double) (-1.0 * gopening);
//	}
//	generatePrimaryLibrary(g, pList, distanceMatrix, score_type_global, GlobalPairwise_Library() );
//
//	if (nSeq < threshold) {
//		// Full triplet...
//		tripletLibraryExtension(g);
//
//		// ... and normal progressive alignment with guide tree
//		Graph<Tree<double> > guideTree;
//		//upgmaTree(distanceMatrix, guideTree, UpgmaMin() );	
//		slowNjTree(distanceMatrix, guideTree);
//		progressiveAlignment(g, guideTree, gOut);
//		clear(guideTree);
//	} else {
//		// Triplet only on groups of sequences
//		Graph<Tree<double> > guideTree;
//		//upgmaTree(distanceMatrix, guideTree);	
//		slowNjTree(distanceMatrix, guideTree);	// More balanced than UPGMA
//		progressiveAlignment(g, guideTree, gOut, threshold);
//		clear(guideTree);
//	}
//	clear(distanceMatrix);
//	clear(g);
//
//	std::cout << "Alignment Score: " << sumOfPairsScore(gOut, score_type_global) << std::endl;
//}
//
//
//
//template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TFileName, typename TAlignmentGraph>
//inline void
//testIslands(StringSet<TString, TSpec> const& seqSet,
//			 StringSet<TName, TSpec2> const& nameSet,
//			 TFileName& matchfile,
//			 TAlignmentGraph& gOut)
//{
//	typedef StringSet<TString, Dependent<> > TStringSet;
//	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
//	typedef typename Size<TGraph>::Type TSize;
//	typedef typename Id<TGraph>::Type TId;
//	TSize nseq = length(seqSet);
//	TSize threshold = 38;
//
//	// Select all pairs
//	TGraph g(seqSet);
//	String<Pair<TId, TId> > pList;
//	selectPairs(g, pList);
//
//
//	TGraph lib1(seqSet);
//	String<double> distanceMatrix;
//	Blosum62 scType(-1,-11);
//	generatePrimaryLibrary(lib1, pList, distanceMatrix, scType, GlobalPairwise_Library() );
//	
//	TGraph lib2(seqSet);
//	generatePrimaryLibrary(lib2, pList, scType, Interleaved_Library() );
//
//	// Weighting of libraries (Signal addition)
//	String<TGraph*> libs;
//	appendValue(libs, &lib2);
//	//appendValue(libs, &lib2);
//	//appendValue(libs, &lib2);
//	combineGraphs(g, libs);
//	//combineGraphs(g, libs, FrequencyCounting() );
//
//	if (nseq < threshold) {
//		// Full triplet...
//		tripletLibraryExtension(g);
//
//		// ... and normal progressive alignment with guide tree
//		Graph<Tree<double> > guideTree;
//		slowNjTree(distanceMatrix, guideTree);
//		progressiveAlignment(g, guideTree, gOut);
//		clear(guideTree);
//	} else {
//		// Triplet only on groups of sequences
//		Graph<Tree<double> > guideTree;
//		slowNjTree(distanceMatrix, guideTree);	// More balanced than UPGMA
//		progressiveAlignment(g, guideTree, gOut, threshold);
//		clear(guideTree);
//	}
//
//	std::cout << gOut << std::endl;
//}
//
//
//
//
//
//template<typename TStringSet, typename TCargo, typename TSpec, typename TString, typename TMatches, typename TSize>
//inline void 
//__alignSeqFragments(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
//					TString& unaligned1,
//					TString& unaligned2,
//					TMatches& matches,
//					TSize seq1,
//					TSize seq2,
//					TSize off1,
//					TSize off2)
//{
//	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
//	typedef typename Size<TGraph>::Type TId;
//	typedef typename Value<TMatches>::Type TFragment;
//	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
//	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
//
//	TStringSet& seqSet = stringSet(g);
//	if ((length(unaligned1) > 0) && (length(unaligned2) > 0)) {
//		// Make a pairwise string-set
//		TStringSet pairSet;
//		TId id1 = positionToId(seqSet, seq1);
//		TId id2 = positionToId(seqSet, seq2);
//		assignValueById(pairSet, unaligned1, id1);
//		assignValueById(pairSet, unaligned2, id2);	
//		TGraph pairAlign(pairSet);
//		
//		Blosum62 st1(-1,-11);
//		TSize from = length(matches);
//		globalAlignment(matches, pairSet, st1, Gotoh() );
//		TSize to = length(matches);
//		
//		for(TSize iter = from; iter < to; ++iter) {
//			fragmentBegin(matches[iter], id1) += off1;
//			fragmentBegin(matches[iter], id2) += off2;
//		}
//	}
//}
//
//
//
//template<typename TStringSet, typename TCargo, typename TSpec, typename TScore, typename TConsensus>
//inline void
//__getAnchorGraph(Graph<Alignment<TStringSet, TCargo, TSpec> >& anchorGraph,
//					TScore const& score_type,
//					TConsensus& consensus)
//{
//	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
//	typedef typename Id<TGraph>::Type TId;
//	typedef typename Size<TGraph>::Type TSize;
//	TSize nSeq = length(stringSet(anchorGraph));
//	
//
//	// String of fragments to combine all pairwise alignments into a multiple alignment
//	typedef Fragment<> TFragment;
//	typedef String<TFragment> TFragmentString;
//	TFragmentString matches;
//
//	TStringSet consSet = stringSet(anchorGraph);
//	appendValue(consSet, consensus);
//	TGraph tmpGraph(consSet);
//	TId idCons = positionToId(consSet, nSeq);
//	TGraph bla(consSet);
//	for(TSize seq = 0; seq < nSeq; ++seq) {
//		TStringSet pairSet;
//		TId id1 = positionToId(consSet, seq);
//		assignValueById(pairSet, consSet, id1);
//		assignValueById(pairSet, consSet, idCons);
//		
//		TGraph bla(pairSet);
//		Blosum62 st(-1,-20);
//		globalAlignment(matches, pairSet, st, AlignConfig<true, false, false, true>(), Gotoh() );	
//		globalAlignment(bla, st, AlignConfig<true, false, false, true>(), Gotoh() );	
//		std::cout << bla << std::endl;
//		clearVertices(bla);
//		globalAlignment(bla, st, AlignConfig<false, true, true, false>(), Gotoh() );	
//		std::cout << bla << std::endl;
//
//	}
//
//	// Refine all matches and create multiple alignment
//	matchRefinement(matches,consSet,const_cast<TScore&>(score_type),tmpGraph);
//
//	String<double> initialDistanceMatrix;
//	getDistanceMatrix(tmpGraph, initialDistanceMatrix, 6, AAGroupsDayhoff(), KmerDistance() );
//	Graph<Tree<double> > guideTree;
//	slowNjTree(initialDistanceMatrix, guideTree);
//
//	tripletLibraryExtension(tmpGraph);
//	TGraph gOut(consSet);
//	progressiveAlignment(tmpGraph, guideTree, gOut);
//
//	clearVertices(anchorGraph);
//	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
//	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
//	TVertexIterator itV(gOut);
//	for(;!atEnd(itV);++itV) {
//		if (sequenceId(gOut, *itV) != idCons) addVertex(anchorGraph, sequenceId(gOut, *itV), fragmentBegin(gOut, *itV), fragmentLength(gOut,*itV));
//	}
//	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
//	TEdgeIterator itE(gOut);
//	for(;!atEnd(itE);++itE) {
//		TVertexDescriptor sV = sourceVertex(itE);
//		TVertexDescriptor tV = targetVertex(itE);
//		if ((sequenceId(gOut, sV) != idCons) && (sequenceId(gOut, tV) != idCons)) {
//			addEdge(anchorGraph, findVertex(anchorGraph, sequenceId(gOut, sV), fragmentBegin(gOut, sV)), findVertex(anchorGraph, sequenceId(gOut, tV), fragmentBegin(gOut, tV)), 1);
//		}
//	}
//}
//
//template<typename TSet1, typename TSet2, typename TSet3>
//inline void
//getCommonSubset(TSet1 const& set1,
//				TSet2 const& set2,
//				TSet3& common)
//{
//	typename TSet1::const_iterator pos1 = set1.begin();
//	typename TSet1::const_iterator pos1End = set1.end();
//	for(;pos1 != pos1End; ++pos1) {
//		if (set2.find(*pos1) != set2.end()) common.insert(*pos1);
//	}
//}
//
//
//
//
//template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TFileName, typename TAlignmentGraph>
//inline void
//testProfiles(StringSet<TString, TSpec> const& seqSet,
//			 StringSet<TName, TSpec2> const& nameSet,
//			 TFileName& matchfile,
//			 TAlignmentGraph& gOut)
//{
//	typedef StringSet<TString, Dependent<> > TStringSet;
//	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
//	typedef typename Size<TGraph>::Type TSize;
//	typedef typename Id<TGraph>::Type TId;
//	TSize nSeq = length(seqSet);
//	TSize threshold = 30;
//
//	// Generate an anchor graph from local alignments
//	TGraph lib1(seqSet);
//	Blosum62 score_type_global(-1,-11);
//	generatePrimaryLibrary(lib1, score_type_global, LocalExtendPairwise_Library() );
//
//	// Build the anchor graph but keep all local alignments
//	TGraph tmpGraph(lib1);
//	cliqueReduction(tmpGraph);
//	String<double> distanceMatrix;
//	getDistanceMatrix(tmpGraph,distanceMatrix,LibraryDistance() );
//	Graph<Tree<double> > guideTree;
//	slowNjTree(distanceMatrix, guideTree);
//	TGraph anchorGraph(seqSet);
//	progressiveAlignment(tmpGraph, guideTree, anchorGraph);
//	clear(guideTree); clear(tmpGraph);
//	reduceGraph(anchorGraph, score_type_global, distanceMatrix);
//	std::cout << anchorGraph << std::endl;
//
//	//// Writing
//	//std::fstream strmBlast;
//	//strmBlast.open(toCString(matchfile), std::ios_base::out | std::ios_base::trunc);
//	//write(strmBlast, anchorGraph, nameSet, BlastLib());
//	//strmBlast.close();
//
//	// Build the new guide tree
//	slowNjTree(distanceMatrix, guideTree);
//	
//	// Make a banded alignment
//	TGraph refinedLib(seqSet);
//	String<char> mat;
//	char gapChar = gapValue<char>();
//	convertAlignment(anchorGraph, mat);
//	TSize len = length(mat) / nSeq;
//	typedef Fragment<> TFragment;
//	typedef String<TFragment> TFragmentString;
//	TFragmentString matches;
//	for(TSize seq1 = 0; seq1 < nSeq - 1; ++seq1) {
//		for(TSize seq2 = seq1 + 1; seq2 < nSeq; ++seq2) {
//			//std::cout << '#' << seq1 << ',' << seq2 << std::endl;
//			TSize offset1 = 0;
//			TSize offset2 = 0;
//			TString unaligned1;
//			TString unaligned2;
//			for(TSize col = 0; col<len; ++col) {
//				if (value(mat, seq1 * len + col) != gapChar) {
//					if (value(mat, seq2 * len + col) != gapChar) {
//						__alignSeqFragments(refinedLib, unaligned1, unaligned2, matches, seq1, seq2, (TSize) (offset1 - length(unaligned1)), (TSize) (offset2 - length(unaligned2)));
//						clear(unaligned1);
//						clear(unaligned2);
//						++offset1;
//						++offset2;
//					} else {
//						++offset1;
//						appendValue(unaligned1, value(mat, seq1 * len + col)); 
//					}
//				} else if (value(mat, seq2 * len + col) != gapChar) {
//					++offset2;
//					appendValue(unaligned2, value(mat, seq2 * len + col)); 
//				}
//			}
//			__alignSeqFragments(refinedLib, unaligned1, unaligned2, matches, seq1, seq2, (TSize) (offset1 - length(unaligned1)), (TSize) (offset2 - length(unaligned2)));
//			clear(unaligned1);
//			clear(unaligned2);
//		}
//	}
//
//		
//	// Refine all matches and create multiple alignment
//	matchRefinement(matches,stringSet(refinedLib),refinedLib);
//
//	// Weighting of libraries (Signal addition)
//	TGraph  g(seqSet);
//	String<TGraph*> libs;
//	appendValue(libs, &anchorGraph);
//	appendValue(libs, &anchorGraph);
//	//appendValue(libs, &anchorGraph);
//	appendValue(libs, &lib1);
//	appendValue(libs, &refinedLib);
//	combineGraphs(g, libs, FrequencyCounting() );
//		
//	// Clear the old libraries
//	clear(anchorGraph);
//	clear(refinedLib);
//	clear(lib1);
//	
//	if (nSeq < threshold) {
//		// Full triplet...
//		tripletLibraryExtension(g);
//		// ... and normal progressive alignment with guide tree
//		progressiveAlignment(g, guideTree, gOut);
//	} else {
//		progressiveAlignment(g, guideTree, gOut, threshold);
//	}
//	//std::cout << gOut << std::endl;
//	clear(guideTree);
//	clear(distanceMatrix);
//	clear(g);
//
//}
//
//
//template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TFileName, typename TAlignmentGraph>
//inline void
//testIslands(StringSet<TString, TSpec> const& seqSet,
//			 StringSet<TName, TSpec2> const& nameSet,
//			 TFileName& matchfile,
//			 TAlignmentGraph& gOut)
//{
//	typedef StringSet<TString, Dependent<> > TStringSet;
//	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
//	typedef typename Size<TGraph>::Type TSize;
//	typedef typename Id<TGraph>::Type TId;
//	TSize nSeq = length(seqSet);
//	TSize threshold = 30;
//
//	// Select all pairs
//	TGraph g(seqSet);
//	String<Pair<TId, TId> > pList;
//	for(TSize i=0; i<nSeq; ++i) {
//		for(TSize j=i+1; j<i+3; ++j) {
//			appendValue(pList, Pair<TId, TId>(positionToId(seqSet, i), positionToId(seqSet, j % nSeq)));
//		}
//	}
//	appendValue(pList, pList[0]);
//
//	typedef Fragment<> TFragment;
//	typedef String<TFragment > TFragmentString;
//	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
//	TFragmentString matches;
//	TFragmentString localMatches;
//	TSize lastStart = 0;
//	Blosum62 scLocalType(-2,-8);
//	TSize k = 0;
//	while(k<length(pList)) {
//		// Make a pairwise string-set
//		TStringSet pairSet;
//		TId id1 = (pList[k]).i1;
//		TId id2 = (pList[k]).i2;
//		assignValueById(pairSet, stringSet(g), id1);
//		assignValueById(pairSet, stringSet(g), id2);
//	
//		// Alignment
//		lastStart = length(localMatches);
//		globalAlignment(localMatches, pairSet, scLocalType, AlignConfig<true, true, false, false>(), Gotoh() );
//
//		if ((k!=0) && (k % 2 == 0)) {
//			TFragmentString copyLocalMatches;
//			for(TSize runThrough = lastStart; runThrough < length(localMatches); ++runThrough) {
//				appendValue(copyLocalMatches, localMatches[runThrough]);
//			}
//			TGraph localGraph(seqSet);
//			matchRefinement(localMatches,stringSet(localGraph), localGraph);
//			cliqueReduction(localGraph);
//			String<unsigned int> components;
//			unsigned int numComps = connected_components(localGraph, components);
//			String<unsigned int> memberCount;
//			fill(memberCount, numComps, 0);
//			String<unsigned int> compLength;
//			fill(compLength, numComps, 0);
//			typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
//			TVertexIterator itVertex(localGraph);
//			for(;!atEnd(itVertex);++itVertex) {
//				memberCount[property(components, *itVertex)] += 1;
//				compLength[property(components, *itVertex)] = fragmentLength(localGraph, *itVertex);
//			}
//			unsigned int maxLength = 0;
//			unsigned int maxComp = 0;
//			for(unsigned int i = 0; i < length(compLength); ++i) {
//				if ((memberCount[i] == 3) && (compLength[i] > maxLength)) {
//					maxLength = compLength[i];
//					maxComp = i;
//				}
//			}
//			if (maxLength > 5) {
//				typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
//				TEdgeIterator itEdge(localGraph);
//				for(;!atEnd(itEdge);++itEdge) {
//					if (property(components, sourceVertex(itEdge)) == maxComp) {
//						appendValue(matches, TFragment(sequenceId(localGraph,sourceVertex(itEdge)), fragmentBegin(localGraph, sourceVertex(itEdge)), sequenceId(localGraph, targetVertex(itEdge)), fragmentBegin(localGraph, targetVertex(itEdge)), fragmentLength(localGraph, targetVertex(itEdge))));
//					}
// 				}
//			}
//			localMatches = copyLocalMatches;
//		}
//		++k;
//	}
//	matchRefinement(matches,stringSet(g),g);
//
//	// Build the anchor graph but keep all local alignments
//	tripletLibraryExtension(g);
//	String<double> distanceMatrix;
//	getDistanceMatrix(g,distanceMatrix,LibraryDistance() );
//	Graph<Tree<double> > guideTree;
//	slowNjTree(distanceMatrix, guideTree);
//	TGraph anchorGraph(seqSet);
//	progressiveAlignment(g, guideTree, anchorGraph);
//	//std::cout << anchorGraph << std::endl;
//
//	// Create an alignment matrix
//	typedef String<char> TAlignMatrix;
//	typedef String<bool> TActive;
//	TSize minBlockLen = 10;
//	TSize minCoverage = 3;
//	TAlignMatrix mat;
//	char gapChar = gapValue<char>();
//	convertAlignment(anchorGraph, mat);
//	TSize lenMat = length(mat) / nSeq;
//	
//	// Mark the gaps and the characters
//	typedef typename Iterator<TAlignMatrix>::Type TAlignIter;
//	typedef typename Iterator<TActive>::Type TActiveIter;
//	TActive active;
//	fill(active, length(mat), true);
//	String<unsigned int> coverage;
//	fill(coverage, lenMat, 0);
//	TAlignIter itA = begin(mat);
//	TAlignIter itAEnd = end(mat);
//	TActiveIter itB = begin(active);
//	TSize col = 0;
//	for(; itA != itAEnd; ++itA, ++itB, col = (++col) % lenMat) {
//		if (*itA == gapChar) value(itB) = false;
//		else value(coverage, col) += 1;
//	}
//
//	// Clear blocks with low coverage
//	itA = begin(mat);
//	itB = begin(active);
//	col = 0;
//	for(;itA != itAEnd; ++itA, ++itB, col = (++col) % lenMat) {
//		if (value(coverage, col) < minCoverage) value(itB) = false;
//	}
//
//	// Collect all blocks of minimum length
//	TSize blockStart = lenMat;
//	TSize blockEnd = 0;
//	typedef std::set<TSize> TSequenceSet;
//	typedef Triple<TSequenceSet, TSize, TSize> TBlock;
//	typedef String<TBlock> TConservedBlock;
//	typedef typename Iterator<TConservedBlock>::Type TConsBlockIter;
//	TConservedBlock conservedBlock;
//	TSequenceSet thisSeqSet;
//	TSequenceSet nextSeqSet;
//	for(TSize column = 0; column < lenMat; ++column) {
//		for(TSize seq = 0; seq < nSeq; ++seq) {
//			if (!value(active, seq * lenMat + column)) continue;
//			nextSeqSet.insert(seq);
//		}
//		if (thisSeqSet != nextSeqSet) {
//			if (blockStart == lenMat) {
//				blockStart = column;
//			} else  {
//				if ((column - blockStart) > minBlockLen) appendValue(conservedBlock, TBlock(thisSeqSet, blockStart, column));
//				if (nextSeqSet.empty()) blockStart = lenMat;
//				else blockStart = column;
//			}
//			thisSeqSet = nextSeqSet;
//		}
//		nextSeqSet.clear();
//	}
//
//	// Get Consensus
//	typedef typename Value<Blosum62>::Type TScoreValue;
//	String<TString> consensus;
//	String<TScoreValue> scores;
//	TConsBlockIter itConsBlock = begin(conservedBlock);
//	TConsBlockIter itConsBlockEnd = end(conservedBlock);
//	for(;itConsBlock != itConsBlockEnd; ++itConsBlock) {
//		TString consString;
//		TScoreValue sumScore = 0;
//		for(TSize column = (value(itConsBlock)).i2; column < (value(itConsBlock)).i3; ++column) {
//			AminoAcid chooseX = AminoAcid(0);
//			TScoreValue maxVal = -1 * _getInfinity<TScoreValue>();
//			for(unsigned int aminoWalker = 0; aminoWalker < 20; ++aminoWalker) {
//				TScoreValue sum = 0;
//				AminoAcid x = AminoAcid(aminoWalker);
//				typename TSequenceSet::const_iterator itCons = (value(itConsBlock)).i1.begin();
//				typename TSequenceSet::const_iterator itConsEnd = (value(itConsBlock)).i1.end();
//				for(;itCons != itConsEnd; ++itCons) {
//					sum += score(scLocalType, value(mat, *itCons * lenMat + column), x);
//				}
//				if (sum > maxVal) {
//					maxVal = sum;
//					chooseX = x;
//				}
//			}
//			appendValue(consString, chooseX);
//			sumScore += maxVal;
//		}
//		appendValue(consensus, consString);
//		appendValue(scores, sumScore);
//	}
//	TScoreValue maxValueCons =  -1 * _getInfinity<TScoreValue>();
//	TSize maxIndexCons = 0;
//	for(TSize i = 0; i<length(scores); ++i) {
//		if (scores[i] > maxValueCons) {
//			maxValueCons = scores[i];
//			maxIndexCons = i;
//		}
//	}
//	
//	// Debug blocks
//	TConsBlockIter itCB = begin(conservedBlock);
//	TConsBlockIter itCBEnd = end(conservedBlock);
//	TSize posCount = 0;
//	for(;itCB != itCBEnd; ++itCB, ++posCount) {
//		typename TSequenceSet::const_iterator itCons = (value(itCB)).i1.begin();
//		typename TSequenceSet::const_iterator itConsEnd = (value(itCB)).i1.end();
//		for(;itCons != itConsEnd; ++itCons) {
//			for(TSize column = (value(itCB)).i2; column < (value(itCB)).i3; ++column) {
//				std::cout << value(mat, *itCons * lenMat + column);
//			}
//			std::cout << std::endl;
//		}
//		std::cout << "Consensus:" << std::endl;
//		std::cout << consensus[posCount] << std::endl;
//		std::cout << scores[posCount] << std::endl;
//		std::cout << std::endl;
//	}
//	
//	TStringSet& strSet = stringSet(anchorGraph);
//	clear(matches);
//	TString myConsensus = consensus[maxIndexCons];
//	TStringSet tmpSet;
//	for(TSize seq = 0; seq < nSeq; ++seq) {
//		TId myId = positionToId(strSet, seq);
//	
//		TStringSet pairSet;
//		assignValueById(tmpSet, strSet[seq], seq);
//		assignValueById(pairSet, strSet[seq], seq);
//		assignValueById(pairSet, myConsensus, nSeq);
//
//		Blosum62 score_type_global(-2,-8);
//
//		typedef String<Pair<TSize, TScoreValue> > TAlignIdAndScore;
//		TAlignIdAndScore alignIdScore;
//		multiLocalAlignment(matches, pairSet, alignIdScore, score_type_global, 3, SmithWatermanIsland() );
//	}
//	assignValueById(tmpSet, myConsensus, nSeq);
//	TGraph tmpGraph(tmpSet);
//	matchRefinement(matches,stringSet(tmpGraph),tmpGraph);
//
//	tripletLibraryExtension(tmpGraph);
//
//	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
//	TEdgeIterator it_tmp(tmpGraph);
//	clear(matches);
//	for(;!atEnd(it_tmp);++it_tmp) {
//		TId id1 = sequenceId(tmpGraph,sourceVertex(it_tmp));
//		TId id2 = sequenceId(tmpGraph,targetVertex(it_tmp));
//		if ((id1 == nSeq) || (id2 ==nSeq)) continue;
//		id1 = positionToId(strSet, id1);
//		id2 = positionToId(strSet, id2);
//		appendValue(matches, TFragment(id1, fragmentBegin(tmpGraph, sourceVertex(it_tmp)), id2, fragmentBegin(tmpGraph, targetVertex(it_tmp)), fragmentLength(tmpGraph, targetVertex(it_tmp))));
//	}
//
//	clearVertices(g);
//	matchRefinement(matches,stringSet(g),g);
//
//	tripletLibraryExtension(g);
//
//	progressiveAlignment(g, guideTree, gOut);
//	
//	//std::cout << gOut << std::endl;
//
//}
////////////////////////////////////////////////////////////////////////////////
//
//
//template<typename TStringSet, typename TCargo, typename TSpec, typename TScore, typename TDistanceMatrix>
//inline void
//reduceGraph(Graph<Alignment<TStringSet, TCargo, TSpec> >& anchorGraph,
//			TScore const& score_type,
//			TDistanceMatrix& distMatrix)
//{
//	typedef String<char> TAlignMatrix;
//	typedef String<bool> TActive;
//	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
//	typedef typename Size<TGraph>::Type TSize;
//	typedef typename Value<TScore>::Type TScoreValue;
//
//	// Anchor-graph parameters
//	TSize minCoverage = 3;
//	TSize minBlockLen = 5;
//	
//	// Create an alignment matrix
//	TAlignMatrix mat;
//	char gapChar = gapValue<char>();
//	convertAlignment(anchorGraph, mat);
//	TStringSet& seqSet = stringSet(anchorGraph);
//	TSize nseq = length(seqSet);
//	TSize len = length(mat) / nseq;
//	
//	// Mark the gaps and the characters
//	typedef typename Iterator<TAlignMatrix>::Type TAlignIter;
//	typedef typename Iterator<TActive>::Type TActiveIter;
//	TActive active;
//	fill(active, length(mat), true);
//	String<unsigned int> coverage;
//	fill(coverage, len, 0);
//	TAlignIter itA = begin(mat);
//	TAlignIter itAEnd = end(mat);
//	TActiveIter itB = begin(active);
//	TSize col = 0;
//	for(; itA != itAEnd; ++itA, ++itB, col = (++col) % len) {
//		if (*itA == gapChar) value(itB) = false;
//		else value(coverage, col) += 1;
//	}
//
//	// Clear blocks with low coverage
//	itA = begin(mat);
//	itB = begin(active);
//	col = 0;
//	for(;itA != itAEnd; ++itA, ++itB, col = (++col) % len) {
//		if (value(coverage, col) < minCoverage) value(itB) = false;
//	}
//
//	// Collect all blocks of minimum length
//	TSize blockStart = len;
//	TSize blockEnd = 0;
//	typedef std::set<TSize> TSequenceSet;
//	typedef Triple<TSequenceSet, TSize, TSize> TBlock;
//	typedef String<TBlock> TConservedBlock;
//	TConservedBlock conservedBlock;
//	TSequenceSet thisSeqSet;
//	TSequenceSet nextSeqSet;
//	for(TSize column = 0; column < len; ++column) {
//		for(TSize seq = 0; seq < nseq; ++seq) {
//			if (!value(active, seq * len + column)) continue;
//			nextSeqSet.insert(seq);
//		}
//		if (thisSeqSet != nextSeqSet) {
//			if (blockStart == len) {
//				blockStart = column;
//			} else  {
//				if ((column - blockStart) > minBlockLen) appendValue(conservedBlock, TBlock(thisSeqSet, blockStart, column));
//				if (nextSeqSet.empty()) blockStart = len;
//				else blockStart = column;
//			}
//			thisSeqSet = nextSeqSet;
//		}
//		nextSeqSet.clear();
//	}
//	clear(active);
//	fill(active, length(mat), false);
//	typedef typename Iterator<TConservedBlock>::Type TConsBlockIter;
//	
//	// Check quality of blocks
//	TConsBlockIter itCB = begin(conservedBlock);
//	TConsBlockIter itCBEnd = end(conservedBlock);
//
//	//typedef std::set<std::pair<TScoreValue, TConsBlockIter> > TBlockQuality;
//	//TBlockQuality bQual;
//	//for(;itCB != itCBEnd; ++itCB) {
//	//	TScoreValue totalScore = 0;
//	//	for(TSize column = (value(itCB)).i2; column < (value(itCB)).i3; ++column) {
//	//		typename TSequenceSet::const_iterator itCons = (value(itCB)).i1.begin();
//	//		typename TSequenceSet::const_iterator itConsEnd = (value(itCB)).i1.end();
//	//		for(;itCons != itConsEnd; ++itCons) {
//	//			typename TSequenceSet::const_iterator pairItCons = itCons;
//	//			++pairItCons;
//	//			for(;pairItCons != itConsEnd; ++pairItCons) {
//	//				totalScore += score(const_cast<TScore&>(score_type), value(mat, *itCons * len + column), value(mat, *pairItCons * len + column));
//	//			}
//	//		}
//	//	}
//	//	totalScore /= (TScoreValue) (value(itCB)).i1.size();
//	//	bQual.insert(std::make_pair(totalScore, itCB));
//	//	//if (totalScore < 0) (value(itCB)).i1.clear();
//	//}
//	//// Outlier detection
//	//TSize n = length(conservedBlock);
//	//typename TBlockQuality::const_iterator scrIt = bQual.begin();
//	//typename TBlockQuality::const_iterator scrItEnd = bQual.end();
//	//TScoreValue upperQuartile = 0;
//	//TScoreValue lowerQuartile = 0;
//	//TSize count = 0;
//	//for(;scrIt != scrItEnd; ++scrIt, ++count) {
//	//	if (n / 4 == count) lowerQuartile = scrIt->first;
//	//	else if ( (3 * n) / 4 == count) upperQuartile = scrIt->first;
//	//}
//	//TScoreValue interQR = upperQuartile - lowerQuartile;
//	//TScoreValue lowerInnerFence = (TScoreValue) ((double) lowerQuartile - 0.5 * (double) interQR);
//	//scrIt = bQual.begin();
//	//scrItEnd = bQual.end();
//	//for(;((scrIt != scrItEnd) && (scrIt->first < lowerInnerFence)); ++scrIt) {
//	//	(value(scrIt->second)).i1.clear();			
//	//}
//
//	//// Check blocks for outliers
//	//typedef std::pair<TScoreValue, TSize> TQualitySeqPair;
//	//typedef std::set<TQualitySeqPair, std::greater<TQualitySeqPair> > TSeqQuality;
//	//itCB = begin(conservedBlock);
//	//itCBEnd = end(conservedBlock);
//	//for(;itCB != itCBEnd; ++itCB) {
//	//	TSeqQuality seqQuality;
//	//	typename TSequenceSet::const_iterator leaveOut = (value(itCB)).i1.begin();
//	//	typename TSequenceSet::const_iterator leaveOutEnd = (value(itCB)).i1.end();
//	//	for(;leaveOut != leaveOutEnd; ++leaveOut) {
//	//		TScoreValue totalScore = 0;
//	//		for(TSize column = (value(itCB)).i2; column < (value(itCB)).i3; ++column) {
//	//			typename TSequenceSet::const_iterator itCons = (value(itCB)).i1.begin();
//	//			typename TSequenceSet::const_iterator itConsEnd = (value(itCB)).i1.end();
//	//			for(;itCons != itConsEnd; ++itCons) {
//	//				if (*itCons == *leaveOut) continue;
//	//				typename TSequenceSet::const_iterator pairItCons = itCons;
//	//				++pairItCons;
//	//				for(;pairItCons != itConsEnd; ++pairItCons) {
//	//					if (*pairItCons == *leaveOut) continue;
//	//					totalScore += score(const_cast<TScore&>(score_type), value(mat, *itCons * len + column), value(mat, *pairItCons * len + column));
//	//				}
//	//			}
//	//		}
//	//		seqQuality.insert(std::make_pair(totalScore, *leaveOut));
//	//	}
//	//	// Outlier detection
//	//	TSize n = (value(itCB)).i1.size();
//	//	typename TSeqQuality::const_iterator screenIt = seqQuality.begin();
//	//	typename TSeqQuality::const_iterator screenItEnd = seqQuality.end();
//	//	TScoreValue upperQuartile = 0;
//	//	TScoreValue lowerQuartile = 0;
//	//	TSize count = 0;
//	//	for(;screenIt != screenItEnd; ++screenIt, ++count) {
//	//		if (n / 4 == count) upperQuartile = screenIt->first;
//	//		else if ( (3 * n) / 4 == count) lowerQuartile = screenIt->first;
//	//	}
//	//	TScoreValue interQR = upperQuartile - lowerQuartile;
//	//	TScoreValue upperInnerFence = upperQuartile + 3 * interQR;
//	//	screenIt = seqQuality.begin();
//	//	screenItEnd = seqQuality.end();
//	//	for(;((screenIt != screenItEnd) && (screenIt->first > upperInnerFence)); ++screenIt) {
//	//		(value(itCB)).i1.erase(screenIt->second);			
//	//	}
//	//}
//
//	
//	
//	// Debug code: All blocks
//	itCB = begin(conservedBlock);
//	itCBEnd = end(conservedBlock);
//	for(;itCB != itCBEnd; ++itCB) {
//		typename TSequenceSet::const_iterator itCons = (value(itCB)).i1.begin();
//		typename TSequenceSet::const_iterator itConsEnd = (value(itCB)).i1.end();
//		for(;itCons != itConsEnd; ++itCons) {
//			for(TSize column = (value(itCB)).i2; column < (value(itCB)).i3; ++column) {
//				std::cout << value(mat, *itCons * len + column);
//			}
//			std::cout << std::endl;
//		}
//		std::cout << std::endl;
//	}
//
//
//	itCB = begin(conservedBlock);
//	itCBEnd = end(conservedBlock);
//	for(;itCB != itCBEnd; ++itCB) {
//		typename TSequenceSet::const_iterator itCons = (value(itCB)).i1.begin();
//		typename TSequenceSet::const_iterator itConsEnd = (value(itCB)).i1.end();
//		for(;itCons != itConsEnd; ++itCons) {
//			for(TSize column = (value(itCB)).i2; column < (value(itCB)).i3; ++column) {
//				value(active, *itCons * len + column) = true;
//			}
//		}
//	}
//
//
//	// Compute a distance matrix from the anchor graph
//	clear(distMatrix);
//	resize(distMatrix, nseq * nseq);
//	for(TSize seq1 = 0; seq1 < nseq-1; ++seq1) {
//		for(TSize seq2 = seq1+1; seq2 < nseq; ++seq2) {
//			typename Value<TDistanceMatrix>::Type diff = 0;
//			for(TSize column = 0; column < len; ++column) {
//				if ((value(active, seq1 * len + column) != value(active, seq2 * len + column)) ||
//					(value(active, seq1 * len + column) == 0))	++diff;
//			}
//			value(distMatrix, seq1 * nseq + seq2) = diff;
//		}
//	}
//
//	// Create the anchor graph
//	typedef Fragment<> TFragment;
//	typedef String<TFragment> TFragmentString;
//	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
//	TFragmentString matches;
//	typedef std::pair<TSize, TSize> TResiduePair;
//	typedef std::set<TResiduePair> TResiduePairSet;
//	String<TResiduePairSet> resPair;
//	resize(resPair, nseq * nseq);	
//	for(TSize seq1 = 0; seq1 < nseq - 1; ++seq1) {
//		for(TSize seq2 = seq1 + 1; seq2 < nseq; ++seq2) {
//			TSize index = seq1 * nseq + seq2;
//			TSize offset1 = 0;
//			TSize offset2 = 0;
//			for(TSize col = 0; col<len; ++col) {
//				if (value(mat, seq1 * len + col) != gapChar) {
//					if (value(mat, seq2 * len + col) != gapChar) {
//						if ((value(active, seq1 * len + col) == true) && (value(active, seq2 * len + col) == true)) {
//							resPair[index].insert(std::make_pair(offset1, offset2));
//						}
//						++offset1;
//						++offset2;
//					} else ++offset1;
//				} else if (value(mat, seq2 * len + col) != gapChar) ++offset2;
//			}
//		}
//	}
//	for(TSize i = 0; i<length(resPair); ++i) {
//		if (resPair[i].empty()) continue;
//		TSize seq1 = i / nseq;
//		TSize seq2 = i % nseq;
//		typename TResiduePairSet::const_iterator pos = resPair[i].begin();
//		typename TResiduePairSet::const_iterator posEnd = resPair[i].end();
//		TSize startMatch1 = pos->first;
//		TSize startMatch2 = pos->second;
//		TSize len = 1;
//		++pos;
//		while(pos != posEnd) {
//			if ((startMatch1 + len == pos->first) && (startMatch2 + len == pos->second)) ++len;
//			else {
//				appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
//				startMatch1 = pos->first;
//				startMatch2 = pos->second;
//				len = 1;
//			}
//			++pos;
//		}
//		appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
//	}
//	clearVertices(anchorGraph);
//	matchRefinement(matches,stringSet(anchorGraph),const_cast<TScore&>(score_type),anchorGraph);
//}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...

