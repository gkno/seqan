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
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_MSA_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_MSA_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Multiple sequence alignment
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Tags
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.MSA_Protein:
	Multiple sequence alignment for protein sequences.
*/
struct MSA_Protein_;
typedef Tag<MSA_Protein_> const MSA_Protein;

/**
.Tag.Global Alignment Algorithms.value.MSA_Dna:
	Multiple sequence alignment for nucleotide sequences.
*/
struct MSA_Dna_;
typedef Tag<MSA_Dna_> const MSA_Dna;


/**
.Tag.Global Alignment Algorithms.value.MSA_Genome:
	Multiple sequence alignment for genomic, closely related sequences.
*/
struct MSA_Genome_;
typedef Tag<MSA_Genome_> const MSA_Genome;


//////////////////////////////////////////////////////////////////////////////
// Global Alignment functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TFileName, typename TAlignmentGraph>
inline void
globalAlignment(StringSet<TString, TSpec> const& seqSet,
				StringSet<TName, TSpec2> const& nameSet,
				TFileName& matchfile,
				TAlignmentGraph& gOut,
				MSA_Protein)
{
	SEQAN_CHECKPOINT

	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	TSize nSeq = length(seqSet);
	TSize threshold = 30;
	
	// Select informative pairs
	TGraph g(seqSet);
	String<Pair<TId, TId> > pList;
	selectPairsForLibraryGeneration(g, pList);

	// Generate a primary library, i.e., all global pairwise alignments
	TGraph lib1(seqSet);
	String<double> distanceMatrix;
	Blosum62 score_type_global(-1,-11);
	generatePrimaryLibrary(lib1, pList, distanceMatrix, score_type_global, GlobalPairwise_Library() );
	
	// Select pairs for end-gaps free and local alignments?
	typedef std::multimap<double, Pair<TId, TId> > TBestPairs;
	TBestPairs bestPairs;
	double dist = 0;
	for(TSize i=0;i<nSeq-1;++i) {
		for(TSize j=i+1;j<nSeq;++j) {
			double d = getValue(distanceMatrix, i*nSeq+j);
			bestPairs.insert(std::make_pair(d, Pair<TId, TId>(positionToId(stringSet(g), i),positionToId(stringSet(g), j)))); 
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

	// Generate a third primary library, i.e., all local pairwise alignments and external matches
	TGraph lib2(seqSet);
	TGraph lib3(seqSet);
	if (length(matchfile)) {
		// Only a selected list of local alignments
		Blosum62 score_type_local(-2,-8);
		generatePrimaryLibrary(lib2, pListLocal, score_type_local, LocalPairwise_Library() );

		// Blast matches
		std::fstream strm_lib;
		strm_lib.open(toCString(matchfile), std::ios_base::in | std::ios_base::binary);
		read(strm_lib, lib3, nameSet, BlastLib());	// Read library
		strm_lib.close();
	} else {
		Blosum62 score_type_local(-2,-8);
		generatePrimaryLibrary(lib2, pList, score_type_local, LocalPairwise_Library() );
	}


	TGraph lib4(seqSet);
	if (dist > 0.75) {
		Blosum30 sT(-4,-20);
		AlignConfig<true,true,true,true> ac;
		generatePrimaryLibrary(lib4, pListLocal, sT, ac, GlobalPairwise_Library() );
	} else if (dist < 0.50) {
		Blosum80 sT(-2, -12);
		AlignConfig<true,true,true,true> ac;
		generatePrimaryLibrary(lib4, pListLocal, sT, ac, GlobalPairwise_Library() );
	} else {
		Blosum62 sT(-3,-14);
		AlignConfig<true,true,true,true> ac;
		generatePrimaryLibrary(lib4, pListLocal, sT, ac, GlobalPairwise_Library() );
	}


	// Weighting of libraries (Signal addition)
	String<TGraph*> libs;
	appendValue(libs, &lib1);
	appendValue(libs, &lib2);
	if (!empty(lib3)) appendValue(libs, &lib3);
	appendValue(libs, &lib4);
	combineGraphs(g, libs, FrequencyCounting() );
	//combineGraphs(g, libs);

	// Clear the old libraries
	clear(lib1);
	clear(lib2);
	clear(lib3);
	clear(lib4);

	if (nSeq < threshold) {
		// Full triplet...
		tripletLibraryExtension(g);

		// ... and normal progressive alignment with guide tree
		Graph<Tree<double> > guideTree;
		//upgmaTree(distanceMatrix, guideTree, UpgmaMin() );	
		slowNjTree(distanceMatrix, guideTree);
		progressiveAlignment(g, guideTree, gOut);
		clear(guideTree);
	} else {
		// Triplet only on groups of sequences
		Graph<Tree<double> > guideTree;
		//upgmaTree(distanceMatrix, guideTree);	
		slowNjTree(distanceMatrix, guideTree);	// More balanced than UPGMA
		progressiveAlignment(g, guideTree, gOut, threshold);
		clear(guideTree);
	}
	clear(distanceMatrix);
	clear(g);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TFileName, typename TAlignmentGraph>
inline void
globalAlignment(StringSet<TString, TSpec> const& seqSet,
				StringSet<TName, TSpec2> const& nameSet,
				TFileName const& matchfile,
				TAlignmentGraph& gOut,
				MSA_Genome)
{
	SEQAN_CHECKPOINT
	
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;

	TGraph g(seqSet);
	TGraph lib1(seqSet);
	if (length(matchfile)) {
		std::fstream strm_lib;
		strm_lib.open(toCString(matchfile), std::ios_base::in | std::ios_base::binary);
		read(strm_lib, lib1, nameSet, BlastLib());	// Read library
		strm_lib.close();
	}

	// Generate a primary library, i.e., all local pairwise alignments
	TGraph lib2(seqSet);
	Score<int> score_type = Score<int>(5,-4,-4,-14);
	generatePrimaryLibrary(lib2, score_type, Lcs_Library() );
	
	
	//String<Pair<TId, TId> > pList;
	//selectPairsForLibraryGeneration(g, pList);
	//Score<int> score_type_global = Score<int>(5,-4,-4,-14);
	//String<double> distanceMatrix;
	//generatePrimaryLibrary(lib1, pList, distanceMatrix, score_type_global, GlobalPairwise_Library() );
	//

	// Weighting of libraries (Signal addition)
	String<TGraph*> libs;
	if (!empty(lib1)) appendValue(libs, &lib1);
	appendValue(libs, &lib2);
	combineGraphs(g, libs, FrequencyCounting() );



	// Calculate a distance matrix using kmer counts
	String<double> distanceMatrix;
	getDistanceMatrix(g, distanceMatrix, KmerDistance());




	// Full triplet...
	tripletLibraryExtension(g);

	// ... and normal progressive alignment with guide tree
	Graph<Tree<double> > guideTree;
	//upgmaTree(distanceMatrix, guideTree);	
	slowNjTree(distanceMatrix, guideTree);
	progressiveAlignment(g, guideTree, gOut);
	clear(guideTree);
	clear(distanceMatrix);
	clear(g);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TFileName, typename TAlignmentGraph>
inline void
globalAlignment(StringSet<TString, TSpec> const& seqSet,
				StringSet<TName, TSpec2> const& nameSet,
				TFileName const& matchfile,
				TAlignmentGraph& gOut,
				MSA_Dna)
{
	SEQAN_CHECKPOINT
	
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;

	// Score objects
	Score<int> score_type_global = Score<int>(5,-4,-4,-14);
	Score<int> score_type_local = Score<int>(5,-4,-4,-14);

	// Select pairs
	TGraph g(seqSet);
	String<Pair<TId, TId> > pList;
	selectPairsForLibraryGeneration(g, pList);

	// Generate a primary library, i.e., all global pairwise alignments
	TGraph lib1(seqSet);
	String<double> distanceMatrix;
	generatePrimaryLibrary(lib1, pList, distanceMatrix, score_type_global, GlobalPairwise_Library() );
	//std::cout << "Library size: " << numVertices(lib1) << " Vertices, " << numEdges(lib1) << " Edges" << std::endl;

	// Generate a primary library, i.e., all local pairwise alignments
	TGraph lib2(seqSet);
	generatePrimaryLibrary(lib2, pList, score_type_local, LocalPairwise_Library() );
	//std::cout << "Library size: " << numVertices(lib2) << " Vertices, " << numEdges(lib2) << " Edges" << std::endl;

	TGraph lib3(seqSet);
	if (length(matchfile)) {
		std::fstream strm_lib;
		strm_lib.open(toCString(matchfile), std::ios_base::in | std::ios_base::binary);
		read(strm_lib, lib3, nameSet, BlastLib());	// Read library
		strm_lib.close();
	}

	// Weighting of libraries (Signal addition)
	String<TGraph*> libs;
	appendValue(libs, &lib1);
	appendValue(libs, &lib2);
	if (!empty(lib3)) appendValue(libs, &lib3);
	//combineGraphs(g, libs);
	//// Only topology
	combineGraphs(g, libs, FrequencyCounting() );
	//std::cout << "Library size: " << numVertices(g) << " Vertices, " << numEdges(g) << " Edges" << std::endl;

	// Clear the old libraries
	clear(lib1);
	clear(lib2);
	clear(lib3);


	TSize nSeq = length(seqSet);
	TSize threshold = 30;
	if (nSeq < threshold) {
		// Full triplet...
		tripletLibraryExtension(g);

		// ... and normal progressive alignment with guide tree
		Graph<Tree<double> > guideTree;
		//upgmaTree(distanceMatrix, guideTree);	
		slowNjTree(distanceMatrix, guideTree);
		progressiveAlignment(g, guideTree, gOut);
		clear(guideTree);
	} else {
		// Triplet only on groups of sequences
		Graph<Tree<double> > guideTree;
		//upgmaTree(distanceMatrix, guideTree);	
		slowNjTree(distanceMatrix, guideTree);	// More balanced than UPGMA
		progressiveAlignment(g, guideTree, gOut, threshold);
		clear(guideTree);
	}
	clear(distanceMatrix);
	clear(g);
}



template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TAlignmentGraph, typename TTag>
inline void
globalAlignment(StringSet<TString, TSpec> const& seqSet,
				StringSet<TName, TSpec2> const& nameSet,
				TAlignmentGraph& gOut,
				TTag)
{
	SEQAN_CHECKPOINT
	globalAlignment(seqSet, nameSet, "", gOut, TTag() );
}


template<typename TString, typename TSpec, typename TDependentSequenceSet, typename TCargo, typename TSpec2, typename TTag>
inline void
globalAlignment(StringSet<TString, TSpec> const& seqSet,
				Graph<Alignment<TDependentSequenceSet, TCargo, TSpec2> >& gOut,
				TTag)
{
	SEQAN_CHECKPOINT
	StringSet<String<char> > names;
	globalAlignment(seqSet, names, gOut, TTag() );
}



















// Testing Stuff



template<typename TStringSet, typename TCargo, typename TSpec, typename TString, typename TMatches, typename TSize>
inline void 
__alignSeqFragments(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					TString& unaligned1,
					TString& unaligned2,
					TMatches& matches,
					TSize seq1,
					TSize seq2,
					TSize off1,
					TSize off2)
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TId;
	typedef typename Value<TMatches>::Type TFragment;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	TStringSet& seqSet = stringSet(g);
	if ((length(unaligned1) > 0) && (length(unaligned2) > 0)) {
		// Make a pairwise string-set
		TStringSet pairSet;
		TId id1 = positionToId(seqSet, seq1);
		TId id2 = positionToId(seqSet, seq2);
		assignValueById(pairSet, unaligned1, id1);
		assignValueById(pairSet, unaligned2, id2);	
		TGraph pairAlign(pairSet);
		
		Blosum62 st1(-1,-11);
		TSize from = length(matches);
		globalAlignment(matches, pairSet, st1, Gotoh() );
		TSize to = length(matches);
		
		for(TSize iter = from; iter < to; ++iter) {
			fragmentBegin(matches[iter], id1) += off1;
			fragmentBegin(matches[iter], id2) += off2;
		}
	}
}



template<typename TStringSet, typename TCargo, typename TSpec, typename TScore, typename TConsensus>
inline void
__getAnchorGraph(Graph<Alignment<TStringSet, TCargo, TSpec> >& anchorGraph,
					TScore const& score_type,
					TConsensus& consensus)
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Size<TGraph>::Type TSize;
	TSize nSeq = length(stringSet(anchorGraph));
	

	// String of fragments to combine all pairwise alignments into a multiple alignment
	typedef Fragment<> TFragment;
	typedef String<TFragment> TFragmentString;
	TFragmentString matches;

	TStringSet consSet = stringSet(anchorGraph);
	appendValue(consSet, consensus);
	TGraph tmpGraph(consSet);
	TId idCons = positionToId(consSet, nSeq);
	TGraph bla(consSet);
	for(TSize seq = 0; seq < nSeq; ++seq) {
		TStringSet pairSet;
		TId id1 = positionToId(consSet, seq);
		assignValueById(pairSet, consSet, id1);
		assignValueById(pairSet, consSet, idCons);
		
		TGraph bla(pairSet);
		Blosum62 st(-1,-20);
		globalAlignment(matches, pairSet, st, AlignConfig<true, false, false, true>(), Gotoh() );	
		globalAlignment(bla, st, AlignConfig<true, false, false, true>(), Gotoh() );	
		std::cout << bla << std::endl;
		clearVertices(bla);
		globalAlignment(bla, st, AlignConfig<false, true, true, false>(), Gotoh() );	
		std::cout << bla << std::endl;

	}

	// Refine all matches and create multiple alignment
	matchRefinement(matches,consSet,const_cast<TScore&>(score_type),tmpGraph);

	String<double> initialDistanceMatrix;
	getDistanceMatrix(tmpGraph, initialDistanceMatrix, 6, AAGroupsDayhoff(), KmerDistance() );
	Graph<Tree<double> > guideTree;
	slowNjTree(initialDistanceMatrix, guideTree);

	tripletLibraryExtension(tmpGraph);
	TGraph gOut(consSet);
	progressiveAlignment(tmpGraph, guideTree, gOut);

	clearVertices(anchorGraph);
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	TVertexIterator itV(gOut);
	for(;!atEnd(itV);++itV) {
		if (sequenceId(gOut, *itV) != idCons) addVertex(anchorGraph, sequenceId(gOut, *itV), fragmentBegin(gOut, *itV), fragmentLength(gOut,*itV));
	}
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TEdgeIterator itE(gOut);
	for(;!atEnd(itE);++itE) {
		TVertexDescriptor sV = sourceVertex(itE);
		TVertexDescriptor tV = targetVertex(itE);
		if ((sequenceId(gOut, sV) != idCons) && (sequenceId(gOut, tV) != idCons)) {
			addEdge(anchorGraph, findVertex(anchorGraph, sequenceId(gOut, sV), fragmentBegin(gOut, sV)), findVertex(anchorGraph, sequenceId(gOut, tV), fragmentBegin(gOut, tV)), 1);
		}
	}
}

template<typename TSet1, typename TSet2, typename TSet3>
inline void
getCommonSubset(TSet1 const& set1,
				TSet2 const& set2,
				TSet3& common)
{
	typename TSet1::const_iterator pos1 = set1.begin();
	typename TSet1::const_iterator pos1End = set1.end();
	for(;pos1 != pos1End; ++pos1) {
		if (set2.find(*pos1) != set2.end()) common.insert(*pos1);
	}
}




template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TFileName, typename TAlignmentGraph>
inline void
testProfiles(StringSet<TString, TSpec> const& seqSet,
			 StringSet<TName, TSpec2> const& nameSet,
			 TFileName& matchfile,
			 TAlignmentGraph& gOut)
{
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	TSize nSeq = length(seqSet);
	TSize threshold = 30;

	// Generate an anchor graph from local alignments
	TGraph lib1(seqSet);
	Blosum62 score_type_global(-1,-11);
	generatePrimaryLibrary(lib1, score_type_global, LocalExtendPairwise_Library() );

	// Build the anchor graph but keep all local alignments
	TGraph tmpGraph(lib1);
	cliqueReduction(tmpGraph);
	String<double> distanceMatrix;
	getDistanceMatrix(tmpGraph,distanceMatrix,LibraryDistance() );
	Graph<Tree<double> > guideTree;
	slowNjTree(distanceMatrix, guideTree);
	TGraph anchorGraph(seqSet);
	progressiveAlignment(tmpGraph, guideTree, anchorGraph);
	clear(guideTree); clear(tmpGraph);
	reduceGraph(anchorGraph, score_type_global, distanceMatrix);
	std::cout << anchorGraph << std::endl;

	// Writing
	std::fstream strmBlast;
	strmBlast.open(toCString(matchfile), std::ios_base::out | std::ios_base::trunc);
	write(strmBlast, anchorGraph, nameSet, BlastLib());
	strmBlast.close();

	// Build the new guide tree
	slowNjTree(distanceMatrix, guideTree);
	
	// Make a banded alignment
	TGraph refinedLib(seqSet);
	String<char> mat;
	char gapChar = gapValue<char>();
	convertAlignment(anchorGraph, mat);
	TSize len = length(mat) / nSeq;
	typedef Fragment<> TFragment;
	typedef String<TFragment> TFragmentString;
	TFragmentString matches;
	for(TSize seq1 = 0; seq1 < nSeq - 1; ++seq1) {
		for(TSize seq2 = seq1 + 1; seq2 < nSeq; ++seq2) {
			//std::cout << '#' << seq1 << ',' << seq2 << std::endl;
			TSize offset1 = 0;
			TSize offset2 = 0;
			TString unaligned1;
			TString unaligned2;
			for(TSize col = 0; col<len; ++col) {
				if (value(mat, seq1 * len + col) != gapChar) {
					if (value(mat, seq2 * len + col) != gapChar) {
						__alignSeqFragments(refinedLib, unaligned1, unaligned2, matches, seq1, seq2, (TSize) (offset1 - length(unaligned1)), (TSize) (offset2 - length(unaligned2)));
						clear(unaligned1);
						clear(unaligned2);
						++offset1;
						++offset2;
					} else {
						++offset1;
						appendValue(unaligned1, value(mat, seq1 * len + col)); 
					}
				} else if (value(mat, seq2 * len + col) != gapChar) {
					++offset2;
					appendValue(unaligned2, value(mat, seq2 * len + col)); 
				}
			}
			__alignSeqFragments(refinedLib, unaligned1, unaligned2, matches, seq1, seq2, (TSize) (offset1 - length(unaligned1)), (TSize) (offset2 - length(unaligned2)));
			clear(unaligned1);
			clear(unaligned2);
		}
	}

		
	// Refine all matches and create multiple alignment
	matchRefinement(matches,stringSet(refinedLib),refinedLib);

	// Weighting of libraries (Signal addition)
	TGraph  g(seqSet);
	String<TGraph*> libs;
	appendValue(libs, &anchorGraph);
	appendValue(libs, &anchorGraph);
	//appendValue(libs, &anchorGraph);
	appendValue(libs, &lib1);
	appendValue(libs, &refinedLib);
	combineGraphs(g, libs, FrequencyCounting() );
		
	// Clear the old libraries
	clear(anchorGraph);
	clear(refinedLib);
	clear(lib1);
	
	if (nSeq < threshold) {
		// Full triplet...
		tripletLibraryExtension(g);
		// ... and normal progressive alignment with guide tree
		progressiveAlignment(g, guideTree, gOut);
	} else {
		progressiveAlignment(g, guideTree, gOut, threshold);
	}
	//std::cout << gOut << std::endl;
	clear(guideTree);
	clear(distanceMatrix);
	clear(g);


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









//////////////////////////////////////////////////////////////////////////////////
//
//template<typename TString, typename TNames, typename TMemeFile, typename TAlignmentGraph>
//inline void
//tCoffeeProteinAlignment(StringSet<TString, Dependent<> > const& strSet,
//						TNames& names,
//						TMemeFile& matchfile,
//						TAlignmentGraph& gOut)
//{
//	typedef StringSet<TString, Dependent<> > TStringSet;
//	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
//	typedef typename Size<TGraph>::Type TSize;
//	typedef typename Id<TGraph>::Type TId;
//	TSize nSeq = length(strSet);
//	TSize threshold = 30;
//	
//	// Select informative pairs
//	TGraph g(strSet);
//	String<Pair<TId, TId> > pList;
//	selectPairsForLibraryGeneration(g, pList);
//
//	// Generate a primary library, i.e., all global pairwise alignments
//	TGraph lib1(strSet);
//	String<double> distanceMatrix;
//	Blosum62 score_type_global(-1,-11);
//	generatePrimaryLibrary(lib1, pList, distanceMatrix, score_type_global, GlobalPairwise_Library() );
//	
//	// Select pairs for end-gaps free and local alignments?
//	typedef std::multimap<double, Pair<TId, TId> > TBestPairs;
//	TBestPairs bestPairs;
//	double dist = 0;
//	for(TSize i=0;i<nSeq-1;++i) {
//		for(TSize j=i+1;j<nSeq;++j) {
//			double d = getValue(distanceMatrix, i*nSeq+j);
//			bestPairs.insert(std::make_pair(d, Pair<TId, TId>(positionToId(stringSet(g), i),positionToId(stringSet(g), j)))); 
//			dist+=d;
//		}
//	}
//	TSize numPairs = nSeq * (nSeq - 1) / 2;
//	dist /= numPairs;
//	String<Pair<TId, TId> > pListLocal;
//	pListLocal = pList;
//	if (nSeq > threshold) {
//		clear(pListLocal);
//		typename TBestPairs::reverse_iterator itBestR = bestPairs.rbegin();
//		typename TBestPairs::iterator itBestF = bestPairs.begin();
//		TSize limit = 4 * nSeq;
//		for(TSize i=0;i<limit;++i, ++itBestR, ++itBestF) {
//			appendValue(pListLocal, itBestF->second);
//			appendValue(pListLocal, itBestR->second);
//		}
//	}
//
//	// Generate a third primary library, i.e., all local pairwise alignments and external matches
//	TGraph lib2(strSet);
//	TGraph lib3(strSet);
//	if (length(matchfile)) {
//		// Only a selected list of local alignments
//		Blosum62 score_type_local(-2,-8);
//		generatePrimaryLibrary(lib2, pListLocal, score_type_local, LocalPairwise_Library() );
//
//		// Blast matches
//		std::fstream strm_lib;
//		strm_lib.open(toCString(matchfile), std::ios_base::in | std::ios_base::binary);
//		read(strm_lib, lib3, names, BlastLib());	// Read library
//		strm_lib.close();
//	} else {
//		Blosum62 score_type_local(-2,-8);
//		generatePrimaryLibrary(lib2, pList, score_type_local, LocalPairwise_Library() );
//	}
//
//
//	TGraph lib4(strSet);
//	if (dist > 0.75) {
//		Blosum30 sT(-4,-20);
//		AlignConfig<true,true,true,true> ac;
//		generatePrimaryLibrary(lib4, pListLocal, sT, ac, GlobalPairwise_Library() );
//	} else if (dist < 0.50) {
//		Blosum80 sT(-2, -12);
//		AlignConfig<true,true,true,true> ac;
//		generatePrimaryLibrary(lib4, pListLocal, sT, ac, GlobalPairwise_Library() );
//	} else {
//		Blosum62 sT(-3,-14);
//		AlignConfig<true,true,true,true> ac;
//		generatePrimaryLibrary(lib4, pListLocal, sT, ac, GlobalPairwise_Library() );
//	}
//
//
//	// Weighting of libraries (Signal addition)
//	String<TGraph*> libs;
//	appendValue(libs, &lib1);
//	appendValue(libs, &lib2);
//	if (!empty(lib3)) appendValue(libs, &lib3);
//	appendValue(libs, &lib4);
//	combineGraphs(g, libs, FrequencyCounting() );
//	//combineGraphs(g, libs);
//
//	// Clear the old libraries
//	clear(lib1);
//	clear(lib2);
//	clear(lib3);
//	clear(lib4);
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
//}
//
////////////////////////////////////////////////////////////////////////////////
//
//template<typename TString, typename TNames, typename TAlignmentGraph>
//inline void
//tCoffeeProteinAlignment(StringSet<TString, Dependent<> > const& strSet,
//						TNames& names,
//						TAlignmentGraph& gOut)
//{
//	tCoffeeProteinAlignment(strSet, names, "", gOut);
//}
//
//
////////////////////////////////////////////////////////////////////////////////
//
//template<typename TString, typename TAlignmentGraph>
//inline void
//tCoffeeProteinAlignment(StringSet<TString, Dependent<> > const& strSet,
//						TAlignmentGraph& gOut)
//{
//	String<String<char> > names;
//	tCoffeeProteinAlignment(strSet, names, gOut);
//}


////////////////////////////////////////////////////////////////////////////////
//
//template<typename TString, typename TNames, typename TFileName, typename TAlignmentGraph>
//inline void
//tCoffeeLongDnaAlignment(StringSet<TString, Dependent<> >& strSet,
//						TNames& names,
//						TFileName& libFile,
//						TAlignmentGraph& gOut)
//{
//	
//	typedef StringSet<TString, Dependent<> > TStringSet;
//	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
//	typedef typename Size<TGraph>::Type TSize;
//	typedef typename Id<TGraph>::Type TId;
//
//	TGraph g(strSet);
//	TGraph lib1(strSet);
//	if (length(libFile)) {
//		std::fstream strm_lib;
//		strm_lib.open(toCString(libFile), std::ios_base::in | std::ios_base::binary);
//		read(strm_lib, lib1, names, BlastLib());	// Read library
//		strm_lib.close();
//	}
//
//	// Generate a primary library, i.e., all local pairwise alignments
//	TGraph lib2(strSet);
//	Score<int> score_type = Score<int>(5,-4,-4,-14);
//	generatePrimaryLibrary(lib2, score_type, Lcs_Library() );
//	
//	
//	//String<Pair<TId, TId> > pList;
//	//selectPairsForLibraryGeneration(g, pList);
//	//Score<int> score_type_global = Score<int>(5,-4,-4,-14);
//	//String<double> distanceMatrix;
//	//generatePrimaryLibrary(lib1, pList, distanceMatrix, score_type_global, GlobalPairwise_Library() );
//	//
//
//	// Weighting of libraries (Signal addition)
//	String<TGraph*> libs;
//	if (!empty(lib1)) appendValue(libs, &lib1);
//	appendValue(libs, &lib2);
//	combineGraphs(g, libs, FrequencyCounting() );
//
//
//
//	// Calculate a distance matrix using kmer counts
//	String<double> distanceMatrix;
//	getDistanceMatrix(g, distanceMatrix, KmerDistance());
//
//
//
//
//	// Full triplet...
//	tripletLibraryExtension(g);
//
//	// ... and normal progressive alignment with guide tree
//	Graph<Tree<double> > guideTree;
//	//upgmaTree(distanceMatrix, guideTree);	
//	slowNjTree(distanceMatrix, guideTree);
//	progressiveAlignment(g, guideTree, gOut);
//	clear(guideTree);
//	clear(distanceMatrix);
//	clear(g);
//}

////////////////////////////////////////////////////////////////////////////////
//
//template<typename TString, typename TAlignmentGraph>
//inline void
//tCoffeeDnaAlignment(StringSet<TString, Dependent<> > const& strSet,
//					TAlignmentGraph& gOut)
//{
//	
//	typedef StringSet<TString, Dependent<> > TStringSet;
//	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
//	typedef typename Size<TGraph>::Type TSize;
//	typedef typename Id<TGraph>::Type TId;
//
//	// Score objects
//	Score<int> score_type_global = Score<int>(5,-4,-4,-14);
//	Score<int> score_type_local = Score<int>(5,-4,-4,-14);
//
//	// Select pairs
//	TGraph g(strSet);
//	String<Pair<TId, TId> > pList;
//	selectPairsForLibraryGeneration(g, pList);
//
//	// Generate a primary library, i.e., all global pairwise alignments
//	TGraph lib1(strSet);
//	String<double> distanceMatrix;
//	generatePrimaryLibrary(lib1, pList, distanceMatrix, score_type_global, GlobalPairwise_Library() );
//	//std::cout << "Library size: " << numVertices(lib1) << " Vertices, " << numEdges(lib1) << " Edges" << std::endl;
//
//	// Generate a primary library, i.e., all local pairwise alignments
//	TGraph lib2(strSet);
//	generatePrimaryLibrary(lib2, pList, score_type_local, LocalPairwise_Library() );
//	//std::cout << "Library size: " << numVertices(lib2) << " Vertices, " << numEdges(lib2) << " Edges" << std::endl;
//
//	// Weighting of libraries (Signal addition)
//	String<TGraph*> libs;
//	appendValue(libs, &lib1);
//	appendValue(libs, &lib2);
//	//combineGraphs(g, libs);
//	//// Only topology
//	combineGraphs(g, libs, FrequencyCounting() );
//	//std::cout << "Library size: " << numVertices(g) << " Vertices, " << numEdges(g) << " Edges" << std::endl;
//
//	// Clear the old libraries
//	clear(lib1);
//	clear(lib2);
//
//
//	TSize nSeq = length(strSet);
//	TSize threshold = 30;
//	if (nSeq < threshold) {
//		// Full triplet...
//		tripletLibraryExtension(g);
//
//		// ... and normal progressive alignment with guide tree
//		Graph<Tree<double> > guideTree;
//		//upgmaTree(distanceMatrix, guideTree);	
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
//}


	




	//// Compute the consensus sequence
	//typedef String<String<unsigned int> > TProfile;
	//TProfile profile;
	//unsigned int alphSize = ValueSize<TAlphabet>::VALUE;
	//fill(profile, len, String<unsigned int>());
	//typedef typename Iterator<TProfile>::Type TProfileIter;
	//TProfileIter itProf = begin(profile);
	//TProfileIter itProfEnd = end(profile);
	//for(;itProf!=itProfEnd;++itProf) fill(*itProf, alphSize, 0);
	//itA = begin(mat);
	//itB = begin(blocks);
	//col = 0;
	//for(; itA != itAEnd; ++itA, ++itB, col = (++col) % len) {
	//	if (value(itB) == true) {
	//		if (*itA != gapChar) value(value(profile, col), (TSize) TAlphabet(*itA)) += 1;
	//	}
	//}
	//itProf = begin(profile);
	//itProfEnd = end(profile);
	//for(;itProf!=itProfEnd;++itProf) {
	//	unsigned int max = 0;
	//	TAlphabet maxLetter;
	//	for(TSize i = 0; i<alphSize; ++i) {
	//		if (value(*itProf, i) > max) {
	//			max = value(*itProf, i);
	//			maxLetter = TAlphabet(i);
	//		}
	//	}
	//	if (max > 0) appendValue(consensus, maxLetter);
	//}

	//__getAnchorGraph(anchorGraph, score_type, consensus);
	//clear(mat);
	//convertAlignment(anchorGraph, mat);
	//len = length(mat) / nseq;
	//clear(blocks);
	//fill(blocks, length(mat), true);
	//clear(coverage);
	//fill(coverage, len, 0);
	//itA = begin(mat);
	//itAEnd = end(mat);
	//itB = begin(blocks);
	//col = 0;
	//for(; itA != itAEnd; ++itA, ++itB, col = (++col) % len) {
	//	if (*itA == gapChar) value(itB) = false;
	//	else value(coverage, col) += 1;
	//}
	//// Clear blocks with low coverage
	//itA = begin(mat);
	//itB = begin(blocks);
	//col = 0;
	//for(; itA != itAEnd; ++itA, ++itB, col = (++col) % len) {
	//	if (value(coverage, col) < minCoverage) value(itB) = false;
	//}




	//// Exclude short blocks
	//TSize blockStart = len;
	//TSize blockEnd = 0;
	//typedef std::set<TSize> TSequenceSet;
	//typedef String<TSequenceSet> TConservedBlock;
	//TConservedBlock conservedBlock;
	//String<TSize> startEndPos;
	//TSequenceSet thisSeqSet;
	//TSequenceSet nextSeqSet;
	//for(TSize column = 0; column < len; ++column) {
	//	for(TSize seq = 0; seq < nseq; ++seq) {
	//		if (!value(blocks, seq * len + column)) continue;
	//		nextSeqSet.insert(seq);
	//	}
	//	if (thisSeqSet != nextSeqSet) {
	//		if (blockStart == len) {
	//			blockStart = column;
	//		} else {
	//			if ((column - blockStart) > minBlockLen) {
	//				appendValue(startEndPos, blockStart);
	//				appendValue(startEndPos, column);
	//				appendValue(conservedBlock, thisSeqSet);
	//			}
	//			if (nextSeqSet.empty()) blockStart = len;
	//			else blockStart = column;
	//		}
	//		thisSeqSet = nextSeqSet;
	//	}
	//	nextSeqSet.clear();
	//}
	








	//// Compute the column scores for all blocks (ignore the gaps)
	//String<TScoreValue> columnScores;
	//fill(columnScores, len, -10000);
	//for(TSize k=0;k<len; ++k) {
	//	TScoreValue totalScore = 0;
	//	bool foundPair = false;
	//	for(TSize i = 0; i<nseq-1; ++i) {
	//		if (value(blocks, i * len + k) == false) continue;
	//		for(TSize j=i+1; j<nseq; ++j) {
	//			if (value(blocks, j * len + k) == false) continue;
	//			foundPair = true;
	//			totalScore += score(const_cast<TScore&>(score_type), value(mat, i*len+k), value(mat, j*len + k));
	//		}
	//	}
	//	if (foundPair) value(columnScores, k) = totalScore;
	//}

	//// Slide a moving average filter over all columns
	//String<TScoreValue> columnScoresAvg;
	//fill(columnScoresAvg, len, 0);
	//for(TSize column=0;column<len; ++column) {
	//	TSize leftBoundary = 0;
	//	TSize rightBoundary = len - 1;
	//	if (column > windowSize / 2) leftBoundary = column - windowSize / 2;
	//	if (column + windowSize / 2 < len) rightBoundary = column + windowSize / 2;
	//	TScoreValue score = 0;
	//	for(TSize point = leftBoundary; point <= rightBoundary; ++point) {
	//		score += value(columnScores, point);
	//	}
	//	value(columnScoresAvg, column) = score / (TScoreValue) (rightBoundary - leftBoundary + 1);
	//}
	//clear(columnScores);
	//// Inactivate columns with negative column score
	//for(TSize column = 0; column < len; ++column) {
	//	if (value(columnScoresAvg, column) <= 0) {
	//		for(TSize seq = 0; seq < nseq; ++seq) {
	//			value(blocks, seq * len + column) = false;
	//		}
	//	}
	//}



//template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TFileName, typename TAlignmentGraph>
//inline void
//testCliques(StringSet<TString, TSpec> const& seqSet,
//			StringSet<TName, TSpec2> const& nameSet,
//			TFileName& matchfile,
//			TAlignmentGraph& gOut)
//{
//	typedef StringSet<TString, Dependent<> > TStringSet;
//	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
//	typedef typename Size<TGraph>::Type TSize;
//	typedef typename Id<TGraph>::Type TId;
//	
//
//	typedef Graph<Tree<double> > TGuideTree;
//	typedef std::set<TSize> TChildrenSet;
//	typedef typename Iterator<TGuideTree, DfsPreorder>::Type TDfsPreorderIterator;
//	typedef typename Iterator<TGuideTree, BfsIterator>::Type TBfsIterator;
//	
//	TGraph anchorGraph(seqSet);
//	String<double> distanceMatrix;
//	getDistanceMatrix(anchorGraph, distanceMatrix, 6, AAGroupsDayhoff(), KmerDistance() );
//	TGuideTree guideTree;
//	slowNjTree(distanceMatrix, guideTree);
//
//
//	// Build groups of sequences
//	TSize seqPerGroup = 10;
//	TGuideTree copy_tree(guideTree);
//	String<TChildrenSet> sequenceGroups;
//	String<TSize> sequenceGroupRoots;
//	while (numChildren(copy_tree, getRoot(copy_tree)) > 0) {
//		TBfsIterator bfsIt(copy_tree, getRoot(copy_tree));
//		for(;!atEnd(bfsIt);goNext(bfsIt)) {
//			TChildrenSet sequenceSet;
//			collectLeaves(copy_tree, *bfsIt, sequenceSet);
//			if ((TSize) sequenceSet.size() <= (TSize) seqPerGroup) {
//				TChildrenSet finalSequenceSet;
//				for(typename TChildrenSet::iterator pos = sequenceSet.begin(); pos != sequenceSet.end(); ++pos) {
//					finalSequenceSet.insert(positionToId(seqSet, *pos));
//				}
//				appendValue(sequenceGroups, finalSequenceSet);
//				// Find the proper root
//				TDfsPreorderIterator dfsItRootFinder(copy_tree, *bfsIt);
//				for(;!atEnd(dfsItRootFinder);goNext(dfsItRootFinder)) {
//					if (numChildren(copy_tree, *dfsItRootFinder) == 1) continue;
//					appendValue(sequenceGroupRoots, *dfsItRootFinder);
//					break;
//				}
//				// Cut the tree
//				if (getRoot(copy_tree) != *bfsIt) {
//					removeChild(copy_tree, parentVertex(guideTree, *bfsIt), *bfsIt);
//				} else {
//					removeAllChildren(copy_tree, *bfsIt);
//				}
//				break;
//			}
//		}
//	}
//	clear(copy_tree);
//
//	// Select all the possible pairs within a group
//	String<Pair<TId, TId> > pList;
//	for(TSize i = 0;i<length(sequenceGroups);++i) {
//		typename TChildrenSet::const_iterator pos1 = (sequenceGroups[i]).begin();
//		while (pos1 != (sequenceGroups[i]).end()) {
//			typename TChildrenSet::const_iterator pos2 = pos1;
//			++pos2;
//			while (pos2 != (sequenceGroups[i]).end()) {
//				appendValue(pList, Pair<TId, TId>(positionToId(seqSet, *pos1), positionToId(seqSet, *pos2)));
//				++pos2;
//			}
//			++pos1;
//		}
//	}
//
//	//// Debug code
//	//std::cout << guideTree << std::endl;
//	//for(TSize i = 0;i<length(sequenceGroups);++i) {
//	//	std::cout << "Sequences: ";
//	//	copy((sequenceGroups[i]).begin(), (sequenceGroups[i]).end(), std::ostream_iterator<int>(std::cout, " "));
//	//	std::cout << std::endl;
//	//}
//	//for(TSize i = 0;i<length(pList);++i) {
//	//	std::cout << '(' << pList[i].i1 << ',' << pList[i].i2 << ')' << ';';
//	//}
//	//std::cout << std::endl;
//
//	// Generate a primary library, i.e., all global pairwise alignments
//	TGraph lib1(seqSet);
//	Blosum62 score_type_global(-1,-14);
//	//generatePrimaryLibrary(lib1, distanceMatrix, score_type_global, AlignConfig<true,true,true,true>(), GlobalPairwise_Library() );
//	generatePrimaryLibrary(lib1, pList, score_type_global, AlignConfig<true,true,true,true>(), GlobalPairwise_Library() );
//	tripletLibraryExtension(lib1);
//	//findMaximumClique(lib1);
//	progressiveAlignment(lib1, guideTree, anchorGraph);
//	std::cout << anchorGraph << std::endl;
//
//
//
//	// Generate a intergroup pairlist
//	String<Pair<TId, TId> > pListInter;
//	for(TSize i = 0;i<length(sequenceGroups)-1;++i) {
//		for(TSize j = i+1;j<length(sequenceGroups);++j) {
//			typename TChildrenSet::const_iterator pos0 = (sequenceGroups[i]).begin();
//			while (pos0 != (sequenceGroups[i]).end()) {
//				typename TChildrenSet::const_iterator pos1 = (sequenceGroups[j]).begin();
//				while (pos1 != (sequenceGroups[j]).end()) {
//					appendValue(pListInter, Pair<TId, TId>(positionToId(seqSet, *pos0), positionToId(seqSet, *pos1)));
//					++pos1;
//				}
//				++pos0;
//			}
//		}
//	}
//
//	TGraph lib2(seqSet);
//	Blosum30 st1(-4,-20);
//	generatePrimaryLibrary(lib2, pListInter, st1, AlignConfig<true,true,true,true>(), GlobalPairwise_Library() );
//
//	TGraph lib3(seqSet);
//	Blosum62 st2(-1,-11);
//	generatePrimaryLibrary(lib3, pListInter, st2, GlobalPairwise_Library() );
//
//	
//	// Weighting of libraries (Signal addition)
//	TGraph g(seqSet);
//	String<TGraph*> libs;
//	appendValue(libs, &anchorGraph);
//	appendValue(libs, &lib2);
//	appendValue(libs, &lib3);
//	combineGraphs(g, libs, FrequencyCounting() );
//	//combineGraphs(g, libs);
//
//	// Clear the old libraries
//	clear(anchorGraph);
//	clear(lib1);
//	clear(lib2);
//	clear(lib3);
//
//	tripletLibraryExtension(g);
//
//	progressiveAlignment(g, guideTree, gOut);
//	clear(guideTree);
//	clear(distanceMatrix);
//	clear(g);
//
//	//std::fstream strmBlast;
//	//strmBlast.open(toCString(matchfile), std::ios_base::out | std::ios_base::trunc);
//	//write(strmBlast, anchorGraph, nameSet, BlastLib());
//	//strmBlast.close();
//
//
//
//
//
//
//}












}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...

