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
// T-Coffee - Main Routines
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TNames, typename TMemeFile, typename TAlignmentGraph>
inline void
tCoffeeProteinAlignment(StringSet<TString, Dependent<> > const& strSet,
						TNames& names,
						TMemeFile& matchfile,
						TAlignmentGraph& gOut)
{
	SEQAN_CHECKPOINT
	
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	TSize nSeq = length(strSet);
	TSize threshold = 30;
	
	// Select informative pairs
	TGraph g(strSet);
	String<Pair<TId, TId> > pList;
	selectPairsForLibraryGeneration(g, pList);

	// Generate a primary library, i.e., all global pairwise alignments
	TGraph lib1(strSet);
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
	TGraph lib2(strSet);
	TGraph lib3(strSet);
	if (length(matchfile)) {
		// Only a selected list of local alignments
		Blosum62 score_type_local(-2,-8);
		generatePrimaryLibrary(lib2, pListLocal, score_type_local, LocalPairwise_Library() );

		// Blast matches
		std::fstream strm_lib;
		strm_lib.open(toCString(matchfile), std::ios_base::in | std::ios_base::binary);
		read(strm_lib, lib3, names, BlastLib());	// Read library
		strm_lib.close();
	} else {
		Blosum62 score_type_local(-2,-8);
		generatePrimaryLibrary(lib2, pList, score_type_local, LocalPairwise_Library() );
	}


	TGraph lib4(strSet);
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
	combineGraphs(g, true, libs);

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

template<typename TString, typename TNames, typename TAlignmentGraph>
inline void
tCoffeeProteinAlignment(StringSet<TString, Dependent<> > const& strSet,
						TNames& names,
						TAlignmentGraph& gOut)
{
	SEQAN_CHECKPOINT
	tCoffeeProteinAlignment(strSet, names, "", gOut);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TAlignmentGraph>
inline void
tCoffeeProteinAlignment(StringSet<TString, Dependent<> > const& strSet,
						TAlignmentGraph& gOut)
{
	SEQAN_CHECKPOINT
	String<String<char> > names;
	tCoffeeProteinAlignment(strSet, names, gOut);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TAlignmentGraph>
inline void
tCoffeeDnaAlignment(StringSet<TString, Dependent<> > const& strSet,
					TAlignmentGraph& gOut)
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
	TGraph g(strSet);
	String<Pair<TId, TId> > pList;
	selectPairsForLibraryGeneration(g, pList);

	// Generate a primary library, i.e., all global pairwise alignments
	TGraph lib1(strSet);
	String<double> distanceMatrix;
	generatePrimaryLibrary(lib1, pList, distanceMatrix, score_type_global, GlobalPairwise_Library() );
	//std::cout << "Library size: " << numVertices(lib1) << " Vertices, " << numEdges(lib1) << " Edges" << std::endl;

	// Generate a primary library, i.e., all local pairwise alignments
	TGraph lib2(strSet);
	generatePrimaryLibrary(lib2, pList, score_type_local, LocalPairwise_Library() );
	//std::cout << "Library size: " << numVertices(lib2) << " Vertices, " << numEdges(lib2) << " Edges" << std::endl;

	// Weighting of libraries (Signal addition)
	String<TGraph*> libs;
	appendValue(libs, &lib1);
	appendValue(libs, &lib2);
	//combineGraphs(g, libs);
	//// Only topology
	combineGraphs(g, true, libs);
	//std::cout << "Library size: " << numVertices(g) << " Vertices, " << numEdges(g) << " Edges" << std::endl;

	// Clear the old libraries
	clear(lib1);
	clear(lib2);


	TSize nSeq = length(strSet);
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

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TNames, typename TFileName, typename TAlignmentGraph>
inline void
tCoffeeLongDnaAlignment(StringSet<TString, Dependent<> >& strSet,
						TNames& names,
						TFileName& libFile,
						TAlignmentGraph& gOut)
{
	SEQAN_CHECKPOINT
	
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;

	TGraph g(strSet);
	TGraph lib1(strSet);
	std::fstream strm_lib;
	strm_lib.open(toCString(libFile), std::ios_base::in | std::ios_base::binary);
	read(strm_lib, lib1, names, BlastLib());	// Read library
	strm_lib.close();

	// Generate a primary library, i.e., all local pairwise alignments
	TGraph lib2(strSet);
	generatePrimaryLibrary(lib2, Lcs_Library() );
	
	
	//String<Pair<TId, TId> > pList;
	//selectPairsForLibraryGeneration(g, pList);
	//Score<int> score_type_global = Score<int>(5,-4,-4,-14);
	//String<double> distanceMatrix;
	//generatePrimaryLibrary(lib1, pList, distanceMatrix, score_type_global, GlobalPairwise_Library() );
	//

	// Weighting of libraries (Signal addition)
	String<TGraph*> libs;
	appendValue(libs, &lib1);
	appendValue(libs, &lib2);
	combineGraphs(g, true, libs);



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

template<typename TPath, typename TFilePrefix, typename TFileSuffix, typename TNames, typename TStringSet>
inline unsigned int
_alignImportSequences(TPath const& in_path, 
					  TFilePrefix const& file_prefix, 
					  TFileSuffix const& file_suffix,
					  TStringSet& origStrSet,
					  TNames& names)
{
	SEQAN_CHECKPOINT

	// Count sequences and read names
	unsigned seqCount = 0;
	std::ifstream file;
	std::stringstream input;
	if (length(file_suffix)) input << in_path << file_prefix << '.' << file_suffix;
	else input << in_path << file_prefix;
	file.open(input.str().c_str(), std::ios_base::in | std::ios_base::binary);
	if (!file.is_open()) return 0;
	while (!_streamEOF(file)) {
		String<char> id;
		readID(file, id, Fasta());
		appendValue(names, id);
		goNext(file, Fasta());
		++seqCount;
	}

	// Import sequences
	file.clear();
	file.seekg(0, std::ios_base::beg);
	resize(origStrSet, seqCount);
	unsigned int count = 0;
	for(unsigned i = 0; (i < seqCount) && !_streamEOF(file); ++i) 	{
		read(file, origStrSet[i], Fasta());
		count += length(origStrSet[i]);
	}
    file.close();

	return count;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TText>
inline void
_alignTiming(std::clock_t& startTime,
			 TText const& text)
{
	std::clock_t endTime=clock();
	double time=((float)(endTime-startTime)/CLOCKS_PER_SEC);
	startTime = endTime;
	std::cout << text << time << " sec" << std::endl;
}



//////////////////////////////////////////////////////////////////////////////
// T-Coffee - Debug stuff
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
		std::cout << std::endl;	

		SEQAN_TASSERT(sequenceId(g,sourceVertex(it_tmp)) != sequenceId(g,targetVertex(it_tmp)))
		SEQAN_TASSERT(fragmentBegin(g,sourceVertex(it_tmp)) < length((stringSet(g))[idToPosition((stringSet(g)), id1)]))
		SEQAN_TASSERT(fragmentBegin(g,sourceVertex(it_tmp)) + fragmentLength(g,sourceVertex(it_tmp)) <= length((stringSet(g))[idToPosition((stringSet(g)), id1)]))
		SEQAN_TASSERT(fragmentBegin(g,targetVertex(it_tmp)) < length((stringSet(g))[idToPosition((stringSet(g)), id2)]))
		SEQAN_TASSERT(fragmentBegin(g,targetVertex(it_tmp)) + fragmentLength(g,targetVertex(it_tmp)) <= length((stringSet(g))[idToPosition((stringSet(g)), id2)]))
		SEQAN_TASSERT(fragmentLength(g,sourceVertex(it_tmp)) == fragmentLength(g,targetVertex(it_tmp)))

	}
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
