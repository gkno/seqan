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

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_BASE_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// T-Coffee - Tags
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Guide Tree Configurator:
..summary:A tag to configure the guide tree construction.
*/

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Guide Tree Configurator.UpgmaMin:
	Uses the min operation in the upgma algorithm
*/

struct UpgmaMin_;
typedef Tag<UpgmaMin_> const UpgmaMin;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Guide Tree Configurator.UpgmaMax:
	Uses the max operation in the upgma algorithm
*/

struct UpgmaMax_;
typedef Tag<UpgmaMax_> const UpgmaMax;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Guide Tree Configurator.UpgmaAvg:
	Uses the average operation in the upgma algorithm
*/

struct UpgmaAvg_;
typedef Tag<UpgmaAvg_> const UpgmaAvg;

/**
.Tag.Distance Calculation:
..summary:A tag to specify how to calculate distance matrices.
*/

/**
.Tag.Distance Calculation.value.LibraryDistance:
	Using the library itself and heaviest common subsequence to determine a distance matrix
*/
struct LibraryDistance_;
typedef Tag<LibraryDistance_> const LibraryDistance;


/**
.Tag.Distance Calculation.value.KmerDistance:
	Using a simple kmer count to determine a distance matrix
*/
struct KmerDistance_;
typedef Tag<KmerDistance_> const KmerDistance;


/**
.Tag.Library Generation:
..summary:A tag to specify how to generate a T-Coffee library.
*/


/**
.Tag.Library Generation.value.GlobalPairwise_Library:
	A primary library of global alignments.
*/

struct GlobalPairwise_Library_;
typedef Tag<GlobalPairwise_Library_> const GlobalPairwise_Library;


/**
.Tag.Library Generation.value.GlobalPairwise_Library:
	A primary library of local alignments.
*/

struct LocalPairwise_Library_;
typedef Tag<LocalPairwise_Library_> const LocalPairwise_Library;


/**
.Tag.Library Generation.value.Overlap_Library:
	A primary library of overlap alignments.
*/

struct Overlap_Library_;
typedef Tag<Overlap_Library_> const Overlap_Library;

/**
.Tag.Library Generation.value.Kmer_Library:
	A primary library of kmer alignments.
*/

struct Kmer_Library_;
typedef Tag<Kmer_Library_> const Kmer_Library;


/**
.Tag.Library Generation.value.MUMPairwise_Library:
	A primary library of mums.
*/

struct MUMPairwise_Library_;
typedef Tag<MUMPairwise_Library_> const MUMPairwise_Library;



//////////////////////////////////////////////////////////////////////////////
// T-Coffee - Library combination
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TLibraries> 
inline void
combineGraphs(Graph<Alignment<TStringSet, TCargo, TSpec> >& outGraph,
			  bool onlyTopology,
			  TLibraries& libs)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	
	// Clear out-library
	clearVertices(outGraph);
	TSize numLibs = length(libs);	// Number of libraries

	// All the matches with score values
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;
	String<TCargo, Block<> > score_values;

	// Max score and index start position for every library
	String<TCargo> max_scores;
	String<TCargo> index_start;
	resize(max_scores, numLibs);
	resize(index_start, numLibs);

	// Get all matches
	TSize count = 0;
	for(TSize i = 0; i<numLibs; ++i) {
		assignValue(index_start, i, count);
		TCargo maxCargoLib = 0;
		TGraph const& lib = *(getValue(libs, i));
		TEdgeIterator it(lib);
		for(;!atEnd(it);++it) {
			TCargo currentCargo = getCargo(*it);
			if (currentCargo > maxCargoLib) maxCargoLib = currentCargo;
			TVertexDescriptor sV = sourceVertex(it);
			TVertexDescriptor tV = targetVertex(it);
			push_back(matches, TFragment( (unsigned int) sequenceId(lib, sV), (unsigned int) fragmentBegin(lib,sV), (unsigned int) sequenceId(lib, tV),  (unsigned int)  fragmentBegin(lib,tV),  (unsigned int)  fragmentLength(lib,tV)));
			push_back(score_values, currentCargo);
			++count;
		}
		assignValue(max_scores, i, maxCargoLib);
	}

	// Match refinement
	TStringSet& str = stringSet(outGraph);	
	matchRefinement(matches,str,outGraph);  // Don't score matches!
	//// Debug code
	//TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	//for(TSize i = 0; i<length(str);++i) {
	//	TId seqId = positionToId(str, i);
	//	TSize j = 0;
	//	TSize len = length(str[i]);
	//	while(j<len) {
	//		TVertexDescriptor nextVertex = findVertex(outGraph, seqId, j);
	//		if (nextVertex == nilVertex) {
	//			std::cout << j << std::endl;
	//			std::cout << findVertex(outGraph, seqId, len - 1) << std::endl;
	//			std::cout << "Nil Vertex!!" << std::endl;
	//			exit(0);
	//		}
	//		j += fragmentLength(outGraph, nextVertex);
	//	}
	//}

	if (!onlyTopology) {
		// Adapt edge weights (fractional weights are used)
		count = 0;
		TSize currentLib = 0;
		double upperBound = 100;
		double scaling = upperBound / (double) getValue(max_scores, currentLib);
		TSize nextLibCounter = length(matches);
		if (currentLib < numLibs - 1) nextLibCounter = (TSize) getValue(index_start, currentLib+1);
		TFragmentStringIter endIt = end(matches);
		for(TFragmentStringIter it = begin(matches); it != endIt; ++it) {
			if (count >= nextLibCounter) {
				++currentLib;
				scaling = upperBound / (double) getValue(max_scores, currentLib);
				if (currentLib < numLibs - 1) nextLibCounter = (TSize) getValue(index_start, currentLib+1);
				else nextLibCounter = length(matches);
			}
			TId id1 = sequenceId(*it,0);
			TId id2 = sequenceId(*it,1);
			TSize pos1 = fragmentBegin(*it, id1);
			TSize pos2 = fragmentBegin(*it, id2);
			TSize end1 = pos1 + fragmentLength(*it, id1);
			while(pos1 < end1) {
				TVertexDescriptor p1 = findVertex(outGraph, id1, pos1);
				TVertexDescriptor p2 = findVertex(outGraph, id2, pos2);
				TEdgeDescriptor e = findEdge(outGraph, p1, p2);
				TSize fragLen = fragmentLength(outGraph, p1); 
				double newVal = (double) fragLen / (double) (end1 - pos1);
				newVal *= scaling;
				newVal *= (double) getValue(score_values, position(it));
				SEQAN_TASSERT(e != 0)  // Refined graph should have all edges (with weight 1)
				cargo(e) += (TCargo) newVal;
				pos1 += fragLen;
				pos2 += fragLen;
			}
			++count;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TLibraries> 
inline void
combineGraphs(Graph<Alignment<TStringSet, TCargo, TSpec> >& outGraph,
			  TLibraries& libs)
{
	SEQAN_CHECKPOINT
	combineGraphs(outGraph, false, libs);
}





//////////////////////////////////////////////////////////////////////////////
// T-Coffee - Triplet extension
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TCargoMap, typename TVertexString, typename TCargoString>
inline void 
_performTripletExtension(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
						 TCargoMap& newCargoMap,
						 TVertexString& edges_vertices,
						 TCargoString& edges_cargo)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	// Remember the old cargo
	resize(newCargoMap, getIdUpperBound(_getEdgeIdManager(g)), Exact());
	TEdgeIterator it(g);
	for(;!atEnd(it);++it) assignProperty(newCargoMap, *it, cargo(*it));

	// Iterate over all vertices
	for(TVertexIterator itVertex(g);!atEnd(itVertex);++itVertex) {
		TOutEdgeIterator outIt1(g, *itVertex);
		while (!atEnd(outIt1)) {
			TOutEdgeIterator outIt2 = outIt1;
			goNext(outIt2);
			// Consider always 2 neighbors
			while (!atEnd(outIt2)) {
				TVertexDescriptor tV1 = targetVertex(outIt1);
				TVertexDescriptor tV2 = targetVertex(outIt2);
				if (sequenceId(g, tV1) != sequenceId(g,tV2)) {
					TEdgeDescriptor e = findEdge(g, tV1, tV2);
					if (e == 0) {
						// New edge
						TCargo val = cargo(*outIt1);
						if (val > cargo(*outIt2)) val = cargo(*outIt2);
						
						// Remember the edge with cargo
						push_back(edges_vertices, tV1);
						push_back(edges_vertices, tV2);
						push_back(edges_cargo, val);
					} else {
						// Increase weight of existing edge
						if (getCargo(*outIt2) > getCargo(*outIt1)) property(newCargoMap, e) += getCargo(*outIt1);
						else property(newCargoMap, e) += getCargo(*outIt2);	
					}
				}
				goNext(outIt2);
			}
			goNext(outIt1);
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TSet, typename TCargoMap, typename TVertexString, typename TCargoString>
inline void 
_performTripletExtension(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
						 TSet const& set1,
						 TCargoMap& newCargoMap,
						 TVertexString& edges_vertices,
						 TCargoString& edges_cargo)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	resize(newCargoMap, getIdUpperBound(_getEdgeIdManager(g)), Exact());
	TEdgeIterator it(g);
	for(;!atEnd(it);++it) assignProperty(newCargoMap, *it, cargo(*it));

	for(TVertexIterator itVertex(g);!atEnd(itVertex);++itVertex) {
		TOutEdgeIterator outIt1(g, *itVertex);
		for(;!atEnd(outIt1);goNext(outIt1)) {
			if (set1.find(sequenceId(g,targetVertex(outIt1))) == set1.end()) continue;
			TOutEdgeIterator outIt2 = outIt1;
			goNext(outIt2);
			for(;!atEnd(outIt2);goNext(outIt2)) {
				if (set1.find(sequenceId(g,targetVertex(outIt2))) == set1.end()) continue;
				TVertexDescriptor tV1 = targetVertex(outIt1);
				TVertexDescriptor tV2 = targetVertex(outIt2);
				if (sequenceId(g, tV1) == sequenceId(g,tV2)) continue;
				TEdgeDescriptor e = findEdge(g, tV1, tV2);
				if (e == 0) {
					// New edge
					TCargo val = cargo(*outIt1);
					if (val > cargo(*outIt2)) val = cargo(*outIt2);

					// Remember the edge with cargo
					push_back(edges_vertices, tV1);
					push_back(edges_vertices, tV2);
					push_back(edges_cargo, val);
				} else {
					// Increase weight of existing edge
					if (getCargo(*outIt2) > getCargo(*outIt1)) property(newCargoMap, e) += getCargo(*outIt1);
					else property(newCargoMap, e) += getCargo(*outIt2);	
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TSet, typename TCargoMap, typename TVertexString, typename TCargoString>
inline void 
_performTripletExtension(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
						 TSet const& set1,
						 TSet const& set2,
						 TCargoMap& newCargoMap,
						 TVertexString& edges_vertices,
						 TCargoString& edges_cargo)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	resize(newCargoMap, getIdUpperBound(_getEdgeIdManager(g)), Exact());
	TEdgeIterator it(g);
	for(;!atEnd(it);++it) assignProperty(newCargoMap, *it, cargo(*it));

	for(TVertexIterator itVertex(g);!atEnd(itVertex);++itVertex) {
		TOutEdgeIterator outIt1(g, *itVertex);
		for(;!atEnd(outIt1);goNext(outIt1)) {
			if ((set1.find(sequenceId(g,targetVertex(outIt1))) == set1.end()) && 
				(set2.find(sequenceId(g,targetVertex(outIt1))) == set2.end())) continue;
			TOutEdgeIterator outIt2 = outIt1;
			goNext(outIt2);
			for(;!atEnd(outIt2);goNext(outIt2)) {
				if ((set1.find(sequenceId(g,targetVertex(outIt2))) == set1.end()) && 
					(set2.find(sequenceId(g,targetVertex(outIt2))) == set2.end())) continue;
				TVertexDescriptor tV1 = targetVertex(outIt1);
				TVertexDescriptor tV2 = targetVertex(outIt2);
				if ((set1.find(sequenceId(g, tV1)) != set1.end()) &&
					(set1.find(sequenceId(g, tV2)) != set1.end())) continue;
				if ((set2.find(sequenceId(g, tV1)) != set2.end()) &&
					(set2.find(sequenceId(g, tV2)) != set2.end())) continue;
				TEdgeDescriptor e = findEdge(g, tV1, tV2);
				if (e == 0) {
					// New edge
					TCargo val = cargo(*outIt1);
					if (val > cargo(*outIt2)) val = cargo(*outIt2);

					// Remember the edge with cargo
					push_back(edges_vertices, tV1);
					push_back(edges_vertices, tV2);
					push_back(edges_cargo, val);
				} else {
					// Increase weight of existing edge
					if (getCargo(*outIt2) > getCargo(*outIt1)) property(newCargoMap, e) += getCargo(*outIt1);
					else property(newCargoMap, e) += getCargo(*outIt2);	
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TCargoMap>
inline void 
_assignNewCargos(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				 TCargoMap& newCargoMap)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;

	// Assign the new weights and clean-up the cargo map
	TEdgeIterator it(g);
	for(;!atEnd(it);++it) cargo(*it) = getProperty(newCargoMap, *it);
	clear(newCargoMap);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexString, typename TCargoString>
inline void 
_addNewEdgesFoundByTriplet(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
						   TVertexString& edges_vertices,
						   TCargoString& edges_cargo)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TVertexString>::Type TVertexStringIter;
	typedef typename Iterator<TCargoString>::Type TCargoStringIter;

	// Finally add the new edges created by the triplet approach
	TVertexStringIter itV = begin(edges_vertices);
	TVertexStringIter endIt = end(edges_vertices);
	TCargoStringIter itC = begin(edges_cargo);
	while(itV != endIt) {
		TVertexStringIter itVNext = itV; ++itVNext;
		// The same edge could have been created multiple times, so check if it exists
		TEdgeDescriptor e = findEdge(g, *itV, *itVNext);
		if (e == 0) addEdge(g, *itV, *itVNext, *itC);
		else cargo(e) += *itC;
		++itV; ++itV;
		++itC;
	}
	SEQAN_TASSERT(itC == end(edges_cargo))
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void 
tripletLibraryExtension(Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	// Two tasks:
	// 1) Add edges for the case that a and c is aligned, b and c is aligned, but a and b are not, give these edges the appropriate weight
	// 2) Augment all existing edges
	String<TCargo> newCargoMap;
	typedef String<TVertexDescriptor, Block<> > TVertexString;
	typedef String<TCargo, Block<> > TCargoString;
	TVertexString edges_vertices;
	TCargoString edges_cargo;
	
	_performTripletExtension(g, newCargoMap, edges_vertices, edges_cargo);
	_assignNewCargos(g, newCargoMap);
	_addNewEdgesFoundByTriplet(g, edges_vertices, edges_cargo);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TSet>
inline void 
tripletLibraryExtension(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
						TSet const& set1)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	// Two tasks:
	// 1) Add edges for the case that a and c is aligned, b and c is aligned, but a and b are not, give these edges the appropriate weight
	// 2) Augment all existing edges
	String<TCargo> newCargoMap;
	typedef String<TVertexDescriptor, Block<> > TVertexString;
	typedef String<TCargo, Block<> > TCargoString;
	TVertexString edges_vertices;
	TCargoString edges_cargo;


	_performTripletExtension(g, set1, newCargoMap, edges_vertices, edges_cargo);
	_assignNewCargos(g, newCargoMap);
	_addNewEdgesFoundByTriplet(g, edges_vertices, edges_cargo);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TSet>
inline void 
tripletLibraryExtension(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
						TSet const& set1,
						TSet const& set2)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	// Two tasks:
	// 1) Add edges for the case that a and c is aligned, b and c is aligned, but a and b are not, give these edges the appropriate weight
	// 2) Augment all existing edges
	String<TCargo> newCargoMap;
	typedef String<TVertexDescriptor, Block<> > TVertexString;
	typedef String<TCargo, Block<> > TCargoString;
	TVertexString edges_vertices;
	TCargoString edges_cargo;

	_performTripletExtension(g, set1, set2, newCargoMap, edges_vertices, edges_cargo);
	_assignNewCargos(g, newCargoMap);
	_addNewEdgesFoundByTriplet(g, edges_vertices, edges_cargo);
}





//////////////////////////////////////////////////////////////////////////////
// T-Coffee - Sum of Pairs Scoring
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScore> 
inline typename Value<TScore>::Type
sumOfPairsScore(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
				TScore const& score_type)
{
	SEQAN_CHECKPOINT
	SEQAN_TASSERT(convertAlignment(g, String<char>()) == true)

	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Infix<TString>::Type TInfix;

	TScoreValue gap = scoreGapExtend(score_type);
	TScoreValue gapOpen = scoreGapOpen(score_type);
	TSize nseq = length(stringSet(g));
	TScoreValue total = 0;
	TScoreValue mismatch = 0;
		
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	TVertexIterator it(g);
	for(;!atEnd(it);++it) {
		typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
		TOutEdgeIterator itEdge(g, *it);
		TSize count = 0;
		for(;!atEnd(itEdge);++itEdge) {
			typedef typename Iterator<TInfix, Rooted>::Type TInfixIter;
			TInfix inf1 = label(g,sourceVertex(itEdge));
			TInfix inf2 = label(g,targetVertex(itEdge));
			TInfixIter sIt1 = begin(inf1);
			TInfixIter sIt2 = begin(inf2);
			while((!atEnd(sIt1)) || (!atEnd(sIt2))) {
				mismatch += score(const_cast<TScore&>(score_type), *sIt1, *sIt2);
				goNext(sIt1); goNext(sIt2);
			}
			++count;
		}
		// How many sequences are left? --> Aligned with gaps
		total += (( (TScoreValue) (nseq - count - 1) ) * (gapOpen + ( (TScoreValue) fragmentLength(g, *it) - 1) * gap));
	}
	return total + (TScoreValue) ((double) 0.5 * (double) mismatch);
}

//////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TSegmentString, typename TScore> 
inline typename Value<TScore>::Type
sumOfPairsScore(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
				TSegmentString const& alignSeq,
				TScore const& score_type)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Infix<TString>::Type TInfix;

	// Initialization
	TScoreValue gap = scoreGapExtend(score_type);
	TScoreValue gapOpen = scoreGapOpen(score_type);
	TScoreValue total = 0;

	// Sum of pair scores
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	TSize alignSeqLen = length(alignSeq);
	for(TSize i = 0; i<alignSeqLen;++i) {
		TSize count = 0;
		TSize vertexSetLen = length(alignSeq[i]);
		TSize fragLen = 0;
		for(TSize j=0; j<vertexSetLen;++j) {
			TVertexDescriptor v1 = getValue(alignSeq[i], j);
			if ((fragLen == 0) && (v1 != nilVertex)) fragLen = fragmentLength(g, v1);
			for(TSize k=j+1; k<vertexSetLen;++k) {
				TVertexDescriptor v2 = getValue(alignSeq[i], k);
				if ((v1 == nilVertex) ||
					(v2 == nilVertex))
				{
					// Count number of pairs where one vertex is a nil vertex
					// If both are nil this pair would be removed in a pairwise alignment, so don't count
					if (!((v1 == nilVertex) &&
						(v2 == nilVertex))) ++count;
					continue;
				}
				typedef typename Iterator<TInfix,Standard>::Type TInfixIter;
				TInfix inf1 = label(g,v1);
				TInfix inf2 = label(g,v2);
				TInfixIter sIt1 = begin(inf1,Standard());
				TInfixIter sIt2 = begin(inf2,Standard());
				TInfixIter sItEnd1 = end(inf1,Standard());
				TInfixIter sItEnd2 = end(inf2,Standard());
				while((sIt1!=sItEnd1) || (sIt2!=sItEnd2)) {
					total += score(const_cast<TScore&>(score_type), *sIt1, *sIt2);
					goNext(sIt1); goNext(sIt2);
				}
			}
		}
		SEQAN_TASSERT(fragLen > 0)
		// How many sequences are left? --> Aligned with gaps
		total += (( (TScoreValue) (count) ) * (gapOpen + ( (TScoreValue) fragLen - 1) * gap));
	}
	return total;
}




}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
