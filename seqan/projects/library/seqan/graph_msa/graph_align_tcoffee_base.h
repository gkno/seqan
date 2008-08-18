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
  $Id: graph_align_tcoffee_base.h 1905 2008-04-29 11:17:24Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_BASE_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Tags
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Alignment Graph Combination:
..summary:A tag to specify how to combine alignment graphs.
*/

/**
.Tag.Alignment Graph Combination.value.FractionalScore:
	Rescores matches with the appropriate fractional score.
*/
struct FractionalScore_;
typedef Tag<FractionalScore_> const FractionalScore;


/**
.Tag.Alignment Graph Combination.value.FrequencyCount:
	Rescores matches with the frequency count for this edge.
*/
struct FrequencyCounting_;
typedef Tag<FrequencyCounting_> const FrequencyCounting;



//////////////////////////////////////////////////////////////////////////////
// Generating an alignment graph from segment matches
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TFragment, typename TSpec1, typename TScoreValue, typename TSpec2, typename TStringSet, typename TCargo, typename TSpec> 
inline void
buildAlignmentGraph(String<TFragment, TSpec1>& matches,
					String<TScoreValue, TSpec2>& scores,
					Graph<Alignment<TStringSet, TCargo, TSpec> >& outGraph,
					FractionalScore)
{
	SEQAN_CHECKPOINT
	typedef String<TFragment, TSpec1> TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	typedef String<TScoreValue, TSpec2> TScoreValues;
	typedef typename Iterator<TScoreValues>::Type TScoreValuesIter;
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TOutGraph;
	typedef typename Size<TFragmentString>::Type TSize;
	typedef typename Id<TOutGraph>::Type TId;
	typedef typename EdgeDescriptor<TOutGraph>::Type TEdgeDescriptor;
	typedef typename VertexDescriptor<TOutGraph>::Type TVertexDescriptor;
	
	// Initialization
	clearVertices(outGraph);
	TStringSet& strSet = stringSet(outGraph);
	
	// Segment-match refinement
	matchRefinement(matches,strSet,outGraph);

	// Find smallest score value
	TScoreValuesIter scoreIt = begin(scores);
	TScoreValuesIter scoreItEnd = end(scores);
	TScoreValue smallest = 0;
	for(;scoreIt != scoreItEnd; ++scoreIt) {
		if (value(scoreIt) < smallest) {
			smallest = value(scoreIt);
		}
	}

	// Clear edge-weights
	typedef typename Iterator<TOutGraph, EdgeIterator>::Type TEdgeIterator;
	TEdgeIterator itE(outGraph);
	for(;!atEnd(itE);goNext(itE)) cargo(value(itE)) = 0;

	// Adapt the scores
	TFragmentStringIter endIt = end(matches);
	scoreIt = begin(scores);
	TScoreValue offset = (-1) * smallest + 1;
	for(TFragmentStringIter it = begin(matches); it != endIt; ++it, ++scoreIt) {
		// Some alignment scores might be negative, offset the scores
		value(scoreIt) += offset;
		TId id1 = sequenceId(*it,0);
		TId id2 = sequenceId(*it,1);
		TSize pos1 = fragmentBegin(*it, id1);
		TSize pos2 = fragmentBegin(*it, id2);
		TSize fragLen = fragmentLength(*it, id1);
		TSize end1 = pos1 + fragmentLength(*it, id1);
		while(pos1 < end1) {
			TVertexDescriptor p1 = findVertex(outGraph, id1, pos1);
			TVertexDescriptor p2 = findVertex(outGraph, id2, pos2);
			TEdgeDescriptor e = findEdge(outGraph, p1, p2);
			cargo(e) += (TCargo) (((double) fragmentLength(outGraph, p1) / (double) fragLen) * (double) value(scoreIt));
			pos1 += fragmentLength(outGraph, p1);
			pos2 += fragmentLength(outGraph, p2);
		}
	}

}


//////////////////////////////////////////////////////////////////////////////

template<typename TFragment, typename TSpec1, typename TScoreValue, typename TSpec2, typename TStringSet, typename TCargo, typename TSpec> 
inline void
buildAlignmentGraph(String<TFragment, TSpec1>& matches,
					String<TScoreValue, TSpec2>&,
					Graph<Alignment<TStringSet, TCargo, TSpec> >& outGraph,
					FrequencyCounting)
{
	SEQAN_CHECKPOINT
	typedef String<TFragment, TSpec1> TFragmentString;
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TOutGraph;
	typedef typename Size<TFragmentString>::Type TSize;
	
	// Initialization
	clearVertices(outGraph);
	TStringSet& strSet = stringSet(outGraph);
	
	// Segment-match refinement
	/*
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	typedef typename Id<TOutGraph>::Type TId;
	typedef typename EdgeDescriptor<TOutGraph>::Type TEdgeDescriptor;
	typedef typename VertexDescriptor<TOutGraph>::Type TVertexDescriptor;
	typedef std::set<std::pair<TId, TId> > TSeqPairs;
	TSize nseq = length(strSet);
	String<TSeqPairs> edgesPerSeqPair;
	resize(edgesPerSeqPair, nseq * nseq);
	clearVertices(outGraph);
	for(TSize i=0;i<nseq;++i) {
		TId id = positionToId(strSet, i);
		for(TSize k=0;k<length(value(strSet,i));++k) {
			addVertex(outGraph, id, k, 1);
		}
	}
	TFragmentStringIter endIt = end(matches);
	for(TFragmentStringIter it = begin(matches); it != endIt; ++it) {
		TId id1 = sequenceId(*it,0);
		TSize pos1 = fragmentBegin(*it, id1);
		TSize len = fragmentLength(*it, id1);
		for(TSize p = 0; p < len; ++p) {
			TSize pos2 = 0;
			TId id2 = 0;
			getProjectedPosition(*it, id1, pos1 + p, id2, pos2);
			TVertexDescriptor v1 = findVertex(outGraph, id1, pos1 + p);
			TVertexDescriptor v2 = findVertex(outGraph, id2, pos2);
			TEdgeDescriptor e = findEdge(outGraph, v1, v2);
			if (e == 0) addEdge(outGraph, v1, v2, 1);
			else cargo(e) += 1;
		}
	}*/
	matchRefinement(matches,strSet,outGraph);
}



//////////////////////////////////////////////////////////////////////////////
// Merging of alignment graphs
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TLibraries, typename TTag> 
inline void
combineGraphs(Graph<Alignment<TStringSet, TCargo, TSpec> >& outGraph,
			  TLibraries& libs,
			  TTag)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	
	// Initialization
	TSize numLibs = length(libs);	
	typedef Fragment<> TFragment;
	typedef String<TFragment> TFragmentString;
	TFragmentString matches;
	String<TCargo> scores;

	// Get all matches
	for(TSize i = 0; i<numLibs; ++i) {
		TGraph const& lib = *(getValue(libs, i));
		TEdgeIterator it(lib);
		for(;!atEnd(it);++it) {
			TVertexDescriptor sV = sourceVertex(it);
			TVertexDescriptor tV = targetVertex(it);
			appendValue(matches, TFragment( (unsigned int) sequenceId(lib, sV), (unsigned int) fragmentBegin(lib,sV), (unsigned int) sequenceId(lib, tV),  (unsigned int)  fragmentBegin(lib,tV),  (unsigned int)  fragmentLength(lib,tV)));
			appendValue(scores, cargo(value(it)));
		}
	}

	// Build graph
	buildAlignmentGraph(matches, scores, outGraph, TTag());
}


//////////////////////////////////////////////////////////////////////////////


/**
.Function.combineGraphs:
..summary:Combines multiple alignment graphs into one single graph.
..cat:Graph
..signature:
combineGraphs(graph, libs, tag)
..param.graph:Out-parameter:The final alignment graph.
...type:Spec.Alignment Graph
..param.libs:String of pointers to alignment graph data structures.
..param.tag:Combination strategy.
...type:Tag.Alignment Graph Combination
..returns:void
*/
template<typename TStringSet, typename TCargo, typename TSpec, typename TLibraries> 
inline void
combineGraphs(Graph<Alignment<TStringSet, TCargo, TSpec> >& outGraph,
			  TLibraries& libs)
{
	SEQAN_CHECKPOINT
	combineGraphs(outGraph, libs, FractionalScore() );
}











//////////////////////////////////////////////////////////////////////////////
// Consistency: Triplet extension
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
						appendValue(edges_vertices, tV1);
						appendValue(edges_vertices, tV2);
						appendValue(edges_cargo, val);
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
					appendValue(edges_vertices, tV1);
					appendValue(edges_vertices, tV2);
					appendValue(edges_cargo, val);
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
					appendValue(edges_vertices, tV1);
					appendValue(edges_vertices, tV2);
					appendValue(edges_cargo, val);
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


/**
.Function.tripletLibraryExtension:
..summary:Performs the consistency extension proposed in T-Coffee.
..cat:Graph
..signature:
tripletLibraryExtension(graph, [,set1] [,set2])
..param.graph:An alignment graph.
...type:Spec.Alignment Graph
..param.set1:A STL set of sequences.
...remarks:If only one set is given the triplet extension is limited to this set of sequences. Edges to other sequences remain unchanged.
..param.set2:A STL set of sequences.
...remarks:If a second set is given, the triplet extension is performed only on edges that have one vertex in set1 and the other one in set2.
..returns:void
*/
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
	typedef String<TVertexDescriptor> TVertexString;
	typedef String<TCargo> TCargoString;
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
	typedef String<TVertexDescriptor> TVertexString;
	typedef String<TCargo> TCargoString;
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
	typedef String<TVertexDescriptor> TVertexString;
	typedef String<TCargo> TCargoString;
	TVertexString edges_vertices;
	TCargoString edges_cargo;

	_performTripletExtension(g, set1, set2, newCargoMap, edges_vertices, edges_cargo);
	_assignNewCargos(g, newCargoMap);
	_addNewEdgesFoundByTriplet(g, edges_vertices, edges_cargo);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void 
reducedTripletLibraryExtension(Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;


	// Just augment existing edges
	String<TCargo> newCargoMap;
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
					if (e != 0) {
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

	// Assign the new weights and clean-up the cargo map
	TEdgeIterator itE(g);
	for(;!atEnd(itE);goNext(itE)) cargo(value(itE)) = property(newCargoMap, value(itE));
	clear(newCargoMap);
}

//////////////////////////////////////////////////////////////////////////////
// Sum of Pairs Scoring
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// This version is sensitive to gap openings

/**
.Function.sumOfPairsScore:
..summary:Given a multiple alignment, this function calculates the sum-of-pairs score.
..cat:Graph
..signature:
sumOfPairsScore(graph, score_type)
..param.graph:An alignment graph.
...type:Spec.Alignment Graph
..param.score_type:A score object.
...type:Class.Score
..remarks:This function does NOT assume independent columns. 
That is, gap openings are properly scored. 
If you want the fast version assuming independ columns use sumOfPairsScoreInd.
..returns:void
*/
template<typename TStringSet, typename TCargo, typename TSpec, typename TScore> 
inline typename Value<TScore>::Type
sumOfPairsScore(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
				TScore const& score_type)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Value<typename Value<TStringSet>::Type>::Type TAlphabet;

	// Convert the graph
	String<char> mat;
	convertAlignment(g, mat);
	char gapChar = gapValue<char>();

	TScoreValue gap = scoreGapExtend(score_type);
	TScoreValue gapOpen = scoreGapOpen(score_type);
	TSize nseq = length(stringSet(g));
	TSize len = length(mat) / nseq;
	
	bool gapOpeni = false;
	bool gapOpenj = false;
	TScoreValue totalScore = 0;
	for(TSize i = 0; i<nseq-1; ++i) {
		for(TSize j=i+1; j<nseq; ++j) {
			for(TSize k=0;k<len; ++k) {
				if (value(mat, i*len+k) != gapChar) {
					if (value(mat, j*len + k) != gapChar) {
						gapOpeni = false;
						gapOpenj = false;
						totalScore += score(const_cast<TScore&>(score_type), TAlphabet(value(mat, i*len+k)), TAlphabet(value(mat, j*len + k)));
					} else {
						if (gapOpenj) {
							totalScore += gap;
						} else {
							gapOpenj = true;
							totalScore += gapOpen;
						}
					}
				} else if (value(mat, j*len + k) != gapChar) {
						if (gapOpeni) {
							totalScore += gap;
						} else {
							gapOpeni = true;
							totalScore += gapOpen;
						}
				}
			}
		}
	}
	return totalScore;
}

//////////////////////////////////////////////////////////////////////////////
// This version is insensitive to gap openings, assumes independent columns
template<typename TStringSet, typename TCargo, typename TSpec, typename TScore> 
inline typename Value<TScore>::Type
sumOfPairsScoreInd(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
				   TScore const& score_type)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Value<typename Value<TStringSet>::Type>::Type TAlphabet;

	// Convert the graph
	String<char> mat;
	convertAlignment(g, mat);
	char gapChar = gapValue<char>();

	TScoreValue gap = scoreGapExtend(score_type);
	TSize nseq = length(stringSet(g));
	TSize len = length(mat) / nseq;
	
	TScoreValue totalScore = 0;
	for(TSize k=0;k<len; ++k) {
		for(TSize i = 0; i<nseq-1; ++i) {
			for(TSize j=i+1; j<nseq; ++j) {
				if (value(mat, i*len+k) != gapChar) {
					if (value(mat, j*len + k) != gapChar) {
						totalScore += score(const_cast<TScore&>(score_type), TAlphabet(value(mat, i*len+k)), TAlphabet(value(mat, j*len + k)));
					} else totalScore += gap;
				} else if (value(mat, j*len + k) != gapChar) {
						totalScore += gap;
				}
			}
		}
	}
	return totalScore;
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.alignmentEvaluation:
..summary:Given a multiple alignment, this function calculates all kinds of alignment statistics.
..cat:Graph
..signature:
alignmentEvaluation(graph, score_type, gapExCount, gapCount, pairCount, numPairs, len)
..param.graph:An alignment graph.
...type:Spec.Alignment Graph
..param.score_type:A score object.
...type:Class.Score
..param.gapExCount:Number of gap extensions.
..param.gapCount:Number of gaps.
..param.pairCount:Number of aligned pairs.
..param.numPairs:Counter for each pair.
..param.len:Alignment length.
..returns:Score of the alignment.
*/
template<typename TStringSet, typename TCargo, typename TSpec, typename TScore, typename TSize> 
inline typename Value<TScore>::Type
alignmentEvaluation(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
					TScore const& score_type,
					TSize& gapExCount,
					TSize& gapCount,
					TSize& pairCount,
					String<TSize>& numPairs,
					TSize& len)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Value<typename Value<TStringSet>::Type>::Type TAlphabet;
	TSize alphSize = ValueSize<TAlphabet>::VALUE;

	// Initialization;
	gapExCount = 0;
	gapCount = 0;
	pairCount = 0;
	clear(numPairs);

	// Convert the graph
	String<char> mat;
	convertAlignment(g, mat);
	char gapChar = gapValue<char>();

	TScoreValue gap = scoreGapExtend(score_type);
	TScoreValue gapOpen = scoreGapOpen(score_type);
	TSize nseq = length(stringSet(g));
	len = length(mat) / nseq;
	
	bool gapOpeni = false;
	bool gapOpenj = false;
	TScoreValue totalScore = 0;
	fill(numPairs, alphSize * alphSize, 0);
	for(TSize i = 0; i<nseq-1; ++i) {
		for(TSize j=i+1; j<nseq; ++j) {
			for(TSize k=0;k<len; ++k) {
				if (value(mat, i*len+k) != gapChar) {
					if (value(mat, j*len + k) != gapChar) {
						gapOpeni = false;
						gapOpenj = false;
						++pairCount;
						TSize index1 = ordValue(TAlphabet(value(mat, i*len+k)));
						TSize index2 = ordValue(TAlphabet(value(mat, j*len + k)));
						value(numPairs, index1 * alphSize + index2) += 1;
						totalScore += score(const_cast<TScore&>(score_type), TAlphabet(value(mat, i*len+k)), TAlphabet(value(mat, j*len + k)));
					} else {
						if (gapOpenj) {
							++gapExCount;
							totalScore += gap;
						} else {
							gapOpenj = true;
							++gapCount;
							totalScore += gapOpen;
						}
					}
				} else if (value(mat, j*len + k) != gapChar) {
						if (gapOpeni) {
							++gapExCount;
							totalScore += gap;
						} else {
							++gapCount;
							gapOpeni = true;
							totalScore += gapOpen;
						}
				}
			}
		}
	}
	return totalScore;
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
