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
// Merging of alignment graphs
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TWeights, typename TLibraries> 
inline void
combineGraphs(Graph<Alignment<TStringSet, TCargo, TSpec> >& outGraph,
			  TLibraries& libs,
			  TWeights& weights,
			  FractionalScore)
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
	typedef String<TFragment> TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;
	String<TCargo> score_values;

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
			appendValue(matches, TFragment( (unsigned int) sequenceId(lib, sV), (unsigned int) fragmentBegin(lib,sV), (unsigned int) sequenceId(lib, tV),  (unsigned int)  fragmentBegin(lib,tV),  (unsigned int)  fragmentLength(lib,tV)));
			appendValue(score_values, currentCargo);
			++count;
		}
		assignValue(max_scores, i, maxCargoLib);
	}

	// Match refinement
	TStringSet& str = stringSet(outGraph);	
	matchRefinement(matches,str,outGraph);  // Don't score matches!

	// Adapt edge weights (fractional weights are used)
	count = 0;
	TSize currentLib = 0;
	double upperBound = 1000 * value(weights, currentLib);
	double scaling = upperBound / (double) getValue(max_scores, currentLib);
	TSize nextLibCounter = length(matches);
	if (currentLib < numLibs - 1) nextLibCounter = (TSize) getValue(index_start, currentLib+1);
	TFragmentStringIter endIt = end(matches);
	TSize posit = 0;
	for(TFragmentStringIter it = begin(matches); it != endIt; ++it, ++posit) {
		if (count >= nextLibCounter) {
			++currentLib;
			upperBound = 1000 * value(weights, currentLib);
			scaling = upperBound / (double) getValue(max_scores, currentLib);
			if (currentLib < numLibs - 1) nextLibCounter = (TSize) getValue(index_start, currentLib+1);
			else nextLibCounter = length(matches);
		}
		TId id1 = sequenceId(*it,0);
		TId id2 = sequenceId(*it,1);
		TSize pos1 = fragmentBegin(*it, id1);
		TSize pos2 = fragmentBegin(*it, id2);
		TSize end1 = pos1 + fragmentLength(*it, id1);
		TSize oldFragLen = (end1 - pos1);
		while(pos1 < end1) {
			TVertexDescriptor p1 = findVertex(outGraph, id1, pos1);
			TVertexDescriptor p2 = findVertex(outGraph, id2, pos2);
			TEdgeDescriptor e = findEdge(outGraph, p1, p2);
			TSize fragLen = fragmentLength(outGraph, p1); 
			double newVal = (double) (fragLen) / (double) (oldFragLen);
			newVal *= scaling;
			newVal *= (double) getValue(score_values, posit);
			SEQAN_TASSERT(e != 0)  // Refined graph should have all edges (with weight 1)
			cargo(e) += (TCargo) newVal;
			pos1 += fragLen;
			pos2 += fragLen;
		}
		++count;
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TLibraries> 
inline void
combineGraphs(Graph<Alignment<TStringSet, TCargo, TSpec> >& outGraph,
			  TLibraries& libs,
			  FrequencyCounting)
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

	// All the matches
	typedef Fragment<> TFragment;
	typedef String<TFragment> TFragmentString;
	TFragmentString matches;

	// Index start position for every library
	String<TCargo> index_start;
	resize(index_start, numLibs);

	// Get all matches
	TSize count = 0;
	for(TSize i = 0; i<numLibs; ++i) {
		assignValue(index_start, i, count);
		TGraph const& lib = *(getValue(libs, i));
		TEdgeIterator it(lib);
		for(;!atEnd(it);++it) {
			TVertexDescriptor sV = sourceVertex(it);
			TVertexDescriptor tV = targetVertex(it);
			appendValue(matches, TFragment( (unsigned int) sequenceId(lib, sV), (unsigned int) fragmentBegin(lib,sV), (unsigned int) sequenceId(lib, tV),  (unsigned int)  fragmentBegin(lib,tV),  (unsigned int)  fragmentLength(lib,tV)));
			++count;
		}
	}

	// Match refinement
	TStringSet& str = stringSet(outGraph);	
	matchRefinement(matches,str,outGraph);  // Don't score matches, just count the number of occurrences
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TLibraries> 
inline void
combineGraphs(Graph<Alignment<TStringSet, TCargo, TSpec> >& outGraph,
			  TLibraries& libs,
			  FractionalScore)
{
	SEQAN_CHECKPOINT
	String<unsigned int> weights;
	fill(weights, length(libs), 1);
	combineGraphs(outGraph, libs, weights, FractionalScore() );
}

//////////////////////////////////////////////////////////////////////////////


/**
.Function.combineGraphs:
..summary:Combines multiple alignment graphs into one single graph.
..cat:Graph
..signature:
combineGraphs(graph, libs [, weights], tag)
..param.graph:Out-parameter:The final alignment graph.
...type:Spec.Alignment Graph
..param.libs:String of pointers to alignment graph data structures.
..param.weights:String of weights.
...remarks:Needs to have the same length as the libs string.
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

/**
.Function.cliqueReduction:
..summary:Performs an edge reduction on an alignment graph producing cliques.
..cat:Graph
..signature:
cliqueReduction(graph)
..param.graph:An alignment graph.
...type:Spec.Alignment Graph
...remarks:
Finds all 3 member cliques and single-edge cliques in an alignment graph.
The found cliques never have two vertices on the same sequence.
..returns:void
*/
template<typename TStringSet, typename TCargo, typename TSpec>
inline void 
cliqueReduction(Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef String<TVertexDescriptor> TVertexString;
	typedef typename Iterator<TVertexString>::Type TVertexStringIter;
	TVertexString keepEdges;

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
						appendValue(keepEdges, *itVertex);
						appendValue(keepEdges, tV1);
						appendValue(keepEdges, *itVertex);
						appendValue(keepEdges, tV2);
						appendValue(keepEdges, tV1);
						appendValue(keepEdges, tV2);
					} 
				}
				goNext(outIt2);
			}
			goNext(outIt1);
		}
	}

	// Remove all edges and re-insert the cliques
	clearEdges(g);
	TVertexStringIter itV = begin(keepEdges);
	TVertexStringIter endIt = end(keepEdges);
	while(itV != endIt) {
		TVertexStringIter itVNext = itV; ++itVNext;
		// The same edge could have been inserted multiple times
		TEdgeDescriptor e = findEdge(g, *itV, *itVNext);
		if (e == 0) addEdge(g, *itV, *itVNext, 1);
		else cargo(e) += 1;
		++itV; ++itV;
	}
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
..param.score_type:A STL set of sequences.
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
						totalScore += score(const_cast<TScore&>(score_type), value(mat, i*len+k), value(mat, j*len + k));
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
						totalScore += score(const_cast<TScore&>(score_type), value(mat, i*len+k), value(mat, j*len + k));
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

template<typename TValue, typename TScore, typename TSize, typename TScoreValue> 
inline void
sumOfPairsScoreInd(String<TValue> const& mat,
				   TScore const& score_type,
				   TSize const nseq,
				   TSize const len,
				   String<TScoreValue>& columnScores)
{
	char gapChar = gapValue<char>();
	TScoreValue gap = scoreGapExtend(score_type);
	for(TSize k=0;k<len; ++k) {
		TScoreValue totalScore = 0;
		for(TSize i = 0; i<nseq-1; ++i) {
			for(TSize j=i+1; j<nseq; ++j) {
				if (value(mat, i*len+k) != gapChar) {
					if (value(mat, j*len + k) != gapChar) {
						totalScore += score(const_cast<TScore&>(score_type), value(mat, i*len+k), value(mat, j*len + k));
					} else totalScore += gap;
				} else if (value(mat, j*len + k) != gapChar) {
						totalScore += gap;
				}
			}
		}
		appendValue(columnScores, totalScore);
	}
}


//////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TScore, typename TDistanceMatrix>
inline void
reduceGraph(Graph<Alignment<TStringSet, TCargo, TSpec> >& anchorGraph,
			TScore const& score_type,
			TDistanceMatrix& distMatrix)
{
	typedef String<char> TAlignMatrix;
	typedef String<bool> TActive;
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TScore>::Type TScoreValue;

	// Anchor-graph parameters
	TSize minCoverage = 3;
	TSize minBlockLen = 3;
	
	// Create an alignment matrix
	TAlignMatrix mat;
	char gapChar = gapValue<char>();
	convertAlignment(anchorGraph, mat);
	TStringSet& seqSet = stringSet(anchorGraph);
	TSize nseq = length(seqSet);
	TSize len = length(mat) / nseq;
	
	// Mark the gaps and the characters
	typedef typename Iterator<TAlignMatrix>::Type TAlignIter;
	typedef typename Iterator<TActive>::Type TActiveIter;
	TActive active;
	fill(active, length(mat), true);
	String<unsigned int> coverage;
	fill(coverage, len, 0);
	TAlignIter itA = begin(mat);
	TAlignIter itAEnd = end(mat);
	TActiveIter itB = begin(active);
	TSize col = 0;
	for(; itA != itAEnd; ++itA, ++itB, col = (++col) % len) {
		if (*itA == gapChar) value(itB) = false;
		else value(coverage, col) += 1;
	}

	// Clear blocks with low coverage
	itA = begin(mat);
	itB = begin(active);
	col = 0;
	for(;itA != itAEnd; ++itA, ++itB, col = (++col) % len) {
		if (value(coverage, col) < minCoverage) value(itB) = false;
	}

	// Collect all blocks of minimum length
	TSize blockStart = len;
	TSize blockEnd = 0;
	typedef std::set<TSize> TSequenceSet;
	typedef Triple<TSequenceSet, TSize, TSize> TBlock;
	typedef String<TBlock> TConservedBlock;
	TConservedBlock conservedBlock;
	TSequenceSet thisSeqSet;
	TSequenceSet nextSeqSet;
	for(TSize column = 0; column < len; ++column) {
		for(TSize seq = 0; seq < nseq; ++seq) {
			if (!value(active, seq * len + column)) continue;
			nextSeqSet.insert(seq);
		}
		if (thisSeqSet != nextSeqSet) {
			if (blockStart == len) {
				blockStart = column;
			} else  {
				if ((column - blockStart) > minBlockLen) appendValue(conservedBlock, TBlock(thisSeqSet, blockStart, column));
				if (nextSeqSet.empty()) blockStart = len;
				else blockStart = column;
			}
			thisSeqSet = nextSeqSet;
		}
		nextSeqSet.clear();
	}
	clear(active);
	fill(active, length(mat), false);
	typedef typename Iterator<TConservedBlock>::Type TConsBlockIter;
	
	// Check quality of blocks
	TConsBlockIter itCB = begin(conservedBlock);
	TConsBlockIter itCBEnd = end(conservedBlock);

	//typedef std::set<std::pair<TScoreValue, TConsBlockIter> > TBlockQuality;
	//TBlockQuality bQual;
	//for(;itCB != itCBEnd; ++itCB) {
	//	TScoreValue totalScore = 0;
	//	for(TSize column = (value(itCB)).i2; column < (value(itCB)).i3; ++column) {
	//		typename TSequenceSet::const_iterator itCons = (value(itCB)).i1.begin();
	//		typename TSequenceSet::const_iterator itConsEnd = (value(itCB)).i1.end();
	//		for(;itCons != itConsEnd; ++itCons) {
	//			typename TSequenceSet::const_iterator pairItCons = itCons;
	//			++pairItCons;
	//			for(;pairItCons != itConsEnd; ++pairItCons) {
	//				totalScore += score(const_cast<TScore&>(score_type), value(mat, *itCons * len + column), value(mat, *pairItCons * len + column));
	//			}
	//		}
	//	}
	//	totalScore /= (TScoreValue) (value(itCB)).i1.size();
	//	bQual.insert(std::make_pair(totalScore, itCB));
	//	//if (totalScore < 0) (value(itCB)).i1.clear();
	//}
	//// Outlier detection
	//TSize n = length(conservedBlock);
	//typename TBlockQuality::const_iterator scrIt = bQual.begin();
	//typename TBlockQuality::const_iterator scrItEnd = bQual.end();
	//TScoreValue upperQuartile = 0;
	//TScoreValue lowerQuartile = 0;
	//TSize count = 0;
	//for(;scrIt != scrItEnd; ++scrIt, ++count) {
	//	if (n / 4 == count) lowerQuartile = scrIt->first;
	//	else if ( (3 * n) / 4 == count) upperQuartile = scrIt->first;
	//}
	//TScoreValue interQR = upperQuartile - lowerQuartile;
	//TScoreValue lowerInnerFence = (TScoreValue) ((double) lowerQuartile - 0.5 * (double) interQR);
	//scrIt = bQual.begin();
	//scrItEnd = bQual.end();
	//for(;((scrIt != scrItEnd) && (scrIt->first < lowerInnerFence)); ++scrIt) {
	//	(value(scrIt->second)).i1.clear();			
	//}

	//// Check blocks for outliers
	//typedef std::pair<TScoreValue, TSize> TQualitySeqPair;
	//typedef std::set<TQualitySeqPair, std::greater<TQualitySeqPair> > TSeqQuality;
	//itCB = begin(conservedBlock);
	//itCBEnd = end(conservedBlock);
	//for(;itCB != itCBEnd; ++itCB) {
	//	TSeqQuality seqQuality;
	//	typename TSequenceSet::const_iterator leaveOut = (value(itCB)).i1.begin();
	//	typename TSequenceSet::const_iterator leaveOutEnd = (value(itCB)).i1.end();
	//	for(;leaveOut != leaveOutEnd; ++leaveOut) {
	//		TScoreValue totalScore = 0;
	//		for(TSize column = (value(itCB)).i2; column < (value(itCB)).i3; ++column) {
	//			typename TSequenceSet::const_iterator itCons = (value(itCB)).i1.begin();
	//			typename TSequenceSet::const_iterator itConsEnd = (value(itCB)).i1.end();
	//			for(;itCons != itConsEnd; ++itCons) {
	//				if (*itCons == *leaveOut) continue;
	//				typename TSequenceSet::const_iterator pairItCons = itCons;
	//				++pairItCons;
	//				for(;pairItCons != itConsEnd; ++pairItCons) {
	//					if (*pairItCons == *leaveOut) continue;
	//					totalScore += score(const_cast<TScore&>(score_type), value(mat, *itCons * len + column), value(mat, *pairItCons * len + column));
	//				}
	//			}
	//		}
	//		seqQuality.insert(std::make_pair(totalScore, *leaveOut));
	//	}
	//	// Outlier detection
	//	TSize n = (value(itCB)).i1.size();
	//	typename TSeqQuality::const_iterator screenIt = seqQuality.begin();
	//	typename TSeqQuality::const_iterator screenItEnd = seqQuality.end();
	//	TScoreValue upperQuartile = 0;
	//	TScoreValue lowerQuartile = 0;
	//	TSize count = 0;
	//	for(;screenIt != screenItEnd; ++screenIt, ++count) {
	//		if (n / 4 == count) upperQuartile = screenIt->first;
	//		else if ( (3 * n) / 4 == count) lowerQuartile = screenIt->first;
	//	}
	//	TScoreValue interQR = upperQuartile - lowerQuartile;
	//	TScoreValue upperInnerFence = upperQuartile + 3 * interQR;
	//	screenIt = seqQuality.begin();
	//	screenItEnd = seqQuality.end();
	//	for(;((screenIt != screenItEnd) && (screenIt->first > upperInnerFence)); ++screenIt) {
	//		(value(itCB)).i1.erase(screenIt->second);			
	//	}
	//}

	
	
	// Debug code: All blocks
	itCB = begin(conservedBlock);
	itCBEnd = end(conservedBlock);
	for(;itCB != itCBEnd; ++itCB) {
		typename TSequenceSet::const_iterator itCons = (value(itCB)).i1.begin();
		typename TSequenceSet::const_iterator itConsEnd = (value(itCB)).i1.end();
		for(;itCons != itConsEnd; ++itCons) {
			for(TSize column = (value(itCB)).i2; column < (value(itCB)).i3; ++column) {
				std::cout << value(mat, *itCons * len + column);
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}


	itCB = begin(conservedBlock);
	itCBEnd = end(conservedBlock);
	for(;itCB != itCBEnd; ++itCB) {
		typename TSequenceSet::const_iterator itCons = (value(itCB)).i1.begin();
		typename TSequenceSet::const_iterator itConsEnd = (value(itCB)).i1.end();
		for(;itCons != itConsEnd; ++itCons) {
			for(TSize column = (value(itCB)).i2; column < (value(itCB)).i3; ++column) {
				value(active, *itCons * len + column) = true;
			}
		}
	}


	// Compute a distance matrix from the anchor graph
	clear(distMatrix);
	resize(distMatrix, nseq * nseq);
	for(TSize seq1 = 0; seq1 < nseq-1; ++seq1) {
		for(TSize seq2 = seq1+1; seq2 < nseq; ++seq2) {
			typename Value<TDistanceMatrix>::Type diff = 0;
			for(TSize column = 0; column < len; ++column) {
				if ((value(active, seq1 * len + column) != value(active, seq2 * len + column)) ||
					(value(active, seq1 * len + column) == 0))	++diff;
			}
			value(distMatrix, seq1 * nseq + seq2) = diff;
		}
	}

	// Create the anchor graph
	typedef Fragment<> TFragment;
	typedef String<TFragment> TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;
	typedef std::pair<TSize, TSize> TResiduePair;
	typedef std::set<TResiduePair> TResiduePairSet;
	String<TResiduePairSet> resPair;
	resize(resPair, nseq * nseq);	
	for(TSize seq1 = 0; seq1 < nseq - 1; ++seq1) {
		for(TSize seq2 = seq1 + 1; seq2 < nseq; ++seq2) {
			TSize index = seq1 * nseq + seq2;
			TSize offset1 = 0;
			TSize offset2 = 0;
			for(TSize col = 0; col<len; ++col) {
				if (value(mat, seq1 * len + col) != gapChar) {
					if (value(mat, seq2 * len + col) != gapChar) {
						if ((value(active, seq1 * len + col) == true) && (value(active, seq2 * len + col) == true)) {
							resPair[index].insert(std::make_pair(offset1, offset2));
						}
						++offset1;
						++offset2;
					} else ++offset1;
				} else if (value(mat, seq2 * len + col) != gapChar) ++offset2;
			}
		}
	}
	for(TSize i = 0; i<length(resPair); ++i) {
		if (resPair[i].empty()) continue;
		TSize seq1 = i / nseq;
		TSize seq2 = i % nseq;
		typename TResiduePairSet::const_iterator pos = resPair[i].begin();
		typename TResiduePairSet::const_iterator posEnd = resPair[i].end();
		TSize startMatch1 = pos->first;
		TSize startMatch2 = pos->second;
		TSize len = 1;
		++pos;
		while(pos != posEnd) {
			if ((startMatch1 + len == pos->first) && (startMatch2 + len == pos->second)) ++len;
			else {
				appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
				startMatch1 = pos->first;
				startMatch2 = pos->second;
				len = 1;
			}
			++pos;
		}
		appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
	}
	clearVertices(anchorGraph);
	matchRefinement(matches,stringSet(anchorGraph),const_cast<TScore&>(score_type),anchorGraph);
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
