#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_BASE_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// T-Coffee - Tags
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.KimuraDistance
..summary:Tag to use the kimura distance correction.
..value.KimuraDistance:Use the Kimura Distance.
*/
struct KimuraDistance_;
typedef Tag<KimuraDistance_> const KimuraDistance;

/**
.Tag.FractionalDistance
..summary:Tag to use the fractional distance correction.
..value.FractionalDistance:Use the Fractional Distance.
*/
struct FractionalDistance_;
typedef Tag<FractionalDistance_> const FractionalDistance;

/**
.Tag.GlobalPairwise_Library
..summary:Tag to specify the type of the library.
..value.GlobalPairwise_Library:Use of a pairwise global library.
*/
struct GlobalPairwise_Library_;
typedef Tag<GlobalPairwise_Library_> const GlobalPairwise_Library;

/**
.Tag.LocalPairwise_Library
..summary:Tag to specify the type of the library.
..value.LocalPairwise_Library:Use of a pairwise local library.
*/
struct LocalPairwise_Library_;
typedef Tag<LocalPairwise_Library_> const LocalPairwise_Library;

/**
.Tag.LocalTriple_Library
..summary:Tag to specify the type of the library.
..value.LocalTriple_Library:Use of pairwise local alignments to find matches present in 3 sequences.
*/
struct LocalTriple_Library_;
typedef Tag<LocalTriple_Library_> const LocalTriple_Library;


/**
.Tag.Lcs_Library
..summary:Tag to specify the type of the library.
..value.Lcs_Library:Use of a lcs library.
*/
struct Lcs_Library_;
typedef Tag<Lcs_Library_> const Lcs_Library;

/**
.Tag.MUM_Library
..summary:Tag to specify the type of the library.
..value.MUM_Library:Use of a maximal unique match library.
*/
struct MUM_Library_;
typedef Tag<MUM_Library_> const MUM_Library;


//////////////////////////////////////////////////////////////////////////////
// T-Coffee - Signal addition, library extension, ...
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TLibraries> 
inline void
combineGraphs(Graph<Alignment<TStringSet, TCargo, TSpec> >& outGraph,
			  bool onlyTopology,
			  TLibraries& libs)
{
	SEQAN_CHECKPOINT
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

template<typename TStringSet, typename TCargo, typename TSpec, typename TSet>
inline void 
restrictedTripletLibraryExtension(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
								  TSet const& set1,
								  TSet const& set2)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	// Two tasks:
	// 1) Add edges for the case that a and c is aligned, b and c is aligned, but a and b are not, give these edges the appropriate weight
	// 2) Augment all existing edges
	String<TCargo> newCargoMap;
	resize(newCargoMap, getIdUpperBound(_getEdgeIdManager(g)), Exact());
	TEdgeIterator it(g);
	for(;!atEnd(it);++it) assignProperty(newCargoMap, *it, cargo(*it));
	
	typedef String<TVertexDescriptor, Block<> > TVertexString;
	typedef String<TCargo, Block<> > TCargoString;
	TVertexString edges_vertices;
	TCargoString edges_cargo;
	TCargo avg = 0;
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
					avg += val;

					// Remember the edge with cargo
					push_back(edges_vertices, tV1);
					push_back(edges_vertices, tV2);
					push_back(edges_cargo, val);
				} else {
					if (getCargo(*outIt2) > getCargo(*outIt1)) property(newCargoMap, e) += getCargo(*outIt1);
					else property(newCargoMap, e) += getCargo(*outIt2);	
				}
			}
		}
	}
	// Assign the new weights and clean-up the cargo map
	goBegin(it);
	for(;!atEnd(it);++it) cargo(*it) = getProperty(newCargoMap, *it);
	clear(newCargoMap);

	// Finally add the new edges created by the triplet approach
	typedef typename Iterator<TVertexString, Rooted>::Type TVertexStringIter;
	typedef typename Iterator<TCargoString, Rooted>::Type TCargoStringIter;
	TVertexStringIter endIt = end(edges_vertices);
	TVertexStringIter itV = begin(edges_vertices);
	TCargoStringIter itC = begin(edges_cargo);
	while(itV != endIt) {
		TVertexStringIter itVNext = itV; ++itVNext;
		// The same edge could have been created multiple times, so check if it exists
		TEdgeDescriptor e = findEdge(g, *itV, *itVNext);
		if (e == 0) {
			addEdge(g, *itV, *itVNext, *itC);
		} else {
			cargo(e) += *itC;
		}
		++itV; ++itV;
		++itC;
	}
	SEQAN_TASSERT(itC == end(edges_cargo))
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TSet>
inline void 
tripletLibraryExtensionOnSet(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
							 TSet const& set1)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	// Two tasks:
	// 1) Add edges for the case that a and c is aligned, b and c is aligned, but a and b are not, give these edges the appropriate weight
	// 2) Augment all existing edges
	String<TCargo> newCargoMap;
	resize(newCargoMap, getIdUpperBound(_getEdgeIdManager(g)), Exact());
	TEdgeIterator it(g);
	for(;!atEnd(it);++it) assignProperty(newCargoMap, *it, cargo(*it));
	
	typedef String<TVertexDescriptor, Block<> > TVertexString;
	typedef String<TCargo, Block<> > TCargoString;
	TVertexString edges_vertices;
	TCargoString edges_cargo;
	TCargo avg = 0;
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
					avg += val;

					// Remember the edge with cargo
					push_back(edges_vertices, tV1);
					push_back(edges_vertices, tV2);
					push_back(edges_cargo, val);
				} else {
					if (getCargo(*outIt2) > getCargo(*outIt1)) property(newCargoMap, e) += getCargo(*outIt1);
					else property(newCargoMap, e) += getCargo(*outIt2);	
				}
			}
		}
	}
	// Assign the new weights and clean-up the cargo map
	goBegin(it);
	for(;!atEnd(it);++it) cargo(*it) = getProperty(newCargoMap, *it);
	clear(newCargoMap);

	// Finally add the new edges created by the triplet approach
	typedef typename Iterator<TVertexString, Rooted>::Type TVertexStringIter;
	typedef typename Iterator<TCargoString, Rooted>::Type TCargoStringIter;
	TVertexStringIter endIt = end(edges_vertices);
	TVertexStringIter itV = begin(edges_vertices);
	TCargoStringIter itC = begin(edges_cargo);
	while(itV != endIt) {
		TVertexStringIter itVNext = itV; ++itVNext;
		// The same edge could have been created multiple times, so check if it exists
		TEdgeDescriptor e = findEdge(g, *itV, *itVNext);
		if (e == 0) {
			addEdge(g, *itV, *itVNext, *itC);
		} else {
			cargo(e) += *itC;
		}
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
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	// Two tasks:
	// 1) Add edges for the case that a and c is aligned, b and c is aligned, but a and b are not, give these edges the appropriate weight
	// 2) Augment all existing edges
	String<TCargo> newCargoMap;
	resize(newCargoMap, getIdUpperBound(_getEdgeIdManager(g)), Exact());
	TEdgeIterator it(g);
	for(;!atEnd(it);++it) assignProperty(newCargoMap, *it, cargo(*it));

	typedef String<TVertexDescriptor, Block<> > TVertexString;
	typedef String<TCargo, Block<> > TCargoString;
	TVertexString edges_vertices;
	TCargoString edges_cargo;
	for(TVertexIterator itVertex(g);!atEnd(itVertex);++itVertex) {
		TOutEdgeIterator outIt1(g, *itVertex);
		while (!atEnd(outIt1)) {
			TOutEdgeIterator outIt2 = outIt1;
			goNext(outIt2);
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
	goBegin(it);
	for(;!atEnd(it);++it) cargo(*it) = getProperty(newCargoMap, *it);
	clear(newCargoMap);

	// Finally add the new edges created by the triplet approach
	typedef typename Iterator<TVertexString, Rooted>::Type TVertexStringIter;
	typedef typename Iterator<TCargoString, Rooted>::Type TCargoStringIter;
	TVertexStringIter endIt = end(edges_vertices);
	TVertexStringIter itV = begin(edges_vertices);
	TCargoStringIter itC = begin(edges_cargo);
	while(itV != endIt) {
		TVertexStringIter itVNext = itV; ++itVNext;
		// The same edge could have been created multiple times, so check if it exists
		TEdgeDescriptor e = findEdge(g, *itV, *itVNext);
		if (e == 0) {
			addEdge(g, *itV, *itVNext, *itC);
		} else {
			cargo(e) += *itC;
		}
		++itV; ++itV;
		++itC;
	}
	SEQAN_TASSERT(itC == end(edges_cargo))
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void 
fullTripletLibraryExtension(Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	// 1) Augment all existing edges according to triplet paradigm
	String<TCargo> newCargoMap;
	resize(newCargoMap, getIdUpperBound(_getEdgeIdManager(g)), Exact());
	TEdgeIterator it(g);
	TSize nEdges = 0;
	for(;!atEnd(it);++it) {
		++nEdges;
		assignProperty(newCargoMap, *it, cargo(*it));
	}
	for(TVertexIterator itVertex(g);!atEnd(itVertex);++itVertex) {
		TOutEdgeIterator outIt1(g, *itVertex);
		while (!atEnd(outIt1)) {
			TOutEdgeIterator outIt2 = outIt1;
			goNext(outIt2);
			while (!atEnd(outIt2)) {
				TVertexDescriptor tV1 = targetVertex(outIt1);
				TVertexDescriptor tV2 = targetVertex(outIt2);
				if (sequenceId(g, tV1) != sequenceId(g,tV2)) {
					TEdgeDescriptor e = findEdge(g, tV1, tV2);
					if (e != 0) {
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
	goBegin(it);
	for(;!atEnd(it);++it) cargo(*it) = getProperty(newCargoMap, *it);
	clear(newCargoMap);


	// 2) Add missing edges according to triplet approach
	typedef String<TVertexDescriptor, Block<> > TVertexString;
	typedef String<TCargo, Block<> > TCargoString;
	TVertexString edges_vertices;
	TCargoString edges_cargo;
	TCargo avg = 0;
	for(TVertexIterator itVertex(g);!atEnd(itVertex);++itVertex) {
		TOutEdgeIterator outIt1(g, *itVertex);
		while (!atEnd(outIt1)) {
			TOutEdgeIterator outIt2 = outIt1;
			goNext(outIt2);
			while (!atEnd(outIt2)) {
				TVertexDescriptor tV1 = targetVertex(outIt1);
				TVertexDescriptor tV2 = targetVertex(outIt2);
				if (sequenceId(g, tV1) != sequenceId(g,tV2)) {
					TEdgeDescriptor e = findEdge(g, tV1, tV2);
					if (e == 0) {
						// New edge
						TCargo val = cargo(*outIt1);
						if (val > cargo(*outIt2)) val = cargo(*outIt2);
						avg += val;
						// Remember the edge with cargo
						push_back(edges_vertices, tV1);
						push_back(edges_vertices, tV2);
						push_back(edges_cargo, val);
					}
				}
				goNext(outIt2);
			}
			goNext(outIt1);
		}
	}
	avg /= length(edges_cargo);
	// Add all triplet edges if the graph is small
	if (nEdges < 1000000) avg = 0;
		
	// Add the newly created edges
	typedef typename Iterator<TVertexString, Rooted>::Type TVertexStringIter;
	typedef typename Iterator<TCargoString, Rooted>::Type TCargoStringIter;
	TVertexStringIter endIt = end(edges_vertices);
	TVertexStringIter itV = begin(edges_vertices);
	TCargoStringIter itC = begin(edges_cargo);
	while(itV != endIt) {
		TVertexStringIter itVNext = itV; ++itVNext;
		// The same edge could have been created multiple times, so check if it exists
		TEdgeDescriptor e = findEdge(g, *itV, *itVNext);
		if (e == 0) {
			if (*itC > avg) addEdge(g, *itV, *itVNext, *itC);
		} else {
			cargo(e) += *itC;
		}
		++itV; ++itV;
		++itC;
	}
}

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
// T-Coffee - Main Routines
//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TAlignmentGraph>
inline void
tCoffeeProteinAlignment(StringSet<TString, Dependent<> > const& strSet,
						TAlignmentGraph& gOut)
{
	SEQAN_CHECKPOINT
	
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	typedef typename Size<TGraph>::Type TSize;

	// Score objects
	Blosum62 score_type_global(-1,-11);
	Blosum62 score_type_local(-2,-8);
	//Score<int, Pam<> > score_type_global(250, -1, -17);
	//Score<int, Pam<> > score_type_local(250, -2, -13);

	// Generate a primary library, i.e., all global pairwise alignments
	TGraph lib1(strSet);
	String<double> distanceMatrix;
	generatePrimaryLibrary(lib1, distanceMatrix, score_type_global, GlobalPairwise_Library() );
	//std::cout << "Library size: " << numVertices(lib1) << " Vertices, " << numEdges(lib1) << " Edges" << std::endl;

	// Generate a primary library, i.e., all local pairwise alignments
	TGraph lib2(strSet);
	generatePrimaryLibrary(lib2, score_type_local, LocalPairwise_Library() );
	//std::cout << "Library size: " << numVertices(lib2) << " Vertices, " << numEdges(lib2) << " Edges" << std::endl;

	// Weighting of libraries (Signal addition)
	TGraph g(strSet);
	String<TGraph*> libs;
	appendValue(libs, &lib1);
	appendValue(libs, &lib2);
	combineGraphs(g, libs);
	//// Only topology
	//combineGraphs(g, true, libs);
	//std::cout << "Library size: " << numVertices(g) << " Vertices, " << numEdges(g) << " Edges" << std::endl;

	// Clear the old libraries
	clear(lib1);
	clear(lib2);


	TSize nSeq = length(strSet);
	TSize firstThreshold = 10;
	TSize secondThreshold = 30;
	if (nSeq < firstThreshold) {
		// Full triplet...
		tripletLibraryExtension(g);
		
		// ...  and alignment according to library information, that is no guide tree
		highestScoreFirstAlignment(g, gOut);

		//// Alternative 1: Pairwise library distances and guide tree
		//pairwiseDistances(g, distanceMatrix);
		//Graph<Tree<double> > guideTree;
		//upgmaTree(distanceMatrix, guideTree);
		//progressiveAlignment(g, guideTree, gOut);
		//clear(guideTree);
	
		//// Alternative 2: Calculate a distance matrix using a compressed alphabet or not
		//getKmerSimilarityMatrix(stringSet(g), distanceMatrix, 6, AAGroupsDayhoff() );
		//similarityToDistanceMatrix(distanceMatrix, KimuraDistance() );
		//Graph<Tree<double> > guideTree;
		//upgmaTree(distanceMatrix, guideTree);
		//progressiveAlignment(g, guideTree, gOut);
		//clear(guideTree);
	} else if (nSeq < secondThreshold) {
		// Full triplet...
		tripletLibraryExtension(g);

		// ... and normal progressive alignment with guide tree
		Graph<Tree<double> > guideTree;
		slowNjTree(distanceMatrix, guideTree);
		progressiveAlignment(g, guideTree, gOut);
		clear(guideTree);
	} else {
		// Triplet only on groups of sequences
		Graph<Tree<double> > guideTree;
		slowNjTree(distanceMatrix, guideTree);	// More balanced than UPGMA
		progressiveAlignmentWithTriplet(g, guideTree, gOut, secondThreshold);
		clear(guideTree);
	}
	clear(distanceMatrix);
	clear(g);
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
