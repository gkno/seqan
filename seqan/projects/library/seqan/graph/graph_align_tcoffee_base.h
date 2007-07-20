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
// T-Coffee - Signal addition and library extension
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TLibraries> 
inline void
combineGraphs(Graph<Alignment<TStringSet, TCargo, TSpec> >& outGraph,
			  TLibraries& libs)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;

	// Clear out-library
	clearVertices(outGraph);
	TSize numLibs = length(libs);

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
		TGraph* lib = getValue(libs, i);
		TEdgeIterator it(*lib);
		for(;!atEnd(it);++it) {
			TCargo currentCargo = getCargo(*it);
			if (currentCargo > maxCargoLib) maxCargoLib = currentCargo;
			TVertexDescriptor sV = sourceVertex(it);
			TVertexDescriptor tV = targetVertex(it);
			push_back(matches, TFragment( (unsigned int) sequenceId(*lib, sV), (unsigned int) fragmentBegin(*lib,sV), (unsigned int) sequenceId(*lib, tV),  (unsigned int)  fragmentBegin(*lib,tV),  (unsigned int)  fragmentLength(*lib,tV)));
			push_back(score_values, currentCargo);
			++count;
		}
		assignValue(max_scores, i, maxCargoLib);
	}

	// Match refinement
	TStringSet& str = stringSet(outGraph);	
	matchRefinement(matches,str,outGraph);

	// Adapt edge weights (fractional weights are used)
	count = 0;
	TSize currentLib = 0;
	TFragmentStringIter endIt = end(matches);
	for(TFragmentStringIter it = begin(matches); it != endIt; ++it) {
		if ((currentLib < numLibs - 1) && (count >= (TSize) getValue(index_start, currentLib+1))) ++currentLib;
		double scaling = (double) 100 / (double) getValue(max_scores, currentLib);
		TId id1 = sequenceId(*it,0);
		TId id2 = sequenceId(*it,1);
		TSize pos1 = fragmentBegin(*it, id1);
		TSize pos2 = fragmentBegin(*it, id2);
		TSize end1 = pos1 + fragmentLength(*it, id1);
		while(pos1 < end1) {
			SEQAN_TASSERT(pos2 < pos2 + fragmentLength(*it, id2))
			TVertexDescriptor p1 = findVertex(outGraph, id1, pos1);
			TVertexDescriptor p2 = findVertex(outGraph, id2, pos2);
			TEdgeDescriptor e = findEdge(outGraph, p1, p2);
			TSize fragLen = fragmentLength(outGraph, p1); 
			double newVal = (double) fragLen / (double) (end1 - pos1);
			newVal *= scaling;
			newVal *= (double) getValue(score_values, position(it));
			if (e != 0) cargo(e) += (TCargo) newVal;
			else addEdge(outGraph, p1, p2, (TCargo) newVal);
			SEQAN_TASSERT(fragLen == fragmentLength(outGraph, p2))
			pos1 += fragLen;
			pos2 += fragLen;
		}
		++count;
	}
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

	if (numEdges(g) > 5000000) return; // Too many edges, no triplet extension, just return

	// Two tasks:
	// 1) Add edges for the case that a and c is aligned, b and c is aligned, but a and b are not, give these edges the appropriate weight
	// 2) Augment all existing edges
	String<TCargo> newCargoMap;
	resize(newCargoMap, getIdUpperBound(_getEdgeIdManager(g)), Exact());
	TEdgeIterator it(g);
	for(;!atEnd(it);++it) assignProperty(newCargoMap, *it, cargo(*it));
	typedef std::map<std::pair<TVertexDescriptor, TVertexDescriptor>, TCargo> TNewEdgeMap;
	TNewEdgeMap edges;
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
						if (tV1 < tV2) {
							typename TNewEdgeMap::iterator pos = edges.find(std::make_pair(tV1, tV2));
							if (pos == edges.end()) {
								edges.insert(std::make_pair(std::make_pair(tV1, tV2), val));
							} else {
								pos->second += val;
							}
						} else {
							typename TNewEdgeMap::iterator pos = edges.find(std::make_pair(tV2, tV1));
							if (pos == edges.end()) {
								edges.insert(std::make_pair(std::make_pair(tV2, tV1), val));
							} else {
								pos->second += val;
							}
						}
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
	for(typename TNewEdgeMap::const_iterator pos = edges.begin();pos!=edges.end();++pos) {
		addEdge(g, pos->first.first, pos->first.second, pos->second);
	}

	// Clean-up the edge map
	edges.clear();
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
