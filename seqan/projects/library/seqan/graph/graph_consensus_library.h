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

#ifndef SEQAN_HEADER_GRAPH_CONSENSUS_LIBRARY_H
#define SEQAN_HEADER_GRAPH_CONSENSUS_LIBRARY_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TRead, typename TSpec2, typename TBegEndPositions>
inline void 
layoutReads(String<TString, TSpec>& names,
			StringSet<TRead, TSpec2>& strSet,
			TBegEndPositions& begEndPos)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TString>::Type TSize;
	typedef typename Value<TBegEndPositions>::Type TPair;
	typedef typename Iterator<String<TString, TSpec> >::Type TNameIter;
	typedef typename Iterator<TString>::Type TIter;

	reserve(begEndPos, length(strSet));
	TNameIter itName = begin(names);
	TNameIter itNameEnd = end(names);
	TSize count = 0;
	for(;itName != itNameEnd; ++itName, ++count) {
		TSize brPoint = 0;
		TIter nIt = begin(value(itName));
		TIter nItEnd = end(value(itName));
		TSize ind = 0;
		for(;nIt!=nItEnd;++ind, ++nIt) {
			if (value(nIt) == ',') {
				brPoint = ind;
				break;
			}
		}
		TSize posI = 0;
		TSize posJ = 0;
		TString inf1 = infix(value(itName), 0, brPoint);
		TString inf2 = infix(value(itName), brPoint+1, length(value(itName)));
		std::stringstream ssStream1(toCString(inf1));
		ssStream1 >> posI; 
		std::stringstream ssStream2(toCString(inf2));
		ssStream2 >> posJ;
		appendValue(begEndPos, TPair(posI, posJ));
		if (posI > posJ) reverseComplementInPlace(value(strSet, count));
	}
}


//////////////////////////////////////////////////////////////////////////////

// Depending on the read offset
template<typename TStringSet, typename TCargo, typename TSpec, typename TBegEndPositions, typename TPairList, typename TReadLength>
inline void 
selectPairsForLibraryGeneration(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
								TBegEndPositions& begEndPos,
								TPairList& pList,
								TReadLength readLen)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TPairList>::Type TPair;
	typedef typename Iterator<TBegEndPositions>::Type TBegEndIter;
	
	TStringSet& str = stringSet(g);
	int threshold = (int) (1.5 * (double) readLen);
	TBegEndIter beIt = begin(begEndPos);
	TBegEndIter beItEnd = end(begEndPos);
	TSize index1 = 0;
	for(;beIt != beItEnd; ++beIt, ++index1) {
		int posI = (value(beIt)).i1;
		if (posI > (int) (value(beIt)).i2) posI = (value(beIt)).i2;
		TBegEndIter beIt2 = beIt;
		++beIt2;
		TSize index2 = index1 + 1;
		for(;beIt2 != beItEnd; ++beIt2, ++index2) {
			int posJ = (value(beIt2)).i1;
			if (posJ > (int) (value(beIt2)).i2) posJ = (value(beIt2)).i2;
			if (((posI > posJ) &&
				(posI - posJ < threshold))) {
					appendValue(pList, TPair(positionToId(str, index2), positionToId(str, index1)));
			} else if (((posI <= posJ) &&
						(posJ - posI < threshold))) {
					appendValue(pList, TPair(positionToId(str, index1), positionToId(str, index2)));
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TFragment, typename TSpec1, typename TStringSet, typename TSize>
inline void 
getAlignmentStatistics(String<TFragment, TSpec1>& matches,
					   TStringSet& str,
					   TSize& from,
					   TSize& matchLength,
					   TSize& overlapLength,
					   TSize& alignLength)
{
	SEQAN_CHECKPOINT
	getAlignmentStatistics(matches, str, (TSize) from, (TSize) length(matches), matchLength, overlapLength, alignLength, typename Value<TStringSet>::Type() );
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSize>
inline void 
__getAlignmentStatistics(Nothing&,
						 TSize,
						 TSize,
						 TSize,
						 double,
						 double)
{
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec, typename TSize>
inline void 
__getAlignmentStatistics(String<TValue, TSpec>& dist,
						 TSize i,
						 TSize j,
						 TSize nseq,
						 double matchLen,
						 double quality)
{
	SEQAN_CHECKPOINT
	TValue sim = matchLen * quality;
	if (sim > 1000) sim = 1000;
	assignValue(dist, i*nseq + j, 1000 - sim);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TSize>
inline void 
__getAlignmentStatistics(Graph<Undirected<TCargo, TSpec> >& dist,
						 TSize i,
						 TSize j,
						 TSize,
						 double matchLen,
						 double quality)
{
	SEQAN_CHECKPOINT
	TCargo sim = matchLen * quality;
	if (sim > 1000) sim = 1000;
	addEdge(dist, i, j, 1000 - sim);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TPair, typename TPairSpec, typename TDistance, typename TScore>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   String<TPair, TPairSpec> const& pList,
					   TDistance& dist,
					   TScore const& score_type,
					   Overlap_Library)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Iterator<String<TPair, TPairSpec> >::Type TPairIter;

	// Clear library
	clearVertices(g);

	// Initialization
	TStringSet& str = stringSet(g);
	TSize nseq = length(str);
	__resizeWithRespectToDistance(dist, nseq);

	// String of fragments to combine all pairwise alignments into a multiple alignment
	typedef Fragment<> TFragment;
	typedef String<TFragment> TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;

	// All pairwise alignments
	TPairIter pairIt = begin(pList);
	TPairIter pairItEnd = end(pList);
	for(;pairIt != pairItEnd; ++pairIt) {
		// Make a pairwise string-set
		TStringSet pairSet;
		TId id1 = (value(pairIt)).i1;
		TId id2 = (value(pairIt)).i2;
		assignValueById(pairSet, str, id1);
		assignValueById(pairSet, str, id2);
		
		// Overlap alignment
		TSize from = length(matches);
		globalAlignment(matches, pairSet, score_type, AlignConfig<true,true,true,true>(), Gotoh() );

		// Determine a sequence weight
		TSize matchLen = 0;
		TSize overlapLen = 0;
		TSize alignLen = 0;
		getAlignmentStatistics(matches, pairSet, from, matchLen, overlapLen, alignLen);
		double quality = (double) matchLen / (double) overlapLen;

		// Get only the good overlap alignments
		if (((quality >= 1) && (matchLen >= 8)) ||
			((quality >= 0.8) && (matchLen >= 15))) {

			// Create a corresponding edge
			TSize i = idToPosition(str, id1);
			TSize j = idToPosition(str, id2);

			if (i<j) __getAlignmentStatistics(dist, i, j, nseq, matchLen, quality);
			else __getAlignmentStatistics(dist, j, i, nseq, matchLen, quality);
		} else {
			resize(matches, from);
			
			//// Debug Code
			//TGraph tmp(pairSet);
			//Score<int> st = Score<int>(5,-2,-4,-14);
			//globalAlignment(tmp, pairSet, st, Gotoh() );
			//std::cout << "Match length: " << matchLen << std::endl;
			//std::cout << "Overlap length: " << overlapLen << std::endl;
			//std::cout << "Align length: " << alignLen << std::endl;
			//std::cout << "Quality: " << quality << std::endl;
			//std::cout << tmp << std::endl;
		}
	}

	// Refine all matches, rescore the matches and create multiple alignment
	matchRefinement(matches,str,const_cast<TScore&>(score_type),g);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TDistance, typename TScoreValue, typename TScoreSpec>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TDistance& dist,
					   Score<TScoreValue, TScoreSpec> const& score_type,
					   Overlap_Library)
{
	SEQAN_CHECKPOINT
	typedef typename Id<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type TId;
	String<Pair<TId, TId> > pList;
	selectPairsForLibraryGeneration(g, pList);
	generatePrimaryLibrary(g, pList, dist, score_type, Overlap_Library());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TId, typename TScoreValue, typename TScoreSpec>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   String<Pair<TId, TId> >& pList,
					   Score<TScoreValue, TScoreSpec> const& score_type,
					   Overlap_Library)
{
	SEQAN_CHECKPOINT
	Nothing noth;
	generatePrimaryLibrary(g, pList, noth, score_type, Overlap_Library());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScore>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TScore const& score_type,
					   Overlap_Library)
{
	SEQAN_CHECKPOINT
	typedef typename Id<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type TId;
	String<Pair<TId, TId> > pList;
	selectPairsForLibraryGeneration(g, pList);
	generatePrimaryLibrary(g, pList, score_type, Overlap_Library());
}





}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
