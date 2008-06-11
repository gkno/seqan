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
  $Id: graph_consensus_library.h 1809 2008-03-31 12:57:59Z rausch@PCPOOL.MI.FU-BERLIN.DE $
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
		TSize brPoint2 = length(value(itName));
		TIter nIt = begin(value(itName));
		TIter nItEnd = end(value(itName));
		TSize ind = 0;
		for(;nIt!=nItEnd;++ind, ++nIt) {
			if (value(nIt) == ',') brPoint = ind;
			if (value(nIt) == '[') {
				brPoint2 = ind;
				break;
			}
		}
		TSize posI = 0;
		TSize posJ = 0;
		TString inf1 = infix(value(itName), 0, brPoint);
		TString inf2 = infix(value(itName), brPoint+1, brPoint2);
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
					   String<TPair, TPairSpec> const& begEndPos,
					   TDistance& dist,
					   TScore const& score_type,
					   Overlap_Library)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Iterator<String<TPair, TPairSpec> >::Type TBegEndIter;

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
	TBegEndIter beIt = begin(begEndPos);
	TBegEndIter beItEnd = end(begEndPos);
	TSize index1 = 0;
	for(;beIt != beItEnd; ++beIt, ++index1) {
		int posIi1 = (value(beIt)).i1;
		int posIi2 = (value(beIt)).i2;
		TBegEndIter beIt2 = beIt;
		++beIt2;
		TSize index2 = index1 + 1;
		for(;beIt2 != beItEnd; ++beIt2, ++index2) {
			int posJi1 = (value(beIt2)).i1;
			int posJi2 = (value(beIt2)).i2;

			// Diagonal boundaries of the band
			// Initialization values are used if one read is contained in the other 
			int diagLow = -1 * (posJi2 - posJi1);
			if (posJi2 < posJi1) diagLow = -1 * (posJi1 - posJi2);
			int diagHigh = posIi2 - posIi1;
			if (posIi2 < posIi1) diagHigh = posIi1 - posIi2;
			int radius = 10;		// Band width

			// Read orientations
			if (posIi1 < posIi2) {
				// 1) Forward - Forward
				if (posJi1 < posJi2) {
					if ((posJi1 >= posIi2) || (posJi2 <= posIi1)) continue;
					if ((posJi2 < posIi2) && (posJi1 < posIi1)) {
						if (-1 * (posIi1 - posJi1) - radius > diagLow) diagLow = -1 * (posIi1 - posJi1) - radius;
						if (-1 * (posIi1 - posJi1) + radius < diagHigh) diagHigh = -1 * (posIi1 - posJi1) + radius;
					} else if ((posJi1 > posIi1) && (posJi2 > posIi2)) {
						if ((posJi1 - posIi1) + radius < diagHigh) diagHigh = (posJi1 - posIi1) + radius;
						if ((posJi1 - posIi1) - radius > diagLow) diagLow = (posJi1 - posIi1) - radius;
					}
				} else { // 2) Forward - Reverse
					if ((posJi2 >= posIi2) || (posJi1 <= posIi1)) continue;
					if ((posJi1 < posIi2) && (posJi2 < posIi1)) {
						if (-1 * (posIi1 - posJi2) - radius > diagLow) diagLow = -1 * (posIi1 - posJi2) - radius;
						if (-1 * (posIi1 - posJi2) + radius < diagHigh) diagHigh = -1 * (posIi1 - posJi2) + radius;
					} else if ((posJi2 > posIi1) && (posJi1 > posIi2)) {
						if ((posJi2 - posIi1) + radius < diagHigh) diagHigh = (posJi2 - posIi1) + radius;
						if ((posJi2 - posIi1) - radius > diagLow) diagLow = (posJi2 - posIi1) - radius;
					} 
				}
			} else { 
				// 3) Reverse - Forward
				if (posJi1 < posJi2) {
					if ((posJi1 >= posIi1) || (posJi2 <= posIi2)) continue;
					if ((posIi1 > posJi2) && (posIi2 > posJi1)) {
						if (-1 * (posIi2 - posJi1) - radius > diagLow) diagLow = -1 * (posIi2 - posJi1) - radius;
						if (-1 * (posIi2 - posJi1) + radius < diagHigh) diagHigh = -1 * (posIi2 - posJi1) + radius;
					} else if ((posJi1 > posIi2) && (posJi2 > posIi1)) {
						if ((posJi1 - posIi2) + radius < diagHigh) diagHigh = (posJi1 - posIi2) + radius;
						if ((posJi1 - posIi2) - radius > diagLow) diagLow = (posJi1 - posIi2) - radius;
					}
				} else { // 4) Reverse - Reverse
					if ((posJi2 >= posIi1) || (posJi1 <= posIi2)) continue;
					if ((posJi1 < posIi1) && (posJi2 < posIi2)) {
						if (-1 * (posIi2 - posJi2) - radius > diagLow) diagLow = -1 * (posIi2 - posJi2) - radius;
						if (-1 * (posIi2 - posJi2) + radius < diagHigh) diagHigh = -1 * (posIi2 - posJi2) + radius;
					} else if ((posJi2 > posIi2) && (posJi1 > posIi1)) {
						if ((posJi2 - posIi2) + radius < diagHigh) diagHigh = (posJi2 - posIi2) + radius;
						if ((posJi2 - posIi2) - radius > diagLow) diagLow = (posJi2 - posIi2) - radius;
					}
				}
			}


			// Make a pairwise string-set
			TStringSet pairSet;
			TId id1 = positionToId(str, index1);
			TId id2 = positionToId(str, index2);
			assignValueById(pairSet, str, id1);
			assignValueById(pairSet, str, id2);

			// Overlap alignment
			TSize from = length(matches);
			globalAlignment(matches, pairSet, score_type, AlignConfig<true,true,true,true>(), diagLow, diagHigh, BandedGotoh() );
			//globalAlignment(matches, pairSet, score_type, AlignConfig<true,true,true,true>(), Gotoh() );

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
