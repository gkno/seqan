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

	TSize nseq = length(strSet);
	reserve(begEndPos, nseq);

	for(TSize i=0; i<nseq; ++i) {
		TSize brPoint = 0;
		for(TSize ind = 0; ind<length(names[i]); ++ind) {
			if ((names[i])[ind] == ',') {
				brPoint = ind;
				break;
			}
		}
		TSize posI = 0;
		TSize posJ = 0;
		TString inf1 = infix(names[i], 0, brPoint);
		TString inf2 = infix(names[i], brPoint+1, length(names[i]));
		std::stringstream ssStream1(toCString(inf1));
		ssStream1 >> posI; 
		std::stringstream ssStream2(toCString(inf2));
		ssStream2 >> posJ;
		appendValue(begEndPos, TPair(posI, posJ));
		if (posI > posJ) reverseComplementInPlace(strSet[i]);
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
	TStringSet& str = stringSet(g);
	TSize nseq = length(str);
	int threshold = (int) (1.5 * (double) readLen);
	for(TSize i=0; i<nseq; ++i) {
		int posI = begEndPos[i].i1;
		if (posI > (int) begEndPos[i].i2) posI = begEndPos[i].i2;
		for(TSize j=i+1; j<nseq; ++j) {
			int posJ = begEndPos[j].i1;
			if (posJ > (int) begEndPos[j].i2) posJ = begEndPos[j].i2;
			if (((posI > posJ) &&
				(posI - posJ < threshold)) ||
				((posI <= posJ) &&
				(posJ - posI < threshold))) {
					appendValue(pList, TPair(positionToId(str, i), positionToId(str, j)));
			}
		}
	}
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

	// Clear library
	clearVertices(g);

	// Initialization
	TStringSet& str = stringSet(g);
	TSize nseq = length(str);
	__resizeWithRespectToDistance(dist, nseq);

	// String of fragments to combine all pairwise alignments into a multiple alignment
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;

	// All pairwise alignments
	TSize amountOfPairs = length(pList);
	for(TSize k=0; k<amountOfPairs; ++k) {
		// Make a pairwise string-set
		TStringSet pairSet;
		TId id1 = (pList[k]).i1;
		TId id2 = (pList[k]).i2;
		assignValueById(pairSet, str, id1);
		assignValueById(pairSet, str, id2);
		
		// Overlap alignment
		TFragmentString pairwise_matches;
		globalAlignment(pairwise_matches, pairSet, score_type, AlignConfig<true,true,true,true>(), Gotoh() );
			
		// Determine a sequence weight
		TSize matchLen = 0;
		TSize overlapLen = 0;
		TSize alignLen = 0;
		getAlignmentStatistics(pairwise_matches, pairSet, matchLen, overlapLen, alignLen);
		double quality = (double) matchLen / (double) overlapLen;

		// Get only the good overlap alignments
		if (((quality >= 1) && (matchLen >= 8)) ||
			((quality >= 0.85) && (matchLen >= 15))) {
			TSize len = length(pairwise_matches);	
			for(TSize z = 0; z<len; ++z) {
				//TGraph tmp(pairSet);
				//globalAlignment(tmp, pairSet, score_type, AlignConfig<true,true,true,true>(), Gotoh() );
				//std::cout << "Match length: " << matchLen << std::endl;
				//std::cout << "Overlap length: " << overlapLen << std::endl;
				//std::cout << "Align length: " << alignLen << std::endl;
				//std::cout << "Quality: " << quality << std::endl;
				//std::cout << tmp << std::endl;

				push_back(matches, pairwise_matches[z]);
			}
			// Create a corresponding edge
			TSize i = idToPosition(str, id1);
			TSize j = idToPosition(str, id2);

			if (i<j) __getAlignmentStatistics(dist, i, j, nseq, matchLen, quality);
			else __getAlignmentStatistics(dist, j, i, nseq, matchLen, quality);
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
