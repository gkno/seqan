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
// Layout-based pair selection
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TBegEndPos, typename TSize, typename TPairList, typename TPos, typename TSpec2>
inline void 
selectPairs(StringSet<TString, TSpec> const& str,
								TBegEndPos const& begEndPos,
								TSize bandwidth,							
								TPairList& pList,
								String<Pair<TPos, TPos>, TSpec2>& dList)
{
	SEQAN_CHECKPOINT
	typedef StringSet<TString, TSpec> TStringSet;
	typedef Pair<TPos, TPos> TDiagPair;
	typedef typename Value<TPairList>::Type TPair;
	typedef typename Iterator<TPairList>::Type TPairIter;
	typedef typename Iterator<TBegEndPos>::Type TBegEndIter;

	// Find all overlapping reads
	TBegEndIter beIt = begin(begEndPos);
	TBegEndIter beItEnd = end(begEndPos);
	TSize index1 = 0;
	for(;beIt != beItEnd; ++beIt, ++index1) {
		TPos posIi1 = (value(beIt)).i1;
		TPos posIi2 = (value(beIt)).i2;
		TBegEndIter beIt2 = beIt;
		++beIt2;
		TSize index2 = index1 + 1;
		for(;beIt2 != beItEnd; ++beIt2, ++index2) {
			TPos posJi1 = (value(beIt2)).i1;
			TPos posJi2 = (value(beIt2)).i2;

			// Diagonal boundaries of the band
			// Initialization values are used if one read is contained in the other 
			TPos diagLow = -1 * (posJi2 - posJi1);
			if (posJi2 < posJi1) diagLow = -1 * (posJi1 - posJi2);
			TPos diagHigh = posIi2 - posIi1;
			if (posIi2 < posIi1) diagHigh = posIi1 - posIi2;

			// If you mix long and short reads and short reads are contained in longer ones you have to increase the band because of long "overlaps" (length of the shorter read)
			TPos radius = (bandwidth + 1) / 2;		// Band width
			TPos lengthDivider = 40;  // Note: overlap / lengthDivider is added to the radius

			// Read orientations
			if (posIi1 < posIi2) {
				// 1) Forward - Forward
				if (posJi1 < posJi2) {
					if ((posJi1 >= posIi2) || (posJi2 <= posIi1)) continue;
					if ((posJi2 < posIi2) && (posJi1 < posIi1)) {
						TPos offset = (posIi1 - posJi1);
						radius += (posJi2 - posJi1 - offset) / lengthDivider;
						if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
						if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
					} else if ((posJi1 > posIi1) && (posJi2 > posIi2)) {
						TPos offset = (posJi1 - posIi1);
						radius += (posIi2 - posIi1 - offset) / lengthDivider;
						if (offset + radius < diagHigh) diagHigh = offset + radius;
						if (offset - radius > diagLow) diagLow = offset - radius;
					} else {
						radius += 40;  // Why???? ToDo???
						if (posIi1 < posJi1) {
							TPos offset = (posJi1 - posIi1);
							radius += (posJi2 - posJi1) / lengthDivider;
							if (offset + radius < diagHigh) diagHigh = offset + radius;
							if (offset - radius > diagLow) diagLow = offset - radius;
						} else {
							TPos offset = (posIi1 - posJi1);
							radius += (posIi2 - posIi1) / lengthDivider;
							if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
							if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
						}
					}
				} else { // 2) Forward - Reverse
					if ((posJi2 >= posIi2) || (posJi1 <= posIi1)) continue;
					if ((posJi1 < posIi2) && (posJi2 < posIi1)) {
						TPos offset = (posIi1 - posJi2);
						radius += (posJi1 - posJi2 - offset) / lengthDivider;
						if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
						if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
					} else if ((posJi2 > posIi1) && (posJi1 > posIi2)) {
						TPos offset = (posJi2 - posIi1);
						radius += (posIi2 - posIi1 - offset) / lengthDivider;
						if (offset + radius < diagHigh) diagHigh = offset + radius;
						if (offset - radius > diagLow) diagLow = offset - radius;
					} else {
						radius += 40;  // Why???? ToDo???
						if (posIi1 < posJi2) {
							TPos offset = (posJi2 - posIi1);
							radius += (posJi1 - posJi2) / lengthDivider;
							if (offset + radius < diagHigh) diagHigh = offset + radius;
							if (offset - radius > diagLow) diagLow = offset - radius;
						} else {
							TPos offset = (posIi1 - posJi2);
							radius += (posIi2 - posIi1) / lengthDivider;
							if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
							if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
						}	
					}
				}
			} else { 
				// 3) Reverse - Forward
				if (posJi1 < posJi2) {
					if ((posJi1 >= posIi1) || (posJi2 <= posIi2)) continue;
					if ((posIi1 > posJi2) && (posIi2 > posJi1)) {
						TPos offset = (posIi2 - posJi1);
						radius += (posJi2 - posJi1 - offset) / lengthDivider;
						if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
						if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
					} else if ((posJi1 > posIi2) && (posJi2 > posIi1)) {
						TPos offset = (posJi1 - posIi2);
						radius += (posIi1 - posIi2 - offset) / lengthDivider;
						if (offset + radius < diagHigh) diagHigh = offset + radius;
						if (offset - radius > diagLow) diagLow = offset - radius;
					} else {
						radius += 40;  // Why???? ToDo???
						if (posIi2 < posJi1) {
							TPos offset = (posJi1 - posIi2);
							radius += (posJi2 - posJi1) / lengthDivider;
							if (offset + radius < diagHigh) diagHigh = offset + radius;
							if (offset - radius > diagLow) diagLow = offset - radius;
						} else {
							TPos offset = (posIi2 - posJi1);
							radius += (posIi1 - posIi2) / lengthDivider;
							if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
							if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
						}
					}
				} else { // 4) Reverse - Reverse
					if ((posJi2 >= posIi1) || (posJi1 <= posIi2)) continue;
					if ((posJi1 < posIi1) && (posJi2 < posIi2)) {
						TPos offset = (posIi2 - posJi2);
						radius += (posJi1 - posJi2 - offset) / lengthDivider;
						if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
						if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
					} else if ((posJi2 > posIi2) && (posJi1 > posIi1)) {
						TPos offset = (posJi2 - posIi2);
						radius += (posIi1 - posIi2 - offset) / lengthDivider;
						if (offset + radius < diagHigh) diagHigh = offset + radius;
						if (offset - radius > diagLow) diagLow = offset - radius;
					} else {
						radius += 40;  // Why???? ToDo???
						if (posIi2 < posJi2) {
							TPos offset = (posJi2 - posIi2);
							radius += (posJi1 - posJi2) / lengthDivider;
							if (offset + radius < diagHigh) diagHigh = offset + radius;
							if (offset - radius > diagLow) diagLow = offset - radius;
						} else {
							TPos offset = (posIi2 - posJi2);
							radius += (posIi1 - posIi2) / lengthDivider;
							if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
							if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
						}
					}
				}
			}

			// Append this pair of reads
			appendValue(pList, TPair(positionToId(str, index1), positionToId(str, index2)));
			appendValue(dList, TDiagPair(diagLow, diagHigh));
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TId, typename TDiagList, typename TScore, typename TSize, typename TSegmentMatches, typename TScoreValues, typename TDistance>
inline void 
appendSegmentMatches(StringSet<TString, TSpec> const& str,
					 String<Pair<TId, TId> > const& pList,
					 TDiagList const& dList,
					 TScore const& score_type,
					 TSize thresholdMatchlength,
					 TSize thresholdQuality,
					 TSegmentMatches& matches,
					 TScoreValues& scores,
					 TDistance& dist,
					 Overlap_Library)
{
	SEQAN_CHECKPOINT
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef String<Pair<TId, TId> > TPairList;
	typedef typename Value<TScoreValues>::Type TScoreValue;
	typedef typename Iterator<TPairList>::Type TPairIter;
	typedef typename Iterator<TDiagList>::Type TDiagIter;

	// Initialization
	double qltThres = (double) thresholdQuality / 100.0;
	TSize nseq = length(str);
	__resizeWithRespectToDistance(dist, nseq);

	// Pairwise alignments
	TPairIter itPair = begin(pList);
	TDiagIter itDiag = begin(dList);
	TPairIter itPairEnd = end(pList);
	for(;itPair != itPairEnd; goNext(itPair), goNext(itDiag)) {
		// Make a pairwise string-set
		TStringSet pairSet;
		TId id1 = (value(itPair)).i1;
		TId id2 = (value(itPair)).i2;
		assignValueById(pairSet, const_cast<StringSet<TString, TSpec>&>(str), id1);
		assignValueById(pairSet, const_cast<StringSet<TString, TSpec>&>(str), id2);

		// Overlap alignment
		TSize from = length(matches);
		TScoreValue myScore = globalAlignment(matches, pairSet, score_type, AlignConfig<true,true,true,true>(), (value(itDiag)).i1, (value(itDiag)).i2, BandedGotoh() );
		TSize to = length(matches);
		
		// Determine a sequence weight
		TSize matchLen = 0;
		TSize overlapLen = 0;
		TSize alignLen = 0;
		getAlignmentStatistics(matches, pairSet, from, matchLen, overlapLen, alignLen);
		double quality = (double) matchLen / (double) overlapLen;

		// Get only the good overlap alignments
		if ((quality >= qltThres) && (matchLen >= thresholdMatchlength)) {

			// Create a corresponding edge
			TSize i = idToPosition(str, id1);
			TSize j = idToPosition(str, id2);
			if (i<j) __getAlignmentStatistics(dist, i, j, nseq, matchLen, quality);
			else __getAlignmentStatistics(dist, j, i, nseq, matchLen, quality);

			// Record the scores
			resize(scores, to);
			typedef typename Iterator<TScoreValues>::Type TScoreIter;
			TScoreIter itScore = begin(scores);
			TScoreIter itScoreEnd = end(scores);
			goFurther(itScore, from);
			for(;itScore != itScoreEnd; ++itScore) value(itScore) = myScore;
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

/*
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
selectPairs(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
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
*/

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
