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
  $Id: graph_align_smith_waterman_clump.h 1764 2008-03-07 10:28:01Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_SMITH_WATERMAN_ISLAND_H
#define SEQAN_HEADER_GRAPH_SMITH_WATERMAN_ISLAND_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Smith Waterman Alignment with Islands, affine gap cost
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TTrace, typename TStringSet, typename TScore, typename TIslandMax, typename TIslandIndex>
inline typename Value<TIslandMax>::Type
_align_smith_waterman_island(TTrace& trace,
							 TStringSet const& str,
							 TScore const & sc,
							 TIslandMax& islandMax,
							 TIslandIndex& islandIndex)
{
	SEQAN_CHECKPOINT
	// TraceBack values for Smith Waterman
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2, Stop = 12};

	// The DP Matrices
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Size<TString>::Type TSize;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef String<TScoreValue> TColumn;
	typedef typename Iterator<TColumn>::Type TColumnIter;
	TColumn mat;	// The DP Matrix for gaps from the left
	TColumn horizontal;	// The DP Matrix for gaps from the top
	TScoreValue vert = 0;

	// The Island Matrices
	typedef typename Value<TIslandIndex>::Type TIndexPair;
	typedef String<TSize> TIslandColumn;
	typedef typename Iterator<TIslandColumn>::Type TIslandColumnIter;
	TIslandColumn islandMat;
	TIslandColumn islandHorizontal;
	TSize islandVert = 0;

	// Initialization
	TString const& str1 = str[0];
	TString const& str2 = str[1];		
	TSize len1 = length(str1);
	TSize len2 = length(str2);
	TScoreValue gap = scoreGapExtend(sc);
	TScoreValue gapOpen = scoreGapOpen(sc);
	fill(mat, (len2 + 1), 0);
	fill(horizontal, len2, gapOpen);
	resize(trace, len1*len2);
	typedef typename Value<TTrace>::Type TTraceValue;
	TTraceValue tvMat=0, tvHorizontal=0, tvVertical=0;
	fill(islandMat, (len2 + 1), 0);
	fill(islandHorizontal, len2, 0);
	appendValue(islandMax, 0); // First value is the number of islands
	appendValue(islandIndex, TIndexPair(0,0));

	// Classical DP
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	TTraceIter it = begin(trace, Standard() );
	for(TSize col = 1; col <= len1; ++col) {
		TScoreValue diagValMat = 0;
		vert = gapOpen;
		islandVert = 0;
		TColumnIter horizontalIt = begin(horizontal);
		TIslandColumnIter islandHorizontalIt = begin(islandHorizontal);
		TColumnIter previousMatIt = begin(mat);
		TIslandColumnIter previousIslandMatIt = begin(islandMat);
		TColumnIter matIt = begin(mat); ++matIt;
		TIslandColumnIter islandMatIt = begin(islandMat); ++ islandMatIt;
		for(TSize row = 1; row <= len2; ++row) {
			// Get the new maximum for vertical
			if (value(previousMatIt) + gapOpen > vert + gap) {
				vert = value(previousMatIt) + gapOpen;
				islandVert = value(previousIslandMatIt);
				tvVertical = (Byte) Diagonal;
			} else {
				vert = vert + gap;
				tvVertical = (Byte) Vertical;
			}
			if (vert < 0) islandVert = 0;

	
			// Get the new maximum for horizontal
			if (value(matIt) + gapOpen > value(horizontalIt) + gap) {
				value(horizontalIt) = value(matIt) + gapOpen;
				value(islandHorizontalIt) = value(previousIslandMatIt);
				tvHorizontal = (Byte) Diagonal;
			} else {
				value(horizontalIt) = value(horizontalIt) + gap;
				tvHorizontal = (Byte) Horizontal;
			}
			if (value(horizontalIt) < 0) value(islandHorizontalIt) = 0;
	
			// Get the new maximum for mat
			
			TScoreValue tmp = diagValMat + score(const_cast<TScore&>(sc), str1[col-1], str2[row-1]);
			tvMat = (Byte) Diagonal;
			if (0 >= tmp) {
				tmp = 0;
				value(islandMatIt) = 0;
				tvMat = (Byte) Stop;
			} else {
				if (value(previousIslandMatIt) == 0) {
					value(islandMatIt) = value(islandMax, 0) + 1;
				} else {
					value(islandMatIt) = value(previousIslandMatIt);
				}
			}
			if (vert > tmp) {
				tmp = vert;
				value(islandMatIt) = islandVert;
				tvMat = (Byte) Vertical;
			}
			if (value(horizontalIt) > tmp) {
				tmp = value(horizontalIt);
				value(islandMatIt) = value(islandHorizontalIt);
				tvMat = (Byte) Horizontal;
			}

			// Assign the new diagonal value
			diagValMat = value(matIt);
			value(matIt) = tmp;

			// Do we have a new Island?
			if (value(islandMatIt) == value(islandMax, 0) + 1) {
				value(islandMax, 0) = value(islandMax, 0) + 1;
				appendValue(islandMax, value(matIt));
				appendValue(islandIndex, TIndexPair(row, col));
			} else if ((value(islandMatIt) != 0) &&
						(value(matIt) > value(islandMax, value(islandMatIt)))) {
				value(islandMax, value(islandMatIt)) = value(matIt);	
				value(islandIndex, value(islandMatIt)).i1 = row;
				value(islandIndex, value(islandMatIt)).i2 = col;
			}

			// Assign the right trace value
			if (tvMat == (Byte) Stop) {
				assignValue(it, 12);
			} else if (tvMat == (Byte) Diagonal) {
				if (tvHorizontal == (Byte) Diagonal) {
					if (tvVertical == (Byte) Diagonal) assignValue(it, 0);
					else assignValue(it, 1);
				} else if (tvHorizontal == (Byte) Horizontal) {
					if (tvVertical == (Byte) Diagonal) assignValue(it, 2);
					else assignValue(it, 3);
				}
			} else if (tvMat == (Byte) Horizontal) {
				if (tvHorizontal == (Byte) Diagonal) {
					if (tvVertical == (Byte) Diagonal) assignValue(it, 4);
					else assignValue(it, 5);
				} else if (tvHorizontal == (Byte) Horizontal) {
					if (tvVertical == (Byte) Diagonal) assignValue(it, 6);
					else assignValue(it, 7);
				}
			} else if (tvMat == (Byte) Vertical) {
				if (tvHorizontal == (Byte) Diagonal) {
					if (tvVertical == (Byte) Diagonal) assignValue(it, 8);
					else assignValue(it, 9);
				} else if (tvHorizontal == (Byte) Horizontal) {
					if (tvVertical == (Byte) Diagonal) assignValue(it, 10);
					else assignValue(it, 11);
				}
			}

			// ForwardPointers
			goNext(matIt); goNext(previousMatIt);
			goNext(islandMatIt); goNext(previousIslandMatIt);
			goNext(horizontalIt);
			goNext(islandHorizontalIt);
			goNext(it);
		}
	}

	return value(islandMax, 0);
}


//////////////////////////////////////////////////////////////////////////////


template<typename TAlign, typename TStringSet, typename TPropertyMap, typename TScore, typename TSize1>
inline void
_localAlignment(TAlign& align,
				TStringSet const& str,
				TPropertyMap& propMap,
				TScore const& sc,
				TSize1 numAlignments,
				SmithWatermanIsland)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Value<TPropertyMap>::Type TRankScorePair;
	typedef Pair<TSize, TSize> TIndexPair;

	// Compute the DP Matrix
	String<bool> forbidden;
	String<TraceBackGotoh> trace;
	String<TScoreValue> islandMax; 
	String<TIndexPair> islandIndex;
	TScoreValue numIslands = _align_smith_waterman_island(trace, str, sc, islandMax, islandIndex);

	// Take the numAlignments best local alignments
	typedef typename Iterator<String<TScoreValue> >::Type TIslandMaxIter;
	TIslandMaxIter islandMaxIt = begin(islandMax);
	++islandMaxIt; // Ignore first value
	TIslandMaxIter islandMaxItEnd = end(islandMax);
	typedef std::multiset<std::pair<TScoreValue, TSize>, std::greater<std::pair<TScoreValue, TSize> > > TScoreValueSet;
	TScoreValueSet scores;
	TSize myPos = 1;
	for(;islandMaxIt != islandMaxItEnd; ++islandMaxIt, ++myPos) scores.insert(std::make_pair(value(islandMaxIt), myPos));
	typename TScoreValueSet::const_iterator scoresIt = scores.begin();
	typename TScoreValueSet::const_iterator scoresItEnd = scores.end();
	for(TSize posit = 0;((scoresIt != scoresItEnd) && (posit < numAlignments)); ++scoresIt, ++posit) {
		TSize from = length(align);
		_align_smith_waterman_trace(align, str, trace, 0, (value(islandIndex, scoresIt->second)).i1, (value(islandIndex, scoresIt->second)).i2, forbidden);
		TSize to = length(align);
		resize(propMap, to);
		for(TSize walk = from; walk < to; ++walk) value(propMap, walk) = TRankScorePair(posit, value(islandMax, scoresIt->second));
	}
}





}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
