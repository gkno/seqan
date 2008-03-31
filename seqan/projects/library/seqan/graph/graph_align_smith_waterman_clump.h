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

#ifndef SEQAN_HEADER_GRAPH_SMITH_WATERMAN_CLUMP_H
#define SEQAN_HEADER_GRAPH_SMITH_WATERMAN_CLUMP_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Smith Waterman Alignment with Clumping, affine gap cost
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TForbidden, typename TScore>
inline typename Value<TScore>::Type
_localAlignment(TAlign& align,
				TStringSet& str,
				TForbidden& forbidden,
				TScore const& sc,
				SmithWatermanClump)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	  
	TScoreValue maxScore;
	TSize best_row = 0;
	TSize best_col = 0;
	
	// Trace
	String<TraceBackGotoh> trace;
	TraceBackGotoh initialDir;

	// Create the trace
	maxScore = _align_smith_waterman(trace, str, sc, initialDir, best_row, best_col, forbidden);	

	//// Debug code
	//for(TSize i= 0; i<length(str[1]);++i) {
	//	for(TSize j= 0; j<length(str[0]);++j) {
	//		std::cout << (TSize) getValue(forbidden, j*length(str[1]) + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;
	
	// Follow the trace and create the alignment
	_align_smith_waterman_trace(align, str, trace, initialDir, best_row, best_col, forbidden);
	
	return maxScore;
}



//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TPropertyMap, typename TScore, typename TSize1>
inline void
_localAlignment(TAlign& align,
				TStringSet const& str,
				TPropertyMap& propMap,
				TScore const& sc,
				TSize1 numAlignments,
				SmithWatermanClump)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TPropertyMap>::Type TRankScorePair;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
  
	// For clumpping remember the used positions
	TSize len0 = length(str[0]);
	TSize len1 = length(str[1]);
	String<bool> forbidden;
	fill(forbidden, len0 * len1, false);

	// Stop looking for local alignments, if there score is too low
	TScoreValue local_score = 0;
	TScoreValue maxScore = 0;
	TSize count = 0;
	do {
		// Create the local alignment
		TSize from = length(align);
		local_score = _localAlignment(align, str, forbidden, sc, SmithWatermanClump());
		TSize to = length(align);
		if (local_score > maxScore) maxScore = local_score;

		resize(propMap, to);
		for(TSize k = from; k<to; ++k) value(propMap, k) = TRankScorePair(count, local_score);
		++count;
	} while ((local_score > 0.5 * maxScore) && (count < (TSize) numAlignments));
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
