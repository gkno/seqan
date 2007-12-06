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
	//for(unsigned int i= 0; i<length(str[1]);++i) {
	//	for(unsigned int j= 0; j<length(str[0]);++j) {
	//		std::cout << (unsigned int) getValue(forbidden, j*length(str[1]) + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;
	
	// Follow the trace and create the alignment
	_align_smith_waterman_trace(align, str, trace, initialDir, best_row, best_col, forbidden);
	
	return maxScore;
}



//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScore, typename TLimitCount>
inline typename Value<TScore>::Type
_localAlignment(TAlign& align,
				TStringSet& str,
				TScore const& sc,
				TLimitCount limit_count,
				SmithWatermanClump)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename VertexDescriptor<TAlign>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TAlign>::Type TEdgeDescriptor;
	typedef typename Id<TAlign>::Type TId;
	typedef typename Cargo<TAlign>::Type TCargo;
	typedef typename Iterator<TAlign, EdgeIterator>::Type TEdgeIterator;
  
	// Score of the best alignment
	TScoreValue maxScore = 0;

	// For clumpping remember the used positions
	TSize len0 = length(str[0]);
	TSize len1 = length(str[1]);
	String<bool> forbidden;
	fill(forbidden, len0 * len1, false);

	// String of fragments
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString, Rooted>::Type TFragmentStringIter;
	TFragmentString matches;
	String<TScoreValue, Block<> > score_values;

	// Stop looking for local alignments, if there score is too low
	TScoreValue local_score = 0;
	TSize count = 0;
	do {
		// Create the local alignment
		local_score = _localAlignment(matches, str, forbidden, sc, SmithWatermanClump());
		if (local_score > maxScore) maxScore = local_score;

		// Remember the confidence in these matches (Score value)
		TSize diff = length(matches) - length(score_values);
		for(TSize k = 0; k<diff; ++k) push_back(score_values, local_score);

		++count;
	} while ((local_score > 0.5 * maxScore) && (count < limit_count));

	// Refine all matches and create multiple alignment
	matchRefinement(matches,str,align);

	// Adapt edge weights
	TFragmentStringIter endIt = end(matches);
	for(TFragmentStringIter it = begin(matches); it != endIt; ++it) {
		TId id1 = sequenceId(*it,0);
		TId id2 = sequenceId(*it,1);
		TSize pos1 = fragmentBegin(*it, id1);
		TSize pos2 = fragmentBegin(*it, id2);
		TSize end1 = pos1 + fragmentLength(*it, id1);
		while(pos1 < end1) {
			TVertexDescriptor p1 = findVertex(align, id1, pos1);
			TVertexDescriptor p2 = findVertex(align, id2, pos2);
			TEdgeDescriptor e = findEdge(align, p1, p2);
			cargo(e) = (TCargo) getValue(score_values, position(it));
			SEQAN_TASSERT(fragmentLength(align, p1) == fragmentLength(align, p2))
			pos1 += fragmentLength(align, p1);
			pos2 += fragmentLength(align, p2);
		}
	}

	return maxScore;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScore>
inline typename Value<TScore>::Type
_localAlignment(TAlign& align,
				TStringSet& str,
				TScore const& sc,
				SmithWatermanClump)
{
	SEQAN_CHECKPOINT
	return _localAlignment(align,str,sc,4,SmithWatermanClump());
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
