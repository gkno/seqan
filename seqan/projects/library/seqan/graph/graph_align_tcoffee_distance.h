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

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_DISTANCE_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_DISTANCE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// T-Coffee - Distance matrix calculation
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix>
inline void 
getDistanceMatrix(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				  TMatrix& distanceMatrix,
				  LibraryDistance)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Value<TMatrix>::Type TValue;

	// Initialization
	clear(distanceMatrix);
	TStringSet& str = stringSet(g);	
	TSize nseq = length(str);
	resize(distanceMatrix, nseq * nseq);

	// All pairwise alignments
	typedef String<String<TVertexDescriptor> > TSegmentString;
	TValue maxScore = 0;
	for(TSize i=0; i<nseq; ++i) {
		TSegmentString seq1;
		TSize len1 = length(str[i]);
		_buildLeafString(g, i, seq1);
		for(TSize j=i+1; j<nseq; ++j) {
			// Align the 2 strings
			TSegmentString seq2;
			TSize len2 = length(str[j]);
			_buildLeafString(g, j, seq2);
			TSegmentString alignSeq;
			TValue score = heaviestCommonSubsequence(g,seq1,seq2,alignSeq);
			
			// Normalize by distance
			if (len1 > len2) score /= len1;
			else score /= len2;
			if (score > maxScore) maxScore = score;
			
			// Remember the value
			assignValue(distanceMatrix, i*nseq+j, score);
		}
	}

	// Normalize values
	for(TSize i=0; i<nseq; ++i) {
		for(TSize j=i+1; j<nseq; ++j) {
			TValue normalizedScore = (TValue) 100 - (TValue) ( ( (double) getValue(distanceMatrix, i*nseq+j) / (double) maxScore ) * 100);
			assignValue(distanceMatrix, i*nseq+j, normalizedScore);
			assignValue(distanceMatrix, j*nseq+i, normalizedScore);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix, typename TSize, typename TAlphabet>
inline void 
getDistanceMatrix(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				  TMatrix& distanceMatrix,
				  TSize ktup,
				  TAlphabet,
				  KmerDistance)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TMatrix>::Type TValue;
	typedef typename Iterator<TMatrix>::Type TMatrixIterator;

	getKmerSimilarityMatrix(stringSet(g), distanceMatrix, ktup, TAlphabet());
	
	// Similarity to distance conversion
	TMatrixIterator matIt = begin(distanceMatrix);
	TMatrixIterator endMatIt = end(distanceMatrix);
	for(;matIt != endMatIt;++matIt) value(matIt) = (1 - (*matIt)) * 100;

	//// Debug code
	//for (TSize row=0;row<nseq;++row) {
	//	for(TSize col=0;col<nseq;++col) {
	//		std::cout << getValue(distanceMatrix, row*nseq+col) << ",";
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix, typename TSize>
inline void 
getDistanceMatrix(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				  TMatrix& distanceMatrix,
				  TSize ktup,
				  KmerDistance)
{
	SEQAN_CHECKPOINT
	getDistanceMatrix(g, distanceMatrix, ktup, typename Value<typename Value<TStringSet>::Type>::Type(), KmerDistance() );
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix>
inline void 
getDistanceMatrix(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				  TMatrix& distanceMatrix,
				  KmerDistance)
{
	SEQAN_CHECKPOINT
	getDistanceMatrix(g, distanceMatrix, 3, KmerDistance() );
}


//////////////////////////////////////////////////////////////////////////////
// T-Coffee - Determining alignment statistics
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TFragmentMatches, typename TStringSet, typename TPos, typename TSize1, typename TAlphabet>
inline void 
getAlignmentStatistics(TFragmentMatches const& matches,
					   TStringSet& str,
					   TPos const from,
					   TPos const to,
					   TSize1& matchLength,	// Number of identical characters
					   TSize1& overlapLength,	// Number of character in overlapping segments (with mismatches and gaps)
					   TSize1& alignLength,	// Length of the alignment
					   TAlphabet)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TFragmentMatches>::Type TSize;
	typedef typename Value<TFragmentMatches>::Type TFragment;
	typedef typename Id<TFragmentMatches>::Type TId;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Infix<TString>::Type TInfix;

	matchLength = 0;
	TSize len1 = length(str[0]);
	TSize len2 = length(str[1]);

	typedef typename Iterator<TInfix, Rooted>::Type TInfixIter;
	TSize minId1 = len1 + len2;
	TSize minId2 = len1 + len2;
	TSize maxId1 = 0;
	TSize maxId2 = 0;
	TSize matchMismatch_length = 0;

	for(TSize i = from;i<to;++i) {
		TId id1 = sequenceId(matches[i], 0);
		TId id2 = sequenceId(matches[i], 1);
		if (fragmentBegin(matches[i], id1) < minId1) minId1 = fragmentBegin(matches[i], id1);
		if (fragmentBegin(matches[i], id2) < minId2) minId2 = fragmentBegin(matches[i], id2);
		if (fragmentBegin(matches[i], id1) + fragmentLength(matches[i], id1) > maxId1) maxId1 = fragmentBegin(matches[i], id1) + fragmentLength(matches[i], id1);
		if (fragmentBegin(matches[i], id2) + fragmentLength(matches[i], id2) > maxId2) maxId2 = fragmentBegin(matches[i], id2) + fragmentLength(matches[i], id2);
		TInfix inf1 = label(matches[i], str, id1);
		TInfix inf2 = label(matches[i], str, id2);
		TInfixIter sIt1 = begin(inf1);
		TInfixIter sIt2 = begin(inf2);
		while((!atEnd(sIt1)) || (!atEnd(sIt2))) {
			if ( (TAlphabet) *sIt1  == (TAlphabet) *sIt2) {
				++matchLength;
			}
			goNext(sIt1); goNext(sIt2);
			++matchMismatch_length;
		}
	}
	alignLength = matchMismatch_length + (len1 - matchMismatch_length) + (len2 - matchMismatch_length);
	overlapLength = alignLength -  minId1 - minId2 - (len1 + len2 - maxId1 - maxId2);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFragmentMatches, typename TStringSet, typename TPos, typename TSize>
inline void 
getAlignmentStatistics(TFragmentMatches& matches,
					   TStringSet& str,
					   TPos from,
					   TPos to,
					   TSize& matchLength,
					   TSize& overlapLength,
					   TSize& alignLength)
{
	SEQAN_CHECKPOINT
	getAlignmentStatistics(matches, str, from, to, matchLength, overlapLength, alignLength, typename Value<typename Value<TStringSet>::Type>::Type() );
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFragmentMatches, typename TStringSet, typename TSize>
inline void 
getAlignmentStatistics(TFragmentMatches& matches,
					   TStringSet& str,
					   TSize& matchLength,
					   TSize& overlapLength,
					   TSize& alignLength)
{
	SEQAN_CHECKPOINT
	getAlignmentStatistics(matches, str, (TSize) 0, (TSize) length(matches), matchLength, overlapLength, alignLength, typename Value<TStringSet>::Type() );
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
