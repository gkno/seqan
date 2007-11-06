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
// T-Coffee - Determining sequence similarity
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
	
template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix>
inline void 
pairwiseDistances(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				  TMatrix& distanceMatrix)
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
			TValue score = hcsPairwiseScore(g,seq1,seq2);
			
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
			assignValue(distanceMatrix, i*nseq+j, 1 - (TValue) getValue(distanceMatrix, i*nseq+j) / (TValue) maxScore);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
	
template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix>
inline void 
pairwiseLibraryDistances(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
						 TMatrix& distanceMatrix)
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
	fill(distanceMatrix, nseq * nseq, 0);

	// All pairwise alignments
	TEdgeIterator itE(g);
	TCargo max_cargo = 0;
	for(;!atEnd(itE);goNext(itE)) if (cargo(*itE) > max_cargo) max_cargo = cargo(*itE);
	goBegin(itE);
	for(;!atEnd(itE);goNext(itE)) {
		TId id1 = sequenceId(g, sourceVertex(itE));
		TId id2 = sequenceId(g, targetVertex(itE));
		TCargo carg = (TCargo) (((double) cargo(*itE)) / ((double) max_cargo) * 100.0);
		value(distanceMatrix, id1 * nseq + id2) += carg;
		value(distanceMatrix, id2 * nseq + id1) += carg;
		cargo(*itE) = carg;
	}

	// Find the maximum value
	TValue max_similarity = 0;
	for(TSize i = 0; i < length(distanceMatrix); ++i) if (max_similarity < distanceMatrix[i]) max_similarity = distanceMatrix[i];
	max_similarity += 1;

	// Build the distance matrix
	TValue max_distance = 0;
	for(TSize i = 0; i < length(distanceMatrix); ++i) {
		value(distanceMatrix, i) = 1.0 - ((double) getValue(distanceMatrix, i) / max_similarity);
		if (getValue(distanceMatrix, i) > max_distance) max_distance = getValue(distanceMatrix, i);
	}

	// Debug code
	for(TSize row = 0;row < nseq; ++row) {
		TSize count = 0;
		for(TSize col = 0;col < nseq; ++col) {
			if (getValue(distanceMatrix, row*nseq + col) == max_distance) {
				++count;
			}
		}
		std::cout << count << ','; 
	}
	std::cout << std::endl;	
}


template<typename TStringSet, typename TCargo, typename TSpec, typename TOutGraph>
inline void 
highestScoreFirstAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
						   TOutGraph& gOut)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;

	// Initialization
	String<double> distanceMatrix;
	TStringSet& str = stringSet(g);	
	TSize nseq = length(str);
	resize(distanceMatrix, nseq * nseq);

	// Build initial distance matrix
	typedef String<String<TVertexDescriptor> > TSegmentString;
	String<TSegmentString> segmentStr;
	String<bool> active;
	resize(segmentStr, nseq);
	fill(active, nseq, true);
	for(TSize i=0; i<nseq; ++i) _buildLeafString(g, i, value(segmentStr,i));
	for(TSize i=0; i<nseq; ++i) {
		TSize len1 = length(value(segmentStr,i));
		for(TSize j=i+1; j<nseq; ++j) {
			// Align the 2 strings
			TSize len2 = length(value(segmentStr,j));
			double score = hcsPairwiseScore(g,value(segmentStr,i),value(segmentStr,j));
			
			// Normalize by distance
			score /= ((len1 + len2) / 2);
						
			// Remember the value
			assignValue(distanceMatrix, i*nseq+j, score);
		}
	}

	bool oneLeft = true;
	do {
		oneLeft = true;

		// Find highest scoring pair
		TSize index_i = 0;
		TSize index_j = 0;
		double maxScore = 0;
		for(TSize i=0; i<nseq; ++i) {
			if (!active[i]) continue;
			for(TSize j=i+1; j<nseq; ++j) {
				if (!active[j]) continue;
				oneLeft = false;
				//std::cout << getValue(distanceMatrix, i*nseq+j) << ',';
				double tmp;
				if ((tmp = getValue(distanceMatrix, i*nseq+j)) > maxScore) {
					maxScore = tmp;
					index_i = i;
					index_j = j;
				}
			}
			//std::cout << std::endl;
		}
		if (oneLeft) break;
		//std::cout << index_i << ',' << index_j << ':' << maxScore << std::endl;

		// Align this sequence pair
		TSegmentString alignSeq;
		heaviestCommonSubsequence(g,segmentStr[index_i],segmentStr[index_j],alignSeq);
		active[index_j] = false;
		clear(value(segmentStr, index_j));
		segmentStr[index_i] = alignSeq;

		// Recalculate distances
		for(TSize i=0; i<nseq; ++i) {
			if ((!active[i]) || (index_i == i)) continue;
			TSize len1 = length(segmentStr[index_i]);
			TSize len2 = length(segmentStr[i]);
			double score = hcsPairwiseScore(g,segmentStr[index_i],segmentStr[i]);
			
			// Normalize by distance
			score /= ((len1 + len2) / 2);
					
			// Remember the value
			if (index_i < i) assignValue(distanceMatrix, index_i*nseq+i, score);
			else assignValue(distanceMatrix, i*nseq+index_i, score);
		}
	} while (!oneLeft);

	// Create alignment graph
	for(TSize i=0; i<nseq; ++i) {
		if (!active[i]) continue;
		_createAlignmentGraph(g, segmentStr[i], gOut);
		break;
	}
}



//////////////////////////////////////////////////////////////////////////////


template<typename TAlignmentGraph, typename TAlphabet>
inline double 
getSequenceSimilarity(TAlignmentGraph& g,
					  TAlphabet)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TAlignmentGraph>::Type TSize;
	typedef typename Host<TAlignmentGraph>::Type TStringSet;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Infix<TString>::Type TInfix;

	TSize len1 = length(stringSet(g)[0]);
	TSize len2 = length(stringSet(g)[1]);

	double sim = 0.0;
	typedef typename Iterator<TAlignmentGraph, EdgeIterator>::Type TEdgeIterator;
	TEdgeIterator it(g);
	for(;!atEnd(it);++it) {
		typedef typename Iterator<TInfix, Rooted>::Type TInfixIter;
		TInfix inf1 = label(g,sourceVertex(it));
		TInfix inf2 = label(g,targetVertex(it));
		TInfixIter sIt1 = begin(inf1);
		TInfixIter sIt2 = begin(inf2);
		while((!atEnd(sIt1)) || (!atEnd(sIt2))) {
			if ( (TAlphabet) *sIt1  == (TAlphabet) *sIt2) ++sim;
			goNext(sIt1); goNext(sIt2);
		}
	}
	
	// Normalize by distance
	if (len1 > len2) return sim / len2;
	else return sim / len1;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlignmentGraph, typename TValue, typename TSpec, typename TAlphabet>
inline void 
getSequenceSimilarity(TAlignmentGraph& g,
					  String<TValue, TSpec>& sim,
					  TAlphabet)
{
	SEQAN_CHECKPOINT
	typedef String<TValue, TSpec> TSimilarityMatrix;
	typedef typename Size<TAlignmentGraph>::Type TSize;
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
	typedef typename Host<TAlignmentGraph>::Type TStringSet;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Infix<TString>::Type TInfix;

	TStringSet& str = stringSet(g);
	TSize numSeqs = length(str);
	clear(sim);
	fill(sim, numSeqs * numSeqs, 0);

	typedef typename Iterator<TAlignmentGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TInfix, Rooted>::Type TInfixIter;
	TEdgeIterator it(g);
	for(;!atEnd(it);++it) {
		TVertexDescriptor sV = sourceVertex(it);
		TVertexDescriptor tV = targetVertex(it);
		TSize pos1 = idToPosition(str, sequenceId(g, sV));
		TSize pos2 = idToPosition(str, sequenceId(g, tV));
		TInfix inf1 = label(g,sV);
		TInfix inf2 = label(g,tV);
		TInfixIter sIt1 = begin(inf1);
		TInfixIter sIt2 = begin(inf2);
		while((!atEnd(sIt1)) || (!atEnd(sIt2))) {
			if ( (TAlphabet) *sIt1  == (TAlphabet) *sIt2) {
				value(sim, pos1 * numSeqs + pos2) += 1;
				value(sim, pos2 * numSeqs + pos1) += 1;
			}
			goNext(sIt1); goNext(sIt2);
		}
	}
	for(TSize i=0; i<numSeqs; ++i) {
		for(TSize j=i+1; j<numSeqs; ++j) {
			TSize len1 = length(str[i]);
			TSize len2 = length(str[j]);

			if (len1 > len2) {
				value(sim, i * numSeqs + j) /= len2;
				value(sim, j * numSeqs + i) /= len2;
				//value(sim, i * numSeqs + j) /= (len2 * len2);
				//value(sim, j * numSeqs + i) /= (len2 * len2);
			} else {
				value(sim, i * numSeqs + j) /= len1;
				value(sim, j * numSeqs + i) /= len1;
				//value(sim, i * numSeqs + j) /= (len1 * len1);
				//value(sim, j * numSeqs + i) /= (len1 * len1);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlignmentGraph, typename TValue, typename TSpec>
inline void 
getSequenceSimilarity(TAlignmentGraph& g,
					  String<TValue, TSpec>& sim)
{
	SEQAN_CHECKPOINT
	typedef typename Host<TAlignmentGraph>::Type TStringSet;
	typedef typename Value<TStringSet>::Type TString;
	getSequenceSimilarity(g, sim, typename Value<TString>::Type());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFragmentMatches, typename TStringSet, typename TSize1,  typename TSize2, typename TSize3, typename TAlphabet>
inline double 
getSequenceSimilarity(TFragmentMatches& matches,
					  TStringSet& str,
					  TSize1 const from,
					  TSize2 const to,
					  TSize3& alignLength,
					  TAlphabet)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TFragmentMatches>::Type TSize;
	typedef typename Id<TFragmentMatches>::Type TId;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Infix<TString>::Type TInfix;

	TSize sim = 0;
	TSize match_length = 0;
	TSize len1 = length(str[0]);
	TSize len2 = length(str[1]);

	typedef typename Iterator<TInfix, Rooted>::Type TInfixIter;
	for(TSize i = from;i<to;++i) {
		TInfix inf1 = label(matches[i], str, sequenceId(matches[i], 0));
		TInfix inf2 = label(matches[i], str, sequenceId(matches[i], 1));
		TInfixIter sIt1 = begin(inf1);
		TInfixIter sIt2 = begin(inf2);
		while((!atEnd(sIt1)) || (!atEnd(sIt2))) {
			if ( (TAlphabet) *sIt1  == (TAlphabet) *sIt2) {
				++sim;
			}
			goNext(sIt1); goNext(sIt2);
			++match_length;
		}
	}

	alignLength = match_length + (len1 - match_length) + (len2 - match_length);
	if (len1 > len2) return (double) sim / (double) len2;
	else return (double) sim / (double) len1;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFragmentMatches, typename TStringSet, typename TSize, typename TAlphabet>
inline TSize 
getOverlapLength(TFragmentMatches& matches,
				 TStringSet& str,
				 TSize& matchLength,
				 TAlphabet)
{
	SEQAN_CHECKPOINT
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
	for(TSize i = 0;i<length(matches);++i) {
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
	//// Debug code
	//std::cout <<  minId1 << std::endl;
	//std::cout <<  minId2 << std::endl;
	//std::cout <<  len1 << std::endl;
	//std::cout <<  len2 << std::endl;
	//std::cout <<  maxId1 << std::endl;
	//std::cout <<  maxId2 << std::endl;
	TSize alignLength = matchMismatch_length + (len1 - matchMismatch_length) + (len2 - matchMismatch_length);
	return alignLength -  minId1 - minId2 - (len1 + len2 - maxId1 - maxId2);
}


//////////////////////////////////////////////////////////////////////////////
// T-Coffee - Similarity to distance matrix
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TMatrix>
inline void
similarityToDistanceMatrix(TMatrix& mat, 
						   KimuraDistance) 
{
	SEQAN_CHECKPOINT
	
	typedef typename Value<TMatrix>::Type TValue;
	typedef typename Size<TMatrix>::Type TSize;

	// Initialize the mat matrix
	TSize nseq = (TSize) sqrt((double)length(mat));

	// Calculate the mat
	for (TSize row=0;row<nseq;++row) {
		for(TSize col=row+1;col<nseq;++col) {
			// Kimura correction
			TValue val = 0.8 * getValue(mat, row*nseq+col) + 0.2;
			//if (val < 0.2) val = 0.2;
			val = (TValue) log((TValue) val - (1.0 - val * val) / 5.0);
			
			// Assign values
			assignValue(mat, row*nseq+col, val);
			assignValue(mat, col*nseq+row, val);
		}
		assignValue(mat, row*nseq+row, 0);
	}

	//// Debug code
	//for (TSize row=0;row<nseq;++row) {
	//	for(TSize col=0;col<nseq;++col) {
	//		std::cout << getValue(mat, row*nseq+col) << ",";
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TMatrix>
inline void
similarityToDistanceMatrix(TMatrix& mat, 
						   FractionalDistance) 
{
	SEQAN_CHECKPOINT
	
	typedef typename Value<TMatrix>::Type TValue;
	typedef typename Size<TMatrix>::Type TSize;

	// Initialize the mat matrix
	TSize nseq = (TSize) sqrt((double)length(mat));

	// Calculate the mat
	for (TSize row=0;row<nseq;++row) {
		for(TSize col=row+1;col<nseq;++col) {
			// Distance must be between 0 and 1
			TValue val = (1 - getValue(mat, row*nseq+col)) * 100;

			// Assign values
			assignValue(mat, row*nseq+col, val);
			assignValue(mat, col*nseq+row, val);
		}
		assignValue(mat, row*nseq+row, 0);
	}

	//// Debug code
	//for (TSize row=0;row<nseq;++row) {
	//	for(TSize col=0;col<nseq;++col) {
	//		std::cout << getValue(mat, row*nseq+col) << ",";
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
