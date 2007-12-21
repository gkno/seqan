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

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_LIBRARY_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_LIBRARY_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// T-Coffee - Library generation
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TScore>
void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TScore const& score_type,
					   LocalPairwise_Library)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Value<TScore>::Type TScoreValue;

	// Clear graph
	clearVertices(g);

	// Pairwise alignments for all pairs of sequences
	TStringSet& str = stringSet(g);	
	TSize nseq = length(str);

	// String of fragments to combine all pairwise alignments into a multiple alignment
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;
	String<TScoreValue, Block<> > score_values;

	for(TSize i=0; i<nseq; ++i) {
		for(TSize j=i+1; j<nseq; ++j) {
			// Pairwise alignment graph
			TStringSet pairSet;
			assignValueById(pairSet, str, positionToId(str, i));
			assignValueById(pairSet, str, positionToId(str, j));
			typedef Graph<Alignment<TStringSet, unsigned int> > TPairGraph;
			typedef typename VertexDescriptor<TPairGraph>::Type TVD;
			typedef typename Iterator<TPairGraph, EdgeIterator>::Type TEI;
			TPairGraph pGraph(pairSet);

			localAlignment(pGraph, score_type, SmithWatermanClump() );
						
			// Remember the matches and their scores
			TEI it(pGraph);
			for(;!atEnd(it);++it) {
				TVD sV = sourceVertex(it);
				TVD tV = targetVertex(it);
				push_back(matches, TFragment( (unsigned int) sequenceId(pGraph, sV), (unsigned int) fragmentBegin(pGraph,sV), (unsigned int) sequenceId(pGraph, tV),  (unsigned int)  fragmentBegin(pGraph,tV),  (unsigned int)  fragmentLength(pGraph,tV)));
				push_back(score_values, (TScoreValue) cargo(*it));
			}
		}
	}

	// Refine all matches, rescore matches and create multiple alignment
	matchRefinement(matches,str,const_cast< TScore&>(score_type),g);

	// Adapt edge weights, scale weights by significance of local match
	TFragmentStringIter endIt = end(matches);
	unsigned int positionIt = 0;
	for(TFragmentStringIter it = begin(matches); it != endIt; ++it, ++positionIt) {
		TId id1 = sequenceId(*it,0);
		TId id2 = sequenceId(*it,1);
		TSize pos1 = fragmentBegin(*it, id1);
		TSize pos2 = fragmentBegin(*it, id2);
		TSize end1 = pos1 + fragmentLength(*it, id1);
		while(pos1 < end1) {
			TVertexDescriptor p1 = findVertex(g, id1, pos1);
			TVertexDescriptor p2 = findVertex(g, id2, pos2);
			TEdgeDescriptor e = findEdge(g, p1, p2);
			if (e != 0) cargo(e) *= (TCargo) getValue(score_values, positionIt);
			else addEdge(g, p1, p2, (TCargo) getValue(score_values, positionIt));
			pos1 += fragmentLength(g, p1);
			pos2 += fragmentLength(g, p2);
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec, typename TSize>
inline void
__resizeWithRespectToDistance(String<TValue, TSpec>& dist, TSize nseq)
{
	SEQAN_CHECKPOINT
	clear(dist);
	TValue infDist = getInfinity<TValue>(); 
	fill(dist, nseq * nseq, infDist);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TSize>
inline void
__resizeWithRespectToDistance(Graph<Undirected<TCargo, TSpec> >& dist, TSize nseq)
{
	SEQAN_CHECKPOINT
	clear(dist);
	for(TSize i=0;i<nseq; ++i) addVertex(dist);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSize>
inline void
__resizeWithRespectToDistance(Nothing&, TSize)
{
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFragmentMatches, typename TStringSet, typename TValue,  typename TSpec, typename TSize>
inline void 
__getAlignmentStatistics(TFragmentMatches& matches,
						 TStringSet& pairSet,
						 String<TValue, TSpec>& dist,
						 TSize i,
						 TSize j,
						 TSize nseq,
						 TSize from)
{
	SEQAN_CHECKPOINT

	// Determine a sequence weight
	TValue matchLen = 0;
	TValue overlapLen = 0;
	TValue alignLen = 0;
	getAlignmentStatistics(matches, pairSet, from, length(matches),  matchLen, overlapLen, alignLen);
			
	// Calculate sequence similarity
	TValue normalizedSimilarity = (matchLen / overlapLen) * (overlapLen / alignLen);

	// Assign the values
	assignValue(dist, i*nseq+j, 1 - normalizedSimilarity);
	assignValue(dist, j*nseq+i, 1 - normalizedSimilarity);	
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFragmentMatches, typename TStringSet, typename TCargo,  typename TSpec, typename TSize>
inline void 
__getAlignmentStatistics(TFragmentMatches& matches,
						 TStringSet& pairSet,
						 Graph<Undirected<TCargo, TSpec> >& dist,
						 TSize i,
						 TSize j,
						 TSize nseq,
						 TSize from)
{
	SEQAN_CHECKPOINT
		
	// Determine a sequence weight
	TCargo matchLen = 0;
	TCargo overlapLen = 0;
	TCargo alignLen = 0;
	getAlignmentStatistics(matches, pairSet, from, length(matches),  matchLen, overlapLen, alignLen);
			
	// Calculate sequence similarity
	TCargo normalizedSimilarity = (matchLen / overlapLen) * (overlapLen / alignLen);

	addEdge(dist, i, j, 1 - normalizedSimilarity);
}



//////////////////////////////////////////////////////////////////////////////

template<typename TFragmentMatches, typename TStringSet, typename TSize>
inline void 
__getAlignmentStatistics(TFragmentMatches&,
						 TStringSet&,
						 Nothing&,
						 TSize,
						 TSize,
						 TSize,
						 TSize)
{
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TDistance, typename TScore>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TDistance& dist,
					   TScore const& score_type,
					   GlobalPairwise_Library)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;

	// Clear graph
	clearVertices(g);

	// Pairwise alignments for all pairs of sequences
	TStringSet& str = stringSet(g);	
	TSize nseq = length(str);
	__resizeWithRespectToDistance(dist, nseq);
	
	
	// String of fragments to combine all pairwise alignments into a multiple alignment
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;

	// All pairwise alignments
	for(TSize i=0; i<nseq; ++i) {
		for(TSize j=i+1; j<nseq; ++j) {
			// Pairwise alignment, get the matches
			TStringSet pairSet;
			assignValueById(pairSet, str, positionToId(str, i));
			assignValueById(pairSet, str, positionToId(str, j));
			TSize from = length(matches);

			// Alignment
			globalAlignment(matches, pairSet, score_type, Gotoh() );
			
			// Get the alignment statistics
			__getAlignmentStatistics(matches, pairSet, dist, i, j, nseq, from);
		}
	}

	// Refine all matches, rescore the matches and create multiple alignment
	matchRefinement(matches,str,const_cast<TScore&>(score_type),g);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue, typename TScoreSpec>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   Score<TScoreValue, TScoreSpec> const& score_type,
					   GlobalPairwise_Library)
{
	SEQAN_CHECKPOINT
	Nothing noth;
	generatePrimaryLibrary(g, noth, score_type, GlobalPairwise_Library());
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

	// Clear library graph and guide tree
	TStringSet& str = stringSet(g);
	TSize nseq = length(str);
	clearVertices(g);
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

template<typename TStringSet, typename TCargo, typename TSpec, typename TPair, typename TPairSpec, typename TScore>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   String<TPair, TPairSpec> const& pList,
					   TScore const& score_type,
					   Overlap_Library)
{
	SEQAN_CHECKPOINT
	Nothing noth;
	generatePrimaryLibrary(g, pList, noth, score_type, Overlap_Library());
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScore, typename TSize>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TScore const& score_type,
					   TSize ktup,
					   bool overlapping,
					   Kmer_Library)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Value<typename Value<TStringSet>::Type>::Type TAlphabet;
	typedef __int64 TWord;
	typedef String<TWord> TTupelString;
	typedef String<TTupelString> TTupelStringSet;
	
	// Initialization
	clearVertices(g);
	TStringSet& str = stringSet(g);	
	TSize nseq = length(str);
	TWord alphabet_size = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

	// All matching kmers
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	TFragmentString matches;

	// Transform the set of strings into a set of strings of k-tupels
	TTupelStringSet tupSet;
	resize(tupSet, nseq);
	for(TSize k=0;k<nseq;++k) {
		if (overlapping) _getTupelString(str[k], tupSet[k], ktup, TAlphabet());
		else _getNonOverlappingTupelString(str[k], tupSet[k], ktup, TAlphabet());
	}

	// Build for each sequence the q-gram Index and count common hits
	String<TWord, External<> > qIndex;
	TWord qIndexSize = 1;
	for (TWord i=0; i<ktup;++i) qIndexSize *= alphabet_size;
	for(TWord i = 0;i < qIndexSize;++i) push_back(qIndex, (TWord) 0);
	for(TSize k=0;k<nseq-1;++k) {
		TId id1 = positionToId(str, k);
		//std::cout << str[k] << std::endl;

		// First pass: Count occurrences
		for(TWord i = 0;i < qIndexSize;++i) qIndex[i] = 0;
		for(TWord i = 0;i < (TWord) length(tupSet[k]);++i) {
			++qIndex[ tupSet[k][i] ];
		}
		// Build incremental sums
		TWord sum = 0;
		for(TWord i = 0;i < qIndexSize;++i) {
			TWord tmp = qIndex[i];
			qIndex[i] = sum;
			sum += tmp;
		}
		// Second pass: Insert hits
		String<TWord> positions;
		resize(positions, sum, Exact());
		for(TSize i = 0;i < (TSize) length(tupSet[k]);++i) {
			positions[qIndex[ tupSet[k][i] ]] = i;
			++qIndex[ tupSet[k][i] ];
		}
		// Reset pointers
		for(TWord i = qIndexSize-1;i > 0; --i) {
			qIndex[i] = qIndex[i-1];
		}
		qIndex[0] = 0;
		for (TSize k2=k+1; k2<nseq; ++k2) {
			TId id2 = positionToId(str, k2);
			//std::cout << str[k2] << std::endl;

			for(TSize i = 0;i < (TSize) length(tupSet[k2]);++i) {
				TWord current_kmer = tupSet[k2][i];
				TWord amount = 0;
				if (current_kmer == qIndexSize - 1) amount = sum - qIndex[current_kmer];
				else amount = qIndex[current_kmer+1] - qIndex[current_kmer];
				for(TWord pos = qIndex[current_kmer]; pos < qIndex[current_kmer] + amount; ++pos) {
					//std::cout << id1 << ',' << positions[pos] << ':' << id2 << ',' << i << std::endl;
					//std::cout << tupSet[k2][i] << ',' << tupSet[k][positions[pos]] << std::endl;
					if (overlapping) push_back(matches, TFragment(id1, positions[pos], id2, i, ktup));
					else {
						TVertexDescriptor v1 = findVertex(g, id1, (unsigned int) positions[pos] * ktup);
						if (v1 == nilVertex) v1 = addVertex(g, id1, (unsigned int) positions[pos] * ktup, (unsigned int) ktup);
						TVertexDescriptor v2 = findVertex(g, id2, (unsigned int) i * ktup);
						if (v2 == nilVertex) v2 = addVertex(g, id2, (unsigned int) i * ktup, (unsigned int) ktup);
						addEdge(g, v1, v2, (unsigned int) ktup);
					}
				}
			}
		}
	}
	
	if (overlapping) matchRefinement(matches,str,const_cast<TScore&>(score_type),g);
	else {
		for(TSize k=0;k<nseq;++k) {
			TId seqId = positionToId(str, k);
			TSize i = 0;
			TSize len = length(str[k]);
			while(i<len) {
				TVertexDescriptor nextVertex = findVertex(g, seqId, i);
				if (nextVertex == nilVertex) {
					TSize j = i + 1;
					while ((j < len) && (findVertex(g, seqId, j) == nilVertex)) ++j;
					nextVertex = addVertex(g, seqId, i, j-i);
				}
				i += fragmentLength(g, nextVertex);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue, typename TScoreSpec, typename TSize>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   Score<TScoreValue, TScoreSpec> const& score_type,
					   TSize ktup,
					   Kmer_Library)
{
	SEQAN_CHECKPOINT
	generatePrimaryLibrary(g,score_type,ktup,true,Kmer_Library());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue, typename TScoreSpec>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   Score<TScoreValue, TScoreSpec> const& score_type,
					   Kmer_Library)
{
	SEQAN_CHECKPOINT
	generatePrimaryLibrary(g,score_type,3,Kmer_Library());
}

/*

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScore>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TScore const& score_type,
					   MUMPairwise_Library)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;

	// Clear graph
	clearVertices(g);

	// Pairwise alignments for all pairs of sequences
	TStringSet& str = stringSet(g);	
	TSize nseq = length(str);
	
	// String of fragments to combine all pairwise alignments into a multiple alignment
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString, Rooted>::Type TFragmentStringIter;
	TFragmentString matches;

	// All pairwise alignments
	for(TSize i=0; i<nseq; ++i) {
		for(TSize j=i+1; j<nseq; ++j) {
			// Pairwise alignment, get the matches
			TStringSet pairSet;
			TId id1 = positionToId(str, i);
			TId id2 = positionToId(str, j);
			assignValueById(pairSet, str, id1);
			assignValueById(pairSet, str, id2);
			
			// Index creation
			typedef Index<TStringSet> TMyIndex;
			TMyIndex myIndex(pairSet);

			typename Iterator< TMyIndex, MUMs >::Type myMUMiterator(myIndex, 10);
			String< typename SAValue<TMyIndex>::Type > occs;
			while (!atEnd(myMUMiterator)) {
					occs = getOccurrences(myMUMiterator);
					if (getSeqNo(occs[0]) == 0)	push_back(matches, TFragment(id1, getSeqOffset(occs[0]), id2, getSeqOffset(occs[1]), repLength(myMUMiterator)));
					else push_back(matches, TFragment(id1, getSeqOffset(occs[1]), id2, getSeqOffset(occs[0]), repLength(myMUMiterator)));
					++myMUMiterator;
			}
		}
	}
	
	// Refine all matches, rescore the matches and create multiple alignment
	matchRefinement(matches,str,const_cast<TScore&>(score_type),g);
}
*/



//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TGuideTree>
inline void 
subtreeMerging(Graph<Undirected<TCargo, TSpec> >& pairGraph,
			   TGuideTree& guideTree)
{
	SEQAN_CHECKPOINT
	typedef Graph<Undirected<TCargo, TSpec> > TPairGraph;
	typedef typename VertexDescriptor<TPairGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TPairGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TPairGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename EdgeDescriptor<TPairGraph>::Type TEdgeDescriptor;
	typedef typename Size<TPairGraph>::Type TSize;
	
	// Copy distance graph
	TPairGraph reference(pairGraph);
	TCargo infCargo = getInfinity<TCargo>();

	// Collect the vertices for each subtree
	typedef String<TVertexDescriptor> TVertexSet;
	String<TVertexSet> setCollection;
	TSize threshold = 100;	// Maximum amount of members per group

	// Start from the best sequence and spread outwards
	TVertexIterator it(pairGraph);
	TVertexDescriptor maxVertex = 0;
	TSize maxDegree = 0;
	for(;!atEnd(it);++it) {
		TSize tmp = degree(pairGraph, *it);
		if (tmp == 0) {
			std::cout << "Not connected!!!" << std::endl;
			exit(0);
		} else if (tmp > maxDegree) {
			maxDegree = tmp;
			maxVertex = *it;
		}
		addVertex(guideTree);	// For each sequence one vertex in the guide tree
	}
	TVertexSet currentSet;
	appendValue(currentSet, maxVertex);
	while (numVertices(pairGraph) > 1) {
		TOutEdgeIterator itE(pairGraph, maxVertex);
		double maxWeight = 0;
		TVertexDescriptor next_maxVertex = 0;
		for(;!atEnd(itE);++itE) {
			if (cargo(*itE) > maxWeight) {
				maxWeight = cargo(*itE);
				next_maxVertex = targetVertex(itE);
			}
		}
		goBegin(itE);
		for(;!atEnd(itE);++itE) {
			if (targetVertex(itE) != next_maxVertex) {
				TEdgeDescriptor e = findEdge(pairGraph, targetVertex(itE), next_maxVertex);
				if (e == 0) addEdge(pairGraph, targetVertex(itE),  next_maxVertex, cargo(*itE));
				else {
					if (cargo(*itE) < infCargo) cargo(e) += cargo(*itE);
					else cargo(e) = infCargo;
				}
			}
		}
		if ((findEdge(reference, maxVertex, next_maxVertex) == 0) ||
			(length(currentSet) > threshold)) {
			appendValue(setCollection, currentSet);
			clear(currentSet);
		} 
		removeVertex(pairGraph, maxVertex);
		maxVertex = next_maxVertex;
		appendValue(currentSet, next_maxVertex);
	}
	TVertexIterator itLast(pairGraph);
	if (findEdge(reference, maxVertex, *itLast) == 0) {
		appendValue(setCollection, currentSet);
		clear(currentSet);
	}
	appendValue(currentSet, *itLast);
	appendValue(setCollection, currentSet);


	// Calculate a guide tree for each subtree
	String<TVertexDescriptor> singletons;
	String<TVertexDescriptor> roots;
	for(TSize z = 0; z<length(setCollection);++z) {
		String<TCargo> dist;
		TSize len = length(setCollection[z]);
		if (len < 3) {
			for(TSize k = 0; k<len;++k) appendValue(singletons, (setCollection[z])[k]);
			continue;
		}
		fill(dist, len * len, infCargo);
		for(TSize i = 0; i<len; ++i) {
			for(TSize j = i+1; j<len;++j) {
				TEdgeDescriptor e = findEdge(reference, (setCollection[z])[i], (setCollection[z])[j]);
				if (e != 0) {
					if (cargo(e) > 1000) assignValue(dist, i * len + j, 0);
					else assignValue(dist, i * len + j, 1000 - cargo(e));
				}
			}
		}
		TGuideTree gTree;
		upgmaTree(dist, gTree);
		std::map<typename VertexDescriptor<TGuideTree>::Type, TVertexDescriptor> mapping;
		typename Iterator<TGuideTree, BfsIterator>::Type itVert(gTree, getRoot(gTree));
		TVertexDescriptor intVertex = addVertex(guideTree);
		appendValue(roots, intVertex);
		mapping.insert(std::make_pair(getRoot(gTree), intVertex));
		++itVert;
		for(;!atEnd(itVert); ++itVert) {
			if (isLeaf(gTree, *itVert)) {
				addEdge(guideTree, (mapping.find(parentVertex(gTree, *itVert)))->second, (setCollection[z])[*itVert]);
			} else {
				intVertex = addVertex(guideTree);
				mapping.insert(std::make_pair(*itVert, intVertex));
				addEdge(guideTree, (mapping.find(parentVertex(gTree, *itVert)))->second, intVertex);
			}
		}
	}
	for(TSize z = 0; z<length(singletons);++z) {
		appendValue(roots, singletons[z]);
	}
	clear(singletons);

	// Combine the subtrees into a full guide tree
	if (length(roots) == 1) {
		assignRoot(guideTree, roots[0]);
	} else {
		TVertexDescriptor internalVertex = addVertex(guideTree);
		addEdge(guideTree, internalVertex, roots[0]);
		addEdge(guideTree, internalVertex, roots[1]);
		for(TSize z = 2; z<length(roots);++z) {
			TVertexDescriptor newInternalVertex = addVertex(guideTree);
			addEdge(guideTree, newInternalVertex, internalVertex);
			addEdge(guideTree, newInternalVertex, roots[z]);
			internalVertex = newInternalVertex;
		}
		assignRoot(guideTree, internalVertex);
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TNames, typename TPairList, typename TReadLength>
inline void 
selectPairsForLibraryGeneration(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
								TNames& names,
								TPairList& pList,
								TReadLength readLen)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TNames>::Type TString;
	typedef typename Value<TPairList>::Type TPair;
	TStringSet& str = stringSet(g);
	TSize nseq = length(str);
	int threshold = (int) (1.5 * (double) readLen);
	//std::cout << "Selected pairs:" << std::endl;
	for(TSize i=0; i<nseq; ++i) {
		int posI = 0;
		std::stringstream ssStream1(toCString(names[i]));
		ssStream1 >> posI; 
		for(TSize j=i+1; j<nseq; ++j) {
			int posJ = 0;
			std::stringstream ssStream2(toCString(names[j]));
			ssStream2 >> posJ;
			if (((posI > posJ) &&
				(posI - posJ < threshold)) ||
				((posI <= posJ) &&
				(posJ - posI < threshold))) {
					//std::cout << i << '(' << posI << ')' << ',' << j << '(' << posJ << ')' << ';';
					appendValue(pList, TPair(positionToId(str, i), positionToId(str, j)));
			}
		}
	}
	//std::cout << std::endl;
}



//////////////////////////////////////////////////////////////////////////////
// T-Coffee Library Reading / Writing
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar, typename TString>
inline void
_parse_readSequenceData(TFile & file,
						TChar & c,
						TString& str)
{
	SEQAN_CHECKPOINT

	append(str, c);

	// Read word
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isLetter(c)) break;
		append(str, c);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec>
inline void
_readLibrary(TFile & file,
			 Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef typename Value<TFile>::Type TValue;
	
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Id<TGraph>::Type TId;

	TValue c;
	bool seq1ToN = false;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	typedef std::pair<unsigned int, unsigned int> TSeqRes;
	typedef std::map<TSeqRes, TVertexDescriptor> TNodeMap;
	TNodeMap node_map;
	TWord seq1 = 0;
	TWord seq2 = 0;
	while (!_streamEOF(file)) {
		_parse_skipWhitespace(file,c);
		if (_streamEOF(file)) break;
		if (c == '#') {
			c = _streamGet(file);
			_parse_skipWhitespace(file,c);
			seq1 = _parse_readNumber(file, c);
			seq2 = _parse_readNumber(file, c);
			if (empty(g))
				if ((seq1 != 0) && (seq2 != 0)) seq1ToN = true;
			if (seq1ToN) {
				--seq1;
				--seq2;
			}
		} else if (c == '!') {
			_parse_skipLine(file, c);
		} else {
			unsigned int res1 = _parse_readNumber(file, c);
			_parse_skipWhitespace(file,c);
			unsigned int res2 = _parse_readNumber(file, c);
			_parse_skipWhitespace(file,c);
			unsigned int weight = _parse_readNumber(file, c);
			_parse_skipLine(file,c);
		
			// Insert new vertex if necessary
			--res1;
			--res2;
			bool newEdge = false;
			TSeqRes key = std::make_pair(seq1, res1);
			typename TNodeMap::iterator nodePos = node_map.find(key);
			TId id1;
			if (nodePos == node_map.end()) {
				//std::cout << seq1 << ',' << res1 << std::endl;
				id1 = addVertex(g, seq1, res1, 1); 
				node_map.insert(std::make_pair(key, id1));
				newEdge = true;
			} else {
				id1 = nodePos->second;
			}

			key = std::make_pair(seq2, res2);
			nodePos = node_map.find(key);
			TId id2;
			if (nodePos == node_map.end()) {
				//std::cout << seq2 << ',' << res2 << std::endl;
				id2 = addVertex(g, seq2, res2, 1); 
				node_map.insert(std::make_pair(key, id2));
				newEdge = true;
			} else {
				id2 = nodePos->second;
			}

			// Insert a new edge or adapt the weight
			if (newEdge) addEdge(g,id1,id2,weight);
			else {
				TEdgeDescriptor e = findEdge(g,id1,id2);
				if( e == 0 ) addEdge(g,id1,id2,weight);
				else cargo(e) += weight;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec>
void 
read(TFile & file,
	 Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
	 TCoffeeLib) 
{
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename Value<TStringSet>::Type TString;

	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	// Ignore first line
	_parse_skipLine(file, c);
	_parse_skipWhitespace(file, c);
	
	// Read number of sequences
	TWord nSeq = (TWord) _parse_readNumber(file, c);

	// Read sequences
	for(TWord i=0; i<nSeq; ++i) {
		_parse_skipWhitespace(file, c);
		//std::cout << _parse_readIdentifier(file, c) << ", ";
		_parse_readIdentifier(file, c);
		_parse_skipWhitespace(file, c);
		//std::cout << _parse_readNumber(file, c) << ", ";
		_parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		TString sequence;
		_parse_readSequenceData(file,c,sequence);
		//std::cout << sequence << std::endl;
	}
	// Reinitialize the graph, because we changed the sequences
	clearVertices(g);

	// Read library
	_readLibrary(file,g);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TString, typename TSpec, typename TNames>
void 
read(TFile & file,
	 StringSet<TString, TSpec>& oriStr,
	 TNames& names,
	 TCoffeeLib) 
{
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;
	
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	// Ignore first line
	_parse_skipLine(file, c);
	_parse_skipWhitespace(file, c);
	
	// Read number of sequences
	TWord nSeq = (TWord) _parse_readNumber(file, c);
	resize(oriStr, nSeq);

	// Read sequences
	for(TWord i=0; i<nSeq; ++i) {
		_parse_skipWhitespace(file, c);
		appendValue(names, _parse_readIdentifier(file, c));
		_parse_skipWhitespace(file, c);
		_parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		_parse_readSequenceData(file,c,oriStr[i]);
	}
}

/*
//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TNames>
void write(TFile & file, 
		   Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
		   TNames& names,
		   TCoffeeLib) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;\
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Size<TStringSet>::Type TSize;

	_streamWrite(file, "! TC_LIB_FORMAT_01\n");
	TSize len = length(getStringSet(g));
	_streamPutInt(file, len);
	_streamPut(file, '\n');
	for(TSize i=0;i<len;++i) {
		_streamWrite(file, names[i]);
		_streamPut(file, ' ');
		TString str = value(getStringSet(g), i);
		_streamPutInt(file, length(str));
		_streamPut(file, ' ');
		_streamWrite(file, str);
		_streamPut(file, '\n');
	}


	typedef std::pair<unsigned int, unsigned int> TSeq;
	typedef Triple<unsigned int, unsigned int, unsigned int> TData;
	typedef std::multimap<TSeq, TData> TMap;
	TMap m;

	typedef typename Iterator<TGraph, EdgeIterator>::Type TIter;
	TIter it(g);
	for(;!atEnd(it);++it) {
		TVertexDescriptor sV = sourceVertex(it);
		TVertexDescriptor tV = targetVertex(it);
		if (sequenceId(g,sV) > sequenceId(g,tV)) {
			TVertexDescriptor tmp = sV;
			sV = tV;
			tV = tmp;
		}
		for(TSize i = 0; i<fragmentLength(g,sV); ++i) {
			m.insert(std::make_pair(TSeq(sequenceId(g,sV), sequenceId(g,tV)), TData(fragmentBegin(g,sV) + i, fragmentBegin(g,tV) + i, getCargo(*it))));
		}
	}
	TSeq old;
	for(TMap::iterator pos = m.begin();pos!=m.end();++pos) {
		if (old != pos->first) {
			old = pos->first; 
			_streamPut(file, '#');
			_streamPutInt(file, pos->first.first + 1);
			_streamPut(file, ' ');
			_streamPutInt(file, pos->first.second + 1);
			_streamPut(file, '\n');		
		}
		_streamPutInt(file, pos->second.i1 + 1);
		_streamPut(file, ' ');
		_streamPutInt(file, pos->second.i2 + 1);
		_streamPut(file, ' ');
		_streamPutInt(file, pos->second.i3);
		_streamPut(file, '\n');	
	}
	_streamWrite(file, "! SEQ_1_TO_N");
	_streamPut(file, '\n');	
}
*/

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TNames>
void write(TFile & file, 
		   Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
		   TNames& names,
		   TCoffeeLib) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;\
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Size<TStringSet>::Type TSize;

	_streamWrite(file, "! TC_LIB_FORMAT_01\n");
	TSize len = length(getStringSet(g));
	_streamPutInt(file, len);
	_streamPut(file, '\n');
	for(TSize i=0;i<len;++i) {
		_streamWrite(file, names[i]);
		_streamPut(file, ' ');
		TString str = value(getStringSet(g), i);
		_streamPutInt(file, length(str));
		_streamPut(file, ' ');
		_streamWrite(file, str);
		_streamPut(file, '\n');
	}

	typedef typename Iterator<TGraph, EdgeIterator>::Type TIter;
	TIter it(g);
	for(;!atEnd(it);++it) {
		TVertexDescriptor sV = sourceVertex(it);
		TVertexDescriptor tV = targetVertex(it);
		if (sequenceId(g,sV) > sequenceId(g,tV)) {
			TVertexDescriptor tmp = sV;
			sV = tV;
			tV = tmp;
		}
		_streamPut(file, '#');
		_streamPutInt(file, sequenceId(g,sV) + 1);
		_streamPut(file, ' ');
		_streamPutInt(file, sequenceId(g,tV) + 1);
		_streamPut(file, '\n');		
		for(TSize i = 0; i<fragmentLength(g,sV); ++i) {
			_streamPutInt(file, fragmentBegin(g,sV) + i + 1);
			_streamPut(file, ' ');
			_streamPutInt(file, fragmentBegin(g,tV) + i + 1);
			_streamPut(file, ' ');
			_streamPutInt(file, getCargo(*it));
			_streamPut(file, '\n');	
		}
	}
	_streamWrite(file, "! SEQ_1_TO_N");
	_streamPut(file, '\n');	
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TCargo, typename TSpec, typename TNames>
void 
read(TFile & file,
	 Graph<Tree<TCargo, TSpec> >& guideTree,
	 TNames& names,
	 NewickFormat) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Tree<TCargo, TSpec> > TGuideTree;
	typedef typename VertexDescriptor<TGuideTree>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGuideTree>::Type TEdgeDescriptor;
	typedef typename Size<TGuideTree>::Type TSize;
	typedef typename Id<TGuideTree>::Type TId;
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;
	typedef typename Value<TNames>::Type TName;
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

	typedef std::map<TName, TId> TNameToId;
	TNameToId nameToId;
	for(TId i=0; i<length(names);++i) {
		addVertex(guideTree);	// Create the sequence vertices
		nameToId.insert(std::make_pair(names[i], i));
	}

	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	TVertexDescriptor lastVertex = nilVertex;
	TVertexDescriptor lastChild = nilVertex;
	while (!_streamEOF(file)) {
		if (c=='(') {
			if (lastVertex == nilVertex) {
				lastVertex = addVertex(guideTree);
				assignRoot(guideTree, lastVertex);
			} else {
				TVertexDescriptor ch = addChild(guideTree, lastVertex);
				lastVertex = ch;
			}
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);
		} else if (c==')') {
			if (!isRoot(guideTree, lastVertex)) {
				lastChild = lastVertex;
				lastVertex = parentVertex(guideTree, lastVertex);
			} else {
				lastChild = lastVertex;
				lastVertex = nilVertex;
			}
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);
		} else if (c==',') {
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);
		} else if (c==':') {
			c = _streamGet(file);
			cargo(findEdge(guideTree, lastVertex, lastChild)) = _parse_readDouble(file,c);
		} else if (c==';') {
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);
		} else {
			TName tmp = _parse_readIdentifier(file, c);
			//std::cout << tmp << std::endl;
			if (lastVertex == nilVertex) {
				// Tree is rooted at a leaf
				// Create artificial root node
				lastVertex = length(names);
				assignRoot(guideTree, addVertex(guideTree));
				addEdge(guideTree, getRoot(guideTree), lastVertex);
				addEdge(guideTree, getRoot(guideTree), nameToId[tmp]);
			} else {
				addEdge(guideTree, lastVertex, nameToId[tmp]);
			}
			lastChild = nameToId[tmp];
		}
	}
	
	// Root the tree if necessary 
	if (outDegree(guideTree, lastChild) > 2) {
		TVertexDescriptor myRoot = addVertex(guideTree);
		assignRoot(guideTree, myRoot);
		typedef typename Iterator<TGuideTree, OutEdgeIterator>::Type TOutEdgeIterator;
		TOutEdgeIterator it(guideTree, lastChild);
		TVertexDescriptor tV = targetVertex(it);
		TCargo c = cargo(*it);
		removeEdge(guideTree, lastChild, tV);
		addEdge(guideTree, myRoot, tV, c);
		addEdge(guideTree, myRoot, lastChild, (TCargo) 0);
	}

	//std::fstream strm1; // Alignment graph as dot
	//strm1.open("D:\\matches\\test\\tree.dot", std::ios_base::out | std::ios_base::trunc);
	//write(strm1,guideTree,DotDrawing());
	//strm1.close();


	//std::cout << guideTree << std::endl;
}

/*
//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TCargo, typename TSpec, typename TNames>
void 
write(TFile & file,
	  Graph<Tree<TCargo, TSpec> >& guideTree,
	  TNames& names,
	  NewickFormat) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Tree<TCargo, TSpec> > TGuideTree;
	typedef Iterator<TGuideTree, DfsPreorder>::Type TDfsPreorder;

	String<unsigned int> predMap;
	String<unsigned int> distMap;
	breadth_first_search(guideTree, getRoot(guideTree), predMap, distMap);
	unsigned int oldLevel = getProperty(distMap, getRoot(guideTree));
	TDfsPreorder dfsIt(guideTree,getRoot(guideTree));
	for(;!atEnd(dfsIt);goNext(dfsIt)) {
		unsigned int newLevel = getProperty(distMap, *dfsIt);
		if (!isLeaf(guideTree, *dfsIt)) {		
			if (newLevel > oldLevel) {
				_streamPut(file, '(');	
				_streamPut(file, '\n');
			} else {
				unsigned int i = newLevel;
				while (i < oldLevel) {
					_streamPut(file, ')');	
					_streamPut(file, '\n');
					++i;
				}
			}
		} else {
			if (oldLevel == newLevel) {
				_streamPut(file, ',');	
				_streamPut(file, '\n');
				_streamWrite(file, names[*dfsIt]);

			} else {
				if (newLevel > oldLevel) {
					_streamPut(file, '(');	
					_streamPut(file, '\n');
				} else {
					unsigned int i = newLevel;
					while (i < oldLevel) {
						_streamPut(file, ')');	
						_streamPut(file, '\n');
						++i;
					}
				}
				_streamWrite(file, names[*dfsIt]);
			}
		}
		oldLevel = newLevel;
	}
	unsigned int i = 0;
	while (i < oldLevel) {
		_streamPut(file, ')');	
		_streamPut(file, '\n');
		++i;
	}
}
*/

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
