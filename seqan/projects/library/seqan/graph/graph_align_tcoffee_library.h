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
	
template<typename TStringSet, typename TCargo, typename TSpec>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   Lcs_Library)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Value<TStringSet>::Type TString;

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

	for(TSize i=0; i<nseq-1; ++i) {
		for(TSize j=i+1; j<nseq; ++j) {
			// Lcs between first and second string
			String<std::pair<unsigned int, unsigned int>, Block<> > pos1;
			longestCommonSubsequence(str[i], str[j], 1000, pos1);

			// Get the significant matches in all 3 sequences
			TSize lenMatch = 1;						
			int last = length(pos1)-1;		
			TSize iBegin = pos1[last].first;
			TSize jBegin = pos1[last].second;
			for(int z = last - 1; z>=0; --z) {
				if ((pos1[z].first == pos1[z+1].first + 1) &&
					(pos1[z].second == pos1[z+1].second + 1)) 
				{
					++lenMatch;
				} else {
					//// Debug code
					//typedef typename Infix<TString>::Type TInfix;
					//TInfix inf1 = infix(str1,iBegin, iBegin + lenMatch);
					//TInfix inf2 = infix(str2,jBegin, jBegin + lenMatch);
					//std::cout << inf1 << std::endl;
					//std::cout << inf2 << std::endl;
						
					push_back(matches, TFragment(positionToId(str, i),iBegin,positionToId(str, j),jBegin,lenMatch));
					lenMatch = 1;
					iBegin = pos1[z].first;
					jBegin = pos1[z].second;
				}
			}
			// Process last match
			push_back(matches, TFragment(positionToId(str, i),iBegin,positionToId(str, j),jBegin,lenMatch));
		}
	}

	// Clear graph
	clearVertices(g);

	// Refine all matches and create multiple alignment
	matchRefinement(matches,str,g);
}

//////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TPairList, typename TScore>
void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TPairList& pList,
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

	// String of fragments to combine all pairwise alignments into a multiple alignment
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;

	// Pairwise alignments
	TSize amountOfPairs = length(pList);
	for(TSize k=0; k<amountOfPairs; ++k) {
		// Make a pairwise string-set
		TStringSet pairSet;
		TId id1 = (pList[k]).i1;
		TId id2 = (pList[k]).i2;
		assignValueById(pairSet, str, id1);
		assignValueById(pairSet, str, id2);

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
		}
	}

	// Refine all matches and create multiple alignment
	matchRefinement(matches,str,g);
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

template<typename TFragment, typename TSpec1, typename TString, typename TSpec2, typename TValue,  typename TSpec, typename TSize>
inline void 
__getAlignmentStatistics(String<TFragment, TSpec1>& matches,
						 StringSet<TString, TSpec2>& pairSet,
						 String<TValue, TSpec>& dist,
						 TSize i,
						 TSize j,
						 TSize nseq,
						 TSize from)
{
	SEQAN_CHECKPOINT
	typedef typename Position<String<TFragment, TSpec1> >::Type TPos;
	
	// Determine a sequence weight
	TValue matchLen = 0;
	TValue overlapLen = 0;
	TValue alignLen = 0;
	getAlignmentStatistics(matches, pairSet, (TPos) from, (TPos) length(matches),  matchLen, overlapLen, alignLen);
			
	// Calculate sequence similarity
	TValue normalizedSimilarity = (matchLen / overlapLen) * (overlapLen / alignLen);

	// Assign the values
	assignValue(dist, i*nseq+j, 1 - normalizedSimilarity);
	assignValue(dist, j*nseq+i, 1 - normalizedSimilarity);	
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFragment, typename TSpec1, typename TString, typename TSpec2, typename TStringSet, typename TCargo,  typename TSpec, typename TSize>
inline void 
__getAlignmentStatistics(String<TFragment, TSpec1>& matches,
						 StringSet<TString, TSpec2>& pairSet,
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

template<typename TFragment, typename TSpec, typename TString, typename TSpec2, typename TSize>
inline void 
__getAlignmentStatistics(String<TFragment, TSpec>&,
						 StringSet<TString, TSpec2>&,
						 Nothing&,
						 TSize,
						 TSize,
						 TSize,
						 TSize)
{
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TPairList, typename TDistance, typename TAlignConfig, typename TScore>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TPairList& pList,
					   TDistance& dist,
					   TScore const& score_type,
					   TAlignConfig const& ac,
					   GlobalPairwise_Library)
{
	SEQAN_CHECKPOINT
	typedef unsigned int TSize;
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Id<TGraph>::Type TId;

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

	// Pairwise alignments
	TSize amountOfPairs = length(pList);
	for(TSize k=0; k<amountOfPairs; ++k) {
		// Make a pairwise string-set
		TStringSet pairSet;
		TId id1 = (pList[k]).i1;
		TId id2 = (pList[k]).i2;
		assignValueById(pairSet, str, id1);
		assignValueById(pairSet, str, id2);
		TSize from = length(matches);
		
		// Alignment
		globalAlignment(matches, pairSet, score_type, ac, Gotoh() );
			
		// Get the alignment statistics
		__getAlignmentStatistics(matches, pairSet, dist, idToPosition(str,id1), idToPosition(str,id2), nseq, from);
	}

	// Refine all matches and create multiple alignment
	matchRefinement(matches,str,g);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TPairList, typename TDistance, typename TScoreValue, typename TScoreSpec>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TPairList& pList,
					   TDistance& dist,
					   Score<TScoreValue, TScoreSpec> const& score_type,
					   GlobalPairwise_Library)
{
	SEQAN_CHECKPOINT
	generatePrimaryLibrary(g, pList, dist, score_type, AlignConfig<>(), GlobalPairwise_Library()); 
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TPairList, typename TScoreValue, typename TScoreSpec>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TPairList& pList,
					   Score<TScoreValue, TScoreSpec> const& score_type,
					   GlobalPairwise_Library)
{
	SEQAN_CHECKPOINT
	Nothing noth;
	generatePrimaryLibrary(g, pList, noth, score_type, GlobalPairwise_Library());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TPairList, typename TScoreValue, typename TScoreSpec, typename TAlignConfig>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TPairList& pList,
					   Score<TScoreValue, TScoreSpec> const& score_type,
					   TAlignConfig const& ac,
					   GlobalPairwise_Library)
{
	SEQAN_CHECKPOINT
	Nothing noth;
	generatePrimaryLibrary(g, pList, noth, score_type, ac, GlobalPairwise_Library());
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

template<typename TStringSet, typename TCargo, typename TSpec, typename TPairList>
inline void 
selectPairsForLibraryGeneration(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
								TPairList& pList)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TPairList>::Type TPair;

	TStringSet& str = stringSet(g);
	TSize nseq = length(str);
	for(TSize i=0; i<nseq; ++i) {
		for(TSize j=i+1; j<nseq; ++j) {
			appendValue(pList, TPair(positionToId(str, i), positionToId(str, j)));
		}
	}
}

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

//template<typename TFile, typename TStringSet, typename TCargo, typename TSpec>
//inline void
//_readLibrary(TFile & file,
//			 Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
//{
//	SEQAN_CHECKPOINT
//	typedef unsigned int TWord;
//	typedef typename Value<TFile>::Type TValue;
//	
//	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
//	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
//	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
//	typedef typename Id<TGraph>::Type TId;
//
//	TValue c;
//	bool seq1ToN = false;
//	if (_streamEOF(file)) return;
//	else c = _streamGet(file);
//
//	typedef std::pair<unsigned int, unsigned int> TSeqRes;
//	typedef std::map<TSeqRes, TVertexDescriptor> TNodeMap;
//	TNodeMap node_map;
//	TWord seq1 = 0;
//	TWord seq2 = 0;
//	while (!_streamEOF(file)) {
//		_parse_skipWhitespace(file,c);
//		if (_streamEOF(file)) break;
//		if (c == '#') {
//			c = _streamGet(file);
//			_parse_skipWhitespace(file,c);
//			seq1 = _parse_readNumber(file, c);
//			seq2 = _parse_readNumber(file, c);
//			if (empty(g))
//				if ((seq1 != 0) && (seq2 != 0)) seq1ToN = true;
//			if (seq1ToN) {
//				--seq1;
//				--seq2;
//			}
//		} else if (c == '!') {
//			_parse_skipLine(file, c);
//		} else {
//			unsigned int res1 = _parse_readNumber(file, c);
//			_parse_skipWhitespace(file,c);
//			unsigned int res2 = _parse_readNumber(file, c);
//			_parse_skipWhitespace(file,c);
//			unsigned int weight = _parse_readNumber(file, c);
//			_parse_skipLine(file,c);
//		
//			// Insert new vertex if necessary
//			--res1;
//			--res2;
//			bool newEdge = false;
//			TSeqRes key = std::make_pair(seq1, res1);
//			typename TNodeMap::iterator nodePos = node_map.find(key);
//			TId id1;
//			if (nodePos == node_map.end()) {
//				//std::cout << seq1 << ',' << res1 << std::endl;
//				id1 = addVertex(g, seq1, res1, 1); 
//				node_map.insert(std::make_pair(key, id1));
//				newEdge = true;
//			} else {
//				id1 = nodePos->second;
//			}
//
//			key = std::make_pair(seq2, res2);
//			nodePos = node_map.find(key);
//			TId id2;
//			if (nodePos == node_map.end()) {
//				//std::cout << seq2 << ',' << res2 << std::endl;
//				id2 = addVertex(g, seq2, res2, 1); 
//				node_map.insert(std::make_pair(key, id2));
//				newEdge = true;
//			} else {
//				id2 = nodePos->second;
//			}
//
//			// Insert a new edge or adapt the weight
//			if (newEdge) addEdge(g,id1,id2,weight);
//			else {
//				TEdgeDescriptor e = findEdge(g,id1,id2);
//				if( e == 0 ) addEdge(g,id1,id2,weight);
//				else cargo(e) += weight;
//			}
//		}
//	}
//}

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

	typedef std::pair<std::pair<TWord, TWord>, TCargo> TResiduePair;
	typedef std::set<TResiduePair> TResiduePairSet;
	TWord nseq = length(stringSet(g));
	String<TResiduePairSet> resPair;
	resize(resPair, nseq * nseq);	
	TWord seq1 = 0;
	TWord seq2 = 0;
	bool firstPass = true;
	while (!_streamEOF(file)) {
		_parse_skipWhitespace(file,c);
		if (_streamEOF(file)) break;
		if (c == '#') {
			c = _streamGet(file);
			_parse_skipWhitespace(file,c);
			seq1 = _parse_readNumber(file, c);
			seq2 = _parse_readNumber(file, c);
			if (firstPass) {
				firstPass = false;
				if ((seq1 != 0) && (seq2 != 0)) seq1ToN = true;
			}
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

			TWord index = 0;
			if (seq1 < seq2) index = seq1 * nseq + seq2;
			else index = seq1 * nseq + seq2;
			resPair[index].insert(std::make_pair(std::make_pair(res1,res2), weight));
		}
	}

	typedef Fragment<TWord, TWord, TWord> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;
	String<TCargo, Block<> > score_values;

	for(unsigned int i = 0; i<length(resPair); ++i) {
		if (resPair[i].empty()) continue;
		TWord seq1 = i / nseq;
		TWord seq2 = i % nseq;
		//std::cout << "#" << seq1 << ',' << seq2 << std::endl;
		typename TResiduePairSet::const_iterator pos = resPair[i].begin();
		typename TResiduePairSet::const_iterator posEnd = resPair[i].end();
		TWord startMatch1 = pos->first.first;
		TWord startMatch2 = pos->first.second;
		TCargo carg = pos->second;
		TWord len = 1;
		++pos;
		while(pos != posEnd) {
			if ((startMatch1 + len == pos->first.first) &&
				(startMatch2 + len == pos->first.second)) {
					carg += pos->second;
					++len;
			} else {
				push_back(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
				push_back(score_values, carg);
				startMatch1 = pos->first.first;
				startMatch2 = pos->first.second;
				carg = pos->second;
				len = 1;
			}
			//std::cout << pos->first.first << ',' << pos->first.second << ',' << pos->second << std::endl;
			++pos;
		}
		push_back(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
		push_back(score_values, carg);
	}

	// Refine all matches, rescore matches and create multiple alignment
	matchRefinement(matches,stringSet(g),g);

	// Adapt edge weights, scale weights by significance of local match
	TFragmentStringIter endIt = end(matches);
	unsigned int positionIt = 0;
	for(TFragmentStringIter it = begin(matches); it != endIt; ++it, ++positionIt) {
		TId id1 = sequenceId(*it,0);
		TId id2 = sequenceId(*it,1);
		TWord pos1 = fragmentBegin(*it, id1);
		TWord pos2 = fragmentBegin(*it, id2);
		TWord origFragLen = fragmentLength(*it, id1);
		TWord end1 = pos1 + origFragLen;
		while(pos1 < end1) {
			TVertexDescriptor p1 = findVertex(g, id1, pos1);
			TVertexDescriptor p2 = findVertex(g, id2, pos2);
			TWord fragLen = fragmentLength(g, p1);
			TEdgeDescriptor e = findEdge(g, p1, p2);
			cargo(e) *= (TCargo) ( (double) fragLen / (double) origFragLen * (double) getValue(score_values, positionIt));
			pos1 += fragLen;
			pos2 += fragLen;
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

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TNames>
void 
read(TFile & file,
	 Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
	 TNames& names, 
	 MemeMotif) 
{
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;
	typedef typename Value<TNames>::Type TName;
	typedef std::map<TName, TId> TNameToId;
	TNameToId nameIdMap;

	// Initialization
	TWord nseq = length(names);
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;
	for(TWord i = 0;i<length(names);++i) {
		nameIdMap.insert(std::make_pair(names[i], i));
	}

	// Read the motives
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);
	_parse_skipLine(file, c);
	
	TName seq;
	String<char> sequence;
	String<TId> ids;
	String<TWord> beginPos;
	reserve(ids, nseq);
	reserve(beginPos, nseq);
	TWord len = 0;
	while (!_streamEOF(file)) {
		if (c == '/') {
			_parse_skipLine(file, c);
			for(TWord i = 0; i<length(ids) - 1;++i) {
				for(TWord j = i+1; j<length(ids);++j) {
					push_back(matches, TFragment(ids[i], beginPos[i], ids[j],  beginPos[j], len));
				}
			}
			len = 0;
			clear(ids);
			clear(beginPos);
			reserve(ids, nseq);
			reserve(beginPos, nseq);
			if (_streamEOF(file)) break;
			_parse_skipLine(file, c);
		}
		clear(seq);
		_parse_skipWhitespace(file, c);
		seq = _parse_readIdentifier(file, c);
		TId seqId = nameIdMap[seq];
		_parse_skipWhitespace(file, c);
		c = _streamGet(file);
		_parse_skipWhitespace(file, c);
		TWord beg = _parse_readNumber(file, c);
		--beg;
		if (len == 0) {
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);
			clear(sequence);
			sequence = _parse_readIdentifier(file, c);
			len = length(sequence);
		}
		_parse_skipLine(file, c);
		appendValue(ids, seqId);
		appendValue(beginPos, beg);
		//std::cout << seqId << ',' << beg << ',' << len << std::endl;
	}

	// Refine all matches, rescore matches and create multiple alignment
	matchRefinement(matches,stringSet(g),g);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TNames>
void 
read(TFile & file,
	 Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
	 TNames& names, 
	 BlastLib) 
{
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;
	typedef typename Value<TNames>::Type TName;
	typedef std::map<TName, TWord> TNameToPosition;
	TNameToPosition namePosMap;

	// Initialization
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;
	for(TWord i = 0;i<length(names);++i) {
		namePosMap.insert(std::make_pair(names[i], i));
	}
	
	// Read the Blast file
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);
	TStringSet& strSet = stringSet(g);

	TName seq1;
	TName seq2;
	TWord window = 1000;
	while (!_streamEOF(file)) {
		clear(seq1);
		clear(seq2);
		_parse_skipWhitespace(file, c);
		seq1 = _parse_readIdentifier(file, c);
		_parse_skipWhitespace(file, c);
		seq2 = _parse_readIdentifier(file, c);
		if (seq1 == seq2) {
			_parse_skipLine(file, c);
			continue;
		}
		TWord seq1Id = namePosMap[seq1];
		TWord seq2Id = namePosMap[seq2];
		_parse_skipWhitespace(file, c);
		_parse_readDouble(file, c);
		_parse_skipWhitespace(file, c);
		_parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		_parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		_parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		TWord beg1 = _parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		TWord end1 = _parse_readNumber(file, c);
		TWord len = end1 - beg1 + 1;
		_parse_skipWhitespace(file, c);
		TWord beg2 = _parse_readNumber(file, c);
		--beg1;
		--beg2;

		if (((beg1 > beg2) && ((beg1 - beg2) > window)) ||
			((beg1 < beg2) && ((beg2 - beg1) > window))) {
				std::cout << "Out of window" << std::endl;
				_parse_skipLine(file, c);
				continue;
		}
		if	(((beg1 + len) > length(strSet[seq1Id])) ||
			((beg2 + len) > length(strSet[seq2Id]))) {
				std::cout << "Wrong match" << std::endl;
				_parse_skipLine(file, c);
				continue;
		}

		// Debug code
		std::cout << seq1Id << ',' << beg1 << ',' << seq2Id << ',' << beg2 << ',' << len << std::endl;
		std::cout << infix(strSet[seq1Id], beg1, beg1+len) << std::endl;
		std::cout << infix(strSet[seq2Id], beg2, beg2+len) << std::endl;
		
		push_back(matches, TFragment(seq1Id, beg1, seq2Id,  beg2, len));
		_parse_skipLine(file, c);
	}

	//_debugMatches(stringSet(g), matches);

	// Refine all matches, rescore matches and create multiple alignment
	matchRefinement(matches,strSet,g);
	
	//std::cout << g << std::endl;
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
	
	typedef std::pair<std::pair<TSize, TSize>, TCargo> TResiduePair;
	typedef std::set<TResiduePair> TResiduePairSet;
	TSize nseq = length(stringSet(g));
	String<TResiduePairSet> resPair;
	resize(resPair, nseq * nseq);
	
	typedef typename Iterator<TGraph, EdgeIterator>::Type TIter;
	TIter it(g);
	for(;!atEnd(it);++it) {
		TVertexDescriptor sV = sourceVertex(it);
		TVertexDescriptor tV = targetVertex(it);
		TSize fragLen = fragmentLength(g,sV);
		TSize fragPos1 = fragmentBegin(g,sV);
		TSize fragPos2 = fragmentBegin(g,tV);
		TSize seq1 = sequenceId(g,sV);
		TSize seq2 = sequenceId(g,tV);
		if (seq1 > seq2) {
			TSize tmp = seq1;
			TSize tmp_pos = fragPos1;
			seq1 = seq2;
			fragPos1 = fragPos2;
			seq2 = tmp;
			fragPos2 = tmp_pos;
		}
		for(TSize i = 0; i<fragLen; ++i) {
			TCargo my_carg = (TCargo) ((double) 1.0 / (double) fragLen * (double) getCargo(*it));
			if (my_carg <= 0) my_carg = 1;
			resPair[seq1 * nseq + seq2].insert(std::make_pair(std::make_pair(fragPos1 + i, fragPos2 + i), my_carg));
		}
	}

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

	for(unsigned int i = 0; i<length(resPair); ++i) {
		if (resPair[i].empty()) continue;
		TSize seq1 = i / nseq;
		TSize seq2 = i % nseq;
		_streamPut(file, '#');
		_streamPutInt(file, seq1 + 1);
		_streamPut(file, ' ');
		_streamPutInt(file, seq2 + 1);
		_streamPut(file, '\n');	
		typename TResiduePairSet::const_iterator pos = resPair[i].begin();
		typename TResiduePairSet::const_iterator posEnd = resPair[i].end();
		while(pos != posEnd) {
			_streamPutInt(file, pos->first.first + 1);
			_streamPut(file, ' ');
			_streamPutInt(file, pos->first.second + 1);
			_streamPut(file, ' ');
			_streamPutInt(file, pos->second);
			_streamPut(file, '\n');	
			++pos;
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
