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
// Alignment graph generation
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Tags
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Library Generation:
..summary:A tag that specifies how to generate alignment graph edges (T-Coffee library).
*/


/**
.Tag.Library Generation.value.GlobalPairwise_Library:
	A primary library of global alignments.
*/

struct GlobalPairwise_Library_;
typedef Tag<GlobalPairwise_Library_> const GlobalPairwise_Library;


/**
.Tag.Library Generation.value.GlobalPairwise_Library:
	A primary library of local alignments.
*/

struct LocalPairwise_Library_;
typedef Tag<LocalPairwise_Library_> const LocalPairwise_Library;


/**
.Tag.Library Generation.value.Overlap_Library:
	A primary library of overlap alignments.
*/

struct Overlap_Library_;
typedef Tag<Overlap_Library_> const Overlap_Library;

/**
.Tag.Library Generation.value.Kmer_Library:
	A primary library of kmer alignments.
*/

struct Kmer_Library_;
typedef Tag<Kmer_Library_> const Kmer_Library;


/**
.Tag.Library Generation.value.Lcs_Library:
	A primary library generated using the longest common subsequence algorithm.
*/

struct Lcs_Library_;
typedef Tag<Lcs_Library_> const Lcs_Library;



//////////////////////////////////////////////////////////////////////////////
// Pair selection to calculate alignment
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

// Dummy function selecting all pairs
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
// Alignment statistics
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TFragment, typename TSpec1, typename TStringSet, typename TPos, typename TSize1, typename TAlphabet>
inline void 
getAlignmentStatistics(String<TFragment, TSpec1> const& matches,
					   TStringSet& str,
					   TPos const from,
					   TPos const to,
					   TSize1& matchLength,	// Number of identical characters
					   TSize1& overlapLength,	// Number of character in overlapping segments (with mismatches and gaps)
					   TSize1& alignLength,	// Length of the alignment
					   TAlphabet)
{
	SEQAN_CHECKPOINT
	typedef String<TFragment, TSpec1> TFragmentMatches;
	typedef typename Size<TFragmentMatches>::Type TSize;
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

template<typename TFragment, typename TSpec1, typename TStringSet, typename TPos, typename TSize>
inline void 
getAlignmentStatistics(String<TFragment, TSpec1>& matches,
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

template<typename TFragment, typename TSpec1, typename TStringSet, typename TSize>
inline void 
getAlignmentStatistics(String<TFragment, TSpec1>& matches,
					   TStringSet& str,
					   TSize& matchLength,
					   TSize& overlapLength,
					   TSize& alignLength)
{
	SEQAN_CHECKPOINT
	getAlignmentStatistics(matches, str, (TSize) 0, (TSize) length(matches), matchLength, overlapLength, alignLength, typename Value<TStringSet>::Type() );
}



//////////////////////////////////////////////////////////////////////////////
// Library Generation
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TScore, typename TPairList>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TPairList& pList,
					   TScore const& score_type,
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

	// String of fragments to combine all pairwise alignments into a multiple alignment
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;

	TSize amountOfPairs = length(pList);
	for(TSize k=0; k<amountOfPairs; ++k) {
		TSize i = idToPosition(str, (pList[k]).i1);
		TSize j = idToPosition(str, (pList[k]).i2);
		
		// Lcs between first and second string
		String<std::pair<unsigned int, unsigned int>, Block<> > pos1;
		longestCommonSubsequence(str[i], str[j], 1000, pos1);

		// Extend the matches as long as possible
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

	// Refine all matches and create multiple alignment
	matchRefinement(matches,str,const_cast<TScore&>(score_type),g);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScore, typename TSize>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TScore const& score_type,
					   TSize ktup,
					   Kmer_Library)
{
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Value<typename Value<TStringSet>::Type>::Type TAlphabet;
	typedef String<TWord> TTupelString;
	typedef String<TTupelString> TTupelStringSet;
	
	// Clear Graph
	clearVertices(g);

	// Initialization
	TStringSet& str = stringSet(g);	
	TSize nseq = length(str);
	TWord alphabet_size = ValueSize<TAlphabet>::VALUE;

	// All matching kmers
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	TFragmentString matches;

	// Transform the set of strings into a set of strings of k-tupels
	TTupelStringSet tupSet;
	resize(tupSet, nseq);
	for(TSize k=0;k<nseq;++k) {
		_getTupelString(str[k], tupSet[k], ktup, TAlphabet());
	}

	// Build for each sequence the q-gram Index and count common hits
	String<TWord, Block<> > qIndex;
	TWord qIndexSize = 1;
	for (TWord i=0; i< (TWord) ktup;++i) qIndexSize *= alphabet_size;
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
			for(TSize i = 0;i < (TSize) length(tupSet[k2]);++i) {
				TWord current_kmer = tupSet[k2][i];
				TWord amount = 0;
				if (current_kmer == qIndexSize - 1) amount = sum - qIndex[current_kmer];
				else amount = qIndex[current_kmer+1] - qIndex[current_kmer];
				for(TWord pos = qIndex[current_kmer]; pos < qIndex[current_kmer] + amount; ++pos) {
					//std::cout << id1 << ',' << positions[pos] << ':' << id2 << ',' << i << std::endl;
					//std::cout << tupSet[k2][i] << ',' << tupSet[k][positions[pos]] << std::endl;
					push_back(matches, TFragment(id1, positions[pos], id2, i, ktup));
				}
			}
		}
	}
	
	matchRefinement(matches,str,const_cast<TScore&>(score_type),g);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScore>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TScore const& score_type,
					   Kmer_Library)
{
	SEQAN_CHECKPOINT
	generatePrimaryLibrary(g, score_type, 3, Kmer_Library());
}


//////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TScore>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TScore const& score_type,
					   Lcs_Library)
{
	SEQAN_CHECKPOINT
	typedef typename Id<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type TId;
	String<Pair<TId, TId> > pList;
	selectPairsForLibraryGeneration(g, pList);
	generatePrimaryLibrary(g, pList, score_type, Lcs_Library() );
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
	matchRefinement(matches,str,const_cast<TScore&>(score_type),g);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScore>
void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TScore const& score_type,
					   LocalPairwise_Library)
{
	SEQAN_CHECKPOINT
	typedef typename Id<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type TId;
	String<Pair<TId, TId> > pList;
	selectPairsForLibraryGeneration(g, pList);
	generatePrimaryLibrary(g, pList, score_type, LocalPairwise_Library() );
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

template<typename TFragment, typename TSpec1, typename TString, typename TSpec2, typename TCargo,  typename TSpec, typename TSize>
inline void 
__getAlignmentStatistics(String<TFragment, TSpec1>& matches,
						 StringSet<TString, TSpec2>& pairSet,
						 Graph<Undirected<TCargo, TSpec> >& dist,
						 TSize i,
						 TSize j,
						 TSize,
						 TSize from)
{
	SEQAN_CHECKPOINT
		
	// Determine a sequence weight
	TCargo matchLen = 0;
	TCargo overlapLen = 0;
	TCargo alignLen = 0;
	getAlignmentStatistics(matches, pairSet, (TSize) from, (TSize) length(matches),  matchLen, overlapLen, alignLen);
			
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
		__getAlignmentStatistics(matches, pairSet, dist, (TSize) idToPosition(str,id1), (TSize) idToPosition(str,id2), (TSize) nseq, (TSize) from);
	}

	// Refine all matches and create multiple alignment
	matchRefinement(matches,str,const_cast<TScore&>(score_type),g);
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

template<typename TStringSet, typename TCargo, typename TSpec, typename TDistance, typename TScoreValue, typename TScoreSpec>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TDistance& dist,
					   Score<TScoreValue, TScoreSpec> const& score_type,
					   GlobalPairwise_Library)
{
	SEQAN_CHECKPOINT
	typedef typename Id<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type TId;
	String<Pair<TId, TId> > pList;
	selectPairsForLibraryGeneration(g, pList);
	generatePrimaryLibrary(g, pList, dist, score_type, GlobalPairwise_Library()); 
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TId, typename TScoreValue, typename TScoreSpec>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   String<Pair<TId, TId> >& pList,
					   Score<TScoreValue, TScoreSpec> const& score_type,
					   GlobalPairwise_Library)
{
	SEQAN_CHECKPOINT
	Nothing noth;
	generatePrimaryLibrary(g, pList, noth, score_type, GlobalPairwise_Library());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue, typename TScoreSpec>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   Score<TScoreValue, TScoreSpec> const& score_type,
					   GlobalPairwise_Library)
{
	SEQAN_CHECKPOINT
	typedef typename Id<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type TId;
	String<Pair<TId, TId> > pList;
	selectPairsForLibraryGeneration(g, pList);
	generatePrimaryLibrary(g, pList, score_type, GlobalPairwise_Library());
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

template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue, typename TScoreSpec, typename TAlignConfig>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   Score<TScoreValue, TScoreSpec> const& score_type,
					   TAlignConfig const& ac,
					   GlobalPairwise_Library)
{
	SEQAN_CHECKPOINT
	typedef typename Id<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type TId;
	String<Pair<TId, TId> > pList;
	selectPairsForLibraryGeneration(g, pList);
	generatePrimaryLibrary(g, pList, score_type, ac, GlobalPairwise_Library());
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
				//std::cout << "Out of window" << std::endl;
				_parse_skipLine(file, c);
				continue;
		}
		if	(((beg1 + len) > length(strSet[seq1Id])) ||
			((beg2 + len) > length(strSet[seq2Id]))) {
				//std::cout << "Wrong match" << std::endl;
				_parse_skipLine(file, c);
				continue;
		}

		//// Debug code
		//std::cout << seq1Id << ',' << beg1 << ',' << seq2Id << ',' << beg2 << ',' << len << std::endl;
		//std::cout << infix(strSet[seq1Id], beg1, beg1+len) << std::endl;
		//std::cout << infix(strSet[seq2Id], beg2, beg2+len) << std::endl;
		
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
