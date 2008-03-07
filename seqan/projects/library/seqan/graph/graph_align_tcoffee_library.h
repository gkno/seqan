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
.Tag.Library Generation.value.LocalPairwise_Library:
	A primary library of local alignments.
*/

struct LocalPairwise_Library_;
typedef Tag<LocalPairwise_Library_> const LocalPairwise_Library;

/**
.Tag.Library Generation.value.LocalExtendPairwise_Library:
	A primary library of extended local alignments.
*/

struct LocalExtendPairwise_Library_;
typedef Tag<LocalExtendPairwise_Library_> const LocalExtendPairwise_Library;


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
	resize(pList, nseq * (nseq - 1) / 2);
	TSize count = 0;
	for(TSize i=0; i<nseq-1; ++i) {
		for(TSize j=i+1; j<nseq; ++j) {
			value(pList, count) = TPair(positionToId(str, i), positionToId(str, j));
			++count;
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
	typedef String<TFragment > TFragmentString;
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
						
				appendValue(matches, TFragment(positionToId(str, i),iBegin,positionToId(str, j),jBegin,lenMatch));
				lenMatch = 1;
				iBegin = pos1[z].first;
				jBegin = pos1[z].second;
			}
		}
		// Process last match
		appendValue(matches, TFragment(positionToId(str, i),iBegin,positionToId(str, j),jBegin,lenMatch));
	}

	// Refine all matches and create multiple alignment
	matchRefinement(matches,str,const_cast<TScore&>(score_type),g);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScore, typename TAlphabet, typename TSize>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TScore const& score_type,
					   TSize ktup,
					   TAlphabet,
					   Kmer_Library)
{
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Id<TGraph>::Type TId;
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
	typedef String<TFragment> TFragmentString;
	TFragmentString matches;

	// Transform the set of strings into a set of strings of k-tupels
	TTupelStringSet tupSet;
	resize(tupSet, nseq);
	for(TSize k=0;k<nseq;++k) {
		_getTupelString(str[k], tupSet[k], ktup, TAlphabet());
	}

	// Build one q-gram Index for all sequences
	typedef std::pair<TWord, TWord> TPosSeqPair;
	typedef std::set<TPosSeqPair> TQGramOcc;
	String<TQGramOcc> qIndex;
	TWord qIndexSize = 1;
	for(TWord i=0; i<(TWord) ktup;++i) qIndexSize *= alphabet_size;
	resize(qIndex, qIndexSize);
	for(TSize k=0;k<nseq;++k) {
		for(TWord i = 0;i < (TWord) length(tupSet[k]);++i) {
			qIndex[ tupSet[k][i] ].insert(std::make_pair(i, k));
		}
	}
	for(TSize q=0;q< (TSize) qIndexSize;++q) {
		typename TQGramOcc::const_iterator pos = qIndex[q].begin();
		typename TQGramOcc::const_iterator posEnd = qIndex[q].end();
		while (pos != posEnd) {
			typename TQGramOcc::const_iterator pos2 = pos;
			++pos2;
			while (pos2 != posEnd) {
				if (pos->second != pos2->second) {
					appendValue(matches, TFragment(positionToId(str, pos->second), pos->first, positionToId(str, pos2->second), pos2->first, ktup));
				}
				++pos2;
			}
			++pos;
		}
	}
	
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
	generatePrimaryLibrary(g, score_type, ktup, typename Value<typename Value<TStringSet>::Type>::Type(), Kmer_Library());
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
	typedef String<TFragment > TFragmentString;
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
			appendValue(matches, TFragment( (unsigned int) sequenceId(pGraph, sV), (unsigned int) fragmentBegin(pGraph,sV), (unsigned int) sequenceId(pGraph, tV),  (unsigned int)  fragmentBegin(pGraph,tV),  (unsigned int)  fragmentLength(pGraph,tV)));
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
	reserve(_getVertexString(dist), nseq);
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
	typedef String<TFragment > TFragmentString;
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

template<typename TStringSet, typename TCargo, typename TSpec, typename TPairList, typename TScore>
void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TPairList& pList,
					   TScore const& score_type,
					   LocalExtendPairwise_Library)
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

	// String of fragments
	typedef Fragment<> TFragment;
	typedef String<TFragment> TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;

	// All local pairwise alignments
	TSize amountOfPairs = length(pList);
	typedef Triple<TSize, TSize, TSize> TTripel;
	typedef std::map<int, TTripel, std::greater<double> > TLocalAlignMap;
	TLocalAlignMap localAlignMap;
	for(TSize k=0; k<amountOfPairs; ++k) {
		// Make a pairwise string-set
		TStringSet pairSet;
		assignValueById(pairSet, str, (pList[k]).i1);
		assignValueById(pairSet, str, (pList[k]).i2);

		// Perform an alignment
		TSize from = length(matches);
		int score = localAlignment(matches, pairSet, score_type, SmithWaterman() );
		TSize to = length(matches);

		localAlignMap.insert(std::make_pair(score, TTripel(k, from, to)));
	}

	// Extend the local alignments to all other sequences
	TSize countMax = 50;
	typename TLocalAlignMap::const_iterator itMap = localAlignMap.begin();
	typename TLocalAlignMap::const_iterator itMapEnd = localAlignMap.end();
	for(TSize count = 0; ((count < countMax) && (itMap != itMapEnd)); ++count, ++itMap) {
		TSize k = itMap->second.i1;
		TSize from = itMap->second.i2;
		TSize to = itMap->second.i3;
		TId id1 = (pList[k]).i1;
		TId id2 = (pList[k]).i2;

		// Get the subsequences taking part in the local alignment
		TSize beg1 = length(getValueById(str, id1));
		TSize beg2 = length(getValueById(str, id2));
		TSize end1 = 0;
		TSize end2 = 0;
		for(TSize iter = from; iter < to; ++iter) {
			if (fragmentBegin(matches[iter], id1) < beg1) beg1 = fragmentBegin(matches[iter], id1);
			if (fragmentBegin(matches[iter], id2) < beg2) beg2 = fragmentBegin(matches[iter], id2); 
			if (fragmentBegin(matches[iter], id1) + fragmentLength(matches[iter], id1) > end1) end1 = fragmentBegin(matches[iter], id1) + fragmentLength(matches[iter], id1);
			if (fragmentBegin(matches[iter], id2) + fragmentLength(matches[iter], id2) > end2) end2 = fragmentBegin(matches[iter], id2) + fragmentLength(matches[iter], id2);
		}
		typename Value<TStringSet>::Type str0 = infix(getValueById(str, id1), beg1, end1);
		typename Value<TStringSet>::Type str1 = infix(getValueById(str, id2), beg2, end2);


		// Align the subsequences with all other sequences
		for(TSize seq = 0; seq < nseq; ++seq) {
			TId id3 = positionToId(str, seq);
			if ((id3 == id1) || (id3 == id2)) continue;

			TStringSet pairSet;
			assignValueById(pairSet, str0, id1);
			assignValueById(pairSet, str, id3);

			from = length(matches);
			localAlignment(matches, pairSet, score_type, SmithWaterman() );
			//globalAlignment(matches, pairSet, score_type, AlignConfig<true, false, false, true>(), Gotoh() );
			to = length(matches);
			
			// Move the matches given our offset
			for(TSize iter = from; iter < to; ++iter) fragmentBegin(matches[iter], id1) += beg1;

			clear(pairSet);
			assignValueById(pairSet, str1, id2);
			assignValueById(pairSet, str, id3);

			from = length(matches);
			//globalAlignment(matches, pairSet, score_type, AlignConfig<true, false, false, true>(), Gotoh() );
			localAlignment(matches, pairSet, score_type, SmithWaterman() );
			to = length(matches);
			
			// Move the matches given our offset
			for(TSize iter = from; iter < to; ++iter) fragmentBegin(matches[iter], id2) += beg2;
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
					   LocalExtendPairwise_Library)
{
	SEQAN_CHECKPOINT
	typedef typename Id<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type TId;
	String<Pair<TId, TId> > pList;
	selectPairsForLibraryGeneration(g, pList);
	generatePrimaryLibrary(g, pList, score_type, LocalExtendPairwise_Library() );
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
