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
  $Id: graph_consensus_base.h 2103 2008-05-23 07:57:13Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_CONSENSUS_BASE_H
#define SEQAN_HEADER_GRAPH_CONSENSUS_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File tags
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.Amos message file:
	Amos message file.
*/
struct TagAmos_;
typedef Tag<TagAmos_> const Amos;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.CeleraFrg message file:
	Celera fragment message file.
*/
struct TagCeleraFrg_;
typedef Tag<TagCeleraFrg_> const CeleraFrg;

/**
.Tag.File Format.tag.CeleraCgb message file:
	Celera cgb file.
*/
struct TagCeleraCgb_;
typedef Tag<TagCeleraCgb_> const CeleraCgb;



//////////////////////////////////////////////////////////////////////////////
// Segment Match Generation tag
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Segment Match Generation.value.Overlap_Library:
	Segment matches from overlap alignments.
*/

struct Overlap_Library_;
typedef Tag<Overlap_Library_> const Overlap_Library;



//////////////////////////////////////////////////////////////////////////////
// Consensus tag
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Consensus Calling:
..summary:A tag that specifies how to call the consensus.
*/


/**
.Tag.Consensus Calling.value.Majority_Vote:
	A consensus based on the most common character.
*/

struct Majority_Vote_;
typedef Tag<Majority_Vote_> const Majority_Vote;

//////////////////////////////////////////////////////////////////////////////
// Read alignment and Consensus Generation
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet, typename TCargo, typename TSpec, typename TValue, typename TBegEndRowPos, typename TSize>
inline void
multireadAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
				   String<TValue>& mat,
				   TBegEndRowPos& readBegEndRowPos,
				   TSize& maxCoverage)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Value<TString>::Type TAlphabet;
	typedef typename Infix<TString>::Type TInfix;
	typedef typename Iterator<TInfix>::Type TInfixIter;
	typedef std::map<std::pair<TIdType, TIdType>, TVertexDescriptor> TPosToVertexMap;
	typedef std::map<TSize, TSize> TComponentLength;
	
	// Strongly Connected Components, topological sort, and length of each component
	String<TSize> component;
	String<TSize> order;
	TComponentLength compLength;
	if (convertAlignment(g, component, order, compLength)) {
		TSize numOfComponents = length(order);
		TStringSet& strSet = stringSet(g);
		TSize nseq = length(strSet);
		clear(readBegEndRowPos);
		resize(readBegEndRowPos, nseq);
		
		// Assign to each sequence the start and end (in terms of component ranks)
		typedef std::map<TSize, TSize> TComponentToRank;
		TComponentToRank compToRank;
		for(TSize compIndex = 0; compIndex < numOfComponents; ++compIndex) {
			compToRank.insert(std::make_pair(order[compIndex], compIndex));
		}
		typedef Pair<TSize, TSize> TRankPair;
		typedef String<TRankPair> TSequenceToRanks;
		TSequenceToRanks seqToRank;
		fill(seqToRank, nseq, TRankPair(0, 0));	
		TVertexIterator itVertex(g);
		for(;!atEnd(itVertex);++itVertex) {
			TVertexDescriptor vert = value(itVertex);
			TSize seq = idToPosition(strSet, sequenceId(g, vert));
			if (fragmentBegin(g, vert) == 0) {
				(value(seqToRank, seq)).i1 = (compToRank.find(getProperty(component,vert)))->second;
			}
			if (fragmentBegin(g, vert) + fragmentLength(g, vert) == length(value(strSet, seq))) {
				(value(seqToRank, seq)).i2 = (compToRank.find(getProperty(component,vert)))->second;
			}
		}
		compToRank.clear();

		// Assign the sequences to rows
		String<TSize> seqToRow;
		resize(seqToRow, nseq);
		maxCoverage = 0;
		typedef std::set<TSize> TLeftOver;
		TLeftOver leftOver;
		for(TSize i=0;i<nseq; ++i) {
			leftOver.insert(i);
		}
		while(!leftOver.empty()) {
			typedef std::set<std::pair<TSize, TSize> > TSeqToBegin;
			TSeqToBegin seqToBegin;
			for(typename TLeftOver::const_iterator pos = leftOver.begin(); pos != leftOver.end(); ++pos) {
				seqToBegin.insert(std::make_pair((seqToRank[*pos]).i1, *pos));
			}
			TSize endPos = 0;
			for(typename TSeqToBegin::const_iterator s = seqToBegin.begin(); s != seqToBegin.end();++s) {
				if (endPos <= (*s).first) {
					TSize currentSeq = (*s).second;
					seqToRow[currentSeq] = maxCoverage;
					endPos = (seqToRank[currentSeq]).i2 + 2;
					leftOver.erase(currentSeq);
				}	
			}
			++maxCoverage;
		}


		// Create the matrix
		TSize len = 0;
		String<TSize> compOffset;
		fill(compOffset, numOfComponents, 0);
		for(TSize compIndex = 0; compIndex < numOfComponents; ++compIndex) {
			compOffset[order[compIndex]] = len;
			len+=compLength[order[compIndex]];
		}
		TValue gapChar = gapValue<TValue>();
		TValue specialGap = '.';
		fill(mat, len * maxCoverage, gapChar);

		// Fill in the segments
		for(typename TPosToVertexMap::const_iterator it = g.data_pvMap.begin();it != g.data_pvMap.end(); ++it) {
			TInfix str = label(g,it->second);
			TSize c = property(component, it->second);
			TSize strPos = idToPosition(strSet, it->first.first);
			TSize row = seqToRow[strPos];
			//if (row == 0) {
			//	std::cout << sequenceId(g, it->second) << ':' << str << ',' << strSet[sequenceId(g, it->second)] << std::endl;
			//	std::cout << getProperty(component, it->second) << ',' << order[compIndex] << std::endl;
			//	std::cout << (seqToRank[sequenceId(g, it->second)]).i1 << ',' << (seqToRank[sequenceId(g, it->second)]).i2 << std::endl;
			//}
			TInfixIter sIt = begin(str);
			TInfixIter sItEnd = end(str);
			TSize i = compOffset[c];
			for(TSize pCol = i;sIt!=sItEnd;goNext(sIt), ++pCol, ++i) assignValue(mat, row * len + pCol, *sIt);
		}
		String<bool> active;
		for(TSize compIndex = 0; compIndex < numOfComponents; ++compIndex) {
			TSize offset = compOffset[order[compIndex]];
			TSize currentCompLength = compLength[order[compIndex]];

			clear(active);
			fill(active, maxCoverage, false);

			// Find the empty rows
			for(TSize i=0;i<nseq; ++i) {
				if (((seqToRank[i]).i1 <= compIndex) && ((seqToRank[i]).i2 >= compIndex)) active[(seqToRow[i])] = true;
			}
			
			// Substitute false gaps with special gap character
			for(TSize i = 0; i < maxCoverage; ++i) {
				if (!(active[i])) {
					for(TSize pCol = offset;pCol < offset + currentCompLength;++pCol) assignValue(mat, i * len + pCol, specialGap);
				}
			}

		}

		// Get the new begin and end positions
		for(TSize i=0;i<nseq; ++i) {
			TVertexDescriptor lastVertex = findVertex(const_cast<TGraph&>(g), positionToId(strSet, i), length(strSet[i]) - 1);
			TSize readBegin = compOffset[getProperty(component, findVertex(const_cast<TGraph&>(g), positionToId(strSet, i), 0))];
			TSize readEnd = compOffset[getProperty(component, lastVertex)] + fragmentLength(const_cast<TGraph&>(g), lastVertex);
			readBegEndRowPos[i].i1 = readBegin;
			readBegEndRowPos[i].i2 = readEnd;
			readBegEndRowPos[i].i3 = seqToRow[i];
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TAlphabet, typename TSpec2, typename TGappedConsensus, typename TCoverage, typename TSize>
inline void
consensusCalling(String<TValue, TSpec> const& mat,
				 String<TAlphabet, TSpec2>& consensus,
				 TGappedConsensus& gappedConsensus,
				 TCoverage& coverage,
				 TSize alignDepth,
				 Majority_Vote)
{
	SEQAN_CHECKPOINT
	typedef String<TValue, TSpec> TMatrix;
	typedef typename Iterator<TMatrix>::Type TMatIter;
	typedef typename Iterator<TCoverage>::Type TCovIter;

	// Initialization
	TSize len = length(mat) / alignDepth;
	TValue gapChar = gapValue<TValue>();
	TValue specialGap = '.';
	TSize alphabetSize = ValueSize<TAlphabet>::VALUE;
	clear(coverage);
	fill(coverage, len, 0);


	// Set-up counter values
	typedef String<TSize> TCounter;
	typedef typename Iterator<TCounter>::Type TCIt;
	typedef String<TCounter> TCounters;
	typedef typename Iterator<TCounters>::Type TCounterIt;
	String<String<TSize> > counterValues;
	resize(counterValues, len);
	for(TSize i=0;i<len; ++i) {
		String<TSize> counter;
		fill(counter, alphabetSize + 1, 0);
		value(counterValues, i) = counter;
	}

	// Count all 
	TMatIter matIt = begin(mat);
	TMatIter matItEnd = end(mat);
	TCounterIt countIt = begin(counterValues);
	TCovIter covIt = begin(coverage);
	TSize pos = 0;
	for(; matIt != matItEnd; goNext(matIt), goNext(countIt), goNext(covIt), ++pos) {
		if (pos % len == 0) {
			countIt = begin(counterValues);
			covIt = begin(coverage);
		}
		TValue c = value(matIt);
		if (c == specialGap) continue;
		else {
			++value(covIt);
			if (c == gapChar) ++value(value(countIt), alphabetSize);
			else ++value(value(countIt), ordValue(TAlphabet(c)));
		}
	}

	// Get the consensus
	TCounterIt countItEnd = end(counterValues);
	countIt = begin(counterValues);
	clear(consensus);
	clear(gappedConsensus);
	for(;countIt != countItEnd; goNext(countIt)) {
		TSize max = 0;
		TValue c = TValue();
		TCIt cIt = begin(value(countIt));
		TCIt cItEnd = end(value(countIt));
		TSize pos = 0;
		for(;cIt != cItEnd; goNext(cIt), ++pos) {
			if (value(cIt) > max) {
				max = value(cIt);
				if (pos == alphabetSize) c = gapChar;
				else c = TAlphabet(pos);
			}
		}
		if (c != gapChar) appendValue(consensus, c);
		appendValue(gappedConsensus, c);
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet, typename TCargo, typename TSpec, typename TBegEndRowPos, typename TNewLibraryGraph>
inline unsigned int
realignLowQualityReads(Graph<Alignment<TStringSet, TCargo, TSpec> > const& gIn,
					   TBegEndRowPos const& readBegEndRowPos,
					   TNewLibraryGraph& gOut)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TBegEndRowPos>::Type TBegEndIter;

	// Initialization
	TStringSet& str = stringSet(gIn);
	clearVertices(gOut);
	
	// Find disrupted reads
	std::set<TId> unalignedRead;
	TBegEndIter beIt = begin(readBegEndRowPos);
	TBegEndIter beItEnd = end(readBegEndRowPos);
	TSize pos = 0;
	for(;beIt != beItEnd; ++beIt, ++pos) {
		TSize lenStr = length(value(str, pos));
		if (((value(beIt)).i2 - (value(beIt)).i1) > (lenStr + lenStr / 3 + 10)) unalignedRead.insert(positionToId(str, pos));
	}

	// Any disrupted reads
	if (unalignedRead.empty()) return 0;
	else exit(0);

	
	//// String of fragments to combine all pairwise alignments into a multiple alignment
	//typedef Fragment<> TFragment;
	//typedef String<TFragment> TFragmentString;
	//typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	//TFragmentString matches;
	//reserve(matches, numEdges(gIn));

	//// Insert all overlaps from the previous alignment
	//TEdgeIterator it_tmp(gIn);
	//for(;!atEnd(it_tmp);++it_tmp) appendValue(matches, TFragment(sequenceId(gIn, sourceVertex(it_tmp)),fragmentBegin(gIn, sourceVertex(it_tmp)), sequenceId(gIn, targetVertex(it_tmp)), fragmentBegin(gIn, targetVertex(it_tmp)), fragmentLength(gIn, sourceVertex(it_tmp))));
	//
	//// Recompute the overlap of interesting pairs
	//Score<int> score_type = Score<int>(2,-1,-4,-6);
	//TPairIter pairIt = begin(pList);
	//TPairIter pairItEnd = end(pList);
	//for(;pairIt != pairItEnd; ++pairIt) {
	//	TId id1 = (value(pairIt)).i1;
	//	TId id2 = (value(pairIt)).i2;

	//	if ((unalignedRead.find(id1) != unalignedRead.end()) || (unalignedRead.find(id2) != unalignedRead.end())) {
	//		// Make a pairwise string-set
	//		TStringSet pairSet;
	//		assignValueById(pairSet, str, id1);
	//		assignValueById(pairSet, str, id2);

	//		// Overlap alignment with a small mismatch score
	//		TSize from = length(matches);
	//		globalAlignment(matches, pairSet, score_type, AlignConfig<true,true,true,true>(), Gotoh() );

	//		// Determine a sequence weight
	//		TSize matchLen = 0;
	//		TSize overlapLen = 0;
	//		TSize alignLen = 0;
	//		getAlignmentStatistics(matches, pairSet, from, matchLen, overlapLen, alignLen);
	//		double quality = (double) matchLen / (double) overlapLen;

	//		// Take all overlaps of good quality
	//		if ((quality < 0.75) || (matchLen < 8)) {
	//			resize(matches, from);
	//		}
	//	}
	//}

	//// Refine all matches
	//matchRefinement(matches,stringSet(gOut),score_type, gOut);

	//return unalignedRead.size();
}


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
