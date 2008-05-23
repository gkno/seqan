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
// Read alignment and Consensus Generation
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet, typename TCargo, typename TSpec, typename TValue, typename TBegEndRowPos, typename TCoverage, typename TGappedConsensus>
inline void
consensusAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
				   String<TValue>& mat,
				   TBegEndRowPos& readBegEndRowPos,
				   TCoverage& coverage,
				   TGappedConsensus& gappedConsensus)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Value<TString>::Type TAlphabet;
	typedef typename Infix<TString>::Type TInfix;
	typedef typename Iterator<TInfix>::Type TInfixIter;
	typedef std::map<std::pair<TIdType, TIdType>, TVertexDescriptor> TPosToVertexMap;
	typedef std::map<unsigned int, unsigned int> TComponentLength;
	
	// Strongly Connected Components, topological sort, and length of each component
	String<unsigned int> component;
	String<unsigned int> order;
	TComponentLength compLength;
	if (convertAlignment(g, component, order, compLength)) {
		unsigned int numOfComponents = length(order);
		TStringSet& strSet = stringSet(g);
		TSize nseq = length(strSet);
		clear(gappedConsensus);
		clear(coverage);
		clear(readBegEndRowPos);
		resize(readBegEndRowPos, nseq);
		
		// Assign to each sequence the start and end (in terms of component ranks)
		typedef std::map<unsigned int, unsigned int> TComponentToRank;
		TComponentToRank compToRank;
		for(unsigned int compIndex = 0; compIndex < numOfComponents; ++compIndex) {
			compToRank.insert(std::make_pair(order[compIndex], compIndex));
		}
		typedef Pair<unsigned int, unsigned int> TRankPair;
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
		String<unsigned int> seqToRow;
		resize(seqToRow, nseq);
		TSize maxCoverage = 0;
		typedef std::set<unsigned int> TLeftOver;
		TLeftOver leftOver;
		for(unsigned int i=0;i<nseq; ++i) {
			leftOver.insert(i);
		}
		while(!leftOver.empty()) {
			typedef std::set<std::pair<unsigned int, unsigned int> > TSeqToBegin;
			TSeqToBegin seqToBegin;
			for(typename TLeftOver::const_iterator pos = leftOver.begin(); pos != leftOver.end(); ++pos) {
				seqToBegin.insert(std::make_pair((seqToRank[*pos]).i1, *pos));
			}
			unsigned int endPos = 0;
			for(typename TSeqToBegin::const_iterator s = seqToBegin.begin(); s != seqToBegin.end();++s) {
				if (endPos <= (*s).first) {
					unsigned int currentSeq = (*s).second;
					seqToRow[currentSeq] = maxCoverage;
					endPos = (seqToRank[currentSeq]).i2 + 2;
					leftOver.erase(currentSeq);
				}	
			}
			++maxCoverage;
		}


		// Create the matrix
		TSize len = 0;
		String<unsigned int> compOffset;
		fill(compOffset, numOfComponents, 0);
		for(unsigned int compIndex = 0; compIndex < numOfComponents; ++compIndex) {
			compOffset[order[compIndex]] = len;
			len+=compLength[order[compIndex]];
		}
		TValue gapChar = gapValue<TValue>();
		TValue specialGap = '.';
		fill(mat, len * maxCoverage, gapChar);
		resize(coverage, len);
		resize(gappedConsensus, len);

		// Fill in the segments
		unsigned int alphabetSize = ValueSize<TAlphabet>::VALUE;
		String<String<unsigned int> > counterValues;
		resize(counterValues, len);
		for(unsigned int i=0;i<len; ++i) {
			String<unsigned int> counter;
			fill(counter, alphabetSize, 0);
			value(counterValues, i) = counter;
		}
		for(typename TPosToVertexMap::const_iterator it = g.data_pvMap.begin();it != g.data_pvMap.end(); ++it) {
			TInfix str = label(g,it->second);
			unsigned int c = property(component, it->second);
			unsigned int strPos = idToPosition(strSet, it->first.first);
			unsigned int row = seqToRow[strPos];
			//if (row == 0) {
			//	std::cout << sequenceId(g, it->second) << ':' << str << ',' << strSet[sequenceId(g, it->second)] << std::endl;
			//	std::cout << getProperty(component, it->second) << ',' << order[compIndex] << std::endl;
			//	std::cout << (seqToRank[sequenceId(g, it->second)]).i1 << ',' << (seqToRank[sequenceId(g, it->second)]).i2 << std::endl;
			//}
			TInfixIter sIt = begin(str);
			TInfixIter sItEnd = end(str);
			unsigned int i = compOffset[c];
			for(unsigned int pCol = i;sIt!=sItEnd;goNext(sIt), ++pCol, ++i) {
				assignValue(mat, row * len + pCol, *sIt);
				++((counterValues[i])[(unsigned int) *sIt]);
			}
		}
		String<bool> active;
		for(unsigned int compIndex = 0; compIndex < numOfComponents; ++compIndex) {
			unsigned int offset = compOffset[order[compIndex]];
			unsigned int currentCompLength = compLength[order[compIndex]];

			clear(active);
			fill(active, maxCoverage, false);

			// Find the empty rows
			unsigned int activeRows = 0;
			for(unsigned int i=0;i<nseq; ++i) {
				if (((seqToRank[i]).i1 <= compIndex) && ((seqToRank[i]).i2 >= compIndex)) {
					++activeRows;
					active[(seqToRow[i])] = true;
				}
			}
			
			// Substitute false gaps with special gap character
			for(unsigned int i = 0; i < maxCoverage; ++i) {
				if (!(active[i])) {
					for(unsigned int pCol = offset;pCol < offset + currentCompLength;++pCol) assignValue(mat, i * len + pCol, specialGap);
				}
			}

			// Build consensus
			for(unsigned int i=offset;i<offset+currentCompLength; ++i) {
				TSize max = 0;
				int index_max = -1;
				TSize max2nd = 0;
				int index_max2nd = -1;
				TSize total_count = 0;
				for(TSize j = 0; j < length(counterValues[i]); ++j) {
					if ((counterValues[i])[j] > max) {
						if (index_max != -1) {
							index_max2nd = index_max;
							max2nd = max;
						}
						max = (counterValues[i])[j];
						index_max = j;
						
					} else if ((counterValues[i])[j] > max2nd) {
						max2nd = (counterValues[i])[j];
						index_max2nd = j;
					}	
					total_count += (counterValues[i])[j];
				}
				// Coverage
				value(coverage, i) = activeRows;
				if ((counterValues[i])[index_max] > (activeRows - total_count)) value(gappedConsensus, i) = TAlphabet((Byte) index_max);
				else value(gappedConsensus, i) = gapChar;
			}
		}

		// Get the new begin and end positions
		for(unsigned int i=0;i<nseq; ++i) {
			TVertexDescriptor lastVertex = findVertex(const_cast<TGraph&>(g), positionToId(strSet, i), length(strSet[i]) - 1);
			unsigned int readBegin = compOffset[getProperty(component, findVertex(const_cast<TGraph&>(g), positionToId(strSet, i), 0))];
			unsigned int readEnd = compOffset[getProperty(component, lastVertex)] + fragmentLength(const_cast<TGraph&>(g), lastVertex);
			readBegEndRowPos[i].i1 = readBegin;
			readBegEndRowPos[i].i2 = readEnd;
			readBegEndRowPos[i].i3 = seqToRow[i];
		}
	}
}



//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet, typename TCargo, typename TSpec, typename TPairList, typename TBegEndRowPos, typename TNewLibraryGraph>
inline unsigned int
realignLowQualityReads(Graph<Alignment<TStringSet, TCargo, TSpec> > const& gIn,
					   TPairList const& pList,
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
	typedef typename Iterator<TPairList>::Type TPairIter;

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
	
	// String of fragments to combine all pairwise alignments into a multiple alignment
	typedef Fragment<> TFragment;
	typedef String<TFragment> TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;
	reserve(matches, numEdges(gIn));

	// Insert all overlaps from the previous alignment
	TEdgeIterator it_tmp(gIn);
	for(;!atEnd(it_tmp);++it_tmp) appendValue(matches, TFragment(sequenceId(gIn, sourceVertex(it_tmp)),fragmentBegin(gIn, sourceVertex(it_tmp)), sequenceId(gIn, targetVertex(it_tmp)), fragmentBegin(gIn, targetVertex(it_tmp)), fragmentLength(gIn, sourceVertex(it_tmp))));
	
	// Recompute the overlap of interesting pairs
	Score<int> score_type = Score<int>(2,-1,-4,-6);
	TPairIter pairIt = begin(pList);
	TPairIter pairItEnd = end(pList);
	for(;pairIt != pairItEnd; ++pairIt) {
		TId id1 = (value(pairIt)).i1;
		TId id2 = (value(pairIt)).i2;

		if ((unalignedRead.find(id1) != unalignedRead.end()) || (unalignedRead.find(id2) != unalignedRead.end())) {
			// Make a pairwise string-set
			TStringSet pairSet;
			assignValueById(pairSet, str, id1);
			assignValueById(pairSet, str, id2);

			// Overlap alignment with a small mismatch score
			TSize from = length(matches);
			globalAlignment(matches, pairSet, score_type, AlignConfig<true,true,true,true>(), Gotoh() );

			// Determine a sequence weight
			TSize matchLen = 0;
			TSize overlapLen = 0;
			TSize alignLen = 0;
			getAlignmentStatistics(matches, pairSet, from, matchLen, overlapLen, alignLen);
			double quality = (double) matchLen / (double) overlapLen;

			// Take all overlaps of good quality
			if ((quality < 0.75) || (matchLen < 8)) {
				resize(matches, from);
			}
		}
	}

	// Refine all matches
	matchRefinement(matches,stringSet(gOut),score_type, gOut);

	return unalignedRead.size();
}


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
