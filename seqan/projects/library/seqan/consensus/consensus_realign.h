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
  $Id: graph_consensus_library.h 1809 2008-03-31 12:57:59Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_CONSENSUS_REALIGN_H
#define SEQAN_HEADER_CONSENSUS_REALIGN_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TGapAnchor, typename TSpec, typename TGapPos>
inline void
removeGap(AlignedReadStoreElement<TPos, TGapAnchor, TSpec>& alignedRead,
		  TGapPos const gapPos) 
{
	if (gapPos < alignedRead.beginPos) {
		--alignedRead.beginPos; --alignedRead.endPos;
	} else if (gapPos < alignedRead.endPos) {
		--alignedRead.endPos;
		typedef String<TGapAnchor> TGaps;
		typedef typename Iterator<TGaps, Standard>::Type TGapIter;
		TGapIter gapIt = upperBoundGapAnchor(alignedRead.gaps, gapPos - alignedRead.beginPos, SortGapPos() );
		TGapIter gapItEnd = end(alignedRead.gaps, Standard());
		// Note: We might create empty gaps here
		for(;gapIt != gapItEnd; goNext(gapIt)) {
			--(gapIt->gapPos);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlignedReads, typename TSpec, typename TGapPos>
inline void
removeGap(String<TAlignedReads, TSpec>& alignedReadStore,
		  TGapPos const gapPos) 
{
	typedef String<TAlignedReads, TSpec> TAlignedReadStore;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignIter;

	TAlignIter alignIt = begin(alignedReadStore, Standard());
	TAlignIter alignItEnd = end(alignedReadStore, Standard());
	for(;alignIt != alignItEnd; goNext(alignIt)) removeGap(value(alignIt), gapPos);
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TGapAnchor, typename TSpec, typename TGapPos>
inline int
insertGap(AlignedReadStoreElement<TPos, TGapAnchor, TSpec>& alignedRead,
		  TGapPos const gapPos) 
{
	if (gapPos <= alignedRead.beginPos) {
		++alignedRead.beginPos; ++alignedRead.endPos;
		return 0;
	} else if (gapPos < alignedRead.endPos) {
		++alignedRead.endPos;
		typedef String<TGapAnchor> TGaps;
		typedef typename Iterator<TGaps, Standard>::Type TGapIter;
		TGapIter gapIt = lowerBoundGapAnchor(alignedRead.gaps, gapPos - alignedRead.beginPos, SortGapPos() );
		TGapIter gapItEnd = end(alignedRead.gaps, Standard());
		TGapPos insertPos = (gapPos - alignedRead.beginPos);
		if (gapIt == gapItEnd) {
			int gapLen = 0;
			if (gapItEnd != begin(alignedRead.gaps)) {
				goPrevious(gapItEnd);
				gapLen = (int) gapItEnd->gapPos - (int) gapItEnd->seqPos;
			}
			appendValue(alignedRead.gaps, TGapAnchor(insertPos - gapLen, insertPos + 1), Generous());
		}
		else {		
			int gapPrev = 0;
			if (gapIt != begin(alignedRead.gaps)) {
				TGapIter gapPrevious = gapIt;
				goPrevious(gapPrevious);
				gapPrev = ((int) gapPrevious->gapPos - (int) gapPrevious->seqPos);
			}
			// If gap is within an existing gap, extend this gap
			if ((gapIt->gapPos - (((int) gapIt->gapPos - (int) gapIt->seqPos) - gapPrev) <= insertPos) && (gapIt->gapPos >= insertPos)) {
				for(;gapIt != gapItEnd; goNext(gapIt)) {
					++(gapIt->gapPos);
				}
			} else {
				// Otherwise, create a new gap
				TGapAnchor tmp = value(gapIt);
				++tmp.gapPos;
				gapIt->gapPos = insertPos + 1;
				gapIt->seqPos = insertPos - gapPrev;
				do {
					goNext(gapIt);
					TGapAnchor newTmp;
					if (gapIt != gapItEnd) {
						newTmp = value(gapIt);
						++newTmp.gapPos;
						value(gapIt) = tmp;
					} else appendValue(alignedRead.gaps, tmp, Generous() );
					tmp = newTmp;
				} while (gapIt != gapItEnd);
			}
		}
		return 1;
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlignedReads, typename TSpec, typename TGapPos>
inline int
insertGap(String<TAlignedReads, TSpec>& alignedReadStore,
		  TGapPos const gapPos) 
{
	typedef String<TAlignedReads, TSpec> TAlignedReadStore;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignIter;

	int numGaps = 0;
	TAlignIter alignIt = begin(alignedReadStore, Standard());
	TAlignIter alignItEnd = end(alignedReadStore, Standard());
	for(;alignIt != alignItEnd; goNext(alignIt)) numGaps += insertGap(value(alignIt), gapPos);
	return numGaps;
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TConsensus>
inline int
scoreConsensus(TConsensus& consensus)
{
	typedef typename Size<TConsensus>::Type TSize;
	typedef typename Value<TConsensus>::Type TProfileAlphabet;
	typedef typename Iterator<TConsensus, Standard>::Type TConsIter;
	TSize alphSize = ValueSize<TProfileAlphabet>::VALUE;

	// Compute the score
	int score = 0;
	TConsIter itCons = begin(consensus, Standard() );
	TConsIter itConsEnd = end(consensus, Standard() );
	for(;itCons != itConsEnd; goNext(itCons)) {
		TSize maxCount = 0;
		TSize sumCount = 0;
		TSize tmp;
		for(TSize i = 0; i<alphSize; ++i) {
			if ((tmp = (value(itCons)).count[i]) > maxCount) maxCount = tmp;
			sumCount += tmp;
		}
		score += (sumCount - maxCount);
	}
	return score;
}



//////////////////////////////////////////////////////////////////////////////////

template<typename TFragSpec, typename TConfig, typename TAlignedRead, typename TSpec, typename TConsensus, typename TScore, typename TBandwidth>
inline void 
reAlign(FragmentStore<TFragSpec, TConfig>& fragStore,
		String<TAlignedRead, TSpec>& contigReads,
		TConsensus& consensus,
		TScore const& consScore,
		TBandwidth const bandwidth)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef String<TAlignedRead, TSpec> TAlignedReadStore;
	typedef typename TFragmentStore::TReadPos TReadPos;
	typedef typename TFragmentStore::TReadSeq TReadSeq;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename TFragmentStore::TReadGapAnchor TGapAnchor;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignedReadIter;
	typedef typename Iterator<TConsensus, Standard>::Type TConsIter;

	// Initialization
	typedef typename Value<TConsensus>::Type TProfileChar;
	TSize gapPos = ValueSize<TProfileChar>::VALUE - 1;

	// Remove each fragment and realign it to the profile
	TAlignedReadIter alignIt = begin(contigReads, Standard() );
	TAlignedReadIter alignItEnd = end(contigReads, Standard() );
	for(;alignIt != alignItEnd; goNext(alignIt)) {
		//// Debug code
		//for(TSize i = 0; i<length(consensus); ++i) {
		//	std::cout << consensus[i] << std::endl;
		//}
		//std::cout << std::endl;
		//std::cout << consensus << std::endl;
		//TAlignedReadIter debugIt = begin(contigReads, Standard() );
		//TAlignedReadIter debugItEnd = end(contigReads, Standard() );
		//for(;debugIt != debugItEnd; goNext(debugIt)) {
		//	std::cout << debugIt->beginPos << ',' << debugIt->endPos << ':';
		//	typedef typename Iterator<String<TGapAnchor> , Standard>::Type TGapIter;
		//	TGapIter gapIt = begin(debugIt->gaps, Standard());
		//	TGapIter gapItEnd = end(debugIt->gaps, Standard());
		//	for(;gapIt != gapItEnd; goNext(gapIt)) {
		//		std::cout << '(' << gapIt->seqPos << ',' << gapIt->gapPos << ')' << ',';
		//	}
		//	std::cout << std::endl;
		//}

		TSize itConsPos = 0;
		TConsIter itCons = begin(consensus, Standard() );
		TConsIter itConsEnd = end(consensus, Standard() );

		// Initialize the consensus of the band
		TConsensus bandConsensus;
		TConsensus myRead;
		fill(myRead, length(value(fragStore.readStore, alignIt->readId).seq), TProfileChar());
		resize(bandConsensus, 2 * bandwidth + (alignIt->endPos - alignIt->beginPos));
		TConsIter bandConsIt = begin(bandConsensus);
		TConsIter myReadIt = begin(myRead);
		TSize bandConsLen = 0;
		TSize myReadLen = 0;
		TReadPos bandOffset = 0;
		if ((TReadPos) bandwidth < (TReadPos) alignIt->beginPos) {
			bandOffset = alignIt->beginPos - bandwidth;
			goFurther(itCons, bandOffset); itConsPos += bandOffset;
		}
		for(TReadPos iPos = bandOffset; iPos < alignIt->beginPos; goNext(itCons), ++itConsPos, ++iPos) {
			value(bandConsIt) = value(itCons); ++bandConsIt; ++bandConsLen;
		}
		alignIt->beginPos = alignIt->endPos = 0; // So this read is discarded in all gap operations


		// Remove sequence from profile and add to the consensus
		typedef typename Iterator<TReadSeq, Standard>::Type TReadIter;
		TReadIter itRead = begin(value(fragStore.readStore, alignIt->readId).seq, Standard() );
		TReadIter itReadEnd = end(value(fragStore.readStore, alignIt->readId).seq, Standard() );
		typedef typename Iterator<String<TGapAnchor>, Standard>::Type TReadGapsIter;
		TReadGapsIter itGaps = begin(alignIt->gaps, Standard() );
		TReadGapsIter itGapsEnd = end(alignIt->gaps, Standard() );
		TReadPos old = 0;
		int diff = 0;
		TReadPos clippedBeginPos = 0;
		TReadPos clippedEndPos = 0;
		if ((itGaps != itGapsEnd) && (itGaps->gapPos == 0)) {
			old = itGaps->seqPos;
			clippedBeginPos = old;
			goFurther(itRead, old);
			diff -= old;
			goNext(itGaps);
		}
		for(;itGaps != itGapsEnd; goNext(itGaps)) {
			TReadPos limit = itGaps->seqPos;
			int newDiff = (itGaps->gapPos - limit);
			if (diff > newDiff) {
				clippedEndPos = diff - newDiff;
				limit -= clippedEndPos;				
			}
			for(;old < limit; ++old, goNext(itRead)) {
				--(value(itCons)).count[ordValue(value(itRead))];
				if (!empty(value(itCons), gapPos)) {
					value(bandConsIt) = value(itCons); ++bandConsIt; ++bandConsLen;
					++itConsPos;
				} else removeGap(contigReads, itConsPos);
				value(myReadIt).count[0] = ordValue(value(itRead)); ++myReadIt; ++myReadLen;
				goNext(itCons);
			}
			for(;diff < newDiff; ++diff) {
				--(value(itCons)).count[gapPos];
				if (!empty(value(itCons), gapPos)) {
					value(bandConsIt) = value(itCons); ++bandConsIt; ++bandConsLen;
					++itConsPos;
				} else removeGap(contigReads, itConsPos);
				goNext(itCons);
			}
		}
		if (!clippedEndPos) {
			for(;itRead!=itReadEnd;goNext(itRead)) {
				--(value(itCons)).count[ordValue(value(itRead))];
				if (!empty(value(itCons), gapPos)) {
					value(bandConsIt) = value(itCons); ++bandConsIt; ++bandConsLen;
					++itConsPos;
				} else removeGap(contigReads, itConsPos);
				(value(myReadIt)).count[0] = ordValue(value(itRead)); ++myReadIt; ++myReadLen;
				goNext(itCons);
			}
		}

		// Go further up to the bandwidth
		for(TReadPos iPos = 0; ((itCons != itConsEnd) && (iPos < (TReadPos) bandwidth)); goNext(itCons), ++iPos) {
			value(bandConsIt) = value(itCons); ++bandConsIt; ++bandConsLen;
		}
		resize(bandConsensus, bandConsLen);
		resize(myRead, myReadLen);


		// ReAlign the consensus with the sequence
		typedef StringSet<TConsensus, Dependent<> > TStringSet;
		TStringSet pairSet;
		appendValue(pairSet, bandConsensus);
		appendValue(pairSet, myRead);

		typedef String<Fragment<> > TFragmentString;
		TFragmentString matches;
		//globalAlignment(matches, pairSet, consScore, AlignConfig<true,false,false,true>(), -5, 10, BandedGotoh() );
		globalAlignment(matches, pairSet, consScore, AlignConfig<true,false,false,true>(), NeedlemanWunsch() );

		//// Debug code
		//std::cout << pairSet[0] << std::endl;
		//std::cout << pairSet[1] << std::endl;
		//Graph<Alignment<TStringSet, void, WithoutEdgeId> > g(pairSet);
		//std::cout << globalAlignment(g, pairSet, consScore, AlignConfig<true,false,false,true>(), NeedlemanWunsch() ) << std::endl;
		//std::cout << g << std::endl;

		// Add the read back to the consensus and build the new consensus
		TConsensus newConsensus = infix(consensus, 0, bandOffset);
		TConsIter bandIt = begin(bandConsensus, Standard());
		TConsIter bandItEnd = end(bandConsensus, Standard());
		typedef typename Iterator<TFragmentString, Standard>::Type TFragIter;
		TFragIter fragIt = end(matches, Standard() );
		TFragIter fragItEnd = begin(matches, Standard() );
		TReadPos consPos = 0;
		TReadPos readPos = 0;
		TReadPos alignPos = 0;
		clear(alignIt->gaps);
		diff = 0;
		if (clippedBeginPos) {
			appendValue(alignIt->gaps, TGapAnchor(clippedBeginPos, 0), Generous() );
			diff -= clippedBeginPos;
		}
		bool firstMatch = true;
		do {
			goPrevious(fragIt);
			int gapLen = fragIt->begin1 - consPos;
			if (firstMatch) gapLen = 0;
			while(consPos < fragIt->begin1) {
				if (!firstMatch) ++value(bandIt).count[gapPos];
				appendValue(newConsensus, value(bandIt));
				goNext(bandIt); ++consPos; ++alignPos;
			}
			while(readPos < fragIt->begin2) {
				if (gapLen) {
					diff += gapLen;
					appendValue(alignIt->gaps, TGapAnchor(clippedBeginPos + readPos, clippedBeginPos + readPos + diff));
					gapLen = 0;
				}
				int numGaps = insertGap(contigReads, bandOffset + alignPos);
				typedef typename Value<TProfileChar>::Type TAlphabet;
				TProfileChar tmpChar;
				++tmpChar.count[myRead[readPos].count[0]];
				tmpChar.count[gapPos] += numGaps;
				appendValue(newConsensus, tmpChar);
				++readPos; ++alignPos;
			}
			for(TSize i = 0; i<fragIt->len; ++i, goNext(bandIt), ++consPos, ++readPos, ++alignPos) {
				if (firstMatch) {
					firstMatch = false;
					alignIt->beginPos = bandOffset + consPos;
				} else if (gapLen) {
					diff += gapLen;
					appendValue(alignIt->gaps, TGapAnchor(clippedBeginPos + readPos, clippedBeginPos + readPos + diff));
					gapLen = 0;
				}
				++value(bandIt).count[myRead[readPos].count[0]];
				appendValue(newConsensus, value(bandIt));
			}
		} while (fragIt != fragItEnd);
		for(; readPos < length(myRead); ++readPos) {
			int numGaps = insertGap(contigReads, bandOffset + alignPos);
			typedef typename Value<TProfileChar>::Type TAlphabet;
			TProfileChar tmpChar;
			++tmpChar.count[myRead[readPos].count[0]];
			tmpChar.count[gapPos] += numGaps;
			appendValue(newConsensus, tmpChar);
			++alignPos;
		}
		alignIt->endPos = alignIt->beginPos + clippedBeginPos + readPos + diff;
		if (clippedEndPos) {
			diff -= clippedEndPos;
			appendValue(alignIt->gaps, TGapAnchor(clippedBeginPos + readPos + clippedEndPos, clippedBeginPos + readPos + clippedEndPos + diff));
		}
		for(;bandIt != bandItEnd; goNext(bandIt)) appendValue(newConsensus, value(bandIt));
		for(;itCons != itConsEnd; goNext(itCons)) appendValue(newConsensus, value(itCons));
		
		// Update the consensus
		consensus = newConsensus;
	}
}


//////////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TConfig, typename TScore, typename TId, typename TBandwidth>
inline void 
reAlign(FragmentStore<TSpec, TConfig>& fragStore,
		TScore const& consScore,
		TId const contigId,
		TBandwidth const bandwidth)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
	typedef typename TFragmentStore::TReadPos TReadPos;
	typedef typename TFragmentStore::TContigPos TContigPos;
	typedef typename TFragmentStore::TReadSeq TReadSeq;
	typedef typename TFragmentStore::TReadGapAnchor TGapAnchor;
	typedef typename Value<typename TFragmentStore::TReadSeq>::Type TAlphabet;
	
	// Sort the reads according to the begin position
	sortAlignedReads(fragStore.alignedReadStore, SortBeginPos());

	// Stable-sort according to contigId
	sortAlignedReads(fragStore.alignedReadStore, SortContigId());

	// Copy all reads belonging to this contig and reverse complement them if necessary
	TAlignedReadStore contigReads;
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore, Standard>::Type TAlignIter;
	TAlignIter alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	TAlignIter alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	TReadPos maxPos = 0;
	TReadPos minPos = SupremumValue<TReadPos>::VALUE;
	for(;alignIt != alignItEnd; goNext(alignIt)) {
		if ((TId) alignIt->contigId == contigId) {
			if (alignIt->beginPos > alignIt->endPos) {
				reverseComplementInPlace(value(fragStore.readStore, alignIt->readId).seq);
				typename Value<TAlignedReadStore>::Type alignedEl = value(alignIt);
				TReadPos tmp = alignedEl.beginPos;
				alignedEl.beginPos = alignedEl.endPos;
				alignedEl.endPos = tmp;
				if (alignedEl.beginPos < minPos) minPos = alignedEl.beginPos;
				if (alignedEl.endPos > maxPos) maxPos = alignedEl.endPos;
				appendValue(contigReads, alignedEl, Generous() );
			} else {
				if (alignIt->beginPos < minPos) minPos = alignIt->beginPos;
				if (alignIt->endPos > maxPos) maxPos = alignIt->endPos;
				appendValue(contigReads, value(alignIt), Generous() );
			}
		}
	}

	// Create the consensus sequence
	typedef ModifiedAlphabet<TAlphabet, ModExpand<'-'> > TProfileChar;
	TSize gapPos = ValueSize<TAlphabet>::VALUE;
	typedef ProfileType<TProfileChar> TProfile;
	typedef String<TProfile> TProfileString;
	typedef typename Iterator<TProfileString, Standard>::Type TConsIter;
	TProfileString consensus;
	fill(consensus, maxPos - minPos, TProfile());
	TConsIter itCons = begin(consensus, Standard() );
	TAlignIter contigReadsIt = begin(contigReads, Standard() );
	TAlignIter contigReadsItEnd = end(contigReads, Standard() );
	for(;contigReadsIt != contigReadsItEnd; goNext(contigReadsIt)) {
		contigReadsIt->beginPos -= minPos;
		contigReadsIt->endPos -= minPos;
		itCons = begin(consensus, Standard() );
		goFurther(itCons, contigReadsIt->beginPos);

		typedef typename Iterator<TReadSeq, Standard>::Type TReadIter;
		TReadIter itRead = begin(value(fragStore.readStore, contigReadsIt->readId).seq, Standard() );
		TReadIter itReadEnd = end(value(fragStore.readStore, contigReadsIt->readId).seq, Standard() );
		typedef typename Iterator<String<typename TFragmentStore::TReadGapAnchor>, Standard>::Type TReadGapsIter;
		TReadGapsIter itGaps = begin(contigReadsIt->gaps, Standard() );
		TReadGapsIter itGapsEnd = end(contigReadsIt->gaps, Standard() );

		TReadPos old = 0;
		int diff = 0;
		bool clippedEnd = false;
		if ((itGaps != itGapsEnd) && (itGaps->gapPos == 0)) {
			old = itGaps->seqPos;
			goFurther(itRead, old);
			diff -= old;
			goNext(itGaps);
		}
		for(;itGaps != itGapsEnd; goNext(itGaps)) {
			TReadPos limit = itGaps->seqPos;
			int newDiff = (itGaps->gapPos - limit);
			if (diff > newDiff) {
				limit -= (diff - newDiff);
				clippedEnd = true;
			}
			for(;old < limit; ++old, goNext(itRead)) ++(value(itCons++)).count[ordValue(value(itRead))];
			for(;diff < newDiff; ++diff) ++(value(itCons++)).count[gapPos];
		}
		if (!clippedEnd) {
			for(;itRead!=itReadEnd;goNext(itRead)) ++(value(itCons++)).count[ordValue(value(itRead))];
		}
	}


	reAlign(fragStore, contigReads, consensus, consScore, bandwidth);
	int score = scoreConsensus(consensus);
	int oldScore = score + 1;
	while(score < oldScore) {
		std::cout << "Score: " << score << std::endl;
		oldScore = score;
		reAlign(fragStore, contigReads, consensus, consScore, bandwidth);
		score = scoreConsensus(consensus);
	}

	// Update all the aligned reads and the new consensus
	alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	TAlignIter contigReadIt = begin(contigReads);
	for(;alignIt != alignItEnd; goNext(alignIt)) {
		if ((TId) alignIt->contigId == contigId) {
			if (alignIt->beginPos > alignIt->endPos) {
				reverseComplementInPlace(value(fragStore.readStore, alignIt->readId).seq);
				alignIt->beginPos = contigReadIt->endPos;
				alignIt->endPos = contigReadIt->beginPos;
			} else {
				alignIt->beginPos = contigReadIt->beginPos;
				alignIt->endPos = contigReadIt->endPos;
			}
			// Remove empty gap anchors
			clear(alignIt->gaps);
			typedef typename Iterator<TGapAnchor, Standard>::Type TGapIter;
			TGapIter gapIt = begin(contigReadIt->gaps, Standard());
			TGapIter gapItEnd = end(contigReadIt->gaps, Standard());
			int diff = 0;
			for(;gapIt != gapItEnd; goNext(gapIt)) {
				if ((int) gapIt->gapPos - (int) gapIt->seqPos != diff) {
					diff = (int) gapIt->gapPos - (int) gapIt->seqPos;
					appendValue(alignIt->gaps, value(gapIt));
				}
			}
			goNext(contigReadIt);
		}
	}
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
	TContigElement& contigEl = value(fragStore.contigStore, contigId);
	typedef typename Iterator<TProfileString, Standard>::Type TConsIter;
	TConsIter itConsensus = begin(consensus, Standard());
	TConsIter itConsensusEnd = end(consensus, Standard());
	char gapChar = gapValue<char>();
	TSize gapLen = 0;
	TContigPos contigPos = 0;
	int diff = 0;
	clear(contigEl.seq);
	clear(contigEl.gaps);
	for(;itConsensus != itConsensusEnd; goNext(itConsensus), ++contigPos) {
		if ((char) value(itConsensus) == gapChar) ++gapLen;
		else {
			if (gapLen) {
				diff += (int) gapLen;
				appendValue(contigEl.gaps, TGapAnchor(contigPos - diff, contigPos));
				gapLen = 0;
			}
			appendValue(contigEl.seq, value(itConsensus));
		}

	}


}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
