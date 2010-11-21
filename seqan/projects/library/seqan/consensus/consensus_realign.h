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
	typedef String<TGapAnchor> TGaps;
	typedef typename Iterator<TGaps, Standard>::Type TGapIter;
	if (gapPos < (TGapPos) alignedRead.beginPos) {
		--alignedRead.beginPos; --alignedRead.endPos;
	} else if (gapPos < (TGapPos) alignedRead.endPos) {
		--alignedRead.endPos;
		TGapIter gapIt = upperBoundGapAnchor(alignedRead.gaps, gapPos - alignedRead.beginPos, SortGapPos() );
		TGapIter gapItEnd = end(alignedRead.gaps, Standard());
		// Note: We might create empty gaps here
		for(;gapIt != gapItEnd; ++gapIt) 
			--(gapIt->gapPos);
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
	for(;alignIt != alignItEnd; ++alignIt) 
		removeGap(*alignIt, gapPos);
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TGapAnchor, typename TSpec, typename TGapPos>
inline int
insertGap(AlignedReadStoreElement<TPos, TGapAnchor, TSpec>& alignedRead,
		  TGapPos const gapPos) 
{
	typedef String<TGapAnchor> TGaps;
	typedef typename Iterator<TGaps, Standard>::Type TGapIter;
	
	if (gapPos <= alignedRead.beginPos) {
		++alignedRead.beginPos; ++alignedRead.endPos;
		return 0;
	} else if (gapPos < alignedRead.endPos) {
		++alignedRead.endPos;
		TGapIter gapIt = lowerBoundGapAnchor(alignedRead.gaps, gapPos - alignedRead.beginPos, SortGapPos() );
		TGapIter gapItEnd = end(alignedRead.gaps, Standard());
		TGapPos insertPos = (gapPos - alignedRead.beginPos);
		if (gapIt == gapItEnd) {
			int gapLen = 0;
			if (gapItEnd != begin(alignedRead.gaps)) {
				--gapItEnd;
				gapLen = (int) gapItEnd->gapPos - (int) gapItEnd->seqPos;
			}
			appendValue(alignedRead.gaps, TGapAnchor(insertPos - gapLen, insertPos + 1), Generous());
		}
		else {		
			int gapPrev = 0;
			if (gapIt != begin(alignedRead.gaps)) {
				TGapIter gapPrevious = gapIt;
				--gapPrevious;
				gapPrev = (int) gapPrevious->gapPos - (int) gapPrevious->seqPos;
			}
			// If gap is within an existing gap, extend this gap
			if ((gapIt->gapPos - (((int) gapIt->gapPos - (int) gapIt->seqPos) - gapPrev) <= insertPos) && (gapIt->gapPos >= insertPos)) {
				for(;gapIt != gapItEnd; ++gapIt) 
					++(gapIt->gapPos);
			} else {
				// Otherwise, create a new gap
				TGapAnchor tmp = value(gapIt);
				++tmp.gapPos;
				gapIt->gapPos = insertPos + 1;
				gapIt->seqPos = insertPos - gapPrev;
				do {
					++gapIt;
					TGapAnchor newTmp;
					if (gapIt != gapItEnd) {
						newTmp = *gapIt;
						++newTmp.gapPos;
						*gapIt = tmp;
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
	for(;alignIt != alignItEnd; ++alignIt) 
		numGaps += insertGap(*alignIt, gapPos);
	return numGaps;
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TConsensus>
inline int
scoreConsensus(TConsensus& consensus)
{
	typedef typename Size<TConsensus>::Type TSize;
	typedef typename Iterator<TConsensus, Standard>::Type TConsIter;

	// Compute the score
	int score = 0;
	TConsIter itCons = begin(consensus, Standard() );
	TConsIter itConsEnd = end(consensus, Standard() );
	TSize maxCount = 0;
	TSize sumCount = 0;
	TSize tmp;
	for(;itCons != itConsEnd; ++itCons) {
		maxCount = 0; sumCount = 0;
		for(TSize i = 0; i < ValueSize<typename Value<TConsensus>::Type>::VALUE; ++i) {
			if ((tmp = (*itCons).count[i]) > maxCount) maxCount = tmp;
			sumCount += tmp;
		}
		score += (sumCount - maxCount);
	}
	return score;
}



//////////////////////////////////////////////////////////////////////////////////

template<typename TFragSpec, typename TConfig, typename TAlignedRead, typename TSpec, typename TConsensus, typename TScore, typename TMethod, typename TBandwidth>
inline void 
reAlign(FragmentStore<TFragSpec, TConfig>& fragStore,
		String<TAlignedRead, TSpec>& contigReads,
		TConsensus& consensus,
		TScore& consScore,
		TMethod const rmethod,
		TBandwidth const bandwidth,
		bool includeReference)
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
	if (includeReference) --alignItEnd;
	TConsensus bandConsensus;
	TConsensus myRead;
	TConsensus newConsensus;
	for(;alignIt != alignItEnd; ++alignIt) {
		//// Debug code
		//for(TSize i = 0; i<length(consensus); ++i) {
		//	std::cout << consensus[i] << std::endl;
		//}
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
		clear(myRead);
		fill(myRead, length(fragStore.readSeqStore[alignIt->readId]), TProfileChar());
		resize(bandConsensus, 2 * bandwidth + (alignIt->endPos - alignIt->beginPos), Generous());
		TConsIter bandConsIt = begin(bandConsensus);
		TConsIter myReadIt = begin(myRead);
		TReadPos bandOffset = 0;
		if (bandwidth < (TBandwidth) alignIt->beginPos) {
			bandOffset = alignIt->beginPos - bandwidth;
			itCons += bandOffset; itConsPos += bandOffset;
			SEQAN_ASSERT_LEQ(itCons, itConsEnd);
		}
		int leftDiag = (alignIt->beginPos - bandOffset) - bandwidth;
		int rightDiag = leftDiag + 2 * bandwidth;
		int increaseBand = 0;
		int removedBeginPos = 0;
		int removedEndPos = 0;
		for(TReadPos iPos = bandOffset; iPos < alignIt->beginPos && itCons != itConsEnd; ++itCons, ++bandConsIt, ++itConsPos, ++iPos)
			*bandConsIt = *itCons; // fill in positions left of readbegin
		TSize itConsPosBegin = itConsPos;  // start position of read basically, right? if(itConsPosBegin != alignIt->beginPos) std::cout <<"nicht unbedingt gleich\n";
		alignIt->beginPos = alignIt->endPos = 0; // So this read is discarded in all gap operations


		// Remove sequence from profile and add to the consensus
		typedef typename Iterator<TReadSeq, Standard>::Type TReadIter;
		TReadIter itRead = begin(fragStore.readSeqStore[alignIt->readId], Standard() );
		TReadIter itReadEnd = end(fragStore.readSeqStore[alignIt->readId], Standard() );
		typedef typename Iterator<String<TGapAnchor>, Standard>::Type TReadGapsIter;
		TReadGapsIter itGaps = begin(alignIt->gaps, Standard() );
		TReadGapsIter itGapsEnd = end(alignIt->gaps, Standard() );
		TReadPos old = 0;
		int diff = 0;
		TReadPos clippedBeginPos = 0;
		TReadPos clippedEndPos = 0;
		SEQAN_ASSERT_LT(itRead, itReadEnd);
		if ((itGaps != itGapsEnd) && (itGaps->gapPos == 0)) {
			old = itGaps->seqPos;
			clippedBeginPos = old;
			itRead += old;
			diff -= old;
			++itGaps;
	        	SEQAN_ASSERT_LT(itRead, itReadEnd);
		}
		for(;itGaps != itGapsEnd && itCons != itConsEnd; ++itGaps) {
			// limit should never be larger than read length 
			TReadPos limit = itGaps->seqPos;
			SEQAN_ASSERT_LT(itGaps->seqPos, (TReadPos)length(fragStore.readSeqStore[alignIt->readId]));
			int newDiff = (itGaps->gapPos - limit);
			SEQAN_ASSERT_LT(itGaps->gapPos, (TReadPos)length(consensus));
			if (diff > newDiff) {
				clippedEndPos = diff - newDiff;
				limit -= clippedEndPos;
			}
			for(;old < limit && itCons != itConsEnd && itRead != itReadEnd; ++old, ++itRead) {
		                //SEQAN_ASSERT_LT(itCons, itConsEnd);
				--(*itCons).count[ordValue(*itRead)];
				if (!empty(*itCons)) {
					*bandConsIt = *itCons; 
					++bandConsIt;
					++itConsPos;
					removedEndPos = 0;
				} else {
					if (itConsPosBegin != itConsPos) {
						++increaseBand;
						++removedEndPos;
					} else ++removedBeginPos;
					removeGap(contigReads, itConsPos);
				}
				(*myReadIt).count[0] = ordValue(*itRead); 
				++myReadIt;
				++itCons;
                		//SEQAN_ASSERT_LT(itRead, itReadEnd);
			}
			for(;diff < newDiff && itCons != itConsEnd; ++diff) {
				++increaseBand;
				//SEQAN_ASSERT_LT(itCons, itConsEnd);
				--(*itCons).count[gapPos];
				if (!empty(*itCons)) {
					*bandConsIt = *itCons; 
					++bandConsIt;
					++itConsPos;
				} else removeGap(contigReads, itConsPos);
				++itCons;
			}
		}
		if (!clippedEndPos) {
			for(;itRead!=itReadEnd && itCons != itConsEnd; ++itRead) {
				//SEQAN_ASSERT_LT(itCons, itConsEnd);
        		        //SEQAN_ASSERT_LT(itRead, itReadEnd);
				--(*itCons).count[ordValue(*itRead)];//subtract the read base to get bandConsensus wo myRead
				if (!empty(*itCons)) {
					*bandConsIt = *itCons; 
					++bandConsIt;
					++itConsPos;
					removedEndPos = 0;
				} else {// only gaps left in this column after removing myRead
					if (itConsPosBegin != itConsPos) {
						++increaseBand;
						++removedEndPos; 
					} else ++removedBeginPos;
					removeGap(contigReads, itConsPos);
				}
				(*myReadIt).count[0] = ordValue(*itRead);   
				++myReadIt;
				++itCons;
			}
		}
		bool singleton = (itConsPosBegin == itConsPos);
		increaseBand -= removedEndPos;

		// Go further up to the bandwidth
		for(TReadPos iPos = 0; ((itCons != itConsEnd) && (iPos < (TReadPos) bandwidth)); ++itCons, ++iPos, ++bandConsIt)
				*bandConsIt = *itCons;
		resize(bandConsensus, bandConsIt - begin(bandConsensus, Standard()), Generous());
		resize(myRead, myReadIt - begin(myRead, Standard()), Generous());

		// ReAlign the consensus with the sequence
		typedef StringSet<TConsensus, Dependent<> > TStringSet;
		TStringSet pairSet;
		appendValue(pairSet, bandConsensus);
		appendValue(pairSet, myRead);

		//for(TSize i = 0; i<length( pairSet[0]); ++i) {
		//	std::cout <<  pairSet[0][i] << std::endl;
		//}
		//std::cout << "_______________" << std::endl;
		//for(TSize i = 0; i<length( pairSet[1]); ++i) {
		//	std::cout <<   pairSet[1][i] << std::endl;
		//}
		//std::cout << "..............." << std::endl;

		typedef String<Fragment<> > TFragmentString;
		TFragmentString matches;
		assignProfile(consScore, bandConsensus);

		
		leftDiag -= removedBeginPos;
		rightDiag -= removedBeginPos;
		if (!singleton) {
			if (rmethod == 0) globalAlignment(matches, pairSet, consScore, AlignConfig<true,false,false,true>(), _max(leftDiag - increaseBand, -1 * (int) length(pairSet[1])), _min(rightDiag + increaseBand, (int) length(pairSet[0])), BandedNeedlemanWunsch());
			else if (rmethod == 1) globalAlignment(matches, pairSet, consScore, AlignConfig<true,false,false,true>(), _max(leftDiag - increaseBand, -1 * (int) length(pairSet[1])), _min(rightDiag + increaseBand, (int) length(pairSet[0])), BandedGotoh());
		}

		//// Debug code
		//Graph<Alignment<TStringSet, void, WithoutEdgeId> > g1(pairSet);
		//int sc1 = globalAlignment(g1, consScore, AlignConfig<true,false,false,true>(), _max(leftDiag - increaseBand, -1 * (int) length(pairSet[1])), _min(rightDiag + increaseBand, (int) length(pairSet[0])), BandedGotoh());
		//std::cout << sc1 << std::endl;
		//std::cout << g1 << std::endl;
		
		// Add the read back to the consensus and build the new consensus
		resize(newConsensus, length(bandConsensus) + length(myRead), Generous());
		TConsIter newConsIt = begin(newConsensus, Standard());
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
		if (fragIt != fragItEnd) {
			do {
				--fragIt;
				int gapLen = fragIt->begin1 - consPos;
				if (firstMatch) gapLen = 0;
				while(consPos < (TReadPos)fragIt->begin1) {
					SEQAN_ASSERT_LT(bandIt, bandItEnd);
					SEQAN_ASSERT_LT(newConsIt, end(newConsensus,Standard()));
					if (!firstMatch) ++(*bandIt).count[gapPos];
					*newConsIt = *bandIt;
					++newConsIt;
					++bandIt; 
					++consPos; 
					++alignPos;
				}
				while(readPos < (TReadPos)fragIt->begin2) {
					SEQAN_ASSERT_LT(readPos, (TReadPos)length(fragStore.readSeqStore[alignIt->readId]));
					SEQAN_ASSERT_LT(newConsIt, end(newConsensus,Standard()));
					if (gapLen) {
						diff += gapLen;
						appendValue(alignIt->gaps, TGapAnchor(clippedBeginPos + readPos, clippedBeginPos + readPos + diff), Generous() );
						gapLen = 0;
					}
					int numGaps = insertGap(contigReads, bandOffset + alignPos);
					TProfileChar tmpChar;
					++tmpChar.count[myRead[readPos].count[0]];
					tmpChar.count[gapPos] += numGaps;
					*newConsIt = tmpChar; ++newConsIt;
					++readPos; ++alignPos;
				}
				for(TSize i = 0; i<fragIt->len; ++i, ++bandIt, ++consPos, ++readPos, ++alignPos, ++newConsIt) {
					SEQAN_ASSERT_LT(bandIt, bandItEnd);
					SEQAN_ASSERT_LT(readPos, (TReadPos)length(fragStore.readSeqStore[alignIt->readId]));
					SEQAN_ASSERT_LT(newConsIt, end(newConsensus,Standard()));
					if (firstMatch) {
						firstMatch = false;
						alignIt->beginPos = bandOffset + consPos;
					} else if (gapLen) {
						diff += gapLen;
						appendValue(alignIt->gaps, TGapAnchor(clippedBeginPos + readPos, clippedBeginPos + readPos + diff), Generous() );
						gapLen = 0;
					}
					SEQAN_ASSERT_LT(bandIt, bandItEnd);
					++(*bandIt).count[myRead[readPos].count[0]];
					*newConsIt = *bandIt; 
				}
			} while (fragIt != fragItEnd);
		}
		for(; readPos < (TReadPos)length(myRead); ++readPos) {
			int numGaps = insertGap(contigReads, bandOffset + alignPos);
			TProfileChar tmpChar;
			++tmpChar.count[myRead[readPos].count[0]];
			tmpChar.count[gapPos] += numGaps;
			SEQAN_ASSERT_LT(newConsIt, end(newConsensus,Standard()));
			*newConsIt = tmpChar; ++newConsIt;
			++alignPos;
		}
		if (singleton) alignIt->beginPos = bandOffset;
		alignIt->endPos = alignIt->beginPos + clippedBeginPos + readPos + diff;
		if (clippedEndPos) {
			diff -= clippedEndPos;
			appendValue(alignIt->gaps, TGapAnchor(clippedBeginPos + readPos + clippedEndPos, clippedBeginPos + readPos + clippedEndPos + diff), Generous() );
		}
		for(;bandIt != bandItEnd; ++bandIt, ++newConsIt) 
		{
			SEQAN_ASSERT_LT(newConsIt, end(newConsensus,Standard()));
			*newConsIt = *bandIt;
		}
		resize(newConsensus, newConsIt - begin(newConsensus, Standard()), Generous());

		infix(consensus, bandOffset, itCons - begin(consensus)) = newConsensus;
	}
}


//////////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TConfig, typename TScore, typename TId, typename TMethod, typename TBandwidth>
inline void 
reAlign(FragmentStore<TSpec, TConfig>& fragStore,
		TScore& consScore,
		TId const contigId,
		TMethod const rmethod,
		TBandwidth const bandwidth,
		bool includeReference)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
	typedef typename TFragmentStore::TReadPos TReadPos;
	
	typedef typename TFragmentStore::TContigStore		TContigStore;
	typedef typename Value<TContigStore>::Type		TContig;
	typedef typename TFragmentStore::TContigPos 		TContigPos;
	typedef typename TFragmentStore::TContigSeq 		TContigSeq;
	typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >	TContigGaps;
	
	typedef typename TFragmentStore::TReadSeq TReadSeq;
	typedef typename TFragmentStore::TReadGapAnchor TGapAnchor;
	typedef typename Value<typename TFragmentStore::TReadSeq>::Type TStoreAlphabet;
	typedef typename BaseAlphabet<TStoreAlphabet>::Type TAlphabet;
	typedef typename Value<TAlignedReadStore>::Type TAlignedElement;
	typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadElement;

	// Sort according to contigId
	sortAlignedReads(fragStore.alignedReadStore, SortContigId());
	
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore, Standard>::Type TAlignIter;
	TAlignIter alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	TAlignIter alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());

	// Sort the reads according to the begin position
	sortAlignedReads(infix(fragStore.alignedReadStore, alignIt - begin(fragStore.alignedReadStore, Standard()), alignItEnd - begin(fragStore.alignedReadStore, Standard())), SortBeginPos());
	alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());

	// Copy all reads belonging to this contig and reverse complement them if necessary
	TAlignedReadStore contigReads;
	TReadPos maxPos = 0;
	TReadPos minPos = SupremumValue<TReadPos>::VALUE;
	for(;alignIt != alignItEnd; ++alignIt) {
		if (alignIt->beginPos > alignIt->endPos) {
			reverseComplementInPlace(fragStore.readSeqStore[alignIt->readId]);
			TAlignedElement alignedEl = *alignIt;
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
	if (includeReference) {
		TId dummyReadId = length(fragStore.readSeqStore);

		appendRead(fragStore, fragStore.contigStore[contigId].seq);
		appendValue(fragStore.readNameStore, fragStore.contigNameStore[contigId], Generous());
		fragStore.contigNameStore[contigId] += "_Consensus";
		TAlignedElement el;
		el.readId = dummyReadId;
		el.contigId = contigId;
		minPos = el.beginPos = 0;
		TContigGaps contigGaps(fragStore.contigStore[contigId].seq, fragStore.contigStore[contigId].gaps);
		maxPos = el.endPos = _max(maxPos,(TReadPos)positionSeqToGap(contigGaps,length(fragStore.contigStore[contigId].seq)-1)+1);
		maxPos = el.endPos = _max(maxPos,(TReadPos)length(fragStore.contigStore[contigId].seq));
		el.gaps = fragStore.contigStore[contigId].gaps;
		appendValue(contigReads, el, Generous() );
		SEQAN_ASSERT_EQ(length(fragStore.readSeqStore[dummyReadId]), length(fragStore.readSeqStore[el.readId]));
	}

	// Create the consensus sequence
	TSize gapPos = ValueSize<TAlphabet>::VALUE;
	typedef ProfileType<TAlphabet> TProfile;
	typedef String<TProfile> TProfileString;
	typedef typename Iterator<TProfileString, Standard>::Type TConsIter;
	TProfileString consensus;
	fill(consensus, maxPos - minPos, TProfile());

	TConsIter itCons = begin(consensus, Standard() );
	TConsIter itConsEnd = end(consensus, Standard());
	TAlignIter contigReadsIt = begin(contigReads, Standard() );
	TAlignIter contigReadsItEnd = end(contigReads, Standard() );
	for(;contigReadsIt != contigReadsItEnd; ++contigReadsIt) {
		contigReadsIt->beginPos -= minPos;
		contigReadsIt->endPos -= minPos;
		itCons = begin(consensus, Standard() );
		itCons += contigReadsIt->beginPos;

		typedef typename Iterator<TReadSeq, Standard>::Type TReadIter;
		TReadIter itRead = begin(fragStore.readSeqStore[contigReadsIt->readId], Standard() );
		TReadIter itReadEnd = end(fragStore.readSeqStore[contigReadsIt->readId], Standard() );
		typedef typename Iterator<String<typename TFragmentStore::TReadGapAnchor>, Standard>::Type TReadGapsIter;
		TReadGapsIter itGaps = begin(contigReadsIt->gaps, Standard() );
		TReadGapsIter itGapsEnd = end(contigReadsIt->gaps, Standard() );

		TReadPos old = 0;
		int diff = 0;
		bool clippedEnd = false;
		if ((itGaps != itGapsEnd) && (itGaps->gapPos == 0)) {
			old = itGaps->seqPos;
			itRead += old;
			diff -= old;
			++itGaps;
		}
		for(;itGaps != itGapsEnd; ++itGaps) {
			TReadPos limit = itGaps->seqPos;
			int newDiff = (itGaps->gapPos - limit);
			SEQAN_ASSERT_LT(itGaps->gapPos, (int)length(consensus));
			if (diff > newDiff) {
				limit -= (diff - newDiff);
				clippedEnd = true;
			}
			for(;old < limit && itRead != itReadEnd && itCons != itConsEnd; ++old, ++itRead) 
			{
				SEQAN_ASSERT_LT(itRead, itReadEnd);
				++(value(itCons++)).count[ordValue(*itRead)];
			}
			for(;diff < newDiff; ++diff) 
				++(value(itCons++)).count[gapPos];
		}
		if (!clippedEnd) {
			for( ; itRead!=itReadEnd && itCons != itConsEnd;++itRead) 
				++(value(itCons++)).count[ordValue(*itRead)];
		}
	}

	
	reAlign(fragStore, contigReads, consensus, consScore, rmethod, bandwidth, includeReference);
	int score = scoreConsensus(consensus);
	int oldScore = score + 1;
	while(score < oldScore) {
	//	std::cout << "Score: " << score << std::endl;
		oldScore = score;
		reAlign(fragStore, contigReads, consensus, consScore, rmethod, bandwidth, includeReference);
		score = scoreConsensus(consensus);
	}
	//std::cout << "FinalScore: " << score << std::endl;

	

	// Update all the aligned reads and the new consensus
	alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	TAlignIter contigReadIt = begin(contigReads, Standard());
	for(;alignIt != alignItEnd; ++alignIt) {
		if (alignIt->beginPos > alignIt->endPos) {
			reverseComplementInPlace(fragStore.readSeqStore[alignIt->readId]);
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
		for(;gapIt != gapItEnd; ++gapIt) {
			if ((int) gapIt->gapPos - (int) gapIt->seqPos != diff) {
				diff = (int) gapIt->gapPos - (int) gapIt->seqPos;
				appendValue(alignIt->gaps, *gapIt, Generous() );
			}
		}
		++contigReadIt;
	}
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
	TContigElement& contigEl = fragStore.contigStore[contigId];
	typedef typename Iterator<TProfileString, Standard>::Type TConsIter;
	TConsIter itConsensus = begin(consensus, Standard());
	TConsIter itConsensusEnd = end(consensus, Standard());
	char gapChar = gapValue<char>();
	TSize gapLen = 0;
	TContigPos contigPos = 0;
	int diff = 0;
	clear(contigEl.seq);
	clear(contigEl.gaps);
	for(;itConsensus != itConsensusEnd; ++itConsensus, ++contigPos) {
		if ((char) *itConsensus == gapChar) ++gapLen;
		else {
			if (gapLen) {
				diff += (int) gapLen;
				appendValue(contigEl.gaps, TGapAnchor(contigPos - diff, contigPos), Generous() );
				gapLen = 0;
			}
			// TODO(weese): Here we convert from ProfileChar<Dna5>->Dna5->Dna5Q
			// instead diverting through Dna5 we could think of directly converting
			// a profile to a quality value, e.g. like the base caller Phred does.
			// Therefore the conversion ProfileChar<Dna5> <-> Dna5Q needs to be
			// defined.
			appendValue(contigEl.seq, (TAlphabet)value(itConsensus), Generous() );
		}

	}
	if (includeReference) 
		appendValue(fragStore.alignedReadStore, contigReads[length(contigReads) - 1], Generous() );
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TConfig, typename TScore, typename TId, typename TBandwidth>
inline void 
reAlign(FragmentStore<TSpec, TConfig>& fragStore,
		TScore& consScore,
		TId const contigId,
		TBandwidth const bandwidth,
		bool includeReference)
{
	reAlign(fragStore, consScore, contigId, 0, bandwidth, includeReference);
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
