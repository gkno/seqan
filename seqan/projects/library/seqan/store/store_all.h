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

#ifndef SEQAN_HEADER_STORE_ALL_H
#define SEQAN_HEADER_STORE_ALL_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Contig Store Configuration
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec = void>
struct FragmentStoreConfig 
{
	typedef String<Dna5Q>	TReadSeq;
	typedef String<Dna5Q>	TContigSeq;
	
	typedef double			TMean;
	typedef double			TStd;
	typedef signed char		TMappingQuality;
		
	typedef void					TReadStoreElementSpec;
	typedef Owner<ConcatDirect<> >	TReadSeqStoreSpec;
	typedef void					TMatePairStoreElementSpec;
	typedef void					TLibraryStoreElementSpec;
	typedef void					TContigStoreElementSpec;
	typedef void					TAlignedReadStoreElementSpec;
	typedef Owner<ConcatDirect<> >	TAlignedReadTagStoreSpec;
	typedef void					TAnnotationStoreElementSpec;
};

//////////////////////////////////////////////////////////////////////////////
// Fragment Store
//////////////////////////////////////////////////////////////////////////////

template <typename TSpec = void, typename TConfig = FragmentStoreConfig<TSpec> >
struct FragmentStore
{
private:
	typedef typename TConfig::TReadStoreElementSpec			TReadStoreElementSpec;
	typedef typename TConfig::TReadSeqStoreSpec				TReadSeqStoreSpec;
	typedef typename TConfig::TMatePairStoreElementSpec		TMatePairStoreElementSpec;
	typedef typename TConfig::TLibraryStoreElementSpec		TLibraryStoreElementSpec;
	typedef typename TConfig::TContigStoreElementSpec		TContigStoreElementSpec;
	typedef typename TConfig::TAlignedReadStoreElementSpec	TAlignedReadStoreElementSpec;
	typedef typename TConfig::TAlignedReadTagStoreSpec		TAlignedReadTagStoreSpec;
	typedef typename TConfig::TAnnotationStoreElementSpec	TAnnotationStoreElementSpec;

public:
	typedef typename TConfig::TMean					TMean;
	typedef typename TConfig::TStd					TStd;
	typedef typename TConfig::TMappingQuality		TMappingQuality;
	
	typedef typename TConfig::TReadSeq				TReadSeq;
	typedef typename TConfig::TContigSeq			TContigSeq;
	typedef typename Position<TReadSeq>::Type		TReadPos;
	typedef typename Position<TContigSeq>::Type		TContigPos;
	
	typedef GapAnchor<TReadPos>						TReadGapAnchor;
	typedef GapAnchor<TContigPos>					TContigGapAnchor;

	typedef String< ReadStoreElement< TReadSeq, TReadPos, TReadStoreElementSpec > >							TReadStore;
	typedef String< MatePairStoreElement< TMatePairStoreElementSpec > >										TMatePairStore;
	typedef String< LibraryStoreElement< TMean, TStd, TLibraryStoreElementSpec > >							TLibraryStore;
	typedef String< ContigStoreElement< TContigSeq, TContigGapAnchor, TContigStoreElementSpec > >			TContigStore;
	typedef String< AlignedReadStoreElement< TContigPos, TReadGapAnchor, TAlignedReadStoreElementSpec > >	TAlignedReadStore;
	typedef String< AlignQualityStoreElement< TMappingQuality >	>											TAlignQualityStore;
	typedef StringSet<CharString, TAlignedReadTagStoreSpec>													TAlignedReadTagStore;
	typedef String< AnnotationStoreElement< TContigPos, TAnnotationStoreElementSpec > >						TAnnotationStore;
	typedef StringSet<TReadSeq, TReadSeqStoreSpec>															TReadSeqStore;
	typedef StringSet<CharString>																			TNameStore;
	
	// main containers
	TReadStore			readStore;			// readId     -> matePairId
	TMatePairStore		matePairStore;		// matePairId -> readId0, readId1, libraryId
	TLibraryStore		libraryStore;		// libraryId  -> libSizeMean, libSizeStd
	TContigStore		contigStore;		// contigId   -> contigSeq, contigGaps
	TAlignedReadStore	alignedReadStore;	//            -> id, readId, contigId, pairMatchId (not matePairId!), beginPos, endPos, gaps
	TAnnotationStore	annotationStore;	// annoId     -> parentId, contigId, beginPos, endPos

											// REMARKS: 
											// 1)
											//    beginPos <= endPos     forward strand
											//    beginPos >  endPos     backward strand (reverse complement)
											// 2) 
											//    The alignedReadStore can arbitrarily be resorted. The unique identifier id should
											//    be used to address additional information for each alignedRead in additional tables.

	// we store the read sequences in a seperate stringset to reduce the memory overhead 
	TReadSeqStore		readSeqStore;

	// extra SAM fields
	TAlignQualityStore		alignQualityStore;
	TAlignedReadTagStore	alignedReadTagStore;

	// retrieve the names of reads, mate-pairs, libraries, contigs, annotations by their ids
	TNameStore	readNameStore;
	TNameStore	matePairNameStore;
	TNameStore	libraryNameStore;
	TNameStore	contigNameStore;
	TNameStore	annotationNameStore;
};

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Read Store Accessors
//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TConfig>
inline void
clearReads(FragmentStore<TSpec, TConfig> &me)
{
	clear(me.readStore);
	clear(me.readSeqStore);
}

template <typename TSpec, typename TConfig, typename TRead, typename TId>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TReadStore>::Type
appendRead(
	FragmentStore<TSpec, TConfig> &me, 
	TRead const &read, 
	TId matePairId)
{
	SEQAN_ASSERT(length(me.readStore) == length(me.readSeqStore))

	typedef typename FragmentStore<TSpec, TConfig>::TReadStore TReadStore;
	typename Value<TReadStore>::Type r;
	r.matePairId = matePairId;

	appendValue(me.readStore, r, Generous());
	appendValue(me.readSeqStore, read, Generous());
	return length(me.readStore) - 1;
}

template <typename TSpec, typename TConfig, typename TRead, typename TId>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TReadStore>::Type
appendRead(
	FragmentStore<TSpec, TConfig> &me, 
	TRead const &read, 
	CharString const &name,
	TId matePairId)
{
	SEQAN_ASSERT(length(me.readStore) == length(me.readSeqStore))

	typedef typename FragmentStore<TSpec, TConfig>::TReadStore TReadStore;
	typename Value<TReadStore>::Type r;
	r.matePairId = matePairId;

	appendValue(me.readStore, r, Generous());
	appendValue(me.readSeqStore, read, Generous());
	appendValue(me.readNameStore, name, Generous());
	return length(me.readStore) - 1;
}

template <typename TSpec, typename TConfig, typename TRead>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TReadStore>::Type
appendRead(
	FragmentStore<TSpec, TConfig> &me, 
	TRead const &read)
{
	typedef typename FragmentStore<TSpec, TConfig>::TReadStore TReadStore;
    typedef typename Value<TReadStore>::Type TReadStoreElement;
    
	return appendRead(me, read, TReadStoreElement::INVALID_ID);
}

template <typename TSpec, typename TConfig, typename TRead>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TReadStore>::Type
appendRead(
	FragmentStore<TSpec, TConfig> &me, 
	TRead const &read,
	CharString const &name)
{
	typedef typename FragmentStore<TSpec, TConfig>::TReadStore TReadStore;
    typedef typename Value<TReadStore>::Type TReadStoreElement;
    
	return appendRead(me, read, name, TReadStoreElement::INVALID_ID);
}

template <typename TSpec, typename TConfig, typename TId>
inline typename Value<typename FragmentStore<TSpec, TConfig>::TReadSeqStore>::Type
getRead(
	FragmentStore<TSpec, TConfig> &me, 
	TId id)
{
	return value(me.readSeqStore, id);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TConfig, typename TRead>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TMatePairStore>::Type
appendMatePair(
	FragmentStore<TSpec, TConfig> &me, 
	TRead const &read1, 
	TRead const &read2)
{
	typedef FragmentStore<TSpec, TConfig>			TFragmentStore;
	typedef typename TFragmentStore::TReadStore		TReadStore;
	typedef typename TFragmentStore::TMatePairStore	TMatePairStore;

	typename Value<TReadStore>::Type r;
	typename Value<TMatePairStore>::Type mp;
	r.matePairId = length(me.matePairStore);
	mp.readId[0] = length(me.readStore);
	mp.readId[1] = length(me.readStore) + 1;

	appendValue(me.readStore, r, Generous());
	appendValue(me.readStore, r, Generous());
	appendValue(me.matePairStore, mp, Generous());
	appendValue(me.readSeqStore, read1, Generous());
	appendValue(me.readSeqStore, read2, Generous());
	return length(me.matePairStore) - 1;
}

template <typename TSpec, typename TConfig, typename TRead>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TMatePairStore>::Type
appendMatePair(
	FragmentStore<TSpec, TConfig> &me, 
	TRead const &read1, 
	TRead const &read2, 
	CharString const &name1,
	CharString const &name2)
{
	SEQAN_ASSERT(length(me.readStore) == length(me.readSeqStore))

	typedef FragmentStore<TSpec, TConfig>			TFragmentStore;
	typedef typename TFragmentStore::TReadStore		TReadStore;
	typedef typename TFragmentStore::TMatePairStore	TMatePairStore;

	typename Value<TReadStore>::Type r;
	typename Value<TMatePairStore>::Type mp;
	r.matePairId = length(me.matePairStore);
	mp.readId[0] = length(me.readStore);
	mp.readId[1] = length(me.readStore) + 1;

	appendValue(me.readStore, r, Generous());
	appendValue(me.readStore, r, Generous());
	appendValue(me.matePairStore, mp, Generous());
	appendValue(me.readSeqStore, read1, Generous());
	appendValue(me.readSeqStore, read2, Generous());
	appendValue(me.readNameStore, name1, Generous());
	appendValue(me.readNameStore, name2, Generous());
	return length(me.matePairStore) - 1;
}

//////////////////////////////////////////////////////////////////////////////

// 1. remove aligned reads with invalid ids
// 2. rename ids beginning with 0
template <typename TSpec, typename TConfig>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TAlignedReadStore>::Type
compactAlignedReads(FragmentStore<TSpec, TConfig> &me)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore				TAlignQualityStore;

	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Id<TAlignedRead>::Type							TId;
	typedef typename Size<TAlignQualityStore>::Type					TAQSize;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignedReadIter;
	typedef typename Iterator<TAlignQualityStore, Standard>::Type	TAlignQualityIter;
	
	sortAlignedReads(me.alignedReadStore, SortId());
	
	TAlignedReadIter itAR = begin(me.alignedReadStore, Standard());
	TAlignedReadIter itARend = end(me.alignedReadStore, Standard());
	TAlignQualityIter itAQ = begin(me.alignQualityStore, Standard());
	TAlignQualityIter itAQbegin = itAQbegin;
	TAQSize aqSize = length(me.alignQualityStore);
	TId newId = 0;
	
	for (; itAR != itARend; ++itAR, ++newId)
	{
		TId id = (*itAR).id;
		if (id == TAlignedRead::INVALID_ID) break;	// we assume that invalid ids are at the end of the AlignedReadStore
		if (id < aqSize)
		{
			*itAQ = *(itAQbegin + id);
			++itAQ;
		}
		(*itAR).id = newId;
	}
	
	resize(me.alignedReadStore, newId, Exact());
	resize(me.alignQualityStore, itAQ - itAQbegin, Exact());
	return newId;
}

// rename pair match ids beginning with 0, returns the number of pair matches
template <typename TSpec, typename TConfig>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TAlignedReadStore>::Type
compactPairMatchIds(FragmentStore<TSpec, TConfig> &me)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;

	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Id<TAlignedRead>::Type							TId;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignedReadIter;
	
	sortAlignedReads(me.alignedReadStore, SortPairMatchId());
	
	TAlignedReadIter itAR = begin(me.alignedReadStore, Standard());
	TAlignedReadIter itARend = end(me.alignedReadStore, Standard());
	if (itAR == itARend) return 0;
	
	TId lastId = (*itAR).pairMatchId;
	TId newId = 0;
	for (; itAR != itARend; ++itAR)
	{
		TId id = (*itAR).pairMatchId;
		if (id == TAlignedRead::INVALID_ID) break;	// we assume that invalid ids are at the end of the AlignedReadStore
		if (lastId < id)
		{
			lastId = id;
			++newId;
		}
		(*itAR).pairMatchId = newId;
	}
	return newId + 1;
}

// calculate outer library size for each pair match
template <typename TLibSizeString, typename TSpec, typename TConfig>
inline void
calculateLibSizes(TLibSizeString &libSize, FragmentStore<TSpec, TConfig> &me)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Id<TAlignedRead>::Type							TId;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignedReadIter;
	typedef typename TFragmentStore::TContigPos						TGPos;

	TAlignedReadIter it = begin(me.alignedReadStore, Standard());
	TAlignedReadIter itEnd = end(me.alignedReadStore, Standard());

	resize(libSize, compactPairMatchIds(me), Exact());
	TId lastId = TAlignedRead::INVALID_ID;
	TGPos leftMatePos = 0;
	for (; it != itEnd; ++it)
	{
		TId id = (*it).pairMatchId;
		if (id == TAlignedRead::INVALID_ID) break;	// we assume that invalid ids are at the end of the AlignedReadStore
		if (id != lastId)
		{
			leftMatePos = (*it).beginPos;
			lastId = id;
		} else
			libSize[id] = (*it).beginPos - leftMatePos;
	}
}

// returns mate number (0..left mate, 1..right mate, -1..no mate pair) for a read
template <typename TSpec, typename TConfig, typename TId>
inline signed char
getMateNo(FragmentStore<TSpec, TConfig> const &me, TId readId)
{
	typedef FragmentStore<TSpec, TConfig>			TFragmentStore;
	typedef typename TFragmentStore::TReadStore		TReadStore;
	typedef typename TFragmentStore::TMatePairStore	TMatePairStore;

	typedef typename Value<TReadStore>::Type		TRead;
	typedef typename Value<TMatePairStore>::Type	TMatePair;
	
	if (readId != TRead::INVALID_ID)
	{
		TRead const &r = me.readStore[readId];
		if (r.matePairId != TRead::INVALID_ID)
		{
			TMatePair const &mp = me.matePairStore[r.matePairId];
			if (mp.readId[0] == readId) return 0;
			if (mp.readId[1] == readId) return 1;
		}
	}
	return -1;
}

// calculate index of the other mate for each pair match
template <typename TMateIndexString, typename TSpec, typename TConfig>
inline void
calculateMateIndices(TMateIndexString &mateIndex, FragmentStore<TSpec, TConfig> &me)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Id<TAlignedRead>::Type							TId;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignedReadIter;

	TAlignedReadIter it = begin(me.alignedReadStore, Standard());
	TAlignedReadIter itEnd = end(me.alignedReadStore, Standard());

	for (TId idx = 0; it != itEnd; ++it, ++idx)
	{
		TId id = (*it).pairMatchId;
		if (id == TAlignedRead::INVALID_ID) continue;
		if (length(mateIndex) < 2*id + 2)
			fill(mateIndex, 2*id + 2, TAlignedRead::INVALID_ID, Generous());
		mateIndex[2*id + 1 - getMateNo(me, (*it).readId)] = idx;
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TConfig, typename TScore>
void convertMatchesToGlobalAlignment(FragmentStore<TSpec, TConfig> &store, TScore &score)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;

	typedef typename TFragmentStore::TReadStore						TReadStore;
	typedef typename TFragmentStore::TReadSeqStore					TReadSeqStore;
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore				TAlignQualityStore;
	typedef typename TFragmentStore::TContigStore					TContigStore;

	typedef typename Value<TReadStore>::Type						TRead;
	typedef typename Value<TContigStore>::Type						TContig;
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;

	typedef typename TFragmentStore::TReadSeq						TReadSeq;
	typedef typename Value<TReadSeqStore>::Type						TReadSeqStored;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignedReadIter;
	typedef typename Id<TAlignedRead>::Type							TId;
	typedef typename GetValue<TAlignQualityStore>::Type				TQuality;
	typedef typename TFragmentStore::TContigPos						TContigPos;
	
	typedef typename TContig::TContigSeq							TContigSeq;
	typedef Align<TReadSeq, ArrayGaps>								TAlign;
	typedef Gaps<TReadSeq, ArrayGaps>								TGaps;

	typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >	TContigGaps;
	typedef Gaps<TReadSeq, AnchorGaps<typename TAlignedRead::TGapAnchors> >	TReadGaps;
	typedef typename Iterator<TContigGaps>::Type							TContigIter;
	typedef typename Iterator<TReadGaps>::Type								TReadIter;
	
	// sort matches by increasing begin positions
	sortAlignedReads(store.alignedReadStore, SortBeginPos());
	sortAlignedReads(store.alignedReadStore, SortContigId());

	TReadSeq readSeq;
	TId lastContigId = TAlignedRead::INVALID_ID;
	TAlignedReadIter it = begin(store.alignedReadStore, Standard());
	TAlignedReadIter itEnd = end(store.alignedReadStore, Standard());
	TAlignedReadIter firstOverlap = begin(store.alignedReadStore, Standard());
	for (; it != itEnd; ++it)
	{
		// 1. Calculate pairwise alignment
		TContigPos left = (*it).beginPos;
		TContigPos right = (*it).endPos;
		TContigPos cBegin = _min(left, right);
		TContigPos cEnd = _max(left, right);
		
		readSeq = store.readSeqStore[(*it).readId];
		if (left > right)
			reverseComplementInPlace(readSeq);
		
		TAlign align;
		resize(rows(align), 2);
		assignSource(row(align, 0), infix(store.contigStore[(*it).contigId].seq, cBegin, cEnd));
		assignSource(row(align, 1), readSeq);
		globalAlignment(align, score);		
		std::cout << align << std::endl;

		// 2. Skip non-overlapping matches
		if (lastContigId != (*it).contigId)
		{
			firstOverlap = it;
			lastContigId = (*it).contigId;
		} else
			while (firstOverlap != itEnd && _max((*firstOverlap).beginPos, (*firstOverlap).endPos) <= cBegin)
				++firstOverlap;

		// 3. Iterate over alignment
		TContigGaps	contigGaps(store.contigStore[(*it).contigId].seq, store.contigStore[(*it).contigId].gaps);
		TReadGaps	readGaps(readSeq, (*it).gaps);
		setBeginPosition(contigGaps, positionSeqToGap(contigGaps, cBegin));
		setEndPosition(contigGaps, positionSeqToGap(contigGaps, cEnd));
		
		TContigIter cIt = begin(contigGaps);
		TReadIter rIt = begin(readGaps);
		typename Iterator<TGaps>::Type it1 = begin(row(align, 0));
		typename Iterator<TGaps>::Type it2 = begin(row(align, 1));
		for (; !atEnd(cIt) && !atEnd(it1); goNext(cIt), goNext(rIt))
		{
			bool isGapContig = isGap(cIt);
			bool isGapLocalContig = isGap(it1);
			if (isGap(it2))
			{
				// copy gaps from alignment
				insertGaps(rIt, 1);
			}
			if (isGapContig != isGapLocalContig)
			{
				if (isGapContig)
				{
					// copy exisiting contig gaps
					insertGaps(rIt, 1);
					goNext(rIt);
					continue;
				} else
				{
					// insert padding gaps in contig and reads
					TContigPos cPos = cIt.current.seqPos;
					insertGaps(cIt, 1);
					for (TAlignedReadIter j = firstOverlap; j != it; ++j)
					{
						TContigPos rBegin = _min((*j).beginPos, (*j).endPos);
						TContigPos rEnd = _max((*j).beginPos, (*j).endPos);
						if (rBegin < cPos && cPos < rEnd)
						{
							TReadGaps gaps(store.readSeqStore[(*j).readId], (*j).gaps);
							insertGap(gaps, positionSeqToGap(gaps, cPos - rBegin));
						}
					}
				}
			}
			goNext(it1);
			goNext(it2);
		}
		{
			TAlignedReadIter j = it + 1;
			TAlignedReadIter jBegin = begin(store.alignedReadStore, Standard());
			TContigGaps	contigGaps(store.contigStore[(*it).contigId].seq, store.contigStore[(*it).contigId].gaps);
			std::cout << contigGaps << std::endl;
			while (j != jBegin)
			{
				goPrevious(j);
				TReadGaps gaps(store.readSeqStore[(*j).readId], (*j).gaps);
				if ((*j).beginPos > (*j).endPos)
					reverseComplementInPlace(source(gaps));
				int rBegin = _min((*j).beginPos, (*j).endPos);
				setBeginPosition(gaps, -rBegin);
				setEndPosition(gaps, _unclippedLength(contigGaps) - rBegin);
				std::cout << gaps << std::endl << std::endl;
			}
		}

/*			if (isGap(it1))
			{
				TContigPos cGapPos = positionSeqToGap(contigGaps, cPos);
				cIt = TContigiter(contigGaps, cGapPos);
			}
			if (!isGap(it1)) ++cPos;
			if (!isGap(it2)) ++rPos;
			insertGap(readGaps, 
			dequeue.pushBack(TDequeueValue(cBegin, cEnd))
		}
*/	}
	it = begin(store.alignedReadStore, Standard());
	for(; it != itEnd; ++it)
	{
		TContigGaps	contigGaps(store.contigStore[(*it).contigId].seq, store.contigStore[(*it).contigId].gaps);
		if ((*it).beginPos < (*it).endPos)
		{
			(*it).beginPos = positionSeqToGap(contigGaps, (*it).beginPos);
			(*it).endPos = positionSeqToGap(contigGaps, (*it).endPos - 1) + 1;
		} else
		{
			(*it).beginPos = positionSeqToGap(contigGaps, (*it).beginPos - 1) + 1;
			(*it).endPos = positionSeqToGap(contigGaps, (*it).endPos);
		}
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
