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

//#include <stdio.h>

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Name Store Cache
//////////////////////////////////////////////////////////////////////////////

	template <typename TNameStore, typename TName>
	struct NameStoreLess
	{
		typedef typename Position<TNameStore>::Type TId;

		TNameStore *store;
		TName *name;
		
		NameStoreLess() {}

		NameStoreLess(TNameStore &_store, TName &_name):
			store(&_store),
			name(&_name) {}
		
		template <typename TId>
		inline bool operator() (TId a, TId b) const
		{
			if (a != supremumValue(a))
			{
				if (b != supremumValue(b))
					return (*store)[a] < (*store)[b];
				else
					return (*store)[a] < *name;
			} else
			{
				if (b != supremumValue(b))
					return *name < (*store)[b];
				else
					return false;
			}
		}
	};
	
	template <typename TNameStore, typename TName>
	struct NameStoreCache
	{
		typedef typename Position<TNameStore>::Type TId;
		typedef NameStoreLess<TNameStore, TName> TLess;
		typedef std::set<TId, TLess> TSet;
		
		TSet nameSet;
		TNameStore &store;
		TName name;

		NameStoreCache(TNameStore &_store):
			nameSet(TLess(_store, name)),
			store(_store)
		{
			for (unsigned i = 0; i < length(store); ++i)
				nameSet.insert(i);
		}
	};


//////////////////////////////////////////////////////////////////////////////
// Contig Store Configuration
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
	typedef void					TContigFileSpec;
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
	typedef typename TConfig::TContigFileSpec				TContigFileSpec;
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
	
	typedef AnnotationStoreElement< TContigPos, TAnnotationStoreElementSpec >	TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TId								TAnnotationStoreElementId;

	typedef String< ReadStoreElement< TReadSeq, TReadPos, TReadStoreElementSpec > >							TReadStore;
	typedef String< MatePairStoreElement< TMatePairStoreElementSpec > >										TMatePairStore;
	typedef String< LibraryStoreElement< TMean, TStd, TLibraryStoreElementSpec > >							TLibraryStore;
	typedef String< ContigStoreElement< TContigSeq, TContigGapAnchor, TContigStoreElementSpec > >			TContigStore;
	typedef String< ContigFile< TContigFileSpec > >															TContigFileStore;
	typedef String< AlignedReadStoreElement< TContigPos, TReadGapAnchor, TAlignedReadStoreElementSpec > >	TAlignedReadStore;
	typedef String< AlignQualityStoreElement< TMappingQuality >	>											TAlignQualityStore;
	typedef StringSet<CharString, TAlignedReadTagStoreSpec>													TAlignedReadTagStore;
	typedef String< TAnnotationStoreElement >																TAnnotationStore;
	typedef String< IntervalTree< TContigPos, TAnnotationStoreElementId > >									TIntervalTreeStore;
	typedef StringSet<TReadSeq, TReadSeqStoreSpec>															TReadSeqStore;
	typedef StringSet<CharString>																			TNameStore;
	
	// main containers
	TReadStore			readStore;			// readId       -> matePairId
	TMatePairStore		matePairStore;		// matePairId   -> readId0, readId1, libraryId
	TLibraryStore		libraryStore;		// libraryId    -> libSizeMean, libSizeStd
	TContigStore		contigStore;		// contigId     -> contigSeq, contigGaps, contigFileId
	TContigFileStore	contigFileStore;	// contigFileId -> fileName, firstContigId
	TAlignedReadStore	alignedReadStore;	//              -> id, readId, contigId, pairMatchId (not matePairId!), beginPos, endPos, gaps
	TAnnotationStore	annotationStore;	// annoId       -> parentId, contigId, beginPos, endPos
	TIntervalTreeStore	intervalTreeStore_F;		// treeId (same as contigId)	-> intervalTree (F: forward strand)
	TIntervalTreeStore	intervalTreeStore_R;		// 						(R: reverse complement strand)

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
	
	NameStoreCache<TNameStore, CharString>	readNameStoreCache;
	NameStoreCache<TNameStore, CharString>	contigNameStoreCache;
	
	FragmentStore():
		readNameStoreCache(readNameStore),
		contigNameStoreCache(contigNameStore) {}
};

//////////////////////////////////////////////////////////////////////////////
// refresh
    
    template <typename TNameStore, typename TName>
    inline void
    refresh(NameStoreCache<TNameStore, TName> &cache)
    {
		cache.nameSet.clear();
		for (unsigned i = 0; i < length(cache.store); ++i)
			cache.nameSet.insert(i);
	}
		

//////////////////////////////////////////////////////////////////////////////
// _getIdByName
    
    template <typename TNameStore, typename TName, typename TPos>
    inline bool 
    _getIdByName(TNameStore &store, TName &name, TPos &pos)
    {
        typedef Iterator<StringSet<CharString> >::Type TNameStoreIter;
        
        // Iterator over read names
        for (TNameStoreIter iter = begin(store); iter != end(store); ++iter)
		{
            // if the element was found
            if (name == getValue(iter))
			{
                // set the ID
                pos = position(iter);
                // and break the loop
                return true;
            }
        }
        return false;
    }
	
    template <typename TNameStore, typename TName, typename TPos, typename TContext>
    inline bool 
    _getIdByName(TNameStore &store, TName &name, TPos &pos, TContext &)
	{
		Nothing nothing;
		return _getIdByName(store, name, pos);
	}
    
    template<typename TNameStore, typename TName, typename TPos>
    inline bool 
    _getIdByName(TNameStore &, TName &name, TPos &pos, NameStoreCache<TNameStore, TName> &context)
    {
        typedef Iterator<StringSet<CharString> >::Type TNameStoreIter;
		typedef typename Position<TNameStore>::Type TId;
		typedef NameStoreCache<TNameStore, TName> TNameStoreCache;
		typedef typename TNameStoreCache::TSet TSet;
		
		TSet const &set = context.nameSet;	
		context.name = name;
		
		typename TSet::const_iterator it = set.find(supremumValue<TId>());
		if (it != set.end())
		{
			pos = *it;
			return true;
		}
		return false;
    }
	
//////////////////////////////////////////////////////////////////////////////

    template<typename TNameStore, typename TName>
    inline void
    appendName(TNameStore &store, TName &name)
    {
		appendValue(store, name);
	}
	
    template<typename TNameStore, typename TName, typename TContext>
    inline void
    appendName(TNameStore &store, TName &name, TContext &)
    {
		appendName(store, name);
	}
	
    template<typename TNameStore, typename TName>
    inline void
    appendName(TNameStore &store, TName &name, NameStoreCache<TNameStore, TName> &context)
    {
		appendValue(store, name);
		context.nameSet.insert(length(store) - 1);
    }
	
    
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
	TAlignQualityIter itAQbegin = itAQ;
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

template <typename TLayoutStringSet, typename TSpec, typename TConfig, typename TContigId>
void layoutAlignment(TLayoutStringSet &layout, FragmentStore<TSpec, TConfig> &store, TContigId contigId)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;

	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignedReadIter;

	typedef typename Id<TAlignedRead>::Type							TId;
	typedef typename TFragmentStore::TContigPos						TContigPos;
	
	typedef typename Value<TLayoutStringSet>::Type					TLayoutString;
	typedef typename Iterator<TLayoutStringSet>::Type				TLayoutStringSetIter;
	
	// sort matches by increasing begin positions
//	sortAlignedReads(store.alignedReadStore, SortBeginPos());
//	sortAlignedReads(store.alignedReadStore, SortContigId());

	clear(layout);
	TAlignedReadIter it = begin(store.alignedReadStore, Standard());
	TAlignedReadIter itEnd = end(store.alignedReadStore, Standard());

	for (TId id = 0; it != itEnd; ++it, ++id)
	{
		if ((*it).contigId != (TId)contigId) continue;
		
		TLayoutStringSetIter lit = begin(layout, Standard());
		TLayoutStringSetIter litEnd = end(layout, Standard());
		
		TContigPos beginPos = _min((*it).beginPos, (*it).endPos);
		
		for (; lit != litEnd; ++lit)
		{
			if (empty(*lit)) break;
			TAlignedRead &align = store.alignedReadStore[back(*lit)];
			if (_max(align.beginPos, align.endPos) < beginPos)			// maybe <= would be better
				break;													// but harder to differ between reads borders
		}
			
		if (lit == litEnd)
		{
			TLayoutString s;
			appendValue(s, id);
			appendValue(layout, s);
		} else
			appendValue(*lit, id);
	}
}
	
template <typename TStream, typename TLayoutStringSet, typename TSpec, typename TConfig, typename TContigId, typename TPos, typename TNum>
void printAlignment(
	TStream &stream, 
	TLayoutStringSet &layout, FragmentStore<TSpec, TConfig> &store, 
	TContigId contigId,
	TPos posBegin, TPos posEnd,
	TNum lineBegin, TNum lineEnd,
	unsigned lastRead)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;

	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename TFragmentStore::TContigStore					TContigStore;

	typedef typename Value<TContigStore>::Type						TContig;
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename TFragmentStore::TReadSeq						TReadSeq;
	typedef typename Id<TAlignedRead>::Type							TId;
	typedef typename TFragmentStore::TContigPos						TContigPos;
	typedef typename TContig::TContigSeq							TContigSeq;

	typedef typename Value<TLayoutStringSet>::Type					TLayoutString;
	typedef typename Size<TLayoutStringSet>::Type					TLayoutStringSize;
	typedef typename Iterator<TLayoutStringSet>::Type				TLayoutStringSetIter;
	typedef typename Iterator<TLayoutString>::Type					TLayoutStringIter;

	typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >	TContigGaps;
	typedef Gaps<CharString, AnchorGaps<typename TAlignedRead::TGapAnchors> >	TReadGaps;
	
	if ((TId)contigId < length(store.contigStore))
	{
		TContigGaps	contigGaps(store.contigStore[contigId].seq, store.contigStore[contigId].gaps);
		setBeginPosition(contigGaps, posBegin);
		setEndPosition(contigGaps, posEnd);
		stream << contigGaps << '\n';
	} else
		stream << '\n';
		
	if ((TLayoutStringSize)lineEnd > length(layout)) lineEnd = length(layout);
	if ((TLayoutStringSize)lineBegin >= (TLayoutStringSize)lineEnd) return;

	TLayoutStringSetIter lit = begin(layout, Standard()) + lineBegin;
	TLayoutStringSetIter litEnd = begin(layout, Standard()) + lineEnd;
	TReadSeq readSeq;
	CharString readSeqString;

	for (; lit < litEnd; ++lit)
	{
		TLayoutStringIter itEnd = end(*lit, Standard());
		TLayoutStringIter left = begin(*lit, Standard());
		TLayoutStringIter right = itEnd;
		TLayoutStringIter mid;

		while (left < right)
		{
			mid = left + (right - left) / 2;
			TAlignedRead &align = store.alignedReadStore[*mid];

			if (align.contigId < contigId || (align.contigId == contigId && (TPos)_max(align.beginPos, align.endPos) <= posBegin))
				left = mid + 1;	// what we search is in the right part
			else
				right = mid;	//            ...           left part
		}
		
		TPos cursor = posBegin;
		for (; mid < itEnd; ++mid)
		{
			if (*mid >= lastRead) continue;
			TAlignedRead &align = store.alignedReadStore[*mid];
			if (align.contigId != contigId) break;

			TReadGaps readGaps(readSeqString, align.gaps);
			TContigPos	left = align.beginPos;
			TContigPos	right = align.endPos;
			TContigPos	cBegin = _min(left, right);
			TContigPos	cEnd = _max(left, right);
			
			if ((TPos)cEnd <= posBegin) continue; // shouldn't occur
			if (posEnd <= (TPos)cBegin) break;
			
			readSeq = store.readSeqStore[align.readId];
			if (left > right)
			{
				reverseComplementInPlace(readSeq);
				readSeqString = readSeq;
				toLowerInPlace(readSeqString);
			} else
				readSeqString = readSeq;
			
			if ((TPos)cBegin < posBegin)
				setBeginPosition(readGaps, posBegin - (TPos)cBegin);
			else
				for (; cursor < (TPos)cBegin; ++cursor)
					stream << ' ';
			
			if (posEnd < (TPos)cEnd)
				setEndPosition(readGaps, posEnd - (TPos)cBegin);
			
			stream << readGaps;
			cursor = cEnd;
		}
		stream << '\n';
	}
}

template <typename TSpec, typename TConfig, typename TScore>
void convertMatchesToGlobalAlignment(FragmentStore<TSpec, TConfig> &store, TScore &score)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;

	typedef typename TFragmentStore::TReadStore						TReadStore;
	typedef typename TFragmentStore::TReadSeqStore					TReadSeqStore;
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename TFragmentStore::TContigStore					TContigStore;

	typedef typename Value<TReadStore>::Type						TRead;
	typedef typename Value<TContigStore>::Type						TContig;
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;

	typedef typename TFragmentStore::TReadSeq						TReadSeq;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignedReadIter;
	typedef typename Id<TAlignedRead>::Type							TId;
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
		TContigPos	left = (*it).beginPos;
		TContigPos	right = (*it).endPos;
		TContigPos	cBegin = _min(left, right);
		TContigPos	cEnd = _max(left, right);
		TContigGaps	contigGaps(store.contigStore[(*it).contigId].seq, store.contigStore[(*it).contigId].gaps);
		TReadGaps	readGaps(readSeq, (*it).gaps);
		
		readSeq = store.readSeqStore[(*it).readId];
		if (left > right)
			reverseComplementInPlace(readSeq);
				
		// 1. Calculate pairwise alignment
		TAlign align;
		resize(rows(align), 2);
		assignSource(row(align, 0), infix(store.contigStore[(*it).contigId].seq, cBegin, cEnd));
		assignSource(row(align, 1), readSeq);
		globalAlignment(align, score);

		// 2. Skip non-overlapping matches
		cBegin = positionSeqToGap(contigGaps, cBegin);
		if (lastContigId != (*it).contigId)
		{
			firstOverlap = it;
			lastContigId = (*it).contigId;
		} else
			while (firstOverlap != itEnd && _max((*firstOverlap).beginPos, (*firstOverlap).endPos) <= cBegin)
				++firstOverlap;

		// 3. Iterate over alignment
		setBeginPosition(contigGaps, cBegin);
		
		TContigIter cIt = begin(contigGaps);
		TReadIter rIt = begin(readGaps);
		typename Iterator<TGaps>::Type it1 = begin(row(align, 0));
		typename Iterator<TGaps>::Type it2 = begin(row(align, 1));

		for (; !atEnd(cIt) && !atEnd(it1); goNext(cIt), goNext(rIt))
		{
			bool isGapContig = isGap(cIt);
			bool isGapLocalContig = isGap(it1);
			if (isGapContig != isGapLocalContig)
			{
				if (isGapContig)
				{
					// *** gap in contig of the global alignment ***
					// copy exisiting contig gap
					insertGaps(rIt, 1);
					continue;
				} else
				{
					// *** gap in contig of the pairwise alignment ***
					// insert padding gaps in contig and reads
					TContigPos insPos = cIt.current.gapPos;
					insertGaps(cIt, 1);
					for (TAlignedReadIter j = firstOverlap; j != it; ++j)
					{
						TContigPos rBegin = _min((*j).beginPos, (*j).endPos);
						TContigPos rEnd = _max((*j).beginPos, (*j).endPos);
						if (rBegin < insPos && insPos < rEnd)
						{
							if (rBegin < insPos)
							{
								TReadGaps gaps(store.readSeqStore[(*j).readId], (*j).gaps);
								insertGap(gaps, insPos - rBegin);
							} else
							{
								// shift beginPos if insertion was at the front of the read
								if ((*j).beginPos < (*j).endPos)
									++(*j).beginPos;
								else
									++(*j).endPos;
							}
							// shift endPos as the alignment was elongated or shifted
							if ((*j).beginPos < (*j).endPos)
								++(*j).endPos;
							else
								++(*j).beginPos;
						}
					}
				}
			}
			if (isGap(it2))
			{
				// *** gap in read of pairwise alignment ***
				// copy gaps from alignment
				insertGaps(rIt, 1);
			}
			goNext(it1);
			goNext(it2);
		}

		// store new gap-space alignment borders
		cEnd = cBegin + length(readGaps);
		if (left < right)
		{
			(*it).beginPos = cBegin;
			(*it).endPos = cEnd;
		} else
		{
			(*it).beginPos = cEnd;
			(*it).endPos = cBegin;
		}
/*		
//		if (interesting)
		{
			String<String<unsigned> > layout;
			layoutAlignment(layout, store, (*it).contigId);
			std::cout << store.readNameStore[(*it).readId] << std::endl;
			std::cout << readGaps << '\t' << cBegin << '\t' << cEnd << std::endl << std::endl;
			printAlignment(std::cout, layout, store, (*it).contigId, (int)cBegin-20, (int)cEnd+20, 0, 40, 1 + (it - begin(store.alignedReadStore, Standard())));
//			getc(stdin);
		}
*/
//		if (store.readNameStore[(*it).readId] == "read3305")
//			return;
	}
}

template <typename TSpec, typename TConfig, typename TContigGapsString>
void convertPairWiseToGlobalAlignment(FragmentStore<TSpec, TConfig> &store, TContigGapsString &gaps)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;

	// stores
	typedef typename TFragmentStore::TReadStore						TReadStore;
	typedef typename TFragmentStore::TReadSeqStore					TReadSeqStore;
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename TFragmentStore::TContigStore					TContigStore;

	// store elements
	typedef typename Value<TReadStore>::Type						TRead;
	typedef typename Value<TReadSeqStore>::Type						TReadSeq;
	typedef typename Value<TContigStore>::Type						TContig;
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;

	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignedReadIter;
	typedef typename Id<TAlignedRead>::Type							TId;
	typedef typename TFragmentStore::TContigPos						TContigPos;
	
	typedef typename TContig::TContigSeq							TContigSeq;
	typedef Align<TReadSeq, ArrayGaps>								TAlign;
	typedef Gaps<TReadSeq, ArrayGaps>								TGaps;

	// gap structures
	typedef Gaps</*TContigSeq*/Nothing, AnchorGaps<typename TContig::TGapAnchors> >			TContigGapsGlobal;
	typedef Gaps</*TContigSeq*/Nothing, AnchorGaps<typename Value<TContigGapsString>::Type> >	TContigGapsPW;
	typedef Gaps<TReadSeq, AnchorGaps<typename TAlignedRead::TGapAnchors> >			TReadGaps;
	
	// gap iterators
	typedef typename Iterator<TContigGapsGlobal>::Type								TContigGlobalIter;	
	typedef typename Iterator<TContigGapsPW>::Type									TContigPWIter;
	typedef typename Iterator<TReadGaps>::Type										TReadIter;
	
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
		TContigPos	left = (*it).beginPos;
		TContigPos	right = (*it).endPos;
		TContigPos	cBegin = _min(left, right);
		TContigPos	cEnd = _max(left, right);
		
		// 1. Initialize gap structures
		TContigGapsGlobal	contigGapsGlobal(/*store.contigStore[(*it).contigId].seq, */store.contigStore[(*it).contigId].gaps);
		TContigGapsPW		contigGapsPW(/*store.contigStore[(*it).contigId].seq, */gaps[(*it).id]);
		TReadGaps			readGaps(store.readSeqStore[(*it).readId], (*it).gaps);
		
		// 2. Skip non-overlapping matches
		cBegin = positionSeqToGap(contigGapsGlobal, cBegin);
		if (lastContigId != (*it).contigId)
		{
			firstOverlap = it;
			lastContigId = (*it).contigId;
		} else
			while (firstOverlap != itEnd && _max((*firstOverlap).beginPos, (*firstOverlap).endPos) <= cBegin)
				++firstOverlap;

		// 3. Iterate over alignment
		setBeginPosition(contigGapsGlobal, cBegin);
		
		TContigGlobalIter cIt = begin(contigGapsGlobal);
		TContigPWIter pIt = begin(contigGapsPW);
		TReadIter rIt = begin(readGaps);
		
		for (; /*!atEnd(cIt) && */!atEnd(rIt); goNext(cIt), goNext(rIt))
		{
			bool isGapContig = isGap(cIt);
			bool isGapLocalContig = isGap(pIt);
			if (isGapContig != isGapLocalContig)
			{
				if (isGapContig)
				{
					// *** gap in contig of the global alignment ***
					// copy exisiting contig gap
					insertGaps(rIt, 1);
					continue;
				} else
				{
					// *** gap in contig of the pairwise alignment ***
					// insert padding gaps in contig and reads
					TContigPos insPos = cIt.current.gapPos;
					insertGaps(cIt, 1);
					for (TAlignedReadIter j = firstOverlap; j != it; ++j)
					{
						TContigPos rBegin = _min((*j).beginPos, (*j).endPos);
						TContigPos rEnd = _max((*j).beginPos, (*j).endPos);
						if (rBegin < insPos && insPos < rEnd)
						{
							if (rBegin < insPos)
							{
								TReadGaps gaps(store.readSeqStore[(*j).readId], (*j).gaps);
								insertGap(gaps, insPos - rBegin);
							} else
							{
								// shift beginPos if insertion was at the front of the read
								if ((*j).beginPos < (*j).endPos)
									++(*j).beginPos;
								else
									++(*j).endPos;
							}
							// shift endPos as the alignment was elongated or shifted
							if ((*j).beginPos < (*j).endPos)
								++(*j).endPos;
							else
								++(*j).beginPos;
						}
					}
				}
			}
			goNext(pIt);
		}

		// store new gap-space alignment borders
		cEnd = cBegin + length(readGaps);
		if (left < right)
		{
			(*it).beginPos = cBegin;
			(*it).endPos = cEnd;
		} else
		{
			(*it).beginPos = cEnd;
			(*it).endPos = cBegin;
		}
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
