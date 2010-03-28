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

/**
.Class.NameStoreCache
..summary:Stores a mapping from names to ids.
..cat:Fragment Store
..signature:FragmentStore<>
..signature:NameStoreCache<TNameStore[, TName]>
..param.TNameStore:The name store to be cached.
...see:Class.FragmentStore
..param.TName:The name type.
...default:$Value<TNameStore>::Type$
...type:Shortcut.CharString

.Memvar.NameStoreCache#NameStoreCache
..summary:Constructor
..signature:NameStoreCache<TNameStore, TName> (nameStore)
..param.nameStore:A name store, e.g. @Memvar.FragmentStore#readNameStore@
...see:Class.FragmentStore
..class:Class.NameStoreCache
*/
	
	template <typename TNameStore, typename TName = typename Value<TNameStore>::Type>
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

/**
.Class.FragmentStore
..summary:Multi-Container to store contigs, reads and a multiple alignment between them.
..cat:Fragment Store
..signature:FragmentStore<>
..signature:FragmentStore<TSpec[, TConfig]>
..param.TSpec:The specializing type.
...default:$void$
..param.TConfig:The configuration struct.
...default:$FragmentStoreConfig<TSpec>$

.Typedef.FragmentStore#TReadStore
..summary:Type of the $readStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TReadSeqStore
..summary:Type of the $readSeqStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TMatePairStore
..summary:Type of the $matePairStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TLibraryStore
..summary:Type of the $libraryStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TContigFileStore
..summary:Type of the $contigFileStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TContigStore
..summary:Type of the $contigStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TAlignedReadStore
..summary:Type of the $alignedReadStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TAnnotationStore
..summary:Type of the $annotationStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TAlignQualityStore
..summary:Type of the $alignQualityStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TAlignedReadTagStore
..summary:Type of the $alignedReadTagStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TReadNameStore
..summary:Type of the $readNameStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TMatePairNameStore
..summary:Type of the $matePairNameStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TLibraryNameStore
..summary:Type of the $libraryNameStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TContigNameStore
..summary:Type of the $contigNameStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TAnnotationNameStore
..summary:Type of the $annotationNameStore$ member.
..class:Class.FragmentStore

.Memvar.FragmentStore#readStore
..summary:String that maps from $readId$ to $<matePairId>$.
..remarks:Value type is @Class.ReadStoreElement@.
..type:Typedef.FragmentStore#TReadStore
..class:Class.FragmentStore
.Memvar.FragmentStore#readSeqStore
..summary:String that maps from $readId$ to $readSeq$.
..type:Typedef.FragmentStore#TReadSeqStore
..class:Class.FragmentStore
.Memvar.FragmentStore#matePairStore
..summary:String that maps from $matePairId$ to $<readId[2], libId>$.
..type:Typedef.FragmentStore#TMatePairStore
..remarks:Value type is @Class.MatePairStoreElement@.
..class:Class.FragmentStore
.Memvar.FragmentStore#libraryStore
..summary:String that maps from $libId$ to $<mean, std>$.
..type:Typedef.FragmentStore#TLibraryStore
..remarks:Value type is @Class.LibraryStoreElement@.
..class:Class.FragmentStore
.Memvar.FragmentStore#contigFileStore
..summary:String that maps from $contigFileId$ to $<fileName, firstContigId>$.
..type:Typedef.FragmentStore#TContigFileStore
..remarks:Value type is @Class.ContigFile@.
..class:Class.FragmentStore
.Memvar.FragmentStore#contigStore
..summary:String that maps from $contigId$ to $<contigSeq, contigGaps, contigFileId>$.
..type:Typedef.FragmentStore#TContigStore
..remarks:Value type is @Class.ContigStoreElement@.
..class:Class.FragmentStore
.Memvar.FragmentStore#alignedReadStore
..summary:String that stores $<alignId, readId, contigId, pairMatchId, beginPos, endPos, gaps>$.
..type:Typedef.FragmentStore#TAlignedReadStore
..remarks:Value type is @Class.AlignedReadStoreElement@.
..class:Class.FragmentStore
.Memvar.FragmentStore#annotationStore
..summary:String that maps from $annoId$ to $<parentId, contigId, beginPos, endPos>$.
..type:Typedef.FragmentStore#TAnnotationStore
..class:Class.FragmentStore
.Memvar.FragmentStore#alignQualityStore
..summary:String that maps from $alignId$ to $<pairScore, score, errors>$.
..type:Typedef.FragmentStore#TAlignQualityStore
..remarks:Value type is @Class.AlignQualityStoreElement@.
..class:Class.FragmentStore
.Memvar.FragmentStore#alignedReadTagStore
..summary:String that maps from $alignId$ to $alignTag$.
..type:Typedef.FragmentStore#TAlignedReadTagStore
..class:Class.FragmentStore
.Memvar.FragmentStore#readNameStore
..summary:String that maps from $readId$ to $readName$.
..type:Typedef.FragmentStore#TReadNameStore
..class:Class.FragmentStore
.Memvar.FragmentStore#matePairNameStore
..summary:String that maps from $contigId$ to $contigName$.
..type:Typedef.FragmentStore#TMatePairNameStore
..class:Class.FragmentStore
.Memvar.FragmentStore#libraryNameStore
..summary:String that maps from $libId$ to $libName$.
..type:Typedef.FragmentStore#TLibraryNameStore
..class:Class.FragmentStore
.Memvar.FragmentStore#contigNameStore
..summary:String that maps from $contigId$ to $contigName$.
..type:Typedef.FragmentStore#TContigNameStore
..class:Class.FragmentStore
.Memvar.FragmentStore#annotationNameStore
..summary:String that maps from $annoId$ to $annoName$.
..type:Typedef.FragmentStore#TAnnotationNameStore
..class:Class.FragmentStore
*/

template <typename TSpec = void, typename TConfig = FragmentStoreConfig<TSpec> >
class FragmentStore
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

	typedef String< ReadStoreElement< TReadStoreElementSpec > >												TReadStore;
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

/**
.Function.refresh:
..summary:Recreate a name store cache.
..cat:Fragment Store
..signature:refresh(cache)
..param.cache:A @Class.NameStoreCache@ object.
...type:Class.NameStoreCache
..see:Function.getIdByName
*/
    
    template <typename TNameStore, typename TName>
    inline void
    refresh(NameStoreCache<TNameStore, TName> &cache)
    {
		cache.nameSet.clear();
		for (unsigned i = 0; i < length(cache.store); ++i)
			cache.nameSet.insert(i);
	}
		

//////////////////////////////////////////////////////////////////////////////
// getIdByName
    
/**
.Function.getIdByName:
..summary:Appends a name to a name store.
..cat:Fragment Store
..signature:getIdByName(nameStore, name, id[, cache])
..param.nameStore:A name store, e.g. @Memvar.FragmentStore#readNameStore@
...see:Class.FragmentStore
..param.name:The name to be searched.
...type:Shortcut.CharString
..param.id:The resulting id.
..param.cache:A structure to efficiently retrieve the id for a given name. If ommited a brute force method is used to search.
...default:Tag.Nothing
...type:Class.NameStoreCache
..returns:$true$ if the name was found and $false$ if not.
..see:Function.getIdByName
*/

    template <typename TNameStore, typename TName, typename TPos>
    inline bool 
    getIdByName(TNameStore &store, TName &name, TPos &pos)
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
    getIdByName(TNameStore &store, TName &name, TPos &pos, TContext &)
	{
		return getIdByName(store, name, pos);
	}
    
    template<typename TNameStore, typename TName, typename TPos>
    inline bool 
    getIdByName(TNameStore &, TName &name, TPos &pos, NameStoreCache<TNameStore, TName> &context)
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

/**
.Function.appendName:
..summary:Appends a name to a name store.
..cat:Fragment Store
..signature:appendRead(nameStore, name[, cache])
..param.nameStore:A name store, e.g. @Memvar.FragmentStore#readNameStore@
...see:Class.FragmentStore
..param.name:The name to be appended.
...type:Shortcut.CharString
..param.cache:A structure to efficiently retrieve the id for a given name. See @Function.getIdByName@.
...default:Tag.Nothing
...type:Class.NameStoreCache
..see:Function.getIdByName
*/

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

/**
.Function.clearReads
..summary:Removes all reads from a fragment store.
..cat:Fragment Store
..signature:clearReads(store)
..param.store:The fragment store.
...type:Class.FragmentStore
..remarks:This function clears the @Memvar.FragmentStore#readStore@, @Memvar.FragmentStore#readSeqStore@ and @Memvar.FragmentStore#readNameStore@.
*/

template <typename TSpec, typename TConfig>
inline void
clearReads(FragmentStore<TSpec, TConfig> &me)
{
	clear(me.readStore);
	clear(me.readSeqStore);
	clear(me.readNameStore);
}

/**
.Function.appendRead:
..summary:Appends a read to a fragment store.
..cat:Fragment Store
..signature:appendRead(store, read[, matePairId])
..signature:appendRead(store, read, name[, matePairId])
..param.store:The fragment store.
...type:Class.FragmentStore
..param.read:The read sequence.
..param.name:The read name.
...type:Shortcut.CharString
..param.matePairId:Id of mate-pair this read is part of.
...default:$INVALID_ID$, which corresponds to an unmated read.
..returns:The $readId$ of the newly appended read.
..remarks:This function appends a single read to the @Memvar.FragmentStore#readStore@ and @Memvar.FragmentStore#readSeqStore@.
If name is given, it is appended to the @Memvar.FragmentStore#readNameStore@.
..see:Function.getRead
*/

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

/**
.Function.getRead
..summary:Returns the read with the given $readId$.
..cat:Fragment Store
..signature:getRead(store, readId)
..param.store:The fragment store.
...type:Class.FragmentStore
..param.readId:The read id.
..returns:The sequence of the read with id $readId$ from the @Memvar.FragmentStore#readSeqStore@.
*/

template <typename TSpec, typename TConfig, typename TId>
inline typename Value<typename FragmentStore<TSpec, TConfig>::TReadSeqStore>::Type
getRead(
	FragmentStore<TSpec, TConfig> &me, 
	TId id)
{
	return value(me.readSeqStore, id);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.appendMatePair
..summary:Appends two paired-end reads to a fragment store.
..cat:Fragment Store
..signature:appendMatePair(store, readId1, readId2)
..signature:appendMatePair(store, readId1, readId2, name1, name2)
..param.store:The fragment store.
...type:Class.FragmentStore
..param.readId1:The read sequence of the first read.
..param.readId2:The read sequence of the second read.
..param.name1:The read name of the first read.
..param.name2:The read name of the second read.
..returns:The $matePairId$ of the newly appended mate-pair.
..remarks:This function appends two reads to the @Memvar.FragmentStore#readStore@ and @Memvar.FragmentStore#readSeqStore@ 
and a mate-pair entry between both of them to the @Memvar.FragmentStore#matePairStore@.
If names are given, they are appended to the @Memvar.FragmentStore#readNameStore@.
..see:Function.appendRead
*/

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

/**
.Function.compactAlignedReads
..summary:Removes invalid aligned reads and rename $alignId$ sequentially beginning with 0.
..cat:Fragment Store
..signature:compactAlignedReads(store)
..param.store:The fragment store.
...type:Class.FragmentStore
..returns:The new size of the @Memvar.FragmentStore#alignedReadStore@.
..remarks:This function removes all entries from @Memvar.FragmentStore#alignedReadStore@ whose $alignId$ equals to $INVALID_ID$ as well as orphan entries
in @Memvar.FragmentStore#alignQualityStore@.
Afterwards the alignIds are renamed sequentially beginning with 0.
This function can be used to remove alignments which are selected by previously setting their id to $INVALID_ID$.
*/

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

/**
.Function.compactPairMatchIds
..summary:Renames $pairMatchId$ sequentially beginning with 0.
..cat:Fragment Store
..signature:compactPairMatchIds(store)
..param.store:The fragment store.
...type:Class.FragmentStore
..returns:The number of pair matches.
..remarks:This function renames the $pairMatchId$ in the @Memvar.FragmentStore#alignedReadStore@ sequentially beginning with 0.
Two read alignments can be identified to be a pair match if they have the same $pairMatchId$.
Please note that paired reads not necessarily have to mapped as a pair match, 
e.g. if they are on different contigs or have the same orientation or a wrong insert size.
*/

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

/**
.Function.calculateInsertSizes
..summary:Calculates a string with insert sizes for each pair match.
..cat:Fragment Store
..signature:compactPairMatchIds(insertSizes, store)
..param.insertSizes:The resulting string of insert sizes.
...remarks:This string is accordingly resized and can be addressed by the $pairMatchId$.
..param.store:The fragment store.
...type:Class.FragmentStore
..remarks:This function calls @Function.compactPairMatchIds@ first and calculate the insert size for every pair match.
The insert size of a pair match is the outer distance between the two matches.
*/

template <typename TLibSizeString, typename TSpec, typename TConfig>
inline void
calculateInsertSizes(TLibSizeString &insertSizes, FragmentStore<TSpec, TConfig> &me)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Id<TAlignedRead>::Type							TId;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignedReadIter;
	typedef typename TFragmentStore::TContigPos						TGPos;

	TAlignedReadIter it = begin(me.alignedReadStore, Standard());
	TAlignedReadIter itEnd = end(me.alignedReadStore, Standard());

	resize(insertSizes, compactPairMatchIds(me), Exact());
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
			insertSizes[id] = (*it).beginPos - leftMatePos;
	}
}

/**
.Function.getMateNo
..summary:Returns the mate number of read for a given $readId$.
..cat:Fragment Store
..signature:getMateNo(store, readId)
..param.store:The fragment store.
...type:Class.FragmentStore
..param.readId:The read id.
..returns:The mate number (0..first mate, 1..second mate) of the read in its mate-pair or -1 if the read is not paired.
*/

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


/**
.Function.calculateMateIndices
..summary:Calculates a string that maps the $readId$ of a read to the $readId$ of its mate.
..cat:Fragment Store
..signature:calculateMateIndices(mateIndices, store)
..param.mateIndices:The resulting string of mate indices.
...remarks:This string is accordingly resized and can be addressed by the $readId$.
..param.store:The fragment store.
...type:Class.FragmentStore
..remarks:Entries of reads without a mate contain $INVALID_ID$.
*/

// calculate index of the other mate for each pair match
template <typename TMateIndexString, typename TSpec, typename TConfig>
inline void
calculateMateIndices(TMateIndexString &mateIndices, FragmentStore<TSpec, TConfig> &me)
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
		if (length(mateIndices) < 2*id + 2)
			fill(mateIndices, 2*id + 2, TAlignedRead::INVALID_ID, Generous());
		mateIndices[2*id + 1 - getMateNo(me, (*it).readId)] = idx;
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.layoutAlignment
..summary:Calculates a visible layout of reads aligned to a contig.
..cat:Fragment Store
..signature:layoutAlignment(layout, store, contigId)
..param.layout:The resulting layout string. Must be a $String<String<unsigned> >$.
...remarks:This string is accordingly resized and can be addressed by the $readId$.
..param.store:The fragment store.
...type:Class.FragmentStore
..param.contigId:The $contigId$ of the affected contig.
..remarks:This function layouts all reads aligned to a contig in rows from up to down reusing empty row spaces.
$layout[row][gappedPos]$ stores the $alignId$ of the alignment beginning at position $gappedPos$ in gap-space.
$row$ is the row of the alignment in the multiple sequence alignment.
..see:Function.printAlignment
*/

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

/**
.Function.printAlignment
..summary:Prints a window of the visible layout of reads into a outstream.
..cat:Fragment Store
..signature:layoutAlignment(stream, layout, store, contigId, posBegin, posEnd, lineBegin, lineEnd)
..param.stream:A C++ outstream, e.g. std::cout.
..param.layout:A layout string created by a previous call of @Function.layoutAlignment@.
..param.store:The fragment store.
...type:Class.FragmentStore
..param.contigId:The $contigId$ of the affected contig.
..param.posBegin:Window begin position in gap-space.
..param.posEnd:Window end position in gap-space.
..param.lineBegin:Begin line of the window.
..param.lineEnd:End line of the window.
..remarks:The window coordinates ($beginPos$, ...) may be chosen bigger than the layout.
The empty space is then filled with whitespaces.
..see:Function.layoutAlignment
*/

template <typename TStream, typename TLayoutStringSet, typename TSpec, typename TConfig, typename TContigId, typename TPos, typename TNum>
void printAlignment(
	TStream &stream, 
	TLayoutStringSet &layout, FragmentStore<TSpec, TConfig> &store, 
	TContigId contigId,
	TPos posBegin, TPos posEnd,
	TNum lineBegin, TNum lineEnd)
//	unsigned lastRead)
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
		TLayoutStringIter mid = itEnd;

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
//			if (*mid >= lastRead) continue;
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

/**
.Function.convertMatchesToGlobalAlignment
..summary:Converts all matches to a multiple global alignment in gap-space.
..cat:Fragment Store
..signature:convertMatchesToGlobalAlignment(store, score)
..param.store:The fragment store.
...type:Class.FragmentStore
..param.score:A score object used by @Function.globalAlignment@ in this function.
..remarks:Before calling this function all $gaps$ structures in @Memvar.FragmentStore#alignedReadStore@ and @Memvar.FragmentStore#contigStore@ must be empty, i.e. there are no gaps in the alignments.
This function iterates over entries in the @Memvar.FragmentStore#alignedReadStore@ and semi-global aligns each read to its contig segments given by begin and end position.
Gaps introduced by these pair-wise alignments are then inserted to the affected contig and reads correspondingly.
..remarks:The invariant that positions in the @Memvar.FragmentStore#alignedReadStore@ are in gap-space holds before (there were no gaps in alignments) and after calling this functions.
*/

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
			while (firstOverlap != it && _max((*firstOverlap).beginPos, (*firstOverlap).endPos) <= cBegin)
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

/**
.Function.convertPairWiseToGlobalAlignment
..summary:Converts pair-wise alignments to a multiple global alignment.
..cat:Fragment Store
..signature:convertPairWiseToGlobalAlignment(store, pairwiseContigGaps)
..param.store:The fragment store.
...type:Class.FragmentStore
..param.pairwiseContigGaps:A string of anchored contig gaps for every pairwise alignment.
..remarks:Before calling this function the $gaps$ structures in the @Memvar.FragmentStore#contigStore@ must be empty, i.e. there are no gaps in the contig.
The pairwise alignment gaps of the reads are stored in the $gaps$ structure in the @Memvar.FragmentStore#alignedReadStore@, whereas the pairwise alignment gaps of the contig are stored in the $pairwiseContigGaps$ string.
..remarks:After calling this functions all positions in the @Memvar.FragmentStore#alignedReadStore@ are in gap-space.
*/

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
			while (firstOverlap != it && _max((*firstOverlap).beginPos, (*firstOverlap).endPos) <= cBegin)
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
