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

#ifndef SEQAN_HEADER_MEMORYMANAGER_H
#define SEQAN_HEADER_MEMORYMANAGER_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Class.MemoryManager:
..cat:Seed Handling
..summary:Saves and manages data using IDs.
..signature:MemoryManager<TValue, Block<size> >
..param.TValue:The value type, that is the type of the items/characters stored in the manager.
...remarks:Use @Metafunction.Value@ to get the value type for a given class.
..param.size:A positive integer that specifies the number of values in each allocated block.
...remarks: Size has to be a power of 2, e.g., 1024.
*/
struct _Free_Pointer;
typedef Tag<_Free_Pointer> const FreeMemoryPointer;	

template<typename TValue, typename TSPec, typename TFree> 
class MemoryManager;


/**
.Spec.MemoryManagerPointer:
..general:Class.MemoryManager
..cat:Seed Handling
..summary: Datasize of at least 4 bytes. Not suitable for many blocks.
..signature:MemoryManager<TValue, Block<size> >
..param.TValue:The value type, that is the type of the items/characters stored in the manager.
...remarks:Use @Metafunction.Value@ to get the value type for a given class.
..param.size:A positive integer that specifies the number of values in each allocated block.
...remarks: Size has to be a power of 2, e.g., 1024.
*/
template<typename TValue, unsigned int SPACE, typename TFree>
class MemoryManager<TValue, Block<SPACE>, TFree> 
{
	typedef String<TValue, Array<SPACE> >				TBlock;
	typedef TBlock*										PBlock;
	typedef Allocator< SinglePool<sizeof(TBlock)> >		TAllocator;

public:

	typedef typename Iterator<TBlock, Standard>::Type	TBlockIter;
	typedef String<PBlock>								TBlockTable;

	TBlockTable		blocks;
	TBlockIter		blockFirst, blockLast;	// current block boundaries
	TBlockIter		lastValue;				// pointer to top value
	TAllocator		alloc;

	
    void *pEmptySpace;
	unsigned int exponent;
	//____________________________________________________________________________
	      
	public:
	MemoryManager():
		blockFirst(TBlockIter()),
		blockLast(TBlockIter()),
		lastValue(TBlockIter()) 
	{
		SEQAN_CHECKPOINT
		unsigned int x = SPACE;
		exponent = 0;
		while (x != 1){
			x >>=1;
			++exponent;
		}
		pEmptySpace = 0;
	}

	template<typename TSource>
	MemoryManager(TSource const& source):
		blockFirst(TBlockIter()),
		blockLast(TBlockIter()),
		lastValue(TBlockIter())
	{
		SEQAN_CHECKPOINT
		unsigned int x = SPACE;
		exponent = 0;
		while (x != 1){
			x >>=1;
			++exponent;
		}
		assign(*this, source);
	} 

	MemoryManager(MemoryManager const & source):
		blockFirst(TBlockIter()),
		blockLast(TBlockIter()),
		lastValue(TBlockIter())
	{
		SEQAN_CHECKPOINT 
		unsigned int x = SPACE;
		exponent = 0;
		while (x != 1){
			x >>=1;
			++exponent;
		}
		assign(*this, source);
	}

	~MemoryManager() 
	{
		clear(*this);
	}

	//____________________________________________________________________________



	public:
		template<typename TPos>
		inline typename Reference<MemoryManager>::Type 
			operator[] (TPos pos) 
		{
		SEQAN_CHECKPOINT	
			return value(*this, pos);
		}

		template<typename TPos>
		inline typename Reference<MemoryManager const>::Type 
			operator[] (TPos pos) const 
		{
		SEQAN_CHECKPOINT
			return value(*this, pos);
		}
	};


template<typename TValue, unsigned int SPACE, typename TFree>
struct DefaultOverflowImplicit< MemoryManager<TValue, Block<SPACE> , TFree> >
{
	typedef Generous Type;
};


///.Metafunction.Value.param.T.type:Class.MemoryManager
template <typename TValue, unsigned int SPACE, typename TFree>
struct Value<MemoryManager<TValue,Block<SPACE>, TFree > >
{
	typedef TValue Type;
};


template<typename TValue, unsigned int SPACE, typename TFree, typename TSource>
inline void 
assign(
	MemoryManager<TValue, Block<SPACE>, TFree>& target, 
	TSource const& source) 
{
	SEQAN_CHECKPOINT
	clear(target);
	_raiseMemory(target,length(source));
	for (unsigned int i = 0; i < length(source);++i){
		target[i] = source[i];
	}

	target.lastValue = &target[length(source)-1];

	target.pEmptySpace = 0;
	const void* pTmpSource = source.pEmptySpace;
	void* pTmpTarget = &target.pEmptySpace;
	while (pTmpSource != 0){
		typename Size<String<TValue, Block<SPACE> > >::Type pos = _getPosition(source,pTmpSource);
		*((void**)pTmpTarget) = &target[pos];
		pTmpSource = (TValue *)source[pos];
		pTmpTarget = &target[pos];
	}
}

template<typename TValue, unsigned int SPACE, typename TFree, typename TPos>
inline typename Reference<String<TValue, Block<SPACE> > >::Type 
value(
	MemoryManager<TValue, Block<SPACE>, TFree >& stack, 
	TPos const pos) 
{
	SEQAN_CHECKPOINT
	return value(*(stack.blocks[pos>>stack.exponent]), pos % SPACE);
}


template<typename TValue, unsigned int SPACE, typename TFree, typename TPos>
inline typename Reference<String<TValue, Block<SPACE> > >::Type 
value(MemoryManager<TValue, Block<SPACE>, TFree > const& stack, 
	  TPos const pos) 
{
	SEQAN_CHECKPOINT
	return value(*(stack.blocks[pos>>stack.exponent]), pos % SPACE);
}

template<typename TValue, unsigned int SPACE, typename TFree>
inline void 
clear(MemoryManager<TValue, Block<SPACE>, TFree >& me)
{
SEQAN_CHECKPOINT
	typedef String<TValue, Block<SPACE>	>			TBlockString;
	typedef typename TBlockString::TBlockTable		TBlockTable;
	typedef typename Iterator<TBlockTable, Standard>::Type	TIter;
	
	TIter it = begin(me.blocks), itEnd = end(me.blocks);
	while (it != itEnd) {
		deallocate(me.alloc, *it, 1);
		++it;
	}
	me.pEmptySpace = 0;
	clear(me.blocks);
	me.lastValue = me.blockLast = typename TBlockString::TBlockIter();
}


template<typename TValue, unsigned int SPACE, typename TFree>
inline typename Size<String<TValue, Block<SPACE> > >::Type
length(MemoryManager<TValue, Block<SPACE>, TFree > const & me) 
{
	SEQAN_CHECKPOINT
	if (length(me.blocks))
		return (length(me.blocks) - 1) * SPACE + (me.lastValue - me.blockFirst) + 1;
	else
		return 0;
}


template<typename TValue, unsigned int SPACE, typename TFree>
inline typename Size<String<TValue, Block<SPACE> > >::Type
capacity(MemoryManager<TValue, Block<SPACE>, TFree > const & me) 
{
SEQAN_CHECKPOINT
	if (length(me.blocks))
		return length(me.blocks) * SPACE;
	else
		return 0;
}


template<typename TValue, unsigned int SPACE, typename TFree> 
void
_raiseMemory(MemoryManager<TValue,Block<SPACE>, TFree > &manager, 
			 typename Size<String<TValue, Block<SPACE> > >::Type leng)
{
	SEQAN_CHECKPOINT
	while (capacity(manager) < leng){
		typename Size< String<TValue, Block<SPACE> > >::Type last = length(manager.blocks);
		resize(manager.blocks, last + 1, typename DefaultOverflowImplicit<String<TValue, Block<SPACE> > >::Type());
		allocate(manager.alloc, manager.blocks[last], 1);
		manager.lastValue = manager.blockFirst = begin(*manager.blocks[last]);
		manager.blockLast = (manager.blockFirst + (SPACE - 1));
	}
}


/**
.Function.releaseID:
..summary: Releases a given id so it can be redistributed later on.
..cat:Seed Handling
..signature:releaseID(manager,pos)
..param.manager: The MemoryManager
...type:Class.MemoryManager
..param.pos: The ID that should be released.
..remarks: The value at this position is not deleted.
*/
template<typename TValue, unsigned int SPACE, typename TFree, typename TPosition> 
void
releaseID(MemoryManager<TValue,Block<SPACE>, TFree > &manager, 
		  TPosition pos)
{
	SEQAN_CHECKPOINT
	void* pTmp = &manager[pos];
	*((void**)pTmp) = manager.pEmptySpace;
	manager.pEmptySpace = &manager[pos];
}


/**
.Function.obtainID:
..summary: btains a new id from the id manager.
..cat:Seed Handling
..signature:obtainID(manager)
..param.manager: The MemoryManager
...type:Class.MemoryManager
..return: A unused ID.
*/
template<typename TValue, unsigned int SPACE, typename TFree> 
inline typename Size<String<TValue, Block<SPACE> > >::Type
obtainID(MemoryManager<TValue,Block<SPACE>,TFree > &manager)
{
	SEQAN_CHECKPOINT
	if (manager.pEmptySpace == 0){
		if (length(manager) == capacity(manager)){
			typename Size< String<TValue, Block<SPACE> > >::Type last = length(manager.blocks);
			resize(manager.blocks, last + 1, typename DefaultOverflowImplicit<String<TValue, Block<SPACE> > >::Type());
			allocate(manager.alloc, manager.blocks[last], 1);
			manager.lastValue = manager.blockFirst = begin(*manager.blocks[last]);
			manager.blockLast = (manager.blockFirst + (SPACE - 1));
		/*	for (int i = 0; i < manager.maps.size(); ++i){
				raiseMemory(manager.maps[i]);
			}*/
		} else
			++manager.lastValue;
		return length(manager)-1;
		
	}else{
		int i = 0;
		long tmp = (TValue*)manager.pEmptySpace - (TValue*)manager.blocks[i];
		while ((tmp >=SPACE)||( tmp < 0)){
			++i;
			tmp = (TValue*)manager.pEmptySpace - (TValue*)manager.blocks[i];
		}
		manager.pEmptySpace = *(void**)manager.pEmptySpace;
		return i * SPACE + tmp;
	}
}

template<typename TValue, unsigned int SPACE, typename TFree> 
inline typename Size<String<TValue, Block<SPACE> > >::Type
_getPosition(MemoryManager<TValue,Block<SPACE>, TFree > const &manager,
			 const void* pointer)
{
	SEQAN_CHECKPOINT
	int i = 0;
	long tmp = (TValue*)pointer - (TValue*)manager.blocks[i];
	while ((tmp >=SPACE)||( tmp < 0)){
		++i;
		tmp = (TValue*)pointer - (TValue*)manager.blocks[i];
	}
	return i * SPACE + tmp;
}


template<typename TValue, unsigned int SPACE, typename TFree, typename TPosition> 
void
assignValue(MemoryManager<TValue,Block<SPACE>, TFree > &manager, 
			TPosition pos, 
			TValue value)
{
	SEQAN_CHECKPOINT
	if (manager.pEmptySpace!=0){
		void* pTmp = &manager[pos];
		manager.pEmptySpace =*((void**)pTmp);
	} 
	manager[pos] = value;
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
