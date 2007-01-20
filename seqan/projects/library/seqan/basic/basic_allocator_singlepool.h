#ifndef SEQAN_HEADER_BASIC_ALLOCATOR_SINGLE_POOL_H
#define SEQAN_HEADER_BASIC_ALLOCATOR_SINGLE_POOL_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// SinglePool Allocator
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Single Pool Allocator:
..cat:Allocators
..general:Class.Allocator
..summary:Allocator that pools memory blocks of specific size.
..signature:Allocator< SinglePool<SIZE, ParentAllocator> >
..param.SIZE:Size of memory blocks that are pooled.
...value:An unsigned integer with $SIZE >= sizeof(void *)$.
..param.ParentAllocator:An allocator that is by the pool allocator used to allocate memory.
...default:@Spec.Simple Allocator@
...note:The single pool allocator only supports @Function.clear@ if this function is also implemented for $ParentAllocator$.
..remarks:A pool allocator allocates several memory blocks at once. 
Freed blocks are not immediately deallocated but recycled in subsequential allocations.
This way, the number of calls to the heap allocator is reduced and that speeds up memory management.
...text:The single pool allocator only pools memory blocks of size $SIZE$.
Blocks of other sizes are allocated and deallocated using an allocator of type $ParentAllocator$.
...text:Using the single pool allocator for blocksizes larger than some KB is not advised.
*/

template <size_t SIZE, typename TParentAllocator = SimpleAllocator>
struct SinglePool;

//////////////////////////////////////////////////////////////////////////////

template <size_t SIZE, typename TParentAllocator>
struct Allocator<SinglePool<SIZE, TParentAllocator> >
{
	enum
	{
		SIZE_PER_ITEM = SIZE,
		ITEMS_PER_BLOCK = (SIZE_PER_ITEM < 0x0100) ? 0x01000 / SIZE_PER_ITEM : 16,
		STORAGE_SIZE = SIZE * ITEMS_PER_BLOCK
	};

	char * data_recycled_blocks;
	char * data_current_begin;
	char * data_current_free;
	Holder<TParentAllocator> data_parent_allocator;

	Allocator():
		data_recycled_blocks(0),
		data_current_begin(0),
		data_current_free(0)
	{
SEQAN_CHECKPOINT
	}

	//Dummy copy
	Allocator(Allocator const &):
		data_recycled_blocks(0),
		data_current_begin(0),
		data_current_free(0)
	{
	}
	inline Allocator &
	operator = (Allocator const &)
	{
		clear(*this);
		return *this;
	}

	~Allocator()
	{
SEQAN_CHECKPOINT
		clear(*this);
	}
};
//////////////////////////////////////////////////////////////////////////////

template <size_t SIZE, typename TParentAllocator>
inline TParentAllocator &
parentAllocator(Allocator<SinglePool<SIZE, TParentAllocator> > & me)
{
	return value(me.data_parent_allocator);
}

//////////////////////////////////////////////////////////////////////////////

template <size_t SIZE, typename TParentAllocator>
inline void
setParentAllocator(Allocator<SinglePool<SIZE, TParentAllocator> > & me,
				   TParentAllocator & alloc_)
{
	setValue(me.data_parent_allocator, alloc_);
}

//////////////////////////////////////////////////////////////////////////////

template <size_t SIZE, typename TParentAllocator>
void
clear(Allocator<SinglePool<SIZE, TParentAllocator> > & me)
{
SEQAN_CHECKPOINT

	me.data_recycled_blocks = 0;
	me.data_current_begin = 0;
	me.data_current_free = 0;

	clear(parentAllocator(me));
}

//////////////////////////////////////////////////////////////////////////////

template <size_t SIZE, typename TParentAllocator, typename TValue, typename TSize, typename TUsage>
inline void
allocate(Allocator<SinglePool<SIZE, TParentAllocator> > & me, 
		 TValue * & data,
		 TSize count,
		 Tag<TUsage> const tag_)
{
SEQAN_CHECKPOINT
	typedef Allocator<SinglePool<SIZE, TParentAllocator> > TAllocator;
	size_t bytes_needed = count * sizeof(TValue);

	if (bytes_needed != TAllocator::SIZE_PER_ITEM)
	{//no blocking
		allocate(parentAllocator(me), data, count, tag_);
		return;
	}

	char * ptr;
	if (me.data_recycled_blocks)
	{//use recycled
		ptr = me.data_recycled_blocks;
		me.data_recycled_blocks = * reinterpret_cast<char **>(ptr);
	}
	else
	{//use new
		ptr = me.data_current_free;
		if (!ptr || (ptr + bytes_needed > me.data_current_begin + TAllocator::STORAGE_SIZE))
		{//not enough free space in current storage: allocate new
			allocate(parentAllocator(me), ptr, (size_t) TAllocator::STORAGE_SIZE, tag_);
			me.data_current_begin = ptr;
		}
		me.data_current_free = ptr + bytes_needed;
	}

	data = reinterpret_cast<TValue *>(ptr);
}

//////////////////////////////////////////////////////////////////////////////

template <size_t SIZE, typename TParentAllocator, typename TValue, typename TSize, typename TUsage>
inline void 
deallocate(Allocator<SinglePool<SIZE, TParentAllocator> > & me,
		   TValue * data, 
		   TSize count,
		   Tag<TUsage> const tag_)
{
SEQAN_CHECKPOINT
	typedef Allocator<SinglePool<SIZE, TParentAllocator> > TAllocator;

	size_t bytes_needed = count * sizeof(TValue);

	if (bytes_needed != TAllocator::SIZE_PER_ITEM)
	{//no blocking
		deallocate(parentAllocator(me), data, count, tag_);
		return;
	}

	//link in recycling list
	*reinterpret_cast<char **>(data) = me.data_recycled_blocks;
	me.data_recycled_blocks = reinterpret_cast<char *>(data);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
