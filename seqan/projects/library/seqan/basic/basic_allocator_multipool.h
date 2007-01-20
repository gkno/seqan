#ifndef SEQAN_HEADER_BASIC_ALLOCATOR_MULTIPOOL_H
#define SEQAN_HEADER_BASIC_ALLOCATOR_MULTIPOOL_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// MultiPool Allocator
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Multi Pool Allocator:
..cat:Allocators
..general:Class.Allocator
..summary:Allocator that pools memory blocks.
..signature:Allocator< MultiPool<ParentAllocator> >
..param.ParentAllocator:An allocator that is by the pool allocator used to allocate memory.
...default:@Spec.Simple Allocator@
...note:The multi pool allocator only supports @Function.clear@ if this function is also implemented for $ParentAllocator$.
..remarks:A pool allocator allocates several memory blocks at once. 
Freed blocks are not immediately deallocated but recycled in subsequential allocations.
This way, the number of calls to the heap allocator is reduced and that speeds up memory management.
...text:Note that memory blocks larger than $Allocator<MultiPool< > >::BLOCKING_LIMIT$ are not pooled 
but immediately allocated and deallocated using $ParentAllocator$.
*/


template <typename TParentAllocator = SimpleAllocator>
struct MultiPool;

//////////////////////////////////////////////////////////////////////////////

typedef Allocator<MultiPool<> > PoolAllocator;

template <typename TParentAllocator>
struct Allocator<MultiPool<TParentAllocator> >
{
	enum
	{
		BLOCKING_LIMIT = 0x100,
		GRANULARITY_BITS = 2,
		BLOCKING_COUNT = BLOCKING_LIMIT >> GRANULARITY_BITS,
		STORAGE_SIZE = 0xf80
	};

	char * data_recycled_blocks [BLOCKING_COUNT];
	char * data_current_begin [BLOCKING_COUNT];
	char * data_current_free [BLOCKING_COUNT];
	Holder<TParentAllocator> data_parent_allocator;

	Allocator()
	{
SEQAN_CHECKPOINT
		::std::memset(data_recycled_blocks, 0, sizeof(data_recycled_blocks));
		::std::memset(data_current_begin, 0, sizeof(data_current_begin));
		::std::memset(data_current_free, 0, sizeof(data_current_free));
	}

	//Dummy copy
	Allocator(Allocator const &)
	{
		::std::memset(data_recycled_blocks, 0, sizeof(data_recycled_blocks));
		::std::memset(data_current_begin, 0, sizeof(data_current_begin));
		::std::memset(data_current_free, 0, sizeof(data_current_free));
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

template <typename TParentAllocator>
inline TParentAllocator &
parentAllocator(Allocator<MultiPool<TParentAllocator> > & me)
{
	return value(me.data_parent_allocator);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TParentAllocator>
inline void
setParentAllocator(Allocator<MultiPool<TParentAllocator> > & me,
				   TParentAllocator & alloc_)
{
	setValue(me.data_parent_allocator, alloc_);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TParentAllocator>
void
clear(Allocator<MultiPool<TParentAllocator> > & me)
{
SEQAN_CHECKPOINT
	::std::memset(me.data_recycled_blocks, 0, sizeof(me.data_recycled_blocks));
	::std::memset(me.data_current_begin, 0, sizeof(me.data_current_begin));
	::std::memset(me.data_current_free, 0, sizeof(me.data_current_free));

	clear(parentAllocator(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TParentAllocator>
inline unsigned int
_allocatorBlockNumber(Allocator<MultiPool<TParentAllocator> > &,
					  size_t size_)
{
SEQAN_CHECKPOINT
	typedef Allocator<MultiPool<TParentAllocator> > TAllocator;

	SEQAN_ASSERT(size_)

	if (size_ < TAllocator::BLOCKING_LIMIT)
	{//blocks
		return size_ >> TAllocator::GRANULARITY_BITS;
	}
	else
	{//no blocking
		return TAllocator::BLOCKING_COUNT;
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TParentAllocator, typename TValue, typename TSize, typename TUsage>
inline void
allocate(Allocator<MultiPool<TParentAllocator> > & me, 
		 TValue * & data,
		 TSize count,
		 Tag<TUsage> const tag_)
{
SEQAN_CHECKPOINT
	typedef Allocator<MultiPool<TParentAllocator> > TAllocator;

	size_t bytes_needed = count * sizeof(TValue);
	char * ptr;

	unsigned int block_number =  _allocatorBlockNumber(me, bytes_needed);
	if (block_number == TAllocator::BLOCKING_COUNT)
	{//no blocking
		return allocate(parentAllocator(me), data, count, tag_);
	}

	bytes_needed = (block_number + 1) << TAllocator::GRANULARITY_BITS;

	if (me.data_recycled_blocks[block_number])
	{//use recycled
		ptr = me.data_recycled_blocks[block_number];
		me.data_recycled_blocks[block_number] = * reinterpret_cast<char **>(ptr);
	}
	else
	{//use new
		ptr = me.data_current_free[block_number];
		if (!ptr || (ptr + bytes_needed > me.data_current_begin[block_number] + TAllocator::STORAGE_SIZE))
		{//not enough free space in current storage: allocate new
			allocate(parentAllocator(me), ptr, (size_t) TAllocator::STORAGE_SIZE, tag_);
			me.data_current_begin[block_number] = ptr;
		}
		me.data_current_free[block_number] = ptr + bytes_needed;
	}

	data = reinterpret_cast<TValue *>(ptr);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TParentAllocator, typename TValue, typename TSize, typename TUsage>
inline void 
deallocate(Allocator<MultiPool<TParentAllocator> > & me,
		   TValue * data, 
		   TSize count,
		   Tag<TUsage> const tag_)
{
SEQAN_CHECKPOINT
	typedef Allocator<MultiPool<TParentAllocator> > TAllocator;

	size_t bytes_needed = count * sizeof(TValue);

	unsigned int block_number = _allocatorBlockNumber(me, bytes_needed);
	if (block_number == TAllocator::BLOCKING_COUNT)
	{//no blocking
		return deallocate(parentAllocator(me), data, count, tag_);
	}

	bytes_needed = (block_number + 1) << TAllocator::GRANULARITY_BITS;

	//link in recycling list
	*reinterpret_cast<char **>(data) = me.data_recycled_blocks[block_number];
	me.data_recycled_blocks[block_number] = reinterpret_cast<char *>(data);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
