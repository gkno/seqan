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

#ifndef SEQAN_HEADER_MISC_DEQUEUE_H
#define SEQAN_HEADER_MISC_DEQUEUE_H

#include <algorithm>
#include <seqan/sequence.h>

//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

/**
.Class.Dequeue:
..cat:Miscellaneous
..summary:A double-ended queue implementation on top of a @Class.String@.
..signature:Dequeue<TValue, TSpec>
..param.TValue:Type of the ungapped sequences.
...metafunction:Metafunction.Value
..param.TSpec:The specializing type of the underlying @Class.String@.
...metafunction:Metafunction.Spec
...default:$Alloc<>$, see @Spec.Alloc String@
..include:seqan/misc.h
*/
template <typename TValue, typename TSpec = Alloc<> >
class Dequeue
{
public:
	typedef String<TValue, TSpec>						TString;
	typedef typename Iterator<TString, Standard>::Type	TIter;

	String<TValue, TSpec> data_string;

	TIter data_begin;	// string beginning
	TIter data_end;		// string end

	TIter data_front;	// front fifo character
	TIter data_back;	// back fifo character
	bool data_empty;	// fifo is empty

//____________________________________________________________________________

public:
	inline Dequeue()
	{
		clear(*this);
	}
};

//////////////////////////////////////////////////////////////////////////////
// Iterators
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Iterator.param.T.type:Class.Dequeue

template<typename TValue, typename TSpec>
struct Iterator<Dequeue<TValue, TSpec>, Standard> 
{
	typedef Iter<Dequeue<TValue, TSpec>, PositionIterator> Type;
};

template<typename TValue, typename TSpec>
struct Iterator<Dequeue<TValue, TSpec> const, Standard> 
{
	typedef Iter<Dequeue<TValue, TSpec> const, PositionIterator> Type;
};

template<typename TValue, typename TSpec>
struct Iterator<Dequeue<TValue, TSpec>, Rooted> 
{
	typedef Iter<Dequeue<TValue, TSpec>, PositionIterator> Type;
};

template<typename TValue, typename TSpec>
struct Iterator<Dequeue<TValue, TSpec> const, Rooted> 
{
	typedef Iter<Dequeue<TValue, TSpec> const, PositionIterator> Type;
};


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline bool
empty(Dequeue<TValue, TSpec> const &me)
{
	return me.data_empty;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline void
clear(Dequeue<TValue, TSpec> &me)
{
	clear(me.data_string);
	me.data_begin = begin(me.data_string, Standard());
	me.data_end = end(me.data_string, Standard());

	me.data_front = me.data_back = me.data_begin;
	me.data_empty = true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TPos>
inline TValue &
value(Dequeue<TValue, TSpec> &me, TPos pos)
{
	typedef typename Size<Dequeue<TValue, TSpec> >::Type TSize;
	TSize wrap = length(me.data_string) - (me.data_front - me.data_begin);
	
	if ((TSize)pos < wrap)
		return value(me.data_front + pos);
	else
		return value(me.data_begin + (pos - wrap));
}

template <typename TValue, typename TSpec, typename TPos>
inline TValue const &
value(Dequeue<TValue, TSpec> const &me, TPos pos)
{
	typedef typename Size<Dequeue<TValue, TSpec> >::Type TSize;
	TSize wrap = length(me.data_string) - (me.data_front - me.data_begin);
	
	if ((TSize)pos < wrap)
		return value(me.data_front + pos);
	else
		return value(me.data_begin + (pos - wrap));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline TValue &
front(Dequeue<TValue, TSpec> &me)
{
	return *me.data_front;
}

template <typename TValue, typename TSpec>
inline TValue const &
front(Dequeue<TValue, TSpec> const &me)
{
	return *me.data_front;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline TValue &
back(Dequeue<TValue, TSpec> &me)
{
	return *me.data_back;
}

template <typename TValue, typename TSpec>
inline TValue const &
back(Dequeue<TValue, TSpec> const &me)
{
	return *me.data_back;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline bool
popFront(Dequeue<TValue, TSpec> &me)
{
	if (me.data_empty) return false;

	if (me.data_front == me.data_back)
		me.data_empty = true;
	else
	{
		if (++me.data_front == me.data_end)
			me.data_front = me.data_begin;
	}

	return true;
}

template <typename TValue, typename TSpec>
inline bool
popBack(Dequeue<TValue, TSpec> &me)
{
	if (me.data_empty) return false;

	if (me.data_front == me.data_back)
		me.data_empty = true;
	else
	{
		if (me.data_back == me.data_begin)
			me.data_back = me.data_end;
		--me.data_back;
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline void
pushFront(Dequeue<TValue, TSpec> &me, TValue const & _value)
{
	typedef typename Dequeue<TValue, TSpec>::TIter TIter;

	if (me.data_empty) 
	{
		if (me.data_begin == me.data_end)
			reserve(me, computeGenerousCapacity(me, length(me.data_string) + 1), Generous());
		me.data_empty = false;
	}
	else 
	{
		TIter new_front = me.data_front;
		if (new_front == me.data_begin)
			new_front = me.data_end;
		--new_front;

		if (new_front == me.data_back)
		{
			reserve(me, computeGenerousCapacity(me, length(me.data_string) + 1), Generous());

			if (me.data_front == me.data_begin)
				me.data_front = me.data_end;
			--me.data_front;
		} else
			me.data_front = new_front;
	}
	assign(*me.data_front, _value);
}

template <typename TValue, typename TSpec>
inline void
pushBack(Dequeue<TValue, TSpec> &me, TValue const & _value)
{
	typedef typename Dequeue<TValue, TSpec>::TIter TIter;

	if (me.data_empty) 
	{
		if (me.data_begin == me.data_end)
			reserve(me, computeGenerousCapacity(me, length(me.data_string) + 1), Generous());
		me.data_empty = false;
	}
	else 
	{
		TIter new_back = me.data_back;
		if (++new_back == me.data_end)
			new_back = me.data_begin;

		if (new_back == me.data_front)
		{
			reserve(me, computeGenerousCapacity(me, length(me.data_string) + 1), Generous());
			// in this case reserve adds new space behind data_back
			++me.data_back;
		} else
			me.data_back = new_back;
	}
	assign(*me.data_back, _value);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline typename Size<Dequeue<TValue, TSpec> >::Type
length(Dequeue<TValue, TSpec> const &me)
{
	if (empty(me)) return 0;

	if (me.data_front <= me.data_back)
		return (me.data_back - me.data_front) + 1;
	else
		return (me.data_end - me.data_begin) - (me.data_front - me.data_back) + 1;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TSize_, typename TExpand>
inline typename Size<Dequeue<TValue, TSpec> >::Type
reserve(Dequeue<TValue, TSpec> &me, TSize_ new_capacity, Tag<TExpand> const tag)
{
	typedef typename Size<Dequeue<TValue, TSpec> >::Type TSize;
//	::std::cout << "resize to "<<new_capacity<<::std::endl;
	TSize len = length(me);
	if (len < new_capacity && length(me.data_string) != new_capacity)
	{
		TSize pos_front = me.data_front - me.data_begin;
		TSize pos_back  = me.data_back  - me.data_begin;
		TSize new_freeSpace = new_capacity - len;

		if (pos_front <= pos_back)
		{
			// |empty|data|empty|
			// 0
			TSize freeSpace = length(me.data_string) - len;
			if (new_freeSpace > freeSpace)
				resize(me.data_string, new_capacity, tag);
			else
			{
				freeSpace -= new_freeSpace;	// reduce the free space by <freeSpace>
				if (pos_front >= freeSpace)
				{
					resizeSpace(me.data_string, pos_front - freeSpace, (TSize)0, pos_front, tag);
					pos_back -= freeSpace;
					pos_front -= freeSpace;
				}
				else
				{
					freeSpace -= pos_front;
					resizeSpace(me.data_string, length(me.data_string) - freeSpace, pos_back + 1, length(me.data_string), tag);
					resizeSpace(me.data_string, (TSize)0, (TSize)0, pos_front, tag);
					pos_back -= pos_front;
					pos_front = 0;
				}
			}
		}
		else
		{
			// |data|empty|data|
			// 0
			resizeSpace(me.data_string, new_freeSpace, pos_back + 1, pos_front, tag);
			pos_front += new_freeSpace;
		}

		me.data_begin = begin(me.data_string, Standard());
		me.data_end = end(me.data_string, Standard());
		me.data_front = me.data_begin + pos_front;
		me.data_back = me.data_begin + pos_back;
	}
	return length(me.data_string);
}

}

#endif

