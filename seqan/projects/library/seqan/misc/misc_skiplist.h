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
  $Id  $
 ==========================================================================*/

#ifndef SEQAN_HEADER_MISC_SKIPLIST_H
#define SEQAN_HEADER_MISC_SKIPLIST_H

//////////////////////////////////////////////////////////////////////////////

#include <seqan/misc/misc_random.h>

//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{



//////////////////////////////////////////////////////////////////////////////
// Skiplist
//////////////////////////////////////////////////////////////////////////////

//forwards

template <typename TKey, typename TValue, typename TSpec = Default>
class Skiplist;

template <typename TKey, typename TValue, typename TSpec>
class SkiplistElement;

template <typename TKey, typename TValue, typename TSpec>
class SkiplistNext;

template <typename TKey, typename TValue, typename TSpec>
class SkiplistPath;

//////////////////////////////////////////////////////////////////////////////
// Tags
struct SkiplistIterator;

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct AllocatorType;

template <typename TKey, typename TValue, typename TSpec>
struct AllocatorType<Skiplist<TKey, TValue, TSpec> >
{
	typedef Allocator<SimpleAlloc<> > Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec>
struct Value<Skiplist<TKey, TValue, TSpec> >
{
	typedef TValue Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec>
struct Key<Skiplist<TKey, TValue, TSpec> >
{
	typedef TKey Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec>
struct Spec<Skiplist<TKey, TValue, TSpec> >
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec, typename TIteratorSpec>
struct Iterator<Skiplist<TKey, TValue, TSpec>, TIteratorSpec >
{
	typedef Skiplist<TKey, TValue, TSpec> TSkiplist;
	typedef Iter<TSkiplist, SkiplistIterator> Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec>
class Skiplist
{
public:
	typedef typename AllocatorType<Skiplist>::Type TAllocator;
	typedef SkiplistElement<TKey, TValue, TSpec> TElement;
	typedef typename Size<Skiplist>::Type TSize;
	typedef typename Value<Skiplist>::Type TValue2;

	enum
	{
		MAX_HEIGHT = 28,
		BLOCK_SIZE = 0x200
	};

	Holder<TAllocator> data_allocator;
	TElement * data_recycle[MAX_HEIGHT];
	unsigned char * data_mem_begin;
	unsigned char * data_mem_end;

	TElement data_border;
	TSize data_length;
	unsigned char data_height;

	Skiplist()
		: data_mem_begin(0)
		, data_mem_end(0)
		, data_length(0)
		, data_height(0)
	{
		for (unsigned char i = 0; i < MAX_HEIGHT; ++i)
		{
			data_recycle[i] = 0;
			valueConstruct(data_border.data_next + i, NonMinimalCtor()); 
		}

		mtRandInit();
	}
	~Skiplist()
	{
	}

	template <typename TKey2>
	inline TValue2 &
	operator [] (TKey2 const & key)
	{
		return value(*this, key);
	}

private:
	Skiplist(Skiplist const &);
	Skiplist const & operator = (Skiplist const &);
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec>
class SkiplistElement
{
public:
	typedef Skiplist<TKey, TValue, TSpec> TSkiplist;
	typedef SkiplistNext<TKey, TValue, TSpec> TNext;

	enum
	{
		MAX_HEIGHT = TSkiplist::MAX_HEIGHT
	};

	struct Data
	{
		TKey key;
		TValue value;
	} data;
	TNext data_next[MAX_HEIGHT]; //note: only parts of this array available

	//indirect constructor in _skiplistConstructElement
	//indirect destructor in _skiplistDestructElement
};

//spec for skiplists without value
template <typename TKey, typename TSpec>
class SkiplistElement<TKey, void, TSpec>
{
public:
	typedef Skiplist<TKey, void, TSpec> TSkiplist;
	typedef SkiplistNext<TKey, void, TSpec> TNext;

	enum
	{
		MAX_HEIGHT = TSkiplist::MAX_HEIGHT
	};

	struct Data
	{
		TKey key;
	} data;
	TNext data_next[MAX_HEIGHT]; //note: only parts of this array available

	//indirect constructor in _skiplistConstructElement
	//indirect destructor in _skiplistDestructElement
};


//////////////////////////////////////////////////////////////////////////////

//representation of "horizontal" pointer in skiplist
//can be overloaded if label is needed
template <typename TKey, typename TValue, typename TSpec>
class SkiplistNext
{
public:
	typedef SkiplistElement<TKey, TValue, TSpec> TElement;

	TElement * data_element;

	SkiplistNext()
	{}
	SkiplistNext(NonMinimalCtor)
		: data_element(0)
	{}
	SkiplistNext(SkiplistNext const & other)
		: data_element(other.data_element)
	{}
	~SkiplistNext()
	{}
	SkiplistNext const & operator = (SkiplistNext const & other)
	{
		data_element = other.data_element;
		return *this;
	}
};

//////////////////////////////////////////////////////////////////////////////

//represents a path from the root to the predecessor of a skiplist element
template <typename TKey, typename TValue, typename TSpec>
class SkiplistPath
{
public:
	typedef Skiplist<TKey, TValue, TSpec> TSkiplist;
	typedef SkiplistElement<TKey, TValue, TSpec> TElement;

	enum
	{
		MAX_HEIGHT = TSkiplist::MAX_HEIGHT
	};

	TElement * data_elements[MAX_HEIGHT];
//	unsigned char data_height; //is stored in Skiplist
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec, typename TSize, typename TUsage>
void allocate(Skiplist<TKey, TValue, TSpec> & me,
			  unsigned char * & buf,
			  TSize count,
			  Tag<TUsage> const tag)
{
	allocate(value(me.data_allocator), buf, count, tag);
}

//////////////////////////////////////////////////////////////////////////////

//create Space for SkiplistElement of given height

template <typename TKey, typename TValue, typename TSpec>
inline SkiplistElement<TKey, TValue, TSpec> &
_skiplistAllocateElement(Skiplist<TKey, TValue, TSpec> & me,
						 unsigned char height)
{
	typedef Skiplist<TKey, TValue, TSpec> TSkiplist;
	typedef SkiplistElement<TKey, TValue, TSpec> TElement;
	typedef SkiplistNext<TKey, TValue, TSpec> TNext;

	TElement * ret;

	if (me.data_recycle[height])
	{//use recycled
		ret = me.data_recycle[height];
		me.data_recycle[height] = * reinterpret_cast<TElement **>(me.data_recycle[height]);
	}
	else
	{
		int need_size = sizeof(typename TElement::Data) + (height+1) * sizeof(TNext);
		int buf_size = me.data_mem_end - me.data_mem_begin;
		if (buf_size < need_size)
		{//need new memory
			if (buf_size >= (sizeof(typename TElement::Data) + sizeof(TNext)))
			{//link rest memory in recycle 
				int rest_height = (buf_size - sizeof(typename TElement::Data)) / sizeof(TNext) - 1; //must be < height, because buf_size < need_size
				* reinterpret_cast<TElement **>(me.data_mem_begin) = me.data_recycle[rest_height];
				me.data_recycle[rest_height] = reinterpret_cast<TElement *>(me.data_mem_begin);
			}
			allocate(me, me.data_mem_begin, TSkiplist::BLOCK_SIZE, TagAllocateStorage());
			me.data_mem_end = me.data_mem_begin + TSkiplist::BLOCK_SIZE;
		}
		ret = reinterpret_cast<TElement *>(me.data_mem_begin);
		me.data_mem_begin += need_size;
	}

	return *ret;
}

//////////////////////////////////////////////////////////////////////////////

//Creates New SkiplistElement
template <typename TKey, typename TValue, typename TSpec, typename TKey2, typename TValue2>
inline SkiplistElement<TKey, TValue, TSpec> &
_skiplistConstructElement(Skiplist<TKey, TValue, TSpec> & me,
						  unsigned char height,
						  TKey2 const & key,
						  TValue2 const & value)
{
	typedef SkiplistElement<TKey, TValue, TSpec> TElement;
	TElement & el = _skiplistAllocateElement(me, height);
	valueConstruct(& (el.data.key), key);
	valueConstruct(& (el.data.value), value);
	//no need to construct the next array

	return el;
}

//variant for TValue == void
template <typename TKey, typename TSpec, typename TKey2>
inline SkiplistElement<TKey, void, TSpec> &
_skiplistConstructElement(Skiplist<TKey, void, TSpec> & me,
						  unsigned char height,
						  TKey2 const & key)
{
	typedef SkiplistElement<TKey, void, TSpec> TElement;
	TElement & el = _skiplistAllocateElement(me, height);
	valueConstruct(& (el.data.key), key);
	//no need to construct the next array

	return el;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec>
inline void
_skiplistDeallocateElement(Skiplist<TKey, TValue, TSpec> & me,
						   SkiplistElement<TKey, TValue, TSpec> & el,
						   unsigned char height)
{
	typedef SkiplistElement<TKey, TValue, TSpec> TElement;
	* reinterpret_cast<TElement *>(&el) = me.data_recycle[height];
	me.data_recycle[height] = reinterpret_cast<TElement *>(&el);
	//the real deallocation is done by the allocator on destruction 
}

//////////////////////////////////////////////////////////////////////////////

//Destroys New SkiplistElement
template <typename TKey, typename TValue, typename TSpec>
inline void
_skiplistDestructElement(Skiplist<TKey, TValue, TSpec> & me,
						 SkiplistElement<TKey, TValue, TSpec> & el,
						 unsigned char height)
{
	valueDestruct(& (el.data.key) );
	valueDestruct(& (el.data.value) );
	//no need to construct the next array
	_skiplistDeallocateElement(me, el, height);
}

//variant for TValue == void
template <typename TKey, typename TSpec>
inline void
_skiplistDestructElement(Skiplist<TKey, void, TSpec> & me,
						 SkiplistElement<TKey, void, TSpec> & el,
						 unsigned char height)
{
	valueDestruct(& (el.data.key) );
	//no need to construct the next array
	_skiplistDeallocateElement(me, el, height);
}

//////////////////////////////////////////////////////////////////////////////

//creates height for new SkiplistElement.
//increases the Skiplist height if necessary
template <typename TKey, typename TValue, typename TSpec>
inline unsigned char
_skiplistCreateHeight(Skiplist<TKey, TValue, TSpec> & me)
{
	typedef Skiplist<TKey, TValue, TSpec> TSkiplist;

	unsigned char height = geomRand<unsigned char>();
	if (height >= TSkiplist::MAX_HEIGHT) height = TSkiplist::MAX_HEIGHT-1;

	if (height > me.data_height) me.data_height = height;

	return height;
}

template <typename TKey, typename TValue, typename TSpec>
inline unsigned char
_skiplistCreateHeight(Skiplist<TKey, TValue, TSpec> & me,
					  SkiplistPath<TKey, TValue, TSpec> & path) //extend path if height is increased
{
	typedef Skiplist<TKey, TValue, TSpec> TSkiplist;

	unsigned char height = geomRand<unsigned char>();
	if (height >= TSkiplist::MAX_HEIGHT) height = TSkiplist::MAX_HEIGHT-1;

	if (height > me.data_height)
	{
		for (unsigned char i = me.data_height + 1; i <= height; ++i)
		{
			path.data_elements[i] = & me.data_border;
		}
		me.data_height = height;
	}

	return height;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec>
inline unsigned char
_skiplistGetHeight(Skiplist<TKey, TValue, TSpec> & me,
				   SkiplistElement<TKey, TValue, TSpec> & el,
				   SkiplistPath<TKey, TValue, TSpec> & path)
{
	int height = me.data_height;
	for (; height > 0 ; --height)
	{
		if (path.elements[height]->data_next[height] == el) break;
	}
	return height;
}

template <typename TKey, typename TValue, typename TSpec>
inline unsigned char
_skiplistGetHeight(Skiplist<TKey, TValue, TSpec> & me, 
				   SkiplistElement<TKey, TValue, TSpec> & el)
{
	typedef SkiplistPath<TKey, TValue, TSpec> TPath;

	TPath path;
	_skiplistFind(me, el, path);

	return _skiplistGetHeight(el, path);
}

//////////////////////////////////////////////////////////////////////////////

//Note: always store paths to the element LEFT of the actually wanted element.
//this element will be the predecessor during insertion.

//Given a key: find path to the elment that is left to:
// - the first element that has the given key, if there is any, or
// - the first element that has the smallest key larger than the given key,
//or a path to the last element, if all keys are smaller than the given key

//Given an element: find path to the element that is left to:
// - the given element, or
// - the first element that has the smallest key larger than the key of the given element,
//or a path to the last element, if all keys are smaller than the key of the given element

template <typename TKey, typename TValue, typename TSpec>
inline bool
_skiplistFindGoNext(SkiplistNext<TKey, TValue, TSpec> & next,
					unsigned char,
					TKey const & key)
{
	return (next.data_element->data.key < key);
}

template <typename TKey, typename TValue, typename TSpec>
inline bool
_skiplistFindGoNext(SkiplistNext<TKey, TValue, TSpec> & next,
					unsigned char,
					SkiplistElement<TKey, TValue, TSpec> const & el)
{
	return (next.data_element->data.key <= el.data.key) && (next.data_element != & el);
}

template <typename TKey, typename TValue, typename TSpec>
inline bool
_skiplistFindGoNext(SkiplistNext<TKey, TValue, TSpec> & next,
					unsigned char height,
					GoEnd)
{
	return next.data_element->data_next[height];
}


template <typename TKey, typename TValue, typename TSpec, typename TFind>
inline void
_skiplistFind(Skiplist<TKey, TValue, TSpec> & me,
			  TFind const & find, //can be a TKey or a SkiplistElement or GoEnd
			  /*OUT*/ SkiplistPath<TKey, TValue, TSpec> & path)
{
	typedef SkiplistElement<TKey, TValue, TSpec> TElement;
	typedef SkiplistNext<TKey, TValue, TSpec> TNext;

	TElement * here = & me.data_border;

	for (int i = me.data_height; i >= 0; --i)
	{
		while (true)
		{
			TNext & next = here->data_next[i];
			if (!next.data_element || !_skiplistFindGoNext(next, i, find)) break;
			here = next.data_element;
		}
		path.data_elements[i] = here;
	}
}



//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec, typename TFind>
inline typename Iterator< Skiplist<TKey, TValue, TSpec> >::Type
findElement(Skiplist<TKey, TValue, TSpec> & me,
			TFind const & find, //can be a TKey or a SkiplistElement or GoEnd
			SkiplistPath<TKey, TValue, TSpec> & path) 
{
	typedef typename Iterator< Skiplist<TKey, TValue, TSpec> >::Type TIterator;

	_skiplistFind(me, find, path);
	return TIterator(path.data_elements[0]->data_next[0].data_element);
}
template <typename TKey, typename TValue, typename TSpec, typename TFind>
inline typename Iterator< Skiplist<TKey, TValue, TSpec> >::Type
findElement(Skiplist<TKey, TValue, TSpec> & me,
			TFind const & find) //can be a TKey or a SkiplistElement or GoEnd
{
	typedef SkiplistPath<TKey, TValue, TSpec> TPath;
	TPath path;
	return findElement(me, find, path);
}

//////////////////////////////////////////////////////////////////////////////

//insert elements after at the position path points to
//Requirements: 
// - height <= height of the skiplist
// - next must be filled at least up to height
// - el.data_next must have space for at least height many entries
template <typename TKey, typename TValue, typename TSpec>
inline void
_skiplistInsertElement(SkiplistElement<TKey, TValue, TSpec> & el,
					   SkiplistPath<TKey, TValue, TSpec> & path,
					   unsigned char height)
{
	for (int i = height; i >= 0; --i)
	{
		el.data_next[i].data_element = path.data_elements[i]->data_next[i].data_element;
		path.data_elements[i]->data_next[i].data_element = & el;
	}
}

template <typename TKey, typename TValue, typename TSpec>
inline void
_skiplistInsertElement(Skiplist<TKey, TValue, TSpec> & me,
					   SkiplistElement<TKey, TValue, TSpec> & el,
					   unsigned char height)
{
	typedef SkiplistElement<TKey, TValue, TSpec> TElement;
	typedef SkiplistPath<TKey, TValue, TSpec> TPath;

	TPath path;
	_skiplistFind(me, el.data.key, path);

	////assert that the height of the first element that has key K is >= the heights of other elements that have the same key K
	//TElement & el2 = path.elements[0]->data_next[0];
	//if (el.data.key == el2.data.key)
	//{
	//	//determine height of el2
	//	int height2 = _skiplistGetHeight(me, el2, path);
	//	if (height2 > height)
	//	{//modify path in a way that el is inserted after el2
	//		for (; height2 > 0; --height2)
	//		{
	//			path.elements[height2] = el2;
	//		}
	//	}
	//}

	_skiplistInsertElement(el, path, height);

	++me.data_length;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec, typename TKey2, typename TValue2>
inline void
insertElement(Skiplist<TKey, TValue, TSpec> & me,
			  TKey2 const & key,
			  TValue2 const & value)
{
	typedef Skiplist<TKey, TValue, TSpec> TSkiplist;
	typedef SkiplistElement<TKey, TValue, TSpec> TElement;

	unsigned char height = _skiplistCreateHeight(me);
	TElement & el = _skiplistConstructElement(me, height, key, value);
	_skiplistInsertElement(me, el, height);
}
template <typename TKey, typename TSpec, typename TKey2>
inline void
insertElement(Skiplist<TKey, void, TSpec> & me,
			  TKey2 const & key)
{
	typedef Skiplist<TKey, void, TSpec> TSkiplist;
	typedef SkiplistElement<TKey, void, TSpec> TElement;

	unsigned char height = _skiplistCreateHeight(me);
	TElement & el = _skiplistConstructElement(me, height, key);
	_skiplistInsertElement(me, el, height);
}

//////////////////////////////////////////////////////////////////////////////

//extract element from Skiplist
template <typename TKey, typename TValue, typename TSpec>
inline void
_skiplistUnlinkElement(Skiplist<TKey, TValue, TSpec> & me,
					   SkiplistElement<TKey, TValue, TSpec> & el)
{
	typedef SkiplistPath<TKey, TValue, TSpec> TPath;

	TPath path;
	_skiplistFind(me, el, path);

	for (int i = me.data_height; i >= 0; --i)
	{
		if (path.data_elements[i]->data_next[i].data_element == & el)
		{
			path.data_elements[i]->data_next[i].data_element = el.data_next[i].data_element;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec, typename TIterator>
inline void
removeElement(Skiplist<TKey, TValue, TSpec> & me,
			  TIterator it)
{
	_skiplistUnlinkElement(me, * it.data_pointer);
	--me.data_length;
}


//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec>
inline typename Size< Skiplist<TKey, TValue, TSpec> >::Type
length(Skiplist<TKey, TValue, TSpec> const & me)
{
	return me.data_length;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec, typename TIteratorSpec>
inline typename Iterator< Skiplist<TKey, TValue, TSpec>, TIteratorSpec>::Type
begin(Skiplist<TKey, TValue, TSpec> & me,
	  TIteratorSpec)
{
	typedef typename Iterator< Skiplist<TKey, TValue, TSpec>, TIteratorSpec>::Type TIterator;
	return TIterator(me);
}
template <typename TKey, typename TValue, typename TSpec>
inline typename Iterator< Skiplist<TKey, TValue, TSpec> >::Type
begin(Skiplist<TKey, TValue, TSpec> & me)
{
	typedef typename Iterator< Skiplist<TKey, TValue, TSpec> >::Type TIterator;
	return TIterator(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec, typename TIteratorSpec>
inline typename Iterator< Skiplist<TKey, TValue, TSpec>, TIteratorSpec>::Type
end(Skiplist<TKey, TValue, TSpec> &,
	TIteratorSpec)
{
	typedef typename Iterator< Skiplist<TKey, TValue, TSpec>, TIteratorSpec>::Type TIterator;
	return TIterator();
}
template <typename TKey, typename TValue, typename TSpec>
inline typename Iterator< Skiplist<TKey, TValue, TSpec> >::Type
end(Skiplist<TKey, TValue, TSpec> &)
{
	typedef typename Iterator< Skiplist<TKey, TValue, TSpec> >::Type TIterator;
	return TIterator();
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec, typename TKey2>
inline typename Value< Skiplist<TKey, TValue, TSpec> >::Type &
value(Skiplist<TKey, TValue, TSpec> & me,
	  TKey2 const & _key)
{
	typedef Skiplist<TKey, TValue, TSpec> TSkiplist;
	typedef SkiplistPath<TKey, TValue, TSpec> TPath;
	typedef SkiplistElement<TKey, TValue, TSpec> TElement;
	typedef typename Iterator<TSkiplist>::Type TIterator;
	typedef typename Value<TSkiplist>::Type TValue2;

	TPath path;
	TIterator it = findElement(me, _key, path);
	if (it && (key(it) == _key))
	{
		return value(it);
	}
	else
	{// insert new value
		unsigned char height = _skiplistCreateHeight(me, path);
		TElement & el = _skiplistConstructElement(me, height, _key, TValue2());
		_skiplistInsertElement(el, path, height);
		++me.data_length;
		return el.data.value;
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec, typename TKey2>
inline bool
hasKey(Skiplist<TKey, TValue, TSpec> & me,
	   TKey2 const & _key)
{
	typedef Skiplist<TKey, TValue, TSpec> TSkiplist;
	typedef typename Iterator<TSkiplist>::Type TIterator;
	TIterator it = findElement(me, _key);

	return (it && (key(it) == _key));
}

//////////////////////////////////////////////////////////////////////////////
// SkiplistIterator: Standard Iterator for Skiplist
//////////////////////////////////////////////////////////////////////////////

template <typename TSkiplist>
class Iter< TSkiplist, SkiplistIterator>
{
public:
	typedef typename Key<TSkiplist>::Type TKey;
	typedef typename Value<TSkiplist>::Type TValue;
	typedef typename Spec<TSkiplist>::Type TSpec;

	typedef SkiplistElement<TKey, TValue, TSpec> TElement;

	TElement * data_pointer;

	Iter()
		: data_pointer(0)
	{
	}
	Iter(Iter const & other)
		: data_pointer(other.data_pointer)
	{
	}
	Iter(TElement * ptr)
		: data_pointer(ptr)
	{
	}
	Iter(TSkiplist & sk)
		: data_pointer(sk.data_border.data_next[0].data_element)
	{
	}
	~Iter()
	{
	}

	Iter const & 
	operator = (Iter const & other)
	{
		data_pointer = other.data_pointer;
		return *this;
	}
	operator bool () const
	{
		return data_pointer;
	}

};

//////////////////////////////////////////////////////////////////////////////

template <typename TSkiplist>
inline bool
atEnd(Iter<TSkiplist, SkiplistIterator> & it)
{
	return (!it.data_pointer);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSkiplist>
inline void
goNext(Iter<TSkiplist, SkiplistIterator> & it)
{
	it.data_pointer = it.data_pointer->data_next[0].data_element;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSkiplist>
inline typename Value<TSkiplist>::Type &
value(Iter<TSkiplist, SkiplistIterator> & it)
{
	return it.data_pointer->data.value;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSkiplist>
inline typename Key<TSkiplist>::Type const &
key(Iter<TSkiplist, SkiplistIterator> & it)
{
	return it.data_pointer->data.key;
}

////////////////////////////////////////////////////////////////////////////////

//____________________________________________________________________________

//////////////////////////////////////////////////////////////////////////////

}

#endif

