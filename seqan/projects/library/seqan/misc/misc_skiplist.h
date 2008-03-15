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

template <typename TKey, typename TValue, typename TSpec>
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


struct _SkiplistFindLast_; //search with _skiplistFind the last element in Skiplist
typedef Tag<_SkiplistFindLast_> _SkiplistFindLast;

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct AllocatorType;

template <typename TKey, typename TValue, typename TSpec>
struct AllocatorType<Skiplist<TKey, TValue, TSpec> >
{
	typedef Allocator<SimpleAlloc<> > Type;
};

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec>
class Skiplist
{
public:
	typedef typename AllocatorType<Skiplist>::Type TAllocator;
	typedef SkiplistElement<TKey, TValue, TSpec> TElement;
	typedef typename Size<Skiplist>::Type TSize;

	enum
	{
		MAX_HEIGHT = 26;
		BLOCK_SIZE = 0x400;
	};


	Holder<TAllocator> data_allocator;
	TElement * data_recycle[TElement::MAX_HEIGHT];
	unsigned char * data_mem_begin;
	unsigned char * data_mem_end;

	SkiplistElement data_border;
	TSize data_length;
	unsigned char data_height;

	Skiplist()
		: data_mem_begin(0)
		, data_mem_end(0)
		, data_length(0)
		, data_height(0)
	{
		arrayFill(data_recycle, data_recycle + MAX_HEIGHT, 0);
		arrayFill(data_border.data_next, data_border.data_next + MAX_HEIGHT, 0);

		mtRandInit();
	}
	~Skiplist()
	{
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
		MAX_HEIGHT = TSkiplist::MAX_HEIGHT;
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
		MAX_HEIGHT = TSkiplist::MAX_HEIGHT;
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

	SkiplistElement * data_element;
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
		MAX_HEIGHT = TSkiplist::MAX_HEIGHT;
	};

	TElement * data_elements[MAX_HEIGHT];
//	unsigned char data_height; //is stored in Skiplist
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec, typename TIteratorSpec>
struct Iterator<Skiplist<TKey, TValue, TSpec>, TIteratorSpec>
{
	typedef Iter<Skiplist<TKey, TValue, TSpec>, SkiplistIterator> Type;
};

//////////////////////////////////////////////////////////////////////////////
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
		me.data_recycle[height] = reinterpret_cast<TElement *>(* me.data_recycle[height]);
	}
	else
	{
		int need_size = sizeof(typename TElement::Data) + height * sizeof(TNext);
		int buf_size = me.data_mem_end - me.data_mem_begin
		if (buf_size < need_size)
		{//need new memory
			if (buf_size >= sizeof(typename TElement::Data))
			{//link rest memory in recycle 
				int rest_height = (buf_size - sizeof(typename TElement::Data)) / sizeof(TNext); //must be < height, because buf_size < need_size
				* reinterpret_cast<TElement *>(me.data_mem_begin) = me.data_recycle[rest_height];
				me.data_recycle[rest_height] = reinterpret_cast<TElement *>(me.data_mem_begin);
			}
			allocate(me), me.data_mem_begin, TSkiplist::BLOCK_SIZE, TagAllocateStorage());
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
						   SkiplistElement<TKey, TValue, TSpec> & el
						   unsigned char height)
{
	* reinterpret_cast<TElement *>(&el) = me.data_recycle[height]
	me.data_recycle[height] = reinterpret_cast<TElement *>(&el);
	//the real deallocation is done by the allocator on destruction 
}

//////////////////////////////////////////////////////////////////////////////

//Destroys New SkiplistElement
template <typename TKey, typename TValue, typename TSpec>
inline void
_skiplistDestructElement(Skiplist<TKey, TValue, TSpec> & me,
						 SkiplistElement<TKey, TValue, TSpec> & el
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
						 SkiplistElement<TKey, void, TSpec> & el
						 unsigned char height)
{
	valueDestruct(& (el.data.key) );
	//no need to construct the next array
	_skiplistDeallocateElement(me, el, height);
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
		if (path.elements[height]->data_next[height2] == el) break;
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

template <typename TKey, typename TValue, typename TSpec, typename TTag>
inline bool
_skiplistFindGoNext(SkiplistNext<TKey, TValue, TSpec> & next,
					unsigned char,
					TKey const & key,
					Tag<TTag> const)
{
	return (next.data_element->data.key < key);
}

template <typename TKey, typename TValue, typename TSpec, typename TTag>
inline bool
_skiplistFindGoNext(SkiplistNext<TKey, TValue, TSpec> & next,
					unsigned char,
					SkiplistElement<TKey, TValue, TSpec> const & el,
					Tag<TTag> const)
{
	return (next.data_element->data.key <= el.data.key) && (next.data_element != & el);
}

template <typename TKey, typename TValue, typename TSpec>
inline bool
_skiplistFindGoNext(SkiplistNext<TKey, TValue, TSpec> & next,
					unsigned char height,
					TFind const &,
					_SkiplistFindLast)
{
	return next.data_element->data_next[height];
}


template <typename TKey, typename TValue, typename TSpec, typename TFind, typename TTag>
inline void
_skiplistFind(Skiplist<TKey, TValue, TSpec> & me,
			  TFind const & find, //can be a TKey or a SkiplistElement
			  /*OUT*/ SkiplistPath<TKey, TValue, TSpec> & path,
			  Tag<TTag> const tag)
{
	typedef SkiplistElement<TKey, TValue, TSpec> TElement;
	typedef SkiplistNext<TKey, TValue, TSpec> TNext;

	TElement * here = me.data_border;

	for (int i = me.data_height; i >= 0; --i)
	{
		TNext & next = here.data_next[i];
		while (next.data_element && _skiplistFindGoNext(next, i, find, tag))
		{
			here = next.data_element;
		}
		path.data_elements[i] = here;
	}
}
template <typename TKey, typename TValue, typename TSpec, typename TFind>
inline void
_skiplistFind(Skiplist<TKey, TValue, TSpec> & me,
			  TFind const & find, //can be a TKey or a SkiplistElement
			  /*OUT*/ SkiplistPath<TKey, TValue, TSpec> & path)
{
	_skiplistFind(me, find, path, Default());
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
		el.data_next[i].element = path.data_elements[i];
		path.data_elements[i]->data_next[i].element = el;
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

//creates height for new SkiplistElement.
//increases the Skiplist height if necessary
template <typename TKey, typename TValue, typename TSpec>
inline unsigned char
_skiplistCreateHeight(Skiplist<TKey, TValue, TSpec> & me)
{
	unsigned char height = geomRand<unsigned char>();
	if (height >= TSkiplist::MAX_HEIGHT) height = TSkiplist::MAX_HEIGHT-1;

	if (height > me.data_height) me.data_height = height;

	return height;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec>
inline void
insertElement(Skiplist<TKey, TValue, TSpec> & me,
			  TKey const & key,
			  TValue const & value)
{
	typedef Skiplist<TKey, TValue, TSpec> TSkiplist;
	typedef SkiplistElement<TKey, TValue, TSpec> TElement;

	unsigned char height = _skiplistCreateHeight(me);
	TElement & el = _skiplistConstructElement(me, height, key, value);
	_skiplistInsertElement(me, el, height);
}
template <typename TKey, typename TSpec>
inline void
insertElement(Skiplist<TKey, void, TSpec> & me,
			  TKey const & key)
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
					   SkiplistPath<TKey, TValue, TSpec> & el)
{
	typedef SkiplistPath<TKey, TValue, TSpec> TPath;

	TPath path;
	_skiplistFind(me, el, path);

	for (int i = me.data_height; i >= 0; --i)
	{
		if (path.elements[i]->data_next[i].element == el)
		{
			path.elements[i]->data_next[i].element = el.data_next[i].element;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// SkiplistIterator
//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec>
class Iter< Skiplist<TKey, TValue, TSpec>, SkiplistIterator>
{
public:
	typedef Skiplist<TKey, TValue, TSpec> TSkiplist;
	typedef SkiplistElement<TKey, TValue, TSpec> TElement;
	typedef SkiplistPath<TKey, TValue, TSpec> TPath;

	TPath data_path;
	TElement * data_pointer;
	TSkiplist * data_container;
	bool data_is_path_valid;

	Iter()
		: data_pointer(0)
		, data_container(0)
		, data_is_path_valid(false)
	{
	}
	Iter(TSkiplist & sk)
		: data_container(& sk)
	{
		goBegin(*this);
	}
	Iter(Iter const & other)
		: data_pointer(other.data_pointer)
		, data_container(other.data_container)
		, data_is_path_valid(other.data_is_path_valid)
	{
		for (int i = data_container->data_height; i > 0 ; --i)
		{
			data_path[i] = other.data_path[i];
		}
	}
	~Iter()
	{
	}

	Iter const & 
	operator = (Iter const & other)
	{
		data_pointer = other.data_pointer;
		data_container = other.data_container;
		data_is_path_valid = other.data_is_path_valid;
		for (int i = data_container->data_height; i > 0 ; --i)
		{
			data_path[i] = other.data_path[i];
		}
	}

};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec>
inline void
goBegin(Iter< Skiplist<TKey, TValue, TSpec>, SkiplistIterator> & it)
{
	typedef Skiplist<TKey, TValue, TSpec> TSkiplist;
	TSkiplist & sk = * (it.data_container);

	it.data_pointer = sk.data_border.data_next[0].data_element;

	for (int i = sk.data_height; i > 0 ; --i)
	{
		it.data_path[i] = & (sk.data_border);
	}

	it.data_is_path_valid = true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec>
inline void
goLast(Iter< Skiplist<TKey, TValue, TSpec>, SkiplistIterator> & it)
{
	_skiplistFind(* it.data_container, 0, it.data_path, _SkiplistFindLast());
	id.data_pointer = id.data_path.data_elements[0]->data_next[0];
	it.data_is_path_valid = true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec>
inline bool
atEnd(Iter< Skiplist<TKey, TValue, TSpec>, SkiplistIterator> & it)
{
	return (!it.data_pointer)
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec, typename TFind>
inline bool
findElement(Iter< Skiplist<TKey, TValue, TSpec>, SkiplistIterator> & it,
			TFind const & find) //find can be a key or an SkiplistElement
{
	_skiplistFind(* it.data_container, key, it.data_path);
	it.data_pointer = id.data_path.data_elements[0]->data_next[0];
	it.data_is_path_valid = true;

	return it.data_pointer;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec, typename TKey2, typename TValue2>
inline bool
insertElement(Iter< Skiplist<TKey, TValue, TSpec>, SkiplistIterator> & it,
			  TKey2 const & key,
			  TValue2 const & value)
{
	typedef SkiplistElement<TKey, TValue, TSpec> TElement;

	if (!it.data_is_path_valid) _skiplistValidatePath(it);

	unsigned char height = _skiplistCreateHeight(me);
	TElement & el = _skiplistConstructElement(*it.data_container, height, key, value);

	_skiplistFind(* it.data_container, key, it.data_path);
	it.data_pointer = id.data_path.data_elements[0]->data_next[0];
	it.data_is_path_valid = true;

	return it.data_pointer;
}
	typedef Skiplist<TKey, TValue, TSpec> TSkiplist;
	typedef SkiplistElement<TKey, TValue, TSpec> TElement;

	unsigned char height = _skiplistCreateHeight(me);
	TElement & el = _skiplistConstructElement(me, height, key, value);
	_skiplistInsertElement(me, el, height);

//////////////////////////////////////////////////////////////////////////////

//____________________________________________________________________________

//////////////////////////////////////////////////////////////////////////////

}

#endif

