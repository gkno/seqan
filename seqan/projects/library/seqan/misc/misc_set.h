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

#ifndef SEQAN_HEADER_MISC_SET_H
#define SEQAN_HEADER_MISC_SET_H

#include <set>
#include "misc_base.h"

namespace SEQAN_NAMESPACE_MAIN
{

	//////////////////////////////////////////////////////////////////////////////
	// 
	//	VectorSet
	// 
	//////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	// forward declaration

	template <typename TElement, typename TSpec>
	class VectorSet;


	//////////////////////////////////////////////////////////////////////////////
	// internal set meta-functions

	template <typename TElement>
	struct _VectorSetKeySize {
		enum { VALUE = ValueSize< typename Key<TElement>::Type >::VALUE };
	};


	template <typename TSet>
	struct _SetSetVector {
		typedef void* Type;
	};
	template <typename TElement, typename TSpec>
	struct _SetSetVector< VectorSet<TElement, TSpec> > {
		typedef String<bool, TSpec> Type;
	};
	template <typename TSet>
	struct _SetSetVector<TSet const> {
		typedef typename _SetSetVector<TSet>::Type const Type;
	};



	template <typename TSet>
	struct _SetObjVector {
		typedef void* Type;
	};
	template <typename TKey, typename TObject, typename TPairSpec, typename TSpec>
	struct _SetObjVector< VectorSet<Pair<TKey, TObject, TPairSpec>, TSpec> > {
		typedef String<TObject, TSpec> Type;
	};
	template <typename TSet>
	struct _SetObjVector<TSet const> {
		typedef typename _SetObjVector<TSet>::Type const Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	// VectorSet class

	template <
		typename TElement = char,
		typename TSpec = Alloc<> /*Array<_VectorSetKeySize<TKey>::VALUE>*/
	>
	class VectorSet {
	public:
		typedef typename _SetSetVector<VectorSet>::Type		TSetVector;
		typedef typename _SetObjVector<VectorSet>::Type		TObjVector;
		typedef typename Size<VectorSet>::Type				TSize;

		TSetVector	vector;
		TObjVector	obj;
		TSize		size;

		VectorSet():
			size(0)	
		{
			_autoSize(*this);
			clear(*this);
		}

		VectorSet(TSize _vectorSize):
			size(0)	
		{
			resize(vector, _vectorSize);
			clear(*this);
		}
		
		template <typename _TSet>
		inline void _autoSize(_TSet &) {}

		template <typename _TElement>
		inline void _autoSize(VectorSet<_TElement, Alloc<> > &) {
			resize(vector, (unsigned)ValueSize<_TElement>::VALUE);
		}

		template <typename _TKey, typename _TObject, typename _TSpec>
		inline void _autoSize(VectorSet<Pair<_TKey, _TObject, _TSpec>, Alloc<> > &) {
			resize(vector, (unsigned)ValueSize<_TKey>::VALUE);
			resize(obj, (unsigned)ValueSize<_TKey>::VALUE);
		}
	};


	//////////////////////////////////////////////////////////////////////////////
	// set meta-functions

	template <typename TElement, typename TSpec>
	struct Value< VectorSet<TElement, TSpec> > {
		typedef TElement Type;
	};
	template <typename TElement, typename TSpec>
	struct Size< VectorSet<TElement, TSpec> >:
		Size< _SetSetVector< VectorSet<TElement, TSpec> > > {};

	template <typename TElement, typename TSpec>
	struct Key< VectorSet<TElement, TSpec> > :
		Key<TElement> {};

	template <typename TElement, typename TSpec>
	struct Object< VectorSet<TElement, TSpec> > :
		Object<TElement> {};

	template <typename TObject, typename TSpec>
	inline typename Size< VectorSet<TObject, TSpec> >::Type 
	length(VectorSet<TObject, TSpec> const &set) {
		return set.size;
	}



	//////////////////////////////////////////////////////////////////////////////
	//

	struct _VectorSetIterator;
	typedef Tag<_VectorSetIterator> VectorSetIterator;

	template <typename TVectorSet>
	class Iter< TVectorSet, VectorSetIterator > 
	{
		typedef typename _SetSetVector<TVectorSet>::Type		TSetVector;
		typedef typename _SetObjVector<TVectorSet>::Type		TObjVector;

		typedef Iter											iterator;
		typedef typename Value<TSetVector>::Type				TValue, value;
		typedef typename Iterator<TSetVector, Rooted>::Type		TSetIter;
		typedef typename Iterator<TObjVector, Standard>::Type	TObjIter;

	public:
		TSetIter	ptr;
		TObjIter	obj;

		Iter() {}

		Iter(TSetIter _ptr, TObjIter _obj):
			ptr(_ptr), obj(_obj)
		{
			while (!(*ptr) && !atEnd(ptr)) {
				++ptr;
				++obj;
			}
		}

//		inline TValue operator*() { return TValue(ptr - begin, *obj); }

		inline iterator operator++() {
			++ptr;
			++obj;
			while (!(*ptr) && !atEnd(ptr)) {
				++ptr;
				++obj;
			}
			return *this;
		}
	};


	template <typename TObject, typename TSpec>
	struct Iterator< VectorSet<TObject, TSpec> > {
		typedef Iter<VectorSet<TObject, TSpec>, VectorSetIterator> Type;
	};
	template <typename TObject, typename TSpec>
	struct Iterator< VectorSet<TObject, TSpec> const> {
		typedef Iter<VectorSet<TObject, TSpec> const, VectorSetIterator> Type;
	};

	template <typename TValue, typename TSpec>
	finline void 
	clear(VectorSet<TValue, TSpec> &set) {
		arrayFill(begin(set.vector), end(set.vector), false);
		set.size = 0;
	}

	template <typename TElement, typename TSetKey, typename TSpec>
	finline void 
	insert(TElement const &element, VectorSet<TSetKey, TSpec> &set) {
		if (!set.vector[(unsigned)(element)]) {
			++set.size;
			set.vector[(unsigned)(element)] = true;
		}
	}
	template <
		typename TElement,
		typename TSetKey, 
		typename TSetObject, 
		typename TPairSpec, 
		typename TSpec>
	finline void 
	insert(TElement const &element, VectorSet< Pair<TSetKey, TSetObject, TPairSpec>, TSpec > &set) {
		if (!set.vector[(unsigned)(keyOf(element))]) {
			++set.size;
			set.vector[(unsigned)(keyOf(element))] = true;
		}
		set.obj[(unsigned)(keyOf(element))] = objectOf(element);
	}

	template <typename TKey, typename TValue, typename TSpec>
	finline void 
	erase(TKey const &key, VectorSet<TValue, TSpec> &set) {
		if (set.vector[(unsigned)(key)]) {
			--set.size;
			set.vector[(unsigned)(key)] = false;
		}
	}

	template <typename TKey, typename TValue, typename TSpec>
	finline bool 
	in(TKey const &key, VectorSet<TValue, TSpec> const &set) {
		return set.vector[(unsigned)(key)];
	}

	template <typename TKey, typename TKey2, typename TSpec>
	inline typename Iterator< VectorSet<TKey2, TSpec> >::Type 
	find(TKey const &key, VectorSet<TKey2, TSpec> &set) {
		if (in(key, set))
			return Iter<VectorSet<TKey2, TSpec>, VectorSetIterator>
				(begin(set.vector, Rooted()) + (unsigned)(key),
				 begin(set.obj, Standard()) + (unsigned)(key));
		else
			return end(set);
	}
	template <typename TKey, typename TKey2, typename TSpec>
	inline typename Iterator< VectorSet<TKey2, TSpec> const>::Type 
	find(TKey const &key, VectorSet<TKey2, TSpec> const &set) {
		if (in(key, set))
			return Iter<VectorSet<TKey2, TSpec> const, VectorSetIterator>
				(begin(set.vector, Rooted()) + (unsigned)(key), 
				 begin(set.obj, Standard()) + (unsigned)(key));
		else
			return end(set);
	}

	template <typename TElement, typename TSpec>
	inline typename Iterator< VectorSet<TElement, TSpec> >::Type 
	begin(VectorSet<TElement, TSpec> &set) {
		return Iter<VectorSet<TElement, TSpec>, VectorSetIterator> 
			(begin(set.vector, Rooted()), begin(set.obj, Standard()));
	}
	template <typename TElement, typename TSpec>
	inline typename Iterator< VectorSet<TElement, TSpec> const>::Type 
	begin(VectorSet<TElement, TSpec> const &set) {
		return Iter<VectorSet<TElement, TSpec> const, VectorSetIterator> 
			(begin(set.vector, Rooted()), begin(set.obj, Standard()));
	}

	template <typename TElement, typename TSpec>
	inline typename Iterator< VectorSet<TElement, TSpec> >::Type 
	end(VectorSet<TElement, TSpec> &set) {
		return Iter<VectorSet<TElement, TSpec>, VectorSetIterator> 
			(end(set.vector, Rooted()), begin(set.obj, Standard()));
	}
	template <typename TElement, typename TSpec>
	inline typename Iterator< VectorSet<TElement, TSpec> const>::Type 
	end(VectorSet<TElement, TSpec> const &set) {
		return Iter<VectorSet<TElement, TSpec> const, VectorSetIterator> 
			(end(set.vector, Rooted()), begin(set.obj, Standard()));
	}

	template <typename TSet>
	inline bool 
	operator==(Iter<TSet, VectorSetIterator> const &a, Iter<TSet, VectorSetIterator> const &b) {
		return a.ptr == b.ptr;
	}
	template <typename TSet>
	inline bool 
	operator!=(Iter<TSet, VectorSetIterator> const &a, Iter<TSet, VectorSetIterator> const &b) {
		return a.ptr != b.ptr;
	}
	template <typename TSet>
	inline bool 
	eof(Iter<TSet, VectorSetIterator> const &a) {
		return atEnd(a.ptr);
	}
	template <typename TSet>
	inline bool 
	atEnd(Iter<TSet, VectorSetIterator> const &a) {
		return atEnd(a.ptr);
	}

	//////////////////////////////////////////////////////////////////////////////

	template <typename TSet>
	inline typename Key<TSet>::Type
	keyOf(Iter<TSet, VectorSetIterator> &it) {
		return position(it.ptr);
	}
	template <typename TSet>
	inline typename Key<TSet>::Type
	keyOf(Iter<TSet, VectorSetIterator> const &it) {
		return position(it.ptr);
	}

	template <typename TSet>
	inline typename Object<TSet>::Type &
	objectOf(Iter<TSet, VectorSetIterator> &it) {
		return *it.obj;
	}
	template <typename TSet>
	inline typename Object<TSet>::Type &
	objectOf(Iter<TSet, VectorSetIterator> const &it) {
		return *it.obj;
	}



//____________________________________________________________________________




	//////////////////////////////////////////////////////////////////////////////
	// 
	//	STL set adaptation
	// 
	//////////////////////////////////////////////////////////////////////////////


	template <typename TElement>
	struct _SetLess : public ::std::binary_function<TElement, TElement, bool>
	{
		// key less
		inline bool operator() (TElement const &a, TElement const &b) {
			return keyOf(a) < keyOf(b);
		}
	};

	template <typename TElement>
	struct Value< ::std::set<TElement> > {
		typedef TElement Type;
	};
	template <typename TElement>
	struct Size< ::std::set<TElement> > {
		typedef typename ::std::set<TElement>::size_type Type;
	};

	template <typename TElement>
	struct Key< ::std::set<TElement> > :
		Key<TElement> {};

	template <typename TElement>
	struct Object< ::std::set<TElement> > :
		Object<TElement> {};



	template <typename TKey>
	inline void 
	clear(::std::set<TKey> &set) {
		set.clear();
	}

	template <typename TElement, typename TSetKey>
	inline void 
	insert(TElement const &element, ::std::set<TSetKey> &set) {
		set.insert(element);
	}

	template <typename TKey, typename TSetKey>
	inline void 
	erase(TKey const &key, ::std::set<TSetKey> &set) {
		set.erase(key);
	}
	template <typename TKey, typename TSetKey, typename TSetObject, typename TPairSpec>
	inline void 
	erase(TKey const &key, ::std::set< Pair<TSetKey, TSetObject, TPairSpec> > &set) {
		set.erase(Pair<TSetKey, TSetObject, TPairSpec>(key, TSetObject()));
	}

	template <typename TKey, typename TSetKey>
	inline bool 
	in(TKey const &key, ::std::set<TSetKey> const &set) {
		return set.count(key) != 0;
	}
	template <typename TKey, typename TSetKey, typename TSetObject, typename TPairSpec>
	inline bool 
	in(TKey const &key, ::std::set<Pair<TSetKey, TSetObject, TPairSpec> > const &set) {
		return set.count(Pair<TSetKey, TSetObject, TPairSpec>(key, TSetObject())) != 0;
	}

	template <typename TKey>
	inline typename Size< ::std::set<TKey> >::Type 
	length(::std::set<TKey> const &set) {
		return set.size();
	}

	template <typename TKey, typename TSetKey>
	inline typename Iterator< ::std::set<TSetKey> >::Type 
	find(TKey const &key, ::std::set<TSetKey> &set) {
		return set.find(key);
	}
	template <typename TKey, typename TSetKey>
	inline typename Iterator< ::std::set<TSetKey> const>::Type 
	find(TKey const &key, ::std::set<TSetKey> const &set) {
		return set.find(key);
	}
	template <typename TKey, typename TSetKey, typename TSetObject, typename TPairSpec>
	inline typename Iterator< ::std::set<Pair<TSetKey, TSetObject, TPairSpec> > >::Type 
	find(TKey const &key, ::std::set<Pair<TSetKey, TSetObject> > &set) {
		return set.find(Pair<TSetKey, TSetObject>(key, TSetObject()));
	}
	template <typename TKey, typename TSetKey, typename TSetObject>
	inline typename Iterator< ::std::set<Pair<TSetKey, TSetObject> > const>::Type 
	find(TKey const &key, ::std::set<Pair<TSetKey, TSetObject> > const &set) {
		return set.find(Pair<TSetKey, TSetObject>(key, TSetObject()));
	}



	//////////////////////////////////////////////////////////////////////////////

	template <typename TObject>
	struct Iterator< ::std::set<TObject> > {
		typedef typename ::std::set<TObject>::iterator Type;
	};
	template <typename TObject>
	struct Iterator< ::std::set<TObject> const > {
		typedef typename ::std::set<TObject>::const_iterator Type;
	};


	template <typename TObject>
	typename Iterator< ::std::set<TObject> >::Type begin(::std::set<TObject> &set) {
		return set.begin();
	}
	template <typename TObject>
	typename Iterator< ::std::set<TObject> const >::Type begin(::std::set<TObject> const &set) {
		return set.begin();
	}
	template <typename TObject>
	typename Iterator< ::std::set<TObject> >::Type end(::std::set<TObject> &set) {
		return set.end();
	}
	template <typename TObject>
	typename Iterator< ::std::set<TObject> const >::Type end(::std::set<TObject> const &set) {
		return set.end();
	}

	//////////////////////////////////////////////////////////////////////////////

	template <typename TElement>
	inline typename Key< ::std::set<TElement> >::Type &
	keyOf(typename ::std::set<TElement>::iterator &it) {
		return keyOf(*it);
	}
	template <typename TElement>
	inline typename Key< ::std::set<TElement> >::Type &
	keyOf(typename ::std::set<TElement>::iterator const &it) {
		return keyOf(*it);
	}
	template <typename TElement>
	inline typename Key< ::std::set<TElement> >::Type const &
	keyOf(typename ::std::set<TElement>::const_iterator &it) {
		return keyOf(*it);
	}
	template <typename TElement>
	inline typename Key< ::std::set<TElement> >::Type const &
	keyOf(typename ::std::set<TElement>::const_iterator const &it) {
		return keyOf(*it);
	}

	template <typename TElement>
	inline typename Object< ::std::set<TElement> >::Type &
	objectOf(typename ::std::set<TElement>::iterator &it) {
		return objectOf(*it);
	}
	template <typename TElement>
	inline typename Object< ::std::set<TElement> >::Type &
	objectOf(typename ::std::set<TElement>::iterator const &it) {
		return objectOf(*it);
	}
	template <typename TElement>
	inline typename Object< ::std::set<TElement> >::Type const &
	objectOf(typename ::std::set<TElement>::const_iterator &it) {
		return objectOf(*it);
	}
	template <typename TElement>
	inline typename Object< ::std::set<TElement> >::Type const &
	objectOf(typename ::std::set<TElement>::const_iterator const &it) {
		return objectOf(*it);
	}



//____________________________________________________________________________




	//////////////////////////////////////////////////////////////////////////////
	// 
	//	Set chooser
	// 
	//////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	// Set meta-function to choose an efficient implementation

	template <typename TKey>
	struct Set {
		typedef ::std::set<TKey> Type;
	};
	template <typename TKey, typename TObject, typename TPairSpec>
	struct Set< Pair<TKey, TObject, TPairSpec> > {
		typedef ::std::set< 
			Pair<TKey, TObject, TPairSpec>, 
			_SetLess< Pair<TKey, TObject, TPairSpec> > > Type;
	};


	template <typename TValue, typename TSpec>
	struct Set< SimpleType<TValue, TSpec> > {
		typedef VectorSet< SimpleType<TValue, TSpec> > Type;
	};
	template <typename TValue, typename TSpec, typename TObject, typename TPairSpec>
	struct Set< Pair<SimpleType<TValue, TSpec>, TObject, TPairSpec> > {
		typedef VectorSet< Pair<SimpleType<TValue, TSpec>, TObject, TPairSpec> > Type;
	};


	template <>
	struct Set<char> {
		typedef VectorSet<char> Type;
	};
	template <typename TObject, typename TPairSpec>
	struct Set< Pair<char, TObject, TPairSpec> > {
		typedef VectorSet< Pair<char, TObject, TPairSpec> > Type;
	};


	template <>
	struct Set<signed char> {
		typedef VectorSet<signed char> Type;
	};
	template <typename TObject, typename TPairSpec>
	struct Set< Pair<signed char, TObject, TPairSpec> > {
		typedef VectorSet< Pair<signed char, TObject, TPairSpec> > Type;
	};

	template <>
	struct Set<unsigned char> {
		typedef VectorSet<unsigned char> Type;
	};
	template <typename TObject, typename TPairSpec>
	struct Set< Pair<unsigned char, TObject, TPairSpec> > {
		typedef VectorSet< Pair<unsigned char, TObject, TPairSpec> > Type;
	};

}

#endif
