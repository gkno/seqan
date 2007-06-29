/*
 *  misc_set.h
 *  SeqAn
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_MISC_SET_H
#define SEQAN_HEADER_MISC_SET_H

#include <set>
#include "misc_base.h"

namespace SEQAN_NAMESPACE_MAIN
{

	//////////////////////////////////////////////////////////////////////////////
	// A SeqAn set expects 2 functions:
	// 1. a unary function to convert an element into an ordinal (for VectorSets)
	// 2. a binary Less-function to compare 2 elements
	//////////////////////////////////////////////////////////////////////////////

	template <typename TKey, typename TOrdinal = unsigned>
	struct SetFunctors : public ::std::binary_function<TKey, TKey, bool>
	{
		// key to ordinal
		inline TOrdinal operator() (TKey const &x) {
			return (TOrdinal)x;
		}
		// key less
		inline bool operator() (TKey const &a, TKey const &b) {
			return a < b;
		}
	};


	//////////////////////////////////////////////////////////////////////////////
	// forward declaration

	template <typename TElement, typename TFunctors, typename TSpec>
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
	template <typename TElement, typename TFunctors, typename TSpec>
	struct _SetSetVector< VectorSet<TElement, TFunctors, TSpec> > {
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
	template <typename TKey, typename TObject, typename TPairSpec, typename TFunctors, typename TSpec>
	struct _SetObjVector< VectorSet<Pair<TKey, TObject, TPairSpec>, TFunctors, TSpec> > {
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
		typename TFunctors = SetFunctors< typename Key<TElement>::Type >,
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

		template <typename _TElement, typename _TFunctors>
		inline void _autoSize(VectorSet<_TElement, _TFunctors, Alloc<> > &) {
			resize(vector, (unsigned)ValueSize<_TElement>::VALUE);
		}

		template <typename _TKey, typename _TObject, typename _TFunctors, typename _TSpec>
		inline void _autoSize(VectorSet<Pair<_TKey, _TObject, _TSpec>, _TFunctors, Alloc<> > &) {
			resize(vector, (unsigned)ValueSize<_TKey>::VALUE);
			resize(obj, (unsigned)ValueSize<_TKey>::VALUE);
		}
	};


	//////////////////////////////////////////////////////////////////////////////
	// set meta-functions

	template <typename TElement, typename TFunctors, typename TSpec>
	struct Value< VectorSet<TElement, TFunctors, TSpec> > {
		typedef TElement Type;
	};
	template <typename TElement, typename TFunctors, typename TSpec>
	struct Size< VectorSet<TElement, TFunctors, TSpec> >:
		Size< _SetSetVector< VectorSet<TElement, TFunctors, TSpec> > > {};

	template <typename TElement, typename TFunctors, typename TSpec>
	struct Key< VectorSet<TElement, TFunctors, TSpec> > :
		Key<TElement> {};

	template <typename TElement, typename TFunctors, typename TSpec>
	struct Object< VectorSet<TElement, TFunctors, TSpec> > :
		Object<TElement> {};

	template <typename TObject, typename TFunctors, typename TSpec>
	inline typename Size< VectorSet<TObject, TFunctors, TSpec> >::Type 
	length(VectorSet<TObject, TFunctors, TSpec> const &set) {
		return set.size;
	}



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
			SetFunctors< Pair<TKey, TObject, TPairSpec> > > Type;
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
/*
	template <typename TVectorSet>
	class Iter< TVectorSet const, VectorSetIterator > {
		typedef typename TVectorSet::TSetVector		TSetVector;
		typedef typename TVectorSet::TObjVector		TObjVector;

		typedef Iter								iterator;
		typedef typename Value<TSetVector>::Type	TValue, value;
		typedef typename Iterator<TSetVector>::Type	TSetIter;
		typedef typename Iterator<TObjVector>::Type	TObjIter;

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
*/

	template <typename TObject, typename TFunctors, typename TSpec>
	struct Iterator< VectorSet<TObject, TFunctors, TSpec> > {
		typedef Iter<VectorSet<TObject, TFunctors, TSpec>, VectorSetIterator> Type;
	};
	template <typename TObject, typename TFunctors, typename TSpec>
	struct Iterator< VectorSet<TObject, TFunctors, TSpec> const> {
		typedef Iter<VectorSet<TObject, TFunctors, TSpec> const, VectorSetIterator> Type;
	};

	template <typename TValue, typename TFunctors, typename TSpec>
	finline void 
	clear(VectorSet<TValue, TFunctors, TSpec> &set) {
		arrayFill(begin(set.vector), end(set.vector), false);
		set.size = 0;
	}

	template <typename TKey, typename TSetKey, typename TFunctors, typename TSpec>
	finline void 
	insert(TKey const &key, VectorSet<TSetKey, TFunctors, TSpec> &set) {
		if (!set.vector[(unsigned)key]) {
			++set.size;
			set.vector[(unsigned)key] = true;
		}
	}
	template <
		typename TPair, 
		typename TSetKey, 
		typename TSetObject, 
		typename TPairSpec, 
		typename TFunctors, 
		typename TSpec>
	finline void 
	insert(TPair const &pair, VectorSet< Pair<TSetKey, TSetObject, TPairSpec>, TFunctors, TSpec > &set) {
		if (!set.vector[(unsigned)pair.i1]) {
			++set.size;
			set.vector[(unsigned)pair.i1] = true;
		}
		set.obj[(unsigned)pair.i1] = pair.i2;
	}

	template <typename TKey, typename TValue, typename TFunctors, typename TSpec>
	finline void 
	erase(TKey const &key, VectorSet<TValue, TFunctors, TSpec> &set) {
		if (set.vector[(unsigned)key]) {
			--set.size;
			set.vector[(unsigned)key] = false;
		}
	}

	template <typename TKey, typename TValue, typename TFunctors, typename TSpec>
	finline bool 
	in(TKey const &key, VectorSet<TValue, TFunctors, TSpec> const &set) {
		return set.vector[(unsigned)key];
	}

	template <typename TKey, typename TSetKey, typename TFunctors, typename TSpec>
	inline typename Iterator< VectorSet<TSetKey, TFunctors, TSpec> >::Type 
	find(TKey const &key, VectorSet<TSetKey, TFunctors, TSpec> &set) {
		if (in(key, set))
			return Iter<VectorSet<TSetKey, TFunctors, TSpec>, VectorSetIterator>
				(begin(set.vector, Rooted()) + (unsigned)key, begin(set.obj, Standard()) + (unsigned)key);
		else
			return end(set);
	}
	template <typename TKey, typename TSetKey, typename TFunctors, typename TSpec>
	inline typename Iterator< VectorSet<TSetKey, TFunctors, TSpec> const>::Type 
	find(TKey const &key, VectorSet<TSetKey, TFunctors, TSpec> const &set) {
		if (in(key, set))
			return Iter<VectorSet<TSetKey, TFunctors, TSpec> const, VectorSetIterator>
				(begin(set.vector, Rooted()) + (unsigned)key, begin(set.obj, Standard()) + (unsigned)key);
		else
			return end(set);
	}
	template <
		typename TKey, 
		typename TSetKey, 
		typename TSetObject, 
		typename TPairSpec, 
		typename TFunctors, 
		typename TSpec>
	inline Iter<VectorSet< Pair<TSetKey, TSetObject, TPairSpec>, TFunctors, TSpec >, VectorSetIterator> 
	find(TKey const &key, VectorSet< Pair<TSetKey, TSetObject, TPairSpec>, TFunctors, TSpec > &set) {
		if (in(key, set))
			return Iter<VectorSet< Pair<TSetKey, TSetObject, TPairSpec>, TFunctors, TSpec >, VectorSetIterator>
				(begin(set.vector, Rooted()) + (unsigned)key, begin(set.obj, Standard()) + (unsigned)key);
		else
			return end(set);
	}
	template <
		typename TKey, 
		typename TSetKey, 
		typename TSetObject, 
		typename TPairSpec, 
		typename TFunctors,
		typename TSpec>
	inline Iter<VectorSet< Pair<TSetKey, TSetObject, TPairSpec>, TFunctors, TSpec > const, VectorSetIterator> 
	find(TKey const &key, VectorSet< Pair<TSetKey, TSetObject, TPairSpec>, TFunctors, TSpec > const &set) {
		if (in(key, set))
			return Iter<VectorSet< Pair<TSetKey, TSetObject, TPairSpec>, TFunctors, TSpec > const, VectorSetIterator>
				(begin(set.vector, Rooted()) + (unsigned)key, begin(set.obj, Standard()) + (unsigned)key);
		else
			return end(set);
	}



	//////////////////////////////////////////////////////////////////////////////
	//

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

	template <typename TKey, typename TSetKey>
	inline void 
	insert(TKey const &key, ::std::set<TSetKey> &set) {
		set.insert(key);
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



	template <typename TElement, typename TFunctors, typename TSpec>
	inline typename Iterator< VectorSet<TElement, TFunctors, TSpec> >::Type 
	begin(VectorSet<TElement, TFunctors, TSpec> &set) {
		return Iter<VectorSet<TElement, TFunctors, TSpec>, VectorSetIterator> 
			(begin(set.vector, Rooted()), begin(set.obj, Standard()));
	}
	template <typename TElement, typename TFunctors, typename TSpec>
	inline typename Iterator< VectorSet<TElement, TFunctors, TSpec> const>::Type 
	begin(VectorSet<TElement, TFunctors, TSpec> const &set) {
		return Iter<VectorSet<TElement, TFunctors, TSpec> const, VectorSetIterator> 
			(begin(set.vector, Rooted()), begin(set.obj, Standard()));
	}
	template <typename TElement, typename TFunctors, typename TSpec>
	inline typename Iterator< VectorSet<TElement, TFunctors, TSpec> >::Type 
	end(VectorSet<TElement, TFunctors, TSpec> &set) {
		return Iter<VectorSet<TElement, TFunctors, TSpec>, VectorSetIterator> 
			(end(set.vector, Rooted()), begin(set.obj, Standard()));
	}
	template <typename TElement, typename TFunctors, typename TSpec>
	inline typename Iterator< VectorSet<TElement, TFunctors, TSpec> const>::Type 
	end(VectorSet<TElement, TFunctors, TSpec> const &set) {
		return Iter<VectorSet<TElement, TFunctors, TSpec> const, VectorSetIterator> 
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
	keyOf(typename ::std::set<TElement>::iterator const &it) {
		return keyOf(*it);
	}
	template <typename TElement>
	inline typename Key< ::std::set<TElement> >::Type const &
	keyOf(typename ::std::set<TElement>::const_iterator const &it) {
		return keyOf(*it);
	}

	template <typename TElement>
	inline typename Object< ::std::set<TElement> >::Type &
	objectOf(typename ::std::set<TElement>::iterator const &it) {
		return objectOf(*it);
	}
	template <typename TElement>
	inline typename Object< ::std::set<TElement> >::Type const &
	objectOf(typename ::std::set<TElement>::const_iterator const &it) {
		return objectOf(*it);
	}

}

#endif
