/*
 *  misc_set.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_MISC_SET_H
#define SEQAN_HEADER_MISC_SET_H

#include <set>

namespace SEQAN_NAMESPACE_MAIN
{

	//////////////////////////////////////////////////////////////////////////////
	//

	template <typename TKey>
	class VectorSet {
	public:
		typedef TKey							TValue;
		typedef bool							TSetEntry;
		typedef TSetEntry						TSetVector[ValueSize<TKey>::VALUE];
		typedef typename Size<TSetVector>::Type	TSize;

		TSetVector	vector;
		TSize		size;

		VectorSet():
			size(0)	{}
	};

	template <typename TKey, typename TObject>
	class VectorSet< Pair<TKey, TObject> > {
	public:
		typedef Pair<TKey, TObject>				TValue;
		typedef bool							TSetEntry;
		typedef TSetEntry						TSetVector[ValueSize<TKey>::VALUE];
		typedef TObject							TObjVector[ValueSize<TKey>::VALUE];
		typedef typename Size<TSetVector>::Type	TSize;

		TSetVector	vector;
		TObjVector	obj;
		TSize		size;

		VectorSet():
			size(0) {}
	};

	template <typename TObject>
	struct Value< VectorSet<TObject> > {
		typedef TObject Type;
	};
	template <typename TObject>
	struct Size< VectorSet<TObject> > {
		typedef typename VectorSet<TObject>::TSize Type;
	};

	template <typename TObject>
	inline typename Size< VectorSet<TObject> >::Type length(VectorSet<TObject> const &set) {
		return set.size;
	}



	template <typename TSet>
	struct _Set_SetVector {
		typedef void* Type;
	};
	
	template <typename TSet>
	struct _Set_ObjVector {
		typedef void* Type;
	};



	template <typename TKey>
	struct _Set_SetVector< VectorSet<TKey> >
	{
		typedef typename VectorSet<TKey>::TSetVector Type;
	};

	template <typename TKey>
	struct _Set_ObjVector< VectorSet<TKey> >
	{
		typedef typename VectorSet<TKey>::TSetVector Type;
	};


	template <typename TKey, typename TObject>
	struct _Set_SetVector< VectorSet< Pair<TKey, TObject> > >
	{
		typedef typename VectorSet< Pair<TKey, TObject> >::TSetVector Type;
	};

	template <typename TKey, typename TObject>
	struct _Set_ObjVector< VectorSet< Pair<TKey, TObject> > >
	{
		typedef typename VectorSet< Pair<TKey, TObject> >::TObjVector Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	//

	template <typename TPair>
	struct SetLess : public ::std::binary_function<TPair, TPair, bool> {
		bool operator() (TPair const &a, TPair const &b) {
			return a.i1 < b.i1;
		}
	};

	template <typename TKey>
	struct Set {
		typedef ::std::set<TKey> Type;
	};
	template <typename TKey, typename TObject>
	struct Set< Pair<TKey, TObject> > {
		typedef ::std::set< Pair<TKey, TObject>, SetLess< Pair<TKey, TObject> > > Type;
	};


	template <typename TValue, typename TSpec>
	struct Set< SimpleType<TValue, TSpec> > {
		typedef VectorSet< SimpleType<TValue, TSpec> > Type;
	};
	template <typename TValue, typename TSpec, typename TObject>
	struct Set< Pair<SimpleType<TValue, TSpec>, TObject> > {
		typedef VectorSet< Pair<SimpleType<TValue, TSpec>, TObject> > Type;
	};


	template <>
	struct Set<char> {
		typedef VectorSet<char> Type;
	};
	template <typename TObject>
	struct Set< Pair<char, TObject> > {
		typedef VectorSet< Pair<char, TObject> > Type;
	};


	template <>
	struct Set<signed char> {
		typedef VectorSet<signed char> Type;
	};
	template <typename TObject>
	struct Set< Pair<signed char, TObject> > {
		typedef VectorSet< Pair<signed char, TObject> > Type;
	};

	template <>
	struct Set<unsigned char> {
		typedef VectorSet<unsigned char> Type;
	};
	template <typename TObject>
	struct Set< Pair<unsigned char, TObject> > {
		typedef VectorSet< Pair<unsigned char, TObject> > Type;
	};



	//////////////////////////////////////////////////////////////////////////////
	//

	struct _VectorSetIterator;
	typedef Tag<_VectorSetIterator> VectorSetIterator;

	template <typename TVectorSet>
	class Iter< TVectorSet, VectorSetIterator > {
		typedef typename _Set_SetVector<TVectorSet>::Type	TSetVector;
		typedef typename _Set_ObjVector<TVectorSet>::Type	TObjVector;

		typedef Iter								iterator;
		typedef typename Value<TSetVector>::Type	TValue, value;
		typedef typename Iterator<TSetVector>::Type	TSetPointer;
		typedef typename Iterator<TObjVector>::Type	TObjPointer;

		TSetPointer	ptr, begin, end;
		TObjPointer	obj;

		Iter(TSetPointer _ptr, TSetPointer _begin, TSetPointer _end, TObjPointer _obj):
			ptr(_ptr), begin(_begin), end(_end), obj(_obj)
		{
			while (!(*ptr) && ptr != end) {
				++ptr;
				++obj;
			}
		}

		inline TValue operator*() { return TValue(ptr - begin, *obj); }

		inline iterator operator++() {
			++ptr;
			while (!(*ptr) && ptr != end) {
				++ptr;
				++obj;
			}
			return *this;
		}
	};

	template <typename TVectorSet>
	class Iter< TVectorSet const, VectorSetIterator > {
		typedef typename _Set_SetVector<TVectorSet>::Type	TSetVector;
		typedef typename _Set_ObjVector<TVectorSet>::Type	TObjVector;

		typedef Iter								iterator;
		typedef typename Value<TSetVector>::Type	TValue, value;
		typedef typename Iterator<TSetVector>::Type	TSetPointer;
		typedef typename Iterator<TObjVector>::Type	TObjPointer;

		TSetPointer	ptr, begin, end;
		TObjPointer	obj;

		Iter(TSetPointer _ptr, TSetPointer _begin, TSetPointer _end, TObjPointer _obj):
			ptr(_ptr), begin(_begin), end(_end), obj(_obj)
		{
			while (!(*ptr) && ptr != end)	{
				++ptr;
				++obj;
			}
		}

		inline TValue operator*() { return TValue(ptr - begin, *obj); }

		inline iterator operator++() {
			++ptr;
			while (!(*ptr) && ptr != end)	{
				++ptr;
				++obj;
			}
			return *this;
		}
	};


	template <typename TObject>
	struct Iterator< VectorSet<TObject> > {
		typedef Iter<VectorSet<TObject>, VectorSetIterator> Type;
	};
	template <typename TObject>
	struct Iterator< VectorSet<TObject> const> {
		typedef Iter<VectorSet<TObject> const, VectorSetIterator> Type;
	};

	template <typename TValue>
	finline void clear(VectorSet<TValue> &set) {
		memset<sizeof(typename VectorSet<TValue>::TSetVector), 0>(set.vector);
		set.size = 0;
	}

	template <typename TKey, typename TSetKey>
	finline void insert(TKey const &key, VectorSet<TSetKey> &set) {
		if (!set.vector[(unsigned)key]) {
			++set.size;
			set.vector[(unsigned)key] = true;
		}
	}
	template <typename TPair, typename TSetKey, typename TSetObject>
	finline void insert(TPair const &pair, VectorSet< Pair<TSetKey, TSetObject> > &set) {
		if (set.vector[(unsigned)pair.i1]) {
			++set.size;
			set.vector[(unsigned)pair.i1] = true;
		}
		set.obj[(unsigned)pair.i1] = pair.i2;
	}

	template <typename TKey, typename TValue>
	finline void erase(TKey const &key, VectorSet<TValue> &set) {
		if (set.vector[(unsigned)key]) {
			--set.size;
			set.vector[(unsigned)key] = false;
		}
	}

	template <typename TKey, typename TValue>
	finline bool in(TKey const &key, VectorSet<TValue> const &set) {
		return set.vector[(unsigned)key];
	}

	template <typename TKey, typename TSetKey>
	inline typename Iterator< VectorSet<TSetKey> >::Type find(TKey const &key, VectorSet<TSetKey> &set) {
		if (in(key, set))
			return Iter<VectorSet<TSetKey>, VectorSetIterator>(&set.vector[(unsigned)key], begin(set.vector), end(set.vector), begin(set.obj));
		else
			return end(set);
	}
	template <typename TKey, typename TSetKey>
	inline typename Iterator< VectorSet<TSetKey> const>::Type find(TKey const &key, VectorSet<TSetKey> const &set) {
		if (in(key, set))
			return Iter<VectorSet<TSetKey> const, VectorSetIterator>(&set.vector[(unsigned)key], begin(set.vector), end(set.vector), begin(set.obj));
		else
			return end(set);
	}
	template <typename TKey, typename TSetKey, typename TSetObject>
	inline Iter<VectorSet< Pair<TSetKey, TSetObject> >, VectorSetIterator> find(TKey const &key, VectorSet< Pair<TSetKey, TSetObject> > &set) {
		if (in(key, set))
			return Iter<VectorSet< Pair<TSetKey, TSetObject> >, VectorSetIterator>(&set.vector[(unsigned)key], begin(set.vector), end(set.vector), begin(set.obj));
		else
			return end(set);
	}
	template <typename TKey, typename TSetKey, typename TSetObject>
	inline Iter<VectorSet< Pair<TSetKey, TSetObject> > const, VectorSetIterator> find(TKey const &key, VectorSet< Pair<TSetKey, TSetObject> > const &set) {
		if (in(key, set))
			return Iter<VectorSet< Pair<TSetKey, TSetObject> > const, VectorSetIterator>(&set.vector[(unsigned)key], begin(set.vector), end(set.vector), begin(set.obj));
		else
			return end(set);
	}



	//////////////////////////////////////////////////////////////////////////////
	//

	template <typename TObject>
	struct Value< ::std::set<TObject> > {
		typedef TObject Type;
	};
	template <typename TObject>
	struct Size< ::std::set<TObject> > {
		typedef typename ::std::set<TObject>::size_type Type;
	};


	template <typename TKey>
	inline void clear(::std::set<TKey> &set) {
		set.clear();
	}

	template <typename TKey, typename TSetKey>
	inline void insert(TKey const &key, ::std::set<TSetKey> &set) {
		set.insert(key);
	}

	template <typename TKey, typename TSetKey>
	inline void erase(TKey const &key, ::std::set<TSetKey> &set) {
		set.erase(key);
	}
	template <typename TKey, typename TSetKey, typename TSetObject>
	inline void erase(TKey const &key, ::std::set< Pair<TSetKey, TSetObject> > &set) {
		set.erase(Pair<TSetKey, TSetObject>(key, TSetObject()));
	}

	template <typename TKey, typename TSetKey>
	inline bool in(TKey const &key, ::std::set<TSetKey> const &set) {
		return set.count(key) != 0;
	}
	template <typename TKey, typename TSetKey, typename TSetObject>
	inline bool in(TKey const &key, ::std::set<Pair<TSetKey, TSetObject> > const &set) {
		return set.count(Pair<TSetKey, TSetObject>(key, TSetObject())) != 0;
	}

	template <typename TKey>
	inline typename Size< ::std::set<TKey> >::Type length(::std::set<TKey> const &set) {
		return set.size();
	}

	template <typename TKey, typename TSetKey>
	inline typename Iterator< ::std::set<TSetKey> >::Type find(TKey const &key, ::std::set<TSetKey> &set) {
		return set.find(key);
	}
	template <typename TKey, typename TSetKey>
	inline typename Iterator< ::std::set<TSetKey> const>::Type find(TKey const &key, ::std::set<TSetKey> const &set) {
		return set.find(key);
	}
	template <typename TKey, typename TSetKey, typename TSetObject>
	inline typename Iterator< ::std::set<Pair<TSetKey, TSetObject> > >::Type find(TKey const &key, ::std::set<Pair<TSetKey, TSetObject> > &set) {
		return set.find(Pair<TSetKey, TSetObject>(key, TSetObject()));
	}
	template <typename TKey, typename TSetKey, typename TSetObject>
	inline typename Iterator< ::std::set<Pair<TSetKey, TSetObject> > const>::Type find(TKey const &key, ::std::set<Pair<TSetKey, TSetObject> > const &set) {
		return set.find(Pair<TSetKey, TSetObject>(key, TSetObject()));
	}



	template <typename TObject>
	inline typename Iterator< VectorSet<TObject> >::Type begin(VectorSet<TObject> &set) {
		return Iter<VectorSet<TObject>, VectorSetIterator> (begin(set.vector), begin(set.vector), end(set.vector), begin(set.obj));
	}
	template <typename TObject>
	inline typename Iterator< VectorSet<TObject> const>::Type begin(VectorSet<TObject> const &set) {
		return Iter<VectorSet<TObject>, VectorSetIterator> (begin(set.vector), begin(set.vector), end(set.vector), begin(set.obj));
	}
	template <typename TObject>
	inline typename Iterator< VectorSet<TObject> >::Type end(VectorSet<TObject> &set) {
		return Iter<VectorSet<TObject>, VectorSetIterator> (begin(set.vector), begin(set.vector), end(set.vector), begin(set.obj));
	}
	template <typename TObject>
	inline typename Iterator< VectorSet<TObject> const>::Type end(VectorSet<TObject> const &set) {
		return Iter<VectorSet<TObject>, VectorSetIterator> (begin(set.vector), begin(set.vector), end(set.vector), begin(set.obj));
	}

	template <typename TObject>
	inline bool operator==(Iter<TObject, VectorSetIterator> const &a, Iter<TObject, VectorSetIterator> const &b) {
		return a.ptr == b.ptr;
	}
	template <typename TObject>
	inline bool operator!=(Iter<TObject, VectorSetIterator> const &a, Iter<TObject, VectorSetIterator> const &b) {
		return a.ptr != b.ptr;
	}
	template <typename TObject>
	inline bool eof(Iter<TObject, VectorSetIterator> const &a) {
		return a.ptr == a.end;
	}

	template <typename TKey>
	inline TKey keyOf(Iter<VectorSet<TKey>, VectorSetIterator> const &it) {
		return it.ptr - it.begin;
	}
	template <typename TKey>
	inline TKey keyOf(Iter<VectorSet<TKey> const, VectorSetIterator> const &it) {
		return it.ptr - it.begin;
	}
	template <typename TKey, typename TObject>
	inline TKey keyOf(Iter<VectorSet< Pair<TKey, TObject> >, VectorSetIterator> const &it) {
		return it.ptr - it.begin;
	}
	template <typename TKey, typename TObject>
	inline TKey keyOf(Iter<VectorSet< Pair<TKey, TObject> const>, VectorSetIterator> const &it) {
		return it.ptr - it.begin;
	}
	template <typename TKey, typename TObject>
	inline TObject& objectOf(Iter<VectorSet< Pair<TKey, TObject> >, VectorSetIterator> &it) {
		return (*it.ptr).i2;
	}
	template <typename TKey, typename TObject>
	inline TObject const & objectOf(Iter<VectorSet< Pair<TKey, TObject> const>, VectorSetIterator> const &it) {
		return (*it.ptr).i2;
	}



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

	template <typename TKey>
	inline TKey& keyOf(typename ::std::set<TKey>::iterator &it) {
		return *it;
	}
	template <typename TKey>
	inline TKey const & keyOf(typename ::std::set<TKey>::const_iterator const &it) {
		return *it;
	}
	template <typename TKey, typename TObject>
	inline TKey& keyOf(typename ::std::set<Pair<TKey, TObject> >::iterator &it) {
		return (*it).i1;
	}
	template <typename TKey, typename TObject>
	inline TKey const & keyOf(typename ::std::set<Pair<TKey, TObject> >::const_iterator const &it) {
		return (*it).i1;
	}
	template <typename TKey, typename TObject>
	inline TObject& objectOf(typename ::std::set<Pair<TKey, TObject> >::iterator &it) {
		return (*it).i2;
	}
	template <typename TKey, typename TObject>
	inline TObject const & objectOf(typename ::std::set<Pair<TKey, TObject> >::const_iterator &it) {
		return (*it).i2;
	}


}

#endif
