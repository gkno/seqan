 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2008

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: $
 ==========================================================================*/

#ifndef SEQAN_HEADER_MAP_ADAPTER_STL_H
#define SEQAN_HEADER_MAP_ADAPTER_STL_H

//////////////////////////////////////////////////////////////////////////////

#include <map>
#include <set>

//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

template <typename TKey, typename TCargo>
struct Value< ::std::map<TKey,TCargo> > 
{
    typedef Pair<TKey,TCargo> Type;
};

template <typename TKey>
struct Value< ::std::set<TKey> > 
{
    typedef bool Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCargo>
struct Key< ::std::map<TKey,TCargo>  >
{
    typedef TKey Type;
};

template <typename TKey>
struct Key< ::std::set<TKey>  >
{
    typedef TKey Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCargo>
struct Cargo< ::std::map<TKey,TCargo> >
{
    typedef TCargo Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey>
struct Size< ::std::set<TKey> >
{
    //typedef ::std::set<TKey>::size_type Type;
    typedef unsigned Type;
};

template <typename TKey, typename TCargo>
struct Size< ::std::map<TKey,TCargo> >
{
    //typedef ::std::map<TKey,TCargo>::size_type Type;
    typedef unsigned Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey>
inline void
assign(::std::set<TKey> & target,
	   ::std::set<TKey> const & source)
{
    target = source;
}

template <typename TKey,typename TCargo>
inline void
assign(::std::map<TKey,TCargo> & target,
	   ::std::map<TKey,TCargo> const & source)
{
    target = source;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TKey2>
inline bool
hasKey(::std::set<TValue> & me, TKey2 const & _key){
    return (me.count(_key) != 0);
}

template <typename TKey, typename TCargo, typename TKey2>
inline bool
hasKey(::std::map<TKey,TCargo> & me, TKey2 const & _key){
    return (me.count(_key) != 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline typename Size< ::std::set<TValue> >::Type
length(::std::set<TValue> const & me)
{
	return me.size();
}

template <typename TKey, typename TCargo>
inline typename Size< ::std::map<TKey,TCargo> >::Type
length(::std::map<TKey,TCargo> const & me)
{
	return me.size();
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TValue2>
inline void
insert(::std::set<TValue> & me,TValue2 const & _value)
{
    me.insert(_value);
}

template <typename TKey, typename TCargo, typename TKey2>
inline void
insert(::std::map<TKey,TCargo> & me,TKey2 const & _key)
{
    me.insert(::std::make_pair(_key,TCargo()));
}

template <typename TKey, typename TCargo, typename TKey2, typename TCargo2, typename TSpec>
inline void
insert(::std::map<TKey,TCargo> & me, Pair<TKey2,TCargo2,TSpec> const & _value)
{
    me.insert(::std::make_pair(_value.i1,_value.i2));
}

template <typename TKey, typename TCargo, typename TKey2, typename TCargo2>
inline void
insert(::std::map<TKey,TCargo> & me, TKey2 const & _key, TCargo2 const & _cargo)
{
    me.insert(::std::make_pair(_key,_cargo));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline void
clear(::std::set<TValue> & me)
{
    me.clear();
}

template <typename TKey, typename TCargo>
inline void
clear(::std::map<TKey,TCargo> & me)
{
    me.clear();
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCargo, typename TKey2>
inline typename Value< ::std::map<TKey,TCargo> >::Type &
value(::std::map<TKey,TCargo> & me,
	  TKey2 const & _key)
{    
    return me[_key];
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCargo, typename TKey2>
inline typename Cargo< ::std::map<TKey,TCargo> >::Type &
cargo(::std::map<TKey,TCargo> & me,
	  TKey2 const & _key)
{
	return me[_key];
}

//////////////////////////////////////////////////////////////////////////////

struct STLSetIterator;
struct STLMapIterator;

//////////////////////////////////////////////////////////////////////////////
//                            MapIterator                                   //
//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCargo, typename TIteratorSpec>
struct Iterator< ::std::map<TKey,TCargo> , TIteratorSpec >
{
    typedef ::std::map<TKey,TCargo> TSTLMap;
	typedef Iter<TSTLMap, STLMapIterator> Type;
};

template <typename TSTLMap>
class Iter< TSTLMap, STLMapIterator>
{
public:    
    typedef typename ::std::map<typename Key<TSTLMap>::Type, typename Cargo<TSTLMap>::Type>::iterator THostIter;
	THostIter _iter;
    Holder<TSTLMap> host_map_holder;

	Iter()
	{
	}

	Iter(Iter const & other)
		: _iter(other._iter)
	{
        host_map_holder = other.host_map_holder;
	}

	Iter(TSTLMap & map)
		: _iter(map.begin())
    {
        setValue(host_map_holder,map);
	}

	~Iter()
	{
	}

	Iter const & 
	operator = (Iter const & other)
	{
        host_map_holder = other.host_map_holder;
        _iter = other._iter;
        return *this;
	}
	operator bool () const
	{
        return (_iter != value(host_map_holder).end());
	}

};

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline bool
operator == (Iter<TSTLMap, STLMapIterator> const & left,
			 Iter<TSTLMap, STLMapIterator> const & right)
{
	return left._iter == right._iter;
}

template <typename TSTLMap>
inline bool
operator != (Iter<TSTLMap, STLMapIterator> const & left,
			 Iter<TSTLMap, STLMapIterator> const & right)
{
    return left._iter != right._iter;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCargo, typename TIteratorSpec>
inline typename Iterator< ::std::map<TKey,TCargo>, TIteratorSpec>::Type
begin(::std::map<TKey,TCargo> & me,
	  TIteratorSpec)
{
    typedef ::std::map<TKey,TCargo> TSTLMap;
    typedef typename Iterator<TSTLMap , TIteratorSpec>::Type TIterator;
    return TIterator(me);
}
template <typename TKey, typename TCargo>
inline typename Iterator< ::std::map<TKey,TCargo> >::Type
begin(::std::map<TKey,TCargo> & me)
{
    typedef ::std::map<TKey,TCargo> TSTLMap;
    typedef typename Iterator<TSTLMap>::Type TIterator;
    return TIterator(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCargo, typename TIteratorSpec>
inline typename Iterator< ::std::map<TKey,TCargo> >::Type
end(::std::map<TKey,TCargo> & me,
	  TIteratorSpec)
{
    typedef ::std::map<TKey,TCargo> TSTLMap;
	typedef typename Iterator<TSTLMap, TIteratorSpec>::Type TIterator;
	TIterator _iter(me);
    _iter._iter = me.end();
    return _iter;
}

template <typename TKey, typename TCargo>
inline typename Iterator< ::std::map<TKey,TCargo> >::Type
end(::std::map<TKey,TCargo> & me)
{
    typedef ::std::map<TKey,TCargo> TSTLMap;
	typedef typename Iterator<TSTLMap>::Type TIterator;
	TIterator _iter(me);
    _iter._iter = me.end();
    return _iter;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline bool
atEnd(Iter<TSTLMap, STLMapIterator> & it)
{
	return (it._iter == value(it.host_map_holder).end());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline void
goNext(Iter<TSTLMap, STLMapIterator> & it)
{
	it._iter++;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline typename Value<TSTLMap>::Type &
value(Iter<TSTLMap, STLMapIterator> & it)
{
    return it._iter->second;
}
template <typename TSTLMap>
inline typename Value<TSTLMap>::Type &
value(Iter<TSTLMap, STLMapIterator> const & it)
{
	return it._iter->second;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline typename Key<TSTLMap>::Type const &
key(Iter<TSTLMap, STLMapIterator> & it)
{
	return it._iter->first;
}
template <typename TSTLMap>
inline typename Key<TSTLMap>::Type const &
key(Iter<TSTLMap, STLMapIterator> const & it)
{
	return it._iter->first;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline typename Cargo<TSTLMap>::Type &
cargo(Iter<TSTLMap, STLMapIterator> & it)
{
	return it._iter->second;
}
template <typename TSTLMap>
inline typename Cargo<TSTLMap>::Type &
cargo(Iter<TSTLMap, STLMapIterator> const & it)
{
	return it._iter->second;
}

////////////////////////////////////////////////////////////////////////////////
//                             SetIterator                                    //
////////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TIteratorSpec>
struct Iterator< ::std::set<TValue> , TIteratorSpec >
{
    typedef ::std::set<TValue> TSTLMap;
	typedef Iter<TSTLMap, STLSetIterator> Type;
};

template <typename TSTLMap>
class Iter< TSTLMap, STLSetIterator>
{
public:    
    typedef typename ::std::set<typename Key<TSTLMap>::Type>::iterator THostIter;
	THostIter _iter;
    Holder<TSTLMap> host_map_holder;

	Iter()
	{
	}

	Iter(Iter const & other)
		: _iter(other._iter)
	{
        host_map_holder = other.host_map_holder;
	}

	Iter(TSTLMap & map)
		: _iter(map.begin())
    {
        setValue(host_map_holder,map);
	}

	~Iter()
	{
	}

	Iter const & 
	operator = (Iter const & other)
	{
        host_map_holder = other.host_map_holder;
        _iter = other._iter;
        return *this;
	}
	operator bool () const
	{
        return (_iter != value(host_map_holder).end());
	}

};

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline bool
operator == (Iter<TSTLMap, STLSetIterator> const & left,
			 Iter<TSTLMap, STLSetIterator> const & right)
{
	return left._iter == right._iter;
}

template <typename TSTLMap>
inline bool
operator != (Iter<TSTLMap, STLSetIterator> const & left,
			 Iter<TSTLMap, STLSetIterator> const & right)
{
    return left._iter != right._iter;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TIteratorSpec>
inline typename Iterator< ::std::set<TValue>, TIteratorSpec>::Type
begin(::std::set<TValue> & me,
	  TIteratorSpec)
{
    typedef ::std::set<TValue> TSTLMap;
    typedef typename Iterator<TSTLMap , TIteratorSpec>::Type TIterator;
    return TIterator(me);
}
template <typename TValue>
inline typename Iterator< ::std::set<TValue> >::Type
begin(::std::set<TValue> & me)
{
    typedef ::std::set<TValue> TSTLMap;
    typedef typename Iterator<TSTLMap>::Type TIterator;
    return TIterator(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TIteratorSpec>
inline typename Iterator< ::std::set<TValue>, TIteratorSpec>::Type
end(::std::set<TValue> & me,
	  TIteratorSpec)
{
    typedef ::std::set<TValue> TSTLMap;
	typedef typename Iterator<TSTLMap, TIteratorSpec>::Type TIterator;
	TIterator _iter(me);
    _iter._iter = me.end();
    return _iter;
}

template <typename TValue>
inline typename Iterator< ::std::set<TValue> >::Type
end(::std::set<TValue> & me)
{
    typedef ::std::set<TValue> TSTLMap;
	typedef typename Iterator<TSTLMap>::Type TIterator;
	TIterator _iter(me);
    _iter._iter = me.end();
    return _iter;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline bool
atEnd(Iter<TSTLMap, STLSetIterator> & it)
{
	return (it._iter == value(it.host_map_holder).end());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline void
goNext(Iter<TSTLMap, STLSetIterator> & it)
{
	it._iter++;
}

//////////////////////////////////////////////////////////////////////////////
//
//template <typename TSTLMap>
//inline typename Value<TSTLMap>::Type &
//value(Iter<TSTLMap, STLSetIterator> & it)
//{
//    return hasKey(*it);
//}
//template <typename TSTLMap>
//inline typename Value<TSTLMap>::Type &
//value(Iter<TSTLMap, STLSetIterator> const & it)
//{
//	return hasKey(*it);
//}

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline typename Key<TSTLMap>::Type const &
key(Iter<TSTLMap, STLSetIterator> & it)
{
    return (*it._iter);
}
template <typename TSTLMap>
inline typename Key<TSTLMap>::Type const &
key(Iter<TSTLMap, STLSetIterator> const & it)
{
	return (*it._iter);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TValue,typename TFind>
inline typename Iterator< ::std::set<TValue> >::Type
find(::std::set<TValue> & me,
	 TFind const & _find)
{
    typedef ::std::set<TValue> TMap;
    typedef typename Iterator< TMap >::Type TMapIterator;  

    TMapIterator _iter(me);
    _iter._iter = me.find(_find);
    if(!_iter){
        _iter._iter = me.upper_bound(_find); 
    }
    return _iter;
}

template <typename TKey, typename TCargo ,typename TFind>
inline typename Iterator< ::std::map<TKey,TCargo> >::Type
find(::std::map<TKey,TCargo> & me,
	 TFind const & _find)
{
    typedef ::std::map<TKey,TCargo> TMap;
    typedef typename Iterator< TMap >::Type TMapIterator;  

    TMapIterator _iter(me);
    _iter._iter = me.find(_find);
    if(!_iter){
        _iter._iter = me.upper_bound(_find); 
    }
    return _iter;
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TMap2>
inline void
erase(::std::set<TValue> & me,
	  Iter<TMap2, STLSetIterator> const & it)
{
    me.erase(it._iter);
}

template <typename TKey, typename TCargo ,typename TMap2>
inline void
erase(::std::map<TKey,TCargo> & me,
	  Iter<TMap2, STLMapIterator> const & it)
{
    me.erase(it._iter);
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TToRemove>
inline void
erase(::std::set<TValue> & me,
	  TToRemove const & to_remove)
{
    me.erase(to_remove);
}

template <typename TKey, typename TCargo ,typename TToRemove>
inline void
erase(::std::map<TKey,TCargo> & me,
	  TToRemove const & to_remove)
{
    me.erase(to_remove);
}

}
#endif // #ifndef SEQAN_HEADER
