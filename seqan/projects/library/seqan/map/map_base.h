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
  $Id: $
 ==========================================================================*/

#ifndef SEQAN_HEADER_MAP_BASE_H
#define SEQAN_HEADER_MAP_BASE_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
//insertion tags

template <typename TSpec = Default>
struct Skiplist;

template <typename TElement, typename TSpec = Skiplist<> >
class Map;


//////////////////////////////////////////////////////////////////////////////
// In SeqAn sets and maps store elements as pairs of (key,cargo) 
// the elements of sets without objects are the keys.
//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TObject, typename TSpec>
struct Key< Pair<TKey, TObject, TSpec> > 
{
	typedef TKey Type;
};

template <typename TKey, typename TCargo, typename TSpec>
struct Cargo< Pair<TKey, TCargo, TSpec> > 
{
	typedef TCargo Type;
};


//////////////////////////////////////////////////////////////////////////////
// Type for mapValue function that implements [] for map types 

template <typename TMap, typename TCargo>
struct _MapValue_Impl
{
	typedef TCargo & Type;
};
template <typename TMap>
struct _MapValue_Impl<TMap, Nothing>
{
	typedef bool Type;
};

template <typename TMap>
struct MapValue:
	_MapValue_Impl< TMap, typename Cargo<TMap>::Type >
{
};



template <typename TCargo>
struct _Impl_mapValue
{
	template <typename TMap, typename TKey2>
	static inline TCargo &
	mapValue_(TMap & me,
		TKey2 const & _key)
	{
		return cargo(me, _key);
	}
};

template <>
struct _Impl_mapValue<Nothing>
{
	template <typename TMap, typename TKey2>
	static inline bool
	mapValue_(TMap & me,
		TKey2 const & _key)
	{
		return hasKey(me, _key);
	}
};

template <typename TMap, typename TKey>
inline typename MapValue<TMap>::Type
mapValue(TMap & me,
		 TKey const & _key)
{
	typedef typename Cargo<TMap>::Type TCargo;
	return _Impl_mapValue<TCargo>::mapValue_(me, _key);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TElement>
inline TElement & 
key(TElement & element) 
{
	return element;
}
template <typename TElement>
inline TElement const & 
key(TElement const & element) 
{
	return element;
}

template <typename TKey, typename TObject, typename TSpec>
inline TKey & 
key(Pair<TKey, TObject, TSpec> & element) 
{
	return element.i1;
}
template <typename TKey, typename TObject, typename TSpec>
inline TKey const &
key(Pair<TKey, TObject, TSpec> const & element) 
{
	return element.i1;
}

//////////////////////////////////////////////////////////////////////////////
//no default cargo function

template <typename TKey, typename TObject, typename TSpec>
inline TObject & 
cargo(Pair<TKey, TObject, TSpec> & element) 
{
	return element.i2;
}
template <typename TKey, typename TObject, typename TSpec>
inline TObject const &
cargo(Pair<TKey, TObject, TSpec> const & element) 
{
	return element.i2;
}

//////////////////////////////////////////////////////////////////////////////

}

#endif

