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

#ifndef SEQAN_HEADER_MISC_BASE_H
#define SEQAN_HEADER_MISC_BASE_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

	//////////////////////////////////////////////////////////////////////////////
	// In SeqAn sets and maps store elements as pairs of (key,object) 
	// the elements of sets without objects are the keys.
	//////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////
	// Key/Object meta-functions
	//////////////////////////////////////////////////////////////////////////////

/* moved to basic_h, see #6
	template <typename TElement>
	struct Key {
		typedef TElement Type;
	};
*/
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

//////////////////////////////////////////////////////////////////////////////

/*trash
	template <typename TElement>
	struct Object {
		typedef Nothing Type;
	};

	template <typename TKey, typename TObject, typename TSpec>
	struct Object< Pair<TKey, TObject, TSpec> > {
		typedef TObject Type;
	};
*/

	//////////////////////////////////////////////////////////////////////////////
	// keyOf function
	//////////////////////////////////////////////////////////////////////////////

	template <typename TElement>
	inline typename Key<TElement>::Type & 
	key(TElement & element) 
	{
		return element;
	}
	template <typename TElement>
	inline typename Key<TElement const>::Type & 
	key(TElement const & element) 
	{
		return element;
	}

	template <typename TKey, typename TObject, typename TSpec>
	inline TKey & 
	key(Pair<TKey, TObject, TSpec> &element) {
		return element.i1;
	}
	template <typename TKey, typename TObject, typename TSpec>
	inline TKey const &
	key(Pair<TKey, TObject, TSpec> const &element) {
		return element.i1;
	}

	//////////////////////////////////////////////////////////////////////////////
	// objectOf
	//////////////////////////////////////////////////////////////////////////////

	template <typename TKey, typename TObject, typename TSpec>
	inline TObject & 
	cargo(Pair<TKey, TObject, TSpec> &element) {
		return element.i2;
	}
	template <typename TKey, typename TObject, typename TSpec>
	inline TObject const &
	cargo(Pair<TKey, TObject, TSpec> const &element) {
		return element.i2;
	}

}

#endif

