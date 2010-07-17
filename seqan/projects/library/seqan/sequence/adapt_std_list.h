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
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
 ==========================================================================*/

// SEQAN_NO_GENERATED_FORWARDS: No forwards are generated for this file.

#ifndef SEQAN_SEQUENCE_ADAPT_STD_LIST_H_
#define SEQAN_SEQUENCE_ADAPT_STD_LIST_H_

#include <list>

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ===========================================================================
// Metafunctions
// ===========================================================================

// TODO(holtgrew): dddoc this

template <typename TValue>
struct Size<std::list<TValue> >
{
	typedef size_t Type;
};

template <typename TValue>
struct Size<std::list<TValue> const>
        : Size<std::list<TValue> > {};

template <typename TValue>
struct Value<std::list<TValue> >
{
	typedef TValue Type;
};

template <typename TValue>
struct Iterator<std::list<TValue> >
{
	typedef typename std::list<TValue>::iterator Type;
};

template <typename TValue>
struct Iterator<std::list<TValue> const>
{
	typedef typename std::list<TValue>::const_iterator Type;
};

template <typename TValue>
struct Reference<std::list<TValue> >
{
	typedef typename std::list<TValue>::iterator _TIterator;
    typedef typename std::iterator_traits<_TIterator>::reference Type;
};

template <typename TValue>
struct Reference<std::list<TValue> const>
{
	typedef typename std::list<TValue>::const_iterator _TIterator;
    typedef typename std::iterator_traits<_TIterator>::reference Type;
};

// ===========================================================================
// Functions
// ===========================================================================

template<typename TValue>
inline typename Iterator<std::list<TValue>, Standard >::Type
begin(std::list<TValue> & list)
{
	return list.begin();
}

template<typename TValue>
inline typename Iterator<std::list<TValue>, Standard >::Type
end(std::list<TValue> & list)
{
	return list.end();
}

template <typename TValue>
inline typename Iterator<std::list<TValue> const, Standard >::Type
begin(const std::list<TValue> & list)
{
	return list.begin();
}

template <typename TValue>
inline typename Iterator<std::list<TValue> const, Standard >::Type
end(std::list<TValue> const & list)
{
	return list.end();
}

template <typename TValue>
inline typename Reference<std::list<TValue> >::Type
front(std::list<TValue> & list)
{
    return list.front();
}

template <typename TValue>
inline typename Reference<std::list<TValue> const>::Type
front(std::list<TValue> const & list)
{
    return list.front();
}

template <typename TValue>
inline typename Reference<std::list<TValue> >::Type
back(std::list<TValue> & list)
{
    return list.back();
}

template <typename TValue>
inline typename Reference<std::list<TValue> const>::Type
back(std::list<TValue> const & list)
{
    return list.back();
}

template <typename TValue>
inline typename Size<std::list<TValue> >::Type
length(std::list<TValue> & list)
{
    return list.size();
}

template <typename TValue>
inline typename Size<std::list<TValue> const>::Type
length(std::list<TValue> const & list)
{
    return list.size();
}

template <typename T, typename T2>
inline void
appendValue(std::list<T> & list, 
			T2 value)
{
	list.push_back(value);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_ADAPT_STD_LIST_H_
