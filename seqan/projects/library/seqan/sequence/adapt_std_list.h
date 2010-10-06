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

/**
.Adaption."std::list"
..summary:Adaption for STL list objects.
*/
    
// ===========================================================================
// Metafunctions
// ===========================================================================

///.Metafunction.IsContiguous.param.T.type:Adaption.std::list
template <typename TChar, typename TAlloc>
struct IsContiguous< ::std::basic_string<TChar, TAlloc> >
{
    enum { VALUE = false };
};

template <typename  TChar, typename TAlloc>
struct IsContiguous< ::std::basic_string<TChar, TAlloc> const>
        : IsContiguous< ::std::basic_string<TChar, TAlloc> > {};

///.Metafunction.Value.param.T.type:Adaption.std::list
template <typename TValue, typename TAlloc>
struct Value< ::std::list<TValue, TAlloc> >
{
	typedef typename ::std::list<TValue, TAlloc>::value_type Type;
};

template <typename TValue, typename TAlloc>
struct Value< ::std::list<TValue, TAlloc> const>
        : Value< ::std::list<TValue, TAlloc> > {};

///.Metafunction.GetValue.param.T.type:Adaption.std::list
template <typename TValue, typename TAlloc>
struct GetValue< ::std::list<TValue, TAlloc> >
{
	typedef typename ::std::list<TValue, TAlloc>::reference Type;
};

template <typename TValue, typename TAlloc>
struct GetValue< ::std::list<TValue, TAlloc> const>
{
	typedef typename ::std::list<TValue, TAlloc>::const_reference Type;
};

///.Metafunction.Iterator.param.T.type:Adaption.std::list
template <typename TValue, typename TAlloc>
struct Iterator< ::std::list<TValue, TAlloc>, Rooted>
{
	typedef ::std::list<TValue, TAlloc> _TString;
	typedef Iter<_TString, StdIteratorAdaptor> _TIterator;
	typedef Iter<_TString, AdaptorIterator<_TIterator> > Type;
};

template <typename TValue, typename TAlloc>
struct Iterator< ::std::list<TValue, TAlloc> const, Rooted>
{
	typedef ::std::list<TValue, TAlloc> const _TString;
	typedef Iter<_TString, StdIteratorAdaptor> _TIterator;
	typedef Iter<_TString, AdaptorIterator<_TIterator> > Type;
};

template <typename TValue, typename TAlloc>
struct Iterator< ::std::list<TValue, TAlloc>, Standard>
{
	typedef Iter< ::std::list<TValue, TAlloc>, StdIteratorAdaptor> Type;
};

template <typename TValue, typename TAlloc>
struct Iterator< ::std::list<TValue, TAlloc> const, Standard>
{
	typedef Iter< ::std::list<TValue, TAlloc> const, StdIteratorAdaptor> Type;
};

///.Metafunction.Position.param.T.type:Adaption.std::list
template <typename TValue, typename TAlloc>
struct Position< ::std::list<TValue, TAlloc> >
{
	typedef typename ::std::list<TValue, TAlloc>::size_type Type;
};

template <typename TValue, typename TAlloc>
struct Position< ::std::list<TValue, TAlloc> const>
        : Position< ::std::list<TValue, TAlloc> > {};

///.Metafunction.Size.param.T.type:Adaption.std::list
template <typename TValue, typename TAlloc>
struct Size< ::std::list<TValue, TAlloc> >
{
	typedef typename ::std::list<TValue, TAlloc>::size_type Type;
};

template <typename TValue, typename TAlloc>
struct Size< ::std::list<TValue, TAlloc> const>
        : Size< ::std::list<TValue, TAlloc> > {};

/**
.Metafunction.StdContainerIterator
..summary:Returns type of the STL container iterator.
..signature:StdContainerIterator<T>::Type
..param.T.type:Adaption.std::list
*/
template <typename TValue, typename TAlloc>
struct StdContainerIterator< ::std::list<TValue, TAlloc> >
{
	typedef ::std::list<TValue, TAlloc> _TContainer;
	typedef typename _TContainer::iterator Type;
};

template <typename TValue, typename TAlloc>
struct StdContainerIterator< ::std::list<TValue, TAlloc> const>
{
	typedef ::std::list<TValue, TAlloc> _TContainer;
	typedef typename _TContainer::const_iterator Type;
};

// ===========================================================================
// Functions
// ===========================================================================

///.Function.begin.param.object.type:Adaption.std::list
template<typename TValue>
inline
typename Iterator<std::list<TValue>, Standard>::Type
begin(std::list<TValue> & list,
      Standard const &)
{
    SEQAN_CHECKPOINT;
	return list.begin();
}

template <typename TValue>
inline
typename Iterator<std::list<TValue> const, Standard>::Type
begin(std::list<TValue> const & list,
      Standard const &)
{
    SEQAN_CHECKPOINT;
	return list.begin();
}

///.Function.end.param.object.type:Adaption.std::list
template<typename TValue>
inline
typename Iterator<std::list<TValue>, Standard>::Type
end(std::list<TValue> & list,
      Standard const &)
{
    SEQAN_CHECKPOINT;
	return list.end();
}

template <typename TValue>
inline
typename Iterator<std::list<TValue> const, Standard>::Type
end(std::list<TValue> const & list,
      Standard const &)
{
    SEQAN_CHECKPOINT;
	return list.end();
}

///.Function.Container#front.param.container.type:Adaption.std::list
template <typename TValue>
inline typename Reference<std::list<TValue> >::Type
front(std::list<TValue> & list)
{
    SEQAN_CHECKPOINT;
    return list.front();
}

template <typename TValue>
inline typename Reference<std::list<TValue> const>::Type
front(std::list<TValue> const & list)
{
    SEQAN_CHECKPOINT;
    return list.front();
}

///.Function.back.param.container.type:Adaption.std::list
template <typename TValue>
inline typename Reference<std::list<TValue> >::Type
back(std::list<TValue> & list)
{
    SEQAN_CHECKPOINT;
    return list.back();
}

template <typename TValue>
inline typename Reference<std::list<TValue> const>::Type
back(std::list<TValue> const & list)
{
    SEQAN_CHECKPOINT;
    return list.back();
}

///.Function.length.param.object.type:Adaption.std::list
template <typename TValue>
inline typename Size<std::list<TValue> >::Type
length(std::list<TValue> & list)
{
    SEQAN_CHECKPOINT;
    return list.size();
}

template <typename TValue>
inline typename Size<std::list<TValue> const>::Type
length(std::list<TValue> const & list)
{
    SEQAN_CHECKPOINT;
    return list.size();
}

/**
.Function.prependValue:
..summary:Prepend a value to a container.
..signature:prependValue(container, value)
..param.container:The container to prepend to.
...type:Adaption.std::list
..param.value:The value to prepend to the container.
*/
template <typename T, typename T2>
inline void
prependValue(std::list<T> & list, 
             T2 value)
{
    SEQAN_CHECKPOINT;
	list.push_front(value);
}

///.Function.appendValue.param.target.type:Adaption.std::list
template <typename T, typename T2>
inline void
appendValue(std::list<T> & list, 
			T2 value)
{
    SEQAN_CHECKPOINT;
	list.push_back(value);
}

///.Function.clear.param.object.type:Adaption.std::list
template <typename T>
inline void
clear(std::list<T> & list)
{
    SEQAN_CHECKPOINT;
    list.clear();
}

/**
.Function.reverse:
..cat:Containers
..signature:reverse(container)
..param.container:The container whose elements to reverse.
...type:Concept.Container
...type:Adaption.std::list
 */
template <typename TContainer>
void reverse(TContainer & container)
{    
    SEQAN_CHECKPOINT;
    std::reverse(begin(container, Standard()), end(container, Standard()));
}

template <typename TElement>
void reverse(std::list<TElement> & container)
{
    SEQAN_CHECKPOINT;
    container.reverse();
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_ADAPT_STD_LIST_H_
