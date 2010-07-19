/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2010

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
 ============================================================================
  A wrapper for an iterator that dereferences the wrapped iterator twice
  on dereferentiation.
 ==========================================================================*/

#ifndef SEQAN_SEEDS2_BASIC_ITER_INDIRECT_H_
#define SEQAN_SEEDS2_BASIC_ITER_INDIRECT_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

template <typename TWrappedIter>
struct Indirect {};

// TODO(holtgrew): Make this an Iter specialization?
//
// Iterator that dereferences twice on derefentiation.
template <typename TContainer, typename TWrappedIter>
class Iter<TContainer, Indirect<TWrappedIter> >
{
public:
    TWrappedIter _wrappedIter;
    
    Iter()
    { SEQAN_CHECKPOINT; }

    Iter(Iter const & other)
            : _wrappedIter(other._wrappedIter)
    { SEQAN_CHECKPOINT; }

    Iter(TWrappedIter const & wrappedIter)
            : _wrappedIter(wrappedIter)
    { SEQAN_CHECKPOINT; }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <typename TContainer, typename TWrappedIter>
struct Reference<Iter<TContainer, Indirect<TWrappedIter> > >
{
    typedef typename Reference<TContainer>::Type Type;
};

// ===========================================================================
// Functions
// ===========================================================================

template <typename TContainer, typename TWrappedIter>
inline bool
operator==(Iter<TContainer, Indirect<TWrappedIter> > const & a, Iter<TContainer, Indirect<TWrappedIter> > const & b)
{
    SEQAN_CHECKPOINT;
    return a._wrappedIter == b._wrappedIter;
}


template <typename TContainer, typename TWrappedIter>
inline bool
operator!=(Iter<TContainer, Indirect<TWrappedIter> > const & a, Iter<TContainer, Indirect<TWrappedIter> > const & b)
{
    SEQAN_CHECKPOINT;
    return a._wrappedIter != b._wrappedIter;
}

template <typename TContainer, typename TWrappedIter, typename TDiff>
Iter<TContainer, Indirect<TWrappedIter> >
operator+(Iter<TContainer, Indirect<TWrappedIter> > const & it, TDiff const diff)
{
    SEQAN_CHECKPOINT;
    typedef Iter<TContainer, Indirect<TWrappedIter> > TIter;
    return TIter(it._wrappedIter + diff);
}


template <typename TContainer, typename TWrappedIter>
inline Iter<TContainer, Indirect<TWrappedIter> > &
operator++(Iter<TContainer, Indirect<TWrappedIter> > & iter)
{
    SEQAN_CHECKPOINT;
    ++iter._wrappedIter;
    return iter;
}


template <typename TContainer, typename TWrappedIter>
inline Iter<TContainer, Indirect<TWrappedIter> >
operator++(Iter<TContainer, Indirect<TWrappedIter> > & iter, int /*postfix*/)
{
    SEQAN_CHECKPOINT;
    typedef Iter<TContainer, Indirect<TWrappedIter> > TIter;
    TIter tmp(iter);
    ++iter;
    return tmp;
}


template <typename TContainer, typename TWrappedIter>
inline typename Value<TWrappedIter>::Type &
operator*(Iter<TContainer, Indirect<TWrappedIter> > & iter)
{
    SEQAN_CHECKPOINT;
    return **iter._wrappedIter;
}


template <typename TContainer, typename TWrappedIter>
inline typename Value<TWrappedIter>::Type &
operator*(Iter<TContainer, Indirect<TWrappedIter> > const & iter)
{
    SEQAN_CHECKPOINT;
    return **iter._wrappedIter;
}


template <typename TContainer, typename TWrappedIter>
inline typename Reference<Iter<TContainer, Indirect<TWrappedIter> > >::Type
value(Iter<TContainer, Indirect<TWrappedIter> > & iter)
{
    SEQAN_CHECKPOINT;
    return **iter._wrappedIter;
}


template <typename TContainer, typename TWrappedIter>
inline typename Reference<Iter<TContainer, Indirect<TWrappedIter> > >::Type
value(Iter<TContainer, Indirect<TWrappedIter> > const & iter)
{
    SEQAN_CHECKPOINT;
    return **iter._wrappedIter;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS2_BASIC_ITER_INDIRECT_H_
