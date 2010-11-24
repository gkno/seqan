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

#ifndef SEQAN_HEADER_SEQUENCE_STREAM_H
#define SEQAN_HEADER_SEQUENCE_STREAM_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec>
inline bool 
_streamEOF(Iter<TContainer, TSpec> const & iter)
{
SEQAN_CHECKPOINT
	return atEnd(iter);
}

//////////////////////////////////////////////////////////////////////////////
 
template <typename TValue, typename TContainer, typename TSpec>
inline ::std::streamsize 
_streamRead(TValue * target,
			Iter<TContainer, TSpec> & source,
			::std::streamsize limit)
{
SEQAN_CHECKPOINT
	if (position(target) + limit > length(container(target)))
		limit = length(container(target)) - position(target);
	Iter<TContainer, TSpec> sourceEnd = source + limit;
	for (; source != sourceEnd; ++source, ++target)
		*target = *source;
	return limit;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec>
inline typename Value<Iter<TContainer, TSpec> >::Type
_streamGet(Iter<TContainer, TSpec> & source)
{
SEQAN_CHECKPOINT
	typename Value<Iter<TContainer, TSpec> >::Type _val = getValue(source);
	goNext(source);
	return _val;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec>
inline typename Value<Iter<TContainer, TSpec> >::Type 
_streamPeek(Iter<TContainer, TSpec> & source)
{
SEQAN_CHECKPOINT
	return getValue(source);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec, typename TChar>
inline void
_streamPut(Iter<TContainer, TSpec> & target,
		   TChar character)
{
SEQAN_CHECKPOINT
	if (atEnd(target))
	{
		typename Container<Iter<TContainer, TSpec> >::Type & container_ = container(target);
		appendValue(container_, character);
		target = begin(container_) + (length(container_) - 1);
	} else
		*target = character;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec>
inline typename Position<Iter<TContainer, TSpec> >::Type
_streamTellG(Iter<TContainer, TSpec> & me)
{
SEQAN_CHECKPOINT
	return position(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec>
inline typename Position<Iter<TContainer, TSpec> >::Type
_streamTellP(Iter<TContainer, TSpec> & me)
{
SEQAN_CHECKPOINT
	return position(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec>
inline void
_streamSeekG(Iter<TContainer, TSpec> & me,
	 typename Position<Iter<TContainer, TSpec> >::Type pos)
{
SEQAN_CHECKPOINT
	me = begin(container(me)) + pos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec>
inline void
_streamSeekP(Iter<TContainer, TSpec> & me,
	 typename Position<Iter<TContainer, TSpec> >::Type pos)
{
SEQAN_CHECKPOINT
	me = begin(container(me)) + pos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec>
inline void
_streamSeek2G(Iter<TContainer, TSpec> & me,
	 int off)
{
SEQAN_CHECKPOINT
	me = begin(container(me)) + (position(me) + off);
}


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
