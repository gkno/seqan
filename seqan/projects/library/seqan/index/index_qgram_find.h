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

#ifndef SEQAN_HEADER_INDEX_QGRAM_FIND_H
#define SEQAN_HEADER_INDEX_QGRAM_FIND_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// QGram finders

struct _Finder_QGramHashLookup; //Finder that simply looks up the q-gram in the hash table

typedef Tag<_Finder_QGramHashLookup> const HashLookup;

//____________________________________________________________________________


template < typename TText, typename TSpecShape >
struct DefaultFinder<Index<TText, Index_QGram<TSpecShape> > > 
{
    typedef HashLookup Type;
};

//////////////////////////////////////////////////////////////////////////////
// Position Spec for Index Finder
template <typename TIndex>
struct Position< Finder<TIndex, HashLookup> > {
	typedef typename SAValue<TIndex>::Type Type;
};


//////////////////////////////////////////////////////////////////////////////
/**
.Spec.HashLookup:
..summary:Finding q-grams in index using a hash table.
..general:Class.Finder
..cat:Index
..signature:Finder<TIndex, HashLookup>
..param.TIndex:The index type.
...type:Spec.Index_QGram
*/

template <typename TIndex>
class Finder <TIndex, HashLookup>
{
public:
	typedef typename Iterator< typename Fibre<TIndex, QGram_SA>::Type const, Standard >::Type TIterator;

	Holder<TIndex>	index;
	Pair<TIterator>	range;
	TIterator		data_iterator;

	Finder() 
	{
		clear(*this);
	}
	Finder(TIndex &_index): index(_index) 
	{
		clear(*this);
	}
	Finder(TIndex const &_index): index(_index)
	{
		clear(*this);
	}
};

//____________________________________________________________________________
template <typename TIndex>
inline typename _Parameter<TIndex>::Type 
host(Finder<TIndex, HashLookup> & me)
{
SEQAN_CHECKPOINT
	return value(me.index);
}

template <typename TIndex>
inline typename _Parameter<TIndex>::Type 
host(Finder<TIndex, HashLookup> const & me)
{
SEQAN_CHECKPOINT
	return value(me.index);
}

template <typename TIndex>
inline typename _Parameter<TIndex>::Type 
container(Finder<TIndex, HashLookup> & me)
{
SEQAN_CHECKPOINT
	return value(me.index);
}

template <typename TIndex>
inline typename _Parameter<TIndex>::Type 
container(Finder<TIndex, HashLookup> const & me)
{
SEQAN_CHECKPOINT
	return value(me.index);
}

//____________________________________________________________________________

template <typename TIndex>
inline void
setHost(Finder<TIndex, HashLookup> & me, typename _Parameter<TIndex>::Type container_)
{
SEQAN_CHECKPOINT
	me.index = container;
}

template <typename TIndex>
inline void
setContainer(Finder<TIndex, HashLookup> & me, typename _Parameter<TIndex>::Type container_)
{
SEQAN_CHECKPOINT
	me.index = container;
}

//____________________________________________________________________________

template <typename TIndex>
inline typename Iterator< typename Fibre<TIndex, QGram_SA>::Type const, Standard >::Type &
hostIterator(Finder<TIndex, HashLookup> & me)
{
SEQAN_CHECKPOINT
	return me.data_iterator;
}

template <typename TIndex>
inline typename Iterator< typename Fibre<TIndex, QGram_SA>::Type const, Standard >::Type const &
hostIterator(Finder<TIndex, HashLookup> const & me)
{
SEQAN_CHECKPOINT
	return me.data_iterator;
}


//____________________________________________________________________________

template <typename TIndex>
inline bool
empty(Finder<TIndex, HashLookup> & me)
{
SEQAN_CHECKPOINT
	return me.range.i1 == me.range.i2;
}

template <typename TIndex>
inline void
clear(Finder<TIndex, HashLookup> & me)
{
SEQAN_CHECKPOINT
	typedef typename Iterator< typename Fibre<TIndex, QGram_SA>::Type const, Standard >::Type TIterator;
	me.range.i1 = me.range.i2 = TIterator();
}

//____________________________________________________________________________

template <typename TIndex>
inline bool
atBegin(Finder<TIndex, HashLookup> & me)
{
SEQAN_CHECKPOINT
	return (empty(me) || hostIterator(me) == me.range.i1);
}

template <typename TIndex>
inline bool
atEnd(Finder<TIndex, HashLookup> & me)
{
SEQAN_CHECKPOINT
	return (empty(me) || hostIterator(me) == me.range.i2);
}

//____________________________________________________________________________

template <typename TIndex>
inline void
goBegin(Finder<TIndex, HashLookup> & me)
{
SEQAN_CHECKPOINT
	hostIterator(me) = me.range.i1;
}

template <typename TIndex>
inline void
goEnd(Finder<TIndex, HashLookup> & me)
{
SEQAN_CHECKPOINT
	hostIterator(me) = me.range.i2;
}

//____________________________________________________________________________
/*
template <typename TPosition>
friend inline void 
setPosition(Finder & me, TPosition pos_)
{
SEQAN_CHECKPOINT
	hostIterator(me) = me.range.i1 + pos_;
}
*/
//____________________________________________________________________________

template <typename TIndex>
inline typename Position<Finder<TIndex, HashLookup> >::Type
position(Finder<TIndex, HashLookup> & me)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(!empty(me))
	return *me.data_iterator;
}

template <typename TIndex>
inline typename Position<Finder<TIndex, HashLookup> >::Type
position(Finder<TIndex, HashLookup> const & me)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(!empty(me))
	return hostIterator(me) - begin(container(me), Rooted());
}


//////////////////////////////////////////////////////////////////////////////
// find implementation

template < typename TFinder, typename TPattern >
inline void _findFirstQGramIndex(TFinder &finder,
								 TPattern const &pattern,
								 HashLookup)
{
	typedef typename Haystack<TFinder>::Type TIndex;

	typedef typename Fibre<TIndex, QGram_SA>::Type TSA;

	typedef typename Fibre<TIndex, QGram_Shape>::Type TShape;
	typedef typename Value<TShape>::Type TShapeValue;

	typedef typename Fibre<TIndex, QGram_Dir>::Type THash;

	typedef typename Iterator<TPattern const, Standard>::Type TPatternIterator;

	TIndex & _index = haystack(finder);
	indexRequire(_index, QGram_SA());
	TSA const & _sa = indexSA(_index);

	//indexRequire(_index, QGram_Dir());

	TShape & _shape = indexShape(_index);

	TPatternIterator _it = begin(pattern, Standard());
	TShapeValue hash_value = hash(_shape, _it);

	THash & _dir = indexDir(_index);
	finder.range.i1 = begin(_sa, Standard()) + _dir[hash_value];
	finder.range.i2 = begin(_sa, Standard()) + _dir[hash_value+1];
}

//////////////////////////////////////////////////////////////////////////////
// find

template < typename TText, typename TSpecShape, typename TSpecFinder, typename TPattern >
inline bool 
find(Finder<Index<TText, Index_QGram<TSpecShape> >, TSpecFinder> & finder,
	 TPattern const &pattern)
{
	if (empty(finder)) 
	{
		_findFirstQGramIndex(finder, needle(pattern), TSpecFinder());
		hostIterator(finder) = finder.range.i1;
	} 
	else
	{
		++hostIterator(finder);
	}	
	return !atEnd(finder);
}

template < typename TText, typename TSpecShape, typename TSpecFinder >
inline bool 
find(Finder<Index<TText, Index_QGram<TSpecShape> >, TSpecFinder> & finder)
{
	if (empty(finder)) return false;
	++hostIterator(finder);
	return !atEnd(finder);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_
