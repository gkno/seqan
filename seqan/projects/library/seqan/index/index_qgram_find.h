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


//____________________________________________________________________________

	friend inline typename _Parameter<TIndex>::Type 
	host(Finder & me)
	{
SEQAN_CHECKPOINT
		return value(me.index);
	}

	friend inline typename _Parameter<TIndex>::Type 
	host(Finder const & me)
	{
SEQAN_CHECKPOINT
		return value(me.index);
	}

	friend inline typename _Parameter<TIndex>::Type 
	container(Finder & me)
	{
SEQAN_CHECKPOINT
		return value(me.index);
	}

	friend inline typename _Parameter<TIndex>::Type 
	container(Finder const & me)
	{
SEQAN_CHECKPOINT
		return value(me.index);
	}

//____________________________________________________________________________

	friend inline void
	setHost(Finder & me, typename _Parameter<TIndex>::Type container_)
	{
SEQAN_CHECKPOINT
		me.index = container;
	}

	friend inline void
	setContainer(Finder & me, typename _Parameter<TIndex>::Type container_)
	{
SEQAN_CHECKPOINT
		me.index = container;
	}

//____________________________________________________________________________

	friend inline TIterator &
	hostIterator(Finder & me)
	{
SEQAN_CHECKPOINT
		return me.data_iterator;
	}

	friend inline TIterator const &
	hostIterator(Finder const & me)
	{
SEQAN_CHECKPOINT
		return me.data_iterator;
	}


//____________________________________________________________________________

	friend inline bool
	empty(Finder & me)
	{
SEQAN_CHECKPOINT
		return me.range.i1 == me.range.i2;
	}

	friend inline void
	clear(Finder & me)
	{
SEQAN_CHECKPOINT
		me.range.i1 = me.range.i2 = TIterator();
	}

//____________________________________________________________________________

	friend inline bool
	atBegin(Finder & me)
	{
SEQAN_CHECKPOINT
		return (empty(me) || hostIterator(me) == me.range.i1);
	}

	friend inline bool
	atEnd(Finder & me)
	{
SEQAN_CHECKPOINT
		return (empty(me) || hostIterator(me) == me.range.i2);
	}

//____________________________________________________________________________

	friend inline void
	goBegin(Finder & me)
	{
SEQAN_CHECKPOINT
		hostIterator(me) = me.range.i1;
	}

	friend inline void
	goEnd(Finder & me)
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

	friend inline typename Position<Finder>::Type
	position(Finder & me)
	{
SEQAN_CHECKPOINT
		if (empty(me)) return 0;
		return *me.data_iterator;
	}

	friend inline typename Position<Finder>::Type
	position(Finder const & me)
	{
SEQAN_CHECKPOINT
		if (empty(me)) return 0;
		return hostIterator(me) - begin(container(me), Rooted());
	}

};


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
