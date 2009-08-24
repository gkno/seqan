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

#ifndef SEQAN_HEADER_STORE_ALIGN_H
#define SEQAN_HEADER_STORE_ALIGN_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Aligned Read Store
//////////////////////////////////////////////////////////////////////////////

template <typename _TPos, typename _TGapAnchor, typename _TSpec = void>
struct AlignedReadStoreElement
{
	typedef typename Id<AlignedReadStoreElement>::Type	TId;
	typedef _TPos										TPos;
	typedef _TGapAnchor									TGapAnchor;
	typedef _TSpec										TSpec;
	typedef String<TGapAnchor>							TGapAnchors;

	static const TId INVALID_ID;
	
	TId			id;
	TId			readId;
	TId			contigId;
	TId			pairMatchId;	// unique id. for multiple mate-pair matches (not matePairId)
	TPos		beginPos;		// begin position of the gapped sequence in gapped contig sequence
	TPos		endPos;			// end position of ..., for reverse aligned reads holds end < begin
	TGapAnchors	gaps;

	AlignedReadStoreElement() : id(INVALID_ID), readId(INVALID_ID), contigId(INVALID_ID), pairMatchId(INVALID_ID), beginPos(0), endPos(0) {}
};

template <typename TScore, typename TSpec = void>
struct AlignQualityStoreElement
{
	TScore				pairScore;		// score of the mate-pair alignment (this read is part of)
	TScore				score;			// score of the single read alignment
	unsigned char		errors;			// absolute number of errors (Hamming or edit distance)
};


//////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TGapAnchor, typename TSpec>
const typename Id<AlignedReadStoreElement<TPos, TGapAnchor, TSpec> >::Type
AlignedReadStoreElement<TPos, TGapAnchor, TSpec>::INVALID_ID = SupremumValue<typename Id<AlignedReadStoreElement<TPos, TGapAnchor, TSpec> >::Type>::VALUE;


//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Sorting tags
//////////////////////////////////////////////////////////////////////////////


struct _SortContigId;
typedef Tag<_SortContigId> const SortContigId;


struct _SortId;
typedef Tag<_SortId> const SortId;

struct _SortBeginPos;
typedef Tag<_SortBeginPos> const SortBeginPos;

struct _SortEndPos;
typedef Tag<_SortEndPos> const SortEndPos;

struct _SortPairMatchId;
typedef Tag<_SortPairMatchId> const SortPairMatchId;

struct _SortReadId;
typedef Tag<_SortReadId> const SortReadId;


//////////////////////////////////////////////////////////////////////////////
// Sorting functors
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


template <typename TAlignedRead, typename TTag>
struct _LessAlignedRead;

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortId> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return (a1.id) < (a2.id);
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortContigId> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return a1.contigId < a2.contigId;
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortBeginPos> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return _min(a1.beginPos, a1.endPos) < _min(a2.beginPos, a2.endPos);
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortEndPos> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return _max(a1.beginPos, a1.endPos) < _max(a2.beginPos, a2.endPos);
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortPairMatchId> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return a1.pairMatchId < a2.pairMatchId;
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortReadId> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return a1.readId < a2.readId;
	}
};

//////////////////////////////////////////////////////////////////////////////
// Sorting function
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSortSpec>
inline void
sortAlignedReads(TAlign& alignStore, Tag<TSortSpec> const &) 
{
	std::stable_sort(
		begin(alignStore, Standard() ), 
		end(alignStore, Standard() ), 
		_LessAlignedRead<typename Value<TAlign>::Type, Tag<TSortSpec> const>() );
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSortSpec>
inline void
sortAlignedReads(TAlign const & alignStore, Tag<TSortSpec> const &) 
{
	std::stable_sort(
		begin(const_cast<TAlign&>(alignStore), Standard() ), 
		end(const_cast<TAlign&>(alignStore), Standard() ), 
		_LessAlignedRead<typename Value<TAlign>::Type, Tag<TSortSpec> const>() );
}

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign & alignStore, TFunctorLess const &less) 
{
	std::stable_sort(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		less);
}

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign const & alignStore, TFunctorLess const &less) 
{
	std::stable_sort(
		begin(const_cast<TAlign&>(alignStore), Standard()), 
		end(const_cast<TAlign&>(alignStore), Standard()), 
		less);
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortContigId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.contigId = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortContigId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortContigId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.contigId = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortContigId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortBeginPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortBeginPos const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.endPos = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortEndPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.endPos = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortEndPos const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortPairMatchId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchId = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortPairMatchId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortPairMatchId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchId = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortPairMatchId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortReadId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.readId = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortReadId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortReadId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.readId = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortReadId const>() );
}

//////////////////////////////////////////////////////////////////////////////


template <typename TGapAnchors = String<GapAnchor<unsigned> > >
struct AnchorGaps;

template <typename TSource, typename TGapAnchors>
class Gaps<TSource, AnchorGaps<TGapAnchors> >
{
public:
	typedef typename Value<TGapAnchors>::Type TGapAnchor;
	typedef typename Position<TGapAnchor>::Type TViewPosition;

public:
	Holder<TSource>		data_source;
	Holder<TGapAnchors>	data_gaps;
	TViewPosition		data_viewBeginPos;
	TViewPosition		data_viewEndPos;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TGapAnchors>
inline Holder<TSource> &
_dataSource(Gaps<TSource, AnchorGaps<TGapAnchors> > & me)
{
SEQAN_CHECKPOINT
	return me.data_source;
}
template <typename TSource, typename TGapAnchors>
inline Holder<TSource> const &
_dataSource(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me)
{
SEQAN_CHECKPOINT
	return me.data_source;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TGapAnchors>
inline TGapAnchors &
_dataAnchors(Gaps<TSource, AnchorGaps<TGapAnchors> > & me)
{
	return value(me.data_gaps);
}

template <typename TSource, typename TGapAnchors>
inline TGapAnchors const &
_dataAnchors(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me)
{
	return value(me.data_gaps);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TGapAnchors, typename TPosition>
inline TPosition
positionGapToSeq(Gaps<TSource, AnchorGaps<TGapAnchors> > & me, TPosition pos)
{
	typedef typename Value<TGapAnchors>::Type TGapAnchor;

	typename Iterator<TGapAnchors, Standard>::Type it;
	TGapAnchors &anchors = _dataAnchors(me);
	TGapAnchor	prevAnchor, nextAnchor;
	TPosition	seqPos;
	TPosition	seqLength = length(source(me));

	if (!empty(anchors))
	{
		// differentiate between 0 (first alignment pos) and the rest because of soft-clipping at the beginning/end
		if (pos == 0)
			it = upperBoundGapAnchor(anchors, 0, SortGapPos());
		else
			it = lowerBoundGapAnchor(anchors, pos, SortGapPos());
		if (it != end(anchors, Standard()))
		{
			if ((*it).gapPos == pos && (*it).seqPos != seqLength)
			{
				prevAnchor = *it;
				++it;
			} else
			{
				if (it != begin(anchors, Standard()))
					prevAnchor = *(it - 1);
				else
					prevAnchor = TGapAnchor(0, 0);
			}
		} else
		{
			if (!empty(anchors))
				prevAnchor = back(anchors);
			else
				prevAnchor = TGapAnchor(0, 0);
		}
		if (it != end(anchors, Standard()))
			nextAnchor = *it;
		else
		{
			nextAnchor.gapPos = prevAnchor.gapPos + (seqLength - prevAnchor.seqPos);
			nextAnchor.seqPos = seqLength;
		}
	} else 
	{
		prevAnchor = TGapAnchor(0, 0);
		nextAnchor = TGapAnchor(seqLength, seqLength);
	}
	seqPos = prevAnchor.seqPos + (pos - prevAnchor.gapPos);
	if (seqPos > nextAnchor.seqPos)
		seqPos = nextAnchor.seqPos;
	return seqPos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource_, typename TGapAnchors_>
class Iter<Gaps<TSource_, AnchorGaps<TGapAnchors_> >, GapsIterator<AnchorGaps<TGapAnchors_> > >
{
public:
	typedef TSource_										TSource;
	typedef TGapAnchors_									TGapAnchors;

	typedef Gaps<TSource, AnchorGaps<TGapAnchors> >			TGaps;
	typedef typename Value<TGapAnchors>::Type				TGapAnchor;
	typedef typename Size<TGapAnchor>::Type					TGapSize;
	typedef typename Iterator<TGapAnchors, Standard>::Type	TAnchorIter;

	TGaps *			data_container;							//the gaps object
	TGapSize		seqLength;
	TGapAnchor		current;
	TGapAnchor		prevAnchor;
	TGapAnchor		nextAnchor;
	TGapAnchor		viewBegin;
	TGapAnchor		viewEnd;
	TAnchorIter		anchorIter;

public:
	Iter() 
	{
SEQAN_CHECKPOINT
		data_container = NULL;
		seqLength = 0;
	}
/*	Iter(Iter const & other_):
		data_container(other_.data_container),
		seqLength(other_.seqLength),
		current(other_.current),
		prevAnchor(other_.prevAnchor),
		nextAnchor(other_.nextAnchor),
		anchorIter(other_.anchorIter)
	{
SEQAN_CHECKPOINT
	}
*/	Iter(TGaps & container_):
		data_container(&container_),
		seqLength(length(source(*data_container)))
	{
SEQAN_CHECKPOINT
		seqLength = length(source(*data_container));
		_goTo_gapAnchorIterator(*this, 0);
		viewBegin = current;
		viewEnd.gapPos = seqLength;
		if (!empty(_dataAnchors(*data_container)))
		{
			viewEnd = back(_dataAnchors(*data_container));
			viewEnd.gapPos += seqLength - viewEnd.seqPos;
		}
		viewEnd.seqPos = positionGapToSeq(*data_container, viewEnd.gapPos);
	}
	Iter(TGaps & container_, TGapSize position):
		data_container(&container_),
		seqLength(length(source(*data_container)))
	{
SEQAN_CHECKPOINT
		_goTo_gapAnchorIterator(*this, position);
		viewBegin.gapPos = 0;
		viewEnd.gapPos = seqLength;
		if (!empty(_dataAnchors(*data_container)))
		{
			viewEnd = back(_dataAnchors(*data_container));
			viewEnd.gapPos += seqLength - viewEnd.seqPos;
		}
		viewBegin.seqPos = positionGapToSeq(*data_container, viewBegin.gapPos);
		viewEnd.seqPos = positionGapToSeq(*data_container, viewEnd.gapPos);
	}
	~Iter()
	{
SEQAN_CHECKPOINT
	}

	Iter const & operator = (Iter const & other_)
	{
SEQAN_CHECKPOINT
		data_container = other_.data_container;
		seqLength = other_.seqLength;
		current = other_.current;
		prevAnchor = other_.prevAnchor;
		nextAnchor = other_.nextAnchor;
		anchorIter = other_.anchorIter;
		viewBegin = other_.viewBegin;
		viewEnd = other_.viewEnd;
		return *this;
	}
};

//____________________________________________________________________________

template <typename TSource, typename TGapAnchors>
inline bool 
isGap(Iter<Gaps<TSource, AnchorGaps<TGapAnchors> >, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
	return me.current.seqPos == me.nextAnchor.seqPos;
}

template <typename TSource, typename TGapAnchors>
inline bool 
isClipped(Iter<Gaps<TSource, AnchorGaps<TGapAnchors> >, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
	return me.current.gapPos == me.nextAnchor.gapPos;
}

//____________________________________________________________________________

template <typename TSource, typename TGapAnchors>
inline bool 
atBegin(Iter<Gaps<TSource, AnchorGaps<TGapAnchors> >, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
//	return me.current.seqPos == 0 && me.current.gapPos == 0;
	return me.current == me.viewBegin;
}

template <typename TSource, typename TGapAnchors>
inline bool 
atBegin(Iter<Gaps<TSource, AnchorGaps<TGapAnchors> >, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
//	return me.current.seqPos == 0 && me.current.gapPos == 0;
	return me.current == me.viewBegin;
}

template <typename TSource, typename TGapAnchors>
inline bool 
atEnd(Iter<Gaps<TSource, AnchorGaps<TGapAnchors> >, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
//	return me.current == me.nextAnchor;
	return me.current == me.viewEnd;
}

template <typename TSource, typename TGapAnchors>
inline bool 
atEnd(Iter<Gaps<TSource, AnchorGaps<TGapAnchors> >, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
//	return me.current == me.nextAnchor;
	return me.current == me.viewEnd;
}

//____________________________________________________________________________

template <typename T>
inline void 
_goNext_gapAnchorIterator(T & me)
{
	if (me.current.seqPos < me.nextAnchor.seqPos)
		++me.current.seqPos;
	if (me.current.gapPos < me.nextAnchor.gapPos)
		++me.current.gapPos;
	if (me.current == me.nextAnchor)
	{
		if (!empty(_dataAnchors(*me.data_container)) && me.anchorIter != end(_dataAnchors(*me.data_container), Standard()))
		{
			++me.anchorIter;
			if (me.anchorIter != end(_dataAnchors(*me.data_container), Standard()))
			{
				me.prevAnchor = me.nextAnchor;
				me.nextAnchor = *me.anchorIter;
				return;
			}
		}
		if (me.seqLength != me.nextAnchor.seqPos)
		{
			me.prevAnchor = me.nextAnchor;
			me.nextAnchor.gapPos += me.seqLength - me.nextAnchor.seqPos;
			me.nextAnchor.seqPos = me.seqLength;
		}
	}
}

template <typename T>
inline void 
_goPrevious_gapAnchorIterator(T & me)
{
	if (me.current == me.prevAnchor)
	{
		if (me.anchorIter != begin(_dataAnchors(*me.data_container), Standard()))
		{
			if (me.anchorIter == end(_dataAnchors(*me.data_container), Standard()))
			{
				me.anchorIter = begin(_dataAnchors(*me.data_container), Standard());
				unsigned len = length(_dataAnchors(*me.data_container));
				if (len)
				{
					me.anchorIter += len - 1;
					if (len > 1 && me.seqLength == (*me.anchorIter).seqPos)
						--me.anchorIter;
				}
			} else
				--me.anchorIter;

			if (me.anchorIter != begin(_dataAnchors(*me.data_container), Standard()))
			{
				me.nextAnchor = me.prevAnchor;
				me.prevAnchor = *(me.anchorIter - 1);
			} else
			{
				me.nextAnchor = me.prevAnchor;
				me.prevAnchor.gapPos = 0;
				me.prevAnchor.seqPos = 0;
			}
		}
	}
	if (me.current.seqPos + me.prevAnchor.gapPos > me.prevAnchor.seqPos + me.current.gapPos)
		--me.current.seqPos;
	else if (me.current.seqPos + me.prevAnchor.gapPos < me.prevAnchor.seqPos + me.current.gapPos)
		--me.current.gapPos;
	else
	{
		--me.current.seqPos;
		--me.current.gapPos;
	}
}

template <typename T, typename TPos>
inline void 
_goTo_gapAnchorIterator(T & me, TPos pos)
{
	typedef typename T::TGapAnchors	TGapAnchors;
	typedef typename T::TGapAnchor	TGapAnchor;

	TGapAnchors & anchors = _dataAnchors(*me.data_container);

	if (!empty(anchors))
	{
		// differentiate between 0 (first alignment pos) and the rest because of soft-clipping at the beginning/end
		if (pos == 0)
			me.anchorIter = upperBoundGapAnchor(anchors, 0, SortGapPos());
		else
			me.anchorIter = lowerBoundGapAnchor(anchors, pos, SortGapPos());
		if (me.anchorIter != end(anchors, Standard()))
		{
			if ((*me.anchorIter).gapPos == pos && (*me.anchorIter).seqPos != me.seqLength)
			{
				me.prevAnchor = *me.anchorIter;
				++me.anchorIter;
			} else
			{
				if (me.anchorIter != begin(anchors, Standard()))
					me.prevAnchor = *(me.anchorIter - 1);
				else
					me.prevAnchor = TGapAnchor(0, 0);
			}
		} else
		{
			if (!empty(anchors))
				me.prevAnchor = back(anchors);
			else
				me.prevAnchor = TGapAnchor(0, 0);
		}
		if (me.anchorIter != end(anchors, Standard()))
			me.nextAnchor = *me.anchorIter;
		else
		{
			me.nextAnchor.gapPos = me.prevAnchor.gapPos + (me.seqLength - me.prevAnchor.seqPos);
			me.nextAnchor.seqPos = me.seqLength;
		}
	} else 
	{
		me.prevAnchor = TGapAnchor(0, 0);
		me.nextAnchor = TGapAnchor(me.seqLength, me.seqLength);
	}
	me.current.seqPos = me.prevAnchor.seqPos + (pos - me.prevAnchor.gapPos);
	if (me.current.seqPos > me.nextAnchor.seqPos)
		me.current.seqPos = me.nextAnchor.seqPos;
	me.current.gapPos = pos;
}

//____________________________________________________________________________

template <typename TSource, typename TGapAnchors>
inline void
goNext(Iter<Gaps<TSource, AnchorGaps<TGapAnchors> >, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
	_goNext_gapAnchorIterator(me);
}

template <typename TSource, typename TGapAnchors>
inline void
goPrevious(Iter<Gaps<TSource, AnchorGaps<TGapAnchors> >, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
	_goPrevious_gapAnchorIterator(me);
}

template <typename TSource, typename TGapAnchors, typename TSize>
inline void
goFurther(Iter<Gaps<TSource, AnchorGaps<TGapAnchors> >, GapsIterator<AnchorGaps<TGapAnchors> > > & me, TSize steps)
{
	_goTo_gapAnchorIterator(me, me.current.gapPos + steps);
}

//____________________________________________________________________________

template <typename TSource, typename TGapAnchors>
inline typename Source<Iter<Gaps<TSource, AnchorGaps<TGapAnchors> >, GapsIterator<AnchorGaps<TGapAnchors> > > >::Type /*returns copy*/
source(Iter<Gaps<TSource, AnchorGaps<TGapAnchors> >, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
SEQAN_CHECKPOINT
	return me.data_source;
}
template <typename TSource, typename TGapAnchors>
inline typename Source<Iter<Gaps<TSource, AnchorGaps<TGapAnchors> >, GapsIterator<AnchorGaps<TGapAnchors> > > const>::Type /*returns copy*/
source(Iter<Gaps<TSource, AnchorGaps<TGapAnchors> >, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
SEQAN_CHECKPOINT
	return me.data_source;
}

//____________________________________________________________________________

template <
	typename TCigar,
	typename TGaps1,
	typename TGaps2,
	typename TThresh>
inline void
getCigarString(
	TCigar &cigar,
	TGaps1 &gaps1,
	TGaps2 &gaps2,
	TThresh splicedGapThresh)
{
	Iterator<TGaps1> it1 = begin(gaps1);
	Iterator<TGaps2> it2 = begin(gaps2);
	clear(cigar);
	char op, lastOp = ' ';
	unsigned numOps = 0;
	while (!atEnd(it1) && !atEnd(it2))
	{
		if (isGap(it1))
		{
			if (isGap(it2))
				op = 'P';
			else if (isClipped(it2))
				op = '?';
			else
				op = 'I';
		} 
		else if (isClipped(it1))
		{
			op = '?';
		}
		else 
		{
			if (isGap(it2))
				op = 'D';
			else if (isClipped(it2))
				op = 'S';
			else
				op = 'M';
		}
		if (lastOp == op || numOps == 0)
			++numOps;
		else
		{
			if (lastOp == 'D' && numOps >= splicedGapThresh)
				lastOp = 'N';
			appendValue(cigar, lastOp);
			std::stringstream num;
			num << numOps;
			append(cigar, num.str());
			numOps = 1;
			lastOp = op;
		}
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
