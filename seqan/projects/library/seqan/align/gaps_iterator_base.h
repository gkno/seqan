#ifndef SEQAN_HEADER_GAPS_ITERATOR_BASE_H
#define SEQAN_HEADER_GAPS_ITERATOR_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Source
//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec>
struct Source<Iter<TGaps, GapsIterator<TSpec> > >
{
	typedef typename Source<TGaps>::Type TSource;
	typedef typename Iterator<TSource>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////
// Value
//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec>
struct Value<Iter<TGaps, GapsIterator<TSpec> > >
{
	typedef typename Source<Iter<TGaps, GapsIterator<TSpec> > >::Type TSource;
	typedef typename Value<TSource>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////
// GetValue
//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec>
struct GetValue<Iter<TGaps, GapsIterator<TSpec> > >:
	Value<Iter<TGaps, GapsIterator<TSpec> > > //no reference
{
};

//////////////////////////////////////////////////////////////////////////////
// Reference
//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec>
struct Reference<Iter<TGaps, GapsIterator<TSpec> > >
{
	typedef Iter<TGaps, GapsIterator<TSpec> > TIterator;
	typedef Proxy<IteratorProxy<TIterator> > Type;
};


//////////////////////////////////////////////////////////////////////////////
//specialization of basic_proxy.h/CompareType implementation

template <typename TGaps, typename TSpec, typename T>
struct CompareType<Proxy<IteratorProxy<Iter<TGaps, GapsIterator<TSpec> > > >, T>
{
	typedef T Type;
};


//???TODO: Symmetrie von CompareType herstellen


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec>
inline Iter<TGaps, GapsIterator<TSpec> > const & 
operator ++(Iter<TGaps, GapsIterator<TSpec> > const & me)
{
SEQAN_CHECKPOINT
	goNext(me);
	return me;
}

template <typename TGaps, typename TSpec>
inline Iter<TGaps, GapsIterator<TSpec> > const 
operator ++(Iter<TGaps, GapsIterator<TSpec> > const & me, int)
{
SEQAN_CHECKPOINT
	Iter<TGaps, GapsIterator<TSpec> > ret = me;
	goNext(me);
	return ret;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec>
inline Iter<TGaps, GapsIterator<TSpec> > const & 
operator --(Iter<TGaps, GapsIterator<TSpec> > const & me)
{
SEQAN_CHECKPOINT
	goPrevious(me);
	return me;
}

template <typename TGaps, typename TSpec>
inline Iter<TGaps, GapsIterator<TSpec> > const 
operator --(Iter<TGaps, GapsIterator<TSpec> > const & me, int)
{
SEQAN_CHECKPOINT
	Iter<TGaps, GapsIterator<TSpec> > ret = me;
	goPrevious(me);
	return ret;
}

//////////////////////////////////////////////////////////////////////////////
//todo???: weitere operatoren


//////////////////////////////////////////////////////////////////////////////

// insert a single gap at current iterator position
template <typename TGaps, typename TSpec>
inline void
insertGap(Iter<TGaps, GapsIterator<TSpec> > const & me)
{
SEQAN_CHECKPOINT
	insertGaps(me, 1);
}

//////////////////////////////////////////////////////////////////////////////

// remove a single gap at current iterator position
template <typename TGaps, typename TSpec>
inline void
removeGap(Iter<TGaps, GapsIterator<TSpec> > const & me)
{
SEQAN_CHECKPOINT
	removeGaps(me, 1);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec>
inline typename Reference< Iter<TGaps, GapsIterator<TSpec> > >::Type
value(Iter<TGaps, GapsIterator<TSpec> > & me)
{
SEQAN_CHECKPOINT
	typedef typename Reference< Iter<TGaps, GapsIterator<TSpec> > >::Type TReference;
	return TReference(me);
}
template <typename TGaps, typename TSpec>
inline typename Reference< Iter<TGaps, GapsIterator<TSpec> > const>::Type
value(Iter<TGaps, GapsIterator<TSpec> > const & me)
{
SEQAN_CHECKPOINT
	typedef typename Reference< Iter<TGaps, GapsIterator<TSpec> > const>::Type TReference;
	return TReference(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec>
inline typename GetValue< Iter<TGaps, GapsIterator<TSpec> > >::Type
getValue(Iter<TGaps, GapsIterator<TSpec> > & me)
{
SEQAN_CHECKPOINT
	typedef typename Value<Iter<TGaps, GapsIterator<TSpec> > >::Type TValue;
	if (isGap(me)) return gapValue<TValue>();
	else return getValue(source(me));
}
template <typename TGaps, typename TSpec>
inline typename GetValue< Iter<TGaps, GapsIterator<TSpec> > const>::Type
getValue(Iter<TGaps, GapsIterator<TSpec> > const & me)
{
SEQAN_CHECKPOINT
	typedef typename Value<Iter<TGaps, GapsIterator<TSpec> > const>::Type TValue;
	if (isGap(me)) return gapValue<TValue>();
	else return getValue(source(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec, typename TValue>
inline void
assignValue(Iter<TGaps, GapsIterator<TSpec> > & me,
			TValue const & val)
{
SEQAN_CHECKPOINT
	if (!isGap(me)) 
	{
		assignValue(source(me), val);
	}
	//else insert??? das waere neu.
}
template <typename TGaps, typename TSpec, typename TValue>
inline void
assignValue(Iter<TGaps, GapsIterator<TSpec> > const & me,
			TValue const & val)
{
SEQAN_CHECKPOINT
	if (!isGap(me)) 
	{
		assignValue(source(me), val);
	}
}
//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec, typename TValue>
inline void
moveValue(Iter<TGaps, GapsIterator<TSpec> > & me,
		  TValue const & val)
{
SEQAN_CHECKPOINT
	if (!isGap(me)) 
	{
		moveValue(source(me), val);
	}
}
template <typename TGaps, typename TSpec, typename TValue>
inline void
moveValue(Iter<TGaps, GapsIterator<TSpec> > const & me,
		  TValue const & val)
{
SEQAN_CHECKPOINT
	if (!isGap(me)) 
	{
		moveValue(source(me), val);
	}
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
