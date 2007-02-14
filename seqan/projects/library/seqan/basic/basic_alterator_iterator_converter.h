#ifndef SEQAN_HEADER_BASIC_ALTERATOR_ITERATOR_CONVERTER_H
#define SEQAN_HEADER_BASIC_ALTERATOR_ITERATOR_CONVERTER_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// METAFUNCTIONS
//////////////////////////////////////////////////////////////////////////////

template <typename TIterator, typename TTargetValue, typename TSpec>
struct _HolderSpec<Alterator<TIterator, Convert<TTargetValue, TSpec> > >
{
	typedef Simple Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T, typename TTargetValue, typename TSpec>
struct Iterator<T, Convert<TTargetValue, TSpec> >
{
	typedef typename Iterator<T, TSpec>::Type TIterator;
	typedef Alterator<TIterator, Convert<TTargetValue> > Type;
};
template <typename T, typename TTargetValue>
struct Iterator<T, Convert<TTargetValue> >
{
	typedef typename Iterator<T>::Type TIterator;
	typedef Alterator<TIterator, Convert<TTargetValue> > Type;
};

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TIterator, typename TTargetValue>
struct Value<Alterator<TIterator, Convert<TTargetValue> > >
{
	typedef TTargetValue Type;
};
template <typename TIterator, typename TTargetValue>
struct Value<Alterator<TIterator, Convert<TTargetValue> > const>
{
	typedef TTargetValue Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TIterator, typename TTargetValue>
struct GetValue<Alterator<TIterator, Convert<TTargetValue> > >:
	Convert<TTargetValue, typename Value<TIterator>::Type> //!
{
};
template <typename TIterator, typename TTargetValue>
struct GetValue<Alterator<TIterator, Convert<TTargetValue> > const>:
	Convert<TTargetValue, typename Value<TIterator>::Type> //!
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename TIterator, typename TTargetValue>
struct Reference<Alterator<TIterator, Convert<TTargetValue> > >
{
	typedef Alterator<TIterator, Convert<TTargetValue> > TAlterator;
	typedef Proxy<IteratorProxy<TAlterator> > Type;
};

template <typename TIterator, typename TTargetValue>
struct Reference<Alterator<TIterator, Convert<TTargetValue> > const>
{
	typedef Alterator<TIterator, Convert<TTargetValue> > const TAlterator;
	typedef Proxy<IteratorProxy<TAlterator> > Type;
};

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// value

template <typename TIterator, typename TTargetValue>
inline typename Reference<Alterator<TIterator, Convert<TTargetValue, void> > >::Type 
value(Alterator<TIterator, Convert<TTargetValue, void> > & me)
{
SEQAN_CHECKPOINT
	return typename Reference<Alterator<TIterator, Convert<TTargetValue, void> > >::Type(me);
}
template <typename TIterator, typename TTargetValue>
inline typename Reference<Alterator<TIterator, Convert<TTargetValue, void> > const>::Type 
value(Alterator<TIterator, Convert<TTargetValue, void> > const & me)
{
SEQAN_CHECKPOINT
	return typename Reference<Alterator<TIterator, Convert<TTargetValue, void> > const>::Type(me);
}

//////////////////////////////////////////////////////////////////////////////
// getValue

template <typename TIterator, typename TTargetValue>
inline typename GetValue<Alterator<TIterator, Convert<TTargetValue, void> > >::Type 
getValue(Alterator<TIterator, Convert<TTargetValue, void> > & me)
{
SEQAN_CHECKPOINT
	return convert<TTargetValue>(getValue(me));
}
template <typename TIterator, typename TTargetValue>
inline typename GetValue<Alterator<TIterator, Convert<TTargetValue, void> > const>::Type 
getValue(Alterator<TIterator, Convert<TTargetValue, void> > const & me)
{
SEQAN_CHECKPOINT
	return convert<TTargetValue>(getValue(me));
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Converting Proxy
//////////////////////////////////////////////////////////////////////////////

template <typename TIteratorLeft, typename TLeft, typename TIteratorRight, typename TRight>
inline bool 
operator == (Proxy<IteratorProxy<Alterator<TIteratorLeft, Convert<TLeft, void> > > > const & left_,
			 Proxy<IteratorProxy<Alterator<TIteratorRight, Convert<TRight, void> > > > const & right_)
{
SEQAN_CHECKPOINT
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(getValue(left_)) == convert<TCompareType>(getValue(right_));
}
template <typename TIteratorLeft, typename TLeft, typename TRight>
inline bool 
operator == (Proxy<IteratorProxy<Alterator<TIteratorLeft, Convert<TLeft, void> > > > const & left_,
			 TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(getValue(left_)) == convert<TCompareType>(right_);
}
template <typename TIteratorLeft, typename TLeft, typename TRightSpec>
inline bool 
operator == (Proxy<IteratorProxy<Alterator<TIteratorLeft, Convert<TLeft, void> > > > const & left_,
			 Proxy<TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef typename Value<Proxy<TRightSpec> >::Type TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(getValue(left_)) == convert<TCompareType>(getValue(right_));
}
template <typename TLeft, typename TIteratorRight, typename TRight>
inline bool 
operator == (TLeft const & left_,
			 Proxy<IteratorProxy<Alterator<TIteratorRight, Convert<TRight, void> > > > const & right_)
{
SEQAN_CHECKPOINT
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(getValue(right_));
}
template <typename TLeftSpec, typename TIteratorRight, typename TRight>
inline bool 
operator == (Proxy<TLeftSpec> const & left_,
			 Proxy<IteratorProxy<Alterator<TIteratorRight, Convert<TRight, void> > > > const & right_)
{
SEQAN_CHECKPOINT
	typedef typename Value<Proxy<TLeftSpec> >::Type TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(getValue(left_)) == convert<TCompareType>(getValue(right_));
}


//TODO???: übrige Vergleichsoperatoren

//____________________________________________________________________________


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
