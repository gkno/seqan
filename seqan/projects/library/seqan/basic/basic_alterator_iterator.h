#ifndef SEQAN_HEADER_BASIC_ALTERATOR_ITERATOR_H
#define SEQAN_HEADER_BASIC_ALTERATOR_ITERATOR_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// container

template <typename THost, typename TSpec>
inline typename Container<Alterator<THost, TSpec> >::Type &
container(Alterator<THost, TSpec> & me)
{
SEQAN_CHECKPOINT
	return container(host(me));
}
template <typename THost, typename TSpec>
inline typename Container<Alterator<THost, TSpec> >::Type &
container(Alterator<THost, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return container(host(me));
}

//////////////////////////////////////////////////////////////////////////////
// setContainer

template <typename THost, typename TSpec, typename TIterator2>
inline void
setContainer(Alterator<THost, TSpec> & me,	
			 TIterator2 & container_)
{
SEQAN_CHECKPOINT
	setContainer(host(me), container_);
}
template <typename THost, typename TSpec, typename TIterator2>
inline void
setContainer(Alterator<THost, TSpec> & me,	
			 TIterator2 const & container_)
{
SEQAN_CHECKPOINT
	setContainer(host(me), container_);
}

//////////////////////////////////////////////////////////////////////////////
// position

template <typename THost, typename TSpec>
inline typename Position<Alterator<THost, TSpec> >::Type 
position(Alterator<THost, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return position(host(me));
}

//////////////////////////////////////////////////////////////////////////////
// setPosition

template <typename THost, typename TSpec, typename TPosition>
inline void 
setPosition(Alterator<THost, TSpec> & me,
			TPosition pos_)
{
SEQAN_CHECKPOINT
	setPosition(host(me));
}

//////////////////////////////////////////////////////////////////////////////
// value, operator *

template <typename THost, typename TSpec>
inline typename Reference<Alterator<THost, TSpec> >::Type 
value(Alterator<THost, TSpec> & me)
{
SEQAN_CHECKPOINT
	return value(host(me));
}
template <typename THost, typename TSpec>
inline typename Reference<Alterator<THost, TSpec> const>::Type 
value(Alterator<THost, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return value(host(me));
}

//____________________________________________________________________________

template <typename THost, typename TSpec>
inline typename Reference<Alterator<THost, TSpec> >::Type 
operator * (Alterator<THost, TSpec> & me)
{
SEQAN_CHECKPOINT
	return value(me);
}
template <typename THost, typename TSpec>
inline typename Reference<Alterator<THost, TSpec> const>::Type 
operator * (Alterator<THost, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return value(me);
}

//////////////////////////////////////////////////////////////////////////////
// getValue

template <typename THost, typename TSpec>
inline typename GetValue<Alterator<THost, TSpec> >::Type 
getValue(Alterator<THost, TSpec> & me)
{
SEQAN_CHECKPOINT
	return getValue(host(me));
}
template <typename THost, typename TSpec>
inline typename GetValue<Alterator<THost, TSpec> >::Type 
getValue(Alterator<THost, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return getValue(host(me));
}

/////////////////////////////////////////////////////////////////////////////
// assignValue

template <typename THost, typename TSpec, typename TValue>
inline void
assignValue(Alterator<THost, TSpec> & me,
			TValue const & _value)
{
SEQAN_CHECKPOINT
	assignValue(host(me), _value);
}
template <typename THost, typename TSpec, typename TValue>
inline void
assignValue(Alterator<THost, TSpec> const & me,
			TValue const & _value)
{
SEQAN_CHECKPOINT
	assignValue(host(me), _value);
}

/////////////////////////////////////////////////////////////////////////////
// moveValue

template <typename THost, typename TSpec, typename TValue>
inline void
moveValue(Alterator<THost, TSpec> & me,
		  TValue const & _value)
{
SEQAN_CHECKPOINT
	moveValue(host(me), _value);
}
template <typename THost, typename TSpec, typename TValue>
inline void
moveValue(Alterator<THost, TSpec> const & me,
		  TValue const & _value)
{
SEQAN_CHECKPOINT
	moveValue(host(me), _value);
}

//////////////////////////////////////////////////////////////////////////////
// operator ==

template <typename THost, typename TSpec>
inline bool 
operator == (Alterator<THost, TSpec> const & left,
			 Alterator<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return host(left) == host(right);
}
template <typename THost, typename TSpec, typename TRight>
inline bool 
operator == (Alterator<THost, TSpec> const & left,
			 TRight const & right)
{
SEQAN_CHECKPOINT
	return host(left) == right;
}
template <typename TLeft, typename THost, typename TSpec>
inline bool 
operator == (TLeft const & left,
			 Alterator<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return left == host(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator !=

template <typename THost, typename TSpec>
inline bool 
operator != (Alterator<THost, TSpec> const & left,
			 Alterator<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return host(left) != host(right);
}
template <typename THost, typename TSpec, typename TRight>
inline bool 
operator != (Alterator<THost, TSpec> const & left,
			 TRight const & right)
{
SEQAN_CHECKPOINT
	return host(left) != right;
}
template <typename TLeft, typename THost, typename TSpec>
inline bool 
operator != (TLeft const & left,
			 Alterator<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return left != host(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator <

template <typename THost, typename TSpec>
inline bool 
operator < (Alterator<THost, TSpec> const & left,
			Alterator<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return host(left) < host(right);
}
template <typename THost, typename TSpec, typename TRight>
inline bool 
operator < (Alterator<THost, TSpec> const & left,
			TRight const & right)
{
SEQAN_CHECKPOINT
	return host(left) < right;
}
template <typename TLeft, typename THost, typename TSpec>
inline bool 
operator < (TLeft const & left,
			Alterator<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return left < host(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator <=

template <typename THost, typename TSpec>
inline bool 
operator <= (Alterator<THost, TSpec> const & left,
			 Alterator<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return host(left) <= host(right);
}
template <typename THost, typename TSpec, typename TRight>
inline bool 
operator <= (Alterator<THost, TSpec> const & left,
			 TRight const & right)
{
SEQAN_CHECKPOINT
	return host(left) <= right;
}
template <typename TLeft, typename THost, typename TSpec>
inline bool 
operator <= (TLeft const & left,
			 Alterator<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return left <= host(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator >

template <typename THost, typename TSpec>
inline bool 
operator > (Alterator<THost, TSpec> const & left,
			Alterator<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return host(left) > host(right);
}
template <typename THost, typename TSpec, typename TRight>
inline bool 
operator > (Alterator<THost, TSpec> const & left,
			TRight const & right)
{
SEQAN_CHECKPOINT
	return host(left) > right;
}
template <typename TLeft, typename THost, typename TSpec>
inline bool 
operator > (TLeft const & left,
			Alterator<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return left > host(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator >=

template <typename THost, typename TSpec>
inline bool 
operator >= (Alterator<THost, TSpec> const & left,
			 Alterator<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return host(left) >= host(right);
}
template <typename THost, typename TSpec, typename TRight>
inline bool 
operator >= (Alterator<THost, TSpec> const & left,
			 TRight const & right)
{
SEQAN_CHECKPOINT
	return host(left) >= right;
}
template <typename TLeft, typename THost, typename TSpec>
inline bool 
operator >= (TLeft const & left,
			 Alterator<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return left >= host(right);
}

//////////////////////////////////////////////////////////////////////////////
// goNext

template <typename THost, typename TSpec>
inline void
goNext(Alterator<THost, TSpec> & me)
{
SEQAN_CHECKPOINT
	goNext(host(me));
}
//____________________________________________________________________________

template <typename THost, typename TSpec>
inline Alterator<THost, TSpec> const &
operator ++ (Alterator<THost, TSpec> & me)
{
SEQAN_CHECKPOINT
	goNext(me);
	return me;
}

template <typename THost, typename TSpec>
inline Alterator<THost, TSpec> const &
operator ++ (Alterator<THost, TSpec> & me, int)
{
SEQAN_CHECKPOINT
	Alterator<THost, TSpec> temp_(me);
	goNext(me);
	return me;
}

//////////////////////////////////////////////////////////////////////////////
// goPrevious

template <typename THost, typename TSpec>
inline void
goPrevious(Alterator<THost, TSpec> & me)
{
SEQAN_CHECKPOINT
	goPrevious(host(me));
}
//____________________________________________________________________________

template <typename THost, typename TSpec>
inline Alterator<THost, TSpec> const &
operator -- (Alterator<THost, TSpec> & me)
{
SEQAN_CHECKPOINT
	goPrevious(me);
	return me;
}

template <typename THost, typename TSpec>
inline Alterator<THost, TSpec> const &
operator -- (Alterator<THost, TSpec> & me, int)
{
SEQAN_CHECKPOINT
	Alterator<THost, TSpec> temp_(me);
	goPrevious(me);
	return me;
}

//////////////////////////////////////////////////////////////////////////////
// operator +

template <typename THost, typename TSpec, typename TIntegral>
inline Alterator<THost, TSpec>  
operator + (Alterator<THost, TSpec> const & left,
			TIntegral right)
{
SEQAN_CHECKPOINT
	return Alterator<THost, TSpec>(host(left) + right);
}
template <typename THost, typename TSpec, typename TIntegral>
inline Alterator<THost, TSpec>  
operator + (TIntegral left,
			Alterator<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return Alterator<THost, TSpec>(host(right) + left);
}

//////////////////////////////////////////////////////////////////////////////
// operator +=

template <typename THost, typename TSpec, typename TIntegral>
inline Alterator<THost, TSpec> &
operator += (Alterator<THost, TSpec> & left,
			 TIntegral right)
{
SEQAN_CHECKPOINT
	host(left) += right;
	return left;
}

//////////////////////////////////////////////////////////////////////////////
// operator -

template <typename THost, typename TSpec, typename TIntegral>
inline Alterator<THost, TSpec>  
operator - (Alterator<THost, TSpec> const & left,
			TIntegral right)
{
SEQAN_CHECKPOINT
	return Alterator<THost, TSpec>(host(left) - right);
}

//____________________________________________________________________________

template <typename THost, typename TSpec>
inline typename Position<Alterator<THost, TSpec> >::Type  
operator - (Alterator<THost, TSpec> const & left,
			Alterator<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return host(left) - host(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator -=

template <typename THost, typename TSpec, typename TIntegral>
inline Alterator<THost, TSpec> &
operator -= (Alterator<THost, TSpec> & left,
			 TIntegral right)
{
SEQAN_CHECKPOINT
	host(left) -= right;
	return left;
}

//////////////////////////////////////////////////////////////////////////////
// atEnd

template <typename THost, typename TSpec>
inline bool
atEnd(Alterator<THost, TSpec> & me)
{
SEQAN_CHECKPOINT
	return atEnd(host(me));
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
