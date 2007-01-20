#ifndef SEQAN_HEADER_BASIC_ITERATOR_ADAPTOR_H
#define SEQAN_HEADER_BASIC_ITERATOR_ADAPTOR_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Tag

//An iterator that adapts a std iterator to a default seqan iterator
template <typename TIterator, typename TSpec = Default>
struct AdaptorIterator;

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct Iterator_Default_Imp<T, Rooted>
{
	typedef typename Iterator<T, Standard>::Type TStandardIterator;
	typedef Iter<T, AdaptorIterator<TStandardIterator> > Type;
};

//////////////////////////////////////////////////////////////////////////////
// Adaptor Iterator
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Adaptor Iterator:
..cat:Iterators
..general:Class.Iter
..summary:Adapts iterators to @Concept.Rooted Iterator@.
..signature:Iter<TContainer, AdaptorIterator<TIterator [, TSpec]> >
..param.TContainer:Type of the container that can be iterated by $TIterator$.
...remarks:Use @Metafunction.Container@ to get the container type for a given iterator.
..param.TIterator:Type of the iterator that is adapted to @Concept.Rooted Iterator@.
..remarks.text:Adaptor iterators can implicitly converted to $TIterator$.
*/

template <typename TContainer, typename TIterator, typename TSpec>
class Iter<TContainer, AdaptorIterator<TIterator, TSpec> >
{
private:
	typename _Pointer<TContainer>::Type data_container;
	TIterator data_iterator;
//____________________________________________________________________________

/**
.Memfunc.AdaptorIterator#Iter:
..class:Spec.Adaptor Iterator
..summary:Constructor
..signature:Iter()
..signature:Iter(iter)
..signature:Iter(container [, iterator])
..param.iter:Another adaptor iterator object.
..param.container:The corresponding container object.
..param.iterator:A iterator of $container$. (optional)
...remarks.text:If this argument is omitted, the adaptor iterator is initialized to the @Function.begin.begin iterator@ of $container$.
*/

public:
	Iter():
		data_container(0)
	{
SEQAN_CHECKPOINT
		data_iterator = TIterator();
	}
	Iter(typename _Parameter<TContainer>::Type container_):
		data_container(_toPointer(container_)),
		data_iterator(begin(container_))
	{
SEQAN_CHECKPOINT
	}
	Iter(typename _Parameter<TContainer>::Type container_, TIterator it_):
		data_container(_toPointer(container_)),
		data_iterator(it_)
	{
SEQAN_CHECKPOINT
	}
	Iter(Iter const & other_):
		data_container(other_.data_container),
		data_iterator(other_.data_iterator)
	{
SEQAN_CHECKPOINT
	}
/*
	template <typename TSource>
	Iter(TSource & source)
	{
SEQAN_CHECKPOINT
		assign(*this, source);
	}
	template <typename TSource>
	Iter(TSource const & source)
	{
SEQAN_CHECKPOINT
		assign(*this, source);
	}
*/

	~Iter()
	{
SEQAN_CHECKPOINT
	}

	Iter const & 
	operator = (Iter const & other_)
	{
SEQAN_CHECKPOINT
		data_container = other_.data_container;
		data_iterator = other_.data_iterator;
		return *this;
	}
/*
	template <typename TSource>
	Iter const & 
	operator = (TSource & source)
	{
SEQAN_CHECKPOINT
		assign(*this, source);
		return *this;
	}
	template <typename TSource>
	Iter const & 
	operator = (TSource const & source)
	{
SEQAN_CHECKPOINT
		assign(*this, source);
		return *this;
	}
*/
//____________________________________________________________________________

	friend inline typename _Parameter<TContainer>::Type 
	container(Iter & me)
	{
SEQAN_CHECKPOINT
		return _toParameter<TContainer>(me.data_container);
	}
	friend inline typename _Parameter<TContainer>::Type 
	container(Iter const & me)
	{
SEQAN_CHECKPOINT
		return _toParameter<TContainer>(me.data_container);
	}
//____________________________________________________________________________

	friend inline void
	setContainer(Iter & me,	typename _Parameter<TContainer>::Type container_)
	{
SEQAN_CHECKPOINT
		if (me.data_container && me.data_iterator != TIterator())
		{
			typename Position<Iter>::Type pos = position(me);
			me.data_container = _toPointer(container_);
			setPosition(me, pos);
		}
		else
		{	
			me.data_container = _toPointer(container_);
		}
	}

//____________________________________________________________________________

	friend inline TIterator &
	hostIterator(Iter & me)
	{
SEQAN_CHECKPOINT
		return me.data_iterator;
	}
	friend inline TIterator const &
	hostIterator(Iter const & me)
	{
SEQAN_CHECKPOINT
		return me.data_iterator;
	}

//____________________________________________________________________________

	operator TIterator () const
	{
SEQAN_CHECKPOINT
		return data_iterator;
	}

//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// position
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec>
inline typename Position<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const>::Type 
position(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me)
{
SEQAN_CHECKPOINT
	return hostIterator(me) - begin(container(me), Standard());
}

//____________________________________________________________________________

template <typename TContainer, typename TIterator, typename TSpec, typename TContainer2>
inline typename Position<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const>::Type 
position(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me,
		 TContainer2 const &)
{
SEQAN_CHECKPOINT
	return hostIterator(me) - begin(container(me), Standard());
}

//////////////////////////////////////////////////////////////////////////////
// setPosition
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec, typename TPosition>
inline void 
setPosition(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me,
			TPosition pos_)
{
SEQAN_CHECKPOINT
	hostIterator(me) = begin(container(me), Standard()) + pos_;
}

//////////////////////////////////////////////////////////////////////////////
// value
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec>
inline typename Reference<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > >::Type 
value(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me)
{
SEQAN_CHECKPOINT
	return value(hostIterator(me));
}
template <typename TContainer, typename TIterator, typename TSpec>
inline typename Reference<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > >::Type 
value(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me)
{
SEQAN_CHECKPOINT
	return value(hostIterator(me));
}

/////////////////////////////////////////////////////////////////////////////
// assignValue
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec, typename TValue>
inline void
assignValue(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me,
			TValue const & _value)
{
SEQAN_CHECKPOINT
	assignValue(hostIterator(me), _value);
}
template <typename TContainer, typename TIterator, typename TSpec, typename TValue>
inline void
assignValue(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me,
			TValue const & _value)
{
SEQAN_CHECKPOINT
	assignValue(hostIterator(me), _value);
}

/////////////////////////////////////////////////////////////////////////////
// moveValue
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec, typename TValue>
inline void
moveValue(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me,
		  TValue const & _value)
{
SEQAN_CHECKPOINT
	moveValue(hostIterator(me), _value);
}
template <typename TContainer, typename TIterator, typename TSpec, typename TValue>
inline void
moveValue(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me,
		  TValue const & _value)
{
SEQAN_CHECKPOINT
	moveValue(hostIterator(me), _value);
}

//////////////////////////////////////////////////////////////////////////////
// operator ==
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec>
inline bool 
operator == (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
			 Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & right)
{
SEQAN_CHECKPOINT
	return hostIterator(left) == hostIterator(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator !=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec>
inline bool 
operator != (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
			 Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & right)
{
SEQAN_CHECKPOINT
	return hostIterator(left) != hostIterator(right);
}

//////////////////////////////////////////////////////////////////////////////
// goNext
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec>
inline void
goNext(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me)
{
SEQAN_CHECKPOINT
	goNext(hostIterator(me));
}

//////////////////////////////////////////////////////////////////////////////
// goPrevious
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec>
inline void
goPrevious(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me)
{
SEQAN_CHECKPOINT
	goPrevious(hostIterator(me));
}

//////////////////////////////////////////////////////////////////////////////
// operator +
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec, typename TIntegral>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> >  
operator + (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
			TIntegral right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, AdaptorIterator<TIterator, TSpec> >(container(left), hostIterator(left) + right);
}
template <typename TContainer, typename TIterator, typename TSpec, typename TIntegral>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> >  
operator + (TIntegral left,
			Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, AdaptorIterator<TIterator, TSpec> >(container(right), hostIterator(right) + left);
}

//////////////////////////////////////////////////////////////////////////////
// operator +=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec, typename TIntegral>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> > &
operator += (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & left,
			 TIntegral right)
{
SEQAN_CHECKPOINT
	hostIterator(left) += right;
	return left;
}

//////////////////////////////////////////////////////////////////////////////
// operator -
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec, typename TIntegral>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> >  
operator - (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
			TIntegral right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, AdaptorIterator<TIterator, TSpec> >(container(left), hostIterator(left) - right);
}

//____________________________________________________________________________

template <typename TContainer, typename TIterator, typename TSpec>
inline typename Position<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > >::Type  
operator - (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
			Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & right)
{
SEQAN_CHECKPOINT
	return hostIterator(left) - hostIterator(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator -=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec, typename TIntegral>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> > &
operator -= (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & left,
			 TIntegral right)
{
SEQAN_CHECKPOINT
	hostIterator(left) -= right;
	return left;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// atEnd
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec>
inline bool
atEnd(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me)
{
SEQAN_CHECKPOINT
	return atEnd(me, container(me));
}


//////////////////////////////////////////////////////////////////////////////
// assign
//////////////////////////////////////////////////////////////////////////////
/* ??? TODO: more sophisticated iterator conversion/assignment

template <typename TContainer, typename TIterator, typename TSpec, typename TSource>
inline void
assign(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me,
	   TSource const & source)
{
SEQAN_CHECKPOINT
}
*/
//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
