#ifndef SEQAN_HEADER_MATRIX_BASE_H
#define SEQAN_HEADER_MATRIX_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec = void>
class Matrix;

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
class Matrix<TValue, void>
{
//____________________________________________________________________________

private:
	typedef typename Size<Matrix>::Type TSize;
	typedef String<TSize> TSizeArr;
	typedef String<TValue> THost;

	TSizeArr data_lengths;
	TSizeArr data_factors;

	Holder<THost> data_host;
//____________________________________________________________________________

public:
	Matrix()
	{
		create(data_host);
	}
	Matrix(Matrix const & other_):
		data_lengths(other_.data_lengths),
		data_factors(other_.data_factors),
		data_host(other_.data_host)
	{
	}
	inline Matrix const &
	operator = (Matrix const & other_)
	{
		data_lengths = other_.data_lengths;
		data_factors = other_.data_factors;
		data_host = other_.data_host;

		return *this;
	}
//____________________________________________________________________________

	friend inline TSizeArr &
	_dataLengths(Matrix & me)
	{
		return me.data_lengths;
	}

	friend inline TSizeArr &
	_dataFactors(Matrix & me)
	{
		return me.data_factors;
	}

//____________________________________________________________________________

	friend inline bool
	dependent(Matrix & me)
	{
		return dependent(me.data_host);
	}

//____________________________________________________________________________

	friend inline void
	setHost(Matrix & me, THost & host_)
	{
		setValue(me.data_host, host_);
	}

//____________________________________________________________________________

	friend inline THost &
	host(Matrix & me)
	{
		return value(me.data_host);
	}
	friend inline THost const &
	host(Matrix const & me)
	{
		return value(me.data_host);
	}

//____________________________________________________________________________

	friend inline void
	assignHost(Matrix & me, THost const & value_)
	{
		assignValue(me.data_host, value_);
	}

//____________________________________________________________________________

	friend inline void
	moveHost(Matrix & me, THost const & value_)
	{
		moveValue(me.data_host, value_);
	}

//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
struct Value< Matrix<TValue, TSpec> >
{
	typedef TValue Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TIteratorSpec>
struct Iterator< Matrix<TValue, TSpec>, TIteratorSpec >
{
	typedef Iter<Matrix<TValue, TSpec>, PositionIterator> Type;
};

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline unsigned int
dimension(Matrix<TValue, TSpec> & me)
{
	return length(_dataLengths(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline void
setDimension(Matrix<TValue, TSpec> & me,
			 unsigned int dim_)
{
	SEQAN_ASSERT(dim_ > 0)

	fill(_dataLengths(me), dim_, 0);

	resize(_dataFactors(me), dim_);
	_dataFactors(me)[0] = 1;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline typename Size<Matrix<TValue, TSpec> >::Type
length(Matrix<TValue, TSpec> & me,
	   unsigned int dim_)
{
	return _dataLengths(me)[dim_];
}

template <typename TValue, typename TSpec>
inline typename Size<Matrix <TValue, TSpec> >::Type
length(Matrix<TValue, TSpec> & me)
{
	return length(host(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TSize>
inline void
setLength(Matrix<TValue, TSpec> & me,
		  unsigned int dim_,
		  TSize length_)
{
	SEQAN_ASSERT(length_ > 0);
	SEQAN_ASSERT(dim_ < dimension(me));

	_dataLengths(me)[dim_] = length_;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline void
resize(Matrix<TValue, TSpec> & me)
{
	typedef Matrix<TValue, TSpec> TMatrix;
	typedef typename Size<TMatrix>::Type TSize;

	unsigned int dimension_ = dimension(me);

	SEQAN_ASSERT(dimension_ > 0);

	TSize factor_ = _dataFactors(me)[0] * length(me, 0);
	for (unsigned int i = 1; (factor_ > 0) && (i < dimension_); ++i)
	{
		_dataFactors(me)[i] = factor_;
		factor_ *= length(me, i);
	}

	if (factor_ > 0)
	{
		resize(host(me), factor_);
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TPosition>
inline typename Position<Matrix <TValue, TSpec> >::Type
nextPosition(Matrix<TValue, TSpec> & me,
			 TPosition position_,
			 unsigned int dimension_)
{
	return position_ + _dataFactors(me)[dimension_];
}

template <typename TValue, typename TSpec, typename TPosition>
inline typename Position<Matrix <TValue, TSpec> >::Type
previousPosition(Matrix<TValue, TSpec> & me,
				 TPosition position_,
				 unsigned int dimension_)
{
	return position_ - _dataFactors(me)[dimension_];
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TPosition>
inline typename Size< Matrix <TValue, TSpec> >::Type
coordinate(Matrix<TValue, TSpec> & me,
		   TPosition position_,
		   unsigned int dimension_)
{
	SEQAN_ASSERT(dimension_ < dimension(me));

	if (dimension_ < dimension(me) - 1)
	{
		return (position_ / _dataFactors(me)[dimension_]) % _dataFactors(me)[dimension_ + 1];
	}
	else
	{
		return position_ / _dataFactors(me)[dimension_];
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TTag>
inline typename Iterator<Matrix <TValue, TSpec> >::Type
begin(Matrix<TValue, TSpec> & me,
	  Tag<TTag> const)
{
	return typename Iterator<Matrix <TValue, TSpec> >::Type(me, 0);
}
template <typename TValue, typename TSpec, typename TTag>
inline typename Iterator<Matrix <TValue, TSpec> >::Type
begin(Matrix<TValue, TSpec> const & me,
	  Tag<TTag> const)
{
	return typename Iterator<Matrix <TValue, TSpec> >::Type(me, 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TTag>
inline typename Iterator<Matrix <TValue, TSpec> >::Type
end(Matrix<TValue, TSpec> & me,
	  Tag<TTag> const)
{
	return typename Iterator<Matrix <TValue, TSpec> >::Type(me, length(host(me)));
}
template <typename TValue, typename TSpec, typename TTag>
inline typename Iterator<Matrix <TValue, TSpec> >::Type
end(Matrix<TValue, TSpec> const & me,
	  Tag<TTag> const)
{
	return typename Iterator<Matrix <TValue, TSpec> >::Type(me, length(host(me)));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TPosition>
inline typename Reference<Matrix<TValue, TSpec> >::Type
value(Matrix<TValue, TSpec> & me,
	  TPosition position_)
{
	return value(host(me), position_);
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Iterator: goNext
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline void
goNext(Iter< Matrix<TValue, TSpec>, PositionIterator > & me,
	   unsigned int dimension_ = 0)
{
	setPosition(me, nextPosition(container(me), position(me), dimension_));
}

//////////////////////////////////////////////////////////////////////////////
// Iterator: goPevious
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline void
goPrevious(Iter< Matrix<TValue, TSpec>, PositionIterator > & me,
		   unsigned int dimension_ = 0)
{
	setPosition(me, previousPosition(container(me), position(me), dimension_));
}

//////////////////////////////////////////////////////////////////////////////
// Iterator: coordinate

template <typename TValue, typename TSpec>
inline typename Size< Matrix<TValue, TSpec> >::Type 
coordinate(Iter< Matrix<TValue, TSpec>, PositionIterator > & me,
		   unsigned int dimension_)
{
	return coordinate(container(me), position(me), dimension_);
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
