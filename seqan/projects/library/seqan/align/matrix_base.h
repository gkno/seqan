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

#ifndef SEQAN_HEADER_MATRIX_BASE_H
#define SEQAN_HEADER_MATRIX_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

struct NDimensional;


template <typename TValue, typename TSpec = NDimensional>
class Matrix;

//////////////////////////////////////////////////////////////////////////////

/**
.Class.Matrix
..summary:A simple n-dimensional matrix type.
..cat:Miscellaneous
*/

template <typename TValue>
class Matrix<TValue, NDimensional>
{
//____________________________________________________________________________

public:
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
	~Matrix()
	{
	}
//____________________________________________________________________________

	friend inline TSizeArr &
	_dataLengths(Matrix & me)
	{
		return me.data_lengths;
	}
	friend inline TSizeArr const &
	_dataLengths(Matrix const & me)
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

	inline TValue & 
	operator () (TSize x1, TSize x2)
	{
		return value(*this, x1, x2);
	}
	inline TValue & 
	operator () (TSize x1, TSize x2, TSize x3)
	{
		return value(*this, x1, x2, x3);
	}
	inline TValue & 
	operator () (TSize x1, TSize x2, TSize x3, TSize x4)
	{
		return value(*this, x1, x2, x3, x4);
	}

//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
struct Value< Matrix<TValue, NDimensional> >
{
	typedef TValue Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TIteratorSpec>
struct Iterator< Matrix<TValue, NDimensional>, TIteratorSpec >
{
	typedef Iter<Matrix<TValue, NDimensional>, PositionIterator> Type;
};

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline unsigned int
dimension(Matrix<TValue, NDimensional> & me)
{
	return length(_dataLengths(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline void
setDimension(Matrix<TValue, NDimensional> & me,
			 unsigned int dim_)
{
	SEQAN_ASSERT(dim_ > 0)

	fill(_dataLengths(me), dim_, 0);

	resize(_dataFactors(me), dim_);
	_dataFactors(me)[0] = 1;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline typename Size<Matrix<TValue, NDimensional> >::Type
length(Matrix<TValue, NDimensional> const & me,
	   unsigned int dim_)
{
	return me.data_lengths[dim_];
}

template <typename TValue>
inline typename Size<Matrix <TValue, NDimensional> >::Type
length(Matrix<TValue, NDimensional> const & me)
{
	return length(host(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSize>
inline void
setLength(Matrix<TValue, NDimensional> & me,
		  unsigned int dim_,
		  TSize length_)
{
	SEQAN_ASSERT(length_ > 0);
	SEQAN_ASSERT(dim_ < dimension(me));

	_dataLengths(me)[dim_] = length_;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline void
resize(Matrix<TValue, NDimensional> & me)
{
	typedef Matrix<TValue, NDimensional> TMatrix;
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

template <typename TValue, typename TPosition>
inline typename Position<Matrix <TValue, NDimensional> >::Type
nextPosition(Matrix<TValue, NDimensional> & me,
			 TPosition position_,
			 unsigned int dimension_)
{
	return position_ + _dataFactors(me)[dimension_];
}

template <typename TValue, typename TPosition>
inline typename Position<Matrix <TValue, NDimensional> >::Type
previousPosition(Matrix<TValue, NDimensional> & me,
				 TPosition position_,
				 unsigned int dimension_)
{
	return position_ - _dataFactors(me)[dimension_];
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TPosition>
inline typename Size< Matrix <TValue, NDimensional> >::Type
coordinate(Matrix<TValue, NDimensional> & me,
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

template <typename TValue, typename TTag>
inline typename Iterator<Matrix <TValue, NDimensional>, Tag<TTag> const>::Type
begin(Matrix<TValue, NDimensional> & me,
	  Tag<TTag> const)
{
	return typename Iterator<Matrix <TValue, NDimensional>, Tag<TTag> const >::Type(me, 0);
}
template <typename TValue, typename TTag>
inline typename Iterator<Matrix <TValue, NDimensional> const, Tag<TTag> const>::Type
begin(Matrix<TValue, NDimensional> const & me,
	  Tag<TTag> const)
{
	return typename Iterator<Matrix <TValue, NDimensional>, Tag<TTag> const >::Type(me, 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TTag>
inline typename Iterator<Matrix <TValue, NDimensional>, Tag<TTag> const >::Type
end(Matrix<TValue, NDimensional> & me,
	  Tag<TTag> const)
{
	return typename Iterator<Matrix <TValue, NDimensional>, Tag<TTag> const >::Type(me, length(host(me)));
}
template <typename TValue, typename TTag>
inline typename Iterator<Matrix <TValue, NDimensional> const, Tag<TTag> const >::Type
end(Matrix<TValue, NDimensional> const & me,
	  Tag<TTag> const)
{
	return typename Iterator<Matrix <TValue, NDimensional>, Tag<TTag> const >::Type(me, length(host(me)));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TPosition>
inline typename Reference<Matrix<TValue, NDimensional> >::Type
value(Matrix<TValue, NDimensional> & me,
	  TPosition position_)
{
	return value(host(me), position_);
}

//____________________________________________________________________________

//two dimensional value access
template <typename TValue, typename TOrdinate1, typename TOrdinate2>
inline typename Reference<Matrix<TValue, NDimensional> >::Type
value(Matrix<TValue, NDimensional> & me,
	  TOrdinate1 i1,
	  TOrdinate2 i2)
{
	return value(host(me), i1 + i2 * _dataFactors(me)[1]);
}

//____________________________________________________________________________

//3 dimensional value access

template <typename TValue, typename TOrdinate1, typename TOrdinate2, typename TOrdinate3>
inline typename Reference<Matrix<TValue, NDimensional> >::Type
value(Matrix<TValue, NDimensional> & me,
	  TOrdinate1 i1, 
	  TOrdinate2 i2,
	  TOrdinate3 i3)
{
	return value(host(me), i1 + i2 * _dataFactors(me)[1] + i3 * _dataFactors(me)[2]);
}

//____________________________________________________________________________

//4 dimensional value access

template <typename TValue, typename TOrdinate1, typename TOrdinate2, typename TOrdinate3, typename TOrdinate4>
inline typename Reference<Matrix<TValue, NDimensional> >::Type
value(Matrix<TValue, NDimensional> & me,
	  TOrdinate1 i1, 
	  TOrdinate2 i2,
	  TOrdinate3 i3,
	  TOrdinate4 i4)
{
	return value(host(me), i1 + i2 * _dataFactors(me)[1] + i3 * _dataFactors(me)[2] + i4 * _dataFactors(me)[3]);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Iterator: goNext
//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline void
goNext(Iter< Matrix<TValue, NDimensional>, PositionIterator > & me,
	   unsigned int dimension_ = 0)
{
	setPosition(me, nextPosition(container(me), position(me), dimension_));
}

//////////////////////////////////////////////////////////////////////////////
// Iterator: goPevious
//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline void
goPrevious(Iter< Matrix<TValue, NDimensional>, PositionIterator > & me,
		   unsigned int dimension_ = 0)
{
	setPosition(me, previousPosition(container(me), position(me), dimension_));
}

//////////////////////////////////////////////////////////////////////////////
// Iterator: coordinate

template <typename TValue>
inline typename Size< Matrix<TValue, NDimensional> >::Type 
coordinate(Iter< Matrix<TValue, NDimensional>, PositionIterator > & me,
		   unsigned int dimension_)
{
	return coordinate(container(me), position(me), dimension_);
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
