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
  Author: Knut Reinert <knut.reinert@fu-berlin.de>
 ============================================================================
  Adaptions for STL vectors to SeqAn strings.
 ==========================================================================*/

#ifndef SEQAN_HEADER_STD_VECTOR_H
#define SEQAN_HEADER_STD_VECTOR_H

#include <vector>

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

/**
.Adaption."std::vector"
..summary:Adaption for STL vector class.
 */

// ===========================================================================
// Metafunctions
// ===========================================================================

    
///.Metafunction.IsContiguous.param.T.type:Adaption.std::vector
template <typename TChar, typename TAlloc>
struct IsContiguous< ::std::vector<TChar, TAlloc> >
{
    enum { VALUE = true };
};

template <typename  TChar, typename TAlloc>
struct IsContiguous< ::std::vector<TChar, TAlloc> const>
        : IsContiguous< ::std::vector<TChar, TAlloc> > {};

///.Metafunction.Value.param.T.type:Adaption.std::vector
template <typename TChar, typename TAlloc>
struct Value< ::std::vector<TChar, TAlloc> >
{
	typedef typename ::std::vector<TChar, TAlloc>::value_type Type;
};

template <typename TChar, typename TAlloc>
struct Value< ::std::vector<TChar, TAlloc> const>
        : Value< ::std::vector<TChar, TAlloc> > {};

///.Metafunction.GetValue.param.T.type:Adaption.std::vector
// TODO(holtgrew): GetValue is a reference?! I thought the reverse was true in respect to Value<>.
template <typename TChar, typename TAlloc>
struct GetValue< ::std::vector<TChar, TAlloc> >
{
	typedef typename ::std::vector<TChar, TAlloc>::reference Type;
};

template <typename TChar, typename TAlloc>
struct GetValue< ::std::vector<TChar,  TAlloc> const>
{
	typedef typename ::std::vector<TChar, TAlloc>::const_reference Type;
};

///.Metafunction.GetValue.param.T.type:Adaption.std::vector
// TODO(holtgrew): GetValue is a reference?! I thought the reverse was true in respect to Value<>.
template <typename TChar, typename TAlloc>
struct Reference< ::std::vector<TChar, TAlloc> >
{
	typedef typename ::std::vector<TChar, TAlloc>::reference Type;
};

template <typename TChar,  typename TAlloc>
struct Reference< ::std::vector<TChar, TAlloc> const>
{
	typedef typename ::std::vector<TChar,  TAlloc>::const_reference Type;
};

///.Metafunction.Iterator.param.T.type:Adaption.std::vector
template <typename TChar, typename TAlloc>
struct Iterator< ::std::vector<TChar, TAlloc>, Rooted>
{
	typedef ::std::vector<TChar, TAlloc> TVector_;
	typedef Iter<TVector_, StdIteratorAdaptor> TIterator_;
	typedef Iter<TVector_, AdaptorIterator<TIterator_> > Type;
};

template <typename TChar, typename TAlloc>
struct Iterator< ::std::vector<TChar, TAlloc> const, Rooted>
{
	typedef ::std::vector<TChar, TAlloc> const TVector_;
	typedef Iter<TVector_, StdIteratorAdaptor> TIterator_;
	typedef Iter<TVector_, AdaptorIterator<TIterator_> > Type;
};

template <typename TChar,  typename TAlloc>
struct Iterator< ::std::vector<TChar, TAlloc>, Standard >
{
	typedef Iter< ::std::vector<TChar,  TAlloc>, StdIteratorAdaptor > Type;
};

template <typename TChar,  typename TAlloc>
struct Iterator< ::std::vector<TChar,  TAlloc> const, Standard>
{
	typedef Iter< ::std::vector<TChar, TAlloc> const, StdIteratorAdaptor > Type;
};

///.Metafunction.Position.param.T.type:Adaption.std::vector
template <typename TChar,  typename TAlloc>
struct Position< ::std::vector<TChar, TAlloc> >
{
	typedef typename ::std::vector<TChar,  TAlloc>::size_type Type;
};

template <typename TChar,  typename TAlloc>
struct Position< ::std::vector<TChar,  TAlloc> const>
        : Position< ::std::vector<TChar,  TAlloc> > {};

///.Metafunction.Size.param.T.type:Adaption.std::vector
template <typename TChar,  typename TAlloc>
struct Size< ::std::vector<TChar, TAlloc> >
{
	typedef typename ::std::vector<TChar, TAlloc>::size_type Type;
};

template <typename TChar, typename TAlloc>
struct Size< ::std::vector<TChar, TAlloc> const>
        : Size< ::std::vector<TChar, TAlloc> > {};

///.Metafunction.Size.param.T.type:Adaption.std::vector
template <typename TChar, typename TAlloc>
struct DefaultOverflowImplicit< ::std::vector<TChar, TAlloc> >
{
	typedef Generous Type;
};

// ===========================================================================
// Functions
// ===========================================================================

///.Function.id.param.object.type:Adaption.std::vector
template <typename TChar,  typename TAlloc>
inline void const * 
id(::std::vector<TChar, TAlloc> const & me)
{
    SEQAN_CHECKPOINT;
	if (me.empty())
		return NULL;
	else
		return (& *(me.end() - 1)) + 1;
}

///.Function.begin.param.object.type:Adaption.std::vector
template <typename TChar,  typename TAlloc>
inline typename Iterator< ::std::vector<TChar,  TAlloc>, Standard>::Type 
begin(::std::vector<TChar,  TAlloc> & me,
	  Standard)
{
    SEQAN_CHECKPOINT;
	return typename Iterator< ::std::vector<TChar,  TAlloc>, Standard>::Type(me.begin());
}
template <typename TChar,  typename TAlloc>
inline typename Iterator< ::std::vector<TChar,  TAlloc> const, Standard>::Type 
begin(::std::vector<TChar, TAlloc> const & me,
	  Standard)
{
    SEQAN_CHECKPOINT;
	return typename Iterator< ::std::vector<TChar,  TAlloc> const, Standard>::Type(me.begin());
}

///.Function.end.param.object.type:Adaption.std::vector
template <typename TChar, typename TAlloc>
inline typename Iterator< ::std::vector<TChar, TAlloc>, Standard>::Type 
end(::std::vector<TChar,  TAlloc> & me,
	Standard)
{
    SEQAN_CHECKPOINT;
	return typename Iterator< ::std::vector<TChar, TAlloc>, Standard>::Type(me.end());
}
template <typename TChar,  typename TAlloc>
inline typename Iterator< ::std::vector<TChar,  TAlloc> const, Standard>::Type 
end(::std::vector<TChar,  TAlloc> const & me,
	Standard)
{
    SEQAN_CHECKPOINT;
	return typename Iterator< ::std::vector<TChar,  TAlloc> const, Standard>::Type(me.end());
}

///.Function.value.param.container.type:Adaption.std::vector
template <typename TChar,  typename TAlloc, typename TPos>
inline typename GetValue< ::std::vector<TChar, TAlloc> >::Type
value(::std::vector<TChar,  TAlloc> & me, 
	  TPos pos)
{
    SEQAN_CHECKPOINT;
	return me[pos];
} 
template <typename TChar,  typename TAlloc, typename TPos>
inline typename GetValue< ::std::vector<TChar,  TAlloc> const>::Type
value(::std::vector<TChar, TAlloc> const & me, 
	  TPos pos)
{
    SEQAN_CHECKPOINT;
	return me[pos];
}

///.Function.length.param.object.type:Adaption.std::vector
template <typename TChar, typename TAlloc>
inline typename Size< ::std::vector<TChar, TAlloc> >::Type
length(::std::vector<TChar, TAlloc> const & me)
{
    SEQAN_CHECKPOINT;
	return me.size();
}

///.Function.capacity.param.object.type:Adaption.std::vector
template <typename TChar, typename TAlloc>
inline typename Size< ::std::vector<TChar, TAlloc> >::Type
capacity(::std::vector<TChar, TAlloc> const & me)
{
    SEQAN_CHECKPOINT;
	return me.capacity();
}

///.Function.empty.param.object.type:Adaption.std::vector
template <typename TChar, typename TAlloc>
inline bool
empty(::std::vector<TChar, TAlloc> const & me)
{
    SEQAN_CHECKPOINT;
	return me.empty();
}

///.Function.clear.param.object.type:Adaption.std::vector
template <typename TChar,  typename TAlloc>
inline void
clear(::std::vector<TChar, TAlloc> & me)
{
    SEQAN_CHECKPOINT;
	me.clear();
}

//////////////////////////////////////////////////////////////////////////////
//assign to ::std::vector

///.Function.assign.param.target.type:Adaption.std::vector
///.Function.assign.param.source.type:Adaption.std::vector

template <typename TChar,  typename TAlloc, typename TSource>
inline void 
assign(::std::vector<TChar,  TAlloc> & target, 
	   TSource & source)
{
    SEQAN_CHECKPOINT;
	assign(target, source, Generous());
}
template <typename TChar,  typename TAlloc, typename TSource>
inline void 
assign(::std::vector<TChar,  TAlloc> & target, 
	   TSource const & source)
{
    SEQAN_CHECKPOINT;
	assign(target, source, Generous());
}

template <typename TChar,  typename TAlloc, typename TSource, typename TSize>
inline void 
assign(::std::vector<TChar,  TAlloc> & target, 
	   TSource & source,
	   TSize limit)
{
    SEQAN_CHECKPOINT;
	assign(target, source, limit, Generous());
}
template <typename TChar,  typename TAlloc, typename TSource, typename TSize>
inline void 
assign(::std::vector<TChar,  TAlloc> & target, 
	   TSource const & source,
	   TSize limit)
{
    SEQAN_CHECKPOINT;
	assign(target, source, limit, Generous());
}

//____________________________________________________________________________

template <typename TChar,  typename TAlloc, typename TSource>
inline void 
assign(::std::vector<TChar, TAlloc> & target, 
	   TSource & source,
	   Generous)
{
    SEQAN_CHECKPOINT;
	target.assign(begin(source, Standard()), end(source, Standard()));
}
template <typename TChar, typename TAlloc, typename TSource>
inline void 
assign(::std::vector<TChar, TAlloc> & target, 
	   TSource const & source,
	   Generous)
{
    SEQAN_CHECKPOINT;
	target.assign(begin(source, Standard()), end(source, Standard()));
}


template <typename TChar,  typename TAlloc, typename TSource>
inline void 
assign_std_vector_Generous_impl(::std::vector<TChar,  TAlloc> & target, 
								TSource & source,
								typename Size< ::std::vector<TChar,  TAlloc> >::Type limit)
{
    SEQAN_CHECKPOINT;
	typename Iterator<TSource const, Standard>::Type source_begin = begin(source, Standard());
	typename Size<TSource const>::Type source_length = length(source);
	if (source_length > limit)
	{
		source_length = limit;
	}
	target.assign(source_begin, source_begin + source_length);
}
template <typename TChar,  typename TAlloc, typename TSource>
inline void 
assign(::std::vector<TChar,  TAlloc> & target, 
	   TSource & source,
	   typename Size< ::std::vector<TChar, TAlloc> >::Type limit,
	   Generous)
{
    SEQAN_CHECKPOINT;
	assign_std_vector_Generous_impl(target, source, limit);
}
template <typename TChar,  typename TAlloc, typename TSource>
inline void 
assign(::std::vector<TChar, TAlloc> & target, 
	   TSource const & source,
	   typename Size< ::std::vector<TChar, TAlloc> >::Type limit,
	   Generous)
{
    SEQAN_CHECKPOINT;
	assign_std_vector_Generous_impl(target, source, limit);
}

//____________________________________________________________________________

template <typename TChar, typename TAlloc, typename TSource>
inline void 
assign(::std::vector<TChar,  TAlloc> & target, 
	   TSource & source,
	   Limit)
{
    SEQAN_CHECKPOINT;
	assign(target, source, target.capacity(), Generous());
}
template <typename TChar, typename TAlloc, typename TSource>
inline void 
assign(::std::vector<TChar,  TAlloc> & target, 
	   TSource const & source,
	   Limit)
{
    SEQAN_CHECKPOINT;
	assign(target, source, target.capacity(), Generous());
}

template <typename TChar, typename TAlloc, typename TSource>
inline void 
assign(::std::vector<TChar, TAlloc> & target, 
	   TSource & source,
	   typename Size< ::std::vector<TChar, TAlloc> >::Type limit,
	   Limit)
{
    SEQAN_CHECKPOINT;
	if (limit > target.capacity()) 
	{
		limit = target.capacity();
	}

	assign(target, source, limit, Generous());
}
template <typename TChar, typename TAlloc, typename TSource>
inline void 
assign(::std::vector<TChar,  TAlloc> & target, 
	   TSource const & source,
	   typename Size< ::std::vector<TChar,  TAlloc> >::Type limit,
	   Limit)
{
    SEQAN_CHECKPOINT;
	if (limit > target.capacity()) 
	{
		limit = target.capacity();
	}

	assign(target, source, limit, Generous());
}

//////////////////////////////////////////////////////////////////////////////
//append to ::std::vector

///.Function.append.param.target.type:Adaption.std::vector
///.Function.append.param.source.type:Adaption.std::vector

template <typename TChar, typename TAlloc, typename TSource>
inline void 
append(::std::vector<TChar,  TAlloc> & target, 
	   TSource const & source,
	   Generous)
{
    SEQAN_CHECKPOINT;
	target.insert(target.end(), begin(source, Standard()), end(source, Standard()));
}

template <typename TChar,  typename TAlloc, typename TSource>
inline void 
append(::std::vector<TChar, TAlloc> & target, 
	   TSource const & source,
	   typename Size< ::std::vector<TChar, TAlloc> >::Type limit,
	   Generous)
{
    SEQAN_CHECKPOINT;
	typename Size< ::std::vector<TChar, TAlloc> >::Type target_length = target.length();
	if (target_length > limit)
	{
		target.resize(limit);
	}
	else
	{
		limit -= target_length;
		typename Iterator<TSource const, Standard>::Type source_begin = begin(source, Standard());
		typename Size<TSource const>::Type source_length = length(source);
		if (source_length > limit)
		{
			source_length = limit;
		}

		target.insert(target.end(), source_begin, source_begin + source_length);
	}
}

//____________________________________________________________________________

template <typename TChar, typename TAlloc, typename TSource>
inline void 
append(::std::vector<TChar,  TAlloc> & target, 
	   TSource const & source,
	   Limit)
{
    SEQAN_CHECKPOINT;
	append(target, source, target.capacity(), Generous());
}

template <typename TChar, typename TAlloc, typename TSource>
inline void 
append(::std::vector<TChar, TAlloc> & target, 
	   TSource const & source,
	   typename Size< ::std::vector<TChar, TAlloc> >::Type limit,
	   Limit)
{
    SEQAN_CHECKPOINT;
	if (limit > target.capacity()) 
	{
		limit = target.capacity();
	}

	append(target, source, limit, Generous());
}

//////////////////////////////////////////////////////////////////////////////
///.Function.appendValue.param.target.type:Adaption.std::vector

template <typename TChar, typename TAlloc, typename TValue, typename TTag>
inline void
appendValue(::std::vector<TChar, TAlloc> & me, 
			TValue const & _value,
			TTag)
{
    SEQAN_CHECKPOINT;
	me.push_back(_value);
} 

template <typename TChar, typename TAlloc, typename TValue>
inline void
appendValue(::std::vector<TChar,  TAlloc> & me, 
			TValue const & _value,
			Limit)
{
    SEQAN_CHECKPOINT;
	if (capacity(me) > length(me)) me.push_back(_value);
} 

//////////////////////////////////////////////////////////////////////////////
//replace to ::std::vector

///.Function.replace.param.target.type:Adaption.std::vector
///.Function.replace.param.source.type:Adaption.std::vector

/*template <typename TChar,  typename TAlloc, typename TSource>
inline void 
replace(::std::vector<TChar, TAlloc> & target,
		typename Position< ::std::vector<TChar, TAlloc> >::Type pos_begin,
		typename Position< ::std::vector<TChar, TAlloc> >::Type pos_end,
		TSource const & source,
		Generous)
{
    SEQAN_CHECKPOINT;
	copy(begin(source, Standard()),begin(source, Standard())+pos_end-pos_begin, pos_begin);
	target.insert(pos_end,begin(source,Standard())+pos_end-pos_begin,end(source,Standard());
	//target.replace(target.begin() + pos_begin, target.begin() + pos_end, begin(source, Standard()), end(source, Standard()));
}

template <typename TChar, typename TAlloc, typename TSource>
inline void 
replace(::std::vector<TChar, TAlloc> & target, 
		typename Position< ::std::vector<TChar, TAlloc> >::Type pos_begin,
		typename Position< ::std::vector<TChar, TAlloc> >::Type pos_end,
		TSource const & source,
		typename Size< ::std::vector<TChar, TAlloc> >::Type limit,
		Generous)
{
    SEQAN_CHECKPOINT;
	if (pos_begin >= limit)
	{
		target.resize(limit);
	}
	else
	{
		typename Iterator<TSource const, Standard>::Type source_begin = begin(source, Standard());
		typename Size<TSource const>::Type source_length = length(source);
		typename Size< ::std::vector<TChar, TAlloc> >::Type pos_mid = pos_begin + source_length;
		if (pos_mid > limit)
		{
			target.replace(target.begin() + pos_begin, target.begin() + limit, source_begin, source_begin + limit - pos_begin);
			target.resize(limit);
		}
		else
		{
			target.replace(target.begin() + pos_begin, target.begin() + pos_end, source_begin, end(source, Standard()));
			if (target.length() > limit)
			{
				target.resize(limit);
			}
		}
	}

}

template <typename TChar,  typename TAlloc, typename TSource>
inline void 
replace(::std::vector<TChar,  TAlloc> & target,
		typename Position< ::std::vector<TChar, TAlloc> >::Type pos_begin,
		typename Position< ::std::vector<TChar, TAlloc> >::Type pos_end,
		TSource const & source,
		Limit)
{
    SEQAN_CHECKPOINT;
//	replace(target, pos_begin, pos_end, source, target.capacity(), Generous());
}

template <typename TChar, typename TAlloc, typename TSource>
inline void 
replace(::std::vector<TChar, TAlloc> & target, 
		typename Position< ::std::vector<TChar,  TAlloc> >::Type pos_begin,
		typename Position< ::std::vector<TChar,  TAlloc> >::Type pos_end,
		TSource const & source,
		typename Size< ::std::vector<TChar, TAlloc> >::Type limit,
		Limit)
{
    SEQAN_CHECKPOINT;
	if (limit > target.capacity()) 
	{
		limit = target.capacity();
	}

//	replace(target, pos_begin, pos_end, source, limit, Generous());
}

//////////////////////////////////////////////////////////////////////////////
// handling of iterators as begin and end

template<typename TChar, typename TCharTraits, typename TAlloc, typename TSource, typename TExpand>
inline void 
replace(::std::vector<TChar, TCharTraits, TAlloc> & target,
		typename Iterator< ::std::vector<TChar, TCharTraits, TAlloc>, Rooted>::Type pos_begin,
		typename Iterator< ::std::vector<TChar, TCharTraits, TAlloc>, Rooted>::Type pos_end,
		TSource & source,
		Tag<TExpand> const tag)
{
	replace(target, position(pos_begin), position(pos_end), source, tag);
}

template<typename TChar, typename TCharTraits, typename TAlloc, typename TSource, typename TExpand>
inline void 
replace(::std::vector<TChar, TCharTraits, TAlloc> & target,
		typename Iterator< ::std::vector<TChar, TCharTraits, TAlloc>, Rooted>::Type pos_begin,
		typename Iterator< ::std::vector<TChar, TCharTraits, TAlloc>, Rooted>::Type pos_end,
		TSource & source,
		typename Size< ::std::vector<TChar, TCharTraits, TAlloc> >::Type limit,
		Tag<TExpand> const tag)
{
	replace(target,  position(pos_begin),  position(pos_end), source, tag);
}
*/

///.Function.reserve.param.object.type:Adaption.std::vector
///.Function.reserve.remarks:For @Adaption.std::vector|STL Adaptions@, $reserve$ is only guaranteed to have the specified behaviour with $Insist$ and $Generous$.
template <typename TChar,  typename TAlloc, typename TSize, typename TExpand>
inline typename Size< ::std::vector<TChar, TAlloc> >::Type 
reserve(
	::std::vector<TChar, TAlloc> & seq, 
	TSize new_capacity,
	Tag<TExpand> const & tag)
{
    SEQAN_CHECKPOINT;
    seq.reserve(new_capacity);
	return _capacityReturned(seq, new_capacity, tag);
}

template <typename TChar, typename TAlloc, typename TSize>
inline typename Size< ::std::vector<TChar, TAlloc> >::Type 
reserve(
	::std::vector<TChar, TAlloc> & seq, 
	TSize new_capacity,
	Insist const &)
{
    SEQAN_CHECKPOINT;
	// do nothing
	return _capacityReturned(seq, new_capacity, Insist());
}

template <typename TChar,  typename TAlloc, typename TSize>
inline typename Size< ::std::vector<TChar, TAlloc> >::Type 
reserve(
	::std::vector<TChar,  TAlloc> & seq, 
	TSize new_capacity,
	Limit const &)
{
    SEQAN_CHECKPOINT;
	// do nothing
	return _capacityReturned(seq, new_capacity, Limit());
}

///.Function.resize.param.object.type:Adaption.std::vector
template <typename TChar,  typename TAlloc, typename TSize, typename TExpand>
inline typename Size< ::std::vector<TChar,  TAlloc> >::Type 
resize(
	::std::vector<TChar, TAlloc> & me,
	TSize new_length,
	Tag<TExpand> const &)
{
    SEQAN_CHECKPOINT;
    me.resize(new_length);
	return me.size();
}

///.Function.fill.param.object.type:Adaption.std::vector
template <typename TChar, typename TAlloc, typename TSize, typename TExpand>
inline typename Size< ::std::vector<TChar,  TAlloc> >::Type 
fill(
	::std::vector<TChar, TAlloc> & me,
	TSize new_length,
	TChar const & val,
	Tag<TExpand> const &)
{
    SEQAN_CHECKPOINT;
    me.resize(new_length, val);
	return me.length();
}

}  // namespace seqan

#endif  // #ifndef SEQAN_HEADER_STD_STRING_H
