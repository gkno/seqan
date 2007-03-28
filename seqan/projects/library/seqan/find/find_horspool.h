#ifndef SEQAN_HEADER_FIND_HORSPOOL_H
#define SEQAN_HEADER_FIND_HORSPOOL_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Horspool
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Horspool:
..summary: Exact string matching using Horspool's algorithm (1980).
..general:Class.Pattern
..cat:Pattern Matching
..signature:Pattern<TNeedle, Horspool>
..param.TNeedle:The needle type.
...type:Class.String
*/

///.Class.Pattern.param.TSpec.type:Spec.Horspool

struct _Horspool;
typedef Tag<_Horspool> Horspool;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, Horspool>
{
//____________________________________________________________________________

public:
	Holder<TNeedle>							data_needle;
	String<typename Size<TNeedle>::Type>	data_map;

//____________________________________________________________________________

public:
	Pattern() {}

	Pattern(Pattern const & other_):
		data_map(other_.data_map) {}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl)
	{
		setHost(*this, ndl);
	}

	Pattern const &
	operator = (Pattern const & other_)
	{
		data_map = other_.data_map;
		return *this;
	}
//____________________________________________________________________________
};


template <typename TNeedle, typename TNeedle2>
void
setHost(Pattern<TNeedle, Horspool> & me, TNeedle2 const & ndl)
{
	typedef typename Value<TNeedle>::Type TValue;

	SEQAN_ASSERT(!empty(ndl));

	typename Size<TNeedle>::Type value_size = ValueSize<TValue>::VALUE;

	//make room for map
	resize(me.data_map, value_size);

	//fill map
	typename Size<TNeedle>::Type jump_width = length(ndl);

	arrayFill(begin(me.data_map), begin(me.data_map) + value_size, jump_width);

	typename Iterator<TNeedle2 const, Standard>::Type it;
	it = begin(ndl, Standard());
	while (jump_width > 1)
	{
		--jump_width;
		unsigned int pos_ = *it; //conversion value type to unsigned int
		me.data_map[pos_] = jump_width;
		++it;
	}

	me.data_needle = ndl;
}

template <typename TNeedle, typename TNeedle2>
void
setHost(Pattern<TNeedle, Horspool> & horsp, TNeedle2 & ndl)
{
	setHost(horsp, reinterpret_cast<TNeedle2 const &>(ndl));
}

//____________________________________________________________________________


template <typename TNeedle>
inline void _finderInit (Pattern<TNeedle, Horspool> & me) 
{
}

//____________________________________________________________________________

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, Horspool>const>::Type & 
host(Pattern<TNeedle, Horspool> & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, Horspool>const>::Type & 
host(Pattern<TNeedle, Horspool> const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

//____________________________________________________________________________

template <typename TFinder, typename TNeedle>
bool
find(TFinder & finder, Pattern<TNeedle, Horspool> & me)
{
SEQAN_CHECKPOINT
	typedef typename Haystack<TFinder>::Type THaystack;
	THaystack & haystack = container(finder);
	typename Position<TFinder>::Type pos = position(finder);
	typename Size<THaystack>::Type hstk_size = length(haystack);
	typename Size<TNeedle>::Type ndl_size = length(host(me));

	if (ndl_size > hstk_size)
	{//needle larger than haystack: nothing to find
		return false;
	}

	typename Position<THaystack>::Type pos_max = hstk_size - ndl_size;
	SEQAN_ASSERT2(pos <= pos_max, "invalid search position")

	//helper variables
	//typename Size<TFinder>::Type window_length = ndl_size - 1;
	typename Position<TFinder>::Type new_pos;

	if (empty(finder))
	{
		_finderInit(me);
		_finderSetNonEmpty(finder);
		goto VALIDATE;
	}

	unsigned int char_i;
	while (true)
	{
		//move to next position
		char_i = *(finder + ndl_size - 1); //conversion to unsigned integer
		new_pos = pos + me.data_map[char_i];
		if ((new_pos > pos_max) || (new_pos < pos))
		{//pos out of range: found nothing
			return false;
		}

		pos = new_pos;
		finder += me.data_map[char_i];

VALIDATE:
		if (infix(haystack, pos, pos + ndl_size) == host(me))
		{//found a hit
			return true;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// Host
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
struct Host< Pattern<TNeedle, Horspool> >
{
	typedef TNeedle Type;
};

template <typename TNeedle>
struct Host< Pattern<TNeedle, Horspool> const>
{
	typedef TNeedle const Type;
};


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
