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
inline void _finderInit (Pattern<TNeedle, Horspool> &) {}

//____________________________________________________________________________

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, Horspool> >::Type & 
host(Pattern<TNeedle, Horspool> & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, Horspool> const>::Type & 
host(Pattern<TNeedle, Horspool> const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

//____________________________________________________________________________

template <typename TFinder, typename TNeedle2>
bool
find_horspool(TFinder & finder, 
			  Pattern<TNeedle2, Horspool> & me,
			  bool find_first)
{
SEQAN_CHECKPOINT
	typedef typename Haystack<TFinder>::Type THaystack;
	THaystack & hayst = haystack(finder);

	typedef Pattern<TNeedle2, Horspool> TPattern;
	typedef typename Needle<TPattern>::Type TNeedle;
	TNeedle & ndl = needle(me);

	typedef typename Size<TNeedle>::Type TNeedleSize;
	TNeedleSize ndl_size = length(ndl);

	typedef typename Iterator<THaystack, Standard>::Type THaystackIterator;
	THaystackIterator haystack_end = end(hayst, Standard());
	THaystackIterator it = begin(hayst, Standard());
	it += position(finder) + ndl_size - 1; //it points to the last character
	THaystackIterator it_next = it;

	typedef typename Iterator<TNeedle, Standard>::Type TNeedleIterator;
	TNeedleIterator nit; //needle iterator
	TNeedleIterator nit_begin = begin(ndl, Standard());
	TNeedleIterator nit_end = end(ndl, Standard()) - 1; //here the verification begins

	unsigned int char_i;

	if (find_first)
	{
		goto VALIDATE;
	}

MOVE_FURTHER:
	//move to next position
	char_i = *it; //conversion to unsigned integer
	it_next = it + me.data_map[char_i];
	if (it_next >= haystack_end)
	{//found nothing
		return false;
	}

	it = it_next;

VALIDATE:
	//validate current position
	for (nit = nit_end; nit >= nit_begin; --nit)
	{
		if (*nit != *it_next)
		{//invalid!
			goto MOVE_FURTHER;
		}
		--it_next;
	}

	//valid! return hit
	setPosition(finder, it - begin(hayst, Standard()) - ndl_size + 1);
	return true;
}

//____________________________________________________________________________
//spec for file reader haystacks

template <typename TFormat, typename TFile, typename TSpec>
struct FileReader;

template <typename TValue, typename TFormat, typename TFile, typename TSpec, typename TNeedle2>
bool
find_horspool(Finder<String<TValue, FileReader<TFormat, TFile, TSpec> > > & finder, 
			  Pattern<TNeedle2, Horspool> & me,
			  bool find_first)
{
SEQAN_CHECKPOINT
	typedef Finder<String<TValue, FileReader<TFormat, TFile, TSpec> > > TFinder;
	typedef typename Haystack<TFinder>::Type THaystack;
	THaystack & hayst = haystack(finder);

	typedef Pattern<TNeedle2, Horspool> TPattern;
	typedef typename Needle<TPattern>::Type TNeedle;
	TNeedle & ndl = needle(me);

	typedef typename Size<TNeedle>::Type TNeedleSize;
	TNeedleSize ndl_size = length(ndl);

	typedef typename Iterator<THaystack, Standard>::Type THaystackIterator;
	THaystackIterator it(hayst, position(finder) + ndl_size - 1); //it points to the last character

	typedef typename Iterator<TNeedle, Standard>::Type TNeedleIterator;
	TNeedleIterator nit; //needle iterator
	TNeedleIterator nit_begin = begin(ndl, Standard());
	TNeedleIterator nit_end = end(ndl, Standard()) - 1; //here the verification begins

	unsigned int char_i;

	if (find_first)
	{
		goto VALIDATE;
	}

MOVE_FURTHER:
	//move to next position
	char_i = *it; //conversion to unsigned integer
	it += me.data_map[char_i];
	if (atEnd(it))
	{//found nothing
		return false;
	}

VALIDATE:
	//validate current position
	for (nit = nit_end; nit >= nit_begin; --nit)
	{
		if (*nit != *it)
		{//invalid!
			it += (nit_end - nit);
			goto MOVE_FURTHER;
		}
		--it;
	}

	//valid! return hit
	setPosition(finder, it - begin(hayst, Standard()) + 1);
	return true;

}

//____________________________________________________________________________
/* groepl variante

template <typename TFinder, typename TNeedle2>
bool
find_horspool(TFinder & finder, 
	Pattern<TNeedle2, Horspool> & me,
	bool find_first)
{
SEQAN_CHECKPOINT
	typedef typename Haystack<TFinder>::Type THaystack;
	THaystack & hayst = haystack(finder);

	typedef Pattern<TNeedle2, Horspool> TPattern;
	typedef typename Needle<TPattern>::Type TNeedle;
	TNeedle & ndl = needle(me);

	typename Size<TNeedle>::Type ndl_size = length(ndl);

	typedef typename Iterator<THaystack, Standard>::Type THaystackIterator;
	THaystackIterator haystack_end = end(hayst, Standard());
	THaystackIterator it = begin(hayst, Standard());
	it += position(finder) + ndl_size - 1; //it points to the last character
	THaystackIterator it2;

	typedef typename Iterator<TNeedle, Standard>::Type TNeedleIterator;
	TNeedleIterator nit; //needle iterator
	TNeedleIterator nit_begin = begin(ndl, Standard());
	TNeedleIterator nit_end = end(ndl, Standard()) - 2; //here the verification begins

	typedef typename Value<TNeedle>::Type TNeedleValue;
	TNeedleValue last_needle_char = value(ndl, ndl_size - 1);

	unsigned int char_i;

	if (!find_first)
	{
		++it;
	}

MOVE_FURTHER:
	//scan for the last character
	while (true)
	{
		if (it >= haystack_end)
		{//found nothing
			return false;
		}
		if (*it == last_needle_char) 
		{
			break;
		}
		++it;
	}

VALIDATE:
	it2 = it;
	//validate current position
	for (nit = nit_end; nit >= nit_begin; --nit)
	{
		--it2;
		if (*nit != *it2)
		{//invalid! skip!
			char_i = *it; //conversion to unsigned integer
			it += me.data_map[char_i];
			goto MOVE_FURTHER;
		}
	}

	//valid! return hit
	setPosition(finder, it - begin(hayst, Standard()) - ndl_size + 1);
	return true;
}
*/

template <typename TFinder, typename TNeedle2>
bool
find(TFinder & finder, Pattern<TNeedle2, Horspool> & me)
{
SEQAN_CHECKPOINT
	bool find_first = empty(finder);
	if (find_first)
	{
		_finderInit(me);
		_finderSetNonEmpty(finder);
	}

	SEQAN_ASSERT(length(needle(me)) > 0)

	return find_horspool(finder, me, find_first);
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
