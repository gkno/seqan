#ifndef SEQAN_HEADER_FIND_SCORE_H
#define SEQAN_HEADER_FIND_SCORE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Score
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Score:
..cat:Pattern Matching
..general:Class.Pattern
..summary:Provides approximate string-matching with a user-definable scoring function.
..signature:Pattern<TNeedle, TScore>
..param.TNeedle:The needle type.
...type:Class.String
..param.TScore:The scoring function.
...type:Class.Score
..remarks.text:The algorithm is based on the Needleman-Wunsch dynamic progamming algorithm. The Pattern object only contains the right-most column of the DP matrix.
*/

///.Class.Pattern.param.TSpec.type:Spec.Score


template <typename TNeedle, typename TScoreValue, typename TScoreSpec>
class Pattern<TNeedle, Score<TScoreValue, TScoreSpec> >
{
public:
	typedef Score<TScoreValue, TScoreSpec> TScore;
	typedef typename Value<TScore>::Type TTabValue;

	Holder<TNeedle>		data_needle;
	Holder<TScore>		data_score;
	TTabValue			data_limit;
	String<TTabValue>	data_tab;

public: 
	Pattern(): 
		data_limit(0)
	{ 
SEQAN_CHECKPOINT
		create(data_score);
	}

	Pattern(TNeedle & _needle, 
			TScore & _score_func, 
			TTabValue _limit = 0): 
		data_score(_score_func),
		data_limit(_limit)
	{ 
SEQAN_CHECKPOINT
		setHost(*this, _needle);
	}

	Pattern(TNeedle & _needle,
			TTabValue _limit = 0): 
		data_limit(_limit)
	{ 
SEQAN_CHECKPOINT
		create(data_score);
		setHost(*this, _needle);
	}

	Pattern(TTabValue _limit): 
		data_limit(_limit)
	{ 
SEQAN_CHECKPOINT
		create(data_score);
	}

	Pattern(Pattern const & other): 
		data_needle( other.data_needle ),
		data_score( other.data_score ), 
		data_limit( other.data_limit ),
		data_tab( other.data_tab )
	{
SEQAN_CHECKPOINT
	}

	inline Pattern & 
	operator = (Pattern const & other) 
	{ 
SEQAN_CHECKPOINT
		this->data_needle = other.data_needle;
		this->data_score = other.data_score;
		this->data_limit = other.data_limit;
		this->data_tab = other.data_tab;

		return *this;
	}

//____________________________________________________________________________

	friend inline typename Host<Pattern const>::Type & 
	host(Pattern & me)
	{
SEQAN_CHECKPOINT
		return value(me.data_needle);
	}

	friend inline typename Host<Pattern const>::Type & 
	host(Pattern const & me)
	{
SEQAN_CHECKPOINT
		return value(me.data_needle);
	}

//____________________________________________________________________________

	friend inline TScore const & 
	scoring(Pattern & me)
	{
SEQAN_CHECKPOINT
		return value(me.data_score);
	}

	friend inline void
	setScoring(Pattern const & me, Score<TScoreValue, TScoreSpec> const & score)
	{
SEQAN_CHECKPOINT
		me.data_score = score;
	}

//____________________________________________________________________________

	friend inline TScoreValue 
	limit(Pattern & me)
	{
SEQAN_CHECKPOINT
		return me.data_limit;
	}
	friend inline void 
	setLimit(Pattern & me, TScoreValue _limit)
	{
SEQAN_CHECKPOINT
		me.data_limit = _limit;
	}

//____________________________________________________________________________

	friend inline typename Host<Pattern>::Type & 
	_dataNeedle(Pattern & me)
	{
SEQAN_CHECKPOINT
		return value(me.data_needle);
	}

	friend inline String<TTabValue> & 
	_dataTab(Pattern & me)
	{
SEQAN_CHECKPOINT
		return me.data_tab;
	}

//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// Host
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TScoreValue, typename TScoreSpec>
struct Host< Pattern<TNeedle, Score<TScoreValue, TScoreSpec> > >
{
	typedef TNeedle Type;
};

template <typename TNeedle, typename TScoreValue, typename TScoreSpec>
struct Host< Pattern<TNeedle, Score<TScoreValue, TScoreSpec> > const>
{
	typedef TNeedle const Type;
};


//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TNeedle2, typename TScoreValue>
void 
setHost(Pattern<TNeedle, Score<TScoreValue, Simple> > & me, TNeedle2 const & ndl)
{
	typedef String<TScoreValue> TTab;
	typedef typename Iterator<TTab, Standard>::Type TIterator;

	TScoreValue score_gap = scoreGapExtend(scoring(me));

	TTab & string_tab = _dataTab(me);

	//allocate enough memory for one column of DP matrix
	/*size_t got_length =*/ resize(string_tab, length(ndl));

//	if (length(ndle) < got_length) throw(0); //???TODO: Throw "not enough memory" exception

	//init matrix
	//note: The column is stored in reverse order
	TIterator tab_end = begin(string_tab, Standard());
	TIterator tab = end(string_tab, Standard());
	
	TScoreValue x = score_gap;

	while (tab > tab_end)
	{
		--tab;
		*tab = x;
		x += score_gap;
	}

	_dataNeedle(me) = ndl;
}

template <typename TNeedle, typename TNeedle2, typename TScoreValue>
void 
setHost(Pattern<TNeedle, Score<TScoreValue, Simple> > & me, TNeedle2 & ndl)
{
	setHost(me, reinterpret_cast<TNeedle const &>(ndl));
}

//////////////////////////////////////////////////////////////////////////////
// returns the score of the last hit position found (note:position = end of occurrence in haystack)

template <typename TNeedle, typename TScoreValue>
inline TScoreValue
getScore(Pattern<TNeedle, Score<TScoreValue, Simple> > & me)
{
	return front(me.data_tab);
}


//////////////////////////////////////////////////////////////////////////////
// find, findNext
//////////////////////////////////////////////////////////////////////////////

//proportional gap cost: Needleman-Wunsch 

//???TODO: Ukkonen trick?

template <typename TFinder, typename TNeedle, typename TScoreValue>
bool 
_find_score_simple_proportional(TFinder & finder, Pattern<TNeedle, Score<TScoreValue, Simple> > & me)
{
	typedef String<TScoreValue> TTab;
	typedef typename Iterator<TTab, Standard>::Type TTabIterator;
	typedef typename Iterator<TNeedle const, Standard>::Type TNeedleIterator;
	typedef typename Value<typename Haystack<TFinder>::Type>::Type THaystackValue;

	String<TScoreValue> & string_tab = _dataTab(me);

	TScoreValue score_gap = scoreGapExtend(scoring(me));
	TScoreValue score_match = scoreMatch(scoring(me));
	TScoreValue score_mismatch = scoreMismatch(scoring(me));

	//init table

	if (empty(finder))
		goBegin(finder);
	else
		goNext(finder);

	//start searching

	TTabIterator tab_begin = end(string_tab, Standard());

	TNeedleIterator it_begin = begin(host(me), Standard());
	TNeedleIterator it_end = end(host(me), Standard());

	//for each character in haystack, do...
	for (; !atEnd(finder); ++finder)
	{
		//get character
		THaystackValue c = *finder;

		//init some variables
		TNeedleIterator it = it_begin;
		TScoreValue * tab = tab_begin;
		TScoreValue h = 0;
		TScoreValue v = 0;

		//fill the column
		while (it < it_end)
		{
			--tab; //note: the column is stored in "reverse order"

			TScoreValue m2 = (c == *it) ? h + score_match : h + score_mismatch;
			h = *tab;
			TScoreValue m1 = (h > v) ? h + score_gap : v + score_gap;

			v = (m1 > m2) ? m1 : m2;
			*tab = v;

			++it;
		}

		if (*tab >= limit(me) )
		{//found a hit
			return true;
		}

	}

	//found nothing
	return false;
}

//____________________________________________________________________________

template <typename TFinder, typename TNeedle, typename TScoreValue>
bool 
_find_score_simple_affine(TFinder & finder, Pattern<TNeedle, Score<TScoreValue, Simple> > & me)
{
//???TODO 
	return false;
}


//////////////////////////////////////////////////////////////////////////////

template <typename TFinder, typename TNeedle, typename TScoreValue>
inline bool 
find(TFinder & finder, Pattern<TNeedle, Score<TScoreValue, Simple> > & me)
{
	if (scoreGapOpen(scoring(me)))
	{
		return _find_score_simple_affine(finder, me);
	}
	else
	{
		return _find_score_simple_proportional(finder, me);
	}
}


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
