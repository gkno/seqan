#ifndef SEQAN_HEADER_FIND_SETHORSPOOL_H
#define SEQAN_HEADER_FIND_SETHORSPOOL_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Set Horspool Algorithm
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.SetHorspool:
..summary: Multiple exact string matching using set horspool algorithm.
..general:Class.Pattern
..cat:Pattern Matching
..signature:Pattern<TNeedle, SetHorspool>
..param.TNeedle:The needle type, a string of keywords.
...type:Class.String
..remarks.text:The types of all keywords in the needle and the haystack have to match.
*/

///.Class.Pattern.param.TSpec.type:Spec.SetHorspool

struct _SetHorspool;
typedef Tag<_SetHorspool> SetHorspool;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, SetHorspool> {
//____________________________________________________________________________
private:
	Pattern(Pattern const& other);
	Pattern const& operator=(Pattern const & other);

//____________________________________________________________________________
public:
	typedef unsigned int TWord;
	typedef typename Size<TNeedle>::Type TSize;
	Holder<TNeedle> data_needle;
	TSize data_keywordIndex;			// Current keyword that produced a hit
	TSize data_needleLength;			// Last length of needle to reposition finder

//____________________________________________________________________________

	Pattern() {
	}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl)
	{
		setHost(*this, ndl);
	}

	~Pattern() {
		SEQAN_CHECKPOINT
	}		
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// Host Metafunctions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
struct Host< Pattern<TNeedle, SetHorspool> >
{
	typedef TNeedle Type;
};

template <typename TNeedle>
struct Host< Pattern<TNeedle, SetHorspool> const>
{
	typedef TNeedle const Type;
};

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, SetHorspool> & me, TNeedle2 const & needle) {
	SEQAN_CHECKPOINT
	me.data_needle = needle;
}

template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, SetHorspool> & me, TNeedle2 & needle)
{
	setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//____________________________________________________________________________

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, SetHorspool>const>::Type & 
host(Pattern<TNeedle, SetHorspool> & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, SetHorspool>const>::Type & 
host(Pattern<TNeedle, SetHorspool> const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

//____________________________________________________________________________


template <typename TNeedle>
inline typename Size<TNeedle>::Type
position(Pattern<TNeedle, SetHorspool> & me)
{
	return me.data_keywordIndex;
}


template <typename TFinder, typename TNeedle>
inline bool find(TFinder & finder, Pattern<TNeedle, SetHorspool> & me) {
	SEQAN_CHECKPOINT
	if (empty(finder))
		goBegin(finder);
	else
		finder += me.data_needleLength;


	return false;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_SETHORSPOOL_H
