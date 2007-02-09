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
	typedef typename Size<TNeedle>::Type TSize;
	typedef typename Value<TNeedle>::Type TValue;
	typedef typename Value<TValue>::Type TAlphabet;
	Holder<TNeedle> data_needle;
	Graph<Automaton<TAlphabet> > data_reverseTrie;  // Search trie
	String<String<TSize> > data_terminalStateMap;
	String<TSize> data_dMap;	// Jump table
	TSize data_lmin;
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
	typedef typename Value<TNeedle>::Type TKeyword;
	typedef typename Size<TKeyword>::Type TSize;
	typedef typename Value<TKeyword>::Type TAlphabet;

	// clean-up
	clear(me.data_reverseTrie);
	clear(me.data_terminalStateMap);
	clear(me.data_dMap);

	// Create Trie
	createTrieOnReverse(me.data_reverseTrie,me.data_terminalStateMap,needle);
	assignRoot(me.data_reverseTrie,0);
	me.data_needle = needle;

	// Create jump map
	TSize alphabet_size = ValueSize<TAlphabet>::VALUE;
	resize(me.data_dMap, alphabet_size);
	me.data_lmin = _get_infinity<TSize>();
	typename Iterator<TNeedle2 const>::Type it = begin(needle);
	for(;!atEnd(it);goNext(it)) {
		TSize tmp = length(*it)-1;
		if (tmp<me.data_lmin) me.data_lmin = tmp;
	}
	for(TSize i=0;i<alphabet_size;++i) {
		me.data_dMap[i]=me.data_lmin;
	}
	goBegin(it);
	for(;!atEnd(it);goNext(it)) {
		for(TSize pos = 0;pos < length(*it) - 1; ++pos) {
			TSize ind = convert<TSize>((*it)[pos]);	
			if ((length(*it)- 1 - pos) < me.data_dMap[ind]) {
				me.data_dMap[ind] = (length(*it) - 1 - pos);
			}
		}
	}

	/*
	fstream strm;
	strm.open(TEST_PATH "my_trie.dot", ios_base::out | ios_base::trunc);
	String<String<char> > nodeMap;
	_createTrieNodeNames(me.data_reverseTrie, me.data_terminalStateMap, nodeMap);
	String<String<char> > edgeMap;
	_createEdgeNames(me.data_reverseTrie,edgeMap);
	write(strm,me.data_reverseTrie,nodeMap,edgeMap,DotDrawing());
	strm.close();
	*/
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
