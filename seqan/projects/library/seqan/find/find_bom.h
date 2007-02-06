#ifndef SEQAN_HEADER_FIND_BOM_H
#define SEQAN_HEADER_FIND_BOM_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// BomAlgo
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.BomAlgo:
..summary: Backward Oracle Matching algorithm. Exact string matching using a factor oracle.
..general:Class.Pattern
..cat:Pattern Matching
..signature:Pattern<TNeedle, BomAlgo>
..param.TNeedle:The needle type.
...type:Class.String
..remarks.text:The types of the needle and the haystack have to match.
*/

///.Class.Pattern.param.TSpec.type:Spec.BomAlgo

struct _BomAlgo;
typedef Tag<_BomAlgo> BomAlgo;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, BomAlgo> {
//____________________________________________________________________________
private:
	Pattern(Pattern const& other);
	Pattern const& operator=(Pattern const & other);

//____________________________________________________________________________
public:
	typedef typename Value<TNeedle>::Type TAlphabet;
	typedef typename Size<TNeedle>::Type TSize;
	Holder<TNeedle> data_needle;
	TSize needleLength;		
	TSize haystackLength;
	TSize step;
	Graph<Automaton<TAlphabet> > oracle;

//____________________________________________________________________________

	Pattern() {	
	}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl)
	{
SEQAN_CHECKPOINT
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
struct Host< Pattern<TNeedle, BomAlgo> >
{
	typedef TNeedle Type;
};

template <typename TNeedle>
struct Host< Pattern<TNeedle, BomAlgo> const>
{
	typedef TNeedle const Type;
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TNeedle2>
inline void 
setHost (Pattern<TNeedle, BomAlgo> & me, TNeedle2 const& needle) 
{
	SEQAN_CHECKPOINT
	me.needleLength = length(needle);
	me.step = 0;
	clear(me.oracle);
	createOracleOnReverse(me.oracle,needle);
	assignRoot(me.oracle,0);
	me.data_needle = needle;
}

template <typename TNeedle, typename TNeedle2>
inline void 
setHost (Pattern<TNeedle, BomAlgo> & me, TNeedle2 & needle)
{
	setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//____________________________________________________________________________

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, BomAlgo>const>::Type & 
host(Pattern<TNeedle, BomAlgo> & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, BomAlgo>const>::Type & 
host(Pattern<TNeedle, BomAlgo> const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

//____________________________________________________________________________


template <typename TFinder, typename TNeedle>
inline bool 
find(TFinder & finder, Pattern<TNeedle, BomAlgo> & me) 
{
	SEQAN_CHECKPOINT
	
	if (empty(finder)) {
		goBegin(finder);
		me.haystackLength = length(container(finder));
	} else
		finder+=me.step;

	if (me.haystackLength < me.needleLength) return false;
	typedef typename Value<TNeedle>::Type TAlphabet;
	typedef Graph<Automaton<TAlphabet> > TOracle;
	typedef typename Size<TNeedle>::Type TSize;
	typedef typename VertexDescriptor<TOracle>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TOracle>::Type TEdgeDescriptor;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	while (position(finder) <= me.haystackLength - me.needleLength) {
		TVertexDescriptor current = getRoot(me.oracle);
		TSize j = me.needleLength;
		while ((j>0) &&
				(current != nilVal))
		{
			TAlphabet c = *(finder+(j-1));
			current = targetVertex(me.oracle,TEdgeDescriptor(current,c));
			--j;
		}
		if (current != nilVal) {
			me.step = j + 1;
			return true;
		}
		finder += j + 1;
	}
	return false;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_SHIFTAND_H
