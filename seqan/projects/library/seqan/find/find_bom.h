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
..cat:Searching
..signature:Pattern<TNeedle, BomAlgo>
..param.TNeedle:The needle type.
...type:Class.String
..remarks.text:The types of the needle and the haystack have to match.
*/

///.Class.Pattern.param.TSpec.type:Spec.BomAlgo

template <typename TSpec>
struct BFA; //backward factor automaton searching

struct Oracle; //Oracle Tag => "BOM"
struct Trie; //Trie Tag => "BTM"

typedef BFA<Oracle> BomAlgo; //for compatibility reasons

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
class Pattern<TNeedle, BFA<TSpec> > {
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
	Graph<Automaton<TAlphabet, void, WithoutEdgeId> > oracle;

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

//BFA<Oracle>: BOM Algorithm
template <typename TNeedle, typename TNeedle2>
inline void 
setHost (Pattern<TNeedle, BFA<Oracle> > & me, TNeedle2 const& needle) 
{
	SEQAN_CHECKPOINT
	me.needleLength = length(needle);
	clear(me.oracle);
	createOracleOnReverse(me.oracle,needle);
	assignRoot(me.oracle,0);
	setValue(me.data_needle, needle);
}
//BFA<Trie>: BTM Algorithm (the same as BOM, but with an trie)
template <typename TNeedle, typename TNeedle2>
inline void 
setHost (Pattern<TNeedle, BFA<Trie> > & me, TNeedle2 const& needle) 
{
	SEQAN_CHECKPOINT
	me.needleLength = length(needle);
	clear(me.oracle);

	String<String<unsigned int> > terminal_state_map; //dummy
	typedef typename Value<TNeedle2 const>::Type TValue;
	String<TValue> reverse_string = needle;
	reverseInPlace(reverse_string);

	setValue(me.data_needle, needle);
}

template <typename TNeedle, typename TNeedle2, typename TSpec>
inline void 
setHost (Pattern<TNeedle, BFA<TSpec> > & me, TNeedle2 & needle)
{
	setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//____________________________________________________________________________


template <typename TNeedle, typename TSpec>
inline void _patternInit (Pattern<TNeedle, BFA<TSpec> > & me) 
{
SEQAN_CHECKPOINT
	me.step = 0;
}


//____________________________________________________________________________

template <typename TNeedle, typename TSpec>
inline typename Host<Pattern<TNeedle, BFA<TSpec> > const>::Type & 
host(Pattern<TNeedle, BFA<TSpec> > & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

template <typename TNeedle, typename TSpec>
inline typename Host<Pattern<TNeedle, BFA<TSpec> > const>::Type & 
host(Pattern<TNeedle, BFA<TSpec> > const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

//____________________________________________________________________________


template <typename TFinder, typename TNeedle, typename TSpec>
inline bool 
find(TFinder & finder, Pattern<TNeedle, BFA<TSpec> > & me) 
{
	SEQAN_CHECKPOINT
	
	if (empty(finder)) {
		_patternInit(me);
		_finderSetNonEmpty(finder);
		me.haystackLength = length(container(finder));
	} else
		finder+=me.step;

	if (me.haystackLength < me.needleLength) return false;
	typedef typename Value<TNeedle>::Type TAlphabet;
	typedef Graph<Automaton<TAlphabet> > TOracle;
	typedef typename Size<TNeedle>::Type TSize;
	typedef typename VertexDescriptor<TOracle>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TOracle>::Type TEdgeDescriptor;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	while (position(finder) <= me.haystackLength - me.needleLength) {
		TVertexDescriptor current = getRoot(me.oracle);
		TSize j = me.needleLength;
		while ((j>0) &&
				(current != nilVal))
		{
			TAlphabet c = *(finder+(j-1));
			current = targetVertex(me.oracle, findEdge(me.oracle, current, c));
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
