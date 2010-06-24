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
  Author: Martin Riese
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_MYERS_UKKONEN2_H
#define SEQAN_HEADER_FIND_MYERS_UKKONEN2_H

namespace SEQAN_NAMESPACE_MAIN 
{

//////////////////////////////////////////////////////////////////////////////
// MyersUkkonen
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Myers:
..cat:Pattern Matching
..general:Class.Pattern
..summary:Provides fast approximate searching of one string in another using Myer's fast bit-parallel algorithm with application of the Ukkonen-trick.
..signature:Pattern<TNeedle, Myers< [TSpec [, TFindBeginPatternSpec] ]> >
..param.TNeedle:The needle type.
...type:Class.String
..param.TSpec:Specialization tag.
...default:Tag.FindInfix
...remarks:This could be @Tag.FindInfix@ for infix search or $FindPrefix$ for prefix search.
..param.TFindBeginPatternSpec:Specialization of @Class.Pattern@ used to find the begin of matches.
...default:@Metafunction.DefaultFindBeginPatternSpec@
...metafunctin:@Metafunction.FindBeginPatternSpec@
...remarks:This must be a finder for prefix search, e.g. @Spec.DPSearch|$DPSearch<TScore, FindPrefix>$@ or @Spec.Myers|$Myers<FindPrefix>$@.
Specify $void$ to suppress prefix searching.
..remarks.text:The needle-length must be smaller than the highest number that can be stored in an unsigned int.
..include:seqan/find.h
*/

///.Class.Pattern.param.TSpec.type:Spec.Myers

template <typename TSpec = FindInfix, 
		  typename THasState = True,
		  typename TFindBeginPatternSpec = typename DefaultFindBeginPatternSpec<EditDistanceScore, THasState>::Type>
struct Myers {};

//FindInfix and FindPrefix are defined int find_base.h
struct AlignTextBanded; // search query in a parallelogram

// TODO(holtgrew): Really deprecated?
//deprecated shortcuts:

/**
.Shortcut.MyersUkkonen:
..cat:Pattern Matching
..summary:Semin-global (query-global, text-local) pattern matching without findBegin() support.
..signature:MyersUkkonen
..shortcutfor:Spec.Myers
...signature:Myers<FindInfix, void>
..see:Spec.Myers
..see:Shortcut.MyersUkkonenGlobal
..see:Shortcut.MyersUkkonenBanded
*/
typedef Myers<FindInfix, True, void> MyersUkkonen;


/**
.Shortcut.MyersUkkonenGlobal:
..cat:Pattern Matching
..summary:Global (query-global, text-global) pattern matching without findBegin() support.
..signature:MyersUkkonen
..shortcutfor:Spec.Myers
...signature:Myers<FindPrefix, void>
..see:Spec.Myers
..see:Shortcut.MyersUkkonen
..see:Shortcut.MyersUkkonenBanded
*/
typedef Myers<FindPrefix, True, void> MyersUkkonenGlobal;


/**
.Shortcut.MyersUkkonenBanded:
..cat:Pattern Matching
..summary:Semin-global (query-global, text-local) pattern matching without findBegin() support.
..signature:MyersUkkonen
..shortcutfor:Spec.Myers
...signature:Myers<AlignedTextBanded, void>
..see:Spec.Myers
..see:Shortcut.MyersUkkonen
..see:Shortcut.MyersUkkonenGlobal
*/
typedef Myers<AlignTextBanded, True, void> MyersUkkonenBanded;


//____________________________________________________________________________
// bit 0 of the HP bit-vector
// 0 for begin-gap-free haystack
// 1 for global alignments of haystack

template <typename T>
struct _MyersUkkonenHP0 {
	enum { VALUE = 0 };
};

template <>
struct _MyersUkkonenHP0<FindPrefix> {
	enum { VALUE = 1 };
};


//////////////////////////////////////////////////////////////////////////////
//overwrite _FindBegin to define host member if find begin is switched on

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TFindBeginPatternSpec2>
struct _FindBegin< Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >, TFindBeginPatternSpec2>
{
private:
	typedef Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > TPattern;
	typedef typename FindBeginPattern<TPattern>::Type TFindBeginPattern;

public:
	TFindBeginPattern data_findBeginPattern;
 	Holder<TNeedle>	data_host;	//defines the 
	typedef True HasHost;
};

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
struct _FindBegin< Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >, void>
{
	typedef False HasHost;
//need no findBegin if FindBeginPatternSpec is void
};

//////////////////////////////////////////////////////////////////////////////

 
// large basic pattern
struct _MyersLargePattern
{
#ifdef SEQAN_SSE2_INT128
	typedef SSE2_int128 TWord;
#else
	typedef unsigned long TWord;
#endif

	unsigned blockCount;		// the number of blocks
	TWord finalScoreMask;		// a mask with a bit set on the position of the last row
};


// large state
struct _MyersLargeState
{
#ifdef SEQAN_SSE2_INT128
	typedef SSE2_int128 TWord;
#else
	typedef unsigned long TWord;
#endif
	unsigned lastBlock;			// the block containing the last active cell
	String<TWord> VP;
	String<TWord> VN;
	TWord scoreMask;			// the mask with a bit set at the position of the last active cell
};

// TODO: should go elsewhere
// template <typename TNeedle, typename TSpec>
// class Pattern{};
// TODO: should go elsewhere
template <typename TSpec>
class _PatternState{};


template<typename TSpec, typename TFindBeginPatternSpec>
class _PatternState<Myers<TSpec, False, TFindBeginPatternSpec> >
{
public:
	  _PatternState() {}
};


template<typename TSpec, typename TFindBeginPatternSpec>
class _PatternState<Myers<TSpec, True, TFindBeginPatternSpec> >{
	
public:
#ifdef SEQAN_SSE2_INT128
	typedef SSE2_int128 TWord;
#else
	typedef unsigned long TWord;
#endif

	unsigned score;				// the current score

	TWord VP0;					// VP[0] (saves one dereferentiation)
	TWord VN0;					// VN[0]

	_MyersLargeState *largeState;

//____________________________________________________________________________

	_PatternState():
		score(0),
		VP0(0),
		VN0(0),
		largeState(NULL)
	{}

	_PatternState(_PatternState const & other):
		score(other.score),
		VP0(other.VP0),
		VN0(other.VN0),
		largeState(NULL)
	{
		if (other.largeState)
		{
			largeState = new _MyersLargeState;
			(*largeState) = *(other.largeState);
		}
	}

	~_PatternState()
	{
		delete largeState;
	}
	
	_PatternState &
	operator = (_PatternState const & other)
	{
		score = other.score;
		VP0 = other.VP0;
		VN0 = other.VN0;
		if (other.largeState)
		{
			if (largeState == NULL)
				largeState = new _MyersLargeState;
			(*largeState) = *(other.largeState);
		} else
			delete largeState;
		return *this;
	}
};


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
class Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >:
	public _FindBegin<Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > >,
	public _PatternState<Myers<TSpec, THasState, TFindBeginPatternSpec> >
{
	
public:
#ifdef SEQAN_SSE2_INT128
	typedef SSE2_int128 TWord;
#else
	typedef unsigned long TWord;
#endif
	typedef typename Iterator<TNeedle, Standard>::Type TIter;
	enum { MACHINE_WORD_SIZE = sizeof(TWord) * 8 };

	typedef _PatternState<Myers<TSpec, THasState, TFindBeginPatternSpec> > TPatternState;

	unsigned needleSize;
	unsigned k;					// the maximal number of differences allowed
	String<TWord> bitMasks;		// encoding the alphabet as bit-masks

	Holder<TNeedle> data_host;

	_MyersLargePattern *largePattern;	// extra preprocessing info for large patterns

//____________________________________________________________________________

	Pattern(int _limit = -1) :
	  TPatternState(),
		needleSize(0),
		k(- _limit),
		largePattern(NULL)
	{}

	Pattern(Pattern const & other) :
	  TPatternState(other),
		needleSize(other.needleSize),
		k(other.k),
		bitMasks(other.bitMasks),
		largePattern(NULL)
	{
		if (other.largePattern)
		{
			largePattern = new _MyersLargePattern;
			(*largePattern) = *(other.largePattern);
		}
		
	}

	~Pattern()
	{
		delete largePattern;
	}

	Pattern &
	operator = (Pattern const & other)
	{
		needleSize = other.needleSize;
		k = other.k;
		bitMasks = other.bitMasks;
		largePattern = NULL;
		if (other.largePattern)
		{
			largePattern = new _MyersLargePattern;
			(*largePattern) = *(other.largePattern);
		}
		return *this;
	}
	
	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl, int _limit = -1):
		needleSize(0),
		k(- _limit),
		largePattern(NULL)
	{
		setHost(*this, ndl);
	}
};


template <typename TNeedle, typename THasState, typename TFindBeginPatternSpec>
class Pattern<TNeedle, Myers<AlignTextBanded, THasState, TFindBeginPatternSpec> > {
//____________________________________________________________________________
public:
#ifdef SEQAN_SSE2_INT128
	typedef SSE2_int128 TWord;
#else
	typedef unsigned long TWord;
#endif

	typedef typename Iterator<TNeedle, Standard>::Type TIter;
	enum { MACHINE_WORD_SIZE = sizeof(TWord) * 8 };

	unsigned needleSize;
	unsigned score;				// the current score
	unsigned blockCount;		// the number of blocks
	unsigned k;					// the maximal number of differences allowed
	unsigned lastBlock;			// the block containing the last active cell

	TWord VP0;					// VP[0] (saves one dereferentiation)
	TWord VN0;					// VN[0]
	int scoreBit;

	Holder<TNeedle>		data_host;

//	String<int> mat;

	String<TWord> VP;
	String<TWord> VN;
	String<TWord> bitMasks;		// encoding the alphabet as bit-masks
	TWord scoreMask;			// the mask with a bit set at the position of the last active cell
	TWord finalScoreMask;		// a mask with a bit set on the position of the last row
	
	TIter ndlIter;				// iterate through the pattern
	
//____________________________________________________________________________

	Pattern(int _limit = -1):
		k(- _limit)
	{}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl, int _limit = -1):
		k(- _limit)
	{
		setHost(*this, ndl);
	}
//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Metafunctions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
struct FindBeginPatternSpec <Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > > 
{
	typedef TFindBeginPatternSpec Type;
};
template <typename TNeedle, typename THasState, typename TFindBeginPatternSpec>
struct FindBeginPatternSpec <Pattern<TNeedle, Myers<FindPrefix, THasState, TFindBeginPatternSpec> > >
{// no find begin for FindPrefix
	typedef void Type;
};


template <typename TPattern>
struct PatternState
{
	
};

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
struct PatternState<Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > >
{
	typedef _PatternState<Myers<TSpec, True, TFindBeginPatternSpec> > Type;
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TNeedle2>
inline void _patternFirstInit(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern, 
							  TNeedle2 & needle)
{
SEQAN_CHECKPOINT

	typedef typename Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >::TWord TWord;
	typedef typename Value<TNeedle>::Type TValue;

	pattern.needleSize = length(needle);
	unsigned blockCount = (pattern.needleSize + pattern.MACHINE_WORD_SIZE - 1) / pattern.MACHINE_WORD_SIZE;

	if (blockCount > 1) 
	{
		if (pattern.largePattern == NULL)
			pattern.largePattern = new _MyersLargePattern();

		pattern.largePattern->blockCount = blockCount;
		pattern.largePattern->finalScoreMask = (TWord)1 << ((pattern.needleSize + pattern.MACHINE_WORD_SIZE - 1) % pattern.MACHINE_WORD_SIZE);
	} 
	else 
	{
		delete pattern.largePattern;
		pattern.largePattern = NULL;
	}

	clear(pattern.bitMasks);
	fill(pattern.bitMasks, (ValueSize<TValue>::VALUE + 1) * blockCount, 0, Exact());

	// encoding the letters as bit-vectors
    for (unsigned j = 0; j < pattern.needleSize; j++)
		pattern.bitMasks[
			blockCount * ordValue((typename Value<TNeedle>::Type) value(needle, j))
			+ j / pattern.MACHINE_WORD_SIZE
		] |= (TWord)1 << (j % pattern.MACHINE_WORD_SIZE);
		//pattern.bitMasks[pattern.blockCount * ordValue((CompareType< Value< TNeedle >::Type, Value< Container< THaystack >::Type >::Type >::Type) needle[j]) + j/pattern.MACHINE_WORD_SIZE] = pattern.bitMasks[pattern.blockCount * ordValue((CompareType< Value< TNeedle >::Type, Value< Container< THaystack >::Type >::Type >::Type) needle[j]) + j/MACHINE_WORD_SIZE] | ((TWord)1 << (j%MACHINE_WORD_SIZE));
		
	_findBeginInit(pattern);
}


template <typename TNeedle, typename THasState, typename TFindBeginPatternSpec, typename TNeedle2>
inline void _patternFirstInit(Pattern<TNeedle, Myers<AlignTextBanded, THasState, TFindBeginPatternSpec> > & me, 
							  TNeedle2 & ndl)
{
SEQAN_CHECKPOINT
	typedef typename Pattern<TNeedle, Myers<AlignTextBanded, THasState, TFindBeginPatternSpec> >::TWord TWord;
	me.needleSize = length(ndl);
	me.finalScoreMask = (TWord)1 << (me.MACHINE_WORD_SIZE - 1);

	_findBeginInit(me);
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline void _patternMatchNOfPatternImpl(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern, 
										bool match)
{
    SEQAN_CHECKPOINT;

	typedef typename Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >::TWord TWord;
	typedef typename Value<TNeedle>::Type TValue;
	unsigned blockCount = (pattern.largePattern == NULL)? 1: pattern.largePattern->blockCount;

	// letters are encoded as bit-vectors
	for (unsigned j = 0; j < pattern.needleSize; j++)
	{
		TWord bit = (TWord)1 << (j % pattern.MACHINE_WORD_SIZE);
		bool allNull = true;
		int idx = j / pattern.MACHINE_WORD_SIZE;

		for (int i = 0; i < 4; ++i, idx += blockCount)
			allNull &= (pattern.bitMasks[idx] & bit) == (TWord)0;

		if (allNull)
		{	// all bits are 0 => this letter must be 'N'
			if (match)
			{
				for (; idx >= 0; idx -= blockCount)
					pattern.bitMasks[idx] |= bit;
			} else
			{
				for (; idx >= 0; idx -= blockCount)
					pattern.bitMasks[idx] &= ~bit;
			}
		}
	}
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
void _patternMatchNOfPattern(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern, 
							 bool match)
{
    SEQAN_CHECKPOINT;
    _patternMatchNOfPatternImpl(pattern, match);
    _patternMatchNOfPatternImpl(pattern.data_findBeginPattern, match);
}


template <typename TNeedle, typename TSpec, typename THasState>
void _patternMatchNOfPattern(Pattern<TNeedle, Myers<TSpec, THasState, void> > & pattern, 
							 bool match)
{
    SEQAN_CHECKPOINT;
    _patternMatchNOfPatternImpl(pattern, match);
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline void _patternMatchNOfFinderImpl(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern, 
									   bool match)
{
	SEQAN_CHECKPOINT;

	typedef typename Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >::TWord TWord;
	unsigned blockCount = (pattern.largePattern == NULL)? 1: pattern.largePattern->blockCount;

	// letters are encoded as bit-vectors
	if (match)
	{
		for (unsigned j = 0; j < pattern.needleSize; j++)
			pattern.bitMasks[blockCount * 4 + j / pattern.MACHINE_WORD_SIZE] |= (TWord)1 << (j % pattern.MACHINE_WORD_SIZE);
	} else {
		for (unsigned j = 0; j < pattern.needleSize; j++)
			pattern.bitMasks[blockCount * 4 + j / pattern.MACHINE_WORD_SIZE] &= ~((TWord)1 << (j % pattern.MACHINE_WORD_SIZE));
	}
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
void _patternMatchNOfFinder(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & me, bool match)
{
    SEQAN_CHECKPOINT;
    _patternMatchNOfFinderImpl(me, match);
    _patternMatchNOfFinderImpl(me.data_findBeginPattern, match);
}


template <typename TNeedle, typename TSpec, typename THasState>
void _patternMatchNOfFinder(Pattern<TNeedle, Myers<TSpec, THasState, void> > & pattern, 
							bool match)
{
    SEQAN_CHECKPOINT;
    _patternMatchNOfFinderImpl(pattern, match);
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TNeedle2>
void _myersSetHost(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern,
				   TNeedle2 const & ndl,
				   True) 
{
	setValue(pattern.data_host, ndl);
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TNeedle2>
void _myersSetHost(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & ,
				   TNeedle2 const & ,
				   False) 
{
	
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TNeedle2>
void setHost(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern,
			 TNeedle2 & ndl)
{
	SEQAN_CHECKPOINT;
	typedef Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > TPattern;
	_myersSetHost(pattern, ndl, typename TPattern::HasHost());
	_patternFirstInit(pattern, ndl);
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TNeedle2>
void setHost(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern,
			 TNeedle2 const & ndl) 
{
	SEQAN_CHECKPOINT;
	typedef Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > TPattern;
	_myersSetHost(pattern, ndl, typename TPattern::HasHost());
	_patternFirstInit(pattern, ndl);
}


//____________________________________________________________________________
/*
template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline TNeedle
host(Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > & me)
{
SEQAN_CHECKPOINT

	typedef typename Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >::TWord TWord;
	typedef typename Value<TNeedle>::Type TValue;

	TNeedle temp;
	resize(temp, me.needleSize, Exact());

	TValue v = TValue();
	for (unsigned i = 0; i < length(me.bitMasks); i += me.blockCount)
	{
		for (unsigned j = 0; j < me.needleSize; j++)
			if ((me.bitMasks[i + j / me.MACHINE_WORD_SIZE] & (TWord)1 << (j % me.MACHINE_WORD_SIZE)) != (TWord)0)
				temp[j] = v;
		++v;
	}
}

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline TNeedle
host(Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >  const & me)
{
SEQAN_CHECKPOINT

	typedef typename Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >::TWord TWord;
	typedef typename Value<TNeedle>::Type TValue;

	TNeedle temp;
	resize(temp, me.needleSize, Exact());

	TValue v = TValue();
	for (unsigned i = 0; i < length(me.bitMasks); i += me.blockCount)
	{
		for (unsigned j = 0; j < me.needleSize; j++)
			if ((me.bitMasks[i + j / me.MACHINE_WORD_SIZE] & (TWord)1 << (j % me.MACHINE_WORD_SIZE)) != (TWord)0)
				temp[j] = v;
		++v;
	}
}
*/


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline typename Host<Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > >::Type & 
host(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern)
{
SEQAN_CHECKPOINT
	return value(pattern.data_host);
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline typename Host<Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > const>::Type & 
host(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > const & pattern)
{
SEQAN_CHECKPOINT
	return value(pattern.data_host);
}


//____________________________________________________________________________

///.Function.scoreLimit.param.pattern.type:Spec.Myers

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline int
scoreLimit(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > const & pattern)
{
SEQAN_CHECKPOINT
	return - (int) pattern.k;
}


//____________________________________________________________________________

///.Function.setScoreLimit.param.pattern.type:Spec.Myers

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TScoreValue>
inline void 
setScoreLimit(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern, 
			  TScoreValue _limit)
{
SEQAN_CHECKPOINT
	pattern.k = (- _limit);
}


//____________________________________________________________________________

///.Function.getScore.param.pattern.type:Spec.Myers

template<typename TSpec, typename TFindBeginPatternSpec>
int getScore(_PatternState<Myers<TSpec, True, TFindBeginPatternSpec> > const & state)
{
	return -(int)state.score;
}


template<typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
int getScore(Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > const & state)
{
	return -(int)state.score;
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TFinder>
void _patternInit(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern, 
				  _PatternState<Myers<TSpec, True, TFindBeginPatternSpec> > & state,
				  TFinder &)
{
SEQAN_CHECKPOINT
	typedef typename Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >::TWord TWord;

	if (pattern.largePattern == NULL)
	{
		state.score = pattern.needleSize;
		state.VP0 = ~(TWord)0;
		state.VN0 = 0;
	} 
	else 
	{
		// TODO: is that good here?
		if (state.largeState == NULL)
			state.largeState = new _MyersLargeState();
		
		_MyersLargePattern &largePattern = *pattern.largePattern;
		_MyersLargeState &largeState = *state.largeState;
		// local_k either stores the score limit (me.k) or the
		// needle size minus one.  It is used for the mask
		// computation and setting the initial score (the
		// minus one is there because of the Ukkonen trick).
		int local_k = _min(pattern.k, pattern.needleSize - 1);
		state.score = local_k + 1;
		largeState.scoreMask = (TWord)1 << (local_k % pattern.MACHINE_WORD_SIZE);
		largeState.lastBlock = local_k / pattern.MACHINE_WORD_SIZE; 
		if (largeState.lastBlock >= largePattern.blockCount)
			largeState.lastBlock = largePattern.blockCount - 1;

		clear(largeState.VP);
		fill(largeState.VP, largePattern.blockCount, ~(TWord)0, Exact());

		clear(largeState.VN);
		fill(largeState.VN, largePattern.blockCount, 0, Exact());
	}
}


template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TFinder>
void _patternInit(Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & pattern, 
				  TFinder & finder)
{
    SEQAN_CHECKPOINT;
    _patternInit(pattern, pattern, finder);
}


template <typename TNeedle, typename TFinder, typename THasState, typename TFindBeginPatternSpec>
void _patternInit(Pattern<TNeedle, Myers<AlignTextBanded, THasState, TFindBeginPatternSpec> > &me, TFinder &finder)
{
SEQAN_CHECKPOINT
	typedef typename Pattern<TNeedle, Myers<AlignTextBanded, THasState, TFindBeginPatternSpec> >::TWord TWord;
	typedef typename Value<TNeedle>::Type TValue;

	me.ndlIter = begin(host(me), Standard());
	unsigned diagWidth = length(container(finder)) - me.needleSize;
	if (diagWidth >= me.needleSize)
		diagWidth = me.needleSize - 1;
	me.blockCount = diagWidth / me.MACHINE_WORD_SIZE + 1;

	clear(me.bitMasks);
	fill(me.bitMasks, ValueSize<TValue>::VALUE * me.blockCount, 0, Exact());

	if (me.blockCount == 1)
	{
		me.score = 0;
		me.scoreMask = 1;
		me.scoreBit = 0;
		me.VP0 = ~0;
		me.VN0 = 0;
	} 
	else 
	{
/*		me.score = me.k + 1;
		me.scoreMask = (TWord)1 << (me.k % me.MACHINE_WORD_SIZE);
		me.lastBlock = me.k / me.MACHINE_WORD_SIZE; 
		if (me.lastBlock >= me.blockCount)
			me.lastBlock = me.blockCount - 1;
*/
		clear(me.VP);
		fill(me.VP, me.blockCount, 0, Exact());

		clear(me.VN);
		fill(me.VN, me.blockCount, 0, Exact());
	}
}


//////////////////////////////////////////////////////////////////////////////
// Myers-Ukkonen for semi-global edit-distance-alignments
// (version for needles longer than one machineword)
//////////////////////////////////////////////////////////////////////////////


template <typename TFinder, typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TSize>
inline bool _findMyersLargePatterns (TFinder & finder, 
									 Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern,
									 _PatternState<Myers<TSpec, True, TFindBeginPatternSpec> > & state,
									 TSize haystack_length) 
{
SEQAN_CHECKPOINT
	typedef typename Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >::TWord TWord;

	TWord X, D0, HN, HP, temp;
	TWord carryD0, carryHP, carryHN;
	unsigned shift, limit, currentBlock;
	_MyersLargePattern &largePattern = *pattern.largePattern;
	_MyersLargeState &largeState = *state.largeState;

	while (position(finder) < haystack_length) 
	{
		carryD0 = carryHN = 0;
		carryHP = (int)_MyersUkkonenHP0<TSpec>::VALUE; // FIXME: replace Noting with TSpec

		// if the active cell is the last of it's block, one additional block has to be calculated
		limit = largeState.lastBlock + (unsigned)(largeState.scoreMask >> (pattern.MACHINE_WORD_SIZE - 1));

		if (limit == largePattern.blockCount)
			limit--;

		shift = largePattern.blockCount * ordValue((typename Value< TNeedle >::Type) *finder);

		// computing the necessary blocks, carries between blocks following one another are stored
		for (currentBlock = 0; currentBlock <= limit; currentBlock++) 
		{
			X = pattern.bitMasks[shift + currentBlock] | largeState.VN[currentBlock];
	
			temp = largeState.VP[currentBlock] + (X & largeState.VP[currentBlock]) + carryD0;
			if (carryD0 != (TWord)0)
				carryD0 = temp <= largeState.VP[currentBlock];
			else
				carryD0 = temp < largeState.VP[currentBlock];
			
			D0 = (temp ^ largeState.VP[currentBlock]) | X;
			HN = largeState.VP[currentBlock] & D0;
			HP = largeState.VN[currentBlock] | ~(largeState.VP[currentBlock] | D0);
			
			X = (HP << 1) | carryHP;
			carryHP = HP >> (pattern.MACHINE_WORD_SIZE - 1);
			
			largeState.VN[currentBlock] = X & D0;

			temp = (HN << 1) | carryHN;
			carryHN = HN >> (pattern.MACHINE_WORD_SIZE - 1);
								
		 	largeState.VP[currentBlock] = temp | ~(X | D0);

			//if the current block is the one containing the last active cell the new score is computed
			if (currentBlock == largeState.lastBlock) {
				if ((HP & largeState.scoreMask) != (TWord)0)
					state.score++;
				else if ((HN & largeState.scoreMask) != (TWord)0)
					state.score--;
			}
		}

		// updating the last active cell
		while (state.score > pattern.k) {
			if ((largeState.VP[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
				state.score--;
			else if ((largeState.VN[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
				state.score++;

			largeState.scoreMask >>= 1;
			if (largeState.scoreMask == (TWord)0) 
			{
				largeState.lastBlock--;
				if (TYPECMP<TSpec, FindPrefix>::VALUE && largeState.lastBlock == (unsigned)-1)
					break;
				largeState.scoreMask = (TWord)1 << (pattern.MACHINE_WORD_SIZE - 1);
			}
		}

		if ((largeState.scoreMask == largePattern.finalScoreMask) && (largeState.lastBlock == largePattern.blockCount - 1))
		{
			_setFinderEnd(finder);
			if (TYPECMP<TSpec, FindPrefix>::VALUE)
			{
				_setFinderLength(finder, endPosition(finder));
			}
			return true;
		}
		else {
			largeState.scoreMask <<= 1;
			if (!largeState.scoreMask) {
				largeState.scoreMask = 1;
				largeState.lastBlock++;
			}
			
			if ((largeState.VP[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
				state.score++;
			else if ((largeState.VN[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
				state.score--;
		}

//		SEQAN_ASSERT (me.score >= 0);

		goNext(finder);
	}

	return false;
}


template <typename TFinder, typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TSize>
inline bool _findMyersSmallPatterns (TFinder & finder, 
									 Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern,
									 _PatternState<Myers<TSpec, True, TFindBeginPatternSpec> > & state,
									 TSize haystack_length) 
{
SEQAN_CHECKPOINT

	typedef typename Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >::TWord TWord;

	TWord X, D0, HN, HP;
	TWord lastBit = (TWord)1 << (pattern.needleSize - 1);

	// computing the blocks
	while (position(finder) < haystack_length) 
	{
		X = pattern.bitMasks[ordValue((typename Value<TNeedle>::Type) *finder)] | state.VN0;
		
		D0 = ((state.VP0 + (X & state.VP0)) ^ state.VP0) | X;
		HN = state.VP0 & D0;
		HP = state.VN0 | ~(state.VP0 | D0);
		X = (HP << 1) | (TWord)(int)_MyersUkkonenHP0<TSpec>::VALUE; // FIXME: replace Nothing by TSpec
		state.VN0 = X & D0;
		state.VP0 = (HN << 1) | ~(X | D0);

		if ((HP & lastBit) != (TWord)0)
			state.score++;
		else if ((HN & lastBit) != (TWord)0)
			state.score--;

		if (state.score <= pattern.k)
		{
			_setFinderEnd(finder);
			if (TYPECMP<TSpec, FindPrefix>::VALUE)
			{
				_setFinderLength(finder, endPosition(finder));
			}
			return true;
		}
		// 
		// if (TYPECMP<TSpec, FindPrefix>::VALUE)
		// {//limit haystack length during prefix search
		// 
		// }
		
		goNext(finder);
	}

	return false;
}


//////////////////////////////////////////////////////////////////////////////
// Myers-Ukkonen as a banded alignment
// the band includes the main diagonal and the diagonals above
// the band width is (blockCount * MACHINE_WORD_SIZE)
//////////////////////////////////////////////////////////////////////////////


template <typename TFinder, typename TNeedle, typename THasState, typename TFindBeginPatternSpec, typename TSize>
inline bool 
_findMyersLargePatterns(
	TFinder & finder, 
	Pattern<TNeedle, Myers<AlignTextBanded, THasState, TFindBeginPatternSpec> > & me,
	TSize) 
{
SEQAN_CHECKPOINT
	typedef typename Pattern<TNeedle, Myers<AlignTextBanded, THasState, TFindBeginPatternSpec> >::TWord TWord;
	typedef typename Value<TNeedle>::Type TValue;

	TWord X, D0, HN, HP, temp;
	TWord carryD0, carryHP, carryHN;
	unsigned shift, limit, currentBlock;

	while (!atEnd(finder)) {
		// shift bitmasks and states
		if (!atEnd(me.ndlIter)) 
		{
			TWord carryVN = 0;
			TWord carryVP = 1;
			for(int j = me.blockCount - 1; j >= 0; --j) 
			{
				TWord newCarryVN = me.VN[j] & 1;
				TWord newCarryVP = me.VP[j] & 1;
				me.VN[j] = (me.VN[j] >> 1) | (carryVN << (me.MACHINE_WORD_SIZE - 1));
				me.VP[j] = (me.VP[j] >> 1) | (carryVP << (me.MACHINE_WORD_SIZE - 1));
				carryVN = newCarryVN;
				carryVP = newCarryVP;
			}
			for(unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i) 
			{
				TWord carry = 0;
				for(int j = me.blockCount - 1; j >= 0; --j)
				{
					unsigned pos = i * ValueSize<TValue>::VALUE + j;
					TWord newCarry = me.bitMasks[pos] & 1;
					me.bitMasks[pos] = (me.bitMasks[pos] >> 1) | (carry << (me.MACHINE_WORD_SIZE - 1));
					carry = newCarry;
				}
			}

			me.bitMasks[me.blockCount * (ordValue(*me.ndlIter) + 1) - 1]
				|= (TWord)1 << (me.MACHINE_WORD_SIZE - 1);

			goNext(me.ndlIter);
		}

		carryD0 = carryHN = 0;
		carryHP = _MyersUkkonenHP0<AlignTextBanded>::VALUE;

		// if the active cell is the last of it's block, one additional block has to be calculated
		limit = me.lastBlock + (me.scoreMask >> (me.MACHINE_WORD_SIZE - 1));

		if (limit == me.blockCount)
			limit--;

		shift = me.blockCount * ordValue((typename Value< TNeedle >::Type) *finder);

		// computing the necessary blocks, carries between blocks following one another are stored
		for (currentBlock = 0; currentBlock <= limit; currentBlock++) {
			X = me.bitMasks[shift + currentBlock] | me.VN[currentBlock];
	
			temp = me.VP[currentBlock] + (X & me.VP[currentBlock]) + carryD0;
			if (carryD0 != (TWord)0)
				carryD0 = temp <= me.VP[currentBlock];
			else
				carryD0 = temp < me.VP[currentBlock];
			
			D0 = (temp ^ me.VP[currentBlock]) | X;
			HN = me.VP[currentBlock] & D0;
			HP = me.VN[currentBlock] | ~(me.VP[currentBlock] | D0);
			
			X = (HP << 1) | carryHP;
			carryHP = HP >> (me.MACHINE_WORD_SIZE - 1);
			
			me.VN[currentBlock] = X & D0;

			temp = (HN << 1) | carryHN;
			carryHN = HN >> (me.MACHINE_WORD_SIZE - 1);
								
		 	me.VP[currentBlock] = temp | ~(X | D0);

			//if the current block is the one containing the last active cell the new score is computed
			if (currentBlock == me.lastBlock) {
				if ((HP & me.scoreMask) != (TWord)0)
					me.score++;
				else if ((HN & me.scoreMask) != (TWord)0)
					me.score--;
			}
		}

		// updating the last active cell
		while (me.score > me.k) {
			if ((me.VP[me.lastBlock] & me.scoreMask) != (TWord)0)
				me.score--;
			else if ((me.VN[me.lastBlock] & me.scoreMask) != (TWord)0)
				me.score++;

			me.scoreMask >>= 1;
			if (me.scoreMask == (TWord)0) {
				me.scoreMask = (TWord)1 << (me.MACHINE_WORD_SIZE - 1);
				me.lastBlock--;
			}
		}

		if ((me.lastBlock == me.blockCount-1) && (me.scoreMask == me.finalScoreMask))
			return true;
		else {
			me.scoreMask <<= 1;
			if (me.scoreMask == (TWord)0) {
				me.scoreMask = 1;
				me.lastBlock++;
			}
			
			if ((me.VP[me.lastBlock] & me.scoreMask) != (TWord)0)
				me.score++;
			else if ((me.VN[me.lastBlock] & me.scoreMask) != (TWord)0)
				me.score--;
		}

//		SEQAN_ASSERT (me.score >= 0);

		goNext(finder);
	}

	return false;
}

// Probably need to put THasState in here too
/* 
template <typename TFinder, typename TNeedle>
inline bool 
_findMyersSmallPatterns(
	TFinder & finder, 
	Pattern<TNeedle, Myers<AlignTextBanded> > & me)
{
SEQAN_CHECKPOINT
	typedef typename Pattern<TNeedle, Myers<AlignTextBanded> >::TWord TWord;

	TWord X, D0, HN, HP;
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
	unsigned SHIFT=me.needleSize;
#endif

	if (!atEnd(me.ndlIter)) 
	{
		// Part 1: Go down the diagonal
		do
		{
			// shift bitmasks and states
			me.VN0 >>= 1;
			me.VP0 = (me.VP0 >> 1) | ((TWord)1 << (me.MACHINE_WORD_SIZE - 1)); // ignore the field left of the leftmost diagonal
			for(unsigned i = 0; i < length(me.bitMasks); ++i)
				me.bitMasks[i] >>= 1;
			me.bitMasks[ordValue(*me.ndlIter)] |= ((TWord)1 << (me.MACHINE_WORD_SIZE - 1));
			
			// Myers core
			X = me.bitMasks[ordValue((typename Value<TNeedle>::Type) *finder)] | me.VN0;
			D0 = ((me.VP0 + (X & me.VP0)) ^ me.VP0) | X;
			HN = me.VP0 & D0;
			HP = me.VN0 | ~(me.VP0 | D0);
			X = (HP << 1);
			me.VN0 = X & D0;
			me.VP0 = (HN << 1) | ~(X | D0);

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
			--SHIFT;
			for(int i=0; i<SHIFT; ++i) 
				::std::cerr << "   ";
			::std::cerr << "dD: ";
			for(int i=me.MACHINE_WORD_SIZE-1; i>=0 ;--i) 
			{
				CharString vd = " 1 ";
				if ((D0 & ((TWord)1 << i)) != (TWord)0) vd = " 0 ";
				::std::cerr << vd;
			}
			::std::cerr << "   ";
			::std::cerr << *finder;
#endif

			if ((D0 & ((TWord)1 << (me.MACHINE_WORD_SIZE - 1))) == (TWord)0)
				++me.score;

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
			::std::cerr << me.score <<::std::endl;
#endif

			goNext(me.ndlIter);
//			if (atEnd(me.ndlIter))
			{
				if (me.score <= me.k)
					return true;
//				break;
			}
			goNext(finder);
		} while (true);
		
		goNext(finder);
	}

	// Part 2: Go to the bottom-right of the parallelogram
	while (!atEnd(finder))
	{
		// Myers core
		X = me.bitMasks[ordValue((typename Value<TNeedle>::Type) *finder)] | me.VN0;
		D0 = ((me.VP0 + (X & me.VP0)) ^ me.VP0) | X;
		HN = me.VP0 & D0;
		HP = me.VN0 | ~(me.VP0 | D0);
		X = (HP << 1) |1;
		me.VN0 = X & D0;
		me.VP0 = (HN << 1) | ~(X | D0);

#ifdef __SEQAN_DEBUG_MYERSBITVECTOR
		::std::cerr << "dH: ";
		for(int i=me.MACHINE_WORD_SIZE-1; i>=0 ;--i) 
		{
			CharString hd = " 0 ";
			if ((HP & ((TWord)1 << i)) != (TWord)0) hd = " 1 ";
			if ((HN & ((TWord)1 << i)) != (TWord)0) hd = "-1 ";
			::std::cerr << hd;
		}
		::std::cerr << "   ";
		::std::cerr << *finder;
#endif

		if ((HP & ((TWord)1 << (me.MACHINE_WORD_SIZE - 1))) != (TWord)0)
			++me.score;
		else if ((HN & ((TWord)1 << (me.MACHINE_WORD_SIZE - 1))) != (TWord)0)
			--me.score;

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
		::std::cerr << me.score <<::std::endl;
#endif

		if (me.score <= me.k)
			return true;
		goNext(finder);
	}
	return false;
}
*/


template <typename TWord, typename TAlignSpec, typename THasState, typename TFindBeginPatternSpec>
inline int
_myersCoreSmall(TWord &VP, TWord &VN, TWord const &bitmap, int scoreBit, Myers<TAlignSpec, THasState, TFindBeginPatternSpec> ) 
{
	// Myers core
	TWord X = bitmap | VN;
	TWord D0 = ((VP + (X & VP)) ^ VP) | X;
	TWord HN = VP & D0;
	TWord HP = VN | ~(VP | D0);
	X = HP << 1;
	VN = X & D0;
	VP = (HN << 1) | ~(X | D0);
	return (int)((HP >> scoreBit) & 1) - (int)((HN >> scoreBit) & 1);
}

template <typename TWord, typename TAlignSpec, typename THasState, typename TFindBeginPatternSpec>
inline int
_myersCoreSmallDiag(TWord &VP, TWord &VN, TWord const &bitmap, int scoreBit, Myers<TAlignSpec, THasState, TFindBeginPatternSpec> ) 
{
	// Myers core
	TWord X = bitmap | VN;
	TWord D0 = ((VP + (X & VP)) ^ VP) | X;
	TWord HN = VP & D0;
	TWord HP = VN | ~(VP | D0);
	X = HP << 1;
	VN = X & D0;
	VP = (HN << 1) | ~(X | D0);
	return (!D0 >> scoreBit) & 1;
}

template <typename TFinder, typename TNeedle, typename THasState, typename TFindBeginPatternSpec, typename TSize>
inline bool 
_findMyersSmallPatterns(
	TFinder & finder, 
	Pattern<TNeedle, Myers<AlignTextBanded, THasState, TFindBeginPatternSpec> > & me,
	TSize)
{
SEQAN_CHECKPOINT
	typedef Pattern<TNeedle, Myers<AlignTextBanded, THasState, TFindBeginPatternSpec> > TPattern;
	typedef typename TPattern::TWord TWord;
	typedef typename TPattern::TIter TIter;

	TWord X, D0, HN, HP;
	//TWord maskMax = (TWord)1 << (length(container(finder)) - me.needleSize);
	int bitMax = length(container(finder)) - me.needleSize;

#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
	int SHIFT = me.needleSize - (me.ndlIter-begin(host(me), Standard()));
#endif

	TIter ndlEnd = end(host(me), Standard());

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
	if (me.ndlIter != ndlEnd && me.scoreBit != bitMax)
#else
	if (me.ndlIter != ndlEnd)
#endif
	{
		// Part 1: left-upper triangle of parallelogram
		do
		{
			me.bitMasks[ordValue(*me.ndlIter)] |= me.scoreMask;
			me.score += _myersCoreSmallDiag(
				me.VP0, me.VN0, me.bitMasks[ordValue((typename Value<TNeedle>::Type) *finder)],
				me.scoreBit, Myers<AlignTextBanded, TFindBeginPatternSpec>());

#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
			--SHIFT;
			::std::cerr << "1D:  ";
			for(int i=0; i<SHIFT; ++i) 
				::std::cerr << "     ";
			for(int i=me.scoreMask; i>0 ;i>>=1) 
			{
/*				char const *vd = "  1  ";
				if ((D0 & i) != (TWord)0) vd = "  0  ";
*/				char const *vd = "  0  ";
				if ((me.VP0 & i) != (TWord)0) vd = "  1  ";
				if ((me.VN0 & i) != (TWord)0) vd = " -1  ";
				::std::cerr << vd;
			}
			::std::cerr << "   ";
			::std::cerr << *finder;
			::std::cerr << me.score <<::std::endl;
#endif
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
			resize(me.VP,1);
			resize(me.VN,1);
			me.VP[0] = me.VP0;
			me.VN[0] = me.VN0;
#endif

			me.scoreMask <<= 1;
			++me.scoreBit;
			goNext(me.ndlIter);

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
			return true;
#endif
			if (me.scoreBit == bitMax || (me.ndlIter == ndlEnd))
			{
				if (me.score <= me.k)
					return true;
				break;
			}
			goNext(finder);
		} while (true);
		
		goNext(finder);
	}


	if (me.ndlIter != ndlEnd)
	{
		// Part 2: go down the parallelogram
		do
		{
			me.bitMasks[ordValue(*me.ndlIter)] |= me.scoreMask;
			me.score += _myersCoreSmallDiag(
				me.VP0, me.VN0, me.bitMasks[ordValue((typename Value<TNeedle>::Type) *finder)],
				me.scoreBit, Myers<AlignTextBanded, TFindBeginPatternSpec>());

#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
			--SHIFT;
			::std::cerr << "2D:  ";
			for(int i=0; i<SHIFT; ++i) 
				::std::cerr << "     ";
			for(int i=me.scoreMask; i>0 ;i>>=1) 
			{
/*				char const *vd = "  0  ";
				if ((me.VP0 & i) != (TWord)0) vd = "  1  ";
				if ((me.VN0 & i) != (TWord)0) vd = " -1  ";
*/				char const *vd = "  1  ";
				if ((D0 & i) != (TWord)0) vd = "  0  ";
				::std::cerr << vd;
			}
			::std::cerr << "   ";
			::std::cerr << *finder;
			::std::cerr << me.score <<::std::endl;
#endif
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
			resize(me.VP,1);
			resize(me.VN,1);
			me.VP[0] = me.VP0;
			me.VN[0] = me.VN0;
#endif

			// shift bitmasks and states
			me.VN0 >>= 1;
			me.VP0 = (me.VP0 >> 1) | ((TWord)1 << (me.MACHINE_WORD_SIZE - 1)); // ignore the field left of the leftmost diagonal
			for(unsigned i = 0; i < length(me.bitMasks); ++i)
				me.bitMasks[i] >>= 1;

			goNext(me.ndlIter);
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
			return true;
#endif
			if (atEnd(me.ndlIter))
			{
				if (me.score <= me.k)
					return true;
				break;
			}
			goNext(finder);
		} while (true);
		
		goNext(finder);
	}

	if (!atEnd(finder))
	{
		// Part 3: go down the parallelogram
		do
		{
			me.score += _myersCoreSmall(
				me.VP0, me.VN0, me.bitMasks[ordValue((typename Value<TNeedle>::Type) *finder)],
				--me.scoreBit, Myers<AlignTextBanded, TFindBeginPatternSpec>());

			// shift bitmasks and states
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
			::std::cerr << "3H:  ";
			for(int i=me.scoreMask; i>0 ;i>>=1) 
			{
				char const *vd = "  0  ";
/*				if ((HP & i) != (TWord)0) vd = "  1  ";
				if ((HN & i) != (TWord)0) vd = " -1  ";
*/				if ((me.VP0 & i) != (TWord)0) vd = "  1  ";
				if ((me.VN0 & i) != (TWord)0) vd = " -1  ";
				::std::cerr << vd;
			}
			::std::cerr << "   ";
			::std::cerr << *finder;
			::std::cerr << me.score <<::std::endl;
#endif
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
			resize(me.VP,1);
			resize(me.VN,1);
			me.VP[0] = me.VP0;
			me.VN[0] = me.VN0;
#endif
			me.VN0 >>= 1;
			me.VP0 = (me.VP0 >> 1) | ((TWord)1 << (me.MACHINE_WORD_SIZE - 1)); // ignore the field left of the leftmost diagonal
			for(unsigned i = 0; i < length(me.bitMasks); ++i)
				me.bitMasks[i] >>= 1;			

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
			return true;
#endif
			if (me.score <= me.k)
				return true;

			goNext(finder);
		} while (!atEnd(finder));
	}
    return false;
}


//////////////////////////////////////////////////////////////////////////////
// find

// First two for AlignTextBanded
template <typename TFinder, typename TNeedle, typename THasState, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder, 
				  Pattern<TNeedle, Myers<AlignTextBanded, THasState, TFindBeginPatternSpec> > & me)
{
SEQAN_CHECKPOINT
	typedef typename Haystack<TFinder>::Type THaystack;
	typedef typename Size<THaystack>::Type TSize;

	int k = scoreLimit(me);

	TSize prefix_begin_position; //for prefix search: the position where the prefix begins

	if (empty(finder))
	{
		// in seqan k is treated as score, here we need it as penalty, that is why it is negated
		me.k = -k;

		_patternInit(me, finder);
		_finderSetNonEmpty(finder);

		prefix_begin_position = position(finder);

		//TODO: adapt myers-ukkonnen to dynamically change k
	}
	else
	{
		if (atEnd(finder)) return false;
		goNext(finder);

		prefix_begin_position = beginPosition(finder);
	}

	//limit search width for prefix search
	TSize haystack_length = length(container(hostIterator(finder)));
	// since this find function is now only for AlignTextBanded it does not need to compare
	// if (TYPECMP<TSpec, FindPrefix>::VALUE)
	// {
	// 	TSize maxlen = prefix_begin_position + me.needleSize - k + 1;
	// 	if (haystack_length > maxlen)
	// 	{
	// 		haystack_length = maxlen;
	// 	}
	// }

	// distinguish between the version for needles not longer than one machineword and the version for longer needles
	if (me.large == NULL) 
		return _findMyersSmallPatterns(finder, me, haystack_length);
	else
		return _findMyersLargePatterns(finder, me, haystack_length);
}


template <typename TFinder, typename TNeedle, typename THasState, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder, 
				  Pattern<TNeedle, Myers<AlignTextBanded, THasState, TFindBeginPatternSpec> > & me, 
				  int const k)
{
	setScoreLimit(me, k);
	return find(finder, me);
}


template <typename TFinder, typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
				  Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern,
				  _PatternState<Myers<TSpec, True, TFindBeginPatternSpec> > & state)
{
SEQAN_CHECKPOINT
	typedef typename Haystack<TFinder>::Type THaystack;
	typedef typename Size<THaystack>::Type TSize;

	int k = scoreLimit(pattern);

	TSize prefix_begin_position; //for prefix search: the position where the prefix begins

	if (empty(finder))
	{
		// in seqan k is treated as score, here we need it as penalty, that is why it is negated
		pattern.k = -k;

		_patternInit(pattern, state, finder);
		_finderSetNonEmpty(finder);

		prefix_begin_position = position(finder);

		//TODO: adapt myers-ukkonnen to dynamically change k
	}
	else
	{
		if (atEnd(finder)) return false;
		goNext(finder);

		prefix_begin_position = beginPosition(finder);
	}

	TSize haystack_length = length(container(hostIterator(finder)));
	// limit search width for prefix search
	if (TYPECMP<TSpec, FindPrefix>::VALUE)
	{
		TSize maxlen = prefix_begin_position + pattern.needleSize - k + 1;
		if (haystack_length > maxlen)
		{
			haystack_length = maxlen;
		}
	}

	// distinguish between the version for needles not longer than one machineword and the version for longer needles
	if (pattern.largePattern == NULL) 
		return _findMyersSmallPatterns(finder, pattern, state, haystack_length);
	else
		return _findMyersLargePatterns(finder, pattern, state, haystack_length);
}


template <typename TFinder, typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
				  Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & pattern)
{
	return find(finder, pattern, pattern);
}


template <typename TFinder, typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
				  Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern,
				  _PatternState<Myers<TSpec, True, TFindBeginPatternSpec> > & state,
				  int const k)
{
	setScoreLimit(pattern, k);
	return find(finder, pattern, state);
}


template <typename TFinder, typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
				  Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & pattern,
				  int const k)
{
	return find(finder, pattern, pattern, k); //static cast
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
