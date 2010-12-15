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
  Author: David Weese <david.weese@fu-berlin.de>
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_MYERS_UKKONEN_H
#define SEQAN_HEADER_FIND_MYERS_UKKONEN_H

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

struct NMatchesNone_;
struct NMatchesN_;
struct NMatchesAll_;

//FindInfix and FindPrefix are defined int find_base.h
template <typename TFinderCharSetPolicy = NMatchesN_, typename TPatternCharSetPolicy = NMatchesN_>
struct AlignTextBanded; // search query in a parallelogram

// TODO(holtgrew): Really deprecated?
//deprecated shortcuts:

/**
.Shortcut.MyersUkkonen:
..status:deprecated, use $Myers<FindInfix>$
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
..status:deprecated, use $Myers<FindPrefix>$
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
..status:deprecated, use $Myers<AlignTextBanded>$
..cat:Pattern Matching
..summary:Semin-global (query-global, text-local) pattern matching without findBegin() support.
..signature:MyersUkkonen
..shortcutfor:Spec.Myers
...signature:Myers<AlignedTextBanded, void>
..see:Spec.Myers
..see:Shortcut.MyersUkkonen
..see:Shortcut.MyersUkkonenGlobal
*/
typedef Myers<AlignTextBanded<NMatchesN_, NMatchesN_>, True, void> MyersUkkonenBanded;


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

template <typename TValue>
struct _MyersSmallAlphabet:
	public Eval<ValueSize<TValue>::VALUE <= 8> {};


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
// 	Holder<TNeedle>	data_host;	//defines the 
	typedef False HasHost;
};

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
struct _FindBegin< Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >, void>
{
	typedef False HasHost;
//need no findBegin if FindBeginPatternSpec is void
};


//////////////////////////////////////////////////////////////////////////////
// State Data
//////////////////////////////////////////////////////////////////////////////

// small state
template <typename TNeedle, typename TSpec>
struct _MyersSmallState
{
#ifdef SEQAN_SSE2_INT128
	typedef SSE2_int128 TWord;
#else
	typedef unsigned long TWord;
#endif

	TWord VP0;					// VP[0] (saves one dereferentiation)
	TWord VN0;					// VN[0]
    unsigned int errors;		// the current number of errors
    unsigned int maxErrors;		// the maximal number of errors allowed
};

template <typename TNeedle, typename TSmallAlphabet>
struct MyersSmallStateBandedShift_ {};
template <typename TNeedle>
struct MyersSmallStateBandedShift_<TNeedle, False> {
    typedef typename Value<TNeedle>::Type TValue;
    unsigned short shift[ValueSize<TValue>::VALUE];
};
template <typename TNeedle, typename TFinderCSP, typename TPatternCSP>
struct _MyersSmallState<TNeedle, AlignTextBanded<TFinderCSP, TPatternCSP> >:
    public MyersSmallStateBandedShift_<TNeedle, typename _MyersSmallAlphabet<typename Value<TNeedle>::Type>::Type>
{
#ifdef SEQAN_SSE2_INT128
	typedef SSE2_int128 TWord;
#else
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
	typedef unsigned char TWord;
#else
	typedef unsigned long TWord;
#endif
#endif
    typedef typename Value<TNeedle>::Type TValue;

    TWord bitMasks[ValueSize<TValue>::VALUE];
	TWord VP0;					// VP[0] (saves one dereferentiation)
	TWord VN0;					// VN[0]
    unsigned short errors;      // the current number of errors
    unsigned short maxErrors;   // the maximal number of errors allowed

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
    String<int> DPMat;
#endif
};

// large state
template <typename TNeedle, typename TSpec>
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
template <typename TNeedle, typename TFinderCSP, typename TPatternCSP>
struct _MyersLargeState<TNeedle, AlignTextBanded<TFinderCSP, TPatternCSP> >
{
#ifdef SEQAN_SSE2_INT128
	typedef SSE2_int128 TWord;
#else
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
	typedef unsigned char TWord;
#else
	typedef unsigned long TWord;
#endif
#endif
	unsigned lastBlock;			// the block containing the last active cell
	unsigned blockCount;		// the number of blocks
	String<TWord> VP;
	String<TWord> VN;
};

// TODO: should go elsewhere
// template <typename TNeedle, typename TSpec>
// class Pattern{};
// TODO: should go elsewhere
template <typename TNeedle, typename TSpec>
class _PatternState{};


template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
class _PatternState<TNeedle, Myers<TSpec, False, TFindBeginPatternSpec> > {};

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
class _PatternState<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> >:
    public _MyersSmallState<TNeedle, TSpec>
{	
public:
    typedef _MyersSmallState<TNeedle, TSpec>    TSmallState;
    typedef _MyersLargeState<TNeedle, TSpec>    TLargeState;
    typedef typename TSmallState::TWord         TWord;

	enum { MACHINE_WORD_SIZE = sizeof(TWord) * 8 };

	TLargeState *largeState;

//____________________________________________________________________________

	_PatternState():
		largeState(NULL) {}

	_PatternState(_PatternState const & other):
        TSmallState(other),
		largeState(NULL)
	{
		if (other.largeState)
			largeState = new TLargeState(*other.largeState);
	}

	~_PatternState()
	{
		delete largeState;
	}
	
	_PatternState &
	operator = (_PatternState const & other)
	{
        TSmallState::operator=(other);
		if (other.largeState)
		{
			if (largeState == NULL)
				largeState = new TLargeState;
			(*largeState) = *(other.largeState);
		} else {
			delete largeState;
            largeState = NULL;
        }
		return *this;
	}
};


//////////////////////////////////////////////////////////////////////////////
// Pattern Data
//////////////////////////////////////////////////////////////////////////////
 
template <typename TNeedle, typename TSpec>
struct _MyersSmallPattern
{
#ifdef SEQAN_SSE2_INT128
	typedef SSE2_int128 TWord;
#else
	typedef unsigned long TWord;
#endif

	String<TWord> bitMasks;		// encode the needle with bitmasks for each alphabet character
	unsigned needleSize;        // needle size

    _MyersSmallPattern():
        needleSize(0) {}
};
template <typename TNeedle, typename TFinderCSP, typename TPatternCSP>
struct _MyersSmallPattern<TNeedle, AlignTextBanded<TFinderCSP, TPatternCSP> >
{
#ifdef SEQAN_SSE2_INT128
	typedef SSE2_int128 TWord;
#else
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
	typedef unsigned char TWord;
#else
	typedef unsigned long TWord;
#endif
#endif

	Holder<TNeedle> data_host;  // needle holder (the banded version needs no preprocessed bitmasks)
};

// large basic pattern
template <typename TNeedle, typename TSpec>
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
template <typename TNeedle, typename TFinderCSP, typename TPatternCSP>
struct _MyersLargePattern<TNeedle, AlignTextBanded<TFinderCSP, TPatternCSP> > {};


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
class Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >:
    public _MyersSmallPattern<TNeedle, TSpec>,
    public _FindBegin<Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > >,
	public _PatternState<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >
{
	
public:
    typedef _MyersSmallPattern<TNeedle, TSpec>      TSmallPattern;
    typedef _MyersLargePattern<TNeedle, TSpec>      TLargePattern;
    typedef typename TSmallPattern::TWord           TWord;
	
	enum { MACHINE_WORD_SIZE = sizeof(TWord) * 8 };

	typedef _PatternState<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > TPatternState;

	TLargePattern *largePattern;	// extra preprocessing info for large patterns

//____________________________________________________________________________

	Pattern():
		largePattern(NULL) {}

	Pattern(int _limit):
		largePattern(NULL) 
	{
		setScoreLimit(*this, _limit);
	}

	Pattern(Pattern const & other) :
        TSmallPattern(other),
        TPatternState(other),
        largePattern(NULL)
	{
		if (other.largePattern)
			largePattern = new TLargePattern(*other.largePattern);
	}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl, int _limit = -1):
		largePattern(NULL)
	{
		setScoreLimit(*this, _limit);
		setHost(*this, ndl);
	}

	~Pattern()
	{
		delete largePattern;
	}

	Pattern &
	operator = (Pattern const & other)
	{
        TSmallPattern::operator=(other);
        TPatternState::operator=(other);
		if (other.largePattern)
        {
            if (largePattern == NULL)
                largePattern = new TLargePattern;
            (*largePattern) = *(other.largePattern);
        } else {
            delete largePattern;
            largePattern = NULL;
        }
		return *this;
	}
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
struct PatternState {};

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
struct PatternState<Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > >
{
	typedef _PatternState<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > Type;
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
			pattern.largePattern = new _MyersLargePattern<TNeedle, TSpec>();

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
			blockCount * ordValue(convert<typename Value<TNeedle>::Type>(getValue(needle, j)))
			+ j / pattern.MACHINE_WORD_SIZE
		] |= (TWord)1 << (j % pattern.MACHINE_WORD_SIZE);
		//pattern.bitMasks[pattern.blockCount * ordValue((CompareType< Value< TNeedle >::Type, Value< Container< THaystack >::Type >::Type >::Type) needle[j]) + j/pattern.MACHINE_WORD_SIZE] = pattern.bitMasks[pattern.blockCount * ordValue((CompareType< Value< TNeedle >::Type, Value< Container< THaystack >::Type >::Type >::Type) needle[j]) + j/MACHINE_WORD_SIZE] | ((TWord)1 << (j%MACHINE_WORD_SIZE));
		
	_findBeginInit(pattern, needle);
}


template <typename TNeedle, typename TFinderCSP, typename TPatternCSP, typename THasState, typename TFindBeginPatternSpec, typename TNeedle2>
inline void _patternFirstInit(Pattern<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, THasState, TFindBeginPatternSpec> > & pattern, 
							  TNeedle2 & ndl)
{
	_findBeginInit(pattern, ndl);
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
inline void 
_patternMatchNOfPattern(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern, bool match)
{
    SEQAN_CHECKPOINT;
    _patternMatchNOfPatternImpl(pattern, match);
    _patternMatchNOfPatternImpl(pattern.data_findBeginPattern, match);
}


template <typename TNeedle, typename TSpec, typename THasState>
inline void 
_patternMatchNOfPattern(Pattern<TNeedle, Myers<TSpec, THasState, void> > & pattern, bool match)
{
    SEQAN_CHECKPOINT;
    _patternMatchNOfPatternImpl(pattern, match);
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline void 
_patternMatchNOfFinderImpl(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern, bool match)
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
inline void 
_patternMatchNOfFinder(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern, bool match)
{
    SEQAN_CHECKPOINT;
    _patternMatchNOfFinderImpl(pattern, match);
    _patternMatchNOfFinderImpl(pattern.data_findBeginPattern, match);
}


template <typename TNeedle, typename TSpec, typename THasState>
inline void 
_patternMatchNOfFinder(Pattern<TNeedle, Myers<TSpec, THasState, void> > & pattern, bool match)
{
    SEQAN_CHECKPOINT;
    _patternMatchNOfFinderImpl(pattern, match);
}


// data_host is not used anymore, the needle can be reconstructed from the bitmasks
template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TNeedle2>
inline void 
_myersSetHost(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > &, TNeedle2 const &) 
{
}

template <typename TNeedle, typename TFinderCSP, typename TPatternCSP, typename THasState, typename TFindBeginPatternSpec, typename TNeedle2>
inline void 
_myersSetHost(Pattern<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, THasState, TFindBeginPatternSpec> > & pattern, TNeedle2 const & ndl) 
{
	setValue(pattern.data_host, ndl);
}

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TNeedle2>
inline void 
setHost(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern, TNeedle2 & ndl)
{
	SEQAN_CHECKPOINT;
	typedef Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > TPattern;
	_myersSetHost(pattern, ndl);
	_patternFirstInit(pattern, ndl);
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TNeedle2>
inline void 
setHost(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern, TNeedle2 const & ndl) 
{
	SEQAN_CHECKPOINT;
	typedef Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > TPattern;
	_myersSetHost(pattern, ndl);
	_patternFirstInit(pattern, ndl);
}


//____________________________________________________________________________


template <typename TNeedle, typename TFinderCSP, typename TPatternCSP, typename THasState, typename TFindBeginPatternSpec>
inline typename Host<Pattern<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, THasState, TFindBeginPatternSpec> > >::Type & 
host(Pattern<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, THasState, TFindBeginPatternSpec> > & pattern)
{
SEQAN_CHECKPOINT
	return value(pattern.data_host);
}


template <typename TNeedle, typename TFinderCSP, typename TPatternCSP, typename THasState, typename TFindBeginPatternSpec>
inline typename Host<Pattern<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, THasState, TFindBeginPatternSpec> > const>::Type & 
host(Pattern<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, THasState, TFindBeginPatternSpec> > const & pattern)
{
SEQAN_CHECKPOINT
	return value(pattern.data_host);
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline TNeedle
host(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > const & pattern)
{
SEQAN_CHECKPOINT

	typedef typename Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >::TWord TWord;
	typedef typename Value<TNeedle>::Type TValue;

	TNeedle temp;
	resize(temp, pattern.needleSize, Exact());

	unsigned blockCount = (pattern.needleSize + pattern.MACHINE_WORD_SIZE - 1) / pattern.MACHINE_WORD_SIZE;
	TValue v = TValue();
	for (unsigned i = 0; i < length(pattern.bitMasks); i += blockCount)
	{
		for (unsigned j = 0; j < pattern.needleSize; j++)
			if ((pattern.bitMasks[i + j / pattern.MACHINE_WORD_SIZE] & (TWord)1 << (j % pattern.MACHINE_WORD_SIZE)) != (TWord)0)
				temp[j] = v;
		++v;
	}
	return temp;
}

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline TNeedle
host(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern)
{
SEQAN_CHECKPOINT
	typedef Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > TPattern;
	return host(const_cast<TPattern const &>(pattern));
}

//____________________________________________________________________________

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline TNeedle
needle(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > const & pattern)
{
SEQAN_CHECKPOINT
	return host(pattern);
}

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline TNeedle
needle(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern)
{
SEQAN_CHECKPOINT
	typedef Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > TPattern;
	return host(const_cast<TPattern const &>(pattern));
}

//____________________________________________________________________________

///.Function.scoreLimit.param.pattern.type:Spec.Myers

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline int 
scoreLimit(_PatternState<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > const & state)
{
SEQAN_CHECKPOINT
	return - (int) state.maxErrors;
}

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline int
scoreLimit(Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > const & pattern)
{
SEQAN_CHECKPOINT
	return - (int) pattern.maxErrors;
}


//____________________________________________________________________________

///.Function.setScoreLimit.param.pattern.type:Spec.Myers

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TScoreValue>
inline void 
setScoreLimit(_PatternState<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & state,
           TScoreValue minScore)
{
SEQAN_CHECKPOINT
    // we need to convert the minimal score into a maximal penalty
    // that is why minScore is negated
    state.maxErrors = -minScore;
}
template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TScoreValue>
inline void 
setScoreLimit(Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & pattern, 
           TScoreValue minScore)
{
SEQAN_CHECKPOINT
    // we need to convert the minimal score into a maximal penalty
    // that is why minScore is negated
    pattern.maxErrors = -minScore;
}



//____________________________________________________________________________

///.Function.getScore.param.pattern.type:Spec.Myers

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline int 
getScore(_PatternState<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > const & state)
{
	return -(int)state.errors;
}
template<typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline int 
getScore(Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > const & state)
{
	return -(int)state.errors;
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TFinder>
inline bool 
_patternInit(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > const & pattern, 
		     _PatternState<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & state,
		     TFinder &)
{
SEQAN_CHECKPOINT
    typedef _MyersLargeState<TNeedle, TSpec> TLargeState;
    typedef typename TLargeState::TWord TWord;

	if (pattern.largePattern == NULL)
	{
		state.errors = pattern.needleSize;
		state.VP0 = ~(TWord)0;
		state.VN0 = 0;
        delete state.largeState;
        state.largeState = NULL;
	} 
	else 
	{
		if (state.largeState == NULL)
			state.largeState = new TLargeState;
		
		TLargeState &largeState = *state.largeState;
		// localMaxErrors either stores the maximal number of errors (me.maxErrors) or the needle size minus one.
		// It is used for the mask computation and setting the initial score (the minus one is there because of the Ukkonen trick).
		int localMaxErrors = _min(state.maxErrors, pattern.needleSize - 1);
		state.errors = localMaxErrors + 1;
		largeState.scoreMask = (TWord)1 << (localMaxErrors % pattern.MACHINE_WORD_SIZE);
		largeState.lastBlock = localMaxErrors / pattern.MACHINE_WORD_SIZE; 
		if (largeState.lastBlock >= pattern.largePattern->blockCount)
			largeState.lastBlock = pattern.largePattern->blockCount - 1;

		clear(largeState.VP);
		fill(largeState.VP, pattern.largePattern->blockCount, ~(TWord)0, Exact());

		clear(largeState.VN);
		fill(largeState.VN, pattern.largePattern->blockCount, 0, Exact());
	}
    return true;
}

//____________________________________________________________________________
// bitmask operations - small alphabet

template <typename TNeedle, typename TSpec>
finline void 
_myersPreInit(_PatternState<TNeedle, TSpec> &state, True) 
{
	typedef typename Value<TNeedle>::Type TValue;
    for (unsigned i = 0; i <= ValueSize<TValue>::VALUE; ++i)
        state.bitMasks[i] = 0;
}

template <typename TNeedle, typename TSpec>
finline void 
_myersPostInit(_PatternState<TNeedle, TSpec> &state, True) 
{
	typedef typename Value<TNeedle>::Type TValue;
	for (unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
		state.bitMasks[i] >>= 1;
}

template <typename TNeedle, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec, typename TValue, typename TShift>
finline void 
_myersAdjustBitmask(_PatternState<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > &state, TValue const value, TShift, True) 
{
	typedef typename _PatternState<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> >::TWord TWord;
    
    // compiler will optimize that
    if (TYPECMP<TPatternCSP, NMatchesAll_>::VALUE && value == unknownValue<TValue>())
    {
        for (unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
            state.bitMasks[i] = (state.bitMasks[i] >> 1) | ((TWord)1 << (BitsPerValue<TWord>::VALUE - 1));
    } 
    else
    {
        for (unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
            state.bitMasks[i] >>= 1;
        if (!(TYPECMP<TPatternCSP, NMatchesNone_>::VALUE && value == unknownValue<TValue>()))
            state.bitMasks[ordValue(value)] |= (TWord)1 << (BitsPerValue<TWord>::VALUE - 1);
    }

    if (TYPECMP<TFinderCSP, NMatchesAll_>::VALUE)
        state.bitMasks[ordValue(unknownValue<TValue>())] |= (TWord)1 << (BitsPerValue<TWord>::VALUE - 1);
    if (TYPECMP<TFinderCSP, NMatchesNone_>::VALUE)
        state.bitMasks[ordValue(unknownValue<TValue>())] &= ~((TWord)1 << (BitsPerValue<TWord>::VALUE - 1));
}

template <typename TNeedle, typename TSpec, typename TValue, typename TShift>
finline typename _PatternState<TNeedle, TSpec>::TWord
_myersGetBitmask(_PatternState<TNeedle, TSpec> &state, TValue const value, TShift, True) 
{
	return state.bitMasks[ordValue(value)];
}



//____________________________________________________________________________
// bitmask operations - large alphabet

template <typename TNeedle, typename TSpec>
finline void 
_myersPreInit(_PatternState<TNeedle, TSpec> &state, False) 
{
	typedef typename Value<TNeedle>::Type TValue;
    memset(state.bitMasks, 0, (ValueSize<TValue>::VALUE + 1) * sizeof(state.bitMasks[0]));
    memset(state.shift, 0, ValueSize<TValue>::VALUE * sizeof(state.shift[0]));
}

template <typename TNeedle, typename TSpec>
finline void 
_myersPostInit(_PatternState<TNeedle, TSpec> &, False) 
{
}

template <typename TNeedle, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec, typename TValue, typename TShift>
finline void 
_myersAdjustBitmask(_PatternState<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > &state, TValue const value, TShift const shift, False) 
{
	typedef typename _PatternState<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> >::TWord TWord;
    
    if (TYPECMP<TPatternCSP, NMatchesNone_>::VALUE && value == unknownValue<TValue>())
        return;
    
    register unsigned ord = ordValue(value);
	register unsigned short x = shift - state.shift[ord];
	if (x < BitsPerValue<TWord>::VALUE)
		state.bitMasks[ord] = (state.bitMasks[ord] >> x) | ((TWord)1 << (BitsPerValue<TWord>::VALUE - 1));
	else
		state.bitMasks[ord] = (TWord)1 << (BitsPerValue<TWord>::VALUE - 1);
	state.shift[ord] = shift;
}

template <typename TNeedle, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec, typename TValue, typename TShift>
finline typename _PatternState<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> >::TWord
_myersGetBitmask(_PatternState<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > &state, TValue const value, TShift const shift, False) 
{
	typedef typename _PatternState<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> >::TWord TWord;

    if (TYPECMP<TFinderCSP, NMatchesNone_>::VALUE && value == unknownValue<TValue>())
        return 0;

    if (TYPECMP<TFinderCSP, NMatchesAll_>::VALUE && value == unknownValue<TValue>())
        return (shift < BitsPerValue<TWord>::VALUE)? -1 << shift: -1;
    
    register unsigned ord = ordValue(value);
    register TWord res;
    register TShift x = shift - state.shift[ord];
	if (x < BitsPerValue<TWord>::VALUE) 
		res = state.bitMasks[ord] >> x;
	else
		res = 0;
    
    if (TYPECMP<TPatternCSP, NMatchesAll_>::VALUE)
    {
        ord = ordValue(unknownValue<TValue>());
        x = shift - state.shift[ord];
        if (x < BitsPerValue<TWord>::VALUE) 
            res |= state.bitMasks[ord] >> x;
    }
    return res;
}


template <typename TFinder, typename TNeedle, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec>
inline bool 
_patternInitSmallStateBanded(
    TFinder &finder,
	TNeedle const & needle, 
    _PatternState<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > & state)
{
	typedef Pattern<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > TPattern;
	typedef typename TPattern::TWord TWord;
	typedef typename Iterator<TNeedle, Standard>::Type TIter;
	typedef typename Value<TNeedle>::Type TValue;

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
    int col = 1;
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
    std::cerr << "     ";
    for (int i = length(needle); i != 0; --i)
        std::cerr << std::setw(5) << needle[i - 1];
    std::cerr << std::endl;
    std::cerr << "     ";
    for (int i = length(needle); i >= 0; --i)
        std::cerr << std::setw(5) << state.DPMat[i];
    std::cerr << std::endl;
#endif
#endif

    _myersPreInit(state, typename _MyersSmallAlphabet<TValue>::Type());

	TIter ndlIter = begin(needle, Standard());
	TIter ndlEnd = end(needle, Standard());

    register TWord VP = -1;
    register TWord VN = 0;
    register unsigned errors = 0;
    register unsigned cutOff = state.maxErrors + (length(container(finder)) - length(needle));

//    std::cerr<<std::hex<<"\t  "<<std::setw(17)<<' '<<"\tVN"<<std::setw(17)<<VN<<"\tVP"<<std::setw(17)<<VP<<std::dec<<std::endl;

    for (register unsigned short shift = 0; ndlIter != ndlEnd; ++ndlIter, goNext(finder), ++shift)
    {
        // PART 1: go down the parallelogram
        
        // adjust bitmask
        _myersAdjustBitmask(state, getValue(ndlIter), shift, typename _MyersSmallAlphabet<TValue>::Type());
        
        // diagonal Myers
        register TWord X = _myersGetBitmask(state, ordValue(*finder), shift, typename _MyersSmallAlphabet<TValue>::Type()) | VN;
        register TWord D0 = ((VP + (X & VP)) ^ VP) | X;
        register TWord HN = VP & D0;
        register TWord HP = VN | ~(VP | D0);
    //    const int PADDING = sizeof(TWord)*2 + 1;
    //    std::cerr << std::hex;
    //    std::cerr << "\tD0"<<std::setw(PADDING)<<(__uint64)D0<<"\tHN"<<std::setw(PADDING)<<(__uint64)HN<<"\tHP"<<std::setw(PADDING)<<(__uint64)HP << std::endl;
        X = D0 >> 1;
        VN = X & HP;
        VP = HN | ~(X | HP);
    //    std::cerr << "\t  "<<std::setw(PADDING)<<' '<<"\tVN"<<std::setw(PADDING)<<(__uint64)VN<<"\tVP"<<std::setw(PADDING)<<(__uint64)VP << std::endl;
    //    std::cerr << std::dec;
        errors += (~D0 >> (BitsPerValue<TWord>::VALUE - 1)) & 1;
        if (errors > cutOff) return false;

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
        std::cerr << "diag ";        
#endif
        int val = errors;
        state.DPMat[col*(length(needle)+1)+col] = val;
        for (int i = length(needle); i >=0; --i)
        {
            if (i > col)
            {
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
                std::cerr << "     ";
#endif
            } else
            {
                int shft = (int)BitsPerValue<TWord>::VALUE-1 - (col-i);
                if (shft >= 0)
                {
                    if (i < col)
                    {
                        TWord mask = (TWord)1 << (shft);
                        val -= ((VP & mask) != (TWord)0)? 1:0;
                        val += ((VN & mask) != (TWord)0)? 1:0;
                    }
                    state.DPMat[col*(length(needle)+1)+i] = val;
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
                    std::cerr << std::setw(5) << val;
                } else
                {
                    std::cerr << "     ";
#endif
                }
            }
        }
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
        std::cerr << std::setw(5) << *finder;
        std::cerr << std::setw(5) << errors << std::endl;
#endif
        ++col;
#endif
    }
    state.VP0 = VP;
    state.VN0 = VN;
    state.errors = errors;
    _myersPostInit(state, typename _MyersSmallAlphabet<TValue>::Type());
    return true;
}


template <typename TFinder, typename TNeedle, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec>
bool _stateInit(
	TFinder &finder,
	TNeedle const & needle, 
	_PatternState<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > & state)
{
SEQAN_CHECKPOINT
	typedef _PatternState<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > TState;
    typedef typename TState::TWord TWord;
    typedef typename TState::TLargeState TLargeState;
	typedef typename Value<TNeedle>::Type TValue;

	unsigned diagWidth = length(container(finder)) - length(needle);
//	if (diagWidth >= length(needle))
//		diagWidth = length() - 1;
	unsigned blockCount = diagWidth / state.MACHINE_WORD_SIZE + 1;
    
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
    clear(state.DPMat);
    fill(state.DPMat, (length(container(finder)) + 1) * (length(needle) + 1), -9);
    for (unsigned i = 0; i <= length(needle); ++i)
        state.DPMat[i] = i;
    for (unsigned i = 0; i <= length(container(finder)); ++i)
        state.DPMat[i * (length(needle) + 1)] = 0;
#endif

	if (blockCount <= 1)
    {
        delete state.largeState;
        state.largeState = NULL;
        return _patternInitSmallStateBanded(finder, needle, state);
	} 
	else 
	{
		// TODO: is that good here?
		if (state.largeState == NULL)
			state.largeState = new TLargeState;
		
		TLargeState &largeState = *state.largeState;
        largeState.blockCount = blockCount;

		clear(largeState.VP);
		fill(largeState.VP, blockCount, ~0, Exact());

		clear(largeState.VN);
		fill(largeState.VN, blockCount, 0, Exact());
        return true;
	}
}

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TFinder>
inline bool 
_patternInit(Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & pattern, TFinder & finder)
{
    SEQAN_CHECKPOINT;
    return _patternInit(pattern, pattern, finder);
}



//////////////////////////////////////////////////////////////////////////////
// Myers-Ukkonen for semi-global edit-distance-alignments
// (version for needles longer than one machineword)
//////////////////////////////////////////////////////////////////////////////


template <typename TFinder, typename TNeedle, typename TSpec, typename THasState, typename THasState2, typename TFindBeginPatternSpec, typename TSize>
inline bool _findMyersLargePatterns (TFinder & finder, 
									 Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern,
									 _PatternState<TNeedle, Myers<TSpec, THasState2, TFindBeginPatternSpec> > & state,
									 TSize haystack_length) 
{
SEQAN_CHECKPOINT
	typedef _MyersLargePattern<TNeedle, TSpec> TLargePattern;
	typedef _MyersLargeState<TNeedle, TSpec> TLargeState;
    typedef typename TLargeState::TWord TWord;
    
	TWord X, D0, HN, HP, temp;
	TWord carryD0, carryHP, carryHN;
	unsigned shift, limit, currentBlock;
    
	TLargePattern &largePattern = *pattern.largePattern;
	TLargeState &largeState = *state.largeState;

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

			// if the current block is the one containing the last active cell
			// the new number of errors is computed
			if (currentBlock == largeState.lastBlock) {
				if ((HP & largeState.scoreMask) != (TWord)0)
					state.errors++;
				else if ((HN & largeState.scoreMask) != (TWord)0)
					state.errors--;
			}
		}

		// updating the last active cell
		while (state.errors > state.maxErrors) {
			if ((largeState.VP[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
				state.errors--;
			else if ((largeState.VN[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
				state.errors++;

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
				state.errors++;
			else if ((largeState.VN[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
				state.errors--;
		}

//		SEQAN_ASSERT (state.errors >= 0);

		goNext(finder);
	}

	return false;
}


template <typename TFinder, typename TNeedle, typename TSpec, typename THasState, typename THasState2, typename TFindBeginPatternSpec, typename TSize>
inline bool 
_findMyersSmallPatterns(
	TFinder & finder, 
	Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern,
	_PatternState<TNeedle, Myers<TSpec, THasState2, TFindBeginPatternSpec> > & state,
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
			state.errors++;
		else if ((HN & lastBit) != (TWord)0)
			state.errors--;

		if (state.errors <= state.maxErrors)
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

template <typename TFinder, typename TNeedle, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec>
finline bool 
_findMyersSmallPatternsBanded(
	TFinder & finder, 
	TNeedle const & needle,
    _PatternState<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > & state)
{
SEQAN_CHECKPOINT
	typedef _PatternState<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > TState;
	typedef typename TState::TWord TWord;
	typedef typename Iterator<TNeedle, Standard>::Type TIter;
	typedef typename Value<TNeedle>::Type TValue;

    register TWord VP = state.VP0;
    register TWord VN = state.VN0;
    register TWord errors = state.errors;
    register TWord const maxErrors = state.maxErrors;

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
    unsigned col = position(finder) + 1;
#endif

	for (; !atEnd(finder); goNext(finder))
    {
		// PART 2: go right

		// normal Myers
		const unsigned short shift = length(needle);
		register TWord X = _myersGetBitmask(state, ordValue(*finder), shift, typename _MyersSmallAlphabet<TValue>::Type()) | VN;
		register TWord D0 = ((VP + (X & VP)) ^ VP) | X;
		register TWord HN = VP & D0;
		register TWord HP = VN | ~(VP | D0);
	//    const int PADDING = sizeof(TWord)*2 + 1;
	//    std::cerr << std::hex;
	//    std::cerr << "\tD0"<<std::setw(PADDING)<<(__uint64)D0<<"\tHN"<<std::setw(PADDING)<<(__uint64)HN<<"\tHP"<<std::setw(PADDING)<<(__uint64)HP<<std::endl;
		X = (HP << 1) | 1;
		VN = X & D0;
		VP = (HN << 1) | ~(X | D0);
	//    std::cerr << "\t  "<<std::setw(PADDING)<<' '<<"\tVN"<<std::setw(PADDING)<<(__uint64)VN<<"\tVP"<<std::setw(PADDING)<<(__uint64)VP<<std::endl;
	//    std::cerr << std::dec;
		errors += (HP >> (BitsPerValue<TWord>::VALUE - 2)) & 1;
        errors -= (HN >> (BitsPerValue<TWord>::VALUE - 2)) & 1;

        // shift bitmasks and states
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
        std::cerr << "horiz";
#endif
        int val = errors;
        state.DPMat[col*(length(needle)+1)+length(needle)] = val;
        for (int i = length(needle); i >= 0; --i)
        {
            int shft = (int)BitsPerValue<TWord>::VALUE-1 - (length(needle)-i);
            if (shft >= 0)
            {
                if (i < (int)length(needle))
                {
                    TWord mask = (TWord)1 << (shft);
                    val -= ((VP & mask) != (TWord)0)? 1:0;
                    val += ((VN & mask) != (TWord)0)? 1:0;
                }
                state.DPMat[col*(length(needle)+1)+i] = val;
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
                std::cerr << std::setw(5) << val;
            } else {
                std::cerr << "     ";
#endif
            }
        }
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
        std::cerr << std::setw(5) << *finder;
        std::cerr << std::setw(5) << errors << std::endl;
#endif
        ++col;
#endif

        if (errors <= maxErrors)
        {
            state.VP0 = VP;
            state.VN0 = VN;
            state.errors = errors;
            return true;
        }
    }
    return false;
}


//////////////////////////////////////////////////////////////////////////////
// find

template <typename TFinder, typename TNeedle, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder, 
				  TNeedle & needle,
                  _PatternState<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > & state)
{
	typedef typename Haystack<TFinder>::Type THaystack;
	typedef typename Size<THaystack>::Type TSize;

	if (empty(finder))
	{
		_finderSetNonEmpty(finder);
		if (!_stateInit(finder, needle, state))
            goEnd(finder);
            
        if (state.errors <= state.maxErrors)
        {
            goPrevious(finder);
            return true;
        }
		//TODO: adapt myers-ukkonnen to dynamically change maxErrors
	}
	else
	{
		if (atEnd(finder)) return false;
		goNext(finder);
	}

	// distinguish between the version for needles not longer than one machineword and the version for longer needles
	if (state.largeState == NULL) 
		return _findMyersSmallPatternsBanded(finder, needle, state);
//	else
//		return _findMyersLargePatterns(finder, needle, state);
	return false;
}

// First two for AlignTextBanded
template <typename TFinder, typename TNeedle, typename TFinderCSP, typename TPatternCSP, typename THasState, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder, 
				  Pattern<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, THasState, TFindBeginPatternSpec> > & pattern,
                  _PatternState<TNeedle, Myers<AlignTextBanded<TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > & state)
{
SEQAN_CHECKPOINT
	return find(finder, host(pattern), state);
}

template <typename TFinder, typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
				  Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern,
				  _PatternState<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & state)
{
SEQAN_CHECKPOINT
	typedef typename Haystack<TFinder>::Type THaystack;
	typedef typename Size<THaystack>::Type TSize;

	TSize prefix_begin_position; //for prefix search: the position where the prefix begins

	if (empty(finder))
	{
		_patternInit(pattern, state, finder);
		_finderSetNonEmpty(finder);

		prefix_begin_position = position(finder);

		//TODO: adapt myers-ukkonnen to dynamically change maxErrors
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
		TSize maxlen = prefix_begin_position + pattern.needleSize - scoreLimit(state) + 1;
		if (haystack_length > maxlen)
			haystack_length = maxlen;
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
				  _PatternState<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & state,
				  int const minScore)
{
	setScoreLimit(state, minScore);
	return find(finder, pattern, state);
}

template <typename TFinder, typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
				  TNeedle & needle,
				  _PatternState<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & state,
				  int const minScore)
{
	setScoreLimit(state, minScore);
	return find(finder, needle, state);
}

template <typename TFinder, typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
				  Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & pattern,
				  int const minScore)
{
	return find(finder, pattern, pattern, minScore); //static cast
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
