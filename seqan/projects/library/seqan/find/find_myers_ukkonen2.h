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
  Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_MYERS_UKKONEN2_H
#define SEQAN_HEADER_FIND_MYERS_UKKONEN2_H

namespace SEQAN_NAMESPACE_MAIN 
{
 
// large basic pattern
struct _MyersLargeBasicPattern
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


// basic pattern
template <typename TNeedle>
class BasicPattern{
	
public:
#ifdef SEQAN_SSE2_INT128
	typedef SSE2_int128 TWord;
#else
	typedef unsigned long TWord;
#endif
	typedef typename Iterator<TNeedle, Standard>::Type TIter;
	enum { MACHINE_WORD_SIZE = sizeof(TWord) * 8 };

	unsigned needleSize;
	unsigned k;					// the maximal number of differences allowed
	String<TWord> bitMasks;		// encoding the alphabet as bit-masks

	Holder<TNeedle> data_host;

	_MyersLargeBasicPattern *largePattern;	// extra preprocessing info for large patterns

//____________________________________________________________________________

	BasicPattern(int _limit = -1):
		needleSize(0),
		k(- _limit),
		largePattern(NULL)
	{}

	BasicPattern(BasicPattern const & other):
		needleSize(other.needleSize),
		k(other.k),
		bitMasks(other.bitMasks),
		largePattern(NULL)
	{
		if (other.largePattern)
		{
			largePattern = new _MyersLargeBasicPattern;
			(*largePattern) = *(other.largePattern);
		}
	}

	~BasicPattern()
	{
		delete largePattern;
	}

	BasicPattern &
	operator = (BasicPattern const & other)
	{
		needleSize = other.needleSize;
		k = other.k;
		bitMasks = other.bitMasks;
		largePattern = NULL;
		if (other.largePattern)
		{
			largePattern = new _MyersLargeBasicPattern;
			(*largePattern) = *(other.largePattern);
		}
		return *this;
	}
	
	template <typename TNeedle2>
	BasicPattern(TNeedle2 const & ndl, int _limit = -1):
		needleSize(0),
		k(- _limit),
		largePattern(NULL)
	{
		setHost(*this, ndl);
	}
};


// state
class PatternState{
	
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

	PatternState():
		score(0),
		VP0(0),
		VN0(0),
		largeState(NULL)
	{}

	PatternState(PatternState const & other):
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

	PatternState &
	operator = (PatternState const & other)
	{
		score = other.score;
		VP0 = other.VP0;
		VN0 = other.VN0;
		if (other.largeState)
		{
			largeState = new _MyersLargeState;
			(*largeState) = *(other.largeState);
		}
		return *this;
	}
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////


template <typename TNeedle, typename TNeedle2>
inline void _patternFirstInit(BasicPattern<TNeedle> & basicPattern, 
							  TNeedle2 & needle)
{
SEQAN_CHECKPOINT

	typedef typename BasicPattern<TNeedle>::TWord TWord;
	typedef typename Value<TNeedle>::Type TValue;

	basicPattern.needleSize = length(needle);
	unsigned blockCount = (basicPattern.needleSize + basicPattern.MACHINE_WORD_SIZE - 1) / basicPattern.MACHINE_WORD_SIZE;

	if (blockCount > 1) 
	{
		if (basicPattern.largePattern == NULL)
			basicPattern.largePattern = new _MyersLargeBasicPattern();

		basicPattern.largePattern->blockCount = blockCount;
		basicPattern.largePattern->finalScoreMask = (TWord)1 << ((basicPattern.needleSize + basicPattern.MACHINE_WORD_SIZE - 1) % basicPattern.MACHINE_WORD_SIZE);
	} 
	else 
	{
		delete basicPattern.largePattern;
		basicPattern.largePattern = NULL;
	}

	clear(basicPattern.bitMasks);
	fill(basicPattern.bitMasks, (ValueSize<TValue>::VALUE + 1) * blockCount, 0, Exact());

	// encoding the letters as bit-vectors
    for (unsigned j = 0; j < basicPattern.needleSize; j++)
		basicPattern.bitMasks[
			blockCount * ordValue((typename Value<TNeedle>::Type) value(needle, j))
			+ j / basicPattern.MACHINE_WORD_SIZE
		] |= (TWord)1 << (j % basicPattern.MACHINE_WORD_SIZE);
		//basicPattern.bitMasks[basicPattern.blockCount * ordValue((CompareType< Value< TNeedle >::Type, Value< Container< THaystack >::Type >::Type >::Type) needle[j]) + j/basicPattern.MACHINE_WORD_SIZE] = basicPattern.bitMasks[basicPattern.blockCount * ordValue((CompareType< Value< TNeedle >::Type, Value< Container< THaystack >::Type >::Type >::Type) needle[j]) + j/MACHINE_WORD_SIZE] | ((TWord)1 << (j%MACHINE_WORD_SIZE));
		
}


template <typename TNeedle, typename TNeedle2>
void _myersSetHost(BasicPattern<TNeedle> & basicPattern, 
				   TNeedle2 const & ndl) 
{
	setValue(basicPattern.data_host, ndl);
}


template <typename TNeedle, typename TNeedle2>
void setHost(BasicPattern<TNeedle> & basicPattern, 
			 TNeedle2 & ndl)
{
	SEQAN_CHECKPOINT;
	// typedef Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > TPattern;
	_myersSetHost(basicPattern, ndl);
	_patternFirstInit(basicPattern, ndl);
}


template <typename TNeedle, typename TNeedle2>
void setHost(BasicPattern<TNeedle> & basicPattern, 
			 TNeedle2 const & ndl) 
{
	SEQAN_CHECKPOINT;
	// typedef Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > TPattern;
	_myersSetHost(basicPattern, ndl);
	_patternFirstInit(basicPattern, ndl);
}

template <typename TNeedle>
inline typename Host<BasicPattern<TNeedle> >::Type & 
(BasicPattern<TNeedle> & basicPattern)
{
SEQAN_CHECKPOINT
	return value(basicPattern.data_host);
}

template <typename TNeedle>
inline typename Host<BasicPattern<TNeedle> const>::Type & 
host(BasicPattern<TNeedle> const & basicPattern)
{
SEQAN_CHECKPOINT
	return value(basicPattern.data_host);
}



template <typename TNeedle>
inline int
scoreLimit(BasicPattern<TNeedle> & basicPattern)
{
SEQAN_CHECKPOINT
	return - (int) basicPattern.k;
}


template <typename TNeedle, typename TScoreValue>
inline void 
setScoreLimit(BasicPattern<TNeedle> & basicPattern, 
			  TScoreValue _limit)
{
SEQAN_CHECKPOINT
	basicPattern.k = (- _limit);
}

int getScore(PatternState const & state)
{
	return -(int)state.score;
}


template <typename TNeedle, typename TFinder>
void _patternInit(BasicPattern<TNeedle> & basicPattern, 
				  PatternState & state,
				  TFinder &)
{
SEQAN_CHECKPOINT
	typedef typename BasicPattern<TNeedle>::TWord TWord;

	if (basicPattern.largePattern == NULL)
	{
		state.score = basicPattern.needleSize;
		state.VP0 = ~(TWord)0;
		state.VN0 = 0;
	} 
	else 
	{
		// TODO: is that good here?
		if (state.largeState == NULL)
			state.largeState = new _MyersLargeState();
		
		_MyersLargeBasicPattern &largePattern = *basicPattern.largePattern;
		_MyersLargeState &largeState = *state.largeState;
		// local_k either stores the score limit (me.k) or the
		// needle size minus one.  It is used for the mask
		// computation and setting the initial score (the
		// minus one is there because of the Ukkonen trick).
		int local_k = _min(basicPattern.k, basicPattern.needleSize - 1);
		state.score = local_k + 1;
		int e = largePattern.blockCount;
		largeState.scoreMask = (TWord)1 << (local_k % basicPattern.MACHINE_WORD_SIZE);
		largeState.lastBlock = local_k / basicPattern.MACHINE_WORD_SIZE; 
		if (largeState.lastBlock >= largePattern.blockCount)
			largeState.lastBlock = largePattern.blockCount - 1;

		clear(largeState.VP);
		fill(largeState.VP, largePattern.blockCount, ~(TWord)0, Exact());

		clear(largeState.VN);
		fill(largeState.VN, largePattern.blockCount, 0, Exact());
	}
}


template <typename TFinder, typename TNeedle, typename TSize>
inline bool _findMyersLargePatterns (TFinder & finder, 
									 PatternState & state,
									 BasicPattern<TNeedle> & basicPattern,
									 TSize haystack_length) 
{
SEQAN_CHECKPOINT
	typedef typename BasicPattern<TNeedle>::TWord TWord;

	TWord X, D0, HN, HP, temp;
	TWord carryD0, carryHP, carryHN;
	unsigned shift, limit, currentBlock;
	_MyersLargeBasicPattern &largePattern = *basicPattern.largePattern;
	_MyersLargeState &largeState = *state.largeState;

	while (position(finder) < haystack_length) 
	{
		carryD0 = carryHN = 0;
		carryHP = (int)_MyersUkkonenHP0<Nothing>::VALUE; // FIXME: replace Noting with TSpec

		// if the active cell is the last of it's block, one additional block has to be calculated
		limit = largeState.lastBlock + (unsigned)(largeState.scoreMask >> (basicPattern.MACHINE_WORD_SIZE - 1));

		if (limit == largePattern.blockCount)
			limit--;

		shift = largePattern.blockCount * ordValue((typename Value< TNeedle >::Type) *finder);

		// computing the necessary blocks, carries between blocks following one another are stored
		for (currentBlock = 0; currentBlock <= limit; currentBlock++) 
		{
			X = basicPattern.bitMasks[shift + currentBlock] | largeState.VN[currentBlock];
	
			temp = largeState.VP[currentBlock] + (X & largeState.VP[currentBlock]) + carryD0;
			if (carryD0 != (TWord)0)
				carryD0 = temp <= largeState.VP[currentBlock];
			else
				carryD0 = temp < largeState.VP[currentBlock];
			
			D0 = (temp ^ largeState.VP[currentBlock]) | X;
			HN = largeState.VP[currentBlock] & D0;
			HP = largeState.VN[currentBlock] | ~(largeState.VP[currentBlock] | D0);
			
			X = (HP << 1) | carryHP;
			carryHP = HP >> (basicPattern.MACHINE_WORD_SIZE - 1);
			
			largeState.VN[currentBlock] = X & D0;

			temp = (HN << 1) | carryHN;
			carryHN = HN >> (basicPattern.MACHINE_WORD_SIZE - 1);
								
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
		while (state.score > basicPattern.k) {
			if ((largeState.VP[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
				state.score--;
			else if ((largeState.VN[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
				state.score++;

			largeState.scoreMask >>= 1;
			if (largeState.scoreMask == (TWord)0) 
			{
				largeState.lastBlock--;
				// if (TYPECMP<TSpec, FindPrefix>::VALUE && largeState.lastBlock == (unsigned)-1)
				// 	break;
				// largeState.scoreMask = (TWord)1 << (me.MACHINE_WORD_SIZE - 1);
			}
		}

		if ((largeState.scoreMask == largePattern.finalScoreMask) && (largeState.lastBlock == largePattern.blockCount - 1))
		{
			_setFinderEnd(finder);
			// if (TYPECMP<TSpec, FindPrefix>::VALUE)
			// {
			// 	_setFinderLength(finder, endPosition(finder));
			// }
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


template <typename TFinder, typename TNeedle, typename TSize>
inline bool _findMyersSmallPatterns (TFinder & finder, 
									 PatternState & state,
									 BasicPattern<TNeedle> & basicPattern,
									 TSize haystack_length) 
{
SEQAN_CHECKPOINT

	typedef typename BasicPattern<TNeedle>::TWord TWord;

	TWord X, D0, HN, HP;
	TWord lastBit = (TWord)1 << (basicPattern.needleSize - 1);

	// computing the blocks
	while (position(finder) < haystack_length) 
	{
		X = basicPattern.bitMasks[ordValue((typename Value<TNeedle>::Type) *finder)] | state.VN0;
		
		D0 = ((state.VP0 + (X & state.VP0)) ^ state.VP0) | X;
		HN = state.VP0 & D0;
		HP = state.VN0 | ~(state.VP0 | D0);
		X = (HP << 1) | (TWord)(int)_MyersUkkonenHP0<Nothing>::VALUE; // FIXME: replace Nothing by TSpec
		state.VN0 = X & D0;
		state.VP0 = (HN << 1) | ~(X | D0);

		if ((HP & lastBit) != (TWord)0)
			state.score++;
		else if ((HN & lastBit) != (TWord)0)
			state.score--;

		if (state.score <= basicPattern.k)
		{
			_setFinderEnd(finder);
			// if (TYPECMP<TSpec, FindPrefix>::VALUE)
			// {
			// 	_setFinderLength(finder, endPosition(finder));
			// }
			return true;
		}
/*
		if (TYPECMP<TSpec, FindPrefix>::VALUE)
		{//limit haystack length during prefix search

		}
*/		
		goNext(finder);
	}

	return false;
}


//////////////////////////////////////////////////////////////////////////////
// find
template <typename TFinder, typename TNeedle>
inline bool find (TFinder & finder,
				  PatternState & state,
				  BasicPattern<TNeedle> & basicPattern)
{
SEQAN_CHECKPOINT
	typedef typename Haystack<TFinder>::Type THaystack;
	typedef typename Size<THaystack>::Type TSize;

	int k = scoreLimit(basicPattern);

	TSize prefix_begin_position; //for prefix search: the position where the prefix begins

	if (empty(finder))
	{
		// in seqan k is treated as score, here we need it as penalty, that is why it is negated
		basicPattern.k = -k;

		_patternInit(basicPattern, state, finder);
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
	//limit search width for prefix search
	// if (TYPECMP<TSpec, FindPrefix>::VALUE)
	// {
	// 	TSize maxlen = prefix_begin_position + basicPattern.needleSize - k + 1;
	// 	if (haystack_length > maxlen)
	// 	{
	// 		haystack_length = maxlen;
	// 	}
	// }

	// distinguish between the version for needles not longer than one machineword and the version for longer needles
	if (basicPattern.largePattern == NULL) 
		return _findMyersSmallPatterns(finder, state, basicPattern, haystack_length);
	else
		return _findMyersLargePatterns(finder, state, basicPattern, haystack_length);
}


template <typename TFinder, typename TNeedle>
inline bool find (TFinder & finder,
				  PatternState & state,
				  BasicPattern<TNeedle> & basicPattern,
				  int const k)
{
	setScoreLimit(basicPattern, k);
	return find(finder, state, basicPattern);
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
