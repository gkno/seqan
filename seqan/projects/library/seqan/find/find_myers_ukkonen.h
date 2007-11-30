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

#ifndef SEQAN_HEADER_FIND_MYERS_UKKONEN_H
#define SEQAN_HEADER_FIND_MYERS_UKKONEN_H

namespace SEQAN_NAMESPACE_MAIN 
{
 
//////////////////////////////////////////////////////////////////////////////
// MyersUkkonen
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.MyersUkkonen:
..cat:Pattern Matching
..general:Class.Pattern
..summary:Provides fast approximate searching of one string in another using Myer's fast bit-parallel algorithm with application of the Ukkonen-trick.
..signature:Pattern<TNeedle, MyersUkkonen>
..param.TNeedle:The needle type.
...type:Class.String
..remarks.text:The needle-length must be smaller than the highest number that can be stored in an unsigned int.
*/

///.Class.Pattern.param.TSpec.type:Spec.MyersUkkonen

struct AlignTextLocal;			// search semi-global (query-global, text-local)
struct AlignTextGlobal;			// search global (query-global, text-global)
struct AlignTextBanded;			// search query in a parallelogram

template <typename TSpec>
struct _MyersUkkonen;

typedef Tag<_MyersUkkonen<AlignTextLocal> >			MyersUkkonen;
typedef Tag<_MyersUkkonen<AlignTextGlobal> >		MyersUkkonenGlobal;
typedef Tag<_MyersUkkonen<AlignTextBanded> >		MyersUkkonenBanded;


//____________________________________________________________________________
// bit 0 of the HP bit-vector
// 0 for begin-gap-free haystack
// 1 for global alignments of haystack

template <typename T>
struct _MyersUkkonenHP0 {
	enum { VALUE = 0 };
};

template <>
struct _MyersUkkonenHP0<AlignTextGlobal> {
	enum { VALUE = 1 };
};


//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
class Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > > {
//____________________________________________________________________________
public:
	typedef unsigned int TWord;
	enum { MACHINE_WORD_SIZE = sizeof(TWord) * 8 };

	unsigned needleSize;
	unsigned score;				// the current score
	unsigned blockCount;		// the number of blocks
	unsigned k;					// the maximal number of differences allowed
	unsigned lastBlock;			// the block containing the last active cell

	Holder<TNeedle>		data_needle;

	TWord VP0;					// VP[0] (saves one dereferentiation)
	TWord VN0;					// VN[0]

	String<TWord> VP;
	String<TWord> VN;
	String<TWord> bitMasks;		// encoding the alphabet as bit-masks
	TWord scoreMask;			// the mask with a bit set at the position of the last active cell
	TWord finalScoreMask;		// a mask with a bit set on the position of the last row
	
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


template <typename TNeedle>
class Pattern<TNeedle, Tag<_MyersUkkonen<AlignTextBanded> > > {
//____________________________________________________________________________
public:
	typedef unsigned int TWord;
	typedef typename Iterator<TNeedle, Rooted>::Type TIter;
	enum { MACHINE_WORD_SIZE = sizeof(TWord) * 8 };

	unsigned needleSize;
	unsigned score;				// the current score
	unsigned blockCount;		// the number of blocks
	unsigned k;					// the maximal number of differences allowed
	unsigned lastBlock;			// the block containing the last active cell

	TWord VP0;					// VP[0] (saves one dereferentiation)
	TWord VN0;					// VN[0]

	Holder<TNeedle>		data_needle;

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
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec, typename TNeedle2>
void _initMyersPattern(Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > > & me, 
					   TNeedle2 & needle) 
{
SEQAN_CHECKPOINT

	typedef typename Value<TNeedle>::Type TValue;

	me.needleSize = length(needle);
	me.blockCount = (me.needleSize + me.MACHINE_WORD_SIZE - 1) / me.MACHINE_WORD_SIZE;
	me.finalScoreMask = 1 << ((me.needleSize + me.MACHINE_WORD_SIZE - 1) % me.MACHINE_WORD_SIZE);

	clear(me.bitMasks);
	fill(me.bitMasks, ValueSize<TValue>::VALUE * me.blockCount, 0, Exact());

	// encoding the letters as bit-vectors
    for (unsigned j = 0; j < me.needleSize; j++)
		me.bitMasks[
			me.blockCount * ordValue((typename Value<TNeedle>::Type) value(needle, j))
			+ j / me.MACHINE_WORD_SIZE
		] |= 1 << (j % me.MACHINE_WORD_SIZE);
		//me.bitMasks[me.blockCount * ordValue((CompareType< Value< TNeedle >::Type, Value< Container< THaystack >::Type >::Type >::Type) needle[j]) + j/me.MACHINE_WORD_SIZE] = me.bitMasks[me.blockCount * ordValue((CompareType< Value< TNeedle >::Type, Value< Container< THaystack >::Type >::Type >::Type) needle[j]) + j/MACHINE_WORD_SIZE] | 1 << (j%MACHINE_WORD_SIZE);
}

template <typename TNeedle, typename TSpec, typename TNeedle2>
void setHost(Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > > & me, 
			 TNeedle2 & ndl)
{
SEQAN_CHECKPOINT
	me.data_needle = ndl;
	_initMyersPattern(me, ndl);
}
template <typename TNeedle, typename TSpec, typename TNeedle2>
void setHost(Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > > & me, 
			 TNeedle2 const & ndl) 
{
	me.data_needle = ndl;
	_initMyersPattern(me, ndl);
}

template <typename TNeedle, typename TNeedle2>
void setHost(Pattern<TNeedle, Tag<_MyersUkkonen<AlignTextBanded> > > & me, 
			 TNeedle2 & ndl) 
{
SEQAN_CHECKPOINT
	me.data_needle = ndl;

	me.needleSize = length(ndl);
	me.ndlIter = begin(ndl, Rooted());
	me.finalScoreMask = 1 << (me.MACHINE_WORD_SIZE - 1);
}

//____________________________________________________________________________

template <typename TNeedle, typename TSpec>
inline typename Host<Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > > >::Type & 
host(Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > > & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

template <typename TNeedle, typename TSpec>
inline typename Host<Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > > const>::Type & 
host(Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > >  const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

//____________________________________________________________________________

///.Function.scoreLimit.param.pattern.type:Spec.MyersUkkonen

template <typename TNeedle, typename TSpec>
inline int 
scoreLimit(Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > > const & me)
{
SEQAN_CHECKPOINT
	return - (int) me.k;
}

//____________________________________________________________________________

///.Function.scoreLimit.param.pattern.type:Spec.MyersUkkonen

template <typename TNeedle, typename TSpec, typename TScoreValue>
inline void 
setScoreLimit(Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > > & me, 
			  TScoreValue _limit)
{
SEQAN_CHECKPOINT
	me.k = (- _limit);
}

//____________________________________________________________________________

///.Function.getScore.param.pattern.type:Spec.MyersUkkonen

template <typename TNeedle, typename TSpec>
int getScore(Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > > & me) 
{
	return -(int)me.score;
}
//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////


//____________________________________________________________________________


template <typename TNeedle, typename TSpec, typename TFinder>
void _patternInit(Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > > &me, TFinder &)
{
SEQAN_CHECKPOINT
	typedef typename Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > >::TWord TWord;
	if (me.blockCount == 1) 
	{
		me.score = me.needleSize;
		me.VP0 = ~0;
		me.VN0 = 0;
	} 
	else 
	{
		me.score = me.k + 1;
		me.scoreMask = 1 << (me.k % me.MACHINE_WORD_SIZE);
		me.lastBlock = me.k / me.MACHINE_WORD_SIZE; 
		if (me.lastBlock >= me.blockCount)
			me.lastBlock = me.blockCount - 1;

		clear(me.VP);
		fill(me.VP, me.blockCount, (TWord) ~0, Exact());

		clear(me.VN);
		fill(me.VN, me.blockCount, 0, Exact());
	}
}


template <typename TNeedle, typename TFinder>
void _patternInit(Pattern<TNeedle, Tag<_MyersUkkonen<AlignTextBanded> > > &me, TFinder &finder)
{
SEQAN_CHECKPOINT
	typedef typename Pattern<TNeedle, Tag<_MyersUkkonen<AlignTextBanded> > >::TWord TWord;
	typedef typename Value<TNeedle>::Type TValue;
/*
	fill(me.mat, (length(container(finder))+1) * (me.needleSize+1), 0, Exact());
	for(int i=1;i<=me.needleSize;++i)
		me.mat[i] = i;
	
	for(int j=1;j<=length(container(finder));++j)
		for(int i=1;i<=me.needleSize;++i) {
			int min = me.mat[(j-1)*(me.needleSize+1)+i-1];
			if (container(finder)[j-1] != container(me.ndlIter)[i-1]) ++min;
			int h = me.mat[(j-1)*(me.needleSize+1)+i]+1;
			int v = me.mat[(j)*(me.needleSize+1)+i-1]+1;
			int d = j-i;
			if (0<=d && d<8)
			{
				if (0<d)
					if (min>h) min=h;
				if (d<7)
					if (min>v) min=v;
				me.mat[j*(me.needleSize+1)+i] = min;
			}
		}

	for(int j=1;j<=length(container(finder));++j) {
		::std::cout << "VD: ";
		for(int i=me.needleSize;i>=1;--i) {
			int diff = me.mat[j*(me.needleSize+1)+i] - me.mat[(j)*(me.needleSize+1)+i-1];
			int d = j-i;
			if (0<=d && d<8) {
				if (diff>=0) ::std::cout<<" ";
				::std::cout << diff<<" ";
			} else
				::std::cout << "   ";
		}
		::std::cout << "   ";
		if (j>0)
			::std::cout << container(finder)[j-1];
		::std::cout << ::std::endl;
	}

	for(int j=1;j<=length(container(finder));++j) {
		::std::cout << "HD: ";
		for(int i=me.needleSize;i>=1;--i) {
			int diff = me.mat[j*(me.needleSize+1)+i];// - me.mat[(j-1)*(me.needleSize+1)+i];
			int d = j-i;
			if (0<=d && d<8) {
				if (diff>=0) ::std::cout<<" ";
				::std::cout << diff<<" ";
			} else
				::std::cout << "   ";
		}
		::std::cout << "   ";
		if (j>0)
			::std::cout << container(finder)[j-1];
		::std::cout << ::std::endl;
	}

	for(int j=1;j<=length(container(finder));++j) {
		::std::cout << "DD: ";
		for(int i=me.needleSize;i>=1;--i) {
			int diff = me.mat[j*(me.needleSize+1)+i] - me.mat[(j-1)*(me.needleSize+1)+i-1];
			int d = j-i;
			if (0<=d && d<8) {
				if (diff>=0) ::std::cout<<" ";
				::std::cout << diff<<" ";
			} else
				::std::cout << "   ";
		}
		::std::cout << "   ";
		if (j>0)
			::std::cout << container(finder)[j-1];
		::std::cout << ::std::endl;
	}

*/

	goBegin(me.ndlIter);
	unsigned diagWidth = length(container(finder)) - me.needleSize;
	if (diagWidth >= me.needleSize)
		diagWidth = me.needleSize - 1;
	me.blockCount = diagWidth / me.MACHINE_WORD_SIZE + 1;

	clear(me.bitMasks);
	fill(me.bitMasks, ValueSize<TValue>::VALUE * me.blockCount, 0, Exact());

	if (me.blockCount == 1) 
	{
		me.score = 0;
		me.VP0 = ~0;
		me.VN0 = 0;
	} 
	else 
	{
/*		me.score = me.k + 1;
		me.scoreMask = 1 << (me.k % me.MACHINE_WORD_SIZE);
		me.lastBlock = me.k / me.MACHINE_WORD_SIZE; 
		if (me.lastBlock >= me.blockCount)
			me.lastBlock = me.blockCount - 1;
*/
		clear(me.VP);
		fill(me.VP, me.blockCount, (TWord) ~0, Exact());

		clear(me.VN);
		fill(me.VN, me.blockCount, 0, Exact());
	}
}




//////////////////////////////////////////////////////////////////////////////
// Myers-Ukkonen for semi-global edit-distance-alignments
// (version for needles longer than one machineword)
//////////////////////////////////////////////////////////////////////////////



template <typename TFinder, typename TNeedle, typename TSpec>
inline bool _findMyersLargePatterns (TFinder & finder, Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > > & me) 
{
SEQAN_CHECKPOINT
	typedef typename Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > >::TWord TWord;

	TWord X, D0, HN, HP, temp;
	TWord carryD0, carryHP, carryHN;
	unsigned shift, limit, currentBlock;

	while (!atEnd(finder)) {
		carryD0 = carryHN = 0;
		carryHP = _MyersUkkonenHP0<TSpec>::VALUE;

		// if the active cell is the last of it's block, one additional block has to be calculated
		limit = me.lastBlock + (me.scoreMask >> (me.MACHINE_WORD_SIZE - 1));

		if (limit == me.blockCount)
			limit--;

		shift = me.blockCount * ordValue((typename Value< TNeedle >::Type) *finder);

		// computing the necessary blocks, carries between blocks following one another are stored
		for (currentBlock = 0; currentBlock <= limit; currentBlock++) {
			X = me.bitMasks[shift + currentBlock] | me.VN[currentBlock];
	
			temp = me.VP[currentBlock] + (X & me.VP[currentBlock]) + carryD0;
			if (carryD0)
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
				if (HP & me.scoreMask)
					me.score++;
				else if (HN & me.scoreMask)
					me.score--;
			}
		}

		// updating the last active cell
		while (me.score > me.k) {
			if (me.VP[me.lastBlock] & me.scoreMask)
				me.score--;
			else if (me.VN[me.lastBlock] & me.scoreMask)
				me.score++;

			me.scoreMask >>= 1;
			if (!me.scoreMask) {
				me.scoreMask = 1 << (me.MACHINE_WORD_SIZE - 1);
				me.lastBlock--;
			}
		}

		if ((me.lastBlock == me.blockCount-1) && (me.scoreMask == me.finalScoreMask))
			return true;
		else {
			me.scoreMask <<= 1;
			if (!me.scoreMask) {
				me.scoreMask = 1;
				me.lastBlock++;
			}
			
			if (me.VP[me.lastBlock] & me.scoreMask)
				me.score++;
			else if (me.VN[me.lastBlock] & me.scoreMask)
				me.score--;
		}

		SEQAN_ASSERT (me.score >= 0);

		goNext(finder);
	}

	return false;
}

//____________________________________________________________________________
// version for needles not longer than one machineword

template <typename TFinder, typename TNeedle, typename TSpec>
inline bool _findMyersSmallPatterns (TFinder & finder, Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > > & me) 
{
SEQAN_CHECKPOINT

	typedef typename Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > >::TWord TWord;

	TWord X, D0, HN, HP;

	// computing the blocks
	while (!atEnd(finder)) {
		X = me.bitMasks[ordValue((typename Value<TNeedle>::Type) *finder)] | me.VN0;
		
		D0 = ((me.VP0 + (X & me.VP0)) ^ me.VP0) | X;
		HN = me.VP0 & D0;
		HP = me.VN0 | ~(me.VP0 | D0);
		X = (HP << 1) | _MyersUkkonenHP0<TSpec>::VALUE;
		me.VN0 = X & D0;
		me.VP0 = (HN << 1) | ~(X | D0);

		if (HP & (1 << (me.needleSize-1)))
			me.score++;
		else if (HN & (1 << (me.needleSize-1)))
			me.score--;

		if (me.score <= me.k)
			return true;
		
		goNext(finder);
	}

	return false;
}




//////////////////////////////////////////////////////////////////////////////
// Myers-Ukkonen as a banded alignment
// the band includes the main diagonal and the diagonals above
// the band width is (blockCount * MACHINE_WORD_SIZE)
//////////////////////////////////////////////////////////////////////////////



template <typename TFinder, typename TNeedle, typename TSpec>
inline bool 
_findMyersLargePatterns(
	TFinder & finder, 
	Pattern<TNeedle, Tag<_MyersUkkonen<AlignTextBanded> > > & me) 
{
SEQAN_CHECKPOINT
	typedef typename Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > >::TWord TWord;
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
				|= 1 << (me.MACHINE_WORD_SIZE - 1);

			goNext(me.ndlIter);
		}

		carryD0 = carryHN = 0;
		carryHP = _MyersUkkonenHP0<TSpec>::VALUE;

		// if the active cell is the last of it's block, one additional block has to be calculated
		limit = me.lastBlock + (me.scoreMask >> (me.MACHINE_WORD_SIZE - 1));

		if (limit == me.blockCount)
			limit--;

		shift = me.blockCount * ordValue((typename Value< TNeedle >::Type) *finder);

		// computing the necessary blocks, carries between blocks following one another are stored
		for (currentBlock = 0; currentBlock <= limit; currentBlock++) {
			X = me.bitMasks[shift + currentBlock] | me.VN[currentBlock];
	
			temp = me.VP[currentBlock] + (X & me.VP[currentBlock]) + carryD0;
			if (carryD0)
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
				if (HP & me.scoreMask)
					me.score++;
				else if (HN & me.scoreMask)
					me.score--;
			}
		}

		// updating the last active cell
		while (me.score > me.k) {
			if (me.VP[me.lastBlock] & me.scoreMask)
				me.score--;
			else if (me.VN[me.lastBlock] & me.scoreMask)
				me.score++;

			me.scoreMask >>= 1;
			if (!me.scoreMask) {
				me.scoreMask = 1 << (me.MACHINE_WORD_SIZE - 1);
				me.lastBlock--;
			}
		}

		if ((me.lastBlock == me.blockCount-1) && (me.scoreMask == me.finalScoreMask))
			return true;
		else {
			me.scoreMask <<= 1;
			if (!me.scoreMask) {
				me.scoreMask = 1;
				me.lastBlock++;
			}
			
			if (me.VP[me.lastBlock] & me.scoreMask)
				me.score++;
			else if (me.VN[me.lastBlock] & me.scoreMask)
				me.score--;
		}

		SEQAN_ASSERT (me.score >= 0);

		goNext(finder);
	}

	return false;
}

template <typename TFinder, typename TNeedle>
inline bool 
_findMyersSmallPatterns(
	TFinder & finder, 
	Pattern<TNeedle, Tag<_MyersUkkonen<AlignTextBanded> > > & me)
{
SEQAN_CHECKPOINT
	typedef typename Pattern<TNeedle, Tag<_MyersUkkonen<AlignTextBanded> > >::TWord TWord;

	TWord X, D0, HN, HP;
//	unsigned SHIFT=me.needleSize;

	unsigned minScore = SupremumValue<unsigned>::VALUE;
	TFinder minPos;

	if (!atEnd(me.ndlIter)) 
	{
		// Part 1: Go down the diagonal
		do
		{
			// shift bitmasks and states
			me.VN0 >>= 1;
			me.VP0 = (me.VP0 >> 1) | (1 << (me.MACHINE_WORD_SIZE - 1)); // ignore the field left of the leftmost diagonal
			for(unsigned i = 0; i < length(me.bitMasks); ++i)
				me.bitMasks[i] >>= 1;
			me.bitMasks[ordValue(*me.ndlIter)] |= (1 << (me.MACHINE_WORD_SIZE - 1));
			
//			--SHIFT;

			// Myers core
			X = me.bitMasks[ordValue((typename Value<TNeedle>::Type) *finder)] | me.VN0;
			D0 = ((me.VP0 + (X & me.VP0)) ^ me.VP0) | X;
			HN = me.VP0 & D0;
			HP = me.VN0 | ~(me.VP0 | D0);
			X = (HP << 1);
			me.VN0 = X & D0;
			me.VP0 = (HN << 1) | ~(X | D0);

/*
					for(int i=0; i<SHIFT; ++i) 
						::std::cout << "   ";
					::std::cout << "DD: ";
					for(int i=me.MACHINE_WORD_SIZE-1; i>=0 ;--i) 
					{
						CharString vd = " 1 ";
						if (D0 & (1 << i)) vd = " 0 ";
						::std::cout << vd;
					}
					::std::cout << "   ";
					::std::cout << *finder;
*/


			if (!(D0 & (1 << (me.MACHINE_WORD_SIZE - 1))))
				++me.score;

//			::std::cout << me.score <<::std::endl;

			goNext(me.ndlIter);
			if (atEnd(me.ndlIter))
			{
//				if (me.score <= me.k)
//					return true;
				minScore = me.score;
				minPos = finder;
				break;
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
		X = (HP << 1);
		me.VN0 = X & D0;
		me.VP0 = (HN << 1) | ~(X | D0);

/*
				for(int i=0; i<SHIFT; ++i) 
					::std::cout << "   ";
				::std::cout << "HD: ";
				for(int i=me.MACHINE_WORD_SIZE-1; i>=0 ;--i) 
				{
					CharString hd = " 0 ";
					if (HP & (1 << i)) hd = " 1 ";
					if (HN & (1 << i)) hd = "-1 ";
					::std::cout << hd;
				}
				::std::cout << "   ";
				::std::cout << *finder;
*/

		if (HP & (1 << (me.MACHINE_WORD_SIZE - 1)))
			++me.score;
		else if (HN & (1 << (me.MACHINE_WORD_SIZE - 1)))
			--me.score;

//		::std::cout << me.score <<::std::endl;
		if (me.score <= minScore) {
			minScore = me.score;
			minPos = finder;
		}
		
		goNext(finder);
	}

	if (minScore <= me.k) {
		finder = minPos;
		me.score = minScore;
		return true;
	}

	return false;
}




//////////////////////////////////////////////////////////////////////////////
// find
template <typename TFinder, typename TNeedle, typename TSpec>
inline bool find (TFinder & finder, 
				  Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > > & me)
{
SEQAN_CHECKPOINT
	int k = scoreLimit(me);

	if (empty(finder))
	{
		// in seqan k is treated as score, here we need it as penalty, that is why it is negated
		me.k = -k;

		_patternInit(me, finder);
		_finderSetNonEmpty(finder);

		//TODO: adapt myers-ukkonnen to dynamically change k

		// distinguish between the version for needles not longer than one machinword and the version for longer needles
		if (me.blockCount == 1) 
			return _findMyersSmallPatterns(finder, me);
		else 
			return _findMyersLargePatterns(finder, me);
	}
	else
	{
		goNext(finder);
		// distinguish between the version for needles not longer than one machineword and the version for longer needles
		if (me.blockCount == 1) 
			return _findMyersSmallPatterns(finder, me);
		else
			return _findMyersLargePatterns(finder, me);
	}
}

template <typename TFinder, typename TNeedle, typename TSpec>
inline bool find (TFinder & finder, 
				  Pattern<TNeedle, Tag<_MyersUkkonen<TSpec> > > & me, 
				  int const k)
{
	setScoreLimit(me, k);
	return find(finder, me);
}

//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
