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

struct _MyersUkkonen;
typedef Tag<_MyersUkkonen> MyersUkkonen;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, MyersUkkonen> {
//____________________________________________________________________________
public:
	typedef unsigned int TWord;
	enum { MACHINE_WORD_SIZE = sizeof(TWord) * 8 };

	TWord needleSize;
	TWord score;				// the current score
	TWord blockCount;			// the number of blocks
	TWord k;					// the maximal number of differences allowed
	TWord lastBlock;			// the block containing the last active cell

	String<TWord> VP;
	String<TWord> VN;
	String<TWord> bitMasks;		// encoding the alphabet as bit-masks
	TWord scoreMask;			// the mask with a bit set at the position of the last active cell
	TWord finalScoreMask;		// a mask with a bit set on the position of the last row

//____________________________________________________________________________

	Pattern() {}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl)
	{
		setHost(*this, ndl);
	}

//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TNeedle2>
void setHost(Pattern<TNeedle, MyersUkkonen> & me, TNeedle2 const & needle) 
{
SEQAN_CHECKPOINT

	typedef typename Value<TNeedle>::Type TValue;

	me.needleSize = length(needle);
	me.blockCount = (me.needleSize + me.MACHINE_WORD_SIZE - 1) / me.MACHINE_WORD_SIZE;

	clear(me.bitMasks);
	fill(me.bitMasks, ValueSize<TValue>::VALUE * me.blockCount, 0, Exact());

	// encoding the letters as bit-vectors
    for (unsigned int j = 0; j < me.needleSize; j++)
		me.bitMasks[me.blockCount * static_cast<unsigned int>((typename Value<TNeedle>::Type) value(needle,j)) + j/me.MACHINE_WORD_SIZE] = me.bitMasks[me.blockCount * static_cast<unsigned int>((typename Value<TNeedle>::Type) value(needle,j)) + j/me.MACHINE_WORD_SIZE] | 1 << (j%me.MACHINE_WORD_SIZE);
		//me.bitMasks[me.blockCount * static_cast<unsigned int>((CompareType< Value< TNeedle >::Type, Value< Container< THaystack >::Type >::Type >::Type) needle[j]) + j/me.MACHINE_WORD_SIZE] = me.bitMasks[me.blockCount * static_cast<unsigned int>((CompareType< Value< TNeedle >::Type, Value< Container< THaystack >::Type >::Type >::Type) needle[j]) + j/MACHINE_WORD_SIZE] | 1 << (j%MACHINE_WORD_SIZE);
}

template <typename TNeedle, typename TNeedle2>
void setHost(Pattern<TNeedle, MyersUkkonen> & me, TNeedle2 & needle)
{
SEQAN_CHECKPOINT
	setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//____________________________________________________________________________


template <typename TNeedle>
void _finderInit(Pattern<TNeedle, MyersUkkonen> & me)
{
SEQAN_CHECKPOINT
	clear(me.VP);
	fill(me.VP, me.blockCount, (unsigned int) ~0, Exact());

	clear(me.VN);
	fill(me.VN, me.blockCount, 0, Exact());
}

//____________________________________________________________________________
// version for needles longer than one machineword

template <typename TFinder, typename TNeedle>
inline bool _findMyersLargePatterns (TFinder & finder, Pattern<TNeedle, MyersUkkonen> & me) 
{
SEQAN_CHECKPOINT

	unsigned int X, D0, HN, HP, temp, shift, limit, currentBlock;
	unsigned int carryD0, carryHP, carryHN;

	while (!atEnd(finder)) {
		carryD0 = carryHP = carryHN = 0;

		// if the active cell is the last of it's block, one additional block has to be calculated
		limit = me.lastBlock + (me.scoreMask >> (me.MACHINE_WORD_SIZE - 1));

		if (limit == me.blockCount)
			limit--;

		shift = me.blockCount * static_cast<unsigned int>((typename Value< TNeedle >::Type) *finder);

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
			carryHP = HP >> (me.MACHINE_WORD_SIZE-1);
			
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

template <typename TFinder, typename TNeedle>
inline bool _findMyersSmallPatterns (TFinder & finder, Pattern<TNeedle, MyersUkkonen> & me) 
{
SEQAN_CHECKPOINT

	unsigned int X, D0, HN, HP;

	// computing the blocks
	while (!atEnd(finder)) {
		X = me.bitMasks[static_cast<unsigned int>((typename Value<TNeedle>::Type) *finder)] | me.VN[0];
		
		D0 = ((me.VP[0] + (X & me.VP[0])) ^ me.VP[0]) | X;
		HN = me.VP[0] & D0;
		HP = me.VN[0] | ~(me.VP[0] | D0);
		X = HP << 1;
		me.VN[0] = X & D0;
		me.VP[0] = (HN << 1) | ~(X | D0);

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

//____________________________________________________________________________

template <typename TFinder, typename TNeedle>
inline bool find (TFinder & finder, Pattern<TNeedle, MyersUkkonen> & me, int const k)
{
SEQAN_CHECKPOINT

	if (empty(finder))
	{
		_finderInit(me);
		_finderSetNonEmpty(finder);

		// in seqan k is treated as score, here we need it as penalty, that is why it is negated
		me.k = -k;
		//TODO: adapt myers-ukkonnen to dynamically change k

		// distinguish between the version for needles not longer than one machinword and the version for longer needles
		if (me.blockCount == 1) 
		{
			me.score = me.needleSize;
			return _findMyersSmallPatterns(finder, me);
		} 
		else 
		{
			me.score = me.k+1;
			me.scoreMask = 1 << (me.k % me.MACHINE_WORD_SIZE);
			me.lastBlock = me.k/me.MACHINE_WORD_SIZE; 
			if (me.lastBlock == me.blockCount)
				me.lastBlock--;
			me.finalScoreMask = 1 << ((me.needleSize + me.MACHINE_WORD_SIZE -1 ) % me.MACHINE_WORD_SIZE);
			return _findMyersLargePatterns(finder, me);
		}
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

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
