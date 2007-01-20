#ifndef SEQAN_HEADER_FIND_MYERS_UKKONEN_H
#define SEQAN_HEADER_FIND_MYERS_UKKONEN_H

namespace SEQAN_NAMESPACE_MAIN {

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

#define MACHINE_WORD_SIZE (sizeof(int) * 8)

private:
	unsigned int needleSize;
	unsigned int score;				// the current score
	unsigned int alphabetSize;
	unsigned int blockCount;		// the number of blocks
	unsigned int k;					// the maximal number of differences allowed
	unsigned int lastBlock;			// the block containing the last active cell

	unsigned int * VP;
	unsigned int * VN;
	unsigned int * bitMasks;		// encoding the alphabet as bit-masks
	unsigned int scoreMask;			// the mask with a bit set at the position of the last active cell
	unsigned int finalScoreMask;	// a mask with a bit set on the position of the last row

//	void * compareClassIdentifier;	// stores the type of compareType. used as a flag whether to init the needle or not

//____________________________________________________________________________

public:
	Pattern () {
SEQAN_CHECKPOINT

		VP = 0;
		VN = 0;
		bitMasks = 0;
//		compareClassIdentifier = 0;
	}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl)
	{
		VP = 0;
		VN = 0;
		bitMasks = 0;
		setHost(*this, ndl);
	}

	~Pattern () {
SEQAN_CHECKPOINT
		
		if (bitMasks != 0) {
			deallocate(this, bitMasks, alphabetSize * blockCount);
			deallocate(this, VP, blockCount);
			deallocate(this, VN, blockCount);
		}
	}
//____________________________________________________________________________

	template <typename TNeedle2>
	friend void
	setHost(Pattern & me, TNeedle2 const & needle) {
SEQAN_CHECKPOINT

		if (me.bitMasks != 0) {
			deallocate(me, me.bitMasks, me.alphabetSize * me.blockCount);
			deallocate(me, me.VP, me.blockCount);
			deallocate(me, me.VN, me.blockCount);
		}

		me.needleSize = length(needle);
		me.alphabetSize = ValueSize< typename Value<TNeedle>::Type >::VALUE;
		me.blockCount = (me.needleSize + MACHINE_WORD_SIZE - 1) / MACHINE_WORD_SIZE;

		allocate (me, me.VP, me.blockCount);
		arrayFill (me.VP, me.blockCount, ~0);

		allocate (me, me.VN, me.blockCount);
		arrayFill (me.VN, me.blockCount, 0);

		allocate (me, me.bitMasks, me.alphabetSize * me.blockCount);
		arrayFill(me.bitMasks, me.alphabetSize * me.blockCount, 0);

		// encoding the letters as bit-vectors
        for (unsigned int j = 0; j < me.needleSize; j++)
			me.bitMasks[me.blockCount * static_cast<unsigned int>((typename Value<TNeedle>::Type) value(needle,j)) + j/MACHINE_WORD_SIZE] = me.bitMasks[me.blockCount * static_cast<unsigned int>((typename Value<TNeedle>::Type) value(needle,j)) + j/MACHINE_WORD_SIZE] | 1 << (j%MACHINE_WORD_SIZE);
			//me.bitMasks[me.blockCount * static_cast<unsigned int>((CompareType< Value< TNeedle >::Type, Value< Container< THaystack >::Type >::Type >::Type) needle[j]) + j/MACHINE_WORD_SIZE] = me.bitMasks[me.blockCount * static_cast<unsigned int>((CompareType< Value< TNeedle >::Type, Value< Container< THaystack >::Type >::Type >::Type) needle[j]) + j/MACHINE_WORD_SIZE] | 1 << (j%MACHINE_WORD_SIZE);
		
//		me.compareClassIdentifier = _ClassIdentifier< CompareType< Value< TNeedle >::Type, Value< Container< THaystack >::Type >::Type >::Type >::getID();

	}
//____________________________________________________________________________

// version for needles longer than one machineword
	template <typename TFinder>
	friend bool
	_findMyersLargePatterns (TFinder & finder, Pattern & me) {
SEQAN_CHECKPOINT

		unsigned int X, D0, HN, HP, temp, shift, limit, currentBlock;
		unsigned int carryD0, carryHP, carryHN;

		while (!atEnd(finder)) {
			carryD0 = carryHP = carryHN = 0;

			// if the active cell is the last of it's block, one additional block has to be calculated
			limit = me.lastBlock + (me.scoreMask >> (MACHINE_WORD_SIZE - 1));

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
				carryHP = HP >> (MACHINE_WORD_SIZE-1);
				
				me.VN[currentBlock] = X & D0;

				temp = (HN << 1) | carryHN;
				carryHN = HN >> (MACHINE_WORD_SIZE - 1);
									
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
					me.scoreMask = 1 << (MACHINE_WORD_SIZE - 1);
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
	template <typename TFinder>
	friend bool
	_findMyersSmallPatterns (TFinder & finder, Pattern & me) {
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

	template <typename TFinder>
	friend bool 
	find (TFinder & finder, Pattern & me, int const k) {
SEQAN_CHECKPOINT

		if (empty(finder))
		{
			goBegin(finder);

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
				me.scoreMask = 1 << (me.k % MACHINE_WORD_SIZE);
				me.lastBlock = me.k/MACHINE_WORD_SIZE; 
				if (me.lastBlock == me.blockCount)
					me.lastBlock--;
				me.finalScoreMask = 1 << ((me.needleSize + MACHINE_WORD_SIZE -1 ) % MACHINE_WORD_SIZE);
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
//____________________________________________________________________________

/*
	friend inline void reset4find(Pattern & me)
	{
		if(me.compareClassIdentifier != 0)
		{
			allocate (me, me.VP, me.blockCount);
			arrayFill (me.VP, me.blockCount, ~0);

			allocate (me, me.VN, me.blockCount);
			arrayFill (me.VN, me.blockCount, 0);
		}
	}
*/
//____________________________________________________________________________
};

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
