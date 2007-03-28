#ifndef SEQAN_HEADER_FIND_BNDMALGO_H
#define SEQAN_HEADER_FIND_BNDMALGO_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// BndmAlgo
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.BndmAlgo:
..summary: Backward Nondeterministic Dawg Matching algorithm. Exact string matching using bit parallelism.
..general:Class.Pattern
..cat:Pattern Matching
..signature:Pattern<TNeedle, BndmAlgo>
..param.TNeedle:The needle type.
...type:Class.String
..remarks.text:The types of the needle and the haystack have to match.
*/

///.Class.Pattern.param.TSpec.type:Spec.BndmAlgo

struct _BndmAlgo;
typedef Tag<_BndmAlgo> BndmAlgo;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, BndmAlgo> {
//____________________________________________________________________________
private:
	Pattern(Pattern const& other);
	Pattern const& operator=(Pattern const & other);

//____________________________________________________________________________
public:
	typedef unsigned int TWord;
	Holder<TNeedle> data_needle;
	TWord* table;			// Look up table for each character in the alphabet (called B in "Navarro")
	TWord* activeFactors;	// The active factors in the pattern (called D in "Navarro")
	TWord alphabetSize;		// e.g., char --> 256
	TWord needleLength;		// e.g., needleLength=33 --> blockCount=2 (iff w=32 bits)
	TWord haystackLength;	// Length of haystack
	TWord blockCount;		// #unsigned ints required to store needle
	TWord last;

//____________________________________________________________________________

	Pattern() {	
		table = 0;
		activeFactors=0;
	}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl)
	{
		SEQAN_CHECKPOINT
		table = 0;
		activeFactors=0;
		setHost(*this, ndl);
	}

	~Pattern() {
		SEQAN_CHECKPOINT
		if (table != 0) {
			deallocate(this, table, alphabetSize * blockCount);
			deallocate(this, activeFactors, blockCount);
		}
	}

//____________________________________________________________________________
};
		
//////////////////////////////////////////////////////////////////////////////
// Host Metafunctions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
struct Host< Pattern<TNeedle, BndmAlgo> >
{
	typedef TNeedle Type;
};

template <typename TNeedle>
struct Host< Pattern<TNeedle, BndmAlgo> const>
{
	typedef TNeedle const Type;
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, BndmAlgo> & me, TNeedle2 const& needle) {
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef typename Value<TNeedle>::Type TValue;
	if (me.table != 0) {
		deallocate(me, me.table, me.alphabetSize * me.blockCount);
	}

	me.needleLength = length(needle);
	me.alphabetSize = ValueSize<TValue>::VALUE;
	if (me.needleLength<1) me.blockCount=1;
	else me.blockCount=((me.needleLength-1) / BitsPerValue<TWord>::VALUE)+1;
			
	allocate (me, me.table, me.blockCount * me.alphabetSize);
	arrayFill (me.table, me.table + me.blockCount * me.alphabetSize, 0);

	for (TWord j = 0; j < me.needleLength; ++j) {
		// Determine character position in array table
		TWord pos = convert<TWord>(getValue(needle,j));
		me.table[me.blockCount*pos + j / BitsPerValue<TWord>::VALUE] |= (1<<(j%BitsPerValue<TWord>::VALUE));
	}

	me.data_needle = needle;
	/*
	// Debug code
	std::cout << "Alphabet size: " << me.alphabetSize << ::std::endl;
	std::cout << "Needle length: " << me.needleLength << ::std::endl;
	std::cout << "Block count: " << me.blockCount << ::std::endl;

	for(unsigned int i=0;i<me.alphabetSize;++i) {
		if ((i<97) || (i>122)) continue;
		std::cout << static_cast<char>(i) << ": ";
		for(int j=0;j<me.blockCount;++j) {
			for(int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
				std::cout << ((me.table[me.blockCount*i+j] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
			}
		}
		std::cout << ::std::endl;
	}
	*/
}

template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, BndmAlgo> & me, TNeedle2 & needle)
{
	setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//____________________________________________________________________________


template <typename TNeedle>
inline void _finderInit (Pattern<TNeedle, BndmAlgo> & me) 
{
SEQAN_CHECKPOINT
	typedef unsigned int TWord;

	if (me.activeFactors != 0) {
		deallocate(me, me.activeFactors, me.blockCount);
	}
	me.last = 0;
	allocate (me, me.activeFactors, me.blockCount);
	arrayFill (me.activeFactors, me.activeFactors + me.blockCount, 0);
}

//____________________________________________________________________________

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, BndmAlgo>const>::Type & 
host(Pattern<TNeedle, BndmAlgo> & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, BndmAlgo>const>::Type & 
host(Pattern<TNeedle, BndmAlgo> const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

//____________________________________________________________________________


template <typename TFinder, typename TNeedle>
bool _findBndm_SmallNeedle(TFinder & finder, Pattern<TNeedle, BndmAlgo> & me) {
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	if (me.haystackLength < me.needleLength) return false;
	while (position(finder) <= me.haystackLength - me.needleLength) {
		me.last=me.needleLength;
		TWord j=me.needleLength;
		me.activeFactors[0]=~0;
		while (me.activeFactors[0]!=0) {
			TWord pos = convert<TWord>(*(finder+j-1));
			me.activeFactors[0] = (me.activeFactors[0] & me.table[me.blockCount*pos]);
			j--;
			if (me.activeFactors[0] & 1 != 0) {
				if (j>0) me.last=j;
				else return true;
			}
			me.activeFactors[0] = me.activeFactors[0] >> 1;
		}
		finder+=me.last;
	}
	return false;
}

template <typename TFinder, typename TNeedle>
bool _findBndm_LargeNeedle(TFinder & finder, Pattern<TNeedle, BndmAlgo> & me) {
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	TWord carryPattern = (1<< (BitsPerValue<TWord>::VALUE - 1));
	while (position(finder) <= me.haystackLength - me.needleLength) {
		me.last=me.needleLength;
		TWord j=me.needleLength;
		for(TWord block=0;block<me.blockCount;++block) me.activeFactors[block]=~0;
		bool zeroPrefSufMatch = false;
		while (!zeroPrefSufMatch) {
			TWord pos = convert<TWord>(*(finder+j-1));

			/*	
			// Debug code
			std::cout << "   ";
			for(int j=0;j<me.blockCount;++j) {
				for(int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
					std::cout << ((me.activeFactors[j] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
				}
			}
			std::cout << ::std::endl;
			*/

			for(TWord block=0;block<me.blockCount;++block) me.activeFactors[block] &= me.table[me.blockCount*pos+block];
			j--;
			if (me.activeFactors[0] & 1 != 0) {
				if (j>0) me.last=j;
				else return true;
			}
			bool carry=0;
			zeroPrefSufMatch=true;
			for(int block=me.blockCount-1;block>=0;--block) {
				bool newCarry=((me.activeFactors[block] & 1)!=0); 
				me.activeFactors[block]>>=1;
				if (carry) me.activeFactors[block]|=carryPattern;
				carry=newCarry;
				if (me.activeFactors[block]!=0) zeroPrefSufMatch=false;
			}
		}
		finder+=me.last;
	}
	return false;
}

template <typename TFinder, typename TNeedle>
bool find(TFinder & finder, Pattern<TNeedle, BndmAlgo> & me) {
	SEQAN_CHECKPOINT

	if (empty(finder)) {
		_finderInit(me);
		_finderSetNonEmpty(finder);
		me.haystackLength = length(container(finder));
	} else
		finder+=me.last;

	// Fast algorithm for needles < machine word?
	if (me.blockCount == 1) {
		return _findBndm_SmallNeedle(finder, me);
	} else {
		return _findBndm_LargeNeedle(finder, me);
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_SHIFTAND_H
