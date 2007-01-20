#ifndef SEQAN_HEADER_FIND_BNDMALGO_H
#define SEQAN_HEADER_FIND_BNDMALGO_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// BndmAlgo
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.BndmAlgo:
..summary: Exact string matching using bit-parallelism. Applicable to small strings and small alphabets.
..general:Class.Pattern
..cat:Pattern Matching
..signature:Pattern<TNeedle, BndmAlgo>
..param.TNeedle:The needle type.
...type:Class.String
..remarks.text:Needle-length and the size of the alphabet have to fit in a machine word.
..remarks.text:Types of needle and haystack have to match.
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
	unsigned int* table;			// Look up table for each character in the alphabet (called B in "Navarro")
	unsigned int* activeFactors;	// The active factors in the pattern (called D in "Navarro")
	unsigned int alphabetSize;		// e.g., char --> 256
	unsigned int needleLength;		// e.g., needleLength=33 --> blockCount=2 (iff w=32 bits)
	unsigned int haystackLength;	// Length of haystack
	unsigned int blockCount;		// #unsigned ints required to store needle
	unsigned int last;

//____________________________________________________________________________

	Pattern() {	
SEQAN_CHECKPOINT
		table = 0;
		activeFactors=0;
	}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl)
	{
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
		

template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, BndmAlgo> & me, TNeedle2 const& needle) {
	SEQAN_CHECKPOINT
	typedef typename Value<TNeedle>::Type TValue;
	if (me.table != 0) {
		deallocate(me, me.table, me.alphabetSize * me.blockCount);
		deallocate(me, me.activeFactors, me.blockCount);
	}

	me.needleLength = length(needle);
	me.alphabetSize = ValueSize<TValue>::VALUE;
	if (me.needleLength<1) me.blockCount=1;
	else me.blockCount=((me.needleLength-1) / BitsPerValue<unsigned int>::VALUE)+1;
			
	allocate (me, me.table, me.blockCount * me.alphabetSize);
	arrayFill (me.table, me.table + me.blockCount * me.alphabetSize, 0);

	allocate (me, me.activeFactors, me.blockCount);
	arrayFill (me.activeFactors, me.activeFactors + me.blockCount, 0);

	for (unsigned int j = 0; j < me.needleLength; ++j) {
		// Determine character position in array table
		unsigned int pos = static_cast<unsigned int>(getValue(needle,j));
		me.table[me.blockCount*pos + j / BitsPerValue<unsigned int>::VALUE] |= (1<<(j%BitsPerValue<unsigned int>::VALUE));
	}

	/*
	// Debug code
	std::cout << "Alphabet size: " << me.alphabetSize << std::endl;
	std::cout << "Needle length: " << me.needleLength << std::endl;
	std::cout << "Block count: " << me.blockCount << std::endl;

	for(unsigned int i=0;i<me.alphabetSize;++i) {
		if ((i<97) || (i>122)) continue;
		std::cout << static_cast<char>(i) << ": ";
		for(int j=0;j<me.blockCount;++j) {
			for(int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
				std::cout << ((me.table[me.blockCount*i+j] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
			}
		}
		std::cout << std::endl;
	}
	*/
}

template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, BndmAlgo> & me, TNeedle2 & needle)
{
	setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}


template <typename TFinder, typename TNeedle>
bool _findBndm_SmallNeedle(TFinder & finder, Pattern<TNeedle, BndmAlgo> & me) {
	SEQAN_CHECKPOINT
	while (position(finder) <= me.haystackLength - me.needleLength) {
		me.last=me.needleLength;
		unsigned int j=me.needleLength;
		me.activeFactors[0]=~0;
		while (me.activeFactors[0]!=0) {
			unsigned int pos = static_cast<unsigned int>(*(finder+j-1));
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
	return 0;
}

template <typename TFinder, typename TNeedle>
bool _findBndm_LargeNeedle(TFinder & finder, Pattern<TNeedle, BndmAlgo> & me) {
	SEQAN_CHECKPOINT
	unsigned int carryPattern = (1<< (BitsPerValue<unsigned int>::VALUE - 1));
	while (position(finder) <= me.haystackLength - me.needleLength) {
		me.last=me.needleLength;
		unsigned int j=me.needleLength;
		for(unsigned int block=0;block<me.blockCount;++block) me.activeFactors[block]=~0;
		bool zeroPrefSufMatch = false;
		while (!zeroPrefSufMatch) {
			unsigned int pos = static_cast<unsigned int>(*(finder+j-1));

			/*	
			// Debug code
			std::cout << "   ";
			for(int j=0;j<me.blockCount;++j) {
				for(int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
					std::cout << ((me.activeFactors[j] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
				}
			}
			std::cout << std::endl;
			*/

			for(unsigned int block=0;block<me.blockCount;++block) me.activeFactors[block] &= me.table[me.blockCount*pos+block];
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
	return 0;
}

template <typename TFinder, typename TNeedle>
bool find(TFinder & finder, Pattern<TNeedle, BndmAlgo> & me) {
	SEQAN_CHECKPOINT

	if (empty(finder)) {
		goBegin(finder);
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
