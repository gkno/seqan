#ifndef SEQAN_HEADER_FIND_SHIFTAND_H
#define SEQAN_HEADER_FIND_SHIFTAND_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// ShiftAnd
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.ShiftAnd:
..summary: Exact string matching using bit-parallelism. Applicable to small strings and small alphabets.
..general:Class.Pattern
..cat:Pattern Matching
..signature:Pattern<TNeedle, ShiftAnd>
..param.TNeedle:The needle type.
...type:Class.String
..remarks.text:Needle-length and the size of the alphabet have to fit in a machine word.
..remarks.text:Types of needle and haystack have to match.
*/

///.Class.Pattern.param.TSpec.type:Spec.ShiftAnd

struct _ShiftAnd;
typedef Tag<_ShiftAnd> ShiftAnd;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, ShiftAnd> {
//____________________________________________________________________________
private:
	Pattern(Pattern const& other);
	Pattern const& operator=(Pattern const & other);

//____________________________________________________________________________
public:
	unsigned int* table;			// Look up table for each character in the alphabet (called B in "Navarro")
	unsigned int* prefSufMatch;		// Set of all the prefixes of needle that match a suffix of haystack (called D in "Navarro")
	unsigned int alphabetSize;		// e.g., char --> 256
	unsigned int needleLength;		// e.g., needleLength=33 --> blockCount=2 (iff w=32 bits)
	unsigned int blockCount;		// #unsigned ints required to store needle	

//____________________________________________________________________________

	Pattern() {
SEQAN_CHECKPOINT
		table = 0;
		prefSufMatch=0;
	}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl)
	{
		table = 0;
		prefSufMatch=0;
		setHost(*this, ndl);
	}

	~Pattern() {
		SEQAN_CHECKPOINT
		if (table != 0) {
			deallocate(this, table, alphabetSize * blockCount);
			deallocate(this, prefSufMatch, blockCount);
		}
	}		
//____________________________________________________________________________
};


template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, ShiftAnd> & me, TNeedle2 const & needle) {
	SEQAN_CHECKPOINT
	typedef typename Value<TNeedle>::Type TValue;
	if (me.table != 0) {
		deallocate(me, me.table, me.alphabetSize * me.blockCount);
		deallocate(me, me.prefSufMatch, me.blockCount);
	}

	
	me.needleLength = length(needle);
	me.alphabetSize = ValueSize<TValue>::VALUE;
	if (me.needleLength<1) me.blockCount=1;
	else me.blockCount=((me.needleLength-1) / BitsPerValue<unsigned int>::VALUE)+1;
			
	allocate (me, me.table, me.blockCount * me.alphabetSize);
	arrayFill (me.table, me.table + me.blockCount * me.alphabetSize, 0);

	allocate (me, me.prefSufMatch, me.blockCount);
	arrayFill (me.prefSufMatch, me.prefSufMatch + me.blockCount, 0);

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
void setHost (Pattern<TNeedle, ShiftAnd> & me, TNeedle2 & needle)
{
	setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}



template <typename TFinder, typename TNeedle>
bool _findShiftAnd_SmallNeedle(TFinder & finder, Pattern<TNeedle, ShiftAnd> & me) {
	SEQAN_CHECKPOINT
	unsigned int compare = (1 << (me.needleLength-1));
	while (!atEnd(finder)) {
		unsigned int pos = static_cast<unsigned int>(*finder);
		me.prefSufMatch[0] = ((me.prefSufMatch[0] << 1) | 1) & me.table[me.blockCount*pos];
		if ((me.prefSufMatch[0] & compare) != 0) {
			finder-=(me.needleLength-1);
			return true; 
		}
		goNext(finder);
	}
	return 0;
}

template <typename TFinder, typename TNeedle>
bool _findShiftAnd_LargeNeedle(TFinder & finder, Pattern<TNeedle, ShiftAnd> & me) {
	SEQAN_CHECKPOINT
	unsigned int compare = (1 << ((me.needleLength-1) % BitsPerValue<unsigned int>::VALUE));
	while (!atEnd(finder)) {
		unsigned int pos = static_cast<unsigned int>(*finder);
		unsigned int carry = 1;
		for(unsigned int block=0;block<me.blockCount;++block) {
			bool newCarry = ((me.prefSufMatch[block] & (1<< (BitsPerValue<unsigned int>::VALUE - 1)))!=0); 
			me.prefSufMatch[block]<<=1;
			me.prefSufMatch[block]|=carry;
			carry = newCarry;
		}
		for(unsigned int block=0;block<me.blockCount;++block) me.prefSufMatch[block] &= me.table[me.blockCount*pos+block];
		if ((me.prefSufMatch[me.blockCount-1] & compare) != 0) {
			finder-=(me.needleLength-1);
			return true; 
		}

		/*
		// Debug code
		std::cout << "   ";
		for(int j=0;j<me.blockCount;++j) {
			for(int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
				std::cout << ((me.prefSufMatch[j] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
			}
		}
		std::cout << std::endl;
		*/
		goNext(finder);
	}
	return 0;
}

template <typename TFinder, typename TNeedle>
inline bool find(TFinder & finder, Pattern<TNeedle, ShiftAnd> & me) {
	SEQAN_CHECKPOINT

	if (empty(finder))
		goBegin(finder);
	else
		finder += me.needleLength;

	// Fast algorithm for needles < machine word?
	if (me.blockCount == 1) {
		return _findShiftAnd_SmallNeedle(finder, me);
	} else {
		return _findShiftAnd_LargeNeedle(finder, me);
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_SHIFTAND_H
