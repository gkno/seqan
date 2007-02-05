#ifndef SEQAN_HEADER_FIND_SHIFTOR_H
#define SEQAN_HEADER_FIND_SHIFTOR_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// ShiftOr
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.ShiftOr:
..summary: Exact string matching using bit parallelism. The Shift-Or algorithm is applicable to search small patterns in texts using a small alphabet.
..general:Class.Pattern
..cat:Pattern Matching
..signature:Pattern<TNeedle, ShiftOr>
..param.TNeedle:The needle type.
...type:Class.String
..remarks.text:The types of the needle and the haystack have to match.
*/

///.Class.Pattern.param.TSpec.type:Spec.ShiftOr

struct _ShiftOr;
typedef Tag<_ShiftOr> ShiftOr;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, ShiftOr> {
//____________________________________________________________________________
private:
	Pattern(Pattern const& other);
	Pattern const& operator=(Pattern const & other);

//____________________________________________________________________________
public:
	typedef typename Size<TNeedle>::Type TSize;
	typedef typename Value<TNeedle>::Type TAlphabet;
	typedef typename Size<TAlphabet>::Type TAlphabetSize;
	TAlphabetSize* table;			// Look up table for each character in the alphabet (called B in "Navarro")
	TAlphabetSize* prefSufMatch;		// Set of all the prefixes of needle that match a suffix of haystack (called D in "Navarro")
	TAlphabetSize alphabetSize;		// e.g., char --> 256
	TSize needleLength;		// e.g., needleLength=33 --> blockCount=2 (iff w=32 bits)
	TSize blockCount;		// #unsigned ints required to store needle	

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
void setHost (Pattern<TNeedle, ShiftOr> & me, TNeedle2 const & needle) {
	SEQAN_CHECKPOINT
	typedef typename Value<TNeedle>::Type TValue;
	typedef typename Size<TNeedle>::Type TSize;
	if (me.table != 0) {
		deallocate(me, me.table, me.alphabetSize * me.blockCount);
		deallocate(me, me.prefSufMatch, me.blockCount);
	}

	
	me.needleLength = length(needle);
	me.alphabetSize = ValueSize<TValue>::VALUE;
	if (me.needleLength<1) me.blockCount=1;
	else me.blockCount=((me.needleLength-1) / BitsPerValue<TSize>::VALUE)+1;
			
	allocate (me, me.table, me.blockCount * me.alphabetSize);
	arrayFill (me.table, me.table + me.blockCount * me.alphabetSize, ~0);

	allocate (me, me.prefSufMatch, me.blockCount);
	arrayFill (me.prefSufMatch, me.prefSufMatch + me.blockCount, ~0);

	for (TSize j = 0; j < me.needleLength; ++j) {
		// Determine character position in array table
		TSize pos = convert<TSize>(getValue(needle,j));
		me.table[me.blockCount*pos + j / BitsPerValue<TSize>::VALUE] ^= (1<<(j%BitsPerValue<TSize>::VALUE));
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
void setHost (Pattern<TNeedle, ShiftOr> & me, TNeedle2 & needle)
{
	setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}



template <typename TFinder, typename TNeedle>
bool _findShiftOr_SmallNeedle(TFinder & finder, Pattern<TNeedle, ShiftOr> & me) {
	SEQAN_CHECKPOINT
	typedef typename Size<TNeedle>::Type TSize;
	TSize compare= (~(1 << (me.needleLength-1)));
	while (!atEnd(finder)) {
		TSize pos = convert<TSize>(*finder);
		me.prefSufMatch[0] = (me.prefSufMatch[0] << 1) | me.table[me.blockCount*pos];
		if ((me.prefSufMatch[0] | compare) != ~0) {
			finder-=(me.needleLength-1);
			return true; 
		}
		goNext(finder);
	}
	return 0;
}

template <typename TFinder, typename TNeedle>
bool _findShiftOr_LargeNeedle(TFinder & finder, Pattern<TNeedle, ShiftOr> & me) {
	SEQAN_CHECKPOINT
	typedef typename Size<TNeedle>::Type TSize;
	TSize compare= (~(1 << ((me.needleLength-1) % BitsPerValue<TSize>::VALUE)));
	while (!atEnd(finder)) {
		TSize pos = convert<TSize>(*finder);
		TSize carry = 0;
		for(TSize block=0;block<me.blockCount;++block) {
			bool newCarry=((me.prefSufMatch[block] & (1<< (BitsPerValue<TSize>::VALUE - 1)))!=0); 
			me.prefSufMatch[block]<<=1;
			me.prefSufMatch[block]|=carry;
			carry=newCarry;
		}
		for(TSize block=0;block<me.blockCount;++block) me.prefSufMatch[block] |= me.table[me.blockCount*pos+block];
		if ((me.prefSufMatch[me.blockCount-1] | compare) != ~0) {
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
inline bool find(TFinder & finder, Pattern<TNeedle, ShiftOr> & me) {
	SEQAN_CHECKPOINT

	if (empty(finder))
		goBegin(finder);
	else
		finder += me.needleLength;

	// Fast algorithm for needles < machine word?
	if (me.blockCount == 1) {
		return _findShiftOr_SmallNeedle(finder, me);
	} else {
		return _findShiftOr_LargeNeedle(finder, me);
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_SHIFTOR_H
