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

#ifndef SEQAN_HEADER_GRAPH_ALIGN_MYERS_H
#define SEQAN_HEADER_GRAPH_ALIGN_MYERS_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Meyer's Bit Vector algorithm
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TString>
typename Size<TString>::Type
_align_myers_bit_vector(TString const & str1,
						TString const & str2)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TString>::Type TAlphabet;
	typedef unsigned int TSize;
	typedef TSize TWord;
	typedef String<TWord> BitVector;
	typedef String<BitVector> TLookupTable;
	TSize alphLen = ValueSize<TAlphabet>::VALUE;

	// Preprocessing
	TLookupTable lT;
	resize(lT, alphLen);

	TSize len1 = length(str1);
	TSize len2 = length(str2);
	TSize blockCount;
	if (len2<1) blockCount=1;
	else blockCount=((len2-1) / BitsPerValue<TWord>::VALUE)+1;
	for(TSize i = 0; i < alphLen; ++i) {
		fill(lT[i],blockCount, 0, Exact());
	}
	for(TSize j = 0; j < len2; ++j) {
		TSize pos = convert<TWord>(getValue(str2,j));
		TSize block = j / BitsPerValue<TWord>::VALUE;
		(lT[pos])[block] |= (1<<(j%BitsPerValue<TWord>::VALUE));
	}
	BitVector VP;
	BitVector VN;
	fill(VP, blockCount, ~0, Exact() );
	fill(VN, blockCount, 0, Exact() );
	TWord err = len2;
	
	//// Debug code
	//std::cout << "Alphabet size: " << alphLen << ::std::endl;
	//std::cout << "Block count: " << blockCount << ::std::endl;
	//for(TSize i=0;i<alphLen;++i) {
	//	if ((i<97) || (i>122)) continue;
	//	std::cout << static_cast<char>(i) << ": ";
	//	for(TSize j=0;j<(TSize)blockCount;++j) {
	//		for(TSize bit_pos=0;bit_pos<BitsPerValue<TWord>::VALUE;++bit_pos) {
	//		  std::cout << ((lT[i][j] & (1<<(bit_pos % BitsPerValue<TSize>::VALUE))) !=0);
	//		}
	//	}
	//	std::cout << ::std::endl;
	//}

	BitVector X;	
	BitVector D0;
	BitVector HN;	
	BitVector HP;
	BitVector HNcopy;	
	BitVector HPcopy;
	resize(X, blockCount);
	resize(D0, blockCount);
	resize(HN, blockCount);
	resize(HP, blockCount);
	resize(HNcopy, blockCount);
	resize(HPcopy, blockCount);

	for(TSize col = 0; col<len1; ++col) {
		TWord pos = convert<TWord>(getValue(str1,col));	  
		// Addition might produce a carry
		bool carry = 0;
		bool newCarry;
		bool HPcarry = 1;  // We want edit distance!!!
		bool HPnewCarry;
		bool HNcarry = 0;
		bool HNnewCarry;
		for(TWord block=0;block<blockCount;++block) {
			X[block] = lT[pos][block] | VN[block];
			D0[block] = X[block] & VP[block];
			if (( (TSize) D0[block] + VP[block] < (TSize) D0[block] ) ||
				( (TSize) D0[block] + VP[block] < (TSize) VP[block] )) newCarry = 1;
			else newCarry = 0;
			D0[block] += VP[block];
			if ((carry) && ( (TSize) D0[block] == (TSize) ~0)) {
				if (newCarry) {
					std::cerr << "Two carries. Error !!!";
					exit(-1);
				} else {
					newCarry = 1;
				}
			}
			D0[block] += carry;
			carry = newCarry;
			D0[block] ^= VP[block];
			D0[block] |= X[block];
			HN[block] = VP[block] & D0[block];
			HP[block] = VN[block] | ~(VP[block] | D0[block]);
			HPcopy[block] = HP[block];
			HNcopy[block] = HN[block];
			HPnewCarry = ((HPcopy[block] & (1<< (BitsPerValue<TWord>::VALUE - 1)))!=0); 
			HPcopy[block] <<= 1;
			HPcopy[block] += HPcarry;
			HPcarry = HPnewCarry;
			X[block] = HPcopy[block];
			VN[block] = X[block] & D0[block];
			HNnewCarry = ((HNcopy[block] & (1<< (BitsPerValue<TWord>::VALUE - 1)))!=0); 
			HNcopy[block] <<= 1;
			HNcopy[block] += HNcarry;
			HNcarry = HNnewCarry;
			VP[block] = HNcopy[block] | ~(X[block] | D0[block]);
		}
		TSize finalBlock = (len2-1) / BitsPerValue<TWord>::VALUE;
		if ((HP[finalBlock] & (1<<((len2-1)%BitsPerValue<TWord>::VALUE))) != 0) ++err;
		else if ((HN[finalBlock] & (1<<((len2-1)%BitsPerValue<TWord>::VALUE))) != 0) --err;

		/*
		std::cout << err << std::endl;
		for(TSize block=0;block<(int)blockCount;++block) {
			for(TSize bit_pos=0;bit_pos<BitsPerValue<TWord>::VALUE;++bit_pos) {
				std::cout << ((D0[block] & (1<<(bit_pos % BitsPerValue<TSize>::VALUE))) !=0);
			}
		}
		std::cout << ::std::endl;
		for(TSize block=0;block<(int)blockCount;++block) {
			for(TSize bit_pos=0;bit_pos<BitsPerValue<TWord>::VALUE;++bit_pos) {
				std::cout << ((HN[block] & (1<<(bit_pos % BitsPerValue<TSize>::VALUE))) !=0);
			}
		}
		std::cout << ::std::endl;
		for(TSize block=0;block<(int)blockCount;++block) {
			for(TSize bit_pos=0;bit_pos<BitsPerValue<TWord>::VALUE;++bit_pos) {
				std::cout << ((HP[block] & (1<<(bit_pos % BitsPerValue<TSize>::VALUE))) !=0);
			}
		}
		std::cout << ::std::endl;
		for(TSize block=0;block<(int)blockCount;++block) {
			for(TSize bit_pos=0;bit_pos<BitsPerValue<TWord>::VALUE;++bit_pos) {
				std::cout << ((VN[block] & (1<<(bit_pos % BitsPerValue<TSize>::VALUE))) !=0);
			}
		}
		std::cout << ::std::endl;
		for(TSize block=0;block<(int)blockCount;++block) {
			for(TSize bit_pos=0;bit_pos<BitsPerValue<TWord>::VALUE;++bit_pos) {
				std::cout << ((VP[block] & (1<<(bit_pos % BitsPerValue<TSize>::VALUE))) !=0);
			}
		}
		std::cout << ::std::endl;
		std::cout << ::std::endl;
		*/

	}
	return err;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScoreValue>
TScoreValue
_globalAlignment(TStringSet const& str,
		 Score<TScoreValue, Simple> const&,
		 MyersBitVector)
{
	SEQAN_CHECKPOINT
	return _align_myers_bit_vector(str[0], str[1]);	
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
