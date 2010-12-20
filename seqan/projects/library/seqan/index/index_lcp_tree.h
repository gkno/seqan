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

#ifndef SEQAN_HEADER_INDEX_LCP_TREE_H
#define SEQAN_HEADER_INDEX_LCP_TREE_H

namespace SEQAN_NAMESPACE_MAIN
{
	
	template <
		class LCPFwdIt,		// lcp table input iterator
		class FlatOutIt >	// flat tree output iterator
	inline FlatOutIt createLcpBinTree(
		LCPFwdIt _First, LCPFwdIt _Last,
		FlatOutIt _Dest)
	{
        typedef typename Value<LCPFwdIt>::Type  TValue;
        typedef typename Size<LCPFwdIt>::Type   TSize;

        TSize size = difference(_First, _Last);
        if (size <= 1) return _Dest;
		--size;

        // calculate the depth of the lcp tree
		unsigned treeLevels = 1;
		TSize _xSize = 1;
		for(; size > _xSize; _xSize *= 2, ++treeLevels) ;
	
		// get output iterators for every level in the flat tree
		FlatOutIt *level = new FlatOutIt[treeLevels];
		for(unsigned i = treeLevels - 1; _xSize; --i, _xSize /= 2) {
			level[i] = _Dest;
			goFurther(_Dest, (size + _xSize - 1) / _xSize);
		}

		// fields to keep track of minimum elements and state
		TValue *minVal = new TValue[treeLevels];
		bool *half = new bool[treeLevels];
		for(unsigned i = 0; i < treeLevels; ++i)
			half[i] = false;

		// it works like a binary counter of half[n]...half[1]half[0]
		for(TSize j = 0; j < size; ++j, ++_First) {
			*(level[0]) = minVal[0] = *_First;
			++(level[0]);
			for(unsigned i = 1; i < treeLevels; ++i) {
				if (half[i]) {
					if (minVal[i-1] < minVal[i]) minVal[i] = minVal[i-1];
					*(level[i]) = minVal[i];	// min[i] is the minimum of last 2 values in min[i-1]
					++(level[i]);
					half[i] = false;
				} else {
					minVal[i] = minVal[i-1];
					half[i] = true;
					break;
				}
			}
		}

		// complete half filled nodes
		bool carry = false;
		for(unsigned i = 1; i < treeLevels; ++i)
			if (half[i] || carry) {
				if (half[i]) {
					if (minVal[i-1] < minVal[i]) minVal[i] = minVal[i-1];
				} else
					minVal[i] = minVal[i-1];
				*(level[i]) = minVal[i];
				++(level[i]);
				carry = true;
			}

		// trailing zero
		*_Dest = 0;
		++_Dest;

		delete[] half;
		delete[] minVal;
		delete[] level;

		return _Dest;
    }



	template < typename TSize >
	inline TSize sizeofLcpe(TSize n)
	{
		if (n < 2) return n;	// 0 -> 0, 1 -> 1, 2 -> 2, 3 -> 4
		--n;
		TSize size = 2;
		for(TSize _xSize = 1; _xSize < n; _xSize *= 2)
			size += (n + _xSize - 1) / _xSize;
		return size;
	}

	template < typename TSize >
	inline TSize sizeofLcph(TSize n)
	{
		return sizeofLcpe(n);
	}

	template <
		class LCPFwdIt,		// lcp table input iterator
		typename TSize >
	inline void sizeofLcpe(LCPFwdIt _First, LCPFwdIt _Last, TSize &_Size)
	{
		_Size = sizeofLcpe(difference(_First, _Last));
	}

	template <
		class LCPFwdIt,		// lcp table input iterator
		typename TSize >
	inline void sizeofLcph(LCPFwdIt _First, LCPFwdIt _Last, TSize &_Size)
	{
		sizeofLcpe(_First, _Last, _Size);
		return;
	}


    template < typename TLCPE, typename TLCP >
    inline void createLcpBinTree(TLCPE &lcp_enhanced, TLCP &lcp) {
        createLcpBinTree(begin(lcp, Standard()), end(lcp, Standard()), begin(lcp_enhanced, Standard()));
    }


	template < typename TSize >
    inline unsigned _treeLevels(TSize lcpSize)
	{
		unsigned treeLevels = 1;
		--lcpSize;
		TSize _xSize = 1;
		for(; lcpSize > _xSize; _xSize *= 2, ++treeLevels) ;
        return treeLevels;
    }

    template < typename TValue, typename TConfig, typename TLCP >
    inline void createLcpBinTree(String<TValue, External<TConfig> > &lcp_enhanced, TLCP &lcp) {
        unsigned writeHeads = _treeLevels(length(lcp)) + 1;   // plus 1 write back buffer
        if (lcp_enhanced.cache.size() < writeHeads)
            lcp_enhanced.resizeCache(writeHeads);
        createLcpBinTree(begin(lcp, Standard()), end(lcp, Standard()), begin(lcp_enhanced));
    }

}

#endif
