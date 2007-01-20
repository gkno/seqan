/*
 *  lcp_tree.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_INDEX_LCP_TREE_H
#define SEQAN_HEADER_INDEX_LCP_TREE_H

namespace SEQAN_NAMESPACE_MAIN
{
	
	template <
		class LCPFwdIt,		// lcp table input iterator
		class FlatOutIt >	// flat tree output iterator
	inline FlatOutIt createLCPBinTree(
		LCPFwdIt _First, LCPFwdIt _Last,
		FlatOutIt _Dest)
	{
        typedef typename Value<LCPFwdIt>::Type  TValue;
        typedef typename Size<LCPFwdIt>::Type   TSize;

        TSize size = distance(_First, _Last);
        if (size <= 1) return _Dest;
		--size;

        // calculate the depth of the lcp tree
		int treeLevels = 1;
		TSize _xSize = 1;
		for(; size > _xSize; _xSize *= 2, ++treeLevels);
	
		// get output iterators for every level in the flat tree
		FlatOutIt *level = new FlatOutIt[treeLevels];
		for(int i = treeLevels - 1; _xSize; --i, _xSize /= 2) {
			level[i] = _Dest;
			goFurther(_Dest, (size + _xSize - 1) / _xSize);
		}

		// fields to keep track of minimum elements and state
		TValue *min = new TValue[treeLevels];
		bool *half = new bool[treeLevels];
		for(int i = 0; i < treeLevels; ++i)
			half[i] = false;

		// it works like a binary counter of half[n]...half[1]half[0]
		for(TSize j = 0; j < size; ++j, ++_First) {
			*(level[0]) = min[0] = *_First;
			++(level[0]);
			for(int i = 1; i < treeLevels; ++i) {
				if (half[i]) {
					if (min[i-1] < min[i]) min[i] = min[i-1];
					*(level[i]) = min[i];	// min[i] is the minimum of last 2 values in min[i-1]
					++(level[i]);
					half[i] = false;
				} else {
					min[i] = min[i-1];
					half[i] = true;
					break;
				}
			}
		}

		// complete half filled nodes
		bool carry = false;
		for(int i = 1; i < treeLevels; ++i)
			if (half[i] || carry) {
				if (half[i]) {
					if (min[i-1] < min[i]) min[i] = min[i-1];
				} else
					min[i] = min[i-1];
				*(level[i]) = min[i];
				++(level[i]);
				carry = true;
			}

		// trailing zero
		*_Dest = 0;
		++_Dest;

		delete[] half;
		delete[] min;
		delete[] level;

		return _Dest;
    }



	// the hybrid tree is an extension of an lcp search tree
	// it contains pairs of lcp and suffix array values
	// which are needed in the enhanced lcp binary search

	template <
		class LCPFwdIt,		// lcp table input iterator
		class SAFwdIt,		// suffix array table input iterator
		class FlatOutIt >	// flat tree output iterator
	inline FlatOutIt createHybridBinTree(
		LCPFwdIt _First, LCPFwdIt _Last,
		SAFwdIt _FirstSA,
		FlatOutIt _Dest)
	{
        typedef Pair<
            typename Value<LCPFwdIt>::Type,
            typename Value<SAFwdIt>::Type,
			Compressed >						TPair;
        typedef typename Size<LCPFwdIt>::Type   TSize;

        TSize size = distance(_First, _Last);
        if (size <= 1) return _Dest;
		--size;

        // if tree contains only no interval,
        // then we need 1 SA elements and no lcp entry
        if (!size) {
            TPair p;
            p.i1 = 0;
            p.i2 = *_FirstSA; ++_FirstSA;
            *_Dest = p;
            return ++_Dest;
        }

        // calculate the depth of the lcp tree
		int treeLevels = 1;
		TSize _xSize = 1;
		for(; size > _xSize; _xSize *= 2, ++treeLevels);
		
		// get output iterators for every level in the flat tree
		FlatOutIt *level = new FlatOutIt[treeLevels];
		for(int i = treeLevels - 1; _xSize; --i, _xSize /= 2) {
			level[i] = _Dest;
            goFurther(_Dest, (size + _xSize - 1) / _xSize);
		}

		// fields to keep track of minimum elements and state
		TPair *out = new TPair[treeLevels];
		bool *half = new bool[treeLevels];
		for(int i = 0; i < treeLevels; ++i)
			half[i] = false;

		// it works like a binary counter of half[n]...half[1]half[0]
		if (_First != _Last)
			do {
				out[0].i1 = *_First; ++_First;
				out[0].i2 = *_FirstSA; ++_FirstSA;
				*(level[0]) = out[0];
				++(level[0]);				
				for(int i = 1; i < treeLevels; ++i) {
					if (half[i]) {
						if (out[i].i1 > out[i-1].i1)
							out[i].i1 = out[i-1].i1;
						*(level[i]) = out[i];	// out[i] is the minimum of last 2 values in out[i-1]
						++(level[i]);
						half[i] = false;
					} else {
						out[i].i1 = out[i-1].i1;
						out[i].i2 = *_FirstSA;
						half[i] = true;
						break;
					}
				}
			} while (_First != _Last);

		// complete half filled nodes
		bool carry = false;
		for(int i = 1; i < treeLevels; ++i)
			if (half[i] || carry) {
				if (half[i]) {
					if (out[i].i1 > out[i-1].i1)
						out[i].i1 = out[i-1].i1;
				} else
					out[i] = out[i-1];
				*(level[i]) = out[i];
				++(level[i]);
				carry = true;
			}

		// push trailing zero and the last SA table entry
		out[0].i1 = 0;
		out[0].i2 = *_FirstSA; ++_FirstSA;
		*_Dest = out[0];
		++_Dest;

		delete[] half;
		delete[] out;
		delete[] level;

		return _Dest;
    }


	template < typename TSize >
	inline TSize sizeofLCPE(TSize n)
	{
		if (n < 2) return n;	// 0 -> 0, 1 -> 1, 2 -> 2, 3 -> 4
		--n;
		TSize size = 2;
		for(TSize _xSize = 1; _xSize < n; _xSize *= 2)
			size += (n + _xSize - 1) / _xSize;
		return size;
	}

	template < typename TSize >
	inline TSize sizeofLCPH(TSize n)
	{
		return sizeofLCPE(n);
	}

	template <
		class LCPFwdIt,		// lcp table input iterator
		typename TSize >
	inline void sizeofLCPE(LCPFwdIt _First, LCPFwdIt _Last, TSize &_Size)
	{
		_Size = sizeofLCPE(distance(_First, _Last));
	}

	template <
		class LCPFwdIt,		// lcp table input iterator
		typename TSize >
	inline void sizeofLCPH(LCPFwdIt _First, LCPFwdIt _Last, TSize &_Size)
	{
		sizeofLCPE(_First, _Last, _Size);
		return;
	}


    template < typename TLCPE, typename TLCP >
    inline void createLCPBinTree(TLCPE &lcp_enhanced, TLCP &lcp) {
        createLCPBinTree(begin(lcp), end(lcp), begin(lcp_enhanced));
    }

    template < typename TLCPH, typename TLCP, typename TSA >
    inline void createHybridBinTree(TLCPH &lcp_hybrid, TLCP &lcp, TSA &sa) {
        createHybridBinTree(begin(lcp), end(lcp), begin(sa), begin(lcp_hybrid));
    }


	template < typename TSize >
    inline unsigned _treeLevels(TSize lcpSize)
	{
		unsigned treeLevels = 1;
		--lcpSize;
		TSize _xSize = 1;
		for(; lcpSize > _xSize; _xSize *= 2, ++treeLevels);
        return treeLevels;
    }

    template < typename TValue, typename TConfig, typename TLCP >
    inline void createLCPBinTree(String<TValue, External<TConfig> > &lcp_enhanced, TLCP &lcp) {
        int writeHeads = _treeLevels(length(lcp)) + 1;   // plus 1 write back buffer
        if (lcp_enhanced.cache.size() < writeHeads)
            lcp_enhanced.resizeCache(writeHeads);
        createLCPBinTree(begin(lcp), end(lcp), begin(lcp_enhanced));
    }

    template < typename TValue, typename TConfig, typename TLCP, typename TSA >
    inline void createHybridBinTree(String<TValue, External<TConfig> > &lcp_hybrid, TLCP &lcp, TSA &sa) {
        int writeHeads = _treeLevels(length(lcp)) + 1;   // plus 1 write back buffer
        if (lcp_hybrid.cache.size() < writeHeads)
            lcp_hybrid.resizeCache(writeHeads);
        createHybridBinTree(begin(lcp), end(lcp), begin(sa), begin(lcp_hybrid));
    }

}

#endif
