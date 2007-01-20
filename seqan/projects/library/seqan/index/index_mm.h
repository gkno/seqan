/*
 *  mm.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_INDEX_MM_H
#define SEQAN_HEADER_INDEX_MM_H

namespace SEQAN_NAMESPACE_MAIN
{

	struct ManberMyers {};


    //////////////////////////////////////////////////////////////////////////////
    // internal Manber Myers algorithm
    //////////////////////////////////////////////////////////////////////////////

    template < typename TSuffixArray,
               typename TText >
    void createSuffixArray(
		TSuffixArray &SA,
		TText &s,
		ManberMyers const &,
		unsigned K = ValueSize< typename Value<TText>::Type >::VALUE,
        unsigned maxdepth = 0)
	{
		typedef typename Value<TSuffixArray>::Type TSize;
		typedef typename Value<TText>::Type TValue;

		TSize n = length(s);

        String<TSize, Alloc<> > ISA;        resize(ISA, n, Exact());
        String<TSize, Alloc<> > count;		resize(count, n, Exact());
        String<bool, Alloc<> > Bh;			resize(Bh, n, Exact());
        String<bool, Alloc<> > B2h;			resize(B2h, n, Exact());
	    
		for(TSize i = 0; i < n; ++i)
			ISA[i] = i;
        SEQAN_PROMARK("Suffix-Array invertiert");

		radixPass(SA, ISA, s, count, K);
		::std::memset(begin(Bh), 0, sizeof(bool) * n);

		TValue c = 0, d;
		for(TSize i = 0; i < n; ++i) {
			if (c != (d = s[SA[i]])) {
				c = d;
				Bh[i] = true;
			}
		}

		TSize j, k, l, *p, Ti, cd = 0;
		for(TSize h = 1; h < n; h<<= 1, ++cd);
        if (maxdepth > 0 && maxdepth < cd) cd = maxdepth;
		for(TSize h = 1; cd > 0; h<<= 1, --cd) {

            #ifdef SEQAN_DEBUG_INDEX
				::std::cout << "[" << cd << "] ";
            #endif
            SEQAN_PROADD(PRODEPTH, 1);
            SEQAN_PROMARK("Beginne Durchlauf");
			::std::memset(begin(count), 0, sizeof(TSize) * n);
			::std::memset(begin(B2h),   0, sizeof(bool)  * n);

			l = 0;
			for(TSize i = 0; i < n; ++i) {
				if (Bh[i]) l = i;
				ISA[SA[i]] = l;
			}

			Ti = n - h;
			p = begin(ISA) + Ti;
			j = count[*p]++;
			*p+= j;
			B2h[*p] = true;

			l = 0;
			for(TSize i = 0; i < n; ++i) {

				if ((Ti = SA[i]) >= h) {
					Ti-= h;
					p = begin(ISA) + Ti;
					j = count[*p]++;
					*p+= j;
					B2h[*p] = true;
				}

				if (i + 1 == n || Bh[i + 1]) {
					for(j = l; j <= i; ++j)
						if ((Ti = SA[j]) >= h) {
							Ti-= h;
							if (B2h[k = ISA[Ti]])
								while (++k < n && !Bh[k] && B2h[k])
									B2h[k] = false;
						}
					l = i + 1;
				}
			}

			for(TSize i = 0; i < n; ++i) {
				SA[ISA[i]] = i;
				Bh[i]|= B2h[i];
			}
		}
        SEQAN_PROSET(PRODEPTH, 0);
        #ifdef SEQAN_DEBUG_INDEX
			::std::cout << ::std::endl;
        #endif
	}

    // creates suffix array sorted by the first maxLCP chars of suffixes
    template < typename TSuffixArray,
               typename TText,
               typename TSize >
    inline void createSuffixArrayPart(
		TSuffixArray &SA,
		TText &s,
		ManberMyers const &_dummy,
        TSize maxLCP,
        unsigned K = ValueSize< typename Value<TText>::Type >::VALUE)
    {
        unsigned depth = 0;
        for(TSize i = 1; i < maxLCP; i*=2) ++depth;
        createSuffixArray(SA, s, _dummy, K, depth);
    }

}

#endif //#ifndef SEQAN_HEADER_...
