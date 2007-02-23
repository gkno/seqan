/*
 *  test_index_creation.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_TEST_INDEX_CREATION_H
#define SEQAN_HEADER_TEST_INDEX_CREATION_H

namespace SEQAN_NAMESPACE_MAIN
{

	template < typename TBuffer >
	void permute(TBuffer &buf) {
		typename Size<TBuffer>::Type i, j, s = length(buf);
//        srand( (unsigned)time( NULL ) );
		for(i = 0; i < s; i++)
			buf[i] = s-i-1;
		for(i = 0; i < s; i++) {
            if (i > 0) {
			    j = i - (rand() % i) - 1;
                assert(0 <= j && j < s);
            } else
                j = 0;
			unsigned tmp = buf[i];
			buf[i] = buf[j];
			buf[j] = tmp;
		}
	}

	template < typename TBuffer >
	void blank(TBuffer &buf) {
		typename Size<TBuffer>::Type i, s = length(buf);
		typename Value<TBuffer>::Type c = typename Value<TBuffer>::Type();
		for(i = 0; i < s; i++)
            buf[i] = c;
	}

	template < typename TBuffer >
	void randomize(TBuffer &buf) {
		typename Size<TBuffer>::Type i, s = length(buf);
		for(i = 0; i < s; i++)
            buf[i] = rand() % s;
	}

	template < typename TBuffer >
	void textRandomize(TBuffer &buf) {
		typename Size<TBuffer>::Type i, s = length(buf);
		for(i = 0; i < s; i++)
			buf[i] = '@' + rand() % 2;//('z'-'@');
	}

	template < typename TValue >
	struct IdentityMap : public ::std::unary_function< TValue, TValue > {
		inline TValue operator() (TValue const i) { return i; }
	};

	template < typename TValue >
	struct SimpleCompare : public ::std::binary_function< TValue const, TValue const, int > {
		inline int operator() (TValue const a, TValue const b) const {
            if (a < b) return -1;
            if (a > b) return 1;
            return 0;
        }
	};
	
    template <typename TSequence>
    bool isPermutation(TSequence const &SA) {
        typedef typename Value<TSequence>::Type TSize;
        TSize n = length(SA);
        bool *seen = new bool[n];
        TSize i;
        for (i = 0;  i < n;  i++) seen[i] = 0;
        for (i = 0;  i < n;  i++)
            if (SA[i] >= 0 && SA[i] < n)
                seen[SA[i]] = 1;
            else
                printf("isPermutation: SA index out of range (n=%d) SA[%d]=%d\n", n, i, SA[i]);

        for (i = 0;  i < n;  i++)
            if (!seen[i]) {
                delete[] seen;
                return 0;
            }
        delete[] seen;
        return 1;
    }

    template <typename TInput, typename TSpec>
    bool isPermutation(Pipe<TInput, TSpec> &SA) {
        typedef typename Value< Pipe<TInput, TSpec> >::Type TSize;
        TSize n = length(SA);
        bool *seen = new bool[n];
        TSize i;
        for (i = 0;  i < n;  i++) seen[i] = 0;
        beginRead(SA);
        for (i = 0;  i < n;  i++, ++SA)
            if (*SA >= 0 && *SA < n) {
                if (seen[*SA])
                    printf("isPermutation: not unique %d->%d\n", i, *SA);
                seen[*SA] = true;
            } else
                printf("isPermutation: SA index out of range (n=%d) SA[%d]=%d\n", n, i, *SA);
        endRead(SA);

        for (i = 0;  i < n;  i++)
            if (!seen[i]) {
                delete[] seen;
                return false;
            }
        delete[] seen;
        return true;
    }



    template <typename TIt, typename ST>
    bool __sleq(TIt s1, TIt s2, ST n1, ST n2) {
        ST n = min(n1, n2);
        for(ST i = 0; i < n; i++, ++s1, ++s2) {
            if (lexLess(*s1,*s2)) return 1;
            if (lexLess(*s2,*s1)) {
                ::std::cout<<(lexLess(*s2,*s1));
                ::std::cout<<*s1;
                ::std::cout<<*s2;
                ::std::cout << "after " << i << " compares not " << (unsigned)*s1 << " leq " << (unsigned)*s2 << ::std::endl;
                return 0;
            }
        }
        return (n1 < n2);
    } 

    // is SA a sorted suffix array for s?
    template <typename TSequence, typename TText>
    bool isSorted(TSequence &SA, TText const &s) {
        typedef typename Value<TSequence>::Type TSize;
        TSize n = length(s);
	    for(TSize i = 1; i < n; ++i) {
		    if (!__sleq(begin(s) + SA[i-1], begin(s) + SA[i], n-SA[i-1], n-SA[i])) {
			    printf("isSorted: sort error s_%d(SA[%d]) >= s_%d(SA[%d])\n",SA[i-1],i-1,SA[i],i);

			String<unsigned, External<> > safile;
			if (!open(safile,"error.sa")) printf("could not open ERROR.SA\n");
			safile = SA;
			    return false;
		    }
	    }
	    return true;  
    }

    template <typename TInput, typename TSpec, typename TText>
    bool isSorted(Pipe<TInput, TSpec> &SA, TText const &s) {
        typedef typename Value< Pipe<TInput, TSpec> >::Type TSize;
        TSize n = length(s);
	    beginRead(SA);
        TSize prev = *SA; ++SA;
	    for(TSize i = 1; !eof(SA); ++SA, ++i) {
		    if (!__sleq(begin(s) + prev, begin(s) + *SA, n-prev, n-*SA)) {
                printf("isSorted: sort error s_%d(SA[%d]) >= s_%d(SA[%d])\n",prev,i-1,*SA,i);
				endRead(SA);

			String<unsigned, External<> > safile;
			if (!open(safile,"error.sa")) printf("could not open ERROR.SA\n");
			safile << SA;

			    return false;
		    }
		    prev = *SA;
	    }
        endRead(SA);
	    return true;  
    }



    template <typename TIt, typename ST>
    bool __sleqLCP(TIt s1, TIt s2, ST n1, ST n2, ST lcp) {
        ST n = min(n1, n2);
        for(ST i = 0; i < n; i++, ++s1, ++s2) {
            if (lexLess(*s1,*s2)) return (i == lcp);
            if (lexLess(*s2,*s1)) {
                ::std::cout << "after " << i << " compares not " << (unsigned)*s1 << " leq " << (unsigned)*s2 << ::std::endl;
                return false;
            }
        }
        return (n1 < n2) && (n == lcp);
    } 

    // is SA a sorted suffix array and LCP the correct LCP-Table for s?
/*    template <
		typename TSize1, typename TSpec1,
		typename TSize2, typename TSpec2,
		typename TText >
    bool isSortedLCP(String<TSize1, TSpec1> &LCP, String<TSize2, TSpec2> &SA, TText const &s) {
        TSize2 n = length(s);
	    for(TSize2 i = 1; i < n; ++i) {
		    if (!__sleq(begin(s) + SA[i-1], begin(s) + SA[i], n-SA[i-1], n-SA[i], LCP[i-1])) {
			    printf("isSorted: sort error s_%d(%d) >= s_%d(%d)\n",i-1,SA[i-1],i,SA[i]);
			    return false;
		    }
	    }
	    return true;  
    }
*/
    template < typename TLCP, typename TSA, typename TText >
    bool isSortedLCP(TLCP &LCP, TSA &SA, TText const &s) {
        typedef typename Value<TSA>::Type TSize;
        typedef typename Iterator<TSA>::Type ISA;
        typedef typename Iterator<TLCP>::Type ILCP;

        TSize n  = length(s);
        ISA  sa  = begin(SA);
        ILCP lcp = begin(LCP);

	    TSize prev = *sa; ++sa;
	    for(TSize i = 1; sa != end(SA); ++sa, ++lcp, ++i) {
		    if (!__sleqLCP(begin(s) + prev, begin(s) + *sa, n-prev, n-*sa, *lcp)) {
                printf("isLCP: sort error s_%d(%d) >= s_%d(%d)\n",i-1,prev,i,*sa);
			    return false;
		    }
		    prev = *sa;
	    }
	    return true;  
    }

    template <typename TSufArray, typename TText>
    bool isSuffixArray(TSufArray &SA, TText const &s) {
        if (length(SA) != length(s)) {
            printf("isSuffixArray: length is bad: SA=%d, s=%d\n", length(SA), length(s));
            return false;
        }
        
        if (!isPermutation(SA)) {
            ::std::cout<<"isSuffixArray: SA is not a permutation!\n";
            return false;
        }

        if (!isSorted(SA, s)) {
            ::std::cout<<"isSuffixArray: SA is not sorted!\n";
			String<unsigned char, External<> > textfile;
			if (!open(textfile,"error.txt")) printf("could not open ERROR.TXT\n");
			textfile=s;
			return false;
        }

//        ::std::cout<<"SATest OK! n="<<length(s)<<std::endl;
        return true;
    }

    template <typename TLCP, typename TSufArray, typename TText>
    bool isLCPTable(TLCP &LCP, TSufArray &SA, TText const &s) {
        if (length(SA) != length(s)) {
            printf("isLCPTable: length is bad: SA=%d, s=%d\n", length(SA), length(s));
            return false;
        }
        
        if (length(LCP) != length(s)) {
            printf("isLCPTable: length is bad: LCP=%d, s=%d\n", length(LCP), length(s));
            return false;
        }
        
        if (!isPermutation(SA)) {
            ::std::cout<<"isLCPTable: SA is not a permutation!\n";
            return false;
        }

        if (!isSortedLCP(LCP, SA, s)) {
            ::std::cout<<"isLCPTable: SA is not sorted!\n";
            return false;
        }

//        ::std::cout<<"LCPTest OK! n="<<length(s)<<std::endl;
        return true;
    }

    template <typename TA, typename TB>
    bool isEqual(TA &_a, TB &_b) {
        typedef typename Iterator<TA>::Type IA;
        typedef typename Iterator<TB>::Type IB;

        IA a = begin(_a), e = end(_a);
        IB b = begin(_b);
        while (a!=e) {
            if (!(*a == *b)) {
                ::std::cout << "isEqual: difference at " << (e-a) << " a=" << *a << "  b=" << *b << ::std::endl;
                return false;
            }
            ++a;
            ++b;
        }

//        ::std::cout<<"EQUALTest OK! n="<<length(_a)<<std::endl;
        return true;
    }

}

#endif

