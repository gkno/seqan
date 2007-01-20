/*
 *  lcp.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_INDEX_LCP_H
#define SEQAN_HEADER_INDEX_LCP_H

namespace SEQAN_NAMESPACE_MAIN
{
	
//namespace SEQAN_NAMESPACE_PIPELINING
//{

	struct Kasai {};
	struct KasaiInPlace {};


    //////////////////////////////////////////////////////////////////////////////
    // external LCP algorithm (modified Kasai et al. for pipelining)
    //////////////////////////////////////////////////////////////////////////////

    template < typename TTextInput, typename TSuffixArrayInput >
    struct Value< Pipe< Bundle2< TTextInput, TSuffixArrayInput >, Kasai > > {
        typedef typename Size<TTextInput>::Type Type;
    };

	template <typename InType, typename Result = typename InType::T2::T>
	struct map_inverse : public ::std::unary_function<InType,Result> {
        Result operator()(const InType& x) const
        { return x.i2[0]; }
    };

    //////////////////////////////////////////////////////////////////////////////
    // lcp class
    template < typename TTextInput, typename TSuffixArrayInput >
    struct Pipe< Bundle2< TTextInput, TSuffixArrayInput >, Kasai >
    {
        // *** SPECIALIZATION ***

        typedef Pipe< TSuffixArrayInput, Echoer<2,false> > TEchoer;
                                        typedef map_inverse<_TypeOf(TEchoer)> map_inverse_t;
		                                typedef typename Size<TTextInput>::Type	TSize;
		typedef Pool< _TypeOf(TEchoer), MapperSpec< MapperConfigSize< map_inverse_t, TSize> > > TInverter;
		                                typedef Pair<TSize,	TSize, Compressed> TCoreType;
		typedef Pool< TCoreType, MapperSpec< MapperConfigSize< getI1<TCoreType>, TSize > > > TLinearMapper;
        typedef Pipe< TLinearMapper, Filter< getI2<TCoreType> > > TFilter;

		TTextInput				*textIn;
        TSuffixArrayInput		*suffixArrayIn;
        TLinearMapper           mapper;
		TFilter					in;
        const LcpConfig			conf;
        
        Pipe():
            in(mapper) {}

        Pipe(const LcpConfig &_conf):
            in(mapper),
			conf(_conf) {}

        Pipe(Bundle2< TTextInput, TSuffixArrayInput > const &_bundleIn):
            textIn(&_bundleIn.in1),
			suffixArrayIn(&_bundleIn.in2),
            in(mapper) {}

        Pipe(Bundle2< TTextInput, TSuffixArrayInput > const &_bundleIn, const LcpConfig &_conf):
            textIn(&_bundleIn.in1),
			suffixArrayIn(&_bundleIn.in2),
            in(mapper),
			conf(_conf) {}
        
        inline void process() {
            process(*textIn, *suffixArrayIn);
        }

		template < typename _TTextInput, typename _TSuffixArrayInput >
        bool process(_TTextInput &textIn, _TSuffixArrayInput &suffixArrayIn) {

            // *** INSTANTIATION ***
			
			TEchoer						echoer(suffixArrayIn);
			TInverter					inverter(echoer);

            #ifdef SEQAN_DEBUG_INDEX
                ::std::cout << "  invert suffix array\n";
            #endif
			inverter << echoer;
			SEQAN_PROMARK("Suffix-Array invertiert");

			lcp_process(textIn, inverter, mapper);
            return true;
        }

        inline typename Value<Pipe>::Type const operator*() const {
            return *in;
        }
        
        inline Pipe& operator++() {
            ++in;
            return *this;
        }
	};

    // not sure which interface is more intuitive, we support both
    // you can call "skew << pipe" or "skew_t skew(pipe); skew.process()"
    // for the first we would need no _in member
	template < typename TInput, typename _TTextInput, typename _TSuffixArrayInput >
    inline bool operator<<(Pipe< TInput, Kasai > &me, Bundle2< _TTextInput, _TSuffixArrayInput > const &bundleIn) {
 	    return me.process(bundleIn.in1, bundleIn.in2);
    }

/**
.Function.createLCPTable:
..summary:Creates a lcp table from a given text and suffix array.
..cat:Index
..signature:createLCPTable(lcp, text, suffixArray[, algo_tag])
..param.lcp:A reference to the resulting lcp table.
..param.text:A given text.
..param.suffixArray:The suffix array of $text$.
..param.algo_tag:A tag that identifies the algorithm which is used for creation.
..remarks:The size of $lcp$ must be at least $length(text)$ before calling this function.
*/


	template < 
        typename TLCPTable,
        typename TValue, 
		typename TConfig,
		typename TObject, 
		typename ConstrSpec>
	void createLCPTable(
		TLCPTable &LCP,
		TObject &s,
		String< TValue, External<TConfig> > &SA,
		Kasai const &spec)
	{
        createLCPTableExt(LCP, s, SA, spec);
	}


    //////////////////////////////////////////////////////////////////////////////
    // internal Kasai algorithm
    //////////////////////////////////////////////////////////////////////////////

    template < typename TLCPTable,
               typename TText,
               typename TSuffixArray >
    void createLCPTable(
		TLCPTable &LCP,
		TText const &s,
		TSuffixArray const &SA,
		Kasai const &)
	{
		typedef typename Value<TSuffixArray>::Type TSize;

		#ifdef SEQAN_DEBUG_INDEX
			if (sizeof(TSize) > 4)
				::std::cout << "WARNING: TSize size is greater 4 (Kasai)" << ::std::endl;
        #endif

		TSize n = length(s);
        if (n < 2) return;

        #ifdef SEQAN_DEBUG_INDEX
            TSize lcpMax = 0, lcpAvrg = 0, lcpNumer = 0, sigma = 1;	// for lcpMax, lcpMean, |Sigma|
        #endif

        String<TSize, Alloc<> > ISA;
		resize(ISA, n, Exact());

		for(TSize i = 0; i < n; ++i)
			ISA[SA[i]] = i;
		
		SEQAN_PROMARK("Suffix-Array invertiert");

		typename Iterator<TText const>::Type Ibegin = begin(s), Iend = end(s);
        typename Iterator<TText const>::Type I = Ibegin, J;
        for(TSize i = 0, h = 0, j, isa; i < n; ++i) {
			if (isa = ISA[i]) {
				J = Ibegin + h + (j = SA[isa - 1]);
                for(TSize hMax = Min(n - i, n - j); h < hMax && *I == *J; ++I, ++J, ++h);
				LCP[isa - 1] = h;
                #ifdef SEQAN_DEBUG_INDEX
                    if ((lcpNumer += h) > n) {
                        lcpNumer -= n;
                        ++lcpAvrg;
                    }
                    if (lcpMax < h) lcpMax = h;
					if (!h) ++sigma;
                #endif
			}
			if (h) --h;
            else ++I;
        }
		LCP[n - 1] = 0;
        #ifdef SEQAN_DEBUG_INDEX
            ::std::cout << "  n: " << n;
            ::std::cout << "  lcpMax: " << lcpMax;
            ::std::cout << "  lcpAvrg: " << (TSize)(lcpAvrg + (lcpNumer + n/2) / n);
            ::std::cout << "  sigma: " << sigma << std::endl;
        #endif
	}

	// HINT:
	// In contrast to the upper functions 
	// createLCPTableInPlace expects the lcp table to be of size n
    template < typename TLCPTable,
               typename TText,
               typename TSuffixArray >
    void createLCPTable(
		TLCPTable &LCP,
		TText const &s,
		TSuffixArray const &SA,
		KasaiInPlace const &)
	{
		typedef typename Value<TSuffixArray>::Type TSize;

		#ifdef SEQAN_DEBUG_INDEX
			if (sizeof(TSize) > 4)
				::std::cout << "WARNING: TSize size is greater 4 (Kasai)" << ::std::endl;
        #endif

		TSize n = length(s);
        if (n < 2) return;

        #ifdef SEQAN_DEBUG_INDEX
            TSize lcpMax = 0, lcpAvrg = 0, lcpNumer = 0, sigma = 1;	// for lcpMax, lcpMean, |Sigma|
        #endif

		TSize mark = ~(~0u>>1);
		TSize mask =   ~0u>>1;

		for(TSize i = 0; i < n; ++i)
			LCP[SA[i]] = i;
		
		SEQAN_PROMARK("Suffix-Array invertiert");
        #ifdef SEQAN_DEBUG_INDEX
			::std::cout << "Suffix-Array invertiert" << ::std::endl;
		#endif

		typename Iterator<TText const>::Type Ibegin = begin(s), Iend = end(s);
        typename Iterator<TText const>::Type I = Ibegin, J;
        for(TSize i = 0, h = 0, j, isa; i < n; ++i) {
			if ((isa = LCP[i] + 1) < n) {
				J = Ibegin + h + (j = SA[isa]);
                for(TSize hMax = Min(n - i, n - j); h < hMax && *I == *J; ++I, ++J, ++h);
				LCP[i] = h | mark;
                #ifdef SEQAN_DEBUG_INDEX
                    if ((lcpNumer += h) > n) {
                        lcpNumer -= n;
                        ++lcpAvrg;
                    }
                    if (lcpMax < h) lcpMax = h;
					if (!h) ++sigma;
                #endif
			}
			if (h) --h;
            else ++I;
        }
		LCP[SA[n - 1]] = mark;
		SEQAN_PROMARK("permutierte LCP-Tabelle erzeugt");
        #ifdef SEQAN_DEBUG_INDEX
			::std::cout << "permutierte LCP-Tabelle erzeugt" << ::std::endl;
		#endif
        for(TSize i = 0, j, tmp; i < n; ++i)
			if (LCP[i] & mark) {
				j = i;
				tmp = LCP[j];
				while (SA[j] != i) {
					LCP[j] = LCP[SA[j]] & mask;
					j = SA[j];
				}
				LCP[j] = tmp & mask;
			}
        #ifdef SEQAN_DEBUG_INDEX
			::std::cout << "LCP-Tabelle erzeugt" << ::std::endl;
		#endif

        #ifdef SEQAN_DEBUG_INDEX
            ::std::cout << "  n: " << n;
            ::std::cout << "  lcpMax: " << lcpMax;
            ::std::cout << "  lcpAvrg: " << (TSize)(lcpAvrg + (lcpNumer + n/2) / n);
            ::std::cout << "  sigma: " << sigma << std::endl;
        #endif
	}

//}

}

#endif
