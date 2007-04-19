/*
 *  bwt.h
 *  SeqAn
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_INDEX_BWT_H
#define SEQAN_HEADER_INDEX_BWT_H

namespace SEQAN_NAMESPACE_MAIN
{
	
//namespace SEQAN_NAMESPACE_PIPELINING
//{

	struct BWT {};


    //////////////////////////////////////////////////////////////////////////////
    // external BWT algorithm
    //////////////////////////////////////////////////////////////////////////////


    template < typename TTextInput, typename TSuffixArrayInput >
    struct Value< Pipe< Bundle2< TTextInput, TSuffixArrayInput >, BWT > > {
        typedef typename Value<TTextInput>::Type Type;
    };

	//////////////////////////////////////////////////////////////////////////////
    // bwt class
    template < typename TTextInput, typename TSuffixArrayInput >
    struct Pipe< Bundle2< TTextInput, TSuffixArrayInput >, BWT >
    {
        // *** SPECIALIZATION ***

        typedef Pipe< TSuffixArrayInput, Counter > TSA;
		                                typedef typename Size<TTextInput>::Type	TSize;
		typedef Pool< _TypeOf(TSA), MapperSpec< MapperConfigSize< filterI1<_TypeOf(TSA)>, TSize> > > TInverter;
        typedef Pipe< TInverter, Filter< filterI2<_TypeOf(TInverter)> > > TCounterFilter;
		typedef Pipe< TTextInput, Shifter< -1, false > > TShiftText;

		typedef Pipe< Bundle2< TCounterFilter, TShiftText >, Joiner > TJoiner;
		typedef Pool< _TypeOf(TJoiner), MapperSpec< MapperConfigSize< filterI1<_TypeOf(TJoiner)>, TSize> > > TLinearMapper;
        typedef Pipe< TLinearMapper, Filter< filterI2<_TypeOf(TLinearMapper)> > > TFilter;

		TTextInput			*textIn;
		TSuffixArrayInput	*suffixArrayIn;
        TLinearMapper		mapper;
		TFilter				in;
        
        Pipe():
            in(mapper) {}

        Pipe(Bundle2< TTextInput, TSuffixArrayInput > const &_bundleIn):
            textIn(&_bundleIn.in1),
			suffixArrayIn(&_bundleIn.in2),
            in(mapper) {}

        inline void process() {
            process(*textIn, *suffixArrayIn);
        }

		template < typename _TTextInput, typename _TSuffixArrayInput >
        bool process(_TTextInput &textIn, _TSuffixArrayInput &suffixArrayIn) {

            // *** INSTANTIATION ***
			
			TSA							sa(suffixArrayIn);
			TInverter					inverter;
			TCounterFilter				filter(inverter);
			
            #ifdef SEQAN_DEBUG_INDEX
				::std::cerr << "  invert suffix array" << ::std::endl;
            #endif
			inverter << sa;
			SEQAN_PROMARK("Suffix-Array invertiert");

			TShiftText					shifter(textIn);
			TJoiner						joiner(bundle2(filter, shifter));
			
            #ifdef SEQAN_DEBUG_INDEX
                ::std::cerr << "  de-invert suffix array" << ::std::endl;
            #endif
			mapper << joiner;
			SEQAN_PROMARK("Suffix-Array linearisiert");

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
    inline bool operator<<(Pipe< TInput, BWT > &me, Bundle2< _TTextInput, _TSuffixArrayInput > const &bundleIn) {
 	    return me.process(bundleIn.in1, bundleIn.in2);
    }



    //////////////////////////////////////////////////////////////////////////////
    // external BWT algorithm (optimized for multiple sequences)
    //////////////////////////////////////////////////////////////////////////////


    template < typename TTextInput, typename TSuffixArrayInput, typename TPair, typename TLimitsString >
    struct Value< Pipe< Bundle2< TTextInput, TSuffixArrayInput >, Multi<BWT, TPair, TLimitsString> > > {
        typedef typename Value<TTextInput>::Type Type;
    };

	template <typename InType, typename TLimitsString, typename Result = typename Value<TLimitsString>::Type>
	struct filter_globalizer : public ::std::unary_function<InType,Result> {
		TLimitsString const &limits;
		filter_globalizer(TLimitsString const &_limits) : limits(_limits) {}
        inline Result operator()(const InType& x) const
        {
			return posGlobalize(x, limits);
		}
    };


	//////////////////////////////////////////////////////////////////////////////
    // bwt class
    template < typename TTextInput, typename TSuffixArrayInput, typename TPair, typename TLimitsString >
    struct Pipe< Bundle2< TTextInput, TSuffixArrayInput >, Multi<BWT, TPair, TLimitsString> >
    {
        // *** SPECIALIZATION ***

										typedef filter_globalizer<_TypeOf(TSuffixArrayInput), TLimitsString, _TSizeOf(TSuffixArrayInput)> filter_globalizer_t;
		typedef Pipe< TSuffixArrayInput, Filter<filter_globalizer_t> > TGlobalizer;
        typedef Pipe< TGlobalizer, Counter > TSA;
		                                typedef typename Size<TTextInput>::Type	TSize;
		typedef Pool< _TypeOf(TSA), MapperSpec< MapperConfigSize< filterI1<_TypeOf(TSA)>, TSize> > > TInverter;
        typedef Pipe< TInverter, Filter< filterI2<_TypeOf(TInverter)> > > TCounterFilter;
		typedef Pipe< TTextInput, Shifter< -1, false > > TShiftText;

		typedef Pipe< Bundle2< TCounterFilter, TShiftText >, Joiner > TJoiner;
		typedef Pool< _TypeOf(TJoiner), MapperSpec< MapperConfigSize< filterI1<_TypeOf(TJoiner)>, TSize> > > TLinearMapper;
        typedef Pipe< TLinearMapper, Filter< filterI2<_TypeOf(TLinearMapper)> > > TFilter;

		TTextInput			*textIn;
		TSuffixArrayInput	*suffixArrayIn;
        TLinearMapper		mapper;
		TFilter				in;

		TLimitsString const	&limits;
        
        Pipe(TLimitsString const &_limits):
            in(mapper),
			limits(_limits)	{}

        Pipe(Bundle2< TTextInput, TSuffixArrayInput > const &_bundleIn, TLimitsString const &_limits):
            textIn(&_bundleIn.in1),
			suffixArrayIn(&_bundleIn.in2),
            in(mapper),
			limits(_limits)	{}

        inline void process() {
            process(*textIn, *suffixArrayIn);
        }

		template < typename _TTextInput, typename _TSuffixArrayInput >
        bool process(_TTextInput &textIn, _TSuffixArrayInput &suffixArrayIn) {

            // *** INSTANTIATION ***
			
			TGlobalizer					globalizer(suffixArrayIn, limits);
			TSA							sa(globalizer);
			TInverter					inverter;
			TCounterFilter				filter(inverter);
			
            #ifdef SEQAN_DEBUG_INDEX
                ::std::cerr << "  invert suffix array" << ::std::endl;
            #endif
			inverter << sa;
			SEQAN_PROMARK("Suffix-Array invertiert");

			TShiftText					shifter(textIn);
			TJoiner						joiner(bundle2(filter, shifter));
			
            #ifdef SEQAN_DEBUG_INDEX
                ::std::cerr << "  de-invert suffix array" << ::std::endl;
            #endif
			mapper << joiner;
			SEQAN_PROMARK("Suffix-Array linearisiert");

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
    // you can call "bwt << pipe" or "bwt_t bwt(pipe); bwt.process()"
    // for the first we would need no _in member
	template < typename TInput, typename _TTextInput, typename _TSuffixArrayInput, typename TPair, typename TLimitsString >
    inline bool operator<<(Pipe< TInput, Multi<BWT, TPair, TLimitsString> > &me, Bundle2< _TTextInput, _TSuffixArrayInput > const &bundleIn) {
 	    return me.process(bundleIn.in1, bundleIn.in2);
    }



    //////////////////////////////////////////////////////////////////////////////
    // internal BWT algorithm
    //////////////////////////////////////////////////////////////////////////////


    template < typename TBWT,
               typename TText,
               typename TSA >
    void createBWTableInt(
		TBWT &bwt,
		TText &s,
		TSA const &SA)
	{
		typedef typename Value<TSA>::Type	TValue;
		typedef typename Size<TSA>::Type	TSize;

		#ifdef SEQAN_DEBUG_INDEX
			if (sizeof(TSize) > 4)
				::std::cerr << "WARNING: TSize size is greater 4 (BWT)" << ::std::endl;
        #endif

		TSize n = length(s);

		for(TSize i = 0; i < n; ++i) {
			TValue sa = SA[i];
			if (sa)
				bwt[i] = s[sa - 1];
			else
				bwt[i] = TSize();
		}
	}

//}

}

#endif
