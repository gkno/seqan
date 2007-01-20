/*
 *  bwt.h
 *  genindex
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

	template <typename TValue, typename Result = TValue>
    struct _scanUndefined : public std::unary_function<TValue,TValue> {
		typedef typename TValue::T1 TSize;
		TSize &undef;

		_scanUndefined(TSize &_undef): undef(_undef) {}
        TValue const & operator()(TValue const &x) const {
			if (x.i1 == 0) undef = x.i2;
			return x; 
		}
    };


	//////////////////////////////////////////////////////////////////////////////
    // Enhanced class (outputs only the childtab (3rd) column of the Enhanced Suffix Array)
    template < typename TTextInput, typename TSuffixArrayInput >
    struct Pipe< Bundle2< TTextInput, TSuffixArrayInput >, BWT >
    {
        // *** SPECIALIZATION ***

        typedef Pipe< TSuffixArrayInput, Counter > TSA;
		typedef Pipe< TSA, Filter<_scanUndefined<_TypeOf(TSA)> > > TScanner;
		                                typedef typename Size<TTextInput>::Type	TSize;
		typedef Pool< _TypeOf(TScanner), MapperSpec< MapperConfigSize< getI1<_TypeOf(TScanner)>, TSize> > > TInverter;
        typedef Pipe< TInverter, Filter< getI2<_TypeOf(TInverter)> > > TCounterFilter;
		typedef Pipe< TTextInput, Shifter< -1, false > > TShiftText;

		typedef Pipe< Bundle2< TCounterFilter, TShiftText >, Joiner > TJoiner;
		typedef Pool< _TypeOf(TJoiner), MapperSpec< MapperConfigSize< getI1<_TypeOf(TJoiner)>, TSize> > > TLinearMapper;
        typedef Pipe< TLinearMapper, Filter< getI2<_TypeOf(TLinearMapper)> > > TFilter;

		TTextInput			*textIn;
		TSuffixArrayInput	*suffixArrayIn;
        TLinearMapper		mapper;
		TFilter				in;
		TSize				undefined;
        
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
			TScanner					scanner(sa, _scanUndefined<_TypeOf(TSA)>(undefined));
			TInverter					inverter;
			TCounterFilter				filter(inverter);
			
            #ifdef SEQAN_DEBUG_INDEX
                ::std::cout << "  invert suffix array\n";
            #endif
			inverter << sa;
			SEQAN_PROMARK("Suffix-Array invertiert");

			TShiftText					shifter(textIn);
			TJoiner						joiner(bundle2(filter, shifter));
			
            #ifdef SEQAN_DEBUG_INDEX
                ::std::cout << "  de-invert suffix array\n";
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


	template < typename TBWTStruct,
               typename TText,
			   typename TSA >
    void createBWTableExt(
		TBWTStruct &bwt,
		TText &s,
		TSA &SA)
	{
		// specialization
		typedef Pipe< TText, Source<> >						srcText_t;
		typedef Pipe< TSA, Source<> >   					srcSA_t;
	    typedef Pipe< Bundle2< srcText_t, srcSA_t >, BWT >	creator_t;

		// instantiation
		srcText_t	srcText(s);
		srcSA_t		srcSA(SA);
		creator_t	creator;

		// processing
	    creator << bundle2(srcText, srcSA);
		bwt.tab << creator;
		bwt.undefined = creator.undefined;
	}

/**
.Function.createBWTable:
..summary:Creates a Burrows-Wheeler table from a given text and suffix array.
..cat:Index
..signature:createBWTable(bwt, text, suffixArray[, algo_tag])
..param.bwt:A reference to the resulting Burrows-Wheeler table.
..param.text:A given text.
..param.suffixArray:The suffix array of $text$.
..param.algo_tag:A tag that identifies the algorithm which is used for creation.
..remarks:The size of $bwt$ must be at least $length(text)$ before calling this function.
*/

	template < typename TValue,
			   typename TConfig,
               typename TText,
			   typename TSA >
    inline void createBWTable(
		BWTStruct< String<TValue, External<TConfig> > > &childtab,
		TText &s,
		TSA &sa)
	{
		createBWTableExt(childtab, s, sa);
	}


    //////////////////////////////////////////////////////////////////////////////
    // internal BWT algorithm
    //////////////////////////////////////////////////////////////////////////////

    template < typename TBWTStruct,
               typename TText,
               typename TSuffixArray >
    void createBWTable(
		TBWTStruct &bwt,
		TText &s,
		TSuffixArray &SA)
	{
		typedef typename Value<TSuffixArray>::Type TSize;

		#ifdef SEQAN_DEBUG_INDEX
			if (sizeof(TSize) > 4)
				::std::cout << "WARNING: TSize size is greater 4 (BWT)" << ::std::endl;
        #endif

		TSize n = length(s);

		for(TSize i = 0, sa; i < n; ++i)
			if (sa = SA[i])
				bwt.tab[i] = s[sa - 1];
			else {
				bwt.tab[i] = TSize();
				bwt.undefined = i;
			}
	}

//}

}

#endif
