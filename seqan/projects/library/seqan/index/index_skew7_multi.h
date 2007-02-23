/*
 *  index_skew7_multi.h
 *  SeqAn
 *
 *  Created by David Weese on 27.01.07.
 *
 */

#ifndef SEQAN_HEADER_INDEX_SKEW7_MULTI_H
#define SEQAN_HEADER_INDEX_SKEW7_MULTI_H

namespace SEQAN_NAMESPACE_MAIN
{

    //////////////////////////////////////////////////////////////////////////////
    // external Skew7 algorithm (optimized for multiple sequences)
	//  - creates a suffix array of pairs (seqNo,seqPos)
    //////////////////////////////////////////////////////////////////////////////


	template < typename TInput, typename TPair, typename TLimitsString >
    struct Value< Pipe< TInput, Multi<Skew7, TPair, TLimitsString> > > {
        typedef TPair Type;
    };


    // *** COMPARATORS & MAPS ***
        
    template <typename InType, typename Result = int>
    struct skew7_ncomp_multi : public ::std::binary_function<InType,InType,Result> {
        inline Result operator()(const InType &a, const InType &b) const
        {
			typedef typename InType::T1 LocalPos;
            typedef typename InType::T2 Septet;
			typedef typename LocalPos::T2 SizeType;
            const typename Septet::T *sa = a.i2.i;
            const typename Septet::T *sb = b.i2.i;

            SizeType n = Septet::size;
			SizeType aLeft = getValueI2(a.i1);
			SizeType bLeft = getValueI2(b.i1);

            if (aLeft < n) n = aLeft;
            if (bLeft < n) n = bLeft;
            for(SizeType i = 0; i < n; i++, ++sa, ++sb) {
                if (*sa == *sb) continue;
                return (*sa < *sb)? -1 : 1;
            }
            if (n < Septet::size) {
				if (aLeft < bLeft) return -1;
				if (aLeft > bLeft) return 1;
				return (getValueI1(a.i1) > getValueI1(b.i1))? -1 : 1;
            }
            return 0;
        }
    };

    // optimized for bitvectors
    template <typename T1, typename TValue, typename Result>
    struct skew7_ncomp_multi< Pair<T1, Tuple<TValue, 7, Compressed>, Compressed >, Result > :
        public ::std::binary_function<
            Pair<T1, Tuple<TValue, 7, Compressed>, Compressed >,
            Pair<T1, Tuple<TValue, 7, Compressed>, Compressed >,
            Result> {       
        inline Result operator()(
            const Pair<T1, Tuple<TValue, 7, Compressed>, Compressed > &a,
            const Pair<T1, Tuple<TValue, 7, Compressed>, Compressed > &b) const
        {
			typedef typename T1::T2 SizeType;
			SizeType aLeft = getValueI2(a.i1);
			SizeType bLeft = getValueI2(b.i1);

			if (aLeft >= 7 && bLeft >= 7) {
				if (a.i2.i < b.i2.i) return -1;
				if (a.i2.i > b.i2.i) return 1;
				return 0;
			}

            SizeType n = 7;

			if (aLeft < n) n = aLeft;
            if (bLeft < n) n = bLeft;
            for(SizeType i = 0; i < n; i++) {
				if (a.i2[i] == b.i2[i]) continue;
				return (a.i2[i] < b.i2[i])? -1 : 1;
			}
			if (aLeft < bLeft) return -1;
			if (aLeft > bLeft) return 1;
			return (getValueI1(a.i1) > getValueI1(b.i1))? -1 : 1;
        }
    };

	template <
		typename InType, typename TLimitsString, typename TResultSize = typename Value<TLimitsString>::Type,
		typename Result = Pair<TResultSize, typename Value<InType,2>::Type, typename Spec<InType>::Type> >
    struct skew7_global_sliced_multi :
		public ::std::unary_function<InType,Result> {

		typedef TResultSize TSize;

		TSize			n4, n2, n1;
		TSize			n24;
		TLimitsString	off[5];

		skew7_global_sliced_multi(TLimitsString const &limits)
		{
			typename Iterator<TLimitsString const>::Type it = begin(limits), itEnd = end(limits);

			n4 = n2 = n1 = n24 = 0;

			if (it == itEnd) return;

			// count the numbers of septets in residue class 1, 2, and 4
			
			TSize size, cur;
			TSize old = *it; ++it;

			while (it != itEnd) {
				cur = *it;
				size = cur - old;
				old = cur;
				
				n4 += (size + 3) / 7;
				n2 += (size + 5) / 7;
				n1 += (size + 6) / 7;

				++it;
			}

			n24 = n2 + n4;

			// precompute the begin positions (off) of septet names 
			// in the sliced string for every sequence and every residue class

			resize(off[1], length(limits) - 1);
			resize(off[2], length(limits) - 1);
			resize(off[4], length(limits) - 1);

			it = begin(limits);
			old = *it; ++it;

			TSize seqNo = 0;
			TSize off4 = 0, off2 = n4, off1 = n24;
			while (it != itEnd) {
				cur = *it;
				size = cur - old;
				old = cur;
				
				off4 += (size + 3) / 7;
				off2 += (size + 5) / 7;
				off1 += (size + 6) / 7;

				off[4][seqNo] = (off4)? off4-1 : 0;
				off[2][seqNo] = (off2)? off2-1 : 0;
				off[1][seqNo] = (off1)? off1-1 : 0;

				++it; ++seqNo;
			}
		}

		inline Result operator()(const InType &x) const {
			TSize seqOfs = getValueI2(x.i1);
			return Result(off[seqOfs % 7][getValueI1(x.i1)] - seqOfs/7, x.i2);
		}
    };



    //////////////////////////////////////////////////////////////////////////////
    // skew7 class
    template < typename TInput, typename TPair, typename TLimitsString >
    struct Pipe< TInput, Multi<Skew7, TPair, TLimitsString> >
    {
  
        // *** SPECIALIZATION ***

        // use compression if lessorequal 16 different values per char
        typedef typename IF< 
            (BitsPerValue<_TypeOf(TInput)>::VALUE > 0) && 
            (BitsPerValue<_TypeOf(TInput)>::VALUE <= 4), 
            Compressed, 
            void>::Type compress;

        // use skew7 for recursion (more I/O-efficient)
        typedef Skew7 recurseSpec;

        // step 1
		typedef Pipe< TInput, Multi<Sampler<7, compress>, TPair, TLimitsString> >  TSamplerDC7;
                                        typedef skew7_ncomp_multi<_TypeOf(TSamplerDC7)> ncomp_t;
        typedef Pool< _TypeOf(TSamplerDC7), SorterSpec< SorterConfigSize<ncomp_t, _TSizeOf(TSamplerDC7) > > > TSortTuples;
		typedef Pipe< TSortTuples, Namer<ncomp_t> > TNamer;
                                        typedef skew7_global_sliced_multi<_TypeOf(TNamer), TLimitsString, _TSizeOf(TNamer)> func_slice_t;

		typedef Pipe< TNamer, Filter<func_slice_t> > TSlicedPos;
        typedef Pool< _TypeOf(TSlicedPos), MapperSpec< MapperConfigSize< filterI1<_TypeOf(TSlicedPos)>, _TSizeOf(TSlicedPos) > > > TNames_Sliced;

        // unique names - shortcut
        typedef Pool< _TypeOf(TNames_Sliced), MapperSpec< MapperConfigSize< func_slice_t, _TSizeOf(TNames_Sliced) > > > TNames_Linear_Unique;

        // non-unique names - recursion
        typedef Pipe< TNames_Sliced, Filter< filterI2<_TypeOf(TNames_Sliced)> > > TFilter;
        typedef Pipe< TFilter, recurseSpec > TRecurse;
		typedef Pipe< TRecurse, Counter > TRenamer;
        typedef Pool< _TypeOf(TRenamer), MapperSpec< MapperConfigSize< filterI1<_TypeOf(TRenamer)>, _TSizeOf(TRenamer) > > > TNames_Linear;

        // step 2
        typedef Pipe< Bundle2< TInput, TNames_Linear >, Extender7Multi<TPair, compress> > TExtender;
                                        typedef skew7_extend_comp<_TypeOf(typename TExtender::Out0),3> extend0_comp_t;
                                        typedef skew7_extend_comp<_TypeOf(typename TExtender::Out6),2> extend6_comp_t;
                                        typedef skew7_extend_comp<_TypeOf(typename TExtender::Out5),1> extend5_comp_t;
                                        typedef skew7_extend_comp<_TypeOf(typename TExtender::Out3),1> extend3_comp_t;
        typedef Pool< _TypeOf(typename TExtender::Out0), SorterSpec< SorterConfigSize< extend0_comp_t, _TSizeOf(typename TExtender::Out0) > > > TSorterS0;
        typedef Pool< _TypeOf(typename TExtender::Out6), SorterSpec< SorterConfigSize< extend6_comp_t, _TSizeOf(typename TExtender::Out6) > > > TSorterS6;
        typedef Pool< _TypeOf(typename TExtender::Out5), SorterSpec< SorterConfigSize< extend5_comp_t, _TSizeOf(typename TExtender::Out5) > > > TSorterS5;
        typedef Pool< _TypeOf(typename TExtender::Out3), SorterSpec< SorterConfigSize< extend3_comp_t, _TSizeOf(typename TExtender::Out3) > > > TSorterS3;

        // step 3
                                        typedef skew7_nmap_extended<_TypeOf(typename TExtender::Out124)> nmap_extended_t;
		typedef Pool< _TypeOf(typename TExtender::Out124), MapperSpec< MapperConfigSize< nmap_extended_t, _TSizeOf(typename TExtender::Out124) > > > TSorterS124;
        typedef Pipe< Bundle5< TSorterS0, TSorterS3, TSorterS5, TSorterS6, TSorterS124 >, Merger7Multi<TLimitsString> > TMerger;

        TInput				*textIn;
        TSorterS0			sortedS0;
        TSorterS3			sortedS3;
        TSorterS5			sortedS5;
        TSorterS6			sortedS6;
        TSorterS124			sortedS124;
        TMerger				in;
		TLimitsString const &limits;
            
        Pipe(TLimitsString const &_limits) :
            textIn(NULL),
			limits(_limits),
			in(bundle5(sortedS0, sortedS3, sortedS5, sortedS6, sortedS124), _limits) {}

		Pipe(TInput& _textIn, TLimitsString const &_limits) :
            textIn(&_textIn),
			limits(_limits),
			in(bundle5(sortedS0, sortedS3, sortedS5, sortedS6, sortedS124), _limits) {}
        
        inline void process(unsigned maxdepth = 0, unsigned depth = 1) {
            process(*textIn, maxdepth, depth);
        }

	    template < typename _TInput >
        bool process(_TInput &textIn, unsigned maxdepth = 0, unsigned depth = 1) {

            SEQAN_PROSET(PRODEPTH, depth);
            SEQAN_PROMARK("Rekursionsabstieg");
            #ifdef SEQAN_DEBUG_INDEX
                ::std::cout << "enter level " << depth << " compression: ";
                ::std::cout << TYPECMP<compress, Compressed>::VALUE << " "<<ValueSize<_TypeOf(TInput)>::VALUE<<"\n";
            #endif
            {


            // *** INSTANTIATION ***

            // step 1
            TSamplerDC7                 sampler(textIn, limits);
//		::std::cout << sampler;
            TSortTuples                 sorter;
            #ifdef SEQAN_DEBUG_INDEX
                ::std::cout << "  sort names (" << length(sampler)<< ")\n";
            #endif
            sorter << sampler;
            SEQAN_PROMARK("Sorter (2) - 7-lets sortieren");

			TNamer                      namer(sorter);
//		::std::cout << namer;
            func_slice_t				func_slice(limits);

			TSlicedPos					slicedPos(namer, func_slice);
            TNames_Sliced               names_sliced;
            #ifdef SEQAN_DEBUG_INDEX
                ::std::cout << "  slice names\n";
            #endif

            names_sliced << slicedPos;
//		::std::cout << names_sliced;

			if (namer.unique() || empty(names_sliced)) {
                // unique names

                clear(sorter);
                SEQAN_PROMARK("Mapper (4) - s124 konstruieren");

				TNames_Linear               names_S1, names_S2, names_S4;

                #ifdef SEQAN_DEBUG_INDEX
                    ::std::cout << "  make names linear\n";
                #endif

				skew7_separate_slices(
					names_sliced, func_slice, 
					names_S1, names_S2, names_S4);

				clear(names_sliced);
                SEQAN_PROMARK("Mapper (10) - ISA124 konstruieren");

                // step 2
                skew7_extend_multi(
					textIn, limits, 
					names_S1, names_S2, names_S4, 
					sortedS0, sortedS3, sortedS5, sortedS6, sortedS124);

            } else {
                // non-unique names

                clear(sorter);
                SEQAN_PROMARK("Mapper (4) - s124 konstruieren");

                TFilter                     filter(names_sliced);
				TRecurse                    recurse(filter);
                recurse.process(maxdepth, depth + 1);					// recursion

                #ifdef SEQAN_TEST_SKEW7
                {
                    String<typename Value<TFilter>::Type, Alloc<> > _text;
                    _text << filter;
                    SEQAN_DO(isSuffixArray(recurse, _text));
                }
                #endif

				clear(filter);
				TRenamer                    renamer(recurse);

				// partition SA by residue classes

				TNames_Linear               names_S1, names_S2, names_S4;

                #ifdef SEQAN_DEBUG_INDEX
                    ::std::cout << "  rename names\n";
                #endif

				skew7_separate_slices(
					renamer, func_slice, 
					names_S1, names_S2, names_S4);

				clear(renamer);
                SEQAN_PROMARK("Mapper (10) - ISA124 konstruieren");
               
                // step 2
                #ifdef SEQAN_DEBUG_INDEX
                    ::std::cout << "  prepare merge\n";
                #endif
                skew7_extend_multi(
					textIn, limits, 
					names_S1, names_S2, names_S4, 
					sortedS0, sortedS3, sortedS5, sortedS6, sortedS124);

                SEQAN_PROMARK("Mapper (12), Sorter (13-16) - SA124, SA3, SA5, SA6, SA0 verschmelzen");
            }
 
            // step 3
            // ... is done on-demand by merger
            }
            #ifdef SEQAN_DEBUG_INDEX
                ::std::cout << "left level " << depth << ::std::endl;
            #endif
            SEQAN_PROMARK("Rekursionsaufstieg");
            SEQAN_PROSUB(PRODEPTH, 1);

            return true;
        }

        inline typename Value<Pipe>::Type const operator*() {
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
	template < typename TInput, typename TObject, typename TPair, typename TLimitsString >
    inline bool operator<<(Pipe< TInput, Multi<Skew7, TPair, TLimitsString> > &me, TObject &textIn) {
        return me.process(textIn);
    }

	template < 
		typename TSA, 
		typename TValue, 
		typename TConfig,
		typename TSpec >
	inline void createSuffixArray(
		TSA &SA,
		StringSet< String<TValue, TConfig/*External<TConfig> */>, TSpec > &s,
		Skew7 const &spec,
		unsigned K = ValueSize<TValue>::VALUE,
        unsigned maxdepth = 0)
	{
        createSuffixArrayExt(SA, s, spec);
	}



}

#endif
