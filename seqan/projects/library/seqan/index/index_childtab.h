/*
 *  childtab.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_INDEX_CHILDTAB_H
#define SEQAN_HEADER_INDEX_CHILDTAB_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

	struct ChildTab {};

    //////////////////////////////////////////////////////////////////////////////
    // external childtab algorithm
    //////////////////////////////////////////////////////////////////////////////

    template < typename TLCPInput >
    struct Value< Pipe< TLCPInput, ChildTab > > {
        typedef typename Size<TLCPInput>::Type Type;
    };


	template < typename TLCPInput, typename TDest >
	inline void childtab_process(TLCPInput &lcpIn, TDest &dest)
	{
		typedef typename Value<TLCPInput>::Type				TValue;
		typedef typename Size<TLCPInput>::Type				TSize;

		typedef Pair<TSize, TValue, Compressed>				TPair;		// (i, lcptab[i])
		typedef ::std::stack<TPair> 						TStack;

		TStack stack_updown;
		TStack stack_nextl;

		stack_updown.push(TPair(0, 0));
		stack_nextl.push(TPair(0, 0));

		dest.undefinedValue = TPair(supremumValue<TSize>(), 0);		// undefined value for unused entries
		resize(dest, length(lcpIn));
		beginRead(lcpIn);
		beginWrite(dest);

		TSize n = length(lcpIn);
		TSize undef = n + 1;
		TSize lastIndex_nextl;
		TPair lastIndex_updown = TPair(undef, 0);

		TValue lcp_i;
		for(TSize i = 0; i < n; ++i) {
			lcp_i = *lcpIn; ++lcpIn;

			//////////////////////////////////////////////////////////////////////////////
			// nextlIndex part
			lastIndex_nextl = undef;
			while (lcp_i < stack_nextl.top().i2)
				stack_nextl.pop();

			if (lcp_i == stack_nextl.top().i2) {
				lastIndex_nextl = stack_nextl.top().i1;
				push(dest, TPair(lastIndex_nextl, i + 1));					// childtab[top].nextl

//				::std::cout << "nextl:\t\t\t" << TPair(lastIndex_nextl, i + 1) << ::std::endl;
				stack_nextl.pop();
			}
			stack_nextl.push(TPair(i + 1, lcp_i));


			//////////////////////////////////////////////////////////////////////////////
			// up/down part
			TPair top;
			while (lcp_i < stack_updown.top().i2) {
				lastIndex_updown = stack_updown.top();
				stack_updown.pop();
				
				top = stack_updown.top();
				if (lcp_i <= top.i2 && top.i2 != lastIndex_updown.i2 && top.i1 != lastIndex_nextl) {
					push(dest, TPair(top.i1, lastIndex_updown.i1));			// childtab[top].down

//					::std::cout << "down:\t\t" << TPair(top.i1, lastIndex_updown.i1) << ::std::endl;			// childtab[top].down
				}
			}
			if (lastIndex_updown.i1 != undef) {
				push(dest, TPair(i, lastIndex_updown.i1));			// childtab[top].up goes to [top - 1]

//				::std::cout << "up:\t" << TPair(i, lastIndex_updown.i1) << ::std::endl;
				lastIndex_updown.i1 = undef;
			}
			stack_updown.push(TPair(i + 1, lcp_i));
		}

		endRead(lcpIn);
		endWrite(dest);
    }


	//////////////////////////////////////////////////////////////////////////////
    // Enhanced class (outputs only the childtab (3rd) column of the Enhanced Suffix Array)
    template < typename TLCPInput >
    struct Pipe< TLCPInput, ChildTab >
    {
        // *** SPECIALIZATION ***

		typedef typename Value<TLCPInput>::Type		TValue;
		typedef typename Size<TLCPInput>::Type		TSize;

		typedef Pair<TSize, TValue, Compressed>		TCoreType;		// (i, lcptab[i])

		typedef Pool< TCoreType, MapperSpec< MapperConfigSize< filterI1<TCoreType>, TSize > > > TLinearMapper;
        typedef Pipe< TLinearMapper, Filter< filterI2<TCoreType> > > TFilter;

		TLCPInput		*lcpIn;
        TLinearMapper   mapper;
		TFilter			in;
        
        Pipe():
            in(mapper) {}

        Pipe(TLCPInput &_in):
            lcpIn(&_in),
            in(mapper) {}

        inline void process() {
            process(*lcpIn);
        }

		template < typename _TLCPInput >
        inline bool process(_TLCPInput &_lcpIn) {

            // *** INSTANTIATION ***
			
			childtab_process(_lcpIn, mapper);
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
	template < typename TInput, typename _TLCPInput >
    inline bool operator<<(Pipe< TInput, ChildTab > &me, _TLCPInput const &in) {
 	    return me.process(in);
    }

	template < typename TLCPTable,
               typename TChildTable >
    void createChildTableExt(
		TChildTable &childtab,
		TLCPTable &lcp)
	{
		typedef Pipe< TLCPTable, Source<ContainerSpec> >	TSource;
		typedef Pipe< TSource, ChildTab >					TESA;

		TSource source(lcp);
		TESA	esa(source);

		esa.process();

		childtab << esa;
	}

/**
.Function.createChildTable:
..summary:Creates a child table from a given lcp table.
..cat:Index
..signature:createChildTable(childTab, lcp[, algo_tag])
..param.childTab:A reference to the resulting child table.
..param.lcp:A given lcp table.
..param.algo_tag:A tag that identifies the algorithm which is used for creation.
..remarks:The size of $childTab$ must be at least $length(text)$ before calling this function.
*/
	template < typename TLCPTable,
               typename TValue,
			   typename TConfig >
    inline void createChildTable(
		String<TValue, External<TConfig> > &childtab,
		TLCPTable &lcp)
	{
		createChildTableExt(childtab, lcp);
	}

	
	
    //////////////////////////////////////////////////////////////////////////////
    // internal childtab algorithm
    //////////////////////////////////////////////////////////////////////////////

	template < typename TLCPInput, typename TDest >
	inline void createChildTable(TDest &dest, TLCPInput const &lcpIn)
	{
		typedef typename Value<TLCPInput>::Type				TValue;
		typedef typename Size<TLCPInput>::Type				TSize;
		typedef typename Iterator<TLCPInput const>::Type	TIter;

		typedef Pair<TSize, TValue>							TPair;		// (i, lcptab[i])
		typedef ::std::stack<TPair> 						TStack;

		TStack stack_updown;
		TStack stack_nextl;

		stack_updown.push(TPair(0, 0));
		stack_nextl.push(TPair(0, 0));

		resize(dest, length(lcpIn));

		TSize n = length(lcpIn);
		TSize undef = n + 1;
		TSize lastIndex_nextl;
		TPair lastIndex_updown = TPair(undef, 0);

		TValue lcp_i;
		TIter lcpI = begin(lcpIn);
/*		for(TSize i = 0; i < n; ++i, ++lcpI) {
			lcp_i = *lcpI;

			//////////////////////////////////////////////////////////////////////////////
			// nextlIndex part
			lastIndex_nextl = undef;
			while (lcp_i < stack_nextl.top().i2)
				stack_nextl.pop();

			if (lcp_i == stack_nextl.top().i2) {
				lastIndex_nextl = stack_nextl.top().i1;
				dest[i + 1] = lastIndex_nextl;						// childtab[top].nextl

				::std::cout << "nextl:\t\t\t" << TPair(lastIndex_nextl, i + 1) << ::std::endl;
				stack_nextl.pop();
			}
			stack_nextl.push(TPair(i + 1, lcp_i));
		}
*/
		for(TSize i = 0; i < n; ++i, ++lcpI) {
			lcp_i = *lcpI;

			//////////////////////////////////////////////////////////////////////////////
			// up/down part
			TPair top;
			while (lcp_i < stack_updown.top().i2) {
				lastIndex_updown = stack_updown.top();
				stack_updown.pop();
				
				top = stack_updown.top();
				if (lcp_i <= top.i2 && top.i2 != lastIndex_updown.i2 /*&& top.i1 != lastIndex_nextl*/) {
					dest[top.i1] = lastIndex_updown.i1;				// childtab[top].down

//					::std::cout << "down:\t\t" << TPair(top.i1, lastIndex_updown.i1) << ::std::endl;			// childtab[top].down
				}
			}
			if (lastIndex_updown.i1 != undef) {
				dest[i] = lastIndex_updown.i1;						// childtab[top].up goes to [top - 1]

//				::std::cout << "up:\t" << TPair(i, lastIndex_updown.i1) << ::std::endl;
				lastIndex_updown.i1 = undef;
			}
			stack_updown.push(TPair(i + 1, lcp_i));
		}

		lcpI = begin(lcpIn);
		for(TSize i = 0; i < n; ++i, ++lcpI) {
			lcp_i = *lcpI;

			//////////////////////////////////////////////////////////////////////////////
			// nextlIndex part
			lastIndex_nextl = undef;
			while (lcp_i < stack_nextl.top().i2)
				stack_nextl.pop();

			if (lcp_i == stack_nextl.top().i2) {
				lastIndex_nextl = stack_nextl.top().i1;
				dest[lastIndex_nextl] = i + 1;						// childtab[top].nextl

//				::std::cout << "nextl:\t\t\t" << TPair(lastIndex_nextl, i + 1) << ::std::endl;
				stack_nextl.pop();
			}
			stack_nextl.push(TPair(i + 1, lcp_i));
		}

    }

//}

}

#endif
