/*
 *  pipe_tupler.h
 *  SeqAn
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_PIPE_TUPLER_H
#define SEQAN_HEADER_PIPE_TUPLER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    struct _ShiftLeftWorker {
        template <typename Arg>
        static inline void body(Arg &arg, unsigned I) {
            arg.i2[I-1] = arg.i2[I];
        }
    };


    template < unsigned tupleLen, bool omitLast = false, typename TCompression = void >
    struct Tupler;

    template < typename TInput, unsigned tupleLen, bool omitLast, typename TCompression >
    struct Value< Pipe< TInput, Tupler< tupleLen, omitLast, TCompression > > > {
        typedef Tuple<typename Value<TInput>::Type, tupleLen, TCompression>	TTuple;
        typedef Pair<typename Size<TInput>::Type, TTuple, Compressed>		Type;
    };


	// output only fully filled tuples
	template < typename TTupler >
	struct _TuplerLastTuples {
		enum { VALUE = 1 };
	};

	// output tupleLen-1 half filled tuples at the end
    template < typename TInput, unsigned tupleLen, typename TCompression >
	struct _TuplerLastTuples< Pipe< TInput, Tupler<tupleLen, false, TCompression> > > {
		enum { VALUE = tupleLen };
	};

/**
.Spec.Tupler:
..cat:Pipelining
..general:Class.Pipe
..summary:Outputs tuples of the $tupleLen$ consecutive elements of the input stream.
..signature:Pipe<TInput, Tupler<tupleLen, omitLast> >
..param.TInput:The type of the pipeline module this module reads from.
..param.tupleLen:The tuple length.
...remarks:The tuples contain elements $in[i]in[i+1]...in[i+(tupleLen-1)]$.
..param.omitLast:Omit half filled tuples.
..param.omitLast:If $true$, the output stream is $tupleLen-1$ elements shorter than the input stream.
..param.omitLast:If $false$, the lengths are identical and the last tuples are filled with blanks (default constructed elements) for undefined entries.
..remarks:The output type is a @Class.Tuple@ of input elements and length $tupleLen$ (i.e. $Tuple<Value<TInput>::Type, tupleLen>$).
..remarks:The tuples are sequences of the form $in[i]in[i-1]in[i-2]..in[i-tupleLen+1]$. For $omitLast=false$ $i$ begins with 0 and for $omitLast=true$ $i$ begins with $tupleLen-1$.
*/

    //////////////////////////////////////////////////////////////////////////////
    // echoer class
    template < typename TInput, unsigned tupleLen, bool omitLast, typename TCompression >
    struct Pipe< TInput, Tupler<tupleLen, omitLast, TCompression> >
    {
		typedef typename Value< typename Value<Pipe>::Type, 2 >::Type	TTuple;
		typedef typename Value<TTuple>::Type							TValue;

		TInput                      &in;
        typename Value<Pipe>::Type	tmp;
		typename Size<TInput>::Type	lastTuples;
        
        Pipe(TInput& _in):
            in(_in) {}

        inline typename Value<Pipe>::Type const & operator*() const {
            return tmp;
        }

        inline Pipe& operator++() {
            if (eof(in)) --lastTuples;
            LOOP<_ShiftLeftWorker, tupleLen - 1>::run(this->tmp);
			++tmp.i1;
			if (lastTuples < _TuplerLastTuples<Pipe>::VALUE)
	            tmp.i2[tupleLen - 1] = TValue();
			else {
				tmp.i2[tupleLen - 1] = *in;
				++in;
			}
            return *this;
        }

        inline void fill() {
            unsigned i;
            for(i = 0; i < tupleLen && !eof(in); ++i, ++in)
                tmp.i2.i[i] = *in;
			lastTuples = eof(in)? 0: _TuplerLastTuples<Pipe>::VALUE;
            for(; i < tupleLen; ++i)
                tmp.i2.i[i] = TValue();
            tmp.i1 = 0;
        }
	};

    template < typename TInput, unsigned tupleLen, bool omitLast >
    struct Pipe< TInput, Tupler<tupleLen, omitLast, Compressed> >
    {
        TInput                      &in;
        typename Value<Pipe>::Type	tmp;
		typename Size<TInput>::Type	lastTuples;
        
        Pipe(TInput& _in):
            in(_in) {}

        inline typename Value<Pipe>::Type const & operator*() const {
            return tmp;
        }

        inline Pipe& operator++() {
            if (eof(in)) --lastTuples;
			tmp.i2 <<= 1;
			++tmp.i1;
			if (lastTuples == _TuplerLastTuples<Pipe>::VALUE) {
				tmp.i2 |= *in;
				++in;
			}
            return *this;
        }

        inline void fill() {
            unsigned i;
            clear(tmp.i2);
			for(i = 0; i < tupleLen && !eof(in); ++i, ++in) {
                tmp.i2 <<= 1;
                tmp.i2 |= *in;
			}
			lastTuples = eof(in)? 0: _TuplerLastTuples<Pipe>::VALUE;
            tmp.i2 <<= (tupleLen - i);
            tmp.i1 = 0;
        }
	};


    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput, unsigned tupleLen, bool omitLast, typename TCompression >
	inline bool 
	control(
		Pipe< TInput, Tupler< tupleLen, omitLast, TCompression > > &me, 
		ControlBeginRead const &command) 
	{
        if (!control(me.in, command)) return false;
		me.fill();
		return true;
	}
    
    template < typename TInput, unsigned tupleLen, bool omitLast, typename TCompression >
	inline bool 
	control(
		Pipe< TInput, Tupler< tupleLen, omitLast, TCompression > > &me, 
		ControlEof const &command) 
	{
		return me.lastTuples == 0;
    }

    template < typename TInput, unsigned tupleLen, bool omitLast, typename TCompression >
    inline typename Size< Pipe< TInput, Tupler< tupleLen, omitLast, TCompression > > >::Type
    length(Pipe< TInput, Tupler< tupleLen, omitLast, TCompression > > const &me) {
		typedef Pipe< TInput, Tupler< tupleLen, omitLast, TCompression > >	TPipe;
		
		if (length(me.in) >= (tupleLen - _TuplerLastTuples<TPipe>::VALUE))
			return length(me.in) - (tupleLen - _TuplerLastTuples<TPipe>::VALUE);
		else
			return 0;
    }

//}

}

#endif
