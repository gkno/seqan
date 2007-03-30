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

//////////////////////////////////////////////////////////////////////////////

    template < unsigned tupleLen, bool omitLast = false, typename TCompression = void >
    struct Tupler;

    template < typename TInput, unsigned tupleLen, bool omitLast, typename TCompression >
    struct Value< Pipe< TInput, Tupler< tupleLen, omitLast, TCompression > > > {
        typedef Tuple<typename Value<TInput>::Type, tupleLen, TCompression>	TTuple;
        typedef Pair<typename Size<TInput>::Type, TTuple, Compressed>		Type;
    };

//////////////////////////////////////////////////////////////////////////////

    template < 
		typename TInput, 
		unsigned tupleLen, 
		bool omitLast, 
		typename TCompression,
		typename TPair, 
		typename TLimitsString >
    struct Value< Pipe< TInput, Multi< Tupler< tupleLen, omitLast, TCompression >, TPair, TLimitsString > > > {
        typedef Tuple<typename Value<TInput>::Type, tupleLen, TCompression>	TTuple;
        typedef Pair<TPair, TTuple, Compressed>								Type;
    };

//////////////////////////////////////////////////////////////////////////////


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

    struct _ShiftLeftWorker {
        template <typename Arg>
        static inline void body(Arg &arg, unsigned I) {
            arg.i2[I-1] = arg.i2[I];
        }
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
    // tupler class
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

//____________________________________________________________________________


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
    // tupler class for multiple sequences
    template < 
		typename TInput, 
		unsigned tupleLen, 
		bool omitLast, 
		typename TCompression, 
		typename TPair, 
		typename TLimitsString >
    struct Pipe< TInput, Multi<Tupler<tupleLen, omitLast, TCompression>, TPair, TLimitsString> >
    {
		typedef typename Value< typename Value<Pipe>::Type, 2 >::Type	TTuple;
		typedef typename Value<TTuple>::Type							TValue;

		typedef _PairIncrementer<TPair, TLimitsString>	Incrementer;

		TInput                      &in;
        Incrementer					localPos;
        typename Value<Pipe>::Type	tmp;
		typename Size<TInput>::Type	seqLength, lastTuples;

		TLimitsString const &limits;
        
        Pipe(TInput& _in, TLimitsString const &_limits):
            in(_in),
			limits(_limits) {}

        inline typename Value<Pipe>::Type const & operator*() const {
            return tmp;
        }

        inline Pipe& operator++() {
			// process next sequence
			if (--lastTuples == 0) {
				fill();
				return *this;
			}

			// shift left 1 character
            LOOP<_ShiftLeftWorker, tupleLen - 1>::run(this->tmp);
			++localPos;
			tmp.i1 = localPos;
			if (lastTuples < _TuplerLastTuples<Pipe>::VALUE) {
	            tmp.i2[tupleLen - 1] = TValue();
			} else {
				tmp.i2[tupleLen - 1] = *in;
				++in;
			}
            return *this;
        }

        inline void fill() {
			unsigned i = 0;
			do {
				while (i > 0)
					++localPos;

				for(; i < tupleLen && !eos(); ++i, ++in) {
					tmp.i2.i[i] = *in;
				}
				lastTuples = _TuplerLastTuples<Pipe>::VALUE;

				// fill up with null chars
				for(; i < tupleLen; ++i)
					tmp.i2.i[i] = TValue();
				
				// eventually, reduce the number of half-filled tuples
				if (lastTuples <= tupleLen - i)
					lastTuples = 0;
				else
					lastTuples -= tupleLen - i;
			} while ((lastTuples == 0) && !eof(in));

			tmp.i1 = localPos;
        }

		inline bool eos() {
			return (getValueI1(localPos) > 0) && (getValueI2(localPos) == 0);
		}
	};

//____________________________________________________________________________


	template < 
		typename TInput, 
		unsigned tupleLen, 
		bool omitLast, 
		typename TPair, 
		typename TLimitsString >
    struct Pipe< TInput, Multi<Tupler<tupleLen, omitLast, Compressed>, TPair, TLimitsString> >
    {
		typedef typename Value< typename Value<Pipe>::Type, 2 >::Type	TTuple;
		typedef typename Value<TTuple>::Type							TValue;

		typedef _PairIncrementer<TPair, TLimitsString>	Incrementer;

		TInput                      &in;
        Incrementer					localPos;
        typename Value<Pipe>::Type	tmp;
		typename Size<TInput>::Type	seqLength, lastTuples;

		TLimitsString const &limits;
        
        Pipe(TInput& _in, TLimitsString const &_limits):
            in(_in),
			limits(_limits) {}

        inline typename Value<Pipe>::Type const & operator*() const {
            return tmp;
        }

        inline Pipe& operator++() {
			// process next sequence
			if (eos())
				if (--lastTuples == 0) {
					assignValueI1(tmp.i1, getValueI1(tmp.i1) + 1);
					fill();
					return *this;
				}

			// shift left 1 character
			tmp.i2 <<= 1;
			assignValueI2(tmp.i1, getValueI2(tmp.i1) + 1);
			if (lastTuples == _TuplerLastTuples<Pipe>::VALUE) {
				tmp.i2 |= *in;
				++localPos;
				++in;
			}
            return *this;
        }

        inline void fill() {
			do {
				unsigned i = 0;
				if (!eof(in))
					do {
						tmp.i2 <<= 1;
						tmp.i2 |= *in;
						++in;
						++i;
						++localPos;
					} while ((i < tupleLen) && !eos());
				lastTuples = _TuplerLastTuples<Pipe>::VALUE;

				// fill up with null chars
	            tmp.i2 <<= (tupleLen - i);
				
				// eventually, reduce the number of half-filled tuples
				if (lastTuples <= tupleLen - i)
					lastTuples = 0;
				else
					lastTuples -= tupleLen - i;

				if (lastTuples == 0)
					assignValueI1(tmp.i1, getValueI1(tmp.i1) + 1);

			} while ((lastTuples == 0) && !eof(in));

			assignValueI2(tmp.i1, 0);
        }

		inline bool eos() {
			return (getValueI1(value(localPos)) > 0) && (getValueI2(value(localPos)) == 0);
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
    
    template < 
		typename TInput,
		unsigned tupleLen,
		bool omitLast,
		typename TCompression,
		typename TPair, 
		typename TLimitsString >
	inline bool 
	control(
		Pipe< TInput, Multi<Tupler< tupleLen, omitLast, TCompression >, TPair, TLimitsString> > &me, 
		ControlBeginRead const &command) 
	{
        if (!control(me.in, command)) return false;
		setHost(me.localPos, me.limits);
		assignValueI1(me.tmp.i1, 0);
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

    template < 
		typename TInput,
		unsigned tupleLen,
		bool omitLast,
		typename TCompression,
		typename TPair, 
		typename TLimitsString >
	inline bool 
	control(
		Pipe< TInput, Multi<Tupler< tupleLen, omitLast, TCompression >, TPair, TLimitsString> > &me, 
		ControlEof const &command) 
	{
		return me.lastTuples == 0;
	}

    template < 
		typename TInput,
		unsigned tupleLen,
		bool omitLast,
		typename TCompression,
		typename TPair, 
		typename TLimitsString >
	inline bool 
	control(
		Pipe< TInput, Multi<Tupler< tupleLen, omitLast, TCompression >, TPair, TLimitsString> > &me, 
		ControlEos const &command) 
	{
		return me.eos();
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
