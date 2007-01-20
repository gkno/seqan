/*
 *  pipe_counter.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_PIPE_COUNTER_H
#define SEQAN_HEADER_PIPE_COUNTER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

	struct Counter;

	template < typename TInput >
    struct Value< Pipe< TInput, Counter > > {
		typedef Pair<
			typename Value<TInput>::Type,
			typename Size<TInput>::Type,
			Compressed
		> Type;
	};


/**
.Spec.Counter:
..cat:Pipelining
..general:Class.Pipe
..summary:Extends the input stream by a second field which enumerates the elements.
..signature:Pipe<TInput, Counter>
..param.TInput:The type of the pipeline module this module reads from.
..remarks:The output type is a @Class.Pair@ of input type and size type (i.e. $Pair<Value<TInput>::Type, Size<TInput>::Type>$).
..remarks:The first output field is the original input stream.
..remarks:The second output field begins with 0 and increases by 1 per element.
*/

    //////////////////////////////////////////////////////////////////////////////
    // counter class
    template < typename TInput >
    struct Pipe< TInput, Counter >
    {
		TInput                      &in;
        typename Value<Pipe>::Type	tmp;
        
        Pipe(TInput& _in):
            in(_in) {}
        
        inline typename Value<Pipe>::Type const & operator*() {
            tmp.i1 = *in;
            return tmp;
        }
        
        inline Pipe& operator++() {
            ++in;
            ++tmp.i2;
            return *this;
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput >
	inline bool control(Pipe< TInput, Counter > &me, ControlBeginRead const &command) {
        if (!control(me.in, command)) return false;
        me.tmp.i2 = 0;
		return true;
	}
    
//}

}

#endif
