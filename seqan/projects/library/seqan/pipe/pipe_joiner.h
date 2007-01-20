/*
 *  pipe_joiner.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_PIPE_JOINER_H
#define SEQAN_HEADER_PIPE_JOINER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

	struct Joiner;

	template < typename TInput1, typename TInput2 >
    struct Value< Pipe< Bundle2< TInput1, TInput2 >, Joiner > > {
		typedef Pair<
			typename Value<TInput1>::Type,
			typename Value<TInput2>::Type,
			Compressed
		> Type;
	};

	template < typename TInput1, typename TInput2, typename TInput3 >
    struct Value< Pipe< Bundle3< TInput1, TInput2, TInput3 >, Joiner > > {
		typedef Triple<
			typename Value<TInput1>::Type,
			typename Value<TInput2>::Type,
			typename Value<TInput3>::Type,
			Compressed
		> Type;
	};


/**
.Spec.Joiner:
..cat:Pipelining
..general:Class.Pipe
..summary:Joins two or three input streams.
..signature:Pipe<Bundle2<TInput1, TInput2>, Joiner>
..signature:Pipe<Bundle3<TInput1, TInput2, TInput3>, Joiner>
..param.TInput1:The type of the first pipeline module this module reads from.
..param.TInput2:The type of the second pipeline module this module reads from.
..param.TInput3:The type of the third pipeline module this module reads from.
..remarks: The output type is a compressed @Class.Pair@ or @Class.Triple@ of the input types $Value<TInputX>::Type$.
*/

    //////////////////////////////////////////////////////////////////////////////
    // joiner class
	template < typename TInput1, typename TInput2 >
    struct Pipe< Bundle2< TInput1, TInput2 >, Joiner >
    {
		Bundle2< TInput1, TInput2 >	in;
        typename Value<Pipe>::Type	tmp;
        
        Pipe(Bundle2< TInput1, TInput2 > _in):
            in(_in) {}
        
        inline typename Value<Pipe>::Type const & operator*() {
            tmp.i1 = *in.in1;
            tmp.i2 = *in.in2;
            return tmp;
        }
        
        inline Pipe& operator++() {
            ++in.in1;
            ++in.in2;
            return *this;
        }
    };

	template < typename TInput1, typename TInput2, typename TInput3 >
    struct Pipe< Bundle3< TInput1, TInput2, TInput3 >, Joiner >
    {
		Bundle3< TInput1, TInput2, TInput3 >	in;
        typename Value<Pipe>::Type				tmp;
        
        Pipe(Bundle3< TInput1, TInput2, TInput3 > _in):
            in(_in) {}
        
        inline typename Value<Pipe>::Type const & operator*() {
            tmp.i1 = *in.in1;
            tmp.i2 = *in.in2;
            tmp.i3 = *in.in2;
            return tmp;
        }
        
        inline Pipe& operator++() {
            ++in.in1;
            ++in.in2;
            ++in.in3;
            return *this;
        }
    };

//}

}

#endif
