/*
 *  pipe_caster.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_PIPE_CASTER_H
#define SEQAN_HEADER_PIPE_CASTER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{
    
    template < typename TValue >
    struct Caster;

	template < typename TInput, typename TValue >
    struct Value< Pipe< TInput, Caster<TValue> > > {
		typedef TValue Type;
	};


/**
.Spec.Caster:
..cat:Pipelining
..general:Class.Pipe
..summary:Casts the input type in a specific output type.
..signature:Pipe<TInput, Caster<TValue> >
..param.TInput:The type of the pipeline module this module reads from.
..param.TValue:The new output type.
..remarks: The input stream is casted using $reinterpret_cast<TValue>$.
*/

    //////////////////////////////////////////////////////////////////////////////
    // filter class
    template <typename TInput, typename TValue >
    struct Pipe< TInput, Caster<TValue> >
    {
		TInput      &in;
        
        Pipe(TInput& _in):
            in(_in) {}
        
        inline TValue const & operator*() const {
            return reinterpret_cast<TValue const &>(*in);
        }

        Pipe& operator++() {
            ++in;
            return *this;
        }
                
    };
    
//}

}

#endif
