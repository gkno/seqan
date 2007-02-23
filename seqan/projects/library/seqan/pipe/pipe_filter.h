/*
 *  pipe_filter.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_PIPE_FILTER_H
#define SEQAN_HEADER_PIPE_FILTER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{
    
    template <typename InType, typename Result = typename InType::T1>
    struct filterI1 : public ::std::unary_function<InType,Result> {
        inline Result operator()(const InType& x) const { return x.i1; }
    };

    template <typename InType, typename Result = typename InType::T2>
    struct filterI2 : public ::std::unary_function<InType,Result> {
        inline Result operator()(const InType& x) const { return x.i2; }
    };

    template <typename InType, typename Result = typename InType::T3>
    struct filterI3 : public ::std::unary_function<InType,Result> {
        inline Result operator()(const InType& x) const { return x.i3; }
    };


    template < typename TFunctor >
    struct Filter;

	template < typename TInput, typename TFunctor >
    struct Value< Pipe< TInput, Filter<TFunctor> > > {
		typedef typename TFunctor::result_type Type;
	};


/**
.Spec.Filter:
..cat:Pipelining
..general:Class.Pipe
..summary:Applies a specific function to the input stream.
..signature:Pipe<TInput, Filter<TFunctor> >
..param.TInput:The type of the pipeline module this module reads from.
..param.TFunctor:A unary function (see STL's $unary_function$).
...remarks:The argument type of $TFunctor$ must be $VALUE<TInput>::Type$.
..remarks: The output type of this pipe is the result type of $TFunctor$.
*/

	//////////////////////////////////////////////////////////////////////////////
    // filter class
    template <typename TInput, typename TFunctor >
    struct Pipe< TInput, Filter<TFunctor> >
    {
		TInput      &in;
        TFunctor    F;
        
/**
.Memfunc.Filter#Pipe:
..class:Spec.Filter
..summary:Constructor
..signature:Pipe<TInput, Filter<TFunctor> > (in)
..signature:Pipe<TInput, Filter<TFunctor> > (in, func)
..param.in:Reference to an input pipe.
..param.func:A TFunctor object (copy constructor).
*/
        Pipe(TInput& _in):
            in(_in) {}
        
        Pipe(TInput& _in, const TFunctor& _F) :
            in(_in),
            F(_F) {}
        
        inline typename Value<Pipe>::Type const operator*() const {
            return F(*in);
        }

        Pipe& operator++() {
            ++in;
            return *this;
        }
                
    };
    
//}

}

#endif
