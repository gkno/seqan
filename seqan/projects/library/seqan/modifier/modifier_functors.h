/*
 *  modifier_functors.h
 *  SeqAn
 *
 *  Created by David Weese on 26.01.06.
 *
 */

#ifndef SEQAN_HEADER_MODIFIER_FUNCTORS_H
#define SEQAN_HEADER_MODIFIER_FUNCTORS_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	// text transformation
	//////////////////////////////////////////////////////////////////////////////

    template <typename InType, typename Result = InType>
	struct FunctorUpcase : public ::std::unary_function<InType,Result> 
	{
        inline Result operator()(InType x) const {
			if (('a' <= x) && (x <= 'z')) return (x + ('A' - 'a'));
			return x; 
		}
    };

    template <typename InType, typename Result = InType>
    struct FunctorDowncase : public ::std::unary_function<InType,Result> 
	{
        inline Result operator()(InType x) const {
			if (('A' <= x) && (x <= 'Z')) return (x + ('a' - 'A'));
			return x; 
		}
    };


	//////////////////////////////////////////////////////////////////////////////
	// alphabet transformation
	//////////////////////////////////////////////////////////////////////////////

    template <typename InType, typename OutType>
    struct FunctorConvert : public ::std::unary_function<InType,OutType> 
	{
        inline OutType operator()(InType x) const {
			return x; 
		}
    };


}

#endif
