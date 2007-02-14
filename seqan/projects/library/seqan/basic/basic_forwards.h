#ifndef SEQAN_HEADER_BASIC_FORWARD2_H
#define SEQAN_HEADER_BASIC_FORWARD2_H

//forward declarations (make GCC 4.x happy)

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// basic_transport.h::assign

template <typename TTarget, typename TSource>
inline void
assign(TTarget & target,
	   TSource & source);

template <typename TTarget, typename TSource>
inline void
assign(TTarget & target,
	   TSource const & source);

//////////////////////////////////////////////////////////////////////////////
// string_pointer.h::assignValue

template <typename TValue, typename TPos>
inline void
assignValue(TValue * me,
			TPos pos, 
			TValue const & _value);

//////////////////////////////////////////////////////////////////////////////
// string_pointer.h::moveValue

template <typename TValue, typename TPos>
inline void
moveValue(TValue * me,
			TPos pos, 
			TValue const & _value);

//////////////////////////////////////////////////////////////////////////////
// string_pointer.h::value

template <typename TValue, typename TPos>
inline TValue &
value(TValue * me,
	  TPos pos);

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
