#ifndef SEQAN_HEADER_BASIC_ALPHABET_TRAIT_BASIC_H
#define SEQAN_HEADER_BASIC_ALPHABET_TRAIT_BASIC_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
//arrayConstruct
//////////////////////////////////////////////////////////////////////////////

template<typename TIterator>
inline void 
_arrayConstruct_Pointer(TIterator begin_, 
						TIterator end_,
						True)
{
SEQAN_CHECKPOINT
	//nothing to do
}
template<typename TIterator>
inline void 
_arrayConstruct_Pointer(TIterator begin_, 
						TIterator end_,
						False)
{
SEQAN_CHECKPOINT
	_arrayConstruct_Default(begin_, end_);
}
template<typename TValue>
inline void 
arrayConstruct(TValue * begin_, 
			   TValue * end_)
{
SEQAN_CHECKPOINT
	_arrayConstruct_Pointer(begin_, end_, typename IsSimple<TValue>::Type() );
}

//____________________________________________________________________________

template<typename TIterator, typename TParam>
inline void 
_arrayConstruct_Pointer(TIterator begin_, 
						TIterator end_, 
						TParam const & param_,
						True)
{
SEQAN_CHECKPOINT
	arrayFill(begin_, end_, param_);
}
template<typename TIterator, typename TParam>
inline void 
_arrayConstruct_Pointer(TIterator begin_, 
						TIterator end_, 
						TParam const & param_,
						False)
{
SEQAN_CHECKPOINT
	_arrayConstruct_Default(begin_, end_, param_);
}
template<typename TValue, typename TParam>
inline void 
arrayConstruct(TValue * begin_, 
			   TValue * end_, 
			   TParam const & param_)
{
SEQAN_CHECKPOINT
	_arrayConstruct_Pointer(begin_, end_, param_, typename IsSimple<TValue>::Type() );
}

//////////////////////////////////////////////////////////////////////////////
//arrayConstructCopy
//////////////////////////////////////////////////////////////////////////////

template<typename TValue>
inline void 
_arrayConstructCopy_Pointer(TValue * source_begin, 
							TValue * source_end, 
							TValue * target_begin,
							True)
{
SEQAN_CHECKPOINT
	arrayCopyForward(source_begin, source_end, target_begin);
}
template<typename TValue>
inline void 
_arrayConstructCopy_Pointer(TValue * source_begin, 
							TValue * source_end, 
							TValue * target_begin,
							False)
{
SEQAN_CHECKPOINT
	_arrayConstructCopy_Default(source_begin, source_end, target_begin);
}
template<typename TValue>
inline void 
arrayConstructCopy(TValue * source_begin, 
				   TValue * source_end, 
				   TValue * target_begin)
{
SEQAN_CHECKPOINT
	_arrayConstructCopy_Pointer(source_begin, source_end, target_begin, typename IsSimple<TValue>::Type() );
}

//////////////////////////////////////////////////////////////////////////////
//arrayConstructMove
//////////////////////////////////////////////////////////////////////////////

template<typename TValue>
inline void 
_arrayConstructMove_Pointer(TValue * source_begin, 
							TValue * source_end, 
							TValue * target_begin,
							True)
{
SEQAN_CHECKPOINT
	arrayMoveForward(source_begin, source_end, target_begin);
}
template<typename TValue>
inline void 
_arrayConstructMove_Pointer(TValue * source_begin, 
							TValue * source_end, 
							TValue * target_begin,
							False)
{
SEQAN_CHECKPOINT
	_arrayConstructMove_Default(source_begin, source_end, target_begin);
}
template<typename TValue>
inline void 
arrayConstructMove(TValue * source_begin, 
				   TValue * source_end, 
				   TValue * target_begin)
{
SEQAN_CHECKPOINT
	_arrayConstructMove_Pointer(source_begin, source_end, target_begin, typename IsSimple<TValue>::Type() );
}

//////////////////////////////////////////////////////////////////////////////
//arrayDestruct
//////////////////////////////////////////////////////////////////////////////

template<typename TValue>
inline void 
_arrayDestruct_Pointer(TValue * begin_, 
					   TValue * end_,
					   True)
{
SEQAN_CHECKPOINT
	//do nothing
}
template<typename TValue>
inline void 
_arrayDestruct_Pointer(TValue * begin_, 
					   TValue * end_,
					   False)
{
SEQAN_CHECKPOINT
	_arrayDestruct_Default(begin_, end_);
}
template<typename TValue>
inline void 
arrayDestruct(TValue * begin_, 
			  TValue * end_)
{
SEQAN_CHECKPOINT
	_arrayDestruct_Pointer(begin_, end_, typename IsSimple<TValue>::Type() );
}

//////////////////////////////////////////////////////////////////////////////
//arrayFill
//////////////////////////////////////////////////////////////////////////////

//no specializiation for pointer to simple

//////////////////////////////////////////////////////////////////////////////
//arrayCopyForward
//////////////////////////////////////////////////////////////////////////////

template<typename TValue>
inline void 
_arrayCopyForward_Pointer(TValue * source_begin, 
						  TValue * source_end, 
						  TValue * target_begin,
						  True)
{
SEQAN_CHECKPOINT
	::std::memcpy(target_begin, source_begin, (source_end - source_begin) * sizeof(TValue));
}
template<typename TValue>
inline void 
_arrayCopyForward_Pointer(TValue * source_begin, 
						  TValue * source_end, 
						  TValue * target_begin,
						  False)
{
SEQAN_CHECKPOINT
	_arrayCopyForward_Default(source_begin, source_end, target_begin);
}
template<typename TValue>
inline void 
arrayCopyForward(TValue * source_begin, 
				 TValue * source_end, 
				 TValue * target_begin)
{
SEQAN_CHECKPOINT
	_arrayCopyForward_Pointer(source_begin, source_end, target_begin, typename IsSimple<TValue>::Type() );
}

//////////////////////////////////////////////////////////////////////////////
//arrayCopyBackward
//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline void 
_arrayCopyBackward_Pointer(TValue * source_begin, 
						   TValue * source_end, 
						   TValue * target_begin,
						   True)
{
SEQAN_CHECKPOINT
	::std::memmove(target_begin, source_begin, (source_end - source_begin) * sizeof(TValue));
}
template <typename TValue>
inline void 
_arrayCopyBackward_Pointer(TValue * source_begin, 
						   TValue * source_end, 
						   TValue * target_begin,
						   False)
{
SEQAN_CHECKPOINT
	_arrayCopyBackward_Default(source_begin, source_end, target_begin); 
}
template<typename TValue>
inline void 
arrayCopyBackward(TValue * source_begin, 
				  TValue * source_end, 
				  TValue * target_begin)
{
SEQAN_CHECKPOINT
	_arrayCopyBackward_Pointer(source_begin, source_end, target_begin, typename IsSimple<TValue>::Type() );
}

//////////////////////////////////////////////////////////////////////////////
//arrayMoveForward
//////////////////////////////////////////////////////////////////////////////

template<typename TValue>
inline void 
_arrayMoveForward_Pointer(TValue * source_begin, 
						  TValue * source_end, 
						  TValue * target_begin,
						  True)
{
SEQAN_CHECKPOINT
	::std::memcpy(target_begin, source_begin, (source_end - source_begin) * sizeof(TValue));
}
template<typename TValue>
inline void 
_arrayMoveForward_Pointer(TValue * source_begin, 
						  TValue * source_end, 
						  TValue * target_begin,
						  False)
{
SEQAN_CHECKPOINT
	_arrayMoveForward_Default(source_begin, source_end, target_begin);
}
template<typename TValue>
inline void 
arrayMoveForward(TValue * source_begin, 
				 TValue * source_end, 
				 TValue * target_begin)
{
SEQAN_CHECKPOINT
	_arrayMoveForward_Pointer(source_begin, source_end, target_begin, typename IsSimple<TValue>::Type() );
}

//////////////////////////////////////////////////////////////////////////////
//arrayMoveBackward
//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline void 
_arrayMoveBackward_Pointer(TValue * source_begin, 
						   TValue * source_end, 
						   TValue * target_begin,
						   True)
{
SEQAN_CHECKPOINT
	::std::memmove(target_begin, source_begin, (source_end - source_begin) * sizeof(TValue));
}
template <typename TValue>
inline void 
_arrayMoveBackward_Pointer(TValue * source_begin, 
						   TValue * source_end, 
						   TValue * target_begin,
						   False)
{
SEQAN_CHECKPOINT
	_arrayMoveBackward_Default(source_begin, source_end, target_begin); 
}
template<typename TValue>
inline void 
arrayMoveBackward(TValue * source_begin, 
				  TValue * source_end, 
				  TValue * target_begin)
{
SEQAN_CHECKPOINT
	_arrayMoveBackward_Pointer(source_begin, source_end, target_begin, typename IsSimple<TValue>::Type() );
}

//////////////////////////////////////////////////////////////////////////////
//arrayClearSpace
//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline void 
_arrayClearSpace_Pointer(TValue * array_begin, 
						size_t array_length, 
						size_t keep_from, 
						size_t move_to,
						True)
{
	if (keep_from == move_to) return;
SEQAN_CHECKPOINT
	arrayCopy(array_begin + keep_from, array_begin + array_length, array_begin + move_to);
}
template <typename TValue>
inline void 
_arrayClearSpace_Pointer(TValue * array_begin, 
						size_t array_length, 
						size_t keep_from, 
						size_t move_to,
						False)
{
	_arrayClearSpace_Default(array_begin, array_length, keep_from, move_to);
}
template <typename TValue>
void arrayClearSpace(TValue * array_begin, 
					 size_t array_length, 
					 size_t keep_from, 
					 size_t move_to)
{
	_arrayClearSpace_Pointer(array_begin, array_length, keep_from, move_to, typename IsSimple<TValue>::Type() );
}


//////////////////////////////////////////////////////////////////////////////
// bool 
//////////////////////////////////////////////////////////////////////////////

template <> struct IsSimple< bool > { typedef True Type; };
template <> struct IsSimple< bool const > { typedef True Type; };

template <> struct BitsPerValue< bool > { enum { VALUE = 1 }; };
template <> struct BitsPerValue< bool const> { enum { VALUE = 1 }; };


//////////////////////////////////////////////////////////////////////////////
// char 
//////////////////////////////////////////////////////////////////////////////

template <> struct IsSimple< char > { typedef True Type; };
template <> struct IsSimple< char const > { typedef True Type; };

inline char const &
supremumValueImpl(char *)
{
SEQAN_CHECKPOINT
	static char const _value = (char) 127;
	return _value;
}
inline char const &
infimumValueImpl(char *)
{
SEQAN_CHECKPOINT
	static char const _value = (char) -128;
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// signed char 
//////////////////////////////////////////////////////////////////////////////

template <> struct IsSimple< signed char > { typedef True Type; };
template <> struct IsSimple< signed char const > { typedef True Type; };

inline signed char const &
supremumValueImpl(signed char *)
{
SEQAN_CHECKPOINT
	static signed char const _value = 127;
	return _value;
}
inline signed char const &
infimumValueImpl(signed char *)
{
SEQAN_CHECKPOINT
	static signed char const _value = -128;
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// unsigned char 
//////////////////////////////////////////////////////////////////////////////

template <> struct IsSimple< unsigned char > { typedef True Type; };
template <> struct IsSimple< unsigned char const > { typedef True Type; };

inline unsigned char const &
supremumValueImpl(unsigned char *)
{
SEQAN_CHECKPOINT
	static unsigned char const _value = 255;
	return _value;
}
inline unsigned char const &
infimumValueImpl(unsigned char *)
{
SEQAN_CHECKPOINT
	static unsigned char const _value = 0;
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// wchar_t
//////////////////////////////////////////////////////////////////////////////

template <> struct IsSimple< wchar_t > { typedef True Type; };
template <> struct IsSimple< wchar_t const > { typedef True Type; };

inline wchar_t const &
supremumValueImpl(wchar_t *)
{
SEQAN_CHECKPOINT
	static wchar_t const _value = 1UL << (BitsPerValue<wchar_t>::VALUE) - 1;
	return _value;
}
inline wchar_t const &
infimumValueImpl(wchar_t *)
{
SEQAN_CHECKPOINT
	static wchar_t const _value = 0;
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// signed short 
//////////////////////////////////////////////////////////////////////////////

template <> struct IsSimple< signed short > { typedef True Type; };
template <> struct IsSimple< signed short const > { typedef True Type; };

inline signed short const &
supremumValueImpl(signed short *)
{
SEQAN_CHECKPOINT
	static signed short const _value = 1 << (BitsPerValue<signed short>::VALUE - 1) - 1;
	return _value;
}
inline signed short const &
infimumValueImpl(signed short *)
{
SEQAN_CHECKPOINT
	static signed short const _value = - (1 << (BitsPerValue<signed short>::VALUE - 1));
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// unsigned short 
//////////////////////////////////////////////////////////////////////////////

template <> struct IsSimple< unsigned short > { typedef True Type; };
template <> struct IsSimple< unsigned short const > { typedef True Type; };

inline unsigned short const &
supremumValueImpl(unsigned short *)
{
SEQAN_CHECKPOINT
	static unsigned short const _value = 1 << BitsPerValue<unsigned short>::VALUE - 1;
	return _value;
}
inline unsigned short const &
infimumValueImpl(unsigned short *)
{
SEQAN_CHECKPOINT
	static unsigned short const _value = 0;
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// signed int 
//////////////////////////////////////////////////////////////////////////////

template <> struct IsSimple< signed int > { typedef True Type; };
template <> struct IsSimple< signed int const > { typedef True Type; };

inline signed int const &
supremumValueImpl(signed int *)
{
SEQAN_CHECKPOINT
	static signed int const _value = 1 << (BitsPerValue<signed int>::VALUE - 1) - 1;
	return _value;
}
inline signed int const &
infimumValueImpl(signed int *)
{
SEQAN_CHECKPOINT
	static signed int const _value = - (1 << (BitsPerValue<signed int>::VALUE - 1));
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// unsigned int 
//////////////////////////////////////////////////////////////////////////////

template <> struct IsSimple< unsigned int > { typedef True Type; };
template <> struct IsSimple< unsigned int const > { typedef True Type; };

inline unsigned int const &
supremumValueImpl(unsigned int *)
{
SEQAN_CHECKPOINT
	static unsigned int const _value = ~0ul;
	return _value;
}
inline unsigned int const &
infimumValueImpl(unsigned int *)
{
SEQAN_CHECKPOINT
	static unsigned int const _value = 0;
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// signed long 
//////////////////////////////////////////////////////////////////////////////

template <> struct IsSimple< signed long > { typedef True Type; };
template <> struct IsSimple< signed long const > { typedef True Type; };

inline signed long const &
supremumValueImpl(signed long *)
{
SEQAN_CHECKPOINT
	static signed long const _value = 1 << (BitsPerValue<signed long>::VALUE - 1) - 1;
	return _value;
}
inline signed long const &
infimumValueImpl(signed long *)
{
SEQAN_CHECKPOINT
	static signed long const _value = - (1 << (BitsPerValue<signed long>::VALUE - 1));
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// unsigned int 
//////////////////////////////////////////////////////////////////////////////

template <> struct IsSimple< unsigned long > { typedef True Type; };
template <> struct IsSimple< unsigned long const > { typedef True Type; };

inline unsigned long const &
supremumValueImpl(unsigned long *)
{
SEQAN_CHECKPOINT
	static unsigned long const _value = ~0ul;
	return _value;
}
inline unsigned long const &
infimumValueImpl(unsigned long *)
{
SEQAN_CHECKPOINT
	static unsigned long const _value = 0;
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// float 
//////////////////////////////////////////////////////////////////////////////

template <> struct IsSimple< float > { typedef True Type; };
template <> struct IsSimple< float const > { typedef True Type; };

inline float const &
supremumValueImpl(float *)
{
SEQAN_CHECKPOINT
static float const _value = ::std::numeric_limits<float>::infinity( );
	return _value;
}
inline float const &
infimumValueImpl(float *)
{
SEQAN_CHECKPOINT
	static float const _value = - ::std::numeric_limits<float>::infinity( );
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// double 
//////////////////////////////////////////////////////////////////////////////

template <> struct IsSimple< double > { typedef True Type; };
template <> struct IsSimple< double const > { typedef True Type; };

inline double const &
supremumValueImpl(double *)
{
SEQAN_CHECKPOINT
static double const _value = ::std::numeric_limits<double>::infinity( );
	return _value;
}
inline double const &
infimumValueImpl(double *)
{
SEQAN_CHECKPOINT
	static double const _value = - ::std::numeric_limits<double>::infinity( );
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// long double 
//////////////////////////////////////////////////////////////////////////////

template <> struct IsSimple< long double > { typedef True Type; };
template <> struct IsSimple< long double const > { typedef True Type; };

inline long double const &
supremumValueImpl(long double *)
{
SEQAN_CHECKPOINT
static long double const _value = ::std::numeric_limits<long double>::infinity( );
	return _value;
}
inline long double const &
infimumValueImpl(long double *)
{
SEQAN_CHECKPOINT
	static long double const _value = - ::std::numeric_limits<long double>::infinity( );
	return _value;
}
//////////////////////////////////////////////////////////////////////////////
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
