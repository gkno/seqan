#ifndef SEQAN_HEADER_BASIC_DEFINITION_H
#define SEQAN_HEADER_BASIC_DEFINITION_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct Tag
{
};

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Default:
..summary:Tag that specifies default behavior.
..value.Default:Use default behavior. 

*/

struct Default_;
typedef Tag<Default_> const Default;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Move Switch:
..summary:Switch to force move.
..value.Move:Move instead of assign. 
..remarks.text:The difference between move constructor and copy constructor
is that the source object is not copied but moved into the target object.
The source object can lose its content and will be empty after
this operation in this case.
A move constructor can sigificantly faster than a copy constructor.
..example.code:String source("hello");
String target(source, Move()); // source is moved to target
std::cout << source; //nothing printed since source lost content
std::cout << target; //"hello"
..see:Function.move
*/

struct Move_;
typedef Tag<Move_> const Move;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Logical Values:
..summary:Tag that represents true and false.
..value.True:The logical value "true".
..value.False:The logical value "false".
*/
struct True {};
struct False {};


//////////////////////////////////////////////////////////////////////////////

///Empty Data Class.
struct Nothing {};


//////////////////////////////////////////////////////////////////////////////

/**
.Internal._RemoveConst:
..signature:_RemoveConst<T>
..returns:$t$ if $T$ is $t const$, otherwise $T$.
*/
template <typename T>
struct _RemoveConst
{
	typedef T Type;
};
template <typename T>
struct _RemoveConst<T const>
{
	typedef T Type;
};


//////////////////////////////////////////////////////////////////////////////
/**
.Internal._MakeUnsigned:
..signature:_MakeUnsigned<T>
..returns:$unsigned t$ if $T$ is not $unsigned t$, otherwise $T$.
*/
template <typename T>
struct _MakeUnsigned
{
	typedef T Type;
};

template <>
struct _MakeUnsigned<char>
{
	typedef unsigned char Type;
};

template <>
struct _MakeUnsigned<const char>
{
	typedef const unsigned char Type;
};

template <>
struct _MakeUnsigned<signed char>
{
	typedef unsigned char Type;
};

template <>
struct _MakeUnsigned<const signed char>
{
	typedef const unsigned char Type;
};

template <>
struct _MakeUnsigned<int>
{
	typedef unsigned int Type;
};

template <>
struct _MakeUnsigned<const int>
{
	typedef const unsigned int Type;
};

template <>
struct _MakeUnsigned<short>
{
	typedef unsigned short Type;
};

template <>
struct _MakeUnsigned<const short>
{
	typedef const unsigned short Type;
};

template <>
struct _MakeUnsigned<long>
{
	typedef unsigned long Type;
};

template <>
struct _MakeUnsigned<const long>
{
	typedef const unsigned long Type;
};
/*
template <>
struct _MakeUnsigned<long long>
{
	typedef unsigned long long Type;
};

template <>
struct _MakeUnsigned<const long long>
{
	typedef const unsigned long long Type;
};
*/

//////////////////////////////////////////////////////////////////////////////
/**
.Internal._MakeSigned:
..signature:_MakeSigned<T>
..returns:$signed t$ if $T$ is not $signed t$, otherwise $T$.
*/
template <typename T>
struct _MakeSigned
{
	typedef T Type;
};

template <>
struct _MakeSigned<char>
{
	typedef signed char Type;
};

template <>
struct _MakeSigned<const char>
{
	typedef const signed char Type;
};

template <>
struct _MakeSigned<signed char>
{
	typedef signed char Type;
};

template <>
struct _MakeSigned<const signed char>
{
	typedef const signed char Type;
};

template <>
struct _MakeSigned<int>
{
	typedef signed int Type;
};

template <>
struct _MakeSigned<const int>
{
	typedef const signed int Type;
};

template <>
struct _MakeSigned<short>
{
	typedef signed short Type;
};

template <>
struct _MakeSigned<const short>
{
	typedef const signed short Type;
};

template <>
struct _MakeSigned<long>
{
	typedef signed long Type;
};

template <>
struct _MakeSigned<const long>
{
	typedef const signed long Type;
};
/*
template <>
struct _MakeSigned<long long>
{
	typedef signed long long Type;
};

template <>
struct _MakeSigned<const long long>
{
	typedef const signed long long Type;
};
*/

//////////////////////////////////////////////////////////////////////////////
/**
.Internal._ClassIdentifier:
..signature:void * _ClassIdentifier<T>::getID()
..returns:A void * that identifies $T$.
...text:The returned values of two calls of $getID$ are equal if and only if
the used type $T$ was the same.
*/
template <typename T>
struct _ClassIdentifier
{
	static inline void *
	getID()
	{
SEQAN_CHECKPOINT
		static bool _id_dummy;
		return &_id_dummy;
	}
};

//////////////////////////////////////////////////////////////////////////////
/**
.Function.log2:
..cat:Miscellaneous
..summary:Computes logarithm of base 2 for integer types
..signature:unsigned int log2(i)
..param.i:An integer type.
..returns:The largest integer smaller or equal than
the logarithm of $i$.
*/

template <int BITS_MAX>
struct _Log2_Impl
{
	template <typename T>
	static inline unsigned int
	log2(T val, unsigned int offset)
	{
		unsigned int val2 = val >> (BITS_MAX / 2);
		if (val2)
		{
			val = val2;
			offset += BITS_MAX / 2;
		}
		return _Log2_Impl<BITS_MAX / 2>::log2(val, offset);
	}
};

template <>
struct _Log2_Impl<1>
{
	template <typename T>
	static inline unsigned int
	log2(T val, unsigned int offset)
	{
		return offset;
	}
};


template <typename T>
inline unsigned int
log2(T val)
{
	enum
	{
//		BITS_PER_VALUE = BitsPerValue<T>::VALUE //TODO???
		BITS_PER_VALUE = sizeof(T) * 8
	};

	return _Log2_Impl<BITS_PER_VALUE>::log2(val, 0);
}

//////////////////////////////////////////////////////////////////////////////
// to avoid conflicts with non-standard macros and namespaces
// we define our own Min/Max functions

template<typename _Tx> inline
const _Tx& Min(const _Tx& _Left, const _Tx& _Right)
{	// return smaller of _Left and _Right
    return (_Right < _Left ? _Right : _Left);
}

template<typename _Tx, typename _Ty> inline
_Tx Min(const _Tx& _Left, const _Ty& _Right)
{	// return smaller of _Left and _Right
    return (_Right < _Left ? _Right : _Left);
}

template<typename _Ty> inline
const _Ty& Max(const _Ty& _Left, const _Ty& _Right)
{	// return larger of _Left and _Right
    return (_Left < _Right ? _Right : _Left);
}

//////////////////////////////////////////////////////////////////////////////

template <typename T1, typename T2>
struct _IsSameType
{
	enum {VALUE = false};
	typedef False Type;
};

template <typename T>
struct _IsSameType<T, T>
{
	enum {VALUE = true};
	typedef True Type;
};

template <typename T1, typename T2>
inline bool 
_isSameType()
{
	return _IsSameType<T1, T2>::VALUE;
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...


