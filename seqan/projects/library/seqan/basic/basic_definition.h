 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id$
 ==========================================================================*/

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

template <typename T>
struct Length;

template <>
struct Length<void>
{
	enum { VALUE = 0 };
};

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.TagList:
..summary:A structure to represent a list of tags.
..signature:TagList<TTag1>
..signature:TagList<TTag1, TagList<TTag2> >
..signature:TagList<TTag1, TagList<TTag2, TagList<TTag3[...]> > >
..param.TTag1:The first tag of the list.
..param.TTag2:The second tag of the list.
..param.TTag3:The third tag of the list.
..include:seqan/basic.h
*/

template <typename TTag = void, typename TSubList = void>
struct TagList
{
	typedef TTag Type;
};

template <typename TTag>
struct Length< TagList<TTag, void> > {
	enum { VALUE = 1 };
};

template <typename TTag, typename TSubList>
struct Length< TagList<TTag, TSubList> > {
	enum { VALUE = Length<TSubList>::VALUE + 1 };
};

template <typename TTagList = void>
struct TagSelector
{
	int tagId;
	
	TagSelector()
	{
		tagId = 0;
	}
};

/**
.Class.TagSelector:
..summary:A structure to select a tag from a @Tag.TagList@.
..signature:TagSelector<TTagList>
..param.TTagList:A tag list.
...type:Tag.TagList
.Memvar.TagSelector#tagId:
..class:Class.TagSelector
..type:nolink:int
..summary:Stores the index of a @Page.Glossary.Tag@ in the tag list.
..include:seqan/basic.h
*/

///

template <typename TTag, typename TSubList>
struct TagSelector< TagList<TTag, TSubList> >:
	TagSelector<TSubList>
{
	typedef TTag					Type;
	typedef TagSelector<TSubList>	Base;
};

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Default:
..summary:Tag that specifies default behavior.
..tag.Default:Use default behavior. 
..include:seqan/basic.h
*/
struct Default_;
typedef Tag<Default_> const Default;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Move Switch:
..summary:Switch to force move.
..tag.Move:Move instead of assign. 
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
..include:seqan/basic.h
*/

struct Move_;
typedef Tag<Move_> const Move;

//////////////////////////////////////////////////////////////////////////////

//Pass to c'tor of iterator to move it to the end
struct GoEnd_;
typedef Tag<GoEnd_> const GoEnd;


//////////////////////////////////////////////////////////////////////////////

//construct without initializing
struct MinimalCtor_;
typedef Tag<MinimalCtor_> const MinimalCtor;

//construct with initializing
struct NonMinimalCtor_;
typedef Tag<NonMinimalCtor_> const NonMinimalCtor;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Nothing:
..summary:Tag that represents an absent parameter or an absent type.
..tag.Nothing:Omit parameter.
..include:seqan/basic.h
*/
///Empty Data Class.
struct Nothing {};



//////////////////////////////////////////////////////////////////////////////
// returns TTo const, if TFrom is const, TTo otherwise

template <typename TFrom, typename TTo>
struct _CopyConst
{
	typedef TTo Type;
};
template <typename TFrom, typename TTo>
struct _CopyConst<TFrom const, TTo>
{
	typedef TTo const Type;
};

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
struct _RemoveConst<T const>:
	public _RemoveConst<T> {};

template <typename T>
struct _RemoveConst<T &>
{
	typedef typename _RemoveConst<T>::Type & Type;
};
template <typename T>
struct _RemoveConst<T *>
{
	typedef typename _RemoveConst<T>::Type * Type;
};
template <typename T, size_t I>
struct _RemoveConst<T const [I]>
{
	typedef T * Type;
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
	typedef
		typename IF< TYPECMP<T, char>::VALUE,         unsigned char,
		typename IF< TYPECMP<T, signed char>::VALUE,  unsigned char,
		typename IF< TYPECMP<T, signed short>::VALUE, unsigned short,
		typename IF< TYPECMP<T, signed int>::VALUE,   unsigned int,
		typename IF< TYPECMP<T, signed long>::VALUE,  unsigned long,
		typename IF< TYPECMP<T, __int64>::VALUE,      __uint64, T
		>::Type>::Type>::Type>::Type>::Type>::Type Type;
};

template <typename T>
struct _MakeUnsigned<T const> {
	typedef typename _MakeUnsigned<T>::Type const Type;
};

//////////////////////////////////////////////////////////////////////////////
/**
.Internal._MakeSigned:
..signature:_MakeSigned<T>
..returns:$signed t$ if $T$ is not $signed t$, otherwise $T$.
*/
template <typename T>
struct _MakeSigned
{
	typedef
		typename IF< TYPECMP<T, char>::VALUE,           signed char,
		typename IF< TYPECMP<T, unsigned char>::VALUE,  signed char,
		typename IF< TYPECMP<T, unsigned short>::VALUE, signed short,
		typename IF< TYPECMP<T, unsigned int>::VALUE,   signed int,
		typename IF< TYPECMP<T, unsigned long>::VALUE,  signed long,
		typename IF< TYPECMP<T, __uint64>::VALUE,       __int64, T
		>::Type>::Type>::Type>::Type>::Type>::Type Type;
};

template <typename T>
struct _MakeSigned<T const> {
	typedef typename _MakeSigned<T>::Type const Type;
};

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
	log2(T /*val*/, unsigned int offset)
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

template <typename TValue, typename TExponent>
inline TValue _intPow(TValue a, TExponent b)
{
SEQAN_CHECKPOINT
	TValue ret = 1;
	while (b != 0)
	{
		if (b & 1) ret *= a;
		a *= a;
		b >>= 1;
	}	
	return ret;
}

//////////////////////////////////////////////////////////////////////////////
// to avoid conflicts with non-standard macros and namespaces
// we define our own Min/Max functions

template<typename _Tx> inline
const _Tx& _min(const _Tx& _Left, const _Tx& _Right)
{	// return smaller of _Left and _Right
	if (_Left < _Right)
		return _Left;
	else
		return _Right;
}

template<typename _Tx, typename _Ty> inline
_Tx _min(const _Tx& _Left, const _Ty& _Right)
{	// return smaller of _Left and _Right
    return (_Right < _Left ? _Right : _Left);
}

template<typename _Ty> inline
const _Ty& _max(const _Ty& _Left, const _Ty& _Right)
{	// return larger of _Left and _Right
	if (_Left < _Right)
		return _Right;
	else
		return _Left;
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

