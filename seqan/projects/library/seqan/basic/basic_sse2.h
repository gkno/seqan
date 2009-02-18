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
  $Id: basic_operator.h 953 2007-07-27 11:48:23Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_SSE2_H
#define SEQAN_HEADER_BASIC_SSE2_H

//#ifdef __SSE2__
#if 1
#include <emmintrin.h>

namespace SEQAN_NAMESPACE_MAIN
{
	
	#define __int128 SSE2_int128

	// may become obsolete when 128-bit integers will be introduced
	struct SSE2_int128
	{
	public:
		__m128i                 data;
		static const __m128i    overflow;

//____________________________________________________________________________

	public:
		SSE2_int128();
		SSE2_int128(SSE2_int128 const &);
		SSE2_int128(__m128i);
		SSE2_int128(__int64, __int64);
		SSE2_int128(int, int, int, int);
		SSE2_int128(short, short, short, short, short, short, short, short);

		SSE2_int128(int);

		inline operator __m128i () { return data; }
	};

	const __m128i SSE2_int128::overflow = SSE2_int128(0, 1, 0, 0);

//////////////////////////////////////////////////////////////////////////////

inline void
clear(SSE2_int128 &me)
{
	me.data = ::_mm_setzero_si128();
}

//////////////////////////////////////////////////////////////////////////////

inline void
assign(SSE2_int128 &me, SSE2_int128 const &source)
{
	me.data = source.data;
}

inline void
assign(SSE2_int128 &me, __m128i source)
{
	me.data = source;
}

//////////////////////////////////////////////////////////////////////////////

inline SSE2_int128::SSE2_int128()
{
	clear(*this);
}
	
inline SSE2_int128::SSE2_int128(SSE2_int128 const &other)
{
	assign(*this, other);
}

inline SSE2_int128::SSE2_int128(__m128i other)
{
	assign(*this, other);
}

// 2x 64bit
inline SSE2_int128::SSE2_int128(__int64 q1, __int64 q0)
{
	data = ::_mm_set_epi32(q1 >> 32, q1, q0 >> 32, q0);
}

// 4x 32bit
inline SSE2_int128::SSE2_int128(int q3, int q2, int q1, int q0)
{
	data = ::_mm_set_epi32(q3, q2, q1, q0);
}

// 8x 16bit
inline SSE2_int128::SSE2_int128(
				   short q7, short q6, short q5, short q4,
				   short q3, short q2, short q1, short q0)
{
	data = ::_mm_set_epi16(q7, q6, q5, q4, q3, q2, q1, q0);
}

// 1x 32bit
inline SSE2_int128::SSE2_int128(int q)
{
	data = ::_mm_set_epi32(0, 0, 0, q);
}

//____________________________________________________________________________
// logical operators

inline SSE2_int128
operator & (SSE2_int128 const &a, SSE2_int128 const &b)
{
	return ::_mm_and_si128(a.data, b.data);
}

inline SSE2_int128
operator | (SSE2_int128 const &a, SSE2_int128 const &b)
{
	return ::_mm_or_si128(a.data, b.data);
}

inline SSE2_int128
operator ^ (SSE2_int128 const &a, SSE2_int128 const &b)
{
	return ::_mm_xor_si128(a.data, b.data);
}

//____________________________________________________________________________
// shift operators

inline SSE2_int128
operator << (SSE2_int128 const &a, const int n)
{
	if (n < 64)
		return ::_mm_or_si128(
			_mm_sll_epi64(a.data, (SSE2_int128)n),
			_mm_srl_epi64(
				_mm_unpacklo_epi64(::_mm_setzero_si128(), a.data),
				(SSE2_int128)(64-n)));
	else
		return _mm_sll_epi64(
			_mm_unpacklo_epi64(::_mm_setzero_si128(), a.data),
			(SSE2_int128)(n-64));
}

inline SSE2_int128
operator >> (SSE2_int128 const &a, const int n)
{
	if (n < 64)
		return ::_mm_or_si128(
			_mm_srl_epi64(a.data, (SSE2_int128)n),
			_mm_sll_epi64(
				_mm_unpackhi_epi64(a.data, ::_mm_setzero_si128()),
				(SSE2_int128)(64-n)));
	else
		return _mm_srl_epi64(
			_mm_unpackhi_epi64(a.data, ::_mm_setzero_si128()),
			(SSE2_int128)(n-64));
}

//____________________________________________________________________________
// artihmetic operators

inline SSE2_int128
operator + (SSE2_int128 const &a, SSE2_int128 const &b)
{
	union {
		__uint64 a[2];
		__m128i v;
	} _sum, _a;
	
	_a.v = a.data;
	_sum.v = ::_mm_add_epi64(a.data, b.data);
	if (_sum.a[0] >= _a.a[0])
		return _sum.v;
	else
		return ::_mm_add_epi64(_sum.v, SSE2_int128::overflow);
}

inline SSE2_int128
operator - (SSE2_int128 const &a, SSE2_int128 const &b)
{
	union {
		__uint64 a[2];
		__m128i v;
	} _diff, _a;
	
	_a.v = a.data;
	_diff.v = ::_mm_sub_epi64(a.data, b.data);
	if (_diff.a[0] <= _a.a[0])
		return _diff.v;
	else
		return ::_mm_sub_epi64(_diff.v, SSE2_int128::overflow);
}


} //namespace SEQAN_NAMESPACE_MAIN

#endif

#endif //#ifndef SEQAN_HEADER_...
