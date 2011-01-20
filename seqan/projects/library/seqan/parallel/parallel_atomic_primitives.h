// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS: No forwards are generated for this file.

#ifndef SEQAN_PARALLEL_PARALLEL_ATOMIC_PRIMITIVES_H_
#define SEQAN_PARALLEL_PARALLEL_ATOMIC_PRIMITIVES_H_

#ifdef PLATFORM_WINDOWS
#include <intrin.h>
#endif  // #ifdef PLATFORM_WINDOWS

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// TODO(holtgrew): One could implement an HasAtomicFunction<T, Function> Metafunction but that would probably be pretty useless.

// ============================================================================
// Functions
// ============================================================================

// TODO(holtgrew): What about correct alignment?!

/**
.Function.atomicInc
..summary:Atomically increment an integer.
..cat:Atomic Operations
..signature:atomicInc(x)
..param.x:Integer, by reference.
..returns:The old value of $x.
..remarks:This is equivalent to an atomic $++x$.
..see:Function.atomicDec
..see:Function.atomicAdd
..see:Function.atomicOr
..see:Function.atomicXor
..see:Function.atomicCas
..include:seqan/parallel.h

.Function.atomicDec
..summary:Atomically decrement an integer.
..cat:Atomic Operations
..signature:atomicDec(x)
..param.x:Integer, by reference.
..returns:The old value of $x$.
..remarks:This is equivalent to an atomic $--x$.
..see:Function.atomicInc
..see:Function.atomicAdd
..see:Function.atomicOr
..see:Function.atomicXor
..see:Function.atomicCas
..include:seqan/parallel.h

.Function.atomicAdd
..summary:Atomically add an integer to another integer.
..cat:Atomic Operations
..signature:atomicAdd(x, y)
..param.x:Integer, by reference.
..param.y:Integer to add to the given value.
..returns:The old value of $x$.
..remarks:This is equivalent to an atomic $x += y$.
..see:Function.atomicInc
..see:Function.atomicDec
..see:Function.atomicOr
..see:Function.atomicXor
..see:Function.atomicCas
..include:seqan/parallel.h

.Function.atomicOr
..summary:Atomically combine two integers with $OR$ operation.
..cat:Atomic Operations
..signature:atomicOr(x, y)
..param.x:Integer, by reference.
..param.y:Integer to combine with $OR$ operation.
..returns:The old value of $x$.
..remarks:This is equivalent to an atomic $x |= y$.
..see:Function.atomicInc
..see:Function.atomicDec
..see:Function.atomicAdd
..see:Function.atomicXor
..see:Function.atomicCas
..include:seqan/parallel.h

.Function.atomicXor
..summary:Atomically combine wto integers with $XOR$ operation.
..cat:Atomic Operations
..signature:atomicXor(x, y)
..param.x:x, by reference.
..param.y:Integer to combine with $XOR$ operation.
..returns:The old value of $x$.
..remarks:This is equivalent to an atomic $x ^= y$.
..see:Function.atomicInc
..see:Function.atomicDec
..see:Function.atomicAdd
..see:Function.atomicOr
..see:Function.atomicCas
..include:seqan/parallel.h

.Function.atomicCas
..summary:Compare-and-Swap operation.
..cat:Atomic Operations
..signature:atomicCas(x, cmp, y)
..param.x:Pointer to the integer to swap.
..param.cmp:Value to compare $x$ with.
..param.y:Value to set $x$ to if it is equal to $cmp$.
..remarks:The pseudo code for this is
...code:
atomic
{
    T val = *(&x);
    if (val == cmp)
        *(&x) = y;
    return val;
}
..see:Function.atomicInc
..see:Function.atomicDec
..see:Function.atomicAdd
..see:Function.atomicOr
..see:Function.atomicXor
..include:seqan/parallel.h
 */

#ifdef PLATFORM_WINDOWS

// ----------------------------------------------------------------------------
// Implementation in MSVC
// ----------------------------------------------------------------------------

// We break the standard code layout here since we only wrap compiler
// intrinsics and it's easier to see things with one glance this way.

inline          long atomicInc(         long volatile & x) { return InterlockedIncrement(&x); }
inline unsigned long atomicInc(unsigned long volatile & x) { return InterlockedIncrement(reinterpret_cast<long volatile *>(&x)); }
#ifdef _WIN64
inline       __int64 atomicInc( __int64 volatile & x) { return _InterlockedIncrement64(&x); }
inline      __uint64 atomicInc(__uint64 volatile & x) { return _InterlockedIncrement64(reinterpret_cast<__int64 volatile *>(&x)); }
#endif  // #ifdef _WIN64

inline          long atomicDec(         long volatile & x) { return _InterlockedDecrement(&x); }
inline unsigned long atomicDec(unsigned long volatile & x) { return _InterlockedDecrement(reinterpret_cast<long volatile *>(&x)); }
#ifdef _WIN64
inline       __int64 atomicDec( __int64 volatile & x) { return _InterlockedDecrement(&x); }
inline      __uint64 atomicInc(__uint64 volatile & x) { return _InterlockedDecrement(reinterpret_cast<__int64 volatile *>(&x)); }
#endif  // #ifdef _WIN64

inline          long atomicAdd(         long volatile & x, long y) { return _InterlockedExchangeAdd(&x, y); }
inline unsigned long atomicAdd(unsigned long volatile & x, long y) { return _InterlockedExchangeAdd(reinterpret_cast<long volatile *>(&x), y); }
#ifdef _WIN64
inline       __int64 atomicAdd( __int64 volatile & x,  _uint64 y) { return _InterlockedExchangeAdd(&x, y); }
inline      __uint64 atomicAdd(__uint64 volatile & x, __uint64 y) { return _InterlockedExchangeAdd(reinterpret_cast<long volatile *>(&x), y); }
#endif  // #ifdef _WIN64

inline           char atomicOr(          char volatile & x,           char y) { return _InterlockedOr8(&x, y); }
inline unsigned  char atomicOr(unsigned  char volatile & x, unsigned  char y) { return _InterlockedOr8(reinterpret_cast<char volatile *>(&x), y); }
inline          short atomicOr(         short volatile & x,          short y) { return _InterlockedOr16(&x, y); }
inline unsigned short atomicOr(unsigned short volatile & x, unsigned short y) { return _InterlockedOr16(reinterpret_cast<short volatile *>(&x), y); }
inline           long atomicOr(          long volatile & x,           long y) { return _InterlockedOr(&x, y); }
inline unsigned  long atomicOr(unsigned  long volatile & x, unsigned  long y) { return _InterlockedOr(reinterpret_cast<long volatile *>(&x), y); }
#ifdef _WIN64
inline        __int64 atomicOr(       __int64 volatile & x,        __int64 y) { return _InterlockedOr64(&x, y); }
inline       __uint64 atomicOr(      __uint64 volatile & x,       __uint64 y) { return _InterlockedOr64(reinterpret_cast<__int64 volatile*>(&x), y); }
#endif  // #ifdef _WIN64

inline           char atomicXor(          char volatile & x,           char y) { return _InterlockedXor8(&x, y); }
inline unsigned  char atomicXor(unsigned  char volatile & x, unsigned  char y) { return _InterlockedXor8(reinterpret_cast<char volatile *>(&x), y); }
inline          short atomicXor(         short volatile & x,          short y) { return _InterlockedXor16(&x, y); }
inline unsigned short atomicXor(unsigned short volatile & x, unsigned short y) { return _InterlockedXor16(reinterpret_cast<short volatile *>(&x), y); }
inline           long atomicXor(          long volatile & x,           long y) { return _InterlockedXor(&x, y); }
inline unsigned  long atomicXor(unsigned  long volatile & x, unsigned  long y) { return _InterlockedXor(reinterpret_cast<long volatile *>(&x), y); }
#ifdef _WIN64
inline        __int64 atomicXor(       __int64 volatile & x,        __int64 y) { return _InterlockedXor64(&x, y); }
inline       __uint64 atomicXor(      __uint64 volatile & x,       __uint64 y) { return _InterlockedXor64(reinterpret_cast<__int64 volatile*>(&x), y); }
#endif  // #ifdef _WIN64

inline          short atomicCas(         short volatile & x,          short cmp,          short y) { return _InterlockedCompareExchange16(&x, y, cmp); }
inline unsigned short atomicCas(unsigned short volatile & x, unsigned short cmp, unsigned short y) { return _InterlockedCompareExchange16(reinterpret_cast<short volatile *>(&x), y, cmp); }
inline           long atomicCas(          long volatile & x,           long cmp,           long y) { return _InterlockedCompareExchange(&x, y, cmp); }
inline unsigned  long atomicCas(unsigned  long volatile & x, unsigned  long cmp, unsigned  long y) { return _InterlockedCompareExchange(reinterpret_cast<long volatile *>(&x), y, cmp); }
#ifdef _WIN64
inline        __int64 atomicCas(       __int64 volatile & x,        __int64 cmp,        __int64 y) { return _InterlockedCompareExchange64(&x, y, cmp); }
inline       __uint64 atomicCas(      __uint64 volatile & x,       __uint64 cmp,       __uint64 y) { return _InterlockedCompareExchange64(reinterpret_cast<__int64 volatile *>(&x), y, cmp); }
#endif  // #ifdef _WIN64

#else  // #ifdef PLATFORM_WINDOWS

// ----------------------------------------------------------------------------
// Implementation in GCC (LLVM is GCC compatible)
// ----------------------------------------------------------------------------

template <typename T>
inline T atomicInc(T volatile & x)
{
    return __sync_fetch_and_add(&x, 1);
}

template <typename T>
inline T atomicDec(T volatile & x)
{
    return __sync_fetch_and_add(&x, -1);
}

template <typename T1, typename T2>
inline T1 atomicAdd(T1 volatile & x, T2 y)
{
    return __sync_fetch_and_add(&x, y);
}

template <typename T>
inline T atomicOr(T volatile & x, T y)
{
    return __sync_fetch_and_or(&x, y);
}

template <typename T>
inline T atomicXor(T volatile & x, T y)
{
    return __sync_fetch_and_xor(&x, y);
}

template <typename T>
inline T atomicCas(T volatile & x, T cmp, T y)
{
    return __sync_val_compare_and_swap(&x, cmp, y);
}

#endif  // #ifdef PLATFORM_WINDOWS
	
} // namespace seqan

#endif  // #define SEQAN_PARALLEL_PARALLEL_ATOMIC_PRIMITIVES_H_
