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

inline __int16  atomicInc(__int16  volatile & x) { return _InterlockedIncrement16(&x); }
inline __int32  atomicInc(__int32  volatile & x) { return _InterlockedIncrement(&x); }
inline __int64  atomicInc(__int64  volatile & x) { return _InterlockedIncrement64(&x); }
inline __uint16 atomicInc(__uint16 volatile & x) { return _InterlockedIncrement16(&x); }
inline __uint32 atomicInc(__uint32 volatile & x) { return _InterlockedIncrement(&x); }
inline __uint64 atomicInc(__uint64 volatile & x) { return _InterlockedIncrement64(&x); }

inline __int16  atomicDec(__int16  volatile & x) { return _InterlockedDecrement16(&x); }
inline __int32  atomicDec(__int32  volatile & x) { return _InterlockedDecrement(&x); }
inline __int64  atomicDec(__int64  volatile & x) { return _InterlockedDecrement64(&x); }
inline __uint16 atomicDec(__uint16 volatile & x) { return _InterlockedDecrement16(&x); }
inline __uint32 atomicDec(__uint32 volatile & x) { return _InterlockedDecrement(&x); }
inline __uint64 atomicDec(__uint64 volatile & x) { return _InterlockedDecrement64(&x); }

inline __int32  atomicAdd(__int32  volatile & x, __int32  y) { return _InterlockedExchangeAdd(&x, y); }
inline __int64  atomicAdd(__int64  volatile & x, __int64  y) { return _InterlockedExchangeAdd64(&x, y); }
inline __uint32 atomicAdd(__uint32 volatile & x, __uint32 y) { return _InterlockedExchangeAdd(&x, y); }
inline __uint64 atomicAdd(__uint64 volatile & x, __uint64 y) { return _InterlockedExchangeAdd64(&x, y); }

inline __int8   atomicOr(__int8   volatile & x, __int8   y) { return _InterlockedOr8(x, y); }
inline __int16  atomicOr(__int16  volatile & x, __int16  y) { return _InterlockedOr16(x, y); }
inline __int32  atomicOr(__int32  volatile & x, __int32  y) { return _InterlockedOr(x, y); }
inline __int64  atomicOr(__int64  volatile & x, __int64  y) { return _InterlockedOr64(x, y); }
inline __uint8  atomicOr(__uint8  volatile & x, __uint8  y) { return _InterlockedOr8(x, y); }
inline __uint16 atomicOr(__uint16 volatile & x, __uint16 y) { return _InterlockedOr16(x, y); }
inline __uint32 atomicOr(__uint32 volatile & x, __uint32 y) { return _InterlockedOr(x, y); }
inline __uint64 atomicOr(__uint64 volatile & x, __uint64 y) { return _InterlockedOr64(x, y); }

inline __int8   atomicXor(__int8   volatile & x, __int8   y) { return _InterlockedXor8(x, y); }
inline __int16  atomicXor(__int16  volatile & x, __int16  y) { return _InterlockedXor16(x, y); }
inline __int32  atomicXor(__int32  volatile & x, __int32  y) { return _InterlockedXor(x, y); }
inline __int64  atomicXor(__int64  volatile & x, __int64  y) { return _InterlockedXor64(x, y); }
inline __uint8  atomicXor(__uint8  volatile & x, __uint8  y) { return _InterlockedXor8(x, y); }
inline __uint16 atomicXor(__uint16 volatile & x, __uint16 y) { return _InterlockedXor16(x, y); }
inline __uint32 atomicXor(__uint32 volatile & x, __uint32 y) { return _InterlockedXor(x, y); }
inline __uint64 atomicXor(__uint64 volatile & x, __uint64 y) { return _InterlockedXor64(x, y); }

inline __int16  atomicCas(__int16  volatile & x, __int16  cmp, __int16  y) { return _InterlockedCompareExchange16(&x, cmp, y); }
inline __int32  atomicCas(__int32  volatile & x, __int32  cmp, __int32  y) { return _InterlockedCompareExchange(&x, cmp, y); }
inline __int64  atomicCas(__int64  volatile & x, __int64  cmp, __int64  y) { return _InterlockedCompareExchange64(&x, cmp, y); }
inline __uint16 atomicCas(__uint16 volatile & x, __uint16 cmp, __uint16 y) { return _InterlockedCompareExchange16(&x, cmp, y); }
inline __uint32 atomicCas(__uint32 volatile & x, __uint32 cmp, __uint32 y) { return _InterlockedCompareExchange(&x, cmp, y); }
inline __uint64 atomicCas(__uint64 volatile & x, __uint64 cmp, __uint64 y) { return _InterlockedCompareExchange64(&x, cmp, y); }

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
