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
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Alphabet concepts stemming from mathematics.
// ==========================================================================

#include <climits>
#include <float.h>

#ifndef SEQAN_BASIC_ALPHABET_MATH_H_
#define SEQAN_BASIC_ALPHABET_MATH_H_

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

// ----------------------------------------------------------------------------
// Metafunction MaxValue
// ----------------------------------------------------------------------------

template <typename T>
struct MaximumValueUnsigned_ { static const T VALUE; };
template <typename T>
struct MaximumValueSigned_ { static const T VALUE; };

template <typename T = void>
struct MaximumValueFloat_ { static const float VALUE; };
template <typename T = void>
struct MaximumValueDouble_ { static const double VALUE; };

template <typename T>
const T MaximumValueUnsigned_<T>::VALUE = ~(T)0;
template <typename T>
const T MaximumValueSigned_<T>::VALUE = ( (((T)1 <<(BitsPerValue<T>::VALUE - 2)) - 1) <<1) + 1;
template <typename T>
const float MaximumValueFloat_<T>::VALUE = FLT_MAX;
template <typename T>
const double MaximumValueDouble_<T>::VALUE = DBL_MAX;

template <>
struct MaximumValueSigned_<bool>
{
    static const bool VALUE = false;
};

// template <>
// const char MaximumValueUnsigned_<char>::VALUE = CHAR_MAX;
// template <>
// const char MaximumValueSigned_<char>::VALUE = CHAR_MAX;

template <
    typename T,
    typename TParent = typename If<
      IsSameType<double, T>::VALUE,
      MaximumValueDouble_<>,
      typename If<
      IsSameType<float, T>::VALUE,
      MaximumValueFloat_<>,
      typename If<
        IsSameType< typename MakeSigned_<T>::Type, T >::VALUE,
        MaximumValueSigned_<T>,
        MaximumValueUnsigned_<T>
        >::Type
      >::Type
    >::Type
  >
struct MaxValue : TParent {};

// ----------------------------------------------------------------------------
// Metafunction MinValue
// ----------------------------------------------------------------------------

template <typename T>
struct MinimumValueUnsigned_ {  static const T VALUE; };
template <typename T>
struct MinimumValueSigned_ { static const T VALUE; };

template <typename T = void>
struct MinimumValueFloat_ { static const float VALUE; };
template <typename T = void>
struct MinimumValueDouble_ { static const double VALUE; };

template <typename T>
const T MinimumValueUnsigned_<T>::VALUE = T(0);
template <typename T>
const T MinimumValueSigned_<T>::VALUE = ~(T)MaximumValueSigned_<T>::VALUE;
template <typename T>
const float MinimumValueFloat_<T>::VALUE = -FLT_MAX;
template <typename T>
const double MinimumValueDouble_<T>::VALUE = -DBL_MAX;

template <>
struct MinimumValueSigned_<bool>
{
    static const bool VALUE = false;
};

// template <>
// const char MinimumValueUnsigned_<char>::VALUE = 0;
// template <>
// const char MinimumValueSigned_<char>::VALUE = 0;

template <
    typename T,
    typename TParent = typename If<
      IsSameType<double, T>::VALUE,
      MinimumValueDouble_<>,
      typename If<
      IsSameType<float, T>::VALUE,
      MinimumValueFloat_<>,
      typename If<
        IsSameType< typename MakeSigned_<T>::Type, T >::VALUE,
        MinimumValueSigned_<T>,
        MinimumValueUnsigned_<T>
        >::Type
      >::Type
    >::Type
  >
struct MinValue : TParent {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function supremumValueImpl
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Function maxValue
// ----------------------------------------------------------------------------

template <typename T>
inline T const &
maxValue()
{
    SEQAN_CHECKPOINT;
    T * _tag = 0;
    return supremumValueImpl(_tag);
}

template <typename T>
inline T const &
maxValue(T)
{
    SEQAN_CHECKPOINT;
    T * _tag = 0;
    return supremumValueImpl(_tag);
}

// ----------------------------------------------------------------------------
// Function infimumValueImpl
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Function minValue
// ----------------------------------------------------------------------------

template <typename T>
inline T const &
minValue()
{
    SEQAN_CHECKPOINT;
    T * _tag = 0;
    return infimumValueImpl(_tag);
}

template <typename T>
inline T const &
minValue(T)
{
    SEQAN_CHECKPOINT;
    T * _tag = 0;
    return infimumValueImpl(_tag);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_ALPHABET_MATH_H_
