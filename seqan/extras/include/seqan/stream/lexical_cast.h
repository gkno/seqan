// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// casts for reading different types from strings
// ==========================================================================+

#ifndef SEQAN_STREAM_LEXICAL_CAST_H
#define SEQAN_STREAM_LEXICAL_CAST_H


namespace seqan {

// template < typename TTarget, typename TSource >
// inline TTarget*
// lexicalCast(TSource *source)
// {
// 
//     std::istringstream str(*source);
// 
//     TTarget* ret = new TTarget();
// 
//     bool success = (str >> *ret);
// 
//     if (success)
//         return ret;
// 
//     delete ret;
//     return 0;
// }


template < typename TTarget, typename TSource >
inline TTarget
lexicalCast(TSource const & source)
{
    std::istringstream str(source);
    TTarget ret;

    str >> ret;
    return ret;
}

template < typename TTarget, typename TValue, typename TSpec>
inline TTarget
lexicalCast(String<TValue, TSpec> const & source)
{
    std::istringstream str(toCString(source));
    TTarget ret;

    str >> ret;
    return ret;
}



template < typename TTarget, typename TSource >
inline bool
lexicalCast2(TTarget & target, TSource const & source)
{
    std::istringstream str(source);
    return (str >> target);
}

template < typename TTarget, typename TValue, typename TSpec>
inline bool
lexicalCast2(TTarget & target, String<TValue, TSpec> const & source)
{
    std::istringstream str(toCString(source));
    return (str >> target);
}


}

#endif //def SEQAN_STREAM_LEXICAL_CAST_H
