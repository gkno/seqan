// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
// Author: Stephan Aiche <stephan.aiche@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_ARG_PARSE_INCLUDE_ARG_PARSE_ARG_PARSE_TYPE_SUPPRT_H_
#define SANDBOX_ARG_PARSE_INCLUDE_ARG_PARSE_ARG_PARSE_TYPE_SUPPRT_H_

namespace seqan {


// ----------------------------------------------------------------------------
// Function _canCast()
// ----------------------------------------------------------------------------
template <typename TTarget, typename TString>
inline TTarget
_canCast(TTarget & dest, TString const s)
{
    std::istringstream stream(toCString(s));
    bool result = (!(stream >> dest).fail()) && (stream.rdbuf()->in_avail() == 0);
    return result;
}

// ----------------------------------------------------------------------------
// Function _cast()
// ----------------------------------------------------------------------------
template <typename TTarget, typename TString>
inline TTarget
_cast(TString const s)
{
    TTarget dst;
    std::istringstream stream(toCString(s));
    bool result = (!(stream >> dst).fail()) && (stream.rdbuf()->in_avail() == 0);
    SEQAN_ASSERT_MSG(result, "could not cast %s", toCString(s));
    return dst;
}

// ----------------------------------------------------------------------------
// Function _isCastable()
// ----------------------------------------------------------------------------
template <typename TTarget, typename TString>
inline bool
_isCastable(TString const s)
{
    TTarget dst;
    std::istringstream stream(toCString(s));
    return (!(stream >> dst).fail()) && (stream.rdbuf()->in_avail() == 0);
}


// ----------------------------------------------------------------------------
// Function _isDouble()
// ----------------------------------------------------------------------------

template <typename TString>
inline bool
_isDouble(TString const s)
{
    return _isCastable<double>(s);
}

// ----------------------------------------------------------------------------
// Function _isInt()
// ----------------------------------------------------------------------------

template <typename TString>
inline bool
_isInt(TString const s)
{
    return _isCastable<int>(s);
}

} // namespace seqan

#endif // SANDBOX_ARG_PARSE_INCLUDE_ARG_PARSE_ARG_PARSE_TYPE_SUPPRT_H_
