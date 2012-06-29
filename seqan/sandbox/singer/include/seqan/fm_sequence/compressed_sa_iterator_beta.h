// ==========================================================================
//                                  FMIndex
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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

#ifndef COMPRESSED_SA_ITERATOR_BETA_H_
#define COMPRESSED_SA_ITERATOR_BETA_H_

namespace seqan {

// ==========================================================================
//Forwards
// ==========================================================================

struct FibreSparseString_;
typedef Tag<FibreSparseString_> const FibreSparseString;    

template <typename TSparseString, typename TLfTable, typename TSpec>
class CompressedSA;

// ==========================================================================
//Metafunctions
// ==========================================================================

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Standard>
{
    typedef Iter<CompressedSA<TSparseString, TLfTable, TSpec> const, PositionIterator> Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Standard>
{
    typedef Iter<CompressedSA<TSparseString, TLfTable, TSpec>, PositionIterator> Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Rooted>:
    Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Standard>{};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Rooted>:
    Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Standard>{};

template <typename TSparseString, typename TLfTable, typename TSpec>
inline Iterator<CompressedSA<TSparseString, TLfTable, TSpec> >
begin(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA)
{
    return Iterator<CompressedSA<TSparseString, TLfTable, TSpec> >(compressedSA, 0);
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline Iterator<CompressedSA<TSparseString, TLfTable, TSpec> >
end(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA)
{
    return Iterator<CompressedSA<TSparseString, TLfTable, TSpec> >(compressedSA, length(compressedSA.compressedSA));
}

}
#endif // COMPRESSED_SA_ITERATOR_BETA_H_

