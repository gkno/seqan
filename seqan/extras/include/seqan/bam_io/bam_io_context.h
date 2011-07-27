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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Class BamIoContext, accessor functions.
// ==========================================================================

// TODO(holtgrew): Test me!

#ifndef EXTRAS_INCLUDE_SEQAN_BAM_IO_BAM_IO_CONTEXT_H_
#define EXTRAS_INCLUDE_SEQAN_BAM_IO_BAM_IO_CONTEXT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TNameStore_, typename TNameStoreCache_ = NameStoreCache<TNameStore_> >
struct BamIOContext
{
    typedef TNameStore_ TNameStore;
    typedef TNameStoreCache_ TNameStoreCache;

    TNameStore * _nameStore;
    TNameStoreCache * _nameStoreCache;

    BamIOContext() : _nameStore(0), _nameStoreCache(0)
    {}

    BamIOContext(TNameStore & nameStore, TNameStoreCache & nameStoreCache) :
            _nameStore(&nameStore), _nameStoreCache(&nameStoreCache)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function nameStore()
// ----------------------------------------------------------------------------

template <typename TNameStore, typename TNameStoreCache>
TNameStore &
nameStore(BamIOContext<TNameStore, TNameStoreCache> & context)
{
    SEQAN_ASSERT(context._nameStore != 0);
    return *context._nameStore;
}

// ----------------------------------------------------------------------------
// Function nameStoreCache()
// ----------------------------------------------------------------------------

template <typename TNameStore, typename TNameStoreCache>
TNameStoreCache &
nameStoreCache(BamIOContext<TNameStore, TNameStoreCache> & context)
{
    SEQAN_ASSERT(context._nameStoreCache != 0);
    return *context._nameStoreCache;
}

}  // namespace seqan

#endif  // #ifndef EXTRAS_INCLUDE_SEQAN_BAM_IO_BAM_IO_CONTEXT_H_
