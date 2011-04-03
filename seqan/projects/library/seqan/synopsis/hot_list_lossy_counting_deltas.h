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
// Implementation of the Lossy Counting (with deltas) hot list algorithm.
// ==========================================================================

#ifndef SEQAN_SYNOPSIS_HOT_LIST_LOSSY_COUNTING_DELTAS_H_
#define SEQAN_SYNOPSIS_HOT_LIST_LOSSY_COUNTING_DELTAS_H_

#include <algorithm>
#include <tr1/unordered_map>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Spec.Lossy Counting Deltas HotList
..cat:Synopsis Data Structures
..general:Class.HotList
..summary:Approximate hot list using the Lossy Counting (with deltas) algorithm.
..signature:HotList<TValue, LossyCountingDeltas>
..param.TValue:Type of the items to keep in hot list.
..remarks:See Gurmeet Singh Manku and Rajeev Motwani.  Approximate Frequency Counts over Data Streams.  VLDB '02.
..include:seqan/synopsis.h

.Memfunc.Lossy Counting Deltas HotList#HotList
..class:Spec.Frequent HotList
..summary:Constructor
..signature:HotList(k)
..param.k:The number of most frequent items to keep.
 */

struct LossyCountingDeltas_;
typedef Tag<LossyCountingDeltas_> LossyCountingDeltas;

template <typename TValue, typename TCallback>
class HotList<TValue, LossyCountingDeltas, TCallback>
{
public:
    typedef HotList<TValue, Frequent> THotList_;
    typedef typename Size<THotList_>::Type TSize_;
    typedef typename std::tr1::unordered_map<TValue, TSize_> TMap_;

    TSize_ _n;   // number of items seen so far, unnecessary?
    unsigned _k;
    TMap_ _map;  // dictionary

    Holder<TCallback> _callbackContext;
    
    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    HotList(unsigned k)
            : _k(k)
    {
        create(_callbackContext);
    }

    HotList(unsigned k, TCallback & callbackContext)
            : _k(k), _callbackContext(callbackContext)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
clear(HotList<TValue, LossyCountingDeltas> & hotList)
{
}

// ----------------------------------------------------------------------------
// Function registerItem()
// ----------------------------------------------------------------------------

template <typename TValue>
inline typename Size<HotList<TValue, Frequent> >::Type
registerItem(HotList<TValue, LossyCountingDeltas> & hotList, TValue const & value)
{
}

// ----------------------------------------------------------------------------
// Function removeItem()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
removeItem(HotList<TValue, LossyCountingDeltas> & hotList, TValue const & value)
{
}

// ----------------------------------------------------------------------------
// Function getItems()
// ----------------------------------------------------------------------------

template <typename TResult, typename TValue>
inline void
getItems(TResult & items,
         HotList<TValue, LossyCountingDeltas> const & hotList)
{
}

}  // namespace seqan

#endif  // SEQAN_SYNOPSIS_HOT_LIST_LOSSY_COUNTING_DELTAS_H_
