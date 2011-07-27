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
// BamHeaderRecord class, supporting types and functions for accessing tags
// in headers.
// ==========================================================================

// TODO(holtgrew): Test me!

#ifndef EXTRAS_INCLUDE_SEQAN_BAM_IO_BAM_HEADER_RECORD_H_
#define EXTRAS_INCLUDE_SEQAN_BAM_IO_BAM_HEADER_RECORD_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

enum BamHeaderRecordType
{
    BAM_HEADER_FIRST       = 0,
    BAM_HEADER_REFERENCE   = 1,
    BAM_HEADER_READ_GROUP  = 2,
    BAM_HEADER_PROGRAM     = 3,
    BAM_HEADER_COMMENT     = 4
};

struct BamHeaderRecord
{
    typedef CharString TTagName;
    typedef CharString TTagValue;
    typedef String<Pair<TTagName, TTagValue> > TTags;

    BamHeaderRecordType type;
    String<Pair<TTagName, TTagValue> > tags;

    BamHeaderRecord() {}
};

struct BamHeader
{
    typedef Pair<CharString, __int32> TSequenceInfo;
    
    String<Pair<CharString, __int32> > sequenceInfos;
    String<BamHeaderRecord> records;
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

inline void
clear(BamHeaderRecord & record)
{
    clear(record.tags);
}

// ----------------------------------------------------------------------------
// Function findKey()
// ----------------------------------------------------------------------------

inline bool
findKey(unsigned & idx, CharString const & key, BamHeaderRecord const & record)
{
    typedef BamHeaderRecord::TTags const TTags;
    typedef Iterator<TTags, Rooted>::Type TIterator;

    idx = 0;
    for (TIterator it = begin(record.tags, Rooted()); !atEnd(it); goNext(it), ++idx)
        if (it->i1 == key)
            return true;

    return false;
}

// ----------------------------------------------------------------------------
// Function getTagValue()
// ----------------------------------------------------------------------------

inline bool
getTagValue(CharString & key, unsigned idx, BamHeaderRecord const & record)
{
    if (idx >= length(record.tags))
        return false;
    key = record.tags[idx].i2;
    return true;
}

inline bool
getTagValue(CharString & value, CharString const & key, BamHeaderRecord const & record)
{
    unsigned idx = 0;
    if (!findKey(idx, key, record))
        return false;
    return getTagValue(value, key, record);
}

}  // namespace seqan

#endif  // #ifndef EXTRAS_INCLUDE_SEQAN_BAM_IO_BAM_HEADER_RECORD_H_
