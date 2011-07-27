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
// The class BamAlignmentRecord, flag checking methods, flag constants.
// ==========================================================================

// TODO(holtgrew): Test me!

#ifndef EXTRAS_INCLUDE_SEQAN_BAM_IO_BAM_RECORD_H_
#define EXTRAS_INCLUDE_SEQAN_BAM_IO_BAM_RECORD_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

enum BAM_FLAGS_
{
    BAM_FLAG_MULTIPLE      = 0x0001,
    BAM_FLAG_ALL_PROPER    = 0x0002,
    BAM_FLAG_UNMAPPED      = 0x0004,
    BAM_FLAG_NEXT_UNMAPPED = 0x0008,
    BAM_FLAG_RC            = 0x0010,
    BAM_FLAG_NEXT_RC       = 0x0020,
    BAM_FLAG_FIRST         = 0x0040,
    BAM_FLAG_LAST          = 0x0080,
    BAM_FLAG_SECONDARY     = 0x0100,
    BAM_FLAG_QC_NO_PASS    = 0x0200,
    BAM_FLAG_DUPLICATE     = 0x0400
};

class BamAlignmentRecord
{
public:
    static __uint32 const INVALID_POS = MaxValue<__uint32>::VALUE;
    static __int32 const INVALID_REFID = -1;
    static __int32 const INVALID_LEN = MaxValue<__int32>::VALUE;
    
    CharString qName;
    __uint16 flag;
    __int32 rId;
    __uint32 pos;
    __uint8 mapQ;
    __uint16 bin;
    String<CigarElement<> > cigar;
    __int32 rNextId;
    __uint32 pNext;
    __int32 tLen;
    CharString seq;
    CharString qual;
    CharString tags;  // raw tags in BAM format
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
clear(BamAlignmentRecord & record)
{
    clear(record.qName);
    clear(record.cigar);
    clear(record.seq);
    clear(record.qual);
    clear(record.tags);
}

// ----------------------------------------------------------------------------
// Function hasFlagMultiple()
// ----------------------------------------------------------------------------

inline bool
hasFlagMultiple(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_MULTIPLE) == BAM_FLAG_MULTIPLE;
}

// ----------------------------------------------------------------------------
// Function hasFlagAllProper()
// ----------------------------------------------------------------------------

inline bool
hasFlagAllProper(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_ALL_PROPER) == BAM_FLAG_ALL_PROPER;
}

// ----------------------------------------------------------------------------
// Function hasFlagUnmapped()
// ----------------------------------------------------------------------------

inline bool
hasFlagUnmapped(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_UNMAPPED) == BAM_FLAG_UNMAPPED;
}

// ----------------------------------------------------------------------------
// Function hasFlagNextUnmapped()
// ----------------------------------------------------------------------------

inline bool
hasFlagNextUnmapped(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_NEXT_UNMAPPED) == BAM_FLAG_NEXT_UNMAPPED;
}

// ----------------------------------------------------------------------------
// Function hasFlagRC()
// ----------------------------------------------------------------------------

inline bool
hasFlagRC(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_RC) == BAM_FLAG_RC;
}

// ----------------------------------------------------------------------------
// Function hasFlagNextRC()
// ----------------------------------------------------------------------------

inline bool
hasFlagNextRC(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_NEXT_RC) == BAM_FLAG_NEXT_RC;
}

// ----------------------------------------------------------------------------
// Function hasFlagFirst()
// ----------------------------------------------------------------------------

inline bool
hasFlagFirst(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_FIRST) == BAM_FLAG_FIRST;
}

// ----------------------------------------------------------------------------
// Function hasFlagLast()
// ----------------------------------------------------------------------------

inline bool
hasFlagLast(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_LAST) == BAM_FLAG_LAST;
}

// ----------------------------------------------------------------------------
// Function hasFlagSecondary()
// ----------------------------------------------------------------------------

inline bool
hasFlagSecondary(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_SECONDARY) == BAM_FLAG_SECONDARY;
}

// ----------------------------------------------------------------------------
// Function hasFlagQcNoPass()
// ----------------------------------------------------------------------------

inline bool
hasFlagQcNoPass(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_QC_NO_PASS) == BAM_FLAG_QC_NO_PASS;
}

// ----------------------------------------------------------------------------
// Function hasFlagDuplicate()
// ----------------------------------------------------------------------------

inline bool
hasFlagDuplicate(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_DUPLICATE) == BAM_FLAG_DUPLICATE;
}

}  // namespace seqan

#endif  // #ifndef EXTRAS_INCLUDE_SEQAN_BAM_IO_BAM_RECORD_H_
