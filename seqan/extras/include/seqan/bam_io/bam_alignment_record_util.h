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

#ifndef EXTRAS_INCLUDE_SEQAN_BAM_IO_BAM_ALIGNMENT_RECORD_UTIL_H_
#define EXTRAS_INCLUDE_SEQAN_BAM_IO_BAM_ALIGNMENT_RECORD_UTIL_H_

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

// ============================================================================
// Functions
// ============================================================================

/**
.Function.bamRecordToAlignment
..signature:bamRecordToAlignment(align, reference, record)
..param.align:The alignment to create.
...type:Class.Align
..param.reference:String of Dna, Dna5, ... characters.
...type:Class.String
..param.record:The alignment record to convert.
...type:Class.BamAlignmentRecord
..returns:$void$
..include:seqan/bam_io.h
..example.code:
StringSet<Dna5String> references;
BamAlignment record;
// Read references and record.
Align<Dna5String> align;
if (record.rId != BamAlignmentRecord::INVALID_REFID)
    bamRecordToAlignment(align, references[record.refId], record);
 */

template <typename TSource, typename TSpec, typename TReference>
void
bamRecordToAlignment(Align<TSource, TSpec> & result, TReference & reference, BamAlignmentRecord & record)
{
    // TODO(holtgrew): Clipping better than copying infix? But is it generic?
    resize(rows(result), 2);

    setSource(row(result, 0), reference, record.pos, record.pos + getAlignmentLengthInRef(record));
    cigarToGapAnchorContig(record.cigar, row(result, 0));
    assignSource(row(result, 1), record.seq);
    cigarToGapAnchorRead(record.cigar, row(result, 1));
}

}  // namespace seqan

#endif  // #ifndef EXTRAS_INCLUDE_SEQAN_BAM_IO_BAM_ALIGNMENT_RECORD_UTIL_H_
