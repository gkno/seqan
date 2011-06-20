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
// Code for writing FASTA and FASTQ files
// ==========================================================================


#ifndef SEQAN_STREAM_WRITE_H_
#define SEQAN_STREAM_WRITE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Enum.FastAQOutputOptions
..cat:Input / Output
..summary:Enum with flags for writing Fasta and Fastq files
..value.NO_FLAGS:No Flags are set.
..value.LINEBREAKS:Break lines at column 70
..include:seqan/stream.h
 */

enum FastAQOutputOptions {
    NO_FLAGS = 0,
    LINEBREAKS = 1,
    CHARBYCHAR = 2 //TODO(h4nn3s): remove once it is decided betweeen line-by-line and char-by-char
//  FOOFOO = 4,
//  BARBAR = 8
};

// ============================================================================
// Metafunctions
// ============================================================================


// ============================================================================
// Functions
// ============================================================================


// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------


/**
.Function.writeRecord
..signature:writeRecord(TStream & stream, TIdString const & meta, TSeqString const & seq, Fasta const &, FastAQOutputOptions const options)
..param.options:default is LINEBREAKS for FastA and NO_FLAG for FastQ
...type:Enum.FastAQOutputOptions
*/

// FASTA
template <typename TSeqString>
inline int
_writeRecordImplLineBreaks(FILE * stream,
            TSeqString const & seq,
            Fasta const & /*tag*/,
            unsigned const options)
{
    (void)options;
    int res = 0;
    // write stream character by character
    typename Iterator<TSeqString>::Type it = begin(seq);
    typename Iterator<TSeqString>::Type it_end = end(seq);
    for (unsigned long l = 0; it < it_end; ++it)
    {
        res = streamPut(stream, *it);
        if (res)
            return res;
        if (++l == 69)
        {
            res = streamPut(stream, '\n');
            l = 0;
            if (res)
                return res;
        }
    }
    return 0;
}

template <typename TStream,
          typename TSeqString>
inline int
_writeRecordImplLineBreaks(TStream & stream,
            TSeqString const & seq,
            Fasta const & /*tag*/,
            unsigned const options)
{
    int res = 0;
    if (options & CHARBYCHAR)
    {
        // write stream character by character
        typename Iterator<TSeqString>::Type it = begin(seq);
        typename Iterator<TSeqString>::Type it_end = end(seq);
        for (unsigned long l = 0; it < it_end; ++it)
        {
            res = streamPut(stream, *it);
            if (res)
                return res;
            if (++l == 69)
            {
                res = streamPut(stream, '\n');
                l = 0;
                if (res)
                    return res;
            }
        }
    } else
    {
//             unsigned long len = length(seq);
//             // write sequence blocks of 70 characters width
//             unsigned long l = 0;
//             while (l < len)
//             {
//                 res = streamPut(stream, infix(seq, l, (l+70 > len) ? len : l+70));
//                 if (res)
//                     return res;
//                 res = streamPut(stream, '\n');
//                 if (res)
//                     return res;
//                 l += 70;
//             }
        unsigned long len = length(seq);
        unsigned long ld = len / 70;
        // write sequence blocks of 70 characters width
        for (unsigned long l = 0; l <ld; ++l)
        {
            res = streamPut(stream, infix(seq, l*70, (l+1)*70));
            if (res)
                return res;
            res = streamPut(stream, '\n');
            if (res)
                return res;
        }
        if (len != ld)
        {
            res = streamPut(stream, infix(seq, ld, len));
            if (res)
                return res;
            res = streamPut(stream, '\n');
            if (res)
                return res;
        }
    }
    return 0;
}

template <typename TStream,
          typename TIdString,
          typename TSeqString>
inline int
writeRecord(TStream & stream,
            TIdString const & meta,
            TSeqString const & seq,
            Fasta const & /*tag*/,
            unsigned const options)
{
    int res = streamPut(stream, '>');
    if (res)
        return res;

    res = streamPut(stream, meta);
    if (res)
        return res;

    res = streamPut(stream, '\n');
    if (res)
        return res;

    if (options & LINEBREAKS)
    {
        // introduce linebreaks
        res = _writeRecordImplLineBreaks(stream, seq, Fasta(), options);
        if (res)
            return res;
    } else
    {
        res = streamPut(stream, seq);
        if (res)
            return res;
    }

    res = streamPut(stream, '\n');
    return res;
}

// FASTA
template <typename TStream,
          typename TIdString,
          typename TSeqString>
inline int
writeRecord(TStream & stream,
            TIdString const & meta,
            TSeqString const & seq,
            Fasta const & /*tag*/)
{
    return writeRecord(stream, meta, seq, Fasta(), LINEBREAKS);
}


/**
.Function.writeRecord
..signature:writeRecord(TStream & stream, TIdString const & meta, TSeqString const & seq, Fastq const &, FastAQOutputOptions const options)
..remarks:Writing to FASTQ without specifying qualities currently will result in all qualities being set to "N"
*/
// FASTQ and we have no qualities
template <typename TStream,
          typename TIdString,
          typename TSeqString>
inline int
writeRecord(TStream & stream,
            TIdString const & meta,
            TSeqString const & seq,
            Fastq const & /*tag*/)
{
    int res = streamPut(stream, '@');
    if (res)
        return res;

    res = streamPut(stream, meta);
    if (res)
        return res;

    res = streamPut(stream, '\n');
    if (res)
        return res;

    res = streamPut(stream, seq);
    if (res)
        return res;

    res = streamPut(stream, "\n+\n");
    if (res)
        return res;

    for (unsigned long i = 0; i < length(seq); ++i)
    {
        //TODO(h4nn3s): is there a better alternative to quality 'N' ?
        res = streamPut(stream, 'N');
        if (res)
            return res;
    }
    //TODO(h4nn3s): are many individual writes wise? One could also construct strings and write those. See comments above.
    res = streamPut(stream, '\n');
    return res;
}

/**
.Function.writeRecord
..signature:writeRecord(TRecordReader & reader, TIdString const & meta, TSeqString const & seq, TQualString const & qual, Fastq const &)
*/
// FASTQ and we have the qualities
template <typename TStream,
          typename TIdString,
          typename TSeqString,
          typename TQualString>
inline int
writeRecord(TStream & stream,
            TIdString const & meta,
            TSeqString const & seq,
            TQualString const & qual,
            Fastq const & /*tag*/)
{
    int res = streamPut(stream, '@');
    if (res)
        return res;

    res = streamPut(stream, meta);
    if (res)
        return res;

    res = streamPut(stream, '\n');
    if (res)
        return res;

    res = streamPut(stream, seq);
    if (res)
        return res;

    res = streamPut(stream, "\n+\n");
    if (res)
        return res;

    res = streamPut(stream, qual);
    if (res)
        return res;

    res = streamPut(stream, '\n');
    return res;
}



// ----------------------------------------------------------------------------
// Function write2()
// ----------------------------------------------------------------------------


/**
.Function.write2
..signature:write2(TStream & stream, StringSet<TIdString, TIdSpec> & sequenceIds, StringSet<TSeqString, TSeqSpec> & sequences, Fasta const &)

*/
// FASTA
template <typename TStream,
          typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec>
int write2(TStream & stream,
         StringSet<TIdString, TIdSpec> & sequenceIds,
         StringSet<TSeqString, TSeqSpec> & sequences,
         Fasta const & /*tag*/)
{
    if (length(sequenceIds) != length(sequences))
        return -1;
    for (unsigned long l = 0; l < length(sequences); ++l)
    {
        int res = writeRecord(stream, sequenceIds[l], sequences[l], Fasta());
        if (res)
            return res;
    }
    return 0;
}

/**
.Function.write2
..signature:write2(TStream & stream, StringSet<TIdString, TIdSpec> & sequenceIds, StringSet<TSeqString, TSeqSpec> & sequences, Fastq const &)

*/
// FASTQ and we have no qualities
template <typename TStream,
          typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec>
int write2(TStream & stream,
         StringSet<TIdString, TIdSpec> & sequenceIds,
         StringSet<TSeqString, TSeqSpec> & sequences,
         Fastq const & /*tag*/)
{
    if (length(sequenceIds) != length(sequences))
        return -1;
    for (unsigned long l = 0; l < length(sequences); ++l)
    {
        int res = writeRecord(stream, sequenceIds[l], sequences[l], Fastq());
        if (res)
            return res;
    }
    return 0;
}

/**
.Function.write2
..signature:write2(TStream & stream, StringSet<TIdString, TIdSpec> & sequenceIds, StringSet<TSeqString, TSeqSpec> & sequences, StringSet<TQualString, TQualSpec> & qualities, Fastq const &)

*/
// FASTQ and we have qualities
template <typename TStream,
          typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec,
          typename TQualString, typename TQualSpec>
int write2(TStream & stream,
         StringSet<TIdString, TIdSpec> & sequenceIds,
         StringSet<TSeqString, TSeqSpec> & sequences,
         StringSet<TQualString, TQualSpec> & qualities,
         Fastq const & /*tag*/)
{
    if (length(sequenceIds) != length(sequences) ||
        length(sequenceIds) != length(qualities))
        return -1;
    for (unsigned long l = 0; l < length(sequences); ++l)
    {
        int res = writeRecord(stream,
                              sequenceIds[l],
                              sequences[l],
                              qualities[l],
                              Fastq());
        if (res)
            return res;
    }
    return 0;
}


}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_WRITE_H_
