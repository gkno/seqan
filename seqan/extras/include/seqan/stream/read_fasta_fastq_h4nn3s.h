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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Record and Document Reading for FASTA and FASTQ files.
// This implementation aims to be more readable than holtgrew's by using
// tokenizing functions from tokenize.h
// ==========================================================================

//TODO(h4nn3s): double-check if we really want to allow EOF inside meta-line
// and also if a fastq file is legal if a record contains no qualities

#ifndef SEQAN_STREAM_READ_FASTA_FASTQ_H_
#define SEQAN_STREAM_READ_FASTA_FASTQ_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Tag.File Format.tag.Fasta:
    FASTA file format for sequences.
..include:seqan/file.h
*/
struct TagFasta_;
typedef Tag<TagFasta_> const Fasta;

/**
.Tag.File Format.tag.Fastq:
    FASTQ file format for sequences.
..include:seqan/file.h
*/
struct TagFastq_;
typedef Tag<TagFastq_> const Fastq;


// ============================================================================
// Metafunctions
// ============================================================================

template <typename TTag>
struct MetaFirstChar_;

template <>
struct MetaFirstChar_<Tag<TagFasta_> >
{
    static const char VALUE = '>';
};

template <>
struct MetaFirstChar_<Tag<TagFastq_> >
{
    static const char VALUE = '@';
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _clearAndReserveMemory() and helpers
// ----------------------------------------------------------------------------

// if target string is CharString, Sequence is any alphabetical
template <typename TRecordReader>
inline int
_countSequenceFastAQ(unsigned int & count,
                     TRecordReader & reader,
                     char const & /* Alphabet type */)
{
    return _countHelper(count, reader, Alpha_(), Whitespace_(), false);
}

// allow fine-grained alphabet-checking for DnaString etc
template <typename TAlph, typename TRecordReader>
inline int
_countSequenceFastAQ(unsigned int & count,
                     TRecordReader & reader,
                     TAlph const & /* Alphabet type */)
{
    return _countHelper(count, reader, Tag<TAlph>(), Whitespace_(), false);
}

template <typename TSeqAlph,
          typename TFile,
          typename TSpec,
          typename TTag>
inline int
_countMetaAndSequence(unsigned int & metaLength,
                      unsigned int & seqLength,
                      RecordReader<TFile, DoublePass<TSpec> > & reader,
                      TTag const & /*tag*/,
                      TSeqAlph const & /* tag*/)
{
    metaLength=0;
    seqLength=0;

    // COUNT META
    if (atEnd(reader) || value(reader) != MetaFirstChar_<TTag>::VALUE)
        return RecordReader<TFile, DoublePass<TSpec> >::INVALID_FORMAT;
    goNext(reader);
    if (resultCode(reader))
        return resultCode(reader);
    if (atEnd(reader)) // empty ID, no sequence, this is legal TODO?
        return 0;

    int res = countLine(metaLength, reader);
    if (res == EOF_BEFORE_SUCCESS)  // EOF half way in ID is legal
        return 0;
    else if (res)
        return res;

    if (atEnd(reader)) // no sequence
        return 0;

    // COUNT SEQUENCE
    res = _countSequenceFastAQ(seqLength, reader, TSeqAlph());
    if (res == EOF_BEFORE_SUCCESS)  // EOF half way in sequence is legal
        return 0;
    else if (res)
        return res;

    return 0;
}

// SINGLE-Pass
template <typename TIdString,
          typename TSeqString,
          typename TFile,
          typename TSpec,
          typename TTag>
inline int
_clearAndReserveMemory(TIdString & meta, TSeqString & seq,
                       RecordReader<TFile, SinglePass<TSpec> > & /**/,
                       TTag const & /*tag*/)
{
    clear(meta);
    clear(seq);
    return 0;
}

// DOUBLE-Pass
template <typename TIdString,
          typename TSeqString,
          typename TFile,
          typename TSpec,
          typename TTag>
inline int
_clearAndReserveMemory(TIdString & meta, TSeqString & seq,
                       RecordReader<TFile, DoublePass<TSpec> > & reader,
                       TTag const & /*tag*/)
{
    clear(meta);
    clear(seq);
    startFirstPass(reader);

    unsigned int metaLength=0;
    unsigned int seqLength=0;
    // COUNT
    int res = _countMetaAndSequence(metaLength,
                                    seqLength,
                                    reader,
                                    TTag(),
                                    typename Value<TSeqString>::Type());
    if ((res != 0) && (res != EOF_BEFORE_SUCCESS))
        return res;

    // RESERVE FOR META
    reserve(meta, metaLength, Exact());

    // RESERVE FOR SEQUENCE
    reserve(seq, seqLength, Exact());

    startSecondPass(reader);
    return 0;
}

// ----------------------------------------------------------------------------
// Function _readMetaAndSequence() and helpers
// ----------------------------------------------------------------------------

// if target string is CharString, Sequence is any alph()
template <typename TSpec, typename TRecordReader>
inline int
_readSequenceFastAQ(String<char, TSpec> & string,
                    TRecordReader & reader)
{
    return _readHelper(string, reader, Alpha_(), Whitespace_(), false);
}

// allow fine-grained alphabet-checking for DnaString etc
template <typename TAlph, typename TSpec, typename TRecordReader>
inline int
_readSequenceFastAQ(String<TAlph, TSpec> & string,
                    TRecordReader & reader)
{
    return _readHelper(string, reader, Tag<TAlph>(), Whitespace_(), false);
}

// if target string is something else, assume alphabet is any alph()
template <typename TString, typename TRecordReader>
inline int
_readSequenceFastAQ(TString & string,
                    TRecordReader & reader)
{
    return _readHelper(string, reader, Alpha_(), Whitespace_(), false);
}

// This reads Meta and Sequence
template <typename TIdString,
          typename TSeqString,
          typename TFile,
          typename TPass,
          typename TTag>
inline int
_readMetaAndSequence(TIdString & meta, TSeqString & seq,
                     RecordReader<TFile, TPass > & reader,
                     TTag const & /*tag*/)
{
    // READ META
    if (atEnd(reader) || value(reader) != MetaFirstChar_<TTag>::VALUE)
        return RecordReader<TFile, TPass>::INVALID_FORMAT;
    goNext(reader);
    if (resultCode(reader))
        return resultCode(reader);
    if (atEnd(reader)) // empty ID, no sequence, this is legal
        return 0;

    int res = readLine(meta, reader);
    if (res == EOF_BEFORE_SUCCESS)  // EOF half way in ID is legal
        return 0;
    else if (res)
        return res;

    if (atEnd(reader)) // no sequence
        return 0;

    // READ SEQUENCE
    res = _readSequenceFastAQ(seq, reader);
    if (res == EOF_BEFORE_SUCCESS)  // EOF half way in sequence is legal
        return 0;
    else if (res)
        return res;

    return 0;
}

// ----------------------------------------------------------------------------
// Function _skipQualityBlock() and _readQualityBlock()
// ----------------------------------------------------------------------------

template <typename TFile, typename TPass>
inline int
_skipQualityBlock(RecordReader<TFile, TPass > & /**/,
                  unsigned const /**/,
                  Fasta const & /*tag*/)
{
    // NOOP for Fasta
    return 0;
}

template <typename TFile, typename TPass>
inline int
_skipQualityBlock(RecordReader<TFile, TPass > & reader,
                  unsigned const seqLength,
                  Fastq const & /*tag*/)
{
    int res = 0;
    // SKIP QUALITIES' META
    skipLine(reader);
    if (res == EOF_BEFORE_SUCCESS)  // EOF half way in ID is legal
        return 0;
    else if (res)
        return res;

    // SKIP QUALITIES
    res = skipNCharsIgnoringWhitespace(reader, seqLength); // there have to be n qualities
    if (res)
        return RecordReader<TFile, TPass >::INVALID_FORMAT;
    skipLine(reader); // goto to next line if it exists, result is unimportant

    return 0;
}

// This reads Meta and Sequence
template <typename TIdString,
          typename TQualString,
          typename TFile,
          typename TPass>
inline int
_readQualityBlock(TQualString & qual,
                  RecordReader<TFile, TPass > & reader,
                  unsigned const seqLength,
                  TIdString const & meta,
                  Fastq const & /*tag*/)
{
    // READ AND CHECK QUALITIES' META
    if (atEnd(reader) || value(reader) != '+')
        return RecordReader<TFile, TPass >::INVALID_FORMAT;
    goNext(reader);
    if (resultCode(reader))
        return resultCode(reader);
    if (atEnd(reader)) // empty ID, no sequence, this is legal? TODO
        return 0;

    CharString qualmeta_buffer;
    int res = readLine(qualmeta_buffer, reader);
    if (res)
        return RecordReader<TFile, TPass >::INVALID_FORMAT;

    // meta string has to be empty or identical to sequence's meta
    if ((qualmeta_buffer != "") && (qualmeta_buffer != meta))
        return RecordReader<TFile, TPass >::INVALID_FORMAT;

    if (atEnd(reader)) // empty qualities, is this legal?
        return 0;

    // READ QUALITIES
    reserve(qual, seqLength, Exact());
    res = readNCharsIgnoringWhitespace(qual, reader, seqLength);
    // there have to be n qualities
    if (res)
        return RecordReader<TFile, TPass >::INVALID_FORMAT;
    skipLine(reader); // goto to next line if it exists, result is unimportant

    return 0;
}

// ----------------------------------------------------------------------------
// Function readRecord()                               [Single and double pass]
// ----------------------------------------------------------------------------

// FASTA or FASTQ, if we don't want the qualities
template <typename TIdString,
          typename TSeqString,
          typename TFile,
          typename TPass,
          typename TTag>
inline int
readRecord(TIdString & meta,
           TSeqString & seq,
           RecordReader<TFile, TPass > & reader,
           TTag const & /*tag*/)
{
    int res = _clearAndReserveMemory(meta, seq, reader, TTag());
    if (res)
        return res;

    res = _readMetaAndSequence(meta, seq, reader, TTag());
    if (res)
        return res;

    return _skipQualityBlock(reader, length(seq), TTag());
}

// FASTQ and we want the qualities
template <typename TIdString,
          typename TSeqString,
          typename TQualString,
          typename TFile,
          typename TPass>
inline int
readRecord(TIdString & meta,
           TSeqString & seq,
           TQualString & qual,
           RecordReader<TFile, TPass > & reader,
           Fastq const & /*tag*/)
{
    int res = _clearAndReserveMemory(meta, seq, reader, Fastq());
    if (res)
        return res;

    res = _readMetaAndSequence(meta, seq, reader, Fastq());
    if (res)
        return res;

    return _readQualityBlock(qual, reader, length(seq), meta, Fastq());
}

// ----------------------------------------------------------------------------
// Function read();  Double-Pass.
// ----------------------------------------------------------------------------

// generic allocation, e.g. allocator string
template <typename TIdString,
          typename TSeqString,
          typename TQualString>
void _readFastAQAllocate(
                StringSet<TIdString, Owner<ConcatDirect<> > > & sequenceIds,
                StringSet<TSeqString, Owner<ConcatDirect<> > > & sequences,
                StringSet<TQualString, Owner<ConcatDirect<> > > & qualities,
                String<unsigned> const & metaLengths,
                String<unsigned> const & seqLengths,
                bool const withQual)
{
    unsigned int len = length(metaLengths);
    resize(sequenceIds, len + 1, Exact());
    resize(sequences, len + 1, Exact());
    if (withQual)
        resize(qualities.limits, len + 1, Exact());

    unsigned long metaLengthsSum = 0;
    unsigned long seqLengthsSum = 0;

    for (unsigned int i = 0; i < len; ++i)
    {
        metaLengthsSum += metaLengths[i];
        seqLengthsSum += seqLengths[i];
    }
    reserve(sequenceIds.concat, metaLengthsSum + 1, Exact());
    reserve(sequences.concat, seqLengthsSum + 1, Exact());
    if (withQual)
        reserve(qualities.concat, seqLengthsSum + 1, Exact());
}

// generic allocation, e.g. allocator string
template <typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec,
          typename TQualString, typename TQualSpec>
void _readFastAQAllocate(StringSet<TIdString, TIdSpec> & sequenceIds,
                StringSet<TSeqString, TSeqSpec> & sequences,
                StringSet<TQualString, TQualSpec> & qualities,
                String<unsigned> const & metaLengths,
                String<unsigned> const & seqLengths,
                bool const withQual)
{
    unsigned int len = length(metaLengths);
    resize(sequenceIds, len + 1, Exact());
    resize(sequences, len + 1, Exact());
    if (withQual)
        resize(qualities, len + 1, Exact());

    for (unsigned int i = 0; i < len; ++i)
    {
        reserve(sequenceIds[i], metaLengths[i], Exact());
        reserve(sequences[i], seqLengths[i], Exact());
        if (withQual)
            reserve(qualities[i], seqLengths[i], Exact());
    }
}

// Reads a whole FASTA/FASTQ file into string sets, optimizing memory usage.
template <typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec,
          typename TQualString, typename TQualSpec,
          typename TFile,
          typename TSpec,
          typename TTag>
int _readFastAQ(StringSet<TIdString, TIdSpec> & sequenceIds,
                StringSet<TSeqString, TSeqSpec> & sequences,
                StringSet<TQualString, TQualSpec> & qualities,
                RecordReader<TFile, DoublePass<TSpec> > & reader,
                bool const withQual,
                TTag const & /*tag*/)
{
    int res = 0;
    String<unsigned> metaLengths;
    String<unsigned> seqLengths;

    // ------------------------------------------------------------------------
    // First Pass: Compute meta and sequence lengths.
    // ------------------------------------------------------------------------
    startFirstPass(reader);
    size_t sequenceCount = 0;
    while (!atEnd(reader))
    {
        sequenceCount += 1;
        unsigned metaLength = 0;
        unsigned seqLength = 0;
        res = _countMetaAndSequence(metaLength,
                                    seqLength,
                                    reader,
                                    TTag(),
                                    typename Value<TSeqString>::Type());
        if (res)
            return res;

        append(metaLengths, metaLength, Generous());
        append(seqLengths, seqLength, Generous());

        res = _skipQualityBlock(reader, seqLength, TTag());
        if (res)
            return res;
    }

    // ------------------------------------------------------------------------
    // Allocate memory.
    // ------------------------------------------------------------------------
    clear(sequenceIds);
    clear(sequences);
    if (withQual)
        clear(qualities);

    _readFastAQAllocate(sequences, sequenceIds, qualities,
                        metaLengths, seqLengths, withQual);

    // ------------------------------------------------------------------------
    // Second Pass: Actually read data.
    // ------------------------------------------------------------------------
    startSecondPass(reader);
    for (unsigned int i = 0; i < sequenceCount; ++i)
    {
        res = _readMetaAndSequence(sequenceIds[i],
                                   sequences[i],
                                   reader,
                                   TTag());
        switch(res)
        {
            case 0:
                break;
            case EOF_BEFORE_SUCCESS:        // file may end without newline
                if (i >= sequenceCount -1)
                    break;
            default:
                return res;
        }
        if (withQual)
            res = _readQualityBlock(qualities[i],
                                    reader,
                                    length(sequences[i]),
                                    sequenceIds[i],
                                    Fastq());
        else
            res = _skipQualityBlock(reader, length(sequences[i]), TTag());
        if (res)
            return res;
    }
    return 0;
}

// FASTA
template <typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec,
          typename TFile,
          typename TSpec>
int read2(StringSet<TIdString, TIdSpec> & sequenceIds,
         StringSet<TSeqString, TSeqSpec> & sequences,
         RecordReader<TFile, DoublePass<TSpec> > & reader,
         Fasta const & /*tag*/)
{
    StringSet<CharString, TSeqSpec> qualities;
    return _readFastAQ(sequenceIds, sequences, qualities, reader, false, Fasta());
}

// FASTQ, if we don't want the qualities
template <typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec,
          typename TFile,
          typename TSpec>
int read2(StringSet<TIdString, TIdSpec> & sequenceIds,
         StringSet<TSeqString, TSeqSpec> & sequences,
         RecordReader<TFile, DoublePass<TSpec> > & reader,
         Fastq const & /*tag*/)
{
    StringSet<CharString, TSeqSpec> qualities;
    return _readFastAQ(sequenceIds, sequences, qualities, reader, false, Fastq());
}

// FASTQ and we want Qualities
template <typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec,
          typename TQualString, typename TQualSpec,
          typename TFile,
          typename TSpec>
int read2(StringSet<TIdString, TIdSpec> & sequenceIds,
         StringSet<TSeqString, TSeqSpec> & sequences,
         StringSet<TQualString, TQualSpec> & qualities,
         RecordReader<TFile, DoublePass<TSpec> > & reader,
         Fastq const & /*tag*/)
{
    return _readFastAQ(sequenceIds, sequences, qualities, reader, true, Fastq());
}



// We can give an especially efficient implementation for ConcatDirect string
// sets.
//TODO

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_READ_FASTA_FASTQ_H_
