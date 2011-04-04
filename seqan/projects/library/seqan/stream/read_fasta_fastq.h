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
// Record and Document Reading for FASTA and FASTQ files.
// ==========================================================================

// TODO(holtgrew): Fasta specializations Fasta<TSpec = WholeMeta>, Fasta<IdOnly>
// TODO(holtgrew): FASTQ reading.

#ifndef SEQAN_STREAM_READ_FASTA_FASTQ_H_
#define SEQAN_STREAM_READ_FASTA_FASTQ_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TIdString, typename TSeqString>
class FastaReaderLambdaContextSinglePass_;

// Used for single-pass and second pass in double-pass.
template <typename TIdString, typename TSeqString>
class FastaReaderLambdaContextSinglePass_
{
public:
    Pair<TIdString, TSeqString> & record;
    FastaReaderLambdaContextSinglePass_(Pair<TIdString, TSeqString> & _record) : record(_record) {}
};

// Used by first pass of double-pass only.
template <typename TIdString, typename TSeqString>
class FastaReaderLambdaContextFirstPass_
{
public:
    unsigned metaLength;
    unsigned sequenceLength;
    FastaReaderLambdaContextFirstPass_() : metaLength(0), sequenceLength(0)
    {}
};

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
// Helper Function _readRecordFastAQMeta
// ----------------------------------------------------------------------------

// Forwards for for _readRecordFastAQMeta and _readRecordFastAQSequence.
template <typename TLambdaContext, typename TRecordReader>
inline void _lambdaFastAQMetaChar(TLambdaContext & lambdaContext, TRecordReader & reader);
template <typename TLambdaContext, typename TRecordReader>
inline void _lambdaFastAQSequenceChar(TLambdaContext & lambdaContext, TRecordReader & reader);

// Read the meta field from a FASTA or FASTQ file.
template <typename TLambdaContext, typename TFile, typename TPass, typename TTag>
int _readRecordFastAQMeta(TLambdaContext & lambdaContext,
                          RecordReader<TFile, TPass> & reader,
                          TTag const & /*tag*/)
{
    typedef RecordReader<TFile, TPass> TRecordReader;

    // Read meta field.
    // std::cerr << __LINE__ << " value(reader) == " << value(reader) << std::endl;
    if (value(reader) != MetaFirstChar_<TTag>::VALUE)
        return TRecordReader::INVALID_FORMAT;  // Did not start with '>'.
    goNext(reader);
    if (atEnd(reader))
        return TRecordReader::INVALID_FORMAT;  // Input stopped after '>'.

    // Read meta field, goes up to the first line break.
    // std::cerr << __LINE__ << " value(reader) == " << value(reader) << std::endl;
    while (value(reader) != '\n' && value(reader) != '\r') {
        _lambdaFastAQMetaChar(lambdaContext, reader);
        goNext(reader);
        if (atEnd(reader))
            return TRecordReader::INVALID_FORMAT;  // Input stopped in meta.
    }
    // Go beyond line break.
    // std::cerr << __LINE__ << " value(reader) == " << value(reader) << std::endl;
    if (value(reader) == '\r') {
        // Skip '\n' in DOS line endings.
        goNext(reader);
        if (atEnd(reader))
            return TRecordReader::INVALID_FORMAT;  // Stopped in meta break.
        // std::cerr << __LINE__ << " value(reader) == " << value(reader) << std::endl;
        if (value(reader) == '\n')
            goNext(reader);
    } else {
        // Skip '\n'.
        goNext(reader);
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Helper Function _readRecordFastAQSequence
// ----------------------------------------------------------------------------

// Read the sequence field from a FASTA or FASTQ file.
template <typename TLambdaContext, typename TFile, typename TPass, typename TTag>
int _readRecordFastAQSequence(TLambdaContext & lambdaContext,
                              RecordReader<TFile, TPass> & reader,
                              TTag const & /*tag*/)
{
    typedef RecordReader<TFile, TPass> TRecordReader;

    // The sequence can be empty.
    while (true) {
        // The file is allowed to end within the sequence.
        if (atEnd(reader))
            return resultCode(reader);
        // Skip line break ('\r', '\n', "\r\n" but also "\n\n" since
        // the code gets simpler this way). The file is allowed to end
        // after the line break.
        if (value(reader) == '\r' || value(reader) == '\n') {
            goNext(reader);
            if (atEnd(reader))
                return resultCode(reader);
            // Skip '\n' after '\r' or first '\n'.
            if (value(reader) == '\n') {
                goNext(reader);
                if (atEnd(reader))
                    return resultCode(reader);
            }
            // '>' after '\r', '\n', '\n\n', or '\r\n'.
            if (value(reader) == MetaFirstChar_<TTag>::VALUE)
                return resultCode(reader);
        }
        _lambdaFastAQSequenceChar(lambdaContext, reader);
        // std::cerr << "i == " << (i++) << " read >>" << value(reader) << "<<" << std::endl;
        goNext(reader);
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function readRecord();  Single-Pass.
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TRecordReader>
inline void
_lambdaFastAQMetaChar(FastaReaderLambdaContextSinglePass_<TIdString, TSeqString> & lambdaContext,
                      TRecordReader & reader)
{
    appendValue(lambdaContext.record.i1, value(reader));
}

template <typename TIdString, typename TSeqString, typename TRecordReader>
inline void
_lambdaFastAQSequenceChar(FastaReaderLambdaContextSinglePass_<TIdString, TSeqString> & lambdaContext,
                          TRecordReader & reader)
{
    appendValue(lambdaContext.record.i2, value(reader));
}

template <typename TIdString, typename TSeqString, typename TFile, typename TSpec>
int readRecord(Pair<TIdString, TSeqString> & record,
               RecordReader<TFile, SinglePass<TSpec> > & reader,
               Fasta const & /*tag*/)
{
    // Nomenclauture: A FASTA file consists of records.  Each record
    // consists of a meta and a sequence field.

    clear(record.i1);
    clear(record.i2);
    FastaReaderLambdaContextSinglePass_<TIdString, TSeqString> lambdaContext(record);

    typedef RecordReader<TFile, SinglePass<TSpec> > TRecordReader;

    // Read meta field.
    int res = _readRecordFastAQMeta(lambdaContext, reader, Fasta());
    if (res)
        return res;
    // std::cerr << "READ " << record.i1 << std::endl;
    // Read sequence field.
    res = _readRecordFastAQSequence(lambdaContext, reader, Fasta());
    if (res)
        return res;

    return 0;
}

// ----------------------------------------------------------------------------
// Function readRecord();  Double-Pass.
// ----------------------------------------------------------------------------
// TODO(holtgrew): This code has to be adjusted to also work for Fastq.

template <typename TIdString, typename TSeqString, typename TFile, typename TSpec>
inline void
_lambdaFastAQMetaChar(FastaReaderLambdaContextFirstPass_<TIdString, TSeqString> & lambdaContext,
                      RecordReader<TFile, DoublePass<TSpec> > & /*reader*/)
{
    // std::cerr << "lambdaContext.metaLength == " << lambdaContext.metaLength << std::endl;
    // std::cerr << "lambdaContext.sequenceLength == " << lambdaContext.sequenceLength << std::endl;
    // std::cerr << "lambdaContext.metaLength += 1;" << std::endl;
    lambdaContext.metaLength += 1;
}

template <typename TIdString, typename TSeqString, typename TFile, typename TSpec>
inline void
_lambdaFastAQSequenceChar(FastaReaderLambdaContextFirstPass_<TIdString, TSeqString> & lambdaContext,
                          RecordReader<TFile, DoublePass<TSpec> > & /*reader*/)
{
    // std::cerr << "lambdaContext.metaLength == " << lambdaContext.metaLength << std::endl;
    // std::cerr << "lambdaContext.sequenceLength == " << lambdaContext.sequenceLength << std::endl;
    // std::cerr << "lambdaContext.sequenceLength += 1;" << std::endl;
    lambdaContext.sequenceLength += 1;
}

template <typename TIdString, typename TSeqString, typename TFile, typename TSpec>
int readRecord(Pair<TIdString, TSeqString> & record,
               RecordReader<TFile, DoublePass<TSpec> > & reader,
               Fasta const & /*tag*/)
{
    // Nomenclauture: A FASTA file consists of records.  Each record
    // consists of a meta and a sequence field.

    // ------------------------------------------------------------------------
    // First Pass: Compute meta and sequence lengths.
    // ------------------------------------------------------------------------
    startFirstPass(reader);
    FastaReaderLambdaContextFirstPass_<TIdString, TSeqString> lambdaContextFirst;
    // Compute meta field length.
    int res = _readRecordFastAQMeta(lambdaContextFirst, reader, Fasta());
    if (res)
        return res;
    // std::cerr << "lambdaContextFirst.metaLength == " << lambdaContextFirst.metaLength << std::endl;
    // std::cerr << "lambdaContextFirst.sequenceLength == " << lambdaContextFirst.sequenceLength << std::endl;
    // std::cerr << "FIRST PASS DONE 1/2" << std::endl;
    // Compute sequence field length.
    res = _readRecordFastAQSequence(lambdaContextFirst, reader, Fasta());
    if (res)
        return res;

    // std::cerr << "lambdaContextFirst.metaLength == " << lambdaContextFirst.metaLength << std::endl;
    // std::cerr << "lambdaContextFirst.sequenceLength == " << lambdaContextFirst.sequenceLength << std::endl;
    // std::cerr << "FIRST PASS DONE 2/2" << std::endl;
    // ------------------------------------------------------------------------
    // Allocate memory.
    // ------------------------------------------------------------------------
    clear(record.i1);
    reserve(record.i1, lambdaContextFirst.metaLength);
    clear(record.i2);
    reserve(record.i2, lambdaContextFirst.sequenceLength);
    // std::cerr << "lambdaContextFirst.metaLength == " << lambdaContextFirst.metaLength << std::endl;
    // std::cerr << "lambdaContextFirst.sequenceLength == " << lambdaContextFirst.sequenceLength << std::endl;

    // ------------------------------------------------------------------------
    // Second Pass: Actually read data.
    // ------------------------------------------------------------------------
    startSecondPass(reader);
    FastaReaderLambdaContextSinglePass_<TIdString, TSeqString> lambdaContextSecond(record);
    res = _readRecordFastAQMeta(lambdaContextSecond, reader, Fasta());
    // std::cerr << "SECOND PASS DONE 1/2 res == " << res << std::endl;
    if (res)
        return res;
    // Compute sequence field length.
    res = _readRecordFastAQSequence(lambdaContextSecond, reader, Fasta());
    // std::cerr << "SECOND PASS DONE 2/2" << std::endl;
    return res;
}

// ----------------------------------------------------------------------------
// Function read();  Double-Pass.
// ----------------------------------------------------------------------------

// Reads a whole FASTA file into string sets, optimizing memory usage.
template <typename TIdString, typename TSeqString, typename TFile, typename TSpec>
int read(StringSet<TIdString> & sequenceIds,
         StringSet<TSeqString> & sequences,
         RecordReader<TFile, DoublePass<TSpec> > & reader,
         Fasta const & /*tag*/)
{
    int res = 0;

    // ------------------------------------------------------------------------
    // First Pass: Compute meta and sequence lengths.
    // ------------------------------------------------------------------------
    startFirstPass(reader);
    FastaReaderLambdaContextFirstPass_<TIdString, TSeqString> lambdaContextFirst;
    size_t sequenceCount = 0;
    while (!atEnd(reader)) {
        sequenceCount += 1;
        // Compute meta field length.
        res = _readRecordFastAQMeta(lambdaContextFirst, reader, Fasta());
        if (res)
            return res;
        // std::cerr << "lambdaContextFirst.metaLength == " << lambdaContextFirst.metaLength << std::endl;
        // std::cerr << "lambdaContextFirst.sequenceLength == " << lambdaContextFirst.sequenceLength << std::endl;
        // std::cerr << "FIRST PASS DONE 1/2" << std::endl;
        // Compute sequence field length.
        res = _readRecordFastAQSequence(lambdaContextFirst, reader, Fasta());
        if (res)
            return res;
    }

    // ------------------------------------------------------------------------
    // Allocate memory.
    // ------------------------------------------------------------------------
    clear(sequenceIds);
    clear(sequences);
    reserve(sequenceIds, sequenceCount + 1, Exact());
    reserve(sequences, sequenceCount + 1, Exact());

    // ------------------------------------------------------------------------
    // Second Pass: Actually read data.
    // ------------------------------------------------------------------------
    startSecondPass(reader);
    Pair<TIdString, TSeqString> record;
    while (!atEnd(reader)) {
        FastaReaderLambdaContextSinglePass_<TIdString, TSeqString> lambdaContextSecond(record);
        res = _readRecordFastAQMeta(lambdaContextSecond, reader, Fasta());
        // std::cerr << "SECOND PASS DONE 1/2 res == " << res << std::endl;
        if (res)
            return res;
        // Compute sequence field length.
        res = _readRecordFastAQSequence(lambdaContextSecond, reader, Fasta());
        // std::cerr << "SECOND PASS DONE 2/2" << std::endl;
        appendValue(sequenceIds, record.i1); // TODO(holtgrew): Move-construction/appending?
        appendValue(sequences, record.i2);
    }
    return res;
}

// We can give an especially efficient implementation for ConcatDirect string
// sets.
template <typename TIdString, typename TSeqString, typename TFile, typename TSpec>
int read(StringSet<TIdString, Owner<ConcatDirect<> > > & sequenceIds,
         StringSet<TSeqString, Owner<ConcatDirect<> > > & sequences,
         RecordReader<TFile, DoublePass<TSpec> > & reader,
         Fasta const & /*tag*/)
{
    typedef typename Concatenator<StringSet<TIdString, Owner<ConcatDirect<> > > >::Type TIdStringConcat;
    typedef typename Concatenator<StringSet<TSeqString, Owner<ConcatDirect<> > > >::Type TSeqStringConcat;
    // std::cerr << reader._string << std::endl;
    int res = 0;

    // ------------------------------------------------------------------------
    // First Pass: Compute meta and sequence lengths.
    // ------------------------------------------------------------------------
    startFirstPass(reader);
    FastaReaderLambdaContextFirstPass_<TIdString, TSeqString> lambdaContextFirst;
    size_t sequenceCount = 0;
    while (!atEnd(reader)) {
        sequenceCount += 1;
        // Compute meta field length.
        res = _readRecordFastAQMeta(lambdaContextFirst, reader, Fasta());
        if (res)
            return res;
        // std::cerr << "lambdaContextFirst.metaLength == " << lambdaContextFirst.metaLength << std::endl;
        // std::cerr << "lambdaContextFirst.sequenceLength == " << lambdaContextFirst.sequenceLength << std::endl;
        // std::cerr << "FIRST PASS DONE 1/2" << std::endl;
        // Compute sequence field length.
        res = _readRecordFastAQSequence(lambdaContextFirst, reader, Fasta());
        if (res)
            return res;
    }

    // ------------------------------------------------------------------------
    // Allocate memory.
    // ------------------------------------------------------------------------
    clear(sequenceIds);
    clear(sequences);
    reserve(sequenceIds.limits, sequenceCount + 1, Exact());
    reserve(sequenceIds.concat, lambdaContextFirst.metaLength, Exact());
    reserve(sequences.limits, sequenceCount + 1, Exact());
    reserve(sequences.concat, lambdaContextFirst.sequenceLength, Exact());

    // ------------------------------------------------------------------------
    // Second Pass: Actually read data.
    // ------------------------------------------------------------------------
    startSecondPass(reader);
    Pair<TIdStringConcat, TSeqStringConcat> record;
    while (!atEnd(reader)) {
        FastaReaderLambdaContextSinglePass_<TIdString, TSeqString> lambdaContextSecond(record);
        res = _readRecordFastAQMeta(lambdaContextSecond, reader, Fasta());
        if (res)
            return res;
        res = _readRecordFastAQSequence(lambdaContextSecond, reader, Fasta());
        appendValue(sequenceIds.limits, length(sequenceIds.concat));
        appendValue(sequences.limits, length(sequences.concat));
    }
    return 0;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_READ_FASTA_FASTQ_H_
