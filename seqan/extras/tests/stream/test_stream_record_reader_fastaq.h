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
// Tests for seqan/stream/read_fasta_fastq.h
// ==========================================================================


#ifndef TEST_STREAM_TEST_STREAM_READ_FASTAQ_H_
#define TEST_STREAM_TEST_STREAM_READ_FASTAQ_H_

// ------------ FASTA -------------

std::fstream* createFastAFile(seqan::CharString &tempFilename)
{
    using namespace seqan;

    tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strncpy(filenameBuffer, toCString(tempFilename), 999);

    std::fstream *file = new std::fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(file->is_open());

    char const * STR =
"> sequenceID_with special chars an irregular linebreaks\n\
AAAACGTGCGGTTGGGCAAAAAACTTTCTTATATTCTATCTATCTTGTAGCTAGCTGTAGCTAGCTAGCATCGTAGC\n\
CCCAGAGTGTCATGCATGTCGA\n\
ACGTGTTTTTGGGGCGGTTATATATATATATATT\n\
\n\
>sequence2... with no linebreaks and no newline at end\n\
ACGTNNNNNNNCGTACTTGCTAGCTAGCTAGCTAGCTAGCATCGTACGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTACTATCATCTACTATCTACTATCATCTACTATCATTCATCGATCGATCGATCGTACGTACGATCGATCGATCGATCGTACGATCGATGCTACGTACGTACG";


    file->write(STR, strlen(STR));
    file->seekg(0);
    file->seekp(0);
    return file;
}

template <typename TRecordReader>
void FASTA_TEST(TRecordReader & reader)
{
    using namespace seqan;

    SEQAN_ASSERT(checkStreamFormat(reader, Fasta()));

    CharString meta;
    Dna5String seq;

    int res = readRecord(meta, seq, reader, Fasta());
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(meta, " sequenceID_with special chars an irregular linebreaks");
    SEQAN_ASSERT_EQ(seq, "AAAACGTGCGGTTGGGCAAAAAACTTTCTTATATTCTATCTATCTTGTAGCTAGCTGTAGCTAGCTAGCATCGTAGCCCCAGAGTGTCATGCATGTCGAACGTGTTTTTGGGGCGGTTATATATATATATATT");
    SEQAN_ASSERT_EQ(value(reader), '>');

    res = readRecord(meta, seq, reader, Fasta());
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(meta, "sequence2... with no linebreaks and no newline at end");
    SEQAN_ASSERT_EQ(seq, "ACGTNNNNNNNCGTACTTGCTAGCTAGCTAGCTAGCTAGCATCGTACGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTACTATCATCTACTATCTACTATCATCTACTATCATTCATCGATCGATCGATCGTACGTACGATCGATCGATCGATCGTACGATCGATGCTACGTACGTACG");
    SEQAN_ASSERT(atEnd(reader));
}

template <typename TRecordReader>
void FASTA_TEST_BATCH(TRecordReader & reader)
{
    using namespace seqan;

    SEQAN_ASSERT(checkStreamFormat(reader, Fasta()));

    StringSet<CharString> metas;
    StringSet<Dna5String> seqs;

    int res = read2(metas, seqs, reader, Fasta());
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(metas[0], " sequenceID_with special chars an irregular linebreaks");
    SEQAN_ASSERT_EQ(seqs[0], "AAAACGTGCGGTTGGGCAAAAAACTTTCTTATATTCTATCTATCTTGTAGCTAGCTGTAGCTAGCTAGCATCGTAGCCCCAGAGTGTCATGCATGTCGAACGTGTTTTTGGGGCGGTTATATATATATATATT");

    SEQAN_ASSERT_EQ(metas[1], "sequence2... with no linebreaks and no newline at end");
    SEQAN_ASSERT_EQ(seqs[1], "ACGTNNNNNNNCGTACTTGCTAGCTAGCTAGCTAGCTAGCATCGTACGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTACTATCATCTACTATCTACTATCATCTACTATCATTCATCGATCGATCGATCGTACGTACGATCGATCGATCGATCGTACGATCGATGCTACGTACGTACG");
    SEQAN_ASSERT(atEnd(reader));
}

template <typename TRecordReader>
void FASTA_TEST_BATCH_CONCAT(TRecordReader & reader)
{
    using namespace seqan;

    SEQAN_ASSERT(checkStreamFormat(reader, Fasta()));

    StringSet<CharString, Owner<ConcatDirect<> > > metas;
    StringSet<Dna5String, Owner<ConcatDirect<> > > seqs;

    int res = read2(metas, seqs, reader, Fasta());

    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(metas[0], " sequenceID_with special chars an irregular linebreaks");
    SEQAN_ASSERT_EQ(seqs[0], "AAAACGTGCGGTTGGGCAAAAAACTTTCTTATATTCTATCTATCTTGTAGCTAGCTGTAGCTAGCTAGCATCGTAGCCCCAGAGTGTCATGCATGTCGAACGTGTTTTTGGGGCGGTTATATATATATATATT");

    SEQAN_ASSERT_EQ(metas[1], "sequence2... with no linebreaks and no newline at end");
    SEQAN_ASSERT_EQ(seqs[1], "ACGTNNNNNNNCGTACTTGCTAGCTAGCTAGCTAGCTAGCATCGTACGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTACTATCATCTACTATCTACTATCATCTACTATCATTCATCGATCGATCGATCGTACGTACGATCGATCGATCGATCGTACGATCGATGCTACGTACGTACG");
    SEQAN_ASSERT(atEnd(reader));
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fasta_single_fstream)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastAFile(filename);

    typedef RecordReader<std::fstream, SinglePass<void> > TRecordReader;
    TRecordReader reader(*file);

    FASTA_TEST(reader);

    file->close();
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fasta_double_fstream)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastAFile(filename);

    typedef RecordReader<std::fstream, DoublePass<void> > TRecordReader;
    TRecordReader reader(*file);

    FASTA_TEST(reader);

    file->close();
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fasta_batch_fstream)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastAFile(filename);

    typedef RecordReader<std::fstream, DoublePass<void> > TRecordReader;
    TRecordReader reader(*file);

    FASTA_TEST_BATCH(reader);

    file->close();
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fasta_single_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastAFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, SinglePass<Mapped> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTA_TEST(reader);

    close(mmapString);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fasta_double_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastAFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, DoublePass<Mapped> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTA_TEST(reader);

    close(mmapString);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fasta_batch_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastAFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, DoublePass<Mapped> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTA_TEST_BATCH(reader);

    close(mmapString);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fasta_batch_concat_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastAFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, DoublePass<Mapped> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTA_TEST_BATCH_CONCAT(reader);

    close(mmapString);
}

// -------------- FASTQ

std::fstream* createFastQFile(seqan::CharString &tempFilename)
{
    using namespace seqan;

    tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strncpy(filenameBuffer, toCString(tempFilename), 999);

    std::fstream *file = new std::fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(file->is_open());

    char const * STR =
"@SEQ_ID\n\
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT\n\
+\n\
!''*((((***+))%%%++)(%%%%).@***-+*''))**55CCF>>>>>>CCCCCCC65\n\
@ 2ndSequence with formatting obscurities\n\
GATTTGGGGTTCAAAGC\n\
AGTATCGATCAAATAGTAAATCCATTT\n\
GTTCAACTCACAGTTT\n\
+ 2ndSequence with formatting obscurities\n\
!''*((((***+))%%%\n\
++)(%%%%).\n\
@***-+*''))**55CCF>>>>>>CCCCCCC65";
    // the second quality-string is line-broken at places to produce lines
    // beginning with '+' and '@'. The new code handles this absolutely fine!

    file->write(STR, strlen(STR));
    file->seekg(0);
    file->seekp(0);
    return file;
}

template <typename TRecordReader>
void FASTQ_TEST(TRecordReader & reader)
{
    using namespace seqan;

    SEQAN_ASSERT(checkStreamFormat(reader, Fastq()));

    CharString meta;
    DnaString seq;
    CharString qual;

    // qualities are not asked for, but should be skipped internally
    int res = readRecord(meta, seq, reader, Fastq());
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(meta, "SEQ_ID");
    SEQAN_ASSERT_EQ(seq,
                "GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT");
    SEQAN_ASSERT_EQ(value(reader), '@');

    // now qualities are also requested
    res = readRecord(meta, seq, qual, reader, Fastq());
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(meta, " 2ndSequence with formatting obscurities");
    SEQAN_ASSERT_EQ(seq,
                "GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT");
    SEQAN_ASSERT_EQ(qual, "!''*((((***+))%%%++)(%%%%).@***-+*''))**55CCF>>>>>>CCCCCCC65");
    SEQAN_ASSERT(atEnd(reader));
}

template <typename TRecordReader>
void FASTQ_TEST_BATCH(TRecordReader & reader)
{
    using namespace seqan;

    SEQAN_ASSERT(checkStreamFormat(reader, Fastq()));

    StringSet<CharString> metas;
    StringSet<DnaString> seqs;
    StringSet<CharString> quals;

    int res = read2(metas, seqs, quals, reader, Fastq());
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(metas[0], "SEQ_ID");
    SEQAN_ASSERT_EQ(seqs[0],
                "GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT");
    SEQAN_ASSERT_EQ(quals[0], "!''*((((***+))%%%++)(%%%%).@***-+*''))**55CCF>>>>>>CCCCCCC65");
    SEQAN_ASSERT_EQ(metas[1], " 2ndSequence with formatting obscurities");
    SEQAN_ASSERT_EQ(seqs[1],
                "GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT");
    SEQAN_ASSERT_EQ(quals[1], "!''*((((***+))%%%++)(%%%%).@***-+*''))**55CCF>>>>>>CCCCCCC65");
    SEQAN_ASSERT(atEnd(reader));
}

template <typename TRecordReader>
void FASTQ_TEST_BATCH_CONCAT(TRecordReader & reader)
{
    using namespace seqan;

    SEQAN_ASSERT(checkStreamFormat(reader, Fastq()));

    StringSet<CharString, Owner<ConcatDirect<> > > metas;
    StringSet<DnaString, Owner<ConcatDirect<> > > seqs;
    StringSet<CharString, Owner<ConcatDirect<> > > quals;

    int res = read2(metas, seqs, quals, reader, Fastq());
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(metas[0], "SEQ_ID");
    SEQAN_ASSERT_EQ(seqs[0],
                "GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT");
    SEQAN_ASSERT_EQ(quals[0], "!''*((((***+))%%%++)(%%%%).@***-+*''))**55CCF>>>>>>CCCCCCC65");
    SEQAN_ASSERT_EQ(metas[1], " 2ndSequence with formatting obscurities");
    SEQAN_ASSERT_EQ(seqs[1],
                "GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT");
    SEQAN_ASSERT_EQ(quals[1], "!''*((((***+))%%%++)(%%%%).@***-+*''))**55CCF>>>>>>CCCCCCC65");
    SEQAN_ASSERT(atEnd(reader));
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fastq_single_fstream)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastQFile(filename);

    typedef RecordReader<std::fstream, SinglePass<void> > TRecordReader;
    TRecordReader reader(*file);

    FASTQ_TEST(reader);

    file->close();
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fastq_double_fstream)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastQFile(filename);

    typedef RecordReader<std::fstream, DoublePass<void> > TRecordReader;
    TRecordReader reader(*file);

    FASTQ_TEST(reader);

    file->close();
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fastq_batch_fstream)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastQFile(filename);

    typedef RecordReader<std::fstream, DoublePass<void> > TRecordReader;
    TRecordReader reader(*file);

    FASTQ_TEST_BATCH(reader);

    file->close();
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fastq_single_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastQFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, SinglePass<Mapped> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTQ_TEST(reader);

    close(mmapString);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fastq_double_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastQFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, DoublePass<Mapped> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTQ_TEST(reader);

    close(mmapString);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fastq_batch_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastQFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, DoublePass<Mapped> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTQ_TEST_BATCH(reader);

    close(mmapString);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fastq_batch_concat_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastQFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, DoublePass<Mapped> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTQ_TEST_BATCH_CONCAT(reader);

    close(mmapString);
}


#endif // def TEST_STREAM_TEST_STREAM_READ_FASTAQ_H_
