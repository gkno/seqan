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
// Tests for seqan/stream/write_fasta_fastq.h
// ==========================================================================


#ifndef TEST_STREAM_TEST_STREAM_WRITE_FASTA_H_
#define TEST_STREAM_TEST_STREAM_WRITE_FASTA_H_

SEQAN_DEFINE_TEST(test_stream_write_record_fasta)
{
    using namespace seqan;

    CharString tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strncpy(filenameBuffer, toCString(tempFilename), 999);

    std::fstream file(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(file.is_open());

    CharString meta1 = "sequence with a nice name !?)=";
    CharString seq1 = "ACGTACGTACGTACGATCACGATACGACTAGACTACGACTAGGGGCTACGACGAGCG\
ACGACGTACGACGACGACGGACGACAGCGATCACGATCTACGAGCTTACTATGGAGCGGGCGATCGAGC";
    CharString meta2 = "           .oO0Oo. .oO0Oo. .oO0Oo. .oO0Oo. .oO0Oo. ";
    CharString seq2 = "GGGGGGGGG";

    int res = writeRecord(file, meta1, seq1, Fasta());
    SEQAN_ASSERT_EQ(res, 0);
    res = writeRecord(file, meta2, seq2, Fasta());
    SEQAN_ASSERT_EQ(res, 0);

    res = writeRecord(file, meta1, seq1, Fasta());
    SEQAN_ASSERT_EQ(res, 0);
    res = writeRecord(file, meta2, seq2, Fasta());
    SEQAN_ASSERT_EQ(res, 0);

    res = writeRecord(file, meta1, seq1, Fasta());
    SEQAN_ASSERT_EQ(res, 0);
    res = writeRecord(file, meta2, seq2, Fasta());
    SEQAN_ASSERT_EQ(res, 0);

    res = writeRecord(file, meta1, seq1, Fasta());
    SEQAN_ASSERT_EQ(res, 0);
    res = writeRecord(file, meta2, seq2, Fasta());
    SEQAN_ASSERT_EQ(res, 0);

    file.seekg(0);
    file.seekp(0);

    typedef RecordReader<std::fstream, SinglePass<void> > TRecordReader;
    TRecordReader reader(file);
    AutoSeqStreamFormat tagSelect;

    SEQAN_ASSERT(checkStreamFormat(reader, tagSelect));
    SEQAN_ASSERT_EQ(tagSelect.tagId, 1);

    for (unsigned int i = 0; i < 8; ++i)
    {
        CharString meta;
        CharString seq;
        res = readRecord(meta, seq, reader, Fasta());
        SEQAN_ASSERT_EQ(res, 0);
        switch(i)
        {
            case 0:
            case 2:
            case 4:
            case 6:
                SEQAN_ASSERT_EQ(meta, meta1);
                SEQAN_ASSERT_EQ(seq, seq1);
                break;
            case 1:
            case 3:
            case 5:
            case 7:
                SEQAN_ASSERT_EQ(meta, meta2);
                SEQAN_ASSERT_EQ(seq, seq2);
                break;
        }
    }

    file.close();
}

#endif // def TEST_STREAM_TEST_STREAM_WRITE_FASTA_H_
