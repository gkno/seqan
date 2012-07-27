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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef EXTRAS_TESTS_BAM_IO_TEST_EASY_BAM_IO_H_

#include <sstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/bam_io.h>

// ---------------------------------------------------------------------------
// Read Header
// ---------------------------------------------------------------------------

void testBamIOBamStreamReadHeader(char const * pathFragment)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, pathFragment);

    seqan::BamStream bamIO(toCString(filePath), seqan::BamStream::READ);
    SEQAN_ASSERT(isGood(bamIO));

    SEQAN_ASSERT_EQ(length(bamIO.header.records), 2u);
    SEQAN_ASSERT_EQ(bamIO.header.records[0].type, seqan::BAM_HEADER_FIRST);
    SEQAN_ASSERT_EQ(bamIO.header.records[0].tags[0].i1, "VN");
    SEQAN_ASSERT_EQ(bamIO.header.records[0].tags[0].i2, "1.3");
    SEQAN_ASSERT_EQ(bamIO.header.records[0].tags[1].i1, "SO");
    SEQAN_ASSERT_EQ(bamIO.header.records[0].tags[1].i2, "coordinate");
    SEQAN_ASSERT_EQ(bamIO.header.records[1].type, seqan::BAM_HEADER_REFERENCE);
    SEQAN_ASSERT_EQ(bamIO.header.records[1].tags[0].i1, "SN");
    SEQAN_ASSERT_EQ(bamIO.header.records[1].tags[0].i2, "REFERENCE");
    SEQAN_ASSERT_EQ(bamIO.header.records[1].tags[1].i1, "LN");
    SEQAN_ASSERT_EQ(bamIO.header.records[1].tags[1].i2, "10000");
    SEQAN_ASSERT_EQ(length(bamIO.header.sequenceInfos), 1u);
    SEQAN_ASSERT_EQ(bamIO.header.sequenceInfos[0].i1, "REFERENCE");
    SEQAN_ASSERT_EQ(bamIO.header.sequenceInfos[0].i2, 10000);
}

SEQAN_DEFINE_TEST(test_bam_io_bam_stream_sam_read_header)
{
    testBamIOBamStreamReadHeader("/extras/tests/bam_io/small.sam");
}

SEQAN_DEFINE_TEST(test_bam_io_bam_stream_bam_read_header)
{
    testBamIOBamStreamReadHeader("/extras/tests/bam_io/small.bam");
}

// ---------------------------------------------------------------------------
// Read Records
// ---------------------------------------------------------------------------

void testBamIOBamStreamReadRecords(char const * pathFragment)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, pathFragment);

    seqan::BamStream bamIO(toCString(filePath), seqan::BamStream::READ);
    SEQAN_ASSERT(isGood(bamIO));

    seqan::BamAlignmentRecord record;

    seqan::String<seqan::BamAlignmentRecord> alignments;
    while (!atEnd(bamIO))
    {
        resize(alignments, length(alignments) + 1);
        SEQAN_ASSERT_EQ(readRecord(back(alignments), bamIO), 0);
        SEQAN_ASSERT(isGood(bamIO));
    }
    SEQAN_ASSERT_EQ(length(alignments), 3u);

    SEQAN_ASSERT_EQ(alignments[0].qName, "READ0");
    SEQAN_ASSERT_EQ(alignments[0].flag, 2);
    SEQAN_ASSERT_EQ(alignments[0].rId, 0);
    SEQAN_ASSERT_EQ(alignments[0].pos, 0);
    SEQAN_ASSERT_EQ(alignments[0].mapQ, 8);
    SEQAN_ASSERT_EQ(length(alignments[0].cigar), 3u);
    SEQAN_ASSERT_EQ(alignments[0].cigar[0].count, 5u);
    SEQAN_ASSERT_EQ(alignments[0].cigar[0].operation, 'M');
    SEQAN_ASSERT_EQ(alignments[0].cigar[1].count, 1u);
    SEQAN_ASSERT_EQ(alignments[0].cigar[1].operation, 'I');
    SEQAN_ASSERT_EQ(alignments[0].cigar[2].count, 4u);
    SEQAN_ASSERT_EQ(alignments[0].cigar[2].operation, 'M');
    SEQAN_ASSERT_EQ(alignments[0].rNextId, 0);
    SEQAN_ASSERT_EQ(alignments[0].pNext, 30);
    SEQAN_ASSERT_EQ(alignments[0].tLen, 40);
    SEQAN_ASSERT_EQ(alignments[0].seq, "AAAAAAAAAA");
    SEQAN_ASSERT_EQ(alignments[0].qual, "!!!!!!!!!!");
    SEQAN_ASSERT_EQ(length(alignments[0].tags), 0u);

    // TODO(holtgrew): Check more alignments?
}

SEQAN_DEFINE_TEST(test_bam_io_bam_stream_sam_read_records)
{
    testBamIOBamStreamReadRecords("/extras/tests/bam_io/small.sam");
}

SEQAN_DEFINE_TEST(test_bam_io_bam_stream_bam_read_records)
{
    testBamIOBamStreamReadRecords("/extras/tests/bam_io/small.bam");
}

// ---------------------------------------------------------------------------
// Write Header
// ---------------------------------------------------------------------------

void testBamIOBamStreamWriteHeader(char const * pathFragmentExpected)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, pathFragmentExpected);

    seqan::CharString tmpPath = SEQAN_TEMP_FILENAME();
    if (seqan::endsWith(pathFragmentExpected, ".bam"))
        append(tmpPath, ".bam");
    else
        append(tmpPath, ".sam");

    // Initialize BamStream, build header.
    seqan::BamStream bamIO(toCString(tmpPath), seqan::BamStream::WRITE);
    resize(bamIO.header.sequenceInfos, 1);
    bamIO.header.sequenceInfos[0].i1 = "REFERENCE";
    bamIO.header.sequenceInfos[0].i2 = 10000;
    resize(bamIO.header.records, 2);
    resize(bamIO.header.records[0].tags, 2);
    bamIO.header.records[0].type = seqan::BAM_HEADER_FIRST;
    bamIO.header.records[0].tags[0].i1 = "VN";
    bamIO.header.records[0].tags[0].i2 = "1.3";
    bamIO.header.records[0].tags[1].i1 = "SO";
    bamIO.header.records[0].tags[1].i2 = "coordinate";
    resize(bamIO.header.records[1].tags, 2);
    bamIO.header.records[1].type = seqan::BAM_HEADER_REFERENCE;
    bamIO.header.records[1].tags[0].i1 = "SN";
    bamIO.header.records[1].tags[0].i2 = "REFERENCE";
    bamIO.header.records[1].tags[1].i1 = "LN";
    bamIO.header.records[1].tags[1].i2 = "10000";

    // Force writing of header on flush.
    flush(bamIO);
    close(bamIO);

    // Compare results.
    if (seqan::endsWith(pathFragmentExpected, ".bam"))
        SEQAN_ASSERT(seqan::_compareBinaryFiles(toCString(tmpPath), toCString(filePath)));
    else
        SEQAN_ASSERT(seqan::_compareTextFiles(toCString(tmpPath), toCString(filePath)));
}

SEQAN_DEFINE_TEST(test_bam_io_bam_stream_sam_write_header)
{
    testBamIOBamStreamWriteHeader("/extras/tests/bam_io/header.sam");
}

SEQAN_DEFINE_TEST(test_bam_io_bam_stream_bam_write_header)
{
    testBamIOBamStreamWriteHeader("/extras/tests/bam_io/header.bam");
}

// ---------------------------------------------------------------------------
// Write Records
// ---------------------------------------------------------------------------

void testBamIOBamStreamWriteRecords(char const * pathFragmentExpected)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, pathFragmentExpected);

    seqan::CharString tmpPath = SEQAN_TEMP_FILENAME();
    if (seqan::endsWith(pathFragmentExpected, ".bam"))
        append(tmpPath, ".bam");
    else
        append(tmpPath, ".sam");

    // Initialize BamStream, build header.
    seqan::BamStream bamIO(toCString(tmpPath), seqan::BamStream::WRITE);
    resize(bamIO.header.sequenceInfos, 1);
    bamIO.header.sequenceInfos[0].i1 = "REFERENCE";
    bamIO.header.sequenceInfos[0].i2 = 10000;
    resize(bamIO.header.records, 2);
    resize(bamIO.header.records[0].tags, 2);
    bamIO.header.records[0].type = seqan::BAM_HEADER_FIRST;
    bamIO.header.records[0].tags[0].i1 = "VN";
    bamIO.header.records[0].tags[0].i2 = "1.3";
    bamIO.header.records[0].tags[1].i1 = "SO";
    bamIO.header.records[0].tags[1].i2 = "coordinate";
    resize(bamIO.header.records[1].tags, 2);
    bamIO.header.records[1].type = seqan::BAM_HEADER_REFERENCE;
    bamIO.header.records[1].tags[0].i1 = "SN";
    bamIO.header.records[1].tags[0].i2 = "REFERENCE";
    bamIO.header.records[1].tags[1].i1 = "LN";
    bamIO.header.records[1].tags[1].i2 = "10000";

    // Construct first records.
    seqan::BamAlignmentRecord record;

    record.qName = "READ0";
    record.flag = 2;
    record.rId = 0;
    record.pos = 0;
    record.mapQ = 8;
    resize(record.cigar, 3);
    record.cigar[0].count = 5;
    record.cigar[0].operation = 'M';
    record.cigar[1].count = 1;
    record.cigar[1].operation = 'I';
    record.cigar[2].count = 4;
    record.cigar[2].operation = 'M';
    record.rNextId = 0;
    record.pNext = 30;
    record.tLen = 40;
    record.seq = "AAAAAAAAAA";
    record.qual = "!!!!!!!!!!";
    SEQAN_ASSERT_EQ(writeRecord(bamIO, record), 0);

    record.qName = "READ0";
    record.flag = 1;
    record.rId = 0;
    record.pos = 1;
    record.mapQ = 8;
    resize(record.cigar, 3);
    record.cigar[0].count = 5;
    record.cigar[0].operation = 'M';
    record.cigar[1].count = 1;
    record.cigar[1].operation = 'I';
    record.cigar[2].count = 4;
    record.cigar[2].operation = 'M';
    record.rNextId = 0;
    record.pNext = 30;
    record.tLen = 40;
    record.seq = "AAAAAAAAAA";
    record.qual = "!!!!!!!!!!";
    SEQAN_ASSERT_EQ(writeRecord(bamIO, record), 0);

    record.qName = "READ0";
    record.flag = 3;
    record.rId = 0;
    record.pos = 2;
    record.mapQ = 8;
    resize(record.cigar, 3);
    record.cigar[0].count = 5;
    record.cigar[0].operation = 'M';
    record.cigar[1].count = 1;
    record.cigar[1].operation = 'I';
    record.cigar[2].count = 4;
    record.cigar[2].operation = 'M';
    record.rNextId = seqan::BamAlignmentRecord::INVALID_REFID;
    record.pNext = seqan::BamAlignmentRecord::INVALID_POS;
    record.tLen = seqan::BamAlignmentRecord::INVALID_LEN;
    record.seq = "AAAAAAAAAA";
    record.qual = "!!!!!!!!!!";
    SEQAN_ASSERT_EQ(writeRecord(bamIO, record), 0);

    // Force writing of everything.
    close(bamIO);

    // Compare results.
    if (seqan::endsWith(pathFragmentExpected, ".bam"))
        SEQAN_ASSERT(seqan::_compareBinaryFiles(toCString(tmpPath), toCString(filePath)));
    else
        SEQAN_ASSERT(seqan::_compareTextFiles(toCString(tmpPath), toCString(filePath)));
}

SEQAN_DEFINE_TEST(test_bam_io_bam_stream_sam_write_records)
{
    testBamIOBamStreamWriteRecords("/extras/tests/bam_io/small.sam");
}

SEQAN_DEFINE_TEST(test_bam_io_bam_stream_bam_write_records)
{
    testBamIOBamStreamWriteRecords("/extras/tests/bam_io/small.bam");
}

#endif  // #ifndef EXTRAS_TESTS_BAM_IO_TEST_EASY_BAM_IO_H_
