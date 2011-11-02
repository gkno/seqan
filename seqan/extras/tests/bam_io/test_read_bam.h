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

#ifndef EXTRAS_TESTS_BAM_IO_TEST_READ_BAM_H_
#define EXTRAS_TESTS_BAM_IO_TEST_READ_BAM_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/bam_io.h>

SEQAN_DEFINE_TEST(test_bam_io_bam_read_header)
{
    using namespace seqan;

    // -----------------------------------------------------------------------
    // Define constant data and input.
    // -----------------------------------------------------------------------

    // File has same contents as in the SAM test.
    CharString bamFilename;
    append(bamFilename, SEQAN_PATH_TO_ROOT());
    append(bamFilename, "/extras/tests/bam_io/small.bam");

    Stream<Bgzf> stream;
    SEQAN_ASSERT(open(stream, toCString(bamFilename), "r"));

    // -----------------------------------------------------------------------
    // Call Code Under Test.
    // -----------------------------------------------------------------------

    StringSet<CharString> referenceNameStore;
    NameStoreCache<StringSet<CharString> > referenceNameStoreCache(referenceNameStore);
    BamIOContext<StringSet<CharString> > bamIOContext(referenceNameStore, referenceNameStoreCache);
    
    BamHeader header;
    SEQAN_ASSERT_EQ(readRecord(header, bamIOContext, stream, Bam()), 0);

    // -----------------------------------------------------------------------
    // Check Results.
    // -----------------------------------------------------------------------

    SEQAN_ASSERT_EQ(length(header.sequenceInfos), 1u);
    SEQAN_ASSERT_EQ(header.sequenceInfos[0].i1, "REFERENCE");
    SEQAN_ASSERT_EQ(header.sequenceInfos[0].i2, 10000);

    SEQAN_ASSERT_EQ(length(referenceNameStore), 1u);
    SEQAN_ASSERT_EQ(referenceNameStore[0], "REFERENCE");

    SEQAN_ASSERT_EQ(length(header.records), 2u);

    SEQAN_ASSERT_EQ(header.records[0].type, BAM_HEADER_FIRST);
    SEQAN_ASSERT_EQ(length(header.records[0].tags), 2u);
    SEQAN_ASSERT_EQ(header.records[0].tags[0].i1, "VN");
    SEQAN_ASSERT_EQ(header.records[0].tags[0].i2, "1.3");
    SEQAN_ASSERT_EQ(header.records[0].tags[1].i1, "SO");
    SEQAN_ASSERT_EQ(header.records[0].tags[1].i2, "coordinate");

    SEQAN_ASSERT_EQ(header.records[1].type, BAM_HEADER_REFERENCE);
    SEQAN_ASSERT_EQ(length(header.records[1].tags), 2u);
    SEQAN_ASSERT_EQ(header.records[1].tags[0].i1, "SN");
    SEQAN_ASSERT_EQ(header.records[1].tags[0].i2, "REFERENCE");
    SEQAN_ASSERT_EQ(header.records[1].tags[1].i1, "LN");
    SEQAN_ASSERT_EQ(header.records[1].tags[1].i2, "10000");
}

SEQAN_DEFINE_TEST(test_bam_io_bam_read_alignment)
{
    using namespace seqan;

    // -----------------------------------------------------------------------
    // Define constant data and input.
    // -----------------------------------------------------------------------

    // File has same contents as in the SAM test.
    CharString bamFilename;
    append(bamFilename, SEQAN_PATH_TO_ROOT());
    append(bamFilename, "/extras/tests/bam_io/small.bam");

    Stream<Bgzf> stream;
    SEQAN_ASSERT(open(stream, toCString(bamFilename), "r"));

    // -----------------------------------------------------------------------
    // Call Code Under Test.
    // -----------------------------------------------------------------------

    StringSet<CharString> referenceNameStore;
    NameStoreCache<StringSet<CharString> > referenceNameStoreCache(referenceNameStore);
    BamIOContext<StringSet<CharString> > bamIOContext(referenceNameStore, referenceNameStoreCache);
    
    BamHeader header;
    SEQAN_ASSERT_EQ(readRecord(header, bamIOContext, stream, Bam()), 0);

    String<BamAlignmentRecord> alignments;
    while (!atEnd(stream))
    {
        resize(alignments, length(alignments) + 1);
        SEQAN_ASSERT_EQ(readRecord(back(alignments), bamIOContext, stream, Bam()), 0);
    }

    // -----------------------------------------------------------------------
    // Check Results.
    // -----------------------------------------------------------------------

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

#endif  // EXTRAS_TESTS_BAM_IO_TEST_READ_BAM_H_
