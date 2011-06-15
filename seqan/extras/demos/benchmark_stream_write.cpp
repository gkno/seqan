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
// Very simple benchmarking tool for stream writing
// ==========================================================================

#include <cstdio>
#include <fstream>
#if SEQAN_HAS_ZLIB
#include <zlib.h>
#endif  // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
#include <bzlib.h>
#endif  // #if SEQAN_HAS_BZIP2

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/stream.h>
#include <seqan/sequence.h>
#include <seqan/misc/misc_cmdparser.h>

// template <typename TMeta, typename TSeq>
// int writeRecordNewIO(std::fstream & stream, TMeta const & meta, TSeq const & seq, Fasta const & /**/)
// {
//     
// 
// 
// }

using namespace seqan;

void constructFastaStrings(StringSet<CharString> & metas, StringSet<CharString> & seqs)
{
    resize(metas, 4);
    resize(seqs, 4);

    reserve(seqs[0], 1024*1024);
    metas[0] = "1MB of as";
    for (int i = 0; i< 1024*1024; ++i)
        appendValue(seqs[0], 'a');

    reserve(seqs[1], 1024*1024*32);
    metas[1] = "32MB of as";
    for (int i = 0; i< 1024*1024*32; ++i)
        appendValue(seqs[1], 'c');

    reserve(seqs[2], 1024*1024*256);
    metas[2] = "256MB of gs";
    for (int i = 0; i< 1024*1024*256; ++i)
        appendValue(seqs[2], 'g');

    reserve(seqs[3], 1024*1024*512);
    metas[3] = "512MB of ts";
    for (int i = 0; i< 1024*1024*512; ++i)
        appendValue(seqs[3], 't');
}

int main(int argc, char const ** argv)
{
    (void)argc;
    (void)argv;

    StringSet<CharString> metas;
    StringSet<CharString> seqs;

    double before = sysTime();
    std::cerr << "Constructing Strings ...." << std::flush;
    constructFastaStrings(metas, seqs);
    double after = sysTime();
    std::cerr << "completed in " << after - before << "s\n" << std::flush;

    {
        CharString tempFilename = SEQAN_TEMP_FILENAME();
        appendValue(tempFilename, '0');
        tempFilename = "0.test";
        char filenameBuffer[1000];
        strncpy(filenameBuffer, toCString(tempFilename), 999);

        std::fstream file(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(file.is_open());

        before = sysTime();
        std::cerr << "Writing with new IO (no linebreaks, new streamPut) ...."<< std::flush;
        for (int i = 0; i < 4; ++i)
            writeRecord0(file, metas[i], seqs[i], Fasta());
        after = sysTime();
        std::cerr << "completed in " << after - before << "s\n"<< std::flush;
        file.close();
    }
    {
        CharString tempFilename = SEQAN_TEMP_FILENAME();
        appendValue(tempFilename, '1');
        tempFilename = "1.test";
        char filenameBuffer[1000];
        strncpy(filenameBuffer, toCString(tempFilename), 999);

        std::fstream file(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(file.is_open());

        before = sysTime();
        std::cerr << "Writing with new IO (line-by-line, new streamPut) ...."<< std::flush;
        for (int i = 0; i < 4; ++i)
            writeRecord(file, metas[i], seqs[i], Fasta());
        after = sysTime();
        std::cerr << "completed in " << after - before << "s\n"<< std::flush;
        file.close();
    }
    {
        CharString tempFilename = SEQAN_TEMP_FILENAME();
        appendValue(tempFilename, '2');
        tempFilename = "2.test";
        char filenameBuffer[1000];
        strncpy(filenameBuffer, toCString(tempFilename), 999);

        std::fstream file(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(file.is_open());

        before = sysTime();
        std::cerr << "Writing with new IO (char-by-char, new streamPut) ...."<< std::flush;
        for (int i = 0; i < 4; ++i)
            writeRecord2(file, metas[i], seqs[i], Fasta());
        after = sysTime();
        std::cerr << "completed in " << after - before << "s\n"<< std::flush;
        file.close();
    }
    {
        CharString tempFilename = SEQAN_TEMP_FILENAME();
        appendValue(tempFilename, '3');
        tempFilename = "3.test";
        char filenameBuffer[1000];
        strncpy(filenameBuffer, toCString(tempFilename), 999);

        std::fstream file(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(file.is_open());

        before = sysTime();
        std::cerr << "Writing with new IO (line-by-line, old streamPut) ...."<< std::flush;
        for (int i = 0; i < 4; ++i)
            writeRecord3(file, metas[i], seqs[i], Fasta());
        after = sysTime();
        std::cerr << "completed in " << after - before << "s\n"<< std::flush;
        file.close();
    }
    {
        CharString tempFilename = SEQAN_TEMP_FILENAME();
        appendValue(tempFilename, '4');
        tempFilename = "4.test";
        char filenameBuffer[1000];
        strncpy(filenameBuffer, toCString(tempFilename), 999);

        std::fstream file(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(file.is_open());

        before = sysTime();
        std::cerr << "Writing with new IO (char-by-char, old streamPut) ...."<< std::flush;
        for (int i = 0; i < 4; ++i)
            writeRecord4(file, metas[i], seqs[i], Fasta());
        after = sysTime();
        std::cerr << "completed in " << after - before << "s\n"<< std::flush;
        file.close();
    }
    {
        CharString tempFilename = SEQAN_TEMP_FILENAME();
        appendValue(tempFilename, '5');
        tempFilename = "5.test";
        char filenameBuffer[1000];
        strncpy(filenameBuffer, toCString(tempFilename), 999);

        std::fstream file(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(file.is_open());

        before = sysTime();
        std::cerr << "Writing with new IO (char-by-char, new streamPut, iterators instead of index) ...."<< std::flush;
        for (int i = 0; i < 4; ++i)
            writeRecord5(file, metas[i], seqs[i], Fasta());
        after = sysTime();
        std::cerr << "completed in " << after - before << "s\n"<< std::flush;
        file.close();
    }
    {
        CharString tempFilename = SEQAN_TEMP_FILENAME();
        appendValue(tempFilename, '6');
        tempFilename = "6.test";
        char filenameBuffer[1000];
        strncpy(filenameBuffer, toCString(tempFilename), 999);

        std::fstream file(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(file.is_open());

        before = sysTime();
        std::cerr << "Writing with new IO (char-by-char, old streamPut, iterators instead of index) ...."<< std::flush;
        for (int i = 0; i < 4; ++i)
            writeRecord6(file, metas[i], seqs[i], Fasta());
        after = sysTime();
        std::cerr << "completed in " << after - before << "s\n"<< std::flush;
        file.close();
    }
    {
        CharString tempFilename = SEQAN_TEMP_FILENAME();
        appendValue(tempFilename, '7');
        tempFilename = "7.test";
        char filenameBuffer[1000];
        strncpy(filenameBuffer, toCString(tempFilename), 999);

        std::fstream file(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(file.is_open());

        before = sysTime();
        std::cerr << "Writing with old IO...."<< std::flush;
        for (int i = 0; i < 4; ++i)
            write(file, seqs[i], metas[i], Fasta());
        after = sysTime();
        std::cerr << "completed in " << after - before << "s\n"<< std::flush;

        file.close();
        return 0;
    }
}
