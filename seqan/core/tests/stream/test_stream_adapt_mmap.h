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
// Tests for mmap stream concept adaption.
// ==========================================================================

#ifndef TEST_STREAM_TEST_STREAM_ADAPT_MMAP_H_
#define TEST_STREAM_TEST_STREAM_ADAPT_MMAP_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "test_stream_generic.h"

// Test of streamPut().
SEQAN_DEFINE_TEST(test_stream_adapt_mmap_streamPut)
{
    using namespace seqan;

    CharString tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strncpy(filenameBuffer, toCString(tempFilename), 999);

    String<char, MMap<> > mmapString;
    open(mmapString, filenameBuffer);

    testStreamPut(mmapString);

    char cmp[] = "c\nseq\nsss\n12\n34\n56\n78\n5.4\n6.5\nA\nACGT\nACGTN\n";

    SEQAN_ASSERT_EQ(CharString(mmapString), CharString(cmp));
    close(mmapString);
}

#endif  // TEST_STREAM_TEST_STREAM_ADAPT_MMAP_H_
