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
// Tests for seqan/stream/tokenize.h
// ==========================================================================

#ifndef TEST_STREAM_TEST_STREAM_TOKENIZING_H_
#define TEST_STREAM_TEST_STREAM_TOKENIZING_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "test_stream_generic.h"

std::fstream* createFile()
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream *file = new std::fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(file->is_open());

    char const * STR = "This is a string ....foobar AAAACCCGGGTTT\n\
    TCG > foo.. SJAUDOF78456KAPLP345LPL AAAAAAAAAAAAA\n\
    \a etstststetststetststetsstetetstetststetstdtetsteststetstedetstet\r\n\
    etstetetstststetststetststetsstetetstetststetstdtetstestst    \r\n\
    etstetetstststetststetstste           tststetstdtetsteststetstedetstet\r\n\
    etstetetstststetststetststetsstetetstetststetstdtetsteststetstedetstet\r\n\
    \n\v\n\r\nACGTACGTACGTACGATACGATCTnn\n\nACGT";
    

    file->write(STR, strlen(STR));
    file->seekg(0);
    file->seekp(0);
    return file;
}


SEQAN_DEFINE_TEST(test_stream_tokenizing)
{
    using namespace seqan;
    std::fstream * file = createFile();
    RecordReader<std::fstream, SinglePass<void> > reader(*file);

    
    CharString buf;
    SEQAN_ASSERT_EQ(readUntilWhitespace(buf,reader), 0);

    SEQAN_ASSERT_EQ(value(reader), ' ');

    SEQAN_ASSERT_EQ(readUntilChar(buf,reader, 'f'), 0);

    SEQAN_ASSERT_EQ(buf, " is a string ....");

    SEQAN_ASSERT_EQ(readNChars(buf,reader, 3), 0);
    
    SEQAN_ASSERT_EQ(buf, "foo");

    SEQAN_ASSERT_EQ(skipUntilBlank(reader), 0);

    SEQAN_ASSERT_EQ(skipUntilWhitespace(reader), 0); // will stay in same place
    
    SEQAN_ASSERT_EQ(readDna5IgnoringWhitespaces(buf, reader), 0);
    SEQAN_ASSERT_EQ(buf, "AAAACCCGGGTTTTCG");
    SEQAN_ASSERT_EQ(value(reader), '>');

    SEQAN_ASSERT_EQ(skipUntilChar(reader, '.'), 0);
    SEQAN_ASSERT_EQ(value(reader), '.');

    SEQAN_ASSERT_EQ(skipUntilChar(reader, 'S'), 0);

    SEQAN_ASSERT_EQ(readLetters(buf, reader), 0);
    SEQAN_ASSERT_EQ(buf, "SJAUDOF");

    SEQAN_ASSERT_EQ(skipUntilString(reader, "KAP"), 0); // -we fail here-
    SEQAN_ASSERT_EQ(readAlphaNums(buf, reader), 0);
    SEQAN_ASSERT_EQ(buf, "LP345LPL");

    SEQAN_ASSERT_EQ(readDna5IgnoringWhitespaces(buf, reader), 0);
    SEQAN_ASSERT_EQ(buf, "AAAAAAAAAAAAA");
    SEQAN_ASSERT_EQ(value(reader), '\a'); // Bell character is not whitespace
    
    SEQAN_ASSERT_EQ(readLine(buf, reader), 0); // do not add \r\n
    SEQAN_ASSERT_EQ(buf, "etstststetststetststetsstetetstetststetstdtetsteststetstedetstet");
    SEQAN_ASSERT_EQ(value(reader), ' ');
    

    //TODO investigate failure
    //TODO investigate compiler warnings
    //TODO(h4nn3s) FINISH


    
}









#endif // ndef TEST_STREAM_TEST_STREAM_TOKENIZING_H_
