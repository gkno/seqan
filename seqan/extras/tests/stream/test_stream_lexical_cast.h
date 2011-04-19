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
// casts for reading different types from strings
// ==========================================================================

template <typename TTest>
void _test1(TTest const & s)
{
    using namespace seqan;
//     TTest s = "12345";

    int i       = lexicalCast<int>(s);
    short sh    = lexicalCast<short>(s);
    long l      = lexicalCast<long>(s);
    unsigned int ui = lexicalCast<unsigned int>(s);

    float f     = lexicalCast<float>(s);
    double d    = lexicalCast<double>(s);

    SEQAN_ASSERT_EQ(i,  12345);
    SEQAN_ASSERT_EQ(sh, 12345);
    SEQAN_ASSERT_EQ(l,  12345l);
    SEQAN_ASSERT_EQ(ui, 12345u);

    SEQAN_ASSERT_EQ(f,  12345.00f);
    SEQAN_ASSERT_EQ(d,  12345.00);
}

template <typename TTest>
void _test2(TTest const & s)
{
    using namespace seqan;
//     TTest s = "12345";

    int i       = lexicalCast<int>(s);
    short sh    = lexicalCast<short>(s);
    long l      = lexicalCast<long>(s);
    unsigned int ui = lexicalCast<unsigned int>(s);

    float f     = lexicalCast<float>(s);
    double d    = lexicalCast<double>(s);

    SEQAN_ASSERT_EQ(i,  -12345);
    SEQAN_ASSERT_EQ(sh, -12345);
    SEQAN_ASSERT_EQ(l,  -12345l);
    SEQAN_ASSERT_EQ(ui, UINT_MAX);

    SEQAN_ASSERT_EQ(f,  -12345.00f);
    SEQAN_ASSERT_EQ(d,  -12345.00);
}

template <typename TTest>
void _test3(TTest const & s)
{
    using namespace seqan;
//     TTest s = "-5.4";

    int i       = lexicalCast<int>(s);
    short sh    = lexicalCast<short>(s);
    long l      = lexicalCast<long>(s);
    unsigned int ui = lexicalCast<unsigned int>(s);

    float f     = lexicalCast<float>(s);
    double d    = lexicalCast<double>(s);

    SEQAN_ASSERT_EQ(i,  -5);
    SEQAN_ASSERT_EQ(sh, -5);
    SEQAN_ASSERT_EQ(l,  -5l);
    SEQAN_ASSERT_EQ(ui,  UINT_MAX);

    SEQAN_ASSERT_EQ(f,  -5.4f);
    SEQAN_ASSERT_EQ(d,  -5.4);
}



SEQAN_DEFINE_TEST(test_stream_lexical_cast_1_stdstring)
{
    using namespace seqan;

    std::string s = "12345";
    _test1(s);

    s = "-12345";
    _test2(s);
    
    s = "-5.4";
    _test3(s);
}

SEQAN_DEFINE_TEST(test_stream_lexical_cast_1_chararray)
{
    using namespace seqan;

    char s[] = "12345";
    _test1(s);

    strcpy(s, "-12345");
    _test2(s);

    strcpy(s, "-5.4");
    _test3(s);
}

SEQAN_DEFINE_TEST(test_stream_lexical_cast_1_seqanstring)
{
    using namespace seqan;

    CharString s = "12345";
    _test1(s);

    s = "-12345";
    _test2(s);

    s = "-5.4";
    _test3(s);
}
/*
SEQAN_DEFINE_TEST(test_stream_lexical_cast_1_chararray)
{

    using namespace seqan;

    char s[] = "12345";

    int i = lexicalCast<int>(s);
    short sh = lexicalCast<short>(s);
    long l = lexicalCast<long>(s);
    unsigned int ui = lexicalCast<unsigned int>(s);

    float f = lexicalCast<float>(s);
    double d = lexicalCast<double>(s);

    SEQAN_ASSERT_EQ(i, 12345);
    SEQAN_ASSERT_EQ(sh, 12345);
    SEQAN_ASSERT_EQ(l, 12345l);
    SEQAN_ASSERT_EQ(ui, 12345u);

    SEQAN_ASSERT_EQ(f, 12345.00);
    SEQAN_ASSERT_EQ(d, 12345.00);

}

SEQAN_DEFINE_TEST(test_stream_lexical_cast_1_seqanstring)
{

    using namespace seqan;

    seqan::CharString s = "12345";

    int i = lexicalCast<int>(s);
    short sh = lexicalCast<short>(s);
    long l = lexicalCast<long>(s);
    unsigned int ui = lexicalCast<unsigned int>(s);

    float f = lexicalCast<float>(s);
    double d = lexicalCast<double>(s);

    SEQAN_ASSERT_EQ(i, 12345);
    SEQAN_ASSERT_EQ(sh, 12345);
    SEQAN_ASSERT_EQ(l, 12345l);
    SEQAN_ASSERT_EQ(ui, 12345u);

    SEQAN_ASSERT_EQ(f, 12345.00);
    SEQAN_ASSERT_EQ(d, 12345.00);

}*/

