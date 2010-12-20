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
// ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Test the find2/find_exact_simple.h header.
// ==========================================================================

#ifndef TESTS_FIND2_TEST_FIND_APPROX_DPSEARCH_H_
#define TESTS_FIND2_TEST_FIND_APPROX_DPSEARCH_H_

#include <seqan/basic.h>
#include <seqan/find2.h>

using namespace seqan;

// Test the basic interface for the DPSearch pattern class.
SEQAN_DEFINE_TEST(test_find2_find_approx_dpsearch_pattern_interface) {
    // TODO(holtgrew): Strings should be const.
    CharString kHaystack = "he is a hero";
    CharString kNeedle = "her";

    typedef Finder<CharString> TFinder;
    typedef Pattern<CharString, DPSearch<EditDistanceScore> > TPattern;

    // Default constructor.
    {
        TPattern pattern;
    }
    // Score limit but not scoring scheme.
    {
        TPattern pattern(kNeedle, -1);
    }

    // Construct pattern and finder that are used below.
    EditDistanceScore scoringScheme;
    TPattern pattern(kNeedle, -1, scoringScheme);
    TFinder finder(kHaystack);

    // Function needle().
    CharString const & patternNeedle = needle(pattern);
    SEQAN_ASSERT_EQ(kNeedle, patternNeedle);
    SEQAN_ASSERT_EQ(3u, length(patternNeedle));

    // Function host().
    CharString const & patternHost = host(pattern);
    SEQAN_ASSERT_EQ(kNeedle, patternHost);
    SEQAN_ASSERT_EQ(3u, length(patternHost));

    // Function find().
    bool ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);

    // Function findBegin().
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);

    // Resulting alignment should be.
    // he-is a hero
    // ||
    // her

    // Function alignment().
    Align<CharString, ArrayGaps> align;
    ret = buildAlignment(finder, pattern, align);
    SEQAN_ASSERT_TRUE(ret);
    // Test the alignment by testing the view-to-source conversion.
    // First row.
    SEQAN_ASSERT_EQ(0u, toViewPosition(row(align, 0), 0));
    SEQAN_ASSERT_EQ(1u, toViewPosition(row(align, 0), 1));
    SEQAN_ASSERT_EQ(3u, toViewPosition(row(align, 0), 2));
    SEQAN_ASSERT_EQ(4u, toViewPosition(row(align, 0), 3));
    SEQAN_ASSERT_EQ(5u, toViewPosition(row(align, 0), 4));
    SEQAN_ASSERT_EQ(6u, toViewPosition(row(align, 0), 5));
    SEQAN_ASSERT_EQ(7u, toViewPosition(row(align, 0), 6));
    SEQAN_ASSERT_EQ(8u, toViewPosition(row(align, 0), 7));
    SEQAN_ASSERT_EQ(9u, toViewPosition(row(align, 0), 8));
    SEQAN_ASSERT_EQ(10u, toViewPosition(row(align, 0), 9));
    SEQAN_ASSERT_EQ(11u, toViewPosition(row(align, 0), 10));
    SEQAN_ASSERT_EQ(12u, toViewPosition(row(align, 0), 11));
    // Second row.
    SEQAN_ASSERT_EQ(0u, toViewPosition(row(align, 1), 0));
    SEQAN_ASSERT_EQ(1u, toViewPosition(row(align, 1), 1));
    SEQAN_ASSERT_EQ(2u, toViewPosition(row(align, 1), 2));

    // Function length().
    SEQAN_ASSERT_EQ(3u, length(pattern));

    // Function infix() for pattern.
    CharString const kExpectedInfixStr = "her";
    CharString const kInfixStr = infix(pattern);
    SEQAN_ASSERT_EQ(kExpectedInfixStr, kInfixStr);

    // Function begin() for pattern.
    {
        typedef Iterator<TPattern::TNeedle const, Standard>::Type TConstIterator;
        TConstIterator kExpectedBegin = begin(kNeedle, Standard());
        TConstIterator kBegin = begin(pattern, Standard());
        SEQAN_ASSERT_EQ(kExpectedBegin, kBegin);
    }

    // Function end() for pattern.
    {
        typedef Iterator<TPattern::TNeedle const, Standard>::Type TConstIterator;
        TConstIterator kExpectedEnd = end(kNeedle, Standard());
        TConstIterator kEnd = end(pattern, Standard());
        SEQAN_ASSERT_EQ(kExpectedEnd, kEnd);
    }

    // Function endPosition(pattern).
    typedef Position<CharString>::Type TPosition;
    TPosition patternEndPosition = endPosition(pattern);
    SEQAN_ASSERT_EQ(3u, patternEndPosition);

    // Function beginPosition(pattern).
    TPosition patternBeginPosition = beginPosition(pattern);
    SEQAN_ASSERT_EQ(0u, patternBeginPosition);
}


// "Easy" (= short, few hits) test of find() using the exact pattern.
SEQAN_DEFINE_TEST(test_find2_find_approx_dpsearch_pattern_find_easy_score_limit_0) {
    typedef Finder<CharString> TFinder;
    typedef Pattern<CharString, DPSearch<EditDistanceScore> > TPattern;

    CharString kHaystack = "he is her hero";
    CharString kNeedle = "he";
    TFinder finder(kHaystack);
    TPattern pattern(kNeedle, 0);

    // We will search for all occurences of the needle using find()
    // and findBegin() in the haystack and check for the correct begin
    // and end position in both finder and needle.  Testing the
    // alignment each time seems overkill.
    bool ret;

    // he is her hero
    // he
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(2u, endPosition(finder));
    SEQAN_ASSERT_EQ(2u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // he is her hero
    //       he
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(8u, endPosition(finder));
    SEQAN_ASSERT_EQ(2u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(6u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // he is her hero
    //           he
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(12u, endPosition(finder));
    SEQAN_ASSERT_EQ(2u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // no more hit
    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// "Easy" (= short, few hits) test of find() using the exact pattern.
// The score limit is set to -1.
SEQAN_DEFINE_TEST(test_find2_find_approx_dpsearch_pattern_find_easy_score_limit_1) {
    typedef Finder<CharString> TFinder;
    typedef Pattern<CharString, DPSearch<EditDistanceScore> > TPattern;

    CharString kHaystack = "he is her hero";
    CharString kNeedle = "her";
    TFinder finder(kHaystack);
    TPattern pattern(kNeedle, -1);

    // We will search for all occurences of the needle using find()
    // and findBegin() in the haystack and check for the correct begin
    // and end position in both finder and needle.  Testing the
    // alignment each time seems overkill.
    bool ret;

    // he- is her hero
    // ||
    // her
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(2u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // he is her hero
    // ||
    // her
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(3u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // he is he-r hero
    //       |||
    //       her
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(8u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(6u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // he is her hero
    //       |||
    //       her
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(9u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(6u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // he is her hero
    //       |||
    //       her-
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(10u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(6u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // he is her he-ro
    //           ||
    //           her
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(12u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // he is her hero
    //           |||
    //           her
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(13u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // he is her hero
    //           |||
    //           her-
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(14u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // no more hit
    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// "Easy" (= short, few hits) test of find() using the exact pattern.
// Prefix search variant.
SEQAN_DEFINE_TEST(test_find2_find_approx_dpsearch_pattern_find_easy_score_limit_0_prefix) {
    typedef Finder<CharString> TFinder;
    typedef Pattern<CharString, DPSearch<EditDistanceScore, FindPrefix> > TPattern;

    {
        CharString kHaystack = "he is her hero";
        CharString kNeedle = "he";
        TFinder finder(kHaystack);
        TPattern pattern(kNeedle, 0);

        bool ret;

        // he is her hero
        // he
        ret = find(finder, pattern);
        SEQAN_ASSERT_TRUE(ret);
        SEQAN_ASSERT_EQ(2u, endPosition(finder));
        SEQAN_ASSERT_EQ(2u, endPosition(pattern));
        ret = findBegin(finder, pattern);
        SEQAN_ASSERT_TRUE(ret);
        SEQAN_ASSERT_EQ(0u, beginPosition(finder));
        SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
        ret = findBegin(finder, pattern);
        SEQAN_ASSERT_NOT(ret);

        ret = find(finder, pattern);
        SEQAN_ASSERT_NOT(ret);
    }

    {
        CharString kHaystack = "Superman is her hero";
        CharString kNeedle = "he";
        TFinder finder(kHaystack);
        TPattern pattern(kNeedle, 0);

        bool ret = find(finder, pattern);
        SEQAN_ASSERT_NOT(ret);
    }
}


// "Easy" (= short, few hits) test of find() using the exact pattern.
// The score limit is set to -1.  Prefix search variant.
SEQAN_DEFINE_TEST(test_find2_find_approx_dpsearch_pattern_find_easy_score_limit_1_prefix) {
    typedef Finder<CharString> TFinder;
    typedef Pattern<CharString, DPSearch<EditDistanceScore, FindPrefix> > TPattern;

    {
        CharString kHaystack = "he is her hero";
        CharString kNeedle = "her";
        TFinder finder(kHaystack);
        TPattern pattern(kNeedle, -1);

        bool ret;

        // he- is her hero
        // ||
        // her
        ret = find(finder, pattern);
        SEQAN_ASSERT_TRUE(ret);
        SEQAN_ASSERT_EQ(-1, score(pattern));
        SEQAN_ASSERT_EQ(2u, endPosition(finder));
        SEQAN_ASSERT_EQ(3u, endPosition(pattern));
        ret = findBegin(finder, pattern);
        SEQAN_ASSERT_TRUE(ret);
        SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
        SEQAN_ASSERT_EQ(0u, beginPosition(finder));
        SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
        ret = findBegin(finder, pattern);
        SEQAN_ASSERT_NOT(ret);

        // he is her hero
        // ||
        // her
        ret = find(finder, pattern);
        SEQAN_ASSERT_TRUE(ret);
        SEQAN_ASSERT_EQ(-1, score(pattern));
        SEQAN_ASSERT_EQ(3u, endPosition(finder));
        SEQAN_ASSERT_EQ(3u, endPosition(pattern));
        ret = findBegin(finder, pattern);
        SEQAN_ASSERT_TRUE(ret);
        SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
        SEQAN_ASSERT_EQ(0u, beginPosition(finder));
        SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
        ret = findBegin(finder, pattern);
        SEQAN_ASSERT_NOT(ret);

        ret = find(finder, pattern);
        SEQAN_ASSERT_NOT(ret);
    }
}


// "Easy" (= short, few hits) test of find() using the exact pattern.
// The score limit is set to -1.  We use the UseScoreLimit variant of FindBegin.
SEQAN_DEFINE_TEST(test_find2_find_approx_dpsearch_pattern_find_easy_score_limit_1_use_score_limit) {
    typedef Finder<CharString> TFinder;
    typedef Pattern<CharString, DPSearch<EditDistanceScore> > TPattern;

    CharString kHaystack = "he is her hero";
    CharString kNeedle = "her";
    TFinder finder(kHaystack);
    TPattern pattern(kNeedle, -1);

    // We will search for all occurences of the needle using find()
    // and findBegin() in the haystack and check for the correct begin
    // and end position in both finder and needle.  Testing the
    // alignment each time seems overkill.
    bool ret;

    // he- is her hero
    // ||
    // her
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(2u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_NOT(ret);

    // he is her hero
    // ||
    // her
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(3u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_NOT(ret);

    // he is he-r hero
    //       |||
    //       her
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(8u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(6u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_NOT(ret);

    // he is her hero
    //       her
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(9u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    // he is her hero
    //        ||
    //        er
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(7u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    // he is her hero
    //       |||
    //       her
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(6u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    // he is her hero
    //       |||
    //      -her
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(5u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_NOT(ret);

    // he is her hero
    //       |||
    //       her-
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(10u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(6u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_NOT(ret);

    // he is her he-ro
    //           ||
    //           her
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(12u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_NOT(ret);

    // he is her hero
    //           her
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(13u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    // he is her hero
    //            ||
    //            er
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(11u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    // he is her hero
    //           |||
    //           her
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    // he is her hero
    //           |||
    //          -her
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(9u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_NOT(ret);

    // he is her hero
    //           |||
    //           her-
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(14u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern, UseScoreLimit());
    SEQAN_ASSERT_NOT(ret);

    // no more hit
    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// "Harder" (= longer, more hits hits) test of find() using the exact
// pattern.  Score limit is set to 0.
SEQAN_DEFINE_TEST(test_find2_find_approx_dpsearch_pattern_find_harder_score_limit_0) {
    typedef Finder<DnaString> TFinder;
    typedef Pattern<DnaString, DPSearch<EditDistanceScore> > TPattern;

    DnaString kHaystack = "AGAAGAAGAGGAAGAAGA";
    DnaString kNeedle = "GAA";
    TFinder finder(kHaystack);
    TPattern pattern(kNeedle, 0);

    // We will search for all occurences of the needle using find()
    // and findBegin() in the haystack and check for the correct begin
    // and end position in both finder and needle.  Testing the
    // alignment each time seems overkill.
    bool ret;

    // AGAAGAAGAGGAAGAAGA
    //  GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // AGAAGAAGAGGAAGAAGA
    //     GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // AGAAGAAGAGGAAGAAGA
    //           GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(13u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // AGAAGAAGAGGAAGAAGA
    //              GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(16u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(13u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // no more hit
    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// "Harder" (= longer, more hits hits) test of find() using the exact
// pattern.  Score limit is set to -1.
SEQAN_DEFINE_TEST(test_find2_find_approx_dpsearch_pattern_find_harder_score_limit_1) {
    typedef Finder<DnaString> TFinder;
    typedef Pattern<DnaString, DPSearch<EditDistanceScore> > TPattern;

    DnaString kHaystack = "AAAGAAGAGGAAGAAGA";
    DnaString kNeedle = "AAA";
    TFinder finder(kHaystack);
    TPattern pattern(kNeedle, -1);

    // We will search for all occurences of the needle using find()
    // and findBegin() in the haystack and check for the correct begin
    // and end position in both finder and needle.  Testing the
    // alignment each time seems overkill.
    bool ret;

    // -AAAGAAGAGGAAGAAGA
    //  ||
    // AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(2u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // AAAGAAGAGGAAGAAGA
    // |||
    // AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(3u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // AAAGAAGAGGAAGAAGA
    //  ||
    //  AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    // AAAGAAGAGGAAGAAGA
    // |||
    // AAA-
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // AAAGAAGAGGAAGAAGA
    //   | |
    //   AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(5u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(2u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    // AAAGAAGAGGAAGAAGA
    //  || |
    //  AA-A
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // AAAGAAGAGGAAGAAGA
    //     ||
    //     AA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(6u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    // AAAGAAGAGGAAGAAGA
    //     ||
    //    AAA
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(3u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    // AAAGAAGAGGAAGAAGA
    //   | ||
    //   A-AA
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(2u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // AAAGAAGAGGAAGAAGA
    //     ||
    //     AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // AAAGAAGAGGAAGAAGA
    //      | |
    //      AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(8u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(5u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    // AAAGAAGAGGAAGAAGA
    //     || |
    //     AA-A
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // AAAGAAGAGGAAGAAGA
    //           ||
    //           AA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(12u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    // AAAGAAGAGGAAGAAGA
    //           ||
    //          AAA
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(9u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // AAAGAAGAGGAAGAAGA
    //           ||
    //           AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(13u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // AAAGAAGAGGAAGAAGA
    //            | |
    //            AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(14u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    // TODO(holtgrew): Why not indel?
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(11u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    // AAAGAAGAGGAAGAAGA
    //           || |
    //           AA-A
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // AAAGAAGAGGAAG-AAGA
    //               ||
    //              AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(15u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(13u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    // AAAGAAGAGGAAGAAGA
    //              ||
    //             AAA
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(12u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    // AAAGAAGAGGAAGAAGA
    //            | ||
    //            A-AA
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(11u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // AAAGAAGAGGAAGAAGA
    //              ||
    //              AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(16u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(13u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // AAAGAAGAGGAAGAAGA
    //               | |
    //               AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(17u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(14u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    // AAAGAAGAGGAAGAAGA
    //              || |
    //              AA-A
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(13u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // no more hit
    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// Tests for setEndPosition() with the DPSearch pattern and a score limit of 0.
SEQAN_DEFINE_TEST(test_find2_find_approx_dpsearch_pattern_set_end_position_score_limit_0) {
    typedef Finder<DnaString> TFinder;
    typedef Pattern<DnaString, DPSearch<EditDistanceScore> > TPattern;

    DnaString kHaystack = "AGAAGAAGAGGAAGAAGA";
    DnaString kNeedle = "GAA";
    TFinder finder(kHaystack);
    TPattern pattern(kNeedle, 0);

    bool ret;

    // Set end position to a hit.
    // AGAAGAAGAGGAAGAAGA
    //  GAA
    ret = setEndPosition(finder, pattern, 4);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // Continue to search from here.
    // AGAAGAAGAGGAAGAAGA
    //     GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // Set end position to a mismatch.
    // AGAAGAAGAGGAAGAAGA
    //  GAA
    ret = setEndPosition(finder, pattern, 5);
    SEQAN_ASSERT_NOT(ret);

    // Continue to search from here.
    // AGAAGAAGAGGAAGAAGA
    //     GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // Set to end.
    // AGAAGAAGAGGAAGAAGA
    //                GAA
    ret = setEndPosition(finder, pattern, 18);
    SEQAN_ASSERT_NOT(ret);

    // No more hit.
    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// Tests for setEndPosition() with the DPSearch pattern and a score limit of -1.
SEQAN_DEFINE_TEST(test_find2_find_approx_dpsearch_pattern_set_end_position_score_limit_1) {
    typedef Finder<DnaString> TFinder;
    typedef Pattern<DnaString, DPSearch<EditDistanceScore> > TPattern;

    {
        // TODO(holtgrew): This does not create the expected output.
//         DnaString kHaystack = "AGAAGAAGAGGAAGAAGA";
//         DnaString kNeedle = "GAA";
//         TFinder finder(kHaystack);
//         TPattern pattern(kNeedle, -1);
//         Align<CharString, ArrayGaps> align;
//         while (find(finder, pattern)) {
//             while (findBegin(finder, pattern)) {
//                 std::cout << beginPosition(finder) << ", " << endPosition(finder) << std::endl;
//                 std::cout << infix(finder) << std::endl;
//                 std::cout << beginPosition(pattern) << ", " << endPosition(pattern) << std::endl;
//                 std::cout << infix(pattern) << std::endl;
//                 buildAlignment(finder, pattern, align);
//                 std::cout << align << std::endl;
//             }
//         }
    }

    DnaString kHaystack = "AGAAGAAGAGGAAGAAGA";
    DnaString kNeedle = "GAA";
    TFinder finder(kHaystack);
    TPattern pattern(kNeedle, -1);

    bool ret;

    // Set end position to a exact hit.
    // AGAAGAAGAGGAAGAAGA
    //  GAA
    ret = setEndPosition(finder, pattern, 4);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // Continue to search from here.
    // AGAAGAAGAGGAAGAAGA
    //  |||
    //  GAA-
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(5u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // Set end position to a mismatch.
    // AGAAGAAGAGGAAGAAGA
    //         GAA
    ret = setEndPosition(finder, pattern, 11);
    SEQAN_ASSERT_NOT(ret);

    // Continue to search from here.
    // AGAAGAAGAGGAAGAAGA
    //           ||
    //           GA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(12u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    // AGAAGAAGAGGAAGAAGA
    //          | |
    //          GAA
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(9u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // Set to end.
    // AGAAGAAGAGGAAGAAG-A
    //                 GAA
    ret = setEndPosition(finder, pattern, 18);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(18u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(16u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // No more hit.
    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// Tests for setEndPosition() with the DPSearch pattern and a score
// limit of 0.  Find prefix variant.
SEQAN_DEFINE_TEST(test_find2_find_approx_dpsearch_pattern_set_end_position_score_limit_0_prefix) {
    typedef Finder<DnaString> TFinder;
    typedef Pattern<DnaString, DPSearch<EditDistanceScore, FindPrefix> > TPattern;

    {
        DnaString kHaystack = "AGAAGAAGAGGAAGAAGA";
        DnaString kNeedle = "GAA";
        TFinder finder(kHaystack);
        TPattern pattern(kNeedle, 0);
        
        bool ret;

        // Set end position to a non-hit
        // AGAAGAAGAGGAAGAAGA
        //   GAA
        ret = setEndPosition(finder, pattern, 4);
        SEQAN_ASSERT_NOT(ret);

        // No more hit.
        ret = find(finder, pattern);
        SEQAN_ASSERT_NOT(ret);

        // Set end position to a non-hit
        // AGAAGAAGAGGAAGAAGA
        // GAA
        ret = setEndPosition(finder, pattern, 3);
        SEQAN_ASSERT_NOT(ret);

        // No more hit.
        ret = find(finder, pattern);
        SEQAN_ASSERT_NOT(ret);
    }

    {
        DnaString kHaystack = "AGAAGAAGAGGAAGAAGA";
        DnaString kNeedle = "AGA";
        TFinder finder(kHaystack);
        TPattern pattern(kNeedle, 0);
        
        bool ret;

        // Set end position to a hit
        // AGAAGAAGAGGAAGAAGA
        // AGA
        ret = setEndPosition(finder, pattern, 3);
        SEQAN_ASSERT_TRUE(ret);
        SEQAN_ASSERT_EQ(0, score(pattern));
        SEQAN_ASSERT_EQ(3u, endPosition(finder));
        SEQAN_ASSERT_EQ(3u, endPosition(pattern));
        ret = findBegin(finder, pattern);
        SEQAN_ASSERT_TRUE(ret);
        SEQAN_ASSERT_EQ(0, findBeginScore(pattern));
        SEQAN_ASSERT_EQ(0u, beginPosition(finder));
        SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

        // No more hit.
        ret = find(finder, pattern);
        SEQAN_ASSERT_NOT(ret);

        // Set end position to a non-hit
        // AGAAGAAGAGGAAGAAGA
        //    AGA
        ret = setEndPosition(finder, pattern, 6);
        SEQAN_ASSERT_NOT(ret);

        // No more hit.
        ret = find(finder, pattern);
        SEQAN_ASSERT_NOT(ret);
    }
}


// Tests for setEndPosition() with the DPSearch pattern and a score
// limit of -1.  Find prefix variant.
SEQAN_DEFINE_TEST(test_find2_find_approx_dpsearch_pattern_set_end_position_score_limit_1_prefix) {
    typedef Finder<DnaString> TFinder;
    typedef Pattern<DnaString, DPSearch<EditDistanceScore, FindPrefix> > TPattern;

    {
        DnaString kHaystack = "AGAAGAAGAGGAAGAAGA";
        DnaString kNeedle = "GAA";
        TFinder finder(kHaystack);
        TPattern pattern(kNeedle, -1);
        
        bool ret;

        // Set end position to a hit.
        // AGAAGAAGAGGAAGAAGA
        // -GAA
        ret = setEndPosition(finder, pattern, 4);
        SEQAN_ASSERT_TRUE(ret);
        SEQAN_ASSERT_EQ(-1, score(pattern));
        SEQAN_ASSERT_EQ(4u, endPosition(finder));
        SEQAN_ASSERT_EQ(3u, endPosition(pattern));
        ret = findBegin(finder, pattern);
        SEQAN_ASSERT_TRUE(ret);
        SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
        SEQAN_ASSERT_EQ(0u, beginPosition(finder));
        SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
        ret = findBegin(finder, pattern);
        SEQAN_ASSERT_NOT(ret);

        // Set end position to a non-hit.
        // -AGAAGAAGAGGAAGAAGA
        // GAA
        ret = setEndPosition(finder, pattern, 2);
        SEQAN_ASSERT_NOT(ret);

        // Search next.
        ret = find(finder, pattern);
        SEQAN_ASSERT_TRUE(ret);
        SEQAN_ASSERT_EQ(-1, score(pattern));
        SEQAN_ASSERT_EQ(4u, endPosition(finder));
        SEQAN_ASSERT_EQ(3u, endPosition(pattern));
        ret = findBegin(finder, pattern);
        SEQAN_ASSERT_TRUE(ret);
        SEQAN_ASSERT_EQ(-1, findBeginScore(pattern));
        SEQAN_ASSERT_EQ(0u, beginPosition(finder));
        SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
        ret = findBegin(finder, pattern);
        SEQAN_ASSERT_NOT(ret);
    }
}


// More or less extensive test for findBegin().  We test find() and
// findBegin() with two "new" alignments, i.e. distinct end positions
// and two "begin alignments" each.
//
// TODO(holtgrew): Move to test_find_approx_find_begin.h?
SEQAN_DEFINE_TEST(test_find2_find_approx_dpsearch_pattern_find_begin) {
    typedef Finder<DnaString> TFinder;
    typedef Pattern<DnaString, DPSearch<EditDistanceScore> > TPattern;

    DnaString kHaystack = "GGCACACA";
    DnaString kNeedle = "CCACA";
    TFinder finder(kHaystack);
    TPattern pattern(kNeedle, -1);

    bool ret;

    //  GGCACACA
    //    ||||
    //    CACA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(6u, endPosition(finder));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(2u, beginPosition(finder));
    //  GGCACACA
    //    ||||
    //   CCACA
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));
    // No further begin match.
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);

    // GGCACACA
    //     ||||
    //     CACA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(8u, endPosition(finder));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
    // GGCACACA
    //     ||||
    //   C-CACA
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(3u, beginPosition(finder));
    // GGCACACA
    //     ||||
    //   C-CACA
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(2u, beginPosition(finder));
    // No further begin match.
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}

#endif  // TESTS_FIND2_TEST_FIND_APPROX_DPSEARCH_H_
