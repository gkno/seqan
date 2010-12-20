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
// Test the find2/find_hamming_simple.h header.
// ==========================================================================

#ifndef TESTS_FIND2_TEST_FIND_HAMMING_SIMPLE_H_
#define TESTS_FIND2_TEST_FIND_HAMMING_SIMPLE_H_

#include <seqan/basic.h>
#include <seqan/find2.h>

using namespace seqan;

// Test the basic interface for the HammingSimple pattern class.
SEQAN_DEFINE_TEST(test_find2_find_hamming_simple_pattern_interface) {
    // TODO(holtgrew): Strings should be const.
    CharString kHaystack = "he is a hero";
    CharString kNeedle = "her";

    typedef Finder<CharString> TFinder;
    typedef Pattern<CharString, HammingSimple> TPattern;

    // Default constructor.
    {
        TPattern pattern;
    }

    // Construct pattern and finder that are used below.
    TPattern pattern(kNeedle, -1);
    TFinder finder(kHaystack);

    // Function getScoreLimit().
    SEQAN_ASSERT_EQ(-1, getScoreLimit(pattern));

    // Function host().
    CharString const & patternHost = host(pattern);
    SEQAN_ASSERT_EQ(kNeedle, patternHost);
    SEQAN_ASSERT_EQ(3u, length(patternHost));

    // Function needle().
    CharString const & patternNeedle = needle(pattern);
    SEQAN_ASSERT_EQ(kNeedle, patternNeedle);
    SEQAN_ASSERT_EQ(3u, length(patternNeedle));

    // Function find().
    bool ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);

    // Function findBegin().
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);

    // Function score().
    SEQAN_ASSERT_EQ(-1, score(pattern));

    // Test found position.
    //
    // he is a hero
    // her
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(finder));

    // Function alignment().
    Align<CharString, ArrayGaps> align;
    ret = buildAlignment(finder, pattern, align);
    SEQAN_ASSERT_TRUE(ret);
    // Test the alignment by testing the view-to-source conversion.
    // First row.
    SEQAN_ASSERT_EQ(0u, toViewPosition(row(align, 0), 0));
    SEQAN_ASSERT_EQ(1u, toViewPosition(row(align, 0), 1));
    SEQAN_ASSERT_EQ(2u, toViewPosition(row(align, 0), 2));
    SEQAN_ASSERT_EQ(3u, toViewPosition(row(align, 0), 3));
    SEQAN_ASSERT_EQ(4u, toViewPosition(row(align, 0), 4));
    SEQAN_ASSERT_EQ(5u, toViewPosition(row(align, 0), 5));
    SEQAN_ASSERT_EQ(6u, toViewPosition(row(align, 0), 6));
    SEQAN_ASSERT_EQ(7u, toViewPosition(row(align, 0), 7));
    SEQAN_ASSERT_EQ(8u, toViewPosition(row(align, 0), 8));
    SEQAN_ASSERT_EQ(9u, toViewPosition(row(align, 0), 9));
    SEQAN_ASSERT_EQ(10u, toViewPosition(row(align, 0), 10));
    SEQAN_ASSERT_EQ(11u, toViewPosition(row(align, 0), 11));
    // Second row.
    SEQAN_ASSERT_EQ(0u, toViewPosition(row(align, 1), 0));
    SEQAN_ASSERT_EQ(1u, toViewPosition(row(align, 1), 1));
    SEQAN_ASSERT_EQ(2u, toViewPosition(row(align, 1), 2));

    // Function find().
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);

    // Function findBegin().
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);

    // Function score().
    SEQAN_ASSERT_EQ(0, score(pattern));

    // Test found position.
    //
    // he is a hero
    //         her
    SEQAN_ASSERT_EQ(8u, beginPosition(finder));
    SEQAN_ASSERT_EQ(11u, endPosition(finder));

    // Function alignment().
    ret = buildAlignment(finder, pattern, align);
    SEQAN_ASSERT_TRUE(ret);
    // Test the alignment by testing the view-to-source conversion.
    // First row.
    SEQAN_ASSERT_EQ(0u, toViewPosition(row(align, 0), 0));
    SEQAN_ASSERT_EQ(1u, toViewPosition(row(align, 0), 1));
    SEQAN_ASSERT_EQ(2u, toViewPosition(row(align, 0), 2));
    SEQAN_ASSERT_EQ(3u, toViewPosition(row(align, 0), 3));
    SEQAN_ASSERT_EQ(4u, toViewPosition(row(align, 0), 4));
    SEQAN_ASSERT_EQ(5u, toViewPosition(row(align, 0), 5));
    SEQAN_ASSERT_EQ(6u, toViewPosition(row(align, 0), 6));
    SEQAN_ASSERT_EQ(7u, toViewPosition(row(align, 0), 7));
    SEQAN_ASSERT_EQ(8u, toViewPosition(row(align, 0), 8));
    SEQAN_ASSERT_EQ(9u, toViewPosition(row(align, 0), 9));
    SEQAN_ASSERT_EQ(10u, toViewPosition(row(align, 0), 10));
    SEQAN_ASSERT_EQ(11u, toViewPosition(row(align, 0), 11));
    // Second row.
    SEQAN_ASSERT_EQ(8u, toViewPosition(row(align, 1), 0));
    SEQAN_ASSERT_EQ(9u, toViewPosition(row(align, 1), 1));
    SEQAN_ASSERT_EQ(10u, toViewPosition(row(align, 1), 2));

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


// "Easy" (= short, few hits) test of find() using the hamming pattern with a
// score limit of 0.
SEQAN_DEFINE_TEST(test_find2_find_hamming_simple_pattern_find_easy_score_limit_0) {
    typedef Finder<CharString> TFinder;
    typedef Pattern<CharString, HammingSimple> TPattern;

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
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(2u, endPosition(finder));
    SEQAN_ASSERT_EQ(2u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // he is her hero
    //       he
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(8u, endPosition(finder));
    SEQAN_ASSERT_EQ(2u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(6u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // he is her hero
    //           he
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(12u, endPosition(finder));
    SEQAN_ASSERT_EQ(2u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // no more hit
    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// "Easy" (= short, few hits) test of find() using the hamming pattern with a
// score limit of -1.
SEQAN_DEFINE_TEST(test_find2_find_hamming_simple_pattern_find_easy_score_limit_1) {
    typedef Finder<CharString> TFinder;
    typedef Pattern<CharString, HammingSimple> TPattern;

    CharString kHaystack = "he is her hero";
    CharString kNeedle = "her";
    TFinder finder(kHaystack);
    TPattern pattern(kNeedle, -1);

    // We will search for all occurences of the needle using find()
    // and findBegin() in the haystack and check for the correct begin
    // and end position in both finder and needle.  Testing the
    // alignment each time seems overkill.
    bool ret;

    // he is her hero
    // her
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(3u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // he is her hero
    //       her
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(9u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(6u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // he is her hero
    //           her
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(13u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // no more hit
    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// "Harder" (= longer, more hits hits) test of find() using the hamming pattern
// and an score limit of 0.
SEQAN_DEFINE_TEST(test_find2_find_hamming_simple_pattern_find_harder_score_limit_0) {
    typedef Finder<DnaString> TFinder;
    typedef Pattern<DnaString, HammingSimple> TPattern;

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
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // AGAAGAAGAGGAAGAAGA
    //     GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // AGAAGAAGAGGAAGAAGA
    //           GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(13u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // AGAAGAAGAGGAAGAAGA
    //              GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(16u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(13u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // no more hit
    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// "Harder" (= longer, more hits hits) test of find() using the hamming pattern
// and an score limit of -1.
SEQAN_DEFINE_TEST(test_find2_find_hamming_simple_pattern_find_harder_score_limit_1) {
    typedef Finder<DnaString> TFinder;
    typedef Pattern<DnaString, HammingSimple> TPattern;

    DnaString kHaystack = "AGAAGAAGAGGAAGAAGA";
    DnaString kNeedle = "GAA";
    TFinder finder(kHaystack);
    TPattern pattern(kNeedle, -1);

    // We will search for all occurences of the needle using find()
    // and findBegin() in the haystack and check for the correct begin
    // and end position in both finder and needle.  Testing the
    // alignment each time seems overkill.
    bool ret;

    // AGAAGAAGAGGAAGAAGA
    //  GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // AGAAGAAGAGGAAGAAGA
    //     GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // AGAAGAAGAGGAAGAAGA
    //        GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(10u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(7u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // AGAAGAAGAGGAAGAAGA
    //          GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(12u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(9u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // AGAAGAAGAGGAAGAAGA
    //           GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(13u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // AGAAGAAGAGGAAGAAGA
    //              GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(16u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(13u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // no more hit
    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// Tests for setEndPosition() with simple pattern with score limit 0.
SEQAN_DEFINE_TEST(test_find2_find_hamming_simple_pattern_set_end_position_score_limit_0) {
    typedef Finder<DnaString> TFinder;
    typedef Pattern<DnaString, HammingSimple> TPattern;

    DnaString kHaystack = "AGAAGAAGAGGAAGAAGA";
    DnaString kNeedle = "GAA";
    TFinder finder(kHaystack);
    TPattern pattern(kNeedle, 0);

    // We will search for all occurences of the needle using find()
    // and findBegin() in the haystack and check for the correct begin
    // and end position in both finder and needle.  Testing the
    // alignment each time seems overkill.  In between, we will call
    // setEndPosition() sometimes.
    bool ret;

    // Set end position to a hit.
    // AGAAGAAGAGGAAGAAGA
    //  GAA
    ret = setEndPosition(finder, pattern, 4);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // Continue to search from here.
    // AGAAGAAGAGGAAGAAGA
    //     GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

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
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // Set to end.
    // AGAAGAAGAGGAAGAAGA
    //                GAA
    ret = setEndPosition(finder, pattern, 18);
    SEQAN_ASSERT_NOT(ret);

    // No more hit.
    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// Tests for setEndPosition() with simple pattern with score limit -1.
SEQAN_DEFINE_TEST(test_find2_find_hamming_simple_pattern_set_end_position_score_limit_1) {
    typedef Finder<DnaString> TFinder;
    typedef Pattern<DnaString, HammingSimple> TPattern;

    DnaString kHaystack = "AGAAGAAGAGGAAGAAGA";
    DnaString kNeedle = "GAA";
    TFinder finder(kHaystack);
    TPattern pattern(kNeedle, -1);

    // We will search for all occurences of the needle using find()
    // and findBegin() in the haystack and check for the correct begin
    // and end position in both finder and needle.  Testing the
    // alignment each time seems overkill.  In between, we will call
    // setEndPosition() sometimes.
    bool ret;

    // Set end position to a hit.
    // AGAAGAAGAGGAAGAAGA
    //  GAA
    ret = setEndPosition(finder, pattern, 4);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // Continue to search from here.
    // AGAAGAAGAGGAAGAAGA
    //     GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // Set end position to a mismatch.
    // AGAAGAAGAGGAAGAAGA
    //   GAA
    ret = setEndPosition(finder, pattern, 5);
    SEQAN_ASSERT_NOT(ret);

    // Continue to search from here.
    // AGAAGAAGAGGAAGAAGA
    //     GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // Set end position to an approximate hit.
    // AGAAGAAGAGGAAGAAGA
    //          GAA
    ret = setEndPosition(finder, pattern, 12);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(12u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(9u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // Continue to search from here.
    // AGAAGAAGAGGAAGAAGA
    //           GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(13u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // Set to end.
    // AGAAGAAGAGGAAGAAGA
    //                GAA
    ret = setEndPosition(finder, pattern, 18);
    SEQAN_ASSERT_NOT(ret);

    // No more hit.
    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// Tests for setBeginPosition() with simple pattern with score limit -1.
SEQAN_DEFINE_TEST(test_find2_find_hamming_simple_pattern_set_begin_position_score_limit_1) {
    typedef Finder<DnaString> TFinder;
    typedef Pattern<DnaString, HammingSimple> TPattern;

    DnaString kHaystack = "AGAAGAAGAGGAAGAAGA";
    DnaString kNeedle = "GAA";
    TFinder finder(kHaystack);
    TPattern pattern(kNeedle, -1);

    // We will search for all occurences of the needle using find()
    // and findBegin() in the haystack and check for the correct begin
    // and end position in both finder and needle.  Testing the
    // alignment each time seems overkill.  In between, we will call
    // setBeginPosition() sometimes.
    bool ret;

    // Set begin position to a hit.
    // AGAAGAAGAGGAAGAAGA
    //  GAA
    ret = setBeginPosition(finder, pattern, 1);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // Continue to search from here.
    // AGAAGAAGAGGAAGAAGA
    //     GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // Set begin position to a mismatch.
    // AGAAGAAGAGGAAGAAGA
    //   GAA
    ret = setBeginPosition(finder, pattern, 2);
    SEQAN_ASSERT_NOT(ret);

    // Continue to search from here.
    // AGAAGAAGAGGAAGAAGA
    //     GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // Set begin position to an approximate hit.
    // AGAAGAAGAGGAAGAAGA
    //          GAA
    ret = setBeginPosition(finder, pattern, 9);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(-1, score(pattern));
    SEQAN_ASSERT_EQ(12u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(9u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // Continue to search from here.
    // AGAAGAAGAGGAAGAAGA
    //           GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0, score(pattern));
    SEQAN_ASSERT_EQ(13u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // Set to end.
    // AGAAGAAGAGGAAGAAGA
    //                GAA
    ret = setBeginPosition(finder, pattern, 15);
    SEQAN_ASSERT_NOT(ret);

    // No more hit.
    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}

#endif  // TESTS_FIND2_TEST_FIND_HAMMING_SIMPLE_H_
