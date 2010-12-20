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
// Test the find2/find_pattern_wild_shiftand.h header.
// ==========================================================================

#ifndef TESTS_FIND2_TEST_FIND_PATTERN_WILD_SHIFTAND_H_
#define TESTS_FIND2_TEST_FIND_PATTERN_WILD_SHIFTAND_H_

#include <seqan/basic.h>
#include <seqan/find2.h>

using namespace seqan;

// Test the _isNumber helper function.
SEQAN_DEFINE_TEST(test_find2_find_pattern_wild_shiftand_is_unsigned) {
    CharString const kUnsignedString = "1234567890";
    CharString const kEmptyString = "";
    CharString const kInvalidString1 = "123456A";
    CharString const kInvalidString2 = "A123456";
    CharString const kInvalidString3 = "_A123456";

    SEQAN_ASSERT_TRUE(_isUnsigned(kUnsignedString));
    SEQAN_ASSERT_NOT(_isUnsigned(kEmptyString));
    SEQAN_ASSERT_NOT(_isUnsigned(kInvalidString1));
    SEQAN_ASSERT_NOT(_isUnsigned(kInvalidString2));
    SEQAN_ASSERT_NOT(_isUnsigned(kInvalidString3));
}


// Tests the _findWildShiftAndIsValid function.
SEQAN_DEFINE_TEST(test_find2_find_pattern_wild_shiftand_is_valid) {
    CharString const kValidString1 = "[A-Z0-9].*";
    CharString const kValidString2 = "x{3,4}";
    CharString const kValidString3 = "3+4?";
    CharString const kValidString4 = ".";

    CharString const kInvalidString1 = "[[A-Z]+";  // Unbalanced parantheses.
    CharString const kInvalidString2 = "[A-Z]]+";  // Unbalanced parantheses 2.
    CharString const kInvalidString3 = "";         // Empty
    CharString const kInvalidString4 = "x{{2,3}";  // Unbalanced parantheses 3.
    CharString const kInvalidString5 = "x{2,3}}";  // Unbalanced parantheses 4.
    CharString const kInvalidString6 = "a\\";      // Wrong usage of escape character.
    CharString const kInvalidString7 = "[A-]+";    // Required chara after dash.

    SEQAN_ASSERT_TRUE(_findWildShiftAndIsValid(kValidString1));
    SEQAN_ASSERT_TRUE(_findWildShiftAndIsValid(kValidString2));
    SEQAN_ASSERT_TRUE(_findWildShiftAndIsValid(kValidString3));
    SEQAN_ASSERT_TRUE(_findWildShiftAndIsValid(kValidString4));

    SEQAN_ASSERT_NOT(_findWildShiftAndIsValid(kInvalidString1));
    SEQAN_ASSERT_NOT(_findWildShiftAndIsValid(kInvalidString2));
    SEQAN_ASSERT_NOT(_findWildShiftAndIsValid(kInvalidString3));
    SEQAN_ASSERT_NOT(_findWildShiftAndIsValid(kInvalidString4));
    SEQAN_ASSERT_NOT(_findWildShiftAndIsValid(kInvalidString5));
    SEQAN_ASSERT_NOT(_findWildShiftAndIsValid(kInvalidString6));
    SEQAN_ASSERT_NOT(_findWildShiftAndIsValid(kInvalidString7));
}


// Tests the _findWildShiftAndLengthWithoutWildcards function.
SEQAN_DEFINE_TEST(test_find2_find_pattern_wild_shiftand_length_without_wildcards) {
    SEQAN_ASSERT_EQ(4u, _findWildShiftAndLengthWithoutWildcards("asdf"));
    SEQAN_ASSERT_EQ(3u, _findWildShiftAndLengthWithoutWildcards("x{1,3}"));
    SEQAN_ASSERT_EQ(3u, _findWildShiftAndLengthWithoutWildcards("x{3,3}"));
    SEQAN_ASSERT_EQ(1u, _findWildShiftAndLengthWithoutWildcards("x*"));
    SEQAN_ASSERT_EQ(2u, _findWildShiftAndLengthWithoutWildcards(".*[A-Z]+"));
    SEQAN_ASSERT_EQ(5u, _findWildShiftAndLengthWithoutWildcards("a?b+c*de"));
}


// Test the find_WildShiftAnd_getCharacterClass function.
SEQAN_DEFINE_TEST(test_find2_find_pattern_wild_shiftand_get_character_class) {
    CharString buffer;

    _findWildShiftAndGetCharacterClass(buffer, "[A-F]", 1, 4);
    SEQAN_ASSERT_EQ("ABCDEF", buffer);

    _findWildShiftAndGetCharacterClass(buffer, "XX[ZA-F]XX", 3, 7);
    SEQAN_ASSERT_EQ("ZABCDEF", buffer);

    _findWildShiftAndGetCharacterClass(buffer, "A-C\\*\\.", 0, 7);
    SEQAN_ASSERT_EQ("ABC*.", buffer);

    _findWildShiftAndGetCharacterClass(buffer, "\\A-\\C\\*\\.", 0, 9);
    SEQAN_ASSERT_EQ("ABC*.", buffer);
}


// Test the basic interface for the Simple pattern class.
SEQAN_DEFINE_TEST(test_find2_find_pattern_wild_shiftand_pattern_interface) {
    // TODO(holtgrew): Strings should be const.
    CharString kHaystack = "he is a hero";
    CharString kNeedle = "her";

    typedef Finder<CharString> TFinder;
    typedef Pattern<CharString, WildShiftAnd> TPattern;

    // Default constructor.
    {
        TPattern pattern;
    }

    // Construt pattern and finder that are used below.
    TPattern pattern(kNeedle);
    TFinder finder(kHaystack);

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
//     ret = findBegin(finder, pattern);
//     SEQAN_ASSERT_TRUE(ret);

    // Test found position.
    SEQAN_ASSERT_EQ(11u, endPosition(finder));
//     SEQAN_ASSERT_EQ(8u, beginPosition(finder));

//     // Function buildAlignment().
//     Align<CharString, ArrayGaps> align;
//     ret = buildAlignment(finder, pattern, align);
//     SEQAN_ASSERT_TRUE(ret);
//     // Test the alignment by testing the view-to-source conversion.
//     // First row.
//     SEQAN_ASSERT_EQ(0u, toViewPosition(row(align, 0), 0));
//     SEQAN_ASSERT_EQ(1u, toViewPosition(row(align, 0), 1));
//     SEQAN_ASSERT_EQ(2u, toViewPosition(row(align, 0), 2));
//     SEQAN_ASSERT_EQ(3u, toViewPosition(row(align, 0), 3));
//     SEQAN_ASSERT_EQ(4u, toViewPosition(row(align, 0), 4));
//     SEQAN_ASSERT_EQ(5u, toViewPosition(row(align, 0), 5));
//     SEQAN_ASSERT_EQ(6u, toViewPosition(row(align, 0), 6));
//     SEQAN_ASSERT_EQ(7u, toViewPosition(row(align, 0), 7));
//     SEQAN_ASSERT_EQ(8u, toViewPosition(row(align, 0), 8));
//     SEQAN_ASSERT_EQ(9u, toViewPosition(row(align, 0), 9));
//     SEQAN_ASSERT_EQ(10u, toViewPosition(row(align, 0), 10));
//     SEQAN_ASSERT_EQ(11u, toViewPosition(row(align, 0), 11));
//     // Second row.
//     SEQAN_ASSERT_EQ(8u, toViewPosition(row(align, 1), 0));
//     SEQAN_ASSERT_EQ(9u, toViewPosition(row(align, 1), 1));
//     SEQAN_ASSERT_EQ(10u, toViewPosition(row(align, 1), 2));

//     // Function length().
//     SEQAN_ASSERT_EQ(3u, length(pattern));

//     // Function infix() for pattern.
//     CharString const kExpectedInfixStr = "her";
//     CharString const kInfixStr = infix(pattern);
//     SEQAN_ASSERT_EQ(kExpectedInfixStr, kInfixStr);

//     // Function begin() for pattern.
//     {
//         typedef Iterator<TPattern::TNeedle const, Standard>::Type TConstIterator;
//         TConstIterator kExpectedBegin = begin(kNeedle, Standard());
//         TConstIterator kBegin = begin(pattern, Standard());
//         SEQAN_ASSERT_EQ(kExpectedBegin, kBegin);
//     }

//     // Function end() for pattern.
//     {
//         typedef Iterator<TPattern::TNeedle const, Standard>::Type TConstIterator;
//         TConstIterator kExpectedEnd = end(kNeedle, Standard());
//         TConstIterator kEnd = end(pattern, Standard());
//         SEQAN_ASSERT_EQ(kExpectedEnd, kEnd);
//     }

//     // Function endPosition(pattern).
//     typedef Position<CharString>::Type TPosition;
//     TPosition patternEndPosition = endPosition(pattern);
//     SEQAN_ASSERT_EQ(3u, patternEndPosition);

//     // Function beginPosition(pattern).
//     TPosition patternBeginPosition = beginPosition(pattern);
//     SEQAN_ASSERT_EQ(0u, patternBeginPosition);
}


// "Easy" (= short, few hits) test of find() using a pattern without
// wildcards.
SEQAN_DEFINE_TEST(test_find2_find_pattern_wild_shiftand_pattern_find_easy) {
    typedef Finder<CharString> TFinder;
    typedef Pattern<CharString, WildShiftAnd> TPattern;

    CharString kHaystack = "he is her hero";
    CharString kNeedle = "he";
    TFinder finder(kHaystack);
    TPattern pattern(kNeedle);

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
//     SEQAN_ASSERT_EQ(0u, beginPosition(finder));
//     SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // he is her hero
    //       he
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(8u, endPosition(finder));
    SEQAN_ASSERT_EQ(2u, endPosition(pattern));
//     SEQAN_ASSERT_EQ(6u, beginPosition(finder));
//     SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // he is her hero
    //           he
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(12u, endPosition(finder));
    SEQAN_ASSERT_EQ(2u, endPosition(pattern));
//     SEQAN_ASSERT_EQ(10u, beginPosition(finder));
//     SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // no more hit
    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// "Harder" (= longer, more hits hits) test of find() using the exact pattern.
SEQAN_DEFINE_TEST(test_find2_find_pattern_wild_shiftand_pattern_find_harder) {
    typedef Finder<DnaString> TFinder;
    typedef Pattern<DnaString, WildShiftAnd> TPattern;

    DnaString kHaystack = "AGAAGAAGAGGAAGAAGA";
    DnaString kNeedle = "GAA";
    TFinder finder(kHaystack);
    TPattern pattern(kNeedle);

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
//     SEQAN_ASSERT_EQ(1u, beginPosition(finder));
//     SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // AGAAGAAGAGGAAGAAGA
    //     GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
//     SEQAN_ASSERT_EQ(4u, beginPosition(finder));
//     SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // AGAAGAAGAGGAAGAAGA
    //           GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(13u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
//     SEQAN_ASSERT_EQ(10u, beginPosition(finder));
//     SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // AGAAGAAGAGGAAGAAGA
    //              GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(16u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
//     SEQAN_ASSERT_EQ(13u, beginPosition(finder));
//     SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // no more hit
    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// Search without any match.
SEQAN_DEFINE_TEST(test_find2_find_pattern_wild_shiftand_pattern_find_nomatch) {
    typedef Finder<DnaString> TFinder;
    typedef Pattern<DnaString, WildShiftAnd> TPattern;

    DnaString kHaystack = "AGAAGAAGAGGAAGAAGA";
    DnaString kNeedle = "CAA";
    TFinder finder(kHaystack);
    TPattern pattern(kNeedle);

    bool ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// Tests for setEndPosition() with simple pattern.
SEQAN_DEFINE_TEST(test_find2_find_pattern_wild_shiftand_pattern_set_end_position) {
    typedef Finder<DnaString> TFinder;
    typedef Pattern<DnaString, WildShiftAnd> TPattern;

    DnaString kHaystack = "AGAAGAAGAGGAAGAAGA";
    DnaString kNeedle = "GAA";
    TFinder finder(kHaystack);
    TPattern pattern(kNeedle);

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
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
//     SEQAN_ASSERT_EQ(1u, beginPosition(finder));
//     SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // Continue to search from here.
    // AGAAGAAGAGGAAGAAGA
    //     GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));
//     SEQAN_ASSERT_EQ(4u, beginPosition(finder));
//     SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

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
//     SEQAN_ASSERT_EQ(4u, beginPosition(finder));
//     SEQAN_ASSERT_EQ(0u, beginPosition(pattern));

    // Set to end.
    // AGAAGAAGAGGAAGAAGA
    //                GAA
    ret = setEndPosition(finder, pattern, 18);
    SEQAN_ASSERT_NOT(ret);

    // No more hit.
    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}

#endif  // TESTS_FIND2_TEST_FIND_PATTERN_WILD_SHIFTAND_H_
