/*
  Test the find2/find_pattern_wild_shiftand.h header.
*/
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


// Tests the _find_WildShiftAnd_isValid function.
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

    SEQAN_ASSERT_TRUE(_find_WildShiftAnd_isValid(kValidString1));
    SEQAN_ASSERT_TRUE(_find_WildShiftAnd_isValid(kValidString2));
    SEQAN_ASSERT_TRUE(_find_WildShiftAnd_isValid(kValidString3));
    SEQAN_ASSERT_TRUE(_find_WildShiftAnd_isValid(kValidString4));

    SEQAN_ASSERT_NOT(_find_WildShiftAnd_isValid(kInvalidString1));
    SEQAN_ASSERT_NOT(_find_WildShiftAnd_isValid(kInvalidString2));
    SEQAN_ASSERT_NOT(_find_WildShiftAnd_isValid(kInvalidString3));
    SEQAN_ASSERT_NOT(_find_WildShiftAnd_isValid(kInvalidString4));
    SEQAN_ASSERT_NOT(_find_WildShiftAnd_isValid(kInvalidString5));
    SEQAN_ASSERT_NOT(_find_WildShiftAnd_isValid(kInvalidString6));
    SEQAN_ASSERT_NOT(_find_WildShiftAnd_isValid(kInvalidString7));
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
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);

    // Test found position.
    SEQAN_ASSERT_EQ(8u, beginPosition(finder));
    SEQAN_ASSERT_EQ(11u, endPosition(finder));

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


// "Easy" (= short, few hits) test of find() using the exact pattern.
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
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));
    SEQAN_ASSERT_EQ(2u, endPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    SEQAN_ASSERT_EQ(2u, endPosition(pattern));

    // he is her hero
    //       he
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(6u, beginPosition(finder));
    SEQAN_ASSERT_EQ(8u, endPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    SEQAN_ASSERT_EQ(2u, endPosition(pattern));

    // he is her hero
    //           he
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(12u, endPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    SEQAN_ASSERT_EQ(2u, endPosition(pattern));

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
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));

    // AGAAGAAGAGGAAGAAGA
    //     GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));

    // AGAAGAAGAGGAAGAAGA
    //           GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(10u, beginPosition(finder));
    SEQAN_ASSERT_EQ(13u, endPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));

    // AGAAGAAGAGGAAGAAGA
    //              GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(13u, beginPosition(finder));
    SEQAN_ASSERT_EQ(16u, endPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));

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
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));

    // Continue to search from here.
    // AGAAGAAGAGGAAGAAGA
    //     GAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));

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
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(0u, beginPosition(pattern));
    SEQAN_ASSERT_EQ(3u, endPosition(pattern));

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
