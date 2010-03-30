/*
  Test the find2/find_multiple_exact_simple.h header.
*/
#ifndef TESTS_FIND2_TEST_FIND_EXACT_SIMPLE_H_
#define TESTS_FIND2_TEST_FIND_EXACT_SIMPLE_H_

#include <seqan/basic.h>
#include <seqan/find2.h>

using namespace seqan;

// "Easy" (short) test of the functionality of the MultipleSimple
// pattern.
SEQAN_DEFINE_TEST(test_find2_find_multiple_exact_simple_pattern_find_easy) {
    SEQAN_ASSERT_FAIL("Write me!");
}


// "Harder" (more exhaustive) test of the functionality of the
// MultipleSimple pattern.
SEQAN_DEFINE_TEST(test_find2_find_multiple_exact_simple_pattern_find_harder) {
    SEQAN_ASSERT_FAIL("Write me!");
}


// Test the interface of the MultipleSimple pattern.
SEQAN_DEFINE_TEST(test_find2_find_multiple_exact_simple_pattern_interface) {
    typedef StringSet<CharString> TCharStringSet;
    typedef Finder<CharString> TFinder;
    typedef Pattern<TCharStringSet, MultipleSimple> TPattern;

    // TODO(holtgrew): Strings should be const.
    CharString kHaystack = "he is a hero";
    CharString kNeedle1 = "her";
    CharString kNeedle2 = " is ";
    TCharStringSet kNeedleSet;
    appendValue(kNeedleSet, kNeedle1);
    appendValue(kNeedleSet, kNeedle2);

    // Default constructor.
    {
        TPattern pattern;
    }

    // Construt pattern and finder that are used below.
    TPattern pattern(kNeedleSet);
    TFinder finder(kHaystack);

    // Function needle().
    StringSet<CharString> const & patternNeedles = needles(pattern);
    SEQAN_ASSERT_EQ(2u, length(patternNeedles));
    SEQAN_ASSERT_EQ(kNeedle1, patternNeedles[0]);
    SEQAN_ASSERT_EQ(kNeedle2, patternNeedles[1]);

    // Function find().
    bool ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);

    // Function needleIndex().
    SEQAN_ASSERT_EQ(1u, needleIndex(pattern));

    // Function findBegin().
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);

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

    /*
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
    */
}


// Test the findEndPosition() functionality of the MultipleSimple pattern.
SEQAN_DEFINE_TEST(test_find2_find_multiple_exact_simple_pattern_set_end_position) {
    SEQAN_ASSERT_FAIL("Write me!");
}

#endif  // TESTS_FIND2_TEST_FIND_EXACT_SIMPLE_H_
