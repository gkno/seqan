/*
  Test the find2/find_finder_default.h header.
*/
#ifndef TESTS_FIND2_TEST_FIND_FINDER_DEFAULT_H_
#define TESTS_FIND2_TEST_FIND_FINDER_DEFAULT_H_

#include <seqan/basic.h>
#include <seqan/find2.h>

using namespace seqan;

// Test the basic interface for the Default finder class.
SEQAN_DEFINE_TEST(test_find2_find_finder_default_interface) {
    // TODO(holtgrew): Strings should be const.
    CharString kHaystack = "he is a hero";
    CharString kNeedle = "her";

    typedef Finder<CharString> TFinder;
    typedef Pattern<CharString, Simple> TPattern;

    // Construt pattern and finder that are used below.
    TPattern pattern(kNeedle);
    TFinder finder(kHaystack);

    // Test needle() and host().
    CharString kExpectedNeedle = "her";
    CharString hst = host(pattern);
    SEQAN_ASSERT_EQ(kExpectedNeedle, hst);
    CharString ndl = needle(pattern);
    SEQAN_ASSERT_EQ(kExpectedNeedle, ndl);

    // Function find().
    bool ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);

    // Function findBegin().
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);

    // Function begin() for haystack.
    {
        typedef Iterator<TFinder::THaystack const, Standard>::Type TConstIterator;
        TConstIterator kExpectedBegin = begin(kHaystack, Standard()) + 8;
        TConstIterator kBegin = begin(finder, Standard());
        SEQAN_ASSERT_EQ(kExpectedBegin, kBegin);
    }

    // Function end() for haystack.
    {
        typedef Iterator<TPattern::TNeedle const, Standard>::Type TConstIterator;
        TConstIterator kExpectedEnd = begin(kHaystack) + 11;
        TConstIterator kEnd = end(finder, Standard());
        SEQAN_ASSERT_EQ(kExpectedEnd, kEnd);
    }

    // Function endPosition(finder).
    typedef Position<CharString>::Type TPosition;
    TPosition finderEndPosition = endPosition(finder);
    SEQAN_ASSERT_EQ(11u, finderEndPosition);

    // Function beginPosition(finder).
    TPosition finderBeginPosition = beginPosition(finder);
    SEQAN_ASSERT_EQ(8u, finderBeginPosition);

    // Function infix() for finder.
    CharString const kExpectedInfixStr = "her";
    CharString const kInfixStr = infix(finder);
    SEQAN_ASSERT_EQ(kExpectedInfixStr, kInfixStr);
}

#endif  // TESTS_FIND2_TEST_FIND_FINDER_DEFAULT_H_
