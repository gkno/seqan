#include <seqan/basic.h>
#include <seqan/find2.h>

using namespace seqan;

// Test the basic interface for the Simple pattern class.
SEQAN_DEFINE_TEST(test_find2_exact_simple_interface) {
    // TODO(holtgrew): Strings should be const.
    CharString kHaystack = "he is a hero";
    CharString kNeedle = "her";

    typedef Finder2<CharString, void> TFinder;
    typedef Pattern2<CharString, Simple> TPattern;

    // Default constructor.
    {
        TPattern pattern;
    }

    // Construct with needle, we will subsequently use this pattern.
    TPattern pattern(kNeedle);
    TFinder finder(kHaystack);

    // Function needle().
    // TODO(holtgrew): Should be needle().
    CharString & patternNeedle = needle2(pattern);
    SEQAN_ASSERT_EQ(kNeedle, patternNeedle);
    SEQAN_ASSERT_EQ(3u, length(patternNeedle));

    // Function find().
    bool ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);

    // Function findBegin().
    ret = findBegin(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);

    // Function endPosition(pattern).
    typedef Position<CharString>::Type TPosition;
    TPosition patternEndPosition = endPosition(pattern);
    SEQAN_ASSERT_EQ(3u, patternEndPosition);

    // Function beginPosition(pattern).
    TPosition patternBeginPosition = beginPosition(pattern);
    SEQAN_ASSERT_EQ(0u, patternBeginPosition);

    // Function endPosition(finder).
    typedef Position<CharString>::Type TPosition;
    TPosition finderEndPosition = endPosition(finder);
    SEQAN_ASSERT_EQ(11u, finderEndPosition);

    // Function beginPosition(finder).
    TPosition finderBeginPosition = beginPosition(finder);
    SEQAN_ASSERT_EQ(8u, finderBeginPosition);
}

SEQAN_BEGIN_TESTSUITE(test_find2) {
    SEQAN_CALL_TEST(test_find2_exact_simple_interface);

    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_approx_dpsearch.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_exact_simple.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_hamming_simple.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/find2/find_multiple_exact_shiftand.h");
}
SEQAN_END_TESTSUITE
