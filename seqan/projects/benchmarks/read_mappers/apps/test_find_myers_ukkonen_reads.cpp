#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "find_myers_ukkonen_reads.h"
#include "find_myers_ukkonen_ext.h"

using namespace seqan;

SEQAN_DEFINE_TEST(test_find_myers_ukkonen_reads_set_end_position) {
    // Setup haystack, needle, finder and pattern.
    CharString hstck = "___AAA___";
    CharString ndl = "AAA";
    Finder<CharString> finder(hstck);
    Pattern<CharString, MyersUkkonenReads> pattern(ndl, -10);

    // First, search for all occurences.
    bool ret;
    // ___AAA___
    // AA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(2u, endPosition(finder));
    SEQAN_ASSERT_EQ(-3, getScore(pattern));
    // ___AAA___
    // AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(3u, endPosition(finder));
    SEQAN_ASSERT_EQ(-3, getScore(pattern));
    // ___AAA___
    //  AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(-2, getScore(pattern));
    // ___AAA___
    //   AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(5u, endPosition(finder));
    SEQAN_ASSERT_EQ(-1, getScore(pattern));
    // ___AAA___
    //    AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(6u, endPosition(finder));
    SEQAN_ASSERT_EQ(0, getScore(pattern));
    // ___AAA___
    //     AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(-1, getScore(pattern));

    // ___AAA___
    // AA
    ret = setEndPosition(finder, pattern, 2u);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(2u, endPosition(finder));
    SEQAN_ASSERT_EQ(-3, getScore(pattern));
    // ___AAA___
    // AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(3u, endPosition(finder));
    SEQAN_ASSERT_EQ(-3, getScore(pattern));

    // ___AAA___
    //    AAA
    ret = setEndPosition(finder, pattern, 6u);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(6u, endPosition(finder));
    SEQAN_ASSERT_EQ(0, getScore(pattern));
    // ___AAA___
    //     AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(-1, getScore(pattern));
}


SEQAN_DEFINE_TEST(test_find_myers_ukkonen_reads_find_begin) {
    // Setup haystack, needle, finder and pattern.
    CharString hstck = "___AAA___";
    CharString ndl = "AAA";
    Finder<CharString> finder(hstck);
    Pattern<CharString, MyersUkkonenReads> pattern(ndl, -10);

    // First, search for all occurences.
    bool ret;
    // ___AAA___
    // AA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(2u, endPosition(finder));
    SEQAN_ASSERT_EQ(-3, getScore(pattern));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));

    // ___AAA___
    //  AA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(3u, endPosition(finder));
    SEQAN_ASSERT_EQ(-3, getScore(pattern));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));

    // ___AAA___
    //   AA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(-2, getScore(pattern));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(2u, beginPosition(finder));

    // ___AAA___
    //    AA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(5u, endPosition(finder));
    SEQAN_ASSERT_EQ(-1, getScore(pattern));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(3u, beginPosition(finder));

    // ___AAA___
    //    AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(6u, endPosition(finder));
    SEQAN_ASSERT_EQ(0, getScore(pattern));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(3u, beginPosition(finder));

    // ___AAA___
    //     AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(-1, getScore(pattern));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
}


SEQAN_DEFINE_TEST(test_find_myers_ukkonen_find_begin) {
    // Setup haystack, needle, finder and pattern.
    CharString hstck = "___AAA___";
    CharString ndl = "AAA";
    Finder<CharString> finder(hstck);
    Pattern<CharString, Myers<FindInfix> > pattern(ndl, -10);

    // First, search for all occurences.
    bool ret;
    // ___AAA___
    // A
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(1u, endPosition(finder));
    SEQAN_ASSERT_EQ(-3, getScore(pattern));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));

    // ___AAA___
    //  A
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(2u, endPosition(finder));
    SEQAN_ASSERT_EQ(-3, getScore(pattern));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));

    // ___AAA___
    //   A
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(3u, endPosition(finder));
    SEQAN_ASSERT_EQ(-3, getScore(pattern));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(2u, beginPosition(finder));

    // ___AAA___
    //    A
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(-2, getScore(pattern));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(3u, beginPosition(finder));

    // ___AAA___
    //    AA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(5u, endPosition(finder));
    SEQAN_ASSERT_EQ(-1, getScore(pattern));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(3u, beginPosition(finder));

    // ___AAA___
    //    AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(6u, endPosition(finder));
    SEQAN_ASSERT_EQ(0, getScore(pattern));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(3u, beginPosition(finder));

    // ___AAA___
    //     AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    SEQAN_ASSERT_EQ(-1, getScore(pattern));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
}


SEQAN_DEFINE_TEST(test_find_myers_ukkonen_set_end_position2) {
    DnaString hstck = "ATTTCGGTCATCAAATAATCATTTATTTTGCCACAACATAAAAATAATTGTCTGAATATGGAATTGTCTGAACCTCACTGAGCTCGTAATAAAATTTCCA";
    Finder<DnaString> finder(hstck);
    DnaString ndl = "CACAACATAAAAAATAATTGTCTGAATATGGAATGT";
    Pattern<DnaString, Myers<FindInfix> > pattern(ndl, -length(ndl));

    bool ret = setEndPosition(finder, pattern, 66u);
    SEQAN_ASSERT_TRUE(ret);
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(31u, beginPosition(finder));
    SEQAN_ASSERT_EQ(66u, endPosition(finder));
//     std::cout << "infix(finder) = " << infix(finder) << std::endl;
//     std::cout << "ndl = " << ndl << std::endl;
    SEQAN_ASSERT_EQ(-3, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT_TRUE(ret);
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_TRUE(ret);
//     std::cout << "infix(finder) = " << infix(finder) << std::endl;
//     std::cout << "ndl = " << ndl << std::endl;
    SEQAN_ASSERT_EQ(31u, beginPosition(finder));
    SEQAN_ASSERT_EQ(67u, endPosition(finder));
    SEQAN_ASSERT_EQ(-2, getScore(pattern));
}


SEQAN_BEGIN_TESTSUITE(test_find_myers_ukkonen_reads) {
    SEQAN_CALL_TEST(test_find_myers_ukkonen_reads_set_end_position);
    SEQAN_CALL_TEST(test_find_myers_ukkonen_reads_find_begin);
    SEQAN_CALL_TEST(test_find_myers_ukkonen_set_end_position2);
    SEQAN_CALL_TEST(test_find_myers_ukkonen_find_begin);
}
SEQAN_END_TESTSUITE
