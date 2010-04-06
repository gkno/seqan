#ifndef SEQAN_TESTS_MISC_TEST_MISC_LONG_WORD_H_
#define SEQAN_TESTS_MISC_TEST_MISC_LONG_WORD_H_

#include <sstream>
#include <string>

#include <seqan/basic.h>
#include <seqan/misc/misc_long_word.h>

using namespace seqan;

// Test the interface of the "Native" specialization.
SEQAN_DEFINE_TEST(test_misc_long_word_native_interface) {
    typedef LongWord<NativeWidth> TLongWord;

    // Default constructor.
    {
        TLongWord word;
        SEQAN_ASSERT_EQ(0u, word);
    }
    // Copy from unsigned.
    {
        TLongWord word(1u);
        SEQAN_ASSERT_EQ(1u, word);
    }
    // Copy from other long word.
    {
        TLongWord word1(1u);
        TLongWord word2(word1);
        SEQAN_ASSERT_EQ(word1, word2);
        SEQAN_ASSERT_EQ(1u, word1);
    }

    // Function length().
    {
        TLongWord word;
        SEQAN_ASSERT_EQ(32u, length(word));
    }

    // Test comparison operators.
    {
        TLongWord word1 = 1u;
        TLongWord word2 = 2u;
        SEQAN_ASSERT_TRUE(word1 != word2);
        SEQAN_ASSERT_NOT(word1 == word2);
        SEQAN_ASSERT_TRUE(word1 < word2);
        SEQAN_ASSERT_TRUE(word1 <= word2);
        SEQAN_ASSERT_NOT(word1 > word2);
        SEQAN_ASSERT_NOT(word1 >= word2);
    }

    // Test shift operators.
    {
        TLongWord word = 1u;
        SEQAN_ASSERT_EQ(0u, word >> 1);
        SEQAN_ASSERT_EQ(1u, word);
        word >>= 1;
        SEQAN_ASSERT_EQ(0u, word);
    }
    {
        TLongWord word = 1u << (sizeof(unsigned) * 8 - 1);
        SEQAN_ASSERT_EQ(0u, word << 1);
        SEQAN_ASSERT_EQ(1u << (sizeof(unsigned) * 8 - 1), word);
        word <<= 1;
        SEQAN_ASSERT_EQ(0u, word);
    }

    // Test bit operators.
    {
        TLongWord word = 0u;
        TLongWord mask = 1u;
        SEQAN_ASSERT_EQ(1u, word | mask);
        word |= mask;
        SEQAN_ASSERT_EQ(1u, word);
    }
    {
        TLongWord word = 0u;
        TLongWord mask = 1u;
        SEQAN_ASSERT_EQ(0u, word & mask);
        word &= mask;
        SEQAN_ASSERT_EQ(0u, word);
        word = 3u;
        SEQAN_ASSERT_EQ(1u, word & mask);
        word &= mask;
        SEQAN_ASSERT_EQ(1u, word);
    }
    {
        TLongWord word = 0u;
        TLongWord mask = 1u;
        SEQAN_ASSERT_EQ(1u, word ^ mask);
        word ^= mask;
        SEQAN_ASSERT_EQ(1u, word);
        word = 3u;
        SEQAN_ASSERT_EQ(2u, word ^ mask);
        word ^= mask;
        SEQAN_ASSERT_EQ(2u, word);
    }

    // Test access to bits.
    {
        // 2147516417d = 10000000000000011000000000000001b.
        TLongWord word = 2147581953u;
        for (unsigned i = 0; i < 32; ++i) {
            if (i == 0u || i == 15u || i == 16u || i == 31u)
                SEQAN_ASSERT_EQ(1u, word[i], "i == %u", i);
            else
                SEQAN_ASSERT_EQ(0u, word[i], "i == %u", i);
        }
        word[0] = 0;
        word[1] = 1;
        for (unsigned i = 0; i < 32; ++i) {
            if (i == 1u || i == 15u || i == 16u || i == 31u)
                SEQAN_ASSERT_EQ(1u, word[i], "i == %u", i);
            else
                SEQAN_ASSERT_EQ(0u, word[i], "i == %u", i);
        }
    }

    // Test const operator[].
    {
        TLongWord const word = 1u;
        SEQAN_ASSERT_EQ(1u, word[0]);
        SEQAN_ASSERT_EQ(0u, word[1]);
    }
}


// Test the interface of the "Static" specialization.  We test the
// type with 129 bits.  This is the first value that will not fit into
// one machine word on today's platforms
SEQAN_DEFINE_TEST(test_misc_long_word_static_interface) {
    typedef LongWord<StaticWidth<129> > TLongWord;

    // Default constructor.
    {
        TLongWord word;
        for (size_t i = 0; i < length(word); ++i)
            SEQAN_ASSERT_EQ(0u, word[i], "i == %u", i);
    }
    // Copy from other long word.
    {
        TLongWord word1;
        word1[0] = 1;
        TLongWord word2(word1);
        SEQAN_ASSERT_EQ(word1, word2);
        SEQAN_ASSERT_EQ(1u, word1[0]);
        for (size_t i = 1; i < length(word2); ++i) {
            SEQAN_ASSERT_EQ(0u, word1[i], "i == %u", i);
        }
    }

    // Test assignment operator.
    {
        TLongWord word1;
        word1[0] = 1u;
        TLongWord word2;
        word2 = word1;
        SEQAN_ASSERT_EQ(1u, word2[0]);
        SEQAN_ASSERT_EQ(0u, word2[1]);
    }

    // Test comparison operators.
    {
        TLongWord word1;
        word1[0] = 1u;
        TLongWord word2;
        word2[1] = 1u;
        SEQAN_ASSERT_TRUE(word1 != word2);
        SEQAN_ASSERT_NOT(word1 == word2);
        SEQAN_ASSERT_TRUE(word1 < word2);
        SEQAN_ASSERT_TRUE(word1 <= word2);
        SEQAN_ASSERT_NOT(word1 > word2);
        SEQAN_ASSERT_NOT(word1 >= word2);
    }

    // Function length().
    {
        TLongWord word;
        SEQAN_ASSERT_EQ(129u, length(word));
    }

    // Test shift operators.
    {
        TLongWord word1;
        word1[0] = 1u;
        TLongWord word2 = word1 >> 1;
        for (size_t i = 0; i < length(word2); ++i)
            SEQAN_ASSERT_EQ(0u, word2[i], "i == %u", i);
        SEQAN_ASSERT_EQ(1u, word1[0]);
        for (size_t i = 1; i < length(word1); ++i)
            SEQAN_ASSERT_EQ(0u, word1[i], "i == %u", i);
        word1 >>= 1;
        for (size_t i = 0; i < length(word1); ++i)
            SEQAN_ASSERT_EQ(0u, word1[i], "i == %u", i);
    }
    {
        TLongWord word1;
        word1[length(word1) - 1] = 1;
        TLongWord word2 = word1 << 1;
        for (size_t i = 0; i < length(word2); ++i)
            SEQAN_ASSERT_EQ(0u, word2[i], "i == %u", i);
        for (size_t i = 0; i < length(word2) - 1; ++i)
            SEQAN_ASSERT_EQ(0u, word2[i], "i == %u", i);
        SEQAN_ASSERT_EQ(1u, word1[length(word1) - 1]);
        word1 <<= 1;
        for (size_t i = 0; i < length(word1); ++i)
            SEQAN_ASSERT_EQ(0u, word1[i], "i == %u", i);
    }

    // Test bit operators.
    {
        TLongWord word;
        TLongWord mask;
        mask[0] = 1;
        TLongWord word2 = word | mask;
        SEQAN_ASSERT_EQ(1u, word2[0]);
        for (size_t i = 1; i < length(word2); ++i)
            SEQAN_ASSERT_EQ(0u, word2[i], "i == %u", i);
        word |= mask;
        SEQAN_ASSERT_EQ(1u, word[0]);
        for (size_t i = 1; i < length(word); ++i)
            SEQAN_ASSERT_EQ(0u, word[i], "i == %u", i);
    }
    {
        TLongWord word;
        TLongWord mask;
        mask[0] = 1;
        TLongWord word2 = word & mask;
        for (size_t i = 0; i < length(word2); ++i)
            SEQAN_ASSERT_EQ(0u, word2[i], "i == %u", i);
        word &= mask;
        for (size_t i = 0; i < length(word); ++i)
            SEQAN_ASSERT_EQ(0u, word[i], "i == %u", i);
        word[0] = 1;
        word[1] = 1;
        TLongWord word3 = word & mask;
        SEQAN_ASSERT_EQ(1u, word3[0]);
        for (size_t i = 1; i < length(word3); ++i)
            SEQAN_ASSERT_EQ(0u, word3[i], "i == %u", i);
        word &= mask;
        SEQAN_ASSERT_EQ(1u, word[0]);
        for (size_t i = 1; i < length(word); ++i)
            SEQAN_ASSERT_EQ(0u, word[i], "i == %u", i);
    }
    {
        TLongWord word;
        TLongWord mask;
        mask[0] = 1;
        TLongWord word2 = word ^ mask;
        SEQAN_ASSERT_EQ(1u, word2[0]);
        for (size_t i = 1; i < length(word2); ++i)
            SEQAN_ASSERT_EQ(0u, word2[i], "i == %u", i);
        word ^= mask;
        SEQAN_ASSERT_EQ(1u, word[0]);
        for (size_t i = 1; i < length(word); ++i)
            SEQAN_ASSERT_EQ(0u, word[i], "i == %u", i);
        word[0] = 1;
        word[1] = 1;
        TLongWord word3 = word ^ mask;
        SEQAN_ASSERT_EQ(0u, word3[0]);
        SEQAN_ASSERT_EQ(1u, word3[1]);
        for (size_t i = 2; i < length(word3); ++i)
            SEQAN_ASSERT_EQ(0u, word3[i], "i == %u", i);
        word ^= mask;
        SEQAN_ASSERT_EQ(0u, word[0]);
        SEQAN_ASSERT_EQ(1u, word[1]);
        for (size_t i = 2; i < length(word); ++i)
            SEQAN_ASSERT_EQ(0u, word[i], "i == %u", i);
    }

    // Test access to bits.
    {
        TLongWord word;
        word[128] = 1;
        word[0] = 1;
        for (unsigned i = 0; i < 129; ++i) {
            if (i == 0u || i == 128u)
                SEQAN_ASSERT_EQ(1u, word[i], "i == %u", i);
            else
                SEQAN_ASSERT_EQ(0u, word[i], "i == %u", i);
        }
        word[0] = 0;
        word[1] = 1;
        for (unsigned i = 0; i < 129; ++i) {
            if (i == 1u || i == 128)
                SEQAN_ASSERT_EQ(1u, word[i], "i == %u", i);
            else
                SEQAN_ASSERT_EQ(0u, word[i], "i == %u", i);
        }
    }

    // Test const operator[].
    {
        TLongWord word;
        word[128] = 1;
        word[0] = 1;
        TLongWord const word2 = word;
        SEQAN_ASSERT_EQ(1u, word2[128]);
        SEQAN_ASSERT_EQ(1u, word2[0]);
    }

    // Test stream shift operator.
    {
        LongWord<StaticWidth<9> > word;
        word[0] = 1;
        word[8] = 1;
        std::stringstream ss(std::stringstream::in | std::stringstream::out);
        ss << word;
        std::string str;
        ss >> str;
        std::string const expectedStr = "100000001";
        SEQAN_ASSERT_EQ(expectedStr, str);
    }
}

#endif  // SEQAN_TESTS_MISC_TEST_MISC_LONG_WORD_H_
