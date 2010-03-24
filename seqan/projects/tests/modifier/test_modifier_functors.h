#ifndef TESTS_MODIFIER_TEST_MODIFIER_FUNCTORS_H_
#define TESTS_MODIFIER_TEST_MODIFIER_FUNCTORS_H_

#include <seqan/basic.h>
#include <seqan/modifier.h>

using namespace seqan;


SEQAN_DEFINE_TEST(test_modifier_functors_functor_upcase) {
    FunctorUpcase<char> func;

    SEQAN_ASSERT_EQ('A', func('a'));
    SEQAN_ASSERT_EQ('A', func('A'));
    SEQAN_ASSERT_EQ('!', func('!'));
}


SEQAN_DEFINE_TEST(test_modifier_functors_functor_lowcase) {
    FunctorLowcase<char> func;

    SEQAN_ASSERT_EQ('a', func('a'));
    SEQAN_ASSERT_EQ('a', func('A'));
    SEQAN_ASSERT_EQ('!', func('!'));
}


SEQAN_DEFINE_TEST(test_modifier_functors_dna_complement) {
    {
        FunctorComplement<Dna> func;

        SEQAN_ASSERT_EQ(Dna('C'), func(Dna('G')));
        SEQAN_ASSERT_EQ(Dna('G'), func(Dna('C')));
        SEQAN_ASSERT_EQ(Dna('A'), func(Dna('T')));
        SEQAN_ASSERT_EQ(Dna('T'), func(Dna('A')));
    }
    {
        FunctorComplement<Dna5> func;

        SEQAN_ASSERT_EQ(Dna5('C'), func(Dna5('G')));
        SEQAN_ASSERT_EQ(Dna5('G'), func(Dna5('C')));
        SEQAN_ASSERT_EQ(Dna5('A'), func(Dna5('T')));
        SEQAN_ASSERT_EQ(Dna5('T'), func(Dna5('A')));
        SEQAN_ASSERT_EQ(Dna5('N'), func(Dna5('N')));
    }
}

#endif  // TESTS_MODIFIER_TEST_MODIFIER_FUNCTORS_H_
