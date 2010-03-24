#ifndef TESTS_MODIFIER_TEST_MODIFIER_VIEW_H_
#define TESTS_MODIFIER_TEST_MODIFIER_VIEW_H_

#include <seqan/basic.h>
#include <seqan/modifier.h>

#include "helpers.h"

using namespace seqan;

// Test the modifier view metafunctions.
SEQAN_DEFINE_TEST(test_modifier_view_iterator_metafunctions) {
    typedef CaesarChiffre<char> TFunctor;
    typedef ModifiedIterator<CharString, ModView<TFunctor> > TModifiedIterator;

    // TODO(holtgrew): We really want a static assertion here.
    {
        typedef ModViewCargo<TFunctor> TExpected;
        typedef Cargo<TModifiedIterator>::Type TResult;
        bool res = TYPECMP<TExpected, TResult>::VALUE;
        SEQAN_ASSERT_TRUE(res);
    }
    // TODO(holtgrew): Should the modified iterator actually have a value function?
    {
        typedef char TExpected;
        typedef Value<TModifiedIterator>::Type TResult;
        bool res = TYPECMP<TExpected, TResult>::VALUE;
        SEQAN_ASSERT_TRUE(res);
    }
    {
        typedef char TExpected;
        typedef GetValue<TModifiedIterator>::Type TResult;
        bool res = TYPECMP<TExpected, TResult>::VALUE;
        SEQAN_ASSERT_TRUE(res);
    }
    {
        typedef char TExpected;
        typedef Reference<TModifiedIterator>::Type TResult;
        bool res = TYPECMP<TExpected, TResult>::VALUE;
        SEQAN_ASSERT_TRUE(res);
    }
}


// Test the modifier view iterator.
SEQAN_DEFINE_TEST(test_modifier_view_iterator) {
    typedef CaesarChiffre<char> TFunctor;
    typedef Iterator<CharString, Rooted> TIterator;
    typedef ModifiedIterator<TIterator, ModView<TFunctor> > TModifiedIterator;

    // The string and functor we will work with.
    CharString myString = "This is a nice string!";
    TFunctor myFunctor(1);

    // Manually shift characters in string by one.
    CharString rotatedString = "This is a nice string!";
    for (size_t i = 0; i < length(rotatedString); ++i) {
        if (rotatedString[i] == ' ' || rotatedString[i] == '!')
            continue;
        rotatedString[i] += 1;
    }

    // Test the various ways to construct the iterator with both a
    // container and a functor.
    {
        TModifiedIterator it;
        assignModViewFunctor(it, myFunctor);
    }
    SEQAN_ASSERT_FAIL("The following does not compile!");
    {
//         TModifiedIterator it(begin(myString, Rooted()));
//         assignModViewFunctor(it, myFunctor);
    }
    {
//         TFunctor const & kMyFunctor = myFunctor;
//         TModifiedIterator it(kMyFunctor);
//         it = begin(myString, Rooted());
    }
    {
//         TModifiedIterator it(myFunctor);
//         it = begin(myString);
    }
    {
//         TModifiedIterator it;
//         TModifiedIterator it2(it);
//         assignModViewFunctor(it, myFunctor);
    }
    {
//         const TModifiedIterator kIt;
//         TModifiedIterator it(kIt);
//         assignModViewFunctor(it, myFunctor);
    }

    // Test value() and getValue().
    {
//         TModifiedIterator it(begin(myString, Rooted()));
//         assignModViewFunctor(it, myFunctor);
//         const TModifiedIterator kIt = it;

        // TODO(holtgrew): The following does not compile.
//         SEQAN_ASSERT_EQ('U', value(it));
//         SEQAN_ASSERT_EQ('U', value(kIt));

//         SEQAN_ASSERT_EQ('U', getValue(it));
//         SEQAN_ASSERT_EQ('U', getValue(kIt));
    }
}


SEQAN_DEFINE_TEST(test_modifier_view_const_iterator) {
    SEQAN_ASSERT_FAIL("Implement me!");
}

// Test the convertInPlace() function.
SEQAN_DEFINE_TEST(test_modifier_convert_in_place) {
    SEQAN_ASSERT_FAIL("Implement me!");
}

#endif  // TESTS_MODIFIER_TEST_MODIFIER_VIEW_H_
