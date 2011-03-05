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
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

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
        bool res = IsSameType<TExpected, TResult>::VALUE;
        SEQAN_ASSERT(res);
    }
    // TODO(holtgrew): Should the modified iterator actually have a value function?
    {
        typedef char TExpected;
        typedef Value<TModifiedIterator>::Type TResult;
        bool res = IsSameType<TExpected, TResult>::VALUE;
        SEQAN_ASSERT(res);
    }
    {
        typedef char TExpected;
        typedef GetValue<TModifiedIterator>::Type TResult;
        bool res = IsSameType<TExpected, TResult>::VALUE;
        SEQAN_ASSERT(res);
    }
    {
        typedef char & TExpected;
        typedef Reference<TModifiedIterator>::Type TResult;
        bool res = IsSameType<TExpected, TResult>::VALUE;
        SEQAN_ASSERT(res);
    }
}


// Test the modifier view iterator.
SEQAN_DEFINE_TEST(test_modifier_view_iterator) {
    typedef CaesarChiffre<char> TFunctor;
    typedef Iterator<CharString, Rooted>::Type TIterator;
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
    {
         TModifiedIterator it(begin(myString, Rooted()));
         assignModViewFunctor(it, myFunctor);
    }
    {
         TFunctor const & kMyFunctor = myFunctor;
         TModifiedIterator it(kMyFunctor);
         it = begin(myString, Rooted());
    }
    {
         TModifiedIterator it(myFunctor);
         it = begin(myString);
    }
    {
         TModifiedIterator it;
         TModifiedIterator it2(it);
         assignModViewFunctor(it, myFunctor);
    }
    {
         const TModifiedIterator kIt;
         TModifiedIterator it(kIt);
         assignModViewFunctor(it, myFunctor);
    }

    // Test value() and getValue().
    {
         TModifiedIterator it(begin(myString, Rooted()));
         assignModViewFunctor(it, myFunctor);
         const TModifiedIterator kIt = it;

        // TODO(holtgrew): The following does not compile.
         SEQAN_ASSERT_EQ('U', value(it));
         SEQAN_ASSERT_EQ('U', value(kIt));

         SEQAN_ASSERT_EQ('U', getValue(it));
         SEQAN_ASSERT_EQ('U', getValue(kIt));
    }
}


SEQAN_DEFINE_TEST(test_modifier_view_const_iterator) {
    SEQAN_ASSERT_FAIL("Implement me!");
}


// Test the modified string class with caesar chiffre.
SEQAN_DEFINE_TEST(test_modifier_view_string_caesar_chiffre) {
    typedef CaesarChiffre<char> TFunctor;
    TFunctor myFunctor(1);

    CharString originalStr = "This is a test!";
    const CharString kExpectedResult = "Uijt jt b uftu!";

    // Test the various ways to initialize a ModifiedString.
    // TODO(holtgrew): Should modified strings not be const to the outside?  Lots of non-const functions are superflous, right?
    {
        ModifiedString<CharString, ModView<TFunctor> > modifiedStr;
        assignModViewFunctor(modifiedStr, myFunctor);
        setHost(modifiedStr, originalStr);

        SEQAN_ASSERT_EQ(kExpectedResult, modifiedStr);
        CharString modifiedStrCopy = modifiedStr;
        SEQAN_ASSERT_EQ(kExpectedResult, modifiedStrCopy);
    }
    {
        ModifiedString<CharString, ModView<TFunctor> > modifiedStr(originalStr);
        assignModViewFunctor(modifiedStr, myFunctor);

        SEQAN_ASSERT_EQ(kExpectedResult, modifiedStr);
        CharString modifiedStrCopy = modifiedStr;
        SEQAN_ASSERT_EQ(kExpectedResult, modifiedStrCopy);
    }
    {
        ModifiedString<CharString, ModView<TFunctor> > modifiedStr(myFunctor);
        setHost(modifiedStr, originalStr);

        SEQAN_ASSERT_EQ(kExpectedResult, modifiedStr);
        CharString modifiedStrCopy = modifiedStr;
        SEQAN_ASSERT_EQ(kExpectedResult, modifiedStrCopy);
    }
    {
        ModifiedString<CharString, ModView<TFunctor> > modifiedStr(originalStr, myFunctor);

        SEQAN_ASSERT_EQ(kExpectedResult, modifiedStr);
        CharString modifiedStrCopy = modifiedStr;
        SEQAN_ASSERT_EQ(kExpectedResult, modifiedStrCopy);
    }
    {
        ModifiedString<CharString, ModView<TFunctor> > modifiedStr(originalStr, myFunctor);

        SEQAN_ASSERT_EQ(kExpectedResult, modifiedStr);
        CharString modifiedStrCopy = modifiedStr;
        SEQAN_ASSERT_EQ(kExpectedResult, modifiedStrCopy);

        ModifiedString<CharString, ModView<TFunctor> > modifiedStr2(modifiedStr);

        SEQAN_ASSERT_EQ(kExpectedResult, modifiedStr2);
        CharString modifiedStrCopy2 = modifiedStr2;
        SEQAN_ASSERT_EQ(kExpectedResult, modifiedStrCopy2);
    }

    // Test operator[].
    {
        ModifiedString<CharString, ModView<TFunctor> > modifiedStr(originalStr, myFunctor);
        SEQAN_ASSERT_EQ('U', modifiedStr[0]);
        SEQAN_ASSERT_EQ('i', modifiedStr[1]);
        SEQAN_ASSERT_EQ('j', modifiedStr[2]);
        SEQAN_ASSERT_EQ('t', modifiedStr[3]);
        SEQAN_ASSERT_EQ(' ', modifiedStr[4]);
    }

    // Test value() and getValue().
    {
        ModifiedString<CharString, ModView<TFunctor> > modifiedStr(originalStr, myFunctor);
        SEQAN_ASSERT_EQ('U', value(modifiedStr, 0));
        SEQAN_ASSERT_EQ('U', getValue(modifiedStr, 0));
    }
}


// Test the modified string class with upper case functor.
SEQAN_DEFINE_TEST(test_modifier_view_string_upper_case) {
    typedef FunctorUpcase<char> TFunctor;
    TFunctor myFunctor;

    CharString originalStr = "This is a test!";
    const CharString kExpectedResult = "THIS IS A TEST!";

    ModifiedString<CharString, ModView<TFunctor> > modifiedStr(originalStr, myFunctor);

    SEQAN_ASSERT_EQ(length(kExpectedResult), length(modifiedStr));
    for (size_t i = 0; i < length(originalStr); ++i)
        SEQAN_ASSERT_EQ_MSG(kExpectedResult[i], modifiedStr[i], "i = %lu", i);

    // TODO(holtgrew): This does not compile.
//     SEQAN_ASSERT_EQ(modifiedStr, kExpectedResult);
    CharString modifiedStrCopy = modifiedStr;
    SEQAN_ASSERT_EQ(modifiedStrCopy, kExpectedResult);

    // We do not test the whole interface as in _caesar_chiffre, so this is
    // more a test of FunctorUpcase.
}


// Test the modified string class with low case functor.
SEQAN_DEFINE_TEST(test_modifier_view_string_low_case) {
    // TODO(holtgrew): Would it make more sense to name this FunctorDowncase?
    typedef FunctorLowcase<char> TFunctor;
    TFunctor myFunctor;

    CharString originalStr = "This is a test!";
    const CharString kExpectedResult = "this is a test!";

    ModifiedString<CharString, ModView<TFunctor> > modifiedStr(originalStr, myFunctor);

    SEQAN_ASSERT_EQ(length(kExpectedResult), length(modifiedStr));
    for (size_t i = 0; i < length(originalStr); ++i)
        SEQAN_ASSERT_EQ_MSG(kExpectedResult[i], modifiedStr[i], "i = %lu", i);

    // TODO(holtgrew): This does not compile.
//     SEQAN_ASSERT_EQ(modifiedStr, kExpectedResult);
    CharString modifiedStrCopy = modifiedStr;
    SEQAN_ASSERT_EQ(modifiedStrCopy, kExpectedResult);

    // We do not test the whole interface as in _caesar_chiffre, so this is
    // more a test of FunctorLowcase
}


// Test the modified string class with alphabet conversion.
SEQAN_DEFINE_TEST(test_modifier_view_string_alphabet_conversion) {
    typedef FunctorConvert<char, Dna5> TFunctor;
    TFunctor myFunctor;

    CharString originalStr = "acgtnACGTN";
    Dna5String const kExpectedResult = "ACGTNACGTN";

    ModifiedString<CharString, ModView<TFunctor> > modifiedStr(originalStr, myFunctor);

    SEQAN_ASSERT_EQ(length(kExpectedResult), length(modifiedStr));
    for (size_t i = 0; i < length(originalStr); ++i)
        SEQAN_ASSERT_EQ_MSG(kExpectedResult[i], modifiedStr[i], "i = %lu", i);

    // TODO(holtgrew): This does not compile.
//     SEQAN_ASSERT_EQ(modifiedStr, kExpectedResult);
    CharString modifiedStrCopy = modifiedStr;
    SEQAN_ASSERT_EQ(modifiedStrCopy, kExpectedResult);

    // We do not test the whole interface as in _caesar_chiffre, so this is
    // more a test of FunctorConvert.
}


// Test the modified string class with nested modifier.
SEQAN_DEFINE_TEST(test_modifier_view_string_nested_modifier) {
    typedef CaesarChiffre<char> TFunctor;
    typedef ModifiedString<ModifiedString<CharString, ModView<TFunctor> >, ModView<TFunctor> > TModifiedString;

    CharString originalStr = "This is a test!";
    CharString const kExpectedResult = "Wklv lv d whvw!";

    TModifiedString modifiedStr(originalStr);
    TFunctor func1(1), func2(2);
    assignModViewFunctor(modifiedStr, func1);
    assignModViewFunctor(host(modifiedStr), func2);

    SEQAN_ASSERT_EQ(modifiedStr, kExpectedResult);

    CharString modifiedStrCopy = modifiedStr;
    SEQAN_ASSERT_EQ(kExpectedResult, modifiedStrCopy);
}


// Test the convert() function.
SEQAN_DEFINE_TEST(test_modifier_convert_in_place) {
    CharString const originalStr = "This is a test!";
    CharString const expectedResult = "Uijt jt b uftu!";

    // Non-const variant on string.
    {
        CharString strCopy = originalStr;
        convert(strCopy, CaesarChiffre<char>(1));
        SEQAN_ASSERT_EQ(expectedResult, strCopy);
    }

    // Non-const variant on segment.
    {
        CharString strCopy = originalStr;
        Segment<CharString, InfixSegment> stringInfix(strCopy);
        convert(strCopy, CaesarChiffre<char>(1));
        SEQAN_ASSERT_EQ(expectedResult, strCopy);
    }

    // Const variant on segment.
    {
        CharString strCopy = originalStr;
        Segment<CharString, InfixSegment> const stringInfix(strCopy);
        convert(strCopy, CaesarChiffre<char>(1));
        SEQAN_ASSERT_EQ(expectedResult, strCopy);
    }
}

#endif  // TESTS_MODIFIER_TEST_MODIFIER_VIEW_H_
