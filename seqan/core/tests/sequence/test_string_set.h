// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================
// This file contains functions to test the functionality of the sequence
// module.
// ==========================================================================

#ifndef CORE_TESTS_SEQUENCE_TEST_STRINGSET_H_
#define CORE_TESTS_SEQUENCE_TEST_STRINGSET_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include "test_sequence.h"

// --------------------------------------------------------------------------
// Generic String Tests
// --------------------------------------------------------------------------

template <typename TReturnStringSet, typename TStringSet>
TStringSet createStringSet(TStringSet & stringSet)
{
    return TReturnStringSet(stringSet);
}

// Test whether PrefixSegments are copy constructible.
template <typename TStringSet>
void testStringSetCopyConstructible(TStringSet & /*Tag*/)
{
    using namespace seqan;

    TStringSet stringSet1;
    resize(stringSet1, 3);
    stringSet1[0] = "AAAA";
    stringSet1[1] = "CC";
    stringSet1[2] = "GGG";

    {
        TStringSet stringSet2(stringSet1);

        SEQAN_ASSERT_EQ(getValue(stringSet2, 0), "AAAA");
        SEQAN_ASSERT_EQ(getValue(stringSet2, 1), "CC");
        SEQAN_ASSERT_EQ(getValue(stringSet2, 2), "GGG");
    }
    {
        TStringSet const stringSet2(stringSet1);

        SEQAN_ASSERT_EQ(getValue(stringSet2, 0), "AAAA");
        SEQAN_ASSERT_EQ(getValue(stringSet2, 1), "CC");
        SEQAN_ASSERT_EQ(getValue(stringSet2, 2), "GGG");
    }
}

// Test whether sequences are default constructible.
template <typename TStringSet>
void testStringSetDefaultConstructible(TStringSet & /*Tag*/)
{
    using namespace seqan;

    TStringSet stringSet;
    SEQAN_ASSERT_EQ(begin(stringSet), end(stringSet));
}

// Test operator<().
template <typename TStringSet>
void testStringSetLessGreaterEqual(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;

    TNonConstStringSet nonConstStringSet1;
    TNonConstStringSet nonConstStringSet2;

    // Nothing is equal to nothing.
    {
        TStringSet stringSet1, stringSet2;
        SEQAN_ASSERT_EQ(stringSet1 < stringSet2, false);
        SEQAN_ASSERT_EQ(stringSet1 <= stringSet2, true);
        SEQAN_ASSERT_EQ(stringSet1 > stringSet2, false);
        SEQAN_ASSERT_EQ(stringSet1 >= stringSet2, true);
        SEQAN_ASSERT_EQ(stringSet1 == stringSet2, true);
        SEQAN_ASSERT_EQ(stringSet1 != stringSet2, false);
    }

    // Something (even uninitialized) is larger than nothing.
    {
        resize(nonConstStringSet1, 3);
        TStringSet stringSet1(nonConstStringSet1);
        TStringSet stringSet2(nonConstStringSet2);
        SEQAN_ASSERT_EQ(stringSet1 < stringSet2, false);
        SEQAN_ASSERT_EQ(stringSet1 <= stringSet2, false);
        SEQAN_ASSERT_EQ(stringSet1 > stringSet2, true);
        SEQAN_ASSERT_EQ(stringSet1 >= stringSet2, true);
        SEQAN_ASSERT_EQ(stringSet1 == stringSet2, false);
        SEQAN_ASSERT_EQ(stringSet1 != stringSet2, true);
    }

    // Something (initialized) is larger than nothing.
    {
        resize(nonConstStringSet1, 3);
        nonConstStringSet1[0] = "AAAA";
        nonConstStringSet1[1] = "CCCC";
        nonConstStringSet1[2] = "GGGG";
        TStringSet stringSet1(nonConstStringSet1);
        TStringSet stringSet2(nonConstStringSet2);
        SEQAN_ASSERT_EQ(stringSet1 < stringSet2, false);
        SEQAN_ASSERT_EQ(stringSet1 <= stringSet2, false);
        SEQAN_ASSERT_EQ(stringSet1 > stringSet2, true);
        SEQAN_ASSERT_EQ(stringSet1 >= stringSet2, true);
        SEQAN_ASSERT_EQ(stringSet1 == stringSet2, false);
        SEQAN_ASSERT_EQ(stringSet1 != stringSet2, true);
    }

    // Something (unintialized) is smaller than something (initialized).
    {
        resize(nonConstStringSet2, 3);
        TStringSet stringSet1(nonConstStringSet1);
        TStringSet stringSet2(nonConstStringSet2);
        SEQAN_ASSERT_EQ(stringSet1 < stringSet2, false);
        SEQAN_ASSERT_EQ(stringSet1 <= stringSet2, false);
        SEQAN_ASSERT_EQ(stringSet1 > stringSet2, true);
        SEQAN_ASSERT_EQ(stringSet1 >= stringSet2, true);
        SEQAN_ASSERT_EQ(stringSet1 == stringSet2, false);
        SEQAN_ASSERT_EQ(stringSet1 != stringSet2, true);
    }
    {/*...*/}
}

// Test of append().
template <typename TStringSet>
void testStringSetAppend(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TStringSet>::Type TString;

    TStringSet stringSet1;
    TStringSet stringSet2;

    // Test the append function
    appendValue(stringSet1, stringSet2);
    SEQAN_ASSERT_EQ(length(stringSet1), 0u);

    // Test the appendValue function
    resize(stringSet2, 3);
    append(stringSet1, stringSet2);
    SEQAN_ASSERT_EQ(length(stringSet1), 3u);
}

// Test of appendValue().
template <typename TStringSet>
void testStringSetAppendValue(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TStringSet>::Type TString;

    TStringSet stringSet;
    TStringSet stringSet2;

    // Test the appendValue function
    TString string = "ACGT";
    appendValue(stringSet, string);
    SEQAN_ASSERT_EQ(length(stringSet), 1u);
    SEQAN_ASSERT_EQ(stringSet[0], "ACGT");

    // Test the appendValue function
    appendValue(stringSet, "CGTA");
    SEQAN_ASSERT_EQ(length(stringSet), 2u);
    SEQAN_ASSERT_EQ(stringSet[0], "ACGT");
    SEQAN_ASSERT_EQ(stringSet[1], "CGTA");
}

// Test of assign().
template <typename TStringSet>
void testStringSetAssign(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TStringSet>::Type TString;

    TStringSet stringSet1;
    TStringSet stringSet2;

    // Test the append function
    assign(stringSet1, stringSet2);
    SEQAN_ASSERT_EQ(length(stringSet1), 0u);
    SEQAN_ASSERT(begin(stringSet1) == end(stringSet1));

    // Test the appendValue function
    resize(stringSet2, 3);
    assign(stringSet1, stringSet2);
    for(unsigned i = 0; i < length(stringSet1); ++i)
        SEQAN_ASSERT(stringSet1[1] == stringSet2[i]);
}

// Test of assignValue().
template <typename TStringSet>
void testStringSetAssignValue(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TStringSet>::Type TString;

    TStringSet stringSet;
    TString string;

    // Test the assignValue function
    resize(stringSet, 3);
    assignValue(stringSet, 1u, string);
    SEQAN_ASSERT_EQ(length(stringSet), 3u);
    SEQAN_ASSERT_EQ(stringSet[1], string);
}

// Test of assignValueById().
template <typename TStringSet>
void testStringSetAssignValueById(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TStringSet>::Type TString;

    TStringSet stringSet1;

    // Assigning a string.
    unsigned id = assignValueById(stringSet1, "ACGT");
    SEQAN_ASSERT_EQ(length(stringSet1), 1u);
    SEQAN_ASSERT_EQ(stringSet1[0], "ACGT");
    SEQAN_ASSERT_EQ(id, 0u);

    // TODO (singer): this is not documented, such that it can be easily understood!
    // Assigning a string from a string set.
    TStringSet stringSet2;
    resize(stringSet2, 3);
    stringSet2[0] = "AAAA";
    stringSet2[1] = "TTTT";
    id = assignValueById(stringSet1, stringSet2, (unsigned)1);

    SEQAN_ASSERT_EQ(length(stringSet1), 2u);
    SEQAN_ASSERT_EQ(stringSet1[0], "ACGT");
    SEQAN_ASSERT_EQ(stringSet1[1], "TTTT");
    SEQAN_ASSERT_EQ(id, 1u);
}

// Test of back() for non const strings.
// TODO (singer): not in docu.
template <typename TStringSet>
void testStringSetBack(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TStringSet>::Type TString;

    TStringSet stringSet;
    resize(stringSet, 3);
    stringSet[2] = "ACGT";

    // val is a reference in contrast to the const version of front()
    TString & val = back(stringSet);
    val = "TTTT";
    SEQAN_ASSERT_EQ(val, stringSet[2]);
}

// Test of back() for const strings.
// TODO (singer): not in docu.
template <typename TStringSet>
void testStringSetBack(TStringSet const & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TStringSet>::Type TString;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;

    TNonConstStringSet nonConstStringSet;
    resize(nonConstStringSet, 3);
    nonConstStringSet[2] = "ACGT";
    TStringSet stringSet(nonConstStringSet);

    // val is a reference in contrast to the const version of front()
    TString val = back(stringSet);
    SEQAN_ASSERT_EQ(val, stringSet[2]);
    val = "TTTT";
    SEQAN_ASSERT_EQ(stringSet[2], "ACGT");
}

// Test of begin().
template <typename TStringSet>
void testStringSetBegin(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TStringSet>::Type TString;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;

    TNonConstStringSet nonConstStringSet;
    resize(nonConstStringSet, 3);
    nonConstStringSet[0] = "ACGT";
    TStringSet stringSet(nonConstStringSet);
    SEQAN_ASSERT_EQ(*begin(stringSet), "ACGT");
    SEQAN_ASSERT_EQ(*begin(stringSet, Standard()), "ACGT");
    SEQAN_ASSERT_EQ(*begin(stringSet, Rooted()), "ACGT");
}

// Test of beginPosition().
template <typename TStringSet>
void testStringSetBeginPosition(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;

    // Test on an empty string.
    TNonConstStringSet nonConstStringSet;
    resize(nonConstStringSet, 3);
    {
        TStringSet stringSet(nonConstStringSet);
        SEQAN_ASSERT_EQ(beginPosition(stringSet), 0u);
    }

    // Test on a non empty string.
    {
        nonConstStringSet[0] = "ACGT";
        TStringSet stringSet(nonConstStringSet);
        SEQAN_ASSERT_EQ(beginPosition(stringSet), 0u);
    }
}

// Test of clear().
template <typename TStringSet>
void testStringSetClear(TStringSet & /*Tag*/)
{
    using namespace seqan;

    // Test on an empty string set.
    {
        TStringSet stringSet;
        clear(stringSet);
        SEQAN_ASSERT(begin(stringSet) == end(stringSet));
    }

    // Test on a non empty string set.
    {
        TStringSet stringSet;
        resize(stringSet,3);
        stringSet[0] = "ACGT";
        clear(stringSet);
        SEQAN_ASSERT(begin(stringSet) == end(stringSet));
    }
}

// Test of concat().
template <typename TStringSet>
void testStringSetConcat(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TStringSet>::Type TString;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;

    // TODO (singer): test fails in non init string sets.
    // Error message: Assertion failed : static_cast<TStringPos>(pos) < static_cast<TStringPos>(length(me)) was: 0 >= 0 (Trying to access an element behind the last one!).
    // Test on an empty string set.
    {
        TStringSet stringSet;
        ConcatenatorManyToOne<TStringSet> concatString = concat(stringSet);
        //SEQAN_ASSERT(begin(concatString) == end(concatString));
    }

    // Test on a non empty string set.
    {
    	TNonConstStringSet nonConstStringSet;
    	resize(nonConstStringSet, 4);
    	nonConstStringSet[0] = "AAAA";
    	nonConstStringSet[1] = "CCCC";
    	nonConstStringSet[2] = "GGGG";
    	nonConstStringSet[3] = "TTTT";
    	TString string = "AAAACCCCGGGGTTTT";
        TStringSet stringSet(nonConstStringSet);
        ConcatenatorManyToOne<TStringSet> concatString = concat(stringSet);
        for (unsigned i = 0; i < length(string); ++i)
        	SEQAN_ASSERT_EQ(string[i], concatString[i]);
    }
}

// Test of end().
template <typename TStringSet>
void testStringSetEnd(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;

    TNonConstStringSet nonConstStringSet;
    resize(nonConstStringSet, 3);
    nonConstStringSet[2] = "ACGT";
    TStringSet stringSet(nonConstStringSet);

     typename Iterator<TStringSet>::Type iter = end(stringSet);
     typename Iterator<TStringSet, Standard>::Type standardIter = end(stringSet, Standard());
     typename Iterator<TStringSet, Rooted>::Type rootedIter = end(stringSet, Rooted());

    --iter;
    --standardIter;
    --rootedIter;

    SEQAN_ASSERT_EQ(*iter, "ACGT");
    SEQAN_ASSERT_EQ(*standardIter, "ACGT");
    SEQAN_ASSERT_EQ(*rootedIter, "ACGT");
}

// Test of endPosition().
// TODO (singer): erase() is not in the docu.
template <typename TStringSet>
void testStringSetEndPosition(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;

    // Test on an empty string.
    TNonConstStringSet nonConstStringSet;
    resize(nonConstStringSet, 3);
    {
        TStringSet stringSet(nonConstStringSet);
        SEQAN_ASSERT_EQ(endPosition(stringSet), 3u);
    }

    // Test on a non empty string.
    {
        nonConstStringSet[2] = "ACGT";
        TStringSet stringSet(nonConstStringSet);
        SEQAN_ASSERT_EQ(endPosition(stringSet), 3u);
    }
}

// Test of erase().
// TODO (singer): erase() is not in the docu.
template <typename TStringSet>
void testStringSetErase(TStringSet & /*Tag*/)
{
    using namespace seqan;

    // Test on an empty string.

    typedef typename Value<TStringSet>::Type TString;

    TStringSet stringSet;
    resize(stringSet, 3);
    stringSet[1] = "ACGTACGTACGT";
    SEQAN_ASSERT_EQ(stringSet[1], "ACGTACGTACGT");
    erase(stringSet, 1);
    SEQAN_ASSERT_EQ(length(stringSet), 2u);
    SEQAN_ASSERT_EQ(stringSet[0], TString());
    SEQAN_ASSERT_EQ(stringSet[1], TString());
}

// Test of eraseBack().
// TODO (singer): eraseBack() is not in the docu.
template <typename TStringSet>
void testStringSetEraseBack(TStringSet & /*Tag*/)
{
    using namespace seqan;


    typedef typename Value<TStringSet>::Type TString;

    TStringSet stringSet;
    // Test on an empty string.
    // TODO (singer): does not work on empty string
    // Assertion failed : length(me) > 0u was: 0 <= 0 (String must have more than 0 characters in eraseBack()!)
    //eraseBack(stringSet);

    // Test on a non empty string.
    resize(stringSet, 3);
    stringSet[2] = "ACGTACGTACGT";
    SEQAN_ASSERT_EQ(stringSet[2], "ACGTACGTACGT");
    eraseBack(stringSet);
    SEQAN_ASSERT_EQ(length(stringSet), 2u);
    SEQAN_ASSERT_EQ(stringSet[0], TString());
    SEQAN_ASSERT_EQ(stringSet[1], TString());
}

// TODO (singer): not in docu.
// Test of front().
template <typename TStringSet>
void testStringSetFront(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TStringSet>::Type TString;

    TStringSet stringSet;
    resize(stringSet, 3);
    stringSet[0] = "ACGT";

    // val is a reference in contrast to the const version of front()
    TString & val = front(stringSet);
    val = "TTTT";
    SEQAN_ASSERT_EQ(val, stringSet[0]);
}

// Test of front() for const strings.
// TODO (singer): not in docu.
template <typename TStringSet>
void testStringSetFront(TStringSet const & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TStringSet>::Type TString;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;

    TNonConstStringSet nonConstStringSet;
    resize(nonConstStringSet, 3);
    nonConstStringSet[0] = "ACGT";
    TStringSet stringSet(nonConstStringSet);

    // val is a reference in contrast to the const version of front()
    TString val = front(stringSet);
    SEQAN_ASSERT_EQ(val, stringSet[0]);
    val = "TTTT";
    SEQAN_ASSERT_EQ(stringSet[0], "ACGT");
}

// Test of getValue().
template <typename TStringSet>
void testStringSetGetValue(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;

     // In contrast to value(), getValue() does not return a reference but a copy.
     // We test this using the variable value_.
    TNonConstStringSet nonConstStringSet;
    resize(nonConstStringSet, 3);
    nonConstStringSet[1]= "ACGT";
    TStringSet stringSet(nonConstStringSet);
    SEQAN_ASSERT_EQ(getValue(stringSet, 1), "ACGT");
}

// TODO (singer): not defined for const string sets.
// Test of getValueById().
template <typename TStringSet>
void testStringSetGetValueById(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TStringSet>::Type TString;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;

    TNonConstStringSet nonConstStringSet;
    resize(nonConstStringSet, 3);
    nonConstStringSet[1]= "ACGT";
    TStringSet stringSet(nonConstStringSet);
    SEQAN_ASSERT_EQ(getValueById(stringSet, 1u), "ACGT");
}

// TODO (singer): define behaviour and adjust test.
// Infix() compiles and does what it is supposed to do?!
// However, it is not very intuitive. For details see comments below.
// There is a need to improve the documentation of this!
// Test of infix()
template <typename TStringSet>
void testStringSetInfix(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TStringSet>::Type TString;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;

    TNonConstStringSet nonConstStringSet;
    resize(nonConstStringSet, 3);
    nonConstStringSet[1]= "ACGT";

    // Only non-const test for this scenario possible.
    TStringSet stringSet(nonConstStringSet);
    TString string = infix(stringSet, 0, 1); // Returns the first character not string!
    // std::cerr << string << std::endl; -> "A"

    // Only non-const test for this scenario possible.
    nonConstStringSet[0] = "TT";

    // Since the infix (should) point to the fist element it should point to "T"
    // std::cerr << string << std::endl; -> "A"
    // Therefore there is a different behaviour to normal strings.
}

// TODO (singer): define behaviour and adjust test.
// Infix() compiles and does what it is supposed to do?!
// However, it is not very intuitive. For details see comments below.
// There is a need to improve the documentation of this!
// Test of infixWithLength()
template <typename TStringSet>
void testStringSetInfixWithLength(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TStringSet>::Type TString;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;

    TNonConstStringSet nonConstStringSet;
    resize(nonConstStringSet, 3);
    nonConstStringSet[1]= "ACGT";

    // Test of const and non-const version
    TStringSet stringSet(nonConstStringSet);
    TString string = infix(stringSet, 0, 1); // Returns the first character not string!
    // std::cerr << string << std::endl; // -> "A"

    // Only non-const test for this scenario possible.
    nonConstStringSet[0] = "TT";

    // Since the infix (should) point to the fist element it should point to "T"
    // std::cerr << string << std::endl; // -> "A"
    // Therefore there is a different behaviour to normal strings.
}

// Test of insert().
// TODO (singer): no insert function implemented.
template <typename TStringSet>
void testStringSetInsert(TStringSet & /*Tag*/)
{
    //using namespace seqan;

    //// Test of inserting an empty string.
    //TStringSet stringSet1;
    //resize(stringSet1, 1);
    //TStringSet stringSet2;
    //insert(stringSet1, 0u, stringSet2);
    //SEQAN_ASSERT_EQ(length(stringSet1), 1u);

    //resize(stringSet2, 3);
    //stringSet2[0] = "ACGT";
    //insert(stringSet1, 0u, stringSet2);
    //SEQAN_ASSERT_EQ(length(stringSet1), 4u);
    //SEQAN_ASSERT_EQ(stringSet1[1], "ACGT");
}

// Test of insertValue().
// TODO (singer): no insertValue function implemented.
template <typename TStringSet>
void testStringSetInsertValue(TStringSet & /*Tag*/)
{
    //using namespace seqan;
    //
    //typedef typename Value<TStringSet>::Type TString;

    //// Test of inserting an empty string.
    //TStringSet stringSet;
    //resize(stringSet, 1);
    //insertValue(stringSet, 0, "ACGT");
    //SEQAN_ASSERT_EQ(length(stringSet), 2u);
    //SEQAN_ASSERT_EQ(stringSet[0], "ACGT");
    //SEQAN_ASSERT_EQ(stringSet[1], TString());
}

// Test of iter().
template <typename TStringSet>
void testStringSetIter(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;
    typedef typename Iterator<TStringSet>::Type TIterator;
    typedef typename Iterator<TStringSet, Standard>::Type TStandardIterator;
    typedef typename Iterator<TStringSet, Rooted>::Type TRootedIterator;

    // Test on an empty string set.
    {
    	TStringSet stringSet;
        TIterator iterator = iter(stringSet, 0);
        TStandardIterator standardIterator = iter(stringSet, 0);
        TRootedIterator rootedIterator = iter(stringSet, 0);
        SEQAN_ASSERT(iterator == begin(stringSet));
        SEQAN_ASSERT(standardIterator == begin(stringSet, Standard()));
        SEQAN_ASSERT(rootedIterator == begin(stringSet, Rooted()));
    }

    // Test on a non empty stringSet.
    {
    	TNonConstStringSet nonConstStringSet;
    	resize(nonConstStringSet, 3);
    	nonConstStringSet[0] = "AAAA";
    	nonConstStringSet[1] = "CCCC";
    	nonConstStringSet[2] = "GGGG";
    	TStringSet stringSet(nonConstStringSet);
        TIterator iterator = iter(stringSet, 0);
        TStandardIterator standardIterator = iter(stringSet, 0);
        TRootedIterator rootedIterator = iter(stringSet, 0);
        SEQAN_ASSERT_EQ(getValue(iterator), "AAAA");
        SEQAN_ASSERT(getValue(iterator) == getValue(stringSet, 0));
        SEQAN_ASSERT(getValue(standardIterator) == getValue(stringSet, 0));
        SEQAN_ASSERT(getValue(rootedIterator) == getValue(stringSet, 0));
    }

    // Test on a non empty stringSet.
    {
    	TNonConstStringSet nonConstStringSet;
    	resize(nonConstStringSet, 4);
    	nonConstStringSet[0] = "AAAA";
    	nonConstStringSet[1] = "CCCC";
    	nonConstStringSet[2] = "GGGG";
    	nonConstStringSet[3] = "TTTT";
    	TStringSet stringSet(nonConstStringSet);
        TIterator iterator = iter(stringSet, 3);
        TStandardIterator standardIterator = iter(stringSet, 3);
        TRootedIterator rootedIterator = iter(stringSet, 3);
        SEQAN_ASSERT_EQ(getValue(iterator), "TTTT");
        SEQAN_ASSERT(getValue(iterator) == getValue(stringSet, 3));
        SEQAN_ASSERT(getValue(standardIterator) == getValue(stringSet, 3));
        SEQAN_ASSERT(getValue(rootedIterator) == getValue(stringSet, 3));
    }
}

// Test of length().
template <typename TStringSet>
void testStringSetLength(TStringSet & /*Tag*/)
{
    using namespace seqan;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;

    // Test on an empty string.
    {
    	TStringSet stringSet;
    	SEQAN_ASSERT_EQ(length(stringSet), 0u);
    }

    // Test on a non empty string.
    {
    	TNonConstStringSet nonConstStringSet;
    	resize(nonConstStringSet, 10);
    	TStringSet stringSet(nonConstStringSet);
    	SEQAN_ASSERT_EQ(length(stringSet), 10u);
    }
}

// Test of moveValue().
template <typename TStringSet>
void testStringSetMoveValue(TStringSet & /*Tag*/)
{
    using namespace seqan;

    TStringSet stringSet;

    resize(stringSet, 2);
    moveValue(stringSet, 1, "ACGT");
    SEQAN_ASSERT_EQ(stringSet[1], "ACGT");
}

// TODO (singer): see infix.
// Test of prefix().
template <typename TStringSet>
void testStringSetPrefix(TStringSet & /*Tag*/)
{
    using namespace seqan;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;
    typedef typename Value<TStringSet>::Type TString;

    TNonConstStringSet nonConstStringSet;
    resize(nonConstStringSet, 10);
    nonConstStringSet[0] = "ACGTACGT";
    TStringSet stringSet(nonConstStringSet);
    TString pref = prefix(stringSet, 3);
    SEQAN_ASSERT_EQ(pref, "ACG");
}

// TODO (singer); replace is not defined for string sets.
//// Test of replace().
//template <typename TStringSet>
//void testStringSetReplace(TStringSet & /*Tag*/)
//{
//    using namespace seqan;
//
//    TStringSet stringSet1;
//    TStringSet stringSet2;
//
//    // This is problematic according to the documentation.
//    // 0 can be a position or an iterator causing compiler errors.
//    replace(stringSet1, 0, 0, stringSet2);
//    SEQAN_ASSERT_EQ(stringSet1, TStringSet());
//
//    appendValue(stringSet1, "AAAA");
//    appendValue(stringSet1, "CCCC");
//    appendValue(stringSet1, "GGGG");
//    appendValue(stringSet1, "TTTT");
//
//    replace(stringSet1, 1, 3, stringSet2);
//    SEQAN_ASSERT_EQ(length(stringSet1), 2u);
//    SEQAN_ASSERT_EQ(stringSet1[0], "AAAA");
//    SEQAN_ASSERT_EQ(stringSet1[1], "TTTT");
//
//    replace(stringSet2, 0, 0, stringSet1);
//    SEQAN_ASSERT_EQ(stringSet2[0], "AAAA");
//    SEQAN_ASSERT_EQ(stringSet2[1], "TTTT");
//
//    clear(stringSet1);
//    clear(stringSet2);
//    appendValue(stringSet1, "AAAA");
//    appendValue(stringSet1, "CCCC");
//    appendValue(stringSet2, "GGGG");
//    appendValue(stringSet2, "TTTT");
//
//    replace(stringSet1, 1, 2, stringSet2);
//    SEQAN_ASSERT_EQ(stringSet1[0], "AAAA");
//    SEQAN_ASSERT_EQ(stringSet2[1], "GGGG");
//    SEQAN_ASSERT_EQ(stringSet2[2], "TTTT");
//
//}

// Test of resize().
template <typename TStringSet>
void testStringSetResize(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TStringSet>::Type TString;

    TStringSet stringSet;

    resize(stringSet, 0);
    SEQAN_ASSERT_EQ(length(stringSet), 0u);

    resize(stringSet, 10);
    SEQAN_ASSERT_EQ(length(stringSet), 10u);
    // TODO (singer): resize should initialize newly allocated memory,
    // which it does not at the moment!
    SEQAN_ASSERT_EQ(stringSet[0], TString());
    SEQAN_ASSERT_EQ(stringSet[9], TString());

    resize(stringSet, 0);
    SEQAN_ASSERT_EQ(length(stringSet), 0u);

    // TODO (singer): resize with a provided value is not possible with string sets.
//    TString string = "ACGT";
//    resize(stringSet, 10, string);
//    SEQAN_ASSERT_EQ(stringSet[0], string);
//    SEQAN_ASSERT_EQ(stringSet[9], string);
}

// TODO (singer): see infix.
// Test of suffix().
template <typename TStringSet>
void testStringSetSuffix(TStringSet & /*Tag*/)
{
    using namespace seqan;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;
    typedef typename Value<TStringSet>::Type TString;

    TNonConstStringSet nonConstStringSet;
    resize(nonConstStringSet, 10);
    nonConstStringSet[9] = "ACGTACGT";
    TStringSet stringSet(nonConstStringSet);
    TString pref = suffix(stringSet, 5);
    SEQAN_ASSERT_EQ(pref, "CGT");
}

// TODO (singer): swap is not working because a constructor string set (StringSet(StringSet, Move)) is missing.
// Test of swap().
template <typename TStringSet>
void testStringSetSwap(TStringSet & /*Tag*/)
{
//    using namespace seqan;
//
//    TStringSet stringSet1;
//    TStringSet stringSet2;
//    TStringSet stringSet3 = stringSet1;
//    TStringSet stringSet4 = stringSet2;
//
//    swap(stringSet1, stringSet2);
//    SEQAN_ASSERT(stringSet1 == stringSet4);
//    SEQAN_ASSERT(stringSet2 == stringSet3);
//
//    appendValue(stringSet1, "ACGT");
//    stringSet3 = stringSet1;
//    stringSet4 = stringSet2;
//
//    swap(stringSet1, stringSet2);
//    SEQAN_ASSERT(stringSet1 == stringSet4);
//    SEQAN_ASSERT(stringSet2 == stringSet3);
//
//    clear(stringSet1);
//    appendValue(stringSet2, "ACGT");
//    stringSet3 = stringSet1;
//    stringSet4 = stringSet2;
//
//    swap(stringSet1, stringSet2);
//    SEQAN_ASSERT(stringSet1 == stringSet4);
//    SEQAN_ASSERT(stringSet2 == stringSet3);
//
//    stringSet1[0] = "ACAC";
//    stringSet2[0] = "GTGT";
//    stringSet3 = stringSet1;
//    stringSet4 = stringSet2;
//
//    swap(stringSet1, stringSet2);
//    SEQAN_ASSERT(stringSet1 == stringSet4);
//    SEQAN_ASSERT(stringSet2 == stringSet3);
}

// Test of value().
template <typename TStringSet>
void testStringSetValue(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TStringSet>::Type TString;

    // In contrast to getValue(), value() does not return a copy but a reference.
    // We test this using the variable value_.
    TStringSet stringSet;
    resize(stringSet, 3);
    stringSet[0] = "ACAC";
    stringSet[1] = "AAAA";
    stringSet[2] = "TTTT";
    TString & value_ = value(stringSet, 0);
    SEQAN_ASSERT_EQ(value_, "ACAC");

    value_ = "GGGG";
    SEQAN_ASSERT_EQ(value_, "GGGG");
    SEQAN_ASSERT_EQ(stringSet[0], "GGGG");
}

// Test of valueById().
template <typename TStringSet>
void testStringSetValueById(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TStringSet>::Type TString;

    // In contrast to getValue(), value() does not return a copy but a reference.
    // We test this using the variable value_.
    TStringSet stringSet;
    resize(stringSet, 3);
    stringSet[0] = "ACAC";
    stringSet[1] = "AAAA";
    stringSet[2] = "TTTT";
    TString & value_ = valueById(stringSet, 0);
    SEQAN_ASSERT_EQ(value_, "ACAC");

    value_ = "GGGG";
    SEQAN_ASSERT_EQ(value_, "GGGG");
    SEQAN_ASSERT_EQ(stringSet[0], "GGGG");
}

// --------------------------------------------------------------------------
// Testing Alloc Strings With Simple Types
// --------------------------------------------------------------------------

// Test whether sequences are copy constructible.
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_copy_constructible)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetCopyConstructible(tag);
}

// Test whether sequences are default constructible.
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_default_constructible)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testStringSetDefaultConstructible(tag);

    String<Dna, Alloc<> > const constTag;
    testStringSetDefaultConstructible(constTag);
}

// Test whether sequences are default constructible.
// TODO (singer): comparison operators are not implemented for string sets.
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_less_greater_equal)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    //testStringSetLessGreaterEqual(tag);

    StringSet<String<Dna, Alloc<> > > const constTag;
    //testStringSetLessGreaterEqual(constTag);
}

// Test of append().
// TODO (singer): append() is not implemented for string sets.
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_append)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    //testStringSetAppend(tag);
}

// Test of appendValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_append_value)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetAppendValue(tag);
}

// Test whether sequences are assignable.
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_assign)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetAssign(tag);
}

// Test of assignValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_assign_value)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetAssignValue(tag);
}

// Test of assignValueById().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_assign_value_by_id)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetAssignValueById(tag);
}

// Test of back().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_back)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetBack(tag);

    StringSet<String<Dna, Alloc<> > > const constTag;
    testStringSetBack(constTag);
}

// Test of begin().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_begin)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetBegin(tag);

    StringSet<String<Dna, Alloc<> > > const constTag;
    testStringSetBegin(constTag);
}

// Test of beginPosition().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_begin_position)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetBeginPosition(tag);

    StringSet<String<Dna, Alloc<> > > const constTag;
    testStringSetBeginPosition(constTag);
}

//// Test of capacity().
//SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_capacity)
//{
//    using namespace seqan;
//
//    String<Dna, Alloc<> > tag;
//    testSequenceCapacity(tag);
//
//    String<Dna, Alloc<> > const constTag;
//    testSequenceCapacity(constTag);
//}
//
// Test of clear().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_clear)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetClear(tag);
}

// Test of end().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_concat)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetConcat(tag);

//    StringSet<String<Dna, Alloc<> > > const constTag;
//    testStringSetConcat(constTag);
}

// Test of end().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_end)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetEnd(tag);

    StringSet<String<Dna, Alloc<> > > const constTag;
    testStringSetEnd(constTag);
}

// Test of endPosition().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_end_position)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetEndPosition(tag);

    StringSet<String<Dna, Alloc<> > > const constTag;
    testStringSetEndPosition(constTag);
}

// Test of erase().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_erase)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetErase(tag);
}

// Test of erase().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_erase_back)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetEraseBack(tag);
}
//// Test of eraseBack().
//SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_erase_back)
//{
//    using namespace seqan;
//
//    String<Dna, Alloc<> > tag;
//    testSequenceEraseBack(tag);
//}
//
// Test of front().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_front)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetFront(tag);

    StringSet<String<Dna, Alloc<> > > const constTag;
    testStringSetFront(constTag);
}

// Test of getValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_get_value)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetGetValue(tag);

    StringSet<String<Dna, Alloc<> > > const constTag;
    testStringSetGetValue(constTag);
}

// Test of getValueById().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_get_value_by_id)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetGetValueById(tag);

    // TODO (singer): not defined for const string sets.
    //StringSet<String<Dna, Alloc<> > > const constTag;
    //testStringSetGetValueById(constTag);
}

// Test of getValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_infix)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetInfix(tag);

    StringSet<String<Dna, Alloc<> > > const constTag;
    testStringSetInfix(constTag);
}

// Test of getValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_infix_with_length)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetInfixWithLength(tag);

    StringSet<String<Dna, Alloc<> > > const constTag;
    testStringSetInfixWithLength(constTag);
}


// Test of insert().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_insert)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetInsert(tag);
}

// Test of insertValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_insert_value)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetInsertValue(tag);
}

// Test of iter().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_iter)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetIter(tag);

    StringSet<String<Dna, Alloc<> > > const constTag;
    testStringSetIter(constTag);
}

// Test of length().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_length)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetLength(tag);

    StringSet<String<Dna, Alloc<> > > const constTag;
    testStringSetLength(constTag);
}

// Test of moveValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_move_value)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetMoveValue(tag);
}

// Test of prefix().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_prefix)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetPrefix(tag);

    StringSet<String<Dna, Alloc<> > > const constTag;
    testStringSetPrefix(constTag);
}

//// Test of replace().
//SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_replace)
//{
//    using namespace seqan;
//
//    StringSet<String<Dna, Alloc<> > > tag;
//    testStringSetReplace(tag);
//}

// Test of resize().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_resize)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetResize(tag);
}

// Test of suffix().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_suffix)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetSuffix(tag);

    StringSet<String<Dna, Alloc<> > > const constTag;
    testStringSetSuffix(constTag);
}

// Test of swap().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_swap)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetSwap(tag);
}

// Test of value().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_value)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetValue(tag);
}

// Test of valueById().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_set_dna_value_by_id)
{
    using namespace seqan;

    StringSet<String<Dna, Alloc<> > > tag;
    testStringSetValueById(tag);
}
#endif // CORE_TESTS_SEQUENCE_TEST_STRINGSET_H_
