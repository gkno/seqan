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

#ifndef CORE_TESTS_SEQUENCE_TEST_SEQUENCE_H_
#define CORE_TESTS_SEQUENCE_TEST_SEQUENCE_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

// --------------------------------------------------------------------------
// CountingChar is used to test sequences of non simple data types.
// --------------------------------------------------------------------------

using namespace seqan;

struct CountingChar
{
    char value;                     // value of the object
    static unsigned numConstruct;   // number of constructor calls
    static unsigned numDeconstruct; // number of destructor calls

    CountingChar()
    {
        numConstruct += 1;
    }

    CountingChar(char const & value) : value(value)
    {
        numConstruct += 1;
    }

    CountingChar(CountingChar const & other) : value(other.value)
    {
        numConstruct += 1;
    }

    ~CountingChar()
    {
        numDeconstruct += 1;
    }

    static void clear()
    {
        numConstruct = 0;
        numDeconstruct = 0;
    }

    bool operator==(CountingChar const & other) const
    {
        return value == other.value;
    }

    bool operator>(CountingChar const & other) const
    {
        return value > other.value;
    }

    bool operator<(CountingChar const & other) const
    {
        return value < other.value;
    }
};

template <typename TStream>
inline TStream & operator<<(TStream & stream, CountingChar const & countingChar)
{
    stream << countingChar.value;

    return stream;
}

template <typename TStream, typename TSpec>
inline TStream & operator<<(TStream & stream, seqan::String<CountingChar, TSpec> const & string)
{
    for (unsigned i = 0; i < length(string); ++i)
        stream << string[i];

    return stream;
}

unsigned CountingChar::numConstruct = 0;
unsigned CountingChar::numDeconstruct = 0;

template <typename TValue>
void testConstructDeconstruct(TValue const & /*tag*/) {}

void testConstructDeconstruct(CountingChar const & /*tag*/)
{
    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct + 1u);
}

template <typename TValue>
void testConstruct(TValue const & /*tag*/) {}

void testConstruct(CountingChar const & /*tag*/)
{
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 1u);
}

template <typename TAlphabetSpecPair>
class StringTest : public seqan::Test
{
public:
    typedef typename seqan::TagListValue<TAlphabetSpecPair, 1>::Type TAlphabet;
    typedef typename seqan::TagListValue<TAlphabetSpecPair, 2>::Type TSpec;
    typedef seqan::String<TAlphabet, TSpec> TString;
};

// --------------------------------------------------------------------------
// Types which are tested for very common functionality
// --------------------------------------------------------------------------

template <typename T>
class StringTestCommon : public StringTest<T>
{};

typedef seqan::TagList<
    seqan::TagList<seqan::Dna5, seqan::TagList<seqan::External<> > >, seqan::TagList<
//     seqan::TagList<char,        seqan::TagList<seqan::External<> > >, seqan::TagList<
    seqan::TagList<int,         seqan::TagList<seqan::External<> > >, seqan::TagList<
    seqan::TagList<CountingChar,seqan::TagList<seqan::External<> > >, seqan::TagList<
    seqan::TagList<seqan::Dna5, seqan::TagList<seqan::MMap<> > >, seqan::TagList<
//     seqan::TagList<char,        seqan::TagList<seqan::MMap<> > >, seqan::TagList<
    seqan::TagList<int,         seqan::TagList<seqan::MMap<> > >, seqan::TagList<
    seqan::TagList<CountingChar,seqan::TagList<seqan::MMap<> > >, seqan::TagList<
    seqan::TagList<seqan::Dna5, seqan::TagList<seqan::Packed<> > >, seqan::TagList<
//     seqan::TagList<char,        seqan::TagList<seqan::Packed<> > >, seqan::TagList<
    seqan::TagList<int,         seqan::TagList<seqan::Packed<> > >, seqan::TagList<
    seqan::TagList<seqan::Dna5, seqan::TagList<seqan::Array<100> > >, seqan::TagList<
//     seqan::TagList<char,        seqan::TagList<seqan::Array<100> > >, seqan::TagList<
    seqan::TagList<int,         seqan::TagList<seqan::Array<100> > >, seqan::TagList<
    seqan::TagList<CountingChar,seqan::TagList<seqan::Array<100> > >, seqan::TagList<
    seqan::TagList<seqan::Dna5, seqan::TagList<seqan::Block<> > >, seqan::TagList<
//     seqan::TagList<char,        seqan::TagList<seqan::Block<> > >, seqan::TagList<
    seqan::TagList<int,         seqan::TagList<seqan::Block<> > >, seqan::TagList<
    seqan::TagList<CountingChar,seqan::TagList<seqan::Block<> > >, seqan::TagList<
    seqan::TagList<seqan::Dna5, seqan::TagList<seqan::Alloc<> > >, seqan::TagList<
//     seqan::TagList<char,        seqan::TagList<seqan::Alloc<> > >, seqan::TagList<
    seqan::TagList<int,         seqan::TagList<seqan::Alloc<> > >, seqan::TagList<
    seqan::TagList<CountingChar,seqan::TagList<seqan::Alloc<> > >//, seqan::TagList<
//     seqan::TagList<char,        seqan::TagList<seqan::CStyle> >
    > > > > > > > > > > > > > > > > > //> > > > > > >
    StringTestCommonTypes;

SEQAN_TYPED_TEST_CASE(StringTestCommon, StringTestCommonTypes);

// Test whether sequences are default constructible.
template <typename TString>
void testSequenceDefaultConstructible(TString & /*Tag*/)
{
    CountingChar::clear();

    TString string;
    SEQAN_ASSERT(begin(string) == end(string));
}

SEQAN_TYPED_TEST(StringTestCommon, DefaultConstructible)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceDefaultConstructible(str);

    typename TestFixture::TString const constStr;
    testSequenceDefaultConstructible(constStr);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test whether sequences are copy constructible.
template <typename TString>
void testSequenceCopyConstructible(TString & /*Tag*/)
{
    using namespace seqan;

    TString string1("ACGCTAGCAT");
    TString string2(string1);
    SEQAN_ASSERT_EQ(string1, string2);
}

SEQAN_TYPED_TEST(StringTestCommon, CopyConstructible)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceCopyConstructible(str);

    typename TestFixture::TString const constStr;
    testSequenceCopyConstructible(constStr);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test operator<().
template <typename TString>
void testSequenceLess(TString & /*Tag*/)
{
    using namespace seqan;

    // Nothing is smaller than something.
    {
        TString string1;
        TString string2("A");
        SEQAN_ASSERT_EQ(string1 < string2, true);
    }

    // Equal is not smaller
    {
        TString string1("A");
        TString string2("A");
        SEQAN_ASSERT_EQ(string1 < string2, false);
    }

    // Sequences with lex. smaller characters and equal length are smaller.
    {
        TString string1("A");
        TString string2("C");
        SEQAN_ASSERT_EQ(string1 < string2, true);
    }

    // Sequences with lex. smaller characters but larger length are smaller.
    {
        TString string1("AA");
        TString string2("C");
        SEQAN_ASSERT_EQ(string1 < string2, true);
    }

    // Sequences with equal characters but smaller length are smaller.
    {
        TString string1("AA");
        TString string2("AAA");
        SEQAN_ASSERT_EQ(string1 < string2, true);
    }
}

SEQAN_TYPED_TEST(StringTestCommon, LessOperator)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceLess(str);

    typename TestFixture::TString const constStr;
    testSequenceLess(constStr);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test operator<=().
template <typename TString>
void testSequenceLessEqual(TString & /*Tag*/)
{
    using namespace seqan;

    // Nothing is equal to nothing.
    {
        TString string1;
        TString string2;
        SEQAN_ASSERT_EQ(string1 <= string2, true);
    }

    // Nothing is smaller than something.
    {
        TString string1;
        TString string2("A");
        SEQAN_ASSERT_EQ(string1 <= string2, true);
    }

    // Sequences of equal characters and length are equal.
    {
        TString string1("A");
        TString string2("A");
        SEQAN_ASSERT_EQ(string1 <= string2, true);
    }

    // Sequences with lex. smaller characters and equal length are smaller.
    {
        TString string1("A");
        TString string2("C");
        SEQAN_ASSERT_EQ(string1 <= string2, true);
    }

    // Sequences with lex. smaller characters but larger length are smaller.
    {
        TString string1("AA");
        TString string2("C");
        SEQAN_ASSERT_EQ(string1 <= string2, true);
    }

    // Sequences with equal characters but smaller length are smaller.
    {
        TString string1("AA");
        TString string2("AAA");
        SEQAN_ASSERT_EQ(string1 <= string2, true);
    }
}

SEQAN_TYPED_TEST(StringTestCommon, LessEqualOperator)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceLessEqual(str);

    typename TestFixture::TString const constStr;
    testSequenceLessEqual(constStr);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test operator>().
template <typename TString>
void testSequenceGreater(TString & /*Tag*/)
{
    using namespace seqan;

    // Something is greater than nothing.
    {
        TString string1("A");
        TString string2;
        SEQAN_ASSERT_EQ(string1 > string2, true);
    }

    // Equal is not greater
    {
        TString string1("A");
        TString string2("A");
        SEQAN_ASSERT_EQ(string1 > string2, false);
    }

    // Sequences with lex. greater characters and equal length are greater.
    {
        TString string1("C");
        TString string2("A");
        SEQAN_ASSERT_EQ(string1 > string2, true);
    }

    // Sequences with equal characters but larger length are larger.
    {
        TString string1("AA");
        TString string2("A");
        SEQAN_ASSERT_EQ(string1 > string2, true);
    }

    // Sequences with lex. greater characters are greater.
    {
        TString string1("C");
        TString string2("AA");
        SEQAN_ASSERT_EQ(string1 > string2, true);
    }
}

SEQAN_TYPED_TEST(StringTestCommon, GreaterOperator)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceGreater(str);

    typename TestFixture::TString const constStr;
    testSequenceGreater(constStr);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test operator>=().
template <typename TString>
void testSequenceGreaterEqual(TString & /*Tag*/)
{
    using namespace seqan;

    // Nothing is equal to nothing.
    {
        TString string1;
        TString string2;
        SEQAN_ASSERT_EQ(string1 >= string2, true);
    }

    // Something is greater than nothing.
    {
        TString string1("A");
        TString string2;
        SEQAN_ASSERT_EQ(string1 >= string2, true);
    }

    // Sequences of equal characters and length are equal.
    {
        TString string1("A");
        TString string2("A");
        SEQAN_ASSERT_EQ(string1 >= string2, true);
    }

    // Sequences with lex. greater characters and equal length are greater.
    {
        TString string1("C");
        TString string2("A");
        SEQAN_ASSERT_EQ(string1 >= string2, true);
    }

    // Sequences with equal characters but larger length are greater.
    {
        TString string1("AA");
        TString string2("A");
        SEQAN_ASSERT_EQ(string1 >= string2, true);
    }

    // Sequences with lex. greater characters are greater.
    {
        TString string1("C");
        TString string2("AA");
        SEQAN_ASSERT_EQ(string1 > string2, true);
    }
}

SEQAN_TYPED_TEST(StringTestCommon, GreaterEqualOperator)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceGreaterEqual(str);

    typename TestFixture::TString const constStr;
    testSequenceGreaterEqual(constStr);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test operator==().
template <typename TString>
void testSequenceEqual(TString & /*Tag*/)
{
    using namespace seqan;

    // Nothing is equal to nothing.
    {
        TString string1;
        TString string2;
        SEQAN_ASSERT_EQ(string1 == string2, true);
    }

    // Something is greater than nothing.
    {
        TString string1("A");
        TString string2;
        SEQAN_ASSERT_EQ(string1 == string2, false);
    }
    {
        TString string1;
        TString string2("A");
        SEQAN_ASSERT_EQ(string1 == string2, false);
    }

    // Sequences of equal characters and length are equal.
    {
        TString string1("A");
        TString string2("A");
        SEQAN_ASSERT_EQ(string1 == string2, true);
    }

    // Sequences of equal characters but different length are not equal.
    {
        TString string1("A");
        TString string2("AA");
        SEQAN_ASSERT_EQ(string1 == string2, false);
    }
    {
        TString string1("AA");
        TString string2("A");
        SEQAN_ASSERT_EQ(string1 == string2, false);
    }

    // Sequences of different characters are not equal.
    {
        TString string1("A");
        TString string2("C");
        SEQAN_ASSERT_EQ(string1 == string2, false);
    }
    {
        TString string1("C");
        TString string2("A");
        SEQAN_ASSERT_EQ(string1 == string2, false);
    }
}

SEQAN_TYPED_TEST(StringTestCommon, EqualOperator)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceEqual(str);

    typename TestFixture::TString const constStr;
    testSequenceEqual(constStr);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test operator!=().
template <typename TString>
void testSequenceUnequal(TString & /*Tag*/)
{
    using namespace seqan;

    // Nothing is equal to nothing.
    {
        TString string1;
            TString string2;
            SEQAN_ASSERT_EQ(string1 != string2, false);
    }

    // Something is greater than nothing.
    {
        TString string1("A");
        TString string2;
        SEQAN_ASSERT_EQ(string1 != string2, true);
    }
    {
        TString string1;
        TString string2("A");
        SEQAN_ASSERT_EQ(string1 != string2, true);
    }

    // Sequences of equal characters and length are equal.
    {
        TString string1("A");
        TString string2("A");
        SEQAN_ASSERT_EQ(string1 != string2, false);
    }

    // Sequences of equal characters but different length are not equal.
    {
        TString string1("A");
        TString string2("AA");
        SEQAN_ASSERT_EQ(string1 != string2, true);
    }
    {
        TString string1("AA");
        TString string2("A");
        SEQAN_ASSERT_EQ(string1 != string2, true);
    }

    // Sequences of different characters are not equal.
    {
        TString string1("A");
        TString string2("C");
        SEQAN_ASSERT_EQ(string1 != string2, true);
    }
    {
        TString string1("C");
        TString string2("A");
        SEQAN_ASSERT_EQ(string1 != string2, true);
    }
}

SEQAN_TYPED_TEST(StringTestCommon, UnequalOperator)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceUnequal(str);

    typename TestFixture::TString const constStr;
    testSequenceUnequal(constStr);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of append().
template <typename TString>
void testSequenceAppend(TString & /*Tag*/)
{
    using namespace seqan;

    // Test the append function on empty strings
    {
        TString string1("");
        TString string2("");
        append(string1, string2);
        SEQAN_ASSERT_EQ(string1, "");
    }

    // Test the append function on one empty and one non empty string
    {
        TString string1("ACGTACGTACGT");
        TString string2("");
        append(string1, string2);
        SEQAN_ASSERT_EQ(string1, "ACGTACGTACGT");
    }

    // Test the append function on one empty and one non empty string
    {
        TString string1("");
        TString string2("TTGGATTAACC");
        append(string1, string2);
        SEQAN_ASSERT_EQ(string1, "TTGGATTAACC");
    }

    // Test the append function on two non empty strings.
    {
        TString string1("ACGTACGTACGT");
        TString string2("TTGGATTAACCC");
        append(string1, string2);
        SEQAN_ASSERT_EQ(string1, "ACGTACGTACGTTTGGATTAACCC");
    }
}

// TODO(singer): no viable overloaded '='
template <typename TValue>
void testSequenceAppend(seqan::String<TValue, seqan::Packed<> > & /*Tag*/) {}

SEQAN_TYPED_TEST(StringTestCommon, Append)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceAppend(str);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of appendValue().
template <typename TString>
void testSequenceAppendValue(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

    TString string("");

    // Test the appendValue function
    TString string2("A");
    TValue value = 'A';
    appendValue(string, value);
    SEQAN_ASSERT(string == string2);

    // Test the appendValue function
    TString string3("AA");
    appendValue(string, 'A');
    SEQAN_ASSERT(string == string3);
}

// TODO(singer): error: no matching function for call to '_setLength'
template <typename TValue>
void testSequenceAppendValue(String<TValue, External<> > & /*Tag*/) {}

SEQAN_TYPED_TEST(StringTestCommon, AppendValue)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceAppendValue(str);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test whether sequences are assignable.
template <typename TString>
void testSequenceAssignable(TString & /*Tag*/)
{
    using namespace seqan;

    {
        // Test on an empty string.
        TString string1("");
        TString string3 = "";

        // Test assign() function.
        {
            TString string2;
            assign(string2, string1);
            SEQAN_ASSERT_EQ(string1, string2);
        }

        // Test operator=().
        {
            TString string2;  // Separate definition and assignment on purpose.
            string2 = string1;
            SEQAN_ASSERT_EQ(string1, string2);
        }
    }
    {
        // Test the basic concept on a non empty string.
        TString string1("ACGTACGTACGT");

        // Test the assign() function.
        {
            TString string2;
            assign(string2, string1);
            SEQAN_ASSERT_EQ(string1, string2);
        }

        // Test operator=().
        {
            TString string2;  // Separate definition and assignment on purpose.
            string2 = string1;
            SEQAN_ASSERT_EQ(string1, string2);
        }
    }
}
 
// TODO(singer): error: no viable conversion from
template <typename TValue>
void testSequenceAssignable(String<TValue, Block<> > & /*Tag*/) {}

SEQAN_TYPED_TEST(StringTestCommon, Assignable)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceAssignable(str);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of assignValue().
template <typename TString>
void testSequenceAssignValue(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

    TString string("AA");

    assignValue(string, 1, TValue('G'));
    SEQAN_ASSERT_EQ(string[1], TValue('G'));
}

template <typename TValue>
void testSequenceAssignValue(String<TValue, Packed<> > & /*Tag*/) {}

SEQAN_TYPED_TEST(StringTestCommon, AssignValue)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceAssignValue(str);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// We need two back() tests, since back() returns a reference or a copy
// (depending on TString). We check whether we can modify the reference.
// Test of back() for non const strings.
template <typename TString>
void testSequenceBack(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    TString string("ACGT");

    // val is a reference in contrast to the const version of back().
    TValue & val = back(string);
    val = 'A';
    SEQAN_ASSERT_EQ(val, string[length(string) - 1]);
}

// Test of back() for const strings.
template <typename TString>
void testSequenceBack(TString const & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    TString const string("ACGT");

    // val is not a reference in contrast to the non const version of back().
    TValue val = back(string);
    SEQAN_ASSERT_EQ(val, string[length(string) - 1]);
}

template <typename TValue>
void testSequenceBack(String<TValue, Packed<> > & /*Tag*/) {}
template <typename TValue>
void testSequenceBack(String<TValue, Packed<> > const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringTestCommon, Back)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceBack(str);

    typename TestFixture::TString const constStr;
    testSequenceBack(constStr);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of begin().
template <typename TString>
void testSequenceBegin(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

    TString string("ACGT");
    SEQAN_ASSERT_EQ(*begin(string), TValue('A'));
    SEQAN_ASSERT_EQ(*begin(string, Standard()), TValue('A'));
    SEQAN_ASSERT_EQ(*begin(string, Rooted()), TValue('A'));
}

template <typename TValue>
void testSequenceBegin(String<TValue, MMap<> > & /*Tag*/) {}
template <typename TValue>
void testSequenceBegin(String<TValue, MMap<> > const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringTestCommon, Begin)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceBegin(str);

    typename TestFixture::TString const constStr;
    testSequenceBegin(constStr);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of beginPosition().
template <typename TString>
void testSequenceBeginPosition(TString & /*Tag*/)
{
    using namespace seqan;

    // Test on an empty string.
    TString string1("");
    SEQAN_ASSERT_EQ(beginPosition(string1), 0u);

    // Test on a non empty string.
    TString string2("ACGT");
    SEQAN_ASSERT_EQ(beginPosition(string2), 0u);
}

SEQAN_TYPED_TEST(StringTestCommon, BeginPosition)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceBeginPosition(str);

    typename TestFixture::TString const constStr;
    testSequenceBeginPosition(constStr);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of capacity().
template <typename TString>
void testSequenceCapacity(TString & /*Tag*/)
{
    using namespace seqan;

    // Test on an empty string.
    TString string1("");
    SEQAN_ASSERT_GEQ(capacity(string1), length(string1));

    // Test on a non empty string.
    TString string2("ACGTACGTACGT");
    SEQAN_ASSERT_GEQ(capacity(string2), length(string2));
}

SEQAN_TYPED_TEST(StringTestCommon, Capacity)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceCapacity(str);

    typename TestFixture::TString const constStr;
    testSequenceCapacity(constStr);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of clear().
template <typename TString>
void testSequenceClear(TString & /*Tag*/)
{
    using namespace seqan;

    {
        // Test on an empty string.
        TString string("");
        clear(string);
        SEQAN_ASSERT(begin(string) == end(string));
        SEQAN_ASSERT_EQ(capacity(string), 0u);
    }
    {
        // Test on a non empty string.
        TString string("ACGTACGTACGT");
        clear(string);
        SEQAN_ASSERT_EQ(string, TString());
    }
}

SEQAN_TYPED_TEST(StringTestCommon, Clear)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceClear(str);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of end().
template <typename TString>
void testSequenceEnd(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

    TString string("ACGT");
    typename Iterator<TString>::Type iter = end(string);
    --iter;
    SEQAN_ASSERT_EQ(*iter, TValue('T'));

    typename Iterator<TString, Standard>::Type standardIter = end(string, Standard());
    --standardIter;
    SEQAN_ASSERT_EQ(*standardIter, TValue('T'));

    typename Iterator<TString, Rooted>::Type rootedIter = end(string, Rooted());
    --rootedIter;
    SEQAN_ASSERT_EQ(*rootedIter, TValue('T'));
}

// TODO(singer): Seg fault
template <typename TValue>
void testSequenceEnd(String<TValue, MMap<> > & /*Tag*/) {}
template <typename TValue>
void testSequenceEnd(String<TValue, MMap<> > const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringTestCommon, End)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceEnd(str);
    
    typename TestFixture::TString const constStr;
    testSequenceEnd(constStr);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of endPosition().
template <typename TString>
void testSequenceEndPosition(TString & /*Tag*/)
{
    using namespace seqan;

    // Test on an empty string.
    TString string1("");
    SEQAN_ASSERT_EQ(endPosition(string1), length(string1));

    // Test on a non empty string.
    TString string2("ACGT");
    SEQAN_ASSERT_EQ(endPosition(string2), length(string2));
}

SEQAN_TYPED_TEST(StringTestCommon, EndPosition)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceEndPosition(str);

    typename TestFixture::TString const constStr;
    testSequenceEndPosition(constStr);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of erase().
template <typename TString>
void testSequenceErase(TString & /*Tag*/)
{
    using namespace seqan;

    // Test on a non empty string.
    TString string("ACGTACGTACGT");
    erase(string, 1);
    SEQAN_ASSERT_EQ(string, "AGTACGTACGT");

    erase(string, 2, 5);
    SEQAN_ASSERT_EQ(string, "AGGTACGT");
}

template <typename TValue>
void testSequenceErase(String<TValue, External<> > & /*Tag*/) {}
template <typename TValue>
void testSequenceErase(String<TValue, MMap<> > & /*Tag*/) {}
template <typename TValue>
void testSequenceErase(String<TValue, Block<> > & /*Tag*/) {}
template <typename TValue>
void testSequenceErase(String<TValue, Packed<> > & /*Tag*/) {}

SEQAN_TYPED_TEST(StringTestCommon, Erase)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceErase(str);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of eraseBack().
template <typename TString>
void testSequenceEraseBack(TString & /*Tag*/)
{
    using namespace seqan;

    {
        // Test on an empty string.
        TString string("");

        // TODO (singer): eraseBack() on an empty string is not possible.
        // Error message: "String must have more than 0 characters in eraseBack()!".
        // If erase() should work on empty strings than eraseBack() as well ???
        // eraseBack(string);
        // SEQAN_ASSERT_EQ(string, "");
    }
    {
        // Test on a non empty string.
        TString string("ACGTACGTACGT");
        TString string2("ACGTACGTACG");
        eraseBack(string);
        SEQAN_ASSERT_EQ(string, string2);

        TString string3("A");
        TString string4;
        eraseBack(string3);
        SEQAN_ASSERT_EQ(string3, string4);
    }
}

SEQAN_TYPED_TEST(StringTestCommon, EraseBack)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceEraseBack(str);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of front() for non const strings.
template <typename TString>
void testSequenceFront(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    TString string("ACGT");

    // val is a reference in contrast to the const version of front()
    TValue & val = front(string);
    val = 'A';
    SEQAN_ASSERT_EQ(val, string[0]);
}

// Test of front() for const strings.
template <typename TString>
void testSequenceFront(TString const & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    TString string("ACGT");

    // val is not a reference in contrast to the non const version of front()
    TValue val = front(string);
    SEQAN_ASSERT_EQ(val, string[0]);
}

template <typename TValue>
void testSequenceFront(String<TValue, Packed<> > & /*Tag*/) {}
template <typename TValue>
void testSequenceFront(String<TValue, Packed<> > const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringTestCommon, Front)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceFront(str);

    typename TestFixture::TString const constStr;
    testSequenceFront(constStr);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of getValue().
template <typename TString>
void testSequenceGetValue(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

     // In contrast to value(), getValue() does not return a reference but a copy.
     // We test this using the variable value_.
    TString string("ACGT");
    TValue dummy_ = 'T';
    TValue & value_ = dummy_;
    SEQAN_ASSERT_EQ(value_, TValue('T'));

    value_ = getValue(string, 0);
    SEQAN_ASSERT_EQ(value_, TValue('A'));

    value_ = 'T';
    SEQAN_ASSERT_EQ(getValue(string, 0), TValue('A'));
}

SEQAN_TYPED_TEST(StringTestCommon, GetValue)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceGetValue(str);

    typename TestFixture::TString const constStr;
    testSequenceGetValue(constStr);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of insert().
template <typename TString>
void testSequenceInsert(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    {
        // Test of inserting an empty string.
        TString string1("");
        TString string2("");
        insert(string1, 0, string2);
        SEQAN_ASSERT_EQ(string1, "");
    }
    {
        // Test of inserting an string.
        TString string1("A");
        TString string2("ACGT");
        insert(string1, 0, string2);
        SEQAN_ASSERT_EQ(string1, "ACGTA");
    }
}

template <typename TValue>
void testSequenceInsert(String<TValue, External<> > & /*Tag*/) {} 
template <typename TValue>
void testSequenceInsert(String<TValue, Block<> > & /*Tag*/) {} 
template <typename TValue>
void testSequenceInsert(String<TValue, Packed<> > & /*Tag*/) {} 

SEQAN_TYPED_TEST(StringTestCommon, Insert)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceInsert(str);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of insertValue().
template <typename TString>
void testSequenceInsertValue(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

    // Test of inserting into an empty string.
    TString string("");
    insert(string, 0, 'A');
    SEQAN_ASSERT_EQ(string, "A");

    // Test of inserting into a non empty string.
    insert(string, 0, 'C');
    SEQAN_ASSERT_EQ(string, "CA");
}

template <typename TValue>
void testSequenceInsertValue(String<TValue, External<> > & /*Tag*/) {} 
template <typename TValue>
void testSequenceInsertValue(String<TValue, Block<> > & /*Tag*/) {} 
template <typename TValue>
void testSequenceInsertValue(String<TValue, Packed<> > & /*Tag*/) {} 

SEQAN_TYPED_TEST(StringTestCommon, InsertValue)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceInsertValue(str);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of iter().
template <typename TString>
void testSequenceIter(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Iterator<TString>::Type TIterator;
    typedef typename Iterator<TString, Standard>::Type TStandardIterator;
    typedef typename Iterator<TString, Rooted>::Type TRootedIterator;

    // Test on an empty string.
    {
        TString string;
        TIterator iterator = iter(string, 0);
        TStandardIterator standardIterator = iter(string, 0);
        TRootedIterator rootedIterator = iter(string, 0);
        SEQAN_ASSERT(iterator == begin(string));
        SEQAN_ASSERT(standardIterator == begin(string));
        SEQAN_ASSERT(rootedIterator == begin(string));
    }

    // Test on a non empty string.
    {
        TString string("A");
        TIterator iterator = iter(string, 0);
        TStandardIterator standardIterator = iter(string, 0);
        TRootedIterator rootedIterator = iter(string, 0);
        SEQAN_ASSERT_EQ(getValue(iterator), getValue(string, 0));
        SEQAN_ASSERT_EQ(getValue(standardIterator), getValue(string, 0));
        SEQAN_ASSERT_EQ(getValue(rootedIterator), getValue(string, 0));
    }

    // Test on a non empty string.
    {
        TString string("ACGT");
        TIterator iterator = iter(string, 3);
        TStandardIterator standardIterator = iter(string, 3);
        TRootedIterator rootedIterator = iter(string, 3);
        SEQAN_ASSERT_EQ(getValue(iterator), getValue(string, 3));
        SEQAN_ASSERT_EQ(getValue(standardIterator), getValue(string, 3));
        SEQAN_ASSERT_EQ(getValue(rootedIterator), getValue(string, 3));
    }
}

SEQAN_TYPED_TEST(StringTestCommon, Iter)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceIter(str);

    typename TestFixture::TString constStr;
    testSequenceIter(constStr);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of length().
template <typename TString>
void testSequenceLength(TString & /*Tag*/)
{
    using namespace seqan;

    // Test on an empty string.
    TString string1;
    SEQAN_ASSERT_EQ(length(string1), 0u);

    // Test on a non empty string.
    TString string2("CGTACGTATC");
    SEQAN_ASSERT_EQ(length(string2), 10u);
}

SEQAN_TYPED_TEST(StringTestCommon, Length)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceLength(str);
 
    typename TestFixture::TString constStr;
    testSequenceLength(constStr);
 
    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of value().
template <typename TString>
void testSequenceMoveValue(TString & /*Tag*/)
{
    using namespace seqan;

    TString string("");

    resize(string, 2);
    moveValue(string, 1, 'G');
    SEQAN_ASSERT_EQ(string[1], 'G');
}

template <typename TValue>
void testSequenceMoveValue(String<TValue, Packed<> > & /*Tag*/) {}

SEQAN_TYPED_TEST(StringTestCommon, MoveValue)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceMoveValue(str);
   
    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of replace().
template <typename TString>
void testSequenceReplace(TString & /*Tag*/)
{
    using namespace seqan;

    {
        TString string1("");
        TString string2("");

        // TODO (singer): This is problematic according to the documentation.
        // 0 can be a position or an iterator causing compiler errors.
        replace(string1, 0, 0, string2);
        SEQAN_ASSERT_EQ(string1, "");
    }
    {
        TString string1("ACACACAC");
        TString string2("");

        replace(string1, 4, 4, string2);
        SEQAN_ASSERT_EQ(string1, "ACACACAC");
    }
    {
        TString string1("");
        TString string2("GTGTGTGT");

        replace(string1, 0, 0, string2);
        SEQAN_ASSERT_EQ(string1, "GTGTGTGT");
    }
    {
        TString string1("ACACACAC");
        TString string2("GTGTGTGT");

        replace(string1, 4, 4, string2);
        SEQAN_ASSERT_EQ(string1, "ACACGTGTGTGTACAC");
    }
}

template <typename TValue>
void testSequenceReplace(String<TValue, External<> > & /*Tag*/) {}
template <typename TValue>
void testSequenceReplace(String<TValue, Block<> > & /*Tag*/) {}
template <typename TValue>
void testSequenceReplace(String<TValue, Packed<> > & /*Tag*/) {}

SEQAN_TYPED_TEST(StringTestCommon, Replace)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceReplace(str);
   
    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of reserve().
template <typename TString>
void testSequenceReserve(TString & /*Tag*/)
{
    using namespace seqan;

    TString string("");

    reserve(string, 0);
    SEQAN_ASSERT_EQ(capacity(string), 0u);

    reserve(string, 1000);
    SEQAN_ASSERT_GEQ(capacity(string), 10u);

    // If the the new capacity is smaller than the current one
    // the new capacity must be larger or equal to the current length.
    reserve(string, 1);
    SEQAN_ASSERT_GEQ(capacity(string), length(string));
}

template <typename TValue>
void testSequenceReserve(String<TValue, MMap<> > & /*Tag*/) {}

SEQAN_TYPED_TEST(StringTestCommon, Reserve)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceReserve(str);
   
    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of resize().
template <typename TString>
void testSequenceResize(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

    TString string("");

    resize(string, 0);
    SEQAN_ASSERT_EQ(length(string), 0u);

    resize(string, 10);
    SEQAN_ASSERT_EQ(length(string), 10u);
    // TODO (singer): resize should initialize newly allocated memory,
    // which it does not at the moment!
    //SEQAN_ASSERT_EQ(string[0], TValue());
    //SEQAN_ASSERT_EQ(string[0], TValue());

    resize(string, 0);
    SEQAN_ASSERT_EQ(length(string), 0u);

    resize(string, 10, TValue('C'));
    SEQAN_ASSERT_EQ(string[0], TValue('C'));
    SEQAN_ASSERT_EQ(string[0], TValue('C'));
}

template <typename TValue>
void testSequenceResize(String<TValue, External<> > & /*Tag*/) {}
template <typename TValue>
void testSequenceResize(String<TValue, Block<> > & /*Tag*/) {}

SEQAN_TYPED_TEST(StringTestCommon, Resize)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceResize(str);
   
    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of swap().
template <typename TString>
void testSequenceSwap(TString & /*Tag*/)
{
    using namespace seqan;
    {
        TString string1("");
        TString string2("");
        TString string3 = string1;
        TString string4 = string2;

        swap(string1, string2);
        SEQAN_ASSERT_EQ(string1, string4);
        SEQAN_ASSERT_EQ(string2, string3);
    }
    {
        TString string1("ACGT");
        TString string2("");
        TString string3 = string1;
        TString string4 = string2;

        swap(string1, string2);
        SEQAN_ASSERT_EQ(string1, string4);
        SEQAN_ASSERT_EQ(string2, string3);
    }
    {
        TString string1("");
        TString string2("ACGT");
        TString string3 = string1;
        TString string4 = string2;

        swap(string1, string2);
        SEQAN_ASSERT_EQ(string1, string4);
        SEQAN_ASSERT_EQ(string2, string3);
    }
    {
        TString string1("ACAC");
        TString string2("GTGT");
        TString string3 = string1;
        TString string4 = string2;

        swap(string1, string2);
        SEQAN_ASSERT_EQ(string1, string4);
        SEQAN_ASSERT_EQ(string2, string3);
    }
}

template <typename TValue>
void testSequenceSwap(String<TValue, External<> > & /*Tag*/) {}
template <typename TValue>
void testSequenceSwap(String<TValue, MMap<> > & /*Tag*/) {}
template <typename TValue>
void testSequenceSwap(String<TValue, Packed<> > & /*Tag*/) {}
template <typename TValue>
void testSequenceSwap(String<TValue, Block<> > & /*Tag*/) {}
template <typename TValue>
void testSequenceSwap(String<TValue, Array<100> > & /*Tag*/) {}

SEQAN_TYPED_TEST(StringTestCommon, Swap)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceSwap(str);
   
    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of value().
template <typename TString>
void testSequenceValue(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

    // In contrast to getValue(), value() does not return a copy but a reference.
    // We test this using the variable value_.
    TString string("ACAC");
    TValue & value_ = value(string, 0);
    SEQAN_ASSERT_EQ(value_, 'A');

    value_ = 'G';
    SEQAN_ASSERT_EQ(value_, 'G');
    SEQAN_ASSERT_EQ(string, "GCAC");
}

// Test of value().
template <typename TString>
void testSequenceValue(TString const & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

    // In contrast to getValue(), value() does not return a copy but a reference.
    // We test this using the variable value_.
    TString const string("ACAC");
    TString const & value_ = value(string, 0);
    SEQAN_ASSERT_EQ(value_, 'A');
}

template <typename TValue>
void testSequenceValue(String<TValue, External<> > & /*Tag*/) {}
template <typename TValue>
void testSequenceValue(String<TValue, External<> > const & /*Tag*/) {}
template <typename TValue>
void testSequenceValue(String<TValue, MMap<> > & /*Tag*/) {}
template <typename TValue>
void testSequenceValue(String<TValue, MMap<> > const & /*Tag*/) {}
template <typename TValue>
void testSequenceValue(String<TValue, Packed<> > & /*Tag*/) {}
template <typename TValue>
void testSequenceValue(String<TValue, Packed<> > const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringTestCommon, Value)
{
    CountingChar::clear();

    typename TestFixture::TString str;
    testSequenceValue(str);
   
    typename TestFixture::TString const constStr;
    testSequenceValue(constStr);

    testConstruct(typename Value<typename TestFixture::TString>::Type());
    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}
#endif  // CORE_TESTS_SEQUENCE_TEST_SEQUENCE_H_
