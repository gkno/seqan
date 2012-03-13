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

// --------------------------------------------------------------------------
// CountingChar is used to test sequences of non simple data types.
// --------------------------------------------------------------------------

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

// --------------------------------------------------------------------------
// Generic String Tests
// --------------------------------------------------------------------------

// Test whether sequences are copy constructible.
template <typename TString>
void testSequenceCopyConstructible(TString & /*Tag*/)
{
    using namespace seqan;

    TString string1 = "ACGCTAGCAT";
    TString string2(string1);
    SEQAN_ASSERT_EQ(string1, string2);
}

// Test whether sequences are default constructible.
template <typename TString>
void testSequenceDefaultConstructible(TString & /*Tag*/)
{
    using namespace seqan;

    TString string;
    SEQAN_ASSERT_EQ(begin(string), end(string));
}

// Test of append().
template <typename TString>
void testSequenceAppend(TString & /*Tag*/)
{
    using namespace seqan;

    // Test the append function on empty strings
    {
        TString string1 = "";
        TString string2 = "";
        append(string1, string2);
        SEQAN_ASSERT_EQ(string1, "");
    }

    // Test the append function on one empty and one non empty string
    {
        TString string1 = "ACGTACGTACGT";
        TString string2 = "";
        append(string1, string2);
        SEQAN_ASSERT_EQ(string1, "ACGTACGTACGT");
    }

    // Test the append function on one empty and one non empty string
    {
        TString string1 = "";
        TString string2 = "TTGGATTAACC";
        append(string1, string2);
        SEQAN_ASSERT_EQ(string1, "TTGGATTAACC");
    }

    // Test the append function on two non empty strings.
    {
        TString string1 = "ACGTACGTACGT";
        TString string2 = "TTGGATTAACCC";
        append(string1, string2);
        SEQAN_ASSERT_EQ(string1, "ACGTACGTACGTTTGGATTAACCC");
    }
}

// Test of appendValue().
template <typename TString>
void testSequenceAppendValue(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

    TString string = "";

    // Test the appendValue function
    TValue value = 'A';
    appendValue(string, value);   
    SEQAN_ASSERT_EQ(string, "A");

    // Test the appendValue function
    appendValue(string, 'A');   
    SEQAN_ASSERT_EQ(string, "AA");
}

// Test whether sequences are assignable.
template <typename TString>
void testSequenceAssignable(TString & /*Tag*/)
{
    using namespace seqan;

    // Test on an empty string.
    TString string1 = "";

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

    // Test the basic concept on a non empty string.
    string1 = "ACGTACGTACGT";

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

// Test of assignValue().
template <typename TString>
void testSequenceAssignValue(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

    TString string = "";

    resize(string, 2);
    assignValue(string, 1, TValue('G'));
    SEQAN_ASSERT_EQ(string[1], TValue('G'));

}

// Test of begin().
template <typename TString>
void testSequenceBegin(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

    TString string = "ACGT";
    SEQAN_ASSERT_EQ(*begin(string), TValue('A'));
    SEQAN_ASSERT_EQ(*begin(string, Standard()), TValue('A'));
    SEQAN_ASSERT_EQ(*begin(string, Rooted()), TValue('A'));
}

// Test of beginPosition().
template <typename TString>
void testSequenceBeginPosition(TString & /*Tag*/)
{
    using namespace seqan;

    // Test on an empty string.
    TString string1 = "";
    SEQAN_ASSERT_EQ(beginPosition(string1), 0u);

    // Test on a non empty string.
    TString string2 = "ACGT";
    SEQAN_ASSERT_EQ(beginPosition(string2), 0u);
}

// Test of capacity().
template <typename TString>
void testSequenceCapacity(TString & /*Tag*/)
{
    using namespace seqan;

    // Test on an empty string.
    TString string1 = "";
    SEQAN_ASSERT_GEQ(capacity(string1), length(string1));

    // Test on a non empty string.
    TString string2 = "ACGTACGTACGT";
    SEQAN_ASSERT_GEQ(capacity(string2), length(string2));
}

// Test of end().
template <typename TString>
void testSequenceEnd(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

    TString string = "ACGT";
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

// Test of endPosition().
template <typename TString>
void testSequenceEndPosition(TString & /*Tag*/)
{
    using namespace seqan;

    // Test on an empty string.
    TString string1 = "";
    SEQAN_ASSERT_EQ(endPosition(string1), length(string1));

    // Test on a non empty string.
    TString string2 = "ACGT";
    SEQAN_ASSERT_EQ(endPosition(string2), length(string2));
}

// Test of getValue().
template <typename TString>
void testSequenceGetValue(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

     // In contrast to value(), getValue() does not return a reference but a copy.
     // We test this using the variable value_.
    TString string = "ACGT";
    TValue dummy_ = 'T';
    TValue & value_ = dummy_; 
    SEQAN_ASSERT_EQ(value_, TValue('T'));

    value_ = getValue(string, 0);   
    SEQAN_ASSERT_EQ(value_, TValue('A'));

    value_ = 'T';
    SEQAN_ASSERT_EQ(getValue(string, 0), TValue('A'));
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
    TString string2 = "CGTACGTATC";
    SEQAN_ASSERT_EQ(length(string2), 10u);
}

// Test of value().
template <typename TString>
void testSequenceMoveValue(TString & /*Tag*/)
{
    using namespace seqan;

    TString string = "";

    resize(string, 2);
    moveValue(string, 1, 'G');  
    SEQAN_ASSERT_EQ(string[1], 'G');
}

// Test of replace().
template <typename TString>
void testSequenceReplace(TString & /*Tag*/)
{
    using namespace seqan;

    TString string1 = "";
    TString string2 = "";

    // This is problematic according to the documentation.
       // 0 can be a position or an iterator causing compiler errors.
    replace(string1, 0, 0, string2);
    SEQAN_ASSERT_EQ(string1, "");

    string1 = "ACACACAC";
    string2 = "";

    replace(string1, 4, 4, string2);
    SEQAN_ASSERT_EQ(string1, "ACACACAC");

    string1 = "";
    string2 = "GTGTGTGT";

    replace(string1, 0, 0, string2);
    SEQAN_ASSERT_EQ(string1, "GTGTGTGT");

    string1 = "ACACACAC";
    string2 = "GTGTGTGT";

    replace(string1, 4, 4, string2);
    SEQAN_ASSERT_EQ(string1, "ACACGTGTGTGTACAC");
}

// Test of append().
template <typename TString>
void testSequenceValue(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

    // In contrast to getValue(), value() does not return a copy but a reference.
    // We test this using the variable value_.
    TString string = "ACAC";
    TValue & value_ = value(string, 0);  
    SEQAN_ASSERT_EQ(value_, 'A');

    value_ = 'G';
    SEQAN_ASSERT_EQ(value_, 'G');
    SEQAN_ASSERT_EQ(string, "GCAC");
}

// --------------------------------------------------------------------------
// Testing Alloc Strings With Simple Types
// --------------------------------------------------------------------------

// Test whether sequences are copy constructible.
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_copy_constructible)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceCopyConstructible(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceCopyConstructible(constTag);
}

// Test whether sequences are default constructible.
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_default_constructible)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceDefaultConstructible(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceDefaultConstructible(constTag);
}

// Test of append().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_append)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceAppend(tag);
}

// Test of appendValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_append_value)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag; 
    testSequenceAppendValue(tag);
}

// Test whether sequences are assignable.
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_assignable)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceAssignable(tag);
}

// Test of assignValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_assign_value)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceAssignValue(tag);
}

// Test of begin().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_begin)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag; 
    testSequenceBegin(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceBegin(constTag);
}

// Test of beginPosition().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_begin_position)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag; 
    testSequenceBeginPosition(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceBeginPosition(constTag);
}

// Test of capacity().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_capacity)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag; 
    testSequenceCapacity(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceCapacity(constTag);
}

// Test of end().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_end)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag; 
    testSequenceEnd(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceEnd(constTag);
}

// Test of endPosition().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_end_position)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag; 
    testSequenceEndPosition(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceEndPosition(constTag);
}

// Test of getValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_get_value)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag; 
    testSequenceGetValue(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceGetValue(constTag);
}

// Test of length().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_length)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag; 
    testSequenceLength(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceLength(constTag);
}

// Test of moveValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_move_value)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceAssignValue(tag);
}

// Test of replace().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_replace)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceReplace(tag);
}

// Test of value().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_value)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceValue(tag);
}

// --------------------------------------------------------------------------
// Testing Alloc Strings With Non Simple Types
// --------------------------------------------------------------------------

// Test whether sequences are copy constructible.
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_copy_constructible)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceCopyConstructible(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceCopyConstructible(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test whether sequences are default constructible.
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_default_constructible)
{
    using namespace seqan;

    CountingChar::clear();

    String<Dna, Alloc<> > tag;
    testSequenceDefaultConstructible(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceDefaultConstructible(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
}

// Test of append().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_append)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceAppend(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of appendValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_append_value)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceAppendValue(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test whether sequences are assignable.
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_assignable)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceAssignable(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of assignValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_assign_value)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceAssignValue(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of begin().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_begin)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag; 
    testSequenceBegin(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceBegin(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of beginPosition().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_begin_position)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag; 
    testSequenceBeginPosition(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceBeginPosition(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of capacity().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_capacity)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag; 
    testSequenceCapacity(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceCapacity(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of end().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_end)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag; 
    testSequenceEnd(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceEnd(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of endPosition().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_end_position)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag; 
    testSequenceEndPosition(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceEndPosition(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of getValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_get_value)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag; 
    testSequenceGetValue(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceGetValue(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test length().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_length)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag; 
    testSequenceLength(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceLength(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of moveValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_move_value)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceMoveValue(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of replace().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_replace)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceReplace(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of value().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_value)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceValue(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}
#endif  // CORE_TESTS_SEQUENCE_TEST_SEQUENCE_H_
