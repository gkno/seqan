/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
  ===========================================================================
  Copyright (C) 2010
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.
  
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  
  ===========================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ===========================================================================
  High-level tests for the sequence_journal module.  In the future, these
  should be broken down and code should be tested more isolated.
  ===========================================================================
*/

#ifndef TEST_SEQUENCE_JOURNAL_TEST_SEQUENCE_JOURNAL_H_
#define TEST_SEQUENCE_JOURNAL_TEST_SEQUENCE_JOURNAL_H_

#include <cstdlib>
#include <sstream>
#include <string>

#include <seqan/basic.h>
#include <seqan/misc/misc_random.h>
#include <seqan/sequence.h>
#include <seqan/sequence_journal.h>

using namespace seqan;

// Test setHost(), host().
template <typename TStringJournalSpec>
void testSequenceJournalHost(TStringJournalSpec const &)
{
    CharString charStr = "test";
    SequenceJournal<CharString, TStringJournalSpec> sequenceJournal;

    setHost(sequenceJournal, charStr);
    SEQAN_ASSERT_EQ(&charStr, &host(sequenceJournal));
    SEQAN_ASSERT_EQ(charStr, host(sequenceJournal));
}


// Test clear().
template <typename TStringJournalSpec>
void testSequenceJournalClear(TStringJournalSpec const &)
{
    CharString charStr = "test";
    SequenceJournal<CharString, TStringJournalSpec> sequenceJournal(charStr);

    insert(sequenceJournal, 4, "!");
    {
        std::stringstream tmp;
        tmp << sequenceJournal;
        SEQAN_ASSERT_EQ("test!", tmp.str());
    }

    clear(sequenceJournal);
    {
        std::stringstream tmp;
        tmp << sequenceJournal;
        SEQAN_ASSERT_EQ("test", tmp.str());
    }
}


// Test erase() with position only
template <typename TStringJournalSpec>
void testSequenceJournalErasePosition(TStringJournalSpec const &)
{
    CharString charStr = "test";
    SequenceJournal<CharString, TStringJournalSpec> sequenceJournal(charStr);

    erase(sequenceJournal, 1);
    std::stringstream tmp;
    tmp << sequenceJournal;
    SEQAN_ASSERT_EQ("tst", tmp.str());
}


// Test erase() with begin/end parameters.
template <typename TStringJournalSpec>
void testSequenceJournalEraseBeginEnd(TStringJournalSpec const &)
{
    CharString charStr = "test";
    SequenceJournal<CharString, TStringJournalSpec> sequenceJournal(charStr);

    erase(sequenceJournal, 1, 3);
    std::stringstream tmp;
    tmp << sequenceJournal;
    SEQAN_ASSERT_EQ("tt", tmp.str());
}


// Test insert().
template <typename TStringJournalSpec>
void testSequenceJournalInsert(TStringJournalSpec const &)
{
    CharString charStr = "test";
    SequenceJournal<CharString, TStringJournalSpec> sequenceJournal(charStr);

    insert(sequenceJournal, 1, "!!");
    std::stringstream tmp;
    tmp << sequenceJournal;
    SEQAN_ASSERT_EQ("t!!est", tmp.str());
}


// Test insertValue().
template <typename TStringJournalSpec>
void testSequenceJournalInsertValue(TStringJournalSpec const &)
{
    CharString charStr = "test";
    SequenceJournal<CharString, TStringJournalSpec> sequenceJournal(charStr);

    insert(sequenceJournal, 2, 'X');
    std::stringstream tmp;
    tmp << sequenceJournal;
    SEQAN_ASSERT_EQ("teXst", tmp.str());
}


// Test assignValue().
template <typename TStringJournalSpec>
void testSequenceJournalAssignValue(TStringJournalSpec const &)
{
    CharString charStr = "test";
    SequenceJournal<CharString, TStringJournalSpec> sequenceJournal(charStr);

    assignValue(sequenceJournal, 2, 'X');
    std::stringstream tmp;
    tmp << sequenceJournal;
    SEQAN_ASSERT_EQ("teXt", tmp.str());
}


// Test assignInfix().
template <typename TStringJournalSpec>
void testSequenceJournalAssignInfix(TStringJournalSpec const &)
{
    CharString charStr = "test";
    SequenceJournal<CharString, TStringJournalSpec> sequenceJournal(charStr);

    assignInfix(sequenceJournal, 2, 4, "quick brown fox");
    std::stringstream tmp;
    tmp << sequenceJournal;
    SEQAN_ASSERT_EQ("tequick brown fox", tmp.str());
}


// Test length().
template <typename TStringJournalSpec>
void testSequenceJournalLength(TStringJournalSpec const &)
{
    CharString charStr = "test";
    SequenceJournal<CharString, TStringJournalSpec> sequenceJournal(charStr);

    assignInfix(sequenceJournal, 2, 3, "XX");
    SEQAN_ASSERT_EQ(length("teXXt"), length(sequenceJournal));
}


template <typename TStringJournalSpec>
void testSequenceJournalCopyConstructor(TStringJournalSpec const &)
{
    CharString charStr = "test";
    SequenceJournal<CharString, TStringJournalSpec> * sequenceJournalPtr = new SequenceJournal<CharString, TStringJournalSpec>(charStr);
    SequenceJournal<CharString, TStringJournalSpec> sequenceJournal(*sequenceJournalPtr);
    delete sequenceJournalPtr;

    insert(sequenceJournal, 1, "XX");

    std::stringstream ss;
    ss << sequenceJournal;
    SEQAN_ASSERT_EQ("tXXest", ss.str());
}


template <typename TStringJournalSpec>
void testSequenceJournalBeginEndIterator(TStringJournalSpec const &)
{
    CharString charStr = "test";
    typedef SequenceJournal<CharString, TStringJournalSpec> TSequenceJournal;
    TSequenceJournal sequenceJournal(charStr);
    insert(sequenceJournal, 2, 'X');

    typedef typename Iterator<TSequenceJournal, Standard>::Type TIterator;

    {  // Test pre-increment iteration.
        CharString buffer;
        for (TIterator it = begin(sequenceJournal), itend = end(sequenceJournal); it != itend; ++it)
            appendValue(buffer, *it);
        SEQAN_ASSERT_EQ("teXst", buffer);
    }
    {  // Test post-increment iteration.
        CharString buffer;
        for (TIterator it = begin(sequenceJournal), itend = end(sequenceJournal); it != itend; it++)
            appendValue(buffer, *it);
        SEQAN_ASSERT_EQ("teXst", buffer);
    }
}


template <typename TStringJournalSpec>
void testSequenceJournalBeginEndConstIterator(TStringJournalSpec const &)
{
    CharString charStr = "test";
    SequenceJournal<CharString, TStringJournalSpec> sequenceJournal(charStr);
    insert(sequenceJournal, 2, 'X');

    typedef typename Iterator<SequenceJournal<CharString, TStringJournalSpec> const, Standard>::Type TIterator;

    CharString buffer;
    for (TIterator it = begin(sequenceJournal, Standard()), itend = end(sequenceJournal, Standard()); it != itend; ++it)
        appendValue(buffer, *it);
    SEQAN_ASSERT_EQ("teXst", buffer);
}


// Perform a random insert/edit/delete of the sequence journal.
template <typename TStringJournalSpec>
void testSequenceJournalFuzzying(TStringJournalSpec const &)
{
    using namespace seqan;

    const unsigned INITIAL_LENGTH = 100;
    const unsigned NUM_CHANGES = 100;
    const unsigned MAX_INSERT = 100;

    const unsigned SEED = 42;
    std::srand(SEED);
    mtRandInit(false);

#define RAND_CHAR() ('A' + mtRand() % ('Z' - 'A'))

    // Build random reference and host string.
    String<char> string;
    reserve(string, INITIAL_LENGTH);
    String<char> host;
    reserve(host, INITIAL_LENGTH);
    for (unsigned i = 0; i < INITIAL_LENGTH; ++i) {
        char c = RAND_CHAR();
        appendValue(string, c);
        appendValue(host, c);
    }

    // Output of initial sequences.
//     std::cout << "reference = " << string << std::endl;
//     std::cout << "host = " << host << std::endl;

    typedef SequenceJournal<String<char>, TStringJournalSpec> TSequenceJournal;

    // Construct sequence journal on host.
    TSequenceJournal sequenceJournal(host);

//     unsigned nextId = 0;
//     std::cerr << "digraph {" << std::endl;

    // We will use a string stream to test the string result of tmp.
    {
        std::stringstream tmp;
        tmp << sequenceJournal;
//         SEQAN_ASSERT_EQ(string, tmp.str());
//         std::cout << "string = " << string << std::endl;
//         std::cout << "jrnld  = " << tmp.str() << std::endl;
//         std::cout << "  tree = " << sequenceJournal._journalTree << std::endl;
//         std::cout << "  orig = " << value(sequenceJournal._host) << std::endl;
//         std::cout << "  buff = " << sequenceJournal._insertionBuffer << std::endl;
//         journalTreeToDot(std::cerr, nextId, sequenceJournal._journalTree);
    }

    size_t expectedLength = length(string);

    for (unsigned i = 0; i < NUM_CHANGES; ++i) {
//         std::cout << "i == " << i << std::endl;
        unsigned changeType = mtRand() % 3;
        if (changeType == 0) {  // edit
            if (length(string) == 0)
                continue;
            unsigned begin = 0;
            unsigned end = 0;
            while (begin == end) {
                begin = mtRand() % length(string);
                end = mtRand() % (length(string) + 1);
            }
            if (begin > end)
                std::swap(begin, end);            
            unsigned len = end - begin;
            String<char> buffer;
            reserve(buffer, len);
            for (unsigned i = 0; i < len; ++i)
                appendValue(buffer, RAND_CHAR());
            // Perform insert.
//             std::cout << "assignInfix(sequenceJournal, " << begin << ", " << end << ", \"" << buffer << "\")" << std::endl;
//             std::cout << "assignInfix(sequenceJournal, " << begin << ", " << end << ", buffer, len(buffer) == " << length(buffer) << ")" << std::endl;
            infix(string, begin, end) = buffer;
//             std::cout << "pre assign infix " << length(sequenceJournal) << std::endl;
            assignInfix(sequenceJournal, begin, end, buffer);
//             std::cout << "post assign infix " << length(sequenceJournal) << std::endl;
        } else if (changeType == 1) {  // insert
            unsigned begin = 0;
            unsigned len = 0;
            while (len == 0) {
                if (length(string) == 0)
                    begin = 0;
                else
                    begin = mtRand() % length(string);
                len = mtRand() % MAX_INSERT + 1;
            }
            expectedLength += len;
            String<char> buffer;
            reserve(buffer, len);
            for (unsigned i = 0; i < len; ++i)
                appendValue(buffer, RAND_CHAR());
            // Perform insert.
//             std::cout << "insert(sequenceJournal, " << begin << ", \"" << buffer << "\")" << std::endl;
//             std::cout << "insert(sequenceJournal, " << begin << ", buffer)" << std::endl;
            infix(string, begin, begin) = buffer;
            insert(sequenceJournal, begin, buffer);
        } else if (changeType == 2) {  // delete
            if (length(string) == 0)
                continue;
            //std::cerr << "length(string) == " << length(string) << std::endl;
            //std::cerr << "string == " << string << std::endl;
            //std::cerr << "journal string== " << sequenceJournal << std::endl;
            //std::cerr << "journal string== " << sequenceJournal._journalTree << std::endl;
            unsigned begin = 0;
            unsigned end = 0;
            while (begin == end) {
                begin = mtRand() % length(string);
                end = mtRand() % (length(string) + 1);
            }
            if (begin > end)
                std::swap(begin, end);
            expectedLength -= (end - begin);
            // Perform erase.
//             std::stringstream tmp;
//             tmp << sequenceJournal;
//             std::cout << ",---" << std::endl;
//             std::cout << "| string = " << string << std::endl;
//             std::cout << "| jrnld  = " << tmp.str() << std::endl;
//             std::cout << "|   tree = " << sequenceJournal._journalTree << std::endl;
//             std::cout << "|   orig = " << value(sequenceJournal._host) << std::endl;
//             std::cout << "|   buff = " << sequenceJournal._insertionBuffer << std::endl;
//             std::cout << "erase(sequenceJournal, " << begin << ", " << end << ")" << std::endl;
//             std::cout << "`---" << std::endl;
//             std::cout << sequenceJournal._journalTree << std::endl;
            erase(string, begin, end);
            erase(sequenceJournal, begin, end);
//             std::cout << sequenceJournal._journalTree << std::endl;
        } else {
            SEQAN_ASSERT_FAIL("Invalid change type.");
        }

        {
            // Check via stream operator<< into stringstream.
            std::stringstream tmp;
            tmp << sequenceJournal;
            SEQAN_ASSERT_EQ(expectedLength, length(tmp.str()));
            SEQAN_ASSERT_EQ(expectedLength, length(string));
            SEQAN_ASSERT_EQ(string, tmp.str());
//             std::cout << "string = " << string << std::endl << "tmp.str() = " << tmp.str() << std::endl;
            // Check via iterator on the journal string.
            std::string buffer;
            reserve(buffer, length(sequenceJournal));
            typedef typename Iterator<TSequenceJournal, Standard>::Type TIterator;
//             std::cout << sequenceJournal._journalTree << std::endl;
            for (TIterator it = begin(sequenceJournal), itend = end(sequenceJournal, Standard()); it != itend; ++it) {
                appendValue(buffer, *it);
            }
//             std::cout << "buffer = " << buffer << std::endl;
            SEQAN_ASSERT_EQ(expectedLength, length(buffer));
            SEQAN_ASSERT_EQ(buffer, tmp.str());
        }

        {
            // Check operator+ and operator+= on sequence journal iterators;
            typedef typename Iterator<TSequenceJournal, Standard>::Type TSequenceJournalIterator;
            typedef typename Iterator<String<char>, Standard>::Type TCharStringIterator;

            TSequenceJournalIterator sjIt = begin(sequenceJournal);
            TCharStringIterator csIt = begin(string);
            size_t remaining = length(sequenceJournal);

            while (remaining > 1) {
                SEQAN_ASSERT_TRUE(csIt != end(string) - 1);
                size_t len = mtRand() % (remaining + 1);
                remaining -= len;
                if (remaining == 0)
                    break;
                SEQAN_ASSERT_EQ(*(sjIt + len), *(csIt + len));
                sjIt += len;
                csIt += len;
                SEQAN_ASSERT_EQ(*sjIt, *csIt);
            }
        }
    }

#undef RAND_CHAR
}


// Tag: Unbalanced()


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_host) {
    testSequenceJournalHost(Unbalanced());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_clear) {
    testSequenceJournalClear(Unbalanced());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_erase_position) {
    testSequenceJournalErasePosition(Unbalanced());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_erase_begin_end) {
    testSequenceJournalEraseBeginEnd(Unbalanced());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_insert) {
    testSequenceJournalInsert(Unbalanced());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_insert_value) {
    testSequenceJournalInsertValue(Unbalanced());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_assign_value) {
    testSequenceJournalAssignValue(Unbalanced());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_assign_infix) {
    testSequenceJournalAssignInfix(Unbalanced());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_length) {
    testSequenceJournalLength(Unbalanced());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_copy_constructor) {
    testSequenceJournalCopyConstructor(Unbalanced());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_begin_end_iterator) {
    testSequenceJournalBeginEndIterator(Unbalanced());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_begin_end_const_iterator) {
    testSequenceJournalBeginEndConstIterator(Unbalanced());
}

SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_fuzzying) {
    testSequenceJournalFuzzying(Unbalanced());
}


// Tag: SortedArray()

SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_host) {
    testSequenceJournalHost(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_clear) {
    testSequenceJournalClear(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_erase_position) {
    testSequenceJournalErasePosition(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_erase_begin_end) {
    testSequenceJournalEraseBeginEnd(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_insert) {
    testSequenceJournalInsert(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_insert_value) {
    testSequenceJournalInsertValue(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_assign_value) {
    testSequenceJournalAssignValue(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_assign_infix) {
    testSequenceJournalAssignInfix(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_length) {
    testSequenceJournalLength(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_copy_constructor) {
    testSequenceJournalCopyConstructor(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_begin_end_iterator) {
    testSequenceJournalBeginEndIterator(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_begin_end_const_iterator) {
    testSequenceJournalBeginEndConstIterator(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_fuzzying) {
    testSequenceJournalFuzzying(SortedArray());
}


#endif  // TEST_SEQUENCE_JOURNAL_TEST_SEQUENCE_JOURNAL_H_
