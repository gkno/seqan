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
  Tests for the sequence_journal module.  Note that we only test journaled
  strings here.  In a perfect world, we would have tests for each atomic
  part of the module but for the moment, this has to suffice.  Instead of
  testing each of the many corner cases for the Journal Entries data
  structures, we rely on "fuzzying", i.e. executing a random set of operations
  and hope that all corner cases occur there.
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


// Test conversion of virtual to host position.
template <typename TStringJournalSpec>
void testSequenceJournalVirtualToHostPosition(TStringJournalSpec const &)
{
    CharString charStr = "test";
    SequenceJournal<CharString, TStringJournalSpec> sequenceJournal(charStr);

    insertValue(sequenceJournal, 0, '!');
    erase(sequenceJournal, 2);
    insertValue(sequenceJournal, 4, '!');

    // sequenceJournal == !tst!
    //
    //                           01234
    // alignment with original  "test"
    //                         "!t-st!"
    //                          01 2345

    SEQAN_ASSERT_EQ(0u, virtualToHostPosition(sequenceJournal, 0));
    SEQAN_ASSERT_EQ(0u, virtualToHostPosition(sequenceJournal, 1));
    SEQAN_ASSERT_EQ(2u, virtualToHostPosition(sequenceJournal, 2));
    SEQAN_ASSERT_EQ(3u, virtualToHostPosition(sequenceJournal, 3));
    SEQAN_ASSERT_EQ(4u, virtualToHostPosition(sequenceJournal, 4));
    SEQAN_ASSERT_EQ(4u, virtualToHostPosition(sequenceJournal, 5));
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

    // Test with non-const iterator.
    {
        typedef typename Iterator<SequenceJournal<CharString, TStringJournalSpec>, Standard>::Type TIterator;
        CharString buffer;
        for (TIterator it = begin(sequenceJournal, Standard()), itend = end(sequenceJournal, Standard()); it != itend; ++it)
            appendValue(buffer, *it);
        SEQAN_ASSERT_EQ("teXst", buffer);
    }
    // Test with const iterator.
    {
        typedef typename Iterator<SequenceJournal<CharString, TStringJournalSpec> const, Standard>::Type TIterator;
        SequenceJournal<CharString, TStringJournalSpec> const & constSJ = sequenceJournal;
        CharString buffer;
        for (TIterator it = begin(constSJ, Standard()), itend = end(constSJ, Standard()); it != itend; ++it)
            appendValue(buffer, *it);
        SEQAN_ASSERT_EQ("teXst", buffer);
    }
    // Test comparison of const and non-const iterators.
//     {
//         typedef typename Iterator<SequenceJournal<CharString, TStringJournalSpec>, Standard>::Type TNonConstIterator;
//         typedef typename Iterator<SequenceJournal<CharString, TStringJournalSpec> const, Standard>::Type TIterator;
//         SequenceJournal<CharString, TStringJournalSpec> const & constSJ = sequenceJournal;

//         SEQAN_ASSERT_TRUE(begin(sequenceJournal, Standard()) == begin(constSJ, Standard()));
//     }
}


template <typename TStringJournalSpec>
void testSequenceJournalSubscriptOperator(TStringJournalSpec const &)
{
    using namespace seqan;

    const unsigned ASSIGN_COUNT = 2;
    const unsigned LENGTH = 1000;
    const unsigned SEED = 42;
    std::srand(SEED);
    mtRandInit(false);

#define RAND_CHAR() ('A' + mtRand() % ('Z' - 'A'))

    // Build random reference and host string.
    String<char> string;
    reserve(string, LENGTH);
    String<char> host;
    reserve(host, LENGTH);
    for (unsigned i = 0; i < LENGTH; ++i) {
        char c = RAND_CHAR();
        appendValue(string, c);
        appendValue(host, c);
    }

    // Create journal string over host and randomly assign infixes in
    // both to fill the tree.
    SequenceJournal<String<char>, TStringJournalSpec> sequenceJournal(host);
    for (unsigned i = 0; i < ASSIGN_COUNT; ++i) {
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
        infix(string, begin, end) = buffer;
        assignInfix(sequenceJournal, begin, end, buffer);
    }

#undef RAND_CHAR

    SEQAN_ASSERT_EQ(length(sequenceJournal), length(string));

    std::stringstream tmp;
    tmp << sequenceJournal;
    CharString str2(tmp.str());
    SEQAN_ASSERT_EQ(str2, string);

    // Now, test the subscript operator.
    for (unsigned i = 0; i < length(sequenceJournal); ++i) {
        SEQAN_ASSERT_EQ_MSG(sequenceJournal[i], string[i], "i = %d", i);
    }
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
//         std::cout << "  tree = " << sequenceJournal._journalEntries << std::endl;
//         std::cout << "  orig = " << value(sequenceJournal._host) << std::endl;
//         std::cout << "  buff = " << sequenceJournal._insertionBuffer << std::endl;
//         journalTreeToDot(std::cerr, nextId, sequenceJournal._journalEntries);
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
            //std::cerr << "journal string== " << sequenceJournal._journalEntries << std::endl;
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
//             std::cout << "|   tree = " << sequenceJournal._journalEntries << std::endl;
//             std::cout << "|   orig = " << value(sequenceJournal._host) << std::endl;
//             std::cout << "|   buff = " << sequenceJournal._insertionBuffer << std::endl;
//             std::cout << "erase(sequenceJournal, " << begin << ", " << end << ")" << std::endl;
//             std::cout << "`---" << std::endl;
//             std::cout << sequenceJournal._journalEntries << std::endl;
            erase(string, begin, end);
            erase(sequenceJournal, begin, end);
//             std::cout << sequenceJournal._journalEntries << std::endl;
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
//             std::cout << sequenceJournal._journalEntries << std::endl;
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


// Tag: UnbalancedTree()


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_host) {
    testSequenceJournalHost(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_clear) {
    testSequenceJournalClear(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_erase_position) {
    testSequenceJournalErasePosition(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_erase_begin_end) {
    testSequenceJournalEraseBeginEnd(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_insert) {
    testSequenceJournalInsert(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_insert_value) {
    testSequenceJournalInsertValue(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_assign_value) {
    testSequenceJournalAssignValue(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_assign_infix) {
    testSequenceJournalAssignInfix(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_length) {
    testSequenceJournalLength(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_virtual_to_host_position) {
    testSequenceJournalVirtualToHostPosition(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_copy_constructor) {
    testSequenceJournalCopyConstructor(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_begin_end_iterator) {
    testSequenceJournalBeginEndIterator(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_begin_end_const_iterator) {
    testSequenceJournalBeginEndConstIterator(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_subscript_operator) {
    testSequenceJournalSubscriptOperator(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_fuzzying) {
    testSequenceJournalFuzzying(UnbalancedTree());
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


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_virtual_to_host_position) {
    testSequenceJournalVirtualToHostPosition(SortedArray());
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


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_subscript_operator) {
    testSequenceJournalSubscriptOperator(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_fuzzying) {
    testSequenceJournalFuzzying(SortedArray());
}


#endif  // TEST_SEQUENCE_JOURNAL_TEST_SEQUENCE_JOURNAL_H_
