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
void testJournaledStringHost(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString;

    setHost(journaledString, charStr);
    SEQAN_ASSERT_EQ(&charStr, &host(journaledString));
    SEQAN_ASSERT_EQ(charStr, host(journaledString));
}


// Test clear().
template <typename TStringJournalSpec>
void testJournaledStringClear(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);

    insert(journaledString, 4, "!");
    {
        std::stringstream tmp;
        tmp << journaledString;
        SEQAN_ASSERT_EQ("test!", tmp.str());
    }

    clear(journaledString);
    {
        std::stringstream tmp;
        tmp << journaledString;
        SEQAN_ASSERT_EQ("test", tmp.str());
    }
}


// Test erase() with position only
template <typename TStringJournalSpec>
void testJournaledStringErasePosition(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);

    erase(journaledString, 1);
    std::stringstream tmp;
    tmp << journaledString;
    SEQAN_ASSERT_EQ("tst", tmp.str());
}


// Test erase() with begin/end parameters.
template <typename TStringJournalSpec>
void testJournaledStringEraseBeginEnd(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);

    erase(journaledString, 1, 3);
    std::stringstream tmp;
    tmp << journaledString;
    SEQAN_ASSERT_EQ("tt", tmp.str());
}


// Test insert().
template <typename TStringJournalSpec>
void testJournaledStringInsert(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);

    insert(journaledString, 1, "!!");
    std::stringstream tmp;
    tmp << journaledString;
    SEQAN_ASSERT_EQ("t!!est", tmp.str());
}


// Test insertValue().
template <typename TStringJournalSpec>
void testJournaledStringInsertValue(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);

    insert(journaledString, 2, 'X');
    std::stringstream tmp;
    tmp << journaledString;
    SEQAN_ASSERT_EQ("teXst", tmp.str());
}


// Test assignValue().
template <typename TStringJournalSpec>
void testJournaledStringAssignValue(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);

    assignValue(journaledString, 2, 'X');
    std::stringstream tmp;
    tmp << journaledString;
    SEQAN_ASSERT_EQ("teXt", tmp.str());
}


// Test assignInfix().
template <typename TStringJournalSpec>
void testJournaledStringAssignInfix(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);

    assignInfix(journaledString, 2, 4, "quick brown fox");
    std::stringstream tmp;
    tmp << journaledString;
    SEQAN_ASSERT_EQ("tequick brown fox", tmp.str());
}


// Test length().
template <typename TStringJournalSpec>
void testJournaledStringLength(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);

    assignInfix(journaledString, 2, 3, "XX");
    SEQAN_ASSERT_EQ(length("teXXt"), length(journaledString));
}


// Test conversion of virtual to host position.
template <typename TStringJournalSpec>
void testJournaledStringVirtualToHostPosition(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);

    insertValue(journaledString, 0, '!');
    erase(journaledString, 2);
    insertValue(journaledString, 4, '!');

    // journaledString == !tst!
    //
    //                           01234
    // alignment with original  "test"
    //                         "!t-st!"
    //                          01 2345

    SEQAN_ASSERT_EQ(0u, virtualToHostPosition(journaledString, 0));
    SEQAN_ASSERT_EQ(0u, virtualToHostPosition(journaledString, 1));
    SEQAN_ASSERT_EQ(2u, virtualToHostPosition(journaledString, 2));
    SEQAN_ASSERT_EQ(3u, virtualToHostPosition(journaledString, 3));
    SEQAN_ASSERT_EQ(4u, virtualToHostPosition(journaledString, 4));
    SEQAN_ASSERT_EQ(4u, virtualToHostPosition(journaledString, 5));
}


template <typename TStringJournalSpec>
void testJournaledStringCopyConstructor(TStringJournalSpec const &)
{
    CharString charStr = "test";
    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;
    TJournaledString * journaledStringPtr = new TJournaledString(charStr);
    TJournaledString journaledString(*journaledStringPtr);
    delete journaledStringPtr;

    insert(journaledString, 1, "XX");

    std::stringstream ss;
    ss << journaledString;
    SEQAN_ASSERT_EQ("tXXest", ss.str());
}


template <typename TStringJournalSpec>
void testJournaledStringBeginEndIterator(TStringJournalSpec const &)
{
    CharString charStr = "test";
    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;
    TJournaledString journaledString(charStr);
    insert(journaledString, 2, 'X');

    typedef typename Iterator<TJournaledString, Standard>::Type TIterator;

    {  // Test pre-increment iteration.
        CharString buffer;
        for (TIterator it = begin(journaledString), itend = end(journaledString); it != itend; ++it)
            appendValue(buffer, *it);
        SEQAN_ASSERT_EQ("teXst", buffer);
    }
    {  // Test post-increment iteration.
        CharString buffer;
        for (TIterator it = begin(journaledString), itend = end(journaledString); it != itend; it++)
            appendValue(buffer, *it);
        SEQAN_ASSERT_EQ("teXst", buffer);
    }
}


template <typename TStringJournalSpec>
void testJournaledStringBeginEndConstIterator(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);
    insert(journaledString, 2, 'X');

    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;

    // Test with non-const iterator.
    {
        typedef typename Iterator<TJournaledString, Standard>::Type TIterator;
        CharString buffer;
        for (TIterator it = begin(journaledString, Standard()), itend = end(journaledString, Standard()); it != itend; ++it)
            appendValue(buffer, *it);
        SEQAN_ASSERT_EQ("teXst", buffer);
    }
    // Test with const iterator.
    {
        typedef typename Iterator<TJournaledString const, Standard>::Type TIterator;
        String<char, Journaled<Alloc<void>, TStringJournalSpec> > const & constSJ = journaledString;
        CharString buffer;
        for (TIterator it = begin(constSJ, Standard()), itend = end(constSJ, Standard()); it != itend; ++it)
            appendValue(buffer, *it);
        SEQAN_ASSERT_EQ("teXst", buffer);
    }
    // Test comparison of const and non-const iterators.
//     {
//         typedef typename Iterator<String<char, Journaled<Alloc<void>, TStringJournalSpec>, Standard>::Type TNonConstIterator;
//         typedef typename Iterator<String<char, Journaled<Alloc<void>, TStringJournalSpec> const, Standard>::Type TIterator;
//         String<char, Journaled<Alloc<void>, TStringJournalSpec> const & constSJ = journaledString;

//         SEQAN_ASSERT_TRUE(begin(journaledString, Standard()) == begin(constSJ, Standard()));
//     }
}


template <typename TStringJournalSpec>
void testJournaledStringSubscriptOperator(TStringJournalSpec const &)
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
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(host);
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
        assignInfix(journaledString, begin, end, buffer);
    }

#undef RAND_CHAR

    SEQAN_ASSERT_EQ(length(journaledString), length(string));

    std::stringstream tmp;
    tmp << journaledString;
    CharString str2(tmp.str());
    SEQAN_ASSERT_EQ(str2, string);

    // Now, test the subscript operator.
    for (unsigned i = 0; i < length(journaledString); ++i) {
        SEQAN_ASSERT_EQ_MSG(journaledString[i], string[i], "i = %d", i);
    }
}


// Perform a random insert/edit/delete of the sequence journal.
template <typename TStringJournalSpec>
void testJournaledStringFuzzying(TStringJournalSpec const &)
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

    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;

    // Construct sequence journal on host.
    TJournaledString journaledString(host);

//     unsigned nextId = 0;
//     std::cerr << "digraph {" << std::endl;

    // We will use a string stream to test the string result of tmp.
    {
        std::stringstream tmp;
        tmp << journaledString;
//         SEQAN_ASSERT_EQ(string, tmp.str());
//         std::cout << "string = " << string << std::endl;
//         std::cout << "jrnld  = " << tmp.str() << std::endl;
//         std::cout << "  tree = " << journaledString._journalEntries << std::endl;
//         std::cout << "  orig = " << value(journaledString._host) << std::endl;
//         std::cout << "  buff = " << journaledString._insertionBuffer << std::endl;
//         journalTreeToDot(std::cerr, nextId, journaledString._journalEntries);
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
//             std::cout << "assignInfix(journaledString, " << begin << ", " << end << ", \"" << buffer << "\")" << std::endl;
//             std::cout << "assignInfix(journaledString, " << begin << ", " << end << ", buffer, len(buffer) == " << length(buffer) << ")" << std::endl;
            infix(string, begin, end) = buffer;
//             std::cout << "pre assign infix " << length(journaledString) << std::endl;
            assignInfix(journaledString, begin, end, buffer);
//             std::cout << "post assign infix " << length(journaledString) << std::endl;
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
//             std::cout << "insert(journaledString, " << begin << ", \"" << buffer << "\")" << std::endl;
//             std::cout << "insert(journaledString, " << begin << ", buffer)" << std::endl;
            infix(string, begin, begin) = buffer;
            insert(journaledString, begin, buffer);
        } else if (changeType == 2) {  // delete
            if (length(string) == 0)
                continue;
            //std::cerr << "length(string) == " << length(string) << std::endl;
            //std::cerr << "string == " << string << std::endl;
            //std::cerr << "journal string== " << journaledString << std::endl;
            //std::cerr << "journal string== " << journaledString._journalEntries << std::endl;
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
//             tmp << journaledString;
//             std::cout << ",---" << std::endl;
//             std::cout << "| string = " << string << std::endl;
//             std::cout << "| jrnld  = " << tmp.str() << std::endl;
//             std::cout << "|   tree = " << journaledString._journalEntries << std::endl;
//             std::cout << "|   orig = " << value(journaledString._host) << std::endl;
//             std::cout << "|   buff = " << journaledString._insertionBuffer << std::endl;
//             std::cout << "erase(journaledString, " << begin << ", " << end << ")" << std::endl;
//             std::cout << "`---" << std::endl;
//             std::cout << journaledString._journalEntries << std::endl;
            erase(string, begin, end);
            erase(journaledString, begin, end);
//             std::cout << journaledString._journalEntries << std::endl;
        } else {
            SEQAN_ASSERT_FAIL("Invalid change type.");
        }

        {
            // Check via stream operator<< into stringstream.
            std::stringstream tmp;
            tmp << journaledString;
            SEQAN_ASSERT_EQ(expectedLength, length(tmp.str()));
            SEQAN_ASSERT_EQ(expectedLength, length(string));
            SEQAN_ASSERT_EQ(string, tmp.str());
//             std::cout << "string = " << string << std::endl << "tmp.str() = " << tmp.str() << std::endl;
            // Check via iterator on the journal string.
            std::string buffer;
            reserve(buffer, length(journaledString));
            typedef typename Iterator<TJournaledString, Standard>::Type TIterator;
//             std::cout << journaledString._journalEntries << std::endl;
            for (TIterator it = begin(journaledString), itend = end(journaledString, Standard()); it != itend; ++it) {
                appendValue(buffer, *it);
            }
//             std::cout << "buffer = " << buffer << std::endl;
            SEQAN_ASSERT_EQ(expectedLength, length(buffer));
            SEQAN_ASSERT_EQ(buffer, tmp.str());
        }

        {
            // Check operator+ and operator+= on sequence journal iterators;
            typedef typename Iterator<TJournaledString, Standard>::Type TJournaledStringIterator;
            typedef typename Iterator<String<char>, Standard>::Type TCharStringIterator;

            TJournaledStringIterator sjIt = begin(journaledString);
            TCharStringIterator csIt = begin(string);
            size_t remaining = length(journaledString);

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
    testJournaledStringHost(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_clear) {
    testJournaledStringClear(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_erase_position) {
    testJournaledStringErasePosition(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_erase_begin_end) {
    testJournaledStringEraseBeginEnd(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_insert) {
    testJournaledStringInsert(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_insert_value) {
    testJournaledStringInsertValue(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_assign_value) {
    testJournaledStringAssignValue(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_assign_infix) {
    testJournaledStringAssignInfix(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_length) {
    testJournaledStringLength(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_virtual_to_host_position) {
    testJournaledStringVirtualToHostPosition(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_copy_constructor) {
    testJournaledStringCopyConstructor(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_begin_end_iterator) {
    testJournaledStringBeginEndIterator(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_begin_end_const_iterator) {
    testJournaledStringBeginEndConstIterator(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_subscript_operator) {
    testJournaledStringSubscriptOperator(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journal_unbalanced_tree_fuzzying) {
    testJournaledStringFuzzying(UnbalancedTree());
}


// Tag: SortedArray()

SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_host) {
    testJournaledStringHost(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_clear) {
    testJournaledStringClear(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_erase_position) {
    testJournaledStringErasePosition(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_erase_begin_end) {
    testJournaledStringEraseBeginEnd(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_insert) {
    testJournaledStringInsert(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_insert_value) {
    testJournaledStringInsertValue(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_assign_value) {
    testJournaledStringAssignValue(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_assign_infix) {
    testJournaledStringAssignInfix(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_length) {
    testJournaledStringLength(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_virtual_to_host_position) {
    testJournaledStringVirtualToHostPosition(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_copy_constructor) {
    testJournaledStringCopyConstructor(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_begin_end_iterator) {
    testJournaledStringBeginEndIterator(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_begin_end_const_iterator) {
    testJournaledStringBeginEndConstIterator(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_subscript_operator) {
    testJournaledStringSubscriptOperator(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journal_sorted_array_fuzzying) {
    testJournaledStringFuzzying(SortedArray());
}


#endif  // TEST_SEQUENCE_JOURNAL_TEST_SEQUENCE_JOURNAL_H_
