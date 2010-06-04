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
SEQAN_DEFINE_TEST(test_sequence_journal_host) {
    CharString charStr = "test";
    SequenceJournal<CharString, Unbalanced> sequenceJournal;

    setHost(sequenceJournal, charStr);
    SEQAN_ASSERT_EQ(&charStr, &host(sequenceJournal));
    SEQAN_ASSERT_EQ(charStr, host(sequenceJournal));
}


// Test clear().
SEQAN_DEFINE_TEST(test_sequence_clear) {
    CharString charStr = "test";
    SequenceJournal<CharString, Unbalanced> sequenceJournal(charStr);

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
SEQAN_DEFINE_TEST(test_sequence_erase_position) {
    CharString charStr = "test";
    SequenceJournal<CharString, Unbalanced> sequenceJournal(charStr);

    erase(sequenceJournal, 1);
    std::stringstream tmp;
    tmp << sequenceJournal;
    SEQAN_ASSERT_EQ("tst", tmp.str());
}


// Test erase() with begin/end parameters.
SEQAN_DEFINE_TEST(test_sequence_erase_begin_end) {
    CharString charStr = "test";
    SequenceJournal<CharString, Unbalanced> sequenceJournal(charStr);

    erase(sequenceJournal, 1, 3);
    std::stringstream tmp;
    tmp << sequenceJournal;
    SEQAN_ASSERT_EQ("tt", tmp.str());
}


// Test insert().
SEQAN_DEFINE_TEST(test_sequence_insert) {
    CharString charStr = "test";
    SequenceJournal<CharString, Unbalanced> sequenceJournal(charStr);

    insert(sequenceJournal, 1, "!!");
    std::stringstream tmp;
    tmp << sequenceJournal;
    SEQAN_ASSERT_EQ("t!!est", tmp.str());
}


// Test insertValue().
SEQAN_DEFINE_TEST(test_sequence_insert_value) {
    CharString charStr = "test";
    SequenceJournal<CharString, Unbalanced> sequenceJournal(charStr);

    insert(sequenceJournal, 2, 'X');
    std::stringstream tmp;
    tmp << sequenceJournal;
    SEQAN_ASSERT_EQ("teXst", tmp.str());
}


// Test assignValue().
SEQAN_DEFINE_TEST(test_sequence_assign_value) {
    CharString charStr = "test";
    SequenceJournal<CharString, Unbalanced> sequenceJournal(charStr);

    assignValue(sequenceJournal, 2, 'X');
    std::stringstream tmp;
    tmp << sequenceJournal;
    SEQAN_ASSERT_EQ("teXt", tmp.str());
}


// Test assignInfix().
SEQAN_DEFINE_TEST(test_sequence_assign_infix) {
    CharString charStr = "test";
    SequenceJournal<CharString, Unbalanced> sequenceJournal(charStr);

    assignInfix(sequenceJournal, 2, 4, "quick brown fox");
    std::stringstream tmp;
    tmp << sequenceJournal;
    SEQAN_ASSERT_EQ("tequick brown fox", tmp.str());
}


// Test length().
SEQAN_DEFINE_TEST(test_sequence_length) {
    CharString charStr = "test";
    SequenceJournal<CharString, Unbalanced> sequenceJournal(charStr);

    assignInfix(sequenceJournal, 2, 3, "XX");
    SEQAN_ASSERT_EQ(length("teXXt"), length(sequenceJournal));
}


SEQAN_DEFINE_TEST(test_sequence_copy_constructor)
{
    CharString charStr = "test";
    SequenceJournal<CharString, Unbalanced> * sequenceJournalPtr = new SequenceJournal<CharString, Unbalanced>(charStr);
    SequenceJournal<CharString, Unbalanced> sequenceJournal(*sequenceJournalPtr);
    delete sequenceJournalPtr;

    insert(sequenceJournal, 1, "XX");

    std::stringstream ss;
    ss << sequenceJournal;
    SEQAN_ASSERT_EQ("tXXest", ss.str());
}


SEQAN_DEFINE_TEST(test_sequence_begin_end_iterator)
{
    CharString charStr = "test";
    SequenceJournal<CharString, Unbalanced> sequenceJournal(charStr);
    insert(sequenceJournal, 2, 'X');

    typedef Iterator<SequenceJournal<CharString, Unbalanced>, Standard>::Type TIterator;

    {  // Test pre-increment iteration.
        CharString buffer;
        for (TIterator it = begin(sequenceJournal, Standard()); it != end(sequenceJournal, Standard()); ++it)
            appendValue(buffer, *it);
        SEQAN_ASSERT_EQ("teXst", buffer);
    }
    {  // Test post-increment iteration.
        CharString buffer;
        for (TIterator it = begin(sequenceJournal, Standard()); it != end(sequenceJournal, Standard()); it++)
            appendValue(buffer, *it);
        SEQAN_ASSERT_EQ("teXst", buffer);
    }
}


// TODO(holtgrew): Broken...
/*
SEQAN_DEFINE_TEST(test_sequence_begin_end_const_iterator)
{
    CharString charStr = "test";
    SequenceJournal<CharString, Unbalanced> sequenceJournal(charStr);
    insert(sequenceJournal, 2, 'X');

    typedef Iterator<SequenceJournal<CharString, Unbalanced> const, Standard>::Type TIterator;

    CharString buffer;
    for (TIterator it = begin(sequenceJournal, Standard()); it != end(sequenceJournal, Standard()); ++it)
        appendValue(buffer, *it);
    SEQAN_ASSERT_EQ("teXst", buffer);
}
*/


// Perform a random insert/edit/delete of the sequence journal.
SEQAN_DEFINE_TEST(test_sequence_journal_fuzzying) {
    using namespace seqan;

    const unsigned INITIAL_LENGTH = 10;
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

    typedef SequenceJournal<String<char>, Unbalanced> TSequenceJournal;

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
            typedef Iterator<TSequenceJournal, Standard>::Type TIterator;
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
            typedef Iterator<TSequenceJournal, Standard>::Type TSequenceJournalIterator;
            typedef Iterator<String<char>, Standard>::Type TCharStringIterator;

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

#endif  // TEST_SEQUENCE_JOURNAL_TEST_SEQUENCE_JOURNAL_H_
