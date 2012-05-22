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
// Tests for the sequence_journaled module.  Note that we only test journaled
// strings here.  In a perfect world, we would have tests for each atomic
// part of the module but for the moment, this has to suffice.  Instead of
// testing each of the many corner cases for the Journal Entries data
// structures, we rely on "fuzzying", i.e. executing a random set of
// operations and hope that all corner cases occur there.
// ==========================================================================

#ifndef TEST_SEQUENCE_JOURNALED_TEST_SEQUENCE_JOURNALED_H_
#define TEST_SEQUENCE_JOURNALED_TEST_SEQUENCE_JOURNALED_H_

#include <cstdlib>
#include <sstream>
#include <string>

#include <seqan/basic.h>
#include <seqan/misc/misc_random.h>
#include <seqan/sequence.h>
#include <seqan/sequence_journaled.h>

using namespace seqan;

// Test assign(), operator=()
template <typename TStringJournalSpec>
void testJournaledStringAssign(TStringJournalSpec const &)
{
    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;
    
    CharString charStr = "test";
    TJournaledString journaledString(charStr);

    CharString charStr2 = "not a test!";

    // Test assignment operator with other string type.
    {
        journaledString = charStr2;
        journaledString = charStr2;

        std::stringstream tmp1;
        tmp1 << journaledString;
        SEQAN_ASSERT_EQ("not a test!", tmp1.str());

        std::stringstream tmp2;
        tmp2 << host(journaledString);
        SEQAN_ASSERT_EQ("test", tmp2.str());
    }

    // Test assign() with other string.
    {
        assign(journaledString, charStr2);

        std::stringstream tmp1;
        tmp1 << journaledString;
        SEQAN_ASSERT_EQ("not a test!", tmp1.str());

        std::stringstream tmp2;
        tmp2 << host(journaledString);
        SEQAN_ASSERT_EQ("test", tmp2.str());
    }

    TJournaledString journaledString2(charStr2);
    
    // Test assignment operator with same journaled string type.
    {
        journaledString = journaledString2;

        std::stringstream tmp1;
        tmp1 << journaledString;
        SEQAN_ASSERT_EQ("not a test!", tmp1.str());

        std::stringstream tmp2;
        tmp2 << host(journaledString);
        SEQAN_ASSERT_EQ("test", tmp2.str());
    }

    // Test assign() with same journaled string type.
    {
        assign(journaledString, journaledString2);

        std::stringstream tmp1;
        tmp1 << journaledString;
        SEQAN_ASSERT_EQ("not a test!", tmp1.str());

        std::stringstream tmp2;
        tmp2 << host(journaledString);
        SEQAN_ASSERT_EQ("test", tmp2.str());
    }
}

// Test set()
template <typename TStringJournalSpec>
void testJournaledStringSet(TStringJournalSpec const &)
{
    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;
    
    CharString charStr = "test";
    TJournaledString journaledString(charStr);

    CharString charStr2 = "not a test!";
    TJournaledString journaledString2(charStr2);

    // Test set() with other string.
    {
        set(journaledString, charStr2);

        std::stringstream tmp1;
        tmp1 << journaledString;
        SEQAN_ASSERT_EQ("not a test!", tmp1.str());

        std::stringstream tmp2;
        tmp2 << host(journaledString);
        SEQAN_ASSERT_EQ("test", tmp2.str());
    }
    
    // Test set() with same journaled string type.
    {
        set(journaledString, journaledString2);

        std::stringstream tmp1;
        tmp1 << journaledString;
        SEQAN_ASSERT_EQ("not a test!", tmp1.str());

        std::stringstream tmp2;
        tmp2 << host(journaledString);
        SEQAN_ASSERT_EQ("not a test!", tmp2.str());
    }
}

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
    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;
    typedef typename Size<TJournaledString>::Type TSize;
    
    CharString charStr = "test";
    TJournaledString journaledString(charStr);

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
        SEQAN_ASSERT_EQ(static_cast<TSize>(4), length(journaledString));
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


// Test operator[]().
template <typename TStringJournalSpec>
void testJournaledStringSubscriptOperator(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);

    SEQAN_ASSERT_EQ(journaledString[0], 't');
    SEQAN_ASSERT_EQ(journaledString[3], 't');

    // static_cast<Nothing>(journaledString[2]);
    journaledString[2] = 'X';

    SEQAN_ASSERT_EQ(journaledString[0], 't');
    SEQAN_ASSERT_EQ(journaledString[2], 'X');
    SEQAN_ASSERT_EQ(journaledString[3], 't');

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

        SEQAN_ASSERT_EQ(0u, virtualToHostPosition(journaledString, 0u));
        SEQAN_ASSERT_EQ(0u, virtualToHostPosition(journaledString, 1u));
        SEQAN_ASSERT_EQ(2u, virtualToHostPosition(journaledString, 2u));
        SEQAN_ASSERT_EQ(3u, virtualToHostPosition(journaledString, 3u));
        SEQAN_ASSERT_EQ(4u, virtualToHostPosition(journaledString, 4u));
        SEQAN_ASSERT_EQ(4u, virtualToHostPosition(journaledString, 5u));
    }
    {
        CharString charStr = "ABCDE";
        String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);
        insert(journaledString, 1, "XX");
        SEQAN_ASSERT_EQ(0u, virtualToHostPosition(journaledString, 0u));
        SEQAN_ASSERT_EQ(1u, virtualToHostPosition(journaledString, 1u));
        SEQAN_ASSERT_EQ(1u, virtualToHostPosition(journaledString, 2u));
        SEQAN_ASSERT_EQ(1u, virtualToHostPosition(journaledString, 3u));
        SEQAN_ASSERT_EQ(2u, virtualToHostPosition(journaledString, 4u));
        SEQAN_ASSERT_EQ(3u, virtualToHostPosition(journaledString, 5u));
        SEQAN_ASSERT_EQ(4u, virtualToHostPosition(journaledString, 6u));
        SEQAN_ASSERT_EQ(5u, virtualToHostPosition(journaledString, 7u));
    }
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
    std::stringstream hss;

    ss << journaledString;
    hss << host(journaledString);

    SEQAN_ASSERT_EQ("tXXest", ss.str());
    SEQAN_ASSERT_EQ("test", hss.str());
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

//         SEQAN_ASSERT(begin(journaledString, Standard()) == begin(constSJ, Standard()));
//     }
}


template <typename TStringJournalSpec>
void testJournaledStringSubscriptOperatorRandomized(TStringJournalSpec const &)
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
                SEQAN_ASSERT(csIt != end(string) - 1);
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

        {
            // Check operator- on sequence journal iterators;
            typedef typename Iterator<TJournaledString, Standard>::Type TJournaledStringIterator;
            typedef typename Iterator<String<char>, Standard>::Type TCharStringIterator;

            TJournaledStringIterator sjIt = begin(journaledString);
            TJournaledStringIterator sjIt2 = sjIt;
            TCharStringIterator csIt = begin(string);
            TCharStringIterator csIt2 = csIt;
            size_t remaining = length(journaledString);

            while (remaining > 1) {
                SEQAN_ASSERT(csIt != end(string) - 1);
                size_t len = mtRand() % (remaining + 1);
                remaining -= len;
                if (remaining == 0)
                    break;
                SEQAN_ASSERT_EQ(*(sjIt + len), *(csIt + len));
                sjIt += len;
                csIt += len;
                SEQAN_ASSERT_EQ(*sjIt, *csIt);

                SEQAN_ASSERT_EQ(csIt2 - csIt, sjIt2 - sjIt);
                SEQAN_ASSERT_EQ(csIt - csIt2, sjIt - sjIt2);
                csIt2 = csIt;
                sjIt2 = sjIt;
            }
        }
    }

#undef RAND_CHAR
}


template <typename TStringJournalSpec>
void testJournaledStringSegmentsReadOnly(TStringJournalSpec const &)
{
    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;
    typedef typename Prefix<TJournaledString>::Type TPrefix;
    typedef typename Suffix<TJournaledString>::Type TSuffix;
    typedef typename Infix<TJournaledString>::Type TInfix;

    CharString charStr = "test";
    TJournaledString journaledString;
    setHost(journaledString, charStr);
    insert(journaledString, 2, "XX");

    // Prefixes.
    {
        TPrefix prefix1 = prefix(journaledString, 3);
        SEQAN_ASSERT(prefix1 == CharString("teX"));
        TPrefix prefix2(journaledString, 3);
        SEQAN_ASSERT(prefix2 == CharString("teX"));
    }
    // Suffixes.
    {   
        TSuffix suffix1 = suffix(journaledString, 3);
        SEQAN_ASSERT(suffix1 == CharString("Xst"));
        TSuffix suffix2(journaledString, 3);
        SEQAN_ASSERT(suffix2 == CharString("Xst"));

    }
    // Infixes.
    {
        TInfix infix1 = infix(journaledString, 1, 5);
        SEQAN_ASSERT(infix1 == CharString("eXXs"));
        TInfix infix2(journaledString, 1, 5);
        SEQAN_ASSERT(infix2 == CharString("eXXs"));
    }
}


template <typename TStringJournalSpec>
void testJournaledStringSegmentsReadWrite(TStringJournalSpec const &)
{
    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;
    typedef typename Prefix<TJournaledString>::Type TPrefix;
    typedef typename Suffix<TJournaledString>::Type TSuffix;
    typedef typename Infix<TJournaledString>::Type TInfix;

    CharString charStr = "test";

    // Prefixes.
    {
        TJournaledString journaledString(charStr);
        insert(journaledString, 2, "XX");

        TPrefix prefix1 = prefix(journaledString, 3);
        SEQAN_ASSERT(prefix1 == CharString("teX"));
        prefix1 = "ABCD";
        SEQAN_ASSERT_EQ(journaledString, "ABCDXst");
        SEQAN_ASSERT_EQ(charStr, "test");
    }
    {
        TJournaledString journaledString(charStr);
        insert(journaledString, 2, "XX");

        TPrefix prefix2(journaledString, 3);
        SEQAN_ASSERT(prefix2 == CharString("teX"));
        prefix2 = "ABCD";
        SEQAN_ASSERT_EQ(journaledString, "ABCDXst");
        SEQAN_ASSERT_EQ(charStr, "test");
    }

    // Suffixes.
    {   
        TJournaledString journaledString(charStr);
        insert(journaledString, 2, "XX");

        TSuffix suffix1 = suffix(journaledString, 3);
        SEQAN_ASSERT(suffix1 == CharString("Xst"));
        suffix1 = "ABCD";
        SEQAN_ASSERT_EQ(journaledString, "teXABCD");
        SEQAN_ASSERT_EQ(charStr, "test");
    }
    {
        TJournaledString journaledString(charStr);
        insert(journaledString, 2, "XX");

        TSuffix suffix2(journaledString, 3);
        SEQAN_ASSERT(suffix2 == CharString("Xst"));
        suffix2 = "ABCD";
        SEQAN_ASSERT_EQ(journaledString, "teXABCD");
        SEQAN_ASSERT_EQ(charStr, "test");
    }

    // Infixes.
    {
        TJournaledString journaledString(charStr);
        insert(journaledString, 2, "XX");

        TInfix infix1 = infix(journaledString, 1, 5);
        SEQAN_ASSERT(infix1 == CharString("eXXs"));
        infix1 = "ABCD";
        SEQAN_ASSERT_EQ(journaledString, "tABCDt");
        SEQAN_ASSERT_EQ(charStr, "test");
    }
    {
        TJournaledString journaledString(charStr);
        insert(journaledString, 2, "XX");

        TInfix infix2(journaledString, 1, 5);
        SEQAN_ASSERT(infix2 == CharString("eXXs"));
        infix2 = "ABCD";
        SEQAN_ASSERT_EQ(journaledString, "tABCDt");
        SEQAN_ASSERT_EQ(charStr, "test");
    }
}


// Tag: UnbalancedTree()

SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_assign) {
    testJournaledStringAssign(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_set) {
    testJournaledStringSet(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_host) {
    testJournaledStringHost(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_clear) {
    testJournaledStringClear(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_erase_position) {
    testJournaledStringErasePosition(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_erase_begin_end) {
    testJournaledStringEraseBeginEnd(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_insert) {
    testJournaledStringInsert(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_insert_value) {
    testJournaledStringInsertValue(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_assign_value) {
    testJournaledStringAssignValue(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_subscript_operator) {
    testJournaledStringSubscriptOperator(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_assign_infix) {
    testJournaledStringAssignInfix(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_length) {
    testJournaledStringLength(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_virtual_to_host_position) {
    testJournaledStringVirtualToHostPosition(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_copy_constructor) {
    testJournaledStringCopyConstructor(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_begin_end_iterator) {
    testJournaledStringBeginEndIterator(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_begin_end_const_iterator) {
    testJournaledStringBeginEndConstIterator(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_subscript_operator_randomized) {
    testJournaledStringSubscriptOperatorRandomized(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_fuzzying) {
    testJournaledStringFuzzying(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_segments_read_only) {
    testJournaledStringSegmentsReadOnly(UnbalancedTree());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_unbalanced_tree_segments_read_write) {
    testJournaledStringSegmentsReadWrite(UnbalancedTree());
}


// Tag: SortedArray()

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_assign) {
    testJournaledStringAssign(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_set) {
    testJournaledStringSet(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_host) {
    testJournaledStringHost(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_clear) {
    testJournaledStringClear(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_erase_position) {
    testJournaledStringErasePosition(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_erase_begin_end) {
    testJournaledStringEraseBeginEnd(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_insert) {
    testJournaledStringInsert(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_insert_value) {
    testJournaledStringInsertValue(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_assign_value) {
    testJournaledStringAssignValue(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_subscript_operator) {
    testJournaledStringSubscriptOperator(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_assign_infix) {
    testJournaledStringAssignInfix(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_length) {
    testJournaledStringLength(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_virtual_to_host_position) {
    testJournaledStringVirtualToHostPosition(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_copy_constructor) {
    testJournaledStringCopyConstructor(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_begin_end_iterator) {
    testJournaledStringBeginEndIterator(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_begin_end_const_iterator) {
    testJournaledStringBeginEndConstIterator(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_subscript_operator_randomized) {
    testJournaledStringSubscriptOperatorRandomized(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_fuzzying) {
    testJournaledStringFuzzying(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_segments_read_only) {
    testJournaledStringSegmentsReadOnly(SortedArray());
}


SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_segments_read_write) {
    testJournaledStringSegmentsReadWrite(SortedArray());
}

#endif  // TEST_SEQUENCE_JOURNALED_TEST_SEQUENCE_JOURNALED_H_
