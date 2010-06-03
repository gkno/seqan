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

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/sequence_journal.h>

// A simple demo of what you can do with journal strings.
SEQAN_DEFINE_TEST(test_sequence_journal_simple_demo) {
    using namespace seqan;

    String<char> charString = "1234567890ABCDEFGHIJKLMNOP";
    String<char, Journal<Alloc<>, Unbalanced> > journalString(charString);
    typedef String<char, Journal<Alloc<>, Unbalanced> > TJournalString;
    clear(journalString);

    erase(journalString, 8, 9);
    erase(journalString, 3, 4);
    erase(journalString, 13, 14);
    std::cout << journalString._journalTree << std::endl;
    typedef JournalTreeIterator<TJournalString::TJournalTree> TIterator;
    /*
    it == begin(journalString._journalTree, STandard);
    it != begin(journalString._journalTree, Standard());
    it == end(journalString._journalTree, Standard());
    it != end(journalString._journalTree, Standard());
    */
    for (TIterator it = begin(journalString._journalTree, Standard()); it != end(journalString._journalTree, Standard()); ++it) {
        std::cout << ">>> " << **it << std::endl;
        for (unsigned i = 0; i < length(_visitedNodes(it)); ++i)
            std::cout << "  " << _visitedNodes(it)[i]->virtualPosition << " ";
        std::cout << std::endl;
    }

//     String<char> charString = "Hello world!";
//     String<char, Journal<Alloc<>, Unbalanced> > journalString(charString);

//     clear(journalString);

//     std::cout << journalString << std::endl;
//     std::cout << "erase(journalString, 0, 6);" << std::endl;
//     erase(journalString, 0, 6);
//     std::cout << journalString << std::endl;
//     std::cout << "clear(journalString);" << std::endl;
//     clear(journalString);
//     std::cout << journalString << std::endl;
//     std::cout << "erase(journalString, 5, 12);" << std::endl;
//     erase(journalString, 5, 12);
//     std::cout << journalString << std::endl;
//     std::cout << "clear(journalString);" << std::endl;
//     clear(journalString);
//     std::cout << journalString << std::endl;
//     std::cout << "erase(journalString, 5, 6);" << std::endl;
//     erase(journalString, 5, 6);
//     std::cout << journalString << std::endl;
//     std::cout << "clear(journalString);" << std::endl;
//     clear(journalString);
//     std::cout << journalString << std::endl;
//     std::cout << "erase(journalString, 5, 6);" << std::endl;
//     erase(journalString, 5, 6);
//     std::cout << journalString << std::endl;
//     std::cout << "erase(journalString, 2, 3);" << std::endl;
//     erase(journalString, 2, 3);
//     std::cout << journalString << std::endl;
//     std::cout << "erase(journalString, 1, 9);" << std::endl;
//     erase(journalString, 1, 9);
//     std::cout << journalString << std::endl;
//     clear(journalString);
//     std::cout << journalString << std::endl;
//     std::cout << "insertValue(journalString, 0, '>');" << std::endl;
//     insertValue(journalString, 0, '>');
//     std::cout << journalString << std::endl;
//     std::cout << "insertValue(journalString, 7, ' ');" << std::endl;
//     insertValue(journalString, 7, ' ');
//     std::cout << journalString << std::endl;
//     std::cout << "insertValue(journalString, 7, 'y');" << std::endl;
//     insertValue(journalString, 7, 'y');
//     std::cout << journalString << std::endl;
//     std::cout << "insertValue(journalString, 7, 'm');" << std::endl;
//     insertValue(journalString, 7, 'm');
//     std::cout << journalString << std::endl;
//     std::cout << "insertValue(journalString, 16, '<');" << std::endl;
//     insertValue(journalString, 16, '<');
//     std::cout << journalString << std::endl;
//     std::cout << "assignValue(journalString, 0, '=');" << std::endl;
//     assignValue(journalString, 0, '=');
//     std::cout << journalString << std::endl;
//     std::cout << "assignInfix(journalString, 0, 4, CharString(\"12345678\"));" << std::endl;
//     assignInfix(journalString, 0, 4, CharString("12345678"));
//     std::cout << journalString << std::endl;
//     std::cout << "front = " << front(journalString) << ", back = " << back(journalString) << std::endl;


//     SEQAN_ASSERT_EQ(journalString, CharString("Hello world!"));
//     erase(journalString, 0, 6);
//     SEQAN_ASSERT_EQ(journalString, CharString("world!"));
}

#endif  // TEST_SEQUENCE_JOURNAL_TEST_SEQUENCE_JOURNAL_H_
