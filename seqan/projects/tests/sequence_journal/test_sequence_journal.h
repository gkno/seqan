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

SEQAN_DEFINE_TEST(test_sequence_journal) {
    using namespace seqan;

    String<char> charString = "Hello world!";
    String<char, Journal<Alloc<>, Unbalanced> > journalString(charString);
    std::cout << journalString._journalTree << std::endl;
    std::cout << journalString << std::endl;
    std::cout << "clear(journalString);" << std::endl;
    clear(journalString);
    std::cout << journalString._journalTree << std::endl;
    std::cout << journalString << std::endl;
    std::cout << "erase(journalString, 0, 6);" << std::endl;
    erase(journalString, 0, 6);
    std::cout << journalString._journalTree << std::endl;
    std::cout << journalString << std::endl;
    std::cout << "clear(journalString);" << std::endl;
    clear(journalString);
    std::cout << journalString._journalTree << std::endl;
    std::cout << journalString << std::endl;
    std::cout << "erase(journalString, 5, 12);" << std::endl;
    erase(journalString, 5, 12);
    std::cout << journalString._journalTree << std::endl;
    std::cout << journalString << std::endl;
    std::cout << "clear(journalString);" << std::endl;
    clear(journalString);
    std::cout << journalString._journalTree << std::endl;
    std::cout << journalString << std::endl;
    std::cout << "erase(journalString, 5, 6);" << std::endl;
    erase(journalString, 5, 6);
    std::cout << journalString._journalTree << std::endl;
    std::cout << journalString << std::endl;
    std::cout << "clear(journalString);" << std::endl;
    clear(journalString);
    std::cout << journalString._journalTree << std::endl;
    std::cout << journalString << std::endl;
    std::cout << "erase(journalString, 5, 6);" << std::endl;
    erase(journalString, 5, 6);
    std::cout << journalString._journalTree << std::endl;
    std::cout << journalString << std::endl;
    std::cout << "erase(journalString, 2, 3);" << std::endl;
    erase(journalString, 2, 3);
    std::cout << journalString._journalTree << std::endl;
    std::cout << journalString << std::endl;
    std::cout << "erase(journalString, 1, 9);" << std::endl;
    erase(journalString, 1, 9);
    std::cout << journalString._journalTree << std::endl;
    std::cout << journalString << std::endl;

//     SEQAN_ASSERT_EQ(journalString, CharString("Hello world!"));
//     erase(journalString, 0, 6);
//     SEQAN_ASSERT_EQ(journalString, CharString("world!"));
}

#endif  // TEST_SEQUENCE_JOURNAL_TEST_SEQUENCE_JOURNAL_H_
