/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
  ===========================================================================
  Copyright (C) 2007
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.
  
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  
  ===========================================================================
  Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de> 
  ===========================================================================
  Tests for the SeqAn module find.
  ===========================================================================*/

#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>
#include <time.h>

#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/basic/basic_testing.h>
#include <seqan/find.h>

using namespace std;
using namespace seqan;


template <typename TAlgorithmSpec>
void Test_OnlineAlg() {
    String<unsigned int> pos;

    //____________________________________________________________________________
    // Test1 - small needle

    String<char> haystack("Dies ist ein Haystack. Ja, das ist wirklich einer!");
    Finder<String<char> > finder(haystack);

    String<char> needle("ist");
    Pattern<String<char>, TAlgorithmSpec> pattern(needle);

    while (find(finder, pattern)) {
        append(pos,position(finder));
        SEQAN_ASSERT_EQ(position(finder), beginPosition(finder));
        SEQAN_ASSERT_EQ(endPosition(finder), beginPosition(finder) + length(finder));
        SEQAN_ASSERT_EQ(length(finder), length(needle));
        SEQAN_ASSERT_EQ(begin(finder), begin(haystack) + beginPosition(finder));
        SEQAN_ASSERT_EQ(end(finder), begin(haystack) + endPosition(finder));
        SEQAN_ASSERT_EQ(infix(finder), needle);
    }

    SEQAN_ASSERT_EQ(host(pattern), needle);
    SEQAN_ASSERT_EQ(host(reinterpret_cast<Pattern<String<char>, TAlgorithmSpec> const &>(pattern)), needle);
    SEQAN_ASSERT_EQ(pos[0], 5);
    SEQAN_ASSERT_EQ(pos[1], 31);
    SEQAN_ASSERT_EQ(length(pos), 2);

    //____________________________________________________________________________
    // Test2 - large needle

    haystack = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefgaabcdef";
    setHost(finder, haystack);
    clear(finder);

    needle = "abcdefghijklmnopqrstuvwxyzabcdefg";
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern)) {
        append(pos,position(finder));
        SEQAN_ASSERT_EQ(position(finder), beginPosition(finder));
        SEQAN_ASSERT_EQ(endPosition(finder), beginPosition(finder) + length(finder));
        SEQAN_ASSERT_EQ(length(finder), length(needle));
        SEQAN_ASSERT_EQ(begin(finder), begin(haystack) + beginPosition(finder));
        SEQAN_ASSERT_EQ(end(finder), begin(haystack) + endPosition(finder));
        SEQAN_ASSERT_EQ(infix(finder), needle);
    }

    SEQAN_ASSERT_EQ(pos[0], 0);
    SEQAN_ASSERT_EQ(pos[1], 26);
    SEQAN_ASSERT_EQ(length(pos), 2);

    //____________________________________________________________________________
    // Test3 - different alphabet, small needle

    String<Dna> hstk = "aaaaaaacaa";
    Finder<String<Dna> > finderDna(hstk);

    String<Dna> ndl = "aa";
    setHost(pattern, ndl);

    clear(pos);
    while (find(finderDna, pattern)) {
        append(pos,position(finderDna));
        SEQAN_ASSERT_EQ(position(finderDna), beginPosition(finderDna));
        SEQAN_ASSERT_EQ(endPosition(finderDna), beginPosition(finderDna) + length(finderDna));
        SEQAN_ASSERT_EQ(length(finderDna), length(ndl));
        SEQAN_ASSERT_EQ(begin(finderDna), begin(hstk) + beginPosition(finderDna));
        SEQAN_ASSERT_EQ(end(finderDna), begin(hstk) + endPosition(finderDna));
        SEQAN_ASSERT_EQ(infix(finderDna), ndl);
    }

    SEQAN_ASSERT_EQ(pos[0], 0);
    SEQAN_ASSERT_EQ(pos[1], 1);
    SEQAN_ASSERT_EQ(pos[2], 2);
    SEQAN_ASSERT_EQ(pos[3], 3);
    SEQAN_ASSERT_EQ(pos[4], 4);
    SEQAN_ASSERT_EQ(pos[5], 5);
    SEQAN_ASSERT_EQ(pos[6], 8);
    SEQAN_ASSERT_EQ(length(pos), 7);

    //____________________________________________________________________________
    // Test3b - different alphabet, small needle, jumping finder

    goBegin(finderDna); // That's a repositioning
    clear(finderDna);       // That's why, clear state
    clear(pos);

    bool firstHit = true;
    while (find(finderDna, pattern)) {
        if (firstHit) {
            firstHit = false;
            finderDna += 2;
            clear(finderDna);  // clear the state of the finder
        } else {
            //unsigned int p = position(finderDna);
            append(pos,position(finderDna));
        }
    }

    SEQAN_ASSERT_EQ(pos[0], 2);
    SEQAN_ASSERT_EQ(pos[1], 3);
    SEQAN_ASSERT_EQ(pos[2], 4);
    SEQAN_ASSERT_EQ(pos[3], 5);
    SEQAN_ASSERT_EQ(pos[4], 8);
    SEQAN_ASSERT_EQ(length(pos), 5);

    //____________________________________________________________________________
    // Test4 - different alphabet, large needle
    String<Dna> text = "taaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaat";
    Finder<String<Dna> > finderText(text);

    String<Dna> query = "taaaataaaataaaataaaataaaataaaataaaataaaataaaat";
    setHost(pattern, query);

    clear(pos);
    while (find(finderText, pattern)) {
        append(pos,position(finderText));
        SEQAN_ASSERT_EQ(position(finderText), beginPosition(finderText));
        SEQAN_ASSERT_EQ(endPosition(finderText), beginPosition(finderText) + length(finderText));
        SEQAN_ASSERT_EQ(length(finderText), length(query));
        SEQAN_ASSERT_EQ(begin(finderText), begin(text) + beginPosition(finderText));
        SEQAN_ASSERT_EQ(end(finderText), begin(text) + endPosition(finderText));
        SEQAN_ASSERT_EQ(infix(finderText), query);
    }

    SEQAN_ASSERT_EQ(pos[0], 0);
    SEQAN_ASSERT_EQ(pos[1], 5);
    SEQAN_ASSERT_EQ(pos[2], 10);
    SEQAN_ASSERT_EQ(pos[3], 15);
    SEQAN_ASSERT_EQ(pos[4], 20);
    SEQAN_ASSERT_EQ(pos[5], 25);
    SEQAN_ASSERT_EQ(length(pos), 6);
}


template <typename TAlgorithmSpec>
void Test_OnlineAlgMulti(bool order_by_begin_position) {
    String<unsigned int> pos;

    //____________________________________________________________________________
    // Test1 - Single keyword
    String<char> haystack("Dies ist ein Haystack. Ja, das ist wirklich einer!");
    Finder<String<char> > finder(haystack);

    typedef String<String<char> > TNeedle;
    TNeedle keywords;
    appendValue(keywords, String<char>("ist"));
    Pattern<TNeedle, TAlgorithmSpec> pattern(keywords);

    while (find(finder, pattern)) {
        append(pos,position(finder));
        SEQAN_ASSERT_EQ(position(finder), beginPosition(finder));
        SEQAN_ASSERT_EQ(endPosition(finder), beginPosition(finder) + length(finder));
        SEQAN_ASSERT_EQ(length(finder), length(keywords[position(pattern)]));
        SEQAN_ASSERT_EQ(begin(finder), begin(haystack) + beginPosition(finder));
        SEQAN_ASSERT_EQ(end(finder), begin(haystack) + endPosition(finder));
        SEQAN_ASSERT_EQ(infix(finder), keywords[position(pattern)]);
    }

    SEQAN_ASSERT_EQ(host(pattern), keywords);
    SEQAN_ASSERT_EQ(host(reinterpret_cast<Pattern<TNeedle, TAlgorithmSpec> const &>(pattern)), keywords);
    SEQAN_ASSERT_EQ(pos[0], 5);
    SEQAN_ASSERT_EQ(pos[1], 31);
    SEQAN_ASSERT_EQ(length(pos), 2);

    haystack = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefgaabcdef";
    setHost(finder, haystack);
    clear(finder);

    clear(keywords);
    appendValue(keywords, String<char>("abcdefghijklmnopqrstuvwxyzabcdefg"));
    setHost(pattern, keywords);
    clear(pos);

    while (find(finder, pattern)) {
        append(pos,position(finder));
        SEQAN_ASSERT_EQ(position(finder), beginPosition(finder));
        SEQAN_ASSERT_EQ(endPosition(finder), beginPosition(finder) + length(finder));
        SEQAN_ASSERT_EQ(length(finder), length(keywords[position(pattern)]));
        SEQAN_ASSERT_EQ(begin(finder), begin(haystack) + beginPosition(finder));
        SEQAN_ASSERT_EQ(end(finder), begin(haystack) + endPosition(finder));
        SEQAN_ASSERT_EQ(infix(finder), keywords[position(pattern)]);
    }

    SEQAN_ASSERT_EQ(pos[0], 0);
    SEQAN_ASSERT_EQ(pos[1], 26);
    SEQAN_ASSERT_EQ(length(pos), 2);


    String<Dna> hstk = "aaaaaaacaa";
    Finder<String<Dna> > finderDna(hstk);

    typedef String<String<Dna> > TDnaNeedle;
    Pattern<TDnaNeedle, TAlgorithmSpec> pattern_dna(keywords);

    TDnaNeedle dna_keywords;
    appendValue(dna_keywords, String<Dna>("aa"));
    setHost(pattern_dna, dna_keywords);

    clear(pos);
    while (find(finderDna, pattern_dna)) {
        append(pos,position(finderDna));
        SEQAN_ASSERT_EQ(position(finderDna), beginPosition(finderDna));
        SEQAN_ASSERT_EQ(endPosition(finderDna), beginPosition(finderDna) + length(finderDna));
        SEQAN_ASSERT_EQ(length(finderDna), length(dna_keywords[position(pattern_dna)]));
        SEQAN_ASSERT_EQ(begin(finderDna), begin(hstk) + beginPosition(finderDna));
        SEQAN_ASSERT_EQ(end(finderDna), begin(hstk) + endPosition(finderDna));
        SEQAN_ASSERT_EQ(infix(finderDna), dna_keywords[position(pattern_dna)]);
    }

    SEQAN_ASSERT_EQ(pos[0], 0);
    SEQAN_ASSERT_EQ(pos[1], 1);
    SEQAN_ASSERT_EQ(pos[2], 2);
    SEQAN_ASSERT_EQ(pos[3], 3);
    SEQAN_ASSERT_EQ(pos[4], 4);
    SEQAN_ASSERT_EQ(pos[5], 5);
    SEQAN_ASSERT_EQ(pos[6], 8);
    SEQAN_ASSERT_EQ(length(pos), 7);

    String<Dna> text = "taaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaat";
    Finder<String<Dna> > finderText(text);

    clear(dna_keywords);
    appendValue(dna_keywords, String<Dna>("taaaataaaataaaataaaataaaataaaataaaataaaataaaat"));
    setHost(pattern_dna, dna_keywords);

    clear(pos);
    while (find(finderText, pattern_dna)) {
        append(pos,position(finderText));
        SEQAN_ASSERT_EQ(position(finderText), beginPosition(finderText));
        SEQAN_ASSERT_EQ(endPosition(finderText), beginPosition(finderText) + length(finderText));
        SEQAN_ASSERT_EQ(length(finderText), length(dna_keywords[position(pattern_dna)]));
        SEQAN_ASSERT_EQ(begin(finderText), begin(text) + beginPosition(finderText));
        SEQAN_ASSERT_EQ(end(finderText), begin(text) + endPosition(finderText));
        SEQAN_ASSERT_EQ(infix(finderText), dna_keywords[position(pattern_dna)]);
    }

    SEQAN_ASSERT_EQ(pos[0], 0);
    SEQAN_ASSERT_EQ(pos[1], 5);
    SEQAN_ASSERT_EQ(pos[2], 10);
    SEQAN_ASSERT_EQ(pos[3], 15);
    SEQAN_ASSERT_EQ(pos[4], 20);
    SEQAN_ASSERT_EQ(pos[5], 25);
    SEQAN_ASSERT_EQ(length(pos), 6);

    //____________________________________________________________________________
    // Test2 - Multiple keywords
    String<char> hst("annual_announce_any_annually");
    Finder<String<char> > fd(hst);

    typedef String<String<char> > TN;
    TN kyw;
    appendValue(kyw, String<char>("announce"));
    appendValue(kyw, String<char>("annual"));
    appendValue(kyw, String<char>("annually"));
    Pattern<TN, TAlgorithmSpec> pt(kyw);

    String<unsigned int> finderPos;
    String<unsigned int> keywordIndex;
    while (find(fd, pt)) {
        append(finderPos,position(fd));
        append(keywordIndex,position(pt));
        SEQAN_ASSERT_EQ(position(fd), beginPosition(fd));
        SEQAN_ASSERT_EQ(endPosition(fd), beginPosition(fd) + length(fd));
        SEQAN_ASSERT_EQ(length(fd), length(kyw[position(pt)]));
        SEQAN_ASSERT_EQ(begin(fd), begin(hst) + beginPosition(fd));
        SEQAN_ASSERT_EQ(end(fd), begin(hst) + endPosition(fd));
        SEQAN_ASSERT_EQ(infix(fd), kyw[position(pt)]);
    }

    SEQAN_ASSERT_EQ(length(finderPos), 4);
    SEQAN_ASSERT_EQ(length(keywordIndex), 4);
    SEQAN_ASSERT_EQ(finderPos[0], 0);
    SEQAN_ASSERT_EQ(keywordIndex[0], 1);
    SEQAN_ASSERT_EQ(finderPos[1], 7);
    SEQAN_ASSERT_EQ(keywordIndex[1], 0);
    SEQAN_ASSERT_EQ(finderPos[2], 20);
    SEQAN_ASSERT_EQ(keywordIndex[2], 1);
    SEQAN_ASSERT_EQ(finderPos[3], 20);
    SEQAN_ASSERT_EQ(keywordIndex[3], 2);

    String<Dna> hstDna("AGATACGATATATAC");
    Finder<String<Dna> > fdDna(hstDna);

    typedef String<String<Dna> > TNDna;
    TNDna kywDna;
    appendValue(kywDna, String<Dna>("ATATATA"));
    appendValue(kywDna, String<Dna>("TATAT"));
    appendValue(kywDna, String<Dna>("ACGATAT"));
    Pattern<TNDna, TAlgorithmSpec> ptDna(kywDna);

    clear(finderPos);
    clear(keywordIndex);
    while (find(fdDna, ptDna)) {
        append(finderPos,position(fdDna));
        append(keywordIndex,position(ptDna));
        SEQAN_ASSERT_EQ(position(fdDna), beginPosition(fdDna));
        SEQAN_ASSERT_EQ(endPosition(fdDna), beginPosition(fdDna) + length(fdDna));
        SEQAN_ASSERT_EQ(length(fdDna), length(kywDna[position(ptDna)]));
        SEQAN_ASSERT_EQ(begin(fdDna), begin(hstDna) + beginPosition(fdDna));
        SEQAN_ASSERT_EQ(end(fdDna), begin(hstDna) + endPosition(fdDna));
        SEQAN_ASSERT_EQ(infix(fdDna), kywDna[position(ptDna)]);
    }

    SEQAN_ASSERT_EQ(length(finderPos), 3);
    SEQAN_ASSERT_EQ(length(keywordIndex), 3);
    if (order_by_begin_position) {
        SEQAN_ASSERT_EQ(finderPos[0], 4);
        SEQAN_ASSERT_EQ(keywordIndex[0], 2);
        SEQAN_ASSERT_EQ(finderPos[1], 7);
        SEQAN_ASSERT_EQ(keywordIndex[1], 0);
        SEQAN_ASSERT_EQ(finderPos[2], 8);
        SEQAN_ASSERT_EQ(keywordIndex[2], 1);
    } else {
        SEQAN_ASSERT_EQ(finderPos[0], 4);
        SEQAN_ASSERT_EQ(keywordIndex[0], 2);
        SEQAN_ASSERT_EQ(finderPos[1], 8);
        SEQAN_ASSERT_EQ(keywordIndex[1], 1);
        SEQAN_ASSERT_EQ(finderPos[2], 7);
        SEQAN_ASSERT_EQ(keywordIndex[2], 0);
    }
    //____________________________________________________________________________
    // Test2 - Multiple keywords that do not fit into a machine word
    String<Dna> my_haystack("AGATACGATATATACAGATACGATATATACAGATACGATATATACAGATACGATATATACAGATACGATATATAC");
    Finder<String<Dna> > my_finder(my_haystack);

    typedef String<String<Dna> > TNeedle_My;
    TNeedle_My my_keywords;
    appendValue(my_keywords, String<Dna>("ATATATA"));
    appendValue(my_keywords, String<Dna>("ACCGATCCAT"));
    appendValue(my_keywords, String<Dna>("TATAT"));
    appendValue(my_keywords, String<Dna>("ACCGAT"));
    appendValue(my_keywords, String<Dna>("ACGATAT"));
    appendValue(my_keywords, String<Dna>("CCAA"));
    Pattern<TNeedle_My, TAlgorithmSpec> my_pattern(my_keywords);

    clear(finderPos);
    clear(keywordIndex);
    while (find(my_finder, my_pattern)) {
        //std::cout << position(my_finder) << "-" << position(my_pattern) << ::std::endl;
        append(finderPos,position(my_finder));
        append(keywordIndex,position(my_pattern));
        SEQAN_ASSERT_EQ(position(my_finder), beginPosition(my_finder));
        SEQAN_ASSERT_EQ(endPosition(my_finder), beginPosition(my_finder) + length(my_finder));
        SEQAN_ASSERT_EQ(length(my_finder), length(my_keywords[position(my_pattern)]));
        SEQAN_ASSERT_EQ(begin(my_finder), begin(my_haystack) + beginPosition(my_finder));
        SEQAN_ASSERT_EQ(end(my_finder), begin(my_haystack) + endPosition(my_finder));
        SEQAN_ASSERT_EQ(infix(my_finder), my_keywords[position(my_pattern)]);
    }

    SEQAN_ASSERT_EQ(length(finderPos), 15);
    SEQAN_ASSERT_EQ(length(keywordIndex), 15);
    if (order_by_begin_position) {
        SEQAN_ASSERT_EQ(finderPos[0], 4);
        SEQAN_ASSERT_EQ(keywordIndex[0], 4);
        SEQAN_ASSERT_EQ(finderPos[1], 7);
        SEQAN_ASSERT_EQ(keywordIndex[1], 0);
        SEQAN_ASSERT_EQ(finderPos[2], 8);
        SEQAN_ASSERT_EQ(keywordIndex[2], 2);
        SEQAN_ASSERT_EQ(finderPos[3], 19);
        SEQAN_ASSERT_EQ(keywordIndex[3], 4);
        SEQAN_ASSERT_EQ(finderPos[4], 22);
        SEQAN_ASSERT_EQ(keywordIndex[4], 0);
        SEQAN_ASSERT_EQ(finderPos[5], 23);
        SEQAN_ASSERT_EQ(keywordIndex[5], 2);
        SEQAN_ASSERT_EQ(finderPos[6], 34);
        SEQAN_ASSERT_EQ(keywordIndex[6], 4);
        SEQAN_ASSERT_EQ(finderPos[7], 37);
        SEQAN_ASSERT_EQ(keywordIndex[7], 0);
        SEQAN_ASSERT_EQ(finderPos[8], 38);
        SEQAN_ASSERT_EQ(keywordIndex[8], 2);
        SEQAN_ASSERT_EQ(finderPos[9], 49);
        SEQAN_ASSERT_EQ(keywordIndex[9], 4);
        SEQAN_ASSERT_EQ(finderPos[10], 52);
        SEQAN_ASSERT_EQ(keywordIndex[10], 0);
        SEQAN_ASSERT_EQ(finderPos[11], 53);
        SEQAN_ASSERT_EQ(keywordIndex[11], 2);
        SEQAN_ASSERT_EQ(finderPos[12], 64);
        SEQAN_ASSERT_EQ(keywordIndex[12], 4);
        SEQAN_ASSERT_EQ(finderPos[13], 67);
        SEQAN_ASSERT_EQ(keywordIndex[13], 0);
        SEQAN_ASSERT_EQ(finderPos[14], 68);
        SEQAN_ASSERT_EQ(keywordIndex[14], 2);
    } else {
        SEQAN_ASSERT_EQ(finderPos[0], 4);
        SEQAN_ASSERT_EQ(keywordIndex[0], 4);
        SEQAN_ASSERT_EQ(finderPos[1], 8);
        SEQAN_ASSERT_EQ(keywordIndex[1], 2);
        SEQAN_ASSERT_EQ(finderPos[2], 7);
        SEQAN_ASSERT_EQ(keywordIndex[2], 0);
        SEQAN_ASSERT_EQ(finderPos[3], 19);
        SEQAN_ASSERT_EQ(keywordIndex[3], 4);
        SEQAN_ASSERT_EQ(finderPos[4], 23);
        SEQAN_ASSERT_EQ(keywordIndex[4], 2);
        SEQAN_ASSERT_EQ(finderPos[5], 22);
        SEQAN_ASSERT_EQ(keywordIndex[5], 0);
        SEQAN_ASSERT_EQ(finderPos[6], 34);
        SEQAN_ASSERT_EQ(keywordIndex[6], 4);
        SEQAN_ASSERT_EQ(finderPos[7], 38);
        SEQAN_ASSERT_EQ(keywordIndex[7], 2);
        SEQAN_ASSERT_EQ(finderPos[8], 37);
        SEQAN_ASSERT_EQ(keywordIndex[8], 0);
        SEQAN_ASSERT_EQ(finderPos[9], 49);
        SEQAN_ASSERT_EQ(keywordIndex[9], 4);
        SEQAN_ASSERT_EQ(finderPos[10], 53);
        SEQAN_ASSERT_EQ(keywordIndex[10], 2);
        SEQAN_ASSERT_EQ(finderPos[11], 52);
        SEQAN_ASSERT_EQ(keywordIndex[11], 0);
        SEQAN_ASSERT_EQ(finderPos[12], 64);
        SEQAN_ASSERT_EQ(keywordIndex[12], 4);
        SEQAN_ASSERT_EQ(finderPos[13], 68);
        SEQAN_ASSERT_EQ(keywordIndex[13], 2);
        SEQAN_ASSERT_EQ(finderPos[14], 67);
        SEQAN_ASSERT_EQ(keywordIndex[14], 0);
    }

    //____________________________________________________________________________
    // Multiple keywords with overlapping matches
    String<Dna> my2_haystack("aaaacaaa");
    Finder<String<Dna> > my2_finder(my2_haystack);

    typedef String<String<Dna> > TNeedle_My2;
    TNeedle_My2 my2_keywords;
    appendValue(my2_keywords, String<Dna>("aa"));
    appendValue(my2_keywords, String<Dna>("aaa"));
    appendValue(my2_keywords, String<Dna>("ac"));
    appendValue(my2_keywords, String<Dna>("aac"));
    appendValue(my2_keywords, String<Dna>("gctccacctgacctagcccatggggcccaaatttccggccttaattcccattt"));
    Pattern<TNeedle_My2, TAlgorithmSpec> my2_pattern(my2_keywords);

    clear(finderPos);
    clear(keywordIndex);
    while (find(my2_finder, my2_pattern)) {
        //std::cout << position(my2_finder) << ":" << position(my2_pattern) << ::std::endl;
        append(finderPos,position(my2_finder));
        append(keywordIndex,position(my2_pattern));
        SEQAN_ASSERT_EQ(position(my2_finder), beginPosition(my2_finder));
        SEQAN_ASSERT_EQ(endPosition(my2_finder), beginPosition(my2_finder) + length(my2_finder));
        SEQAN_ASSERT_EQ(length(my2_finder), length(my2_keywords[position(my2_pattern)]));
        SEQAN_ASSERT_EQ(begin(my2_finder), begin(my2_haystack) + beginPosition(my2_finder));
        SEQAN_ASSERT_EQ(end(my2_finder), begin(my2_haystack) + endPosition(my2_finder));
        SEQAN_ASSERT_EQ(infix(my2_finder), my2_keywords[position(my2_pattern)]);
    }

    if (order_by_begin_position) {
        SEQAN_ASSERT_EQ(finderPos[0], 0);
        SEQAN_ASSERT_EQ(keywordIndex[0], 0);
        SEQAN_ASSERT_EQ(finderPos[1], 0);
        SEQAN_ASSERT_EQ(keywordIndex[1], 1);
        SEQAN_ASSERT_EQ(finderPos[2], 1);
        SEQAN_ASSERT_EQ(keywordIndex[2], 0);
        SEQAN_ASSERT_EQ(finderPos[3], 1);
        SEQAN_ASSERT_EQ(keywordIndex[3], 1);
        SEQAN_ASSERT_EQ(finderPos[4], 2);
        SEQAN_ASSERT_EQ(keywordIndex[4], 0);
        SEQAN_ASSERT_EQ(finderPos[5], 2);
        SEQAN_ASSERT_EQ(keywordIndex[5], 3);
        SEQAN_ASSERT_EQ(finderPos[6], 3);
        SEQAN_ASSERT_EQ(keywordIndex[6], 2);
        SEQAN_ASSERT_EQ(finderPos[7], 5);
        SEQAN_ASSERT_EQ(keywordIndex[7], 0);
        SEQAN_ASSERT_EQ(finderPos[8], 5);
        SEQAN_ASSERT_EQ(keywordIndex[8], 1);
        SEQAN_ASSERT_EQ(finderPos[9], 6);
        SEQAN_ASSERT_EQ(keywordIndex[9], 0);
    } else{
        SEQAN_ASSERT_EQ(finderPos[0], 0);
        SEQAN_ASSERT_EQ(keywordIndex[0], 0);
        SEQAN_ASSERT_EQ(finderPos[1], 1);
        SEQAN_ASSERT_EQ(keywordIndex[1], 0);
        SEQAN_ASSERT_EQ(finderPos[2], 0);
        SEQAN_ASSERT_EQ(keywordIndex[2], 1);
        SEQAN_ASSERT_EQ(finderPos[3], 2);
        SEQAN_ASSERT_EQ(keywordIndex[3], 0);
        SEQAN_ASSERT_EQ(finderPos[4], 1);
        SEQAN_ASSERT_EQ(keywordIndex[4], 1);
        SEQAN_ASSERT_EQ(finderPos[5], 3);
        SEQAN_ASSERT_EQ(keywordIndex[5], 2);
        SEQAN_ASSERT_EQ(finderPos[6], 2);
        SEQAN_ASSERT_EQ(keywordIndex[6], 3);
        SEQAN_ASSERT_EQ(finderPos[7], 5);
        SEQAN_ASSERT_EQ(keywordIndex[7], 0);
        SEQAN_ASSERT_EQ(finderPos[8], 6);
        SEQAN_ASSERT_EQ(keywordIndex[8], 0);
        SEQAN_ASSERT_EQ(finderPos[9], 5);
        SEQAN_ASSERT_EQ(keywordIndex[9], 1);
    }

    //____________________________________________________________________________
    // Multiple duplicated keywords with overlapping matches, jumping finder
    goBegin(my2_finder); // That's a repositioning
    clear(my2_finder);      // That's why, clear state
    clear(finderPos);
    clear(keywordIndex);

    unsigned int hits = 0;
    while (find(my2_finder, my2_pattern)) {
        if (hits, 2) break;
        //std::cout << position(my2_finder) << ":" << position(my2_pattern) << ::std::endl;
        append(finderPos,position(my2_finder));
        append(keywordIndex,position(my2_pattern));
        ++hits;
    }
    goBegin(my2_finder);
    my2_finder+=5;
    clear(my2_finder);
    clear(finderPos);
    clear(keywordIndex);
    while (find(my2_finder, my2_pattern)) {
        //std::cout << position(my2_finder) << ":" << position(my2_pattern) << ::std::endl;
        append(finderPos,position(my2_finder));
        append(keywordIndex,position(my2_pattern));
    }

    if (order_by_begin_position) {
        SEQAN_ASSERT_EQ(finderPos[0], 5);
        SEQAN_ASSERT_EQ(keywordIndex[0], 0);
        SEQAN_ASSERT_EQ(finderPos[1], 5);
        SEQAN_ASSERT_EQ(keywordIndex[1], 1);
        SEQAN_ASSERT_EQ(finderPos[2], 6);
        SEQAN_ASSERT_EQ(keywordIndex[2], 0);
    } else {
        SEQAN_ASSERT_EQ(finderPos[0], 5);
        SEQAN_ASSERT_EQ(keywordIndex[0], 0);
        SEQAN_ASSERT_EQ(finderPos[1], 6);
        SEQAN_ASSERT_EQ(keywordIndex[1], 0);
        SEQAN_ASSERT_EQ(finderPos[2], 5);
        SEQAN_ASSERT_EQ(keywordIndex[2], 1);
    }
}


template <typename TAlgorithmSpec>
void Test_OnlineAlgWildcards() {
    String<unsigned int> pos;

    //____________________________________________________________________________
    // Test1 - simple find wo wildcards
    String<char> haystack("Dies ist ein Haystack. Ja, das ist wirklich einer!");
    Finder<String<char> > finder(haystack);

    String<char> needle("ist");
    Pattern<String<char>, TAlgorithmSpec> pattern(needle);
    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));

    SEQAN_ASSERT_EQ(pos[0], 7);
    SEQAN_ASSERT_EQ(pos[1], 33);
    SEQAN_ASSERT_EQ(host(pattern), needle);
    SEQAN_ASSERT_EQ(host(reinterpret_cast<Pattern<String<char>, TAlgorithmSpec> const &>(pattern)), needle);

    //____________________________________________________________________________
    // Test - validation of patterns
    needle = "ist";
    setHost(pattern, needle);
    SEQAN_ASSERT_EQ(valid(pattern), true);
    SEQAN_ASSERT_EQ(valid(reinterpret_cast<Pattern<String<char>, TAlgorithmSpec> const &>(pattern)), true);

    needle = "i[a-z]s{3,4}t?a*a+c..\\a";
    setHost(pattern, needle);
    SEQAN_ASSERT_EQ(valid(pattern), true);

    needle = "i[st";
    setHost(pattern, needle);
    SEQAN_ASSERT_EQ(valid(pattern), false);

    needle = "ist\\";
    setHost(pattern, needle);
    SEQAN_ASSERT_EQ(valid(pattern), false);

    needle = "ist?*";
    setHost(pattern, needle);
    SEQAN_ASSERT_EQ(valid(pattern), false);

    needle = "ist{5,4}";
    setHost(pattern, needle);
    SEQAN_ASSERT_EQ(valid(pattern), false);

    needle = "ist{5,a}";
    setHost(pattern, needle);
    SEQAN_ASSERT_EQ(valid(pattern), false);

    needle = "ist{5,";
    setHost(pattern, needle);
    SEQAN_ASSERT_EQ(valid(pattern), false);

    //____________________________________________________________________________
    // Test - searching with invalid needles
    haystack = "Dies i[st ein Haystack. Ja, das i[st wirklich einer!";
    setHost(finder, haystack);
    clear(finder);

    needle = "i[st";
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));
        
    SEQAN_ASSERT_EQ(length(pos), 0);

    //____________________________________________________________________________
    // Test - handle needles with wildcards
    // to produce two \ in the pattern you need to escape both of them
    needle = "aa+c*[a-z]xx?aa\\\\";
    SEQAN_ASSERT_EQ(_length_wo_wild(needle), 9);
        
    //____________________________________________________________________________
    // Test - optional characters (?)
    haystack = "abc__ac";
    setHost(finder, haystack);
    clear(finder);

    needle = "ab?c";
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));
        
    SEQAN_ASSERT_EQ(pos[0], 2);
    SEQAN_ASSERT_EQ(pos[1], 6);
    SEQAN_ASSERT_EQ(length(pos), 2);

    //____________________________________________________________________________
    // Test - repeatable characters (+)
    haystack = "abc__abbbbbc_ac";
    setHost(finder, haystack);
    clear(finder);

    needle = "ab+c";
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));
        
    SEQAN_ASSERT_EQ(pos[0], 2);
    SEQAN_ASSERT_EQ(pos[1], 11);
    SEQAN_ASSERT_EQ(length(pos), 2);

    //____________________________________________________________________________
    // Test - repeatable characters (*)
    haystack = "abc__abbbbbc_ac";
    setHost(finder, haystack);
    clear(finder);

    needle = "ab*c";
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));
        
    SEQAN_ASSERT_EQ(pos[0], 2);
    SEQAN_ASSERT_EQ(pos[1], 11);
    SEQAN_ASSERT_EQ(pos[2], 14);
        
    SEQAN_ASSERT_EQ(length(pos), 3);

    //____________________________________________________________________________
    // Test - wildcard matching
    haystack = "acccdfabdeeef";
    setHost(finder, haystack);
    clear(finder);

    needle = "ab?c*de+f";
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));
        
    SEQAN_ASSERT_EQ(pos[0], 12);
    SEQAN_ASSERT_EQ(length(pos), 1);

    //____________________________________________________________________________
    // Test - wildcard matching (hard case)
    haystack = "aacccdfacccdeeef";
    setHost(finder, haystack);
    clear(finder);

    needle = "a*c*de+f";
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));
        
    SEQAN_ASSERT_EQ(pos[0], 15);
    SEQAN_ASSERT_EQ(length(pos), 1);


    //____________________________________________________________________________
    // Test - character classes matching
    haystack = "annual_Annual_znnual_Znnual";
    setHost(finder, haystack);
    clear(finder);

    needle = "[a-zA]nnual"; 
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));
        
    SEQAN_ASSERT_EQ(pos[0], 5);
    SEQAN_ASSERT_EQ(pos[1], 12);
    SEQAN_ASSERT_EQ(pos[2], 19);
    SEQAN_ASSERT_EQ(length(pos), 3);

    //____________________________________________________________________________
    // Test - long needles
    haystack = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefgaabcdef";
    setHost(finder, haystack);
    clear(finder);

    needle = "abcdefghijklmnopqrstuvwxyzabcdefg";   
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));

    SEQAN_ASSERT_EQ(pos[0], 32);
    SEQAN_ASSERT_EQ(pos[1], 58);
    SEQAN_ASSERT_EQ(length(pos), 2);

    //____________________________________________________________________________
    // Test - long needles with character classes
    //              abcdefghijklmnopqrstuvwxyzabcdefghijkl
    //                                        abcdefghijklmnopqrstuvwxyzabcdefghijkl
    haystack = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuzwxyzabcdefghzjklaabcdefhijkl";
    setHost(finder, haystack);
    clear(finder);

    needle = "abcdefghijklmnopqrstu[vz]wxyzabcdefgh[iz]jkl";        
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));

    SEQAN_ASSERT_EQ(pos[0], 37);
    SEQAN_ASSERT_EQ(pos[1], 63);
    SEQAN_ASSERT_EQ(length(pos), 2);


    //____________________________________________________________________________
    // Test - long needles with repeating characters
    //              abcdefghijklmnopqrstuvwxyzabcdefghijkl
    //                                                                                                        abcdefghijklmnopqrstuvwxyzabcdefghijkl
    haystack = "abcdefghijklmnopqrstuvwxyzabcdefghiiiiijkl____aaaaabcdefghijklmnopqrstuvwxyzabcdeghijkl__aaaaabcdefghijklmnopqrstuvwxyzabcdefghjkl";
    setHost(finder, haystack);
    clear(finder);

    needle = "aa*bcdefghijklmnopqrstuvwxyzabcdef?g?hi+jkl"; 
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));

    SEQAN_ASSERT_EQ(pos[0], 41);
    SEQAN_ASSERT_EQ(pos[1], 86);
    SEQAN_ASSERT_EQ(length(pos), 2);

    //____________________________________________________________________________
    // Test - handle .
    haystack = "annual_Annual_znnual";
    setHost(finder, haystack);
    clear(finder);

    needle = ".nnual";      
    setHost(pattern, needle);
    clear(pos);
    while (find(finder, pattern)) {
        append(pos,position(finder));
    }
    SEQAN_ASSERT_EQ(pos[0], 5);
    SEQAN_ASSERT_EQ(pos[1], 12);
    SEQAN_ASSERT_EQ(pos[2], 19);
    SEQAN_ASSERT_EQ(length(pos), 3);

    //____________________________________________________________________________
    // Test - handle backslash 
    haystack = "annual_Annual_.nnual";
    setHost(finder, haystack);
    clear(finder);

    needle = "\\.nnual";    
    setHost(pattern, needle);
    clear(pos);
    while (find(finder, pattern)){
        append(pos,position(finder));
    }
    SEQAN_ASSERT_EQ(pos[0], 19);
    SEQAN_ASSERT_EQ(length(pos), 1);



    //____________________________________________________________________________
    // Test - handle bounded length repeats {n,m}
    haystack = "aannual_aaaannual_annual";
    setHost(finder, haystack);
    clear(finder);

    needle = "a{2,5}n{2}ual";       
        
    SEQAN_ASSERT_EQ(_length_wo_wild(needle), 10);
        
    setHost(pattern, needle);
    clear(pos);
    while (find(finder, pattern)) {
        append(pos,position(finder));
    }
    SEQAN_ASSERT_EQ(pos[0], 6);
    SEQAN_ASSERT_EQ(pos[1], 16);
    SEQAN_ASSERT_EQ(length(pos), 2);


    //____________________________________________________________________________
    // Test - handle different types of Pattern and Needle
    String<Dna> dna_haystack("AAACCTATGGGTTTAAAACCCTGAAACCCC");
    Finder<String<Dna> > dna_finder(dna_haystack);

    String<char> char_needle("a{3}c+t[ag].");
    Pattern<String<Dna>, TAlgorithmSpec> dna_pattern(char_needle);
    clear(pos);

    while (find(dna_finder, dna_pattern))
        append(pos,position(dna_finder));

    SEQAN_ASSERT_EQ(pos[0], 7);
    SEQAN_ASSERT_EQ(pos[1], 23);
    SEQAN_ASSERT_EQ(length(pos), 2);

}


template <typename TPatternSpec>
void Test_Approx_EditDist() {
    //test DPSearch
    String<char> hstk("any_annealing");
    String<char> nl("annual");

    Finder<String<char> > fd(hstk);

    Pattern<String<char>, TPatternSpec> pt(nl, -2);
    SEQAN_ASSERT_TRUE(find(fd, pt));
    SEQAN_ASSERT_EQ(position(fd), 8);
    SEQAN_ASSERT_EQ(getScore(pt), -2);
    SEQAN_ASSERT_TRUE(find(fd, pt));
    SEQAN_ASSERT_EQ(position(fd), 9);
    SEQAN_ASSERT_EQ(getScore(pt), -1);
    SEQAN_ASSERT_TRUE(find(fd, pt));
    SEQAN_ASSERT_EQ(position(fd), 10);
    SEQAN_ASSERT_EQ(getScore(pt), -2);

    SEQAN_ASSERT_NOT(find(fd,pt));

    String<char> haystk("Dies ist der Haystack des Tests. Ja, das ist er wirklich!");
    String<char> ndl("des");

    Finder<String<char> > fnd(haystk);

    Pattern<String<char>, TPatternSpec> pat(ndl, -2);
    SEQAN_ASSERT_EQ(host(pat), ndl);
    SEQAN_ASSERT_EQ(host(reinterpret_cast<Pattern<String<char>, TPatternSpec> const &>(pat)), ndl);

    SEQAN_ASSERT_EQ(scoreLimit(pat), -2);
    setScoreLimit(pat, -1);;
    SEQAN_ASSERT_EQ(scoreLimit(pat), -1);

    SEQAN_ASSERT_TRUE(find(fnd, pat));
    SEQAN_ASSERT_EQ(position(fnd), 3);
    SEQAN_ASSERT_EQ(getScore(pat), -1);

    SEQAN_ASSERT_TRUE(find(fnd, pat));
    SEQAN_ASSERT_EQ(position(fnd), 10);
    SEQAN_ASSERT_EQ(getScore(pat), -1);

    SEQAN_ASSERT_TRUE(find(fnd, pat));
    SEQAN_ASSERT_EQ(position(fnd), 11);
    SEQAN_ASSERT_EQ(getScore(pat), -1);

    SEQAN_ASSERT_TRUE(find(fnd, pat));
    SEQAN_ASSERT_EQ(position(fnd), 23);
    SEQAN_ASSERT_EQ(getScore(pat), -1);

    SEQAN_ASSERT_TRUE(find(fnd, pat));
    SEQAN_ASSERT_EQ(position(fnd), 24);
    SEQAN_ASSERT_EQ(getScore(pat), 0);

    SEQAN_ASSERT_TRUE(find(fnd, pat));
    SEQAN_ASSERT_EQ(position(fnd), 25);
    SEQAN_ASSERT_EQ(getScore(pat), -1);

    SEQAN_ASSERT_TRUE(find(fnd, pat));
    SEQAN_ASSERT_EQ(position(fnd), 28);
    SEQAN_ASSERT_EQ(getScore(pat), -1);

    SEQAN_ASSERT_TRUE(find(fnd, pat));
    SEQAN_ASSERT_EQ(position(fnd), 39);
    SEQAN_ASSERT_EQ(getScore(pat), -1);

    SEQAN_ASSERT_NOT(find(fnd, pat));

    // Test with long needles and a Dna Alphabet
    String<Dna> long_haystk("taaaataaaatacaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaat");
    String<Dna> long_ndl("taaaataaaatacaataaaataaaatataataaaataaaataaaat");

    Finder<String<Dna> > long_fnd(long_haystk);

    Pattern<String<Dna>, TPatternSpec> long_pat(long_ndl, -2);

    SEQAN_ASSERT_TRUE(find(long_fnd,long_pat));
    SEQAN_ASSERT_EQ(position(long_fnd), 44);
    SEQAN_ASSERT_EQ(getScore(long_pat), -2);

    SEQAN_ASSERT_TRUE(find(long_fnd,long_pat));
    SEQAN_ASSERT_EQ(position(long_fnd), 45);
    SEQAN_ASSERT_EQ(getScore(long_pat), -1);

    SEQAN_ASSERT_TRUE(find(long_fnd,long_pat));
    SEQAN_ASSERT_EQ(position(long_fnd), 46);
    SEQAN_ASSERT_EQ(getScore(long_pat), -2);

    SEQAN_ASSERT_TRUE(find(long_fnd,long_pat));
    SEQAN_ASSERT_EQ(position(long_fnd), 60);
    SEQAN_ASSERT_EQ(getScore(long_pat), -2);

    SEQAN_ASSERT_TRUE(find(long_fnd,long_pat));
    SEQAN_ASSERT_EQ(position(long_fnd), 65);
    SEQAN_ASSERT_EQ(getScore(long_pat), -2);

    SEQAN_ASSERT_TRUE(find(long_fnd,long_pat));
    SEQAN_ASSERT_EQ(position(long_fnd), 70);
    SEQAN_ASSERT_EQ(getScore(long_pat), -2);

    SEQAN_ASSERT_NOT(find(long_fnd,long_pat));

    //____________________________________________________________________________

    String<char> haystack_1 = "123XXXabaXXX45aba123";
    String<char> needle_1 = "XXXaba";
    Finder<String<char> > finder_1(haystack_1);
    Pattern<String<char>, TPatternSpec> pattern_1(needle_1, -2);

    SEQAN_ASSERT_TRUE(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(endPosition(finder_1), 7);
    SEQAN_ASSERT_EQ(getScore(pattern_1), -2);
    SEQAN_ASSERT_TRUE(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXXa");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT_NOT(findBegin(finder_1, pattern_1));

    SEQAN_ASSERT_TRUE(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(endPosition(finder_1), 8);
    SEQAN_ASSERT_EQ(getScore(pattern_1), -1);
    SEQAN_ASSERT_TRUE(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXab");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT_TRUE(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXXab");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -1);
    SEQAN_ASSERT_TRUE(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "3XXXab");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT_NOT(findBegin(finder_1, pattern_1));

    SEQAN_ASSERT_TRUE(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(endPosition(finder_1), 9);
    SEQAN_ASSERT_EQ(getScore(pattern_1), 0);
    SEQAN_ASSERT_TRUE(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "Xaba");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT_TRUE(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXaba");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -1);
    SEQAN_ASSERT_TRUE(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXXaba");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), 0);
    SEQAN_ASSERT_TRUE(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "3XXXaba");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -1);
    SEQAN_ASSERT_TRUE(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "23XXXaba");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT_NOT(findBegin(finder_1, pattern_1));

    SEQAN_ASSERT_TRUE(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(endPosition(finder_1), 10);
    SEQAN_ASSERT_EQ(getScore(pattern_1), -1);
    SEQAN_ASSERT_TRUE(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXabaX");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT_TRUE(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXXabaX");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -1);
    SEQAN_ASSERT_TRUE(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "3XXXabaX");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT_NOT(findBegin(finder_1, pattern_1));

    SEQAN_ASSERT_TRUE(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(endPosition(finder_1), 11);
    SEQAN_ASSERT_EQ(getScore(pattern_1), -2);
    SEQAN_ASSERT_TRUE(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXXabaXX");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT_NOT(findBegin(finder_1, pattern_1));

    SEQAN_ASSERT_TRUE(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(endPosition(finder_1), 15);
    SEQAN_ASSERT_EQ(getScore(pattern_1), -2);
    SEQAN_ASSERT_TRUE(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXX45a");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT_NOT(findBegin(finder_1, pattern_1));

    SEQAN_ASSERT_TRUE(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(endPosition(finder_1), 17);
    SEQAN_ASSERT_EQ(getScore(pattern_1), -2);
    SEQAN_ASSERT_TRUE(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "X45aba");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT_TRUE(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XX45aba");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT_TRUE(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXX45aba");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT_NOT(findBegin(finder_1, pattern_1));

    SEQAN_ASSERT_NOT(find(finder_1, pattern_1));
}


// Test prefix search.
template <typename TPatternSpec>
void Test_Approx_Prefix_EditDist() {
    String<char> haystack_1 = "mississippi";
    String<char> needle_1 = "misssi";
    Finder<String<char> > finder_1(haystack_1);
    Pattern<String<char>, TPatternSpec> pattern_1(needle_1, -2);
    SEQAN_ASSERT_TRUE(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(position(finder_1), 3);
    SEQAN_ASSERT_EQ(length(finder_1), 4);
    SEQAN_ASSERT_EQ(beginPosition(finder_1), 0);
    SEQAN_ASSERT_EQ(endPosition(finder_1), 4);
    SEQAN_ASSERT_EQ(infix(finder_1), "miss");
    SEQAN_ASSERT_TRUE(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(position(finder_1), 4);
    SEQAN_ASSERT_EQ(infix(finder_1), "missi");
    SEQAN_ASSERT_TRUE(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(position(finder_1), 5);
    SEQAN_ASSERT_EQ(infix(finder_1), "missis");
    SEQAN_ASSERT_TRUE(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(position(finder_1), 6);
    SEQAN_ASSERT_TRUE(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(position(finder_1), 7);
    SEQAN_ASSERT_NOT(find(finder_1, pattern_1));


    String<char> haystack_2 = "yyyXXaba";
    String<char> needle_2 = "yyyaba";
    Finder<String<char> > finder_2(haystack_2);
    Pattern<String<char>, TPatternSpec> pattern_2(needle_2, -2);
    SEQAN_ASSERT_TRUE(find(finder_2, pattern_2));
    SEQAN_ASSERT_EQ(position(finder_2), 5);
    SEQAN_ASSERT_EQ(infix(finder_2), "yyyXXa");
    SEQAN_ASSERT_TRUE(find(finder_2, pattern_2));
    SEQAN_ASSERT_EQ(position(finder_2), 7);
    SEQAN_ASSERT_EQ(infix(finder_2), "yyyXXaba");
    SEQAN_ASSERT_NOT(find(finder_2, pattern_2));


    String<char> haystack_3 = "testtexttext";
    String<char> needle_3 = "mismatch";
    Finder<String<char> > finder_3(haystack_3);
    Pattern<String<char>, TPatternSpec> pattern_3(needle_3, -2);
    SEQAN_ASSERT_NOT(find(finder_3, pattern_3));


    String<char> haystack_4 = "testtext";
    String<char> needle_4 = "a longer mismatch";
    Finder<String<char> > finder_4(haystack_4);
    Pattern<String<char>, TPatternSpec> pattern_4(needle_4, -2);
    SEQAN_ASSERT_NOT(find(finder_4, pattern_4));


    String<char> haystack_5 = "exactmatching";
    String<char> needle_5 = "exact";
    Finder<String<char> > finder_5(haystack_5);
    Pattern<String<char>, TPatternSpec> pattern_5(needle_5, 0);
    SEQAN_ASSERT_TRUE(find(finder_5, pattern_5));
    SEQAN_ASSERT_EQ(position(finder_5), 4);
    SEQAN_ASSERT_NOT(find(finder_5, pattern_5));


    String<char> haystack_6 = "this is a text that is a bit longer than one machine word of 32 or 64 bits. AAXYX";
    String<char> needle_6 =   "this is a text that is a bit longer than one machine word of 32 or 64 bits. XYX";
    Finder<String<char> > finder_6(haystack_6);
    Pattern<String<char>, TPatternSpec> pattern_6(needle_6, -2);
    SEQAN_ASSERT_TRUE(find(finder_6, pattern_6));
    SEQAN_ASSERT_EQ(infix(finder_6), "this is a text that is a bit longer than one machine word of 32 or 64 bits. AAX");
    SEQAN_ASSERT_TRUE(find(finder_6, pattern_6));
    SEQAN_ASSERT_EQ(infix(finder_6), "this is a text that is a bit longer than one machine word of 32 or 64 bits. AAXYX");
    SEQAN_ASSERT_NOT(find(finder_6, pattern_6));

}


SEQAN_DEFINE_TEST(test_find_online_Simple) {
    Test_OnlineAlg<Simple>();   
}


SEQAN_DEFINE_TEST(test_find_online_Horspool) {
    Test_OnlineAlg<Horspool>();   
}


SEQAN_DEFINE_TEST(test_find_online_ShiftAnd) {
    Test_OnlineAlg<ShiftAnd>();   
}


SEQAN_DEFINE_TEST(test_find_online_ShiftOr) {
    Test_OnlineAlg<ShiftOr>();   
}


SEQAN_DEFINE_TEST(test_find_online_BndmAlgo) {
    Test_OnlineAlg<BndmAlgo>();   
}


SEQAN_DEFINE_TEST(test_find_online_BFAM_Oracle) {
    Test_OnlineAlg<BFAM<Oracle> >();   
}


SEQAN_DEFINE_TEST(test_find_online_BFAM_Trie) {
    Test_OnlineAlg<BFAM<Trie> >();   
}


SEQAN_DEFINE_TEST(test_find_online_wildcards) {
    Test_OnlineAlgWildcards<WildShiftAnd>();
}


SEQAN_DEFINE_TEST(test_find_online_multi_AhoCorasick) {
    Test_OnlineAlgMulti<AhoCorasick>(false);
}


SEQAN_DEFINE_TEST(test_find_online_multi_MultipleShiftAnd) {
    // TODO(holtgrew): Original comment: "leaks".
    // TODO(holtgrew): Fails, but was commented out in original code.
    // Test_OnlineAlgMulti<MultipleShiftAnd>(false);
}


SEQAN_DEFINE_TEST(test_find_online_multi_SetHorspool) {
    // TODO(holtgrew): Original comment: Does not compile.
    // TODO(holtgrew): Crashes, but was commented out in original code.
    // Test_OnlineAlgMulti<SetHorspool>(false);
}


SEQAN_DEFINE_TEST(test_find_online_multi_WuManber) {
    Test_OnlineAlgMulti<WuManber>(true);
}


SEQAN_DEFINE_TEST(test_find_online_multi_MultiBFAM_Oracle) {
    Test_OnlineAlgMulti<MultiBFAM<Oracle> >(true);
}


SEQAN_DEFINE_TEST(test_find_online_multi_MultiBFAM_Trie) {
    Test_OnlineAlgMulti<MultiBFAM<Trie> >(true);
}


SEQAN_DEFINE_TEST(test_find_approx_prefix_edit_dist_dpsearch) {
    Test_Approx_Prefix_EditDist<DPSearch<Score<>, FindPrefix> >();
}


SEQAN_DEFINE_TEST(test_approx_prefix_edit_dist_myers) {
    Test_Approx_Prefix_EditDist<Myers<FindPrefix> >();
}


SEQAN_DEFINE_TEST(test_approx_edit_dist_dp_search_simple_score) {
    Test_Approx_EditDist<DPSearch<SimpleScore> >();
}


SEQAN_DEFINE_TEST(test_approx_edit_dp_search_simple_score_legacy_case) {
    // TODO(holtgrew): This was written out like this in Test_Approx() in original code.
    // Test DPSearch.
    Pattern<String<char>, DPSearch<SimpleScore> > pat1;

    SimpleScore sc;
    scoreGap(sc) = -10;
    setScoringScheme(pat1, sc);
    SEQAN_ASSERT_EQ(scoreGap(scoringScheme(pat1)), -10);

    scoreGap(sc) = -1;
    SEQAN_ASSERT_EQ(scoreGap(scoringScheme(pat1)), -10);
    setScoringScheme(pat1, sc);
    SEQAN_ASSERT_EQ(scoreGap(scoringScheme(pat1)), -1);
}


SEQAN_DEFINE_TEST(test_approx_edit_dist_myers) {
    Test_Approx_EditDist<Myers<> >();
}


SEQAN_DEFINE_TEST(test_approx_edit_dist_abndm_algo) {
    Test_Approx_EditDist<AbndmAlgo >();
}


SEQAN_DEFINE_TEST(test_approx_edit_dist_pex_non_hierarchical) {
    Test_Approx_EditDist<PexNonHierarchical>();
}


SEQAN_DEFINE_TEST(test_approx_edit_dist_pex_hierarchical) {
    Test_Approx_EditDist<PexHierarchical>();
}


SEQAN_DEFINE_TEST(test_approx_edit_dist_pex_non_hierarchical_aho_corasick) {
    Test_Approx_EditDist< Pex<NonHierarchical,AhoCorasick > >();
}


SEQAN_DEFINE_TEST(test_approx_edit_dist_pex_non_hierarchical_multi_bfam) {
    Test_Approx_EditDist< Pex<NonHierarchical,MultiBFAM<> > >();
}


SEQAN_DEFINE_TEST(test_regression_rmbench) {
    // The data to test with:  Needle, haystack and score limit.
    const char *kCharNeedle = "CCATATGCTTGTGTCGCGGGTTTATTTGCATTCGACCCAGTTGACTCGGAAGTCGAAATGTTCCTGCCCCGTTTCTGCGTTCCGTGCAGTTGCGCGGTCTGGTTGGGCGGGTCCCCCCCTGA";
    const char *kCharHaystack = "ATTCCATATGCTTGTGTCGCGGGTTTATTTGCATTCGACCCAGTTGACTCGGAAGTCGAAATGTTCCTGCCCCGTTTCTGCGTTCCGTGCAGTTGCGCGGTCTGGTTGGGCGGGTCCCCCGCCTACGGATGCACGTTCTCCCGGGCTCGTAAATCC";
    const int kScoreLimit = -100000;

    // Build actual DNA needle and haystack, the pattern and the
    // finder.
    DnaString needle(kCharNeedle);
    DnaString haystack(kCharHaystack);
    Finder<DnaString> finder(haystack);
    Pattern<DnaString, MyersUkkonen> pattern(needle);
    setScoreLimit(pattern, kScoreLimit);

    // The following invariant should always hold: The pattern score
    // should be smaller than the needle size.
    SEQAN_ASSERT_LEQ(pattern.score, pattern.needleSize);
    while (find(finder, pattern)) {
        SEQAN_ASSERT_LEQ(pattern.score, pattern.needleSize);
    }
    SEQAN_ASSERT_LEQ(pattern.score, pattern.needleSize);
}


SEQAN_BEGIN_TESTSUITE(test_find) {
    // Testing MyersUkkonen with large needle and manual score limit.
    SEQAN_CALL_TEST(test_regression_rmbench);

    // Call all tests.
    SEQAN_CALL_TEST(test_find_online_Simple);
    SEQAN_CALL_TEST(test_find_online_Horspool);
    SEQAN_CALL_TEST(test_find_online_ShiftAnd);
    SEQAN_CALL_TEST(test_find_online_ShiftOr);
    SEQAN_CALL_TEST(test_find_online_BndmAlgo);
    SEQAN_CALL_TEST(test_find_online_BFAM_Oracle);
    SEQAN_CALL_TEST(test_find_online_BFAM_Trie);
    SEQAN_CALL_TEST(test_find_online_wildcards);
    SEQAN_CALL_TEST(test_find_online_multi_AhoCorasick);
    SEQAN_CALL_TEST(test_find_online_multi_MultipleShiftAnd);
    SEQAN_CALL_TEST(test_find_online_multi_SetHorspool);
    SEQAN_CALL_TEST(test_find_online_multi_WuManber);
    SEQAN_CALL_TEST(test_find_online_multi_MultiBFAM_Oracle);
    SEQAN_CALL_TEST(test_find_online_multi_MultiBFAM_Trie);
    SEQAN_CALL_TEST(test_find_approx_prefix_edit_dist_dpsearch);
    SEQAN_CALL_TEST(test_approx_prefix_edit_dist_myers);
    SEQAN_CALL_TEST(test_approx_edit_dist_dp_search_simple_score);
    SEQAN_CALL_TEST(test_approx_edit_dist_myers);
    // Test DP Serach.
    SEQAN_CALL_TEST(test_approx_edit_dp_search_simple_score_legacy_case);
    // Test other approximate search algorithms.
    SEQAN_CALL_TEST(test_approx_edit_dist_abndm_algo);
    SEQAN_CALL_TEST(test_approx_edit_dist_pex_non_hierarchical);
    SEQAN_CALL_TEST(test_approx_edit_dist_pex_hierarchical);
    // Tests with different multifinder.
    SEQAN_CALL_TEST(test_approx_edit_dist_pex_non_hierarchical_aho_corasick);
    SEQAN_CALL_TEST(test_approx_edit_dist_pex_non_hierarchical_multi_bfam);


    // TODO(holtgrew): Enable checkpoints again.
    /*
    //      debug::verifyCheckpoints("projects/library/seqan/find/find_myers_ukkonen.h");

    debug::verifyCheckpoints("projects/library/seqan/find/find_wild_shiftand.h");
    debug::verifyCheckpoints("projects/library/seqan/find/find_horspool.h");
    debug::verifyCheckpoints("projects/library/seqan/find/find_base.h");
    debug::verifyCheckpoints("projects/library/seqan/find/find_shiftand.h");
    debug::verifyCheckpoints("projects/library/seqan/find/find_shiftor.h");
    debug::verifyCheckpoints("projects/library/seqan/find/find_bndm.h");
    debug::verifyCheckpoints("projects/library/seqan/find/find_bom.h");
    //debug::verifyCheckpoints("projects/library/seqan/find/find_quasar.h");
    debug::verifyCheckpoints("projects/library/seqan/find/find_ahocorasick.h");
    debug::verifyCheckpoints("projects/library/seqan/find/find_multiple_shiftand.h");
    debug::verifyCheckpoints("projects/library/seqan/find/find_set_horspool.h");
    debug::verifyCheckpoints("projects/library/seqan/find/find_wumanber.h");
    debug::verifyCheckpoints("projects/library/seqan/find/find_abndm.h");
    debug::verifyCheckpoints("projects/library/seqan/find/find_pex.h");
    */
}
SEQAN_END_TESTSUITE
