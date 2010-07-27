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
  Author: Carsten Kemena <carsten.kemena@crg.es>
 ===========================================================================
  Tests for the banded chain alignment.  Contains code that is based on
  the test code by Carsten Kemena.
 ===========================================================================
*/

#ifndef TEST_SEEDS_TEST_ALIGN_CHAIN_BANDED_H_
#define TEST_SEEDS_TEST_ALIGN_CHAIN_BANDED_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.

// Test banded alignment algorithm around a chain.  Linear gap cost
// case.
SEQAN_DEFINE_TEST(test_align_chain_banded_align_linear)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    
    // Test on whole strings.
    {
        CharString query = "ACGTCCTCGTACACCGTCTTAA";
        CharString database = "TACGATCCACACCGCGTCT";
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(1, 2, 6, 8));
        appendValue(seedChain, TSeed(12, 10, 19, 19));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), query);
        assignSource(row(alignment, 1), database);

        //cout << "Score: " << bandedChainAlignment(seedChain1, 2, alignment1, scoreMatrix) << endl;
        int result = bandedChainAlignment(seedChain, 2, alignment, scoringScheme);
        SEQAN_ASSERT_EQ(result, 11);

        std::cout << alignment << std::endl;
        SEQAN_ASSERT_TRUE(row(alignment, 0) == "ACGTCCTCGTACACCGTCTTAA");
        SEQAN_ASSERT_TRUE(row(alignment, 1) == "TACGATC-C--ACACCG-CGTCT");
    }
    // Test on infixes.
    {
        CharString query = "ACGTCCTCGTACACCGTCTTAA";
        CharString database = "TACGATCCACACCGCGTCT";
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed > seedChain;
        appendValue(seedChain, TSeed(1, 2, 6, 8));
        appendValue(seedChain, TSeed(12, 10, 19, 19));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        // TODO(holtgrew): This infix assignment for alignments is a bit creepy, maybe one of the too many shortcuts?
        assignSource(row(alignment, 0), query, 1, length(query));
        assignSource(row(alignment, 1), database, 2, length(database));

        //cout << "Score: " << bandedChainAlignment(seedChain1, 2, alignment2, scoreMatrix) << endl;
        int result = bandedChainAlignment(seedChain, 2, alignment, scoringScheme);
        SEQAN_ASSERT_EQ(result, 11);

        //cout << alignment2 << endl;
        SEQAN_ASSERT_TRUE(row(alignment, 0) == "CGTCCTCGTACACCGTCTTAA" );
        SEQAN_ASSERT_TRUE(row(alignment, 1) == "CGATC-C--ACACCG-CGTCT");
    }
}


// Test banded alignment algorithm around a chain.  Linear gap cost
// case.
SEQAN_DEFINE_TEST(test_align_chain_banded_align_affine)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    
    // Test on whole strings.
    {
        CharString query = "ACGTCCTCGTACACCGTCTTAA";
        CharString database = "TACGATCCACACCGCGTCT";
        Score<int, Simple> scoringScheme(3, -2, -1, -3);
        
        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(1, 2, 6, 8));
        appendValue(seedChain, TSeed(12, 10, 19, 19));
        
        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), query);
        assignSource(row(alignment, 1), database);
        
        //cout << "Score: " << bandedChainAlignment(seedChain2, 2, alignment3, scoreMatrix2) << endl;
        int result = bandedChainAlignment(seedChain, 2, alignment, scoringScheme);
        SEQAN_ASSERT_EQ(result, 24);
        
        //cout << alignment3 << endl;
        SEQAN_ASSERT_TRUE(row(alignment, 0) == "ACG-TCCTCGTACAC--CGTCTTAA");
        SEQAN_ASSERT_TRUE(row(alignment, 1) == "TACGATCC----ACACCGCGTCT");
    }
    // Test on infixes.
    {
        CharString query = "ACGTCCTCGTACACCGTCTTAA";
        CharString database = "TACGATCCACACCGCGTCT";
        Score<int, Simple> scoringScheme(3, -2, -1, -3);
        
        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(1, 2, 6, 8));
        appendValue(seedChain, TSeed(12, 10, 19, 19));
        
        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), query, 1, length(query));
        assignSource(row(alignment, 1), database, 2, length(database));
        
        //cout << "Score: " << bandedChainAlignment(seedChain2, 2, alignment4, scoreMatrix2) << endl;
        int result = bandedChainAlignment(seedChain, 2, alignment, scoringScheme);
        SEQAN_ASSERT_EQ(result, 21);
        
        //cout << alignment4 << endl;
        SEQAN_ASSERT_TRUE(row(alignment, 0) == "CG-TCCTCGTACAC--CGTCTTAA");
        SEQAN_ASSERT_TRUE(row(alignment, 1) == "CGATCC----ACACCGCGTCT");
    }
}

#endif  // TEST_SEEDS_TEST_ALIGN_CHAIN_BANDED_H_
