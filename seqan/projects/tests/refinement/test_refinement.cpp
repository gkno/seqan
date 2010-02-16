/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de
  ===========================================================================
  Copyright (C) 2007-2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  ===========================================================================
  Author: Anne-Katrin Emde <emde@fu-berlin.de>
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ===========================================================================
  Tests for the SeqAn moduel "refinement".
  ===========================================================================*/

#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_VERBOSE

// Test path
#define TEST_PATH "projects/tests/refinement/"
#define LIB_PATH "projects/library/seqan/refinement/"

//#define LIB_PATH "/home/bude2/emde/seqan/projects/library/seqan/refinement/"
//#define TEST_PATH "/home/bude2/emde/seqan/projects/tests/refinement/"

// External / STL
#include <iostream>
#include <fstream>
#include <string>

// SeqAn Includes
#include <seqan/refinement.h>
#include <seqan/basic/basic_testing.h>

// Test files
#include "test_graph_impl_align.h"
#include "test_graph_match_refinement.h"
#include "test_graph_interval_tree.h"

using namespace seqan;


//////////////////////////////////////////////////////////////////////////////

SEQAN_BEGIN_TESTSUITE(test_refinement)
{
    // Redirect std::cout
//     std::ofstream file(TEST_PATH "redirect.txt");
//     std::streambuf* strm_puffer = std::cout.rdbuf();
//     std::cout.rdbuf(file.rdbuf());

//      Test_AlignmentGraph();
//      Test_GraphMatchRefinement();// Test Match Refinement

    // Test AlignmentGraph.
//     SEQAN_CALL_TEST(AlignmentGraphFunctions);
//     SEQAN_CALL_TEST(HeaviestCommonSubsequence);
//     SEQAN_CALL_TEST(OutEdgeIteratorAlignment);

    // Test Match Refinement.
    SEQAN_CALL_TEST(GraphMatchRefine);
    SEQAN_CALL_TEST(RefineAlign);

    // Test IntervalTree class.
    SEQAN_CALL_TEST(Graph_Interval_Tree__IntervalTreeTest_IntervalTree__int);
    SEQAN_CALL_TEST(Graph_Interval_Tree__IntervalTreeTest_IntervalTreeFromIterator__int);
    SEQAN_CALL_TEST(Graph_Interval_Tree__IntervalTreeTest_NonFullLength__int);
    SEQAN_CALL_TEST(Graph_Interval_Tree__IntervalTreeTest_AddInterval__int);
    SEQAN_CALL_TEST(Graph_Interval_Tree__IntervalTreeTest_TreeStructure__int);
    SEQAN_CALL_TEST(Graph_Interval_Tree__IntervalTreeTest_FindIntervalExcludeTouching__int);
    SEQAN_CALL_TEST(Graph_Interval_Tree__IntervalTreeTest_FindNoInterval__int);
    SEQAN_CALL_TEST(Graph_Interval_Tree__IntervalTreeTest_GraphMap__int_ComputeCenter_StoreIntervals);
    SEQAN_CALL_TEST(Graph_Interval_Tree__IntervalTreeTest_Random__int_RandomCenter_StorePointsOnly);

    // Restore std::cout
//      std::cout.rdbuf(strm_puffer);

    // TODO(holtgrew): Enable checkpoint system again.
//      debug::verifyCheckpoints(LIB_PATH "graph_impl_interval_types.h");
//      debug::verifyCheckpoints(LIB_PATH "graph_impl_interval_tree.h");
//      debug::verifyCheckpoints(LIB_PATH "graph_algorithm_refine_scoring.h");
//      debug::verifyCheckpoints(LIB_PATH "graph_algorithm_refine_fragment.h");
//      debug::verifyCheckpoints(LIB_PATH "graph_algorithm_refine_exact.h");
//  debug::verifyCheckpoints(LIB_PATH "graph_algorithm_refine_aligngraph.h");
//  debug::verifyCheckpoints(LIB_PATH "graph_algorithm_refine_align.h");
//  debug::verifyCheckpoints(LIB_PATH "graph_algorithm_refine_annotation.h");
//  debug::verifyCheckpoints(LIB_PATH "graph_algorithm_refine_inexact.h");
}
SEQAN_END_TESTSUITE
