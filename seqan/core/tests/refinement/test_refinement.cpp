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
// Author: Anne-Katrin Emde <emde@fu-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Tests for the SeqAn moduel "refinement".
// ==========================================================================

#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_VERBOSE

// Test path
#define TEST_PATH "projects/tests/refinement/"
#define LIB_PATH "core/include/seqan/refinement/"

//#define LIB_PATH "/home/bude2/emde/seqan/core/include/seqan/refinement/"
//#define TEST_PATH "/home/bude2/emde/seqan/projects/tests/refinement/"

// External / STL
#include <iostream>
#include <fstream>
#include <string>

// SeqAn Includes
#include <seqan/refinement.h>
#include <seqan/basic.h>

// Test files
#include "test_graph_impl_align.h"
#include "test_graph_match_refinement.h"

using namespace seqan;


//////////////////////////////////////////////////////////////////////////////

SEQAN_BEGIN_TESTSUITE(test_refinement)
{
    // Redirect std::cout
//     std::ofstream file(TEST_PATH "redirect.txt");
//     std::streambuf* strm_puffer = std::cout.rdbuf();
//     std::cout.rdbuf(file.rdbuf());

//      Test_AlignmentGraph();
     // Test_GraphMatchRefinement();// Test Match Refinement

    // Test AlignmentGraph.
//     SEQAN_CALL_TEST(AlignmentGraphFunctions);
  //   SEQAN_CALL_TEST(HeaviestCommonSubsequence);
    // SEQAN_CALL_TEST(OutEdgeIteratorAlignment);

    // Test Match Refinement.
    SEQAN_CALL_TEST(RefineMatchesSelfEdges);

    //SEQAN_CALL_TEST(GraphMatchRefine);
    SEQAN_CALL_TEST(RefineAlign);
    SEQAN_CALL_TEST(RefineInexactFragment);

    // Restore std::cout
//      std::cout.rdbuf(strm_puffer);

    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/refinement/graph_algorithm_refine_align.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/refinement/graph_algorithm_refine_aligngraph.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/refinement/graph_algorithm_refine_annotation.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/refinement/graph_algorithm_refine_exact.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/refinement/graph_algorithm_refine_exact_iterative.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/refinement/graph_algorithm_refine_fragment.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/refinement/graph_algorithm_refine_inexact.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/refinement/graph_algorithm_refine_scoring.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/refinement/graph_impl_align.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/refinement/graph_impl_align_adapt.h");

}
SEQAN_END_TESTSUITE
