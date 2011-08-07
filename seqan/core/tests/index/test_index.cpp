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

#include <iostream>
#include <fstream>
#include <functional>
#include <typeinfo>

#define SEQAN_DEBUG 
#define SEQAN_TEST 
#define SEQAN_ENABLE_CHECKPOINTS 0

#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/find.h>
#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/pipe.h>
#include <seqan/misc/misc_random.h>

#include "test_index_helpers.h"
#include "test_crosscompare.h"
#include "test_index_creation.h"
#include "test_qgram_index.h"
#include "test_sa_bwtwalk.h"
#include "test_shapes.h"
#include "test_stree_iterators.h"

using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

SEQAN_BEGIN_TESTSUITE(test_index)
{
	// test_shapes.h
	SEQAN_CALL_TEST(testShapes);
    
	// test_qgram_index.h
	SEQAN_CALL_TEST(testGappedShapes);
	SEQAN_CALL_TEST(testUngappedShapes);
	SEQAN_CALL_TEST(testUngappedQGramIndex);
	SEQAN_CALL_TEST(testQGramFind);
    
	// test_stree_iterators.h
	SEQAN_CALL_TEST(testSTreeIterators_Wotd);
	SEQAN_CALL_TEST(testSTreeIterators_WotdOriginal);
	SEQAN_CALL_TEST(testSTreeIterators_Esa);
	SEQAN_CALL_TEST(testFind_Esa_Mlr);
	SEQAN_CALL_TEST(testCompareIndices_Esa_Wotd);
	SEQAN_CALL_TEST(testMultiIndex);
	SEQAN_CALL_TEST(testMUMs);
	SEQAN_CALL_TEST(testMaxRepeats);
	SEQAN_CALL_TEST(testSuperMaxRepeats);
	SEQAN_CALL_TEST(testSuperMaxRepeatsFast);

	// test_crosscompare.h
	SEQAN_CALL_TEST(testIndexCrossCompare);

	// test_index_creation.h
	SEQAN_CALL_TEST(testIndexCreation);

	// test_sa_bwtwalk.h
	SEQAN_CALL_TEST(testBWTWalk);

//	testBuild();

    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/find_index_esa.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/find_index_qgram.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/find_index.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/find_quasar.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/find_swift.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_base.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_bwt.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_childtab.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_dfi.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_esa_algs_multi.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_esa_algs.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_esa_base.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_esa_drawing.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_esa_stree.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_lcp_tree.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_lcp.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_pizzachili_find.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_pizzachili_string.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_pizzachili.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_qgram_openaddressing.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_qgram.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_sa_btree.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_sa_bwtwalk.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_sa_lss.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_sa_mm.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_sa_qsort.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_shawarma.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_shims.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_skew3.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_skew7_multi.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_skew7.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/index_wotd.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/pipe_merger3.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/pipe_merger7.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/pizzachili_api.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/pump_extender3.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/pump_extender7.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/pump_lcp_core.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/pump_separator7.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/radix.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/repeat_base.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/shape_base.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/shape_gapped.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/shape_onegapped.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/shape_predefined.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/index/shape_threshold.h");
}
SEQAN_END_TESTSUITE
