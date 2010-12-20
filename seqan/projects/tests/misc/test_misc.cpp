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
// ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <ctime>
#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/basic.h>

#include <seqan/map.h>

#include <seqan/misc/edit_environment.h>
#include <seqan/misc/misc_base.h>
#include <seqan/misc/misc_cmdparser.h>
#include <seqan/misc/misc_dequeue.h>
#include <seqan/misc/misc_long_word.h>
#include <seqan/misc/misc_map.h>
#include <seqan/misc/misc_random.h>
#include <seqan/misc/misc_set.h>
#include <seqan/misc/misc_interval_tree.h>
#include <seqan/misc/priority_type_base.h>
#include <seqan/misc/priority_type_heap.h>

#include "test_misc_long_word.h"
#include "test_misc_interval_tree.h"

using namespace std;
using namespace seqan;


SEQAN_DEFINE_TEST(test_misc_random) {
	mtRandInit();

	for (unsigned int i=0; i<100; ++i)
	{
		cout << mtRand() << ", ";
	}
	cout << "\n\n";

	for (unsigned int i=0; i<100; ++i)
	{
		cout << geomRand<int>() << ", ";
	}

}


SEQAN_BEGIN_TESTSUITE(test_misc) {
    //     SEQAN_CALL_TEST(test_misc_random);

    SEQAN_CALL_TEST(test_misc_long_word_native_interface);
    SEQAN_CALL_TEST(test_misc_long_word_static_interface);
    SEQAN_CALL_TEST(test_misc_long_word_dynamic_interface);


    // Test IntervalTree class
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_IntervalTree__int);
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_IntervalTreeFromIterator__int);
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_NonFullLength__int);
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_AddInterval__int);
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_TreeStructure__int);
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_FindIntervalExcludeTouching__int);
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_FindNoInterval__int);
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_GraphMap__int_ComputeCenter_StoreIntervals);
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_FindIntervalsIntervals__int_ComputeCenter);
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_Random__int_RandomCenter_StorePointsOnly);



    // Verify checkpoints.
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/edit_environment.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/misc_base.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/misc_cmdparser.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/misc_dequeue.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/misc_map.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/misc_long_word.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/misc_parsing.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/misc_random.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/misc_set.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/priority_type_base.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/priority_type_heap.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/misc_interval_tree.h");
}
SEQAN_END_TESTSUITE

