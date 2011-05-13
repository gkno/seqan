// ==========================================================================
//                                  parse_lm
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_HOLTGREW_TESTS_PARSE_LM_TEST_PARSE_LM_H_
#define SANDBOX_HOLTGREW_TESTS_PARSE_LM_TEST_PARSE_LM_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/parse_lm.h>

SEQAN_DEFINE_TEST(test_parse_lm_local_match_constructor)
{
    using namespace seqan;

    typedef LocalMatch<unsigned, unsigned> TLocalMatch;

    // Default constructor.
    {
        TLocalMatch localMatch;
        unsigned const maxU = MaxValue<unsigned>::VALUE;

        SEQAN_ASSERT_EQ(maxU, localMatch.id);
        SEQAN_ASSERT_EQ(maxU, localMatch.subjectId);
        SEQAN_ASSERT_EQ(maxU, localMatch.subjectBeginPos);
        SEQAN_ASSERT_EQ(maxU, localMatch.subjectEndPos);
        SEQAN_ASSERT_EQ(maxU, localMatch.queryId);
        SEQAN_ASSERT_EQ(maxU, localMatch.queryBeginPos);
        SEQAN_ASSERT_EQ(maxU, localMatch.queryEndPos);
    }
    // Constructor with options.
    {
        TLocalMatch localMatch(0, 1, 2, 3, 4, 5, 6);

        SEQAN_ASSERT_EQ(0u, localMatch.id);
        SEQAN_ASSERT_EQ(1u, localMatch.subjectId);
        SEQAN_ASSERT_EQ(2u, localMatch.subjectBeginPos);
        SEQAN_ASSERT_EQ(3u, localMatch.subjectEndPos);
        SEQAN_ASSERT_EQ(4u, localMatch.queryId);
        SEQAN_ASSERT_EQ(5u, localMatch.queryBeginPos);
        SEQAN_ASSERT_EQ(6u, localMatch.queryEndPos);
    }
}

SEQAN_DEFINE_TEST(test_parse_lm_local_match_store_constructor)
{
    using namespace seqan;

    // Default constructor.
    LocalMatchStore<> store;
}

SEQAN_DEFINE_TEST(test_parse_lm_local_match_store_append_local_match)
{
    using namespace seqan;

    LocalMatchStore<> store;

    // Append with sequence names.
    {
        appendLocalMatch(store, "seq0", 1, 2, "seq1", 3, 4);

        SEQAN_ASSERT_EQ(length(store.sequenceNameStore), 2u);
        SEQAN_ASSERT_EQ(store.sequenceNameStore[0], CharString("seq0"));
        SEQAN_ASSERT_EQ(store.sequenceNameStore[1], CharString("seq1"));

        SEQAN_ASSERT_EQ(length(store.matchStore), 1u);
        SEQAN_ASSERT_EQ(back(store.matchStore).id, 0u);
        SEQAN_ASSERT_EQ(back(store.matchStore).subjectId, 0u);
        SEQAN_ASSERT_EQ(back(store.matchStore).subjectBeginPos, 1u);
        SEQAN_ASSERT_EQ(back(store.matchStore).subjectEndPos, 2u);
        SEQAN_ASSERT_EQ(back(store.matchStore).queryId, 1u);
        SEQAN_ASSERT_EQ(back(store.matchStore).queryBeginPos, 3u);
        SEQAN_ASSERT_EQ(back(store.matchStore).queryEndPos, 4u);
    }

    // Append with sequence ids.
    {
        appendLocalMatch(store, 0, 5, 6, 1, 7, 8);

        SEQAN_ASSERT_EQ(length(store.sequenceNameStore), 2u);
        SEQAN_ASSERT_EQ(store.sequenceNameStore[0], CharString("seq0"));
        SEQAN_ASSERT_EQ(store.sequenceNameStore[1], CharString("seq1"));

        SEQAN_ASSERT_EQ(length(store.matchStore), 2u);
        SEQAN_ASSERT_EQ(back(store.matchStore).id, 1u);
        SEQAN_ASSERT_EQ(back(store.matchStore).subjectId, 0u);
        SEQAN_ASSERT_EQ(back(store.matchStore).subjectBeginPos, 5u);
        SEQAN_ASSERT_EQ(back(store.matchStore).subjectEndPos, 6u);
        SEQAN_ASSERT_EQ(back(store.matchStore).queryId, 1u);
        SEQAN_ASSERT_EQ(back(store.matchStore).queryBeginPos, 7u);
        SEQAN_ASSERT_EQ(back(store.matchStore).queryEndPos, 8u);
    }
}

#endif  // SANDBOX_HOLTGREW_TESTS_PARSE_LM_TEST_PARSE_LM_H_
