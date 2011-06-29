// ==========================================================================
//                                  ext_lh3
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

#include <seqan/basic.h>
#include <seqan/file.h>

#include <seqan/ext_lh3.h>

#include "test_stream_bgzf.h"


SEQAN_BEGIN_TESTSUITE(test_ext_lh3)
{
    // Test Stream<Bgzf>.
    SEQAN_CALL_TEST(test_stream_bgzf_metafunctions);
    SEQAN_CALL_TEST(test_stream_bgzf_read_simple_usage);
    SEQAN_CALL_TEST(test_stream_bgzf_read_complex_usage);
    SEQAN_CALL_TEST(test_stream_bgzf_write_simple_usage);
    SEQAN_CALL_TEST(test_stream_bgzf_write_complex_usage);
    SEQAN_CALL_TEST(test_stream_bgzf_eof);
    SEQAN_CALL_TEST(test_stream_bgzf_peek);
    SEQAN_CALL_TEST(test_stream_bgzf_read_char);
    SEQAN_CALL_TEST(test_stream_bgzf_read_block);
    SEQAN_CALL_TEST(test_stream_bgzf_write_block);
    SEQAN_CALL_TEST(test_stream_bgzf_streamPut);
    SEQAN_CALL_TEST(test_stream_bgzf_write_char);
    SEQAN_CALL_TEST(test_stream_bgzf_flush);
    SEQAN_CALL_TEST(test_stream_bgzf_seek);
    SEQAN_CALL_TEST(test_stream_bgzf_tell);
}
SEQAN_END_TESTSUITE
