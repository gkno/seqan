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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Tests for the stream module.
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/stream.h>

#include "test_stream_char_array.h"
#if SEQAN_HAS_ZLIB
#include "test_stream_gz_file.h"
#endif  // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
#include "test_stream_bz2_file.h"
#endif  // #if SEQAN_HAS_BZIP2
#include "test_stream_adapt_cstdio.h"
#include "test_stream_adapt_fstream.h"
#include "test_stream_tokenize.h"
#include "test_stream_lexical_cast.h"
#include "test_stream_record_reader_fastaq.h"

SEQAN_BEGIN_TESTSUITE(test_stream)
{
    // Tests for Char Array Stream.
    SEQAN_CALL_TEST(test_stream_char_array_metafunctions);
    SEQAN_CALL_TEST(test_stream_char_array_read_simple_usage);
    SEQAN_CALL_TEST(test_stream_char_array_read_complex_usage);
    SEQAN_CALL_TEST(test_stream_char_array_write_simple_usage);
    SEQAN_CALL_TEST(test_stream_char_array_write_complex_usage);
    SEQAN_CALL_TEST(test_stream_char_array_eof);
    SEQAN_CALL_TEST(test_stream_char_array_peek);
    SEQAN_CALL_TEST(test_stream_char_array_read_char);
    SEQAN_CALL_TEST(test_stream_char_array_read_block);
    SEQAN_CALL_TEST(test_stream_char_array_write_block);
    SEQAN_CALL_TEST(test_stream_char_array_streamPut);
    SEQAN_CALL_TEST(test_stream_char_array_write_char);
    SEQAN_CALL_TEST(test_stream_char_array_flush);
    SEQAN_CALL_TEST(test_stream_char_array_seek);
    SEQAN_CALL_TEST(test_stream_char_array_tell);

#if SEQAN_HAS_ZLIB  // Enable tests for Stream<GZFile> if available.
    // Tests for BZ2 File Stream.
    SEQAN_CALL_TEST(test_stream_gz_file_metafunctions);
    SEQAN_CALL_TEST(test_stream_gz_file_read_simple_usage);
    SEQAN_CALL_TEST(test_stream_gz_file_read_complex_usage);
    SEQAN_CALL_TEST(test_stream_gz_file_write_simple_usage);
    SEQAN_CALL_TEST(test_stream_gz_file_write_complex_usage);
    SEQAN_CALL_TEST(test_stream_gz_file_eof);
    SEQAN_CALL_TEST(test_stream_gz_file_peek);
    SEQAN_CALL_TEST(test_stream_gz_file_read_char);
    SEQAN_CALL_TEST(test_stream_gz_file_read_block);
    SEQAN_CALL_TEST(test_stream_gz_file_write_block);
    SEQAN_CALL_TEST(test_stream_gz_file_streamPut);
    SEQAN_CALL_TEST(test_stream_gz_file_write_char);
    SEQAN_CALL_TEST(test_stream_gz_file_flush);
    SEQAN_CALL_TEST(test_stream_gz_file_seek);
    SEQAN_CALL_TEST(test_stream_gz_file_tell);
#endif  // #if SEQAN_HAS_ZLIB

#if SEQAN_HAS_BZIP2  // Enable tests for Stream<BZ2File> if available.
    // Tests for BZ2 File Stream.
    SEQAN_CALL_TEST(test_stream_bz2_file_metafunctions);
    SEQAN_CALL_TEST(test_stream_bz2_file_read_simple_usage);
    SEQAN_CALL_TEST(test_stream_bz2_file_read_complex_usage);
    SEQAN_CALL_TEST(test_stream_bz2_file_write_simple_usage);
    SEQAN_CALL_TEST(test_stream_bz2_file_write_complex_usage);
    SEQAN_CALL_TEST(test_stream_bz2_file_eof);
    SEQAN_CALL_TEST(test_stream_bz2_file_read_char);
    SEQAN_CALL_TEST(test_stream_bz2_file_read_block);
    SEQAN_CALL_TEST(test_stream_bz2_file_write_block);
    SEQAN_CALL_TEST(test_stream_bz2_file_write_char);
    SEQAN_CALL_TEST(test_stream_bz2_file_streamPut);
    SEQAN_CALL_TEST(test_stream_bz2_file_flush);
#endif  // #if SEQAN_HAS_BZIP2

    // Tests for cstdio.
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_metafunctions);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_read_simple_usage);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_read_complex_usage);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_write_simple_usage);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_write_complex_usage);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_eof);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_peek);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_read_char);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_read_block);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_write_block);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_write_char);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_streamPut);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_flush);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_seek);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_tell);

    // Tests for std::fstream adaptions.
    SEQAN_CALL_TEST(test_stream_adapt_fstream_metafunctions);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_read_simple_usage);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_read_complex_usage);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_write_simple_usage);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_write_complex_usage);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_eof);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_peek);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_read_char);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_read_block);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_write_block);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_write_char);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_streamPut);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_flush);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_seek);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_tell);

    // Tests for std::ifstream adaptions.
    SEQAN_CALL_TEST(test_stream_adapt_ifstream_metafunctions);
    SEQAN_CALL_TEST(test_stream_adapt_ifstream_read_simple_usage);
    SEQAN_CALL_TEST(test_stream_adapt_ifstream_read_complex_usage);
    SEQAN_CALL_TEST(test_stream_adapt_ifstream_eof);
    SEQAN_CALL_TEST(test_stream_adapt_ifstream_peek);
    SEQAN_CALL_TEST(test_stream_adapt_ifstream_read_char);
    SEQAN_CALL_TEST(test_stream_adapt_ifstream_read_block);
    SEQAN_CALL_TEST(test_stream_adapt_ifstream_seek);
    SEQAN_CALL_TEST(test_stream_adapt_ifstream_tell);

    // Tests for std::ofstream adaptions.
    SEQAN_CALL_TEST(test_stream_adapt_ofstream_metafunctions);
    SEQAN_CALL_TEST(test_stream_adapt_ofstream_write_simple_usage);
    SEQAN_CALL_TEST(test_stream_adapt_ofstream_write_complex_usage);
    SEQAN_CALL_TEST(test_stream_adapt_ofstream_write_block);
    SEQAN_CALL_TEST(test_stream_adapt_ofstream_write_char);
    SEQAN_CALL_TEST(test_stream_adapt_ofstream_streamPut);
    SEQAN_CALL_TEST(test_stream_adapt_ofstream_flush);

    // TODO(h4nn3s) finished tests for tokenize.h
    SEQAN_CALL_TEST(test_stream_tokenizing_readUntil);
    SEQAN_CALL_TEST(test_stream_tokenizing_readNChars);
    SEQAN_CALL_TEST(test_stream_tokenizing_readIgnoring);
    SEQAN_CALL_TEST(test_stream_tokenizing_readLine);
    SEQAN_CALL_TEST(test_stream_tokenizing_skipUntil);
    SEQAN_CALL_TEST(test_stream_tokenizing_skipWhile);
    SEQAN_CALL_TEST(test_stream_tokenizing_skipLine);
    SEQAN_CALL_TEST(test_stream_tokenizing_skipUntilString);
    SEQAN_CALL_TEST(test_stream_tokenizing_skipUntilLineBeginsWithChar);
    SEQAN_CALL_TEST(test_stream_tokenizing_skipUntilLineBeginsWithStr);
    SEQAN_CALL_TEST(test_stream_tokenizing_skipUntilLineBeginsWithOneCharOfStr);
    

    // Tests for lexical_cast
    SEQAN_CALL_TEST(test_stream_lexical_cast_1_stdstring);
    SEQAN_CALL_TEST(test_stream_lexical_cast_1_chararray);
    SEQAN_CALL_TEST(test_stream_lexical_cast_1_seqanstring);
    SEQAN_CALL_TEST(test_stream_lexical_cast_2_stdstring);
    SEQAN_CALL_TEST(test_stream_lexical_cast_2_chararray);
    SEQAN_CALL_TEST(test_stream_lexical_cast_2_seqanstring);

    
    // TODO(holtgrew): Tests for record reader class hierarchy. It's enough to test on fstream as a representative of the stream concept and the and memory mapped string specialization. Each test should be run with one and two records, as representatives for one/many case. Note that the buffer sizes could/should be set small enough or data large enough such that the fist buffer actually runs full in the two records variant. We test the FASTA reader as a representative.

    
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_single_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_double_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_batch_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_single_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_double_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_batch_mmap);

    SEQAN_CALL_TEST(test_stream_record_reader_fastq_single_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_double_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_batch_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_single_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_double_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_batch_mmap);

    
/*    SEQAN_CALL_TEST(test_stream_record_reader_fasta_single_fstream_two_records);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_double_fstream_one_record);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_double_fstream_two_records);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_single_mapped_one_record);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_single_mapped_two_records);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_double_mapped_one_record);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_double_mapped_two_records);
    */
    
    // TODO(holtgrew): Tests for FASTA and FASTQ record reader. Here, we test the single and double reader on the FASTA (above) and FASTQ (here) variant, both the one and two records case but only data from fstream.

    /*
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_single_fstream_one_record);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_single_fstream_two_records);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_double_fstream_one_record);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_double_fstream_two_records);
    */

    // TODO(holtgrew): Tests for the FASTA and FASTQ document reading functions, once they are in place. We test single/double pass for FASTA/FASTQ with one/two records, memory mapped file only, anything else is forbiddingly slow.
    // TODO(holtgrew): Probably, mapped only is sufficient.
    /*
    SEQAN_CALL_TEST(test_stream_read_document_fasta_single_mapped_one_record);
    SEQAN_CALL_TEST(test_stream_read_document_fasta_single_mapped_two_records);
    SEQAN_CALL_TEST(test_stream_read_document_fasta_double_mapped_one_record);
    SEQAN_CALL_TEST(test_stream_read_document_fasta_double_mapped_two_records);
    SEQAN_CALL_TEST(test_stream_read_document_fastq_single_mapped_one_record);
    SEQAN_CALL_TEST(test_stream_read_document_fastq_single_mapped_two_records);
    SEQAN_CALL_TEST(test_stream_read_document_fastq_double_mapped_one_record);
    SEQAN_CALL_TEST(test_stream_read_document_fastq_double_mapped_two_records);
    */
}
SEQAN_END_TESTSUITE

