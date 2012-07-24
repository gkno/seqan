cmake_minimum_required (VERSION 2.6)
project (sandbox_ez_seq_io_tests_seq_io)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_seq_io
                          test_seq_io.cpp
                          test_sequence_stream.h
                          test_fai_index.h
                          test_genomic_region.h)
