cmake_minimum_required (VERSION 2.6)
project (extras_tests_ext_lh3)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_ext_lh3 test_ext_lh3.cpp test_stream_bgzf.h)
