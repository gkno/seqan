cmake_minimum_required (VERSION 2.6)
project (extras_tests_stream)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_stream test_stream.cpp test_stream.h)
