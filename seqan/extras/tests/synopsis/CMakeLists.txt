cmake_minimum_required (VERSION 2.6)
project (extras_tests_synopsis)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_synopsis test_synopsis.cpp test_synopsis_counter_buckets.h test_synopsis_hot_list.h test_synopsis_histogram_set.h)
