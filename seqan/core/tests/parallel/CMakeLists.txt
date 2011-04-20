cmake_minimum_required (VERSION 2.6)
project (core_tests_parallel)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_parallel test_parallel.cpp
                          test_parallel_atomic_misc.h
					      test_parallel_atomic_primitives.h
					      test_parallel_splitting.h)
