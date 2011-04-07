cmake_minimum_required (VERSION 2.6)
project (core_tests_index)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_index test_index.cpp
				     test_crosscompare.h
				     test_index_creation.h
				     test_index_helpers.h
				     test_qgram_index.h
				     test_sa_bwtwalk.h
				     test_shapes.h
				     test_stree_iterators.h)
