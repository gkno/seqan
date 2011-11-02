cmake_minimum_required (VERSION 2.6)
project (core_tests_index)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_index_creation
                          test_index_creation.cpp
                          test_index_creation.h)

seqan_add_test_executable(test_index_crosscompare_char
                          test_index_crosscompare_char.cpp
                          test_crosscompare.h)

seqan_add_test_executable(test_index_crosscompare_dna
                          test_index_crosscompare_dna.cpp
                          test_crosscompare.h)

seqan_add_test_executable(test_index_qgram
                          test_index_qgram.cpp
                          test_qgram_index.h)

seqan_add_test_executable(test_index_sa_bwtwalk
                          test_index_sa_bwtwalk.cpp
                          test_sa_bwtwalk.h)

seqan_add_test_executable(test_index_shapes
                          test_index_shapes.cpp
                          test_shapes.h)

seqan_add_test_executable(test_index_stree_iterators
                          test_index_stree_iterators.cpp
                          test_stree_iterators.h)
