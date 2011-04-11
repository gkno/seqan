cmake_minimum_required (VERSION 2.6)
project (core_tests_basic)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_basic test_basic.cpp
                                     test_basic_aggregates.h
                                     test_basic_allocator.h
                                     test_basic_alphabet.h
                                     test_basic_construct_destruct.h
                                     test_basic_holder.h
                                     test_basic_iterator.h
                                     test_basic_math.h
                                     test_basic_metaprogramming.h
                                     test_basic_proxy.h
                                     test_basic_tag.h
                                     test_basic_transport.h
                                     test_basic_type.h)
