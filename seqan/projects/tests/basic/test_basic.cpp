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

#include <iostream>
#include <cstdio>
#include <vector>

#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "test_basic_alphabet.h"
#include "test_basic_aggregates.h"
#include "test_basic_allocator.h"
#include "test_basic_alphabet.h"
#include "test_basic_common.h"
#include "test_basic_construct_destruct.h"
#include "test_basic_holder.h"
#include "test_basic_iterator.h"
#include "test_basic_math.h"
#include "test_basic_metaprogramming.h"
#include "test_basic_proxy.h"
#include "test_basic_tag.h"
#include "test_basic_transport.h"
#include "test_basic_type.h"

SEQAN_BEGIN_TESTSUITE(test_basic)
{
    // =======================================================================
    // Tests for Miscalleneous Code
    // =======================================================================

    SEQAN_CALL_TEST(test_basic_metaprogramming_true);
    SEQAN_CALL_TEST(test_basic_metaprogramming_false);
    SEQAN_CALL_TEST(test_basic_metaprogramming_eval);
    SEQAN_CALL_TEST(test_basic_metaprogramming_or);
    SEQAN_CALL_TEST(test_basic_metaprogramming_and);
    SEQAN_CALL_TEST(test_basic_metaprogramming_if);
    SEQAN_CALL_TEST(test_basic_metaprogramming_is_same_type);
    SEQAN_CALL_TEST(test_basic_metaprogramming_switch);
    SEQAN_CALL_TEST(test_basic_metaprogramming_loop);
    SEQAN_CALL_TEST(test_basic_metaprogramming_loop_reverse);
    SEQAN_CALL_TEST(test_basic_metaprogramming_log2);
    SEQAN_CALL_TEST(test_basic_metaprogramming_log2_floor);
    SEQAN_CALL_TEST(test_basic_metaprogramming_power);
    SEQAN_CALL_TEST(test_basic_metaprogramming_make_unsigned);
    SEQAN_CALL_TEST(test_basic_metaprogramming_make_signed);
    SEQAN_CALL_TEST(test_basic_metaprogramming_make_remove_const);
    SEQAN_CALL_TEST(test_basic_metaprogramming_make_copy_const);
    SEQAN_CALL_TEST(test_basic_metaprogramming_make_is_const);
    SEQAN_CALL_TEST(test_basic_metaprogramming_class_identifier);
    SEQAN_CALL_TEST(test_basic_metaprogramming_memset);

    SEQAN_CALL_TEST(test_basic_math_int_pow);
    SEQAN_CALL_TEST(test_basic_math_log2);
    SEQAN_CALL_TEST(test_basic_math_min);
    SEQAN_CALL_TEST(test_basic_math_max);
    SEQAN_CALL_TEST(test_basic_math_abs);

    SEQAN_CALL_TEST(test_basic_tag_tag_struct);
    SEQAN_CALL_TEST(test_basic_tag_tag_basic_tags);
    SEQAN_CALL_TEST(test_basic_tag_move);
    SEQAN_CALL_TEST(test_basic_tag_misc_tags1);
    SEQAN_CALL_TEST(test_basic_tag_misc_tags2);
    SEQAN_CALL_TEST(test_basic_tag_tag_list_selector);
    SEQAN_CALL_TEST(test_basic_tag_misc_tags3);
    SEQAN_CALL_TEST(test_basic_tag_misc_tags4);

    SEQAN_CALL_TEST(seqan_basic_type_metafunction_value);
    SEQAN_CALL_TEST(seqan_basic_type_metafunction_get_value);
    SEQAN_CALL_TEST(seqan_basic_type_metafunction_reference);
    SEQAN_CALL_TEST(seqan_basic_type_metafunction_size);
    SEQAN_CALL_TEST(seqan_basic_type_metafunction_difference);
    SEQAN_CALL_TEST(seqan_basic_type_metafunction_position);
    SEQAN_CALL_TEST(seqan_basic_type_metafunction_host);
    SEQAN_CALL_TEST(seqan_basic_type_metafunction_spec);
    SEQAN_CALL_TEST(seqan_basic_type_metafunction_deepest_spec);
    SEQAN_CALL_TEST(seqan_basic_type_metafunction_cargo);
    SEQAN_CALL_TEST(seqan_basic_type_metafunction_vertex_descriptor);
    SEQAN_CALL_TEST(seqan_basic_type_metafunction_id);
    SEQAN_CALL_TEST(seqan_basic_type_metafunction_key);
    SEQAN_CALL_TEST(seqan_basic_type_metafunction_object);
    SEQAN_CALL_TEST(seqan_basic_type_metafunction_source);
    SEQAN_CALL_TEST(seqan_basic_type_metafunction_parameter);
    SEQAN_CALL_TEST(seqan_basic_type_to_parameter);
    SEQAN_CALL_TEST(seqan_basic_type_metafunction_const_parameter);
    SEQAN_CALL_TEST(seqan_basic_type_metafunction_length);
    SEQAN_CALL_TEST(seqan_basic_type_metafunction_is_integral);

    // =======================================================================
    // Tests for Aggregates
    // =======================================================================

    // -----------------------------------------------------------------------
    // Tests for Pairs
    // -----------------------------------------------------------------------
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_metafunctions);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_constructors);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_assign);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_set);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_move);
    // SEQAN_CALL_TEST(test_basic_aggregates_pair_base_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_get_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_assign_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_set_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_move_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_comparison_same_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_comparison_different_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_stream_output);

    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_metafunctions);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_constructors);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_assign);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_set);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_move);
    // SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_get_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_assign_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_set_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_move_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_comparison_same_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_comparison_different_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_stream_output);

    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_compressed_metafunctions);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_compressed_constructors);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_compressed_assign);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_compressed_set);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_compressed_move);
    // SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_compressed_value);  // TODO(holtgrew): Need proxy for this.
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_compressed_get_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_compressed_assign_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_compressed_set_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_compressed_move_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_compressed_comparison_same_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_compressed_comparison_different_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_compressed_stream_output);

    // -----------------------------------------------------------------------
    // Tests for Triples
    // -----------------------------------------------------------------------
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_metafunctions);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_constructors);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_assign);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_set);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_move);
    // SEQAN_CALL_TEST(test_basic_aggregates_triple_base_value);  // TODO(holtgrew): Need proxy for this.
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_get_value);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_assign_value);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_set_value);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_move_value);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_comparison_same_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_comparison_different_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_stream_output);

    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_metafunctions);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_constructors);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_assign);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_set);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_move);
    // SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_value);  // TODO(holtgrew): Need proxy for this.
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_get_value);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_assign_value);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_set_value);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_move_value);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_comparison_same_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_comparison_different_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_stream_output);

    // -----------------------------------------------------------------------
    // Tests for Tuples
    // -----------------------------------------------------------------------
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_metafunctions);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_constructors);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_assign);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_set);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_move);
    // SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_value);  // TODO(holtgrew): Need proxy for this.
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_get_value);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_assign_value);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_set_value);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_move_value);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_shift_left);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_shift_right);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_clear);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_length);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_comparison_same_spec);
    // SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_comparison_different_spec);  // TODO(holtgrew): Could be added for completeness case, not supported right now.
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_stream_output);

    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_compressed_metafunctions);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_compressed_constructors);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_compressed_assign);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_compressed_set);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_compressed_move);
    // SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_compressed_value);  // TODO(holtgrew): Need proxy for this.
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_compressed_get_value);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_compressed_assign_value);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_compressed_set_value);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_compressed_move_value);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_compressed_shift_left);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_compressed_shift_right);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_compressed_clear);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_compressed_length);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_compressed_comparison_same_spec);
    // SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_compressed_comparison_different_spec);  // TODO(holtgrew): Could be added for completeness case, not supported right now.
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_compressed_stream_output);

    // -----------------------------------------------------------------------
    // Tests for Allocators
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_allocator_simple);
    SEQAN_CALL_TEST(test_basic_allocator_pool);
    SEQAN_CALL_TEST(test_basic_allocator_multi_pool);
    SEQAN_CALL_TEST(test_basic_allocator_chunk_pool);

    // -----------------------------------------------------------------------
    // Tests for Construction / Destruction
    // -----------------------------------------------------------------------

    // TODO(holtgrew): Test on non-pointers.
    SEQAN_CALL_TEST(test_basic_construct_destruct_construct_value_pointer);
    SEQAN_CALL_TEST(test_basic_construct_destruct_destruct_value_pointer);
    SEQAN_CALL_TEST(test_basic_construct_destruct_array_construct_pointer);
    SEQAN_CALL_TEST(test_basic_construct_destruct_array_construct_copy_pointer);
    SEQAN_CALL_TEST(test_basic_construct_destruct_array_construct_move_pointer);
    SEQAN_CALL_TEST(test_basic_construct_destruct_array_destruct_pointer);
    SEQAN_CALL_TEST(test_basic_construct_destruct_array_fill_pointer);
    SEQAN_CALL_TEST(test_basic_construct_destruct_array_copy_forward_pointer);
    SEQAN_CALL_TEST(test_basic_construct_destruct_array_copy_backward_pointer);
    SEQAN_CALL_TEST(test_basic_construct_destruct_array_copy_pointer);
    SEQAN_CALL_TEST(test_basic_construct_destruct_array_move_forward_pointer);
    SEQAN_CALL_TEST(test_basic_construct_destruct_array_move_backward_pointer);
    SEQAN_CALL_TEST(test_basic_construct_destruct_array_move_pointer);
    SEQAN_CALL_TEST(test_basic_construct_destruct_array_clear_space_pointer);

    // -----------------------------------------------------------------------
    // Tests for Holders
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_holder);

    SEQAN_CALL_TEST(test_basic_alphabet_interface);
    // SEQAN_CALL_TEST(test_basic_conversions);
    SEQAN_CALL_TEST(test_basic_alphabet_extreme_values);
    // SEQAN_CALL_TEST(test_basic_simple_types);
    SEQAN_CALL_TEST(test_basic_array_functions);
    SEQAN_CALL_TEST(test_basic_suprema_infima);
    SEQAN_CALL_TEST(test_basic_common_definition);
    SEQAN_CALL_TEST(test_basic_common_type);
    SEQAN_CALL_TEST(test_basic_common_iterator_adapt_std);
    SEQAN_CALL_TEST(test_basic_iterator_basic);
    SEQAN_CALL_TEST(test_basic_iterator_adaptor);
    SEQAN_CALL_TEST(test_basic_iterator_position);
    SEQAN_CALL_TEST(test_basic_proxy_iterator);
    SEQAN_CALL_TEST(test_basic_transport);

}
SEQAN_END_TESTSUITE

