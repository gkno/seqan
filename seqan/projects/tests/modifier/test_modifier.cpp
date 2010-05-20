#include <seqan/file.h>
#include <seqan/modifier.h>

#include "helpers.h"
#include "test_modifier_alphabet.h"
#include "test_modifier_view.h"
#include "test_modifier_functors.h"
#include "test_modifier_shortcuts.h"

using namespace std;
using namespace seqan;


SEQAN_BEGIN_TESTSUITE(test_modifier) {
    // Tests for modifier_alphabet.h and modifier_alphabet_expansion.h.
    SEQAN_CALL_TEST(test_modifier_alphabet_size_metafunctions);
    SEQAN_CALL_TEST(test_modifier_alphabet_convert);
	SEQAN_CALL_TEST(test_modifier_DnaQ);
    SEQAN_CALL_TEST(test_modifier_alphabet_enumerate);
    SEQAN_CALL_TEST(test_modifier_alphabet_ord_value);
    SEQAN_CALL_TEST(test_modifier_alphabet_operator_eq);
    SEQAN_CALL_TEST(test_modifier_alphabet_operator_neq);
    SEQAN_CALL_TEST(test_modifier_alphabet_operator_lt);
    SEQAN_CALL_TEST(test_modifier_alphabet_operator_gt);
    SEQAN_CALL_TEST(test_modifier_alphabet_operator_leq);
    SEQAN_CALL_TEST(test_modifier_alphabet_operator_geq);

    // Tests for modifier_functors.h.
    SEQAN_CALL_TEST(test_modifier_functors_functor_upcase);
    SEQAN_CALL_TEST(test_modifier_functors_functor_lowcase);
    SEQAN_CALL_TEST(test_modifier_functors_dna_complement);

    // Tests for modifier_iterator.h.
    // TODO(holtgrew): Write me!

    // Tests for modifier_reverse.h.
    // TODO(holtgrew): Write me!

    // Tests for modifier_shortcuts.h.
    SEQAN_CALL_TEST(test_modifer_shortcuts_dna_string_reverse);
    SEQAN_CALL_TEST(test_modifer_shortcuts_dna5_string_reverse);
    SEQAN_CALL_TEST(test_modifer_shortcuts_dna_string_complement);
    SEQAN_CALL_TEST(test_modifer_shortcuts_dna5_string_complement);
    SEQAN_CALL_TEST(test_modifer_shortcuts_dna_string_reverse_complement);
    SEQAN_CALL_TEST(test_modifer_shortcuts_dna5_string_reverse_complement);
    SEQAN_CALL_TEST(test_modifer_shortcuts_complement_in_place_string);
    SEQAN_CALL_TEST(test_modifer_shortcuts_complement_in_place_string_set);
    SEQAN_CALL_TEST(test_modifer_shortcuts_reverse_in_place_string);
    SEQAN_CALL_TEST(test_modifer_shortcuts_reverse_in_place_string_set);
    SEQAN_CALL_TEST(test_modifer_shortcuts_reverse_complement_in_place_string);
    SEQAN_CALL_TEST(test_modifer_shortcuts_reverse_complement_in_place_string_set);
    SEQAN_CALL_TEST(test_modifer_shortcuts_to_lower_in_place_string);
    SEQAN_CALL_TEST(test_modifer_shortcuts_to_lower_in_place_string_set);
    SEQAN_CALL_TEST(test_modifer_shortcuts_to_upper_in_place_string);
    SEQAN_CALL_TEST(test_modifer_shortcuts_to_upper_in_place_string_set);

    // Tests for modifier_string.h.
    // TODO(holtgrew): Write me!

    // Tests for modifier_view.h.
    SEQAN_CALL_TEST(test_modifier_view_iterator_metafunctions);
    SEQAN_CALL_TEST(test_modifier_view_iterator);
    SEQAN_CALL_TEST(test_modifier_view_const_iterator);
    SEQAN_CALL_TEST(test_modifier_convert_in_place);

    SEQAN_CALL_TEST(test_modifier_view_string_caesar_chiffre);
    SEQAN_CALL_TEST(test_modifier_view_string_upper_case);
    SEQAN_CALL_TEST(test_modifier_view_string_low_case);
    SEQAN_CALL_TEST(test_modifier_view_string_alphabet_conversion);
    SEQAN_CALL_TEST(test_modifier_view_string_alphabet_conversion);
    SEQAN_CALL_TEST(test_modifier_view_string_nested_modifier);

    // Verify check points for all headers in the module modifier.
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/modifier/modifier_alphabet.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/modifier/modifier_alphabet_expansion.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/modifier/modifier_functors.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/modifier/modifier_iterator.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/modifier/modifier_reverse.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/modifier/modifier_shortcuts.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/modifier/modifier_string.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/modifier/modifier_view.h");
}
SEQAN_END_TESTSUITE
