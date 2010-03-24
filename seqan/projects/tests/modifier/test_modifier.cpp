#include <seqan/file.h>
#include <seqan/modifier.h>

#include "helpers.h"
#include "test_modifier_alphabet.h"
#include "test_modifier_view.h"
#include "test_modifier_functors.h"

using namespace std;
using namespace seqan;


SEQAN_DEFINE_TEST(test_modifier_reverse_string) {
		String<char> origin = "A man, a plan, a canal-Panama";

	//____________________________________________________________________________
	// Test1 - reverse string

		typedef ModifiedString< String<char>, ModReverse> TModString;

		TModString reverse(origin);

		cout << "*** Test1: Reverse String ***" << endl;
		cout << "origin:  " << origin << endl;
		cout << "reverse: " << reverse << endl << endl;

	//____________________________________________________________________________
	// Test2 - complement string

		DnaString dna = "attacgg";

		cout << "*** Test2: DNA symmetry ***" << endl;
		cout << "origin:             " << dna << endl;
		cout << "reverse:            " << DnaStringReverse(dna) << endl;
		cout << "complement:         " << DnaStringComplement(dna) << endl;
		cout << "reverse complement: " << DnaStringReverseComplement(dna) << endl << endl;

	//____________________________________________________________________________
	// Test3 - in-place conversions

		DnaString dna2 = dna;
		cout << "*** Test3: in-place conversions ***" << endl;
		cout << "origin:             " << dna2 << endl;

		reverseInPlace(dna2);
		cout << "reverse:            " << dna2 << endl;

		dna2 = dna; 
		complementInPlace(dna2);
		cout << "complement:         " << dna2 << endl;

		dna2 = dna; 
		reverseComplementInPlace(dna2);
		cout << "reverse complement: " << dna2 << endl << endl;

}


SEQAN_BEGIN_TESTSUITE(test_modifier) {
    // Tests for modifier_alphabet.h and modifier_alphabet_expansion.h.
    SEQAN_CALL_TEST(test_modifier_alphabet_size_metafunctions);
    SEQAN_CALL_TEST(test_modifier_alphabet_convert);
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
    // TODO(holtgrew): Write me!

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

    // Older tests...
    // TODO(holtgrew): Remove.
    SEQAN_CALL_TEST(test_modifier_reverse_string);

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
