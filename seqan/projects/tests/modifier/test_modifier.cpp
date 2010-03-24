#include <iostream>
#include <fstream>
#include <functional>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>
#include <time.h>

#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/file.h>
#include <seqan/modifier.h>

#include "helpers.h"
#include "test_modifier_alphabet.h"
#include "test_modifier_view.h"

using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(test_modifier_view_iteratorXX) {
		String<char> origin = "Vjku ku qwt qtkikpcn uvtkpi";

	//____________________________________________________________________________
	// Test1 - no modification (default)
	{

		typedef ModifiedIterator< Iterator<String<char>, Rooted>::Type > TModIterDefault;

		TModIterDefault it, itEnd(end(origin));

		setContainer(it, container(itEnd));		// test setContainer

		cout << "*** Iterator Test: Caesar chiffre ***" << endl;
		cout << "chiffre:  ";
		
		it = begin(origin);
		while (it != itEnd) {
			cout << *it;
			++it;
		}
		cout << endl;

	}
	//____________________________________________________________________________
	// Test2 - Caesar chiffre
	{

		typedef CaesarChiffre<char> TEncode;
		typedef ModifiedIterator< Iterator<String<char>, Rooted>::Type, ModView<TEncode> > TModIterCaesar;

		TEncode encode(-2);
		TModIterCaesar it(encode), itEnd(end(origin));

		setContainer(it, container(itEnd));		// test setContainer

		cout << "original: ";

		it = begin(origin);
		while (it != itEnd) {
			cout << *it;
			it = it + 1;
		}
		cout << endl << endl;

	}
}


SEQAN_DEFINE_TEST(test_modifier_view_string) {
		String<char> origin = "This is our original string";

	//____________________________________________________________________________
	// Test1 - no modification (default)

		typedef ModifiedString< String<char> > TModStringDefault;

		TModStringDefault nomod(origin);

		cout << "*** Test1/2: Caesar chiffre ***" << endl;
		cout << "origin:  " << nomod << endl;
		
	//____________________________________________________________________________
	// Test2 - Caesar chiffre

		typedef CaesarChiffre<char> TEncode;
		typedef ModifiedString< String<char>, ModView<TEncode> > TModStringCaesar;

		TEncode encode(2);
		TModStringCaesar chiffre(origin);
		assignModViewFunctor(chiffre, encode);

		cout << "chiffre: " << chiffre << endl << endl;

	//____________________________________________________________________________
	// Test3 - upcase/lowcase

		typedef ModifiedString< String<char>, ModView< FunctorUpcase<char> > > TModStringUp;
		typedef ModifiedString< String<char>, ModView< FunctorLowcase<char> > > TModStringLow;

		TModStringUp	up(origin);
		TModStringLow	low(origin);

		cout << "*** Test3: upcase/lowcase ***" << endl;
		cout << "upcase:   " << up << endl;
		cout << "lowcase: " << low << endl << endl;

	//____________________________________________________________________________
	// Test4 - alphabet conversion

		String<char> originDNA = "acgtnACGTN";
		typedef ModifiedString< String<char>, ModView< FunctorConvert<char, Dna5> > > TModStringDNA;

		TModStringDNA	dna(originDNA);

		cout << "*** Test4: alphabet conversion ***" << endl;
		cout << "origin: " << originDNA << endl;
		cout << "Dna5:   " << dna << endl << endl;

	//____________________________________________________________________________
	// Test5 - nested modifiers

		typedef CaesarChiffre<char> TEncode;
		typedef ModifiedString< 
		//			ModifiedString< 
						ModifiedString< 
							String<char>, 
							ModView<TEncode> 
						>, 
						ModView<TEncode> 
		//			>, 
		//			ModView<TEncode>
				> TModStringNested;

		TModStringNested nested(origin);

		TEncode enc2(2), enc3(-2), enc_5(-5);

		assignModViewFunctor(nested, enc2);
		assignModViewFunctor(host(nested), enc3);
//		assignModViewFunctor(host(host(nested)), enc_5);
		
		cout << (int) (begin(nested)).data_cargo.func._delta << "  ";
		cout << (int) host((begin(nested))).data_cargo.func._delta << "  ";
//		cout << (int) host(host(nested)).data_cargo.func.delta << "  ";
		
		cout << "*** Test5: nested modifiers ***" << endl;
		cout << "nested: " << nested << endl << endl;

}


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

//____________________________________________________________________________

SEQAN_DEFINE_TEST(test_modifier_alphabet_modifier) {
	typedef ModifiedAlphabet<Dna, ModExpand<'-'> > TDnaGap;
	typedef String<TDnaGap> TString;

	TString str = "aCgT-AcGt";
	cout << str << endl;

	SEQAN_TASSERT(str == "aCgT-AcGt");
	SEQAN_TASSERT(str == "AcGt-aCgT");
}


SEQAN_BEGIN_TESTSUITE(test_modifier) {
    // Tests for modifier_alphabet.h.
    SEQAN_CALL_TEST(test_modifier_alphabet_size_metafunctions);
    SEQAN_CALL_TEST(test_modifier_alphabet_convert);
    SEQAN_CALL_TEST(test_modifier_alphabet_ord_value);
    SEQAN_CALL_TEST(test_modifier_alphabet_operator_eq);
    SEQAN_CALL_TEST(test_modifier_alphabet_operator_neq);
    SEQAN_CALL_TEST(test_modifier_alphabet_operator_lt);
    SEQAN_CALL_TEST(test_modifier_alphabet_operator_gt);
    SEQAN_CALL_TEST(test_modifier_alphabet_operator_leq);
    SEQAN_CALL_TEST(test_modifier_alphabet_operator_geq);

    // Tests for modifier_view.h.
    SEQAN_CALL_TEST(test_modifier_view_iterator_metafunctions);
    SEQAN_CALL_TEST(test_modifier_view_iterator);
    SEQAN_CALL_TEST(test_modifier_view_const_iterator);
    SEQAN_CALL_TEST(test_modifier_convert_in_place);

    // Older tests...
    SEQAN_CALL_TEST(test_modifier_view_string);
    SEQAN_CALL_TEST(test_modifier_view_iterator);
    SEQAN_CALL_TEST(test_modifier_reverse_string);
    SEQAN_CALL_TEST(test_modifier_alphabet_modifier);

    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/modifier/modifier_alphabet.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/modifier/modifier_alphabet_expansion.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/modifier/modifier_view.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/modifier/modifier_functors.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/modifier/modifier_generated_forwards.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/modifier/modifier_iterator.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/modifier/modifier_reverse.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/modifier/modifier_shortcuts.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/modifier/modifier_string.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/modifier/modifier_view.h");
}
SEQAN_END_TESTSUITE
